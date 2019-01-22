package org.broadinstitute.hellbender.tools.walkers.mutect.filtering;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.util.OverlapDetector;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import org.apache.commons.lang3.mutable.MutableDouble;
import org.apache.commons.lang3.mutable.MutableInt;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.walkers.contamination.ContaminationRecord;
import org.broadinstitute.hellbender.tools.walkers.contamination.MinorAlleleFractionRecord;
import org.broadinstitute.hellbender.tools.walkers.mutect.Mutect2Engine;
import org.broadinstitute.hellbender.tools.walkers.mutect.MutectStats;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.io.File;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Simple class to hold all information needed by Mutect2 filters, including: the M2FilteringArgumentCollection,
 * the set of normal samples, the samples' contamination, and the tumor sample CNV segments
 */
public class Mutect2FilteringInfo {
    private static final double FIRST_PASS_THRESHOLD = 0.5;
    private static final double EPSILON = 1.0e-10;

    private final M2FiltersArgumentCollection MTFAC;
    private final Set<String> normalSamples;
    private final Map<String, Double> contaminationBySample;
    private final Map<String, OverlapDetector<MinorAlleleFractionRecord>> tumorSegments;

    private double artifactProbabilityThreshold = FIRST_PASS_THRESHOLD;

    private double log10PriorOfSomaticSNV;
    private double log10PriorOfSomaticIndel;

    private double priorProbOfArtifactVersusVariant;

    private final MutableDouble realVariantCount = new MutableDouble(0);
    private final MutableDouble realSNVCount = new MutableDouble(0);
    private final MutableDouble realIndelCount = new MutableDouble(0);
    private final MutableDouble technicalArtifactCount = new MutableDouble(0);



    private List<Mutect2VariantFilter> filters;

    private MutableInt passingVariants = new MutableInt(0);
    private Map<String, MutableDouble> expectedFalsePositivesPerFilter = new HashMap<>();
    private Map<String, MutableDouble> expectedFalseNegativesPerFilter = new HashMap<>();
    private final MutableDouble expectedFalsePositives = new MutableDouble(0);
    private final MutableDouble expectedTruePositives = new MutableDouble(0);
    private final MutableDouble expectedFalseNegatives = new MutableDouble(0);

    private OptionalLong totalCallableSites = OptionalLong.empty();

    // for each PID, the positions with PGTs of filtered genotypes
    private final Map<String, ImmutablePair<Integer, Set<String>>> filteredPhasedCalls;

    // artifact posteriors from first pass
    List<Double> firstPassArtifactProbabilities = new ArrayList<>();

    public Mutect2FilteringInfo(M2FiltersArgumentCollection MTFAC, final VCFHeader vcfHeader) {
        this.MTFAC = MTFAC;
        normalSamples = vcfHeader.getMetaDataInInputOrder().stream()
                .filter(line -> line.getKey().equals(Mutect2Engine.NORMAL_SAMPLE_KEY_IN_VCF_HEADER))
                .map(VCFHeaderLine::getValue)
                .collect(Collectors.toSet());

        contaminationBySample = MTFAC.contaminationTables.stream()
                .map(file -> ContaminationRecord.readFromFile(file).get(0))
                .collect(Collectors.toMap(rec -> rec.getSample(), rec -> rec.getContamination()));

        for (final String sample : vcfHeader.getSampleNamesInOrder()) {
            if (!contaminationBySample.containsKey(sample)) {
                contaminationBySample.put(sample, MTFAC.contaminationEstimate);
            }
        }

        tumorSegments = MTFAC.tumorSegmentationTables.stream()
                .map(MinorAlleleFractionRecord::readFromFile)
                .collect(Collectors.toMap(ImmutablePair::getLeft, p -> OverlapDetector.create(p.getRight())));

        filteredPhasedCalls = new HashMap<>();

        log10PriorOfSomaticSNV = MTFAC.log10PriorProbOfSomaticSNV;
        log10PriorOfSomaticIndel = MTFAC.log10PriorProbOfSomaticIndel;

        priorProbOfArtifactVersusVariant = MTFAC.initialPriorOfArtifactVersusVariant;

        filters = new ArrayList<>();
        filters.add(new TumorEvidenceFilter());
        filters.add(new BaseQualityFilter());
        filters.add(new MappingQualityFilter());
        filters.add(new DuplicatedAltReadFilter());
        filters.add(new StrandArtifactFilter());
        filters.add(new ContaminationFilter());
        filters.add(new PanelOfNormalsFilter());
        filters.add(new NormalArtifactFilter());
        filters.add(new ReadOrientationFilter());
        filters.add(new NRatioFilter());
        filters.add(new StrictStrandBiasFilter());
        filters.add(new ReadPositionFilter());

        if (MTFAC.mitochondria) {
            filters.add(new LogOddsOverDepthFilter());
            filters.add(new ChimericOriginalAlignmentFilter());
        } else {
            filters.add(new ClusteredEventsFilter());
            filters.add(new MultiallelicFilter());
            filters.add(new FragmentLengthFilter());
            filters.add(new PolymeraseSlippageFilter());
            filters.add(new FilteredHaplotypeFilter());
            filters.add(new GermlineFilter());
        }

        filters.forEach(filter -> expectedFalsePositivesPerFilter.put(filter.filterName(), new MutableDouble(0)));
        filters.forEach(filter -> expectedFalseNegativesPerFilter.put(filter.filterName(), new MutableDouble(0)));
    }

    public void inputMutectStats(final File mutectStatsTable) {
        final Map<String, Double> stats = MutectStats.readFromFile(mutectStatsTable).stream()
                .collect(Collectors.toMap(MutectStats::getStatistic, MutectStats::getValue));

        totalCallableSites = OptionalLong.of(Math.round(stats.get(Mutect2Engine.CALLABLE_SITES_NAME)));
    }

    public M2FiltersArgumentCollection getMTFAC() {
        return MTFAC;
    }

    public Set<String> getNormalSamples() {
        return normalSamples;
    }

    public Map<String, Double> getContaminationBySample() {
        return contaminationBySample;
    }

    public Map<String, OverlapDetector<MinorAlleleFractionRecord>> getTumorSegments() {
        return tumorSegments;
    }

    public Map<String, ImmutablePair<Integer, Set<String>>> getFilteredPhasedCalls() {
        return filteredPhasedCalls;
    }

    public double getArtifactProbabilityThreshold() {
        return artifactProbabilityThreshold;
    }

    public double getLog10PriorOfSomaticVariant(final VariantContext vc) {
        return vc.isSNP() ? (MathUtils.LOG10_ONE_THIRD + log10PriorOfSomaticSNV) : log10PriorOfSomaticIndel;
    }

    public double getPriorProbOfArtifactVersusVariant() { return priorProbOfArtifactVersusVariant; }

    public void recordFilteredHaplotypes(final VariantContext vc) {
        final Map<String, Set<String>> phasedGTsForEachPhaseID = vc.getGenotypes().stream()
                .filter(gt -> !normalSamples.contains(gt.getSampleName()))
                .filter(Mutect2FilteringInfo::hasPhaseInfo)
                .collect(Collectors.groupingBy(g -> (String) g.getExtendedAttribute(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_ID_KEY, ""),
                        Collectors.mapping(g -> (String) g.getExtendedAttribute(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_GT_KEY, ""), Collectors.toSet())));

        for (final Map.Entry<String, Set<String>> pidAndPgts : phasedGTsForEachPhaseID.entrySet()) {
            filteredPhasedCalls.put(pidAndPgts.getKey(), new ImmutablePair<>(vc.getStart(), pidAndPgts.getValue()));
        }
    }

    public void addFirstPassArtifactProbability(final double prob) {
        firstPassArtifactProbabilities.add(prob);
    }

    public void addRealVariantCount(final double x, final boolean isSNV) {
        realVariantCount.add(x);
        (isSNV ? realSNVCount : realIndelCount).add(x);
    }

    public void addTechnicalArtifactCount(final double x) { technicalArtifactCount.add(x); }

    public void learnPriorProbOfArtifactVersusVariant() {
        priorProbOfArtifactVersusVariant = (technicalArtifactCount.getValue() + 1) / (realVariantCount.getValue() + technicalArtifactCount.getValue() + 2);
    }

    public void learnPriorProbOfVariant() {
        if (totalCallableSites.isPresent()) {
            log10PriorOfSomaticSNV = Math.log10(realSNVCount.getValue() / totalCallableSites.getAsLong());
            log10PriorOfSomaticIndel = Math.log10(realIndelCount.getValue() / totalCallableSites.getAsLong());
        }
    }

    public void adjustThreshold() {
        final double threshold;
        switch (MTFAC.thresholdStrategy) {
            case CONSTANT:
                threshold = MTFAC.posteriorThreshold;
                break;
            case FALSE_DISCOVERY_RATE:
                threshold = calculateThresholdBasedOnFalseDiscoveryRate(firstPassArtifactProbabilities, MTFAC.maxFalsePositiveRate);
                break;
            case OPTIMAL_F_SCORE:
                threshold = calculateThresholdBasedOnOptimalFScore(firstPassArtifactProbabilities, MTFAC.fScoreBeta);
                break;
            default:
                throw new GATKException.ShouldNeverReachHereException("Invalid threshold strategy type: " + MTFAC.thresholdStrategy + ".");
        }

        firstPassArtifactProbabilities.clear();

        artifactProbabilityThreshold = threshold;
    }

    /**
     *
     * Compute the filtering threshold that ensures that the false positive rate among the resulting pass variants
     * will not exceed the requested false positive rate
     *
     * @param posteriors A list of posterior probabilities, which gets sorted
     * @param requestedFPR We set the filtering threshold such that the FPR doesn't exceed this value
     * @return
     */
    @VisibleForTesting
    static double calculateThresholdBasedOnFalseDiscoveryRate(final List<Double> posteriors, final double requestedFPR){
        ParamUtils.isPositiveOrZero(requestedFPR, "requested FPR must be non-negative");
        final double thresholdForFilteringNone = 1.0;
        final double thresholdForFilteringAll = 0.0;

        Collections.sort(posteriors);

        final int numPassingVariants = posteriors.size();
        double cumulativeExpectedFPs = 0.0;

        for (int i = 0; i < numPassingVariants; i++){
            final double posterior = posteriors.get(i);

            // One can show that the cumulative error rate is monotonically increasing in i
            final double expectedFPR = (cumulativeExpectedFPs + posterior) / (i + 1);
            if (expectedFPR > requestedFPR){
                return i > 0 ? posteriors.get(i-1) : thresholdForFilteringAll;
            }

            cumulativeExpectedFPs += posterior;
        }

        // If the expected FP rate never exceeded the max tolerable value, then we can let everything pass
        return thresholdForFilteringNone;
    }

    /**
     * Compute the filtering threshold that maximizes the F_beta score
     *
     * @param posteriors A list of posterior probabilities, which gets sorted
     * @param beta relative weight of recall to precision
     */
    @VisibleForTesting
    static double calculateThresholdBasedOnOptimalFScore(final List<Double> posteriors, final double beta){
        ParamUtils.isPositiveOrZero(beta, "requested F-score beta must be non-negative");
        final double thresholdForFilteringNone = 1.0;
        final double thresholdForFilteringAll = 0.0;

        Collections.sort(posteriors);

        final double expectedTruePositives = posteriors.stream()
                .mapToDouble(prob -> 1 - prob).sum();


        // starting from filtering everything (threshold = 0) increase the threshold to maximize the F score
        final MutableDouble truePositives = new MutableDouble(0);
        final MutableDouble falsePositives = new MutableDouble(0);
        final MutableDouble falseNegatives = new MutableDouble(expectedTruePositives);
        int optimalIndexInclusive = -1; // include all indices up to and including this. -1 mean filter all.
        double optimalFScore = 0;   // if you exclude everything, recall is zero

        final int N = posteriors.size();

        for (int n = 0; n < N; n++){
            truePositives.add(1 - posteriors.get(n));
            falsePositives.add(posteriors.get(n));
            falseNegatives.subtract(1 - posteriors.get(n));
            final double F = (1+beta*beta)*truePositives.getValue() /
                    ((1+beta*beta)*truePositives.getValue() + beta*beta*falseNegatives.getValue() + falsePositives.getValue());
            if (F >= optimalFScore) {
                optimalIndexInclusive = n;
                optimalFScore = F;
            }
        }

        return optimalIndexInclusive == -1 ? 0 : (optimalIndexInclusive == N - 1 ? 1 : posteriors.get(optimalIndexInclusive));
    }

    public static boolean hasPhaseInfo(final Genotype genotype) {
        return genotype.hasExtendedAttribute(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_GT_KEY) && genotype.hasExtendedAttribute(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_ID_KEY);
    }

    public void accumulateDataForLearning(final VariantContext vc) {
        filters.forEach(f -> f.accumulateDataForLearning(vc, this));
    }

    public void accumulateRealAndArtifactCounts(final VariantContext vc) {
        double artifactProbability = 0;
        double technicalArtifactProbability = 0;

        for (final Mutect2VariantFilter filter : filters) {
            final double prob = filter.artifactProbability(vc, this);
            artifactProbability = Math.max(artifactProbability, prob);
            technicalArtifactProbability = filter.isTechnicalArtifact() ? Math.max(technicalArtifactProbability, prob) : technicalArtifactProbability;
        }

        final Pair<Double, Double> overallAndTechnicalArtifactProbs = overallAndTechnicalOnlyArtifactProbabilities(vc);

        addRealVariantCount(1 - overallAndTechnicalArtifactProbs.getLeft(), vc.isSNP());
        addTechnicalArtifactCount(overallAndTechnicalArtifactProbs.getRight());
    }

    // TODO: dupes from FilteRMutectCalls
    private Pair<Double, Double> overallAndTechnicalOnlyArtifactProbabilities(final VariantContext vc) {
        return overallAndTechnicalOnlyArtifactProbabilities(filters.stream()
                .collect(Collectors.toMap(f -> f, f -> f.calculateArtifactProbability(vc, this))));
    }

    private Pair<Double, Double> overallAndTechnicalOnlyArtifactProbabilities(final Map<Mutect2VariantFilter, Double> map) {
        final double technicalArtifactProbability = map.entrySet().stream().filter(entry -> entry.getKey().isTechnicalArtifact())
                .mapToDouble(entry -> entry.getValue()).max().orElseGet(() -> 0);
        final double nonTechnicalArtifactProbability = map.entrySet().stream().filter(entry -> !entry.getKey().isTechnicalArtifact())
                .mapToDouble(entry -> entry.getValue()).max().orElseGet(() -> 0);

        final double overallProbability = 1 - (1 - technicalArtifactProbability) * (1 - nonTechnicalArtifactProbability);

        return ImmutablePair.of(overallProbability, technicalArtifactProbability);
    }

    public void learnFilterParameters() {
        filters.forEach(f -> f.learnParameters());
    }

    public void writeFilteringStats(final File filteringStatsFile) {
        final int totalCalls = passingVariants.getValue();
        final List<FilterStats> filterStats = filters.stream().map(Mutect2VariantFilter::filterName)
                .map(filter -> {
                    final double falseNegatives = expectedFalseNegativesPerFilter.get(filter).getValue();
                    final double falsePositives = expectedFalsePositivesPerFilter.get(filter).getValue();
                    final double fdr = falsePositives / totalCalls;
                    final double totalTrueVariants = expectedTruePositives.getValue() + expectedFalseNegatives.getValue();

                    return new FilterStats(filter, falsePositives, fdr, falseNegatives, falseNegatives / totalTrueVariants);
                })
                .filter(stats -> stats.getFalsePositiveCount() > 0 || stats.getFalseNegativeCount() > 0)
                .collect(Collectors.toList());


        FilterStats.writeM2FilterSummary(filterStats, filteringStatsFile, getArtifactProbabilityThreshold(), totalCalls,
                expectedTruePositives.getValue(), expectedFalsePositives.getValue(), expectedFalseNegatives.getValue());
    }

    public VariantContext makeFilteredVariant(final VariantContext vc) {
        final VariantContextBuilder vcb = new VariantContextBuilder(vc);
        vcb.filters(new HashSet<>());

        final Map<Mutect2VariantFilter, Double> artifactProbabilities = filters.stream()
                .collect(Collectors.toMap(f -> f, f -> f.artifactProbability(vc, this)));

        final double overallArtifactProbability = overallAndTechnicalOnlyArtifactProbabilities(artifactProbabilities).getLeft();

        final boolean filtered = overallArtifactProbability > getArtifactProbabilityThreshold() - EPSILON;

        if (filtered) {
            expectedFalseNegatives.add(1 - overallArtifactProbability);
        } else {
            passingVariants.increment();
            expectedFalsePositives.add(overallArtifactProbability);
            expectedTruePositives.add(1 - overallArtifactProbability);
        }

        for (final Map.Entry<Mutect2VariantFilter, Double> entry : artifactProbabilities.entrySet()) {
            final String filter = entry.getKey().filterName();
            final double artifactProbability = entry.getValue();

            final Optional<String> posteriorAnnotation = entry.getKey().phredScaledPosteriorAnnotationName();
            if (posteriorAnnotation.isPresent() && entry.getKey().requiredAnnotations().stream().allMatch(vc::hasAttribute)) {
                vcb.attribute(posteriorAnnotation.get(), QualityUtils.errorProbToQual(artifactProbability));
            }

            if (artifactProbability > EPSILON && artifactProbability > getArtifactProbabilityThreshold() - EPSILON) {
                vcb.filter(filter);
                expectedFalseNegativesPerFilter.get(filter).add(1 - overallArtifactProbability);
            } else if (!filtered) {
                expectedFalsePositivesPerFilter.get(filter).add(artifactProbability);
            }
        }
        return vcb.make();
    }

    public void accumulateArtifactPosteriorsAndBadHaplotypes(final VariantContext vc) {
        final double artifactProbability = overallAndTechnicalOnlyArtifactProbabilities(vc).getLeft();

        if (artifactProbability > getArtifactProbabilityThreshold() - EPSILON) {
            recordFilteredHaplotypes(vc);
        }

        addFirstPassArtifactProbability(artifactProbability);
    }

}
