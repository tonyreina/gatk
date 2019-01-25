package org.broadinstitute.hellbender.tools.walkers.mutect.filtering;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import org.apache.commons.lang3.mutable.MutableDouble;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.math3.util.MathArrays;
import org.broadinstitute.hellbender.tools.walkers.mutect.Mutect2Engine;
import org.broadinstitute.hellbender.tools.walkers.mutect.MutectStats;
import org.broadinstitute.hellbender.utils.GATKProtectedVariantContextUtils;
import org.broadinstitute.hellbender.utils.IndexRange;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.io.File;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Simple class to hold all information needed by Mutect2 filters, including: the M2FilteringArgumentCollection,
 * the set of normal samples, the samples' contamination, and the tumor sample CNV segments
 */
public class Mutect2FilteringEngine {
    public static final double EPSILON = 1.0e-10;

    private final List<Mutect2VariantFilter> filters = new ArrayList<>();
    private final Set<String> normalSamples;

    // TODO: make this private to BadHaplotypeFilter
    // TODO: make it probabilistic and eliminate the FIRST_PASS_THRESHOLD
    // for each PID, the positions with PGTs of filtered genotypes
    private final Map<String, ImmutablePair<Integer, Set<String>>> filteredPhasedCalls;

    /**
     * DATA ACCUMULATED AND LEARNED ON EACH PASS OF {@link FilterMutectCalls}
     */
    private final ThresholdCalculator thresholdCalculator;
    private final FilteringOutputStats filteringOutputStats;
    private final SomaticPriorModel somaticPriorModel;

    public Mutect2FilteringEngine(M2FiltersArgumentCollection MTFAC, final VCFHeader vcfHeader, final File mutectStatsTable) {
        thresholdCalculator = new ThresholdCalculator(MTFAC.thresholdStrategy, MTFAC.initialPosteriorThreshold, MTFAC.maxFalsePositiveRate, MTFAC.fScoreBeta);

        somaticPriorModel = new SomaticPriorModel(MTFAC, mutectStatsTable.exists() ? MutectStats.readFromFile(mutectStatsTable) : Collections.emptyList());

        normalSamples = vcfHeader.getMetaDataInInputOrder().stream()
                .filter(line -> line.getKey().equals(Mutect2Engine.NORMAL_SAMPLE_KEY_IN_VCF_HEADER))
                .map(VCFHeaderLine::getValue)
                .collect(Collectors.toSet());

        filteredPhasedCalls = new HashMap<>();

        buildFiltersList(MTFAC);
        filteringOutputStats = new FilteringOutputStats(filters);
    }

    //THE FOLLOWING ARE HELPER METHODS FOR FILTERS THAT IMPLEMENT {@link Mutect2VariantFilter}
    public boolean isNormal(final Genotype genotype) { return normalSamples.contains(genotype.getSampleName()); }

    public boolean isTumor(final Genotype genotype) { return !isNormal(genotype); }

    public Map<String, ImmutablePair<Integer, Set<String>>> getFilteredPhasedCalls() {
        return filteredPhasedCalls;
    }

    public double getArtifactProbabilityThreshold() { return thresholdCalculator.getThreshold(); }

    public double getLog10PriorOfSomaticVariant(final VariantContext vc) {
        return somaticPriorModel.getLog10PriorOfSomaticVariant(vc);
    }

    public double getPriorProbOfArtifactVersusVariant() {
        return somaticPriorModel.getPriorProbOfArtifactVersusVariant();
    }

    public static boolean hasPhaseInfo(final Genotype genotype) {
        return genotype.hasExtendedAttribute(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_GT_KEY) && genotype.hasExtendedAttribute(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_ID_KEY);
    }

    public int[] sumADsOverSamples(final VariantContext vc, final boolean includeTumor, final boolean includeNormal) {
        final int[] ADs = new int[vc.getNAlleles()];
        vc.getGenotypes().stream().filter(g -> (includeTumor && isTumor(g)) || (includeNormal && isNormal(g)))
                .map(Genotype::getAD).forEach(ad -> new IndexRange(0, vc.getNAlleles()).forEach(n -> ADs[n] += ad[n]));
        return ADs;
    }

    public double[] weightedAverageOfTumorAFs(final VariantContext vc) {
        final MutableDouble totalWeight = new MutableDouble(0);
        final double[] AFs = new double[vc.getNAlleles() - 1];
        vc.getGenotypes().stream().filter(this::isTumor).forEach(g ->  {
            final double weight = MathUtils.sum(g.getAD());
            totalWeight.add(weight);
            final double[] sampleAFs = GATKProtectedVariantContextUtils.getAttributeAsDoubleArray(g, VCFConstants.ALLELE_FREQUENCY_KEY,
                    () -> new double[] {0.0}, 0.0);
            MathArrays.scaleInPlace(weight, sampleAFs);
            MathUtils.addToArrayInPlace(AFs, sampleAFs);
        });
        MathArrays.scaleInPlace(1/totalWeight.getValue(), AFs);
        return AFs;
    }
    // END HELPER METHODS

    /**
     * record data from a potential variant in a non-final pass of {@link FilterMutectCalls}
     */
    public void accumulateData(final VariantContext vc) {
        filters.forEach(f -> f.accumulateDataForLearning(vc, this));

        final Pair<Double, Double> overallAndTechnicalArtifactProbs = overallAndTechnicalOnlyArtifactProbabilities(vc);
        somaticPriorModel.record(vc, overallAndTechnicalArtifactProbs.getLeft(), overallAndTechnicalArtifactProbs.getRight());

        // bad haplotypes and artifact posteriors
        if (overallAndTechnicalArtifactProbs.getLeft() > getArtifactProbabilityThreshold() - EPSILON) {
            recordFilteredHaplotypes(vc);
        }

        thresholdCalculator.addArtifactProbability(overallAndTechnicalArtifactProbs.getLeft());
    }

    /**
     * Refine model parameters based on data acquired in a non-final pass of {@link FilterMutectCalls}
     */
    public void learnParameters() {
        filters.forEach(Mutect2VariantFilter::learnParametersAndClearAccumulatedData);
        somaticPriorModel.learnAndClearAccumulatedData();
        thresholdCalculator.relearnThresholdAndClearAcumulatedProbabilities();

        filteringOutputStats.clear();
    }

    /**
     * Create a filtered variant and record statistics for the final pass of {@link FilterMutectCalls}
     */
    public VariantContext applyFiltersAndAccumulateOutputStats(final VariantContext vc) {
        final VariantContextBuilder vcb = new VariantContextBuilder(vc).filters(new HashSet<>());

        final Map<Mutect2VariantFilter, Double> artifactProbabilities = filters.stream()
                .collect(Collectors.toMap(f -> f, f -> f.errorProbability(vc, this)));

        final double overallArtifactProbability = overallAndTechnicalOnlyArtifactProbabilities(artifactProbabilities).getLeft();

        final boolean filtered = overallArtifactProbability > getArtifactProbabilityThreshold() - EPSILON;

        filteringOutputStats.recordCall(filtered, overallArtifactProbability, artifactProbabilities, getArtifactProbabilityThreshold());

        for (final Map.Entry<Mutect2VariantFilter, Double> entry : artifactProbabilities.entrySet()) {
            final double artifactProbability = entry.getValue();

            entry.getKey().phredScaledPosteriorAnnotationName().ifPresent(annotation -> {
                if (entry.getKey().requiredAnnotations().stream().allMatch(vc::hasAttribute)) {
                    vcb.attribute(annotation, QualityUtils.errorProbToQual(artifactProbability));
                }
            });

            // TODO: clarify this logic
            if (artifactProbability > EPSILON && artifactProbability > getArtifactProbabilityThreshold() - EPSILON) {
                vcb.filter(entry.getKey().filterName());
            }
        }

        return vcb.make();
    }

    /**
     * Write statistics collected in the final pass of {@link FilterMutectCalls}
     * @param filteringStatsFile
     */
    public void writeFilteringStats(final File filteringStatsFile) {
        filteringOutputStats.writeFilteringStats(filteringStatsFile, getArtifactProbabilityThreshold());
    }

    private void buildFiltersList(final M2FiltersArgumentCollection MTFAC) {
        filters.add(new TumorEvidenceFilter());
        filters.add(new BaseQualityFilter(MTFAC.minMedianBaseQuality));
        filters.add(new MappingQualityFilter(MTFAC.minMedianMappingQuality, MTFAC.longIndelLength));
        filters.add(new DuplicatedAltReadFilter(MTFAC.uniqueAltReadCount));
        filters.add(new StrandArtifactFilter());
        filters.add(new ContaminationFilter(MTFAC.contaminationTables, MTFAC.contaminationEstimate));
        filters.add(new PanelOfNormalsFilter());
        filters.add(new NormalArtifactFilter());
        filters.add(new ReadOrientationFilter());
        filters.add(new NRatioFilter(MTFAC.nRatio));
        filters.add(new StrictStrandBiasFilter(MTFAC.minReadsOnEachStrand));
        filters.add(new ReadPositionFilter(MTFAC.minMedianReadPosition));

        if (MTFAC.mitochondria) {
            filters.add(new LogOddsOverDepthFilter(MTFAC.minLog10OddsDividedByDepth));
            filters.add(new ChimericOriginalAlignmentFilter(MTFAC.maxNuMTFraction));
        } else {
            filters.add(new ClusteredEventsFilter(MTFAC.maxEventsInRegion));
            filters.add(new MultiallelicFilter(MTFAC.numAltAllelesThreshold));
            filters.add(new FragmentLengthFilter(MTFAC.maxMedianFragmentLengthDifference));
            filters.add(new PolymeraseSlippageFilter(MTFAC.minSlippageLength, MTFAC.slippageRate));
            filters.add(new FilteredHaplotypeFilter(MTFAC.maxDistanceToFilteredCallOnSameHaplotype));
            filters.add(new GermlineFilter(MTFAC.tumorSegmentationTables));
        }
    }

    private void recordFilteredHaplotypes(final VariantContext vc) {
        final Map<String, Set<String>> phasedGTsForEachPhaseID = vc.getGenotypes().stream()
                .filter(gt -> !normalSamples.contains(gt.getSampleName()))
                .filter(Mutect2FilteringEngine::hasPhaseInfo)
                .collect(Collectors.groupingBy(g -> (String) g.getExtendedAttribute(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_ID_KEY, ""),
                        Collectors.mapping(g -> (String) g.getExtendedAttribute(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_GT_KEY, ""), Collectors.toSet())));

        for (final Map.Entry<String, Set<String>> pidAndPgts : phasedGTsForEachPhaseID.entrySet()) {
            filteredPhasedCalls.put(pidAndPgts.getKey(), new ImmutablePair<>(vc.getStart(), pidAndPgts.getValue()));
        }
    }

    private Pair<Double, Double> overallAndTechnicalOnlyArtifactProbabilities(final VariantContext vc) {
        return overallAndTechnicalOnlyArtifactProbabilities(filters.stream()
                .collect(Collectors.toMap(f -> f, f -> f.calculateErrorProbability(vc, this))));
    }

    private Pair<Double, Double> overallAndTechnicalOnlyArtifactProbabilities(final Map<Mutect2VariantFilter, Double> map) {
        final double technicalArtifactProbability = map.entrySet().stream().filter(entry -> entry.getKey().errorType())
                .mapToDouble(entry -> entry.getValue()).max().orElseGet(() -> 0);
        final double nonTechnicalArtifactProbability = map.entrySet().stream().filter(entry -> !entry.getKey().errorType())
                .mapToDouble(entry -> entry.getValue()).max().orElseGet(() -> 0);

        final double overallProbability = 1 - (1 - technicalArtifactProbability) * (1 - nonTechnicalArtifactProbability);

        return ImmutablePair.of(overallProbability, technicalArtifactProbability);
    }

    private final class ErrorProbabilities {
        private final Map<Mutect2VariantFilter, Double> probabilitiesByFilter;
        private final EnumMap<ErrorType, Double> probabilitiesByType;


        public ErrorProbabilities(final VariantContext vc) {
            probabilitiesByFilter = filters.stream().collect(Collectors.toMap(f -> f, f -> f.errorProbability(vc, this)));

            probabilitiesByType = new EnumMap<>(ErrorType.class);

            final Map<ErrorType, List<Map.Entry<Mutect2VariantFilter, Double>>> doogle = probabilitiesByFilter.entrySet().stream().collect(Collectors.groupingBy(entry -> entry.getKey().errorType()));


        }
    }


}
