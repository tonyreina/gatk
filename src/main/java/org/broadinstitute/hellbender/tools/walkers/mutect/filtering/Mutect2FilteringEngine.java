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
    private static final double EPSILON = 1.0e-10;

    private final List<Mutect2VariantFilter> filters = new ArrayList<>();

    /**
     * CONSTANT PARAMETERS
     */
    private final Set<String> normalSamples;
    private OptionalLong totalCallableSites = OptionalLong.empty();


    // TODO: make this private to BadHaplotypeFilter
    // TODO: make it probabilistic and eliminate the FIRST_PASS_THRESHOLD
    // for each PID, the positions with PGTs of filtered genotypes
    private final Map<String, ImmutablePair<Integer, Set<String>>> filteredPhasedCalls;

    /**
     * DATA ACCUMULATED AND LEARNED ON EACH PASS OF {@link FilterMutectCalls}
     */
    private final ThresholdCalculator thresholdCalculator;
    private final OutputStats outputStats = new OutputStats();
    private final SomaticPriorModel somaticPriorModel;

    public Mutect2FilteringEngine(M2FiltersArgumentCollection MTFAC, final VCFHeader vcfHeader) {
        thresholdCalculator = new ThresholdCalculator(MTFAC.thresholdStrategy, MTFAC.initialPosteriorThreshold, MTFAC.maxFalsePositiveRate, MTFAC.fScoreBeta);
        somaticPriorModel = new SomaticPriorModel(MTFAC.log10PriorProbOfSomaticSNV, MTFAC.log10PriorProbOfSomaticIndel, MTFAC.initialPriorOfArtifactVersusVariant);

        normalSamples = vcfHeader.getMetaDataInInputOrder().stream()
                .filter(line -> line.getKey().equals(Mutect2Engine.NORMAL_SAMPLE_KEY_IN_VCF_HEADER))
                .map(VCFHeaderLine::getValue)
                .collect(Collectors.toSet());

        filteredPhasedCalls = new HashMap<>();

        buildFiltersList(MTFAC);
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

    public void inputMutectStats(final File mutectStatsTable) {
        final Map<String, Double> stats = MutectStats.readFromFile(mutectStatsTable).stream()
                .collect(Collectors.toMap(MutectStats::getStatistic, MutectStats::getValue));

        totalCallableSites = OptionalLong.of(Math.round(stats.get(Mutect2Engine.CALLABLE_SITES_NAME)));
    }

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

    public static boolean hasPhaseInfo(final Genotype genotype) {
        return genotype.hasExtendedAttribute(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_GT_KEY) && genotype.hasExtendedAttribute(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_ID_KEY);
    }

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

    public void learnParameters() {
        filters.forEach(Mutect2VariantFilter::learnParametersAndClearAccumulatedData);
        somaticPriorModel.learnAndClearAccumulatedData();
        thresholdCalculator.relearnThresholdAndClearAcumulatedProbabilities();

        outputStats.clear();
    }

    public VariantContext applyFiltersAndAccumulateOutputStats(final VariantContext vc) {
        final VariantContextBuilder vcb = new VariantContextBuilder(vc).filters(new HashSet<>());

        final Map<Mutect2VariantFilter, Double> artifactProbabilities = filters.stream()
                .collect(Collectors.toMap(f -> f, f -> f.artifactProbability(vc, this)));

        final double overallArtifactProbability = overallAndTechnicalOnlyArtifactProbabilities(artifactProbabilities).getLeft();

        final boolean filtered = overallArtifactProbability > getArtifactProbabilityThreshold() - EPSILON;

        outputStats.recordCall(filtered, overallArtifactProbability, artifactProbabilities);

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

    public void writeFilteringStats(final File filteringStatsFile) {
        outputStats.writeFilteringStats(filteringStatsFile);
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

    /**
     * Helper class used on the final pass of {@link FilterMutectCalls} to record total expected true positives, false positives,
     * and false negatives, as well as false positives and false negatives attributable to each filter
     */
    private class OutputStats {
        private int pass = 0;
        private double TPs = 0;
        private double FPs = 0;
        private double FNs = 0;

        private Map<Mutect2VariantFilter, MutableDouble> filterFPs = makeEmptyFilterCounts();
        private Map<Mutect2VariantFilter, MutableDouble> filterFNs = makeEmptyFilterCounts();

        public void recordCall(final boolean filtered, final double artifactProbability,
                               final Map<Mutect2VariantFilter, Double> artifactProbabilities) {
            if (filtered) {
                FNs += 1 - artifactProbability;
            } else {
                pass++;
                FPs += artifactProbability;
                TPs += 1 - artifactProbability;
            }

            for (final Map.Entry<Mutect2VariantFilter, Double> entry : artifactProbabilities.entrySet()) {
                final double filterArtifactProbability = entry.getValue();
                if (filterArtifactProbability > EPSILON && filterArtifactProbability > getArtifactProbabilityThreshold() - EPSILON) {
                    filterFNs.get(entry.getKey()).add(1 - artifactProbability);
                } else if (!filtered) {
                    filterFPs.get(entry.getKey()).add(filterArtifactProbability);
                }
            }
        }

        public void writeFilteringStats(final File filteringStatsFile) {
            final double totalTrueVariants = TPs + FNs;

            final List<FilterStats> filterStats = filters.stream()
                    .map(f -> new FilterStats(f.filterName(), filterFPs.get(f).getValue(), filterFPs.get(f).getValue() / pass,
                            filterFNs.get(f).getValue(), filterFNs.get(f).getValue() / totalTrueVariants))
                    .filter(stats -> stats.getFalsePositiveCount() > 0 || stats.getFalseNegativeCount() > 0)
                    .collect(Collectors.toList());

            FilterStats.writeM2FilterSummary(filterStats, filteringStatsFile, getArtifactProbabilityThreshold(), pass, TPs, FPs, FNs);
        }

        private Map<Mutect2VariantFilter, MutableDouble> makeEmptyFilterCounts() {
            return filters.stream().collect(Collectors.toMap(f -> f, f -> new MutableDouble(0)));
        }


        public void clear() {
            pass = 0;
            TPs = 0;
            FPs = 0;
            FNs = 0;

            filterFPs = makeEmptyFilterCounts();
            filterFNs = makeEmptyFilterCounts();
        }
    }

}
