package org.broadinstitute.hellbender.tools.walkers.mutect.filtering;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import org.apache.commons.lang3.mutable.MutableDouble;
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

        buildFiltersList(MTFAC);
        filteringOutputStats = new FilteringOutputStats(filters);
    }

    //THE FOLLOWING ARE HELPER METHODS FOR FILTERS THAT IMPLEMENT {@link Mutect2VariantFilter}
    public boolean isNormal(final Genotype genotype) { return normalSamples.contains(genotype.getSampleName()); }

    public boolean isTumor(final Genotype genotype) { return !isNormal(genotype); }

    public double getArtifactProbabilityThreshold() { return thresholdCalculator.getThreshold(); }

    public double getLog10PriorOfSomaticVariant(final VariantContext vc) {
        return somaticPriorModel.getLog10PriorOfSomaticVariant(vc);
    }

    /**
     * Correct the tumor log-10 odds (ie the TLOD) from Mutect2 to account for allele fraction clustering
     * @param tumorLog10Odds the original tumor log-10 odds of Mutect2 that assumes a flat prior on allele fractions
     * @return the log-10 odds corrected for the allele fraction clustering learned by the {@link SomaticPriorModel}
     */
    public double clusteringCorrectedLog10Odds(final double tumorLog10Odds, final int altCount, final int refCount) {
        return somaticPriorModel.clusteringCorrectedLog10Odds(tumorLog10Odds, altCount, refCount);
    }

    //TODO: need a method to renormalize log10Lieklihoods from SomaticPriorodel clustering

    public double getPriorProbOfArtifactVersusVariant() {
        return somaticPriorModel.getPriorProbOfArtifactVersusVariant();
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
        final ErrorProbabilities errorProbabilities = new ErrorProbabilities(filters, vc, this);
        filters.forEach(f -> f.accumulateDataForLearning(vc, errorProbabilities, this));
        final int[] tumorADs = sumADsOverSamples(vc, true, false);
        final double[] tumorLog10Odds = GATKProtectedVariantContextUtils.getAttributeAsDoubleArray(vc, GATKVCFConstants.TUMOR_LOD_KEY);

        // TODO: this needs to get both technical and non-somatic error probabilities
        somaticPriorModel.record(tumorADs, tumorLog10Odds, errorProbabilities.getTechnicalArtifactProbability(),
                errorProbabilities.getNonSomaticProbability(), vc.getType());
        thresholdCalculator.addArtifactProbability(errorProbabilities.getErrorProbability());
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

        final ErrorProbabilities errorProbabilities = new ErrorProbabilities(filters, vc, this);
        filteringOutputStats.recordCall(errorProbabilities, getArtifactProbabilityThreshold() - EPSILON);

        for (final Map.Entry<Mutect2VariantFilter, Double> entry : errorProbabilities.getProbabilitiesByFilter().entrySet()) {
            final double errorProbability = entry.getValue();

            entry.getKey().phredScaledPosteriorAnnotationName().ifPresent(annotation -> {
                if (entry.getKey().requiredAnnotations().stream().allMatch(vc::hasAttribute)) {
                    vcb.attribute(annotation, QualityUtils.errorProbToQual(errorProbability));
                }
            });

            // TODO: clarify this logic
            if (errorProbability > EPSILON && errorProbability > getArtifactProbabilityThreshold() - EPSILON) {
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
        filteringOutputStats.writeFilteringStats(filteringStatsFile, getArtifactProbabilityThreshold(), somaticPriorModel.clusteringMetadata());
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
}
