package org.broadinstitute.hellbender.tools.walkers.mutect.filtering;

import htsjdk.samtools.util.OverlapDetector;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.lang3.mutable.MutableDouble;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.broadinstitute.hellbender.tools.walkers.contamination.MinorAlleleFractionRecord;
import org.broadinstitute.hellbender.tools.walkers.mutect.GermlineProbabilityCalculator;
import org.broadinstitute.hellbender.utils.GATKProtectedVariantContextUtils;
import org.broadinstitute.hellbender.utils.IndexRange;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.io.File;
import java.util.*;
import java.util.stream.Collectors;

public class GermlineFilter extends Mutect2VariantFilter {
    private static final double MIN_ALLELE_FRACTION_FOR_GERMLINE_HOM_ALT = 0.9;

    private final Map<String, OverlapDetector<MinorAlleleFractionRecord>> tumorSegments;

    public GermlineFilter(final List<File> tumorSegmentationTables) {
        tumorSegments = tumorSegmentationTables.stream()
                .map(MinorAlleleFractionRecord::readFromFile)
                .collect(Collectors.toMap(ImmutablePair::getLeft, p -> OverlapDetector.create(p.getRight())));
    }

    @Override
    public ErrorType errorType() { return ErrorType.NON_SOMATIC; }

    @Override
    public double calculateErrorProbability(final VariantContext vc, final Mutect2FilteringEngine filteringInfo) {
        final double[] tumorLog10OddsIfSomatic = GATKProtectedVariantContextUtils.getAttributeAsDoubleArray(vc, GATKVCFConstants.TUMOR_LOD_KEY);
        final Optional<double[]> normalLods = vc.hasAttribute(GATKVCFConstants.NORMAL_LOD_KEY) ?
                Optional.of(GATKProtectedVariantContextUtils.getAttributeAsDoubleArray(vc, GATKVCFConstants.NORMAL_LOD_KEY)) : Optional.empty();
        final double[] negativeLog10AlleleFrequencies = GATKProtectedVariantContextUtils.getAttributeAsDoubleArray(vc, GATKVCFConstants.POPULATION_AF_VCF_ATTRIBUTE);
        final double[] populationAlleleFrequencies = MathUtils.applyToArray(negativeLog10AlleleFrequencies, x -> Math.pow(10,-x));

        final MutableDouble weightedSumOfMafs = new MutableDouble(0);
        vc.getGenotypes().stream().filter(filteringInfo::isTumor).forEach(tumorGenotype -> {
            final String sample = tumorGenotype.getSampleName();
            final List<MinorAlleleFractionRecord> segments = tumorSegments.containsKey(sample) ? tumorSegments.get(sample).getOverlaps(vc).stream().collect(Collectors.toList())
                    : Collections.emptyList();

            // minor allele fraction -- we abbreviate the name to make the formulas below less cumbersome
            final double maf = segments.isEmpty() ? 0.5 : segments.get(0).getMinorAlleleFraction();

            weightedSumOfMafs.add(maf * MathUtils.sum(tumorGenotype.getAD()));
        });

        final double[] altAlleleFractions = filteringInfo.weightedAverageOfTumorAFs(vc);

        // note that this includes the ref
        final int[] alleleCounts = filteringInfo.sumADsOverSamples(vc, true, false);

        // weighted average of sample minor allele fractions.  This is the expected allele fraction of a germline het in the aggregated read counts
        final double maf = weightedSumOfMafs.getValue() / MathUtils.sum(alleleCounts);

        // exclude the ref
        final int[] altCounts = Arrays.copyOfRange(alleleCounts, 1, alleleCounts.length);

        final int refCount = alleleCounts[0];

        // this is \chi in the docs, the correction factor for tumor likelihoods if forced to have maf or 1 - maf
        // as the allele fraction
        final double[] log10OddsOfGermlineHetVsSomatic = new IndexRange(0, altAlleleFractions.length).mapToDouble(n -> {
            final double log10GermlineAltMinorLikelihood = refCount * Math.log10(1 - maf) + altCounts[n] * Math.log10(maf);
            final double log10GermlineAltMajorLikelihood = refCount * Math.log10(maf) + altCounts[n] * Math.log10(1 - maf);
            final double log10GermlineLikelihood = MathUtils.LOG10_ONE_HALF + MathUtils.log10SumLog10(log10GermlineAltMinorLikelihood, log10GermlineAltMajorLikelihood);

            final double log10SomaticLikelihood = refCount * Math.log10(1 - altAlleleFractions[n]) + altCounts[n] * Math.log10(altAlleleFractions[n]);
            return log10GermlineLikelihood - log10SomaticLikelihood;
        });

        // see docs -- basically the tumor likelihood for a germline hom alt is approximately equal to the somatic likelihood
        // as long as the allele fraction is high
        final double[] log10OddsOfGermlineHomAltVsSomatic = MathUtils.applyToArray(altAlleleFractions, x-> x < MIN_ALLELE_FRACTION_FOR_GERMLINE_HOM_ALT ? Double.NEGATIVE_INFINITY : 0);

        final double[] log10GermlinePosteriors = GermlineProbabilityCalculator.calculateGermlineProbabilities(
                populationAlleleFrequencies, log10OddsOfGermlineHetVsSomatic, log10OddsOfGermlineHomAltVsSomatic, normalLods, filteringInfo.getLog10PriorOfSomaticVariant(vc));

        final int indexOfMaxTumorLod = MathUtils.maxElementIndex(tumorLog10OddsIfSomatic);
        return Math.pow(10.0, log10GermlinePosteriors[indexOfMaxTumorLod]);
    }

    @Override
    public String filterName() {
        return GATKVCFConstants.GERMLINE_RISK_FILTER_NAME;
    }

    @Override
    public Optional<String> phredScaledPosteriorAnnotationName() {
        return Optional.of(GATKVCFConstants.GERMLINE_QUAL_VCF_ATTRIBUTE);
    }

    @Override
    protected List<String> requiredAnnotations() {
        return Arrays.asList(GATKVCFConstants.TUMOR_LOD_KEY, GATKVCFConstants.POPULATION_AF_VCF_ATTRIBUTE);
    }
}
