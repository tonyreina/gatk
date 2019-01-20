package org.broadinstitute.hellbender.tools.walkers.mutect.filtering;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.broadinstitute.hellbender.utils.GATKProtectedVariantContextUtils;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Optional;

public class ContaminationFilter extends Mutect2VariantFilter {
    @Override
    public double calculateArtifactProbability(final VariantContext vc, final Mutect2FilteringInfo filteringInfo) {
        final double somaticPriorProb = Math.pow(10, filteringInfo.getLog10PriorOfSomaticVariant(vc));
        final List<ImmutablePair<Integer, Double>> depthsAndPosteriors = new ArrayList<>();

        for (final Genotype tumorGenotype : vc.getGenotypes()) {
            if (filteringInfo.getNormalSamples().contains(tumorGenotype.getSampleName())) {
                continue;
            }

            final double contamination = filteringInfo.getContaminationBySample().get(tumorGenotype.getSampleName());

            final double[] alleleFractions = GATKProtectedVariantContextUtils.getAttributeAsDoubleArray(tumorGenotype, VCFConstants.ALLELE_FREQUENCY_KEY,
                    () -> new double[] {1.0}, 1.0);
            final int maxFractionIndex = MathUtils.maxElementIndex(alleleFractions);
            final int[] ADs = tumorGenotype.getAD();
            final int altCount = ADs[maxFractionIndex + 1];   // AD is all alleles, while AF is alts only, hence the +1 offset
            final int depth = (int) MathUtils.sum(ADs);
            final double[] negativeLog10AlleleFrequencies = GATKProtectedVariantContextUtils.getAttributeAsDoubleArray(vc,
                    GATKVCFConstants.POPULATION_AF_VCF_ATTRIBUTE, () -> new double[]{Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY}, Double.POSITIVE_INFINITY);
            final double alleleFrequency = MathUtils.applyToArray(negativeLog10AlleleFrequencies, x -> Math.pow(10,-x))[maxFractionIndex];

            final double somaticLikelihood = 1.0 / (depth + 1);

            final double singleContaminantLikelihood = 2 * alleleFrequency * (1 - alleleFrequency) * MathUtils.binomialProbability(depth, altCount, contamination /2)
                    + MathUtils.square(alleleFrequency) * MathUtils.binomialProbability(depth, altCount, contamination);
            final double manyContaminantLikelihood = MathUtils.binomialProbability(depth, altCount, contamination * alleleFrequency);
            final double contaminantLikelihood = Math.max(singleContaminantLikelihood, manyContaminantLikelihood);
            final double posteriorProbOfContamination = (1 - somaticPriorProb) * contaminantLikelihood / ((1 - somaticPriorProb) * contaminantLikelihood + somaticPriorProb * somaticLikelihood);

            depthsAndPosteriors.add(ImmutablePair.of(altCount, posteriorProbOfContamination));
        }

        return weightedMedianPosteriorProbability(depthsAndPosteriors);
    }

    @Override
    public boolean isTechnicalArtifact() { return false; }

    @Override
    public String filterName() {
        return GATKVCFConstants.CONTAMINATION_FILTER_NAME;
    }

    @Override
    public Optional<String> phredScaledPosteriorAnnotationName() {
        return Optional.of(GATKVCFConstants.CONTAMINATION_QUAL_ATTRIBUTE);
    }

    @Override
    protected List<String> requiredAnnotations() { return Collections.singletonList(GATKVCFConstants.POPULATION_AF_VCF_ATTRIBUTE); }
}
