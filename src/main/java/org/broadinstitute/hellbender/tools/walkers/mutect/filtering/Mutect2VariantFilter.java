package org.broadinstitute.hellbender.tools.walkers.mutect.filtering;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.math3.util.MathArrays;
import org.broadinstitute.hellbender.utils.GATKProtectedVariantContextUtils;
import org.broadinstitute.hellbender.utils.IndexRange;
import org.broadinstitute.hellbender.utils.MathUtils;

import java.util.Comparator;
import java.util.List;
import java.util.Optional;
import java.util.Set;

public abstract class Mutect2VariantFilter {
    public Mutect2VariantFilter() { }

    public double artifactProbability(final VariantContext vc, final Mutect2FilteringInfo filteringInfo) {
        return requiredAnnotations().stream().allMatch(vc::hasAttribute) ? calculateArtifactProbability(vc, filteringInfo) : 0;
    }

    protected abstract double calculateArtifactProbability(final VariantContext vc, final Mutect2FilteringInfo filteringInfo);

    // by default do nothing, but we may override to allow some filters to learn their parameters in the first pass of {@link FilterMutectCalls}
    protected void accumulateDataForLearning(final VariantContext vc, final Mutect2FilteringInfo filteringInfo) { }
    protected void learnParameters() { }

    public abstract String filterName();

    public abstract Optional<String> phredScaledPosteriorAnnotationName();

    protected abstract List<String> requiredAnnotations();

    protected static int[] sumADsOverSamples(final VariantContext vc, final Set<String> normalSamples, final boolean includeTumor, final boolean includeNormal) {
        final int[] ADs = new int[vc.getNAlleles()];
        vc.getGenotypes().stream()
                .filter(g -> includeTumor || normalSamples.contains(g.getSampleName()))
                .filter(g -> includeNormal || !normalSamples.contains(g.getSampleName()))
                .map(Genotype::getAD).forEach(ad -> new IndexRange(0, vc.getNAlleles()).forEach(n -> ADs[n] += ad[n]));
        return ADs;
    }

    protected static double[] weightedAverageOfTumorAFs(final VariantContext vc, final Set<String> normalSamples) {
        double totalWeight = 0.0;
        final double[] AFs = new double[vc.getNAlleles() - 1];
        for (final Genotype g : vc.getGenotypes()) {
            if (normalSamples.contains(g.getSampleName())) {
                continue;
            } else {
                final double weight = MathUtils.sum(g.getAD());
                totalWeight += weight;
                final double[] sampleAFs = GATKProtectedVariantContextUtils.getAttributeAsDoubleArray(g, VCFConstants.ALLELE_FREQUENCY_KEY,
                        () -> new double[] {0.0}, 0.0);
                MathArrays.scaleInPlace(weight, sampleAFs);
                MathUtils.addToArrayInPlace(AFs, sampleAFs);
            }
        }
        MathArrays.scaleInPlace(1/totalWeight, AFs);
        return AFs;
    }

    // weighted median -- what's the lowest posterior probability that accounts for samples with half of the total alt depth
    protected static double weightedMedianPosteriorProbability(List<ImmutablePair<Integer, Double>> depthsAndPosteriors) {
        final int totalAltDepth = depthsAndPosteriors.stream().mapToInt(ImmutablePair::getLeft).sum();

        // sort from lowest to highest posterior probability of artifact
        depthsAndPosteriors.sort(Comparator.comparingDouble(p -> p.getRight()));

        int cumulativeAltCount = 0;

        for (final ImmutablePair<Integer, Double> pair : depthsAndPosteriors) {
            cumulativeAltCount += pair.getLeft();
            if (cumulativeAltCount * 2 >= totalAltDepth) {
                return pair.getRight();
            }
        }
        return 0;
    }
}
