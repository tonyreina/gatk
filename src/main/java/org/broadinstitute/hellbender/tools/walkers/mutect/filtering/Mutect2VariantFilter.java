package org.broadinstitute.hellbender.tools.walkers.mutect.filtering;

import com.google.common.annotations.VisibleForTesting;
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
    protected void clearAccumulatedData() { }
    protected void learnParameters() { }
    protected void learnParametersAndClearAccumulatedData() {
        learnParameters();
        clearAccumulatedData();
    }

    // by default assume that anything filtered is a technical artifact, but some filters, for example the germline and
    // contamination filters involve evidence of real, non-somatic variation.  These should not inform our models of
    // technical artifacts.
    public boolean isTechnicalArtifact() { return true; }

    public abstract String filterName();

    public abstract Optional<String> phredScaledPosteriorAnnotationName();

    protected abstract List<String> requiredAnnotations();

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

    @VisibleForTesting
    static double posteriorProbabilityOfError(final double log10OddsOfRealVersusError, final double log10PriorOfReal) {
        final double[] unweightedPosteriorOfRealAndError = new double[] {log10OddsOfRealVersusError + log10PriorOfReal,
                MathUtils.log10OneMinusPow10(log10PriorOfReal)};

        final double[] posteriorOfRealAndError = MathUtils.normalizeFromLog10ToLinearSpace(unweightedPosteriorOfRealAndError);

        return posteriorOfRealAndError[1];
    }
}
