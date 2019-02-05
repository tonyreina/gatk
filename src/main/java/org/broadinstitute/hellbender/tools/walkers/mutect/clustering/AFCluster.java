package org.broadinstitute.hellbender.tools.walkers.mutect.clustering;

import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.math3.special.Gamma;
import org.broadinstitute.hellbender.tools.walkers.readorientation.BetaDistributionShape;
import org.broadinstitute.hellbender.utils.MathUtils;

import java.util.HashSet;
import java.util.Set;

// an allele fraction cluster is a beta prior on allele fraction with a method for
// adjusting a tumor log-likelihood to account for this beta prior as opposed to the flat prior used for Mutect2's
// initial calls
public class AFCluster {
    private BetaDistributionShape betaShape;
    private static double STD_DEV_OVER_MEAN = 0.01;
    private final Set<Datum> members = new HashSet<>();    // the indices of data in this cluster
    private final boolean isBackgroundCluster;

    private static int totalAssignments = 0;

    public static AFCluster makeBackgroundCuster(final BetaDistributionShape betaShape) {
        return new AFCluster(betaShape, true);
    }

    public static AFCluster makeCluster(final double alleleFraction) {
        return new AFCluster(getFuzzyBinomial(alleleFraction, STD_DEV_OVER_MEAN), false);
    }

    private AFCluster(final BetaDistributionShape betaShape, final boolean isBackgroundCluster) {
        this.betaShape = betaShape;
        this.isBackgroundCluster = isBackgroundCluster;
    }

    public void remove(final Datum datum) {
        members.remove(datum);
        totalAssignments--;
    }

    public static void clearAssignmentCount() {
        totalAssignments = 0;
    }

    public void add(final Datum datum) {
        members.add(datum);
        totalAssignments++;
    }

    public int size() { return members.size(); }

    public int SNVCount() { return (int) members.stream().filter(datum -> datum.getType() == VariantContext.Type.SNP).count(); }

    public boolean isEmpty() { return members.isEmpty(); }

    public void relearn() {
        final double altCount = members.stream().mapToDouble(Datum::getAltCount).sum();
        final long totalCount = members.stream().mapToInt(Datum::getTotalCount).sum();
        if (totalCount == 0) {
            members.clear();
            return;
        }
        final double mean = Math.min(altCount / totalCount, 1 - STD_DEV_OVER_MEAN);


        if (isBackgroundCluster) {
            //TODO: learn background cluster here
            final double rate = 0.01;
            final double maxStep = 0.1;

            double alpha = betaShape.getAlpha();
            double beta = betaShape.getBeta();

            for (int epoch = 0; epoch < 10; epoch++) {
                for (final Datum datum : members) {
                    final int alt = datum.getAltCount();
                    final int total = datum.getTotalCount();
                    final int ref = total - alt;

                    //TODO re-use terms
                    final double alphaGradient = Gamma.digamma(alpha + alt) - Gamma.digamma(total + alpha + beta) - Gamma.digamma(alpha) + Gamma.digamma(alpha + beta);
                    final double betaGradient = Gamma.digamma(beta + ref) - Gamma.digamma(total + alpha + beta) - Gamma.digamma(beta) + Gamma.digamma(alpha + beta);

                    alpha = alpha + rate * alphaGradient;
                    beta = beta + rate * betaGradient;
                }
                int j = 10;
            }

            betaShape = new BetaDistributionShape(alpha, beta);


            // TODO: THIS IS EMPTY!!!!!
        } else {
            betaShape = getFuzzyBinomial(mean, STD_DEV_OVER_MEAN);
        }
    }

    public double getLog10ChineseRestaurantFactor(final double concentration) {
        return Math.log10((double) members.size() / (totalAssignments + concentration));
    }

    public static double getNewClusterLog10ChineseRestaurantFactor(final double concentration) {
        return Math.log10(concentration / (totalAssignments + concentration));
    }

    public double log10Likelihood(final Datum datum) {
        final int altCount = datum.getAltCount();
        final int refCount = datum.getTotalCount() - altCount;
        return datum.getTumorLog10Odds() + SomaticClusteringModel.log10OddsCorrection(SomaticClusteringModel.FLAT_BETA, betaShape, altCount, refCount);
    }

    private static BetaDistributionShape getFuzzyBinomial(final double mean, final double stdDevOverMean) {
        final double alphaPlusBeta = ((1 - mean) / (mean * MathUtils.square(stdDevOverMean))) - 1;
        final double alpha = mean * alphaPlusBeta;
        final double beta = alphaPlusBeta - alpha;
        return new BetaDistributionShape(alpha, beta);
    }
}
