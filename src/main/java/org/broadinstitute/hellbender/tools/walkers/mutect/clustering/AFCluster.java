package org.broadinstitute.hellbender.tools.walkers.mutect.clustering;

import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.math3.special.Gamma;
import org.broadinstitute.hellbender.tools.walkers.readorientation.BetaDistributionShape;
import org.broadinstitute.hellbender.utils.MathUtils;

import java.util.HashSet;
import java.util.Set;

public class AFCluster {
    private BetaDistributionShape betaShape;
    private static double STD_DEV_OVER_MEAN = 0.01;
    private final boolean isBackgroundCluster;

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

    private static BetaDistributionShape getFuzzyBinomial(final double mean, final double stdDevOverMean) {
        final double alphaPlusBeta = ((1 - mean) / (mean * MathUtils.square(stdDevOverMean))) - 1;
        final double alpha = mean * alphaPlusBeta;
        final double beta = alphaPlusBeta - alpha;
        return new BetaDistributionShape(alpha, beta);
    }
}
