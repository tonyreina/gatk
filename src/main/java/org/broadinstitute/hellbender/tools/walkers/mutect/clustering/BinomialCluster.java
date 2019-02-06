package org.broadinstitute.hellbender.tools.walkers.mutect.clustering;

import org.broadinstitute.hellbender.tools.walkers.readorientation.BetaDistributionShape;
import org.broadinstitute.hellbender.utils.MathUtils;

import java.util.List;

public class BinomialCluster implements AlleleFractionCluster {

    private static final double STD_DEV_OVER_MEAN = 0.01;

    private BetaDistributionShape betaDistributionShape;

    public BinomialCluster(final double mean) {
        betaDistributionShape = getFuzzyBinomial(mean);
    }

    @Override
    public double log10Likelihood(final Datum datum) {
        return BetaBinomialCluster.log10Likelihood(datum, betaDistributionShape);
    }

    @Override
    public void learn(final List<Datum> data) {
        final double altCount = data.stream().mapToInt(Datum::getAltCount).sum() + 0.0001;
        final double totalCount = data.stream().mapToInt(Datum::getTotalCount).sum() + 0.0001;
        betaDistributionShape = getFuzzyBinomial((double) altCount / totalCount);

    }

    private static BetaDistributionShape getFuzzyBinomial(final double mean) {
        final double alphaPlusBeta = ((1 - mean) / (mean * MathUtils.square(STD_DEV_OVER_MEAN))) - 1;
        final double alpha = mean * alphaPlusBeta;
        final double beta = alphaPlusBeta - alpha;
        return new BetaDistributionShape(alpha, beta);
    }

    @Override
    public String toString() {
        return String.format("mean = %.3f", betaDistributionShape.getAlpha() / (betaDistributionShape.getAlpha() + betaDistributionShape.getBeta()));
    }
}
