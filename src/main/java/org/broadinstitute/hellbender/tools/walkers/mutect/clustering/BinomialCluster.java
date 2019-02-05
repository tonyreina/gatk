package org.broadinstitute.hellbender.tools.walkers.mutect.clustering;

import org.broadinstitute.hellbender.tools.walkers.readorientation.BetaDistributionShape;
import org.broadinstitute.hellbender.utils.MathUtils;

import java.util.List;

public class BinomialCluster implements AlleleFractionCluster {

    private static final double STD_DEV_OVER_MEAN = 0.01;




    @Override
    public void learn(final List<Datum> data) {
        final long altCount = data.stream().mapToInt(Datum::getAltCount).sum();
        final long totalCount = data.stream().mapToInt(Datum::getTotalCount).sum();
        
        betaDistributionShape = new BetaDistributionShape(alpha, beta);

    }




    private static BetaDistributionShape getFuzzyBinomial(final double mean) {
        final double alphaPlusBeta = ((1 - mean) / (mean * MathUtils.square(STD_DEV_OVER_MEAN))) - 1;
        final double alpha = mean * alphaPlusBeta;
        final double beta = alphaPlusBeta - alpha;
        return new BetaDistributionShape(alpha, beta);
    }
}
