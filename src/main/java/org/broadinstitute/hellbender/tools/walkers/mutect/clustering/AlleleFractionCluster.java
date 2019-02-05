package org.broadinstitute.hellbender.tools.walkers.mutect.clustering;

import java.util.List;

public interface AlleleFractionCluster {
    double log10Likelihood(final Datum datum);

    void learn(final List<Datum> data);
}
