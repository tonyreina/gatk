package org.broadinstitute.hellbender.tools.walkers.mutect.clustering;

import java.util.List;

public class SequencingError implements AlleleFractionCluster {

    public SequencingError() { }

    @Override
    public double log10Likelihood(final Datum datum) {
        return 0;
    }

    @Override
    public void learn(final List<Datum> data) { }

    @Override
    public String toString() {
        return "sequencing error";
    }
}
