package org.broadinstitute.hellbender.tools.walkers.mutect.clustering;

import htsjdk.variant.variantcontext.VariantContext;

public class Datum {
    private final double tumorLog10Odds;
    private final double artifactProb;
    private final double nonSequencingErrorProb;
    private final int altCount;
    private final int totalCount;
    private final VariantContext.Type type;

    private AFCluster cluster = null;

    public Datum(final double tumorLog10Odds, final double artifactProb, final double nonSomaticProb, final int altCount, final int totalCount, final VariantContext.Type type) {
        this.tumorLog10Odds = tumorLog10Odds;
        this.artifactProb = artifactProb;
        this.altCount = altCount;
        this.totalCount = totalCount;
        this.type = type;
        nonSequencingErrorProb = 1 - (1 - artifactProb) * (1 - nonSomaticProb);
    }

    public void unassign() {
        if (cluster != null) {
            cluster.remove(this);
            cluster = null;
        }
    }

    public void assign(final AFCluster cluster) {
        this.cluster = cluster;
        cluster.add(this);
    }

    public double getTumorLog10Odds() { return tumorLog10Odds; }

    public double getArtifactProb() { return artifactProb; }

    public double getNonSequencingErrorProb() { return nonSequencingErrorProb; }

    public int getAltCount() { return altCount; }

    public int getTotalCount() { return totalCount; }

    public VariantContext.Type getType() { return type; }

    public AFCluster getCluster() { return cluster; }
}
