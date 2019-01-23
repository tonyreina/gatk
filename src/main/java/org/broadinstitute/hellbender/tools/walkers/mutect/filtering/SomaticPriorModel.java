package org.broadinstitute.hellbender.tools.walkers.mutect.filtering;

public class SomaticPriorModel {
    private double log10SNVPrior;
    private double log10IndelPrior;
    private double artifactVsVariantPrior;

    public SomaticPriorModel(final double initialLog10SNVPrior, final double initialLog10IndelPrior,
                             final double initialArtifactVsVariantPrior) {
        log10SNVPrior = initialLog10SNVPrior;
        log10IndelPrior = initialLog10IndelPrior;
        artifactVsVariantPrior = initialArtifactVsVariantPrior;
    }


    public void clear() {

    }
}
