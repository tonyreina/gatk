package org.broadinstitute.hellbender.tools.walkers.mutect.filtering;

import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.lang3.mutable.MutableDouble;
import org.broadinstitute.hellbender.utils.MathUtils;

public class SomaticPriorModel {
    private double log10SNVPrior;
    private double log10IndelPrior;
    private double artifactVsVariantPrior;

    private final MutableDouble realVariantCount = new MutableDouble(0);
    private final MutableDouble realSNVCount = new MutableDouble(0);
    private final MutableDouble realIndelCount = new MutableDouble(0);
    private final MutableDouble technicalArtifactCount = new MutableDouble(0);

    public SomaticPriorModel(final double initialLog10SNVPrior, final double initialLog10IndelPrior,
                             final double initialArtifactVsVariantPrior, final ) {
        log10SNVPrior = initialLog10SNVPrior;
        log10IndelPrior = initialLog10IndelPrior;
        artifactVsVariantPrior = initialArtifactVsVariantPrior;
    }

    public double getLog10PriorOfSomaticVariant(final VariantContext vc) {
        return vc.isSNP() ? (MathUtils.LOG10_ONE_THIRD + log10SNVPrior) : log10IndelPrior;
    }

    public double getPriorProbOfArtifactVersusVariant() { return artifactVsVariantPrior; }

    public void record(final VariantContext vc, final double artifactProbability, final double technicalArtifactProbability) {
        realVariantCount.add(1 - artifactProbability);
        (vc.isSNP() ? realSNVCount : realIndelCount).add(1 - artifactProbability);
        technicalArtifactCount.add(technicalArtifactProbability);
    }

    // by default, clear accumulated data after learning
    public void learnAndClearAccumulatedData() {
        learn();
        clear();
    }

    public void learn() {
        // global parameters
        artifactVsVariantPrior = (technicalArtifactCount.getValue() + 1) / (realVariantCount.getValue() + technicalArtifactCount.getValue() + 2);
        if (totalCallableSites.isPresent()) {
            log10SNVPrior = Math.log10(realSNVCount.getValue() / totalCallableSites.getAsLong());
            log10IndelPrior = Math.log10(realIndelCount.getValue() / totalCallableSites.getAsLong());
        }
    }

    public void clear() {
        realVariantCount.setValue(0);
        realSNVCount.setValue(0);
        realIndelCount.setValue(0);
        technicalArtifactCount.setValue(0);
    }
}
