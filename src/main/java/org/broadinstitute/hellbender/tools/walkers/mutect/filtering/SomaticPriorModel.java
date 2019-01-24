package org.broadinstitute.hellbender.tools.walkers.mutect.filtering;

import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.lang3.mutable.MutableDouble;
import org.broadinstitute.hellbender.tools.walkers.mutect.Mutect2Engine;
import org.broadinstitute.hellbender.tools.walkers.mutect.MutectStats;
import org.broadinstitute.hellbender.utils.MathUtils;

import java.util.List;
import java.util.OptionalDouble;

public class SomaticPriorModel {
    private double log10SNVPrior;
    private double log10IndelPrior;
    private double artifactVsVariantPrior;
    private final OptionalDouble callableSites;

    private final MutableDouble realVariantCount = new MutableDouble(0);
    private final MutableDouble realSNVCount = new MutableDouble(0);
    private final MutableDouble realIndelCount = new MutableDouble(0);
    private final MutableDouble technicalArtifactCount = new MutableDouble(0);

    public SomaticPriorModel(final M2FiltersArgumentCollection MTFAC, final List<MutectStats> mutectStats) {
        log10SNVPrior = MTFAC.log10SNVPrior;
        log10IndelPrior = MTFAC.log10IndelPrior;
        artifactVsVariantPrior = MTFAC.initialPriorOfArtifactVersusVariant;

        callableSites = mutectStats.stream().filter(stat -> stat.getStatistic().equals(Mutect2Engine.CALLABLE_SITES_NAME))
                .mapToDouble(MutectStats::getValue).findFirst();
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
        if (callableSites.isPresent()) {
            log10SNVPrior = Math.log10(realSNVCount.getValue() / callableSites.getAsDouble());
            log10IndelPrior = Math.log10(realIndelCount.getValue() / callableSites.getAsDouble());
        }
    }

    public void clear() {
        realVariantCount.setValue(0);
        realSNVCount.setValue(0);
        realIndelCount.setValue(0);
        technicalArtifactCount.setValue(0);
    }
}
