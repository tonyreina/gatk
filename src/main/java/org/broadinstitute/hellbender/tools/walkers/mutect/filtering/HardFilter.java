package org.broadinstitute.hellbender.tools.walkers.mutect.filtering;

import htsjdk.variant.variantcontext.VariantContext;

import java.util.Optional;

public abstract class HardFilter extends Mutect2VariantFilter {
    @Override
    public double calculateArtifactProbability(final VariantContext vc, final Mutect2FilteringEngine filteringInfo) {
        return isArtifact(vc, filteringInfo) ? 1 : 0;
    }

    public abstract boolean isArtifact(final VariantContext vc, final Mutect2FilteringEngine filteringInfo);

    // the posterior of a hard filter is 0 or 1, hence there's no reason to annotate it
    @Override
    public Optional<String> phredScaledPosteriorAnnotationName() {
        return Optional.empty();
    }
}
