package org.broadinstitute.hellbender.tools.walkers.mutect.filtering;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.utils.GATKProtectedVariantContextUtils;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.Collections;
import java.util.List;
import java.util.Optional;

public class TumorEvidenceFilter extends Mutect2VariantFilter {
    @Override
    public double calculateArtifactProbability(final VariantContext vc, final Mutect2FilteringInfo filteringInfo) {
        final double[] tumorLods = GATKProtectedVariantContextUtils.getAttributeAsDoubleArray(vc, GATKVCFConstants.TUMOR_LOD_KEY);
        return posteriorProbabilityOfError(MathUtils.arrayMax(tumorLods), filteringInfo.getLog10PriorOfSomaticVariant());
    }

    // a lack of evidence means that the model of independent reads with error rates given by the base qualities explains
    // the apparent variant.  This is not an artifact.
    @Override
    public boolean isTechnicalArtifact() { return false; }

    @Override
    public Optional<String> phredScaledPosteriorAnnotationName() {
        return Optional.of(GATKVCFConstants.SEQUENCING_QUAL_VCF_ATTRIBUTE);
    }

    @Override
    public String filterName() {
        return GATKVCFConstants.TUMOR_LOD_FILTER_NAME;
    }

    @Override
    protected List<String> requiredAnnotations() { return Collections.singletonList(GATKVCFConstants.TUMOR_LOD_KEY); }

}
