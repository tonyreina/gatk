package org.broadinstitute.hellbender.tools.walkers.mutect.filtering;

import com.google.common.annotations.VisibleForTesting;
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
        return posteriorProbabilityOfError(MathUtils.arrayMax(tumorLods), filteringInfo.getMTFAC().log10PriorProbOfSomaticEvent);
    }

    public Optional<String> phredScaledPosteriorAnnotationName() {
        return Optional.of(GATKVCFConstants.SEQUENCING_QUAL_VCF_ATTRIBUTE);
    }

    public String filterName() {
        return GATKVCFConstants.TUMOR_LOD_FILTER_NAME;
    }

    protected List<String> requiredAnnotations() { return Collections.singletonList(GATKVCFConstants.TUMOR_LOD_KEY); }

    @VisibleForTesting
    static double posteriorProbabilityOfError(final double log10OddsOfRealVersusError, final double log10PriorOfReal) {
        final double[] unweightedPosteriorOfRealAndError = new double[] {log10OddsOfRealVersusError + log10PriorOfReal,
                MathUtils.log10OneMinusPow10(log10PriorOfReal)};

        final double[] posteriorOfRealAndError = MathUtils.normalizeFromLog10ToLinearSpace(unweightedPosteriorOfRealAndError);

        return posteriorOfRealAndError[1];
    }

}
