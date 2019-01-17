package org.broadinstitute.hellbender.tools.walkers.mutect.filtering;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.Collections;
import java.util.List;

public class LogOddsOverDepthFilter extends HardFilter {
    @Override
    public boolean isArtifact(final VariantContext vc, final Mutect2FilteringInfo filteringInfo) {
        if(!vc.isBiallelic()) {
            return false;
        }

        final Double lod = vc.getAttributeAsDouble(GATKVCFConstants.TUMOR_LOD_KEY, 1);
        final Double depth = vc.getAttributeAsDouble(VCFConstants.DEPTH_KEY, 1);
        return lod / depth < filteringInfo.getMTFAC().minLog10OddsDividedByDepth;
    }

    public String filterName() {
        return GATKVCFConstants.LOW_AVG_ALT_QUALITY_FILTER_NAME;
    }

    protected List<String> requiredAnnotations() { return Collections.singletonList(GATKVCFConstants.TUMOR_LOD_KEY); }
}
