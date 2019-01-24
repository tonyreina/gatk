package org.broadinstitute.hellbender.tools.walkers.mutect.filtering;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.utils.GATKProtectedVariantContextUtils;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.Collections;
import java.util.List;

public class BaseQualityFilter extends HardFilter {
    private final double minMedianBaseQuality;

    public BaseQualityFilter(final double minMedianBaseQuality) {
        this.minMedianBaseQuality = minMedianBaseQuality;
    }

    @Override
    public boolean isArtifact(final VariantContext vc, final Mutect2FilteringEngine filteringInfo) {
        final List<Integer> baseQualityByAllele = vc.getAttributeAsIntList(GATKVCFConstants.MEDIAN_BASE_QUALITY_KEY, 0);
        final double[] tumorLods = GATKProtectedVariantContextUtils.getAttributeAsDoubleArray(vc, GATKVCFConstants.TUMOR_LOD_KEY);
        final int indexOfMaxTumorLod = MathUtils.maxElementIndex(tumorLods);

        return baseQualityByAllele.get(indexOfMaxTumorLod + 1) < minMedianBaseQuality;
    }

    @Override
    public String filterName() {
        return GATKVCFConstants.MEDIAN_BASE_QUALITY_FILTER_NAME;
    }

    @Override
    protected List<String> requiredAnnotations() { return Collections.singletonList(GATKVCFConstants.MEDIAN_BASE_QUALITY_KEY); }
}
