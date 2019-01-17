package org.broadinstitute.hellbender.tools.walkers.mutect.filtering;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.tools.walkers.annotator.ReadPosition;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.Collections;
import java.util.List;

public class ReadPositionFilter extends HardFilter {
    @Override
    public boolean isArtifact(final VariantContext vc, final Mutect2FilteringInfo filteringInfo) {
        final List<Integer> readPositionByAllele = vc.getAttributeAsIntList(ReadPosition.KEY, 0);

        // a negative value is possible due to a bug: https://github.com/broadinstitute/gatk/issues/5492
        return readPositionByAllele.get(0) > -1 && readPositionByAllele.get(0) < filteringInfo.getMTFAC().minMedianReadPosition;
    }

    @Override
    public String filterName() {
        return GATKVCFConstants.READ_POSITION_FILTER_NAME;
    }

    @Override
    protected List<String> requiredAnnotations() { return Collections.singletonList(ReadPosition.KEY); }
}
