package org.broadinstitute.hellbender.tools.walkers.mutect.filtering;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.lang3.mutable.MutableInt;
import org.broadinstitute.hellbender.tools.walkers.annotator.StrandBiasBySample;
import org.broadinstitute.hellbender.utils.GATKProtectedVariantContextUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.Collections;
import java.util.List;

public class StrictStrandBiasFilter extends HardFilter {
    @Override
    public boolean isArtifact(final VariantContext vc, final Mutect2FilteringInfo filteringInfo) {
        if (filteringInfo.getMTFAC().minReadsOnEachStrand == 0) {
            return false;
        }

        final MutableInt altForwardCount = new MutableInt(0);
        final MutableInt altReverseCount = new MutableInt(0);

        for (final Genotype g : vc.getGenotypes()) {
            if (filteringInfo.getNormalSamples().contains(g.getSampleName())) {
                continue;
            } else if (!g.hasExtendedAttribute(GATKVCFConstants.STRAND_BIAS_BY_SAMPLE_KEY)) {
                return false;
            } else {
                final int[] strandBiasCounts = GATKProtectedVariantContextUtils.getAttributeAsIntArray(g, GATKVCFConstants.STRAND_BIAS_BY_SAMPLE_KEY, () -> null, 0);
                altForwardCount.add(StrandBiasBySample.getAltForwardCountFromFlattenedContingencyTable(strandBiasCounts));
                altReverseCount.add(StrandBiasBySample.getAltReverseCountFromFlattenedContingencyTable(strandBiasCounts));
            }
        }

        // filter if there is no alt evidence in the forward or reverse strand
        return Math.min(altForwardCount.getValue(), altReverseCount.getValue()) >= filteringInfo.getMTFAC().minReadsOnEachStrand;
    }

    @Override
    public String filterName() {
        return GATKVCFConstants.STRICT_STRAND_BIAS_FILTER_NAME;
    }

    @Override
    protected List<String> requiredAnnotations() { return Collections.emptyList(); }
}
