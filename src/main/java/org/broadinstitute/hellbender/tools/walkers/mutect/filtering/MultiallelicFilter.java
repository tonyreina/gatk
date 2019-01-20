package org.broadinstitute.hellbender.tools.walkers.mutect.filtering;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.utils.GATKProtectedVariantContextUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public class MultiallelicFilter extends HardFilter {
    @Override
    public boolean isArtifact(final VariantContext vc, final Mutect2FilteringInfo filteringInfo) {
        final double[] tumorLods = GATKProtectedVariantContextUtils.getAttributeAsDoubleArray(vc, GATKVCFConstants.TUMOR_LOD_KEY);

        final long numPassingAltAlleles = Arrays.stream(tumorLods)
                .map(lod -> posteriorProbabilityOfError(lod, filteringInfo.getLog10PriorOfSomaticVariant(vc)))
                .filter(prob -> prob < filteringInfo.getArtifactProbabilityThreshold())
                .count();

        return numPassingAltAlleles > filteringInfo.getMTFAC().numAltAllelesThreshold;
    }

    public String filterName() {
        return GATKVCFConstants.MULTIALLELIC_FILTER_NAME;
    }

    protected List<String> requiredAnnotations() { return Collections.singletonList(GATKVCFConstants.TUMOR_LOD_KEY); }
}
