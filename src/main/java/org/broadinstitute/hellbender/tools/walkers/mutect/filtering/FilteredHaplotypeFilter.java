package org.broadinstitute.hellbender.tools.walkers.mutect.filtering;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.utils.GATKProtectedVariantContextUtils;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Set;

public class FilteredHaplotypeFilter extends HardFilter {
    private final double maxIntraHaplotypeDistance;

    public FilteredHaplotypeFilter(final double maxIntraHaplotypeDistance) {
        this.maxIntraHaplotypeDistance = maxIntraHaplotypeDistance;
    }

    @Override
    public ErrorType errorType() { return ErrorType.ARTIFACT; }

    @Override
    public boolean isArtifact(final VariantContext vc, final Mutect2FilteringEngine filteringInfo) {
        // use phasing of tumor genotype with greatest allele fraction
        final Genotype tumorGenotype = vc.getGenotypes().stream().filter(filteringInfo::isTumor)
                .max(Comparator.comparingDouble(g -> MathUtils.arrayMax(GATKProtectedVariantContextUtils.getAttributeAsDoubleArray(g, VCFConstants.ALLELE_FREQUENCY_KEY,
                        () -> new double[] {0.0}, 0.0)))).get();

        if (!Mutect2FilteringEngine.hasPhaseInfo(tumorGenotype)) {
            return false;
        }

        final String pgt = (String) tumorGenotype.getExtendedAttribute(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_GT_KEY, "");
        final String pid = (String) tumorGenotype.getExtendedAttribute(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_ID_KEY, "");
        final int position = vc.getStart();

        final Pair<Integer, Set<String>> filteredCall = filteringInfo.getFilteredPhasedCalls().get(pid);
        if (filteredCall == null) {
            return false;
        }

        // Check that vc occurs on a filtered haplotype
        return filteredCall.getRight().contains(pgt) && Math.abs(filteredCall.getLeft() - position) <= maxIntraHaplotypeDistance;
    }

    @Override
    public String filterName() {
        return GATKVCFConstants.BAD_HAPLOTYPE_FILTER_NAME;
    }

    @Override
    protected List<String> requiredAnnotations() { return Collections.emptyList(); }
}
