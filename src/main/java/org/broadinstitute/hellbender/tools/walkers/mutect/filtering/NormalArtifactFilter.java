package org.broadinstitute.hellbender.tools.walkers.mutect.filtering;

import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.math3.distribution.BinomialDistribution;
import org.broadinstitute.hellbender.utils.GATKProtectedVariantContextUtils;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.Arrays;
import java.util.List;

public class NormalArtifactFilter extends HardFilter {
    private static final double MIN_NORMAL_ARTIFACT_RATIO = 0.1;    // don't call normal artifact if allele fraction in normal is much smaller than allele fraction in tumor
    private static final int IMPUTED_NORMAL_BASE_QUALITY = 30;  // only used if normal base quality annotation fails somehow

    @Override
    public boolean isArtifact(final VariantContext vc, final Mutect2FilteringInfo filteringInfo) {
        final double[] tumorLods = GATKProtectedVariantContextUtils.getAttributeAsDoubleArray(vc, GATKVCFConstants.TUMOR_LOD_KEY);
        final int indexOfMaxTumorLod = MathUtils.maxElementIndex(tumorLods);

        final int[] tumorAlleleDepths = sumADsOverSamples(vc, filteringInfo.getNormalSamples(), true, false);
        final int tumorDepth = (int) MathUtils.sum(tumorAlleleDepths);
        final int tumorAltDepth = tumorAlleleDepths[indexOfMaxTumorLod + 1];

        final int[] normalAlleleDepths = sumADsOverSamples(vc, filteringInfo.getNormalSamples(), false, true);
        final int normalDepth = (int) MathUtils.sum(normalAlleleDepths);
        final int normalAltDepth = normalAlleleDepths[indexOfMaxTumorLod + 1];

        // if normal AF << tumor AF, don't filter regardless of LOD
        final double tumorAlleleFraction = (double) tumorAltDepth / tumorDepth;
        final double normalAlleleFraction = normalDepth == 0 ? 0 : (double) normalAltDepth / normalDepth;

        if (normalAlleleFraction < MIN_NORMAL_ARTIFACT_RATIO * tumorAlleleFraction)  {
            return false;
        }

        // negative because vcf shows log odds of not artifact / artifact (in order to have bigger positive --> good)
        final double[] normalArtifactLods = GATKProtectedVariantContextUtils.getAttributeAsDoubleArray(vc, GATKVCFConstants.NORMAL_ARTIFACT_LOD_ATTRIBUTE);
        if (-normalArtifactLods[indexOfMaxTumorLod] > filteringInfo.getMTFAC().NORMAL_ARTIFACT_LOD_THRESHOLD) {
            return true;
        }

        // the above filter misses artifacts whose support in the normal consists entirely of low base quality reads
        // Since a lot of low-BQ reads is itself evidence of an artifact, we filter these by hand via an estimated LOD
        // that uses the average base quality of *ref* reads in the normal
        final int medianRefBaseQuality = vc.getAttributeAsIntList(GATKVCFConstants.MEDIAN_BASE_QUALITY_KEY, IMPUTED_NORMAL_BASE_QUALITY).get(0);
        final double normalPValue = 1 - new BinomialDistribution(null, normalDepth, QualityUtils.qualToErrorProb(medianRefBaseQuality))
                .cumulativeProbability(normalAltDepth - 1);

        return normalPValue < M2FiltersArgumentCollection.normalPileupPValueThreshold;
    }

    public String filterName() {
        return GATKVCFConstants.ARTIFACT_IN_NORMAL_FILTER_NAME;
    }

    protected List<String> requiredAnnotations() {
        return Arrays.asList(GATKVCFConstants.NORMAL_ARTIFACT_LOD_ATTRIBUTE, GATKVCFConstants.TUMOR_LOD_KEY);
    }
}
