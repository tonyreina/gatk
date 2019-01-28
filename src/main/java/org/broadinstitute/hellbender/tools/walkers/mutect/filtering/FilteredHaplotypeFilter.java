package org.broadinstitute.hellbender.tools.walkers.mutect.filtering;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.utils.GATKProtectedVariantContextUtils;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.*;

public class FilteredHaplotypeFilter extends Mutect2VariantFilter {
    private final double maxIntraHaplotypeDistance;

    // for each pgt + pid phasing string, a list of loci-error probability pairs
    private final Map<String, List<Pair<Integer, Double>>> phasedProbabilities = new HashMap<>();

    public FilteredHaplotypeFilter(final double maxIntraHaplotypeDistance) {
        this.maxIntraHaplotypeDistance = maxIntraHaplotypeDistance;
    }

    @Override
    public ErrorType errorType() { return ErrorType.ARTIFACT; }

    @Override
    public double calculateErrorProbability(final VariantContext vc, final Mutect2FilteringEngine filteringInfo) {
        // use phasing of tumor genotype with greatest allele fraction
        final Genotype tumorGenotype = vc.getGenotypes().stream().filter(filteringInfo::isTumor)
                .max(Comparator.comparingDouble(g -> MathUtils.arrayMax(GATKProtectedVariantContextUtils.getAttributeAsDoubleArray(g, VCFConstants.ALLELE_FREQUENCY_KEY,
                        () -> new double[] {0.0}, 0.0)))).get();

        final Optional<String> phasingString = makePhasingString(tumorGenotype);
        if (!phasingString.isPresent()) {
            return 0.0;
        }

        final List<Pair<Integer, Double>> phasedProbs = phasedProbabilities.get(phasingString.get());

        if (phasedProbs == null) {
            return 0.0;
        }

        return phasedProbs.stream()
                .filter(pair -> Math.abs(pair.getLeft() - vc.getStart()) <= maxIntraHaplotypeDistance)
                .mapToDouble(Pair::getRight)
                .max().orElse(0.0);
    }

    @Override
    protected void accumulateDataForLearning(final VariantContext vc, final ErrorProbabilities errorProbabilities, final Mutect2FilteringEngine filteringInfo) {
        final double technicalArtifactProbability = errorProbabilities.getTechnicalArtifactProbability();

        for (final Genotype tumorGenotype : vc.getGenotypes()) {
            if (!filteringInfo.isTumor(tumorGenotype)) {
                continue;
            }

            final Optional<String> phasingString = makePhasingString(tumorGenotype);

            if (!phasingString.isPresent()) {
                continue;
            }

            if (!phasedProbabilities.containsKey(phasingString.get())) {
                phasedProbabilities.put(phasingString.get(), new ArrayList<>());
            }

            phasedProbabilities.get(phasingString.get()).add(ImmutablePair.of(vc.getStart(), technicalArtifactProbability));
        }
    }

    // although we don't have to do so, it's worth explicitly overriding to make clear that although there is a
    // non-trivial accumulate method, nonetheless we don't do anything in the clear and learn methods
    @Override
    protected void clearAccumulatedData() { }

    @Override
    protected void learnParameters() { }

    @Override
    public String filterName() {
        return GATKVCFConstants.BAD_HAPLOTYPE_FILTER_NAME;
    }

    @Override
    protected List<String> requiredAnnotations() { return Collections.emptyList(); }

    @Override
    public Optional<String> phredScaledPosteriorAnnotationName() { return Optional.empty(); }

    private static boolean hasPhaseInfo(final Genotype genotype) {
        return genotype.hasExtendedAttribute(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_GT_KEY) && genotype.hasExtendedAttribute(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_ID_KEY);
    }

    // concatenate the PGT and PID strings, if present
    private static Optional<String> makePhasingString(final Genotype genotype) {
        final String pgt = (String) genotype.getExtendedAttribute(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_GT_KEY, null);
        final String pid = (String) genotype.getExtendedAttribute(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_ID_KEY, null);
        return (pgt == null || pid == null) ? Optional.empty() : Optional.of(pgt + pid);
    }
}
