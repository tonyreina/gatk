package org.broadinstitute.hellbender.tools.walkers.mutect.filtering;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.util.OverlapDetector;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import org.apache.commons.lang3.mutable.MutableDouble;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.walkers.contamination.ContaminationRecord;
import org.broadinstitute.hellbender.tools.walkers.contamination.MinorAlleleFractionRecord;
import org.broadinstitute.hellbender.tools.walkers.mutect.Mutect2Engine;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.*;
import java.util.stream.Collectors;

/**
 * Simple class to hold all information needed by Mutect2 filters, including: the M2FilteringArgumentCollection,
 * the set of normal samples, the samples' contamination, and the tumor sample CNV segments
 */
public class Mutect2FilteringInfo {
    private static final double FIRST_PASS_THRESHOLD = 0.5;

    private final M2FiltersArgumentCollection MTFAC;
    private final Set<String> normalSamples;
    private final Map<String, Double> contaminationBySample;
    private final Map<String, OverlapDetector<MinorAlleleFractionRecord>> tumorSegments;

    private double artifactProbabilityThreshold = FIRST_PASS_THRESHOLD;

    private double log10PriorOfSomaticVariant;

    // for each PID, the positions with PGTs of filtered genotypes
    private final Map<String, ImmutablePair<Integer, Set<String>>> filteredPhasedCalls;

    // artifact posteriors from first pass
    List<Double> firstPassArtifactProbabilities = new ArrayList<>();

    public Mutect2FilteringInfo(M2FiltersArgumentCollection MTFAC, final VCFHeader vcfHeader) {
        this.MTFAC = MTFAC;
        normalSamples = vcfHeader.getMetaDataInInputOrder().stream()
                .filter(line -> line.getKey().equals(Mutect2Engine.NORMAL_SAMPLE_KEY_IN_VCF_HEADER))
                .map(VCFHeaderLine::getValue)
                .collect(Collectors.toSet());

        contaminationBySample = MTFAC.contaminationTables.stream()
                .map(file -> ContaminationRecord.readFromFile(file).get(0))
                .collect(Collectors.toMap(rec -> rec.getSample(), rec -> rec.getContamination()));

        for (final String sample : vcfHeader.getSampleNamesInOrder()) {
            if (!contaminationBySample.containsKey(sample)) {
                contaminationBySample.put(sample, MTFAC.contaminationEstimate);
            }
        }

        tumorSegments = MTFAC.tumorSegmentationTables.stream()
                .map(MinorAlleleFractionRecord::readFromFile)
                .collect(Collectors.toMap(ImmutablePair::getLeft, p -> OverlapDetector.create(p.getRight())));

        filteredPhasedCalls = new HashMap<>();

        log10PriorOfSomaticVariant = MTFAC.log10PriorProbOfSomaticEvent;
    }

    public M2FiltersArgumentCollection getMTFAC() {
        return MTFAC;
    }

    public Set<String> getNormalSamples() {
        return normalSamples;
    }

    public Map<String, Double> getContaminationBySample() {
        return contaminationBySample;
    }

    public Map<String, OverlapDetector<MinorAlleleFractionRecord>> getTumorSegments() {
        return tumorSegments;
    }

    public Map<String, ImmutablePair<Integer, Set<String>>> getFilteredPhasedCalls() {
        return filteredPhasedCalls;
    }

    public double getArtifactProbabilityThreshold() {
        return artifactProbabilityThreshold;
    }

    public double getLog10PriorOfSomaticVariant() { return log10PriorOfSomaticVariant; }

    public void recordFilteredHaplotypes(final VariantContext vc) {
        final Map<String, Set<String>> phasedGTsForEachPhaseID = vc.getGenotypes().stream()
                .filter(gt -> !normalSamples.contains(gt.getSampleName()))
                .filter(Mutect2FilteringInfo::hasPhaseInfo)
                .collect(Collectors.groupingBy(g -> (String) g.getExtendedAttribute(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_ID_KEY, ""),
                        Collectors.mapping(g -> (String) g.getExtendedAttribute(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_GT_KEY, ""), Collectors.toSet())));

        for (final Map.Entry<String, Set<String>> pidAndPgts : phasedGTsForEachPhaseID.entrySet()) {
            filteredPhasedCalls.put(pidAndPgts.getKey(), new ImmutablePair<>(vc.getStart(), pidAndPgts.getValue()));
        }
    }

    public void addFirstPassArtifactProbability(final double prob) {
        firstPassArtifactProbabilities.add(prob);
    }

    public void adjustThreshold() {
        final double threshold;
        switch (MTFAC.thresholdStrategy) {
            case CONSTANT:
                threshold = MTFAC.posteriorThreshold;
                break;
            case FALSE_DISCOVERY_RATE:
                threshold = calculateThresholdBasedOnFalseDiscoveryRate(firstPassArtifactProbabilities, MTFAC.maxFalsePositiveRate);
                break;
            case OPTIMAL_F_SCORE:
                threshold = calculateThresholdBasedOnOptimalFScore(firstPassArtifactProbabilities, MTFAC.fScoreBeta);
                break;
            default:
                throw new GATKException.ShouldNeverReachHereException("Invalid threshold strategy type: " + MTFAC.thresholdStrategy + ".");
        }

        firstPassArtifactProbabilities.clear();

        artifactProbabilityThreshold = threshold;
    }

    /**
     *
     * Compute the filtering threshold that ensures that the false positive rate among the resulting pass variants
     * will not exceed the requested false positive rate
     *
     * @param posteriors A list of posterior probabilities, which gets sorted
     * @param requestedFPR We set the filtering threshold such that the FPR doesn't exceed this value
     * @return
     */
    @VisibleForTesting
    static double calculateThresholdBasedOnFalseDiscoveryRate(final List<Double> posteriors, final double requestedFPR){
        ParamUtils.isPositiveOrZero(requestedFPR, "requested FPR must be non-negative");
        final double thresholdForFilteringNone = 1.0;
        final double thresholdForFilteringAll = 0.0;

        Collections.sort(posteriors);

        final int numPassingVariants = posteriors.size();
        double cumulativeExpectedFPs = 0.0;

        for (int i = 0; i < numPassingVariants; i++){
            final double posterior = posteriors.get(i);

            // One can show that the cumulative error rate is monotonically increasing in i
            final double expectedFPR = (cumulativeExpectedFPs + posterior) / (i + 1);
            if (expectedFPR > requestedFPR){
                return i > 0 ? posteriors.get(i-1) : thresholdForFilteringAll;
            }

            cumulativeExpectedFPs += posterior;
        }

        // If the expected FP rate never exceeded the max tolerable value, then we can let everything pass
        return thresholdForFilteringNone;
    }

    /**
     * Compute the filtering threshold that maximizes the F_beta score
     *
     * @param posteriors A list of posterior probabilities, which gets sorted
     * @param beta relative weight of recall to precision
     */
    @VisibleForTesting
    static double calculateThresholdBasedOnOptimalFScore(final List<Double> posteriors, final double beta){
        ParamUtils.isPositiveOrZero(beta, "requested F-score beta must be non-negative");
        final double thresholdForFilteringNone = 1.0;
        final double thresholdForFilteringAll = 0.0;

        Collections.sort(posteriors);

        final double expectedTruePositives = posteriors.stream()
                .mapToDouble(prob -> 1 - prob).sum();


        // starting from filtering everything (threshold = 0) increase the threshold to maximize the F score
        final MutableDouble truePositives = new MutableDouble(0);
        final MutableDouble falsePositives = new MutableDouble(0);
        final MutableDouble falseNegatives = new MutableDouble(expectedTruePositives);
        int optimalIndexInclusive = -1; // include all indices up to and including this. -1 mean filter all.
        double optimalFScore = 0;   // if you exclude everything, recall is zero

        final int N = posteriors.size();

        for (int n = 0; n < N; n++){
            truePositives.add(1 - posteriors.get(n));
            falsePositives.add(posteriors.get(n));
            falseNegatives.subtract(1 - posteriors.get(n));
            final double F = (1+beta*beta)*truePositives.getValue() /
                    ((1+beta*beta)*truePositives.getValue() + beta*beta*falseNegatives.getValue() + falsePositives.getValue());
            if (F >= optimalFScore) {
                optimalIndexInclusive = n;
                optimalFScore = F;
            }
        }

        return optimalIndexInclusive == -1 ? 0 : (optimalIndexInclusive == N - 1 ? 1 : posteriors.get(optimalIndexInclusive));
    }

    public static boolean hasPhaseInfo(final Genotype genotype) {
        return genotype.hasExtendedAttribute(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_GT_KEY) && genotype.hasExtendedAttribute(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_ID_KEY);
    }

}
