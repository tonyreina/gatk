package org.broadinstitute.hellbender.tools.walkers.mutect.clustering;

import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.lang3.mutable.MutableInt;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.math3.distribution.BetaDistribution;
import org.apache.commons.math3.distribution.EnumeratedIntegerDistribution;
import org.apache.commons.math3.random.RandomDataGenerator;
import org.broadinstitute.hellbender.tools.walkers.mutect.Mutect2Engine;
import org.broadinstitute.hellbender.tools.walkers.mutect.MutectStats;
import org.broadinstitute.hellbender.tools.walkers.mutect.SomaticLikelihoodsEngine;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.M2FiltersArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.readorientation.BetaDistributionShape;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class SomaticClusteringModel {
    //TODO: might need to move
    public static final BetaDistributionShape FLAT_BETA = new BetaDistributionShape(1,1);

    private double log10SNVPrior;
    private double log10IndelPrior;
    private double log10NoVariantPrior;
    private double log10VariantVsArtifactPrior;
    private final OptionalDouble callableSites;

    private static final double CONCENTRATION = 0.5;
    private static final int NUM_ITERATIONS = 5;
    final List<Datum> data = new ArrayList<>();

    List<Pair<Double, BetaDistributionShape>> alleleFractionClusters;

    public SomaticClusteringModel(final M2FiltersArgumentCollection MTFAC, final List<MutectStats> mutectStats) {
        log10SNVPrior = MTFAC.log10SNVPrior;
        log10IndelPrior = MTFAC.log10IndelPrior;
        log10NoVariantPrior = MathUtils.log10OneMinusPow10(MathUtils.log10SumLog10(log10SNVPrior, log10IndelPrior));
        log10VariantVsArtifactPrior = MTFAC.initialLog10PriorOfVariantVersusArtifact;
        callableSites = mutectStats.stream().filter(stat -> stat.getStatistic().equals(Mutect2Engine.CALLABLE_SITES_NAME))
                .mapToDouble(MutectStats::getValue).findFirst();

        // initialize clusters so that all weight (that is log-10(weight = 1) = 0) is in the flat beta prior used in the Mutect2 tumor LOD
        alleleFractionClusters = new ArrayList<>();
        alleleFractionClusters.add(ImmutablePair.of(0.0, FLAT_BETA));
    }

    public double getLog10PriorOfSomaticVariant(final VariantContext vc) {
        return vc.isSNP() ? (MathUtils.LOG10_ONE_THIRD + log10SNVPrior) : log10IndelPrior;
    }

    public double getLog10PriorProbOfVariantVersusArtifact() { return log10VariantVsArtifactPrior; }

    /**
     * Correct the tumor log-10 odds (ie the TLOD) from Mutect2 to account for allele fraction clustering
     * @param tumorLog10Odds the original tumor log-10 odds of Mutect2 that assumes a flat prior on allele fractions
     * @return the log-10 odds corrected for the allele fraction clustering learned by this class
     */
    public double clusteringCorrectedLog10Odds(final double tumorLog10Odds, final int altCount, final int refCount) {
        return tumorLog10Odds + MathUtils.log10SumLog10(alleleFractionClusters.stream()
                .mapToDouble(pair -> pair.getLeft() + log10OddsCorrection(FLAT_BETA, pair.getRight(), altCount, refCount))
                .toArray());
    }

    public void record(final int[] tumorADs, final double[] tumorLog10Odds, final double artifactProbability, final double nonSomaticProbability, final VariantContext.Type type) {
        final int totalAD = (int) MathUtils.sum(tumorADs);
        // split into one-vs-all biallelics for clustering
        for (int i = 0; i < tumorLog10Odds.length; i++) {
            data.add(new Datum(tumorLog10Odds[i], artifactProbability, nonSomaticProbability, tumorADs[i+1], totalAD, type));
        }
    }

    public void learnAndClearAccumulatedData() {
        Utils.resetRandomGenerator();
        final RandomDataGenerator rng = Utils.getRandomDataGenerator();
        AFCluster.clearAssignmentCount();

        Set<AFCluster> clusters = new HashSet<>();

        //initialize flat prior background cluster
        final AFCluster background = AFCluster.makeBackgroundCuster(FLAT_BETA);
        data.stream().limit(1).forEach(background::add);
        clusters.add(background);

        for (int iteration = 0; iteration < NUM_ITERATIONS; iteration++) {
            int k = 0;
            int[] indices = IntStream.range(0, clusters.size() + 2).toArray();
            for (final Datum datum : data) {
                k++;
                datum.unassign();
                //stochastically assign to error (other than sequencing error, which is handled within the clustering)
                if (rng.nextUniform(0,1) < datum.getNonSequencingErrorProb()) {
                    continue;
                }

                final double log10VariantPrior = datum.getType() == VariantContext.Type.SNP ? (MathUtils.LOG10_ONE_THIRD + log10SNVPrior) : log10IndelPrior;

                // the posterior contains a prior that it's a variant, a CRP factor within the set of variant clusters,
                // and the cluster likelihood
                final List<Pair<AFCluster, Double>> clusterLog10RelativePosteriors = clusters.stream().map(cluster -> {
                    final double posterior = log10VariantPrior + cluster.log10Likelihood(datum) + cluster.getLog10ChineseRestaurantFactor(CONCENTRATION);
                        return ImmutablePair.of(cluster, posterior);
                }).collect(Collectors.toList());

                // we'll make a list of all relative posteriors -- the existing clusters first, then a hypothetical new cluster,
                // then sequencing error -- and sample from the normalized posteriors
                final double[] posteriors = new double[clusters.size() + 2];
                for (int i = 0; i < clusters.size(); i++) {
                    posteriors[i] = clusterLog10RelativePosteriors.get(i).getValue();
                }
                // the new cluster likelihood integrates over all AFs with a flat prior, which is exactly the original Mutect2 tumor LOD!
                posteriors[clusters.size()] = datum.getTumorLog10Odds() + log10VariantPrior +
                        AFCluster.getNewClusterLog10ChineseRestaurantFactor(CONCENTRATION);

                // note that the log-10 likelihood of no variant is WLOG zero because the tumor LOD is a log likelihood *ratio*
                // we may arbitrarily set the variant likelihood to the LOD and the non-variant likelihood to 0
                posteriors[clusters.size() + 1] = log10NoVariantPrior;

                MathUtils.normalizeLog10(posteriors, false, true);   //normalize in-place

                final int index = new EnumeratedIntegerDistribution(rng.getRandomGenerator(), indices, posteriors).sample();
                if (index < clusters.size()) {  // existing cluster
                    datum.assign(clusterLog10RelativePosteriors.get(index).getLeft());
                } else if (index == clusters.size()) { // new cluster
                    final double newClusterAlleleFraction = new BetaDistribution(rng.getRandomGenerator(), datum.getAltCount() + 1, datum.getTotalCount() - datum.getAltCount() + 1).sample();
                    final AFCluster newCluster = AFCluster.makeCluster(newClusterAlleleFraction);
                    datum.assign(newCluster);
                    clusters.add(newCluster);
                    indices = IntStream.range(0, clusters.size() + 2).toArray();
                }
            }

            clusters = clusters.stream().filter(cluster -> !cluster.isEmpty()).collect(Collectors.toSet());
            clusters.forEach(AFCluster::relearn);

            final double technicalArtifactCount = data.stream().mapToDouble(Datum::getArtifactProb).sum();
            //final long realSNVCount = clusters.stream().mapToInt(AFCluster::SNVCount).sum();
            //final long realIndelCount = clusters.stream().mapToInt(AFCluster::size).sum() - realSNVCount;

            final double realSNVCount = data.stream().filter(data1 -> data1.getType() == VariantContext.Type.SNP).mapToDouble(datum -> {
                final double tumorLog10Odds = clusteringCorrectedLog10Odds(datum.getTumorLog10Odds(), datum.getAltCount(), datum.getTotalCount() - datum.getAltCount());
                final double[] variantVsSequencingErrorProbs = MathUtils.normalizeFromLog10ToLinearSpace(new double[] { log10SNVPrior + tumorLog10Odds, log10NoVariantPrior});
                return (1 - datum.getNonSequencingErrorProb()) * variantVsSequencingErrorProbs[0];
            }).sum();

            final double realIndelCount = data.stream().filter(data1 -> data1.getType() != VariantContext.Type.SNP).mapToDouble(datum -> {
                final double tumorLog10Odds = clusteringCorrectedLog10Odds(datum.getTumorLog10Odds(), datum.getAltCount(), datum.getTotalCount() - datum.getAltCount());
                final double[] variantVsSequencingErrorProbs = MathUtils.normalizeFromLog10ToLinearSpace(new double[] { log10IndelPrior + tumorLog10Odds, log10NoVariantPrior});
                return (1 - datum.getNonSequencingErrorProb()) * variantVsSequencingErrorProbs[0];
            }).sum();

            if (callableSites.isPresent()) {
                log10SNVPrior = Math.log10(Math.max(realSNVCount / callableSites.getAsDouble(), 1.0e-8));
                log10IndelPrior = Math.log10(Math.max(realIndelCount / callableSites.getAsDouble(), 1.0e-8));
                log10NoVariantPrior = MathUtils.log10OneMinusPow10(MathUtils.log10SumLog10(log10SNVPrior, log10IndelPrior));
            }
            log10VariantVsArtifactPrior = Math.log10((realSNVCount + realIndelCount + 1) / (realSNVCount + realIndelCount + technicalArtifactCount + 2));
        }
        alleleFractionClusters = clusters.stream()
                .map(cluster -> ImmutablePair.of(cluster.getLog10ChineseRestaurantFactor(CONCENTRATION), cluster.betaShape))
                .collect(Collectors.toList());
        data.clear();
    }

    public List<Pair<String, String>> clusteringMetadata() {
        final List<Pair<String, String>> result = new ArrayList<>();
        result.add(ImmutablePair.of("log10 SNV prior", Double.toString(log10SNVPrior)));
        result.add(ImmutablePair.of("log10 Indel prior", Double.toString(log10IndelPrior)));

        final MutableInt clusterIndex = new MutableInt(1);
        alleleFractionClusters.stream().sorted(Comparator.comparingDouble(pair-> -pair.getLeft())).forEach(log10WeightAndShape -> {
            final double weight = Math.pow(10, log10WeightAndShape.getLeft());
            final double alpha = log10WeightAndShape.getRight().getAlpha();
            final double beta = log10WeightAndShape.getRight().getBeta();
            result.add(ImmutablePair.of("cluster " + clusterIndex.getValue().toString(),
                    String.format("weight = %.3f, alpha = %.2f, beta = %.2f", weight, alpha, beta)));
            clusterIndex.increment();
        });

        return result;
    }

    // TODO: should be abe to move stuff and make private
    public static double log10OddsCorrection(final BetaDistributionShape originalBeta, final BetaDistributionShape newBeta, final int altCount, final int refCount) {
        return g(newBeta.getAlpha(), newBeta.getBeta()) - g(newBeta.getAlpha() + altCount, newBeta.getBeta() + refCount)
                - g(originalBeta.getAlpha(), originalBeta.getBeta()) + g(originalBeta.getAlpha() + altCount, originalBeta.getBeta() + refCount);
    }

    protected static double g(final double... omega) {
        return SomaticLikelihoodsEngine.log10DirichletNormalization(omega);
    }
}
