package org.broadinstitute.hellbender.tools.walkers.mutect.filtering;

import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.lang3.mutable.MutableDouble;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.math3.distribution.BetaDistribution;
import org.apache.commons.math3.distribution.EnumeratedIntegerDistribution;
import org.apache.commons.math3.random.RandomDataGenerator;
import org.broadinstitute.hellbender.tools.walkers.mutect.MutectStats;
import org.broadinstitute.hellbender.tools.walkers.mutect.SomaticLikelihoodsEngine;
import org.broadinstitute.hellbender.tools.walkers.readorientation.BetaDistributionShape;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class SomaticPriorModel {
    private double log10SNVPrior;
    private double log10IndelPrior;
    private double log10NoVariantPrior;
    private double artifactVsVariantPrior;
    private final OptionalDouble callableSites;
    private static final double CONCENTRATION = 0.5;
    private static final int NUM_ITERATIONS = 5;

    final List<Datum> data = new ArrayList<>();

    private final MutableDouble realVariantCount = new MutableDouble(0);
    private final MutableDouble realSNVCount = new MutableDouble(0);
    private final MutableDouble realIndelCount = new MutableDouble(0);
    private final MutableDouble technicalArtifactCount = new MutableDouble(0);

    public SomaticPriorModel(final M2FiltersArgumentCollection MTFAC, final List<MutectStats> mutectStats) {
        log10SNVPrior = MTFAC.log10SNVPrior;
        log10IndelPrior = MTFAC.log10IndelPrior;
        log10NoVariantPrior = MathUtils.log10OneMinusPow10(MathUtils.log10SumLog10(log10SNVPrior, log10IndelPrior));
        artifactVsVariantPrior = MTFAC.initialPriorOfArtifactVersusVariant;

        //callableSites = mutectStats.stream().filter(stat -> stat.getStatistic().equals(Mutect2Engine.CALLABLE_SITES_NAME))
        //        .mapToDouble(MutectStats::getValue).findFirst();
        callableSites = OptionalDouble.empty();
    }

    public double getLog10PriorOfSomaticVariant(final VariantContext vc) {
        return vc.isSNP() ? (MathUtils.LOG10_ONE_THIRD + log10SNVPrior) : log10IndelPrior;
    }

    public double getPriorProbOfArtifactVersusVariant() { return artifactVsVariantPrior; }

    // TODO: original -- build the new and delete this
    public void record(final VariantContext vc, final ErrorProbabilities errorProbabilities) {
        realVariantCount.add(1 - errorProbabilities.getErrorProbability());
        (vc.isSNP() ? realSNVCount : realIndelCount).add(1 - errorProbabilities.getErrorProbability());
        technicalArtifactCount.add(errorProbabilities.getTechnicalArtifactProbability());
    }

    public void recordDatum(final int[] tumorADs, final double[] tumorLog10Odds, final double artifactProbability, final VariantContext.Type type) {
        final int totalAD = (int) MathUtils.sum(tumorADs);
        for (int i = 0; i < tumorLog10Odds.length; i++) {
            data.add(new Datum(tumorLog10Odds[i], artifactProbability, tumorADs[i+1], totalAD, type));
        }

    }

    // by default, clear accumulated data after learning
    public void learnAndClearAccumulatedData() {
        learn();
        data.clear();
        clear();
    }

    public void learn() {
        Utils.resetRandomGenerator();
        final RandomDataGenerator rng = Utils.getRandomDataGenerator();

        Set<AFCluster> clusters = new HashSet<>();

        //initialize flat prior background cluster
        final AFCluster background = AFCluster.makeBackgroundCuster(new BetaDistributionShape(1,1));
        data.stream().limit(100).forEach(background::add);
        clusters.add(background);

        for (int iteration = 0; iteration < NUM_ITERATIONS; iteration++) {
            int[] indices = IntStream.range(0, clusters.size() + 2).toArray();
            for (final Datum datum : data) {
                datum.unassign();
                //stochastically assign to artifact
                if (rng.nextUniform(0,1) < datum.artifactProb) {
                    continue;
                }

                final double log10VariantPrior = datum.getType() == VariantContext.Type.SNP ? log10SNVPrior : log10IndelPrior;

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
                posteriors[clusters.size()] = datum.getTumorLog10Odds() + log10VariantPrior;
                        AFCluster.getNewClusterLog10ChineseRestaurantFactor(CONCENTRATION);

                // note that the log-10 likelihood of no variant is WLOG zero because the tumor LOD is a log likelihood *ratio*
                // we may arbitrarily set the variant likelihood to the LOD and the non-variant likelihood to 0
                posteriors[clusters.size() + 1] = log10NoVariantPrior;

                MathUtils.normalizeLog10(posteriors);   //normalize in-place

                final int index = new EnumeratedIntegerDistribution(rng.getRandomGenerator(), indices, posteriors).sample();
                if (index < clusters.size()) {  // existing cluster
                    clusterLog10RelativePosteriors.get(index).getLeft().add(datum);
                } else if (index == clusters.size()) { // new cluster
                    final double newClusterAlleleFraction = new BetaDistribution(rng.getRandomGenerator(), datum.getAltCount() + 1, datum.getTotalCount() - datum.getAltCount() + 1).sample();
                    final AFCluster newCluster = AFCluster.makeCluster(newClusterAlleleFraction);
                    newCluster.add(datum);
                    clusters.add(newCluster);
                    indices = IntStream.range(0, clusters.size() + 2).toArray();
                }
            }

            clusters = clusters.stream().filter(cluster -> !cluster.isEmpty()).collect(Collectors.toSet());
            clusters.forEach(AFCluster::relearn);
        }





        // global parameters
        artifactVsVariantPrior = (technicalArtifactCount.getValue() + 1) / (realVariantCount.getValue() + technicalArtifactCount.getValue() + 2);
        if (callableSites.isPresent()) {
            log10SNVPrior = Math.log10(realSNVCount.getValue() / callableSites.getAsDouble());
            log10IndelPrior = Math.log10(realIndelCount.getValue() / callableSites.getAsDouble());
        }
    }

    public void clear() {
        realVariantCount.setValue(0);
        realSNVCount.setValue(0);
        realIndelCount.setValue(0);
        technicalArtifactCount.setValue(0);
    }

    private static class Datum {
        private final double tumorLog10Odds;
        private final double artifactProb;
        private final int altCount;
        private final int totalCount;
        private final VariantContext.Type type;

        private AFCluster cluster = null;

        public Datum(final double tumorLog10Odds, final double artifactProb, final int altCount, final int totalCount, final VariantContext.Type type) {
            this.tumorLog10Odds = tumorLog10Odds;
            this.artifactProb = artifactProb;
            this.altCount = altCount;
            this.totalCount = totalCount;
            this.type = type;
        }

        public void unassign() {
            if (cluster != null) {
                cluster.remove(this);
                cluster = null;
            }
        }

        public void assign(final AFCluster cluster) {
            this.cluster = cluster;
            cluster.add(this);
        }

        public double getTumorLog10Odds() { return tumorLog10Odds; }

        public double getArtifactProb() { return artifactProb; }

        public int getAltCount() { return altCount; }

        public int getTotalCount() { return totalCount; }

        public VariantContext.Type getType() { return type; }

        public AFCluster getCluster() { return cluster; }
    }

    // an allele fraction cluster is a beta prior on allele fraction with a method for
    // adjusting a tumor log-likelihood to account for this beta prior as opposed to the flat prior used for Mutect2's
    // initial calls
    private static class AFCluster {
        private BetaDistributionShape betaShape;
        private static double STD_DEV_OVER_MEAN = 0.01;
        private final Set<Datum> members = new HashSet<>();    // the indices of data in this cluster
        private final boolean isBackgroundCluster;

        private static int totalAssignments = 0;

        public static AFCluster makeBackgroundCuster(final BetaDistributionShape betaShape) {
            return new AFCluster(betaShape, true);
        }

        public static AFCluster makeCluster(final double alleleFraction) {
            return new AFCluster(getFuzzyBinomial(alleleFraction, STD_DEV_OVER_MEAN), false);
        }

        private AFCluster(final BetaDistributionShape betaShape, final boolean isBackgroundCluster) {
            this.betaShape = betaShape;
            this.isBackgroundCluster = isBackgroundCluster;
        }

        public void remove(final Datum datum) {
            members.remove(datum);
            totalAssignments--;
        }

        public void add(final Datum datum) {
            members.add(datum);
            totalAssignments++;
        }

        public boolean isEmpty() { return members.isEmpty(); }

        public void relearn() {
            final double altCount = members.stream().mapToDouble(Datum::getAltCount).sum();
            final long totalCount = members.stream().mapToInt(Datum::getTotalCount).sum();
            if (totalCount == 0) {
                members.clear();
                return;
            }
            final double mean = altCount / totalCount;


            if (isBackgroundCluster) {
                //TODO: learn background cluster here



                // TODO: THIS IS EMPTY!!!!!
            } else {
                betaShape = getFuzzyBinomial(mean, STD_DEV_OVER_MEAN);
            }
        }

        public double getLog10ChineseRestaurantFactor(final double concentration) {
            return Math.log10((double) members.size() / (totalAssignments + concentration));
        }

        public static double getNewClusterLog10ChineseRestaurantFactor(final double concentration) {
            return Math.log10(concentration / (totalAssignments + concentration));
        }

        public double log10Likelihood(final Datum datum) {
            final int altCount = datum.getAltCount();
            final int refCount = datum.getTotalCount() - altCount;
            final double alpha = betaShape.getAlpha();
            final double beta = betaShape.getBeta();
            return datum.getTumorLog10Odds() + g(alpha, beta) - g(alpha + altCount, beta + refCount)
                    + g(1,1) - g(1 + altCount, 1 + refCount);
        }

        private static double g(final double... omega) {
            return SomaticLikelihoodsEngine.log10DirichletNormalization(omega);
        }

        private static BetaDistributionShape getFuzzyBinomial(final double mean, final double stdDevOverMean) {
            final double alphaPlusBeta = ((1 - mean) / (mean * MathUtils.square(stdDevOverMean))) - 1;
            final double alpha = mean * alphaPlusBeta;
            final double beta = alphaPlusBeta - alpha;
            return new BetaDistributionShape(alpha, beta);
        }
    }
}
