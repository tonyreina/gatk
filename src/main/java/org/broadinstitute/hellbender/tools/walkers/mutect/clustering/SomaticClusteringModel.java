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
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.M2FiltersArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.readorientation.BetaDistributionShape;
import org.broadinstitute.hellbender.utils.IndexRange;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class SomaticClusteringModel {

    private double log10SNVPrior;
    private double log10IndelPrior;
    private double log10NoVariantPrior;
    private double log10VariantVsArtifactPrior;
    private final OptionalDouble callableSites;

    private static final double INITIAL_HIGH_AF_WEIGHT = 0.01;
    private static final double INITIAL_BACKGROUND_WEIGHT = 0.1;

    private double log10HighAFWeight = Math.log10(INITIAL_HIGH_AF_WEIGHT);
    private double log10BackgroundWeight = Math.log10(INITIAL_BACKGROUND_WEIGHT);
    private double log10SparseClustersWeight = MathUtils.log10OneMinusX(MathUtils.log10SumLog10(log10HighAFWeight, log10BackgroundWeight));

    private static final double CONCENTRATION = 0.5;
    private static final int NUM_ITERATIONS = 5;

    private static final int SEQUENCING_ERROR_INDEX = 0;
    private static final int HIGH_AF_INDEX = 1;
    private static final int BACKGROUND_INDEX = 2;
    private static final int OFFSET = 3;

    private static final BetaDistributionShape INITIAL_HIGH_AF_BETA = new BetaDistributionShape(10, 1);
    private static final BetaDistributionShape INITIAL_BACKGROUND_BETA = BetaDistributionShape.FLAT_BETA;
    private static final AlleleFractionCluster NEW_CLUSTER = new BetaBinomialCluster(BetaDistributionShape.FLAT_BETA);

    private final RandomDataGenerator rng = Utils.getRandomDataGenerator();


    final List<Datum> data = new ArrayList<>();
    List<AlleleFractionCluster> clusters;

    List<OptionalInt> clusterAssignments = new ArrayList<>();
    final List<MutableInt> clusterCounts = new ArrayList<>();
    final MutableInt totalSparseClusterCount = new MutableInt(0);

    final MutableInt snvCount = new MutableInt(0);
    final MutableInt indelCount = new MutableInt(0);

    private boolean firstPass = true;

    public SomaticClusteringModel(final M2FiltersArgumentCollection MTFAC, final List<MutectStats> mutectStats) {
        log10SNVPrior = MTFAC.log10SNVPrior;
        log10IndelPrior = MTFAC.log10IndelPrior;
        log10NoVariantPrior = MathUtils.log10OneMinusPow10(MathUtils.log10SumLog10(log10SNVPrior, log10IndelPrior));
        log10VariantVsArtifactPrior = MTFAC.initialLog10PriorOfVariantVersusArtifact;
        callableSites = mutectStats.stream().filter(stat -> stat.getStatistic().equals(Mutect2Engine.CALLABLE_SITES_NAME))
                .mapToDouble(MutectStats::getValue).findFirst();

        clusters = new ArrayList<>();
        clusters.add(SEQUENCING_ERROR_INDEX, new SequencingError());
        clusters.add(HIGH_AF_INDEX, new BetaBinomialCluster(INITIAL_HIGH_AF_BETA));
        clusters.add(BACKGROUND_INDEX, new BetaBinomialCluster(INITIAL_BACKGROUND_BETA));
    }

    public double getLog10PriorOfSomaticVariant(final VariantContext vc) {
        return vc.isSNP() ? (MathUtils.LOG10_ONE_THIRD + log10SNVPrior) : log10IndelPrior;
    }

    public double getLog10PriorProbOfVariantVersusArtifact() { return log10VariantVsArtifactPrior; }

    public double probabilityOfSequencingError(final Datum datum) {
        return clusterProbabilities(datum)[SEQUENCING_ERROR_INDEX];
    }

    public void record(final int[] tumorADs, final double[] tumorLog10Odds, final double artifactProbability, final double nonSomaticProbability, final VariantContext.Type type) {
        final int totalAD = (int) MathUtils.sum(tumorADs);
        // split into one-vs-all biallelics for clustering
        for (int i = 0; i < tumorLog10Odds.length; i++) {
            data.add(new Datum(tumorLog10Odds[i], artifactProbability, nonSomaticProbability, tumorADs[i+1], totalAD, type));
        }
    }

    public void learnAndClearAccumulatedData() {
        if (firstPass) {
            clusterAssignments.addAll(Collections.nCopies(data.size(), OptionalInt.empty()));
            clusterCounts.addAll(Collections.nCopies(clusters.size(), new MutableInt(0)));
        }

        for (int iteration = 0; iteration < NUM_ITERATIONS; iteration++) {
            for (int datumIndex = 0; datumIndex < data.size(); datumIndex++) {
                final Datum datum = popDatum(datumIndex);

                //stochastically assign to non-sequencing error (sequencing error is handled within the clustering)
                if (rng.nextUniform(0,1) < datum.getNonSequencingErrorProb()) {
                    continue;
                }

                final double[] clusterPosteriors = clusterProbabilities(datum);

                final int[] indices = IntStream.range(0, clusters.size() + 1).toArray();
                final int clusterIndex = new EnumeratedIntegerDistribution(rng.getRandomGenerator(), indices, clusterPosteriors).sample();
                assignDatum(datumIndex, clusterIndex);

            }

            pruneEmptyClusters();

            final List<List<Datum>> dataByCluster = clusters.stream().map(c -> new ArrayList<Datum>()).collect(Collectors.toList());
            for (final MutableInt datumIndex = new MutableInt(0); datumIndex.getValue() < clusterAssignments.size(); datumIndex.increment()) {
                clusterAssignments.get(datumIndex.getValue()).ifPresent(c -> dataByCluster.get(c).add(data.get(datumIndex.getValue())));
            }

            for (int clusterIndex = 0; clusterIndex < clusters.size(); clusterIndex++) {
                clusters.get(clusterIndex).learn(dataByCluster.get(clusterIndex));
            }

            learnWeightsAndPriors();
        }

        firstPass = false;
        data.clear();
    }

    private void pruneEmptyClusters() {
        final Map<Integer, Integer> oldToNewClusterIndices = new TreeMap<>();

        int newIndex = OFFSET;
        for (int oldIndex = OFFSET; oldIndex < clusters.size(); oldIndex++) {
            if (clusterCounts.get(oldIndex).getValue() > 0) {
                oldToNewClusterIndices.put(oldIndex, newIndex);

                if (newIndex != oldIndex) {
                    clusters.set(newIndex, clusters.get(oldIndex));
                    clusterCounts.set(newIndex, clusterCounts.get(oldIndex));
                }
                newIndex++;
            }
        }

        clusterAssignments = clusterAssignments.stream()
                .map(a -> a.isPresent() ? OptionalInt.of(oldToNewClusterIndices.get(a)) : a)
                .collect(Collectors.toList());
    }

    private double[] clusterProbabilities(final Datum datum) {
        final double log10VariantPrior = datum.getType() == VariantContext.Type.SNP ? (MathUtils.LOG10_ONE_THIRD + log10SNVPrior) : log10IndelPrior;

        final double[] log10ClusterPosteriors = new IndexRange(0, clusters.size() + 1).mapToDouble(c -> {
            final double log10Likelihood = c < clusters.size() ? clusters.get(c).log10Likelihood(datum) :
                    NEW_CLUSTER.log10Likelihood(datum);
            if (c == SEQUENCING_ERROR_INDEX) {
                return log10NoVariantPrior + log10Likelihood;
            } else if (c == HIGH_AF_INDEX) {
                return log10VariantPrior + log10HighAFWeight + log10Likelihood;
            } else if (c == BACKGROUND_INDEX) {
                return log10VariantPrior + log10BackgroundWeight + log10Likelihood;
            } else if (c < clusters.size()) {   // existing sparse cluster
                return log10VariantPrior + log10SparseClustersWeight + log10CRPWeight(c)
                        + log10Likelihood;
            } else {    // new sparse cluster
                return log10VariantPrior + log10SparseClustersWeight + log10CRPWeight(c)
                        + log10Likelihood;
            }
        });

        return MathUtils.normalizeLog10(log10ClusterPosteriors, false, false);
    }

    private double log10CRPWeight(final int clusterIndex) {
        Utils.validate(clusterIndex >= OFFSET, "Chinese restaurant process does not apply to error, high-AF, and backgorund clusters");
        final double numerator = clusterIndex == clusters.size() ? CONCENTRATION : clusterCounts.get(clusterIndex).getValue();
        final double denominator = totalSparseClusterCount.getValue() + CONCENTRATION;
        return Math.log10(numerator / denominator);

    }

    private Datum popDatum(final int datumIndex) {
        // pop datum off its cluster and decrement the sparse cluster count if appropriate
        clusterAssignments.get(datumIndex).ifPresent(c -> {
            clusterCounts.get(c).decrement();
            if (OFFSET <= c) {
                totalSparseClusterCount.decrement();
            }
            if (c != SEQUENCING_ERROR_INDEX) {
                (data.get(datumIndex).getType() == VariantContext.Type.SNP ? snvCount : indelCount).decrement();
            }
        });
        clusterAssignments.set(datumIndex, OptionalInt.empty());
        return data.get(datumIndex);
    }

    private void assignDatum(int datumIndex, int clusterIndex) {
        final Datum datum = data.get(datumIndex);

        // make new cluster
        if (clusterIndex == clusters.size()) {
            final double newClusterAlleleFraction = new BetaDistribution(rng.getRandomGenerator(), datum.getAltCount() + 1, datum.getTotalCount() - datum.getAltCount() + 1).sample();
            clusters.add(new BinomialCluster(newClusterAlleleFraction));
        }

        if (OFFSET <= clusterIndex) {
            totalSparseClusterCount.increment();
        }

        if (clusterIndex != SEQUENCING_ERROR_INDEX) {
            (datum.getType() == VariantContext.Type.SNP ? snvCount : indelCount).increment();
        }

        clusterAssignments.set(datumIndex, OptionalInt.of(clusterIndex));
        clusterCounts.get(clusterIndex).increment();
    }

    private void learnWeightsAndPriors() {
        final double totalVariants = clusterCounts.get(HIGH_AF_INDEX).getValue() + clusterCounts.get(BACKGROUND_INDEX).getValue()
                + totalSparseClusterCount.getValue() + 0.0001;

        log10HighAFWeight = Math.log10((double) clusterCounts.get(HIGH_AF_INDEX).getValue() / totalVariants);
        log10BackgroundWeight = Math.log10((double) clusterCounts.get(BACKGROUND_INDEX).getValue() / totalVariants);
        log10SparseClustersWeight = MathUtils.log10OneMinusX(MathUtils.log10SumLog10(log10HighAFWeight, log10BackgroundWeight));

        final double technicalArtifactCount = data.stream().mapToDouble(Datum::getArtifactProb).sum();

        if (callableSites.isPresent()) {
            log10SNVPrior = Math.log10(Math.max(snvCount.getValue() / callableSites.getAsDouble(), 1.0e-8));
            log10IndelPrior = Math.log10(Math.max(indelCount.getValue() / callableSites.getAsDouble(), 1.0e-8));
            log10NoVariantPrior = MathUtils.log10OneMinusPow10(MathUtils.log10SumLog10(log10SNVPrior, log10IndelPrior));
        }
        final int variantCount = snvCount.getValue() + indelCount.getValue();
        log10VariantVsArtifactPrior = Math.log10((variantCount + 1) / (variantCount + technicalArtifactCount + 2));
    }

    public List<Pair<String, String>> clusteringMetadata() {
        final List<Pair<String, String>> result = new ArrayList<>();
        result.add(ImmutablePair.of("log10 SNV prior", Double.toString(log10SNVPrior)));
        result.add(ImmutablePair.of("log10 Indel prior", Double.toString(log10IndelPrior)));
        result.add(ImmutablePair.of("High-AF beta-binomial cluster",
                String.format("weight = %.4f, %s", Math.pow(10, log10HighAFWeight), clusters.get(HIGH_AF_INDEX).toString())));
        result.add(ImmutablePair.of("Background beta-binomial cluster",
                String.format("weight = %.4f, %s", Math.pow(10, log10BackgroundWeight), clusters.get(BACKGROUND_INDEX).toString())));

        final MutableInt clusterIndex = new MutableInt(1);
        IntStream.range(OFFSET, clusters.size()).boxed()
                .sorted(Comparator.comparingDouble(c -> -log10CRPWeight(c)))
                .forEach(c -> result.add(ImmutablePair.of("Binomial cluster " + clusterIndex.toString(),
                        String.format("weight = %.4f, %s", Math.pow(10, log10CRPWeight(c)), clusters.get(c).toString()))));
        return result;
    }
}
