package org.broadinstitute.hellbender.tools.walkers.mutect.filtering;

import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.lang3.mutable.MutableDouble;
import org.apache.commons.lang3.mutable.MutableInt;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.math3.random.RandomDataGenerator;
import org.bouncycastle.crypto.prng.RandomGenerator;
import org.broadinstitute.hellbender.tools.walkers.mutect.Mutect2Engine;
import org.broadinstitute.hellbender.tools.walkers.mutect.MutectStats;
import org.broadinstitute.hellbender.tools.walkers.readorientation.BetaDistributionShape;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;

public class SomaticPriorModel {
    private double log10SNVPrior;
    private double log10IndelPrior;
    private double artifactVsVariantPrior;
    private final OptionalDouble callableSites;

    final List<Datum> data = new ArrayList<>();

    private final MutableDouble realVariantCount = new MutableDouble(0);
    private final MutableDouble realSNVCount = new MutableDouble(0);
    private final MutableDouble realIndelCount = new MutableDouble(0);
    private final MutableDouble technicalArtifactCount = new MutableDouble(0);

    public SomaticPriorModel(final M2FiltersArgumentCollection MTFAC, final List<MutectStats> mutectStats) {
        log10SNVPrior = MTFAC.log10SNVPrior;
        log10IndelPrior = MTFAC.log10IndelPrior;
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

    public void recordDatum(final int[] tumorADs, final double[] tumorLog10Odds, final double artifactProbability) {
        final int totalAD = (int) MathUtils.sum(tumorADs);
        for (int i = 0; i < tumorLog10Odds.length; i++) {
            data.add(new Datum(tumorLog10Odds[i], artifactProbability, tumorADs[i+1], totalAD));
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
        final RandomDataGenerator rng = Utils.getRandomDataGenerator()

        final int noAssignment = -1;
        final int artifactTable = 0;
        final int totalPatrons = data.size();

        // Chinese Restaurant tables:
        // table 0: artifacts
        // table 1: Beta background
        // tables 2, 3 . . . single-AF tables

        final List<Integer> tableAssignments = Collections.nCopies(data.size(), noAssignment);
        final List<Pair<BetaDistributionShape, MutableInt>> betaShapesAndCounts = new ArrayList<>();
        betaShapesAndCounts.add(ImmutablePair.of(null, new MutableInt(0)));

        for (int iteration = 0; iteration < 5; iteration++) {
            for (int n = 0; n < totalPatrons; n++) {
                final Datum datum = data.get(n);
                // stochastically assign to artifact
                if (rng.nextUniform(0,1) < datum.artifactProb) {
                    tableAssignments.set(n, artifactTable);
                    betaShapesAndCounts.get(artifactTable).getRight().increment();
                } else {
                    final int table = 9;
                    final double[] likelihoods

                            tableAssignments.set(n, table);
                    betaShapesAndCounts.get(table).getRight().increment();
                }
            }
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
        final double tumorLog10Odds;
        final double artifactProb;
        final int altCount;
        final int totalCount;

        public Datum(final double tumorLog10Odds, final double artifactProb, final int altCount, final int totalCount) {
            this.tumorLog10Odds = tumorLog10Odds;
            this.artifactProb = artifactProb;
            this.altCount = altCount;
            this.totalCount = totalCount;
        }
    }
}
