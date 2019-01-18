package org.broadinstitute.hellbender.tools.walkers.mutect.filtering;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import org.apache.commons.lang3.mutable.MutableDouble;
import org.apache.commons.lang3.mutable.MutableInt;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.walkers.mutect.Mutect2;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.QualityUtils;
import picard.cmdline.programgroups.VariantFilteringProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.tools.exome.FilterByOrientationBias;
import org.broadinstitute.hellbender.tools.walkers.contamination.CalculateContamination;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;

import java.io.File;
import java.util.*;
import java.util.stream.Collectors;

/**
 * <p>Filter variants in a Mutect2 VCF callset.</p>
 *
 * <p>
 *     FilterMutectCalls encapsulates GATK3 MuTect2's filtering functionality and adds additional filters.
 *     Thresholds for filters are contained in {@link M2FiltersArgumentCollection} and described in
 *     <a href='https://github.com/broadinstitute/gatk/tree/master/docs/mutect/mutect.pdf' target='_blank'>https://github.com/broadinstitute/gatk/tree/master/docs/mutect/mutect.pdf</a>.
 *     To filter based on sequence context artifacts, see {@link FilterByOrientationBias}.
 * </p>
 * <p>
 *     Filtering thresholds for both normal-artifact-lod (default threshold 0.0) and tumor-lod (default threshold 5.3) can be set in this tool.
 *     If the normal artifact log odds is larger than the threshold, then FilterMutectCalls applies the artifact-in-normal filter.
 *     For matched normal analyses with tumor contamination in the normal, consider increasing the normal-artifact-lod threshold.
 *     If the tumor log odds is smaller than the threshold, then FilterMutectCalls filters the variant.
 * </p>
 * <p>
 *     If given a --contamination-table file, e.g. results from
 *     {@link CalculateContamination}, the tool will additionally
 *     filter on contamination fractions. Alternatively, provide a numerical fraction to filter with the --contamination argument.
 * </p>
 * <p>
 *     This tool is featured in the Somatic Short Mutation calling Best Practice Workflow.
 *     See <a href="https://software.broadinstitute.org/gatk/documentation/article?id=11136">Tutorial#11136</a> for a
 *     step-by-step description of the workflow and <a href="https://software.broadinstitute.org/gatk/documentation/article?id=11127">Article#11127</a>
 *     for an overview of what traditional somatic calling entails. For the latest pipeline scripts, see the
 *     <a href="https://github.com/broadinstitute/gatk/tree/master/scripts/mutect2_wdl">Mutect2 WDL scripts directory</a>.
 * </p>
 *
 * <h3>Usage example</h3>
 * <pre>
 * gatk FilterMutectCalls \
 *   -V somatic.vcf.gz \
 *   --contamination-table contamination.table \
 *   -O filtered.vcf.gz
 * </pre>
 *
 */
@CommandLineProgramProperties(
        summary = "Filter somatic SNVs and indels called by Mutect2",
        oneLineSummary = "Filter somatic SNVs and indels called by Mutect2",
        programGroup = VariantFilteringProgramGroup.class
)
@DocumentedFeature
public final class FilterMutectCalls extends MultiplePassVariantWalker {

    public static final String FILTERING_STATS_LONG_NAME = "filtering-stats";

    private static final double EPSILON = 1.0e-10;

    public static final String FILTERING_STATUS_VCF_KEY = "filtering_status";

    public static final String FILTERING_STATS_EXTENSION = ".filteringStats.tsv";

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName =StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc="The output filtered VCF file", optional=false)
    private final String outputVcf = null;

    @Argument(fullName = Mutect2.MUTECT_STATS_SHORT_NAME, doc="The Mutect stats file output by Mutect2", optional=true)
    private final String statsTable = null;

    @Argument(fullName = FILTERING_STATS_LONG_NAME, doc="The output filtering stats file", optional=true)
    private final String filteringStatsOutput = null;

    @ArgumentCollection
    protected M2FiltersArgumentCollection MTFAC = new M2FiltersArgumentCollection();

    private VariantContextWriter vcfWriter;

    private List<Mutect2VariantFilter> filters;

    private Mutect2FilteringInfo filteringInfo;

    private Map<String, MutableDouble> expectedFalsePositivesPerFilter = new HashMap<>();
    private Map<String, MutableDouble> expectedFalseNegativesPerFilter = new HashMap<>();
    private final MutableDouble expectedFalsePositives = new MutableDouble(0);
    private final MutableDouble expectedTruePositives = new MutableDouble(0);
    private final MutableDouble expectedFalseNegatives = new MutableDouble(0);

    private MutableInt passingVariants = new MutableInt(0);

    private Map<String, String> filterPhredPosteriorAnnotations = new HashMap<>();

    @Override
    protected int numberOfPasses() { return 3; }

    @Override
    public void onTraversalStart() {
        final VCFHeader inputHeader = getHeaderForVariants();
        final Set<VCFHeaderLine> headerLines = inputHeader.getMetaDataInSortedOrder().stream()
                .filter(line -> !line.getKey().equals(FILTERING_STATUS_VCF_KEY)) //remove header line from Mutect2 stating that calls are unfiltered.
                .collect(Collectors.toSet());
        headerLines.add(new VCFHeaderLine(FILTERING_STATUS_VCF_KEY, "These calls have been filtered by " + FilterMutectCalls.class.getSimpleName() + " to label false positives with a list of failed filters and true positives with PASS."));

        GATKVCFConstants.MUTECT_FILTER_NAMES.stream().map(GATKVCFHeaderLines::getFilterLine).forEach(headerLines::add);

        headerLines.addAll(getDefaultToolVCFHeaderLines());

        final VCFHeader vcfHeader = new VCFHeader(headerLines, inputHeader.getGenotypeSamples());
        vcfWriter = createVCFWriter(new File(outputVcf));
        vcfWriter.writeHeader(vcfHeader);


        final File mutect2StatsTable = new File(statsTable == null ? drivingVariantFile + Mutect2.DEFAULT_STATS_EXTENSION : statsTable);

        filteringInfo = new Mutect2FilteringInfo(MTFAC, vcfHeader);
        if (mutect2StatsTable.exists()) {
            filteringInfo.inputMutectStats(mutect2StatsTable);
        } else {
            logger.warn("Mutect stats table " + mutect2StatsTable + " not found.  Filtering will proceed without this information.");
        }
        filters = new ArrayList<>();
        filters.add(new TumorEvidenceFilter());
        filters.add(new BaseQualityFilter());
        filters.add(new MappingQualityFilter());
        filters.add(new DuplicatedAltReadFilter());
        filters.add(new StrandArtifactFilter());
        filters.add(new ContaminationFilter());
        filters.add(new PanelOfNormalsFilter());
        filters.add(new NormalArtifactFilter());
        filters.add(new ReadOrientationFilter());
        filters.add(new NRatioFilter());
        filters.add(new StrictStrandBiasFilter());
        filters.add(new ReadPositionFilter());

        if (MTFAC.mitochondria) {
            filters.add(new LogOddsOverDepthFilter());
            filters.add(new ChimericOriginalAlignmentFilter());
        } else {
            filters.add(new ClusteredEventsFilter());
            filters.add(new MultiallelicFilter());
            filters.add(new FragmentLengthFilter());
            filters.add(new PolymeraseSlippageFilter());
            filters.add(new FilteredHaplotypeFilter());
            filters.add(new GermlineFilter());
        }

        filters.forEach(filter -> expectedFalsePositivesPerFilter.put(filter.filterName(), new MutableDouble(0)));
        filters.forEach(filter -> expectedFalseNegativesPerFilter.put(filter.filterName(), new MutableDouble(0)));

        for (final Mutect2VariantFilter filter : filters) {
            filter.phredScaledPosteriorAnnotationName().ifPresent(ann -> filterPhredPosteriorAnnotations.put(filter.filterName(), ann));
        }
    }

    @Override
    protected void nthPassApply(final VariantContext variant,
                                final ReadsContext readsContext,
                                final ReferenceContext referenceContext,
                                final FeatureContext featureContext,
                                final int n) {
        if (n == 0) {
            filters.forEach(f -> f.accumulateDataForLearning(variant, filteringInfo));

            double artifactProbability = 0;
            double technicalArtifactProbability = 0;

            for (final Mutect2VariantFilter filter : filters) {
                final double prob = filter.artifactProbability(variant, filteringInfo);
                artifactProbability = Math.max(artifactProbability, prob);
                technicalArtifactProbability = filter.isTechnicalArtifact() ? Math.max(technicalArtifactProbability, prob) : technicalArtifactProbability;
            }

            filteringInfo.addRealVariantCount(1 - artifactProbability);
            filteringInfo.addTechnicalArtifactCount(technicalArtifactProbability);

            // TODO: if real variant probability is significant, send allele fraction info the filteringInfo
        } else if (n == 1) {
            secondPassOptimizeThresholdApply(variant, readsContext, referenceContext, featureContext);
        } else if (n == 2) {
            thirdPassMakeFilteringCallsApply(variant, readsContext, referenceContext, featureContext);
        } else {
            throw new GATKException.ShouldNeverReachHereException("This two-pass walker should never reach (zero-indexed) pass " + n);
        }
    }

    @Override
    protected void afterNthPass(final int n) {
        if (n == 0) {
            filteringInfo.learnPriorProbOfArtifactVersusVariant();
            filteringInfo.learnPriorProbOfVariant();
            filters.forEach(f -> f.learnParameters());
        } else if (n == 1) {
            filteringInfo.adjustThreshold();
        } else if (n == 2) {
            final int totalCalls = passingVariants.getValue();
            final List<FilterStats> filterStats = filters.stream().map(Mutect2VariantFilter::filterName)
                    .map(filter -> {
                        final double falseNegatives = expectedFalseNegativesPerFilter.get(filter).getValue();
                        final double falsePositives = expectedFalsePositivesPerFilter.get(filter).getValue();
                        final double fdr = falsePositives / totalCalls;
                        final double totalTrueVariants = expectedTruePositives.getValue() + expectedFalseNegatives.getValue();

                        return new FilterStats(filter, falsePositives, fdr, falseNegatives, falseNegatives / totalTrueVariants);
                    })
                    .filter(stats -> stats.getFalsePositiveCount() > 0 || stats.getFalseNegativeCount() > 0)
                    .collect(Collectors.toList());

            final File filteringStatsFile = new File(filteringStatsOutput != null ? filteringStatsOutput : outputVcf + FILTERING_STATS_EXTENSION);

            FilterStats.writeM2FilterSummary(filterStats, filteringStatsFile, filteringInfo.getArtifactProbabilityThreshold(), totalCalls,
                    expectedTruePositives.getValue(), expectedFalsePositives.getValue(), expectedFalseNegatives.getValue());
        } else {
            throw new GATKException.ShouldNeverReachHereException("This three-pass walker should never reach (zero-indexed) pass " + n);
        }
    }

    private void secondPassOptimizeThresholdApply(final VariantContext vc, final ReadsContext readsContext, final ReferenceContext refContext, final FeatureContext fc) {

        final double[] artifactProbabilities = filters.stream().mapToDouble(filter -> filter.artifactProbability(vc, filteringInfo)).toArray();
        final double artifactProbability = MathUtils.arrayMax(artifactProbabilities);

        if (artifactProbability > filteringInfo.getArtifactProbabilityThreshold() - EPSILON) {
            filteringInfo.recordFilteredHaplotypes(vc);
        }

        filteringInfo.addFirstPassArtifactProbability(artifactProbability);
    }

    public void thirdPassMakeFilteringCallsApply(final VariantContext vc, final ReadsContext readsContext, final ReferenceContext refContext, final FeatureContext fc) {
        final VariantContextBuilder vcb = new VariantContextBuilder(vc);
        vcb.filters(new HashSet<>());

        final Map<Mutect2VariantFilter, Double> artifactProbabilities = filters.stream()
                .collect(Collectors.toMap(f -> f, f -> f.artifactProbability(vc, filteringInfo)));

        final double maxArtifactProbability = artifactProbabilities.values().stream()
                .mapToDouble(x->x).max().orElseGet(() -> 0);

        final boolean filtered = maxArtifactProbability > filteringInfo.getArtifactProbabilityThreshold() - EPSILON;

        if (filtered) {
            expectedFalseNegatives.add(1 - maxArtifactProbability);
        } else {
            passingVariants.increment();
            expectedFalsePositives.add(maxArtifactProbability);
            expectedTruePositives.add(1 - maxArtifactProbability);
        }

        for (final Map.Entry<Mutect2VariantFilter, Double> entry : artifactProbabilities.entrySet()) {
            final String filter = entry.getKey().filterName();
            final double artifactProbability = entry.getValue();

            final String posteriorAnnotation = filterPhredPosteriorAnnotations.get(filter);
            if (posteriorAnnotation != null && entry.getKey().requiredAnnotations().stream().allMatch(vc::hasAttribute)) {
                vcb.attribute(posteriorAnnotation, QualityUtils.errorProbToQual(artifactProbability));
            }

            if (artifactProbability > EPSILON && artifactProbability > filteringInfo.getArtifactProbabilityThreshold() - EPSILON) {
                vcb.filter(filter);
                expectedFalseNegativesPerFilter.get(filter).add(1 - maxArtifactProbability);
            } else if (!filtered) {
                expectedFalsePositivesPerFilter.get(filter).add(artifactProbability);
            }
        }

        vcfWriter.add(vcb.make());
    }

    @Override
    public void closeTool() {
        if ( vcfWriter != null ) {
            vcfWriter.close();
        }
    }

}
