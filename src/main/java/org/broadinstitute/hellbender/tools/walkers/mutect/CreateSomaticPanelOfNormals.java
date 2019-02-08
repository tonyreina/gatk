package org.broadinstitute.hellbender.tools.walkers.mutect;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.apache.commons.math3.special.Gamma;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.VariantWalker;
import org.broadinstitute.hellbender.tools.walkers.readorientation.BetaDistributionShape;
import org.broadinstitute.hellbender.tools.walkers.validation.basicshortmutpileup.BetaBinomialDistribution;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.OptimizationUtils;
import picard.cmdline.programgroups.VariantFilteringProgramGroup;

import java.io.File;
import java.util.List;
import java.util.OptionalInt;
import java.util.function.DoubleUnaryOperator;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Create a panel of normals (PoN) containing germline and artifactual sites for use with Mutect2.
 *
 * <p>
 *     The tool takes multiple normal sample callsets produced by {@link Mutect2}'s tumor-only mode and collates sites present in two or more samples
 *     into a sites-only VCF. The PoN captures common artifactual and germline variant sites.
 *     Mutect2 then uses the PoN to filter variants at the site-level.
 * </p>
 * <p>
 *     This tool is featured in the Somatic Short Mutation calling Best Practice Workflow.
 *     See <a href="https://software.broadinstitute.org/gatk/documentation/article?id=11136">Tutorial#11136</a> for a
 *     step-by-step description of the workflow and <a href="https://software.broadinstitute.org/gatk/documentation/article?id=11127">Article#11127</a>
 *     for an overview of what traditional somatic calling entails. For the latest pipeline scripts, see the
 *     <a href="https://github.com/broadinstitute/gatk/tree/master/scripts/mutect2_wdl">Mutect2 WDL scripts directory</a>.
 * </p>
 * <h3>Example workflow</h3>
 *
 * <h4>Step 1. Run Mutect2 in tumor-only mode for each normal sample.</h4>
 * <pre>
 * gatk Mutect2 -R reference.fasta -I normal1.bam -O normal1.vcf.gz
 * </pre>
 *
 * <h4>Step 2. Create a GenomicsDB from the normal Mutect2 calls.</h4>
 *
 *  <pre>
 *    gatk GenomicsDBImport \
 *       --genomicsdb-workspace-path pon_db \
 *       -V normal1.vcf.gz \
 *       -V normal2.vcf.gz \
 *       -V normal3.vcf.gz
 *  </pre>
 *
 * <h4>Step 3. Combine the normal calls using CreateSomaticPanelOfNormals.</h4>
 *
 * <pre>
 * gatk CreateSomaticPanelOfNormals -V gendb:///pon_db -O pon.vcf.gz
 * </pre>
 */
@CommandLineProgramProperties(
        summary = "Make a panel of normals (PoN) for use with Mutect2",
        oneLineSummary = "Make a panel of normals for use with Mutect2",
        programGroup = VariantFilteringProgramGroup.class
)
@DocumentedFeature
@BetaFeature
public class CreateSomaticPanelOfNormals extends VariantWalker {

    public static final String MIN_SAMPLE_COUNT_LONG_NAME = "min-sample-count";
    public static final int DEFAULT_MIN_SAMPLE_COUNT = 2;

    public static final String FRACTION_INFO_FIELD = "FRACTION";
    public static final String BETA_SHAPE_INFO_FIELD = "BETA";

    @Argument(fullName = MIN_SAMPLE_COUNT_LONG_NAME, doc="Number of samples containing a variant site required to include it in the panel of normals.", optional = true)
    private int minSampleCount = DEFAULT_MIN_SAMPLE_COUNT;

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc="Output vcf")
    private String outputVcf;

    private VariantContextWriter vcfWriter;

    private int numSamples;

    @Override
    public void onTraversalStart() {
        final VCFHeader outputHeader = new VCFHeader();
        getDefaultToolVCFHeaderLines().forEach(outputHeader::addMetaDataLine);
        outputHeader.setSequenceDictionary(getHeaderForVariants().getSequenceDictionary());

        outputHeader.addMetaDataLine(new VCFInfoHeaderLine(FRACTION_INFO_FIELD, 1, VCFHeaderLineType.Float, "Fraction of samples exhibiting artifact"));
        outputHeader.addMetaDataLine(new VCFInfoHeaderLine(BETA_SHAPE_INFO_FIELD, 2, VCFHeaderLineType.Float, "Beta distribution parameters to fit artifact allele fractions"));

        vcfWriter = createVCFWriter(new File(outputVcf));
        vcfWriter.writeHeader(outputHeader);

        numSamples = getHeaderForVariants().getNGenotypeSamples();
    }

    @Override
    public void apply(final VariantContext vc, final ReadsContext rc, final ReferenceContext ref, final FeatureContext fc) {

        final OptionalInt spanDelAltIndex = IntStream.range(0, vc.getNAlleles() - 1)
                .filter(n -> vc.getAlternateAllele(n).basesMatch(Allele.SPAN_DEL))
                .findFirst();

        final List<Genotype> variantGenotypes = vc.getGenotypes().stream()
                .filter(g -> altCount(g, spanDelAltIndex) > 0).collect(Collectors.toList());

        if (variantGenotypes.size() < minSampleCount) {
            return;
        }

        final double fraction = (double) variantGenotypes.size() / numSamples;

        final List<int[]> altAndRefCounts = variantGenotypes.stream()
                .map(g -> new int[] {altCount(g, spanDelAltIndex), g.getAD()[0]})
                .collect(Collectors.toList());

        final BetaDistributionShape betaDistributionShape = fitBeta(altAndRefCounts);

        final VariantContext outputVc = new VariantContextBuilder(vc.getSource(), vc.getContig(), vc.getStart(), vc.getEnd(), vc.getAlleles())
                .attribute(FRACTION_INFO_FIELD, fraction)
                .attribute(BETA_SHAPE_INFO_FIELD, new double[] { betaDistributionShape.getAlpha(), betaDistributionShape.getBeta()})
                .make();

        vcfWriter.add(outputVc);
    }

    private static final int altCount(final Genotype g, final OptionalInt spanDelAltIndex) {
        if (!g.hasAD()) {
            return 0;
        }
        final int[] AD = g.getAD();
        return (int) MathUtils.sum(AD) - AD[0] - (spanDelAltIndex.isPresent() ? AD[spanDelAltIndex.getAsInt() + 1] : 0);
    }

    private BetaDistributionShape fitBeta(final List<int[]> altAndRefCounts) {
        final int totalAltCount = altAndRefCounts.stream().mapToInt(pair -> pair[0]).sum();
        final int totalRefCount = altAndRefCounts.stream().mapToInt(pair -> pair[1]).sum();
        final int min = Math.min(totalAltCount, totalRefCount);

        // keeping the ratio of alpha and beta equal to the ratio of baseAlpha and baseBeta gives the empirical mean
        final double baseAlpha = (totalAltCount + 1.0) / (min + 1);
        final double baseBeta = (totalRefCount + 1.0) / (min + 1);

        final DoubleUnaryOperator logLikelihood = s -> {
            final double alpha = baseAlpha * s;
            final double beta = baseBeta * s;

            return altAndRefCounts.stream().mapToDouble(pair -> {
                final int n = pair[0] + pair[1];
                final int k = pair[0];
                return new BetaBinomialDistribution(null, alpha, beta, n).logProbability(k);
            }).sum();
        };

        final double scale = OptimizationUtils.max(logLikelihood, 0.01, 100, 1, 0.01, 0.1, 100).getPoint();

        return new BetaDistributionShape(baseAlpha * scale, baseBeta * scale);
    }

    @Override
    public Object onTraversalSuccess() {
        return "SUCCESS";
    }

    @Override
    public void closeTool() {
        if ( vcfWriter != null ) {
            vcfWriter.close();
        }
    }
}
