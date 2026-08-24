package com.synapticfour.gatkrs.parity;

import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeCalculationArgumentCollection;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.hellbender.tools.walkers.genotyper.afcalc.AFCalculationResult;
import org.broadinstitute.hellbender.tools.walkers.genotyper.afcalc.AlleleFrequencyCalculator;

import java.io.BufferedReader;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.Collections;
import java.util.Locale;

/** Parity dump for GATK {@link AlleleFrequencyCalculator} (G-D01 / g2-af). */
public final class HcParityAfEm {

    private HcParityAfEm() {}

    /**
     * Site QUAL (Phred) from log10 GLs after HTSJDK PL round-trip + biallelic AF EM
     * ({@code use-posteriors-to-calculate-qual} false). Used by {@code HcParityRegionVcf}
     * synthetic VCF dumps.
     */
    public static double siteQualPhredFromLog10Gl(final double[] log10Gl) {
        final GenotypeCalculationArgumentCollection genotypeArgs =
                new GenotypeCalculationArgumentCollection();
        final Genotype genotype =
                new GenotypeBuilder("sample")
                        .alleles(
                                Arrays.asList(
                                        Allele.create("A", true), Allele.create("C", false)))
                        .PL(log10Gl)
                        .make();
        final double[] glAfterPl = genotype.getLikelihoods().getAsVector();
        final HcParityAfEmIterations.EmMetrics metrics =
                HcParityAfEmIterations.computeBiallelicMetrics(glAfterPl, genotypeArgs);
        return -10.0 * metrics.log10PNoVariant;
    }

    public static void dumpAfEmFixture(final Path fixture) throws Exception {
        final GenotypeCalculationArgumentCollection genotypeArgs =
                new GenotypeCalculationArgumentCollection();
        final AlleleFrequencyCalculator calculator =
                AlleleFrequencyCalculator.makeCalculator(genotypeArgs);
        final Allele ref = Allele.create("A", true);
        final Allele alt = Allele.create("C", false);
        try (BufferedReader br = Files.newBufferedReader(fixture, StandardCharsets.UTF_8)) {
            String line;
            while ((line = br.readLine()) != null) {
                final String t = line.trim();
                if (t.isEmpty() || t.startsWith("#")) {
                    continue;
                }
                final String[] cols = t.split("\t");
                if (cols.length < 2) {
                    throw new IllegalArgumentException("af-em row needs gl: " + t);
                }
                final String caseId = cols[0];
                final double[] gl = parseCsvDoubles(cols[1]);
                final Genotype genotype =
                        new GenotypeBuilder("sample")
                                .alleles(Arrays.asList(ref, alt))
                                .PL(gl)
                                .make();
                final double[] glAfterPl = genotype.getLikelihoods().getAsVector();
                final HcParityAfEmIterations.EmMetrics metrics =
                        HcParityAfEmIterations.computeBiallelicMetrics(glAfterPl, genotypeArgs);
                calculator.calculate(
                        new VariantContextBuilder("parity", "1", 1, 1, Arrays.asList(ref, alt))
                                .genotypes(Collections.singletonList(genotype))
                                .make());
                System.out.println(caseId + "\talt_ac\t" + metrics.altAc);
                System.out.printf(Locale.ROOT, "%s\taf\t%.6f%n", caseId, metrics.af);
                System.out.printf(
                        Locale.ROOT,
                        "%s\tlog10_p_no_variant\t%.6f%n",
                        caseId,
                        metrics.log10PNoVariant);
                System.out.println(caseId + "\tem_iterations\t" + metrics.iterations);
            }
        }
    }

    private static double[] parseCsvDoubles(final String csv) {
        final String[] parts = csv.split(",");
        final double[] out = new double[parts.length];
        for (int i = 0; i < parts.length; i++) {
            out[i] = Double.parseDouble(parts[i].trim());
        }
        return out;
    }
}
