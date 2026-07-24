package org.broadinstitute.hellbender.tools.walkers.genotyper;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import java.io.BufferedReader;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.List;
import java.util.Locale;
import java.util.stream.Collectors;

/** Parity dump for GATK {@link AlleleSubsettingUtils#subsetAlleles}. */
public final class HcParitySubsetAlleles {

    private HcParitySubsetAlleles() {}

    public static void dumpSubsetAllelesVcFixture(final Path fixture) throws Exception {
        final Allele aref = Allele.create("A", true);
        final Allele c = Allele.create("C", false);
        final Allele g = Allele.create("G", false);
        try (BufferedReader br = Files.newBufferedReader(fixture, StandardCharsets.UTF_8)) {
            String line;
            while ((line = br.readLine()) != null) {
                final String t = line.trim();
                if (t.isEmpty() || t.startsWith("#")) {
                    continue;
                }
                final String[] cols = t.split("\t");
                final String caseId = cols[0];
                final double[] log10Pl = parseCsvDoubles(cols[1]);
                final int[] ad = parseCsvInts(cols[2]);
                final int[] sac = parseCsvInts(cols[3]);
                final Genotype gt =
                        new GenotypeBuilder("sample")
                                .alleles(Arrays.asList(aref, c))
                                .PL(log10Pl)
                                .AD(ad)
                                .attribute(GATKVCFConstants.STRAND_COUNT_BY_SAMPLE_KEY, sac)
                                .GQ(200)
                                .make();
                final VariantContext original =
                        new VariantContextBuilder("parity", "20", 10, 10, Arrays.asList(aref, c, g))
                                .genotypes(GenotypesContext.create(gt))
                                .make();
                final List<Allele> keep = Arrays.asList(aref, c);
                final GenotypesContext subset =
                        AlleleSubsettingUtils.subsetAlleles(
                                original.getGenotypes(),
                                2,
                                original.getAlleles(),
                                keep,
                                null,
                                GenotypeAssignmentMethod.USE_PLS_TO_ASSIGN);
                final Genotype out = subset.get("sample");
                System.out.println(caseId + "\tallele_count_before\t3");
                System.out.println(caseId + "\tallele_count_after\t2");
                printSubsetVcRows(caseId, out);
            }
        }
    }

    public static void dumpSubsetAllelesFixture(final Path fixture) throws Exception {
        final Allele aref = Allele.create("A", true);
        final Allele c = Allele.create("C", false);
        final Allele g = Allele.create("G", false);
        try (BufferedReader br = Files.newBufferedReader(fixture, StandardCharsets.UTF_8)) {
            String line;
            while ((line = br.readLine()) != null) {
                final String t = line.trim();
                if (t.isEmpty() || t.startsWith("#")) {
                    continue;
                }
                final String[] cols = t.split("\t");
                final String caseId = cols[0];
                final double[] log10Pl = parseCsvDoubles(cols[1]);
                final int[] ad = parseCsvInts(cols[2]);
                final List<Allele> keep;
                if ("trim_acg_to_ag".equals(caseId)) {
                    keep = Arrays.asList(aref, g);
                } else {
                    keep = Arrays.asList(aref, c);
                }
                final Genotype gt =
                        new GenotypeBuilder("sample")
                                .alleles(keep)
                                .PL(log10Pl)
                                .AD(ad)
                                .GQ(200)
                                .make();
                final VariantContext original =
                        new VariantContextBuilder("parity", "20", 10, 10, Arrays.asList(aref, c, g))
                                .genotypes(GenotypesContext.create(gt))
                                .make();
                final GenotypesContext subset =
                        AlleleSubsettingUtils.subsetAlleles(
                                original.getGenotypes(),
                                2,
                                original.getAlleles(),
                                keep,
                                null,
                                GenotypeAssignmentMethod.USE_PLS_TO_ASSIGN);
                final Genotype out = subset.get("sample");
                System.out.println(caseId + "\tallele_count_before\t3");
                System.out.println(caseId + "\tallele_count_after\t2");
                printSubsetVcRows(caseId, out);
            }
        }
    }

    private static void printSubsetVcRows(final String caseId, final Genotype out) {
        System.out.println(
                caseId + "\tpl_length\t" + (out.hasPL() ? out.getPL().length : 0));
        if (out.hasPL()) {
            for (int i = 0; i < out.getPL().length; i++) {
                System.out.printf(Locale.ROOT, "%s\tpl_%d\t%d%n", caseId, i, out.getPL()[i]);
            }
        }
        if (out.hasAD()) {
            System.out.println(
                    caseId
                            + "\tad\t"
                            + Arrays.stream(out.getAD())
                                    .mapToObj(Integer::toString)
                                    .collect(Collectors.joining(",")));
        }
        if (out.hasGQ()) {
            System.out.println(caseId + "\tgq\t" + out.getGQ());
        }
        if (out.hasExtendedAttribute(GATKVCFConstants.STRAND_COUNT_BY_SAMPLE_KEY)) {
            final int[] newSac =
                    (int[]) out.getExtendedAttribute(GATKVCFConstants.STRAND_COUNT_BY_SAMPLE_KEY);
            System.out.println(
                    caseId
                            + "\tsac\t"
                            + Arrays.stream(newSac)
                                    .mapToObj(Integer::toString)
                                    .collect(Collectors.joining(",")));
        }
    }

    private static int[] parseCsvInts(final String csv) {
        final String[] parts = csv.split(",");
        final int[] out = new int[parts.length];
        for (int i = 0; i < parts.length; i++) {
            out[i] = Integer.parseInt(parts[i].trim());
        }
        return out;
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
