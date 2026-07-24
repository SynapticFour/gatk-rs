package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Locale;
import java.util.stream.Collectors;
import org.broadinstitute.hellbender.tools.walkers.ReferenceConfidenceVariantContextMerger;
import org.broadinstitute.hellbender.tools.walkers.annotator.VariantAnnotatorEngine;
import org.broadinstitute.hellbender.utils.SimpleInterval;

/**
 * H-D02: production {@link ReferenceConfidenceVariantContextMerger#merge} parity dumps.
 */
public final class HcParityGvcfMerger {

    private HcParityGvcfMerger() {}

    public static void dumpMergeCase(final String caseId) {
        final List<VariantContext> toMerge = buildUnitTestMergeInputs(caseId);
        if (toMerge == null) {
            throw new IllegalArgumentException("unknown merge case: " + caseId);
        }
        final int start = 10;
        final SimpleInterval loc = new SimpleInterval("20", start, start);
        final VariantAnnotatorEngine engine =
                new VariantAnnotatorEngine(
                        Collections.emptyList(),
                        null,
                        Collections.emptyList(),
                        false,
                        false);
        final ReferenceConfidenceVariantContextMerger merger =
                new ReferenceConfidenceVariantContextMerger(
                        engine, new htsjdk.variant.vcf.VCFHeader());
        final VariantContext merged =
                merger.merge(toMerge, loc, null, true, false);
        if (merged == null) {
            System.out.println("case_id\t" + caseId);
            System.out.println("contig\t" + loc.getContig());
            System.out.println("pos\t" + loc.getStart());
            System.out.println("merged_null\ttrue");
            return;
        }
        printMerged(caseId, merged);
    }

    /** @deprecated use {@link #dumpMergeCase(String)} */
    public static void dumpMergeSnpsRef() {
        dumpMergeCase("merge_snps_ref");
    }

    private static void printMerged(final String caseId, final VariantContext merged) {
        System.out.println("case_id\t" + caseId);
        System.out.println("contig\t" + merged.getContig());
        System.out.println("pos\t" + merged.getStart());
        System.out.println("merged_null\tfalse");
        System.out.println("allele_count\t" + merged.getAlleles().size());
        for (int i = 0; i < merged.getAlleles().size(); i++) {
            System.out.println("allele_" + i + "\t" + merged.getAlleles().get(i).getDisplayString());
        }
        System.out.println("sample_count\t" + merged.getNSamples());
        int gi = 0;
        for (final Genotype g : merged.getGenotypes()) {
            System.out.println("genotype_" + gi + "_name\t" + g.getSampleName());
            if (g.hasPL()) {
                System.out.println(
                        "genotype_"
                                + gi
                                + "_pl\t"
                                + Arrays.stream(g.getPL())
                                        .mapToObj(Integer::toString)
                                        .collect(Collectors.joining(",")));
            }
            if (g.hasGQ()) {
                System.out.printf(Locale.ROOT, "genotype_%d_gq\t%d%n", gi, g.getGQ());
            }
            if (g.hasAD()) {
                System.out.println(
                        "genotype_"
                                + gi
                                + "_ad\t"
                                + Arrays.stream(g.getAD())
                                        .mapToObj(Integer::toString)
                                        .collect(Collectors.joining(",")));
            }
            gi++;
        }
        System.out.println("has_non_ref\t" + merged.hasAllele(Allele.NON_REF_ALLELE));
    }

    private static List<VariantContext> buildUnitTestMergeInputs(final String caseId) {
        final Allele aref = Allele.create("A", true);
        final Allele c = Allele.create("C", false);
        final Allele g = Allele.create("G", false);
        final Allele atc = Allele.create("ATC", false);
        final Allele aaref = Allele.create("AA", true);
        final Allele aAlt = Allele.create("A", false);
        final int start = 10;
        final VariantContext vcBase =
                new VariantContextBuilder("test", "20", start, start, Arrays.asList(aref)).make();
        final VariantContext vcBase2 =
                new VariantContextBuilder("test2", "20", start, start, Arrays.asList(aref)).make();
        final VariantContext vcPrevBase =
                new VariantContextBuilder("test", "20", start - 1, start - 1, Arrays.asList(aref))
                        .make();
        final List<Allele> noCalls = Arrays.asList(Allele.NO_CALL, Allele.NO_CALL);
        final int[] standardPLs = new int[] {30, 20, 10, 71, 72, 73};

        switch (caseId) {
            case "merge_single_vc":
            case "test00":
                return Arrays.asList(
                        vcWith(
                                vcBase,
                                Arrays.asList(aref, c, Allele.NON_REF_ALLELE),
                                g("A_C", standardPLs, noCalls)));
            case "merge_two_snps":
            case "test01":
                return Arrays.asList(
                        vcWith(
                                vcBase,
                                Arrays.asList(aref, c, Allele.NON_REF_ALLELE),
                                g("A_C", standardPLs, noCalls)),
                        vcWith(
                                vcBase,
                                Arrays.asList(aref, g, Allele.NON_REF_ALLELE),
                                g("A_G", standardPLs, noCalls)));
            case "merge_snp_indel":
            case "test02":
                return Arrays.asList(
                        vcWith(
                                vcBase,
                                Arrays.asList(aref, c, Allele.NON_REF_ALLELE),
                                g("A_C", standardPLs, noCalls)),
                        vcWith(
                                vcBase,
                                Arrays.asList(aref, atc, Allele.NON_REF_ALLELE),
                                g("A_ATC", standardPLs, noCalls)));
            case "merge_snp_three_alleles":
            case "test03":
                return Arrays.asList(
                        vcWith(
                                vcBase,
                                Arrays.asList(aref, c, Allele.NON_REF_ALLELE),
                                g("A_C", standardPLs, noCalls)),
                        vcWith(
                                vcBase,
                                Arrays.asList(aref, c, g, Allele.NON_REF_ALLELE),
                                g(
                                        "A_C_G",
                                        new int[] {40, 20, 30, 20, 10, 30, 71, 72, 73, 74},
                                        noCalls)));
            case "merge_snps_ref":
            case "test04":
                return Arrays.asList(
                        vcWith(
                                vcBase,
                                Arrays.asList(aref, c, Allele.NON_REF_ALLELE),
                                g("A_C", standardPLs, noCalls)),
                        vcWith(
                                vcBase,
                                Arrays.asList(aref, Allele.NON_REF_ALLELE),
                                g("A", new int[] {0, 100, 1000}, noCalls)));
            case "merge_spanning_del":
            case "test06":
                return Arrays.asList(
                        vcWith(
                                vcBase,
                                Arrays.asList(aref, c, Allele.NON_REF_ALLELE),
                                g("A_C", standardPLs, noCalls)),
                        vcWith(
                                vcPrevBase,
                                Arrays.asList(aaref, aAlt, Allele.NON_REF_ALLELE),
                                g("AA_A", standardPLs, noCalls)));
            case "merge_all_combined":
            case "test07":
                return Arrays.asList(
                        vcWith(vcBase, Arrays.asList(aref, c, Allele.NON_REF_ALLELE), g("A_C", standardPLs, noCalls)),
                        vcWith(vcBase, Arrays.asList(aref, g, Allele.NON_REF_ALLELE), g("A_G", standardPLs, noCalls)),
                        vcWith(vcBase, Arrays.asList(aref, atc, Allele.NON_REF_ALLELE), g("A_ATC", standardPLs, noCalls)),
                        vcWith(
                                vcBase,
                                Arrays.asList(aref, c, g, Allele.NON_REF_ALLELE),
                                g("A_C_G", new int[] {40, 20, 30, 20, 10, 30, 71, 72, 73, 74}, noCalls)),
                        vcWith(vcBase, Arrays.asList(aref, Allele.NON_REF_ALLELE), g("A", new int[] {0, 100, 1000}, noCalls)),
                        vcWith(
                                vcPrevBase,
                                Arrays.asList(aaref, Allele.NON_REF_ALLELE),
                                g("AA", new int[] {0, 80, 800}, noCalls)),
                        vcWith(
                                vcPrevBase,
                                Arrays.asList(aaref, aAlt, Allele.NON_REF_ALLELE),
                                g("AA_A", standardPLs, noCalls)));
            case "merge_spanning_ref_only":
            case "test08":
                return Arrays.asList(
                        vcWith(
                                vcPrevBase,
                                Arrays.asList(aaref, Allele.NON_REF_ALLELE),
                                g("AA", new int[] {0, 80, 800}, noCalls)));
            case "merge_ad_pl_mix":
            case "test12":
                return Arrays.asList(
                        vcWith(
                                vcBase,
                                Arrays.asList(aref, atc, Allele.NON_REF_ALLELE),
                                gAdPl("A_ATC", standardPLs, new int[] {20, 10}, noCalls)),
                        vcWith(
                                vcBase,
                                Arrays.asList(aref, c, g, Allele.NON_REF_ALLELE),
                                gAdPl(
                                        "A_C_G",
                                        new int[] {40, 20, 30, 20, 10, 30, 71, 72, 73, 74},
                                        new int[] {30, 0, 8},
                                        noCalls)));
            case "merge_ad_only_overlap":
            case "test13":
                return Arrays.asList(
                        vcWith(
                                vcBase,
                                Arrays.asList(aref, c, g, Allele.NON_REF_ALLELE),
                                gAd("A_C_G", new int[] {60, 9, 20}, noCalls)),
                        vcWith(
                                vcBase,
                                Arrays.asList(aref, c, g, Allele.NON_REF_ALLELE),
                                gAd("A_C_G", new int[] {60, 9, 20}, noCalls)));
            case "merge_ad_only_distinct":
            case "test14":
                final Allele aa = Allele.create("AA", false);
                return Arrays.asList(
                        vcWith(
                                vcBase,
                                Arrays.asList(aref, c, g, Allele.NON_REF_ALLELE),
                                gAd("A_C_G", new int[] {60, 9, 20}, noCalls)),
                        vcWith(
                                vcBase,
                                Arrays.asList(aref, atc, aa, Allele.NON_REF_ALLELE),
                                gAd("A_ATC_AA", new int[] {30, 8, 40}, noCalls)));
            default:
                return null;
        }
    }

    private static Genotype gAd(final String name, final int[] ad, final List<Allele> alleles) {
        return new GenotypeBuilder(name).AD(ad).alleles(alleles).make();
    }

    private static Genotype gAdPl(
            final String name,
            final int[] pl,
            final int[] ad,
            final List<Allele> alleles) {
        return new GenotypeBuilder(name).PL(pl).AD(ad).alleles(alleles).make();
    }

    private static Genotype g(
            final String name, final int[] pl, final List<Allele> alleles) {
        return new GenotypeBuilder(name).PL(pl).alleles(alleles).make();
    }

    private static VariantContext vcWith(
            final VariantContext base,
            final List<Allele> alleles,
            final Genotype genotype) {
        return new VariantContextBuilder(base).alleles(alleles).genotypes(genotype).make();
    }
}
