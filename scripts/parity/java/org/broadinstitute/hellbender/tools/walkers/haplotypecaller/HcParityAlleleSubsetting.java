package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;

import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;

/**
 * Parity dump for GATK {@link HaplotypeCallerGenotypingEngine#whichAllelesToKeepBasedonHapScores}.
 */
public final class HcParityAlleleSubsetting {

    private HcParityAlleleSubsetting() {}

    public static void dumpAlleleSubsetting(
            final String sumsCsv, final String isRefCsv, final int maxAlleles) {
        final double[] sums = parseCsvDoubles(sumsCsv);
        final String[] refTok = isRefCsv.split(",");
        final int n = sums.length;
        final boolean[] isRef = parseIsRefFlags(refTok, n);
        final List<Integer> kept = keptIndices(sums, isRef, maxAlleles);
        int refIdx = 0;
        for (int i = 0; i < refTok.length; i++) {
            if ("1".equals(refTok[i].trim()) || "true".equalsIgnoreCase(refTok[i].trim())) {
                refIdx = i;
                break;
            }
        }
        int altIdx = refIdx;
        double altSum = Double.NEGATIVE_INFINITY;
        for (int i = 0; i < n; i++) {
            if (i != refIdx && sums[i] > altSum) {
                altSum = sums[i];
                altIdx = i;
            }
        }
        System.out.println("haplotype_count\t" + n);
        System.out.println("kept_indices\t" + joinInts(kept, ','));
        System.out.println("ref_haplotype_index\t" + refIdx);
        System.out.println("alt_haplotype_index\t" + altIdx);
    }

    /** Live assembly-region dump: per-hap scores + production trim indices (G-D05). */
    public static void dumpLiveSubsetExtension(
            final double[] sums,
            final boolean[] isRef,
            final List<Haplotype> haps,
            final int maxAlleles) {
        for (int i = 0; i < sums.length; i++) {
            System.out.println("haplotype_" + i + "_log10_sum\t" + sums[i]);
        }
        for (int i = 0; i < isRef.length; i++) {
            System.out.println("haplotype_" + i + "_is_reference\t" + isRef[i]);
        }
        for (int i = 0; i < haps.size(); i++) {
            System.out.println(
                    "haplotype_"
                            + i
                            + "_bases\t"
                            + new String(haps.get(i).getBases(), java.nio.charset.StandardCharsets.US_ASCII));
        }
        System.out.println("max_allele_count\t" + maxAlleles);
        final List<Integer> kept = keptIndices(sums, isRef, maxAlleles);
        System.out.println("kept_indices\t" + joinInts(kept, ','));
        System.out.println("trim_triggered\t" + (sums.length > maxAlleles));
    }

    private static List<Integer> keptIndices(
            final double[] sums, final boolean[] isRef, final int maxAlleles) {
        final int n = sums.length;
        final Allele[] alleles = new Allele[n];
        final LinkedHashMap<Allele, List<Haplotype>> alleleMapper = new LinkedHashMap<>();
        for (int i = 0; i < n; i++) {
            final char base = "ACGT".charAt(i % 4);
            alleles[i] = Allele.create(String.valueOf(base), isRef[i]);
            final Haplotype hap = new Haplotype(new byte[] {'A'});
            hap.setScore(sums[i]);
            alleleMapper.put(alleles[i], Collections.singletonList(hap));
        }
        final List<Allele> keptAlleles =
                HaplotypeCallerGenotypingEngine.whichAllelesToKeepBasedonHapScores(
                        alleleMapper, maxAlleles);
        final List<Integer> kept = new ArrayList<>();
        for (int i = 0; i < n; i++) {
            if (keptAlleles.contains(alleles[i])) {
                kept.add(i);
            }
        }
        Collections.sort(kept);
        return kept;
    }

    private static boolean[] parseIsRefFlags(final String[] refTok, final int n) {
        final boolean[] isRef = new boolean[n];
        for (int i = 0; i < n; i++) {
            isRef[i] =
                    i < refTok.length
                            && ("1".equals(refTok[i].trim())
                                    || "true".equalsIgnoreCase(refTok[i].trim()));
        }
        return isRef;
    }

    private static double[] parseCsvDoubles(final String csv) {
        if (csv == null || csv.isEmpty()) {
            return new double[0];
        }
        final String[] parts = csv.split(",");
        final double[] out = new double[parts.length];
        for (int i = 0; i < parts.length; i++) {
            out[i] = Double.parseDouble(parts[i].trim());
        }
        return out;
    }

    private static String joinInts(final List<Integer> values, final char sep) {
        final StringBuilder sb = new StringBuilder();
        for (int i = 0; i < values.size(); i++) {
            if (i > 0) {
                sb.append(sep);
            }
            sb.append(values.get(i));
        }
        return sb.toString();
    }
}
