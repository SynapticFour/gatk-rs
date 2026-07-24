import org.broadinstitute.hellbender.tools.walkers.genotyper.HcParityAfEm;

import java.util.List;
import java.util.Locale;

/**
 * Rust {@code region_vcf_emit} / Phase J VCF identity + FORMAT dumps.
 */
public final class HcParityRegionVcf {

    private HcParityRegionVcf() {}

    static void dumpVariantVcf(final boolean emitted, final String chrom, final long pos, final String id,
            final String ref, final String alt, final String qual, final String filter) {
        if (!emitted) {
            System.out.println("variant_emitted\tfalse");
            return;
        }
        System.out.println("variant_emitted\ttrue");
        System.out.println("chrom\t" + chrom);
        System.out.println("pos\t" + pos);
        System.out.println("id\t" + id);
        System.out.println("ref\t" + ref);
        System.out.println("alt\t" + alt);
        System.out.println("qual\t" + qual);
        System.out.println("filter\t" + filter);
    }

    static void dumpVariantFormat(final boolean emitted, final String gt, final int gq, final int dp,
            final String ad, final String pl) {
        if (!emitted) {
            System.out.println("format_emitted\tfalse");
            return;
        }
        System.out.println("format_emitted\ttrue");
        System.out.println("gt\t" + gt);
        System.out.println("gq\t" + gq);
        System.out.println("dp\t" + dp);
        System.out.println("ad\t" + ad);
        System.out.println("pl\t" + pl);
    }

    static boolean shouldEmitSynthetic(final double[] gl, final int[] ad) throws Exception {
        final P7GenotypeFieldsDump.Fields f = P7GenotypeFieldsDump.emit(gl, ad);
        return bestPlIndex(f.pl) != 0 && f.gq >= 10;
    }

    static P7GenotypeFieldsDump.Fields emitSynthetic(final double[] gl, final int[] ad) throws Exception {
        return P7GenotypeFieldsDump.emit(gl, ad);
    }

    private static int bestPlIndex(final List<Integer> pl) {
        int best = 0;
        for (int i = 1; i < pl.size(); i++) {
            if (pl.get(i) < pl.get(best)) {
                best = i;
            }
        }
        return best;
    }

    static void dumpSyntheticVcf(
            final String contig,
            final long pos,
            final String refAllele,
            final String altAllele,
            final double[] gl,
            final int[] ad) throws Exception {
        final P7GenotypeFieldsDump.Fields f = emitSynthetic(gl, ad);
        final List<Integer> pl = f.pl;
        if (!shouldEmitSynthetic(gl, ad)) {
            dumpVariantVcf(false, "", 0, "", "", "", "", "");
            return;
        }
        final double qualPhred = HcParityAfEm.siteQualPhredFromLog10Gl(gl);
        final String qual = String.format(Locale.ROOT, "%.6f", qualPhred);
        dumpVariantVcf(true, contig, pos, ".", refAllele, altAllele, qual, ".");
    }

    static void dumpSyntheticFormat(
            final String contig,
            final long pos,
            final String refAllele,
            final String altAllele,
            final double[] gl,
            final int[] ad) throws Exception {
        final P7GenotypeFieldsDump.Fields f = emitSynthetic(gl, ad);
        final List<Integer> pl = f.pl;
        if (!shouldEmitSynthetic(gl, ad)) {
            dumpVariantFormat(false, "", 0, 0, "", "");
            return;
        }
        final int best = bestPlIndex(pl);
        final String gt = best == 1 ? "0/1" : "1/1";
        final String adStr = ad[0] + "," + ad[1];
        final String plStr = pl.get(0) + "," + pl.get(1) + "," + pl.get(2);
        dumpVariantFormat(true, gt, f.gq, f.dp, adStr, plStr);
    }
}
