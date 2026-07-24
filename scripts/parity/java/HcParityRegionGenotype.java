import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import java.util.Locale;

/**
 * Rust {@code hc_genotyping_engine} + {@code call_region} genotype slice for hc-full-parity G.2.
 *
 * <p>Uses {@link HcParityNativePairHmm} (aligned with Rust {@code log10_pairhmm_likelihood} defaults:
 * INS/DEL/GCP 45/45/10). PCR error adjustment is not applied in this parity helper yet.
 */
public final class HcParityRegionGenotype {

    private HcParityRegionGenotype() {}

    static final class GenotypeDump {
        final int haplotypeCount;
        final int readCount;
        final boolean genotyped;
        final int refHapIdx;
        final int altHapIdx;
        final int bestHapIdx;
        final double[] genotypeLog10;
        final P7GenotypeFieldsDump.Fields format;

        GenotypeDump(
                final int haplotypeCount,
                final int readCount,
                final boolean genotyped,
                final int refHapIdx,
                final int altHapIdx,
                final int bestHapIdx,
                final double[] genotypeLog10,
                final P7GenotypeFieldsDump.Fields format) {
            this.haplotypeCount = haplotypeCount;
            this.readCount = readCount;
            this.genotyped = genotyped;
            this.refHapIdx = refHapIdx;
            this.altHapIdx = altHapIdx;
            this.bestHapIdx = bestHapIdx;
            this.genotypeLog10 = genotypeLog10;
            this.format = format;
        }
    }

    static GenotypeDump genotypeFromLikelihoodMatrix(
            final List<double[]> readRows, final boolean[] isReference) {
        final int hapCount = isReference.length;
        final int readCount = readRows.size();
        if (hapCount == 0 || readCount == 0) {
            return new GenotypeDump(hapCount, readCount, false, 0, 0, 0, new double[0], null);
        }
        final double[] sums = new double[hapCount];
        for (final double[] row : readRows) {
            for (int i = 0; i < hapCount; i++) {
                sums[i] += row[i];
            }
        }
        int best = 0;
        for (int i = 1; i < hapCount; i++) {
            if (sums[i] > sums[best]) {
                best = i;
            }
        }
        int refIdx = 0;
        for (int i = 0; i < hapCount; i++) {
            if (isReference[i]) {
                refIdx = i;
                break;
            }
        }
        int altIdx = refIdx;
        double altSum = Double.NEGATIVE_INFINITY;
        for (int i = 0; i < hapCount; i++) {
            if (i != refIdx && !isReference[i] && sums[i] > altSum) {
                altSum = sums[i];
                altIdx = i;
            }
        }
        double g0 = 0.0;
        double g1 = 0.0;
        double g2 = 0.0;
        int refAd = 0;
        int altAd = 0;
        final double log10Ploidy = Math.log10(2.0);
        final double denominator = readCount * log10Ploidy;
        for (final double[] row : readRows) {
            final double lr = row[refIdx];
            final double la = row[altIdx];
            g0 += lr + log10Ploidy;
            g2 += la + log10Ploidy;
            g1 += org.broadinstitute.hellbender.utils.MathUtils.log10SumLog10(lr, la);
            if (lr >= la) {
                refAd++;
            } else {
                altAd++;
            }
        }
        if (refAd == 0 && altAd == 0) {
            refAd = 1;
        }
        final double[] gls =
                new double[] {g0 - denominator, g1 - denominator, g2 - denominator};
        final P7GenotypeFieldsDump.Fields format =
                P7GenotypeFieldsDump.emit(gls, new int[] {refAd, altAd});
        return new GenotypeDump(
                hapCount, readCount, true, refIdx, altIdx, best, gls, format);
    }

    static void printDump(final GenotypeDump d) {
        System.out.println("haplotype_count\t" + d.haplotypeCount);
        System.out.println("read_count\t" + d.readCount);
        System.out.println("genotyped\t" + d.genotyped);
        if (!d.genotyped || d.format == null) {
            return;
        }
        System.out.println("ref_haplotype_index\t" + d.refHapIdx);
        System.out.println("alt_haplotype_index\t" + d.altHapIdx);
        System.out.println("best_haplotype_index\t" + d.bestHapIdx);
        for (int i = 0; i < d.genotypeLog10.length; i++) {
            System.out.printf(
                    Locale.ROOT, "genotype_%d_log10\t%s%n", i, Double.toString(d.genotypeLog10[i]));
        }
        System.out.println("pl\t" + P7GenotypeFieldsDump.joinInts(d.format.pl, ','));
        System.out.println("gq\t" + d.format.gq);
        System.out.println("ad\t" + P7GenotypeFieldsDump.joinInts(d.format.ad, ','));
        System.out.println("dp\t" + d.format.dp);
    }

    static void sortReads(final List<org.broadinstitute.hellbender.utils.read.GATKRead> reads) {
        reads.sort(
                Comparator.comparing(
                                (org.broadinstitute.hellbender.utils.read.GATKRead read) ->
                                        read.getName())
                        .thenComparingInt(
                                org.broadinstitute.hellbender.utils.read.GATKRead::getStart));
    }
}
