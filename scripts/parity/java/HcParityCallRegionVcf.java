import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Locale;
import org.broadinstitute.hellbender.engine.AssemblyRegion;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AssemblyResultSet;

/**
 * Emit VCF identity dump from one assembled active region (mirrors Rust {@code try_emit_call_region_variant}).
 */
public final class HcParityCallRegionVcf {

    private static final double STAND_EMIT_CONFIDENCE = 10.0;

    private HcParityCallRegionVcf() {}

    static void emitFromRegion(final AssemblyRegion r, final AssemblyResultSet ars) throws Exception {
        final List<GATKRead> reads = new ArrayList<>(r.getReads());
        HcParityRegionGenotype.sortReads(reads);
        final List<Haplotype> haps = new ArrayList<>(ars.getHaplotypeList());
        if (haps.isEmpty()) {
            HcParityRegionVcf.dumpVariantVcf(false, "", 0, "", "", "", "", "");
            return;
        }
        final List<double[]> readRows = new ArrayList<>();
        final boolean[] isRef = new boolean[haps.size()];
        for (int hi = 0; hi < haps.size(); hi++) {
            isRef[hi] = haps.get(hi).isReference();
        }
        for (final GATKRead read : reads) {
            final byte[] quals = read.getBaseQualities();
            final String readBases =
                    new String(read.getBases(), StandardCharsets.US_ASCII);
            final double[] row = new double[haps.size()];
            Arrays.fill(row, Double.NEGATIVE_INFINITY);
            for (int hi = 0; hi < haps.size(); hi++) {
                final String hapBases =
                        new String(haps.get(hi).getBases(), StandardCharsets.US_ASCII);
                row[hi] =
                        HcParityNativePairHmm.pairhmmLog10Likelihood(
                                readBases, quals, read.getMappingQuality(), hapBases);
            }
            readRows.add(row);
        }
        final HcParityRegionGenotype.GenotypeDump dump =
                HcParityRegionGenotype.genotypeFromLikelihoodMatrix(readRows, isRef);
        if (!dump.genotyped || dump.format == null) {
            HcParityRegionVcf.dumpVariantVcf(false, "", 0, "", "", "", "", "");
            return;
        }
        if (dump.format.gq < STAND_EMIT_CONFIDENCE) {
            HcParityRegionVcf.dumpVariantVcf(false, "", 0, "", "", "", "", "");
            return;
        }
        final int refIdx = dump.refHapIdx;
        final int altIdx = dump.altHapIdx;
        final byte[] refBases = haps.get(refIdx).getBases();
        final byte[] altBases = haps.get(altIdx).getBases();
        if (Arrays.equals(refBases, altBases)) {
            HcParityRegionVcf.dumpVariantVcf(false, "", 0, "", "", "", "", "");
            return;
        }
        final Long pos = firstDifferingPosition(r.getStart(), refBases, altBases);
        if (pos == null) {
            HcParityRegionVcf.dumpVariantVcf(false, "", 0, "", "", "", "", "");
            return;
        }
        final String refAllele = alleleAt(refBases, pos, r.getStart());
        final String altAllele = alleleAt(altBases, pos, r.getStart());
        if (refAllele.equals(altAllele)) {
            HcParityRegionVcf.dumpVariantVcf(false, "", 0, "", "", "", "", "");
            return;
        }
        final int bestPl = bestPlIndex(dump.format.pl);
        if (bestPl == 0) {
            HcParityRegionVcf.dumpVariantVcf(false, "", 0, "", "", "", "", "");
            return;
        }
        final String qual =
                dump.format.pl.get(bestPl) <= 0
                        ? "99.000000"
                        : String.format(
                                Locale.ROOT,
                                "%.6f",
                                (double) Math.min(99, dump.format.pl.get(bestPl)));
        HcParityRegionVcf.dumpVariantVcf(
                true, r.getContig(), pos, ".", refAllele, altAllele, qual, ".");
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

    private static Long firstDifferingPosition(
            final int regionStart1, final byte[] refBases, final byte[] altBases) {
        final int max = Math.max(refBases.length, altBases.length);
        for (int i = 0; i < max; i++) {
            final byte rb = i < refBases.length ? refBases[i] : (byte) 'N';
            final byte ab = i < altBases.length ? altBases[i] : (byte) 'N';
            if (rb != ab) {
                return (long) regionStart1 + i;
            }
        }
        return null;
    }

    private static String alleleAt(
            final byte[] hapBases, final long pos1, final int regionStart1) {
        final int off = (int) (pos1 - regionStart1);
        if (off < 0 || off >= hapBases.length) {
            return "N";
        }
        return new String(new byte[] {hapBases[off]}, StandardCharsets.US_ASCII);
    }
}
