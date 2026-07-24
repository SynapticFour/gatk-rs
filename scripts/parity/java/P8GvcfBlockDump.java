import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;

/**
 * Standalone reference implementation matching {@code build_gvcf_blocks_with_semantics}
 * (phase 8 steps 101-103): band bucketization, adjacency merge, RGQ delta guard.
 *
 * Not invoked by production Rust - used only as an independent parity oracle under Docker.
 */
public final class P8GvcfBlockDump {

    private static final int MAX_RGQ_DELTA_WITHIN_BLOCK = 10;

    static final class Locus {
        final int position1Based;
        final int gq;
        final int dp;

        Locus(int position1Based, int gq, int dp) {
            this.position1Based = position1Based;
            this.gq = gq;
            this.dp = dp;
        }
    }

    static final class Block {
        int start1Based;
        int end1Based;
        int gqBandUpper;
        int minRgq;
        int minDp;
        int maxDp;
    }

    private static int gqBandUpper(int gq, int[] bands) {
        for (int b : bands) {
            if (gq <= b) {
                return b;
            }
        }
        return bands[bands.length - 1];
    }

    public static List<Block> buildBlocks(List<Locus> loci, int[] bands, int maxRgqDelta) {
        List<Block> out = new ArrayList<>();
        if (loci.isEmpty()) {
            return out;
        }
        Locus first = loci.get(0);
        Block cur = new Block();
        cur.start1Based = first.position1Based;
        cur.end1Based = first.position1Based;
        cur.gqBandUpper = gqBandUpper(first.gq, bands);
        cur.minRgq = first.gq;
        cur.minDp = first.dp;
        cur.maxDp = first.dp;

        for (int i = 1; i < loci.size(); i++) {
            Locus locus = loci.get(i);
            int band = gqBandUpper(locus.gq, bands);
            boolean contiguous = locus.position1Based == cur.end1Based + 1;
            boolean sameBand = band == cur.gqBandUpper;
            boolean rgqCompatible =
                    Math.abs(locus.gq - cur.minRgq) <= maxRgqDelta;
            if (contiguous && sameBand && rgqCompatible) {
                cur.end1Based = locus.position1Based;
                cur.minRgq = Math.min(cur.minRgq, locus.gq);
                cur.minDp = Math.min(cur.minDp, locus.dp);
                cur.maxDp = Math.max(cur.maxDp, locus.dp);
            } else {
                out.add(cur);
                cur = new Block();
                cur.start1Based = locus.position1Based;
                cur.end1Based = locus.position1Based;
                cur.gqBandUpper = band;
                cur.minRgq = locus.gq;
                cur.minDp = locus.dp;
                cur.maxDp = locus.dp;
            }
        }
        out.add(cur);
        return out;
    }

    /** Mirrors {@code gvcf_block_to_record_fields}: END == inclusive block end. */
    static String row(String caseId, Block b) {
        return String.format(
                Locale.ROOT,
                "%s\t%d\t%d\t%d\t%d\t%d\t%d",
                caseId,
                b.start1Based,
                b.end1Based,
                b.gqBandUpper,
                b.minRgq,
                b.minDp,
                b.maxDp);
    }

    public static void main(String[] args) throws IOException {
        if (args.length != 1) {
            System.err.println("usage: P8GvcfBlockDump <fixture.tsv>");
            System.exit(2);
        }
        Path fixture = Path.of(args[0]);
        List<String> lines = Files.readAllLines(fixture, StandardCharsets.UTF_8);

        Map<String, List<Locus>> byCase = new LinkedHashMap<>();
        for (String raw : lines) {
            String line = raw.trim();
            if (line.isEmpty() || line.startsWith("#")) {
                continue;
            }
            String[] c = line.split("\t");
            if (c.length != 4) {
                throw new IllegalArgumentException("bad fixture row: " + line);
            }
            String caseId = c[0];
            int pos = Integer.parseInt(c[1]);
            int gq = Integer.parseInt(c[2]);
            int dp = Integer.parseInt(c[3]);
            byCase.computeIfAbsent(caseId, k -> new ArrayList<>()).add(new Locus(pos, gq, dp));
        }

        int[] bands = new int[] {9, 19, 29, 99};

        System.out.println("# case_id\tstart_1based\tend_info\tgq_band_upper\tmin_rgq\tmin_dp\tmax_dp");
        for (Map.Entry<String, List<Locus>> e : byCase.entrySet()) {
            List<Locus> loci = e.getValue();
            for (int i = 1; i < loci.size(); i++) {
                if (loci.get(i - 1).position1Based >= loci.get(i).position1Based) {
                    throw new IllegalArgumentException(
                            "case " + e.getKey() + ": positions must be strictly increasing");
                }
            }
            List<Block> blocks = buildBlocks(loci, bands, MAX_RGQ_DELTA_WITHIN_BLOCK);
            for (Block b : blocks) {
                System.out.println(row(e.getKey(), b));
            }
        }
    }
}
