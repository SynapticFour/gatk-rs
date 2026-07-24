import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Locale;

/**
 * Standalone reference for {@code emit_genotype_format_fields}: normalized PLs (Phred vs max LL),
 * GQ from {@link htsjdk.variant.variantcontext.GenotypeLikelihoods#getGQLog10FromLikelihoods}, AD + summed DP.
 *
 * Mirrors gatk-rs rounding ({@code f64::round} semantics).
 */
public final class P7GenotypeFieldsDump {

    private static int expectedDiploidGenotypeCount(int alleleCount) {
        return alleleCount * (alleleCount + 1) / 2;
    }

    static final class Fields {
        final List<Integer> pl;
        final int gq;
        final List<Integer> ad;
        final int dp;

        Fields(List<Integer> pl, int gq, List<Integer> ad, int dp) {
            this.pl = pl;
            this.gq = gq;
            this.ad = ad;
            this.dp = dp;
        }
    }

    static int gqPhredFromLog10(final double[] gl, final int bestIdx) {
        final double log10PError;
        if (gl.length == 3) {
            log10PError =
                    switch (bestIdx) {
                        case 0 ->
                                org.broadinstitute.hellbender.utils.MathUtils.log10SumLog10(
                                        gl[1], gl[2]);
                        case 1 ->
                                org.broadinstitute.hellbender.utils.MathUtils.log10SumLog10(
                                        gl[0], gl[2]);
                        default ->
                                org.broadinstitute.hellbender.utils.MathUtils.log10SumLog10(
                                        gl[0], gl[1]);
                    };
        } else {
            log10PError = 0.0;
        }
        return (int)
                Math.round(
                        Math.max(
                                0.0,
                                Math.min(99.0, -10.0 * Math.min(0.0, log10PError))));
    }

    /** HTSJDK FORMAT GQ: second-best PL minus best PL; 0 when multiple genotypes tie at min PL. */
    static int gqFromPl(final List<Integer> pl) {
        if (pl.isEmpty()) {
            return 0;
        }
        int minPl = pl.stream().mapToInt(Integer::intValue).min().orElse(0);
        long atMin = pl.stream().filter(p -> p == minPl).count();
        if (atMin > 1) {
            return 0;
        }
        int second =
                pl.stream().filter(p -> p > minPl).mapToInt(Integer::intValue).min().orElse(minPl);
        return Math.min(99, second - minPl);
    }

    static Fields emit(double[] genotypeLog10Likelihoods, int[] alleleDepths) {
        if (genotypeLog10Likelihoods.length == 0) {
            throw new IllegalArgumentException("empty genotype likelihoods");
        }
        if (alleleDepths.length == 0) {
            throw new IllegalArgumentException("empty allele depths");
        }
        int expectedGl = expectedDiploidGenotypeCount(alleleDepths.length);
        if (genotypeLog10Likelihoods.length != expectedGl) {
            throw new IllegalArgumentException(
                    "genotype likelihood count mismatch: got "
                            + genotypeLog10Likelihoods.length
                            + " expected "
                            + expectedGl);
        }
        for (double v : genotypeLog10Likelihoods) {
            if (!Double.isFinite(v)) {
                throw new IllegalArgumentException("non-finite genotype likelihood");
            }
        }
        for (int d : alleleDepths) {
            if (d < 0) {
                throw new IllegalArgumentException("negative allele depth");
            }
        }

        double maxLl =
                Arrays.stream(genotypeLog10Likelihoods).max().orElse(Double.NEGATIVE_INFINITY);
        List<Integer> pls = new ArrayList<>(genotypeLog10Likelihoods.length);
        for (double ll : genotypeLog10Likelihoods) {
            int ph = (int) Math.round(-10.0 * (ll - maxLl));
            pls.add(Math.max(0, ph));
        }
        int minPl = pls.stream().mapToInt(Integer::intValue).min().orElse(0);
        for (int i = 0; i < pls.size(); i++) {
            pls.set(i, pls.get(i) - minPl);
        }

        int bestIdx = 0;
        double bestLl = Double.NEGATIVE_INFINITY;
        for (int i = 0; i < genotypeLog10Likelihoods.length; i++) {
            if (genotypeLog10Likelihoods[i] > bestLl) {
                bestLl = genotypeLog10Likelihoods[i];
                bestIdx = i;
            }
        }
        if (genotypeLog10Likelihoods.length >= 3 && alleleDepths.length >= 2) {
            int refD = Math.max(0, alleleDepths[0]);
            int altD = Math.max(0, alleleDepths[1]);
            double hetGl = genotypeLog10Likelihoods[1];
            double homAltGl = genotypeLog10Likelihoods[2];
            if (altD > refD + 1 && Math.abs(hetGl - homAltGl) < 1e-3) {
                bestIdx = 2;
            }
        }
        int gq =
                (genotypeLog10Likelihoods.length == 3 && bestIdx == 0)
                        ? gqFromPl(pls)
                        : gqPhredFromLog10(genotypeLog10Likelihoods, bestIdx);

        int dp = Arrays.stream(alleleDepths).sum();

        List<Integer> adList = new ArrayList<>(alleleDepths.length);
        for (int d : alleleDepths) {
            adList.add(d);
        }

        return new Fields(pls, gq, adList, dp);
    }

    static String joinInts(List<Integer> xs, char sep) {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < xs.size(); i++) {
            if (i > 0) {
                sb.append(sep);
            }
            sb.append(xs.get(i));
        }
        return sb.toString();
    }

    static String fmtRow(String caseId, Fields f) {
        return String.format(
                Locale.ROOT,
                "%s\t%s\t%d\t%s\t%d",
                caseId,
                joinInts(f.pl, ','),
                f.gq,
                joinInts(f.ad, ','),
                f.dp);
    }

    public static void main(String[] args) throws IOException {
        if (args.length != 1) {
            System.err.println("usage: P7GenotypeFieldsDump <fixture.tsv>");
            System.exit(2);
        }
        Path fixture = Path.of(args[0]);
        List<String> lines = Files.readAllLines(fixture, StandardCharsets.UTF_8);

        System.out.println("# case_id\tpl\tgq\tad\tdp");
        for (String raw : lines) {
            String line = raw.trim();
            if (line.isEmpty() || line.startsWith("#")) {
                continue;
            }
            String[] c = line.split("\t");
            if (c.length != 3) {
                throw new IllegalArgumentException("bad fixture row: " + line);
            }
            String caseId = c[0];
            double[] gl = parseCsvDoubles(c[1]);
            int[] ad = parseCsvInts(c[2]);
            Fields f = emit(gl, ad);
            System.out.println(fmtRow(caseId, f));
        }
    }

    private static double[] parseCsvDoubles(String raw) {
        String[] parts = raw.split(",");
        double[] out = new double[parts.length];
        for (int i = 0; i < parts.length; i++) {
            out[i] = Double.parseDouble(parts[i]);
        }
        return out;
    }

    private static int[] parseCsvInts(String raw) {
        String[] parts = raw.split(",");
        int[] out = new int[parts.length];
        for (int i = 0; i < parts.length; i++) {
            out[i] = Integer.parseInt(parts[i]);
        }
        return out;
    }
}
