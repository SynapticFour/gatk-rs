import java.util.Arrays;

/**
 * Scalar PairHMM matching {@code gatk_haplotypecaller::pairhmm::pairhmm_log10_likelihood} (default
 * {@code PairHmmParams}) for hc-full-parity L2. Not GATK {@code Log10PairHMM}; numerically aligned to
 * Rust reference implementation.
 */
public final class HcParityScalarPairHmm {

    private static final double LOG10_NEG_INF = Double.NEGATIVE_INFINITY;

    private static final double GAP_OPEN = 1e-2;
    private static final double GAP_EXTEND = 1e-1;
    private static final double INS_EMIT = 0.25;

    private HcParityScalarPairHmm() {}

    private static double probToLog10(final double p) {
        if (p <= 0.0) {
            return LOG10_NEG_INF;
        }
        return Math.log10(p);
    }

    private static double log10Sum(final double a, final double b) {
        if (Double.isInfinite(a) && a < 0) {
            return b;
        }
        if (Double.isInfinite(b) && b < 0) {
            return a;
        }
        final double m = Math.max(a, b);
        return m + Math.log10(Math.pow(10.0, a - m) + Math.pow(10.0, b - m));
    }

    private static double log10Sum3(final double a, final double b, final double c) {
        return log10Sum(log10Sum(a, b), c);
    }

    private static int effectiveQual(final int baseQual, final int mapq) {
        return Math.min(baseQual, mapq);
    }

    private static double emissionLog10(final byte readBase, final byte hapBase, final int effQ) {
        if (readBase == 'N' || readBase == 'n' || hapBase == 'N' || hapBase == 'n') {
            return probToLog10(0.25);
        }
        final double errorProb = Math.pow(10.0, -(double) effQ / 10.0);
        final char r = Character.toUpperCase((char) (readBase & 0xFF));
        final char h = Character.toUpperCase((char) (hapBase & 0xFF));
        if (r == h) {
            return probToLog10(1.0 - errorProb);
        }
        return probToLog10(errorProb / 3.0);
    }

    /** Log10 P(read | haplotype) with default Rust scaffold parameters. */
    public static double pairhmmLog10Likelihood(
            final String readBases,
            final byte[] quals,
            final int mapq,
            final String haplotypeBases) {
        final byte[] r = readBases.getBytes(java.nio.charset.StandardCharsets.US_ASCII);
        final byte[] h = haplotypeBases.getBytes(java.nio.charset.StandardCharsets.US_ASCII);
        final int rn = r.length;
        final int hn = h.length;
        if (rn == 0) {
            return 0.0;
        }
        if (quals.length != rn) {
            throw new IllegalArgumentException("quals length must match read bases");
        }
        if (hn == 0) {
            throw new IllegalArgumentException("haplotype must be non-empty");
        }

        final double[][] m = new double[rn + 1][hn + 1];
        final double[][] ins = new double[rn + 1][hn + 1];
        final double[][] del = new double[rn + 1][hn + 1];
        for (int i = 0; i <= rn; i++) {
            Arrays.fill(m[i], LOG10_NEG_INF);
            Arrays.fill(ins[i], LOG10_NEG_INF);
            Arrays.fill(del[i], LOG10_NEG_INF);
        }
        m[0][0] = 0.0;
        final double initDel = -Math.log10(hn);
        for (int j = 1; j <= hn; j++) {
            del[0][j] = initDel;
        }

        final double pGo = probToLog10(GAP_OPEN);
        final double pGe = probToLog10(GAP_EXTEND);
        final double pStay = probToLog10(Math.max(1e-12, 1.0 - 2.0 * GAP_OPEN));
        final double pInsEmit = probToLog10(Math.max(1e-12, INS_EMIT));

        for (int i = 1; i <= rn; i++) {
            final int q = effectiveQual(quals[i - 1] & 0xFF, mapq);
            for (int j = 1; j <= hn; j++) {
                final double e = emissionLog10(r[i - 1], h[j - 1], q);
                final double fromMatch = m[i - 1][j - 1] + pStay;
                final double fromIns = ins[i - 1][j - 1] + pGo;
                final double fromDel = del[i - 1][j - 1] + pGo;
                m[i][j] = log10Sum3(fromMatch, fromIns, fromDel) + e;

                final double insFromM = m[i - 1][j] + pGo;
                final double insFromI = ins[i - 1][j] + pGe;
                ins[i][j] = log10Sum(insFromM, insFromI) + pInsEmit;

                final double delFromM = m[i][j - 1] + pGo;
                final double delFromD = del[i][j - 1] + pGe;
                del[i][j] = log10Sum(delFromM, delFromD);
            }
        }

        double terminal = LOG10_NEG_INF;
        for (int j = 1; j <= hn; j++) {
            terminal = log10Sum(terminal, m[rn][j]);
            terminal = log10Sum(terminal, ins[rn][j]);
        }
        return terminal;
    }
}
