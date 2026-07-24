import java.lang.reflect.Method;
import java.util.Arrays;
import java.util.Locale;
import org.broadinstitute.hellbender.utils.pairhmm.Log10PairHMM;
import org.broadinstitute.hellbender.utils.pairhmm.PairHMM;

/**
 * GATK production {@link Log10PairHMM} (GKL/native when available) for hc-full-parity Phase F.2.
 */
public final class HcParityNativePairHmm {

    private static final Method COMPUTE_LOG10;

    static {
        try {
            COMPUTE_LOG10 =
                    PairHMM.class.getDeclaredMethod(
                            "computeReadLikelihoodGivenHaplotypeLog10",
                            byte[].class,
                            byte[].class,
                            byte[].class,
                            byte[].class,
                            byte[].class,
                            byte[].class,
                            boolean.class,
                            byte[].class);
            COMPUTE_LOG10.setAccessible(true);
        } catch (final ReflectiveOperationException e) {
            throw new ExceptionInInitializerError(e);
        }
    }

    private HcParityNativePairHmm() {}

    private static byte[] fill(final int n, final int value) {
        final byte[] out = new byte[n];
        Arrays.fill(out, (byte) value);
        return out;
    }

    public static void dumpCases(final java.nio.file.Path casesPath, final Appendable out)
            throws Exception {
        final Log10PairHMM hmm = new Log10PairHMM(true);
        try (java.io.BufferedReader br = java.nio.file.Files.newBufferedReader(casesPath)) {
            out.append("case_id\tlog10_likelihood\n");
            String line;
            while ((line = br.readLine()) != null) {
                final String t = line.trim();
                if (t.isEmpty() || t.startsWith("#")) {
                    continue;
                }
                final String[] c = t.split("\t");
                if (c.length < 5) {
                    throw new IllegalArgumentException("pairhmm cases row needs 5 cols: " + line);
                }
                final String caseId = c[0];
                final byte[] read = c[1].getBytes();
                final String[] qparts = c[2].split(",");
                final byte[] quals = new byte[qparts.length];
                for (int i = 0; i < qparts.length; i++) {
                    quals[i] = (byte) Integer.parseInt(qparts[i].trim());
                }
                final byte[] hap = c[4].getBytes();
                final byte[] insQuals = fill(read.length, 45);
                final byte[] delQuals = fill(read.length, 45);
                final byte[] gcps = fill(read.length, 10);
                hmm.initialize(read.length, hap.length);
                final double ll =
                        (double)
                                COMPUTE_LOG10.invoke(
                                        hmm,
                                        hap,
                                        read,
                                        quals,
                                        insQuals,
                                        delQuals,
                                        gcps,
                                        true,
                                        null);
                out.append(caseId)
                        .append('\t')
                        .append(formatLog10(ll))
                        .append('\n');
            }
        } finally {
            hmm.close();
        }
    }

    public static double pairhmmLog10Likelihood(
            final String readBases,
            final byte[] readBaseQuals,
            final int readMapq,
            final String haplotypeBases) {
        final byte[] read = readBases.getBytes();
        final byte[] hap = haplotypeBases.getBytes();
        final byte[] insQuals = fill(read.length, 45);
        final byte[] delQuals = fill(read.length, 45);
        final byte[] gcps = fill(read.length, 10);
        final Log10PairHMM hmm = new Log10PairHMM(true);
        try {
            hmm.initialize(read.length, hap.length);
            return (double)
                    COMPUTE_LOG10.invoke(
                            hmm,
                            hap,
                            read,
                            readBaseQuals,
                            insQuals,
                            delQuals,
                            gcps,
                            true,
                            null);
        } catch (final ReflectiveOperationException e) {
            throw new IllegalStateException("native PairHMM invoke failed", e);
        } finally {
            hmm.close();
        }
    }

    public static String formatLog10(final double ll) {
        return String.format(Locale.ROOT, "%.17g", ll);
    }
}
