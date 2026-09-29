import java.io.BufferedReader;
import java.io.FileReader;
import java.lang.reflect.Field;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;
import org.broadinstitute.hellbender.utils.MathUtils;

/** 6R.301: frozen GATK 4.4.0.0 heterozygote combine oracle. */
public final class Gl301 {
    public static void main(final String[] args) throws Exception {
        if (args.length != 1) {
            throw new IllegalArgumentException("usage: Gl301 6r297_java_hap_ll.tsv");
        }
        final Class<?> jlt =
                Class.forName("org.broadinstitute.hellbender.utils.MathUtils$JacobianLogTable");
        final Field cacheField = jlt.getDeclaredField("cache");
        cacheField.setAccessible(true);
        final double[] cache = (double[]) cacheField.get(null);
        final Field inv = jlt.getDeclaredField("INV_STEP");
        inv.setAccessible(true);
        kv("cache_len", Integer.toString(cache.length));
        kv("inv_step", fmt(inv.getDouble(null)));
        long xor = 0L;
        for (double value : cache) {
            xor ^= Double.doubleToRawLongBits(value);
        }
        kv("cache_xor", String.format("0x%016x", xor));
        for (int k : new int[] {0, 1, 2, 3, 10000, 79999, 80000}) {
            kv("cache_k", k + "\t" + fmt(cache[k]));
        }
        final String[] genotypes = new String[] {"T/TTTG", "T/TGTTTG", "TTTG/TGTTTG"};
        final int[][] pairs = new int[][] {{0, 1}, {0, 2}, {1, 2}};
        try (BufferedReader br = new BufferedReader(new FileReader(args[0]))) {
            String line;
            while ((line = br.readLine()) != null) {
                if (!line.startsWith("6R297\tallele_ll\t")) {
                    continue;
                }
                final String[] p = line.split("\t");
                final double[] alleles =
                        new double[] {bits(p[4]), bits(p[5]), bits(p[6])};
                for (int g = 0; g < pairs.length; g++) {
                    final double a = alleles[pairs[g][0]];
                    final double b = alleles[pairs[g][1]];
                    kv(
                            "pair",
                            genotypes[g]
                                    + "\t"
                                    + p[2]
                                    + "\t"
                                    + p[3]
                                    + "\t"
                                    + fmt(a)
                                    + "\t"
                                    + fmt(b)
                                    + "\t"
                                    + describe(a, b, cache));
                }
            }
        }
        final String[] labels =
                new String[] {
                    "0",
                    "1e-15",
                    "1e-12",
                    "1e-9",
                    "1",
                    "4",
                    "7.999999999999",
                    "8.0",
                    "8.000000000001",
                    "10"
                };
        for (String label : labels) {
            final double diff = Double.parseDouble(label);
            kv("threshold", label + "\t" + fmt(diff) + "\t" + describe(0.0, diff, cache));
            kv("threshold_swapped", label + "\t" + fmt(diff) + "\t" + describe(diff, 0.0, cache));
        }
    }

    private static String describe(final double a, final double b, final double[] cache) {
        final double result = MathUtils.approximateLog10SumLog10(a, b);
        final double smaller = a > b ? b : a;
        final double larger = a > b ? a : b;
        final String branch;
        final String index;
        if (smaller == Double.NEGATIVE_INFINITY) {
            branch = "neg_inf";
            index = "NA";
        } else {
            final double diff = larger - smaller;
            if (diff < 8.0) {
                branch = "table";
                final int idx = MathUtils.fastRound(diff * 10000.0);
                index = idx + "\t" + fmt(cache[idx]);
            } else {
                branch = "ignore";
                index = "NA";
            }
        }
        return fmt(result) + "\t" + branch + "\t" + index;
    }

    private static double bits(final String cell) {
        final String hex = cell.substring(cell.indexOf("bits=") + 5).trim();
        return Double.longBitsToDouble(Long.parseUnsignedLong(hex.substring(2), 16));
    }

    private static String fmt(final double x) {
        return String.format(Locale.US, "%.17g bits=0x%016x", x, Double.doubleToRawLongBits(x));
    }

    private static void kv(final String key, final String value) {
        System.out.println("6R301\t" + key + "\t" + value);
    }
}
