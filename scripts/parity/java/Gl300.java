import htsjdk.variant.variantcontext.Allele;
import java.io.BufferedReader;
import java.io.FileReader;
import java.lang.reflect.Method;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Locale;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeAlleleCounts;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeLikelihoodCalculator;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeLikelihoodCalculators;
import org.broadinstitute.hellbender.utils.genotyper.AlleleList;
import org.broadinstitute.hellbender.utils.genotyper.LikelihoodMatrix;

/**
 * 6R.300: per-read genotype-likelihood contributions from the canonical
 * LOGLESS allele matrix, using GATK's GenotypeLikelihoodCalculator.
 */
public final class Gl300 {
    public static void main(final String[] args) throws Exception {
        if (args.length != 1) {
            throw new IllegalArgumentException("usage: Gl300 6r297_java_hap_ll.tsv");
        }
        final List<String> names = new ArrayList<String>();
        final List<Integer> flags = new ArrayList<Integer>();
        final List<double[]> rows = new ArrayList<double[]>();
        try (BufferedReader br = new BufferedReader(new FileReader(args[0]))) {
            String line;
            while ((line = br.readLine()) != null) {
                if (!line.startsWith("6R297\tallele_ll\t")) {
                    continue;
                }
                final String[] p = line.split("\t");
                names.add(p[2]);
                flags.add(Integer.valueOf(p[3]));
                rows.add(new double[] {bits(p[4]), bits(p[5]), bits(p[6])});
            }
        }
        final Allele[] alleles =
                new Allele[] {
                    Allele.create("T".getBytes(), true),
                    Allele.create("TTTG".getBytes(), false),
                    Allele.create("TGTTTG".getBytes(), false)
                };
        final Matrix matrix = new Matrix(alleles, rows);
        final GenotypeLikelihoodCalculator calc =
                new GenotypeLikelihoodCalculators().getInstance(2, 3);
        calc.ensureReadCapacity(rows.size());
        final Method components =
                GenotypeLikelihoodCalculator.class.getDeclaredMethod(
                        "readLikelihoodComponentsByAlleleCount", LikelihoodMatrix.class);
        components.setAccessible(true);
        final Method byRead =
                GenotypeLikelihoodCalculator.class.getDeclaredMethod(
                        "genotypeLikelihoodByRead", double[].class, int.class);
        byRead.setAccessible(true);
        final double[] component =
                (double[]) components.invoke(calc, matrix);
        final double[][] perRead =
                (double[][]) byRead.invoke(calc, component, Integer.valueOf(rows.size()));
        final StringBuilder order = new StringBuilder();
        for (int g = 0; g < calc.genotypeCount(); g++) {
            if (g > 0) {
                order.append(',');
            }
            final GenotypeAlleleCounts counts = calc.genotypeAlleleCountsAt(g);
            order.append(counts);
        }
        kv("genotype_order", order.toString());
        kv("n_reads", Integer.toString(rows.size()));
        kv("n_genotypes", Integer.toString(calc.genotypeCount()));
        for (int r = 0; r < rows.size(); r++) {
            final StringBuilder row = new StringBuilder();
            row.append(names.get(r)).append('\t').append(flags.get(r));
            for (int g = 0; g < calc.genotypeCount(); g++) {
                row.append('\t').append(fmt(perRead[g][r]));
            }
            kv("read_gl", row.toString());
        }
    }

    private static double bits(final String cell) {
        final String hex = cell.substring(cell.indexOf("bits=") + 5).trim();
        return Double.longBitsToDouble(Long.parseUnsignedLong(hex.substring(2), 16));
    }

    private static String fmt(final double x) {
        return String.format(Locale.US, "%.17g bits=0x%016x", x, Double.doubleToRawLongBits(x));
    }

    private static void kv(final String key, final String value) {
        System.out.println("6R300\t" + key + "\t" + value);
    }

    private static final class Matrix implements LikelihoodMatrix<String, Allele>, AlleleList<Allele> {
        private final Allele[] alleles;
        private final List<double[]> rows;

        Matrix(final Allele[] alleles, final List<double[]> rows) {
            this.alleles = alleles;
            this.rows = rows;
        }

        public List<String> evidence() {
            return null;
        }

        public List<Allele> alleles() {
            return Arrays.asList(alleles);
        }

        public void set(final int alleleIndex, final int evidenceIndex, final double value) {
            rows.get(evidenceIndex)[alleleIndex] = value;
        }

        public double get(final int alleleIndex, final int evidenceIndex) {
            return rows.get(evidenceIndex)[alleleIndex];
        }

        public int indexOfAllele(final Allele allele) {
            for (int i = 0; i < alleles.length; i++) {
                if (alleles[i].equals(allele)) {
                    return i;
                }
            }
            return -1;
        }

        public int indexOfEvidence(final String evidence) {
            return -1;
        }

        public int numberOfAlleles() {
            return alleles.length;
        }

        public int evidenceCount() {
            return rows.size();
        }

        public Allele getAllele(final int index) {
            return alleles[index];
        }

        public String getEvidence(final int index) {
            return Integer.toString(index);
        }

        public void copyAlleleLikelihoods(final int alleleIndex, final double[] dest, final int offset) {
            for (int r = 0; r < rows.size(); r++) {
                dest[offset + r] = rows.get(r)[alleleIndex];
            }
        }
    }
}
