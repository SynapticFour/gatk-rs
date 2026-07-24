import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;

/**
 * Rust {@code compute_core_variant_annotations} / Phase I {@code i1-core} gate.
 */
public final class HcParityCoreAnnotations {

    private HcParityCoreAnnotations() {}

    static final class CoreResult {
        final List<Integer> ac;
        final int an;
        final List<Double> af;
        final int ns;
        final int dp;

        CoreResult(final List<Integer> ac, final int an, final List<Double> af, final int ns, final int dp) {
            this.ac = ac;
            this.an = an;
            this.af = af;
            this.ns = ns;
            this.dp = dp;
        }
    }

    static final class SampleRow {
        final List<Integer> genotypeAlleles;
        final Integer dp;

        SampleRow(final List<Integer> genotypeAlleles, final Integer dp) {
            this.genotypeAlleles = genotypeAlleles;
            this.dp = dp;
        }
    }

    static CoreResult compute(final int altAlleleCount, final List<SampleRow> samples) {
        if (altAlleleCount <= 0) {
            throw new IllegalArgumentException("core annotations require at least one ALT allele");
        }
        final List<Integer> ac = new ArrayList<>();
        for (int i = 0; i < altAlleleCount; i++) {
            ac.add(0);
        }
        int an = 0;
        int ns = 0;
        int dp = 0;
        for (final SampleRow sample : samples) {
            if (sample.dp != null) {
                if (sample.dp < 0) {
                    throw new IllegalArgumentException("core annotations do not allow negative sample DP");
                }
                dp += sample.dp;
            }
            boolean sampleHasCalled = false;
            for (final int allele : sample.genotypeAlleles) {
                if (allele < 0) {
                    continue;
                }
                sampleHasCalled = true;
                an += 1;
                if (allele > 0) {
                    final int altIdx = allele - 1;
                    if (altIdx >= altAlleleCount) {
                        throw new IllegalArgumentException(
                                "allele index " + allele + " exceeds ALT allele count " + altAlleleCount);
                    }
                    ac.set(altIdx, ac.get(altIdx) + 1);
                }
            }
            if (sampleHasCalled) {
                ns += 1;
            }
        }
        final List<Double> af = new ArrayList<>();
        if (an > 0) {
            for (final int c : ac) {
                af.add(c / (double) an);
            }
        } else {
            for (int i = 0; i < altAlleleCount; i++) {
                af.add(0.0);
            }
        }
        return new CoreResult(ac, an, af, ns, dp);
    }

    static List<SampleRow> parseSamplesTsv(final Path path) throws IOException {
        final List<SampleRow> out = new ArrayList<>();
        for (final String line : Files.readAllLines(path)) {
            final String t = line.trim();
            if (t.isEmpty() || t.startsWith("#")) {
                continue;
            }
            final String[] cols = t.split("\t", -1);
            if (cols.length < 2) {
                throw new IllegalArgumentException("samples row needs >=2 cols: " + line);
            }
            final String[] gParts = cols[1].split(",");
            final List<Integer> gts = new ArrayList<>();
            for (final String p : gParts) {
                gts.add(Integer.parseInt(p.trim()));
            }
            Integer dp = null;
            if (cols.length >= 3 && !cols[2].isEmpty() && !"-".equals(cols[2])) {
                dp = Integer.parseInt(cols[2].trim());
            }
            out.add(new SampleRow(gts, dp));
        }
        return out;
    }

    static void dumpAnnotatedSite(final CoreResult site) {
        final String[] infoKeys = {"AC", "AN", "AF", "NS", "DP"};
        final String[] formatKeys = {"GT", "GQ", "AD", "DP", "PL"};
        System.out.println("info_key_count\t" + infoKeys.length);
        for (int i = 0; i < infoKeys.length; i++) {
            System.out.println("info_key_" + i + "\t" + infoKeys[i]);
        }
        System.out.println("format_key_count\t" + formatKeys.length);
        for (int i = 0; i < formatKeys.length; i++) {
            System.out.println("format_key_" + i + "\t" + formatKeys[i]);
        }
        System.out.println("ac\t" + formatIntList(site.ac));
        System.out.println("an\t" + site.an);
        System.out.println("af\t" + formatFloatList(site.af));
        System.out.println("ns\t" + site.ns);
        System.out.println("dp\t" + site.dp);
    }

    static void dumpAnnotationManifest() {
        final String[] infoKeys = {"AC", "AN", "AF", "NS", "DP"};
        final String[] formatKeys = {"GT", "GQ", "AD", "DP", "PL"};
        System.out.println("parity_v1_info_count\t" + infoKeys.length);
        for (int i = 0; i < infoKeys.length; i++) {
            System.out.println("parity_v1_info_" + i + "\t" + infoKeys[i]);
        }
        System.out.println("parity_v1_format_count\t" + formatKeys.length);
        for (int i = 0; i < formatKeys.length; i++) {
            System.out.println("parity_v1_format_" + i + "\t" + formatKeys[i]);
        }
    }

    private static String formatIntList(final List<Integer> v) {
        final StringBuilder sb = new StringBuilder();
        for (int i = 0; i < v.size(); i++) {
            if (i > 0) {
                sb.append(',');
            }
            sb.append(v.get(i));
        }
        return sb.toString();
    }

    private static String formatFloatList(final List<Double> v) {
        final StringBuilder sb = new StringBuilder();
        for (int i = 0; i < v.size(); i++) {
            if (i > 0) {
                sb.append(',');
            }
            sb.append(String.format(Locale.ROOT, "%.6f", v.get(i)));
        }
        return sb.toString();
    }
}
