import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.tools.walkers.annotator.TandemRepeat;
import org.broadinstitute.hellbender.utils.SimpleInterval;

/**
 * 6R.292: arguments and result of {@code TandemRepeat.getNumTandemRepeatUnits}
 * for A&gt;AT at 20:29455644. Stops at the repeat count.
 */
public final class Repeat292 {
    public static void main(final String[] args) throws Exception {
        if (args.length != 1) {
            throw new IllegalArgumentException("usage: Repeat292 ref");
        }
        final Path ref = Paths.get(args[0]);
        final ReferenceDataSource reference = ReferenceDataSource.of(ref);
        final SimpleInterval active = new SimpleInterval("20", 29455560, 29455744);
        final ReferenceContext refCtx = new ReferenceContext(reference, active, 100, 100);
        final VariantContext vc =
                new VariantContextBuilder("6R292", "20", 29455644, 29455644, Arrays.asList(
                                Allele.create("A", true), Allele.create("AT", false)))
                        .make();
        final SimpleInterval window = refCtx.getWindow();
        kv("window", window.getContig() + ":" + window.getStart() + "-" + window.getEnd());
        final byte[] refBases = refCtx.getBases();
        final int anchor = 29455644 - window.getStart();
        kv("anchor_base", String.valueOf((char) refBases[anchor]));
        kv("ref_at_event_and_after", new String(refBases, anchor, 16));
        final int startIndex = vc.getStart() + 1 - window.getStart();
        kv("start_index", Integer.toString(startIndex));
        kv("java_context_start", Integer.toString(window.getStart() + startIndex));
        final byte[] remaining = Arrays.copyOfRange(refBases, startIndex, Math.min(refBases.length, startIndex + 16));
        kv("java_remaining_prefix", new String(remaining));
        final byte[] refAlleleBases = Arrays.copyOfRange(vc.getReference().getBases(), 1, vc.getReference().length());
        final byte[] altBases = Arrays.copyOfRange(vc.getAlternateAllele(0).getBases(), 1, vc.getAlternateAllele(0).length());
        kv("ref_allele", vc.getReference().getDisplayString());
        kv("alt_allele", vc.getAlternateAllele(0).getDisplayString());
        kv("ref_allele_after_anchor", new String(refAlleleBases));
        kv("alt_after_anchor", new String(altBases));
        final Pair<java.util.List<Integer>, byte[]> rep = TandemRepeat.getNumTandemRepeatUnits(refCtx, vc);
        if (rep == null || rep.getRight() == null) {
            kv("repeat_unit", "NONE");
            kv("repeat_counts", "NONE");
            return;
        }
        kv("repeat_unit", new String(rep.getRight()));
        kv("repeat_unit_len", Integer.toString(rep.getRight().length));
        kv("repeat_counts", rep.getLeft().toString());
        int most = 0;
        for (final int n : rep.getLeft()) {
            most = Math.max(most, n);
        }
        kv("repeat_count_max", Integer.toString(most));
        kv("str_addition", Integer.toString(most * rep.getRight().length));
    }

    private static void kv(final String key, final String value) {
        System.out.println("6R292\t" + key + "\t" + value);
    }
}
