import htsjdk.samtools.SAMFileHeader;
import java.lang.reflect.Field;
import java.lang.reflect.Method;
import java.util.List;
import java.util.Map;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCallerEngine;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.PairHMMLikelihoodCalculationEngine;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.ReadLikelihoodCalculationEngine;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.genotyper.IndexedAlleleList;
import org.broadinstitute.hellbender.utils.genotyper.LikelihoodMatrix;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.pairhmm.PairHMM;
import org.broadinstitute.hellbender.utils.read.GATKRead;

/**
 * TEST-ONLY wrapper: dump {@code filterPoorlyModeledEvidence} inputs at the Java call, then
 * delegate the production filter. Does not change keep/drop arithmetic.
 */
public final class HcParityFilterPoorlyModeledDump implements ReadLikelihoodCalculationEngine {

    private static int callSeq = 0;
    /** Default 6R.96 post-kernel dump. 6R.98 sets {@code 6R98}. */
    public static String DUMP_PREFIX = "6R96";
    /** 6R.98: GATK native PairHMM double-precision requested before engine construction. */
    public static boolean USE_DOUBLE_PRECISION = false;
    private final ReadLikelihoodCalculationEngine inner;

    public static void installOn(final HaplotypeCallerEngine engine) throws Exception {
        final Field f =
                HaplotypeCallerEngine.class.getDeclaredField("likelihoodCalculationEngine");
        f.setAccessible(true);
        final ReadLikelihoodCalculationEngine inner =
                (ReadLikelihoodCalculationEngine) f.get(engine);
        f.set(engine, new HcParityFilterPoorlyModeledDump(inner));
    }

    private HcParityFilterPoorlyModeledDump(final ReadLikelihoodCalculationEngine inner) {
        this.inner = inner;
    }

    @Override
    public void close() {
        inner.close();
    }

    @Override
    @SuppressWarnings("unchecked")
    public AlleleLikelihoods<GATKRead, Haplotype> computeReadLikelihoods(
            final List<Haplotype> haplotypeList,
            final SAMFileHeader hdr,
            final SampleList samples,
            final Map<String, List<GATKRead>> perSampleReadList,
            final boolean filterPoorly) {
        if (!(inner instanceof PairHMMLikelihoodCalculationEngine)) {
            return inner.computeReadLikelihoods(
                    haplotypeList, hdr, samples, perSampleReadList, filterPoorly);
        }
        try {
            return computeAndDump((PairHMMLikelihoodCalculationEngine) inner, haplotypeList, samples, perSampleReadList);
        } catch (final Exception e) {
            throw new RuntimeException("6R.93 filter-poorly-modeled dump failed", e);
        }
    }

    private static AlleleLikelihoods<GATKRead, Haplotype> computeAndDump(
            final PairHMMLikelihoodCalculationEngine phmm,
            final List<Haplotype> haplotypeList,
            final SampleList samples,
            final Map<String, List<GATKRead>> perSampleReadList)
            throws Exception {
        final Method init =
                PairHMMLikelihoodCalculationEngine.class.getDeclaredMethod(
                        "initializePairHMM", List.class, Map.class);
        init.setAccessible(true);
        init.invoke(phmm, haplotypeList, perSampleReadList);

        final AlleleLikelihoods<GATKRead, Haplotype> result =
                new AlleleLikelihoods<>(
                        samples, new IndexedAlleleList<>(haplotypeList), perSampleReadList);
        final Method computeOne =
                PairHMMLikelihoodCalculationEngine.class.getDeclaredMethod(
                        "computeReadLikelihoods", LikelihoodMatrix.class);
        computeOne.setAccessible(true);
        for (int i = 0; i < result.numberOfSamples(); i++) {
            computeOne.invoke(phmm, result.sampleMatrix(i));
        }

        // 6R.97: primitive GKL/PairHMM double[] vs AlleleLikelihoods after SampleMatrix.set.
        dumpKernelResultBuffer(phmm, result);

        // 6R.96: materialize immediately after PairHMM.computeLog10Likelihoods,
        // before normalizeLikelihoods / filterPoorlyModeledEvidence.
        dumpLikelihoodMatrix(result, "post_kernel");

        final double log10global = fieldDouble(phmm, "log10globalReadMismappingRate");
        final boolean sym = fieldBoolean(phmm, "symmetricallyNormalizeAllelesToReference");
        result.normalizeLikelihoods(log10global, sym);

        final boolean dynamic = fieldBoolean(phmm, "dynamicDisqualification");
        final double expectedError = fieldDouble(phmm, "expectedErrorRatePerBase");
        final double scale = fieldDouble(phmm, "readDisqualificationScale");
        dumpBeforeFilter(result, dynamic, expectedError, scale);
        phmm.filterPoorlyModeledEvidence(result, dynamic, expectedError, scale);
        return result;
    }

    /**
     * 6R.97: compare PairHMM {@code mLogLikelihoodArray} (GKL {@code double[]}) to
     * {@code AlleleLikelihoods} after {@code SampleMatrix.set}. Assignment copy only:
     * no numeric transform. Hap-list order may permute cells; compare sorted bit
     * multisets.
     */
    private static void dumpKernelResultBuffer(
            final PairHMMLikelihoodCalculationEngine phmm,
            final AlleleLikelihoods<GATKRead, Haplotype> ll)
            throws Exception {
        final int seq = callSeq + 1;
        final Field pf =
                PairHMMLikelihoodCalculationEngine.class.getDeclaredField("pairHMM");
        pf.setAccessible(true);
        final PairHMM hmm = (PairHMM) pf.get(phmm);
        final double[] buf = hmm.getLogLikelihoodArray();
        int matrixN = 0;
        int matrixF32 = 0;
        final java.util.ArrayList<Long> matrixBits = new java.util.ArrayList<>();
        for (int s = 0; s < ll.numberOfSamples(); s++) {
            final LikelihoodMatrix<GATKRead, Haplotype> mx = ll.sampleMatrix(s);
            for (int r = 0; r < mx.evidenceCount(); r++) {
                for (int a = 0; a < mx.numberOfAlleles(); a++) {
                    final double v = mx.get(a, r);
                    matrixN++;
                    if (Double.doubleToRawLongBits((double) (float) v)
                            == Double.doubleToRawLongBits(v)) {
                        matrixF32++;
                    }
                    matrixBits.add(Double.doubleToRawLongBits(v));
                }
            }
        }
        int bufN = 0;
        int bufF32 = 0;
        final java.util.ArrayList<Long> bufBits = new java.util.ArrayList<>();
        if (buf != null) {
            bufN = buf.length;
            for (final double v : buf) {
                if (Double.doubleToRawLongBits((double) (float) v)
                        == Double.doubleToRawLongBits(v)) {
                    bufF32++;
                }
                bufBits.add(Double.doubleToRawLongBits(v));
            }
        }
        java.util.Collections.sort(matrixBits);
        java.util.Collections.sort(bufBits);
        final String prefix = USE_DOUBLE_PRECISION ? "6R98" : "6R97";
        kvLine(
                prefix,
                seq,
                "kernel_buffer_type",
                "double[] mLogLikelihoodArray");
        kvLine(prefix, seq, "kernel_buffer_n", Integer.toString(bufN));
        kvLine(prefix, seq, "kernel_buffer_f32_wide", Integer.toString(bufF32));
        kvLine(prefix, seq, "matrix_n", Integer.toString(matrixN));
        kvLine(prefix, seq, "matrix_f32_wide", Integer.toString(matrixF32));
        kvLine(
                prefix,
                seq,
                "buffer_matrix_sorted_bits_equal",
                Boolean.toString(matrixBits.equals(bufBits)));
        kvLine(prefix, seq, "likelihood_set", "SampleMatrix.set assignment");
        kvLine(
                prefix,
                seq,
                "use_double_precision",
                Boolean.toString(USE_DOUBLE_PRECISION));
        kvLine(prefix, seq, "pairhmm_class", hmm.getClass().getSimpleName());
        kvLine(
                prefix,
                seq,
                "gkl_double_path",
                "g_use_double ? g_compute_full_prob_double + log10 - LOG10_INITIAL_CONSTANT");
    }

    /**
     * Dump haplotype FNV + exact f64 bits keyed by (qname, flags). Used for the post-kernel
     * object (6R.96) before {@code normalizeLikelihoods}.
     */
    private static void dumpLikelihoodMatrix(
            final AlleleLikelihoods<GATKRead, Haplotype> ll, final String stage) {
        final int seq = callSeq + 1;
        kvDump(seq, "stage", stage);
        kvDump(seq, "n_alleles", Integer.toString(ll.numberOfAlleles()));
        kvDump(seq, "n_samples", Integer.toString(ll.numberOfSamples()));
        kvDump(seq, "evidence_count", Integer.toString(ll.evidenceCount()));
        kvDump(seq, "likelihood_object", "AlleleLikelihoods<GATKRead,Haplotype>");
        if (ll.numberOfSamples() > 0) {
            final LikelihoodMatrix<GATKRead, Haplotype> mx = ll.sampleMatrix(0);
            kvDump(seq, "column_count", Integer.toString(mx.numberOfAlleles()));
            for (int a = 0; a < mx.numberOfAlleles(); a++) {
                final Haplotype hap = mx.getAllele(a);
                final byte[] bases = hap.getBases();
                kvDump(
                        seq,
                        "hap_" + a,
                        "index="
                                + a
                                + "\tis_ref="
                                + hap.isReference()
                                + "\tlen="
                                + bases.length
                                + "\tfnv="
                                + Long.toHexString(fnv1a64(bases)));
            }
        }
        for (int s = 0; s < ll.numberOfSamples(); s++) {
            final List<GATKRead> ev = ll.sampleEvidence(s);
            kvDump(seq, "sample_evidence_n_" + s, Integer.toString(ev.size()));
            final LikelihoodMatrix<GATKRead, Haplotype> mx = ll.sampleMatrix(s);
            for (int r = 0; r < ev.size(); r++) {
                final GATKRead read = ev.get(r);
                kvDump(
                        seq,
                        "row_" + s + "_" + r,
                        "row="
                                + r
                                + "\tqname="
                                + read.getName()
                                + "\tflags="
                                + read.getFlags()
                                + "\tn="
                                + mx.numberOfAlleles());
                final StringBuilder bits = new StringBuilder();
                for (int a = 0; a < mx.numberOfAlleles(); a++) {
                    if (a > 0) {
                        bits.append(',');
                    }
                    bits.append(
                            Long.toHexString(Double.doubleToRawLongBits(mx.get(a, r))));
                }
                kvDump(
                        seq,
                        "rowbits_" + s + "_" + r,
                        "qname="
                                + read.getName()
                                + "\tflags="
                                + read.getFlags()
                                + "\tn="
                                + mx.numberOfAlleles()
                                + "\tbits="
                                + bits);
            }
        }
    }

    private static void dumpBeforeFilter(
            final AlleleLikelihoods<GATKRead, Haplotype> ll,
            final boolean dynamic,
            final double expectedErrorRatePerBase,
            final double scale) {
        final int seq = ++callSeq;
        kv(seq, "n_alleles", Integer.toString(ll.numberOfAlleles()));
        kv(seq, "n_samples", Integer.toString(ll.numberOfSamples()));
        kv(seq, "evidence_count", Integer.toString(ll.evidenceCount()));
        kv(seq, "dynamic_disqualification", Boolean.toString(dynamic));
        kv(seq, "expected_error_rate_per_base", Double.toString(expectedErrorRatePerBase));
        kv(seq, "read_disqualification_scale", Double.toString(scale));
        kv(seq, "cap_likelihoods", Boolean.toString(!dynamic));
        kv(seq, "likelihood_object", "AlleleLikelihoods<GATKRead,Haplotype>");
        kv(seq, "likelihood_matrix_class", "LikelihoodMatrix / sampleMatrix");
        kv(seq, "evidence_object", "sampleEvidence GATKRead");
        kv(seq, "column_object", "IndexedAlleleList haplotype columns");
        kv(seq, "max_ll_producer", "maximumLikelihoodOverAllAlleles");
        if (ll.numberOfSamples() > 0) {
            final LikelihoodMatrix<GATKRead, Haplotype> mx = ll.sampleMatrix(0);
            kv(seq, "column_count", Integer.toString(mx.numberOfAlleles()));
            for (int a = 0; a < mx.numberOfAlleles(); a++) {
                final Haplotype hap = mx.getAllele(a);
                final byte[] bases = hap.getBases();
                kv94(
                        seq,
                        "hap_" + a,
                        "index="
                                + a
                                + "\tis_ref="
                                + hap.isReference()
                                + "\tlen="
                                + bases.length
                                + "\tfnv="
                                + Long.toHexString(fnv1a64(bases)));
            }
        }
        for (int s = 0; s < ll.numberOfSamples(); s++) {
            final String sample = ll.getSample(s);
            final List<GATKRead> ev = ll.sampleEvidence(s);
            kv(seq, "sample_" + s, sample);
            kv(seq, "sample_evidence_n_" + s, Integer.toString(ev.size()));
            final LikelihoodMatrix<GATKRead, Haplotype> mx = ll.sampleMatrix(s);
            for (int r = 0; r < ev.size(); r++) {
                final GATKRead read = ev.get(r);
                final int readLen = read.getLength();
                final int bqLen = read.getBaseQualities().length;
                final Object hmmTag =
                        read.getTransientAttribute(
                                PairHMMLikelihoodCalculationEngine.HMM_BASE_QUALITIES_TAG);
                final boolean hmmPresent = hmmTag instanceof byte[];
                final int hmmBqLen = hmmPresent ? ((byte[]) hmmTag).length : -1;
                final int qualifiedLen = hmmPresent ? hmmBqLen : readLen;
                final double maxErrors =
                        Math.min(2.0, Math.ceil(qualifiedLen * expectedErrorRatePerBase));
                final double threshold = maxErrors * -4.0;
                final int argmax = argmaxAllele(mx, r);
                final double maxLl = argmax >= 0 ? mx.get(argmax, r) : Double.NEGATIVE_INFINITY;
                final boolean keep = !(maxLl < threshold);
                final Haplotype argHap = argmax >= 0 ? mx.getAllele(argmax) : null;
                final String argFnv =
                        argHap == null ? "none" : Long.toHexString(fnv1a64(argHap.getBases()));
                kv(
                        seq,
                        "read_" + s + "_" + r,
                        "sample="
                                + sample
                                + "\tqname="
                                + read.getName()
                                + "\tflags="
                                + read.getFlags()
                                + "\tread_len="
                                + readLen
                                + "\tbq_len="
                                + bqLen
                                + "\thmm_bq_present="
                                + hmmPresent
                                + "\thmm_bq_len="
                                + hmmBqLen
                                + "\tqualifiedLen="
                                + qualifiedLen
                                + "\tmax_errors="
                                + Double.toString(maxErrors)
                                + "\tthreshold="
                                + Double.toString(threshold)
                                + "\tmax_ll="
                                + Double.toString(maxLl)
                                + "\tkeep="
                                + keep
                                + "\tn_haps="
                                + ll.numberOfAlleles());
                kv94(
                        seq,
                        "row_" + s + "_" + r,
                        "sample="
                                + sample
                                + "\trow="
                                + r
                                + "\tqname="
                                + read.getName()
                                + "\tflags="
                                + read.getFlags()
                                + "\tstart="
                                + read.getStart()
                                + "\tend="
                                + read.getEnd()
                                + "\tcigar="
                                + read.getCigar().toString()
                                + "\tmax_ll="
                                + Double.toString(maxLl)
                                + "\targmax_col="
                                + argmax
                                + "\targmax_fnv="
                                + argFnv
                                + "\targmax_is_ref="
                                + (argHap != null && argHap.isReference())
                                + "\targmax_len="
                                + (argHap == null ? -1 : argHap.getBases().length)
                                + "\tcolumn_count="
                                + mx.numberOfAlleles()
                                + "\tkeep="
                                + keep);
                final StringBuilder bits = new StringBuilder();
                for (int a = 0; a < mx.numberOfAlleles(); a++) {
                    if (a > 0) {
                        bits.append(',');
                    }
                    bits.append(
                            Long.toHexString(
                                    Double.doubleToRawLongBits(mx.get(a, r))));
                }
                kv95(
                        seq,
                        "rowbits_" + s + "_" + r,
                        "qname="
                                + read.getName()
                                + "\tflags="
                                + read.getFlags()
                                + "\tn="
                                + mx.numberOfAlleles()
                                + "\tbits="
                                + bits);
            }
        }
    }

    /** FNV-1a 64 of haplotype bases; language-independent column identity. */
    static long fnv1a64(final byte[] data) {
        long h = 0xcbf29ce484222325L;
        for (int i = 0; i < data.length; i++) {
            h ^= (data[i] & 0xffL);
            h *= 0x100000001b3L;
        }
        return h;
    }

    /** Java {@code maximumLikelihoodOverAllAlleles}: first index with strict {@code >}. */
    private static int argmaxAllele(
            final LikelihoodMatrix<GATKRead, Haplotype> mx, final int evidenceIndex) {
        double result = Double.NEGATIVE_INFINITY;
        int arg = -1;
        for (int a = 0; a < mx.numberOfAlleles(); a++) {
            final double v = mx.get(a, evidenceIndex);
            if (v > result) {
                result = v;
                arg = a;
            }
        }
        return arg;
    }

    private static double fieldDouble(final Object obj, final String name) throws Exception {
        final Field f = obj.getClass().getDeclaredField(name);
        f.setAccessible(true);
        return f.getDouble(obj);
    }

    private static boolean fieldBoolean(final Object obj, final String name) throws Exception {
        final Field f = obj.getClass().getDeclaredField(name);
        f.setAccessible(true);
        return f.getBoolean(obj);
    }

    static void kv(final int seq, final String key, final String value) {
        System.out.println("6R93\t" + seq + "\t" + key + "\t" + value);
    }

    static void kv94(final int seq, final String key, final String value) {
        System.out.println("6R94\t" + seq + "\t" + key + "\t" + value);
    }

    static void kv95(final int seq, final String key, final String value) {
        System.out.println("6R95\t" + seq + "\t" + key + "\t" + value);
    }

    static void kv96(final int seq, final String key, final String value) {
        System.out.println("6R96\t" + seq + "\t" + key + "\t" + value);
    }

    static void kvDump(final int seq, final String key, final String value) {
        System.out.println(DUMP_PREFIX + "\t" + seq + "\t" + key + "\t" + value);
    }

    static void kvLine(
            final String prefix, final int seq, final String key, final String value) {
        System.out.println(prefix + "\t" + seq + "\t" + key + "\t" + value);
    }

    static void kv97(final int seq, final String key, final String value) {
        System.out.println("6R97\t" + seq + "\t" + key + "\t" + value);
    }
}
