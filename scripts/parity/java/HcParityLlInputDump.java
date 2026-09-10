import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMUtils;
import htsjdk.variant.variantcontext.VariantContext;
import java.lang.reflect.Field;
import java.lang.reflect.Method;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.TreeSet;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AssemblyBasedCallerUtils;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCallerEngine;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.PairHMMLikelihoodCalculationEngine;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.ReadLikelihoodCalculationEngine;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.genotyper.IndexedAlleleList;
import org.broadinstitute.hellbender.utils.genotyper.LikelihoodMatrix;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;
import org.broadinstitute.hellbender.utils.haplotype.EventMap;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.pairhmm.PairHMM;
import org.broadinstitute.hellbender.utils.pairhmm.PairHMMInputScoreImputation;
import org.broadinstitute.hellbender.utils.pairhmm.PairHMMInputScoreImputator;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

/**
 * TEST-ONLY 6R.124: dump likelihood-engine input objects and primitive vs stored
 * haplotype likelihoods. Delegates production arithmetic (initialize, PairHMM,
 * normalize, filter). Does not inspect PairHMM DP.
 */
public final class HcParityLlInputDump implements ReadLikelihoodCalculationEngine {

    private final ReadLikelihoodCalculationEngine inner;
    private static String dumpPrefix = "6R124";
    private static int vcLoc = -1;

    public static void installOn(final HaplotypeCallerEngine engine) throws Exception {
        installOn(engine, "6R124", -1);
    }

    public static void installOn(
            final HaplotypeCallerEngine engine, final String prefix, final int loc)
            throws Exception {
        dumpPrefix = prefix;
        vcLoc = loc;
        final Field f =
                HaplotypeCallerEngine.class.getDeclaredField("likelihoodCalculationEngine");
        f.setAccessible(true);
        final ReadLikelihoodCalculationEngine inner =
                (ReadLikelihoodCalculationEngine) f.get(engine);
        f.set(engine, new HcParityLlInputDump(inner));
    }

    private HcParityLlInputDump(final ReadLikelihoodCalculationEngine inner) {
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
            return computeAndDump(
                    (PairHMMLikelihoodCalculationEngine) inner,
                    haplotypeList,
                    samples,
                    perSampleReadList);
        } catch (final Exception e) {
            throw new RuntimeException("6R.124 likelihood-input dump failed", e);
        }
    }

    private static AlleleLikelihoods<GATKRead, Haplotype> computeAndDump(
            final PairHMMLikelihoodCalculationEngine phmm,
            final List<Haplotype> haplotypeList,
            final SampleList samples,
            final Map<String, List<GATKRead>> perSampleReadList)
            throws Exception {
        kv("pairhmm_class", pairHmmClass(phmm));
        kv("constant_gcp", Byte.toString(fieldByte(phmm, "constantGCP")));
        kv("pcr_error_model", String.valueOf(fieldObject(phmm, "pcrErrorModel")));
        kv("modify_softclipped_bases", Boolean.toString(fieldBoolean(phmm, "modifySoftclippedBases")));
        kv("disable_cap_mapq", Boolean.toString(fieldBoolean(phmm, "disableCapReadQualitiesToMapQ")));
        kv("bq_threshold", Byte.toString(fieldByte(phmm, "baseQualityScoreThreshold")));
        kv("hap_count", Integer.toString(haplotypeList.size()));
        int origN = 0;
        for (final List<GATKRead> reads : perSampleReadList.values()) {
            origN += reads.size();
        }
        kv("orig_evidence_count", Integer.toString(origN));

        for (int i = 0; i < haplotypeList.size(); i++) {
            dumpHap(i, haplotypeList.get(i));
        }
        dumpEventMapsAndVcs(haplotypeList);
        int oi = 0;
        for (final List<GATKRead> reads : perSampleReadList.values()) {
            for (final GATKRead r : reads) {
                dumpRead("orig", oi++, r, null);
            }
        }

        final Method init =
                PairHMMLikelihoodCalculationEngine.class.getDeclaredMethod(
                        "initializePairHMM", List.class, Map.class);
        init.setAccessible(true);
        init.invoke(phmm, haplotypeList, perSampleReadList);

        final AlleleLikelihoods<GATKRead, Haplotype> result =
                new AlleleLikelihoods<>(
                        samples, new IndexedAlleleList<>(haplotypeList), perSampleReadList);
        final Method modify =
                PairHMMLikelihoodCalculationEngine.class.getDeclaredMethod(
                        "modifyReadQualities", List.class);
        modify.setAccessible(true);
        final Method computeOne =
                PairHMMLikelihoodCalculationEngine.class.getDeclaredMethod(
                        "computeReadLikelihoods", LikelihoodMatrix.class);
        computeOne.setAccessible(true);
        final PairHMMInputScoreImputator imputator =
                (PairHMMInputScoreImputator) fieldObject(phmm, "inputScoreImputator");

        int procN = 0;
        for (int s = 0; s < result.numberOfSamples(); s++) {
            @SuppressWarnings("unchecked")
            final List<GATKRead> processed =
                    (List<GATKRead>) modify.invoke(phmm, result.sampleEvidence(s));
            for (final GATKRead r : processed) {
                dumpRead("proc", procN++, r, imputator);
            }
            computeOne.invoke(phmm, result.sampleMatrix(s));
        }
        kv("proc_evidence_count", Integer.toString(procN));

        dumpMatrix(result, "prim");
        dumpTransfer(phmm, result);

        final double log10global = fieldDouble(phmm, "log10globalReadMismappingRate");
        final boolean sym = fieldBoolean(phmm, "symmetricallyNormalizeAllelesToReference");
        result.normalizeLikelihoods(log10global, sym);
        dumpMatrix(result, "norm");

        final boolean dynamic = fieldBoolean(phmm, "dynamicDisqualification");
        final double expectedError = fieldDouble(phmm, "expectedErrorRatePerBase");
        final double scale = fieldDouble(phmm, "readDisqualificationScale");
        phmm.filterPoorlyModeledEvidence(result, dynamic, expectedError, scale);
        dumpMatrix(result, "stored");
        kv("stored_evidence_count", Integer.toString(result.evidenceCount()));
        return result;
    }

    private static void dumpEventMapsAndVcs(final List<Haplotype> haplotypeList) {
        final TreeSet<String> union = new TreeSet<>();
        int nWithMap = 0;
        int nRefEvents = 0;
        for (int i = 0; i < haplotypeList.size(); i++) {
            final Haplotype h = haplotypeList.get(i);
            final EventMap em = h.getEventMap();
            if (em == null) {
                kv("em", "idx=" + i + "\thash=" + fnv1a64Hex(h.getBases()) + "\tevent_n=0\tmap=null");
                continue;
            }
            nWithMap++;
            final ArrayList<VariantContext> evs = new ArrayList<>();
            for (final Object o : em.values()) {
                evs.add((VariantContext) o);
            }
            kv(
                    "em",
                    "idx="
                            + i
                            + "\thash="
                            + fnv1a64Hex(h.getBases())
                            + "\tisRef="
                            + h.isReference()
                            + "\tevent_n="
                            + evs.size());
            for (final VariantContext vc : evs) {
                final String alt =
                        vc.getAlternateAlleles().isEmpty()
                                ? "."
                                : vc.getAlternateAlleles().get(0).getBaseString();
                final String key =
                        vc.getStart() + ":" + vc.getReference().getBaseString() + ">" + alt;
                union.add(key);
                if (h.isReference()) {
                    nRefEvents++;
                }
                kv(
                        "em_event",
                        "idx="
                                + i
                                + "\thash="
                                + fnv1a64Hex(h.getBases())
                                + "\tstart="
                                + vc.getStart()
                                + "\tend="
                                + vc.getEnd()
                                + "\tref="
                                + vc.getReference().getBaseString()
                                + "\talt="
                                + alt);
            }
        }
        kv("em_haps_with_map", Integer.toString(nWithMap));
        kv("em_union_n", Integer.toString(union.size()));
        kv("em_ref_hap_event_rows", Integer.toString(nRefEvents));
        for (final String key : union) {
            kv("em_union", key);
        }
        if (vcLoc <= 0) {
            return;
        }
        try {
            final List events =
                    AssemblyBasedCallerUtils.getVariantContextsFromActiveHaplotypes(
                            vcLoc, haplotypeList, true);
            kv("vc_loc", Integer.toString(vcLoc));
            kv("vc_n", Integer.toString(events.size()));
            for (int i = 0; i < events.size(); i++) {
                final VariantContext vc = (VariantContext) events.get(i);
                kv(
                        "vc",
                        "i="
                                + i
                                + "\tstart="
                                + vc.getStart()
                                + "\tend="
                                + vc.getEnd()
                                + "\tref="
                                + vc.getReference().getBaseString()
                                + "\talts="
                                + altList(vc)
                                + "\tnAlleles="
                                + vc.getNAlleles());
            }
            final List eventsNoSpan =
                    AssemblyBasedCallerUtils.getVariantContextsFromActiveHaplotypes(
                            vcLoc, haplotypeList, false);
            kv("vc_n_no_span", Integer.toString(eventsNoSpan.size()));
        } catch (final Exception e) {
            kv("vc_error", e.getClass().getSimpleName() + "\t" + String.valueOf(e.getMessage()));
        }
    }

    private static String altList(final VariantContext vc) {
        final StringBuilder sb = new StringBuilder();
        for (int i = 0; i < vc.getAlternateAlleles().size(); i++) {
            if (i > 0) {
                sb.append(",");
            }
            sb.append(vc.getAlternateAlleles().get(i).getBaseString());
        }
        return sb.length() == 0 ? "." : sb.toString();
    }

    private static void dumpHap(final int idx, final Haplotype h) {
        final byte[] bases = h.getBases();
        final String loc;
        if (h.getGenomeLocation() == null) {
            loc = ".";
        } else {
            loc = h.getGenomeLocation().toString();
        }
        final String cigar = h.getCigar() == null ? "." : h.getCigar().toString();
        kv(
                "hap",
                idx
                        + "\t"
                        + fnv1a64Hex(bases)
                        + "\tlen="
                        + bases.length
                        + "\tisRef="
                        + h.isReference()
                        + "\tloc="
                        + loc
                        + "\talignStart="
                        + h.getAlignmentStartHapwrtRef()
                        + "\tcigar="
                        + cigar
                        + "\tbases="
                        + new String(bases));
    }

    private static void dumpRead(
            final String kind,
            final int idx,
            final GATKRead r,
            final PairHMMInputScoreImputator imputator) {
        final byte[] bases = r.getBasesNoCopy();
        final byte[] bq = r.getBaseQualitiesNoCopy();
        byte[] iq;
        byte[] dq;
        byte[] gcp = new byte[0];
        if (imputator != null) {
            final PairHMMInputScoreImputation imp = imputator.impute(r);
            iq = imp.insOpenPenalties();
            dq = imp.delOpenPenalties();
            gcp = imp.gapContinuationPenalties();
        } else {
            iq = ReadUtils.getBaseInsertionQualities(r);
            dq = ReadUtils.getBaseDeletionQualities(r);
        }
        final Object hmm =
                r.getTransientAttribute(PairHMMLikelihoodCalculationEngine.HMM_BASE_QUALITIES_TAG);
        final String hmmFq =
                hmm instanceof byte[] ? SAMUtils.phredToFastq((byte[]) hmm) : ".";
        final String cigar = r.getCigar() == null ? "." : r.getCigar().toString();
        kv(
                kind,
                idx
                        + "\tqname="
                        + r.getName()
                        + "\tflags="
                        + r.getFlags()
                        + "\tstart="
                        + r.getStart()
                        + "\tend="
                        + r.getEnd()
                        + "\tuStart="
                        + r.getUnclippedStart()
                        + "\tuEnd="
                        + r.getUnclippedEnd()
                        + "\tmapq="
                        + r.getMappingQuality()
                        + "\tlen="
                        + r.getLength()
                        + "\tcigar="
                        + cigar
                        + "\tbasesHash="
                        + fnv1a64Hex(bases)
                        + "\tbqHash="
                        + fnv1a64Hex(bq)
                        + "\tiqHash="
                        + fnv1a64Hex(iq)
                        + "\tdqHash="
                        + fnv1a64Hex(dq)
                        + "\tgcpHash="
                        + fnv1a64Hex(gcp)
                        + "\thmmHash="
                        + (hmm instanceof byte[] ? fnv1a64Hex((byte[]) hmm) : ".")
                        + "\tbases="
                        + new String(bases)
                        + "\tbq="
                        + SAMUtils.phredToFastq(bq)
                        + "\tiq="
                        + SAMUtils.phredToFastq(iq)
                        + "\tdq="
                        + SAMUtils.phredToFastq(dq)
                        + "\tgcp="
                        + (gcp.length == 0 ? "." : SAMUtils.phredToFastq(gcp))
                        + "\thmmBq="
                        + hmmFq);
    }

    private static void dumpMatrix(
            final AlleleLikelihoods<GATKRead, Haplotype> ll, final String stage) {
        if (ll.numberOfSamples() < 1) {
            return;
        }
        final LikelihoodMatrix<GATKRead, Haplotype> mx = ll.sampleMatrix(0);
        kv(stage + "_n_hap", Integer.toString(mx.numberOfAlleles()));
        kv(stage + "_n_ev", Integer.toString(mx.evidenceCount()));
        final String[] hashes = new String[mx.numberOfAlleles()];
        for (int a = 0; a < mx.numberOfAlleles(); a++) {
            hashes[a] = fnv1a64Hex(mx.getAllele(a).getBases());
        }
        for (int r = 0; r < mx.evidenceCount(); r++) {
            final GATKRead ev = mx.getEvidence(r);
            for (int a = 0; a < mx.numberOfAlleles(); a++) {
                final double v = mx.get(a, r);
                kv(
                        stage,
                        ev.getName()
                                + "\tflags="
                                + ev.getFlags()
                                + "\tstart="
                                + ev.getStart()
                                + "\thap="
                                + hashes[a]
                                + "\tll="
                                + String.format(Locale.US, "%.12f", v)
                                + "\tbits="
                                + Long.toHexString(Double.doubleToRawLongBits(v))
                                + "\tf32wide="
                                + Boolean.toString(
                                        Double.doubleToRawLongBits((double) (float) v)
                                                == Double.doubleToRawLongBits(v))
                                + "\tfinite="
                                + Boolean.toString(Double.isFinite(v)));
            }
        }
    }

    private static void dumpTransfer(
            final PairHMMLikelihoodCalculationEngine phmm,
            final AlleleLikelihoods<GATKRead, Haplotype> ll)
            throws Exception {
        final Field pf = PairHMMLikelihoodCalculationEngine.class.getDeclaredField("pairHMM");
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
        kv("kernel_buffer_type", "double[] mLogLikelihoodArray");
        kv("kernel_buffer_n", Integer.toString(bufN));
        kv("kernel_buffer_f32_wide", Integer.toString(bufF32));
        kv("matrix_n", Integer.toString(matrixN));
        kv("matrix_f32_wide", Integer.toString(matrixF32));
        kv("buffer_matrix_sorted_bits_equal", Boolean.toString(matrixBits.equals(bufBits)));
        kv("likelihood_set", "SampleMatrix.set assignment");
    }

    private static String pairHmmClass(final PairHMMLikelihoodCalculationEngine phmm)
            throws Exception {
        final Field pf = PairHMMLikelihoodCalculationEngine.class.getDeclaredField("pairHMM");
        pf.setAccessible(true);
        final PairHMM hmm = (PairHMM) pf.get(phmm);
        return hmm.getClass().getSimpleName();
    }

    private static Object fieldObject(final Object obj, final String name) throws Exception {
        final Field f = obj.getClass().getDeclaredField(name);
        f.setAccessible(true);
        return f.get(obj);
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

    private static byte fieldByte(final Object obj, final String name) throws Exception {
        final Field f = obj.getClass().getDeclaredField(name);
        f.setAccessible(true);
        return f.getByte(obj);
    }

    private static String fnv1a64Hex(final byte[] data) {
        if (data == null) {
            return ".";
        }
        long h = 0xcbf29ce484222325L;
        for (int i = 0; i < data.length; i++) {
            h ^= (data[i] & 0xffL);
            h *= 0x100000001b3L;
        }
        return String.format("%016x", h);
    }

    private static void kv(final String key, final String value) {
        System.out.println(dumpPrefix + "\t" + key + "\t" + value);
    }
}
