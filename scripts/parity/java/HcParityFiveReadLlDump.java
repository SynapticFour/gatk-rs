import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMUtils;
import java.lang.reflect.Field;
import java.lang.reflect.Method;
import java.util.List;
import java.util.Locale;
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
import org.broadinstitute.hellbender.utils.pairhmm.PairHMMInputScoreImputation;
import org.broadinstitute.hellbender.utils.pairhmm.PairHMMInputScoreImputator;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

/**
 * TEST-ONLY 6R.230: dump PairHMM inputs + prim/norm cells for the frozen five
 * reads. Delegates production arithmetic. Does not change keep/drop.
 */
public final class HcParityFiveReadLlDump implements ReadLikelihoodCalculationEngine {

    static final String PFX = "6R230";
    static final String[] QNAMES = HcParityHapLlMembershipDump.QNAMES;
    static final int[] FLAGS = HcParityHapLlMembershipDump.FLAGS;

    private final ReadLikelihoodCalculationEngine inner;

    public static void installOn(final HaplotypeCallerEngine engine) throws Exception {
        final Field f =
                HaplotypeCallerEngine.class.getDeclaredField("likelihoodCalculationEngine");
        f.setAccessible(true);
        final ReadLikelihoodCalculationEngine inner =
                (ReadLikelihoodCalculationEngine) f.get(engine);
        f.set(engine, new HcParityFiveReadLlDump(inner));
    }

    private HcParityFiveReadLlDump(final ReadLikelihoodCalculationEngine inner) {
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
            throw new RuntimeException("6R.230 five-read ll dump failed", e);
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
        kv("hap_count", Integer.toString(haplotypeList.size()));
        int origN = 0;
        for (final List<GATKRead> reads : perSampleReadList.values()) {
            origN += reads.size();
        }
        kv("orig_evidence_count", Integer.toString(origN));
        for (int i = 0; i < haplotypeList.size(); i++) {
            dumpHap(i, haplotypeList.get(i));
        }
        for (final List<GATKRead> reads : perSampleReadList.values()) {
            for (final GATKRead r : reads) {
                if (isFive(r)) {
                    dumpRead("orig", r, null);
                }
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
        for (int s = 0; s < result.numberOfSamples(); s++) {
            modify.invoke(phmm, result.sampleEvidence(s));
            @SuppressWarnings("unchecked")
            final List<GATKRead> evidence = (List<GATKRead>) result.sampleEvidence(s);
            for (final GATKRead r : evidence) {
                if (isFive(r)) {
                    dumpRead("proc", r, imputator);
                }
            }
            computeOne.invoke(phmm, result.sampleMatrix(s));
        }
        dumpFiveMatrix(result, "prim");
        final double log10global = fieldDouble(phmm, "log10globalReadMismappingRate");
        final boolean sym = fieldBoolean(phmm, "symmetricallyNormalizeAllelesToReference");
        result.normalizeLikelihoods(log10global, sym);
        dumpFiveMatrix(result, "norm");
        final boolean dynamic = fieldBoolean(phmm, "dynamicDisqualification");
        final double expectedError = fieldDouble(phmm, "expectedErrorRatePerBase");
        final double scale = fieldDouble(phmm, "readDisqualificationScale");
        phmm.filterPoorlyModeledEvidence(result, dynamic, expectedError, scale);
        kv("stored_evidence_count", Integer.toString(result.evidenceCount()));
        return result;
    }

    private static boolean isFive(final GATKRead r) {
        for (int i = 0; i < QNAMES.length; i++) {
            if (QNAMES[i].equals(r.getName()) && r.getFlags() == FLAGS[i]) {
                return true;
            }
        }
        return false;
    }

    private static void dumpHap(final int idx, final Haplotype h) {
        final byte[] bases = h.getBases();
        final String loc = h.getGenomeLocation() == null ? "." : h.getGenomeLocation().toString();
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
                        + cigar);
    }

    private static void dumpRead(
            final String kind, final GATKRead r, final PairHMMInputScoreImputator imputator) {
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
        kv(
                kind,
                "qname="
                        + r.getName()
                        + "\tflags="
                        + r.getFlags()
                        + "\tstart="
                        + r.getStart()
                        + "\tend="
                        + r.getEnd()
                        + "\tmapq="
                        + r.getMappingQuality()
                        + "\tlen="
                        + r.getLength()
                        + "\tcigar="
                        + (r.getCigar() == null ? "." : r.getCigar().toString())
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
                        + (gcp.length == 0 ? "." : SAMUtils.phredToFastq(gcp)));
    }

    private static void dumpFiveMatrix(
            final AlleleLikelihoods<GATKRead, Haplotype> ll, final String stage) {
        if (ll.numberOfSamples() < 1) {
            return;
        }
        final LikelihoodMatrix<GATKRead, Haplotype> mx = ll.sampleMatrix(0);
        kv(stage + "_n_hap", Integer.toString(mx.numberOfAlleles()));
        kv(stage + "_n_ev", Integer.toString(mx.evidenceCount()));
        final String[] hashes = new String[mx.numberOfAlleles()];
        final boolean[] isRef = new boolean[mx.numberOfAlleles()];
        for (int a = 0; a < mx.numberOfAlleles(); a++) {
            hashes[a] = fnv1a64Hex(mx.getAllele(a).getBases());
            isRef[a] = mx.getAllele(a).isReference();
        }
        for (int r = 0; r < mx.evidenceCount(); r++) {
            final GATKRead ev = mx.getEvidence(r);
            if (!isFive(ev)) {
                continue;
            }
            double maxLl = Double.NEGATIVE_INFINITY;
            int win = -1;
            for (int a = 0; a < mx.numberOfAlleles(); a++) {
                final double v = mx.get(a, r);
                if (v > maxLl) {
                    maxLl = v;
                    win = a;
                }
                kv(
                        stage,
                        ev.getName()
                                + "\tflags="
                                + ev.getFlags()
                                + "\thap="
                                + hashes[a]
                                + "\tisRef="
                                + isRef[a]
                                + "\tll="
                                + String.format(Locale.US, "%.12f", v)
                                + "\tf32wide="
                                + Boolean.toString(
                                        Double.doubleToRawLongBits((double) (float) v)
                                                == Double.doubleToRawLongBits(v)));
            }
            kv(
                    stage + "_max",
                    ev.getName()
                            + "\tflags="
                            + ev.getFlags()
                            + "\tmax_ll="
                            + String.format(Locale.US, "%.12f", maxLl)
                            + "\twin="
                            + (win < 0 ? "." : hashes[win])
                            + "\twinIsRef="
                            + (win >= 0 && isRef[win]));
        }
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
        System.out.println(PFX + "\t" + key + "\t" + value);
    }
}
