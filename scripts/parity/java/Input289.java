import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import java.lang.reflect.Field;
import java.lang.reflect.Method;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import org.broadinstitute.hellbender.engine.AssemblyRegion;
import org.broadinstitute.hellbender.engine.AssemblyRegionIterator;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.MultiIntervalLocalReadShard;
import org.broadinstitute.hellbender.engine.ReadsPathDataSource;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.spark.AssemblyRegionArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.annotator.VariantAnnotatorEngine;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCallerArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCallerEngine;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.PairHMMLikelihoodCalculationEngine;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.ReadLikelihoodCalculationEngine;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.pairhmm.PairHMM;
import org.broadinstitute.hellbender.utils.pairhmm.PairHMMInputScoreImputation;
import org.broadinstitute.hellbender.utils.pairhmm.PairHMMInputScoreImputator;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.PairHMMNativeArgumentCollection;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

/**
 * 6R.289: PairHMM inputs for one read and haplotype 0, after
 * {@code modifyReadQualities} and before the PairHMM recurrence.
 */
public final class Input289 implements ReadLikelihoodCalculationEngine {
    private static final String QNAME = "HISEQ1:11:H8GV6ADXX:1:2116:18670:99941";
    private static final int FLAGS = 99;
    private static final int LOC = 29455649;
    private static final int PADDING = 100;

    private final ReadLikelihoodCalculationEngine inner;
    private static int calls = 0;

    public static void main(final String[] args) throws Exception {
        if (args.length != 2) {
            throw new IllegalArgumentException("usage: Input289 ref bam");
        }
        Utils.resetRandomGenerator();
        final Path bam = Paths.get(args[1]);
        final ReadsPathDataSource readsSource = new ReadsPathDataSource(bam);
        final SAMFileHeader header = readsSource.getHeader();
        final Path ref = Paths.get(args[0]);
        final CachingIndexedFastaSequenceFile refReader = new CachingIndexedFastaSequenceFile(ref);
        final ReferenceDataSource reference = ReferenceDataSource.of(ref);
        final AssemblyRegionArgumentCollection asmArgs = new AssemblyRegionArgumentCollection();
        asmArgs.assemblyRegionPadding = PADDING;
        final HaplotypeCallerArgumentCollection hcArgs = new HaplotypeCallerArgumentCollection();
        hcArgs.likelihoodArgs.pairHMM = PairHMM.Implementation.LOGLESS_CACHING;
        final Field threads =
                PairHMMNativeArgumentCollection.class.getDeclaredField("pairHmmNativeThreads");
        threads.setAccessible(true);
        threads.setInt(hcArgs.likelihoodArgs.pairHMMNativeArgs, 1);
        final Field useDouble =
                PairHMMNativeArgumentCollection.class.getDeclaredField("useDoublePrecision");
        useDouble.setAccessible(true);
        useDouble.setBoolean(hcArgs.likelihoodArgs.pairHMMNativeArgs, false);
        final VariantAnnotatorEngine annotationEngine =
                new VariantAnnotatorEngine(
                        Collections.emptyList(), null, Collections.emptyList(), false, false);
        final HaplotypeCallerEngine engine =
                new HaplotypeCallerEngine(
                        hcArgs, asmArgs, false, false, header, refReader, annotationEngine);
        install(engine);
        final List<ReadFilter> readFilters = HaplotypeCallerEngine.makeStandardHCReadFilters();
        for (final ReadFilter f : readFilters) {
            f.setHeader(header);
        }
        final SAMSequenceDictionary dict = header.getSequenceDictionary();
        final SimpleInterval interval = new SimpleInterval("20", 29455000, 29456500);
        final MultiIntervalLocalReadShard shard =
                new MultiIntervalLocalReadShard(
                        Collections.singletonList(interval), PADDING, readsSource);
        shard.setPreReadFilterTransformer(HaplotypeCallerEngine.makeStandardHCReadTransformer());
        shard.setReadFilter(ReadFilter.fromList(readFilters, header));
        shard.setPostReadFilterTransformer(
                org.broadinstitute.hellbender.transformers.ReadTransformer.identity());
        final AssemblyRegionIterator iter =
                new AssemblyRegionIterator(shard, header, reference, null, engine, asmArgs, false);
        boolean found = false;
        while (iter.hasNext()) {
            final AssemblyRegion region = iter.next();
            if (!region.isActive() || LOC < region.getStart() || LOC > region.getEnd()) {
                continue;
            }
            kv("region", region.getContig() + ":" + region.getStart() + "-" + region.getEnd());
            final ReferenceContext refCtx =
                    new ReferenceContext(reference, region.getSpan(), PADDING, PADDING);
            engine.callRegion(region, new FeatureContext(), refCtx);
            found = true;
            break;
        }
        engine.shutdown();
        readsSource.close();
        if (!found) {
            throw new IllegalStateException("no active region covering " + LOC);
        }
        if (calls == 0) {
            throw new IllegalStateException("target read was not passed to PairHMM");
        }
    }

    private static void install(final HaplotypeCallerEngine engine) throws Exception {
        final Field f = HaplotypeCallerEngine.class.getDeclaredField("likelihoodCalculationEngine");
        f.setAccessible(true);
        final ReadLikelihoodCalculationEngine inner =
                (ReadLikelihoodCalculationEngine) f.get(engine);
        f.set(engine, new Input289(inner));
    }

    private Input289(final ReadLikelihoodCalculationEngine inner) {
        this.inner = inner;
    }

    @Override
    public void close() {
        inner.close();
    }

    @Override
    public AlleleLikelihoods<GATKRead, Haplotype> computeReadLikelihoods(
            final List<Haplotype> haplotypeList,
            final SAMFileHeader hdr,
            final SampleList samples,
            final Map<String, List<GATKRead>> perSampleReadList,
            final boolean filterPoorly) {
        if (inner instanceof PairHMMLikelihoodCalculationEngine) {
            try {
                dumpInputs((PairHMMLikelihoodCalculationEngine) inner, haplotypeList, perSampleReadList);
            } catch (final Exception e) {
                throw new RuntimeException("6R.289 input dump failed", e);
            }
        }
        return inner.computeReadLikelihoods(
                haplotypeList, hdr, samples, perSampleReadList, filterPoorly);
    }

    private static void dumpInputs(
            final PairHMMLikelihoodCalculationEngine phmm,
            final List<Haplotype> haplotypeList,
            final Map<String, List<GATKRead>> perSampleReadList)
            throws Exception {
        GATKRead orig = null;
        for (final List<GATKRead> reads : perSampleReadList.values()) {
            for (final GATKRead r : reads) {
                if (QNAME.equals(r.getName()) && r.getFlags() == FLAGS) {
                    orig = r;
                }
            }
        }
        if (orig == null) {
            return;
        }
        final int call = ++calls;
        if (call == 1) {
            kv("constant_gcp", Byte.toString(fieldByte(phmm, "constantGCP")));
            kv("pcr_error_model", String.valueOf(fieldObject(phmm, "pcrErrorModel")));
            kv("bq_threshold", Byte.toString(fieldByte(phmm, "baseQualityScoreThreshold")));
            kv(
                    "disable_cap_mapq",
                    Boolean.toString(fieldBoolean(phmm, "disableCapReadQualitiesToMapQ")));
            kv(
                    "modify_softclipped_bases",
                    Boolean.toString(fieldBoolean(phmm, "modifySoftclippedBases")));
        }
        final Haplotype hap0 = haplotypeList.get(0);
        kv("call", Integer.toString(call));
        kv("n_haps", Integer.toString(haplotypeList.size()));
        kv("hap0_len", Integer.toString(hap0.getBases().length));
        kv("hap0_is_ref", Boolean.toString(hap0.isReference()));
        kv("hap0_align_start", Integer.toString(hap0.getAlignmentStartHapwrtRef()));
        kv("hap0_bases", new String(hap0.getBases()));
        dumpRead("orig", call, orig, null);

        final Method init =
                PairHMMLikelihoodCalculationEngine.class.getDeclaredMethod(
                        "initializePairHMM", List.class, Map.class);
        init.setAccessible(true);
        init.invoke(phmm, haplotypeList, perSampleReadList);
        final Method modify =
                PairHMMLikelihoodCalculationEngine.class.getDeclaredMethod(
                        "modifyReadQualities", List.class);
        modify.setAccessible(true);
        final PairHMMInputScoreImputator imputator =
                (PairHMMInputScoreImputator) fieldObject(phmm, "inputScoreImputator");
        boolean foundProc = false;
        for (final List<GATKRead> reads : perSampleReadList.values()) {
            int target = -1;
            for (int i = 0; i < reads.size(); i++) {
                final GATKRead r = reads.get(i);
                if (QNAME.equals(r.getName()) && r.getFlags() == FLAGS) {
                    target = i;
                }
            }
            if (target < 0) {
                continue;
            }
            @SuppressWarnings("unchecked")
            final List<GATKRead> processed = (List<GATKRead>) modify.invoke(phmm, reads);
            kv("imputator_null", Boolean.toString(imputator == null));
            kv("proc_name", processed.get(target).getName());
            dumpRead("proc", call, processed.get(target), imputator);
            foundProc = true;
        }
        if (foundProc) {
            System.out.flush();
            System.exit(0);
        }
    }

    private static void dumpRead(
            final String kind,
            final int call,
            final GATKRead r,
            final PairHMMInputScoreImputator imputator) {
        final byte[] bases = r.getBases();
        final byte[] bq = r.getBaseQualities();
        byte[] iq;
        byte[] dq;
        byte[] gcp;
        if (imputator != null) {
            final PairHMMInputScoreImputation imp = imputator.impute(r);
            iq = imp.insOpenPenalties();
            dq = imp.delOpenPenalties();
            gcp = imp.gapContinuationPenalties();
        } else {
            iq = ReadUtils.getBaseInsertionQualities(r);
            dq = ReadUtils.getBaseDeletionQualities(r);
            gcp = new byte[0];
        }
        final String cigar = r.getCigar() == null ? "." : r.getCigar().toString();
        kv(
                kind,
                "call="
                        + call
                        + "\tqname="
                        + r.getName()
                        + "\tflags="
                        + r.getFlags()
                        + "\tmapq="
                        + r.getMappingQuality()
                        + "\tstart="
                        + r.getStart()
                        + "\tend="
                        + r.getEnd()
                        + "\tuStart="
                        + r.getUnclippedStart()
                        + "\tuEnd="
                        + r.getUnclippedEnd()
                        + "\tlen="
                        + bases.length
                        + "\treverse="
                        + r.isReverseStrand()
                        + "\tcigar="
                        + cigar
                        + "\tbases="
                        + new String(bases)
                        + "\tbq="
                        + csv(bq)
                        + "\tiq="
                        + csv(iq)
                        + "\tdq="
                        + csv(dq)
                        + "\tgcp="
                        + (gcp.length == 0 ? "." : csv(gcp)));
    }

    private static String csv(final byte[] q) {
        final StringBuilder b = new StringBuilder();
        for (int i = 0; i < q.length; i++) {
            if (i > 0) {
                b.append(',');
            }
            b.append(q[i] & 0xff);
        }
        return b.toString();
    }

    private static void kv(final String key, final String value) {
        System.out.println("6R289\t" + key + "\t" + value);
    }

    private static Object fieldObject(final Object o, final String name) throws Exception {
        final Field f = o.getClass().getDeclaredField(name);
        f.setAccessible(true);
        return f.get(o);
    }

    private static byte fieldByte(final Object o, final String name) throws Exception {
        return ((Byte) fieldObject(o, name)).byteValue();
    }

    private static boolean fieldBoolean(final Object o, final String name) throws Exception {
        return ((Boolean) fieldObject(o, name)).booleanValue();
    }
}
