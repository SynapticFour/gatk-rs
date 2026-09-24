import htsjdk.samtools.SAMFileHeader;
import java.lang.reflect.Field;
import java.lang.reflect.Method;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.engine.AssemblyRegion;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AssemblyBasedCallerUtils;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AssemblyRegionTrimmer;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AssemblyResultSet;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCallerArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCallerEngine;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.PairHMMLikelihoodCalculationEngine;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.ReadLikelihoodCalculationEngine;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading.ReadThreadingAssembler;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.clipping.ReadClipper;
import org.broadinstitute.hellbender.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.genotyper.IndexedAlleleList;
import org.broadinstitute.hellbender.utils.genotyper.LikelihoodMatrix;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.read.AlignmentUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanAligner;
import htsjdk.variant.variantcontext.VariantContext;
import java.util.Collections;
import java.util.Comparator;
import java.util.TreeSet;

/**
 * TEST-ONLY 6R.229: dump frozen-five membership through Java hap_ll construction.
 * Production arithmetic is delegated. Does not add/remove reads.
 */
public final class HcParityHapLlMembershipDump implements ReadLikelihoodCalculationEngine {

    static final String PFX = "6R229";
    static final String[] QNAMES = {
        "HISEQ1:13:H8G92ADXX:1:2102:7192:18079",
        "HISEQ1:9:H8962ADXX:1:1116:1789:43193",
        "HWI-D00360:6:H81VLADXX:1:1111:8050:25694",
        "HWI-D00360:6:H81VLADXX:1:2102:1733:39463",
        "HWI-D00360:7:H88WKADXX:2:1201:4043:96748"
    };
    static final int[] FLAGS = {99, 163, 163, 163, 163};

    private final ReadLikelihoodCalculationEngine inner;

    public static void installOn(final HaplotypeCallerEngine engine) throws Exception {
        final Field f =
                HaplotypeCallerEngine.class.getDeclaredField("likelihoodCalculationEngine");
        f.setAccessible(true);
        final ReadLikelihoodCalculationEngine inner =
                (ReadLikelihoodCalculationEngine) f.get(engine);
        f.set(engine, new HcParityHapLlMembershipDump(inner));
    }

    private HcParityHapLlMembershipDump(final ReadLikelihoodCalculationEngine inner) {
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
            throw new RuntimeException("6R.229 hap-ll membership dump failed", e);
        }
    }

    private static AlleleLikelihoods<GATKRead, Haplotype> computeAndDump(
            final PairHMMLikelihoodCalculationEngine phmm,
            final List<Haplotype> haplotypeList,
            final SampleList samples,
            final Map<String, List<GATKRead>> perSampleReadList)
            throws Exception {
        final List<GATKRead> orig = new ArrayList<>();
        for (final List<GATKRead> reads : perSampleReadList.values()) {
            orig.addAll(reads);
        }
        dumpStage("pairhmm_input", orig);

        final Method init =
                PairHMMLikelihoodCalculationEngine.class.getDeclaredMethod(
                        "initializePairHMM", List.class, Map.class);
        init.setAccessible(true);
        init.invoke(phmm, haplotypeList, perSampleReadList);

        final AlleleLikelihoods<GATKRead, Haplotype> result =
                new AlleleLikelihoods<>(
                        samples, new IndexedAlleleList<>(haplotypeList), perSampleReadList);
        dumpLl("initial_hap_ll", result);

        final Method modify =
                PairHMMLikelihoodCalculationEngine.class.getDeclaredMethod(
                        "modifyReadQualities", List.class);
        modify.setAccessible(true);
        final Method computeOne =
                PairHMMLikelihoodCalculationEngine.class.getDeclaredMethod(
                        "computeReadLikelihoods", LikelihoodMatrix.class);
        computeOne.setAccessible(true);
        for (int s = 0; s < result.numberOfSamples(); s++) {
            modify.invoke(phmm, result.sampleEvidence(s));
            computeOne.invoke(phmm, result.sampleMatrix(s));
        }
        dumpLl("after_pairhmm_kernel", result);

        final double log10global = fieldDouble(phmm, "log10globalReadMismappingRate");
        final boolean sym = fieldBoolean(phmm, "symmetricallyNormalizeAllelesToReference");
        result.normalizeLikelihoods(log10global, sym);
        dumpLl("after_normalize", result);

        final boolean dynamic = fieldBoolean(phmm, "dynamicDisqualification");
        final double expectedError = fieldDouble(phmm, "expectedErrorRatePerBase");
        final double scale = fieldDouble(phmm, "readDisqualificationScale");
        dumpPoorlyModeledMetrics(result, expectedError, dynamic, scale);
        phmm.filterPoorlyModeledEvidence(result, dynamic, expectedError, scale);
        dumpLl("after_poorly_modeled", result);
        kv("stored_evidence_count", Integer.toString(result.evidenceCount()));
        return result;
    }

    static void dumpLl(final String stage, final AlleleLikelihoods<?, ?> ll) {
        if (ll.numberOfSamples() < 1) {
            dumpStage(stage, Collections.emptyList());
            return;
        }
        @SuppressWarnings("unchecked")
        final List<GATKRead> ev = (List<GATKRead>) ll.sampleEvidence(0);
        dumpStage(stage, ev);
    }

    static void dumpStage(final String stage, final Collection<? extends GATKRead> reads) {
        kv(stage + "_count", Integer.toString(reads.size()));
        kv(stage + "_five", bitmap(reads));
        dumpFiveDetail(stage, reads);
    }

    static String bitmap(final Collection<? extends GATKRead> reads) {
        final StringBuilder sb = new StringBuilder(5);
        for (int i = 0; i < QNAMES.length; i++) {
            sb.append(find(reads, QNAMES[i], FLAGS[i]) != null ? '1' : '0');
        }
        return sb.toString();
    }

    static GATKRead find(
            final Collection<? extends GATKRead> reads, final String qname, final int flags) {
        for (final GATKRead r : reads) {
            if (qname.equals(r.getName()) && r.getFlags() == flags) {
                return r;
            }
        }
        return null;
    }

    static void dumpFiveDetail(final String stage, final Collection<? extends GATKRead> reads) {
        for (int i = 0; i < QNAMES.length; i++) {
            final GATKRead r = find(reads, QNAMES[i], FLAGS[i]);
            if (r == null) {
                kv(
                        stage + "_read",
                        "i="
                                + i
                                + "\tqname="
                                + QNAMES[i]
                                + "\tflags="
                                + FLAGS[i]
                                + "\tpresent=0");
                continue;
            }
            kv(
                    stage + "_read",
                    "i="
                            + i
                            + "\tqname="
                            + r.getName()
                            + "\tflags="
                            + r.getFlags()
                            + "\tpresent=1"
                            + "\tstart="
                            + r.getStart()
                            + "\tend="
                            + r.getEnd()
                            + "\tuStart="
                            + r.getUnclippedStart()
                            + "\tuEnd="
                            + r.getUnclippedEnd()
                            + "\tcigar="
                            + String.valueOf(r.getCigar())
                            + "\tlen="
                            + r.getLength()
                            + "\tunclippedLen="
                            + AlignmentUtils.unclippedReadLength(r)
                            + "\tmapq="
                            + r.getMappingQuality()
                            + "\tmateContig="
                            + String.valueOf(r.getMateContig())
                            + "\tmateStart="
                            + r.getMateStart()
                            + "\tpaired="
                            + r.isPaired()
                            + "\tmateUnmapped="
                            + r.mateIsUnmapped()
                            + "\tunmapped="
                            + r.isUnmapped()
                            + "\tempty="
                            + r.isEmpty()
                            + "\tmatePred="
                            + (ReadFilterLibrary.MATE_ON_SAME_CONTIG_OR_NO_MAPPED_MATE.test(r)
                                    ? "PASS"
                                    : "FAIL"));
        }
    }

    static void dumpFilterPredicates(
            final String stage,
            final Collection<? extends GATKRead> reads,
            final int mqThreshold,
            final String keepRg) {
        for (int i = 0; i < QNAMES.length; i++) {
            final GATKRead r = find(reads, QNAMES[i], FLAGS[i]);
            if (r == null) {
                kv(stage + "_pred", "i=" + i + "\tqname=" + QNAMES[i] + "\tabsent=1");
                continue;
            }
            final int unclipped = AlignmentUtils.unclippedReadLength(r);
            final boolean lenOk = unclipped >= 10;
            final boolean mqOk = r.getMappingQuality() >= mqThreshold;
            final boolean mateOk =
                    ReadFilterLibrary.MATE_ON_SAME_CONTIG_OR_NO_MAPPED_MATE.test(r);
            final boolean rgOk =
                    keepRg == null
                            || (r.getReadGroup() != null && keepRg.equals(r.getReadGroup()));
            final boolean keep = lenOk && mqOk && mateOk && rgOk;
            kv(
                    stage + "_pred",
                    "i="
                            + i
                            + "\tqname="
                            + r.getName()
                            + "\tflags="
                            + r.getFlags()
                            + "\tunclippedLen="
                            + unclipped
                            + "\tlenOk="
                            + lenOk
                            + "\tmapq="
                            + r.getMappingQuality()
                            + "\tmqThresh="
                            + mqThreshold
                            + "\tmqOk="
                            + mqOk
                            + "\tmatePred="
                            + (mateOk ? "PASS" : "FAIL")
                            + "\trgOk="
                            + rgOk
                            + "\tkeep="
                            + keep);
        }
    }

    static void dumpClipAgainstSpan(
            final String stage,
            final Collection<? extends GATKRead> reads,
            final SimpleInterval paddedSpan) {
        kv(
                stage + "_span",
                paddedSpan.getContig() + ":" + paddedSpan.getStart() + "-" + paddedSpan.getEnd());
        for (int i = 0; i < QNAMES.length; i++) {
            final GATKRead r = find(reads, QNAMES[i], FLAGS[i]);
            if (r == null) {
                kv(stage + "_clip", "i=" + i + "\tqname=" + QNAMES[i] + "\tabsent_before_clip=1");
                continue;
            }
            final boolean origOverlap = r.overlaps(paddedSpan);
            final GATKRead clipped =
                    ReadClipper.hardClipToRegion(
                            r, paddedSpan.getStart(), paddedSpan.getEnd());
            final boolean empty = clipped.isEmpty();
            final boolean overlapAfter =
                    !empty && clipped.overlaps(paddedSpan);
            final boolean keep = !empty && overlapAfter;
            kv(
                    stage + "_clip",
                    "i="
                            + i
                            + "\tqname="
                            + r.getName()
                            + "\tflags="
                            + r.getFlags()
                            + "\torigStart="
                            + r.getStart()
                            + "\torigEnd="
                            + r.getEnd()
                            + "\torigCigar="
                            + String.valueOf(r.getCigar())
                            + "\torigOverlap="
                            + origOverlap
                            + "\tclipStart="
                            + clipped.getStart()
                            + "\tclipEnd="
                            + clipped.getEnd()
                            + "\tclipCigar="
                            + String.valueOf(clipped.getCigar())
                            + "\tclipLen="
                            + clipped.getLength()
                            + "\tempty="
                            + empty
                            + "\toverlapAfter="
                            + overlapAfter
                            + "\tkeep="
                            + keep);
        }
    }

    private static void dumpPoorlyModeledMetrics(
            final AlleleLikelihoods<GATKRead, Haplotype> ll,
            final double expectedError,
            final boolean dynamic,
            final double scale) {
        kv("poorly_dynamic", Boolean.toString(dynamic));
        kv("poorly_expected_error", Double.toString(expectedError));
        kv("poorly_scale", Double.toString(scale));
        if (ll.numberOfSamples() < 1) {
            return;
        }
        final LikelihoodMatrix<GATKRead, Haplotype> mx = ll.sampleMatrix(0);
        final List<GATKRead> ev = ll.sampleEvidence(0);
        for (int i = 0; i < QNAMES.length; i++) {
            int idx = -1;
            for (int r = 0; r < ev.size(); r++) {
                if (QNAMES[i].equals(ev.get(r).getName()) && ev.get(r).getFlags() == FLAGS[i]) {
                    idx = r;
                    break;
                }
            }
            if (idx < 0) {
                kv(
                        "poorly_read",
                        "i="
                                + i
                                + "\tqname="
                                + QNAMES[i]
                                + "\tflags="
                                + FLAGS[i]
                                + "\tpresent=0\tkeep=N/A");
                continue;
            }
            final GATKRead read = ev.get(idx);
            final Object hmmTag =
                    read.getTransientAttribute(
                            PairHMMLikelihoodCalculationEngine.HMM_BASE_QUALITIES_TAG);
            final int qlen = hmmTag instanceof byte[] ? ((byte[]) hmmTag).length : read.getLength();
            final double maxErrors = Math.min(2.0, Math.ceil(qlen * expectedError));
            final double threshold = maxErrors * -4.0;
            double maxLl = Double.NEGATIVE_INFINITY;
            for (int a = 0; a < mx.numberOfAlleles(); a++) {
                maxLl = Math.max(maxLl, mx.get(a, idx));
            }
            final boolean keep = !(maxLl < threshold);
            kv(
                    "poorly_read",
                    "i="
                            + i
                            + "\tqname="
                            + QNAMES[i]
                            + "\tflags="
                            + FLAGS[i]
                            + "\tpresent=1"
                            + "\tqlen="
                            + qlen
                            + "\tmax_ll="
                            + String.format(Locale.US, "%.6f", maxLl)
                            + "\tthresh="
                            + String.format(Locale.US, "%.6f", threshold)
                            + "\tkeep="
                            + (keep ? "KEEP" : "DROP"));
        }
    }

    static AssemblyRegion copyRegion(final AssemblyRegion src) {
        final AssemblyRegion dst =
                new AssemblyRegion(
                        src.getSpan(), src.getPaddedSpan(), src.isActive(), src.getHeader());
        for (final GATKRead read : src.getReads()) {
            dst.add(read.copy());
        }
        return dst;
    }

    @SuppressWarnings("unchecked")
    static void dumpPrefixStages(
            final HaplotypeCallerEngine engine,
            final AssemblyRegion regionCopy,
            final HaplotypeCallerArgumentCollection hcArgs,
            final SampleList samplesList,
            final Logger logger,
            final CachingIndexedFastaSequenceFile refReader,
            final ReadThreadingAssembler assembler,
            final SmithWatermanAligner aligner,
            final AssemblyRegionTrimmer trimmer,
            final ReferenceContext refCtx)
            throws Exception {
        dumpStage("iterator_region_reads", regionCopy.getReads());
        kv(
                "iterator_span",
                regionCopy.getContig()
                        + ":"
                        + regionCopy.getStart()
                        + "-"
                        + regionCopy.getEnd());
        kv(
                "iterator_padded",
                regionCopy.getPaddedSpan().getContig()
                        + ":"
                        + regionCopy.getPaddedSpan().getStart()
                        + "-"
                        + regionCopy.getPaddedSpan().getEnd());
        kv("mq_threshold", Integer.toString(hcArgs.mappingQualityThreshold));
        kv("keep_rg", hcArgs.keepRG == null ? "null" : hcArgs.keepRG);

        final AssemblyResultSet ars =
                AssemblyBasedCallerUtils.assembleReads(
                        regionCopy,
                        Collections.emptyList(),
                        hcArgs,
                        regionCopy.getHeader(),
                        samplesList,
                        logger,
                        refReader,
                        assembler,
                        aligner,
                        !hcArgs.doNotCorrectOverlappingBaseQualities,
                        hcArgs.fbargs,
                        false);
        dumpStage("after_assemble_finalize", regionCopy.getReads());

        final TreeSet<VariantContext> sortedUnion =
                new TreeSet<>(
                        Comparator.comparingInt(VariantContext::getStart)
                                .thenComparingInt(VariantContext::getEnd)
                                .thenComparing(vc -> vc.getReference().toString())
                                .thenComparing(vc -> vc.getAlternateAlleles().toString()));
        sortedUnion.addAll(ars.getVariationEvents(hcArgs.maxMnpDistance));
        final AssemblyRegionTrimmer.Result trimmingResult =
                trimmer.trim(regionCopy, sortedUnion, refCtx);
        kv("trim_variation_present", Boolean.toString(trimmingResult.isVariationPresent()));
        if (!trimmingResult.isVariationPresent()) {
            kv("trim", "skipped_no_variation");
            return;
        }
        final AssemblyRegion variantRegion = trimmingResult.getVariantRegion();
        kv(
                "variant_span",
                variantRegion.getContig()
                        + ":"
                        + variantRegion.getStart()
                        + "-"
                        + variantRegion.getEnd());
        kv(
                "variant_padded",
                variantRegion.getPaddedSpan().getContig()
                        + ":"
                        + variantRegion.getPaddedSpan().getStart()
                        + "-"
                        + variantRegion.getPaddedSpan().getEnd());
        dumpClipAgainstSpan(
                "trim_hardclip", regionCopy.getReads(), variantRegion.getPaddedSpan());
        dumpStage("after_trim", variantRegion.getReads());

        final List<GATKRead> stubs = new ArrayList<>();
        for (final GATKRead r : variantRegion.getReads()) {
            if (AlignmentUtils.unclippedReadLength(r)
                    < AssemblyBasedCallerUtils.MINIMUM_READ_LENGTH_AFTER_TRIMMING) {
                stubs.add(r);
            }
        }
        kv("stub_n", Integer.toString(stubs.size()));
        dumpStage("stubs", stubs);
        variantRegion.removeAll(stubs);
        dumpStage("after_stub", variantRegion.getReads());
        dumpFilterPredicates(
                "after_stub",
                variantRegion.getReads(),
                hcArgs.mappingQualityThreshold,
                hcArgs.keepRG);

        final Method fnpr =
                HaplotypeCallerEngine.class.getDeclaredMethod(
                        "filterNonPassingReads", AssemblyRegion.class);
        fnpr.setAccessible(true);
        final Collection<?> removed = (Collection<?>) fnpr.invoke(engine, variantRegion);
        kv("filter_non_passing_removed_n", Integer.toString(removed.size()));
        dumpStage("filter_non_passing_removed", (Collection<GATKRead>) removed);
        dumpStage("after_filter_non_passing", variantRegion.getReads());
        dumpFilterPredicates(
                "after_filter_non_passing",
                variantRegion.getReads(),
                hcArgs.mappingQualityThreshold,
                hcArgs.keepRG);
    }

    public static void kv(final String key, final String value) {
        System.out.println(PFX + "\t" + key + "\t" + value);
    }

    private static double fieldDouble(final Object o, final String name) throws Exception {
        final Field f = o.getClass().getDeclaredField(name);
        f.setAccessible(true);
        return f.getDouble(o);
    }

    private static boolean fieldBoolean(final Object o, final String name) throws Exception {
        final Field f = o.getClass().getDeclaredField(name);
        f.setAccessible(true);
        return f.getBoolean(o);
    }
}
