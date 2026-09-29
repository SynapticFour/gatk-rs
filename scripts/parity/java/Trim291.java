import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import java.lang.reflect.Field;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.TreeSet;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
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
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AssemblyBasedCallerUtils;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AssemblyRegionTrimmer;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AssemblyResultSet;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCallerArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCallerEngine;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.PairHMMNativeArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading.ReadThreadingAssembler;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;
import org.broadinstitute.hellbender.utils.pairhmm.PairHMM;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanAligner;
import org.broadinstitute.hellbender.tools.walkers.annotator.TandemRepeat;
import org.apache.commons.lang3.tuple.Pair;
import htsjdk.variant.variantcontext.VariantContext;

/**
 * 6R.291: trim-span inputs for one read, from the raw BAM record through
 * {@code AssemblyRegion.trim}'s {@code hardClipToRegion}. Stops before PairHMM.
 */
public final class Trim291 {
    private static final String QNAME = "HISEQ1:11:H8GV6ADXX:1:2116:18670:99941";
    private static final int FLAGS = 99;
    private static final int LOC = 29455649;
    private static final int PADDING = 100;

    public static void main(final String[] args) throws Exception {
        if (args.length != 2) {
            throw new IllegalArgumentException("usage: Trim291 ref bam");
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
        final VariantAnnotatorEngine annotationEngine =
                new VariantAnnotatorEngine(
                        Collections.emptyList(), null, Collections.emptyList(), false, false);
        final HaplotypeCallerEngine engine =
                new HaplotypeCallerEngine(
                        hcArgs, asmArgs, false, false, header, refReader, annotationEngine);
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
            dump(engine, region, refReader, reference);
            found = true;
            break;
        }
        engine.shutdown();
        readsSource.close();
        if (!found) {
            throw new IllegalStateException("no active region covering " + LOC);
        }
    }

    private static void dump(
            final HaplotypeCallerEngine engine,
            final AssemblyRegion region,
            final CachingIndexedFastaSequenceFile refReader,
            final ReferenceDataSource reference)
            throws Exception {
        final HaplotypeCallerArgumentCollection hcArgs = field(engine, "hcArgs");
        final SampleList samples = field(engine, "samplesList");
        final AssemblyRegionTrimmer trimmer = field(engine, "trimmer");
        final ReadThreadingAssembler assembler = field(engine, "assemblyEngine");
        final SmithWatermanAligner aligner = field(engine, "aligner");
        final Logger logger = LogManager.getLogger(Trim291.class);

        kv("active_span", region.getContig() + ":" + region.getStart() + "-" + region.getEnd());
        kv(
                "padded_span",
                region.getPaddedSpan().getContig()
                        + ":"
                        + region.getPaddedSpan().getStart()
                        + "-"
                        + region.getPaddedSpan().getEnd());
        kv("min_base_qual", Integer.toString(hcArgs.minBaseQualityScore));
        kv("dont_use_soft_clipped_bases", Boolean.toString(hcArgs.dontUseSoftClippedBases));
        kv("soft_clip_low_qual_ends", Boolean.toString(hcArgs.softClipLowQualityEnds));

        final GATKRead raw = find(region.getReads());
        if (raw == null) {
            throw new IllegalStateException("target read absent from iterator region");
        }
        dumpRead("iterator", raw);

        final AssemblyRegion finalized = copyRegion(region);
        AssemblyBasedCallerUtils.finalizeRegion(
                finalized,
                hcArgs.assemblerArgs.errorCorrectReads,
                hcArgs.dontUseSoftClippedBases,
                (byte) (hcArgs.minBaseQualityScore - 1),
                region.getHeader(),
                samples,
                !hcArgs.doNotCorrectOverlappingBaseQualities,
                hcArgs.softClipLowQualityEnds,
                hcArgs.overrideSoftclipFragmentCheck,
                hcArgs.fbargs,
                false);
        final GATKRead afterFinalize = find(finalized.getReads());
        if (afterFinalize == null) {
            kv("after_finalize", "ABSENT");
        } else {
            dumpRead("after_finalize", afterFinalize);
        }

        final AssemblyRegion forAssembly = copyRegion(region);
        final AssemblyResultSet ars =
                AssemblyBasedCallerUtils.assembleReads(
                        forAssembly,
                        Collections.emptyList(),
                        hcArgs,
                        region.getHeader(),
                        samples,
                        logger,
                        refReader,
                        assembler,
                        aligner,
                        !hcArgs.doNotCorrectOverlappingBaseQualities,
                        hcArgs.fbargs,
                        false);
        final GATKRead afterAssemble = find(forAssembly.getReads());
        if (afterAssemble == null) {
            kv("after_assemble_finalize", "ABSENT");
        } else {
            dumpRead("after_assemble_finalize", afterAssemble);
        }

        final TreeSet<VariantContext> events =
                new TreeSet<>(
                        Comparator.comparingInt(VariantContext::getStart)
                                .thenComparingInt(VariantContext::getEnd)
                                .thenComparing(vc -> vc.getReference().toString())
                                .thenComparing(vc -> vc.getAlternateAlleles().toString()));
        events.addAll(ars.getVariationEvents(hcArgs.maxMnpDistance));
        kv("n_events", Integer.toString(events.size()));
        if (!events.isEmpty()) {
            final VariantContext left = events.first();
            kv(
                    "leftmost_event",
                    left.getContig()
                            + ":"
                            + left.getStart()
                            + "-"
                            + left.getEnd()
                            + " "
                            + left.getReference().getDisplayString()
                            + ">"
                            + left.getAlternateAllele(0).getDisplayString()
                            + " indel="
                            + left.isIndel());
        }
        final ReferenceContext refCtx =
                new ReferenceContext(reference, region.getSpan(), PADDING, PADDING);
        final AssemblyRegionTrimmer.Result trimming = trimmer.trim(forAssembly, events, refCtx);
        final Object asmArgsObj = field(trimmer, "assemblyRegionArgs");
        final AssemblyRegionArgumentCollection liveArgs = (AssemblyRegionArgumentCollection) asmArgsObj;
        kv("snp_padding", Integer.toString(liveArgs.snpPaddingForGenotyping));
        kv("indel_padding", Integer.toString(liveArgs.indelPaddingForGenotyping));
        kv("str_padding", Integer.toString(liveArgs.strPaddingForGenotyping));
        kv("legacy", Boolean.toString(liveArgs.enableLegacyAssemblyRegionTrimming));
        int minStart = Integer.MAX_VALUE;
        int maxEnd = Integer.MIN_VALUE;
        int nOverlap = 0;
        for (final VariantContext vc : events) {
            final boolean overlap = vc.getStart() <= region.getEnd() && vc.getEnd() >= region.getStart();
            if (!overlap) {
                kv("event_outside", vc.getContig()+":"+vc.getStart()+"-"+vc.getEnd()+" "+vc.getReference().getDisplayString()+">"+vc.getAlternateAllele(0).getDisplayString());
                continue;
            }
            nOverlap++;
            minStart = Math.min(minStart, vc.getStart());
            maxEnd = Math.max(maxEnd, vc.getEnd());
        }
        kv("n_overlap", Integer.toString(nOverlap));
        kv("raw_min_start", Integer.toString(minStart));
        kv("raw_max_end", Integer.toString(maxEnd));
        int runMin = minStart;
        int runMax = maxEnd;
        int i = 0;
        for (final VariantContext vc : events) {
            if (!(vc.getStart() <= region.getEnd() && vc.getEnd() >= region.getStart())) {
                continue;
            }
            int padding = liveArgs.snpPaddingForGenotyping;
            String str = ".";
            if (vc.isIndel()) {
                padding = liveArgs.indelPaddingForGenotyping;
                final Pair<java.util.List<Integer>, byte[]> rep = TandemRepeat.getNumTandemRepeatUnits(refCtx, vc);
                if (rep != null && rep.getRight() != null) {
                    final int repeatLength = rep.getRight().length;
                    int most = 0;
                    for (final int n : rep.getLeft()) {
                        most = Math.max(most, n);
                    }
                    final int longest = most * repeatLength;
                    padding = liveArgs.strPaddingForGenotyping + longest;
                    str = "unit="+repeatLength+" repeats="+most+" longest="+longest;
                }
            }
            final int left = Math.max(vc.getStart() - padding, 1);
            final int right = vc.getEnd() + padding;
            runMin = Math.min(runMin, left);
            runMax = Math.max(runMax, right);
            kv("event", "i="+i
                +"\tstart="+vc.getStart()
                +"\tend="+vc.getEnd()
                +"\tref="+vc.getReference().getDisplayString()
                +"\talt="+vc.getAlternateAllele(0).getDisplayString()
                +"\tindel="+vc.isIndel()
                +"\tpadding="+padding
                +"\tstr="+str
                +"\tleft="+left
                +"\tright="+right
                +"\trun_min="+runMin
                +"\trun_max="+runMax);
            i++;
        }
        kv("padded_before_intersect", runMin+"-"+runMax);

        kv("trim_variation_present", Boolean.toString(trimming.isVariationPresent()));
        if (!trimming.isVariationPresent()) {
            throw new IllegalStateException("trim found no variation");
        }
        final AssemblyRegion variantRegion = trimming.getVariantRegion();
        kv(
                "variant_span",
                variantRegion.getContig()
                        + ":"
                        + variantRegion.getStart()
                        + "-"
                        + variantRegion.getEnd());
        kv(
                "variant_padded_span",
                variantRegion.getPaddedSpan().getContig()
                        + ":"
                        + variantRegion.getPaddedSpan().getStart()
                        + "-"
                        + variantRegion.getPaddedSpan().getEnd());
        final GATKRead afterTrim = find(variantRegion.getReads());
        if (afterTrim == null) {
            kv("after_trim_hardclip", "ABSENT");
        } else {
            dumpRead("after_trim_hardclip", afterTrim);
        }
    }

    private static AssemblyRegion copyRegion(final AssemblyRegion src) {
        final AssemblyRegion dst =
                new AssemblyRegion(
                        src.getSpan(), src.getPaddedSpan(), src.isActive(), src.getHeader());
        for (final GATKRead read : src.getReads()) {
            dst.add(read.copy());
        }
        return dst;
    }

    private static GATKRead find(final Iterable<GATKRead> reads) {
        for (final GATKRead r : reads) {
            if (QNAME.equals(r.getName()) && r.getFlags() == FLAGS) {
                return r;
            }
        }
        return null;
    }

    private static void dumpRead(final String stage, final GATKRead r) {
        final String bases = new String(r.getBases());
        kv(
                stage,
                "qname="
                        + r.getName()
                        + "\tflags="
                        + r.getFlags()
                        + "\tcigar="
                        + r.getCigar()
                        + "\tstart="
                        + r.getStart()
                        + "\tend="
                        + r.getEnd()
                        + "\tlen="
                        + r.getLength()
                        + "\treverse="
                        + r.isReverseStrand()
                        + "\tbases="
                        + bases);
    }

    @SuppressWarnings("unchecked")
    private static <T> T field(final Object o, final String name) throws Exception {
        final Field f = o.getClass().getDeclaredField(name);
        f.setAccessible(true);
        return (T) f.get(o);
    }

    private static void kv(final String key, final String value) {
        System.out.println("6R291\t" + key + "\t" + value);
    }
}
