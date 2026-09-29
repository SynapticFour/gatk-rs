import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import java.lang.reflect.Method;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Collections;
import java.util.List;
import org.broadinstitute.hellbender.engine.AssemblyRegion;
import org.broadinstitute.hellbender.engine.AssemblyRegionIterator;
import org.broadinstitute.hellbender.engine.MultiIntervalLocalReadShard;
import org.broadinstitute.hellbender.engine.ReadsPathDataSource;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.engine.spark.AssemblyRegionArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.annotator.VariantAnnotatorEngine;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCallerArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCallerEngine;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.hellbender.utils.read.GATKRead;

/**
 * 6R.305: trace one read through Java HC evidence selection up to
 * filterNonPassingReads. Does not compute PairHMM or genotype likelihoods.
 */
public final class Membership305 {
    private static final String QNAME = "HISEQ1:11:H8GV6ADXX:2:1103:14252:55237";
    private static final int FLAGS = 97;
    private static final int LOC = 29455649;

    public static void main(final String[] args) throws Exception {
        if (args.length != 2) {
            throw new IllegalArgumentException("usage: Membership305 ref bam");
        }
        final Path bam = Paths.get(args[1]);
        final ReadsPathDataSource readsSource = new ReadsPathDataSource(bam);
        final SAMFileHeader header = readsSource.getHeader();
        final Path ref = Paths.get(args[0]);
        final CachingIndexedFastaSequenceFile refReader = new CachingIndexedFastaSequenceFile(ref);
        final ReferenceDataSource reference = ReferenceDataSource.of(ref);
        final AssemblyRegionArgumentCollection asmArgs = new AssemblyRegionArgumentCollection();
        asmArgs.assemblyRegionPadding = 100;
        final HaplotypeCallerArgumentCollection hcArgs = new HaplotypeCallerArgumentCollection();
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
                        Collections.singletonList(interval), 100, readsSource);
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
            kv("iterator_n", Integer.toString(region.getReads().size()));
            final GATKRead hit = find(region.getReads());
            kv("iterator_present", Boolean.toString(hit != null));
            if (hit == null) {
                break;
            }
            kv("name", hit.getName());
            kv("flags", Integer.toString(hit.getFlags()));
            kv("contig", String.valueOf(hit.getContig()));
            kv("start", Integer.toString(hit.getStart()));
            kv("end", Integer.toString(hit.getEnd()));
            kv("cigar", String.valueOf(hit.getCigar()));
            kv("mapq", Integer.toString(hit.getMappingQuality()));
            kv("mate_contig", String.valueOf(hit.getMateContig()));
            kv("mate_start", Integer.toString(hit.getMateStart()));
            kv("paired", Boolean.toString(hit.isPaired()));
            kv("proper_pair", Boolean.toString(hit.isProperlyPaired()));
            kv("unmapped", Boolean.toString(hit.isUnmapped()));
            kv("mate_unmapped", Boolean.toString(hit.mateIsUnmapped()));
            kv("reverse", Boolean.toString(hit.isReverseStrand()));
            kv("mate_reverse", Boolean.toString(hit.mateIsReverseStrand()));
            kv("secondary", Boolean.toString(hit.isSecondaryAlignment()));
            kv("supplementary", Boolean.toString(hit.isSupplementaryAlignment()));
            kv("duplicate", Boolean.toString(hit.isDuplicate()));
            kv("vendor_fail", Boolean.toString(hit.failsVendorQualityCheck()));
            kv("read1", Boolean.toString(hit.isFirstOfPair()));
            for (final ReadFilter f : readFilters) {
                kv(
                        "standard_filter",
                        f.getClass().getSimpleName() + "\t" + f.test(hit));
            }
            final boolean mate =
                    ReadFilterLibrary.MATE_ON_SAME_CONTIG_OR_NO_MAPPED_MATE.test(hit);
            kv("mate_on_same_contig", Boolean.toString(mate));
            kv("mq_threshold", Integer.toString(hcArgs.mappingQualityThreshold));
            final int before = region.getReads().size();
            final Method fnpr =
                    HaplotypeCallerEngine.class.getDeclaredMethod(
                            "filterNonPassingReads", AssemblyRegion.class);
            fnpr.setAccessible(true);
            final java.util.Collection<?> removed =
                    (java.util.Collection<?>) fnpr.invoke(engine, region);
            boolean removedHit = false;
            for (final Object o : removed) {
                final GATKRead r = (GATKRead) o;
                if (QNAME.equals(r.getName()) && r.getFlags() == FLAGS) {
                    removedHit = true;
                }
            }
            kv("filter_before", Integer.toString(before));
            kv("filter_removed_n", Integer.toString(removed.size()));
            kv("filter_after", Integer.toString(region.getReads().size()));
            kv("filter_removed_this_read", Boolean.toString(removedHit));
            kv("filter_still_present", Boolean.toString(find(region.getReads()) != null));
            found = true;
            break;
        }
        engine.shutdown();
        readsSource.close();
        if (!found) {
            throw new IllegalStateException("read not traced");
        }
    }

    private static GATKRead find(final List<GATKRead> reads) {
        for (final GATKRead r : reads) {
            if (QNAME.equals(r.getName()) && r.getFlags() == FLAGS) {
                return r;
            }
        }
        return null;
    }

    private static void kv(final String key, final String value) {
        System.out.println("6R305\t" + key + "\t" + value);
    }
}
