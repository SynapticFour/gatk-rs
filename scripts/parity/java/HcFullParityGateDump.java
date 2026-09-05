import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.utils.GenomeLoc;
import org.broadinstitute.hellbender.utils.GenomeLocParser;
import org.broadinstitute.hellbender.engine.AssemblyRegion;
import org.broadinstitute.hellbender.engine.AlignmentAndReferenceContext;
import org.broadinstitute.hellbender.engine.AlignmentContext;
import org.broadinstitute.hellbender.engine.AssemblyRegionIterator;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.FeatureInput;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.engine.FeatureManager;
import org.broadinstitute.hellbender.engine.MultiIntervalLocalReadShard;
import org.broadinstitute.hellbender.engine.ReadsPathDataSource;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.engine.filters.CountingReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.MappingQualityReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.engine.filters.WellformedReadFilter;
import org.broadinstitute.hellbender.engine.spark.AssemblyRegionArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.annotator.VariantAnnotatorEngine;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCallerArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.PairHMMNativeArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AssemblyRegionTrimmer;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AssemblyBasedCallerUtils;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AssemblyResultSet;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading.ReadThreadingAssembler;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCallerEngine;
import com.synapticfour.gatkrs.parity.HcParityAlleleSubsetting;
import com.synapticfour.gatkrs.parity.HcParityGvcfMerger;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.ReferenceConfidenceModel;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.RefVsAnyResult;
import org.broadinstitute.hellbender.tools.walkers.genotyper.MinimalGenotypingEngine;
import org.broadinstitute.hellbender.utils.MathUtils;
import java.lang.reflect.Field;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.PileupReadErrorCorrector;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.ReadThreadingAssemblerArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.BaseEdge;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.SeqGraph;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.SeqVertex;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.AdaptiveChainPruner;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.ChainPruner;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.GraphBasedKBestHaplotypeFinder;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.JunctionTreeKBestHaplotypeFinder;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.KBestHaplotype;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.LowWeightChainPruner;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.MultiSampleEdge;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading.MultiDeBruijnVertex;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading.AbstractReadThreadingGraph;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading.JunctionTreeLinkedDeBruijnGraph;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading.ReadThreadingGraph;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWParameters;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanAligner;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.activityprofile.ActivityProfile;
import org.broadinstitute.hellbender.utils.activityprofile.ActivityProfileState;
import org.broadinstitute.hellbender.utils.activityprofile.BandPassActivityProfile;
import org.broadinstitute.hellbender.transformers.ReadTransformer;
import org.broadinstitute.hellbender.utils.downsampling.AlleleBiasedDownsamplingUtils;
import org.broadinstitute.hellbender.utils.downsampling.DownsamplingMethod;
import org.broadinstitute.hellbender.utils.downsampling.PositionalDownsampler;
import org.broadinstitute.hellbender.utils.downsampling.ReadsDownsampler;
import org.broadinstitute.hellbender.utils.pileup.PileupElement;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;
import org.broadinstitute.hellbender.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.hellbender.utils.iterators.IntervalLocusIterator;
import org.broadinstitute.hellbender.utils.iterators.ReadCachingIterator;
import org.broadinstitute.hellbender.utils.locusiterator.IntervalAlignmentContextIterator;
import org.broadinstitute.hellbender.utils.locusiterator.LocusIteratorByState;
import org.broadinstitute.hellbender.utils.clipping.ReadClipper;
import org.broadinstitute.hellbender.utils.fragments.FragmentCollection;
import org.broadinstitute.hellbender.utils.read.AlignmentUtils;
import org.broadinstitute.hellbender.utils.read.CigarUtils;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadCoordinateComparator;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;
import org.broadinstitute.hellbender.utils.genotyper.IndexedSampleList;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.tribble.Feature;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.lang.reflect.Method;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Deque;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.OptionalInt;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.stream.Collectors;

/**
 * Java L2 emitter for hc-full-parity (pinned GATK 4.4.0.0). Subcommands mirror
 * {@code gatk-haplotypecaller/examples/hc_full_parity_gate_dump.rs}.
 */
public final class HcFullParityGateDump {

    private static final String HC_READ_FILTER_COUNT_SECTION = "---HC_READ_FILTER_COUNTS---";

    private static final int DEFAULT_PADDING = 100;

    private HcFullParityGateDump() {}

    public static void main(String[] args) throws Exception {
        if (args.length < 1) {
            usage();
            return;
        }
        final String cmd = args[0];
        final String[] rest = Arrays.copyOfRange(args, 1, args.length);
        switch (cmd) {
            case "read-shards":
                readShards(rest);
                break;
            case "read-filters":
                readFilters(rest);
                break;
            case "read-shard-pipeline":
                readShardPipeline(rest);
                break;
            case "read-pre-softclip":
                readPreSoftclip(rest);
                break;
            case "read-pre-len":
                readPreLen(rest);
                break;
            case "read-pre-mq":
                readPreMq(rest);
                break;
            case "read-pre-overlap":
                readPreOverlap(rest);
                break;
            case "assembly-regions":
                assemblyRegions(rest, false);
                break;
            case "assembly-regions-force-active":
                assemblyRegions(rest, true);
                break;
            case "locus-pileup":
                locusPileup(rest);
                break;
            case "assembly-region-reads":
                assemblyRegionReads(rest);
                break;
            case "assembly-region-reference":
                assemblyRegionReference(rest);
                break;
            case "assembly-region-features":
                assemblyRegionFeatures(rest);
                break;
            case "assembly-region-trim":
                assemblyRegionTrim(rest);
                break;
            case "assembly-region-pileup-track":
                assemblyRegionPileupTrack(rest);
                break;
            case "assembly-region-haplotypes":
                assemblyRegionHaplotypes(rest);
                break;
            case "assembly-region-kmer-probe":
                assemblyRegionKmerProbe(rest);
                break;
            case "assembly-region-assembly-stages":
                assemblyRegionAssemblyStages(rest);
                break;
            case "assembly-region-assembly-stages-finalize":
                assemblyRegionAssemblyStagesFinalize(rest);
                break;
            case "assembly-region-finalize-reads":
                assemblyRegionFinalizeReads(rest);
                break;
            case "assembly-region-kbest-paths":
                assemblyRegionKbestPaths(rest);
                break;
            case "assembly-region-seqgraph-edges":
                assemblyRegionSeqgraphEdges(rest);
                break;
            case "assembly-region-pairhmm-likelihoods":
                assemblyRegionPairhmmLikelihoods(rest);
                break;
            case "pairhmm-bq-cap":
                pairhmmBqCap(rest);
                break;
            case "pairhmm-haplotype-filter":
                pairhmmHaplotypeFilter(rest);
                break;
            case "apply-summary":
            case "walker-traversal-summary":
                applySummary(rest);
                break;
            case "raw-activity":
                rawActivity(rest);
                break;
            case "raw-activity-force":
                rawActivityForce(rest);
                break;
            case "smoothed-activity":
                smoothedActivity(rest);
                break;
            case "active-locus":
                activeLocus(rest);
                break;
            case "genotype-likelihoods":
                genotypeLikelihoodActivity(rest);
                break;
            case "assembly-graph":
                assemblyGraph(rest);
                break;
            case "assembly-graph-multi":
                assemblyGraphMulti(rest);
                break;
            case "assembly-graph-summary":
                assemblyGraphSummary(rest);
                break;
            case "assembly-graph-dangling-summary":
                assemblyGraphDanglingSummary(rest);
                break;
            case "assembly-graph-non-unique-summary":
                assemblyGraphNonUniqueSummary(rest);
                break;
            case "assembly-haplotype-cigars":
                assemblyHaplotypeCigars(rest);
                break;
            case "assembly-haplotypes":
                assemblyHaplotypes(rest);
                break;
            case "assembly-kbest-paths":
                assemblyKbestPaths(rest);
                break;
            case "assembly-haplotypes-cap":
                assemblyHaplotypesCap(rest);
                break;
            case "assembly-haplotypes-production":
                assemblyHaplotypesProduction(rest);
                break;
            case "assembly-junction-haplotypes":
                assemblyJunctionHaplotypes(rest);
                break;
            case "read-error-correction":
                readErrorCorrection(rest);
                break;
            case "assembly-seqgraph-summary":
                assemblySeqGraphSummary(rest);
                break;
            case "assembly-assemble":
                assemblyAssemble(rest);
                break;
            case "pairhmm-likelihoods":
                pairhmmLikelihoods(rest);
                break;
            case "pairhmm-native-likelihoods":
                pairhmmNativeLikelihoods(rest);
                break;
            case "genotyping-aggregate":
                genotypingAggregate(rest);
                break;
            case "genotype-format":
                genotypeFormat(rest);
                break;
            case "annotate-core":
                annotateCore(rest);
                break;
            case "annotation-manifest":
                HcParityCoreAnnotations.dumpAnnotationManifest();
                break;
            case "variant-vcf-from-gl-ad":
                variantVcfFromGlAd(rest);
                break;
            case "variant-format-from-gl-ad":
                variantFormatFromGlAd(rest);
                break;
            case "call-region-vcf":
                callRegionVcf(rest);
                break;
            case "ad-annotation-call":
                adAnnotationCall(rest);
                break;
            case "filter-poorly-modeled-call":
                filterPoorlyModeledCall(rest);
                break;
            case "filter-poorly-modeled-call-double":
                filterPoorlyModeledCallDouble(rest);
                break;
            case "call-region-format":
                System.err.println(
                        "call-region-format: use call-region-vcf + variant-format-from-gl-ad for parity");
                System.exit(2);
                break;
            case "af-em":
                HcParityDeferredGates.dumpAfEmFixture(Paths.get(rest[0]));
                break;
            case "genotype-limits":
                HcParityDeferredGates.dumpGenotypeLimits(
                        Integer.parseInt(rest[0]), Integer.parseInt(rest[1]));
                break;
            case "genotype-phasing":
                HcParityDeferredGates.dumpGenotypePhasing(
                        rest[0],
                        "1".equals(rest[1]),
                        "-".equals(rest[2]) ? null : Integer.parseInt(rest[2]));
                break;
            case "force-calling-genotype":
                HcParityDeferredGates.dumpForceCallingGenotype(
                        Paths.get(rest[0]), rest[1], Long.parseLong(rest[2]), "1".equals(rest[3]));
                break;
            case "allele-subsetting":
                HcParityDeferredGates.dumpAlleleSubsetting(
                        rest[0], rest[1], Integer.parseInt(rest[2]));
                break;
            case "subset-alleles-pl":
                HcParityDeferredGates.dumpSubsetAllelesPl(java.nio.file.Paths.get(rest[0]));
                break;
            case "gvcf-header":
                HcParityDeferredGates.dumpGvcfHeader(rest[0], Long.parseLong(rest[1]));
                break;
            case "gvcf-writer-blocks":
                HcParityDeferredGates.dumpGvcfWriterBlocks(java.nio.file.Paths.get(rest[0]));
                break;
            case "subset-alleles-vc":
                HcParityDeferredGates.dumpSubsetAllelesVc(java.nio.file.Paths.get(rest[0]));
                break;
            case "subset-alleles-integration":
                HcParityDeferredGates.dumpSubsetAllelesIntegration(
                        rest[0], rest[1], Integer.parseInt(rest[2]), java.nio.file.Paths.get(rest[3]));
                break;
            case "standard-annotations":
                HcParityDeferredGates.dumpStandardAnnotations(
                        Integer.parseInt(rest[0]),
                        Integer.parseInt(rest[1]),
                        Integer.parseInt(rest[2]),
                        Integer.parseInt(rest[3]),
                        Double.parseDouble(rest[4]),
                        Integer.parseInt(rest[5]),
                        rest[6],
                        rest[7],
                        rest.length > 8 ? rest[8] : "-",
                        rest.length > 9 ? rest[9] : "-",
                        rest.length > 10 ? rest[10] : "-",
                        rest.length > 11 ? rest[11] : "-");
                break;
            case "as-annotations":
                HcParityDeferredGates.dumpAsAnnotations(
                        Double.parseDouble(rest[0]), Double.parseDouble(rest[1]));
                break;
            case "excess-het":
                HcParityDeferredGates.dumpExcessHet(
                        Integer.parseInt(rest[0]),
                        Integer.parseInt(rest[1]),
                        Integer.parseInt(rest[2]));
                break;
            case "depth-per-sample-hc":
                HcParityDeferredGates.dumpDepthPerSampleHc(rest[0]);
                break;
            case "annotation-plugin":
                HcParityDeferredGates.dumpAnnotationPlugin(
                        rest[0],
                        Integer.parseInt(rest[1]),
                        Integer.parseInt(rest[2]),
                        Integer.parseInt(rest[3]),
                        Integer.parseInt(rest[4]),
                        Double.parseDouble(rest[5]),
                        Integer.parseInt(rest[6]),
                        rest[7],
                        rest[8]);
                break;
            case "emit-mode-decision":
                HcParityDeferredGates.dumpEmitModeDecision(
                        rest[0], "1".equals(rest[1]), Integer.parseInt(rest[2]));
                break;
            case "bamout-stub":
                HcParityDeferredGates.dumpBamoutStub("1".equals(rest[0]), Integer.parseInt(rest[1]));
                break;
            case "dragen-mode-branch":
                HcParityDeferredGates.dumpDragenModeBranch();
                break;
            case "gvcf-l5-merged":
                HcParityDeferredGates.dumpGvcfL5Merged(Paths.get(rest[0]));
                break;
            case "dragstr-calibration":
                HcParityDeferredGates.dumpDragstrCalibration("1".equals(rest[0]));
                break;
            case "assembly-debug-stub":
                // args: failure_bam graph_dot
                System.out.println("assembly_failure_bam_enabled\t" + "1".equals(rest[0]));
                System.out.println("graph_dot_enabled\t" + "1".equals(rest[1]));
                System.out.println("dump_ready\tfalse");
                break;
            case "assembly-region-genotype":
                assemblyRegionGenotype(rest);
                break;
            case "assembly-region-genotype-subset":
                assemblyRegionGenotypeSubset(rest);
                break;
            case "gvcf-merge-ref-confidence":
                if (rest.length < 1) {
                    usage();
                }
                HcParityGvcfMerger.dumpMergeCase(rest[0]);
                break;
            case "reference-confidence-locus":
                referenceConfidenceLocus(rest);
                break;
            case "inactive-reference-model":
                inactiveReferenceModel(rest);
                break;
            case "ploidy-resolution":
                ploidyResolution(rest);
                break;
            case "downsample-positional":
                downsamplePositional(rest);
                break;
            case "allele-biased-target-counts":
                alleleBiasedTargetCounts(rest);
                break;
            case "allele-biased-evidence":
                alleleBiasedEvidence(rest);
                break;
            case "raw-activity-contam":
                rawActivityContam(rest);
                break;
            case "soft-clip-mean":
                softClipMean(rest);
                break;
            default:
                usage();
        }
    }

    private static void usage() {
        System.err.println(
                "Usage: HcFullParityGateDump <read-shards|read-filters|read-shard-pipeline|assembly-regions|assembly-regions-force-active|"
                        + "assembly-region-reads|assembly-region-reference|assembly-region-features|"
                        + "assembly-region-trim|assembly-region-pileup-track|assembly-region-haplotypes|"
                        + "assembly-region-kmer-probe|assembly-region-assembly-stages|assembly-region-assembly-stages-finalize|"
                        + "assembly-region-finalize-reads|assembly-region-kbest-paths|assembly-region-seqgraph-edges|"
                        + "assembly-region-pairhmm-likelihoods|pairhmm-bq-cap|pairhmm-haplotype-filter|"
                        + "locus-pileup|apply-summary|walker-traversal-summary|"
                        + "raw-activity|raw-activity-force|smoothed-activity|active-locus|genotype-likelihoods|"
                        + "assembly-graph|assembly-graph-multi|assembly-graph-summary|"
                        + "assembly-graph-dangling-summary|assembly-graph-non-unique-summary|"
                        + "assembly-haplotype-cigars|assembly-haplotypes|assembly-kbest-paths|"
                        + "assembly-haplotypes-cap|assembly-haplotypes-production|"
                        + "assembly-junction-haplotypes|read-error-correction|"
                        + "assembly-seqgraph-summary|assembly-assemble|pairhmm-likelihoods|pairhmm-native-likelihoods|"
                        + "pairhmm-bq-cap|pairhmm-haplotype-filter|"
                        + "genotyping-aggregate|genotype-format|annotate-core|annotation-manifest|"
                        + "call-region-vcf|call-region-format|ad-annotation-call|filter-poorly-modeled-call|filter-poorly-modeled-call-double|variant-vcf-from-gl-ad|variant-format-from-gl-ad|"
                        + "af-em|subset-alleles-pl|subset-alleles-vc|subset-alleles-integration|"
                        + "gvcf-header|gvcf-writer-blocks|"
                        + "assembly-region-genotype|assembly-region-genotype-subset|"
                        + "gvcf-merge-ref-confidence|reference-confidence-locus|inactive-reference-model|"
                        + "ploidy-resolution|"
                        + "downsample-positional|allele-biased-target-counts|allele-biased-evidence|"
                        + "raw-activity-contam|soft-clip-mean|read-shard-pipeline|read-pre-softclip|read-pre-len|read-pre-mq|read-pre-overlap> ...");
        System.exit(2);
    }

    private static final class AssemblyReadRow {
        final String sequence;
        final byte qual;

        AssemblyReadRow(final String sequence, final byte qual) {
            this.sequence = sequence;
            this.qual = qual;
        }
    }

    private static final class AlignedAssemblyReadRow {
        final byte[] bases;
        final byte[] quals;
        final int start1;

        AlignedAssemblyReadRow(final byte[] bases, final byte[] quals, final int start1) {
            this.bases = bases;
            this.quals = quals;
            this.start1 = start1;
        }
    }

    private static List<AssemblyReadRow> loadAssemblyReadsTsv(final Path path) throws Exception {
        final List<AssemblyReadRow> out = new ArrayList<>();
        try (BufferedReader br = Files.newBufferedReader(path, StandardCharsets.UTF_8)) {
            String line;
            while ((line = br.readLine()) != null) {
                final String t = line.trim();
                if (t.isEmpty() || t.startsWith("#")) {
                    continue;
                }
                final String[] parts = t.split("\\s+");
                if (parts.length < 2) {
                    throw new IllegalArgumentException("reads tsv: missing columns in " + path);
                }
                final int q = Integer.parseInt(parts[1]);
                if (q < 0 || q > 255) {
                    throw new IllegalArgumentException("reads tsv: invalid qual in " + path);
                }
                out.add(new AssemblyReadRow(parts[0], (byte) q));
            }
        }
        return out;
    }

    private static List<AlignedAssemblyReadRow> loadAlignedAssemblyReadsTsv(final Path path)
            throws Exception {
        final List<AlignedAssemblyReadRow> out = new ArrayList<>();
        try (BufferedReader br = Files.newBufferedReader(path, StandardCharsets.UTF_8)) {
            String line;
            while ((line = br.readLine()) != null) {
                final String t = line.trim();
                if (t.isEmpty() || t.startsWith("#")) {
                    continue;
                }
                final String[] parts = t.split("\\s+");
                if (parts.length < 3) {
                    throw new IllegalArgumentException("reads tsv: missing columns in " + path);
                }
                final byte[] bases = parts[0].getBytes(StandardCharsets.US_ASCII);
                final int q = Integer.parseInt(parts[1]);
                if (q < 0 || q > 255) {
                    throw new IllegalArgumentException("reads tsv: invalid qual in " + path);
                }
                final byte[] quals = new byte[bases.length];
                Arrays.fill(quals, (byte) q);
                final int start1 = Integer.parseInt(parts[2]);
                out.add(new AlignedAssemblyReadRow(bases, quals, start1));
            }
        }
        return out;
    }

    /** GATK {@code AbstractReadThreadingGraph.addRead}: high-qual runs only, count = 1. */
    private static void addReadThreadingSequence(
            final AbstractReadThreadingGraph graph,
            final String name,
            final byte[] bases,
            final byte[] quals,
            final int minQual) {
        addReadThreadingSequence(graph, name, bases, quals, minQual, false);
    }

    private static void addReadThreadingSequence(
            final AbstractReadThreadingGraph graph,
            final String name,
            final byte[] bases,
            final byte[] quals,
            final int minQual,
            final boolean isRef) {
        int lastGood = -1;
        for (int end = 0; end <= bases.length; end++) {
            final boolean unusable =
                    end == bases.length
                            || quals[end] < minQual
                            || !isBaseUsableForAssembly(bases[end]);
            if (unusable) {
                if (lastGood != -1) {
                    final int start = lastGood;
                    final int len = end - start;
                    if (len >= graph.getKmerSize()) {
                        graph.addSequence(
                                name + "_" + start + "_" + end,
                                Arrays.copyOfRange(bases, start, end),
                                1,
                                isRef);
                    }
                }
                lastGood = -1;
            } else if (lastGood == -1) {
                lastGood = end;
            }
        }
    }

    private static boolean isBaseUsableForAssembly(final byte base) {
        switch (base) {
            case 'A':
            case 'C':
            case 'G':
            case 'T':
            case 'N':
                return true;
            default:
                return false;
        }
    }

    private static ReadThreadingGraph buildReadThreadingGraph(
            final Path readsPath, final int kmerSize, final int minQual) throws Exception {
        final ReadThreadingGraph graph = new ReadThreadingGraph(kmerSize);
        final List<AssemblyReadRow> rows = loadAssemblyReadsTsv(readsPath);
        int i = 0;
        for (final AssemblyReadRow row : rows) {
            final byte[] bases = row.sequence.getBytes(StandardCharsets.US_ASCII);
            final byte[] quals = new byte[bases.length];
            Arrays.fill(quals, row.qual);
            addReadThreadingSequence(graph, "r" + (i++), bases, quals, minQual);
        }
        graph.buildGraphIfNecessary();
        return graph;
    }

    private static ReadThreadingGraph buildReadThreadingGraphFromRefAndReads(
            final Path refPath, final Path readsPath, final int kmerSize, final int minQual)
            throws Exception {
        final ReadThreadingGraph graph = new ReadThreadingGraph(kmerSize);
        final List<AssemblyReadRow> refRows = loadAssemblyReadsTsv(refPath);
        if (refRows.isEmpty()) {
            throw new IllegalArgumentException("ref tsv: no sequence row in " + refPath);
        }
        final AssemblyReadRow ref = refRows.get(0);
        final byte[] refBases = ref.sequence.getBytes(StandardCharsets.US_ASCII);
        final byte[] refQuals = new byte[refBases.length];
        Arrays.fill(refQuals, ref.qual);
        addReadThreadingSequence(graph, "ref", refBases, refQuals, minQual, true);
        int i = 0;
        for (final AssemblyReadRow row : loadAssemblyReadsTsv(readsPath)) {
            final byte[] bases = row.sequence.getBytes(StandardCharsets.US_ASCII);
            final byte[] quals = new byte[bases.length];
            Arrays.fill(quals, row.qual);
            addReadThreadingSequence(graph, "r" + (i++), bases, quals, minQual, false);
        }
        graph.buildGraphIfNecessary();
        return graph;
    }

    private static SWParameters danglingEndSwParameters() {
        return new SWParameters(25, -50, -110, -6);
    }

    private static boolean isRefSink(final ReadThreadingGraph graph, final MultiDeBruijnVertex v) {
        if (graph.outDegreeOf(v) != 0) {
            return false;
        }
        for (final MultiSampleEdge edge : graph.incomingEdgesOf(v)) {
            if (edge.isRef()) {
                return true;
            }
        }
        return false;
    }

    private static int countDanglingTailCandidates(final ReadThreadingGraph graph) {
        int n = 0;
        for (final MultiDeBruijnVertex v : graph.vertexSet()) {
            if (graph.outDegreeOf(v) == 0 && graph.inDegreeOf(v) > 0 && !isRefSink(graph, v)) {
                n++;
            }
        }
        return n;
    }

    private static String vertexKmer(final MultiDeBruijnVertex v) {
        return v.getSequenceString();
    }

    private static final class EdgeRow implements Comparable<EdgeRow> {
        final String fromKmer;
        final String toKmer;
        final int support;
        final int kmerSize;

        EdgeRow(final String fromKmer, final String toKmer, final int support, final int kmerSize) {
            this.fromKmer = fromKmer;
            this.toKmer = toKmer;
            this.support = support;
            this.kmerSize = kmerSize;
        }

        @Override
        public int compareTo(final EdgeRow o) {
            int c = Integer.compare(kmerSize, o.kmerSize);
            if (c != 0) {
                return c;
            }
            c = fromKmer.compareTo(o.fromKmer);
            if (c != 0) {
                return c;
            }
            return toKmer.compareTo(o.toKmer);
        }
    }

    private static List<EdgeRow> collectEdgeRows(final ReadThreadingGraph graph) {
        final List<EdgeRow> rows = new ArrayList<>();
        for (final MultiSampleEdge edge : graph.edgeSet()) {
            final MultiDeBruijnVertex from = graph.getEdgeSource(edge);
            final MultiDeBruijnVertex to = graph.getEdgeTarget(edge);
            rows.add(
                    new EdgeRow(
                            vertexKmer(from),
                            vertexKmer(to),
                            edge.getPruningMultiplicity(),
                            graph.getKmerSize()));
        }
        Collections.sort(rows);
        return rows;
    }

    private static ChainPruner<MultiDeBruijnVertex, MultiSampleEdge> makeChainPruner(
            final int minPruneFactor, final boolean adaptive) {
        if (adaptive) {
            return new AdaptiveChainPruner<>(
                    0.001,
                    ReadThreadingAssemblerArgumentCollection.DEFAULT_PRUNING_LOG_ODDS_THRESHOLD,
                    ReadThreadingAssemblerArgumentCollection.DEFAULT_PRUNING_SEEDING_LOG_ODDS_THRESHOLD,
                    100);
        }
        return new LowWeightChainPruner<>(minPruneFactor);
    }

    private static String formatSummaryDouble(final double v) {
        if (Double.isInfinite(v) && v < 0) {
            return "-inf";
        }
        if (Double.isInfinite(v) && v > 0) {
            return "inf";
        }
        return String.format(Locale.ROOT, "%.8f", v);
    }

    private static void assemblyGraph(final String[] args) throws Exception {
        if (args.length < 3) {
            usage();
        }
        final Path readsPath = Paths.get(args[0]);
        final int kmerSize = Integer.parseInt(args[1]);
        final int minQual = Integer.parseInt(args[2]);
        final ReadThreadingGraph graph = buildReadThreadingGraph(readsPath, kmerSize, minQual);
        System.out.println("from_kmer\tto_kmer\tsupport");
        for (final EdgeRow row : collectEdgeRows(graph)) {
            System.out.printf("%s\t%s\t%d%n", row.fromKmer, row.toKmer, row.support);
        }
    }

    private static void assemblyGraphMulti(final String[] args) throws Exception {
        if (args.length < 3) {
            usage();
        }
        final Path readsPath = Paths.get(args[0]);
        final String[] kmerParts = args[1].split(",");
        final int minQual = Integer.parseInt(args[2]);
        System.out.println("kmer_size\tfrom_kmer\tto_kmer\tsupport");
        for (final String kmerToken : kmerParts) {
            final int kmerSize = Integer.parseInt(kmerToken.trim());
            final ReadThreadingGraph graph = buildReadThreadingGraph(readsPath, kmerSize, minQual);
            for (final EdgeRow row : collectEdgeRows(graph)) {
                System.out.printf(
                        "%d\t%s\t%s\t%d%n",
                        row.kmerSize, row.fromKmer, row.toKmer, row.support);
            }
        }
    }

    private static void assemblyGraphSummary(final String[] args) throws Exception {
        if (args.length < 5) {
            usage();
        }
        final Path readsPath = Paths.get(args[0]);
        final int kmerSize = Integer.parseInt(args[1]);
        final int minQual = Integer.parseInt(args[2]);
        final int minPrune = Integer.parseInt(args[3]);
        final boolean adaptive = "1".equals(args[4]) || "true".equalsIgnoreCase(args[4]);
        final ReadThreadingGraph graph = buildReadThreadingGraph(readsPath, kmerSize, minQual);
        final int edgesBefore = graph.edgeSet().size();
        final ChainPruner<MultiDeBruijnVertex, MultiSampleEdge> pruner =
                makeChainPruner(minPrune, adaptive);
        pruner.pruneLowWeightChains(graph);
        graph.removeSingletonOrphanVertices();
        final int edgesAfter = graph.edgeSet().size();
        final int edgesPruned = Math.max(0, edgesBefore - edgesAfter);
        long sumSupport = 0;
        int maxSupport = 0;
        for (final MultiSampleEdge edge : graph.edgeSet()) {
            final int m = edge.getPruningMultiplicity();
            sumSupport += m;
            maxSupport = Math.max(maxSupport, m);
        }
        final double log10Max =
                maxSupport <= 0 ? Double.NEGATIVE_INFINITY : Math.log10(Math.max(1, maxSupport));
        final double log10Sum =
                sumSupport <= 0 ? Double.NEGATIVE_INFINITY : Math.log10(sumSupport);
        final double lodThreshold =
                ReadThreadingAssemblerArgumentCollection.DEFAULT_PRUNING_LOG_ODDS_THRESHOLD;
        System.out.println("metric\tvalue");
        System.out.printf("node_count\t%d%n", graph.vertexSet().size());
        System.out.printf("edge_count\t%d%n", edgesAfter);
        System.out.printf("log10_max_edge_support\t%s%n", formatSummaryDouble(log10Max));
        System.out.printf("log10_sum_edge_support\t%s%n", formatSummaryDouble(log10Sum));
        System.out.printf("pruning_lod_threshold_ln\t%s%n", formatSummaryDouble(lodThreshold));
        System.out.printf("adaptive_pruning\t%s%n", adaptive ? "true" : "false");
        System.out.printf("min_prune_factor\t%d%n", minPrune);
        System.out.printf("edges_pruned\t%d%n", edgesPruned);
    }

    private static void assemblyGraphDanglingSummary(final String[] args) throws Exception {
        if (args.length < 7) {
            usage();
        }
        final Path refPath = Paths.get(args[0]);
        final Path readsPath = Paths.get(args[1]);
        final int kmerSize = Integer.parseInt(args[2]);
        final int minQual = Integer.parseInt(args[3]);
        final int minPrune = Integer.parseInt(args[4]);
        final int minDangling = Integer.parseInt(args[5]);
        final boolean recoverHeads = "1".equals(args[6]) || "true".equalsIgnoreCase(args[6]);
        final boolean recoverAll = false;
        final ReadThreadingGraph graph =
                buildReadThreadingGraphFromRefAndReads(refPath, readsPath, kmerSize, minQual);
        final ChainPruner<MultiDeBruijnVertex, MultiSampleEdge> pruner =
                makeChainPruner(minPrune, false);
        pruner.pruneLowWeightChains(graph);
        graph.removeSingletonOrphanVertices();
        final int edgesBefore = graph.edgeSet().size();
        final int tailsAttempted = countDanglingTailCandidates(graph);
        final SmithWatermanAligner aligner =
                SmithWatermanAligner.getAligner(SmithWatermanAligner.Implementation.JAVA);
        final SWParameters swParams = danglingEndSwParameters();
        graph.recoverDanglingTails(
                minPrune, minDangling, recoverAll, aligner, swParams);
        graph.removeSingletonOrphanVertices();
        int headsAttempted = 0;
        int headsRecovered = 0;
        if (recoverHeads) {
            for (final MultiDeBruijnVertex v : graph.vertexSet()) {
                if (graph.inDegreeOf(v) == 0 && graph.outDegreeOf(v) > 0) {
                    headsAttempted++;
                }
            }
            final int edgesBeforeHeads = graph.edgeSet().size();
            graph.recoverDanglingHeads(
                    minPrune, minDangling, recoverAll, aligner, swParams);
            graph.removeSingletonOrphanVertices();
            headsRecovered = Math.max(0, graph.edgeSet().size() - edgesBeforeHeads);
        }
        final int edgesAfter = graph.edgeSet().size();
        final int edgesMerged = Math.max(0, edgesAfter - edgesBefore);
        final int tailsRecovered = Math.max(0, edgesMerged - headsRecovered);
        System.out.println("metric\tvalue");
        System.out.printf("edges_before\t%d%n", edgesBefore);
        System.out.printf("edges_after\t%d%n", edgesAfter);
        System.out.printf("tails_attempted\t%d%n", tailsAttempted);
        System.out.printf("tails_recovered\t%d%n", tailsRecovered);
        System.out.printf("heads_attempted\t%d%n", headsAttempted);
        System.out.printf("heads_recovered\t%d%n", headsRecovered);
        System.out.printf("edges_merged\t%d%n", edgesMerged);
    }

    private static Field findThreadingField(final ReadThreadingGraph graph, final String name)
            throws Exception {
        for (Class<?> c = graph.getClass(); c != null; c = c.getSuperclass()) {
            try {
                final Field f = c.getDeclaredField(name);
                f.setAccessible(true);
                return f;
            } catch (final NoSuchFieldException ignored) {
                // try superclass (fields live on AbstractReadThreadingGraph)
            }
        }
        throw new NoSuchFieldException(name);
    }

    @SuppressWarnings("unchecked")
    private static int readThreadingUniqueKmerCount(final ReadThreadingGraph graph) throws Exception {
        // GATK 4.4: `kmerToVertexMap` on AbstractReadThreadingGraph (was `uniqueKmers` in older decomp).
        return ((Map<?, ?>) findThreadingField(graph, "kmerToVertexMap").get(graph)).size();
    }

    @SuppressWarnings("unchecked")
    private static int readThreadingNonUniqueKmerCount(final ReadThreadingGraph graph) throws Exception {
        final Object set = findThreadingField(graph, "nonUniqueKmers").get(graph);
        return set == null ? 0 : ((Set<?>) set).size();
    }

    private static int maxKmerMultiplicity(final ReadThreadingGraph graph) {
        final Map<String, Integer> mult = new HashMap<>();
        for (final MultiDeBruijnVertex v : graph.vertexSet()) {
            final String kmer = new String(v.getSequence(), StandardCharsets.US_ASCII);
            mult.merge(kmer, 1, Integer::sum);
        }
        int max = 0;
        for (final int c : mult.values()) {
            max = Math.max(max, c);
        }
        return max;
    }

    private static void assemblyGraphNonUniqueSummary(final String[] args) throws Exception {
        if (args.length < 4) {
            usage();
        }
        final String refToken = args[0];
        final Path readsPath = Paths.get(args[1]);
        final int kmerSize = Integer.parseInt(args[2]);
        final int minQual = Integer.parseInt(args[3]);
        final ReadThreadingGraph graph;
        if ("-".equals(refToken)) {
            graph = buildReadThreadingGraph(readsPath, kmerSize, minQual);
        } else {
            graph = buildReadThreadingGraphFromRefAndReads(Paths.get(refToken), readsPath, kmerSize, minQual);
        }
        System.out.println("metric\tvalue");
        System.out.printf("node_count\t%d%n", graph.vertexSet().size());
        System.out.printf("edge_count\t%d%n", graph.edgeSet().size());
        System.out.printf("unique_kmer_count\t%d%n", readThreadingUniqueKmerCount(graph));
        System.out.printf("non_unique_kmer_count\t%d%n", readThreadingNonUniqueKmerCount(graph));
        System.out.printf(
                "is_low_complexity\t%s%n", graph.isLowQualityGraph() ? "true" : "false");
        System.out.printf("max_kmer_multiplicity\t%d%n", maxKmerMultiplicity(graph));
    }

  private static String formatScore(final double v) {
        if (v == 0.0) {
            return "0";
        }
        return String.format(Locale.ROOT, "%.8f", v);
    }

    private static ReadThreadingGraph prepareReadThreadingGraphForHaplotypeDump(
            final Path refPath,
            final Path readsPath,
            final int kmerSize,
            final int minQual,
            final int minPrune,
            final int minDangling,
            final boolean recoverHeads)
            throws Exception {
        final ReadThreadingGraph graph =
                buildReadThreadingGraphFromRefAndReads(refPath, readsPath, kmerSize, minQual);
        final ChainPruner<MultiDeBruijnVertex, MultiSampleEdge> pruner =
                makeChainPruner(minPrune, false);
        pruner.pruneLowWeightChains(graph);
        graph.removeSingletonOrphanVertices();
        final SmithWatermanAligner aligner =
                SmithWatermanAligner.getAligner(SmithWatermanAligner.Implementation.JAVA);
        final SWParameters swParams = danglingEndSwParameters();
        graph.recoverDanglingTails(minPrune, minDangling, false, aligner, swParams);
        if (recoverHeads) {
            graph.recoverDanglingHeads(minPrune, minDangling, false, aligner, swParams);
        }
        graph.removeSingletonOrphanVertices();
        graph.removePathsNotConnectedToRef();
        return graph;
    }

    private static List<KBestHaplotype<MultiDeBruijnVertex, MultiSampleEdge>> rankKbestPaths(
            final ReadThreadingGraph graph, final int maxHaplotypes) {
        final MultiDeBruijnVertex source = graph.getReferenceSourceVertex();
        final MultiDeBruijnVertex sink = graph.getReferenceSinkVertex();
        if (source == null || sink == null) {
            return Collections.emptyList();
        }
        final List<KBestHaplotype<MultiDeBruijnVertex, MultiSampleEdge>> ranked =
                new ArrayList<>(
                        new GraphBasedKBestHaplotypeFinder<>(graph, source, sink)
                                .findBestHaplotypes(maxHaplotypes));
        ranked.sort(
                Comparator.comparingDouble((KBestHaplotype<MultiDeBruijnVertex, MultiSampleEdge> p) -> p.score())
                        .reversed()
                        .thenComparing(
                                p -> new String(p.getBases()),
                                Comparator.reverseOrder()));
        return ranked;
    }

    private static void assemblyHaplotypes(final String[] args) throws Exception {
        if (args.length < 7) {
            usage();
        }
        emitAssemblyHaplotypesWithCap(args, 128, false);
    }

    private static void assemblyHaplotypesCap(final String[] args) throws Exception {
        if (args.length < 8) {
            usage();
        }
        emitAssemblyHaplotypesWithCap(args, Integer.parseInt(args[7]), false);
    }

    private static void emitAssemblyHaplotypesWithCap(
            final String[] args, final int maxHaplotypes, final boolean productionRefTagging)
            throws Exception {
        final Path refPath = Paths.get(args[0]);
        final Path readsPath = Paths.get(args[1]);
        final int kmerSize = Integer.parseInt(args[2]);
        final int minQual = Integer.parseInt(args[3]);
        final int minPrune = Integer.parseInt(args[4]);
        final int minDangling = Integer.parseInt(args[5]);
        final boolean recoverHeads = "1".equals(args[6]) || "true".equalsIgnoreCase(args[6]);
        final ReadThreadingGraph graph =
                prepareReadThreadingGraphForHaplotypeDump(
                        refPath, readsPath, kmerSize, minQual, minPrune, minDangling, recoverHeads);
        final byte[] refBases =
                loadAssemblyReadsTsv(refPath).get(0).sequence.getBytes(StandardCharsets.US_ASCII);
        final SmithWatermanAligner aligner =
                SmithWatermanAligner.getAligner(SmithWatermanAligner.Implementation.JAVA);
        final SWParameters hapSw =
                org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanAlignmentConstants
                        .NEW_SW_PARAMETERS;
        final List<KBestHaplotype<MultiDeBruijnVertex, MultiSampleEdge>> ranked =
                rankKbestPaths(graph, maxHaplotypes);
        System.out.println("rank\tsequence\tscore\tis_reference\tcigar");
        int rank = 0;
        boolean wroteRef = false;
        for (final KBestHaplotype<MultiDeBruijnVertex, MultiSampleEdge> path : ranked) {
            final Haplotype h = path.haplotype();
            final htsjdk.samtools.Cigar cigar =
                    CigarUtils.calculateCigar(
                            refBases,
                            h.getBases(),
                            aligner,
                            hapSw,
                            org.broadinstitute.gatk.nativebindings.smithwaterman.SWOverhangStrategy
                                    .SOFTCLIP);
            if (cigar == null) {
                continue;
            }
            if (cigar.getReferenceLength() != refBases.length
                    && refBases.length >= 30
                    && cigar.getReferenceLength() < 30) {
                continue;
            }
            final boolean isRefPath =
                    productionRefTagging
                            ? Arrays.equals(h.getBases(), refBases)
                            : path.isReference();
            if (Arrays.equals(h.getBases(), refBases)) {
                wroteRef = true;
            }
            System.out.printf(
                    "%d\t%s\t%s\t%s\t%s%n",
                    rank++,
                    new String(h.getBases()),
                    formatScore(path.score()),
                    isRefPath,
                    cigar);
        }
        if (!wroteRef) {
            final htsjdk.samtools.Cigar refCigar =
                    CigarUtils.calculateCigar(
                            refBases,
                            refBases,
                            aligner,
                            hapSw,
                            org.broadinstitute.gatk.nativebindings.smithwaterman.SWOverhangStrategy
                                    .SOFTCLIP);
            System.out.printf(
                    "%d\t%s\t0\ttrue\t%s%n",
                    rank,
                    new String(refBases),
                    refCigar == null ? "" : refCigar.toString());
        }
    }

    private static void assemblyHaplotypesProduction(final String[] args) throws Exception {
        if (args.length < 7) {
            usage();
        }
        emitAssemblyHaplotypesWithCap(args, 128, true);
    }

    private static void assemblyKbestPaths(final String[] args) throws Exception {
        if (args.length < 8) {
            usage();
        }
        final Path refPath = Paths.get(args[0]);
        final Path readsPath = Paths.get(args[1]);
        final int kmerSize = Integer.parseInt(args[2]);
        final int minQual = Integer.parseInt(args[3]);
        final int minPrune = Integer.parseInt(args[4]);
        final int minDangling = Integer.parseInt(args[5]);
        final boolean recoverHeads = "1".equals(args[6]) || "true".equalsIgnoreCase(args[6]);
        final int maxHaplotypes = Integer.parseInt(args[7]);
        final ReadThreadingGraph graph =
                prepareReadThreadingGraphForHaplotypeDump(
                        refPath, readsPath, kmerSize, minQual, minPrune, minDangling, recoverHeads);
        final List<KBestHaplotype<MultiDeBruijnVertex, MultiSampleEdge>> ranked =
                rankKbestPaths(graph, maxHaplotypes);
        System.out.println("rank\tsequence\tscore\tis_reference");
        int rank = 0;
        for (final KBestHaplotype<MultiDeBruijnVertex, MultiSampleEdge> path : ranked) {
            System.out.printf(
                    "%d\t%s\t%s\t%s%n",
                    rank++,
                    new String(path.getBases()),
                    formatScore(path.score()),
                    path.isReference());
        }
    }

    private static JunctionTreeLinkedDeBruijnGraph buildJunctionGraphFromRefAndReads(
            final Path refPath, final Path readsPath, final int kmerSize, final int minQual)
            throws Exception {
        final JunctionTreeLinkedDeBruijnGraph graph = new JunctionTreeLinkedDeBruijnGraph(kmerSize);
        final List<AssemblyReadRow> refRows = loadAssemblyReadsTsv(refPath);
        if (refRows.isEmpty()) {
            throw new IllegalArgumentException("ref tsv: no sequence row in " + refPath);
        }
        final AssemblyReadRow ref = refRows.get(0);
        final byte[] refBases = ref.sequence.getBytes(StandardCharsets.US_ASCII);
        final byte[] refQuals = new byte[refBases.length];
        Arrays.fill(refQuals, ref.qual);
        addReadThreadingSequence(graph, "ref", refBases, refQuals, minQual, true);
        int i = 0;
        for (final AssemblyReadRow row : loadAssemblyReadsTsv(readsPath)) {
            final byte[] bases = row.sequence.getBytes(StandardCharsets.US_ASCII);
            final byte[] quals = new byte[bases.length];
            Arrays.fill(quals, row.qual);
            addReadThreadingSequence(graph, "r" + (i++), bases, quals, minQual, false);
        }
        graph.buildGraphIfNecessary();
        graph.generateJunctionTrees();
        return graph;
    }

    private static void assemblyJunctionHaplotypes(final String[] args) throws Exception {
        if (args.length < 6) {
            usage();
        }
        final Path refPath = Paths.get(args[0]);
        final Path readsPath = Paths.get(args[1]);
        final int kmerSize = Integer.parseInt(args[2]);
        final int minQual = Integer.parseInt(args[3]);
        final boolean recoverEdges =
                "1".equals(args[4]) || "true".equalsIgnoreCase(args[4]);
        final int maxHaplotypes = Integer.parseInt(args[5]);
        final JunctionTreeLinkedDeBruijnGraph graph =
                buildJunctionGraphFromRefAndReads(refPath, readsPath, kmerSize, minQual);
        final MultiDeBruijnVertex source = graph.getReferenceSourceVertex();
        final MultiDeBruijnVertex sink = graph.getReferenceSinkVertex();
        if (source == null || sink == null) {
            System.out.println("rank\tsequence\tscore\tis_reference");
            return;
        }
        final List<KBestHaplotype<MultiDeBruijnVertex, MultiSampleEdge>> ranked =
                new ArrayList<>(
                        new JunctionTreeKBestHaplotypeFinder<>(
                                        graph,
                                        source,
                                        sink,
                                        JunctionTreeKBestHaplotypeFinder
                                                .DEFAULT_OUTGOING_JT_EVIDENCE_THRESHOLD_TO_BELEIVE,
                                        recoverEdges)
                                .setJunctionTreeEvidenceWeightThreshold(1)
                                .findBestHaplotypes(maxHaplotypes));
        ranked.sort(
                Comparator.comparingDouble((KBestHaplotype<MultiDeBruijnVertex, MultiSampleEdge> p) -> p.score())
                        .reversed()
                        .thenComparing(
                                p -> new String(p.getBases()),
                                Comparator.reverseOrder()));
        System.out.println("rank\tsequence\tscore\tis_reference");
        int rank = 0;
        for (final KBestHaplotype<MultiDeBruijnVertex, MultiSampleEdge> path : ranked) {
            System.out.printf(
                    "%d\t%s\t%s\t%s%n",
                    rank++,
                    new String(path.getBases()),
                    formatScore(path.score()),
                    path.isReference());
        }
    }

    private static void readErrorCorrection(final String[] args) throws Exception {
        if (args.length < 2) {
            usage();
        }
        final Path readsPath = Paths.get(args[0]);
        final double logOdds = Double.parseDouble(args[1]);
        final List<AlignedAssemblyReadRow> rows = loadAlignedAssemblyReadsTsv(readsPath);
        final SAMReadGroupRecord rg = new SAMReadGroupRecord("rg");
        rg.setSample("sample");
        rg.setPlatform("ILLUMINA");
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeaderWithReadGroup(rg);
        final String contig =
                header.getSequenceDictionary().getSequences().get(0).getSequenceName();
        final List<GATKRead> reads = new ArrayList<>();
        int idx = 0;
        for (final AlignedAssemblyReadRow row : rows) {
            final GATKRead r =
                    ArtificialReadUtils.createArtificialRead(
                            header,
                            "r" + (idx++),
                            contig,
                            row.start1,
                            row.bases,
                            row.quals);
            r.setReadGroup(rg.getReadGroupId());
            reads.add(r);
        }
        final PileupReadErrorCorrector corrector = new PileupReadErrorCorrector(logOdds, header);
        final List<GATKRead> corrected = corrector.correctReads(reads);
        System.out.println("read_index\tsequence\tmean_qual");
        for (int i = 0; i < corrected.size(); i++) {
            final GATKRead r = corrected.get(i);
            final byte[] quals = r.getBaseQualities();
            int mean = 0;
            if (quals.length > 0) {
                int sum = 0;
                for (final byte q : quals) {
                    sum += q & 0xff;
                }
                mean = sum / quals.length;
            }
            System.out.printf(
                    "%d\t%s\t%d%n", i, new String(r.getBases(), StandardCharsets.US_ASCII), mean);
        }
    }

    private static void printSeqGraphMetric(final String metric, final String value) {
        System.out.printf("%s\t%s%n", metric, value);
    }

    private static void rustParitySimplifySeqGraph(final SeqGraph seqGraph) {
        seqGraph.zipLinearChains();
        for (int i = 0; i < 100; i++) {
            if (!seqGraph.zipLinearChains()) {
                break;
            }
        }
    }

    /** Rust {@code SeqGraph::cleanup_seq_graph} aligned with production {@code simplifyGraph} (GAP-E-03). */
    private static String cleanupSeqGraphRustParity(final SeqGraph seqGraph) {
        seqGraph.zipLinearChains();
        seqGraph.removeSingletonOrphanVertices();
        seqGraph.removeVerticesNotConnectedToRefRegardlessOfEdgeDirection();
        seqGraph.simplifyGraph();
        if (seqGraph.getReferenceSourceVertex() == null || seqGraph.getReferenceSinkVertex() == null) {
            return "just_assembled_reference";
        }
        seqGraph.removePathsNotConnectedToRef();
        seqGraph.simplifyGraph();
        if (seqGraph.vertexSet().size() == 1) {
            final SeqVertex complete = seqGraph.vertexSet().iterator().next();
            final SeqVertex dummy = new SeqVertex("");
            seqGraph.addVertex(dummy);
            seqGraph.addEdge(complete, dummy, new BaseEdge(true, 0));
        }
        return "assembled_some_variation";
    }

    private static void assemblySeqGraphSummary(final String[] args) throws Exception {
        if (args.length < 7) {
            usage();
        }
        final Path refPath = Paths.get(args[0]);
        final Path readsPath = Paths.get(args[1]);
        final int kmerSize = Integer.parseInt(args[2]);
        final int minQual = Integer.parseInt(args[3]);
        final int minPrune = Integer.parseInt(args[4]);
        final int minDangling = Integer.parseInt(args[5]);
        final boolean recoverHeads = "1".equals(args[6]) || "true".equalsIgnoreCase(args[6]);
        System.out.println("metric\tvalue");
        final ReadThreadingGraph graph =
                prepareReadThreadingGraphForHaplotypeDump(
                        refPath, readsPath, kmerSize, minQual, minPrune, minDangling, recoverHeads);
        if (graph.getReferenceSourceVertex() == null || graph.getReferenceSinkVertex() == null) {
            printSeqGraphMetric("status", "no_graph");
            printSeqGraphMetric("node_count", "0");
            printSeqGraphMetric("edge_count", "0");
            printSeqGraphMetric("ref_path_len", "0");
            printSeqGraphMetric("ref_path_sequence", "");
            return;
        }
        final SeqGraph seqGraph = graph.toSequenceGraph();
        seqGraph.cleanNonRefPaths();
        final String status = cleanupSeqGraphRustParity(seqGraph);
        final byte[] refPathBytes;
        if (seqGraph.getReferenceSourceVertex() != null && seqGraph.getReferenceSinkVertex() != null) {
            refPathBytes =
                    seqGraph.getReferenceBytes(
                            seqGraph.getReferenceSourceVertex(),
                            seqGraph.getReferenceSinkVertex(),
                            true,
                            true);
        } else {
            refPathBytes = new byte[0];
        }
        printSeqGraphMetric("status", status);
        printSeqGraphMetric("node_count", Integer.toString(seqGraph.vertexSet().size()));
        printSeqGraphMetric("edge_count", Integer.toString(seqGraph.edgeSet().size()));
        printSeqGraphMetric("ref_path_len", Integer.toString(refPathBytes.length));
        printSeqGraphMetric("ref_path_sequence", new String(refPathBytes, StandardCharsets.US_ASCII));
    }

    private static final class AssembleAttempt {
        final String status;
        final int kmerSize;
        final List<KBestHaplotype<MultiDeBruijnVertex, MultiSampleEdge>> ranked;
        final byte[] refBases;

        AssembleAttempt(
                final String status,
                final int kmerSize,
                final List<KBestHaplotype<MultiDeBruijnVertex, MultiSampleEdge>> ranked,
                final byte[] refBases) {
            this.status = status;
            this.kmerSize = kmerSize;
            this.ranked = ranked;
            this.refBases = refBases;
        }
    }

    private static AssembleAttempt tryAssembleKmer(
            final Path refPath,
            final Path readsPath,
            final int kmerSize,
            final int minQual,
            final int minPrune,
            final int minDangling,
            final boolean recoverHeads)
            throws Exception {
        final byte[] refBases =
                loadAssemblyReadsTsv(refPath).get(0).sequence.getBytes(StandardCharsets.US_ASCII);
        if (refBases.length < kmerSize) {
            return new AssembleAttempt("failed", kmerSize, Collections.emptyList(), refBases);
        }
        final ReadThreadingGraph graph =
                prepareReadThreadingGraphForHaplotypeDump(
                        refPath, readsPath, kmerSize, minQual, minPrune, minDangling, recoverHeads);
        if (graph.getReferenceSourceVertex() == null || graph.getReferenceSinkVertex() == null) {
            return new AssembleAttempt(
                    "just_assembled_reference", kmerSize, Collections.emptyList(), refBases);
        }
        final List<KBestHaplotype<MultiDeBruijnVertex, MultiSampleEdge>> ranked =
                rankKbestPaths(graph, 128);
        final String status =
                ranked.size() <= 1 ? "just_assembled_reference" : "assembled_some_variation";
        return new AssembleAttempt(status, kmerSize, ranked, refBases);
    }

    private static void emitAssemblyAssembleResult(final AssembleAttempt attempt) throws Exception {
        System.out.println("status\t" + attempt.status);
        System.out.println("kmer_size\t" + attempt.kmerSize);
        final SmithWatermanAligner aligner =
                SmithWatermanAligner.getAligner(SmithWatermanAligner.Implementation.JAVA);
        final SWParameters hapSw =
                org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanAlignmentConstants
                        .NEW_SW_PARAMETERS;
        System.out.println("rank\tsequence\tscore\tis_reference\tcigar");
        int rank = 0;
        boolean wroteRef = false;
        for (final KBestHaplotype<MultiDeBruijnVertex, MultiSampleEdge> path : attempt.ranked) {
            final Haplotype h = path.haplotype();
            final htsjdk.samtools.Cigar cigar =
                    CigarUtils.calculateCigar(
                            attempt.refBases,
                            h.getBases(),
                            aligner,
                            hapSw,
                            org.broadinstitute.gatk.nativebindings.smithwaterman.SWOverhangStrategy
                                    .SOFTCLIP);
            if (cigar == null) {
                continue;
            }
            if (cigar.getReferenceLength() != attempt.refBases.length
                    && attempt.refBases.length >= 30
                    && cigar.getReferenceLength() < 30) {
                continue;
            }
            if (Arrays.equals(h.getBases(), attempt.refBases)) {
                wroteRef = true;
            }
            System.out.printf(
                    "%d\t%s\t%s\t%s\t%s%n",
                    rank++,
                    new String(h.getBases(), StandardCharsets.US_ASCII),
                    formatScore(path.score()),
                    path.isReference(),
                    cigar);
        }
        if (!wroteRef) {
            final htsjdk.samtools.Cigar refCigar =
                    CigarUtils.calculateCigar(
                            attempt.refBases,
                            attempt.refBases,
                            aligner,
                            hapSw,
                            org.broadinstitute.gatk.nativebindings.smithwaterman.SWOverhangStrategy
                                    .SOFTCLIP);
            System.out.printf(
                    "%d\t%s\t0\ttrue\t%s%n",
                    rank,
                    new String(attempt.refBases, StandardCharsets.US_ASCII),
                    refCigar == null ? "" : refCigar.toString());
        }
    }

    private static void assemblyAssemble(final String[] args) throws Exception {
        if (args.length < 2) {
            usage();
        }
        final Path refPath = Paths.get(args[0]);
        final Path readsPath = Paths.get(args[1]);
        final int minQual = 10;
        final int minPrune = 2;
        final int minDangling = 4;
        final boolean recoverHeads = true;
        final int[] kmers = {10, 25};
        AssembleAttempt lastFail = null;
        for (int idx = 0; idx < kmers.length; idx++) {
            final AssembleAttempt attempt =
                    tryAssembleKmer(
                            refPath, readsPath, kmers[idx], minQual, minPrune, minDangling, recoverHeads);
            if (!"failed".equals(attempt.status)) {
                emitAssemblyAssembleResult(attempt);
                return;
            }
            lastFail = attempt;
        }
        emitAssemblyAssembleResult(
                lastFail != null
                        ? lastFail
                        : new AssembleAttempt(
                                "failed",
                                kmers[0],
                                Collections.emptyList(),
                                loadAssemblyReadsTsv(refPath)
                                        .get(0)
                                        .sequence
                                        .getBytes(StandardCharsets.US_ASCII)));
    }

    private static void assemblyHaplotypeCigars(final String[] args) throws Exception {
        if (args.length < 2) {
            usage();
        }
        final Path refPath = Paths.get(args[0]);
        final Path hapsPath = Paths.get(args[1]);
        final List<AssemblyReadRow> refRows = loadAssemblyReadsTsv(refPath);
        if (refRows.isEmpty()) {
            throw new IllegalArgumentException("ref tsv: no sequence row in " + refPath);
        }
        final byte[] refBases = refRows.get(0).sequence.getBytes(StandardCharsets.US_ASCII);
        final SmithWatermanAligner aligner =
                SmithWatermanAligner.getAligner(SmithWatermanAligner.Implementation.JAVA);
        final SWParameters swParams =
                org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanAlignmentConstants
                        .NEW_SW_PARAMETERS;
        System.out.println("haplotype_idx\tsequence\tcigar");
        int idx = 0;
        for (final AssemblyReadRow row : loadAssemblyReadsTsv(hapsPath)) {
            final byte[] hapBases = row.sequence.getBytes(StandardCharsets.US_ASCII);
            final htsjdk.samtools.Cigar cigar =
                    CigarUtils.calculateCigar(
                            refBases,
                            hapBases,
                            aligner,
                            swParams,
                            org.broadinstitute.gatk.nativebindings.smithwaterman.SWOverhangStrategy
                                    .SOFTCLIP);
            final String cigarStr = cigar == null ? "" : cigar.toString();
            System.out.printf("%d\t%s\t%s%n", idx++, row.sequence, cigarStr);
        }
    }

    private static int parsePadding(String token) {
        if (token == null || token.isEmpty() || "-".equals(token)) {
            return DEFAULT_PADDING;
        }
        return Integer.parseInt(token);
    }

    /** BAM header sample names for {@link AssemblyBasedCallerUtils#assembleReads} (not hard-coded {@code s1}). */
    private static SampleList sampleListFromHeader(final SAMFileHeader header) {
        final java.util.Set<String> samples = ReadUtils.getSamplesFromHeader(header);
        if (samples.isEmpty()) {
            return SampleList.singletonSampleList("unknown");
        }
        if (samples.size() == 1) {
            return SampleList.singletonSampleList(samples.iterator().next());
        }
        return new IndexedSampleList(samples);
    }

    private static List<String> splitIntervals(String intervalCli) {
        final List<String> out = new ArrayList<>();
        for (final String part : intervalCli.split(";")) {
            if (!part.trim().isEmpty()) {
                out.add(part.trim());
            }
        }
        return out;
    }

    private static List<SimpleInterval> parseIntervals(SAMSequenceDictionary dict, String intervalCli) {
        final GenomeLocParser parser = new GenomeLocParser(dict);
        final List<GenomeLoc> locs =
                IntervalUtils.parseIntervalArguments(parser, splitIntervals(intervalCli));
        final List<SimpleInterval> intervals = new ArrayList<>(locs.size());
        for (final GenomeLoc gl : locs) {
            intervals.add(new SimpleInterval(gl));
        }
        return intervals;
    }

    @SuppressWarnings("unchecked")
    private static List<List<Locatable>> groupByContig(final List<SimpleInterval> intervals) {
        return (List<List<Locatable>>) (List<?>) IntervalUtils.groupIntervalsByContig(intervals);
    }

    private static double[] cappedGenotypeLikelihoods(final RefVsAnyResult refVsAny)
            throws Exception {
        final Method m =
                RefVsAnyResult.class.getDeclaredMethod(
                        "getGenotypeLikelihoodsCappedByHomRefLikelihood");
        m.setAccessible(true);
        return (double[]) m.invoke(refVsAny);
    }

    private static String formatProb(double v) {
        if (v == 0.0) {
            return "0";
        }
        return String.format(Locale.ROOT, "%.8f", v);
    }

    private static String formatKind(ActivityProfileState.Type t) {
        if (t == ActivityProfileState.Type.HIGH_QUALITY_SOFT_CLIPS) {
            return "hq_soft_clips";
        }
        return "none";
    }

    private static FeatureContext featureContextForLocus(final HcContext ctx, final SimpleInterval pileupInterval) {
        if (ctx.forceAllelesInput == null) {
            return new FeatureContext((FeatureManager) null, pileupInterval);
        }
        final Map<FeatureInput<? extends Feature>, Class<? extends Feature>> m =
                Collections.singletonMap(ctx.forceAllelesInput, VariantContext.class);
        return FeatureContext.createFeatureContextForTesting(
                m,
                "HcFullParityGateDump",
                pileupInterval,
                0,
                0,
                0,
                Paths.get(ctx.refPathForFeatures));
    }

    private static final class HcContext implements AutoCloseable {
        final String refPathForFeatures;
        final SAMFileHeader header;
        final ReadsPathDataSource readsSource;
        final ReferenceDataSource reference;
        final HaplotypeCallerEngine engine;
        final AssemblyRegionArgumentCollection asmArgs;
        final List<ReadFilter> readFilters;
        final FeatureInput<VariantContext> forceAllelesInput;

        HcContext(final String refPath, final String bamPath, final int padding) throws Exception {
            this(refPath, bamPath, padding, null, false);
        }

        HcContext(
                final String refPath,
                final String bamPath,
                final int padding,
                final String forceAllelesVcfOrNull)
                throws Exception {
            this(refPath, bamPath, padding, forceAllelesVcfOrNull, false);
        }

        HcContext(
                final String refPath,
                final String bamPath,
                final int padding,
                final String forceAllelesVcfOrNull,
                final boolean nativePairHmmUseDoublePrecision)
                throws Exception {
            refPathForFeatures = refPath;
            Utils.resetRandomGenerator();
            final Path bam = Paths.get(bamPath);
            readsSource = new ReadsPathDataSource(bam);
            header = readsSource.getHeader();
            final Path ref = Paths.get(refPath);
            final CachingIndexedFastaSequenceFile refReader =
                    new CachingIndexedFastaSequenceFile(ref);
            reference = ReferenceDataSource.of(ref);
            asmArgs = new AssemblyRegionArgumentCollection();
            asmArgs.assemblyRegionPadding = padding;
            final HaplotypeCallerArgumentCollection hcArgs = new HaplotypeCallerArgumentCollection();
            FeatureInput<VariantContext> forceIn = null;
            if (forceAllelesVcfOrNull != null) {
                forceIn = new FeatureInput<>(new GATKPath(forceAllelesVcfOrNull));
                hcArgs.alleles = forceIn;
            }
            forceAllelesInput = forceIn;
            if (nativePairHmmUseDoublePrecision) {
                final Field useDbl =
                        PairHMMNativeArgumentCollection.class.getDeclaredField(
                                "useDoublePrecision");
                useDbl.setAccessible(true);
                useDbl.setBoolean(hcArgs.likelihoodArgs.pairHMMNativeArgs, true);
            }
            final VariantAnnotatorEngine annotationEngine =
                    new VariantAnnotatorEngine(
                            Collections.emptyList(),
                            null,
                            Collections.emptyList(),
                            false,
                            false);
            engine =
                    new HaplotypeCallerEngine(
                            hcArgs,
                            asmArgs,
                            false,
                            false,
                            header,
                            refReader,
                            annotationEngine);
            readFilters = HaplotypeCallerEngine.makeStandardHCReadFilters();
            for (final ReadFilter f : readFilters) {
                f.setHeader(header);
            }
        }

        ReadFilter combinedFilter() {
            return ReadFilter.fromList(readFilters, header);
        }

        @Override
        public void close() {
            if (engine != null) {
                engine.shutdown();
            }
            if (readsSource != null) {
                readsSource.close();
            }
        }
    }

    /** Same shard wiring as {@link org.broadinstitute.hellbender.engine.AssemblyRegionWalker#traverse()}. */
    private static void configureHcProductionReadShard(
            final MultiIntervalLocalReadShard shard, final HcContext ctx) {
        shard.setPreReadFilterTransformer(HaplotypeCallerEngine.makeStandardHCReadTransformer());
        shard.setReadFilter(ctx.combinedFilter());
        shard.setPostReadFilterTransformer(ReadTransformer.identity());
        final ReadsDownsampler downsampler = createHcDownsampler(ctx);
        if (downsampler != null) {
            shard.setDownsampler(downsampler);
        }
    }

    private static ReadsDownsampler createHcDownsampler(final HcContext ctx) {
        if (ctx.asmArgs.maxReadsPerAlignmentStart <= 0) {
            return null;
        }
        return new PositionalDownsampler(
                ctx.asmArgs.maxReadsPerAlignmentStart, ctx.header, false);
    }

    /** Rust {@code ReadFilterParams::default()}-style HC gates (Phase D.2 allele / D.3 soft-clip dumps). */
    private static void configureSoftClipReadShard(
            final MultiIntervalLocalReadShard shard, final HcContext ctx) {
        shard.setPreReadFilterTransformer(ReadTransformer.identity());
        final List<ReadFilter> minimal =
                Arrays.asList(
                        new WellformedReadFilter(ctx.header),
                        ReadFilterLibrary.MAPPING_QUALITY_AVAILABLE,
                        new MappingQualityReadFilter(20),
                        ReadFilterLibrary.NOT_DUPLICATE,
                        ReadFilterLibrary.NOT_SECONDARY_ALIGNMENT,
                        ReadFilterLibrary.NOT_SUPPLEMENTARY_ALIGNMENT,
                        new ReadFilterLibrary.MappedReadFilter());
        for (final ReadFilter f : minimal) {
            f.setHeader(ctx.header);
        }
        shard.setReadFilter(ReadFilter.fromList(minimal, ctx.header));
        shard.setPostReadFilterTransformer(ReadTransformer.identity());
        shard.setDownsampler(null);
    }

    private static Iterator<AlignmentContext> makeLocusIteratorSoftClip(
            final MultiIntervalLocalReadShard shard, final HcContext ctx) {
        configureSoftClipReadShard(shard, ctx);
        final SAMFileHeader header = ctx.header;
        final ReadCachingIterator cache = new ReadCachingIterator(shard.iterator());
        final LocusIteratorByState libs =
                new LocusIteratorByState(
                        cache,
                        DownsamplingMethod.NONE,
                        ReadUtils.getSamplesFromHeader(header),
                        header,
                        true);
        final IntervalLocusIterator intervalIt =
                new IntervalLocusIterator(shard.getIntervals().iterator());
        return new IntervalAlignmentContextIterator(
                libs, intervalIt, header.getSequenceDictionary());
    }

    private static void pairhmmNativeLikelihoods(final String[] args) throws Exception {
        if (args.length < 1) {
            usage();
        }
        HcParityNativePairHmm.dumpCases(Paths.get(args[0]), System.out);
    }

    private static final int MINIMUM_PUTATIVE_PLOIDY_FOR_ACTIVE_REGION_DISCOVERY = 2;

    private static void ploidyResolution(final String[] args) throws Exception {
        if (args.length < 1) {
            usage();
        }
        final int samplePloidy = Integer.parseInt(args[0].trim());
        final int activityEval =
                Math.max(MINIMUM_PUTATIVE_PLOIDY_FOR_ACTIVE_REGION_DISCOVERY, samplePloidy);
        System.out.println("sample_ploidy\t" + samplePloidy);
        System.out.println("activity_eval_ploidy\t" + activityEval);
        System.out.println("genotyping_ploidy\t" + samplePloidy);
    }

    private static void genotypingAggregate(final String[] args) throws Exception {
        if (args.length < 1) {
            usage();
        }
        final Path casesPath = Paths.get(args[0]);
        final List<String> hapLabels = new ArrayList<>();
        final List<double[]> readRows = new ArrayList<>();
        try (BufferedReader br = Files.newBufferedReader(casesPath, StandardCharsets.UTF_8)) {
            String line;
            while ((line = br.readLine()) != null) {
                final String t = line.trim();
                if (t.isEmpty() || t.startsWith("#")) {
                    continue;
                }
                final String[] c = t.split("\t");
                if (c.length < 5) {
                    throw new IllegalArgumentException("pairhmm cases row needs 5 cols: " + line);
                }
                final String readBases = c[1];
                final String[] qparts = c[2].split(",");
                final byte[] quals = new byte[qparts.length];
                for (int i = 0; i < qparts.length; i++) {
                    quals[i] = (byte) Integer.parseInt(qparts[i].trim());
                }
                final int mapq = Integer.parseInt(c[3].trim());
                final String hap = c[4];
                int hapIdx = hapLabels.indexOf(hap);
                if (hapIdx < 0) {
                    hapLabels.add(hap);
                    hapIdx = hapLabels.size() - 1;
                }
                final double ll =
                        HcParityScalarPairHmm.pairhmmLog10Likelihood(readBases, quals, mapq, hap);
                double[] row = null;
                for (final double[] existing : readRows) {
                    row = existing;
                    break;
                }
                if (row == null) {
                    row = new double[hapLabels.size()];
                    Arrays.fill(row, Double.NEGATIVE_INFINITY);
                    readRows.add(row);
                }
                if (row.length < hapLabels.size()) {
                    row = Arrays.copyOf(row, hapLabels.size());
                    readRows.set(0, row);
                }
                row[hapIdx] = ll;
            }
        }
        final int hapCount = hapLabels.size();
        final double[] sums = new double[hapCount];
        for (final double[] row : readRows) {
            for (int i = 0; i < hapCount; i++) {
                sums[i] += row[i];
            }
        }
        int best = 0;
        for (int i = 1; i < hapCount; i++) {
            if (sums[i] > sums[best]) {
                best = i;
            }
        }
        System.out.println("haplotype_count\t" + hapCount);
        System.out.println("read_count\t" + readRows.size());
        for (int i = 0; i < hapCount; i++) {
            System.out.printf(Locale.ROOT, "haplotype_%d_log10_sum\t%s%n", i, Double.toString(sums[i]));
        }
        System.out.println("best_haplotype_index\t" + best);
    }

    private static void genotypeFormat(final String[] args) throws Exception {
        if (args.length < 1) {
            usage();
        }
        final Path fixture = Paths.get(args[0]);
        final List<String> lines = Files.readAllLines(fixture, StandardCharsets.UTF_8);
        System.out.println("# case_id\tpl\tgq\tad\tdp");
        for (final String raw : lines) {
            final String line = raw.trim();
            if (line.isEmpty() || line.startsWith("#")) {
                continue;
            }
            final String[] c = line.split("\t");
            if (c.length != 3) {
                throw new IllegalArgumentException("genotype-format fixture row needs 3 cols: " + line);
            }
            final double[] gl = parseCsvDoubles(c[1]);
            final int[] ad = parseCsvInts(c[2]);
            final P7GenotypeFieldsDump.Fields f = P7GenotypeFieldsDump.emit(gl, ad);
            System.out.println(P7GenotypeFieldsDump.fmtRow(c[0], f));
        }
    }

    private static void variantVcfFromGlAd(final String[] args) throws Exception {
        if (args.length < 6) {
            usage();
        }
        final String contig = args[0];
        final long pos = Long.parseLong(args[1]);
        final String ref = args[2];
        final String alt = args[3];
        final double[] gl = parseCsvDoubles(args[4]);
        final int[] ad = parseCsvInts(args[5]);
        HcParityRegionVcf.dumpSyntheticVcf(contig, pos, ref, alt, gl, ad);
    }

    private static void variantFormatFromGlAd(final String[] args) throws Exception {
        if (args.length < 6) {
            usage();
        }
        final String contig = args[0];
        final long pos = Long.parseLong(args[1]);
        final String ref = args[2];
        final String alt = args[3];
        final double[] gl = parseCsvDoubles(args[4]);
        final int[] ad = parseCsvInts(args[5]);
        HcParityRegionVcf.dumpSyntheticFormat(contig, pos, ref, alt, gl, ad);
    }

    private static void annotateCore(final String[] args) throws Exception {
        if (args.length < 2) {
            usage();
        }
        final int altCount = Integer.parseInt(args[0]);
        final List<HcParityCoreAnnotations.SampleRow> samples =
                HcParityCoreAnnotations.parseSamplesTsv(Paths.get(args[1]));
        final HcParityCoreAnnotations.CoreResult site =
                HcParityCoreAnnotations.compute(altCount, samples);
        HcParityCoreAnnotations.dumpAnnotatedSite(site);
    }

    private static double[] parseCsvDoubles(final String raw) {
        final String[] parts = raw.split(",");
        final double[] out = new double[parts.length];
        for (int i = 0; i < parts.length; i++) {
            out[i] = Double.parseDouble(parts[i]);
        }
        return out;
    }

    private static int[] parseCsvInts(final String raw) {
        final String[] parts = raw.split(",");
        final int[] out = new int[parts.length];
        for (int i = 0; i < parts.length; i++) {
            out[i] = Integer.parseInt(parts[i]);
        }
        return out;
    }

    private static void assemblyRegionGenotype(final String[] args) throws Exception {
        if (args.length < 3) {
            usage();
        }
        final String refPath = args[0];
        final String bamPath = args[1];
        final String intervalCli = args[2];
        final boolean[] wantInactive = {false};
        final int padding = parsePaddingAndTargetIndex(args, 3, wantInactive);
        final boolean dumpInactive = wantInactive[0];
        try (HcContext ctx = new HcContext(refPath, bamPath, padding);
                CachingIndexedFastaSequenceFile refReader =
                        new CachingIndexedFastaSequenceFile(Paths.get(refPath))) {
            final Field hcArgsField = HaplotypeCallerEngine.class.getDeclaredField("hcArgs");
            hcArgsField.setAccessible(true);
            final HaplotypeCallerArgumentCollection hcArgs =
                    (HaplotypeCallerArgumentCollection) hcArgsField.get(ctx.engine);
            final Field assemblerField =
                    HaplotypeCallerEngine.class.getDeclaredField("assemblyEngine");
            assemblerField.setAccessible(true);
            final ReadThreadingAssembler assembler =
                    (ReadThreadingAssembler) assemblerField.get(ctx.engine);
            final Field alignerField = HaplotypeCallerEngine.class.getDeclaredField("aligner");
            alignerField.setAccessible(true);
            final SmithWatermanAligner aligner =
                    (SmithWatermanAligner) alignerField.get(ctx.engine);
            final Logger logger = LogManager.getLogger(HcFullParityGateDump.class);
            final SampleList samplesList = SampleList.singletonSampleList("s1");

            final List<SimpleInterval> intervals =
                    parseIntervals(ctx.header.getSequenceDictionary(), intervalCli);
            final List<List<Locatable>> byContig = groupByContig(intervals);
            for (final List<Locatable> contigIntervals : byContig) {
                final List<SimpleInterval> contigSimple =
                        contigIntervals.stream()
                                .map(SimpleInterval::new)
                                .collect(Collectors.toList());
                final MultiIntervalLocalReadShard shard =
                        new MultiIntervalLocalReadShard(
                                contigSimple, padding, ctx.readsSource);
                configureHcProductionReadShard(shard, ctx);
                final AssemblyRegionIterator iter =
                        new AssemblyRegionIterator(
                                shard,
                                ctx.header,
                                ctx.reference,
                                null,
                                ctx.engine,
                                ctx.asmArgs,
                                false);
                while (iter.hasNext()) {
                    final AssemblyRegion r = iter.next();
                    if (dumpInactive) {
                        if (r.isActive()) {
                            continue;
                        }
                        System.out.println("region_contig\t" + r.getContig());
                        System.out.println("region_start\t" + r.getStart());
                        System.out.println("region_end\t" + r.getEnd());
                        System.out.println("is_active\tfalse");
                        System.out.println("genotyped\tfalse");
                        System.out.println("haplotype_count\t0");
                        System.out.println("read_count\t0");
                        return;
                    }
                    if (!r.isActive()) {
                        continue;
                    }
                    final AssemblyResultSet ars =
                            AssemblyBasedCallerUtils.assembleReads(
                                    r,
                                    Collections.emptyList(),
                                    hcArgs,
                                    ctx.header,
                                    samplesList,
                                    logger,
                                    refReader,
                                    assembler,
                                    aligner,
                                    !hcArgs.doNotCorrectOverlappingBaseQualities,
                                    hcArgs.fbargs,
                                    false);
                    System.out.println("region_contig\t" + r.getContig());
                    System.out.println("region_start\t" + r.getStart());
                    System.out.println("region_end\t" + r.getEnd());
                    System.out.println("is_active\ttrue");
                    final List<org.broadinstitute.hellbender.utils.read.GATKRead> reads =
                            new ArrayList<>(r.getReads());
                    HcParityRegionGenotype.sortReads(reads);
                    final List<Haplotype> haps = new ArrayList<>(ars.getHaplotypeList());
                    final List<double[]> readRows = new ArrayList<>();
                    final boolean[] isRef = new boolean[haps.size()];
                    for (int hi = 0; hi < haps.size(); hi++) {
                        isRef[hi] = haps.get(hi).isReference();
                    }
                    for (final org.broadinstitute.hellbender.utils.read.GATKRead read : reads) {
                        final byte[] quals = read.getBaseQualities();
                        final String readBases =
                                new String(read.getBases(), StandardCharsets.US_ASCII);
                        final double[] row = new double[haps.size()];
                        Arrays.fill(row, Double.NEGATIVE_INFINITY);
                        for (int hi = 0; hi < haps.size(); hi++) {
                            final String hapBases =
                                    new String(
                                            haps.get(hi).getBases(), StandardCharsets.US_ASCII);
                            row[hi] =
                                    HcParityNativePairHmm.pairhmmLog10Likelihood(
                                            readBases, quals, read.getMappingQuality(), hapBases);
                        }
                        readRows.add(row);
                    }
                    final HcParityRegionGenotype.GenotypeDump dump =
                            HcParityRegionGenotype.genotypeFromLikelihoodMatrix(readRows, isRef);
                    HcParityRegionGenotype.printDump(dump);
                    return;
                }
            }
            throw new IllegalArgumentException("no active assembly region in interval");
        }
    }

    private static void applyAssemblyProfile(
            final ReadThreadingAssembler assembler, final String profile) throws Exception {
        if (!"sensitive".equals(profile)) {
            return;
        }
        final Field kmerField = ReadThreadingAssembler.class.getDeclaredField("kmerSizes");
        kmerField.setAccessible(true);
        kmerField.set(assembler, new ArrayList<>(Arrays.asList(3, 5, 10)));
        final Field recoverAllField =
                ReadThreadingAssembler.class.getDeclaredField("recoverAllDanglingBranches");
        recoverAllField.setAccessible(true);
        recoverAllField.setBoolean(assembler, true);
        final Field dontIncreaseField =
                ReadThreadingAssembler.class.getDeclaredField("dontIncreaseKmerSizesForCycles");
        dontIncreaseField.setAccessible(true);
        dontIncreaseField.setBoolean(assembler, false);
        final Field minDanglingField =
                ReadThreadingAssembler.class.getDeclaredField("minDanglingBranchLength");
        minDanglingField.setAccessible(true);
        minDanglingField.setInt(assembler, 2);
    }

    private static void assemblyRegionGenotypeSubset(final String[] args) throws Exception {
        if (args.length < 5) {
            usage();
        }
        final int maxAlleles = Integer.parseInt(args[args.length - 1]);
        final String assemblyProfile = args[args.length - 2];
        final String[] base = Arrays.copyOf(args, args.length - 2);
        if (base.length < 3) {
            usage();
        }
        final String refPath = base[0];
        final String bamPath = base[1];
        final String intervalCli = base[2];
        final boolean[] wantInactive = {false};
        final int padding = parsePaddingAndTargetIndex(base, 3, wantInactive);
        final boolean dumpInactive = wantInactive[0];
        try (HcContext ctx = new HcContext(refPath, bamPath, padding);
                CachingIndexedFastaSequenceFile refReader =
                        new CachingIndexedFastaSequenceFile(Paths.get(refPath))) {
            final Field hcArgsField = HaplotypeCallerEngine.class.getDeclaredField("hcArgs");
            hcArgsField.setAccessible(true);
            final HaplotypeCallerArgumentCollection hcArgs =
                    (HaplotypeCallerArgumentCollection) hcArgsField.get(ctx.engine);
            final Field assemblerField =
                    HaplotypeCallerEngine.class.getDeclaredField("assemblyEngine");
            assemblerField.setAccessible(true);
            final ReadThreadingAssembler assembler =
                    (ReadThreadingAssembler) assemblerField.get(ctx.engine);
            applyAssemblyProfile(assembler, assemblyProfile);
            final Field alignerField = HaplotypeCallerEngine.class.getDeclaredField("aligner");
            alignerField.setAccessible(true);
            final SmithWatermanAligner aligner =
                    (SmithWatermanAligner) alignerField.get(ctx.engine);
            final Logger logger = LogManager.getLogger(HcFullParityGateDump.class);
            final SampleList samplesList = SampleList.singletonSampleList("s1");

            final List<SimpleInterval> intervals =
                    parseIntervals(ctx.header.getSequenceDictionary(), intervalCli);
            final List<List<Locatable>> byContig = groupByContig(intervals);
            for (final List<Locatable> contigIntervals : byContig) {
                final List<SimpleInterval> contigSimple =
                        contigIntervals.stream()
                                .map(SimpleInterval::new)
                                .collect(Collectors.toList());
                final MultiIntervalLocalReadShard shard =
                        new MultiIntervalLocalReadShard(
                                contigSimple, padding, ctx.readsSource);
                configureHcProductionReadShard(shard, ctx);
                final AssemblyRegionIterator iter =
                        new AssemblyRegionIterator(
                                shard,
                                ctx.header,
                                ctx.reference,
                                null,
                                ctx.engine,
                                ctx.asmArgs,
                                false);
                while (iter.hasNext()) {
                    final AssemblyRegion r = iter.next();
                    if (dumpInactive) {
                        if (r.isActive()) {
                            continue;
                        }
                        System.out.println("region_contig\t" + r.getContig());
                        System.out.println("region_start\t" + r.getStart());
                        System.out.println("region_end\t" + r.getEnd());
                        System.out.println("is_active\tfalse");
                        System.out.println("genotyped\tfalse");
                        System.out.println("haplotype_count\t0");
                        System.out.println("read_count\t0");
                        return;
                    }
                    if (!r.isActive()) {
                        continue;
                    }
                    final AssemblyResultSet ars =
                            AssemblyBasedCallerUtils.assembleReads(
                                    r,
                                    Collections.emptyList(),
                                    hcArgs,
                                    ctx.header,
                                    samplesList,
                                    logger,
                                    refReader,
                                    assembler,
                                    aligner,
                                    !hcArgs.doNotCorrectOverlappingBaseQualities,
                                    hcArgs.fbargs,
                                    false);
                    System.out.println("region_contig\t" + r.getContig());
                    System.out.println("region_start\t" + r.getStart());
                    System.out.println("region_end\t" + r.getEnd());
                    System.out.println("is_active\ttrue");
                    System.out.println("assembly_profile\t" + assemblyProfile);
                    final List<org.broadinstitute.hellbender.utils.read.GATKRead> reads =
                            new ArrayList<>(r.getReads());
                    HcParityRegionGenotype.sortReads(reads);
                    final List<Haplotype> haps = new ArrayList<>(ars.getHaplotypeList());
                    final List<double[]> readRows = new ArrayList<>();
                    final boolean[] isRef = new boolean[haps.size()];
                    for (int hi = 0; hi < haps.size(); hi++) {
                        isRef[hi] = haps.get(hi).isReference();
                    }
                    for (final org.broadinstitute.hellbender.utils.read.GATKRead read : reads) {
                        final byte[] quals = read.getBaseQualities();
                        final String readBases =
                                new String(read.getBases(), StandardCharsets.US_ASCII);
                        final double[] row = new double[haps.size()];
                        Arrays.fill(row, Double.NEGATIVE_INFINITY);
                        for (int hi = 0; hi < haps.size(); hi++) {
                            final String hapBases =
                                    new String(
                                            haps.get(hi).getBases(), StandardCharsets.US_ASCII);
                            row[hi] =
                                    HcParityNativePairHmm.pairhmmLog10Likelihood(
                                            readBases, quals, read.getMappingQuality(), hapBases);
                        }
                        readRows.add(row);
                    }
                    final HcParityRegionGenotype.GenotypeDump dump =
                            HcParityRegionGenotype.genotypeFromLikelihoodMatrix(readRows, isRef);
                    HcParityRegionGenotype.printDump(dump);
                    final double[] sums = new double[haps.size()];
                    for (final double[] row : readRows) {
                        for (int i = 0; i < haps.size(); i++) {
                            sums[i] += row[i];
                        }
                    }
                    HcParityAlleleSubsetting.dumpLiveSubsetExtension(sums, isRef, haps, maxAlleles);
                    return;
                }
            }
            throw new IllegalArgumentException("no active assembly region in interval");
        }
    }

    private static void callRegionVcf(final String[] args) throws Exception {
        if (args.length < 3) {
            usage();
        }
        final String refPath = args[0];
        final String bamPath = args[1];
        final String intervalCli = args[2];
        final int padding = parsePaddingAndTargetIndex(args, 3, new boolean[] {false});
        try (HcContext ctx = new HcContext(refPath, bamPath, padding);
                CachingIndexedFastaSequenceFile refReader =
                        new CachingIndexedFastaSequenceFile(Paths.get(refPath))) {
            final Field hcArgsField = HaplotypeCallerEngine.class.getDeclaredField("hcArgs");
            hcArgsField.setAccessible(true);
            final HaplotypeCallerArgumentCollection hcArgs =
                    (HaplotypeCallerArgumentCollection) hcArgsField.get(ctx.engine);
            final Field assemblerField =
                    HaplotypeCallerEngine.class.getDeclaredField("assemblyEngine");
            assemblerField.setAccessible(true);
            final ReadThreadingAssembler assembler =
                    (ReadThreadingAssembler) assemblerField.get(ctx.engine);
            final Field alignerField = HaplotypeCallerEngine.class.getDeclaredField("aligner");
            alignerField.setAccessible(true);
            final SmithWatermanAligner aligner =
                    (SmithWatermanAligner) alignerField.get(ctx.engine);
            final Logger logger = LogManager.getLogger(HcFullParityGateDump.class);
            final SampleList samplesList = SampleList.singletonSampleList("s1");

            final List<SimpleInterval> intervals =
                    parseIntervals(ctx.header.getSequenceDictionary(), intervalCli);
            final List<List<Locatable>> byContig = groupByContig(intervals);
            for (final List<Locatable> contigIntervals : byContig) {
                final List<SimpleInterval> contigSimple =
                        contigIntervals.stream()
                                .map(SimpleInterval::new)
                                .collect(Collectors.toList());
                final MultiIntervalLocalReadShard shard =
                        new MultiIntervalLocalReadShard(
                                contigSimple, padding, ctx.readsSource);
                configureHcProductionReadShard(shard, ctx);
                final AssemblyRegionIterator iter =
                        new AssemblyRegionIterator(
                                shard,
                                ctx.header,
                                ctx.reference,
                                null,
                                ctx.engine,
                                ctx.asmArgs,
                                false);
                while (iter.hasNext()) {
                    final AssemblyRegion r = iter.next();
                    if (!r.isActive()) {
                        continue;
                    }
                    final AssemblyResultSet ars =
                            AssemblyBasedCallerUtils.assembleReads(
                                    r,
                                    Collections.emptyList(),
                                    hcArgs,
                                    ctx.header,
                                    samplesList,
                                    logger,
                                    refReader,
                                    assembler,
                                    aligner,
                                    !hcArgs.doNotCorrectOverlappingBaseQualities,
                                    hcArgs.fbargs,
                                    false);
                    HcParityCallRegionVcf.emitFromRegion(r, ars);
                    return;
                }
            }
        }
        HcParityRegionVcf.dumpVariantVcf(false, "", 0, "", "", "", "", "");
    }

    /**
     * Live GATK 4.4 {@code HaplotypeCallerEngine.callRegion} with TEST-ONLY dump of
     * {@code DepthPerAlleleBySample.annotate} inputs (6R.91). Args: ref bam interval [padding].
     */
    private static void adAnnotationCall(final String[] args) throws Exception {
        if (args.length < 3) {
            usage();
        }
        final String refPath = args[0];
        final String bamPath = args[1];
        final String intervalCli = args[2];
        final int padding = parsePaddingAndTargetIndex(args, 3, new boolean[] {false});
        try (HcContext ctx = new HcContext(refPath, bamPath, padding)) {
            HcParityAdAnnotationDump.installOn(ctx.engine);
            final List<SimpleInterval> intervals =
                    parseIntervals(ctx.header.getSequenceDictionary(), intervalCli);
            final List<List<Locatable>> byContig = groupByContig(intervals);
            for (final List<Locatable> contigIntervals : byContig) {
                final List<SimpleInterval> contigSimple =
                        contigIntervals.stream()
                                .map(SimpleInterval::new)
                                .collect(Collectors.toList());
                final MultiIntervalLocalReadShard shard =
                        new MultiIntervalLocalReadShard(
                                contigSimple, padding, ctx.readsSource);
                configureHcProductionReadShard(shard, ctx);
                final AssemblyRegionIterator iter =
                        new AssemblyRegionIterator(
                                shard,
                                ctx.header,
                                ctx.reference,
                                null,
                                ctx.engine,
                                ctx.asmArgs,
                                false);
                while (iter.hasNext()) {
                    final AssemblyRegion r = iter.next();
                    if (!r.isActive()) {
                        continue;
                    }
                    System.out.println(
                            "6R91\t0\tregion\t"
                                    + r.getContig()
                                    + ":"
                                    + r.getStart()
                                    + "-"
                                    + r.getEnd()
                                    + "\treads="
                                    + r.getReads().size());
                    final FeatureContext features = new FeatureContext();
                    final ReferenceContext refCtx =
                            new ReferenceContext(
                                    ctx.reference, r.getSpan(), padding, padding);
                    final List<VariantContext> calls =
                            ctx.engine.callRegion(r, features, refCtx);
                    System.out.println("6R91\t0\tcall_region_n\t" + calls.size());
                    for (int i = 0; i < calls.size(); i++) {
                        final VariantContext vc = calls.get(i);
                        final Genotype gt =
                                vc.getGenotypes().isEmpty() ? null : vc.getGenotype(0);
                        System.out.println(
                                "6R91\t0\tcall_region_vc_"
                                        + i
                                        + "\t"
                                        + vc.getContig()
                                        + ":"
                                        + vc.getStart()
                                        + "-"
                                        + vc.getEnd()
                                        + "\talleles="
                                        + HcParityAdAnnotationDump.alleleList(vc.getAlleles())
                                        + "\tad="
                                        + (gt != null && gt.hasAD()
                                                ? HcParityAdAnnotationDump.ints(gt.getAD())
                                                : "ABSENT")
                                        + "\tgt="
                                        + (gt == null ? "." : gt.getGenotypeString()));
                    }
                }
            }
        }
        System.out.println("6R91\t0\tno_active_region\ttrue");
    }

    /**
     * Live GATK 4.4 {@code HaplotypeCallerEngine.callRegion} with TEST-ONLY dump of
     * {@code filterPoorlyModeledEvidence} inputs (6R.93). Args: ref bam interval [padding].
     */
    private static void filterPoorlyModeledCall(final String[] args) throws Exception {
        if (args.length < 3) {
            usage();
        }
        final String refPath = args[0];
        final String bamPath = args[1];
        final String intervalCli = args[2];
        final int padding = parsePaddingAndTargetIndex(args, 3, new boolean[] {false});
        try (HcContext ctx = new HcContext(refPath, bamPath, padding)) {
            HcParityFilterPoorlyModeledDump.installOn(ctx.engine);
            final List<SimpleInterval> intervals =
                    parseIntervals(ctx.header.getSequenceDictionary(), intervalCli);
            final List<List<Locatable>> byContig = groupByContig(intervals);
            for (final List<Locatable> contigIntervals : byContig) {
                final List<SimpleInterval> contigSimple =
                        contigIntervals.stream()
                                .map(SimpleInterval::new)
                                .collect(Collectors.toList());
                final MultiIntervalLocalReadShard shard =
                        new MultiIntervalLocalReadShard(
                                contigSimple, padding, ctx.readsSource);
                configureHcProductionReadShard(shard, ctx);
                final AssemblyRegionIterator iter =
                        new AssemblyRegionIterator(
                                shard,
                                ctx.header,
                                ctx.reference,
                                null,
                                ctx.engine,
                                ctx.asmArgs,
                                false);
                while (iter.hasNext()) {
                    final AssemblyRegion r = iter.next();
                    if (!r.isActive()) {
                        continue;
                    }
                    System.out.println(
                            "6R93\t0\tregion\t"
                                    + r.getContig()
                                    + ":"
                                    + r.getStart()
                                    + "-"
                                    + r.getEnd()
                                    + "\treads="
                                    + r.getReads().size());
                    final FeatureContext features = new FeatureContext();
                    final ReferenceContext refCtx =
                            new ReferenceContext(
                                    ctx.reference, r.getSpan(), padding, padding);
                    final List<VariantContext> calls =
                            ctx.engine.callRegion(r, features, refCtx);
                    System.out.println("6R93\t0\tcall_region_n\t" + calls.size());
                    for (int i = 0; i < calls.size(); i++) {
                        final VariantContext vc = calls.get(i);
                        final Genotype gt =
                                vc.getGenotypes().isEmpty() ? null : vc.getGenotype(0);
                        System.out.println(
                                "6R93\t0\tcall_region_vc_"
                                        + i
                                        + "\t"
                                        + vc.getContig()
                                        + ":"
                                        + vc.getStart()
                                        + "-"
                                        + vc.getEnd()
                                        + "\talleles="
                                        + HcParityAdAnnotationDump.alleleList(vc.getAlleles())
                                        + "\tad="
                                        + (gt != null && gt.hasAD()
                                                ? HcParityAdAnnotationDump.ints(gt.getAD())
                                                : "ABSENT")
                                        + "\tgt="
                                        + (gt == null ? "." : gt.getGenotypeString())
                                        + "\tqual="
                                        + vc.getPhredScaledQual());
                    }
                }
            }
        }
        System.out.println("6R93\t0\tno_active_region\ttrue");
    }

    /**
     * Same live dump as {@link #filterPoorlyModeledCall} with GATK
     * {@code --native-pair-hmm-use-double-precision} equivalent set before engine
     * construction. TEST-ONLY 6R.98 oracle; does not change default dump path.
     */
    private static void filterPoorlyModeledCallDouble(final String[] args) throws Exception {
        if (args.length < 3) {
            usage();
        }
        final String refPath = args[0];
        final String bamPath = args[1];
        final String intervalCli = args[2];
        final int padding = parsePaddingAndTargetIndex(args, 3, new boolean[] {false});
        HcParityFilterPoorlyModeledDump.DUMP_PREFIX = "6R98";
        HcParityFilterPoorlyModeledDump.USE_DOUBLE_PRECISION = true;
        try (HcContext ctx = new HcContext(refPath, bamPath, padding, null, true)) {
            HcParityFilterPoorlyModeledDump.installOn(ctx.engine);
            final List<SimpleInterval> intervals =
                    parseIntervals(ctx.header.getSequenceDictionary(), intervalCli);
            final List<List<Locatable>> byContig = groupByContig(intervals);
            for (final List<Locatable> contigIntervals : byContig) {
                final List<SimpleInterval> contigSimple =
                        contigIntervals.stream()
                                .map(SimpleInterval::new)
                                .collect(Collectors.toList());
                final MultiIntervalLocalReadShard shard =
                        new MultiIntervalLocalReadShard(
                                contigSimple, padding, ctx.readsSource);
                configureHcProductionReadShard(shard, ctx);
                final AssemblyRegionIterator iter =
                        new AssemblyRegionIterator(
                                shard,
                                ctx.header,
                                ctx.reference,
                                null,
                                ctx.engine,
                                ctx.asmArgs,
                                false);
                while (iter.hasNext()) {
                    final AssemblyRegion r = iter.next();
                    if (!r.isActive()) {
                        continue;
                    }
                    System.out.println(
                            "6R93\t0\tregion\t"
                                    + r.getContig()
                                    + ":"
                                    + r.getStart()
                                    + "-"
                                    + r.getEnd()
                                    + "\treads="
                                    + r.getReads().size());
                    final FeatureContext features = new FeatureContext();
                    final ReferenceContext refCtx =
                            new ReferenceContext(
                                    ctx.reference, r.getSpan(), padding, padding);
                    final List<VariantContext> calls =
                            ctx.engine.callRegion(r, features, refCtx);
                    System.out.println("6R93\t0\tcall_region_n\t" + calls.size());
                    for (int i = 0; i < calls.size(); i++) {
                        final VariantContext vc = calls.get(i);
                        final Genotype gt =
                                vc.getGenotypes().isEmpty() ? null : vc.getGenotype(0);
                        System.out.println(
                                "6R93\t0\tcall_region_vc_"
                                        + i
                                        + "\t"
                                        + vc.getContig()
                                        + ":"
                                        + vc.getStart()
                                        + "-"
                                        + vc.getEnd()
                                        + "\talleles="
                                        + HcParityAdAnnotationDump.alleleList(vc.getAlleles())
                                        + "\tad="
                                        + (gt != null && gt.hasAD()
                                                ? HcParityAdAnnotationDump.ints(gt.getAD())
                                                : "ABSENT")
                                        + "\tgt="
                                        + (gt == null ? "." : gt.getGenotypeString())
                                        + "\tqual="
                                        + vc.getPhredScaledQual());
                    }
                }
            }
        }
        System.out.println("6R93\t0\tno_active_region\ttrue");
    }

    private static void pairhmmLikelihoods(final String[] args) throws Exception {
        if (args.length < 1) {
            usage();
        }
        final Path casesPath = Paths.get(args[0]);
        System.out.println("case_id\tlog10_likelihood");
        try (BufferedReader br = Files.newBufferedReader(casesPath, StandardCharsets.UTF_8)) {
            String line;
            while ((line = br.readLine()) != null) {
                final String t = line.trim();
                if (t.isEmpty() || t.startsWith("#")) {
                    continue;
                }
                final String[] c = t.split("\t");
                if (c.length < 5) {
                    throw new IllegalArgumentException("pairhmm cases row needs 5 cols: " + line);
                }
                final String caseId = c[0];
                final String readBases = c[1];
                final String[] qparts = c[2].split(",");
                final byte[] quals = new byte[qparts.length];
                for (int i = 0; i < qparts.length; i++) {
                    quals[i] = (byte) Integer.parseInt(qparts[i].trim());
                }
                final int mapq = Integer.parseInt(c[3].trim());
                final String hap = c[4];
                final double ll =
                        HcParityScalarPairHmm.pairhmmLog10Likelihood(readBases, quals, mapq, hap);
                System.out.printf(Locale.ROOT, "%s\t%s%n", caseId, Double.toString(ll));
            }
        }
    }

    private static void downsamplePositional(final String[] args) throws Exception {
        if (args.length < 2) {
            usage();
        }
        final Path bamPath = Paths.get(args[0]);
        final int cap = Integer.parseInt(args[1]);
        final boolean nonRandom =
                args.length < 3 || !args[2].equalsIgnoreCase("random");
        System.out.println("alignment_start\tkept_count\tkept_qnames");
        try (SamReader reader = SamReaderFactory.makeDefault().open(bamPath)) {
            Utils.resetRandomGenerator();
            final SAMFileHeader header = reader.getFileHeader();
            final TreeMap<Integer, List<SAMRecord>> byStart = new TreeMap<>();
            for (final SAMRecord rec : reader) {
                byStart.computeIfAbsent(rec.getAlignmentStart(), k -> new ArrayList<>()).add(rec);
            }
            for (final Map.Entry<Integer, List<SAMRecord>> e : byStart.entrySet()) {
                final int start1Based = e.getKey();
                // Match rust_htslib Record::pos(): 0-based; htsjdk SAMRecord.getAlignmentStart() is 1-based when mapped.
                final int alignmentStart0 = start1Based <= 0 ? -1 : start1Based - 1;
                final List<GATKRead> greads = new ArrayList<>();
                for (final SAMRecord r : e.getValue()) {
                    greads.add(new SAMRecordToGATKReadAdapter(r));
                }
                final PositionalDownsampler pd =
                        new PositionalDownsampler(cap, header, nonRandom);
                for (final GATKRead g : greads) {
                    pd.submit(g);
                }
                pd.signalEndOfInput();
                final List<GATKRead> kept = pd.consumeFinalizedItems();
                final List<String> names = new ArrayList<>();
                for (final GATKRead g : kept) {
                    names.add(g.getName());
                }
                Collections.sort(names);
                System.out.printf(
                        Locale.ROOT,
                        "%d\t%d\t%s%n",
                        alignmentStart0,
                        names.size(),
                        String.join(",", names));
            }
        }
    }

    private enum AlleleRoundRobinClass {
        REF,
        ALT
    }

    private static boolean passesRustDefaultReadFilter(final SAMRecord rec) {
        if (rec.getMappingQuality() < 20) {
            return false;
        }
        if (rec.isSecondaryAlignment()) {
            return false;
        }
        if (rec.getSupplementaryAlignmentFlag()) {
            return false;
        }
        if (rec.getDuplicateReadFlag()) {
            return false;
        }
        return true;
    }

    /** Match Rust {@code query_index_at_reference_position} (0-based alignment start and ref pos). */
    private static Byte queryBaseAtRef0(final SAMRecord rec, final int ref0) {
        if (rec.getReadUnmappedFlag()) {
            return null;
        }
        final int start1 = rec.getAlignmentStart();
        if (start1 <= 0) {
            return null;
        }
        final int start0 = start1 - 1;
        if (ref0 < start0) {
            return null;
        }
        int r = start0;
        int q = 0;
        for (final htsjdk.samtools.CigarElement el : rec.getCigar().getCigarElements()) {
            final int n = el.getLength();
            switch (el.getOperator()) {
                case M:
                case EQ:
                case X:
                    if (ref0 < r + n) {
                        final int qi = q + (ref0 - r);
                        final byte[] bases = rec.getReadBases();
                        if (qi < 0 || qi >= bases.length) {
                            return null;
                        }
                        return bases[qi];
                    }
                    r += n;
                    q += n;
                    break;
                case I:
                case S:
                    q += n;
                    break;
                case D:
                case N:
                    if (ref0 < r + n) {
                        return null;
                    }
                    r += n;
                    break;
                case H:
                case P:
                    break;
                default:
                    break;
            }
        }
        return null;
    }

    /** Mirror of package-private {@link AlleleBiasedDownsamplingUtils#targetAlleleCounts} for parity dumps. */
    private static int[] gateTargetAlleleCounts(final int[] alleleCounts, final int numReadsToRemove) {
        final int numAlleles = alleleCounts.length;
        int maxScore = gateScoreAlleleCounts(alleleCounts);
        int[] alleleCountsOfMax = alleleCounts;
        final int numReadsToRemovePerAllele = numReadsToRemove / 2;
        for (int i = 0; i < numAlleles; i++) {
            for (int j = i; j < numAlleles; j++) {
                final int[] newCounts = alleleCounts.clone();
                if (i == j) {
                    newCounts[i] = Math.max(0, newCounts[i] - numReadsToRemove);
                } else {
                    newCounts[i] = Math.max(0, newCounts[i] - numReadsToRemovePerAllele);
                    newCounts[j] = Math.max(0, newCounts[j] - numReadsToRemovePerAllele);
                }
                final int score = gateScoreAlleleCounts(newCounts);
                if (score < maxScore) {
                    maxScore = score;
                    alleleCountsOfMax = newCounts;
                }
            }
        }
        return alleleCountsOfMax;
    }

    private static int gateScoreAlleleCounts(final int[] alleleCounts) {
        if (alleleCounts.length < 2) {
            return 0;
        }
        final int[] sorted = alleleCounts.clone();
        Arrays.sort(sorted);
        final int maxCount = sorted[sorted.length - 1];
        final int nextBestCount = sorted[sorted.length - 2];
        final int remainderCount =
                Arrays.stream(sorted).sum() - maxCount - nextBestCount;
        return Math.min(
                maxCount - nextBestCount + remainderCount,
                Math.abs(nextBestCount + remainderCount));
    }

    private static void alleleBiasedTargetCounts(final String[] args) {
        if (args.length < 2) {
            usage();
        }
        final String[] parts = args[0].split(",");
        final int[] counts = new int[parts.length];
        for (int i = 0; i < parts.length; i++) {
            counts[i] = Integer.parseInt(parts[i].trim());
        }
        final int numRemove = Integer.parseInt(args[1]);
        final int[] target = gateTargetAlleleCounts(counts, numRemove);
        final StringBuilder sb = new StringBuilder();
        for (int i = 0; i < target.length; i++) {
            if (i > 0) {
                sb.append(',');
            }
            sb.append(target[i]);
        }
        System.out.println("allele_counts\tnum_reads_to_remove\ttarget_counts");
        System.out.printf(Locale.ROOT, "%s\t%d\t%s%n", args[0], numRemove, sb);
    }

    private static void alleleBiasedEvidence(final String[] args) throws Exception {
        if (args.length < 5) {
            usage();
        }
        final String refPath = args[0];
        final Path bamPath = Paths.get(args[1]);
        final String contig = args[2];
        final int pos1 = Integer.parseInt(args[3]);
        final double contamination = Double.parseDouble(args[4]);
        final int ref0 = pos1 - 1;
        Utils.resetRandomGenerator();
        byte refBase;
        try (CachingIndexedFastaSequenceFile refReader =
                new CachingIndexedFastaSequenceFile(Paths.get(refPath))) {
            refBase = refReader.getSubsequenceAt(contig, pos1, pos1).getBases()[0];
        }
        final SAMFileHeader header;
        final List<GATKRead> atLocus = new ArrayList<>();
        try (SamReader reader = SamReaderFactory.makeDefault().open(bamPath)) {
            header = reader.getFileHeader();
            final int tid = header.getSequenceIndex(contig);
            for (final SAMRecord rec : reader) {
                if (rec.getReferenceIndex() != tid) {
                    continue;
                }
                if (!passesRustDefaultReadFilter(rec)) {
                    continue;
                }
                final Byte qb = queryBaseAtRef0(rec, ref0);
                if (qb == null) {
                    continue;
                }
                atLocus.add(new SAMRecordToGATKReadAdapter(rec));
            }
        }
        final Map<Allele, List<GATKRead>> readMap = new LinkedHashMap<>();
        for (final GATKRead read : atLocus) {
            final Byte qb = queryBaseAtRef0(read.convertToSAMRecord(header), ref0);
            if (qb == null) {
                continue;
            }
            final Allele allele = Allele.create(qb);
            readMap.computeIfAbsent(allele, k -> new ArrayList<>()).add(read);
        }
        final List<GATKRead> removed =
                AlleleBiasedDownsamplingUtils.selectAlleleBiasedEvidence(readMap, contamination);
        final List<String> names = new ArrayList<>();
        for (final GATKRead r : removed) {
            names.add(r.getName());
        }
        Collections.sort(names);
        System.out.println("contig\tpos\tcontamination_fraction\tremoved_count\tremoved_qnames");
        System.out.printf(
                Locale.ROOT,
                "%s\t%d\t%s\t%d\t%s%n",
                contig,
                pos1,
                contamination,
                names.size(),
                String.join(",", names));
    }

    private static void rawActivityContam(final String[] args) throws Exception {
        if (args.length < 4) {
            usage();
        }
        final String refPath = args[0];
        final String bamPath = args[1];
        final String intervalCli = args[2];
        final double contamination = Double.parseDouble(args[3]);
        final int padding = args.length > 4 ? parsePadding(args[4]) : DEFAULT_PADDING;
        System.out.println("contig\tpos\tactive_prob\toriginal_active_prob\tkind");
        Utils.resetRandomGenerator();
        try (HcContext ctx = new HcContext(refPath, bamPath, padding)) {
            walkRawActivityWithContamination(ctx, intervalCli, contamination);
        }
    }

    private static void walkRawActivityWithContamination(
            final HcContext ctx, final String intervalCli, final double contaminationFraction)
            throws Exception {
        final List<SimpleInterval> intervals =
                parseIntervals(ctx.header.getSequenceDictionary(), intervalCli);
        for (final SimpleInterval interval : intervals) {
            Utils.resetRandomGenerator();
            final MultiIntervalLocalReadShard shard =
                    new MultiIntervalLocalReadShard(
                            Collections.singletonList(interval),
                            ctx.asmArgs.assemblyRegionPadding,
                            ctx.readsSource);
            final Iterator<AlignmentContext> locusIt = makeLocusIterator(shard, ctx);
            while (locusIt.hasNext()) {
                final AlignmentContext pileup = locusIt.next();
                final SimpleInterval loc = new SimpleInterval(pileup);
                final ReferenceContext refCtx = new ReferenceContext(ctx.reference, loc);
                final ReadPileup filtered =
                        contaminationFilteredPileup(
                                pileup.getBasePileup(), loc.getStart() - 1, contaminationFraction);
                final AlignmentContext contaminated = new AlignmentContext(loc, filtered);
                final FeatureContext featCtx = featureContextForLocus(ctx, loc);
                final ActivityProfileState raw =
                        ctx.engine.isActive(contaminated, refCtx, featCtx);
                final String prob = formatProb(raw.isActiveProb());
                System.out.printf(
                        "%s\t%d\t%s\t%s\t%s%n",
                        loc.getContig(), loc.getStart(), prob, prob, formatKind(raw.getResultState()));
            }
        }
    }

    private static ReadPileup contaminationFilteredPileup(
            final ReadPileup pileup, final int ref0, final double contaminationFraction) {
        if (pileup.size() == 0 || contaminationFraction <= 0.0) {
            return pileup;
        }
        final Map<Allele, List<GATKRead>> readMap = new LinkedHashMap<>();
        for (final PileupElement element : pileup) {
            final GATKRead read = element.getRead();
            final Allele allele = Allele.create(element.getBase());
            readMap.computeIfAbsent(allele, k -> new ArrayList<>()).add(read);
        }
        final List<GATKRead> toRemove =
                AlleleBiasedDownsamplingUtils.selectAlleleBiasedEvidence(readMap, contaminationFraction);
        final Set<String> removeNames = new HashSet<>();
        for (final GATKRead r : toRemove) {
            removeNames.add(r.getName());
        }
        final List<PileupElement> kept = new ArrayList<>();
        for (final PileupElement element : pileup) {
            if (!removeNames.contains(element.getRead().getName())) {
                kept.add(element);
            }
        }
        return new ReadPileup(pileup.getLocation(), kept);
    }

    private static void softClipMean(final String[] args) throws Exception {
        if (args.length < 3) {
            usage();
        }
        final String refPath = args[0];
        final String bamPath = args[1];
        final String intervalCli = args[2];
        System.out.println("contig\tpos\thq_soft_clip_mean");
        Utils.resetRandomGenerator();
        try (HcContext ctx = new HcContext(refPath, bamPath, 0)) {
            final Field rcmField = HaplotypeCallerEngine.class.getDeclaredField("referenceConfidenceModel");
            rcmField.setAccessible(true);
            final ReferenceConfidenceModel rcm =
                    (ReferenceConfidenceModel) rcmField.get(ctx.engine);
            final Field genoField =
                    HaplotypeCallerEngine.class.getDeclaredField("activeRegionEvaluationGenotyperEngine");
            genoField.setAccessible(true);
            final MinimalGenotypingEngine genoEngine =
                    (MinimalGenotypingEngine) genoField.get(ctx.engine);
            final int ploidy = genoEngine.getConfiguration().genotypeArgs.samplePloidy;
            final Field hcArgsField = HaplotypeCallerEngine.class.getDeclaredField("hcArgs");
            hcArgsField.setAccessible(true);
            final HaplotypeCallerArgumentCollection hcArgs =
                    (HaplotypeCallerArgumentCollection) hcArgsField.get(ctx.engine);
            final int minBaseQuality = hcArgs.minBaseQualityScore;

            final List<SimpleInterval> intervals =
                    parseIntervals(ctx.header.getSequenceDictionary(), intervalCli);
            for (final SimpleInterval interval : intervals) {
                final MultiIntervalLocalReadShard shard =
                        new MultiIntervalLocalReadShard(
                                Collections.singletonList(interval), 0, ctx.readsSource);
                final Iterator<AlignmentContext> locusIt = makeLocusIteratorSoftClip(shard, ctx);
                while (locusIt.hasNext()) {
                    final AlignmentContext pileup = locusIt.next();
                    final SimpleInterval loc = new SimpleInterval(pileup);
                    final ReferenceContext refCtx =
                            new ReferenceContext(ctx.reference, loc);
                    final MathUtils.RunningAverage hqAvg = new MathUtils.RunningAverage();
                    rcm.calcGenotypeLikelihoodsOfRefVsAny(
                            ploidy,
                            pileup.getBasePileup(),
                            refCtx.getBase(),
                            (byte) minBaseQuality,
                            hqAvg,
                            false);
                    final double mean = hqAvg.mean();
                    final String formatted =
                            mean == 0.0 ? "0" : String.format(Locale.ROOT, "%.8f", mean);
                    System.out.printf(
                            Locale.ROOT,
                        "%s\t%d\t%s%n",
                            loc.getContig(),
                            loc.getStart(),
                            formatted);
                }
            }
        }
    }

    private static void readShards(final String[] args) throws Exception {
        if (args.length < 2) {
            usage();
        }
        final String refPath = args[0];
        final String intervalCli = args[1];
        final int padding = args.length > 2 ? parsePadding(args[2]) : DEFAULT_PADDING;
        final CachingIndexedFastaSequenceFile refReader =
                new CachingIndexedFastaSequenceFile(Paths.get(refPath));
        try {
            final SAMSequenceDictionary dict = refReader.getSequenceDictionary();
            final List<SimpleInterval> intervals = parseIntervals(dict, intervalCli);
            final List<List<Locatable>> byContig = groupByContig(intervals);
            for (final List<Locatable> contigIntervals : byContig) {
                final List<SimpleInterval> contigSimple =
                        contigIntervals.stream()
                                .map(SimpleInterval::new)
                                .collect(Collectors.toList());
                final List<SimpleInterval> padded =
                        IntervalUtils.getIntervalsWithFlanks(contigSimple, padding, dict);
                for (final SimpleInterval iv : padded) {
                    System.out.printf(
                            "%s\t%d\t%d%n", iv.getContig(), iv.getStart(), iv.getEnd());
                }
            }
        } finally {
            refReader.close();
        }
    }

    /** Phase PRE.1: HC soft-clip policy before assembly (first step of finalizeRegion). */
    private static void readPreSoftclip(final String[] args) throws Exception {
        if (args.length < 1) {
            usage();
        }
        final String bamPath = args[0];
        final boolean dontUseSoftClippedBases =
                args.length > 1 && "1".equals(args[1]);
        final boolean overrideSoftclipFragmentCheck =
                args.length > 2 && "1".equals(args[2]);
        System.out.println(
                "qname\tflags\tfragment_length\tcigar_in\tcigar_out\tseq_len_in\tseq_len_out\taction\tos\toe");
        try (SamReader reader =
                SamReaderFactory.makeDefault().open(Paths.get(bamPath))) {
            final SAMFileHeader header = reader.getFileHeader();
            for (final SAMRecord rec : reader) {
                final GATKRead in = new SAMRecordToGATKReadAdapter(rec);
                final String cigarIn = in.getCigar().toString();
                final int lenIn = in.getLength();
                final boolean hardClip =
                        dontUseSoftClippedBases
                                || !(overrideSoftclipFragmentCheck
                                        || ReadUtils.hasWellDefinedFragmentSize(in));
                final GATKRead out;
                final String action;
                if (hardClip) {
                    out = ReadClipper.hardClipSoftClippedBases(in);
                    action = "hard_clip";
                } else {
                    final int softStart = in.getStart();
                    final int softEnd = in.getEnd();
                    out = ReadClipper.revertSoftClippedBases(in);
                    if (!out.isUnmapped()) {
                        out.setAttribute(
                                ReferenceConfidenceModel.ORIGINAL_SOFTCLIP_START_TAG, softStart);
                        out.setAttribute(
                                ReferenceConfidenceModel.ORIGINAL_SOFTCLIP_END_TAG, softEnd);
                    }
                    action = "revert";
                }
                final String os =
                        out.hasAttribute(ReferenceConfidenceModel.ORIGINAL_SOFTCLIP_START_TAG)
                                ? String.valueOf(
                                        out.getAttributeAsInteger(
                                                ReferenceConfidenceModel
                                                        .ORIGINAL_SOFTCLIP_START_TAG))
                                : "";
                final String oe =
                        out.hasAttribute(ReferenceConfidenceModel.ORIGINAL_SOFTCLIP_END_TAG)
                                ? String.valueOf(
                                        out.getAttributeAsInteger(
                                                ReferenceConfidenceModel
                                                        .ORIGINAL_SOFTCLIP_END_TAG))
                                : "";
                System.out.printf(
                        Locale.ROOT,
                        "%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s\t%s%n",
                        rec.getReadName(),
                        rec.getFlags(),
                        rec.getInferredInsertSize(),
                        cigarIn,
                        out.getCigar().toString(),
                        lenIn,
                        out.getLength(),
                        action,
                        os,
                        oe);
            }
        }
    }

    /** Phase PRE.2: {@code filterNonPassingReads} read-length gate ({@code unclippedReadLength >= 10}). */
    private static void readPreLen(final String[] args) throws Exception {
        if (args.length < 1) {
            usage();
        }
        final String bamPath = args[0];
        final int threshold = 10; // HaplotypeCallerEngine.READ_LENGTH_FILTER_THRESHOLD
        System.out.println("qname\tcigar\tread_length\tunclipped_length\tpasses_len_filter");
        try (SamReader reader =
                SamReaderFactory.makeDefault().open(Paths.get(bamPath))) {
            for (final SAMRecord rec : reader) {
                final GATKRead in = new SAMRecordToGATKReadAdapter(rec);
                final int unclipped = AlignmentUtils.unclippedReadLength(in);
                final boolean pass = unclipped >= threshold;
                System.out.printf(
                        Locale.ROOT,
                        "%s\t%s\t%d\t%d\t%b%n",
                        rec.getReadName(),
                        in.getCigar().toString(),
                        in.getLength(),
                        unclipped,
                        pass);
            }
        }
    }

    /** Phase PRE.3: {@code filterNonPassingReads} MQ gate ({@code mapq >= threshold}, default 20). */
    private static void readPreMq(final String[] args) throws Exception {
        if (args.length < 1) {
            usage();
        }
        final String bamPath = args[0];
        final int threshold =
                args.length > 1
                        ? Integer.parseInt(args[1])
                        : HaplotypeCallerEngine.DEFAULT_READ_QUALITY_FILTER_THRESHOLD;
        System.out.println("qname\tmapq\tmq_threshold\tpasses_mq_filter");
        try (SamReader reader =
                SamReaderFactory.makeDefault().open(Paths.get(bamPath))) {
            for (final SAMRecord rec : reader) {
                final GATKRead in = new SAMRecordToGATKReadAdapter(rec);
                final int mapq = in.getMappingQuality();
                final boolean pass = mapq >= threshold;
                System.out.printf(
                        Locale.ROOT,
                        "%s\t%d\t%d\t%b%n",
                        rec.getReadName(), mapq, threshold, pass);
            }
        }
    }

    /** Phase PRE.4: {@code cleanOverlappingReadPairs} / {@code FragmentUtils.adjustQualsOfOverlappingPairedFragments}. */
    private static void readPreOverlap(final String[] args) throws Exception {
        if (args.length < 1) {
            usage();
        }
        final String bamPath = args[0];
        System.out.println("qname\tqual_in\tqual_out\toverlap_pair");
        try (SamReader reader =
                SamReaderFactory.makeDefault().open(Paths.get(bamPath))) {
            final SAMFileHeader header = reader.getFileHeader();
            final List<GATKRead> reads = new ArrayList<>();
            final Map<String, String> qualIn = new LinkedHashMap<>();
            for (final SAMRecord rec : reader) {
                final GATKRead g = new SAMRecordToGATKReadAdapter(rec);
                reads.add(g);
                qualIn.put(g.getName(), formatQualCsv(g.getBaseQualities()));
            }
            reads.sort(new ReadCoordinateComparator(header));
            final SampleList samplesList = SampleList.singletonSampleList("s1");
            AssemblyBasedCallerUtils.cleanOverlappingReadPairs(
                    reads,
                    samplesList,
                    header,
                    true,
                    OptionalInt.empty(),
                    OptionalInt.empty());
            final Set<String> inPair = new LinkedHashSet<>();
            final FragmentCollection<GATKRead> fragments = FragmentCollection.create(reads);
            for (final org.apache.commons.lang3.tuple.Pair<GATKRead, GATKRead> p :
                    fragments.getOverlappingPairs()) {
                inPair.add(p.getLeft().getName());
                inPair.add(p.getRight().getName());
            }
            for (final GATKRead g : reads) {
                System.out.printf(
                        Locale.ROOT,
                        "%s\t%s\t%s\t%b%n",
                        g.getName(),
                        qualIn.get(g.getName()),
                        formatQualCsv(g.getBaseQualities()),
                        inPair.contains(g.getName()));
            }
        }
    }

    private static String formatQualCsv(final byte[] quals) {
        final StringBuilder sb = new StringBuilder();
        for (int i = 0; i < quals.length; i++) {
            if (i > 0) {
                sb.append(',');
            }
            sb.append(quals[i] & 0xff);
        }
        return sb.toString();
    }

    /** Phase D.4: IUPAC pre, HC filters, identity post (no downsampler). */
    private static void readShardPipeline(final String[] args) throws Exception {
        if (args.length < 1) {
            usage();
        }
        final String bamPath = args[0];
        System.out.println(
                "qname\tflags\tmapq\tseq_raw\tseq_after_pre\tpasses_hc_filter\tseq_after_post");
        try (SamReader reader =
                SamReaderFactory.makeDefault().open(Paths.get(bamPath))) {
            final SAMFileHeader samHeader = reader.getFileHeader();
            final ReadTransformer pre = HaplotypeCallerEngine.makeStandardHCReadTransformer();
            final ReadTransformer post = ReadTransformer.identity();
            final List<ReadFilter> filters = HaplotypeCallerEngine.makeStandardHCReadFilters();
            for (final ReadFilter f : filters) {
                f.setHeader(samHeader);
            }
            final ReadFilter combined = ReadFilter.fromList(filters, samHeader);
            for (final SAMRecord rec : reader) {
                final byte[] rawBases = rec.getReadBases().clone();
                final String seqRaw = new String(rawBases);
                final GATKRead read = new SAMRecordToGATKReadAdapter(rec);
                pre.apply(read);
                final String seqAfterPre = new String(read.getBases());
                final boolean pass = combined.test(read);
                post.apply(read);
                final String seqAfterPost = new String(read.getBases());
                System.out.printf(
                        Locale.ROOT,
                        "%s\t%d\t%d\t%s\t%s\t%s\t%s%n",
                        rec.getReadName(),
                        rec.getFlags(),
                        rec.getMappingQuality(),
                        seqRaw,
                        seqAfterPre,
                        pass ? "true" : "false",
                        seqAfterPost);
            }
        }
    }

    private static void readFilters(final String[] args) throws Exception {
        if (args.length < 1) {
            usage();
        }
        final String bamPath = args[0];
        System.out.println("qname\tflags\tmapq\tpasses_hc_filter");
        try (SamReader reader =
                SamReaderFactory.makeDefault().open(Paths.get(bamPath))) {
            final SAMFileHeader samHeader = reader.getFileHeader();
            final List<ReadFilter> filters = HaplotypeCallerEngine.makeStandardHCReadFilters();
            for (final ReadFilter f : filters) {
                f.setHeader(samHeader);
            }
            final List<CountingReadFilter> wrapped = new ArrayList<>();
            for (final ReadFilter f : filters) {
                wrapped.add(new CountingReadFilter(f));
            }
            CountingReadFilter combined = wrapped.get(0);
            for (int i = 1; i < wrapped.size(); i++) {
                combined = combined.and(wrapped.get(i));
            }
            for (final SAMRecord rec : reader) {
                final GATKRead gatkRead = new SAMRecordToGATKReadAdapter(rec);
                final boolean pass = combined.test(gatkRead);
                System.out.printf(
                        "%s\t%d\t%d\t%s%n",
                        rec.getReadName(),
                        rec.getFlags(),
                        rec.getMappingQuality(),
                        pass ? "true" : "false");
            }
            System.out.println(HC_READ_FILTER_COUNT_SECTION);
            System.out.println("filter\tfiltered_count");
            for (final CountingReadFilter w : wrapped) {
                System.out.printf("%s\t%d%n", w.getName(), w.getFilteredCount());
            }
        }
    }

    private static Iterator<AlignmentContext> makeLocusIterator(
            final MultiIntervalLocalReadShard shard, final HcContext ctx) {
        configureHcProductionReadShard(shard, ctx);
        final SAMFileHeader header = ctx.header;
        final ReadCachingIterator cache = new ReadCachingIterator(shard.iterator());
        final LocusIteratorByState libs =
                new LocusIteratorByState(
                        cache,
                        DownsamplingMethod.NONE,
                        ReadUtils.getSamplesFromHeader(header),
                        header,
                        true);
        final IntervalLocusIterator intervalIt =
                new IntervalLocusIterator(shard.getIntervals().iterator());
        return new IntervalAlignmentContextIterator(
                libs, intervalIt, header.getSequenceDictionary());
    }

    @SuppressWarnings("unchecked")
    private static List<ActivityProfileState> profileStates(final BandPassActivityProfile profile)
            throws Exception {
        final Method m = ActivityProfile.class.getDeclaredMethod("getStateList");
        m.setAccessible(true);
        return (List<ActivityProfileState>) m.invoke(profile);
    }

    private static void walkActivity(
            final HcContext ctx,
            final String intervalCli,
            final ActivityRowWriter writer)
            throws Exception {
        final List<SimpleInterval> intervals =
                parseIntervals(ctx.header.getSequenceDictionary(), intervalCli);
        for (final SimpleInterval interval : intervals) {
            final MultiIntervalLocalReadShard shard =
                    new MultiIntervalLocalReadShard(
                            Collections.singletonList(interval),
                            ctx.asmArgs.assemblyRegionPadding,
                            ctx.readsSource);
            final BandPassActivityProfile profile =
                    new BandPassActivityProfile(
                            ctx.asmArgs.maxProbPropagationDistance,
                            ctx.asmArgs.activeProbThreshold,
                            BandPassActivityProfile.MAX_FILTER_SIZE,
                            BandPassActivityProfile.DEFAULT_SIGMA,
                            ctx.header);
            final Iterator<AlignmentContext> locusIt = makeLocusIterator(shard, ctx);
            while (locusIt.hasNext()) {
                final AlignmentContext pileup = locusIt.next();
                final SimpleInterval pileupInterval = new SimpleInterval(pileup);
                final ReferenceContext refCtx =
                        new ReferenceContext(ctx.reference, pileupInterval);
                final FeatureContext featCtx = featureContextForLocus(ctx, pileupInterval);
                final ActivityProfileState raw =
                        ctx.engine.isActive(pileup, refCtx, featCtx);
                writer.onRaw(pileupInterval, raw);
                profile.add(raw);
            }
            final List<ActivityProfileState> states = profileStates(profile);
            if (states.isEmpty()) {
                continue;
            }
            final int start = states.get(0).getLoc().getStart();
            for (int i = 0; i < states.size(); i++) {
                final ActivityProfileState st = states.get(i);
                final int pos = start + i;
                if (pos < interval.getStart() || pos > interval.getEnd()) {
                    continue;
                }
                writer.onSmoothed(interval.getContig(), pos, st);
            }
        }
    }

    private interface ActivityRowWriter {
        void onRaw(SimpleInterval loc, ActivityProfileState raw) throws Exception;

        void onSmoothed(String contig, int pos, ActivityProfileState st) throws Exception;
    }

    private static void rawActivity(final String[] args) throws Exception {
        if (args.length < 3) {
            usage();
        }
        final String refPath = args[0];
        final String bamPath = args[1];
        final String intervalCli = args[2];
        final int padding = args.length > 3 ? parsePadding(args[3]) : DEFAULT_PADDING;
        System.out.println("contig\tpos\tactive_prob\toriginal_active_prob\tkind");
        try (HcContext ctx = new HcContext(refPath, bamPath, padding)) {
            walkActivity(
                    ctx,
                    intervalCli,
                    new ActivityRowWriter() {
                        @Override
                        public void onRaw(final SimpleInterval loc, final ActivityProfileState raw) {
                            System.out.printf(
                                    "%s\t%d\t%s\t%s\t%s%n",
                                    loc.getContig(),
                                    loc.getStart(),
                                    formatProb(raw.isActiveProb()),
                                    formatProb(raw.isActiveProb()),
                                    formatKind(raw.getResultState()));
                        }

                        @Override
                        public void onSmoothed(
                                final String contig, final int pos, final ActivityProfileState st) {}
                    });
        }
    }

    /** Tier B C.5.3: {@code --force-calling-alleles-present} with a real {@link FeatureContext}. */
    private static void rawActivityForce(final String[] args) throws Exception {
        if (args.length < 4) {
            usage();
        }
        final String refPath = args[0];
        final String bamPath = args[1];
        final String intervalCli = args[2];
        final String vcfPath = args[3];
        final int padding = args.length > 4 ? parsePadding(args[4]) : DEFAULT_PADDING;
        System.out.println("contig\tpos\tactive_prob\toriginal_active_prob\tkind");
        try (HcContext ctx = new HcContext(refPath, bamPath, padding, vcfPath)) {
            walkActivity(
                    ctx,
                    intervalCli,
                    new ActivityRowWriter() {
                        @Override
                        public void onRaw(final SimpleInterval loc, final ActivityProfileState raw) {
                            System.out.printf(
                                    "%s\t%d\t%s\t%s\t%s%n",
                                    loc.getContig(),
                                    loc.getStart(),
                                    formatProb(raw.isActiveProb()),
                                    formatProb(raw.isActiveProb()),
                                    formatKind(raw.getResultState()));
                        }

                        @Override
                        public void onSmoothed(
                                final String contig, final int pos, final ActivityProfileState st) {}
                    });
        }
    }

    private static void smoothedActivity(final String[] args) throws Exception {
        if (args.length < 3) {
            usage();
        }
        final String refPath = args[0];
        final String bamPath = args[1];
        final String intervalCli = args[2];
        final int padding = args.length > 3 ? parsePadding(args[3]) : DEFAULT_PADDING;
        System.out.println("contig\tpos\tsmoothed_active_prob");
        try (HcContext ctx = new HcContext(refPath, bamPath, padding)) {
            walkActivity(
                    ctx,
                    intervalCli,
                    new ActivityRowWriter() {
                        @Override
                        public void onRaw(final SimpleInterval loc, final ActivityProfileState raw) {}

                        @Override
                        public void onSmoothed(
                                final String contig, final int pos, final ActivityProfileState st) {
                            System.out.printf(
                                    "%s\t%d\t%s%n",
                                    contig, pos, formatProb(st.isActiveProb()));
                        }
                    });
        }
    }

    private static final class RefConfidenceEngines {
        final ReferenceConfidenceModel rcm;
        final MinimalGenotypingEngine genoEngine;
        final int ploidy;
        final byte minBaseQuality;

        RefConfidenceEngines(final HcContext ctx) throws Exception {
            final Field rcmField =
                    HaplotypeCallerEngine.class.getDeclaredField("referenceConfidenceModel");
            rcmField.setAccessible(true);
            rcm = (ReferenceConfidenceModel) rcmField.get(ctx.engine);
            final Field genoField =
                    HaplotypeCallerEngine.class.getDeclaredField(
                            "activeRegionEvaluationGenotyperEngine");
            genoField.setAccessible(true);
            genoEngine = (MinimalGenotypingEngine) genoField.get(ctx.engine);
            ploidy = genoEngine.getConfiguration().genotypeArgs.samplePloidy;
            final Field hcArgsField = HaplotypeCallerEngine.class.getDeclaredField("hcArgs");
            hcArgsField.setAccessible(true);
            final HaplotypeCallerArgumentCollection hcArgs =
                    (HaplotypeCallerArgumentCollection) hcArgsField.get(ctx.engine);
            minBaseQuality = (byte) hcArgs.minBaseQualityScore;
        }
    }

    private static final class RefConfidenceLocusRow {
        final double[] gl;
        final int gq;
        final int dp;
        final double activeProb;

        RefConfidenceLocusRow(
                final double[] gl, final int gq, final int dp, final double activeProb) {
            this.gl = gl;
            this.gq = gq;
            this.dp = dp;
            this.activeProb = activeProb;
        }
    }

    /** Rust {@code reference_gq_from_log10_gl} (max vs second-best, phred-scaled). */
    private static int referenceGqFromLog10Gl(final double[] gl) {
        if (gl.length == 0) {
            return 0;
        }
        int bestIdx = 0;
        for (int i = 1; i < gl.length; i++) {
            if (gl[i] > gl[bestIdx]) {
                bestIdx = i;
            }
        }
        double second = Double.NEGATIVE_INFINITY;
        for (int i = 0; i < gl.length; i++) {
            if (i != bestIdx && gl[i] > second) {
                second = gl[i];
            }
        }
        final double best = gl[bestIdx];
        if (!Double.isFinite(best) || !Double.isFinite(second)) {
            return 0;
        }
        return (int)
                Math.round(Math.max(0.0, Math.min(99.0, -10.0 * (second - best))));
    }

    private static RefConfidenceLocusRow refConfidenceAtPileup(
            final RefConfidenceEngines engines,
            final AlignmentContext pileup,
            final ReferenceContext refCtx)
            throws Exception {
        final MathUtils.RunningAverage hqAvg = new MathUtils.RunningAverage();
        final RefVsAnyResult refVsAny =
                (RefVsAnyResult)
                        engines.rcm.calcGenotypeLikelihoodsOfRefVsAny(
                                engines.ploidy,
                                pileup.getBasePileup(),
                                refCtx.getBase(),
                                engines.minBaseQuality,
                                hqAvg,
                                false);
        final double[] gl = cappedGenotypeLikelihoods(refVsAny);
        final double activeProb =
                engines.genoEngine.calculateSingleSampleRefVsAnyActiveStateProfileValue(gl);
        final int gq = referenceGqFromLog10Gl(gl);
        final int dp = pileup.getBasePileup().size();
        return new RefConfidenceLocusRow(gl, gq, dp, activeProb);
    }

    /** Phase H.1: per-locus RCM GL + GQ + DP + activity (Rust ref_confidence_dump). */
    private static void referenceConfidenceLocus(final String[] args) throws Exception {
        if (args.length < 3) {
            usage();
        }
        final String refPath = args[0];
        final String bamPath = args[1];
        final String intervalCli = args[2];
        final int padding = args.length > 3 ? parsePadding(args[3]) : DEFAULT_PADDING;
        System.out.println("contig\tpos\tgl0\tgl1\tgl2\tgq\tdp\tactive_prob");
        try (HcContext ctx = new HcContext(refPath, bamPath, padding)) {
            final RefConfidenceEngines engines = new RefConfidenceEngines(ctx);
            final List<SimpleInterval> intervals =
                    parseIntervals(ctx.header.getSequenceDictionary(), intervalCli);
            for (final SimpleInterval interval : intervals) {
                final MultiIntervalLocalReadShard shard =
                        new MultiIntervalLocalReadShard(
                                Collections.singletonList(interval),
                                ctx.asmArgs.assemblyRegionPadding,
                                ctx.readsSource);
                final Iterator<AlignmentContext> locusIt = makeLocusIterator(shard, ctx);
                while (locusIt.hasNext()) {
                    final AlignmentContext pileup = locusIt.next();
                    final SimpleInterval loc = new SimpleInterval(pileup);
                    final ReferenceContext refCtx =
                            new ReferenceContext(ctx.reference, loc);
                    final RefConfidenceLocusRow row =
                            refConfidenceAtPileup(engines, pileup, refCtx);
                    System.out.printf(
                            Locale.ROOT,
                            "%s\t%d\t%s\t%s\t%s\t%d\t%d\t%s%n",
                            loc.getContig(),
                            loc.getStart(),
                            formatProb(row.gl[0]),
                            row.gl.length > 1 ? formatProb(row.gl[1]) : formatProb(0),
                            row.gl.length > 2 ? formatProb(row.gl[2]) : formatProb(0),
                            row.gq,
                            row.dp,
                            formatProb(row.activeProb));
                }
            }
        }
    }

    /** Phase H.1.2: first inactive region referenceModelForNoVariation summary. */
    private static void inactiveReferenceModel(final String[] args) throws Exception {
        if (args.length < 3) {
            usage();
        }
        final String refPath = args[0];
        final String bamPath = args[1];
        final String intervalCli = args[2];
        final int padding = args.length > 3 ? parsePadding(args[3]) : DEFAULT_PADDING;
        try (HcContext ctx = new HcContext(refPath, bamPath, padding)) {
            final RefConfidenceEngines engines = new RefConfidenceEngines(ctx);
            final List<SimpleInterval> intervals =
                    parseIntervals(ctx.header.getSequenceDictionary(), intervalCli);
            final List<List<Locatable>> byContig = groupByContig(intervals);
            AssemblyRegion inactive = null;
            for (final List<Locatable> contigIntervals : byContig) {
                final List<SimpleInterval> contigSimple =
                        contigIntervals.stream()
                                .map(SimpleInterval::new)
                                .collect(Collectors.toList());
                final MultiIntervalLocalReadShard shard =
                        new MultiIntervalLocalReadShard(
                                contigSimple, padding, ctx.readsSource);
                configureHcProductionReadShard(shard, ctx);
                final AssemblyRegionIterator iter =
                        new AssemblyRegionIterator(
                                shard,
                                ctx.header,
                                ctx.reference,
                                null,
                                ctx.engine,
                                ctx.asmArgs,
                                true);
                while (iter.hasNext()) {
                    final AssemblyRegion r = iter.next();
                    if (!r.isActive()) {
                        inactive = r;
                        break;
                    }
                }
                if (inactive != null) {
                    break;
                }
            }
            if (inactive == null) {
                throw new IllegalStateException("no inactive assembly region in interval");
            }
            final List<AlignmentAndReferenceContext> pileups =
                    new ArrayList<>(inactive.getAlignmentData());
            pileups.sort(
                    Comparator.comparingInt(
                            p -> p.getAlignmentContext().getStart()));
            final List<RefConfidenceLocusRow> loci = new ArrayList<>();
            for (final AlignmentAndReferenceContext ar : pileups) {
                loci.add(
                        refConfidenceAtPileup(
                                engines, ar.getAlignmentContext(), ar.getReferenceContext()));
            }
            System.out.println("region_contig\t" + inactive.getContig());
            System.out.println("region_start\t" + inactive.getStart());
            System.out.println("region_end\t" + inactive.getEnd());
            System.out.println("is_active\tfalse");
            System.out.println("path\treferenceModelForNoVariation");
            System.out.println("emit_mode\tGVCF");
            System.out.println("locus_count\t" + loci.size());
            System.out.println("reference_blocks_emitted\t" + loci.size());
            System.out.println("reference_sites_emitted\t0");
            if (!loci.isEmpty()) {
                System.out.println("first_gq\t" + loci.get(0).gq);
                System.out.println("first_dp\t" + loci.get(0).dp);
                final RefConfidenceLocusRow last = loci.get(loci.size() - 1);
                System.out.println("last_gq\t" + last.gq);
                System.out.println("last_dp\t" + last.dp);
            }
            emitInactiveGvcfBlocksFromPileups(pileups, loci);
        }
    }

    /** L5.1 inactive dump extension: per-locus + gVCF block rows (Rust ref_confidence_dump parity). */
    private static void emitInactiveGvcfBlocksFromPileups(
            final List<AlignmentAndReferenceContext> pileups,
            final List<RefConfidenceLocusRow> loci) {
        if (loci.isEmpty()) {
            return;
        }
        final int[] gqBands = {
            1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23,
            24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44,
            45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 70, 80, 90, 99
        };
        final List<Integer> positions = new ArrayList<>(pileups.size());
        for (final AlignmentAndReferenceContext ar : pileups) {
            positions.add(ar.getAlignmentContext().getStart());
        }
        final List<InactiveGvcfBlockRow> blocks = composeInactiveGvcfBlocks(positions, loci, gqBands);
        System.out.println("gvcf_block_count\t" + blocks.size());
        System.out.println("locus_header\tpos\tgq\tdp");
        for (int i = 0; i < loci.size(); i++) {
            final RefConfidenceLocusRow row = loci.get(i);
            System.out.println(
                    "locus\t" + positions.get(i) + "\t" + row.gq + "\t" + row.dp);
        }
        System.out.println("block_header\tstart\tend\tmin_dp\tmax_dp\tgq_band\tmin_rgq");
        for (final InactiveGvcfBlockRow block : blocks) {
            System.out.println(
                    "block\t"
                            + block.start
                            + "\t"
                            + block.end
                            + "\t"
                            + block.minDp
                            + "\t"
                            + block.maxDp
                            + "\t"
                            + block.gqBandUpper
                            + "\t"
                            + block.minRgq);
        }
    }

    private static final class InactiveGvcfBlockRow {
        final int start;
        final int end;
        final int minDp;
        final int maxDp;
        final int gqBandUpper;
        final int minRgq;

        InactiveGvcfBlockRow(
                final int start,
                final int end,
                final int minDp,
                final int maxDp,
                final int gqBandUpper,
                final int minRgq) {
            this.start = start;
            this.end = end;
            this.minDp = minDp;
            this.maxDp = maxDp;
            this.gqBandUpper = gqBandUpper;
            this.minRgq = minRgq;
        }
    }

    private static List<InactiveGvcfBlockRow> composeInactiveGvcfBlocks(
            final List<Integer> positions,
            final List<RefConfidenceLocusRow> loci,
            final int[] gqBands) {
        final List<InactiveGvcfBlockRow> out = new ArrayList<>();
        int blockStart = positions.get(0);
        int blockEnd = positions.get(0);
        int minDp = loci.get(0).dp;
        int maxDp = loci.get(0).dp;
        int minRgq = loci.get(0).gq;
        int partLo = gvcfPartitionLower(loci.get(0).gq, gqBands);
        int partHi = gvcfPartitionUpper(loci.get(0).gq, gqBands);
        for (int i = 1; i < loci.size(); i++) {
            final int pos = positions.get(i);
            final RefConfidenceLocusRow row = loci.get(i);
            final boolean contiguous = pos == blockEnd + 1;
            final boolean samePart =
                    row.gq >= partLo
                            && row.gq < partHi
                            && gvcfPartitionLower(row.gq, gqBands) == partLo
                            && gvcfPartitionUpper(row.gq, gqBands) == partHi;
            if (contiguous && samePart) {
                blockEnd = pos;
                minDp = Math.min(minDp, row.dp);
                maxDp = Math.max(maxDp, row.dp);
                minRgq = Math.min(minRgq, row.gq);
            } else {
                out.add(
                        new InactiveGvcfBlockRow(
                                blockStart, blockEnd, minDp, maxDp, partHi, minRgq));
                blockStart = pos;
                blockEnd = pos;
                minDp = row.dp;
                maxDp = row.dp;
                minRgq = row.gq;
                partLo = gvcfPartitionLower(row.gq, gqBands);
                partHi = gvcfPartitionUpper(row.gq, gqBands);
            }
        }
        out.add(new InactiveGvcfBlockRow(blockStart, blockEnd, minDp, maxDp, partHi, minRgq));
        return out;
    }

    private static int gvcfPartitionUpper(final int gq, final int[] gqBands) {
        final int capped = Math.min(Math.max(gq, 0), 99);
        int lower = 0;
        for (final int upper : gqBands) {
            if (capped < upper) {
                return upper;
            }
            lower = upper;
        }
        return 100;
    }

    private static int gvcfPartitionLower(final int gq, final int[] gqBands) {
        final int capped = Math.min(Math.max(gq, 0), 99);
        int lower = 0;
        for (final int upper : gqBands) {
            if (capped < upper) {
                return lower;
            }
            lower = upper;
        }
        return lower;
    }

    /** Phase C.4: per-locus log10 GL vector + activity (MinimalGenotypingEngine path). */
    private static void genotypeLikelihoodActivity(final String[] args) throws Exception {
        if (args.length < 3) {
            usage();
        }
        final String refPath = args[0];
        final String bamPath = args[1];
        final String intervalCli = args[2];
        final int padding = args.length > 3 ? parsePadding(args[3]) : DEFAULT_PADDING;
        System.out.println("contig\tpos\tgl0\tgl1\tgl2\tactive_prob\toriginal_active_prob");
        try (HcContext ctx = new HcContext(refPath, bamPath, padding)) {
            walkGenotypeLikelihood(ctx, intervalCli);
        }
    }

    private static void walkGenotypeLikelihood(
            final HcContext ctx, final String intervalCli) throws Exception {
        final Field rcmField = HaplotypeCallerEngine.class.getDeclaredField("referenceConfidenceModel");
        rcmField.setAccessible(true);
        final ReferenceConfidenceModel rcm = (ReferenceConfidenceModel) rcmField.get(ctx.engine);
        final Field genoField =
                HaplotypeCallerEngine.class.getDeclaredField("activeRegionEvaluationGenotyperEngine");
        genoField.setAccessible(true);
        final MinimalGenotypingEngine genoEngine = (MinimalGenotypingEngine) genoField.get(ctx.engine);
        final int ploidy = genoEngine.getConfiguration().genotypeArgs.samplePloidy;
        final Field hcArgsField = HaplotypeCallerEngine.class.getDeclaredField("hcArgs");
        hcArgsField.setAccessible(true);
        final HaplotypeCallerArgumentCollection hcArgs =
                (HaplotypeCallerArgumentCollection) hcArgsField.get(ctx.engine);
        final int minBaseQuality = hcArgs.minBaseQualityScore;

        final List<SimpleInterval> intervals =
                parseIntervals(ctx.header.getSequenceDictionary(), intervalCli);
        for (final SimpleInterval interval : intervals) {
            final MultiIntervalLocalReadShard shard =
                    new MultiIntervalLocalReadShard(
                            Collections.singletonList(interval),
                            ctx.asmArgs.assemblyRegionPadding,
                            ctx.readsSource);
            final Iterator<AlignmentContext> locusIt = makeLocusIterator(shard, ctx);
            while (locusIt.hasNext()) {
                final AlignmentContext pileup = locusIt.next();
                final SimpleInterval loc = new SimpleInterval(pileup);
                final ReferenceContext refCtx =
                        new ReferenceContext(ctx.reference, loc);
                final MathUtils.RunningAverage hqAvg = new MathUtils.RunningAverage();
                final RefVsAnyResult refVsAny =
                        (RefVsAnyResult)
                                rcm.calcGenotypeLikelihoodsOfRefVsAny(
                                        ploidy,
                                        pileup.getBasePileup(),
                                        refCtx.getBase(),
                                        (byte) minBaseQuality,
                                        hqAvg,
                                        false);
                final double[] gl = cappedGenotypeLikelihoods(refVsAny);
                final double activeProb =
                        genoEngine.calculateSingleSampleRefVsAnyActiveStateProfileValue(gl);
                double maxGl = gl[0];
                for (int i = 1; i < gl.length; i++) {
                    if (gl[i] > maxGl) {
                        maxGl = gl[i];
                    }
                }
                final double original = maxGl - gl[0];
                System.out.printf(
                        "%s\t%d\t%s\t%s\t%s\t%s\t%s%n",
                        loc.getContig(),
                        loc.getStart(),
                        formatProb(gl[0]),
                        gl.length > 1 ? formatProb(gl[1]) : formatProb(0),
                        gl.length > 2 ? formatProb(gl[2]) : formatProb(0),
                        formatProb(activeProb),
                        formatProb(original));
            }
        }
    }

    private static void activeLocus(final String[] args) throws Exception {
        if (args.length < 3) {
            usage();
        }
        final String refPath = args[0];
        final String bamPath = args[1];
        final String intervalCli = args[2];
        final int padding = args.length > 3 ? parsePadding(args[3]) : DEFAULT_PADDING;
        final double threshold;
        try (HcContext ctx = new HcContext(refPath, bamPath, padding)) {
            threshold = ctx.asmArgs.activeProbThreshold;
            System.out.println("contig\tpos\tis_active");
            walkActivity(
                    ctx,
                    intervalCli,
                    new ActivityRowWriter() {
                        @Override
                        public void onRaw(final SimpleInterval loc, final ActivityProfileState raw) {}

                        @Override
                        public void onSmoothed(
                                final String contig, final int pos, final ActivityProfileState st) {
                            final boolean active = st.isActiveProb() > threshold;
                            System.out.printf(
                                    "%s\t%d\t%s%n", contig, pos, active ? "true" : "false");
                        }
                    });
        }
    }

    private static void locusPileup(final String[] args) throws Exception {
        if (args.length < 3) {
            usage();
        }
        final String refPath = args[0];
        final String bamPath = args[1];
        final String intervalCli = args[2];
        final int padding = args.length > 3 ? parsePadding(args[3]) : DEFAULT_PADDING;
        System.out.println("contig\tpos\tpileup_depth");
        try (HcContext ctx = new HcContext(refPath, bamPath, padding)) {
            final List<SimpleInterval> intervals =
                    parseIntervals(ctx.header.getSequenceDictionary(), intervalCli);
            final List<List<Locatable>> byContig = groupByContig(intervals);
            for (final List<Locatable> contigIntervals : byContig) {
                final List<SimpleInterval> contigSimple =
                        contigIntervals.stream()
                                .map(SimpleInterval::new)
                                .collect(Collectors.toList());
                final MultiIntervalLocalReadShard shard =
                        new MultiIntervalLocalReadShard(
                                contigSimple, padding, ctx.readsSource);
                final Iterator<AlignmentContext> locusIt = makeLocusIterator(shard, ctx);
                while (locusIt.hasNext()) {
                    final AlignmentContext pileup = locusIt.next();
                    final SimpleInterval loc = new SimpleInterval(pileup);
                    System.out.printf(
                            "%s\t%d\t%d%n",
                            loc.getContig(), loc.getStart(), pileup.size());
                }
            }
        }
    }

    private static void assemblyRegionReads(final String[] args) throws Exception {
        if (args.length < 3) {
            usage();
        }
        final String refPath = args[0];
        final String bamPath = args[1];
        final String intervalCli = args[2];
        final int padding = args.length > 3 ? parsePadding(args[3]) : DEFAULT_PADDING;
        System.out.println("contig\tstart\tend\tis_active\tread_count");
        try (HcContext ctx = new HcContext(refPath, bamPath, padding)) {
            final List<SimpleInterval> intervals =
                    parseIntervals(ctx.header.getSequenceDictionary(), intervalCli);
            final List<List<Locatable>> byContig = groupByContig(intervals);
            for (final List<Locatable> contigIntervals : byContig) {
                final List<SimpleInterval> contigSimple =
                        contigIntervals.stream()
                                .map(SimpleInterval::new)
                                .collect(Collectors.toList());
                final MultiIntervalLocalReadShard shard =
                        new MultiIntervalLocalReadShard(
                                contigSimple, padding, ctx.readsSource);
                configureHcProductionReadShard(shard, ctx);
                final AssemblyRegionIterator iter =
                        new AssemblyRegionIterator(
                                shard,
                                ctx.header,
                                ctx.reference,
                                null,
                                ctx.engine,
                                ctx.asmArgs,
                                false);
                while (iter.hasNext()) {
                    final AssemblyRegion r = iter.next();
                    System.out.printf(
                            "%s\t%d\t%d\t%s\t%d%n",
                            r.getContig(),
                            r.getStart(),
                            r.getEnd(),
                            r.isActive() ? "true" : "false",
                            r.getReads().size());
                    final List<org.broadinstitute.hellbender.utils.read.GATKRead> reads =
                            new ArrayList<>(r.getReads());
                    reads.sort(
                            Comparator.comparing(
                                            (org.broadinstitute.hellbender.utils.read.GATKRead read) ->
                                                    read.getName())
                                    .thenComparingInt(
                                            org.broadinstitute.hellbender.utils.read.GATKRead
                                                    ::getStart));
                    for (final org.broadinstitute.hellbender.utils.read.GATKRead read : reads) {
                        System.out.printf(
                                "read\t%s\t%d\t%d%n",
                                read.getName(), read.getStart(), read.getEnd());
                    }
                }
            }
        }
    }

    private static void assemblyRegionReference(final String[] args) throws Exception {
        if (args.length < 3) {
            usage();
        }
        final String refPath = args[0];
        final String bamPath = args[1];
        final String intervalCli = args[2];
        final int padding = args.length > 3 ? parsePadding(args[3]) : DEFAULT_PADDING;
        System.out.println(
                "contig\tstart\tend\tis_active\textended_start\textended_end\tref_window_start\tref_window_end\tref_len");
        try (HcContext ctx = new HcContext(refPath, bamPath, padding)) {
            final List<SimpleInterval> intervals =
                    parseIntervals(ctx.header.getSequenceDictionary(), intervalCli);
            final List<List<Locatable>> byContig = groupByContig(intervals);
            for (final List<Locatable> contigIntervals : byContig) {
                final List<SimpleInterval> contigSimple =
                        contigIntervals.stream()
                                .map(SimpleInterval::new)
                                .collect(Collectors.toList());
                final MultiIntervalLocalReadShard shard =
                        new MultiIntervalLocalReadShard(
                                contigSimple, padding, ctx.readsSource);
                configureHcProductionReadShard(shard, ctx);
                final AssemblyRegionIterator iter =
                        new AssemblyRegionIterator(
                                shard,
                                ctx.header,
                                ctx.reference,
                                null,
                                ctx.engine,
                                ctx.asmArgs,
                                false);
                while (iter.hasNext()) {
                    final AssemblyRegion r = iter.next();
                    final SimpleInterval padded = r.getPaddedSpan();
                    final ReferenceContext refCtx =
                            new ReferenceContext(ctx.reference, padded);
                    final byte[] bases = refCtx.getBases();
                    final String refBases =
                            new String(bases, java.nio.charset.StandardCharsets.US_ASCII);
                    System.out.printf(
                            "%s\t%d\t%d\t%s\t%d\t%d\t%d\t%d\t%d%n",
                            r.getContig(),
                            r.getStart(),
                            r.getEnd(),
                            r.isActive() ? "true" : "false",
                            padded.getStart(),
                            padded.getEnd(),
                            refCtx.getStart(),
                            refCtx.getEnd(),
                            bases.length);
                    System.out.printf("ref_bases\t%s%n", refBases);
                }
            }
        }
    }

    private static boolean parseRegionTargetInactive(final String[] args, final int startIdx) {
        if (args.length <= startIdx) {
            return false;
        }
        return "inactive".equalsIgnoreCase(args[startIdx]);
    }

    private static int parsePaddingAndTargetIndex(
            final String[] args, final int startIdx, final boolean[] wantInactiveOut) {
        wantInactiveOut[0] = false;
        if (args.length <= startIdx) {
            return DEFAULT_PADDING;
        }
        if ("inactive".equalsIgnoreCase(args[startIdx])) {
            wantInactiveOut[0] = true;
            return DEFAULT_PADDING;
        }
        if ("active".equalsIgnoreCase(args[startIdx])) {
            return DEFAULT_PADDING;
        }
        if ("-".equals(args[startIdx])) {
            if (args.length > startIdx + 1) {
                wantInactiveOut[0] = parseRegionTargetInactive(args, startIdx + 1);
            }
            return DEFAULT_PADDING;
        }
        final int padding = parsePadding(args[startIdx]);
        if (args.length > startIdx + 1) {
            wantInactiveOut[0] = parseRegionTargetInactive(args, startIdx + 1);
        }
        return padding;
    }

    /**
     * ASM-1 materialize path: iterator reads (no {@code finalizeRegion}), {@code runLocalAssembly} only.
     */
    private static void emitAssemblyRegionHaplotypesFromMaterial(
            final RegionAssemblyMaterial material,
            final HcContext ctx,
            final int padding,
            final HaplotypeCallerArgumentCollection hcArgs,
            final ReadThreadingAssembler assembler,
            final SmithWatermanAligner aligner,
            final CachingIndexedFastaSequenceFile refReader)
            throws Exception {
        final SimpleInterval activeSpan =
                new SimpleInterval(material.contig, material.start, material.end);
        final AssemblyRegion region = new AssemblyRegion(activeSpan, padding, ctx.header);
        final String readGroupId = ctx.header.getReadGroups().get(0).getReadGroupId();
        for (int r = 0; r < material.readBases.size(); r++) {
            final GATKRead read =
                    ArtificialReadUtils.createArtificialRead(
                            ctx.header,
                            "mat_" + r,
                            material.contig,
                            material.start,
                            material.readBases.get(r),
                            material.readQuals.get(r));
            read.setReadGroup(readGroupId);
            region.add(read);
        }
        final SimpleInterval paddedReferenceLoc =
                AssemblyBasedCallerUtils.getPaddedReferenceLoc(
                        region,
                        AssemblyBasedCallerUtils.REFERENCE_PADDING_FOR_ASSEMBLY,
                        refReader);
        final Haplotype refHaplotype =
                AssemblyBasedCallerUtils.createReferenceHaplotype(
                        region, paddedReferenceLoc, refReader);
        final AssemblyResultSet ars =
                assembler.runLocalAssembly(
                        region,
                        refHaplotype,
                        material.refBases,
                        paddedReferenceLoc,
                        null,
                        ctx.header,
                        aligner,
                        null,
                        hcArgs.getDanglingEndSWParameters(),
                        hcArgs.getHaplotypeToReferenceSWParameters());
        int kmerSize = 0;
        try {
            kmerSize = ars.getMinimumKmerSize();
        } catch (final IllegalStateException ignored) {
            kmerSize = 0;
        }
        final String status =
                ars.isVariationPresent()
                        ? "assembled_some_variation"
                        : "just_assembled_reference";
        System.out.println("region_contig\t" + material.contig);
        System.out.println("region_start\t" + material.start);
        System.out.println("region_end\t" + material.end);
        System.out.println("is_active\ttrue");
        System.out.println("status\t" + status);
        System.out.println("kmer_size\t" + kmerSize);
        System.out.println("rank\tsequence\tscore\tis_reference\tcigar");
        final List<Haplotype> haps = new ArrayList<>(ars.getHaplotypeList());
        haps.sort(
                Comparator.comparingDouble((Haplotype h) -> 0.0)
                        .reversed()
                        .thenComparing(
                                h -> new String(h.getBases()),
                                Comparator.reverseOrder()));
        final byte[] refBases = ars.getReferenceHaplotype().getBases();
        final SWParameters hapSw =
                org.broadinstitute.hellbender.utils.smithwaterman
                        .SmithWatermanAlignmentConstants.NEW_SW_PARAMETERS;
        int rank = 0;
        boolean wroteRef = false;
        for (final Haplotype h : haps) {
            final htsjdk.samtools.Cigar cigar =
                    CigarUtils.calculateCigar(
                            refBases,
                            h.getBases(),
                            aligner,
                            hapSw,
                            org.broadinstitute.gatk.nativebindings.smithwaterman.SWOverhangStrategy
                                    .SOFTCLIP);
            if (cigar == null) {
                continue;
            }
            if (Arrays.equals(h.getBases(), refBases)) {
                wroteRef = true;
            }
            System.out.printf(
                    "%d\t%s\t0\t%s\t%s%n",
                    rank++,
                    new String(h.getBases(), StandardCharsets.US_ASCII),
                    h.isReference(),
                    cigar);
        }
        if (!wroteRef) {
            final htsjdk.samtools.Cigar refCigar =
                    CigarUtils.calculateCigar(
                            refBases,
                            refBases,
                            aligner,
                            hapSw,
                            org.broadinstitute.gatk.nativebindings.smithwaterman.SWOverhangStrategy
                                    .SOFTCLIP);
            System.out.printf(
                    "%d\t%s\t0\ttrue\t%s%n",
                    rank,
                    new String(refBases, StandardCharsets.US_ASCII),
                    refCigar == null ? "" : refCigar.toString());
        }
    }

    private static void assemblyRegionHaplotypes(final String[] args) throws Exception {
        if (args.length < 3) {
            usage();
        }
        final String refPath = args[0];
        final String bamPath = args[1];
        final String intervalCli = args[2];
        final boolean[] wantInactive = {false};
        final int padding = parsePaddingAndTargetIndex(args, 3, wantInactive);
        final boolean dumpInactive = wantInactive[0];
        if (dumpInactive) {
            try (HcContext ctx = new HcContext(refPath, bamPath, padding)) {
                final List<SimpleInterval> intervals =
                        parseIntervals(ctx.header.getSequenceDictionary(), intervalCli);
                for (final List<Locatable> contigIntervals : groupByContig(intervals)) {
                    final List<SimpleInterval> contigSimple =
                            contigIntervals.stream()
                                    .map(SimpleInterval::new)
                                    .collect(Collectors.toList());
                    final MultiIntervalLocalReadShard shard =
                            new MultiIntervalLocalReadShard(
                                    contigSimple, padding, ctx.readsSource);
                    configureHcProductionReadShard(shard, ctx);
                    final AssemblyRegionIterator iter =
                            new AssemblyRegionIterator(
                                    shard,
                                    ctx.header,
                                    ctx.reference,
                                    null,
                                    ctx.engine,
                                    ctx.asmArgs,
                                    false);
                    while (iter.hasNext()) {
                        final AssemblyRegion r = iter.next();
                        if (r.isActive()) {
                            continue;
                        }
                        System.out.println("region_contig\t" + r.getContig());
                        System.out.println("region_start\t" + r.getStart());
                        System.out.println("region_end\t" + r.getEnd());
                        System.out.println("is_active\tfalse");
                        System.out.println("status\tinactive_skip");
                        System.out.println("kmer_size\t0");
                        return;
                    }
                }
            }
            throw new IllegalArgumentException("no inactive assembly region in interval");
        }
        final RegionAssemblyMaterial material =
                materializeFirstActiveRegion(refPath, bamPath, intervalCli, padding);
        try (HcContext ctx = new HcContext(refPath, bamPath, padding);
                CachingIndexedFastaSequenceFile refReader =
                        new CachingIndexedFastaSequenceFile(Paths.get(refPath))) {
            final Field hcArgsField = HaplotypeCallerEngine.class.getDeclaredField("hcArgs");
            hcArgsField.setAccessible(true);
            final HaplotypeCallerArgumentCollection hcArgs =
                    (HaplotypeCallerArgumentCollection) hcArgsField.get(ctx.engine);
            final ReadThreadingAssembler assembler = productionAssemblerFromContext(ctx);
            final Field alignerField = HaplotypeCallerEngine.class.getDeclaredField("aligner");
            alignerField.setAccessible(true);
            final SmithWatermanAligner aligner =
                    (SmithWatermanAligner) alignerField.get(ctx.engine);
            emitAssemblyRegionHaplotypesFromMaterial(
                    material, ctx, padding, hcArgs, assembler, aligner, refReader);
        }
    }

    private static void assemblyRegionPairhmmLikelihoods(final String[] args) throws Exception {
        if (args.length < 3) {
            usage();
        }
        final String refPath = args[0];
        final String bamPath = args[1];
        final String intervalCli = args[2];
        final boolean[] wantInactive = {false};
        final int padding = parsePaddingAndTargetIndex(args, 3, wantInactive);
        final boolean dumpInactive = wantInactive[0];
        final byte bqThreshold = 18;
        try (HcContext ctx = new HcContext(refPath, bamPath, padding);
                CachingIndexedFastaSequenceFile refReader =
                        new CachingIndexedFastaSequenceFile(Paths.get(refPath))) {
            final Field hcArgsField = HaplotypeCallerEngine.class.getDeclaredField("hcArgs");
            hcArgsField.setAccessible(true);
            final HaplotypeCallerArgumentCollection hcArgs =
                    (HaplotypeCallerArgumentCollection) hcArgsField.get(ctx.engine);
            final Field assemblerField =
                    HaplotypeCallerEngine.class.getDeclaredField("assemblyEngine");
            assemblerField.setAccessible(true);
            final ReadThreadingAssembler assembler =
                    (ReadThreadingAssembler) assemblerField.get(ctx.engine);
            final Field alignerField = HaplotypeCallerEngine.class.getDeclaredField("aligner");
            alignerField.setAccessible(true);
            final SmithWatermanAligner aligner =
                    (SmithWatermanAligner) alignerField.get(ctx.engine);
            final Logger logger = LogManager.getLogger(HcFullParityGateDump.class);
            final SampleList samplesList = SampleList.singletonSampleList("s1");

            final List<SimpleInterval> intervals =
                    parseIntervals(ctx.header.getSequenceDictionary(), intervalCli);
            final List<List<Locatable>> byContig = groupByContig(intervals);
            for (final List<Locatable> contigIntervals : byContig) {
                final List<SimpleInterval> contigSimple =
                        contigIntervals.stream()
                                .map(SimpleInterval::new)
                                .collect(Collectors.toList());
                final MultiIntervalLocalReadShard shard =
                        new MultiIntervalLocalReadShard(
                                contigSimple, padding, ctx.readsSource);
                configureHcProductionReadShard(shard, ctx);
                final AssemblyRegionIterator iter =
                        new AssemblyRegionIterator(
                                shard,
                                ctx.header,
                                ctx.reference,
                                null,
                                ctx.engine,
                                ctx.asmArgs,
                                false);
                while (iter.hasNext()) {
                    final AssemblyRegion r = iter.next();
                    if (dumpInactive) {
                        if (r.isActive()) {
                            continue;
                        }
                        System.out.println("region_contig\t" + r.getContig());
                        System.out.println("region_start\t" + r.getStart());
                        System.out.println("region_end\t" + r.getEnd());
                        System.out.println("is_active\tfalse");
                        System.out.println("read_count\t0");
                        System.out.println("haplotype_count\t0");
                        return;
                    }
                    if (!r.isActive()) {
                        continue;
                    }
                    final AssemblyResultSet ars =
                            AssemblyBasedCallerUtils.assembleReads(
                                    r,
                                    Collections.emptyList(),
                                    hcArgs,
                                    ctx.header,
                                    samplesList,
                                    logger,
                                    refReader,
                                    assembler,
                                    aligner,
                                    !hcArgs.doNotCorrectOverlappingBaseQualities,
                                    hcArgs.fbargs,
                                    false);
                    System.out.println("region_contig\t" + r.getContig());
                    System.out.println("region_start\t" + r.getStart());
                    System.out.println("region_end\t" + r.getEnd());
                    System.out.println("is_active\ttrue");
                    final List<org.broadinstitute.hellbender.utils.read.GATKRead> reads =
                            new ArrayList<>(r.getReads());
                    reads.sort(
                            Comparator.comparing(
                                            (org.broadinstitute.hellbender.utils.read.GATKRead read) ->
                                                    read.getName())
                                    .thenComparingInt(
                                            org.broadinstitute.hellbender.utils.read.GATKRead
                                                    ::getStart));
                    final List<Haplotype> haps = new ArrayList<>(ars.getHaplotypeList());
                    System.out.println("read_count\t" + reads.size());
                    System.out.println("haplotype_count\t" + haps.size());
                    System.out.println("read_index\thaplotype_index\tlog10_likelihood");
                    for (int ri = 0; ri < reads.size(); ri++) {
                        final org.broadinstitute.hellbender.utils.read.GATKRead read = reads.get(ri);
                        final byte[] quals =
                                HcParityPairHmmQual.capBaseQualities(
                                        read.getBaseQualities(),
                                        read.getMappingQuality(),
                                        bqThreshold,
                                        false);
                        final String readBases =
                                new String(read.getBases(), StandardCharsets.US_ASCII);
                        for (int hi = 0; hi < haps.size(); hi++) {
                            final String hapBases =
                                    new String(
                                            haps.get(hi).getBases(), StandardCharsets.US_ASCII);
                            final double ll =
                                    HcParityNativePairHmm.pairhmmLog10Likelihood(
                                            readBases, quals, read.getMappingQuality(), hapBases);
                            System.out.printf(
                                    Locale.ROOT,
                                    "%d\t%d\t%s%n",
                                    ri,
                                    hi,
                                    HcParityNativePairHmm.formatLog10(ll));
                        }
                    }
                    return;
                }
            }
            throw new IllegalArgumentException("no active assembly region in interval");
        }
    }

    private static void pairhmmBqCap(final String[] args) throws Exception {
        if (args.length < 1) {
            usage();
        }
        final Path casesPath = Paths.get(args[0]);
        System.out.println("case_id\tcapped_quals");
        try (BufferedReader br = Files.newBufferedReader(casesPath, StandardCharsets.UTF_8)) {
            String line;
            while ((line = br.readLine()) != null) {
                final String t = line.trim();
                if (t.isEmpty() || t.startsWith("#")) {
                    continue;
                }
                final String[] c = t.split("\t");
                if (c.length < 4) {
                    throw new IllegalArgumentException("bq-cap row needs 4 cols: " + line);
                }
                final String[] qparts = c[1].split(",");
                final byte[] quals = new byte[qparts.length];
                for (int i = 0; i < qparts.length; i++) {
                    quals[i] = (byte) Integer.parseInt(qparts[i].trim());
                }
                final byte threshold = (byte) Integer.parseInt(c[2].trim());
                final int mapq = Integer.parseInt(c[3].trim());
                final byte[] capped =
                        HcParityPairHmmQual.capBaseQualities(
                                quals, mapq, threshold, true);
                final StringBuilder sb = new StringBuilder();
                for (int i = 0; i < capped.length; i++) {
                    if (i > 0) {
                        sb.append(',');
                    }
                    sb.append(capped[i] & 0xff);
                }
                System.out.println(c[0] + "\t" + sb);
            }
        }
    }

    private static void pairhmmHaplotypeFilter(final String[] args) throws Exception {
        if (args.length < 1) {
            usage();
        }
        final Path casesPath = Paths.get(args[0]);
        System.out.println("case_id\thaplotype_index\tkept\tmax_log10_likelihood");
        try (BufferedReader br = Files.newBufferedReader(casesPath, StandardCharsets.UTF_8)) {
            String line;
            while ((line = br.readLine()) != null) {
                final String t = line.trim();
                if (t.isEmpty() || t.startsWith("#")) {
                    continue;
                }
                final String[] c = t.split("\t");
                if (c.length < 5) {
                    throw new IllegalArgumentException("hap-filter row needs 5 cols: " + line);
                }
                final String caseId = c[0];
                final String readBases = c[1];
                final String[] qparts = c[2].split(",");
                final byte[] quals = new byte[qparts.length];
                for (int i = 0; i < qparts.length; i++) {
                    quals[i] = (byte) Integer.parseInt(qparts[i].trim());
                }
                final double threshold = Double.parseDouble(c[3]);
                for (int hi = 4; hi < c.length; hi++) {
                    final double ll =
                            HcParityNativePairHmm.pairhmmLog10Likelihood(
                                    readBases, quals, 60, c[hi]);
                    final String kept = ll > threshold ? "true" : "false";
                    System.out.printf(
                            Locale.ROOT,
                            "%s\t%d\t%s\t%s%n",
                            caseId,
                            hi - 4,
                            kept,
                            HcParityNativePairHmm.formatLog10(ll));
                }
            }
        }
    }

    private static void assemblyRegionFeatures(final String[] args) throws Exception {
        if (args.length < 3) {
            usage();
        }
        final String refPath = args[0];
        final String bamPath = args[1];
        final String intervalCli = args[2];
        final int padding = args.length > 3 ? parsePadding(args[3]) : DEFAULT_PADDING;
        final String featuresVcf =
                args.length > 4 && !"-".equals(args[4]) ? args[4] : null;
        System.out.println(
                "contig\tstart\tend\tis_active\textended_start\textended_end\tfeat_has_backing\tfeat_count");
        try (HcContext ctx = new HcContext(refPath, bamPath, padding)) {
            final List<SimpleInterval> intervals =
                    parseIntervals(ctx.header.getSequenceDictionary(), intervalCli);
            final List<List<Locatable>> byContig = groupByContig(intervals);
            for (final List<Locatable> contigIntervals : byContig) {
                final List<SimpleInterval> contigSimple =
                        contigIntervals.stream()
                                .map(SimpleInterval::new)
                                .collect(Collectors.toList());
                final MultiIntervalLocalReadShard shard =
                        new MultiIntervalLocalReadShard(
                                contigSimple, padding, ctx.readsSource);
                configureHcProductionReadShard(shard, ctx);
                final AssemblyRegionIterator iter =
                        new AssemblyRegionIterator(
                                shard,
                                ctx.header,
                                ctx.reference,
                                null,
                                ctx.engine,
                                ctx.asmArgs,
                                false);
                while (iter.hasNext()) {
                    final AssemblyRegion r = iter.next();
                    final SimpleInterval padded = r.getPaddedSpan();
                    final boolean hasBacking = featuresVcf != null;
                    final List<VariantContext> overlapping = new ArrayList<>();
                    if (featuresVcf != null) {
                        final File vcfFile = new File(featuresVcf);
                        try (VCFFileReader reader = new VCFFileReader(vcfFile, false)) {
                            try (CloseableIterator<VariantContext> vcIt = reader.iterator()) {
                                while (vcIt.hasNext()) {
                                    final VariantContext vc = vcIt.next();
                                    if (vc.getContig().equals(padded.getContig())
                                            && vc.getStart() <= padded.getEnd()
                                            && vc.getEnd() >= padded.getStart()) {
                                        overlapping.add(vc);
                                    }
                                }
                            }
                        }
                    }
                    overlapping.sort(
                            Comparator.comparingInt(VariantContext::getStart)
                                    .thenComparingInt(VariantContext::getEnd)
                                    .thenComparing(v -> v.getReference().getDisplayString()));
                    System.out.printf(
                            "%s\t%d\t%d\t%s\t%d\t%d\t%s\t%d%n",
                            r.getContig(),
                            r.getStart(),
                            r.getEnd(),
                            r.isActive() ? "true" : "false",
                            padded.getStart(),
                            padded.getEnd(),
                            hasBacking ? "true" : "false",
                            overlapping.size());
                    for (final VariantContext vc : overlapping) {
                        final String alts =
                                vc.getAlternateAlleles().stream()
                                        .map(a -> a.getDisplayString())
                                        .collect(Collectors.joining(","));
                        System.out.printf(
                                "feature\talleles\t%d\t%d\t%s\t%s\t%s%n",
                                vc.getStart(),
                                vc.getEnd(),
                                vc.getContig(),
                                vc.getReference().getDisplayString(),
                                alts);
                    }
                }
            }
        }
    }

    private static void assemblyRegionPileupTrack(final String[] args) throws Exception {
        if (args.length < 3) {
            usage();
        }
        final String refPath = args[0];
        final String bamPath = args[1];
        final String intervalCli = args[2];
        final int padding = args.length > 3 ? parsePadding(args[3]) : DEFAULT_PADDING;
        final boolean track =
                args.length > 4
                        && ("1".equals(args[4])
                                || "true".equalsIgnoreCase(args[4])
                                || "track".equalsIgnoreCase(args[4]));
        System.out.println("track_enabled\t" + track);
        System.out.println("contig\tstart\tend\tis_active\tpileup_count");
        try (HcContext ctx = new HcContext(refPath, bamPath, padding)) {
            final List<SimpleInterval> intervals =
                    parseIntervals(ctx.header.getSequenceDictionary(), intervalCli);
            final List<List<Locatable>> byContig = groupByContig(intervals);
            for (final List<Locatable> contigIntervals : byContig) {
                final List<SimpleInterval> contigSimple =
                        contigIntervals.stream()
                                .map(SimpleInterval::new)
                                .collect(Collectors.toList());
                final MultiIntervalLocalReadShard shard =
                        new MultiIntervalLocalReadShard(
                                contigSimple, padding, ctx.readsSource);
                configureHcProductionReadShard(shard, ctx);
                final AssemblyRegionIterator iter =
                        new AssemblyRegionIterator(
                                shard,
                                ctx.header,
                                ctx.reference,
                                null,
                                ctx.engine,
                                ctx.asmArgs,
                                track);
                while (iter.hasNext()) {
                    final AssemblyRegion r = iter.next();
                    final List<AlignmentAndReferenceContext> pileups = r.getAlignmentData();
                    System.out.printf(
                            "%s\t%d\t%d\t%s\t%d%n",
                            r.getContig(),
                            r.getStart(),
                            r.getEnd(),
                            r.isActive() ? "true" : "false",
                            pileups.size());
                    pileups.sort(
                            Comparator.comparingInt(
                                            (AlignmentAndReferenceContext p) ->
                                                    p.getAlignmentContext().getStart())
                                    .thenComparingInt(
                                            p -> p.getAlignmentContext().getEnd()));
                    for (final AlignmentAndReferenceContext p : pileups) {
                        final AlignmentContext ac = p.getAlignmentContext();
                        System.out.printf(
                                "pileup\t%s\t%d\t%d%n",
                                ac.getContig(), ac.getStart(), ac.size());
                    }
                }
            }
        }
    }

    private static void assemblyRegionTrim(final String[] args) throws Exception {
        if (args.length < 6) {
            usage();
        }
        final String refPath = args[0];
        final String contig = args[1];
        final int start = Integer.parseInt(args[2]);
        final int end = Integer.parseInt(args[3]);
        final int extStart = Integer.parseInt(args[4]);
        final int extEnd = Integer.parseInt(args[5]);
        final String variantsPath =
                args.length > 6 && !"-".equals(args[6]) ? args[6] : null;
        final boolean legacy =
                args.length > 7
                        && ("1".equals(args[7])
                                || "true".equalsIgnoreCase(args[7])
                                || "legacy".equalsIgnoreCase(args[7]));
        System.out.println(
                "contig\torig_start\torig_end\torig_ext_start\torig_ext_end\tvariation_present\tvariant_start\tvariant_end\ttrim_start\ttrim_end\ttrim_ext_start\ttrim_ext_end");
        final CachingIndexedFastaSequenceFile refReader =
                new CachingIndexedFastaSequenceFile(Paths.get(refPath));
        try {
            final SAMSequenceDictionary dict = refReader.getSequenceDictionary();
            final SAMFileHeader header = new SAMFileHeader(dict);
            final AssemblyRegionArgumentCollection asmArgs = new AssemblyRegionArgumentCollection();
            asmArgs.enableLegacyAssemblyRegionTrimming = legacy;
            final AssemblyRegionTrimmer trimmer = new AssemblyRegionTrimmer(asmArgs, dict);
            final SimpleInterval active = new SimpleInterval(contig, start, end);
            final SimpleInterval padded = new SimpleInterval(contig, extStart, extEnd);
            final AssemblyRegion region = new AssemblyRegion(active, padded, true, header);
            final ReferenceContext refCtx =
                    new ReferenceContext(ReferenceDataSource.of(Paths.get(refPath)), padded);
            final TreeSet<VariantContext> variants = loadTrimVariantsTsv(variantsPath);
            final AssemblyRegionTrimmer.Result result =
                    trimmer.trim(region, variants, refCtx);
            final AssemblyRegion trimmed =
                    result.isVariationPresent() ? result.getVariantRegion() : region;
            final String vs;
            final String ve;
            if (result.isVariationPresent()) {
                final AssemblyRegion variantOnly = result.getVariantRegion();
                vs = Integer.toString(variantOnly.getStart());
                ve = Integer.toString(variantOnly.getEnd());
            } else {
                vs = "-";
                ve = "-";
            }
            System.out.printf(
                    "%s\t%d\t%d\t%d\t%d\t%s\t%s\t%s\t%d\t%d\t%d\t%d%n",
                    region.getContig(),
                    region.getStart(),
                    region.getEnd(),
                    region.getPaddedSpan().getStart(),
                    region.getPaddedSpan().getEnd(),
                    result.isVariationPresent() ? "true" : "false",
                    vs,
                    ve,
                    trimmed.getStart(),
                    trimmed.getEnd(),
                    trimmed.getPaddedSpan().getStart(),
                    trimmed.getPaddedSpan().getEnd());
        } finally {
            refReader.close();
        }
    }

    private static TreeSet<VariantContext> loadTrimVariantsTsv(final String path)
            throws Exception {
        final TreeSet<VariantContext> out = new TreeSet<>(Comparator.comparingInt(VariantContext::getStart));
        if (path == null) {
            return out;
        }
        try (BufferedReader br = new BufferedReader(new FileReader(path))) {
            String line;
            while ((line = br.readLine()) != null) {
                line = line.trim();
                if (line.isEmpty() || line.startsWith("#")) {
                    continue;
                }
                final String[] cols = line.split("\t");
                if (cols.length < 4) {
                    throw new IllegalArgumentException("expected contig\\tstart\\tend\\tis_indel: " + line);
                }
                final String contig = cols[0];
                final int pos = Integer.parseInt(cols[1]);
                final int stop = Integer.parseInt(cols[2]);
                final boolean isIndel =
                        "true".equalsIgnoreCase(cols[3])
                                || "1".equals(cols[3])
                                || "yes".equalsIgnoreCase(cols[3]);
                final List<Allele> alleles;
                if (isIndel) {
                    alleles =
                            Arrays.asList(
                                    Allele.create("ACG", true), Allele.create("A", false));
                } else {
                    alleles =
                            Arrays.asList(
                                    Allele.create("A", true), Allele.create("G", false));
                }
                final VariantContext vc =
                        new VariantContextBuilder("parity", contig, pos, stop, alleles).make();
                out.add(vc);
            }
        }
        return out;
    }

    private static void assemblyRegions(final String[] args, final boolean forceActive)
            throws Exception {
        if (args.length < 3) {
            usage();
        }
        final String refPath = args[0];
        final String bamPath = args[1];
        final String intervalCli = args[2];
        final int padding = args.length > 3 ? parsePadding(args[3]) : DEFAULT_PADDING;
        System.out.println(
                "contig\tstart\tend\tis_active\textended_start\textended_end\textension");
        try (HcContext ctx = new HcContext(refPath, bamPath, padding)) {
            final List<SimpleInterval> intervals =
                    parseIntervals(ctx.header.getSequenceDictionary(), intervalCli);
            final List<List<Locatable>> byContig =
                    groupByContig(intervals);
            for (final List<Locatable> contigIntervals : byContig) {
                final List<SimpleInterval> contigSimple =
                        contigIntervals.stream()
                                .map(SimpleInterval::new)
                                .collect(Collectors.toList());
                final MultiIntervalLocalReadShard shard =
                        new MultiIntervalLocalReadShard(
                                contigSimple, padding, ctx.readsSource);
                configureHcProductionReadShard(shard, ctx);
                final AssemblyRegionIterator iter =
                        new AssemblyRegionIterator(
                                shard,
                                ctx.header,
                                ctx.reference,
                                null,
                                ctx.engine,
                                ctx.asmArgs,
                                false);
                while (iter.hasNext()) {
                    final AssemblyRegion r = iter.next();
                    final boolean active = forceActive || r.isActive();
                    System.out.printf(
                            "%s\t%d\t%d\t%s\t%d\t%d\t%d%n",
                            r.getContig(),
                            r.getStart(),
                            r.getEnd(),
                            active ? "true" : "false",
                            r.getPaddedSpan().getStart(),
                            r.getPaddedSpan().getEnd(),
                            ctx.asmArgs.assemblyRegionPadding);
                }
            }
        }
    }

    private static void applySummary(final String[] args) throws Exception {
        if (args.length < 3) {
            usage();
        }
        final String refPath = args[0];
        final String bamPath = args[1];
        final String intervalCli = args[2];
        final int padding = args.length > 3 ? parsePadding(args[3]) : DEFAULT_PADDING;
        int total = 0;
        int inactive = 0;
        int active = 0;
        try (HcContext ctx = new HcContext(refPath, bamPath, padding)) {
            final List<SimpleInterval> intervals =
                    parseIntervals(ctx.header.getSequenceDictionary(), intervalCli);
            final List<List<Locatable>> byContig =
                    groupByContig(intervals);
            for (final List<Locatable> contigIntervals : byContig) {
                final List<SimpleInterval> contigSimple =
                        contigIntervals.stream()
                                .map(SimpleInterval::new)
                                .collect(Collectors.toList());
                final MultiIntervalLocalReadShard shard =
                        new MultiIntervalLocalReadShard(
                                contigSimple, padding, ctx.readsSource);
                configureHcProductionReadShard(shard, ctx);
                final AssemblyRegionIterator iter =
                        new AssemblyRegionIterator(
                                shard,
                                ctx.header,
                                ctx.reference,
                                null,
                                ctx.engine,
                                ctx.asmArgs,
                                false);
                while (iter.hasNext()) {
                    total++;
                    if (iter.next().isActive()) {
                        active++;
                    } else {
                        inactive++;
                    }
                }
            }
        }
        System.out.println("total_apply\tinactive_fast_path\tactive_full");
        System.out.printf("%d\t%d\t%d%n", total, inactive, active);
    }

    /** Padded reference + reads for first active region (ASM-1 Java/Rust side-by-side). */
    private static final class RegionAssemblyMaterial {
        final String contig;
        final int start;
        final int end;
        final byte[] refBases;
        final byte[] refQuals;
        final List<byte[]> readBases;
        final List<byte[]> readQuals;

        RegionAssemblyMaterial(
                final String contig,
                final int start,
                final int end,
                final byte[] refBases,
                final byte[] refQuals,
                final List<byte[]> readBases,
                final List<byte[]> readQuals) {
            this.contig = contig;
            this.start = start;
            this.end = end;
            this.refBases = refBases;
            this.refQuals = refQuals;
            this.readBases = readBases;
            this.readQuals = readQuals;
        }
    }

    private static RegionAssemblyMaterial materializeFirstActiveRegion(
            final String refPath,
            final String bamPath,
            final String intervalCli,
            final int padding)
            throws Exception {
        try (HcContext ctx = new HcContext(refPath, bamPath, padding);
                CachingIndexedFastaSequenceFile refReader =
                        new CachingIndexedFastaSequenceFile(Paths.get(refPath))) {
            final List<SimpleInterval> intervals =
                    parseIntervals(ctx.header.getSequenceDictionary(), intervalCli);
            for (final List<Locatable> contigIntervals : groupByContig(intervals)) {
                final List<SimpleInterval> contigSimple =
                        contigIntervals.stream()
                                .map(SimpleInterval::new)
                                .collect(Collectors.toList());
                final MultiIntervalLocalReadShard shard =
                        new MultiIntervalLocalReadShard(
                                contigSimple, padding, ctx.readsSource);
                configureHcProductionReadShard(shard, ctx);
                final AssemblyRegionIterator iter =
                        new AssemblyRegionIterator(
                                shard,
                                ctx.header,
                                ctx.reference,
                                null,
                                ctx.engine,
                                ctx.asmArgs,
                                false);
                RegionAssemblyMaterial earliest = null;
                while (iter.hasNext()) {
                    final AssemblyRegion r = iter.next();
                    if (!r.isActive()) {
                        continue;
                    }
                    if (earliest != null && r.getStart() >= earliest.start) {
                        continue;
                    }
                    final byte[] refBases =
                            r.getAssemblyRegionReference(
                                    refReader,
                                    AssemblyBasedCallerUtils.REFERENCE_PADDING_FOR_ASSEMBLY);
                    final byte[] refQuals = new byte[refBases.length];
                    Arrays.fill(refQuals, (byte) 30);
                    final List<byte[]> readBases = new ArrayList<>();
                    final List<byte[]> readQuals = new ArrayList<>();
                    for (final GATKRead read : r.getReads()) {
                        readBases.add(Arrays.copyOf(read.getBases(), read.getLength()));
                        readQuals.add(Arrays.copyOf(read.getBaseQualities(), read.getLength()));
                    }
                    earliest =
                            new RegionAssemblyMaterial(
                                    r.getContig(),
                                    r.getStart(),
                                    r.getEnd(),
                                    refBases,
                                    refQuals,
                                    readBases,
                                    readQuals);
                }
                if (earliest != null) {
                    return earliest;
                }
            }
        }
        throw new IllegalArgumentException("no active assembly region in interval");
    }

    private static ReadThreadingGraph buildReadThreadingGraphFromRegion(
            final RegionAssemblyMaterial material, final int kmerSize, final int minQual) {
        final ReadThreadingGraph graph = new ReadThreadingGraph(kmerSize);
        addReadThreadingSequence(
                graph, "ref", material.refBases, material.refQuals, minQual, true);
        int i = 0;
        for (int r = 0; r < material.readBases.size(); r++) {
            addReadThreadingSequence(
                    graph,
                    "r" + (i++),
                    material.readBases.get(r),
                    material.readQuals.get(r),
                    minQual,
                    false);
        }
        graph.buildGraphIfNecessary();
        return graph;
    }

    private static boolean referenceHasNonUniqueKmers(
            final byte[] refBases, final int kmerSize) {
        if (refBases.length < kmerSize) {
            return false;
        }
        final Set<String> seen = new HashSet<>();
        final int stop = refBases.length - kmerSize;
        for (int i = 0; i <= stop; i++) {
            final String kmer =
                    new String(refBases, i, kmerSize, StandardCharsets.US_ASCII);
            if (!seen.add(kmer)) {
                return true;
            }
        }
        return false;
    }

    @SuppressWarnings("unchecked")
    private static List<Integer> assemblerKmerSizes(final ReadThreadingAssembler assembler)
            throws Exception {
        final Field f = ReadThreadingAssembler.class.getDeclaredField("kmerSizes");
        f.setAccessible(true);
        return new ArrayList<>((List<Integer>) f.get(assembler));
    }

    private static boolean assemblerDontIncreaseKmerSizes(final ReadThreadingAssembler assembler)
            throws Exception {
        final Field f =
                ReadThreadingAssembler.class.getDeclaredField("dontIncreaseKmerSizesForCycles");
        f.setAccessible(true);
        return f.getBoolean(assembler);
    }

    private static boolean assemblerAllowNonUniqueInRef(final ReadThreadingAssembler assembler)
            throws Exception {
        final Field f =
                ReadThreadingAssembler.class.getDeclaredField("allowNonUniqueKmersInRef");
        f.setAccessible(true);
        return f.getBoolean(assembler);
    }

    private static ReadThreadingAssembler productionAssemblerFromContext(final HcContext ctx)
            throws Exception {
        final Field assemblerField =
                HaplotypeCallerEngine.class.getDeclaredField("assemblyEngine");
        assemblerField.setAccessible(true);
        return (ReadThreadingAssembler) assemblerField.get(ctx.engine);
    }

    private static int assemblerPruneFactor(final ReadThreadingAssembler assembler)
            throws Exception {
        final Field f = ReadThreadingAssembler.class.getDeclaredField("pruneFactor");
        f.setAccessible(true);
        return f.getInt(assembler);
    }

    private static int assemblerMinDanglingBranchLength(final ReadThreadingAssembler assembler)
            throws Exception {
        final Field f =
                ReadThreadingAssembler.class.getDeclaredField("minDanglingBranchLength");
        f.setAccessible(true);
        return f.getInt(assembler);
    }

    private static Set<MultiDeBruijnVertex> rtRefSpineVertices(final ReadThreadingGraph graph) {
        final MultiDeBruijnVertex source = graph.getReferenceSourceVertex();
        final MultiDeBruijnVertex sink = graph.getReferenceSinkVertex();
        if (source == null || sink == null) {
            return Collections.emptySet();
        }
        final Set<MultiDeBruijnVertex> fromSource = new HashSet<>();
        final Deque<MultiDeBruijnVertex> stack = new ArrayDeque<>();
        stack.push(source);
        fromSource.add(source);
        while (!stack.isEmpty()) {
            final MultiDeBruijnVertex v = stack.pop();
            for (final MultiDeBruijnVertex t : graph.outgoingVerticesOf(v)) {
                if (fromSource.add(t)) {
                    stack.push(t);
                }
            }
        }
        final Set<MultiDeBruijnVertex> fromSink = new HashSet<>();
        stack.push(sink);
        fromSink.add(sink);
        while (!stack.isEmpty()) {
            final MultiDeBruijnVertex v = stack.pop();
            for (final MultiDeBruijnVertex p : graph.incomingVerticesOf(v)) {
                if (fromSink.add(p)) {
                    stack.push(p);
                }
            }
        }
        fromSource.retainAll(fromSink);
        return fromSource;
    }

    private static int countRtBranches(final ReadThreadingGraph graph, final Set<MultiDeBruijnVertex> spine) {
        int n = 0;
        for (final MultiDeBruijnVertex v : spine) {
            if (graph.outDegreeOf(v) > 1) {
                n++;
            }
        }
        return n;
    }

    private static int countRtNonRefEdgesOnSpine(
            final ReadThreadingGraph graph, final Set<MultiDeBruijnVertex> spine) {
        int n = 0;
        for (final MultiSampleEdge e : graph.edgeSet()) {
            final MultiDeBruijnVertex from = graph.getEdgeSource(e);
            final MultiDeBruijnVertex to = graph.getEdgeTarget(e);
            if (spine.contains(from) && spine.contains(to) && !e.isRef()) {
                n++;
            }
        }
        return n;
    }

    private static int countRtBranchesAll(final ReadThreadingGraph graph) {
        int n = 0;
        for (final MultiDeBruijnVertex v : graph.vertexSet()) {
            if (graph.outDegreeOf(v) > 1) {
                n++;
            }
        }
        return n;
    }

    private static int countRtNonRefEdgesAll(final ReadThreadingGraph graph) {
        int n = 0;
        for (final MultiSampleEdge e : graph.edgeSet()) {
            if (!e.isRef()) {
                n++;
            }
        }
        return n;
    }

    private static final class AssemblyStageRow {
        final String graphKind;
        final String stage;
        final int nodes;
        final int edges;
        final int spineVerts;
        final int branchSpine;
        final int branchAll;
        final int nonRefSpine;
        final int nonRefAll;
        final int kbestPaths;
        final int extractedHaps;
        final int nonRefHaps;
        final int topPathLen;
        final boolean topPathEqRef;

        AssemblyStageRow(
                final String graphKind,
                final String stage,
                final int nodes,
                final int edges,
                final int spineVerts,
                final int branchSpine,
                final int branchAll,
                final int nonRefSpine,
                final int nonRefAll,
                final int kbestPaths,
                final int extractedHaps,
                final int nonRefHaps,
                final int topPathLen,
                final boolean topPathEqRef) {
            this.graphKind = graphKind;
            this.stage = stage;
            this.nodes = nodes;
            this.edges = edges;
            this.spineVerts = spineVerts;
            this.branchSpine = branchSpine;
            this.branchAll = branchAll;
            this.nonRefSpine = nonRefSpine;
            this.nonRefAll = nonRefAll;
            this.kbestPaths = kbestPaths;
            this.extractedHaps = extractedHaps;
            this.nonRefHaps = nonRefHaps;
            this.topPathLen = topPathLen;
            this.topPathEqRef = topPathEqRef;
        }
    }

    private static void printAssemblyStageRow(final AssemblyStageRow r) {
        System.out.printf(
                "%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s%n",
                r.graphKind,
                r.stage,
                r.nodes,
                r.edges,
                r.spineVerts,
                r.branchSpine,
                r.branchAll,
                r.nonRefSpine,
                r.nonRefAll,
                r.kbestPaths,
                r.extractedHaps,
                r.nonRefHaps,
                r.topPathLen,
                r.topPathEqRef ? "true" : "false");
    }

    private static AssemblyStageRow rtStageRow(
            final String stage,
            final ReadThreadingGraph graph,
            final byte[] refBases,
            final int maxHaps)
            throws Exception {
        final Set<MultiDeBruijnVertex> spine = rtRefSpineVertices(graph);
        final List<KBestHaplotype<MultiDeBruijnVertex, MultiSampleEdge>> paths =
                rankKbestPaths(graph, maxHaps);
        int nonRefHaps = 0;
        for (final KBestHaplotype<MultiDeBruijnVertex, MultiSampleEdge> p : paths) {
            if (!p.isReference()) {
                nonRefHaps++;
            }
        }
        final int topLen =
                paths.isEmpty() ? 0 : paths.get(0).getBases().length;
        final boolean topEq =
                !paths.isEmpty()
                        && Arrays.equals(paths.get(0).getBases(), refBases);
        return new AssemblyStageRow(
                "rt",
                stage,
                graph.vertexSet().size(),
                graph.edgeSet().size(),
                spine.size(),
                countRtBranches(graph, spine),
                countRtBranchesAll(graph),
                countRtNonRefEdgesOnSpine(graph, spine),
                countRtNonRefEdgesAll(graph),
                paths.size(),
                paths.size(),
                nonRefHaps,
                topLen,
                topEq);
    }

    private static Set<SeqVertex> seqRefSpineVertices(final SeqGraph graph) {
        final SeqVertex source = graph.getReferenceSourceVertex();
        final SeqVertex sink = graph.getReferenceSinkVertex();
        if (source == null || sink == null) {
            return Collections.emptySet();
        }
        final Set<SeqVertex> fromSource = new HashSet<>();
        final Deque<SeqVertex> stack = new ArrayDeque<>();
        stack.push(source);
        fromSource.add(source);
        while (!stack.isEmpty()) {
            final SeqVertex v = stack.pop();
            for (final SeqVertex t : graph.outgoingVerticesOf(v)) {
                if (fromSource.add(t)) {
                    stack.push(t);
                }
            }
        }
        final Set<SeqVertex> fromSink = new HashSet<>();
        stack.push(sink);
        fromSink.add(sink);
        while (!stack.isEmpty()) {
            final SeqVertex v = stack.pop();
            for (final SeqVertex p : graph.incomingVerticesOf(v)) {
                if (fromSink.add(p)) {
                    stack.push(p);
                }
            }
        }
        fromSource.retainAll(fromSink);
        return fromSource;
    }

    private static AssemblyStageRow seqStageRow(
            final String stage,
            final SeqGraph graph,
            final byte[] refBases,
            final int maxHaps)
            throws Exception {
        final Set<SeqVertex> spine = seqRefSpineVertices(graph);
        int branchSpine = 0;
        int branchAll = 0;
        int nonRefSpine = 0;
        int nonRefAll = 0;
        for (final SeqVertex v : graph.vertexSet()) {
            if (graph.outDegreeOf(v) > 1) {
                branchAll++;
                if (spine.contains(v)) {
                    branchSpine++;
                }
            }
        }
        for (final BaseEdge e : graph.edgeSet()) {
            if (!e.isRef()) {
                nonRefAll++;
                final SeqVertex from = graph.getEdgeSource(e);
                final SeqVertex to = graph.getEdgeTarget(e);
                if (spine.contains(from) && spine.contains(to)) {
                    nonRefSpine++;
                }
            }
        }
        final SeqVertex source = graph.getReferenceSourceVertex();
        final SeqVertex sink = graph.getReferenceSinkVertex();
        final List<KBestHaplotype<SeqVertex, BaseEdge>> paths;
        if (source != null && sink != null) {
            paths =
                    new ArrayList<>(
                            new GraphBasedKBestHaplotypeFinder<>(graph, source, sink)
                                    .findBestHaplotypes(maxHaps));
        } else {
            paths = Collections.emptyList();
        }
        int nonRefHaps = 0;
        for (final KBestHaplotype<SeqVertex, BaseEdge> p : paths) {
            if (!p.isReference()) {
                nonRefHaps++;
            }
        }
        final int topLen = paths.isEmpty() ? 0 : paths.get(0).getBases().length;
        final boolean topEq =
                !paths.isEmpty()
                        && Arrays.equals(paths.get(0).getBases(), refBases);
        return new AssemblyStageRow(
                "seq",
                stage,
                graph.vertexSet().size(),
                graph.edgeSet().size(),
                spine.size(),
                branchSpine,
                branchAll,
                nonRefSpine,
                nonRefAll,
                paths.size(),
                paths.size(),
                nonRefHaps,
                topLen,
                topEq);
    }

    private static void assemblyRegionKmerProbe(final String[] args) throws Exception {
        if (args.length < 3) {
            usage();
        }
        final String refPath = args[0];
        final String bamPath = args[1];
        final String intervalCli = args[2];
        final int padding = args.length > 3 ? parsePadding(args[3]) : DEFAULT_PADDING;
        final RegionAssemblyMaterial material =
                materializeFirstActiveRegion(refPath, bamPath, intervalCli, padding);
        try (HcContext ctx = new HcContext(refPath, bamPath, padding)) {
            final ReadThreadingAssembler assembler = productionAssemblerFromContext(ctx);
            final int minQual = assembler.getMinBaseQualityToUseInAssembly();
            final int minPrune = assemblerPruneFactor(assembler);
            final int minDangling = assemblerMinDanglingBranchLength(assembler);
            final boolean recoverHeads = assembler.isRecoverDanglingBranches();
            final boolean dontIncrease = assemblerDontIncreaseKmerSizes(assembler);
            final boolean allowNonUniqueRef = assemblerAllowNonUniqueInRef(assembler);
            final List<Integer> kmerSizes = assemblerKmerSizes(assembler);
            Collections.sort(kmerSizes);

            System.out.println("region_contig\t" + material.contig);
            System.out.println("region_start\t" + material.start);
            System.out.println("region_end\t" + material.end);
            System.out.println("padded_ref_len\t" + material.refBases.length);
            System.out.println("read_count\t" + material.readBases.size());
            System.out.println(
                    "phase\tkmer\tallow_low_complexity\tallow_non_unique\toutcome\tthread_nodes\tthread_edges\tcleanup_status\thas_ref_source\thas_ref_sink\tref_path_matches\tkbest_paths\textracted_haps\tnon_ref_haps\tpath_bases_len\tpath_eq_ref");

            for (final int kmerSize : kmerSizes) {
                emitKmerProbeRow(
                        material,
                        kmerSize,
                        "configured",
                        dontIncrease,
                        allowNonUniqueRef,
                        minQual,
                        minPrune,
                        minDangling,
                        recoverHeads);
            }
            if (!dontIncrease) {
                int kmerSize = kmerSizes.isEmpty() ? 35 : kmerSizes.get(kmerSizes.size() - 1) + 10;
                for (int iter = 1; iter <= 6; iter++) {
                    final boolean last = iter == 6;
                    emitKmerProbeRow(
                            material,
                            kmerSize,
                            "expanded",
                            last || dontIncrease,
                            last || allowNonUniqueRef,
                            minQual,
                            minPrune,
                            minDangling,
                            recoverHeads);
                    kmerSize += 10;
                }
            }
        }
    }

    private static void emitKmerProbeRow(
            final RegionAssemblyMaterial material,
            final int kmerSize,
            final String phase,
            final boolean allowLowComplexity,
            final boolean allowNonUnique,
            final int minQual,
            final int minPrune,
            final int minDangling,
            final boolean recoverHeads)
            throws Exception {
        String outcome = "";
        int threadNodes = 0;
        int threadEdges = 0;
        String cleanupStatus = "";
        boolean hasRefSource = false;
        boolean hasRefSink = false;
        boolean refPathMatches = false;
        int kbestPaths = 0;
        int extractedHaps = 0;
        int nonRefHaps = 0;
        int pathBasesLen = 0;
        boolean pathEqRef = false;

        if (material.refBases.length < kmerSize) {
            outcome = "ref_shorter_than_kmer";
        } else if (!allowNonUnique
                && referenceHasNonUniqueKmers(material.refBases, kmerSize)) {
            outcome = "skip_non_unique_ref_kmers";
        } else {
            final ReadThreadingGraph rt =
                    prepareReadThreadingGraphForHaplotypeDumpFromRegion(
                            material,
                            kmerSize,
                            minQual,
                            minPrune,
                            minDangling,
                            recoverHeads);
            if (rt == null) {
                outcome = "no_threading_graph";
            } else {
                threadNodes = rt.vertexSet().size();
                threadEdges = rt.edgeSet().size();
                final SeqGraph seqGraph = rt.toSequenceGraph();
                    seqGraph.cleanNonRefPaths();
                    cleanupStatus = cleanupSeqGraphRustParity(seqGraph);
                    hasRefSource = seqGraph.getReferenceSourceVertex() != null;
                    hasRefSink = seqGraph.getReferenceSinkVertex() != null;
                    if (hasRefSource && hasRefSink) {
                        final byte[] refPathBytes =
                                seqGraph.getReferenceBytes(
                                        seqGraph.getReferenceSourceVertex(),
                                        seqGraph.getReferenceSinkVertex(),
                                        true,
                                        true);
                        refPathMatches = Arrays.equals(refPathBytes, material.refBases);
                    }
                    if ("just_assembled_reference".equals(cleanupStatus)) {
                        outcome = "cleanup_just_assembled_reference";
                    } else if (!hasRefSource || !hasRefSink) {
                        outcome = "dropped_no_ref_endpoints";
                    } else {
                        final List<KBestHaplotype<SeqVertex, BaseEdge>> paths =
                                new ArrayList<>(
                                        new GraphBasedKBestHaplotypeFinder<>(
                                                        seqGraph,
                                                        seqGraph.getReferenceSourceVertex(),
                                                        seqGraph.getReferenceSinkVertex())
                                                .findBestHaplotypes(128));
                        kbestPaths = paths.size();
                        extractedHaps = paths.size();
                        for (final KBestHaplotype<SeqVertex, BaseEdge> p : paths) {
                            if (!p.isReference()) {
                                nonRefHaps++;
                            }
                        }
                        if (!paths.isEmpty()) {
                            pathBasesLen = paths.get(0).getBases().length;
                            pathEqRef =
                                    Arrays.equals(
                                            paths.get(0).getBases(), material.refBases);
                        }
                        outcome =
                                pathEqRef
                                        ? "variation_graph_kbest_path_is_ref_bases"
                                        : (nonRefHaps > 0
                                                ? "variation_graph_with_alt_haps"
                                                : (kbestPaths > 1
                                                        ? "variation_graph_kbest_no_extracted_alts"
                                                        : "variation_graph_ref_only_kbest"));
                }
            }
        }
        System.out.printf(
                "%s\t%d\t%s\t%s\t%s\t%d\t%d\t%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%s%n",
                phase,
                kmerSize,
                allowLowComplexity ? "true" : "false",
                allowNonUnique ? "true" : "false",
                outcome,
                threadNodes,
                threadEdges,
                cleanupStatus,
                hasRefSource ? "true" : "false",
                hasRefSink ? "true" : "false",
                refPathMatches ? "true" : "false",
                kbestPaths,
                extractedHaps,
                nonRefHaps,
                pathBasesLen,
                pathEqRef ? "true" : "false");
    }

    /** Rust {@code build_threading_graph_for_haplotype_dump} order on region material. */
    private static ReadThreadingGraph prepareReadThreadingGraphForHaplotypeDumpFromRegion(
            final RegionAssemblyMaterial material,
            final int kmerSize,
            final int minQual,
            final int minPrune,
            final int minDangling,
            final boolean recoverHeads)
            throws Exception {
        if (material.refBases.length < kmerSize) {
            return null;
        }
        final ReadThreadingGraph graph =
                buildReadThreadingGraphFromRegion(material, kmerSize, minQual);
        final ChainPruner<MultiDeBruijnVertex, MultiSampleEdge> pruner =
                makeChainPruner(minPrune, false);
        pruner.pruneLowWeightChains(graph);
        final SmithWatermanAligner aligner =
                SmithWatermanAligner.getAligner(SmithWatermanAligner.Implementation.JAVA);
        final SWParameters swParams = danglingEndSwParameters();
        if (graph.getReferenceSourceVertex() != null) {
            graph.recoverDanglingTails(minPrune, minDangling, false, aligner, swParams);
            if (recoverHeads) {
                graph.recoverDanglingHeads(minPrune, minDangling, false, aligner, swParams);
            }
        }
        if (graph.getReferenceSourceVertex() == null || graph.getReferenceSinkVertex() == null) {
            return null;
        }
        graph.removePathsNotConnectedToRef();
        if (graph.getReferenceSourceVertex() == null || graph.getReferenceSinkVertex() == null) {
            return null;
        }
        return graph;
    }

    private static void printRtRegionKbestRows(
            final String stage,
            final boolean stripCycles,
            final byte[] refBases,
            final List<KBestHaplotype<MultiDeBruijnVertex, MultiSampleEdge>> paths) {
        int rank = 0;
        for (final KBestHaplotype<MultiDeBruijnVertex, MultiSampleEdge> path : paths) {
            final byte[] bases = path.getBases();
            System.out.printf(
                    "rt\t%s\t%s\t%d\t%s\t%s\t%d\t%s\t%s%n",
                    stage,
                    stripCycles,
                    rank++,
                    formatScore(path.score()),
                    path.isReference(),
                    bases.length,
                    Arrays.equals(bases, refBases),
                    new String(bases));
        }
    }

    private static void printSeqRegionKbestRows(
            final String stage,
            final byte[] refBases,
            final SeqGraph seq,
            final int maxHaps)
            throws Exception {
        final SeqVertex source = seq.getReferenceSourceVertex();
        final SeqVertex sink = seq.getReferenceSinkVertex();
        final List<KBestHaplotype<SeqVertex, BaseEdge>> paths;
        if (source != null && sink != null) {
            paths =
                    new ArrayList<>(
                            new GraphBasedKBestHaplotypeFinder<>(seq, source, sink)
                                    .findBestHaplotypes(maxHaps));
        } else {
            paths = Collections.emptyList();
        }
        int rank = 0;
        for (final KBestHaplotype<SeqVertex, BaseEdge> path : paths) {
            final byte[] bases = path.getBases();
            System.out.printf(
                    "seq\t%s\tfalse\t%d\t%s\t%s\t%d\t%s\t%s%n",
                    stage,
                    rank++,
                    formatScore(path.score()),
                    path.isReference(),
                    bases.length,
                    Arrays.equals(bases, refBases),
                    new String(bases));
        }
    }

    private static void assemblyRegionKbestPaths(final String[] args) throws Exception {
        if (args.length < 3) {
            usage();
        }
        final String refPath = args[0];
        final String bamPath = args[1];
        final String intervalCli = args[2];
        final int padding = args.length > 3 ? parsePadding(args[3]) : DEFAULT_PADDING;
        final int maxHaps = 128;
        final RegionAssemblyMaterial material =
                materializeFirstActiveRegion(refPath, bamPath, intervalCli, padding);
        try (HcContext ctx = new HcContext(refPath, bamPath, padding)) {
            final ReadThreadingAssembler assembler = productionAssemblerFromContext(ctx);
            final int minQual = assembler.getMinBaseQualityToUseInAssembly();
            final int minPrune = assemblerPruneFactor(assembler);
            final int minDangling = assemblerMinDanglingBranchLength(assembler);
            final boolean recoverHeads = assembler.isRecoverDanglingBranches();
            final int kmer = 85;
            System.out.println("region_contig\t" + material.contig);
            System.out.println("region_start\t" + material.start);
            System.out.println("region_end\t" + material.end);
            System.out.println("padded_ref_len\t" + material.refBases.length);
            if (material.refBases.length < kmer) {
                System.out.println(
                        "graph\tstage\tstrip_cycles\trank\tscore\tis_reference\tpath_len\teq_ref\tsequence");
                return;
            }
            ReadThreadingGraph graph =
                    buildReadThreadingGraphFromRegion(material, kmer, minQual);
            if (graph.hasCycles()) {
                System.out.println("warn\tgraph_has_cycles\ttrue");
            }
            final ChainPruner<MultiDeBruijnVertex, MultiSampleEdge> pruner =
                    makeChainPruner(minPrune, false);
            pruner.pruneLowWeightChains(graph);
            final SmithWatermanAligner aligner =
                    SmithWatermanAligner.getAligner(SmithWatermanAligner.Implementation.JAVA);
            final SWParameters swParams = danglingEndSwParameters();
            if (graph.getReferenceSourceVertex() != null) {
                graph.recoverDanglingTails(minPrune, minDangling, false, aligner, swParams);
                if (recoverHeads) {
                    graph.recoverDanglingHeads(
                            minPrune, minDangling, false, aligner, swParams);
                }
            }
            System.out.println(
                    "graph\tstage\tstrip_cycles\trank\tscore\tis_reference\tpath_len\teq_ref\tsequence");
            printRtRegionKbestRows(
                    "threading_after_dangling_pre_remove_paths",
                    true,
                    material.refBases,
                    rankKbestPaths(graph, maxHaps));
            if (graph.getReferenceSourceVertex() != null
                    && graph.getReferenceSinkVertex() != null) {
                graph.removePathsNotConnectedToRef();
            }
            if (graph.getReferenceSourceVertex() != null
                    && graph.getReferenceSinkVertex() != null) {
                printRtRegionKbestRows(
                        "threading_after_remove_paths",
                        true,
                        material.refBases,
                        rankKbestPaths(graph, maxHaps));
                SeqGraph seq = graph.toSequenceGraph();
                seq.cleanNonRefPaths();
                printSeqRegionKbestRows(
                        "seq_after_to_sequence_graph",
                        material.refBases,
                        seq,
                        maxHaps);
                final SeqGraph seqCleanup = seq.clone();
                seqCleanup.zipLinearChains();
                seqCleanup.removeSingletonOrphanVertices();
                seqCleanup.removeVerticesNotConnectedToRefRegardlessOfEdgeDirection();
                seqCleanup.simplifyGraph();
                if (seqCleanup.getReferenceSourceVertex() != null
                        && seqCleanup.getReferenceSinkVertex() != null) {
                    seqCleanup.removePathsNotConnectedToRef();
                    seqCleanup.simplifyGraph();
                    if (seqCleanup.vertexSet().size() == 1) {
                        final SeqVertex complete =
                                seqCleanup.vertexSet().iterator().next();
                        final SeqVertex dummy = new SeqVertex("");
                        seqCleanup.addVertex(dummy);
                        seqCleanup.addEdge(complete, dummy, new BaseEdge(true, 0));
                    }
                    printSeqRegionKbestRows(
                            "seq_after_cleanup_seq_graph",
                            material.refBases,
                            seqCleanup,
                            maxHaps);
                }
            }
        }
    }

    /**
     * 6R.79 forensic: production k=25 SeqGraph vertices/edges/multiplicities after
     * {@code cleanupSeqGraph}, plus k-best scores. Graph reference is padding-0
     * extended span ({@code getAssemblyRegionReference(reader, 0)}).
     *
     * <p>Args: ref bam interval [coverLocus] [padding]
     */
    private static void assemblyRegionSeqgraphEdges(final String[] args) throws Exception {
        if (args.length < 3) {
            usage();
        }
        final String refPath = args[0];
        final String bamPath = args[1];
        final String intervalCli = args[2];
        final int coverLocus = args.length > 3 ? Integer.parseInt(args[3]) : -1;
        final int padding = args.length > 4 ? parsePadding(args[4]) : DEFAULT_PADDING;
        final String j0 =
                "CATGGAGCCTGACCTTATTTGAAGTAGGGCATTTGCAGATGTATTTAAGATATTTGAGGCTGGGCACAGTGGCTCACGTCTGTAATCCCAGCACTTTGAAAGGCCGAGGCAGGTGGATTCACCTGAGGTCAGGAGTTTGAGACCAGCCTGTCCCACATGGTGAAAAGCCCGTATCTACCAAAAATACAAACGTTAGCTGTGTGTGGTGGTGGCGGCACCTGTAATCCCAGCTACTCGAGAGCCAGAG";
        final String j1 =
                "CATGGAGCCTGACTTTATTTGAAGTAGGGCATTTGCAGATGTATTTAAGATATTTGAGGCCGGGCACAGTGGCTCATGTCTGTAATCCCAGCACTTTAAAAGGCTGAGGCAGGTGTATTCACCTGAGGTCAGGAGTTCGAGACCAGCCTGGCCAACATGGTGAAAAGCCCGTATCTACCAAAAATACAAAAGTTAGCTGGGTGTGGTGGCGGCACCTGTAATCCCAGCTACTCGAGAGCCAGAG";
        final String r0 =
                "CATGGAGCCTGACCTTATTTGAAGTAGGGCATTTGCAGATGTATTTAAGATATTTGAGGCTGGGCACAGTGGCTCACGTCTGTAATCCCAGCACTTTGAAAGGCCGAGGCAGGTGGATTCACCTGAGGTCAGGAGTTTGAGACCAGCCTGTCCCACATGGTGAAAAGCCCGTATCTACCAAAAATACAAACGTTAGCTGTGTGTGGTGGTGGCGGCACCTGTAATCGCAGCTACTCGAGAGCCAGAG";
        final String k25c = "ACCTGTAATCCCAGCTACTCGAGAG";
        final String k25g = "ACCTGTAATCGCAGCTACTCGAGAG";
        try (HcContext ctx = new HcContext(refPath, bamPath, padding);
                CachingIndexedFastaSequenceFile refReader =
                        new CachingIndexedFastaSequenceFile(Paths.get(refPath))) {
            final FinalizedActiveRegion finalized =
                    loadCoveringActiveRegionWithFinalize(
                            ctx, refReader, intervalCli, padding, coverLocus);
            final ReadThreadingAssembler assembler = productionAssemblerFromContext(ctx);
            final int minQual = assembler.getMinBaseQualityToUseInAssembly();
            final int minPrune = assemblerPruneFactor(assembler);
            final int minDangling = assemblerMinDanglingBranchLength(assembler);
            final boolean recoverHeads = assembler.isRecoverDanglingBranches();
            final byte[] graphRef =
                    sliceGraphReference(refReader, finalized.contig, finalized.paddedStart, finalized.paddedEnd);
            final byte[] graphQuals = new byte[graphRef.length];
            Arrays.fill(graphQuals, (byte) 30);
            final RegionAssemblyMaterial graphMaterial =
                    new RegionAssemblyMaterial(
                            finalized.contig,
                            finalized.start,
                            finalized.end,
                            graphRef,
                            graphQuals,
                            finalized.material.readBases,
                            finalized.material.readQuals);
            System.out.println("metric\tvalue");
            System.out.println("region_contig\t" + finalized.contig);
            System.out.println("region_start\t" + finalized.start);
            System.out.println("region_end\t" + finalized.end);
            System.out.println("extended_start\t" + finalized.paddedStart);
            System.out.println("extended_end\t" + finalized.paddedEnd);
            System.out.println("graph_ref_len\t" + graphRef.length);
            System.out.println("finalized_reads\t" + finalized.material.readBases.size());
            final int kmer = 25;
            if (referenceHasNonUniqueKmers(graphRef, kmer)) {
                System.out.println("k25_non_unique_ref\ttrue");
            } else {
                System.out.println("k25_non_unique_ref\tfalse");
            }
            final ReadThreadingGraph rt =
                    prepareReadThreadingGraphForHaplotypeDumpFromRegion(
                            graphMaterial, kmer, minQual, minPrune, minDangling, recoverHeads);
            if (rt == null) {
                System.out.println("rt_status\tnull");
                return;
            }
            System.out.println("rt_nodes\t" + rt.vertexSet().size());
            System.out.println("rt_edges\t" + rt.edgeSet().size());
            MultiDeBruijnVertex cV = null;
            MultiDeBruijnVertex gV = null;
            for (final MultiDeBruijnVertex v : rt.vertexSet()) {
                final String seq = v.getSequenceString();
                if (k25c.equals(seq)) {
                    cV = v;
                } else if (k25g.equals(seq)) {
                    gV = v;
                }
            }
            System.out.println("rt_c25\t" + (cV != null));
            System.out.println("rt_g25\t" + (gV != null));
            dumpRtStar("C25", rt, cV);
            dumpRtStar("G25", rt, gV);
            final SeqGraph seqGraph = rt.toSequenceGraph();
            seqGraph.cleanNonRefPaths();
            dumpJavaSeqMotif("seq_after_to_sequence_graph", seqGraph);
            final String status = cleanupSeqGraphRustParity(seqGraph);
            System.out.println("cleanup_status\t" + status);
            dumpJavaSeqMotif("seq_after_cleanup", seqGraph);
            System.out.println(
                    "seq_nodes\t" + seqGraph.vertexSet().size());
            System.out.println("seq_edges\t" + seqGraph.edgeSet().size());
            System.out.println(
                    "seq_edge\tfrom_seq\tto_seq\tmultiplicity\tis_ref\tout_mult");
            for (final SeqVertex v : seqGraph.vertexSet()) {
                int outMult = 0;
                for (final BaseEdge e : seqGraph.outgoingEdgesOf(v)) {
                    outMult += e.getMultiplicity();
                }
                for (final BaseEdge e : seqGraph.outgoingEdgesOf(v)) {
                    final SeqVertex t = seqGraph.getEdgeTarget(e);
                    if (seqGraph.outDegreeOf(v) < 2) {
                        continue;
                    }
                    System.out.printf(
                            "seq_edge\t%s\t%s\t%d\t%s\t%d%n",
                            v.getSequenceString(),
                            t.getSequenceString(),
                            e.getMultiplicity(),
                            e.isRef(),
                            outMult);
                }
            }
            final SeqVertex source = seqGraph.getReferenceSourceVertex();
            final SeqVertex sink = seqGraph.getReferenceSinkVertex();
            if (source == null || sink == null) {
                System.out.println("kbest\tno_source_or_sink");
                return;
            }
            final List<KBestHaplotype<SeqVertex, BaseEdge>> paths =
                    new ArrayList<>(
                            new GraphBasedKBestHaplotypeFinder<>(seqGraph, source, sink)
                                    .findBestHaplotypes(128));
            System.out.println(
                    "kbest\trank\tscore\tis_ref\tlen\thas_j0\thas_j1\thas_r0");
            int rank = 0;
            Integer j0Rank = null;
            Integer j1Rank = null;
            Integer r0Rank = null;
            Double cutoff = null;
            for (final KBestHaplotype<SeqVertex, BaseEdge> path : paths) {
                final String bases = new String(path.getBases(), StandardCharsets.US_ASCII);
                final boolean hasJ0 = bases.contains(j0);
                final boolean hasJ1 = bases.contains(j1);
                final boolean hasR0 = bases.contains(r0);
                if (hasJ0 && j0Rank == null) {
                    j0Rank = rank;
                }
                if (hasJ1 && j1Rank == null) {
                    j1Rank = rank;
                }
                if (hasR0 && r0Rank == null) {
                    r0Rank = rank;
                }
                if (rank == 127) {
                    cutoff = path.score();
                }
                System.out.printf(
                        "kbest\t%d\t%s\t%s\t%d\t%s\t%s\t%s%n",
                        rank,
                        formatScore(path.score()),
                        path.isReference(),
                        bases.length(),
                        hasJ0,
                        hasJ1,
                        hasR0);
                rank++;
            }
            System.out.println("java_j0_rank\t" + j0Rank);
            System.out.println("java_j1_rank\t" + j1Rank);
            System.out.println("java_r0_rank\t" + r0Rank);
            System.out.println("java_k128_cutoff\t" + cutoff);
        }
    }

    private static byte[] sliceGraphReference(
            final CachingIndexedFastaSequenceFile refReader,
            final String contig,
            final int start1,
            final int end1) {
        return refReader.getSubsequenceAt(contig, start1, end1).getBases();
    }

    private static void dumpRtStar(
            final String tag, final ReadThreadingGraph rt, final MultiDeBruijnVertex v) {
        if (v == null) {
            System.out.println("rt_star\t" + tag + "\tmissing");
            return;
        }
        int outMult = 0;
        for (final MultiSampleEdge e : rt.outgoingEdgesOf(v)) {
            outMult += e.getMultiplicity();
        }
        for (final MultiSampleEdge e : rt.outgoingEdgesOf(v)) {
            final MultiDeBruijnVertex t = rt.getEdgeTarget(e);
            System.out.printf(
                    "rt_star\t%s\t%s\t%s\tmult=%d\tprune=%d\tref=%s\toutMult=%d%n",
                    tag,
                    v.getSequenceString(),
                    t.getSequenceString(),
                    e.getMultiplicity(),
                    e.getPruningMultiplicity(),
                    e.isRef(),
                    outMult);
        }
    }

    private static void dumpJavaSeqMotif(final String stage, final SeqGraph g) {
        System.out.println("motif_stage\t" + stage);
        for (final SeqVertex v : g.vertexSet()) {
            if (g.outDegreeOf(v) < 2) {
                continue;
            }
            boolean hasC = false;
            boolean hasG = false;
            boolean hasA = false;
            boolean hasT = false;
            for (final SeqVertex t : g.outgoingVerticesOf(v)) {
                final String s = t.getSequenceString();
                if (s.startsWith("C")) {
                    hasC = true;
                }
                if (s.startsWith("G")) {
                    hasG = true;
                }
                if (s.startsWith("A")) {
                    hasA = true;
                }
                if (s.startsWith("T")) {
                    hasT = true;
                }
            }
            final String seq = v.getSequenceString();
            final boolean snp = hasC && hasG;
            final boolean at = seq.contains("CCAGCTACTCGAGAG") && hasA && hasT;
            if (!snp && !at) {
                continue;
            }
            int outMult = 0;
            for (final BaseEdge e : g.outgoingEdgesOf(v)) {
                outMult += e.getMultiplicity();
            }
            System.out.printf(
                    "motif\t%s\t%s\tseq=%s\tout=%d\toutMult=%d%n",
                    stage, snp ? "CG_BUBBLE" : "AT_FORK", seq, g.outDegreeOf(v), outMult);
            for (final BaseEdge e : g.outgoingEdgesOf(v)) {
                final SeqVertex t = g.getEdgeTarget(e);
                System.out.printf(
                        "motif_edge\t%s\t%s\t%d\tref=%s\toutMult=%d%n",
                        seq,
                        t.getSequenceString(),
                        e.getMultiplicity(),
                        e.isRef(),
                        outMult);
            }
        }
    }

    private static FinalizedActiveRegion loadCoveringActiveRegionWithFinalize(
            final HcContext ctx,
            final CachingIndexedFastaSequenceFile refReader,
            final String intervalCli,
            final int padding,
            final int coverLocus)
            throws Exception {
        final HaplotypeCallerArgumentCollection hcArgs = hcArgsFromContext(ctx);
        final SampleList samplesList = sampleListFromHeader(ctx.header);
        final List<SimpleInterval> intervals =
                parseIntervals(ctx.header.getSequenceDictionary(), intervalCli);
        for (final List<Locatable> contigIntervals : groupByContig(intervals)) {
            final List<SimpleInterval> contigSimple =
                    contigIntervals.stream()
                            .map(SimpleInterval::new)
                            .collect(Collectors.toList());
            final MultiIntervalLocalReadShard shard =
                    new MultiIntervalLocalReadShard(contigSimple, padding, ctx.readsSource);
            configureHcProductionReadShard(shard, ctx);
            final AssemblyRegionIterator iter =
                    new AssemblyRegionIterator(
                            shard,
                            ctx.header,
                            ctx.reference,
                            null,
                            ctx.engine,
                            ctx.asmArgs,
                            false);
            FinalizedActiveRegion chosen = null;
            while (iter.hasNext()) {
                final AssemblyRegion r = iter.next();
                if (!r.isActive()) {
                    continue;
                }
                if (coverLocus > 0
                        && (r.getStart() > coverLocus || r.getEnd() < coverLocus)) {
                    continue;
                }
                if (chosen != null && coverLocus <= 0 && r.getStart() >= chosen.start) {
                    continue;
                }
                chosen = finalizeActiveRegionSnapshot(hcArgs, samplesList, ctx, refReader, r);
                if (coverLocus > 0) {
                    return chosen;
                }
            }
            if (chosen != null) {
                return chosen;
            }
        }
        throw new IllegalArgumentException("no covering active assembly region in interval");
    }

    private static FinalizedActiveRegion finalizeActiveRegionSnapshot(
            final HaplotypeCallerArgumentCollection hcArgs,
            final SampleList samplesList,
            final HcContext ctx,
            final CachingIndexedFastaSequenceFile refReader,
            final AssemblyRegion r)
            throws Exception {
        final List<GATKRead> rawReads = new ArrayList<>();
        for (final GATKRead read : r.getReads()) {
            rawReads.add(read.copy());
        }
        AssemblyBasedCallerUtils.finalizeRegion(
                r,
                hcArgs.assemblerArgs.errorCorrectReads,
                hcArgs.dontUseSoftClippedBases,
                (byte) (hcArgs.minBaseQualityScore - 1),
                ctx.header,
                samplesList,
                !hcArgs.doNotCorrectOverlappingBaseQualities,
                hcArgs.softClipLowQualityEnds,
                hcArgs.overrideSoftclipFragmentCheck,
                hcArgs.fbargs,
                true);
        final byte[] refBases =
                r.getAssemblyRegionReference(
                        refReader, AssemblyBasedCallerUtils.REFERENCE_PADDING_FOR_ASSEMBLY);
        final byte[] refQuals = new byte[refBases.length];
        Arrays.fill(refQuals, (byte) 30);
        final List<byte[]> readBases = new ArrayList<>();
        final List<byte[]> readQuals = new ArrayList<>();
        final List<GATKRead> finalizedReads = new ArrayList<>();
        for (final GATKRead read : r.getReads()) {
            finalizedReads.add(read);
            readBases.add(Arrays.copyOf(read.getBases(), read.getLength()));
            readQuals.add(Arrays.copyOf(read.getBaseQualities(), read.getLength()));
        }
        final RegionAssemblyMaterial material =
                new RegionAssemblyMaterial(
                        r.getContig(),
                        r.getStart(),
                        r.getEnd(),
                        refBases,
                        refQuals,
                        readBases,
                        readQuals);
        return new FinalizedActiveRegion(
                r.getContig(),
                r.getStart(),
                r.getEnd(),
                r.getPaddedSpan().getStart(),
                r.getPaddedSpan().getEnd(),
                rawReads,
                finalizedReads,
                material);
    }

    private static void assemblyRegionAssemblyStages(final String[] args) throws Exception {
        if (args.length < 3) {
            usage();
        }
        final String refPath = args[0];
        final String bamPath = args[1];
        final String intervalCli = args[2];
        final int padding = args.length > 3 ? parsePadding(args[3]) : DEFAULT_PADDING;
        final RegionAssemblyMaterial material =
                materializeFirstActiveRegion(refPath, bamPath, intervalCli, padding);
        try (HcContext ctx = new HcContext(refPath, bamPath, padding)) {
            emitAssemblyRegionAssemblyStages(material, ctx, padding, "materialize");
        }
    }

    private static void assemblyRegionAssemblyStagesFinalize(final String[] args) throws Exception {
        if (args.length < 3) {
            usage();
        }
        final String refPath = args[0];
        final String bamPath = args[1];
        final String intervalCli = args[2];
        final int padding = args.length > 3 ? parsePadding(args[3]) : DEFAULT_PADDING;
        try (HcContext ctx = new HcContext(refPath, bamPath, padding);
                CachingIndexedFastaSequenceFile refReader =
                        new CachingIndexedFastaSequenceFile(Paths.get(refPath))) {
            final FinalizedActiveRegion finalized =
                    loadFirstActiveRegionWithFinalize(ctx, refReader, intervalCli, padding);
            emitAssemblyRegionAssemblyStages(finalized.material, ctx, padding, "finalize");
        }
    }

    private static void assemblyRegionFinalizeReads(final String[] args) throws Exception {
        if (args.length < 3) {
            usage();
        }
        final String refPath = args[0];
        final String bamPath = args[1];
        final String intervalCli = args[2];
        final int padding = args.length > 3 ? parsePadding(args[3]) : DEFAULT_PADDING;
        try (HcContext ctx = new HcContext(refPath, bamPath, padding);
                CachingIndexedFastaSequenceFile refReader =
                        new CachingIndexedFastaSequenceFile(Paths.get(refPath))) {
            final FinalizedActiveRegion finalized =
                    loadFirstActiveRegionWithFinalize(ctx, refReader, intervalCli, padding);
            System.out.println("region_contig\t" + finalized.contig);
            System.out.println("region_start\t" + finalized.start);
            System.out.println("region_end\t" + finalized.end);
            System.out.println("extended_start\t" + finalized.paddedStart);
            System.out.println("extended_end\t" + finalized.paddedEnd);
            System.out.println("read_path\tfinalize");
            System.out.println("raw_read_count\t" + finalized.rawReads.size());
            System.out.println("finalize_read_count\t" + finalized.finalizedReads.size());
            System.out.println(
                    "read\tqname\tphase\tflags\tmapq\tstart\tend\tcigar\tseq_len\tunclipped_len\tfragment_length\tunmapped");
            final List<GATKRead> rawSorted = new ArrayList<>(finalized.rawReads);
            rawSorted.sort(new ReadCoordinateComparator(ctx.header));
            for (final GATKRead read : rawSorted) {
                printAssemblyFinalizeReadRow(read, "raw");
            }
            final List<GATKRead> finSorted = new ArrayList<>(finalized.finalizedReads);
            finSorted.sort(new ReadCoordinateComparator(ctx.header));
            for (final GATKRead read : finSorted) {
                printAssemblyFinalizeReadRow(read, "finalize");
            }
        }
    }

    /** Active region with raw shard reads and post-{@code finalizeRegion} reads (production assembly input). */
    private static final class FinalizedActiveRegion {
        final String contig;
        final int start;
        final int end;
        final int paddedStart;
        final int paddedEnd;
        final List<GATKRead> rawReads;
        final List<GATKRead> finalizedReads;
        final RegionAssemblyMaterial material;

        FinalizedActiveRegion(
                final String contig,
                final int start,
                final int end,
                final int paddedStart,
                final int paddedEnd,
                final List<GATKRead> rawReads,
                final List<GATKRead> finalizedReads,
                final RegionAssemblyMaterial material) {
            this.contig = contig;
            this.start = start;
            this.end = end;
            this.paddedStart = paddedStart;
            this.paddedEnd = paddedEnd;
            this.rawReads = rawReads;
            this.finalizedReads = finalizedReads;
            this.material = material;
        }
    }

    private static HaplotypeCallerArgumentCollection hcArgsFromContext(final HcContext ctx)
            throws Exception {
        final Field hcArgsField = HaplotypeCallerEngine.class.getDeclaredField("hcArgs");
        hcArgsField.setAccessible(true);
        return (HaplotypeCallerArgumentCollection) hcArgsField.get(ctx.engine);
    }

    private static void printAssemblyFinalizeReadRow(final GATKRead read, final String phase) {
        final boolean unmapped = read.isUnmapped() || read.isEmpty();
        final int start = unmapped ? 0 : read.getStart();
        final int end = unmapped ? 0 : read.getEnd();
        final String cigar = unmapped ? "" : read.getCigar().toString();
        System.out.printf(
                "read\t%s\t%s\t%d\t%d\t%d\t%d\t%s\t%d\t%d\t%d\t%s%n",
                read.getName(),
                phase,
                read.getFlags(),
                read.getMappingQuality(),
                start,
                end,
                cigar,
                read.getLength(),
                AlignmentUtils.unclippedReadLength(read),
                read.getFragmentLength(),
                unmapped);
    }

    private static FinalizedActiveRegion loadFirstActiveRegionWithFinalize(
            final HcContext ctx,
            final CachingIndexedFastaSequenceFile refReader,
            final String intervalCli,
            final int padding)
            throws Exception {
        final HaplotypeCallerArgumentCollection hcArgs = hcArgsFromContext(ctx);
        final SampleList samplesList = sampleListFromHeader(ctx.header);
        final List<SimpleInterval> intervals =
                parseIntervals(ctx.header.getSequenceDictionary(), intervalCli);
        for (final List<Locatable> contigIntervals : groupByContig(intervals)) {
            final List<SimpleInterval> contigSimple =
                    contigIntervals.stream()
                            .map(SimpleInterval::new)
                            .collect(Collectors.toList());
            final MultiIntervalLocalReadShard shard =
                    new MultiIntervalLocalReadShard(contigSimple, padding, ctx.readsSource);
            configureHcProductionReadShard(shard, ctx);
            final AssemblyRegionIterator iter =
                    new AssemblyRegionIterator(
                            shard,
                            ctx.header,
                            ctx.reference,
                            null,
                            ctx.engine,
                            ctx.asmArgs,
                            false);
            FinalizedActiveRegion earliest = null;
            while (iter.hasNext()) {
                final AssemblyRegion r = iter.next();
                if (!r.isActive()) {
                    continue;
                }
                if (earliest != null && r.getStart() >= earliest.start) {
                    continue;
                }
                final List<GATKRead> rawReads = new ArrayList<>();
                for (final GATKRead read : r.getReads()) {
                    rawReads.add(read.copy());
                }
                AssemblyBasedCallerUtils.finalizeRegion(
                        r,
                        hcArgs.assemblerArgs.errorCorrectReads,
                        hcArgs.dontUseSoftClippedBases,
                        (byte) (hcArgs.minBaseQualityScore - 1),
                        ctx.header,
                        samplesList,
                        !hcArgs.doNotCorrectOverlappingBaseQualities,
                        hcArgs.softClipLowQualityEnds,
                        hcArgs.overrideSoftclipFragmentCheck,
                        hcArgs.fbargs,
                        true);
                final byte[] refBases =
                        r.getAssemblyRegionReference(
                                refReader,
                                AssemblyBasedCallerUtils.REFERENCE_PADDING_FOR_ASSEMBLY);
                final byte[] refQuals = new byte[refBases.length];
                Arrays.fill(refQuals, (byte) 30);
                final List<byte[]> readBases = new ArrayList<>();
                final List<byte[]> readQuals = new ArrayList<>();
                final List<GATKRead> finalizedReads = new ArrayList<>();
                for (final GATKRead read : r.getReads()) {
                    finalizedReads.add(read);
                    readBases.add(Arrays.copyOf(read.getBases(), read.getLength()));
                    readQuals.add(Arrays.copyOf(read.getBaseQualities(), read.getLength()));
                }
                final RegionAssemblyMaterial material =
                        new RegionAssemblyMaterial(
                                r.getContig(),
                                r.getStart(),
                                r.getEnd(),
                                refBases,
                                refQuals,
                                readBases,
                                readQuals);
                earliest =
                        new FinalizedActiveRegion(
                                r.getContig(),
                                r.getStart(),
                                r.getEnd(),
                                r.getPaddedSpan().getStart(),
                                r.getPaddedSpan().getEnd(),
                                rawReads,
                                finalizedReads,
                                material);
            }
            if (earliest != null) {
                return earliest;
            }
        }
        throw new IllegalArgumentException("no active assembly region in interval");
    }

    private static void emitAssemblyRegionAssemblyStages(
            final RegionAssemblyMaterial material,
            final HcContext ctx,
            final int padding,
            final String readPath)
            throws Exception {
        final int kmer = 85;
        final ReadThreadingAssembler assembler = productionAssemblerFromContext(ctx);
        final int minQual = assembler.getMinBaseQualityToUseInAssembly();
        final int minPrune = assemblerPruneFactor(assembler);
        final int minDangling = assemblerMinDanglingBranchLength(assembler);
        final boolean recoverHeads = assembler.isRecoverDanglingBranches();
        System.out.println("read_path\t" + readPath);
        System.out.println("region_contig\t" + material.contig);
        System.out.println("region_start\t" + material.start);
        System.out.println("region_end\t" + material.end);
        System.out.println("padded_ref_len\t" + material.refBases.length);
        System.out.println("read_count\t" + material.readBases.size());

        final List<AssemblyStageRow> rows = new ArrayList<>();
        if (material.refBases.length >= kmer) {
            ReadThreadingGraph graph =
                    buildReadThreadingGraphFromRegion(material, kmer, minQual);
            if (graph.hasCycles()) {
                System.out.println("warn\tgraph_has_cycles\ttrue");
            }
            rows.add(
                    rtStageRow("threading_after_build", graph, material.refBases, 128));
            final ChainPruner<MultiDeBruijnVertex, MultiSampleEdge> pruner =
                    makeChainPruner(minPrune, false);
            pruner.pruneLowWeightChains(graph);
            rows.add(
                    rtStageRow(
                            "threading_after_prune_before_dangling",
                            graph,
                            material.refBases,
                            128));
            int tailsAttempted = countDanglingTailCandidates(graph);
            int headsAttempted = 0;
            if (recoverHeads) {
                for (final MultiDeBruijnVertex v : graph.vertexSet()) {
                    if (graph.inDegreeOf(v) == 0 && graph.outDegreeOf(v) > 0) {
                        headsAttempted++;
                    }
                }
            }
            final int edgesBefore = graph.edgeSet().size();
            final SmithWatermanAligner aligner =
                    SmithWatermanAligner.getAligner(SmithWatermanAligner.Implementation.JAVA);
            final SWParameters swParams = danglingEndSwParameters();
            if (graph.getReferenceSourceVertex() != null) {
                graph.recoverDanglingTails(minPrune, minDangling, false, aligner, swParams);
                if (recoverHeads) {
                    graph.recoverDanglingHeads(minPrune, minDangling, false, aligner, swParams);
                }
            }
            final int edgesAfter = graph.edgeSet().size();
            System.out.printf(
                    "dangling_recovery\ttails_attempted=%d\ttails_recovered=unknown\theads_attempted=%d\theads_recovered=unknown\tedges_before=%d\tedges_after=%d%n",
                    tailsAttempted,
                    headsAttempted,
                    edgesBefore,
                    edgesAfter);
            rows.add(
                    rtStageRow("threading_after_dangling", graph, material.refBases, 128));
            if (graph.getReferenceSourceVertex() != null
                    && graph.getReferenceSinkVertex() != null) {
                graph.removePathsNotConnectedToRef();
            }
            if (graph.getReferenceSourceVertex() != null
                    && graph.getReferenceSinkVertex() != null) {
                rows.add(
                        rtStageRow(
                                "threading_after_prune_dangling",
                                graph,
                                material.refBases,
                                128));
                SeqGraph seq = graph.toSequenceGraph();
                seq.cleanNonRefPaths();
                rows.add(
                        seqStageRow(
                                "after_to_sequence_graph", seq, material.refBases, 128));
                seq.zipLinearChains();
                seq.removeSingletonOrphanVertices();
                seq.removeVerticesNotConnectedToRefRegardlessOfEdgeDirection();
                rows.add(
                        seqStageRow(
                                "after_zip_orphans_prune", seq, material.refBases, 128));
                seq.simplifyGraph();
                rows.add(
                        seqStageRow("after_first_simplify", seq, material.refBases, 128));
                if (seq.getReferenceSourceVertex() != null
                        && seq.getReferenceSinkVertex() != null) {
                    seq.removePathsNotConnectedToRef();
                    seq.simplifyGraph();
                    rows.add(
                            seqStageRow(
                                    "after_remove_paths_final_simplify",
                                    seq,
                                    material.refBases,
                                    128));
                }
            }
        }

        System.out.println(
                "graph\tstage\tnodes\tedges\tref_spine_vertices\tbranch_vertices\tbranch_vertices_all\tnon_ref_edges_on_spine\tnon_ref_edges_all\tkbest_paths\textracted_haps\tnon_ref_haps\ttop_path_len\ttop_path_eq_ref");
        for (final AssemblyStageRow r : rows) {
            printAssemblyStageRow(r);
        }
    }
}
