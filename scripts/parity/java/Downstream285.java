import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.utils.fasta.CachingIndexedFastaSequenceFile;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypeLikelihoods;
import htsjdk.variant.variantcontext.VariantContext;
import java.io.BufferedReader;
import java.io.FileReader;
import java.lang.reflect.Field;
import java.lang.reflect.Method;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Locale;
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
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCallerGenotypingEngine;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.PairHMMNativeArgumentCollection;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.genotyper.IndexedAlleleList;
import org.broadinstitute.hellbender.utils.genotyper.LikelihoodMatrix;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.pairhmm.PairHMM;
import org.broadinstitute.hellbender.utils.read.GATKRead;

/**
 * 6R.285: one live GATK 4.4 genotyping pass at 20:29455649. Prints the
 * continuous genotype-likelihood vector before {@code getAsPLs}.
 */
public final class Downstream285 {
    private static final int LOC = 29455649;
    private static final int PADDING = 100;

    public static void main(final String[] args) throws Exception {
        if (args.length != 3) {
            throw new IllegalArgumentException("usage: Downstream285 ref bam frozen_inputs.tsv");
        }
        final List<String> frozenHaps = loadHaps(args[2]);
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
        hcArgs.likelihoodArgs.pairHMM = PairHMM.Implementation.AVX_LOGLESS_CACHING;
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
        install(engine, frozenHaps);
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
    }

    private static void install(final HaplotypeCallerEngine engine, final List<String> frozenHaps)
            throws Exception {
        final Field hcArgsF = HaplotypeCallerEngine.class.getDeclaredField("hcArgs");
        hcArgsF.setAccessible(true);
        final HaplotypeCallerArgumentCollection hcArgs =
                (HaplotypeCallerArgumentCollection) hcArgsF.get(engine);
        final Field samplesF = HaplotypeCallerEngine.class.getDeclaredField("samplesList");
        samplesF.setAccessible(true);
        final SampleList samples = (SampleList) samplesF.get(engine);
        final Dump dump =
                new Dump(
                        hcArgs,
                        samples,
                        !hcArgs.doNotRunPhysicalPhasing,
                        hcArgs.applyBQD,
                        frozenHaps);
        final Field ae = HaplotypeCallerEngine.class.getDeclaredField("annotationEngine");
        ae.setAccessible(true);
        dump.setAnnotationEngine(
                (VariantAnnotatorEngine) ae.get(engine));
        final Field ge = HaplotypeCallerEngine.class.getDeclaredField("genotypingEngine");
        ge.setAccessible(true);
        ge.set(engine, dump);
    }

    private static List<String> loadHaps(final String path) throws Exception {
        final List<String> haps = new ArrayList<String>();
        try (BufferedReader br = new BufferedReader(new FileReader(path))) {
            String line;
            while ((line = br.readLine()) != null) {
                if (line.startsWith("HAP\t")) {
                    final String[] f = line.split("\t", -1);
                    haps.add(f[3]);
                }
            }
        }
        if (haps.size() != 5) {
            throw new IllegalStateException("expected 5 frozen haplotypes, got " + haps.size());
        }
        return haps;
    }

    private static void kv(final String key, final String value) {
        System.out.println("6R285\t" + key + "\t" + value);
    }

    private static String bits(final double x) {
        return String.format(
                Locale.US, "%.17g bits=0x%016x", x, Double.doubleToRawLongBits(x));
    }

    private static final class Dump extends HaplotypeCallerGenotypingEngine {
        private final List<String> frozenHaps;

        Dump(
                final HaplotypeCallerArgumentCollection configuration,
                final SampleList samples,
                final boolean doPhysicalPhasing,
                final boolean applyBQD,
                final List<String> frozenHaps) {
            super(configuration, samples, doPhysicalPhasing, applyBQD);
            this.frozenHaps = frozenHaps;
        }

        @Override
        @SuppressWarnings({"rawtypes", "unchecked"})
        public org.broadinstitute.hellbender.tools.walkers.haplotypecaller.CalledHaplotypes
                assignGenotypeLikelihoods(
                        final List haplotypes,
                        final AlleleLikelihoods readLikelihoods,
                        final java.util.Map perSampleFilteredReadList,
                        final byte[] ref,
                        final SimpleInterval refLoc,
                        final SimpleInterval activeRegionWindow,
                        final FeatureContext tracker,
                        final List givenAlleles,
                        final boolean emitReferenceConfidence,
                        final int maxMnpDistance,
                        final htsjdk.samtools.SAMFileHeader header,
                        final boolean withBamOut,
                        final java.util.Set suspiciousLocations,
                        final AlleleLikelihoods preFilteringAlleleLikelihoods) {
            dumpFrozenColumns(haplotypes, readLikelihoods);
            return super.assignGenotypeLikelihoods(
                    haplotypes,
                    readLikelihoods,
                    perSampleFilteredReadList,
                    ref,
                    refLoc,
                    activeRegionWindow,
                    tracker,
                    givenAlleles,
                    emitReferenceConfidence,
                    maxMnpDistance,
                    header,
                    withBamOut,
                    suspiciousLocations,
                    preFilteringAlleleLikelihoods);
        }

        @SuppressWarnings({"rawtypes", "unchecked"})
        private void dumpFrozenColumns(final List haplotypes, final AlleleLikelihoods readLikelihoods) {
            final int[] index = new int[frozenHaps.size()];
            Arrays.fill(index, -1);
            for (int h = 0; h < haplotypes.size(); h++) {
                final byte[] bases = ((Haplotype) haplotypes.get(h)).getBases();
                final String seq = new String(bases);
                for (int k = 0; k < frozenHaps.size(); k++) {
                    if (seq.equals(frozenHaps.get(k))) {
                        index[k] = h;
                    }
                }
            }
            kv("frozen_hap_index", Arrays.toString(index));
            kv("n_haplotypes", Integer.toString(haplotypes.size()));
            if (readLikelihoods.numberOfSamples() < 1) {
                return;
            }
            final LikelihoodMatrix matrix = readLikelihoods.sampleMatrix(0);
            int matched = 0;
            for (int k = 0; k < index.length; k++) {
                if (index[k] >= 0) {
                    matched++;
                }
            }
            kv("frozen_haps_found", Integer.toString(matched));
            if (matched == 0) {
                return;
            }
            for (int r = 0; r < matrix.evidenceCount(); r++) {
                final GATKRead ev = (GATKRead) matrix.getEvidence(r);
                final StringBuilder row = new StringBuilder();
                row.append(ev.getName()).append('\t').append(ev.getFlags());
                for (int k = 0; k < index.length; k++) {
                    row.append('\t');
                    if (index[k] < 0) {
                        row.append("ABSENT");
                    } else {
                        row.append(bits(matrix.get(index[k], r)));
                    }
                }
                kv("hap_ll", row.toString());
            }
        }

        @Override
        @SuppressWarnings({"rawtypes", "unchecked"})
        protected htsjdk.variant.variantcontext.GenotypesContext calculateGLsForThisEvent(
                final AlleleLikelihoods readLikelihoods,
                final VariantContext mergedVC,
                final List noCallAlleles,
                final byte[] paddedReference,
                final int offsetForRefIntoEvent,
                final org.broadinstitute.hellbender.utils.dragstr.DragstrReferenceAnalyzer dragstrs) {
            if (mergedVC.getStart() == LOC) {
                dumpEvent(readLikelihoods, mergedVC, paddedReference, offsetForRefIntoEvent, dragstrs);
            }
            return super.calculateGLsForThisEvent(
                    readLikelihoods,
                    mergedVC,
                    noCallAlleles,
                    paddedReference,
                    offsetForRefIntoEvent,
                    dragstrs);
        }

        @SuppressWarnings({"rawtypes", "unchecked"})
        private void dumpEvent(
                final AlleleLikelihoods readLikelihoods,
                final VariantContext mergedVC,
                final byte[] paddedReference,
                final int offsetForRefIntoEvent,
                final org.broadinstitute.hellbender.utils.dragstr.DragstrReferenceAnalyzer dragstrs) {
            try {
                kv(
                        "event",
                        mergedVC.getContig()
                                + ":"
                                + mergedVC.getStart()
                                + " "
                                + mergedVC.getReference().getBaseString()
                                + " "
                                + mergedVC.getAlternateAlleles());
                final LikelihoodMatrix matrix = readLikelihoods.sampleMatrix(0);
                final StringBuilder cols = new StringBuilder();
                for (int a = 0; a < matrix.numberOfAlleles(); a++) {
                    if (a > 0) {
                        cols.append(',');
                    }
                    cols.append(((Allele) matrix.getAllele(a)).getBaseString());
                }
                kv("allele_columns", cols.toString());
                kv("n_evidence", Integer.toString(matrix.evidenceCount()));
                final double[] sum = new double[matrix.numberOfAlleles()];
                for (int r = 0; r < matrix.evidenceCount(); r++) {
                    final GATKRead ev = (GATKRead) matrix.getEvidence(r);
                    final StringBuilder row = new StringBuilder();
                    row.append(ev.getName()).append('\t').append(ev.getFlags());
                    for (int a = 0; a < matrix.numberOfAlleles(); a++) {
                        final double v = matrix.get(a, r);
                        sum[a] += v;
                        row.append('\t').append(bits(v));
                    }
                    kv("allele_ll", row.toString());
                }
                for (int a = 0; a < sum.length; a++) {
                    kv("allele_sum", a + "\t" + cols.toString().split(",")[a] + "\t" + bits(sum[a]));
                }
                final Field modelF =
                        HaplotypeCallerGenotypingEngine.class.getDeclaredField("genotypingModel");
                modelF.setAccessible(true);
                final Object model = modelF.get(this);
                final Field ploidyF =
                        HaplotypeCallerGenotypingEngine.class.getDeclaredField("ploidyModel");
                ploidyF.setAccessible(true);
                final Object ploidy = ploidyF.get(this);
                final List vcAlleles = mergedVC.getAlleles();
                final Object alleleList =
                        readLikelihoods.numberOfAlleles() == vcAlleles.size()
                                ? readLikelihoods
                                : new IndexedAlleleList<Allele>(vcAlleles);
                final Class<?> dataClass =
                        Class.forName(
                                "org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypingData");
                final Object data =
                        dataClass
                                .getConstructor(
                                        Class.forName(
                                                "org.broadinstitute.hellbender.tools.walkers.genotyper.PloidyModel"),
                                        AlleleLikelihoods.class)
                                .newInstance(ploidy, readLikelihoods);
                final Method calc =
                        model.getClass()
                                .getMethod(
                                        "calculateLikelihoods",
                                        Class.forName(
                                                "org.broadinstitute.hellbender.utils.genotyper.AlleleList"),
                                        dataClass,
                                        byte[].class,
                                        int.class,
                                        org.broadinstitute.hellbender.utils.dragstr.DragstrReferenceAnalyzer
                                                .class);
                final Object likelihoods =
                        calc.invoke(
                                model,
                                alleleList,
                                data,
                                paddedReference,
                                Integer.valueOf(offsetForRefIntoEvent),
                                dragstrs);
                final Method sampleLikelihoods =
                        likelihoods.getClass().getMethod("sampleLikelihoods", int.class);
                final GenotypeLikelihoods gls =
                        (GenotypeLikelihoods) sampleLikelihoods.invoke(likelihoods, Integer.valueOf(0));
                final double[] vector = gls.getAsVector();
                double max = Double.NEGATIVE_INFINITY;
                for (double v : vector) {
                    max = Math.max(max, v);
                }
                final StringBuilder gl = new StringBuilder();
                final StringBuilder cont = new StringBuilder();
                final int[] pl = gls.getAsPLs();
                final StringBuilder pli = new StringBuilder();
                for (int i = 0; i < vector.length; i++) {
                    if (i > 0) {
                        gl.append('\t');
                        cont.append('\t');
                        pli.append(',');
                    }
                    gl.append(bits(vector[i]));
                    cont.append(bits(-10.0 * (vector[i] - max)));
                    pli.append(pl[i]);
                }
                kv("java_gl_vector", gl.toString());
                kv("java_gl_max", bits(max));
                kv("java_continuous_pl_vector", cont.toString());
                kv("java_integer_pl", pli.toString());
                kv(
                        "pl_transform",
                        "GLsToPLs: adjust=max(log10 GL); PL[i]=Math.round(min(-10*(GL[i]-adjust), Integer.MAX_VALUE)); double in, int out; no other clamp");
            } catch (final Exception e) {
                throw new IllegalStateException(e);
            }
        }
    }
}
