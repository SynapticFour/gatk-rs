import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeLikelihoods;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import java.lang.reflect.Field;
import java.lang.reflect.Method;
import java.util.Arrays;
import java.util.Collections;
import java.util.IdentityHashMap;
import java.util.List;
import java.util.Map;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.tools.walkers.annotator.VariantAnnotatorEngine;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypingEngine;
import org.broadinstitute.hellbender.tools.walkers.genotyper.OutputMode;
import org.broadinstitute.hellbender.tools.walkers.genotyper.afcalc.AFCalculationResult;
import org.broadinstitute.hellbender.tools.walkers.genotyper.afcalc.AlleleFrequencyCalculator;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AssemblyBasedCallerUtils;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.CalledHaplotypes;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCallerArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCallerEngine;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCallerGenotypingEngine;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.genotyper.LikelihoodMatrix;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;
import org.broadinstitute.hellbender.utils.haplotype.EventMap;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;

/**
 * TEST-ONLY: dump {@code calculateGLsForThisEvent} + {@code calculateGenotypes} inputs at one loc.
 * Does not change genotyping arithmetic: delegates to production {@code assignGenotypeLikelihoods}.
 */
public final class HcParityGenotypeEmitDump extends HaplotypeCallerGenotypingEngine {

    private static final double AFC_PASSES_THRESHOLD_EPSILON = 1.0e-10;

    private final int targetStart;

    public static void installOn(final HaplotypeCallerEngine engine, final int targetStart)
            throws Exception {
        final Field hcArgsF = HaplotypeCallerEngine.class.getDeclaredField("hcArgs");
        hcArgsF.setAccessible(true);
        final HaplotypeCallerArgumentCollection hcArgs =
                (HaplotypeCallerArgumentCollection) hcArgsF.get(engine);
        final Field samplesF = HaplotypeCallerEngine.class.getDeclaredField("samplesList");
        samplesF.setAccessible(true);
        final SampleList samples = (SampleList) samplesF.get(engine);
        final HcParityGenotypeEmitDump dump =
                new HcParityGenotypeEmitDump(
                        hcArgs,
                        samples,
                        !hcArgs.doNotRunPhysicalPhasing,
                        hcArgs.applyBQD,
                        targetStart);
        final Field ae = HaplotypeCallerEngine.class.getDeclaredField("annotationEngine");
        ae.setAccessible(true);
        dump.setAnnotationEngine((VariantAnnotatorEngine) ae.get(engine));
        final Field ge = HaplotypeCallerEngine.class.getDeclaredField("genotypingEngine");
        ge.setAccessible(true);
        ge.set(engine, dump);
    }

    private HcParityGenotypeEmitDump(
            final HaplotypeCallerArgumentCollection configuration,
            final SampleList samples,
            final boolean doPhysicalPhasing,
            final boolean applyBQD,
            final int targetStart) {
        super(configuration, samples, doPhysicalPhasing, applyBQD);
        this.targetStart = targetStart;
    }

    @Override
    @SuppressWarnings({"rawtypes", "unchecked"})
    public CalledHaplotypes assignGenotypeLikelihoods(
            final List haplotypes,
            final AlleleLikelihoods readLikelihoods,
            final Map perSampleFilteredReadList,
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
        dumpAtLoc(haplotypes, readLikelihoods, ref, refLoc, activeRegionWindow, header, maxMnpDistance);
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
    private void dumpAtLoc(
            final List haplotypes,
            final AlleleLikelihoods readLikelihoods,
            final byte[] ref,
            final SimpleInterval refLoc,
            final SimpleInterval activeRegionWindow,
            final htsjdk.samtools.SAMFileHeader header,
            final int maxMnpDistance) {
        try {
            EventMap.buildEventMapsForHaplotypes(
                    haplotypes, ref, refLoc, false, maxMnpDistance);
            kv("target", Integer.toString(targetStart));
            kv(
                    "active_window",
                    activeRegionWindow.getContig()
                            + ":"
                            + activeRegionWindow.getStart()
                            + "-"
                            + activeRegionWindow.getEnd());
            kv("hap_count", Integer.toString(haplotypes.size()));
            kv(
                    "stand_call_conf",
                    Double.toString(configuration.genotypeArgs.standardConfidenceForCalling));
            kv("output_mode", configuration.outputMode.toString());
            kv(
                    "genotype_assignment_method",
                    configuration.genotypeArgs.genotypeAssignmentMethod.toString());

            if (targetStart < activeRegionWindow.getStart()
                    || targetStart > activeRegionWindow.getEnd()) {
                kv("in_active_window", "false");
                return;
            }
            kv("in_active_window", "true");

            final List events =
                    AssemblyBasedCallerUtils.getVariantContextsFromActiveHaplotypes(
                            targetStart, haplotypes, true);
            kv("events_with_spanning_n", Integer.toString(events.size()));
            for (int i = 0; i < events.size(); i++) {
                final VariantContext vc = (VariantContext) events.get(i);
                kv(
                        "event_" + i,
                        vc.getStart()
                                + "\t"
                                + vc.getReference()
                                + "\t"
                                + vc.getAlternateAlleles()
                                + "\tend="
                                + vc.getEnd());
            }
            final Allele refAllele =
                    Allele.create(new byte[] {ref[targetStart - refLoc.getStart()]}, true);
            final List replaced =
                    HaplotypeCallerGenotypingEngine.replaceSpanDels(
                            events, refAllele, targetStart);
            kv("after_replace_span_dels_n", Integer.toString(replaced.size()));
            for (int i = 0; i < replaced.size(); i++) {
                final VariantContext vc = (VariantContext) replaced.get(i);
                kv(
                        "replaced_" + i,
                        vc.getStart() + "\t" + vc.getReference() + "\t" + vc.getAlternateAlleles());
            }
            VariantContext merged =
                    AssemblyBasedCallerUtils.makeMergedVariantContext(replaced);
            if (merged == null) {
                kv("merged_vc", "null");
                return;
            }
            kv("merged_alleles", merged.getAlleles().toString());
            kv("merged_n_alleles", Integer.toString(merged.getNAlleles()));

            kv(
                    "hap_ll_n_alleles",
                    Integer.toString(readLikelihoods.numberOfAlleles()));
            kv(
                    "hap_ll_n_evidence",
                    Integer.toString(readLikelihoods.evidenceCount()));
            kv(
                    "hap_ll_sample0_evidence",
                    Integer.toString(readLikelihoods.sampleEvidenceCount(0)));

            final Map alleleMapper =
                    AssemblyBasedCallerUtils.createAlleleMapper(
                            merged, targetStart, haplotypes, true);
            dumpAlleleMapper(alleleMapper, haplotypes, targetStart);

            AlleleLikelihoods readAlleleLikelihoods = readLikelihoods.marginalize(alleleMapper);
            final SimpleInterval relevant =
                    new SimpleInterval(merged)
                            .expandWithinContig(
                                    2, header.getSequenceDictionary());
            readAlleleLikelihoods.retainEvidence(
                    (Object ev) -> relevant.overlaps((GATKRead) ev));
            dumpReadAlleleLikelihoods(readAlleleLikelihoods, merged);

            final int ploidy = configuration.genotypeArgs.samplePloidy;
            final List noCallAlleles = GATKVariantContextUtils.noCallAlleles(ploidy);
            final Method calcGLs =
                    HaplotypeCallerGenotypingEngine.class.getDeclaredMethod(
                            "calculateGLsForThisEvent",
                            AlleleLikelihoods.class,
                            VariantContext.class,
                            List.class,
                            byte[].class,
                            int.class,
                            org.broadinstitute.hellbender.utils.dragstr.DragstrReferenceAnalyzer
                                    .class);
            calcGLs.setAccessible(true);
            final htsjdk.variant.variantcontext.GenotypesContext genotypes =
                    (htsjdk.variant.variantcontext.GenotypesContext)
                            calcGLs.invoke(
                                    this,
                                    readAlleleLikelihoods,
                                    merged,
                                    noCallAlleles,
                                    ref,
                                    targetStart - refLoc.getStart(),
                                    null);
            final Genotype g0 = genotypes.get(0);
            kv("gl_sample", g0.getSampleName());
            kv("gl_called_before_calculateGenotypes", Boolean.toString(g0.isCalled()));
            kv("gl_alleles_before", g0.getAlleles().toString());
            kv("gl_has_pl", Boolean.toString(g0.hasLikelihoods()));
            final int[] pl = g0.hasLikelihoods() ? g0.getPL() : new int[0];
            kv("pl", ints(pl));
            kv("n_genotypes", Integer.toString(pl.length));
            final int bestPl = minIndex(pl);
            kv("best_pl_index", Integer.toString(bestPl));
            kv("assigned_gt_from_pl", diploidGtFromPlIndex(merged.getNAlleles(), bestPl));
            if (g0.hasLikelihoods()) {
                final GenotypeLikelihoods gls = g0.getLikelihoods();
                kv("log10_gl", doubles(gls.getAsVector()));
            }
            kv("gq_before", g0.hasGQ() ? Integer.toString(g0.getGQ()) : "ABSENT");

            final VariantContext withGls =
                    new VariantContextBuilder(merged).genotypes(genotypes).make();

            final Field afcField =
                    GenotypingEngine.class.getDeclaredField("alleleFrequencyCalculator");
            afcField.setAccessible(true);
            final AlleleFrequencyCalculator afc =
                    (AlleleFrequencyCalculator) afcField.get(this);
            final AFCalculationResult af =
                    afc.calculate(withGls, configuration.genotypeArgs.samplePloidy);
            kv("log10_prob_only_ref", Double.toString(af.log10ProbOnlyRefAlleleExists()));
            kv("log10_prob_variant_present", Double.toString(af.log10ProbVariantPresent()));
            final double standCall = configuration.genotypeArgs.standardConfidenceForCalling;
            final double callConfLog10 = QualityUtils.qualToErrorProbLog10(standCall);
            kv("call_conf_log10", Double.toString(callConfLog10));

            boolean siteIsMonomorphic = true;
            final StringBuilder subset = new StringBuilder();
            for (final Allele allele : (List<Allele>) af.getAllelesUsedInGenotyping()) {
                if (allele.isReference()) {
                    continue;
                }
                final double log10Absent = af.getLog10PosteriorOfAlleleAbsent(allele);
                final boolean isPlausible =
                        log10Absent + AFC_PASSES_THRESHOLD_EPSILON < callConfLog10;
                kv(
                        "alt_threshold",
                        allele
                                + "\tlog10_absent="
                                + log10Absent
                                + "\tpassesThreshold30="
                                + isPlausible);
                siteIsMonomorphic &= !isPlausible;
                if (isPlausible) {
                    if (subset.length() > 0) {
                        subset.append(",");
                    }
                    subset.append(allele);
                }
            }
            kv("calculateOutputAlleleSubset_alts", subset.length() == 0 ? "[]" : subset.toString());
            kv("siteIsMonomorphic", Boolean.toString(siteIsMonomorphic));

            final double log10Confidence =
                    !siteIsMonomorphic
                            ? af.log10ProbOnlyRefAlleleExists() + 0.0
                            : af.log10ProbVariantPresent() + 0.0;
            final double phredScaled = (-10.0 * log10Confidence) + 0.0;
            kv("phred_scaled_confidence", Double.toString(phredScaled));
            final boolean emitAllConfident =
                    configuration.outputMode == OutputMode.EMIT_ALL_CONFIDENT_SITES;
            final boolean passesCall = phredScaled >= standCall;
            final boolean passesEmit =
                    (emitAllConfident || !siteIsMonomorphic) && passesCall;
            kv("passesCallThreshold", Boolean.toString(passesCall));
            kv("passesEmitThreshold", Boolean.toString(passesEmit));

            final VariantContext call = calculateGenotypes(withGls);
            if (call == null) {
                kv("calculateGenotypes", "null");
            } else {
                kv("calculateGenotypes_alleles", call.getAlleles().toString());
                kv("calculateGenotypes_qual", Double.toString(call.getPhredScaledQual()));
                final Genotype cg = call.getGenotype(0);
                kv("calculateGenotypes_gt", cg.getGenotypeString());
                kv("calculateGenotypes_pl", cg.hasPL() ? ints(cg.getPL()) : "ABSENT");
                kv("calculateGenotypes_gq", cg.hasGQ() ? Integer.toString(cg.getGQ()) : "ABSENT");
                kv("calculateGenotypes_ad", cg.hasAD() ? ints(cg.getAD()) : "ABSENT");
            }
        } catch (final Exception e) {
            kv("dump_error", e.toString());
            e.printStackTrace(System.err);
        }
    }

    @SuppressWarnings({"rawtypes", "unchecked"})
    private static void dumpAlleleMapper(
            final Map alleleMapper, final List haplotypes, final int loc) {
        kv106("mapper_n_keys", Integer.toString(alleleMapper.size()));
        final StringBuilder keys = new StringBuilder();
        for (final Object k : alleleMapper.keySet()) {
            final Allele a = (Allele) k;
            if (keys.length() > 0) {
                keys.append(",");
            }
            keys.append(alleleLabel(a));
        }
        kv106("mapper_keys", keys.toString());

        final IdentityHashMap mapped = new IdentityHashMap<>();
        for (final Object e : alleleMapper.entrySet()) {
            final Map.Entry entry = (Map.Entry) e;
            final Allele a = (Allele) entry.getKey();
            final String label = alleleLabel(a);
            final List haps = (List) entry.getValue();
            kv106("mapper_pool", label + "\tn=" + haps.size());
            for (final Object ho : haps) {
                final Haplotype h = (Haplotype) ho;
                final String prev = (String) mapped.get(h);
                mapped.put(h, prev == null ? label : prev + "," + label);
            }
        }

        int nRef = 0;
        int nT = 0;
        int nStar = 0;
        int nUnmapped = 0;
        int nMulti = 0;
        for (int i = 0; i < haplotypes.size(); i++) {
            final Haplotype h = (Haplotype) haplotypes.get(i);
            final String mappedAllele = (String) mapped.get(h);
            final String assignment = mappedAllele == null ? "unmapped" : mappedAllele;
            if (mappedAllele == null) {
                nUnmapped++;
            } else {
                final String[] parts = mappedAllele.split(",");
                if (parts.length > 1) {
                    nMulti++;
                }
                for (final String p : parts) {
                    if (p.endsWith("*") && !p.equals("*")) {
                        nRef++;
                    } else if (p.equals("T")) {
                        nT++;
                    } else if (p.equals("*")) {
                        nStar++;
                    }
                }
            }
            final List spanning;
            if (h.getEventMap() == null) {
                spanning = Collections.emptyList();
            } else {
                spanning = h.getEventMap().getOverlappingEvents(loc);
            }
            final StringBuilder evs = new StringBuilder();
            for (int j = 0; j < spanning.size(); j++) {
                final VariantContext vc = (VariantContext) spanning.get(j);
                if (j > 0) {
                    evs.append(";");
                }
                evs.append(vc.getStart())
                        .append(":")
                        .append(vc.getReference().getBaseString())
                        .append("/")
                        .append(vc.getAlternateAlleles());
            }
            kv106(
                    "hap",
                    i
                            + "\t"
                            + fnv1a64Hex(h.getBases())
                            + "\tlen="
                            + h.getBases().length
                            + "\tisRef="
                            + h.isReference()
                            + "\tmapped="
                            + assignment
                            + "\tspanning_n="
                            + spanning.size()
                            + "\tspanning="
                            + (evs.length() == 0 ? "." : evs.toString()));
        }
        kv106("pool_C_ref", Integer.toString(nRef));
        kv106("pool_T_alt", Integer.toString(nT));
        kv106("pool_star", Integer.toString(nStar));
        kv106("pool_unmapped", Integer.toString(nUnmapped));
        kv106("pool_multi", Integer.toString(nMulti));
        kv106("hap_count", Integer.toString(haplotypes.size()));
    }

    private static String alleleLabel(final Allele a) {
        if (a.isSymbolic() && a.getDisplayString().equals("*")) {
            return "*";
        }
        return a.getBaseString() + (a.isReference() ? "*" : "");
    }

    @SuppressWarnings({"rawtypes", "unchecked"})
    private static void dumpReadAlleleLikelihoods(
            final AlleleLikelihoods likelihoods, final VariantContext merged) {
        kv106("allele_ll_n_alleles", Integer.toString(likelihoods.numberOfAlleles()));
        kv106("allele_ll_n_evidence", Integer.toString(likelihoods.evidenceCount()));
        kv106("allele_ll_sample0_evidence", Integer.toString(likelihoods.sampleEvidenceCount(0)));
        kv106("allele_ll_sample_count", Integer.toString(likelihoods.numberOfSamples()));
        kv106("merged_alleles_for_gls", merged.getAlleles().toString());
        kv106(
                "gls_uses_likelihood_allele_list",
                Boolean.toString(likelihoods.numberOfAlleles() == merged.getNAlleles()));
        final StringBuilder cols = new StringBuilder();
        for (final Object a : likelihoods.alleles()) {
            if (cols.length() > 0) {
                cols.append(",");
            }
            cols.append(alleleLabel((Allele) a));
        }
        kv106("allele_ll_columns", cols.toString());
        if (likelihoods.numberOfSamples() < 1) {
            return;
        }
        final LikelihoodMatrix matrix = likelihoods.sampleMatrix(0);
        final int nEv = matrix.evidenceCount();
        final int nAl = matrix.numberOfAlleles();
        kv106("matrix_evidence", Integer.toString(nEv));
        kv106("matrix_alleles", Integer.toString(nAl));
        for (int r = 0; r < nEv; r++) {
            final GATKRead ev = (GATKRead) matrix.getEvidence(r);
            final StringBuilder ll = new StringBuilder();
            for (int a = 0; a < nAl; a++) {
                if (a > 0) {
                    ll.append(",");
                }
                final Allele al = (Allele) matrix.getAllele(a);
                ll.append(alleleLabel(al))
                        .append("=")
                        .append(String.format(java.util.Locale.US, "%.12f", matrix.get(a, r)));
            }
            kv106(
                    "read_ll",
                    ev.getName()
                            + "\tflags="
                            + ev.getFlags()
                            + "\tstart="
                            + ev.getStart()
                            + "\t"
                            + ll.toString());
        }
    }

    private static String fnv1a64Hex(final byte[] data) {
        long h = 0xcbf29ce484222325L;
        for (int i = 0; i < data.length; i++) {
            h ^= (data[i] & 0xffL);
            h *= 0x100000001b3L;
        }
        return String.format("%016x", h);
    }

    private static void kv106(final String key, final String value) {
        System.out.println("6R106\t" + key + "\t" + value);
    }

    private static int minIndex(final int[] pl) {
        if (pl.length == 0) {
            return -1;
        }
        int best = 0;
        for (int i = 1; i < pl.length; i++) {
            if (pl[i] < pl[best]) {
                best = i;
            }
        }
        return best;
    }

    /** GATK diploid genotype index: for each j, all (i,j) with i <= j. */
    private static String diploidGtFromPlIndex(final int nAlleles, final int plIndex) {
        if (plIndex < 0) {
            return "./.";
        }
        int idx = 0;
        for (int j = 0; j < nAlleles; j++) {
            for (int i = 0; i <= j; i++) {
                if (idx == plIndex) {
                    return i + "/" + j;
                }
                idx++;
            }
        }
        return "?";
    }

    private static String ints(final int[] v) {
        final StringBuilder sb = new StringBuilder();
        for (int i = 0; i < v.length; i++) {
            if (i > 0) {
                sb.append(",");
            }
            sb.append(v[i]);
        }
        return sb.toString();
    }

    private static String doubles(final double[] v) {
        final StringBuilder sb = new StringBuilder();
        for (int i = 0; i < v.length; i++) {
            if (i > 0) {
                sb.append(",");
            }
            sb.append(String.format(java.util.Locale.US, "%.12f", v[i]));
        }
        return sb.toString();
    }

    private static void kv(final String key, final String value) {
        System.out.println("6R105\t" + key + "\t" + value);
    }
}
