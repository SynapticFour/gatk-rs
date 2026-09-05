import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import java.lang.reflect.Field;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.tools.walkers.annotator.DepthPerAlleleBySample;
import org.broadinstitute.hellbender.tools.walkers.annotator.GenotypeAnnotation;
import org.broadinstitute.hellbender.tools.walkers.annotator.VariantAnnotatorEngine;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCallerEngine;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCallerGenotypingEngine;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.genotyper.LikelihoodMatrix;
import org.broadinstitute.hellbender.utils.read.GATKRead;

/**
 * TEST-ONLY live dump of {@link DepthPerAlleleBySample#annotate} inputs/outputs.
 * Installed as the sole genotype annotation on a live {@link HaplotypeCallerEngine}.
 * Does not change annotation arithmetic: delegates to production {@link DepthPerAlleleBySample}.
 */
public final class HcParityAdAnnotationDump implements GenotypeAnnotation {

    private static final DepthPerAlleleBySample DELEGATE = new DepthPerAlleleBySample();
    private static int callSeq = 0;

    public static void installOn(final HaplotypeCallerEngine engine) throws Exception {
        final VariantAnnotatorEngine dumping =
                new VariantAnnotatorEngine(
                        Collections.singletonList(new HcParityAdAnnotationDump()),
                        null,
                        Collections.emptyList(),
                        false,
                        false);
        final Field ae = HaplotypeCallerEngine.class.getDeclaredField("annotationEngine");
        ae.setAccessible(true);
        ae.set(engine, dumping);
        final Field ge = HaplotypeCallerEngine.class.getDeclaredField("genotypingEngine");
        ge.setAccessible(true);
        final HaplotypeCallerGenotypingEngine genotyping =
                (HaplotypeCallerGenotypingEngine) ge.get(engine);
        genotyping.setAnnotationEngine(dumping);
    }

    @Override
    public List<String> getKeyNames() {
        return DELEGATE.getKeyNames();
    }

    @Override
    public void annotate(
            final ReferenceContext ref,
            final VariantContext vc,
            final Genotype g,
            final GenotypeBuilder gb,
            final AlleleLikelihoods<GATKRead, Allele> likelihoods) {
        final int seq = ++callSeq;
        dumpBefore(seq, vc, g, likelihoods);
        DELEGATE.annotate(ref, vc, g, gb, likelihoods);
        dumpAfter(seq, vc, g, gb);
    }

    private static void dumpBefore(
            final int seq,
            final VariantContext vc,
            final Genotype g,
            final AlleleLikelihoods<GATKRead, Allele> likelihoods) {
        kv(seq, "call_class", vc == null ? "null" : vc.getClass().getName());
        if (vc != null) {
            kv(seq, "call_contig", vc.getContig());
            kv(seq, "call_start", Integer.toString(vc.getStart()));
            kv(seq, "call_end", Integer.toString(vc.getEnd()));
            kv(seq, "call_id", vc.getID());
            kv(seq, "call_n_alleles", Integer.toString(vc.getNAlleles()));
            kv(seq, "call_identity", Integer.toString(System.identityHashCode(vc)));
            dumpAlleles(seq, "call", vc.getAlleles());
        }
        if (g != null) {
            kv(seq, "genotype_sample", g.getSampleName());
            kv(seq, "genotype_called", Boolean.toString(g.isCalled()));
            kv(seq, "genotype_has_ad_before", Boolean.toString(g.hasAD()));
            kv(seq, "genotype_ad_before", g.hasAD() ? ints(g.getAD()) : "ABSENT");
            kv(seq, "genotype_alleles", alleleList(g.getAlleles()));
        } else {
            kv(seq, "genotype", "null");
        }
        if (likelihoods == null) {
            kv(seq, "likelihoods", "null");
            return;
        }
        kv(seq, "ll_class", likelihoods.getClass().getName());
        kv(seq, "ll_identity", Integer.toString(System.identityHashCode(likelihoods)));
        kv(seq, "ll_n_alleles", Integer.toString(likelihoods.numberOfAlleles()));
        kv(seq, "ll_evidence_count", Integer.toString(likelihoods.evidenceCount()));
        kv(seq, "ll_n_samples", Integer.toString(likelihoods.numberOfSamples()));
        final SimpleInterval subset = likelihoods.getVariantCallingSubsetApplied();
        kv(seq, "ll_variant_calling_subset", subset == null ? "null" : subset.toString());
        dumpAlleles(seq, "ll", likelihoods.alleles());
        compareCallVsLikelihoodsAlleles(seq, vc, likelihoods);
        for (int s = 0; s < likelihoods.numberOfSamples(); s++) {
            kv(seq, "sample_index_" + s, likelihoods.getSample(s));
            kv(
                    seq,
                    "sample_evidence_count_" + s,
                    Integer.toString(likelihoods.sampleEvidenceCount(s)));
            kv(
                    seq,
                    "sample_filtered_evidence_count_" + s,
                    Integer.toString(likelihoods.filteredSampleEvidence(s).size()));
        }
        dumpEvidenceAndMatrix(seq, likelihoods);
        dumpFilteredEvidence(seq, likelihoods);
        if (vc != null && g != null && g.isCalled()) {
            dumpIndependentReconstruction(seq, vc, g, likelihoods);
        }
    }

    private static void dumpAfter(
            final int seq, final VariantContext vc, final Genotype g, final GenotypeBuilder gb) {
        final Genotype after = gb.make();
        kv(seq, "ad_after_has", Boolean.toString(after.hasAD()));
        kv(seq, "ad_after", after.hasAD() ? ints(after.getAD()) : "ABSENT");
        kv(seq, "ad_first_write", Boolean.toString(g != null && !g.hasAD() && after.hasAD()));
        if (vc != null) {
            kv(
                    seq,
                    "ad_site",
                    vc.getContig()
                            + ":"
                            + vc.getStart()
                            + "\t"
                            + alleleList(vc.getAlleles())
                            + "\t"
                            + (after.hasAD() ? ints(after.getAD()) : "ABSENT"));
        }
    }

    private static void dumpAlleles(final int seq, final String tag, final List<Allele> alleles) {
        kv(seq, tag + "_allele_count", Integer.toString(alleles.size()));
        kv(seq, tag + "_allele_display", alleleList(alleles));
        for (int i = 0; i < alleles.size(); i++) {
            final Allele a = alleles.get(i);
            kv(
                    seq,
                    tag + "_allele_" + i,
                    i
                            + "\t"
                            + a.getBaseString()
                            + "\t"
                            + a.getDisplayString()
                            + "\t"
                            + "isRef="
                            + a.isReference()
                            + "\t"
                            + "isNonRef="
                            + a.isNonReference()
                            + "\t"
                            + "isSymbolic="
                            + a.isSymbolic()
                            + "\t"
                            + "len="
                            + a.length()
                            + "\t"
                            + "id="
                            + System.identityHashCode(a));
        }
    }

    private static void compareCallVsLikelihoodsAlleles(
            final int seq, final VariantContext vc, final AlleleLikelihoods<GATKRead, Allele> ll) {
        if (vc == null) {
            return;
        }
        final List<Allele> call = vc.getAlleles();
        final List<Allele> llAll = ll.alleles();
        kv(seq, "alleles_equal_lists", Boolean.toString(call.equals(llAll)));
        kv(seq, "alleles_same_order_identity", Boolean.toString(sameOrderIdentity(call, llAll)));
        kv(seq, "ll_containsAll_call", Boolean.toString(llAll.containsAll(call)));
        kv(seq, "call_containsAll_ll", Boolean.toString(call.containsAll(llAll)));
        final int n = Math.max(call.size(), llAll.size());
        for (int i = 0; i < n; i++) {
            final Allele c = i < call.size() ? call.get(i) : null;
            final Allele l = i < llAll.size() ? llAll.get(i) : null;
            kv(
                    seq,
                    "allele_side_by_side_" + i,
                    "idx="
                            + i
                            + "\tcall="
                            + (c == null ? "MISSING" : alleleId(c))
                            + "\tll="
                            + (l == null ? "MISSING" : alleleId(l))
                            + "\tequals="
                            + (c != null && l != null && c.equals(l))
                            + "\tbases_eq="
                            + (c != null
                                    && l != null
                                    && c.getBaseString().equals(l.getBaseString()))
                            + "\tref_eq="
                            + (c != null && l != null && c.isReference() == l.isReference())
                            + "\tsym_eq="
                            + (c != null && l != null && c.isSymbolic() == l.isSymbolic())
                            + "\tsame_object="
                            + (c != null && l != null && c == l));
        }
        for (int i = 0; i < call.size(); i++) {
            final Allele c = call.get(i);
            int found = -1;
            boolean sameObj = false;
            for (int j = 0; j < llAll.size(); j++) {
                if (c.equals(llAll.get(j))) {
                    found = j;
                    sameObj = c == llAll.get(j);
                    break;
                }
            }
            kv(
                    seq,
                    "call_in_ll_" + i,
                    alleleId(c)
                            + "\tll_index="
                            + found
                            + "\tsame_object="
                            + sameObj);
        }
    }

    private static boolean sameOrderIdentity(final List<Allele> a, final List<Allele> b) {
        if (a.size() != b.size()) {
            return false;
        }
        for (int i = 0; i < a.size(); i++) {
            if (a.get(i) != b.get(i)) {
                return false;
            }
        }
        return true;
    }

    private static void dumpEvidenceAndMatrix(
            final int seq, final AlleleLikelihoods<GATKRead, Allele> ll) {
        for (int s = 0; s < ll.numberOfSamples(); s++) {
            final String sample = ll.getSample(s);
            final List<GATKRead> ev = ll.sampleEvidence(s);
            final LikelihoodMatrix<GATKRead, Allele> mx = ll.sampleMatrix(s);
            for (int r = 0; r < ev.size(); r++) {
                final GATKRead read = ev.get(r);
                final StringBuilder llCols = new StringBuilder();
                for (int a = 0; a < mx.numberOfAlleles(); a++) {
                    if (a > 0) {
                        llCols.append(',');
                    }
                    llCols.append(Double.toString(mx.get(a, r)));
                }
                kv(
                        seq,
                        "evidence_" + s + "_" + r,
                        "sample="
                                + sample
                                + "\tqname="
                                + read.getName()
                                + "\tcontig="
                                + read.getContig()
                                + "\tstart="
                                + read.getStart()
                                + "\tend="
                                + read.getEnd()
                                + "\tflags="
                                + read.getFlags()
                                + "\tcigar="
                                + (read.getCigar() == null ? "." : read.getCigar().toString())
                                + "\trg="
                                + nullToDot(read.getReadGroup())
                                + "\tfp="
                                + fingerprint(read)
                                + "\tll="
                                + llCols);
            }
        }
    }

    private static void dumpIndependentReconstruction(
            final int seq,
            final VariantContext vc,
            final Genotype g,
            final AlleleLikelihoods<GATKRead, Allele> likelihoods) {
        final Set<Allele> remaining = new LinkedHashSet<>(vc.getAlleles());
        kv(seq, "remaining_alleles", alleleList(new ArrayList<>(remaining)));
        final boolean containsAll = likelihoods.alleles().containsAll(remaining);
        kv(seq, "remaining_subset_of_ll", Boolean.toString(containsAll));
        if (!containsAll) {
            kv(seq, "independent_ad", "SKIP_NOT_SUBSET");
            return;
        }
        final Map<Allele, List<Allele>> alleleSubset =
                remaining.stream().collect(Collectors.toMap(a -> a, Arrays::asList));
        for (final Map.Entry<Allele, List<Allele>> e : alleleSubset.entrySet()) {
            kv(
                    seq,
                    "mapping_" + e.getKey().getDisplayString(),
                    alleleList(e.getValue()));
        }
        final AlleleLikelihoods<GATKRead, Allele> subsetted = likelihoods.marginalize(alleleSubset);
        kv(seq, "subsetted_identity", Integer.toString(System.identityHashCode(subsetted)));
        kv(seq, "subsetted_n_alleles", Integer.toString(subsetted.numberOfAlleles()));
        kv(seq, "subsetted_evidence_count", Integer.toString(subsetted.evidenceCount()));
        kv(seq, "subsetted_alleles", alleleList(subsetted.alleles()));
        dumpAlleles(seq, "subsetted", subsetted.alleles());

        final Map<Allele, Integer> alleleCounts = new LinkedHashMap<>();
        for (final Allele allele : vc.getAlleles()) {
            alleleCounts.put(allele, 0);
        }
        int nInformative = 0;
        int nUninformative = 0;
        int evIdx = 0;
        for (final AlleleLikelihoods<GATKRead, Allele>.BestAllele ba :
                subsetted.bestAllelesBreakingTies(g.getSampleName())) {
            final GATKRead read = (GATKRead) ba.evidence;
            kv(
                    seq,
                    "best_" + evIdx,
                    "fp="
                            + fingerprint(read)
                            + "\tqname="
                            + read.getName()
                            + "\tsample="
                            + ba.sample
                            + "\tbest="
                            + (ba.allele == null ? "null" : alleleId(ba.allele))
                            + "\tsecond="
                            + (ba.second_best_allele == null
                                    ? "null"
                                    : alleleId(ba.second_best_allele))
                            + "\tbest_ll="
                            + ba.likelihood
                            + "\tsecond_ll="
                            + ba.secondBestLikelihood
                            + "\tconf="
                            + ba.confidence
                            + "\tinformative="
                            + ba.isInformative());
            if (ba.isInformative()) {
                nInformative++;
                alleleCounts.compute(ba.allele, (allele, prev) -> prev + 1);
            } else {
                nUninformative++;
            }
            evIdx++;
        }
        final int[] counts = new int[alleleCounts.size()];
        counts[0] = alleleCounts.get(vc.getReference());
        for (int i = 0; i < vc.getNAlleles() - 1; i++) {
            counts[i + 1] = alleleCounts.get(vc.getAlternateAllele(i));
        }
        kv(seq, "independent_ad", ints(counts));
        kv(seq, "independent_n_informative", Integer.toString(nInformative));
        kv(seq, "independent_n_uninformative", Integer.toString(nUninformative));
    }

    private static void dumpFilteredEvidence(
            final int seq, final AlleleLikelihoods<GATKRead, Allele> ll) {
        for (int s = 0; s < ll.numberOfSamples(); s++) {
            final String sample = ll.getSample(s);
            final List<GATKRead> filtered = ll.filteredSampleEvidence(s);
            kv(seq, "filtered_n_" + s, Integer.toString(filtered.size()));
            for (int r = 0; r < filtered.size(); r++) {
                final GATKRead read = filtered.get(r);
                kv(
                        seq,
                        "filtered_" + s + "_" + r,
                        "sample="
                                + sample
                                + "\tqname="
                                + read.getName()
                                + "\tcontig="
                                + read.getContig()
                                + "\tstart="
                                + read.getStart()
                                + "\tend="
                                + read.getEnd()
                                + "\tflags="
                                + read.getFlags()
                                + "\tcigar="
                                + (read.getCigar() == null ? "." : read.getCigar().toString())
                                + "\trg="
                                + nullToDot(read.getReadGroup())
                                + "\tfp="
                                + fingerprint(read));
            }
        }
    }

    static String fingerprint(final GATKRead read) {
        final String bases = new String(read.getBases(), StandardCharsets.US_ASCII);
        final byte[] q = read.getBaseQualities();
        final StringBuilder bq = new StringBuilder(q.length);
        for (final byte b : q) {
            bq.append((char) (Byte.toUnsignedInt(b) + 33));
        }
        return bases
                + "|"
                + bq
                + "|"
                + read.getFlags()
                + "|"
                + (read.getCigar() == null ? "." : read.getCigar().toString());
    }

    static String alleleId(final Allele a) {
        return a.getDisplayString()
                + "/isRef="
                + a.isReference()
                + "/sym="
                + a.isSymbolic()
                + "/len="
                + a.length();
    }

    static String alleleList(final List<Allele> alleles) {
        return alleles.stream().map(Allele::getDisplayString).collect(Collectors.joining(","));
    }

    static String ints(final int[] v) {
        final StringBuilder sb = new StringBuilder();
        for (int i = 0; i < v.length; i++) {
            if (i > 0) {
                sb.append(',');
            }
            sb.append(v[i]);
        }
        return sb.toString();
    }

    static String nullToDot(final String s) {
        return s == null || s.isEmpty() ? "." : s;
    }

    static void kv(final int seq, final String key, final String value) {
        System.out.println("6R91\t" + seq + "\t" + key + "\t" + value);
    }
}
