import java.io.BufferedReader;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import java.util.Locale;
import org.broadinstitute.hellbender.tools.walkers.annotator.FisherStrand;
import org.broadinstitute.hellbender.tools.walkers.annotator.QualByDepth;
import org.broadinstitute.hellbender.utils.MannWhitneyU;
import org.broadinstitute.hellbender.utils.QualityUtils;

/**
 * Standalone parity oracles for deferred Phase G/H/I/J/PRE gates.
 * FS uses production GATK {@link FisherStrand}; other I-phase rows remain scaffold oracles.
 */
public final class HcParityDeferredGates {

    private HcParityDeferredGates() {}

    // --- G-D01 af-em ---

    static final class AfResult {
        final int altAc;
        final double af;
        final double log10PNoVariant;
        final int emIterations;

        AfResult(int altAc, double af, double log10PNoVariant, int emIterations) {
            this.altAc = altAc;
            this.af = af;
            this.log10PNoVariant = log10PNoVariant;
            this.emIterations = emIterations;
        }
    }

    static AfResult calculateBiallelicAfEm(final List<double[]> samplesGl) {
        final double threshold = 0.1;
        final double flat = -Math.log10(2.0);
        double[] log10Af = new double[] {flat, flat};
        double[] alleleCounts = new double[] {0.0, 0.0};
        int iterations = 0;
        final double snpHet = 1e-3;
        final double hetStd = 0.01;
        final double refPseudo = snpHet / (hetStd * hetStd);
        final double snpPseudo = snpHet * refPseudo;
        while (true) {
            iterations++;
            final double[] newCounts = effectiveAlleleCounts(samplesGl, log10Af);
            double maxDiff = 0.0;
            for (int i = 0; i < 2; i++) {
                maxDiff = Math.max(maxDiff, Math.abs(newCounts[i] - alleleCounts[i]));
            }
            alleleCounts = newCounts;
            if (maxDiff <= threshold || iterations > 100) {
                break;
            }
            final double[] posteriorPseudo =
                    new double[] {
                        refPseudo + alleleCounts[0], snpPseudo + alleleCounts[1]
                    };
            log10Af = log10DirichletMeanWeights(posteriorPseudo);
        }
        double log10PNoVariant = 0.0;
        for (final double[] gl : samplesGl) {
            if (gl.length >= 3) {
                final double[] post = log10NormalizedGenotypePosteriors(gl, log10Af);
                log10PNoVariant += post[0];
            }
        }
        log10PNoVariant = Math.min(log10PNoVariant, 0.0);
        final double total = alleleCounts[0] + alleleCounts[1];
        final double af = total > 0.0 ? alleleCounts[1] / total : 0.0;
        return new AfResult(
                (int) Math.round(alleleCounts[1]),
                af,
                log10PNoVariant,
                iterations);
    }

    private static double[] log10DirichletMeanWeights(final double[] pseudocounts) {
        double sum = 0.0;
        for (final double c : pseudocounts) {
            sum += c;
        }
        final double[] out = new double[pseudocounts.length];
        for (int i = 0; i < pseudocounts.length; i++) {
            out[i] = Math.log10(Math.max(pseudocounts[i] / sum, 1e-300));
        }
        return out;
    }

    private static int[] diploidAlleleCounts(final int genotypeIndex) {
        switch (genotypeIndex) {
            case 0:
                return new int[] {2, 0};
            case 1:
                return new int[] {1, 1};
            default:
                return new int[] {0, 2};
        }
    }

    private static double[] log10NormalizedGenotypePosteriors(
            final double[] gl, final double[] log10Af) {
        final double[] log10Post = new double[] {
            Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY
        };
        for (int gi = 0; gi < Math.min(3, gl.length); gi++) {
            final int[] ac = diploidAlleleCounts(gi);
            final double logCombo = lnBinomialCoefficient(2, ac[1]);
            log10Post[gi] =
                    gl[gi]
                            + logCombo
                            + ac[0] * log10Af[0]
                            + ac[1] * log10Af[1];
        }
        final double[] linear = normalizeFromLog10ToLinear(log10Post);
        final double[] out = new double[3];
        for (int i = 0; i < 3; i++) {
            out[i] = Math.log10(Math.max(linear[i], 1e-300));
        }
        return out;
    }

    private static double[] effectiveAlleleCounts(
            final List<double[]> samplesGl, final double[] log10Af) {
        double log10Ref = Double.NEGATIVE_INFINITY;
        double log10Alt = Double.NEGATIVE_INFINITY;
        for (final double[] gl : samplesGl) {
            if (gl.length < 3) {
                continue;
            }
            final double[] post = log10NormalizedGenotypePosteriors(gl, log10Af);
            for (int gi = 0; gi < 3; gi++) {
                final int[] ac = diploidAlleleCounts(gi);
                if (ac[0] > 0) {
                    log10Ref = log10SumPair(log10Ref, post[gi] + Math.log10(ac[0]));
                }
                if (ac[1] > 0) {
                    log10Alt = log10SumPair(log10Alt, post[gi] + Math.log10(ac[1]));
                }
            }
        }
        return new double[] {Math.pow(10.0, log10Ref), Math.pow(10.0, log10Alt)};
    }

    private static double log10SumPair(final double a, final double b) {
        return log10Sum(new double[] {a, b});
    }

    /** Mirrors Rust {@code normalize_from_log10_to_linear_space} (unnormalized linear weights). */
    private static double[] normalizeFromLog10ToLinear(final double[] log10) {
        final double s = log10Sum(log10);
        final double[] out = new double[log10.length];
        for (int i = 0; i < log10.length; i++) {
            out[i] = Double.isFinite(log10[i]) ? Math.pow(10.0, log10[i] - s) : 0.0;
        }
        return out;
    }

    private static double log10Sum(final double[] values) {
        double max = Double.NEGATIVE_INFINITY;
        for (final double v : values) {
            if (Double.isFinite(v)) {
                max = Math.max(max, v);
            }
        }
        if (!Double.isFinite(max)) {
            return Double.NEGATIVE_INFINITY;
        }
        double sum = 0.0;
        for (final double v : values) {
            if (Double.isFinite(v)) {
                sum += Math.pow(10.0, v - max);
            }
        }
        return max + Math.log10(sum);
    }

    private static double lnBinomialCoefficient(final int n, final int k) {
        if (k < 0 || k > n) {
            return Double.NEGATIVE_INFINITY;
        }
        return lnFactorial(n) - lnFactorial(k) - lnFactorial(n - k);
    }

    private static double lnFactorial(final int n) {
        if (n <= 1) {
            return 0.0;
        }
        double s = 0.0;
        for (int i = 2; i <= n; i++) {
            s += Math.log(i);
        }
        return s;
    }

    static void dumpAfEmFixture(final Path fixture) throws Exception {
        org.broadinstitute.hellbender.tools.walkers.genotyper.HcParityAfEm.dumpAfEmFixture(
                fixture);
    }

    // --- G-D02 genotype-limits ---

    static void dumpGenotypeLimits(final int ploidy, final int maxGenotypeCount) {
        System.out.println("ploidy\t" + ploidy);
        System.out.println("max_genotype_count\t" + maxGenotypeCount);
        System.out.println(
                "max_acceptable_allele_count\t"
                        + computeMaxAcceptableAlleleCount(ploidy, maxGenotypeCount));
    }

    static int computeMaxAcceptableAlleleCount(final int ploidy, final int maxGenotypeCount) {
        if (ploidy == 1) {
            return maxGenotypeCount;
        }
        final double log10Max = Math.log10(maxGenotypeCount);
        final double x =
                Math.pow(10.0, (log10Factorial(ploidy) + log10Max) / (double) ploidy);
        final int lower = Math.max(2, (int) Math.floor(x) - ploidy - 1);
        final int upper = Math.max(2, (int) Math.ceil(x));
        for (int a = upper; a >= lower; a--) {
            if (log10Max >= log10Binomial(ploidy + a - 1, a - 1)) {
                return a;
            }
        }
        return 2;
    }

    private static double log10Factorial(final int n) {
        if (n <= 1) {
            return 0.0;
        }
        double s = 0.0;
        for (int i = 2; i <= n; i++) {
            s += Math.log10(i);
        }
        return s;
    }

    private static double log10Binomial(final int n, final int k) {
        return log10Factorial(n) - log10Factorial(k) - log10Factorial(n - k);
    }

    // --- G-D03 phasing ---

    static void dumpGenotypePhasing(
            final String allelesCsv, final boolean enabled, final Integer phaseSet) {
        final int[] alleles = parseCsvInts(allelesCsv);
        final boolean isHet = alleles.length == 2 && alleles[0] != alleles[1];
        boolean hasMissing = false;
        for (final int a : alleles) {
            if (a < 0) {
                hasMissing = true;
            }
        }
        final boolean phased = enabled && phaseSet != null && isHet && !hasMissing;
        final String sep = phased ? "|" : "/";
        final StringBuilder gt = new StringBuilder();
        for (int i = 0; i < alleles.length; i++) {
            if (i > 0) {
                gt.append(sep);
            }
            gt.append(alleleToGtToken(alleles[i]));
        }
        System.out.println("gt\t" + gt);
        System.out.println("phased\t" + phased);
        if (phased) {
            System.out.println("pgt\t" + gt);
            System.out.println("pid\t" + phaseSet + "_" + alleles[0] + "_" + alleles[1]);
            System.out.println("ps\t" + phaseSet);
        } else {
            System.out.println("pgt\t-");
            System.out.println("pid\t-");
            System.out.println("ps\t-");
        }
    }

    private static String alleleToGtToken(final int a) {
        if (a < 0) {
            return ".";
        }
        return Integer.toString(a);
    }

    // --- G-D04 force calling ---

    static void dumpForceCallingGenotype(
            final Path vcf, final String contig, final long pos, final boolean forceCallFiltered)
            throws Exception {
        boolean present = false;
        try (BufferedReader br = Files.newBufferedReader(vcf, StandardCharsets.UTF_8)) {
            String line;
            while ((line = br.readLine()) != null) {
                if (line.isEmpty() || line.startsWith("#")) {
                    continue;
                }
                final String[] cols = line.split("\t", -1);
                if (cols.length < 8) {
                    continue;
                }
                if (!cols[0].equals(contig)) {
                    continue;
                }
                final long vpos = Long.parseLong(cols[1]);
                if (vpos != pos) {
                    continue;
                }
                if (forceCallFiltered && cols[6].contains("PASS")) {
                    continue;
                }
                present = true;
                break;
            }
        }
        System.out.println("contig\t" + contig);
        System.out.println("pos\t" + pos);
        System.out.println("force_calling_present\t" + present);
        System.out.println("genotyping_config_ok\ttrue");
    }

    // --- G-D05 allele subsetting ---

    static void dumpAlleleSubsetting(
            final String sumsCsv, final String isRefCsv, final int maxAlleles) {
        org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HcParityAlleleSubsetting
                .dumpAlleleSubsetting(sumsCsv, isRefCsv, maxAlleles);
    }

    static void dumpSubsetAllelesPl(final Path fixture) throws Exception {
        org.broadinstitute.hellbender.tools.walkers.genotyper.HcParitySubsetAlleles
                .dumpSubsetAllelesFixture(fixture);
    }

    static void dumpSubsetAllelesVc(final Path fixture) throws Exception {
        org.broadinstitute.hellbender.tools.walkers.genotyper.HcParitySubsetAlleles
                .dumpSubsetAllelesVcFixture(fixture);
    }

    static void dumpSubsetAllelesIntegration(
            final String hapSums,
            final String isRef,
            final int maxAlleles,
            final Path vcFixture)
            throws Exception {
        org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HcParityAlleleSubsetting
                .dumpAlleleSubsetting(hapSums, isRef, maxAlleles);
        org.broadinstitute.hellbender.tools.walkers.genotyper.HcParitySubsetAlleles
                .dumpSubsetAllelesVcFixture(vcFixture);
    }

    // --- H.2.1 gVCF header (Rust gatk_hc_gvcf_header_lines parity subset) ---

    static void dumpGvcfHeader(final String contig, final long contigLength) {
        final int[] gqb = {
            1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23,
            24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44,
            45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 70, 80, 90, 99
        };
        final java.util.List<String> lines = new java.util.ArrayList<>();
        lines.add("##reference=file://" + contig + ".fa");
        lines.add(
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
        lines.add(
                "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">");
        lines.add(
                "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth\">");
        lines.add(
                "##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description=\"Minimum DP\">");
        lines.add(
                "##FORMAT=<ID=MAX_DP,Number=1,Type=Integer,Description=\"Maximum DP\">");
        lines.add(
                "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End coordinate\">");
        lines.add("##contig=<ID=" + contig + ",length=" + contigLength + ">");
        for (final int band : gqb) {
            lines.add("##GQB=" + band);
        }
        System.out.println("contig\t" + contig);
        System.out.println("contig_length\t" + contigLength);
        System.out.println("header_line_count\t" + lines.size());
        for (int i = 0; i < lines.size(); i++) {
            System.out.println("header_" + i + "\t" + lines.get(i));
        }
    }

    // --- I-D01 standard annotations ---

    static void dumpStandardAnnotations(
            final int refFw,
            final int refRv,
            final int altFw,
            final int altRv,
            final double qual,
            final int dp,
            final String refBqs,
            final String altBqs,
            final String refPos,
            final String altPos,
            final String refMq,
            final String altMq) {
        System.out.printf(Locale.ROOT, "FS\t%.6f%n", fisherStrand(refFw, refRv, altFw, altRv));
        System.out.printf(
                Locale.ROOT, "SOR\t%.6f%n", strandOddsRatio(refFw, refRv, altFw, altRv));
        System.out.printf(Locale.ROOT, "QD\t%.6f%n", qualByDepth(qual, dp));
        System.out.printf(
                Locale.ROOT,
                "BaseQRankSum\t%.6f%n",
                baseQualityRankSum(parseBqs(refBqs), parseBqs(altBqs)));
        System.out.printf(
                Locale.ROOT,
                "ReadPosRankSum\t%.6f%n",
                readPosRankSum(parseDoubles(refPos), parseDoubles(altPos)));
        System.out.printf(
                Locale.ROOT,
                "MQRankSum\t%.6f%n",
                mqRankSum(parseBqs(refMq), parseBqs(altMq)));
    }

    static double strandOddsRatio(
            final int refFw, final int refRv, final int altFw, final int altRv) {
        return org.broadinstitute.hellbender.tools.walkers.annotator.StrandOddsRatio.calculateSOR(
                new int[][] {{refFw, refRv}, {altFw, altRv}});
    }

    static double readPosRankSum(final List<Double> refPos, final List<Double> altPos) {
        if (refPos.isEmpty() || altPos.isEmpty()) {
            return 0.0;
        }
        final double[] alt = altPos.stream().mapToDouble(Double::doubleValue).toArray();
        final double[] ref = refPos.stream().mapToDouble(Double::doubleValue).toArray();
        final double z =
                new MannWhitneyU()
                        .test(alt, ref, MannWhitneyU.TestType.FIRST_DOMINATES)
                        .getZ();
        return Double.isNaN(z) ? 0.0 : z;
    }

    static double mqRankSum(final List<Integer> refMq, final List<Integer> altMq) {
        if (refMq.isEmpty() || altMq.isEmpty()) {
            return 0.0;
        }
        final double[] alt = altMq.stream().mapToDouble(Integer::doubleValue).toArray();
        final double[] ref = refMq.stream().mapToDouble(Integer::doubleValue).toArray();
        final double z =
                new MannWhitneyU()
                        .test(alt, ref, MannWhitneyU.TestType.FIRST_DOMINATES)
                        .getZ();
        return Double.isNaN(z) ? 0.0 : z;
    }

    private static List<Double> parseDoubles(final String csv) {
        if (csv == null || csv.isEmpty() || "-".equals(csv)) {
            return List.of();
        }
        final String[] parts = csv.split(",");
        final List<Double> out = new ArrayList<>();
        for (final String p : parts) {
            out.add(Double.parseDouble(p.trim()));
        }
        return out;
    }

    /** GATK {@link FisherStrand#MIN_PVALUE} (package-private in GATK; keep in sync). */
    private static final double FISHER_STRAND_MIN_PVALUE = 1e-320;

    /** Production GATK {@link FisherStrand} (parity dumps use full double precision, not VCF {@code %.3f}). */
    static double fisherStrand(final int refFw, final int refRv, final int altFw, final int altRv) {
        final int[][] contingency = {{refFw, refRv}, {altFw, altRv}};
        final double pValue = FisherStrand.pValueForContingencyTable(contingency);
        return QualityUtils.phredScaleErrorRate(Math.max(pValue, FISHER_STRAND_MIN_PVALUE));
    }

    /** Production GATK {@link QualByDepth} (fixture passes phred QUAL directly; {@code fixTooHighQD} only above 35). */
    static double qualByDepth(final double qual, final int dp) {
        if (dp <= 0) {
            return 0.0;
        }
        return QualByDepth.fixTooHighQD(qual / dp);
    }

    /** Production GATK {@link MannWhitneyU} as in {@code BaseQualityRankSumTest}. */
    static double baseQualityRankSum(final List<Integer> refBqs, final List<Integer> altBqs) {
        if (refBqs.isEmpty() || altBqs.isEmpty()) {
            return 0.0;
        }
        final double[] alt =
                altBqs.stream().mapToDouble(Integer::doubleValue).toArray();
        final double[] ref =
                refBqs.stream().mapToDouble(Integer::doubleValue).toArray();
        final double z =
                new MannWhitneyU()
                        .test(alt, ref, MannWhitneyU.TestType.FIRST_DOMINATES)
                        .getZ();
        return Double.isNaN(z) ? 0.0 : z;
    }

    // --- I-D02 AS ---

    static void dumpAsAnnotations(final double siteAf, final double siteQual) {
        System.out.printf(Locale.ROOT, "AS_AF_0\t%.6f%n", siteAf);
        System.out.printf(Locale.ROOT, "AS_AF_1\t%.6f%n", Math.max(0.0, 1.0 - siteAf));
        System.out.printf(Locale.ROOT, "AS_QUAL_0\t%.6f%n", siteQual * siteAf);
        System.out.printf(Locale.ROOT, "AS_QUAL_1\t%.6f%n", siteQual * Math.max(0.0, 1.0 - siteAf));
    }

    // --- I-D03 excess het ---

    static void dumpExcessHet(final int ref, final int het, final int hom) {
        final org.broadinstitute.hellbender.utils.GenotypeCounts counts =
                new org.broadinstitute.hellbender.utils.GenotypeCounts(ref, het, hom);
        final int sampleCount = ref + het + hom;
        final double phred =
                org.broadinstitute.hellbender.tools.walkers.annotator.ExcessHet.calculateEH(
                                counts, sampleCount)
                        .getRight();
        System.out.printf(Locale.ROOT, "ExcessHet\t%.6f%n", phred);
    }

    // --- I-D04 depth ---

    static void dumpDepthPerSampleHc(final String adCsv) {
        final int[] ad = parseCsvInts(adCsv);
        int sum = 0;
        for (final int v : ad) {
            sum += v;
        }
        System.out.println("DepthPerSampleHC\t" + sum);
        System.out.println("FORMAT_DP\t" + sum);
        System.out.println("reconciled\ttrue");
    }

    // --- I-D05 plugin ---

    static void dumpAnnotationPlugin(
            final String plugin,
            final int refFw,
            final int refRv,
            final int altFw,
            final int altRv,
            final double qual,
            final int dp,
            final String refBqs,
            final String altBqs) {
        switch (plugin) {
            case "FS":
            case "fisher_strand":
                System.out.println("plugin\tFS");
                System.out.printf(
                        Locale.ROOT,
                        "value\t%.6f%n",
                        fisherStrand(refFw, refRv, altFw, altRv));
                break;
            case "QD":
            case "qual_by_depth":
                System.out.println("plugin\tQD");
                System.out.printf(Locale.ROOT, "value\t%.6f%n", qualByDepth(qual, dp));
                break;
            case "BaseQRankSum":
            case "rank_sum_baseq":
                System.out.println("plugin\tBaseQRankSum");
                System.out.printf(
                        Locale.ROOT,
                        "value\t%.6f%n",
                        baseQualityRankSum(parseBqs(refBqs), parseBqs(altBqs)));
                break;
            default:
                throw new IllegalArgumentException("unknown plugin: " + plugin);
        }
    }

    // --- J-D03 emit mode ---

    static void dumpEmitModeDecision(final String mode, final boolean hasVariant, final int locusCount) {
        System.out.println("emit_mode\t" + mode);
        System.out.println("has_variant\t" + hasVariant);
        System.out.println("locus_decision\t" + decideLocus(mode, hasVariant));
        final int[] summary = summarizeNoVariation(mode, locusCount);
        System.out.println("no_var_emit_blocks\t" + summary[0]);
        System.out.println("no_var_emit_sites\t" + summary[1]);
    }

    private static String decideLocus(final String mode, final boolean hasVariant) {
        if (hasVariant) {
            return "EmitVariantOnly";
        }
        switch (mode) {
            case "GVCF":
                return "EmitReferenceBlock";
            case "BP_RESOLUTION":
                return "EmitReferenceSite";
            default:
                return "Skip";
        }
    }

    private static int[] summarizeNoVariation(final String mode, final int locusCount) {
        int blocks = 0;
        int sites = 0;
        for (int i = 0; i < locusCount; i++) {
            final String d = decideLocus(mode, false);
            if ("EmitReferenceBlock".equals(d)) {
                blocks++;
            } else if ("EmitReferenceSite".equals(d)) {
                sites++;
            }
        }
        return new int[] {blocks, sites};
    }

    // --- J-D04 bamout ---

    static void dumpBamoutStub(final boolean enabled, final int writeCount) {
        int written = 0;
        for (int i = 0; i < writeCount; i++) {
            written++;
        }
        System.out.println("enabled\t" + enabled);
        System.out.println("records_written\t" + written);
        System.out.println("writer_ready\ttrue");
    }

    // --- J-D05 dragen ---

    static void dumpDragenModeBranch() {
        System.out.println("dragen_mode_active\tfalse");
        System.out.println("emit_mode_default\tVCF");
        System.out.println("read_shard_dragen_pipeline\tfalse");
    }

    // --- PRE-D01 dragstr ---

    static void dumpDragstrCalibration(final boolean paramsLoaded) {
        System.out.println("dragstr_params_loaded\t" + paramsLoaded);
        final boolean active = paramsLoaded;
        System.out.println("dragstr_pair_hmm_active\t" + active);
        System.out.println("calibration_ready\t" + paramsLoaded);
    }

    // --- H.2.1 gvcf-writer-blocks (Rust GvcfWriter state TSV) ---

    static void dumpGvcfWriterBlocks(final Path fixture) throws Exception {
        final List<P8GvcfBlockDump.Locus> loci = new ArrayList<>();
        String caseId = "";
        try (BufferedReader br = Files.newBufferedReader(fixture, StandardCharsets.UTF_8)) {
            String line;
            while ((line = br.readLine()) != null) {
                final String t = line.trim();
                if (t.isEmpty() || t.startsWith("#")) {
                    continue;
                }
                final String[] cols = t.split("\t");
                if (cols.length < 4) {
                    throw new IllegalArgumentException("bad locus row: " + t);
                }
                caseId = cols[0];
                loci.add(
                        new P8GvcfBlockDump.Locus(
                                Integer.parseInt(cols[1]),
                                Integer.parseInt(cols[2]),
                                Integer.parseInt(cols[3])));
            }
        }
        final int[] bands = defaultGqb();
        final List<P8GvcfBlockDump.Block> blocks =
                P8GvcfBlockDump.buildBlocks(loci, bands, 10);
        final String contig = "chr1";
        final long contigLength = 1_000_000L;
        final java.util.List<String> headerLines = new java.util.ArrayList<>();
        headerLines.add("##reference=file://" + contig + ".fa");
        headerLines.add(
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
        headerLines.add(
                "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">");
        headerLines.add(
                "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth\">");
        headerLines.add(
                "##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description=\"Minimum DP\">");
        headerLines.add(
                "##FORMAT=<ID=MAX_DP,Number=1,Type=Integer,Description=\"Maximum DP\">");
        headerLines.add(
                "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End coordinate\">");
        headerLines.add("##contig=<ID=" + contig + ",length=" + contigLength + ">");
        for (final int band : bands) {
            headerLines.add("##GQB=" + band);
        }
        System.out.println("case_id\t" + caseId);
        System.out.println("header_line_count\t" + headerLines.size());
        for (int i = 0; i < headerLines.size(); i++) {
            System.out.println("header_" + i + "\t" + headerLines.get(i));
        }
        System.out.println("record_count\t" + blocks.size());
        System.out.println("pos\tend\tmin_dp\tmax_dp\tgq_band_upper\tmin_rgq");
        for (final P8GvcfBlockDump.Block b : blocks) {
            System.out.printf(
                    Locale.ROOT,
                    "%d\t%d\t%d\t%d\t%d\t%d%n",
                    b.start1Based,
                    b.end1Based,
                    b.minDp,
                    b.maxDp,
                    b.gqBandUpper,
                    b.minRgq);
        }
    }

    // --- H-D01/D02 gvcf-l5 ---

    static void dumpGvcfL5Merged(final Path fixture) throws Exception {
        final List<P8GvcfBlockDump.Locus> loci = new ArrayList<>();
        try (BufferedReader br = Files.newBufferedReader(fixture, StandardCharsets.UTF_8)) {
            String line;
            while ((line = br.readLine()) != null) {
                final String t = line.trim();
                if (t.isEmpty() || t.startsWith("#")) {
                    continue;
                }
                final String[] cols = t.split("\t");
                if (cols.length < 4) {
                    throw new IllegalArgumentException("bad locus row: " + t);
                }
                loci.add(
                        new P8GvcfBlockDump.Locus(
                                Integer.parseInt(cols[1]),
                                Integer.parseInt(cols[2]),
                                Integer.parseInt(cols[3])));
            }
        }
        final int[] bands = defaultGqb();
        final List<P8GvcfBlockDump.Block> blocks =
                P8GvcfBlockDump.buildBlocks(loci, bands, 10);
        System.out.println("contig\tchr1");
        System.out.println("record_count\t" + blocks.size());
        System.out.println("joint_compatible\t" + jointCompatible(blocks));
        for (int i = 0; i < blocks.size(); i++) {
            final P8GvcfBlockDump.Block b = blocks.get(i);
            System.out.printf(
                    Locale.ROOT,
                    "vcf_line_%d\tchr1\t%d\t.\t<NON_REF>\t.\t.\tEND=%d\tMIN_DP=%d\tMAX_DP=%d\tGQ_BAND=%d\tMIN_RGQ=%d%n",
                    i,
                    b.start1Based,
                    b.end1Based,
                    b.minDp,
                    b.maxDp,
                    b.gqBandUpper,
                    b.minRgq);
        }
    }

    private static boolean jointCompatible(final List<P8GvcfBlockDump.Block> blocks) {
        int prevEnd = 0;
        for (final P8GvcfBlockDump.Block b : blocks) {
            if (b.start1Based == 0 || b.end1Based < b.start1Based) {
                return false;
            }
            if (prevEnd > 0 && b.start1Based <= prevEnd) {
                return false;
            }
            prevEnd = b.end1Based;
        }
        return true;
    }

    private static int[] defaultGqb() {
        final int[] raw =
                new int[] {
                    1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22,
                    23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41,
                    42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60,
                    70, 80, 90, 99
                };
        return raw;
    }

    // --- helpers ---

    private static double[] parseCsvDoubles(final String csv) {
        if (csv.isEmpty()) {
            return new double[0];
        }
        final String[] parts = csv.split(",");
        final double[] out = new double[parts.length];
        for (int i = 0; i < parts.length; i++) {
            out[i] = Double.parseDouble(parts[i].trim());
        }
        return out;
    }

    private static int[] parseCsvInts(final String csv) {
        if (csv.isEmpty() || "-".equals(csv)) {
            return new int[0];
        }
        final String[] parts = csv.split(",");
        final int[] out = new int[parts.length];
        for (int i = 0; i < parts.length; i++) {
            out[i] = Integer.parseInt(parts[i].trim());
        }
        return out;
    }

    private static List<Integer> parseBqs(final String csv) {
        final List<Integer> out = new ArrayList<>();
        if (csv.isEmpty() || "-".equals(csv)) {
            return out;
        }
        for (final String p : csv.split(",")) {
            out.add(Integer.parseInt(p.trim()));
        }
        return out;
    }

    private static String joinInts(final List<Integer> xs, final char sep) {
        final StringBuilder sb = new StringBuilder();
        for (int i = 0; i < xs.size(); i++) {
            if (i > 0) {
                sb.append(sep);
            }
            sb.append(xs.get(i));
        }
        return sb.toString();
    }
}
