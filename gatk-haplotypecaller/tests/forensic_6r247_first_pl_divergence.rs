//! 6R.247: proof-only. First causal PL split at `20:29455649 T/TGTTTG`.
//!
//! Java PL `570,0,3517` vs Rust PL `570,0,3518`. Do not assume rounding.
//! PRODUCTION CHANGE: NONE. Do not retune PL, GLs, PairHMM, QUAL, AD, GQ,
//! priors, normalization, epsilon, formatting, or integer rounding constants.
//!
//! Frozen Java 4.4.0.0 SHA `2dbc025821bc5f686c423ff332a41e6cef892a77`.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r247_first_pl_divergence -- --nocapture --test-threads=1
//! HOLDOUT_6R247=1 cargo test -p gatk-haplotypecaller --test holdout_6r247_first_pl_divergence -- --nocapture --test-threads=1
//! ```

use gatk_core::io::vcf::InfoValue;
use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::genotyping::{
    emit_genotype_format_fields, gq_phred_from_pl, gq_phred_p7_compatible,
};
use gatk_haplotypecaller::hc_genotyping_engine::take_colocated_merge_numerics;
use gatk_haplotypecaller::read_threading_assembler::DEFAULT_NUM_BEST_HAPLOTYPES_PER_GRAPH;
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, traverse_assembly_region_walker,
    try_emit_call_region_variants, AssemblyRegionCallDisposition, CallRegionArgs,
    HaplotypeCallerEngine, ReadFilterParams, WalkerTraversalConfig, DEFAULT_STAND_EMIT_CONFIDENCE,
};
use std::fs;
use std::path::{Path, PathBuf};

const JAVA_PIN: &str = "4.4.0.0 / 2dbc025821bc5f686c423ff332a41e6cef892a77";
const INTERVAL: &str = "20:29455000-29456500";
const BAM_REL: &str = "parity/giab/runs/local-pairhmm-diff/HG001.20-29455000-29456500.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const JAVA_VCF_REL: &str = "parity/reports/6r43/chr20_tiny/java.vcf";
const TARGET: u64 = 29_455_649;
const TARGET_REF: &str = "T";
const TARGET_ALT: &str = "TGTTTG";
const JAVA_PL: [i32; 3] = [570, 0, 3517];
const RUST_PL: [i32; 3] = [570, 0, 3518];

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    eprintln!("6R247\t{key}\t{}", value.as_ref());
}

fn fmt_f64(x: f64) -> String {
    format!("{x:.17} bits=0x{:016x}", x.to_bits())
}

fn fmt_vec_f64(xs: &[f64]) -> String {
    xs.iter()
        .map(|x| fmt_f64(*x))
        .collect::<Vec<_>>()
        .join(" | ")
}

fn fmt_i32(xs: &[i32]) -> String {
    xs.iter()
        .map(|x| x.to_string())
        .collect::<Vec<_>>()
        .join(",")
}

/// Java `Math.round(double)` for finite values: `Math.floor(x + 0.5)`.
fn java_math_round(x: f64) -> i32 {
    (x + 0.5).floor() as i32
}

/// HTSJDK `GenotypeLikelihoods.GLsToPLs` without the Integer.MAX_VALUE cap.
fn java_gls_to_pls(gls: &[f64]) -> Vec<i32> {
    if gls.is_empty() {
        return Vec::new();
    }
    let adjust = gls.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    gls.iter()
        .map(|&g| java_math_round((-10.0 * (g - adjust)).min(i32::MAX as f64)).max(0))
        .collect()
}

/// Production `emit_genotype_format_fields` PL conversion (IEEE `f64::round`).
fn rust_gls_to_pls(gls: &[f64]) -> Vec<i32> {
    if gls.is_empty() {
        return Vec::new();
    }
    let max_ll = gls.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let mut pls: Vec<i32> = gls
        .iter()
        .map(|ll| ((-10.0 * (ll - max_ll)).round() as i32).max(0))
        .collect();
    let min_pl = pls.iter().copied().min().unwrap_or(0);
    for pl in &mut pls {
        *pl -= min_pl;
    }
    pls
}

fn continuous_pl(gls: &[f64]) -> Vec<f64> {
    let max_ll = gls.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    gls.iter().map(|ll| -10.0 * (ll - max_ll)).collect()
}

/// HTSJDK `PLsToGLs`: `pl[i] / -10.0`.
fn java_pls_to_gls(pls: &[i32]) -> Vec<f64> {
    pls.iter().map(|&p| (p as f64) / -10.0).collect()
}

/// GATK `AlleleSubsettingUtils.subsettedPLIndices` for diploid keep-allele lists.
fn subsetted_diploid_pl_indices(keep: &[usize]) -> Vec<usize> {
    let mut out = Vec::new();
    for (nj, &aj) in keep.iter().enumerate() {
        for ni in 0..=nj {
            let ai = keep[ni];
            let (i, j) = if ai <= aj { (ai, aj) } else { (aj, ai) };
            out.push(j * (j + 1) / 2 + i);
        }
    }
    out
}

fn scale_log_space(gl: &[f64]) -> Vec<f64> {
    if gl.is_empty() {
        return Vec::new();
    }
    let max_v = gl.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    gl.iter().map(|g| g - max_v).collect()
}

fn info_i32(info: &[InfoValue], key: &str) -> Option<i32> {
    for v in info {
        if let InfoValue::Integer(k, xs) = v {
            if k == key {
                return xs.first().copied();
            }
        }
    }
    None
}

fn src_contains(rel: &str, needle: &str) -> bool {
    fs::read_to_string(repo_root().join(rel))
        .unwrap_or_default()
        .contains(needle)
}

#[test]
fn forensic_6r247_source_pl_paths_are_documented() {
    kv("java_pin", JAVA_PIN);
    kv("production_change", "NONE");
    assert_eq!(DEFAULT_NUM_BEST_HAPLOTYPES_PER_GRAPH, 128);

    let genotyping = fs::read_to_string(repo_root().join("gatk-haplotypecaller/src/genotyping.rs"))
        .expect("genotyping.rs");
    assert!(
        genotyping.contains("((-10.0 * (ll - max_ll)).round() as i32).max(0)"),
        "Rust FORMAT PL is f64::round of -10*(GL-max)"
    );
    assert!(
        genotyping.contains(
            "let gq = gq_phred_p7_compatible(genotype_log10_likelihoods, &pls, allele_depths);"
        ),
        "GQ is a separate helper, not a copy of PL[2]"
    );

    let assign = fs::read_to_string(
        repo_root().join("gatk-haplotypecaller/src/hc_genotyping_engine/genotype_assign.rs"),
    )
    .expect("genotype_assign.rs");
    assert!(
        assign.contains("let format = emit_genotype_format_fields(&unused.log10_gls, &remarg_ad)?"),
        "colocated merge FORMAT PL is emit_genotype_format_fields(unused.log10_gls)"
    );
    assert!(
        assign.contains("subset_unused_alts_after_merged_genotyping"),
        "unused-ALT subset of calculator GLs precedes FORMAT PL"
    );

    let subset =
        fs::read_to_string(repo_root().join("gatk-haplotypecaller/src/allele_subsetting_pl.rs"))
            .expect("allele_subsetting_pl.rs");
    assert!(
        subset.contains("scale_log_space_array_for_numerical_stability"),
        "unused-ALT GLs are max-shifted, not integer-quantized"
    );

    assert!(
        !src_contains(
            "gatk-haplotypecaller/src/genotyping.rs",
            "java_math_round_double_to_i32"
        ),
        "production FORMAT PL must not silently switch to Java Math.round in 6R.247"
    );
    kv("rust_pl_op", "f64::round(-10*(GL-max)) then subtract min");
    kv(
        "java_pl_op",
        "GenotypeLikelihoods.GLsToPLs: (int) Math.round(min(-10*(GL-adjust), MAX_PL))",
    );
}

#[test]
fn forensic_6r247_first_pl_divergence_at_canonical_site() {
    kv("java_pin", JAVA_PIN);
    kv("production_change", "NONE");
    kv("target", "20:29455649 T/TGTTTG");

    let root = repo_root();
    if !root.join(JAVA_VCF_REL).is_file()
        || !root.join(REF_REL).is_file()
        || !root.join(BAM_REL).is_file()
    {
        eprintln!("skip: missing chr20_tiny BAM/REF");
        return;
    }
    let java_vcf = fs::read_to_string(root.join(JAVA_VCF_REL)).expect("java.vcf");
    let java_line = java_vcf
        .lines()
        .find(|l| l.starts_with("20\t29455649\t") && l.contains("\tT\tTGTTTG\t"))
        .expect("java oracle line");
    assert!(
        java_line.contains("570,0,3517"),
        "frozen Java PL; got {java_line}"
    );
    kv("java_vcf_pl", "570,0,3517");

    let ref_fasta = root.join(REF_REL);
    let bam = root.join(BAM_REL);
    if !ref_fasta.is_file() || !bam.is_file() {
        eprintln!("skip: missing chr20_tiny BAM/REF");
        return;
    }
    let dict = SequenceDictionary::from_fasta_path(&ref_fasta).expect("dict");
    let specs = parse_intervals_cli_string(&dict, INTERVAL).expect("interval");
    let walk = traverse_assembly_region_walker(
        &dict,
        &specs,
        &ref_fasta,
        &bam,
        &ReadFilterParams::gatk_standard_hc(),
        &WalkerTraversalConfig::gatk_haplotype_caller_production(100),
    )
    .expect("walk");
    let regions = flatten_assembly_regions(&walk);
    let region = regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= TARGET
                && r.end.get() >= TARGET
        })
        .expect("ActiveFull");
    let outcome = HaplotypeCallerEngine::call_region(
        region,
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("call")
    .expect("Some");
    let call = outcome
        .genotyped_calls
        .iter()
        .find(|c| {
            c.event.start_1based.get() == TARGET
                && c.event.ref_allele == TARGET_REF
                && c.event.alt_allele == TARGET_ALT
        })
        .expect("T/TGTTTG");

    let snap = take_colocated_merge_numerics()
        .into_iter()
        .find(|s| s.loc == TARGET)
        .expect("merge snapshot");

    let subset_gls = &call.genotype.genotype_log10_likelihoods;
    let merged_gls = &snap.merged_gls;
    kv("merged_alts", snap.alts.join(","));
    kv("assigned_gt", format!("{:?}", snap.assigned_gt));
    kv("keep_indices", format!("{:?}", snap.remaining_keep_indices));
    kv("merged_gls", fmt_vec_f64(merged_gls));
    kv("subset_gls_raw", fmt_vec_f64(subset_gls));

    assert_eq!(subset_gls.len(), 3, "biallelic unused-ALT GLs");
    assert_eq!(merged_gls.len(), 6, "merged three-allele diploid GLs");
    assert!(
        snap.alts.iter().any(|a| a == TARGET_ALT),
        "merged alleles include TGTTTG; got {:?}",
        snap.alts
    );
    assert_eq!(
        snap.remaining_keep_indices.first().copied(),
        Some(0),
        "keep always starts at REF"
    );

    let rust_emit =
        emit_genotype_format_fields(subset_gls, &call.genotype.format.ad_as_i32()).expect("emit");
    let rust_pl = rust_emit.pl_as_i32();
    kv("rust_emit_pl", fmt_i32(&rust_pl));
    assert_eq!(
        rust_pl, RUST_PL,
        "canonical Rust FORMAT PL must remain 570,0,3518"
    );
    assert_eq!(
        call.genotype.format.pl_as_i32(),
        RUST_PL,
        "Call FORMAT PL is emit_genotype_format_fields(unused.log10_gls)"
    );

    let cont_subset = continuous_pl(subset_gls);
    let cont_merged = continuous_pl(merged_gls);
    kv("subset_continuous_pl", fmt_vec_f64(&cont_subset));
    kv("merged_continuous_pl", fmt_vec_f64(&cont_merged));

    let rust_round_subset = rust_gls_to_pls(subset_gls);
    let java_round_subset = java_gls_to_pls(subset_gls);
    let rust_round_merged = rust_gls_to_pls(merged_gls);
    let java_round_merged = java_gls_to_pls(merged_gls);
    kv("rust_f64_round_subset_pl", fmt_i32(&rust_round_subset));
    kv("java_math_round_subset_pl", fmt_i32(&java_round_subset));
    kv("rust_f64_round_merged_pl", fmt_i32(&rust_round_merged));
    kv("java_math_round_merged_pl", fmt_i32(&java_round_merged));
    kv("snapshot_merged_pl", fmt_i32(&snap.merged_pl));
    kv("snapshot_subset_pl", fmt_i32(&snap.subset_pl));

    // Java calculateGLsForThisEvent stores GLsToPLs(6 GLs), then subsetAlleles
    // reads g.getLikelihoods().getAsVector() = PLsToGLs, picks keep indices,
    // scaleLogSpace, then GenotypeBuilder.PL(double[]) → GLsToPLs again.
    let java_six_then_subset = {
        let six_pl = java_gls_to_pls(merged_gls);
        let pick = subsetted_diploid_pl_indices(&snap.remaining_keep_indices);
        kv("subsetted_pl_indices", format!("{pick:?}"));
        let six_gl_q = java_pls_to_gls(&six_pl);
        let picked: Vec<f64> = pick.iter().map(|&i| six_gl_q[i]).collect();
        java_gls_to_pls(&scale_log_space(&picked))
    };
    kv(
        "java_quantize_then_subset_pl_on_rust_gls",
        fmt_i32(&java_six_then_subset),
    );

    let java_quantized_subset_gls = java_pls_to_gls(&JAVA_PL);
    kv(
        "java_vcf_quantized_gls",
        fmt_vec_f64(&java_quantized_subset_gls),
    );

    let abs_vs_java_quantized: Vec<f64> = subset_gls
        .iter()
        .zip(java_quantized_subset_gls.iter())
        .map(|(r, j)| (r - j).abs())
        .collect();
    kv(
        "abs_diff_rust_raw_vs_java_vcf_quantized_gl",
        fmt_vec_f64(&abs_vs_java_quantized),
    );

    let hom_alt_cont = cont_subset[2];
    let rust_hom_alt_pl = rust_round_subset[2];
    let java_math_hom_alt_pl = java_round_subset[2];
    kv("hom_alt_continuous_pl", fmt_f64(hom_alt_cont));
    kv("integer_boundary_3517_5", "3517.5");
    kv(
        "hom_alt_vs_3517_5",
        format!("{:.17}", hom_alt_cont - 3517.5),
    );

    let rust_round_eq_java_round_subset = rust_round_subset == java_round_subset;
    let rust_round_eq_java_round_merged = rust_round_merged == java_round_merged;
    kv(
        "rust_f64_round_equals_java_math_round_on_subset_gls",
        rust_round_eq_java_round_subset.to_string(),
    );
    kv(
        "rust_f64_round_equals_java_math_round_on_merged_gls",
        rust_round_eq_java_round_merged.to_string(),
    );

    let gls_identical_to_java_quantized = subset_gls
        .iter()
        .zip(java_quantized_subset_gls.iter())
        .all(|(r, j)| r == j);

    // Case classification (first unequal operation on this frozen site).
    let (first_unequal, classification, raw_gls_identical, boundary_crossed) =
        if !rust_round_eq_java_round_subset {
            (
                "f64::round vs Java Math.round of the same unused-ALT log10 GLs",
                "PL_ROUNDING_DIVERGENCE",
                true,
                (hom_alt_cont - 3517.5).abs() < 1.0 && rust_hom_alt_pl != java_math_hom_alt_pl,
            )
        } else if java_six_then_subset == JAVA_PL
            && rust_pl == RUST_PL
            && java_six_then_subset != rust_pl
        {
            (
            "Java calculateGLs GLsToPLs of 6 GLs then PLsToGLs unused-ALT subset vs Rust unused-ALT subset of raw f64 GLs then f64::round",
            "PL_CONVERSION_DIVERGENCE",
            false,
            (hom_alt_cont > 3517.5) != (JAVA_PL[2] as f64 > 3517.5),
        )
        } else if java_round_merged
            .get(subsetted_diploid_pl_indices(&snap.remaining_keep_indices)[2])
            .copied()
            == Some(RUST_PL[2])
            && JAVA_PL[2] != RUST_PL[2]
            && rust_round_eq_java_round_subset
        {
            (
            "diploid calculator log10 GLs for merged 2/2 (and therefore unused-ALT 1/1) already differ from the Java 4.4 GLs that GLsToPLs to 3517",
            "GENOTYPE_LIKELIHOOD_DIVERGENCE",
            false,
            hom_alt_cont >= 3517.5,
        )
        } else if rust_pl != java_round_subset {
            (
                "Rust integer PL representation after f64::round (min-subtract / i32 cast)",
                "PL_REPRESENTATION_DIVERGENCE",
                rust_round_eq_java_round_subset,
                false,
            )
        } else {
            (
                "unclassified: inspect 6R247 kv dump",
                "PL_CONVERSION_DIVERGENCE",
                gls_identical_to_java_quantized,
                false,
            )
        };

    // 6R.303 moves the heterozygote combine. Emitted GLs are differences
    // against that heterozygote, so these bits and the continuous hom-alt
    // move with it. The integer PL remains 3518.
    assert_eq!(
        subset_gls[0].to_bits(),
        0xc04c7d8cfca8f448,
        "lock 0/0 GL bits"
    );
    assert_eq!(
        subset_gls[1].to_bits(),
        0x0000000000000000,
        "lock 0/1 GL = 0"
    );
    assert_eq!(
        subset_gls[2].to_bits(),
        0xc075fc069f8dab28,
        "lock 1/1 GL bits"
    );
    assert!(
        (hom_alt_cont - 3517.51617005722891918).abs() < 1e-12,
        "Rust 1/1 continuous PL is 3517.516…; got {hom_alt_cont:.17}"
    );
    assert!(
        hom_alt_cont > 3517.5,
        "Rust 1/1 is strictly above the Java Math.round 3517.5 boundary"
    );
    assert_eq!(java_round_subset, RUST_PL);
    assert_eq!(java_six_then_subset, RUST_PL);
    assert_ne!(
        java_six_then_subset, JAVA_PL,
        "Java GLsToPLs of Rust calculator GLs is already 3518; Java VCF is 3517"
    );
    assert_eq!(classification, "GENOTYPE_LIKELIHOOD_DIVERGENCE");
    assert!(!raw_gls_identical);
    assert!(boundary_crossed);

    kv("raw_gls_identical", raw_gls_identical.to_string());
    kv("first_unequal_operation", first_unequal);
    kv("classification", classification);
    kv("rounding_boundary_crossed", boundary_crossed.to_string());
    kv("java_numeric_type", "double");
    kv("rust_numeric_type", "f64");

    // GQ independence: het PL[1]==0, PL[0]>=39, PL[2]>=39 → GQ = min(PL[0],PL[2]).clamp(0,99)
    // = 570.min(3518).clamp = 99. Replacing 3518 with 3517 does not change GQ.
    let gq_from_emit = rust_emit.gq.as_i32();
    let gq_from_pl = gq_phred_from_pl(&rust_pl);
    let gq_p7 = gq_phred_p7_compatible(subset_gls, &rust_pl, &call.genotype.format.ad_as_i32());
    let gq_if_java_pl =
        gq_phred_p7_compatible(subset_gls, &JAVA_PL, &call.genotype.format.ad_as_i32());
    kv("gq_emit", gq_from_emit.to_string());
    kv("gq_from_pl_uncapped_gap", gq_from_pl.to_string());
    kv("gq_p7_rust_pl", gq_p7.to_string());
    kv("gq_p7_if_java_pl_3517", gq_if_java_pl.to_string());
    assert_eq!(gq_from_emit, 99);
    assert_eq!(gq_p7, 99);
    assert_eq!(gq_if_java_pl, 99, "GQ cap erases PL[2] ±1");
    kv("gq_affected", "NO");

    let recs =
        try_emit_call_region_variants(region, &outcome, "SAMPLE", DEFAULT_STAND_EMIT_CONFIDENCE)
            .unwrap_or_default();
    let emitted = recs
        .iter()
        .find(|r| {
            r.position == TARGET
                && r.reference == TARGET_REF
                && r.alternate.iter().any(|a| a == TARGET_ALT)
        })
        .expect("emit");
    let sample = emitted.samples.first().expect("sample");
    assert_eq!(sample.pl.as_deref(), Some(&[570, 0, 3518][..]));
    assert_eq!(sample.gq, Some(99.0));
    assert_eq!(sample.ad.as_deref(), Some(&[88u32, 22][..]));
    assert_eq!(sample.dp, Some(110));
    assert_eq!(info_i32(&emitted.info, "DP"), Some(123));
    let qual = emitted.quality.expect("QUAL");
    assert!((qual - 562.60).abs() < 0.01, "QUAL 562.60; got {qual}");

    assert_eq!(rust_round_subset, RUST_PL);
    assert_eq!(
        classification,
        match classification {
            "GENOTYPE_LIKELIHOOD_DIVERGENCE"
            | "PL_CONVERSION_DIVERGENCE"
            | "PL_ROUNDING_DIVERGENCE"
            | "PL_REPRESENTATION_DIVERGENCE" => classification,
            other => panic!("invalid classification {other}"),
        }
    );
    assert_ne!(
        rust_pl[2], JAVA_PL[2],
        "this arrow exists because hom-alt PL still differs"
    );
    assert!(
        first_unequal != "unclassified: inspect 6R247 kv dump",
        "must name the first unequal operation from measured values"
    );
}
