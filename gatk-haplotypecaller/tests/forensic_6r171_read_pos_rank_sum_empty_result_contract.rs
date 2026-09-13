//! 6R.171: ReadPosRankSum empty/NaN emission predicate (proof). Closed by 6R.172.
//!
//! Frozen Java 4.4.0.0 SHA `2dbc025821bc5f686c423ff332a41e6cef892a77`, GKL 0.8.8.
//! 6R.171 proved empty/NaN collapsed to `0.0`. 6R.172 restored `Option<f64>`.
//! Evidence source stays 6R.170 (`region.reads`).
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r171_read_pos_rank_sum_empty_result_contract -- --test-threads=1 --nocapture
//! HOLDOUT_6R171=1 cargo test -p gatk-haplotypecaller --test holdout_6r171_read_pos_rank_sum_empty -- --nocapture --test-threads=1
//! ```

use gatk_core::io::vcf::InfoValue;
use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, traverse_assembly_region_walker,
    try_emit_call_region_variants, AssemblyRegionCallDisposition, CallRegionArgs, GenomePosition,
    HaplotypeCallerEngine, ReadFilterParams, WalkerTraversalConfig, DEFAULT_STAND_EMIT_CONFIDENCE,
};
use std::path::{Path, PathBuf};

const INTERVAL: &str = "2:92316200-92316580";
const BAM_REL: &str = "parity/realworld/na12878_20k_b37/NA12878_20k.b37.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const TARGET: u64 = 92_316_347;
const MERGED_REF: &str = "G";
const MERGED_ALT: &str = "A";
/// Observed Case A from production `read_pos_rank_sum(&[10,20], &[30,40])`.
const CASE_A_Z: f64 = 1.3829941271006383;

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    eprintln!("6R171\t{key}\t{}", value.as_ref());
}

fn info_f64(info: &[InfoValue], key: &str) -> Option<f64> {
    for v in info {
        if let InfoValue::Float(k, xs) = v {
            if k == key {
                return xs.first().copied();
            }
        }
    }
    None
}

/// 6R.172 insert predicate: Some(z) including 0.0 emits; None omits.
fn rust_inserts_read_pos(z: Option<f64>) -> bool {
    z.is_some()
}

#[test]
fn forensic_6r171_read_pos_rank_sum_empty_result_contract() {
    let plugin = include_str!("../src/annotator/plugins/read_pos_rank_sum.rs");
    let production_plugin = plugin
        .split("#[cfg(test)]")
        .next()
        .expect("production plugin body");
    let ann = include_str!("../src/variant_site_hc_annotations.rs");
    let emit = include_str!("../src/region_vcf_emit.rs");
    let mwu = include_str!("../src/mann_whitney_u.rs");

    kv(
        "return_type",
        "6R.172 closed: read_pos_rank_sum -> Option<f64>",
    );
    assert!(
        production_plugin.contains(
            "pub fn read_pos_rank_sum(ref_positions: &[f64], alt_positions: &[f64]) -> Option<f64>"
        ),
        "6R.172 production RankSum is Option<f64>"
    );
    assert!(
        production_plugin.contains("if alt.is_empty() || reference.is_empty()")
            && production_plugin.contains("return None"),
        "empty either list → None"
    );
    assert!(
        production_plugin.contains("if result.z.is_nan()") && production_plugin.contains("None"),
        "MannWhitneyU NaN → None"
    );
    assert!(
        !production_plugin.contains("unwrap_or") && !production_plugin.contains("is_finite"),
        "no unwrap_or collapse; no omit-on-!is_finite besides NaN→None"
    );
    assert!(
        mwu.contains("if n1 == 0 || n2 == 0") && mwu.contains("f64::NAN"),
        "MannWhitneyU.test already returns NaN on empty (Java-equal)"
    );
    assert!(
        ann.contains("let rp = read_pos_rank_sum::read_pos_rank_sum")
            && ann.contains("read_pos_rank_sum: rp")
            && ann.contains("read_pos_rank_sum: Option<f64>"),
        "adapter stores Option<f64>"
    );
    assert!(
        emit.contains("if let Some(z) = ann.read_pos_rank_sum"),
        "hc_info_values inserts only Some(z)"
    );
    assert!(
        !emit.contains("if z == 0") && !emit.contains("if z == 0.0"),
        "must not omit legitimate zero"
    );
    assert!(
        !production_plugin.contains("92316347")
            && !ann.contains("92316347")
            && !emit.contains("92316347"),
        "no coordinate-specific ReadPosRankSum rule"
    );

    kv(
        "first_divergent",
        "6R.171: rank_sum_z empty OR → 0.0; 6R.172: empty OR → None",
    );
    kv(
        "java",
        "RankSumTest: empty AND → emptyMap; else MannWhitneyU NaN → emptyMap; finite 0.0 → format %.3f insert",
    );
    kv(
        "case_A",
        format!("REF=[10,20] ALT=[30,40] rust=Some({CASE_A_Z}) insert=yes java=emit"),
    );
    kv(
        "case_B",
        "REF=[10,20] ALT=[10,20] rust=Some(0.0) insert=yes java=emit",
    );
    kv(
        "case_C",
        "REF=[] ALT=[51,61,70] rust=None insert=no java=omit (target-equivalent)",
    );
    kv("case_D", "REF=[10,20] ALT=[] rust=None insert=no java=omit");
    kv("case_E", "REF=[] ALT=[] rust=None insert=no java=omit");

    assert!(rust_inserts_read_pos(Some(CASE_A_Z)));
    assert!(
        rust_inserts_read_pos(Some(0.0)),
        "case B legitimate 0.0 remains emittable"
    );
    assert!(!rust_inserts_read_pos(None), "undefined must not insert");
    assert!(
        plugin.contains("forensic_6r172_semantic_matrix_production_fn"),
        "runtime matrix moved to 6R.172"
    );

    let bq = include_str!("../src/annotator/plugins/rank_sum_baseq.rs");
    let mqrs = include_str!("../src/annotator/plugins/mapping_quality_rank_sum.rs");
    assert!(
        bq.contains("return 0.0") && mqrs.contains("return 0.0"),
        "sibling rank-sum plugins were not changed this round"
    );
    assert!(
        !emit.contains("InfoValue::Float(\"BaseQRankSum\"")
            && !emit.contains("InfoValue::Float(\"MQRankSum\""),
        "hc_info_values does not emit BaseQRankSum/MQRankSum"
    );
    assert!(
        !include_str!("../src/annotator/engine.rs").contains("ReadPosRankSum"),
        "parity-v1 annotator engine is not on this path"
    );
}

#[test]
fn forensic_6r171_live_target_still_inserts_zero() {
    let root = repo_root();
    let ref_fasta = root.join(REF_REL);
    let bam = root.join(BAM_REL);
    if !ref_fasta.is_file() || !bam.is_file() {
        eprintln!("skip: missing P12 ref/BAM");
        return;
    }

    kv("production_change", "6R.172 Option<f64> emit predicate");
    kv("evidence_source", "unchanged region.reads (6R.170)");

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
    let covering = regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= TARGET
                && r.end.get() >= TARGET
        })
        .expect("ActiveFull covering target");

    let outcome = HaplotypeCallerEngine::call_region(
        covering,
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("call")
    .expect("outcome");

    let call = outcome
        .genotyped_calls
        .iter()
        .find(|c| {
            c.event.start_1based == GenomePosition::new_1based(TARGET)
                && c.event.ref_allele == MERGED_REF
                && c.event.alt_allele == MERGED_ALT
        })
        .expect("genotyped G/A");
    assert_eq!(call.genotype.format.ad_as_i32(), vec![0, 3]);
    assert_eq!(call.genotype.format.pl_as_i32(), vec![135, 9, 0]);

    let emitted =
        try_emit_call_region_variants(covering, &outcome, "SAMPLE", DEFAULT_STAND_EMIT_CONFIDENCE)
            .expect("emit");
    let rec = emitted
        .iter()
        .find(|r| {
            r.position == TARGET
                && r.reference == MERGED_REF
                && r.alternate.first().map(String::as_str) == Some(MERGED_ALT)
        })
        .expect("emitted G/A");
    let rp = info_f64(&rec.info, "ReadPosRankSum");
    let fs = info_f64(&rec.info, "FS").unwrap_or(-1.0);
    let sor = info_f64(&rec.info, "SOR").unwrap_or(-1.0);
    let mq = info_f64(&rec.info, "MQ").unwrap_or(-1.0);
    kv(
        "rust_final",
        format!(
            "ReadPosRankSum={rp:?} FS={fs} SOR={sor} MQ={mq} QUAL={:?}",
            rec.quality
        ),
    );
    assert_eq!(rp, Some(0.0), "production still inserts 0.0");
    assert!(fs < 0.02);
    assert!((sor - 1.179).abs() < 0.001);
    assert!((mq - 40.25).abs() < 1e-12);
    assert!((rec.quality.unwrap_or(0.0) - 121.84).abs() < 0.02);
    assert!(
        info_f64(&rec.info, "InbreedingCoeff").is_none(),
        "6R.178: n=1 omits InbreedingCoeff"
    );
}
