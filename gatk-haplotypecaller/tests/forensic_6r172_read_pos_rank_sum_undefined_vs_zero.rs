//! 6R.172: ReadPosRankSum undefined-vs-zero emit predicate.
//!
//! Frozen Java 4.4.0.0 SHA `2dbc025821bc5f686c423ff332a41e6cef892a77`, GKL 0.8.8.
//! Production: `read_pos_rank_sum -> Option<f64>`; HC adapter inserts only `Some(z)`
//! including `Some(0.0)`. Evidence source stays 6R.170 (`region.reads` pileup).
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --lib forensic_6r172_semantic_matrix_production_fn -- --test-threads=1 --nocapture
//! cargo test -p gatk-haplotypecaller --test forensic_6r172_read_pos_rank_sum_undefined_vs_zero -- --test-threads=1 --nocapture
//! HOLDOUT_6R172=1 cargo test -p gatk-haplotypecaller --test holdout_6r172_read_pos_rank_sum_undefined -- --nocapture --test-threads=1
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
const CASE_A_Z: f64 = 1.3829941271006383;

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    eprintln!("6R172\t{key}\t{}", value.as_ref());
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

fn info_has(info: &[InfoValue], key: &str) -> bool {
    info.iter().any(|v| match v {
        InfoValue::Flag(k)
        | InfoValue::Integer(k, _)
        | InfoValue::Float(k, _)
        | InfoValue::String(k, _)
        | InfoValue::Character(k, _) => k == key,
    })
}

/// Adapter contract: Some(z) including 0.0 emits; None omits. Never omit-on-zero.
fn rust_inserts_read_pos(z: Option<f64>) -> bool {
    z.is_some()
}

#[test]
fn forensic_6r172_read_pos_rank_sum_undefined_vs_zero() {
    let plugin = include_str!("../src/annotator/plugins/read_pos_rank_sum.rs");
    let production_plugin = plugin
        .split("#[cfg(test)]")
        .next()
        .expect("production plugin body");
    let ann = include_str!("../src/variant_site_hc_annotations.rs");
    let emit = include_str!("../src/region_vcf_emit.rs");

    kv(
        "return_type",
        "read_pos_rank_sum -> Option<f64> (Some(0.0) != None)",
    );
    assert!(
        production_plugin.contains(
            "pub fn read_pos_rank_sum(ref_positions: &[f64], alt_positions: &[f64]) -> Option<f64>"
        ),
        "production RankSum is Option<f64>"
    );
    assert!(
        production_plugin.contains("if alt.is_empty() || reference.is_empty()")
            && production_plugin.contains("return None"),
        "empty REF or ALT → None"
    );
    assert!(
        production_plugin.contains("if result.z.is_nan()") && production_plugin.contains("None"),
        "MannWhitneyU NaN → None"
    );
    assert!(
        !production_plugin.contains("return 0.0"),
        "must not use 0.0 as the undefined sentinel"
    );
    assert!(
        emit.contains("if let Some(z) = ann.read_pos_rank_sum"),
        "ReadPosRankSum insert is Some-only"
    );
    assert!(
        !emit.contains("if z == 0") && !emit.contains("if z == 0.0") && !emit.contains("z.abs() <"),
        "must not omit legitimate zero in hc_info_values"
    );
    assert!(
        emit.contains("InfoValue::Float(\"FS\".to_string(), vec![ann.fs])"),
        "FS insert stays unconditional (FS=0 must remain)"
    );
    assert!(
        ann.contains("read_offset_evidence_at_site")
            && ann.contains("read_pos_rank_sum: Option<f64>"),
        "evidence source unchanged; field is Option"
    );
    assert!(
        !production_plugin.contains("92316347")
            && !ann.contains("92316347")
            && !emit.contains("92316347"),
        "no coordinate-specific ReadPosRankSum rule"
    );
    assert!(
        plugin.contains("forensic_6r172_semantic_matrix_production_fn"),
        "runtime A–E matrix must call the production function"
    );

    assert!(rust_inserts_read_pos(Some(CASE_A_Z)), "A emit");
    assert!(rust_inserts_read_pos(Some(0.0)), "B emit finite zero");
    assert!(!rust_inserts_read_pos(None), "C/D/E omit");
    kv("case_A", "REF=[10,20] ALT=[30,40] Some(nonzero) insert=yes");
    kv("case_B", "REF=[10,20] ALT=[10,20] Some(0.0) insert=yes");
    kv(
        "case_C",
        "REF=[] ALT=[51,61,70] None insert=no (Java-equivalent target)",
    );
    kv("case_D", "REF=[10,20] ALT=[] None insert=no");
    kv("case_E", "REF=[] ALT=[] None insert=no");
}

#[test]
fn forensic_6r172_live_target_format_and_fs_zero_preserved() {
    let root = repo_root();
    let ref_fasta = root.join(REF_REL);
    let bam = root.join(BAM_REL);
    if !ref_fasta.is_file() || !bam.is_file() {
        eprintln!("skip: missing P12 ref/BAM");
        return;
    }

    kv(
        "java_pin",
        "4.4.0.0 / 2dbc025821bc5f686c423ff332a41e6cef892a77",
    );
    kv("variant", format!("2:{TARGET} G/A"));
    kv(
        "production_change",
        "read_pos_rank_sum Option<f64> + Some-only INFO insert",
    );
    kv("evidence_source", "unchanged region.reads pileup (6R.170)");
    kv(
        "java_equivalent",
        "REF=0 ALT=3 → None → no ReadPosRankSum (plugin Case C; not live pileup)",
    );

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
    assert_eq!(call.genotype.format.gq.as_i32(), 9);

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
    let fs = info_f64(&rec.info, "FS");
    let sor = info_f64(&rec.info, "SOR").unwrap_or(-1.0);
    let mq = info_f64(&rec.info, "MQ").unwrap_or(-1.0);
    kv(
        "rust_final",
        format!(
            "ReadPosRankSum={rp:?} FS={fs:?} SOR={sor} MQ={mq} QUAL={:?} keys={:?}",
            rec.quality,
            rec.info
                .iter()
                .map(|v| match v {
                    InfoValue::Flag(k)
                    | InfoValue::Integer(k, _)
                    | InfoValue::Float(k, _)
                    | InfoValue::String(k, _)
                    | InfoValue::Character(k, _) => k.as_str(),
                })
                .collect::<Vec<_>>()
        ),
    );
    assert!(
        info_has(&rec.info, "FS"),
        "FS=0 must still be inserted; a global skip-zero would be a regression"
    );
    assert!(fs.unwrap_or(-1.0) < 0.02);
    assert!((sor - 1.179).abs() < 0.001);
    assert!((mq - 40.25).abs() < 1e-12);
    assert!((rec.quality.unwrap_or(0.0) - 121.84).abs() < 0.02);
    assert!(
        !info_has(&rec.info, "InbreedingCoeff"),
        "6R.178: n=1 omits InbreedingCoeff"
    );
    kv(
        "live_pileup_note",
        "live ReadPosRankSum follows pileup (6R.170), not Java-equivalent REF=0/ALT=3",
    );
}
