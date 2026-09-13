//! 6R.174 live: INFO DP is Java Coverage evidenceCount (3) while FORMAT DP stays 2.
//! Skipped unless `HOLDOUT_6R174=1`.
//!
//! ```text
//! HOLDOUT_6R174=1 cargo test -p gatk-haplotypecaller --test holdout_6r174_info_dp -- --nocapture --test-threads=1
//! ```

use gatk_core::io::vcf::InfoValue;
use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, traverse_assembly_region_walker,
    try_emit_call_region_variants, AssemblyRegionCallDisposition, CallRegionArgs, GenomePosition,
    HaplotypeCallerEngine, ReadFilterParams, WalkerTraversalConfig, DEFAULT_STAND_EMIT_CONFIDENCE,
};
use std::path::{Path, PathBuf};

const INTERVAL: &str = "2:92305500-92305850";
const BAM_REL: &str = "parity/realworld/na12878_20k_b37/NA12878_20k.b37.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const TARGET: u64 = 92_305_634;
const MERGED_REF: &str = "G";
const MERGED_ALT: &str = "T";

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    println!("6R174\t{key}\t{}", value.as_ref());
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

#[test]
fn holdout_6r174_info_dp() {
    if std::env::var("HOLDOUT_6R174").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R174=1");
        return;
    }
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
    kv("target", "2:92305634 G/T");
    kv(
        "production_change",
        "INFO DP ← Coverage.evidenceCount (retainEvidence unique); FORMAT DP unchanged",
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
        .expect("covering");
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
        .expect("call");
    assert_eq!(call.genotype.format.dp.as_i32(), 2);
    assert_eq!(call.genotype.format.ad_as_i32(), vec![0, 2]);
    let emitted =
        try_emit_call_region_variants(covering, &outcome, "NA12878", DEFAULT_STAND_EMIT_CONFIDENCE)
            .expect("emit");
    let rec = emitted
        .iter()
        .find(|r| r.position == TARGET)
        .expect("record");
    kv(
        "info_dp",
        info_i32(&rec.info, "DP").unwrap_or(-1).to_string(),
    );
    kv(
        "format_dp",
        rec.samples
            .first()
            .and_then(|s| s.dp)
            .map(|d| d.to_string())
            .unwrap_or_default(),
    );
    assert_eq!(info_i32(&rec.info, "DP"), Some(3));
    assert_eq!(
        rec.samples.first().and_then(|s| s.dp).map(|d| d as i32),
        Some(2)
    );
    let sor = info_f64(&rec.info, "SOR").unwrap_or(0.0);
    kv("sor_after_6r176", format!("{sor}"));
    assert!(
        (sor - 0.6931471805599453).abs() < 1e-6 || (sor - 0.693).abs() < 0.002,
        "after 6R.176, SOR is Java 0.693, got {sor}"
    );
}
