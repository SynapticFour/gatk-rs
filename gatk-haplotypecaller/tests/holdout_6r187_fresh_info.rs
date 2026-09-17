//! 6R.187 live: first remaining INFO DP at `2:92305635 A/G` after 6R.186.
//! Skipped unless `HOLDOUT_6R187=1`.
//!
//! ```text
//! HOLDOUT_6R187=1 cargo test -p gatk-haplotypecaller --test holdout_6r187_fresh_info -- --nocapture --test-threads=1
//! ```

use gatk_core::io::vcf::InfoValue;
use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, traverse_assembly_region_walker,
    try_emit_call_region_variants, AssemblyRegionCallDisposition, CallRegionArgs, GenomePosition,
    HaplotypeCallerEngine, ReadFilterParams, WalkerTraversalConfig, DEFAULT_STAND_EMIT_CONFIDENCE,
};
use std::collections::BTreeSet;
use std::path::{Path, PathBuf};

const GAP_INTERVAL: &str = "2:92305500-92305850";
const TG_INTERVAL: &str = "2:92307200-92307550";
const BAM_REL: &str = "parity/realworld/na12878_20k_b37/NA12878_20k.b37.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const TARGET: u64 = 92_305_635;
const CLOSED_SNP: u64 = 92_305_634;
const CLOSED_TG: u64 = 92_307_333;

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    println!("6R187\t{key}\t{}", value.as_ref());
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

fn info_has(info: &[InfoValue], key: &str) -> bool {
    info.iter().any(|v| match v {
        InfoValue::Flag(k)
        | InfoValue::Integer(k, _)
        | InfoValue::Float(k, _)
        | InfoValue::String(k, _)
        | InfoValue::Character(k, _) => k == key,
    })
}

#[test]
fn holdout_6r187_fresh_info() {
    if std::env::var("HOLDOUT_6R187").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R187=1");
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
    kv("production_change", "NONE");
    kv("target", "2:92305635 A/G");
    kv(
        "fresh_pass",
        "6R.187 is a fresh post-6R.186 reconnaissance pass. It does not assume that the remaining annotation differences share the 6R.181/6R.186 cause.",
    );

    let dict = SequenceDictionary::from_fasta_path(&ref_fasta).expect("dict");
    let specs = parse_intervals_cli_string(&dict, GAP_INTERVAL).expect("interval");
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

    let closed = outcome
        .genotyped_calls
        .iter()
        .find(|c| {
            c.event.start_1based == GenomePosition::new_1based(CLOSED_SNP)
                && c.event.ref_allele == "G"
                && c.event.alt_allele == "T"
        })
        .expect("closed");
    let target = outcome
        .genotyped_calls
        .iter()
        .find(|c| {
            c.event.start_1based == GenomePosition::new_1based(TARGET)
                && c.event.ref_allele == "A"
                && c.event.alt_allele == "G"
        })
        .expect("target");
    assert_eq!(closed.genotype.format.dp.as_i32(), 2);
    assert_eq!(target.genotype.format.dp.as_i32(), 2);
    assert_eq!(target.genotype.format.ad_as_i32(), vec![0, 2]);
    assert_eq!(target.genotype.format.pl_as_i32(), vec![90, 6, 0]);
    let closed_ann: BTreeSet<usize> = closed
        .annotation_likelihoods
        .iter()
        .map(|c| c.read_index.get())
        .collect();
    let target_ann: BTreeSet<usize> = target
        .annotation_likelihoods
        .iter()
        .map(|c| c.read_index.get())
        .collect();
    assert!(
        closed.annotation_likelihoods.is_empty() || closed_ann.len() == 3,
        "closed 92305634 is not the n=1 object"
    );
    assert_eq!(target_ann.len(), 3, "6R.189 loc-loop n=3");

    let emitted =
        try_emit_call_region_variants(covering, &outcome, "NA12878", DEFAULT_STAND_EMIT_CONFIDENCE)
            .expect("emit");
    let closed_rec = emitted
        .iter()
        .find(|r| r.position == CLOSED_SNP && r.reference == "G")
        .expect("closed emit");
    let target_rec = emitted
        .iter()
        .find(|r| r.position == TARGET && r.reference == "A")
        .expect("target emit");
    let closed_sample = closed_rec.samples.first().expect("sample");
    let target_sample = target_rec.samples.first().expect("sample");
    assert_eq!(closed_sample.dp.map(|d| d as i32), Some(2));
    assert_eq!(target_sample.dp.map(|d| d as i32), Some(2));
    assert_eq!(info_i32(&closed_rec.info, "DP"), Some(3));
    assert_eq!(info_i32(&target_rec.info, "DP"), Some(3));
    let closed_mq = info_f64(&closed_rec.info, "MQ").unwrap_or(-1.0);
    let target_mq = info_f64(&target_rec.info, "MQ").unwrap_or(-1.0);
    let closed_sor = info_f64(&closed_rec.info, "SOR").unwrap_or(-1.0);
    let target_sor = info_f64(&target_rec.info, "SOR").unwrap_or(-1.0);
    assert!((closed_mq - 41.96).abs() < 0.005);
    assert!((target_mq - 41.96).abs() < 0.005);
    assert!((closed_sor - 0.693).abs() < 0.002);
    assert!((target_sor - 0.693).abs() < 0.002);
    assert!(!info_has(&closed_rec.info, "InbreedingCoeff"));
    kv(
        "info",
        format!("closed DP=3 MQ={closed_mq} SOR={closed_sor}; target DP=3 MQ={target_mq} SOR={target_sor}"),
    );

    let tg_specs = parse_intervals_cli_string(&dict, TG_INTERVAL).expect("tg");
    let tg_walk = traverse_assembly_region_walker(
        &dict,
        &tg_specs,
        &ref_fasta,
        &bam,
        &ReadFilterParams::gatk_standard_hc(),
        &WalkerTraversalConfig::gatk_haplotype_caller_production(100),
    )
    .expect("tg walk");
    let tg_regions = flatten_assembly_regions(&tg_walk);
    let tg_covering = tg_regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= CLOSED_TG
                && r.end.get() >= CLOSED_TG
        })
        .expect("tg covering");
    let tg_outcome = HaplotypeCallerEngine::call_region(
        tg_covering,
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("tg call")
    .expect("tg outcome");
    let tg_call = tg_outcome
        .genotyped_calls
        .iter()
        .find(|c| {
            c.event.start_1based == GenomePosition::new_1based(CLOSED_TG)
                && c.event.ref_allele == "T"
                && c.event.alt_allele == "G"
        })
        .expect("T/G");
    let tg_ann: BTreeSet<usize> = tg_call
        .annotation_likelihoods
        .iter()
        .map(|c| c.read_index.get())
        .collect();
    assert_eq!(tg_ann.len(), 1, "6R.186 n=1 stays");
    let tg_emitted = try_emit_call_region_variants(
        tg_covering,
        &tg_outcome,
        "NA12878",
        DEFAULT_STAND_EMIT_CONFIDENCE,
    )
    .expect("tg emit");
    let tg_rec = tg_emitted
        .iter()
        .find(|r| r.position == CLOSED_TG && r.reference == "T")
        .expect("tg record");
    assert_eq!(info_i32(&tg_rec.info, "DP"), Some(1));
    let tg_mq = info_f64(&tg_rec.info, "MQ").unwrap_or(-1.0);
    let tg_sor = info_f64(&tg_rec.info, "SOR").unwrap_or(-1.0);
    assert!((tg_mq - 44.0).abs() < 0.005);
    assert!((tg_sor - 1.609).abs() < 0.002);
    kv(
        "closed_92307333",
        format!("annotation_n=1 DP=1 MQ={tg_mq} SOR={tg_sor}"),
    );
}
