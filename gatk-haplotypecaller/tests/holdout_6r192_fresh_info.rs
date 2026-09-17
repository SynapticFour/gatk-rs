//! 6R.192 live: first remaining INFO is Java-only BaseQRankSum/MQRankSum at
//! `2:92307359 CT/C`. Closed 6R.189/6R.191/6R.186 sites stay closed.
//! Skipped unless `HOLDOUT_6R192=1`.
//!
//! ```text
//! HOLDOUT_6R192=1 cargo test -p gatk-haplotypecaller --test holdout_6r192_fresh_info -- --nocapture --test-threads=1
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
const TARGET: u64 = 92_307_359;
const CLOSED_AG: u64 = 92_305_635;
const CLOSED_AC: u64 = 92_305_716;
const CLOSED_TG: u64 = 92_307_333;

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    println!("6R192\t{key}\t{}", value.as_ref());
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
fn holdout_6r192_fresh_info() {
    if std::env::var("HOLDOUT_6R192").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R192=1");
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
    kv(
        "production_change",
        "NONE at 6R.192; BaseQRankSum emit closed by 6R.193; MQRankSum closed by 6R.194",
    );
    kv("target", "2:92307359 CT/C");
    kv(
        "fresh_pass",
        "6R.192 proved Java-only BaseQ at 2:92307359; 6R.193 closed BaseQ emission. 6R.194 closed MQRankSum.",
    );

    let dict = SequenceDictionary::from_fasta_path(&ref_fasta).expect("dict");
    let gap_specs = parse_intervals_cli_string(&dict, GAP_INTERVAL).expect("gap");
    let tg_specs = parse_intervals_cli_string(&dict, TG_INTERVAL).expect("tg");
    let gap_walk = traverse_assembly_region_walker(
        &dict,
        &gap_specs,
        &ref_fasta,
        &bam,
        &ReadFilterParams::gatk_standard_hc(),
        &WalkerTraversalConfig::gatk_haplotype_caller_production(100),
    )
    .expect("gap walk");
    let tg_walk = traverse_assembly_region_walker(
        &dict,
        &tg_specs,
        &ref_fasta,
        &bam,
        &ReadFilterParams::gatk_standard_hc(),
        &WalkerTraversalConfig::gatk_haplotype_caller_production(100),
    )
    .expect("tg walk");
    let gap_regions = flatten_assembly_regions(&gap_walk);
    let tg_regions = flatten_assembly_regions(&tg_walk);
    let covering_ag = gap_regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= CLOSED_AG
                && r.end.get() >= CLOSED_AG
        })
        .expect("covering A/G");
    let covering_ac = gap_regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= CLOSED_AC
                && r.end.get() >= CLOSED_AC
        })
        .expect("covering A/C");
    let covering_tg = tg_regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= TARGET
        })
        .expect("tg covering");
    let ag_outcome = HaplotypeCallerEngine::call_region(
        covering_ag,
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("ag call")
    .expect("ag outcome");
    let ac_outcome = HaplotypeCallerEngine::call_region(
        covering_ac,
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("ac call")
    .expect("ac outcome");
    let tg_outcome = HaplotypeCallerEngine::call_region(
        covering_tg,
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("tg call")
    .expect("tg outcome");

    let closed_ag = ag_outcome
        .genotyped_calls
        .iter()
        .find(|c| {
            c.event.start_1based == GenomePosition::new_1based(CLOSED_AG)
                && c.event.ref_allele == "A"
                && c.event.alt_allele == "G"
        })
        .expect("A/G");
    let target = tg_outcome
        .genotyped_calls
        .iter()
        .find(|c| {
            c.event.start_1based == GenomePosition::new_1based(TARGET)
                && c.event.ref_allele == "CT"
                && c.event.alt_allele == "C"
        })
        .expect("CT/C");
    let tg_call = tg_outcome
        .genotyped_calls
        .iter()
        .find(|c| {
            c.event.start_1based == GenomePosition::new_1based(CLOSED_TG)
                && c.event.ref_allele == "T"
                && c.event.alt_allele == "G"
        })
        .expect("T/G");
    let ag_n: BTreeSet<usize> = closed_ag
        .annotation_likelihoods
        .iter()
        .map(|c| c.read_index.get())
        .collect();
    let tg_n: BTreeSet<usize> = tg_call
        .annotation_likelihoods
        .iter()
        .map(|c| c.read_index.get())
        .collect();
    assert_eq!(ag_n.len(), 3);
    assert_eq!(tg_n.len(), 1);
    assert_eq!(target.genotype.format.ad_as_i32(), vec![1, 1]);
    assert_eq!(target.genotype.format.dp.as_i32(), 2);

    let ag_emitted = try_emit_call_region_variants(
        covering_ag,
        &ag_outcome,
        "NA12878",
        DEFAULT_STAND_EMIT_CONFIDENCE,
    )
    .expect("ag emit");
    let ac_emitted = try_emit_call_region_variants(
        covering_ac,
        &ac_outcome,
        "NA12878",
        DEFAULT_STAND_EMIT_CONFIDENCE,
    )
    .expect("ac emit");
    let tg_emitted = try_emit_call_region_variants(
        covering_tg,
        &tg_outcome,
        "NA12878",
        DEFAULT_STAND_EMIT_CONFIDENCE,
    )
    .expect("tg emit");
    let ag_rec = ag_emitted
        .iter()
        .find(|r| r.position == CLOSED_AG && r.reference == "A")
        .expect("A/G emit");
    let ac_rec = ac_emitted
        .iter()
        .find(|r| r.position == CLOSED_AC && r.reference == "A")
        .expect("A/C emit");
    let target_rec = tg_emitted
        .iter()
        .find(|r| r.position == TARGET && r.reference == "CT")
        .expect("CT/C emit");
    let sample = target_rec.samples.first().expect("sample");
    assert_eq!(
        sample.gt.as_ref().map(|g| g.to_string()).as_deref(),
        Some("0/1")
    );
    assert_eq!(sample.ad.as_deref(), Some(&[1u32, 1][..]));
    assert_eq!(info_i32(&ag_rec.info, "DP"), Some(3));
    assert!((info_f64(&ag_rec.info, "MQ").unwrap_or(-1.0) - 41.96).abs() < 0.005);
    assert!(!info_has(&ac_rec.info, "ReadPosRankSum"));
    assert!(info_has(&target_rec.info, "ReadPosRankSum"));
    let bq = info_f64(&target_rec.info, "BaseQRankSum").unwrap_or(-1.0);
    assert!(
        info_has(&target_rec.info, "BaseQRankSum") && bq.abs() < 0.0005,
        "6R.193 emits BaseQRankSum finite zero"
    );
    let mqrs = info_f64(&target_rec.info, "MQRankSum").expect("6R.194 MQRankSum");
    assert_eq!((mqrs * 1000.0).round(), -674.0);
    kv("format", "GT=0/1 AD=1,1 DP=2");
    kv(
        "classification",
        "F — WRONG EMISSION PREDICATE (BaseQ closed by 6R.193; MQ closed by 6R.194)",
    );
}
