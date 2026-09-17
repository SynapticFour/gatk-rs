//! 6R.190 live: first remaining INFO is rust-only ReadPosRankSum at `2:92305716 A/C`.
//! Closed `2:92305635 A/G` stays FORMAT/QUAL/DP/MQ/SOR/annotation n=3.
//! Skipped unless `HOLDOUT_6R190=1`.
//!
//! ```text
//! HOLDOUT_6R190=1 cargo test -p gatk-haplotypecaller --test holdout_6r190_fresh_info -- --nocapture --test-threads=1
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
const BAM_REL: &str = "parity/realworld/na12878_20k_b37/NA12878_20k.b37.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const TARGET: u64 = 92_305_716;
const CLOSED_SNP: u64 = 92_305_634;
const CLOSED_AG: u64 = 92_305_635;

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    println!("6R190\t{key}\t{}", value.as_ref());
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
fn holdout_6r190_fresh_info() {
    if std::env::var("HOLDOUT_6R190").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R190=1");
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
    kv("target", "2:92305716 A/C");
    kv(
        "fresh_pass",
        "6R.190 is a fresh post-6R.189 reconnaissance pass. It does not assume another FORMAT-subset or 6R.181–6R.189 cause.",
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
    let covering_closed = regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= CLOSED_AG
                && r.end.get() >= CLOSED_SNP
        })
        .expect("covering closed");
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
    let closed_outcome = HaplotypeCallerEngine::call_region(
        covering_closed,
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("closed call")
    .expect("closed outcome");
    let outcome = HaplotypeCallerEngine::call_region(
        covering,
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("call")
    .expect("outcome");

    let closed_gt = closed_outcome
        .genotyped_calls
        .iter()
        .find(|c| {
            c.event.start_1based == GenomePosition::new_1based(CLOSED_SNP)
                && c.event.ref_allele == "G"
                && c.event.alt_allele == "T"
        })
        .expect("closed G/T");
    let closed_ag = closed_outcome
        .genotyped_calls
        .iter()
        .find(|c| {
            c.event.start_1based == GenomePosition::new_1based(CLOSED_AG)
                && c.event.ref_allele == "A"
                && c.event.alt_allele == "G"
        })
        .expect("closed A/G");
    let target = outcome
        .genotyped_calls
        .iter()
        .find(|c| {
            c.event.start_1based == GenomePosition::new_1based(TARGET)
                && c.event.ref_allele == "A"
                && c.event.alt_allele == "C"
        })
        .expect("target A/C");

    let closed_ag_ann: BTreeSet<usize> = closed_ag
        .annotation_likelihoods
        .iter()
        .map(|c| c.read_index.get())
        .collect();
    assert_eq!(closed_ag_ann.len(), 3, "6R.189 annotation n=3");
    assert_eq!(closed_ag.genotype.format.ad_as_i32(), vec![0, 2]);
    assert_eq!(closed_ag.genotype.format.dp.as_i32(), 2);
    assert_eq!(closed_ag.genotype.format.gq.as_i32(), 6);
    assert_eq!(closed_ag.genotype.format.pl_as_i32(), vec![90, 6, 0]);
    assert_eq!(target.genotype.format.ad_as_i32(), vec![0, 3]);
    assert_eq!(target.genotype.format.dp.as_i32(), 3);
    assert_eq!(closed_gt.genotype.format.dp.as_i32(), 2);

    let closed_emitted = try_emit_call_region_variants(
        covering_closed,
        &closed_outcome,
        "NA12878",
        DEFAULT_STAND_EMIT_CONFIDENCE,
    )
    .expect("closed emit");
    let emitted =
        try_emit_call_region_variants(covering, &outcome, "NA12878", DEFAULT_STAND_EMIT_CONFIDENCE)
            .expect("emit");
    let closed_gt_rec = closed_emitted
        .iter()
        .find(|r| r.position == CLOSED_SNP && r.reference == "G")
        .expect("G/T emit");
    let closed_ag_rec = closed_emitted
        .iter()
        .find(|r| r.position == CLOSED_AG && r.reference == "A")
        .expect("A/G emit");
    let target_rec = emitted
        .iter()
        .find(|r| r.position == TARGET && r.reference == "A")
        .expect("A/C emit");
    let ag_sample = closed_ag_rec.samples.first().expect("sample");
    assert_eq!(
        ag_sample.gt.as_ref().map(|g| g.to_string()).as_deref(),
        Some("1/1")
    );
    assert_eq!(ag_sample.ad.as_deref(), Some(&[0u32, 2][..]));
    assert_eq!(ag_sample.dp.map(|d| d as i32), Some(2));
    assert!((closed_ag_rec.quality.unwrap_or(0.0) - 78.32).abs() < 0.05);
    assert_eq!(info_i32(&closed_gt_rec.info, "DP"), Some(3));
    assert_eq!(info_i32(&closed_ag_rec.info, "DP"), Some(3));
    assert_eq!(info_i32(&target_rec.info, "DP"), Some(4));
    let ag_mq = info_f64(&closed_ag_rec.info, "MQ").unwrap_or(-1.0);
    let ag_sor = info_f64(&closed_ag_rec.info, "SOR").unwrap_or(-1.0);
    assert!((ag_mq - 41.96).abs() < 0.005);
    assert!((ag_sor - 0.693).abs() < 0.002);
    assert!(!info_has(&target_rec.info, "ReadPosRankSum"));
    assert!(!info_has(&target_rec.info, "BaseQRankSum"));
    assert!(!info_has(&target_rec.info, "MQRankSum"));
    kv(
        "info",
        format!(
            "closed A/G DP=3 MQ={ag_mq} SOR={ag_sor} annotation_n=3; target A/C DP=4 ReadPosRankSum={:?}",
            info_f64(&target_rec.info, "ReadPosRankSum")
        ),
    );
    kv("classification", "A — WRONG SOURCE OBJECT");
}
