//! 6R.197 live: after 6R.212, `2:92325268 C/T` INFO DP/MQ/SOR match Java.
//! After 6R.218, first remaining genuine split is `20:29455379 G/A` FORMAT/AD.
//! PRODUCTION CHANGE: NONE.
//! Skipped unless `HOLDOUT_6R197=1`.
//!
//! ```text
//! HOLDOUT_6R197=1 cargo test -p gatk-haplotypecaller --test holdout_6r197_fresh_info -- --nocapture --test-threads=1
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
const TARGET: u64 = 92_307_403;
const CLOSED_AG: u64 = 92_305_635;
const CLOSED_QUAL: u64 = 92_307_359;
const CLOSED_INDEL: u64 = 92_307_324;

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    println!("6R197\t{key}\t{}", value.as_ref());
}

fn unique_indices(likelihoods: &[gatk_haplotypecaller::RegionReadLikelihood]) -> BTreeSet<usize> {
    likelihoods.iter().map(|c| c.read_index.get()).collect()
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
fn holdout_6r197_fresh_info() {
    if std::env::var("HOLDOUT_6R197").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R197=1");
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
    kv("target", "2:92307403 C/A");
    kv(
        "fresh_pass",
        "6R.197 live pin at 2:92307403; 6R.218 closed 20:29455015; first remaining is 20:29455379 FORMAT/AD (record only)",
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
    let covering_for = |pos: u64| {
        tg_regions
            .iter()
            .find(|r| {
                matches!(
                    call_disposition(r),
                    AssemblyRegionCallDisposition::ActiveFull
                ) && r.start.get() <= pos
                    && r.end.get() >= pos
            })
            .unwrap_or_else(|| panic!("covering {pos}"))
    };
    let covering_qual = covering_for(CLOSED_QUAL);
    let covering_indel = covering_for(CLOSED_INDEL);
    let covering_tg = covering_for(TARGET);
    let ag_outcome = HaplotypeCallerEngine::call_region(
        covering_ag,
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("ag call")
    .expect("ag outcome");
    let qual_outcome = HaplotypeCallerEngine::call_region(
        covering_qual,
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("qual call")
    .expect("qual outcome");
    let indel_outcome = HaplotypeCallerEngine::call_region(
        covering_indel,
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("indel call")
    .expect("indel outcome");
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
    let closed_indel = indel_outcome
        .genotyped_calls
        .iter()
        .find(|c| {
            c.event.start_1based == GenomePosition::new_1based(CLOSED_INDEL)
                && c.event.ref_allele == "TTC"
                && c.event.alt_allele == "T"
        })
        .expect("TTC/T");
    let closed_qual = qual_outcome
        .genotyped_calls
        .iter()
        .find(|c| {
            c.event.start_1based == GenomePosition::new_1based(CLOSED_QUAL)
                && c.event.ref_allele == "CT"
                && c.event.alt_allele == "C"
        })
        .expect("CT/C");
    let target = tg_outcome
        .genotyped_calls
        .iter()
        .find(|c| {
            c.event.start_1based == GenomePosition::new_1based(TARGET)
                && c.event.ref_allele == "C"
                && c.event.alt_allele == "A"
        })
        .expect("C/A");

    assert_eq!(unique_indices(&closed_ag.annotation_likelihoods).len(), 3);
    assert_eq!(
        unique_indices(&closed_indel.annotation_likelihoods).len(),
        1
    );
    assert_eq!(unique_indices(&closed_qual.annotation_likelihoods).len(), 0);
    assert_eq!(unique_indices(&target.annotation_likelihoods).len(), 6);
    assert_eq!(target.genotype.format.ad_as_i32(), vec![2, 4]);
    assert_eq!(target.genotype.format.dp.as_i32(), 6);
    assert_eq!(closed_qual.genotype.format.pl_as_i32(), vec![39, 0, 39]);

    let ag_emitted = try_emit_call_region_variants(
        covering_ag,
        &ag_outcome,
        "NA12878",
        DEFAULT_STAND_EMIT_CONFIDENCE,
    )
    .expect("ag emit");
    let qual_emitted = try_emit_call_region_variants(
        covering_qual,
        &qual_outcome,
        "NA12878",
        DEFAULT_STAND_EMIT_CONFIDENCE,
    )
    .expect("qual emit");
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
    let pin_rec = qual_emitted
        .iter()
        .find(|r| r.position == CLOSED_QUAL && r.reference == "CT")
        .expect("CT/C emit");
    let target_rec = tg_emitted
        .iter()
        .find(|r| r.position == TARGET && r.reference == "C")
        .expect("C/A emit");
    assert_eq!(info_i32(&ag_rec.info, "DP"), Some(3));
    let pin_q = pin_rec.quality.expect("pin QUAL");
    assert!((pin_q - 31.60).abs() < 0.005);
    assert_eq!(info_i32(&pin_rec.info, "DP"), Some(2));
    assert!(info_has(&pin_rec.info, "BaseQRankSum"));
    assert!(info_has(&pin_rec.info, "MQRankSum"));
    assert!(info_has(&pin_rec.info, "ReadPosRankSum"));
    assert_eq!(info_i32(&target_rec.info, "DP"), Some(6));
    let bq = info_f64(&target_rec.info, "BaseQRankSum").expect("6R.203 BaseQ");
    assert_eq!((bq * 1000.0).round(), -1834.0);
    let rprs = info_f64(&target_rec.info, "ReadPosRankSum").expect("6R.203 ReadPos");
    assert_eq!((rprs * 1000.0).round(), 1282.0);
    let mqrs = info_f64(&target_rec.info, "MQRankSum").expect("MQRankSum");
    assert_eq!((mqrs * 1000.0).round(), 1834.0);
    kv("format", "GT=0/1 AD=2,4 DP=6 GQ=72 PL=162,0,72");
    kv("java_info_dp", "6");
    kv("rust_info_dp", "6");
    kv(
        "classification",
        "6R.218 closed 20:29455015; first remaining 20:29455379 FORMAT/AD recorded only",
    );
}
