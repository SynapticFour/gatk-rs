//! 6R.185 live: after the last P12 refresh, stored hap likelihoods at
//! `2:92307333 T/G` are Java-order filtered n=2. Skipped unless `HOLDOUT_6R185=1`.
//!
//! ```text
//! HOLDOUT_6R185=1 cargo test -p gatk-haplotypecaller --test holdout_6r185_post_refresh_filter -- --nocapture --test-threads=1
//! ```

use gatk_core::io::vcf::InfoValue;
use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::{
    begin_likelihood_pipeline_observe, call_disposition, flatten_assembly_regions,
    take_likelihood_pipeline_snaps, traverse_assembly_region_walker, try_emit_call_region_variants,
    AssemblyRegionCallDisposition, CallRegionArgs, GenomePosition, HaplotypeCallerEngine,
    ReadFilterParams, WalkerTraversalConfig, DEFAULT_STAND_EMIT_CONFIDENCE,
};
use std::collections::BTreeSet;
use std::path::{Path, PathBuf};

const INTERVAL: &str = "2:92307200-92307550";
const BAM_REL: &str = "parity/realworld/na12878_20k_b37/NA12878_20k.b37.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const TARGET: u64 = 92_307_333;
const MERGED_REF: &str = "T";
const MERGED_ALT: &str = "G";
const CLOSED_INDEL: u64 = 92_307_324;
const CLOSED_SNP: u64 = 92_305_634;

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    println!("6R185\t{key}\t{}", value.as_ref());
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
fn holdout_6r185_post_refresh_filter() {
    if std::env::var("HOLDOUT_6R185").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R185=1");
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
    kv("target", "2:92307333 T/G");
    kv(
        "production_change",
        "Java-order normalize+filter after last P12 refresh",
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
    begin_likelihood_pipeline_observe();
    let outcome = HaplotypeCallerEngine::call_region(
        covering,
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("call")
    .expect("outcome");
    let snaps = take_likelihood_pipeline_snaps();
    let last = snaps.last().expect("snap");
    kv("last_snap_stage", last.stage);
    assert_eq!(last.stage, "filter");

    let stored: BTreeSet<usize> = outcome
        .read_likelihoods
        .iter()
        .map(|c| c.read_index.get())
        .collect();
    kv("stored_n", stored.len().to_string());
    assert_eq!(stored.len(), 2);

    let call = outcome
        .genotyped_calls
        .iter()
        .find(|c| {
            c.event.start_1based == GenomePosition::new_1based(TARGET)
                && c.event.ref_allele == MERGED_REF
                && c.event.alt_allele == MERGED_ALT
        })
        .expect("call");
    assert_eq!(call.genotype.format.dp.as_i32(), 1);
    assert_eq!(call.genotype.format.ad_as_i32(), vec![0, 1]);
    assert_eq!(call.genotype.format.pl_as_i32(), vec![45, 3, 0]);
    assert!(
        call.annotation_likelihoods.is_empty(),
        "6R.181 B not implemented this round"
    );

    let emitted =
        try_emit_call_region_variants(covering, &outcome, "NA12878", DEFAULT_STAND_EMIT_CONFIDENCE)
            .expect("emit");
    let rec = emitted
        .iter()
        .find(|r| {
            r.position == TARGET
                && r.reference == MERGED_REF
                && r.alternate.first().map(String::as_str) == Some(MERGED_ALT)
        })
        .expect("record");
    let sample = rec.samples.first().expect("sample");
    assert_eq!(
        sample.gt.as_ref().map(|g| g.to_string()).as_deref(),
        Some("1/1")
    );
    assert_eq!(sample.ad.as_deref(), Some(&[0u32, 1][..]));
    assert_eq!(sample.dp.map(|d| d as i32), Some(1));
    assert_eq!(sample.gq.map(|g| g as i32), Some(3));
    assert_eq!(sample.pl.as_deref(), Some(&[45u32, 3, 0][..]));
    assert!((rec.quality.unwrap_or(0.0) - 35.48).abs() < 0.05);
    kv(
        "target_info_fallback",
        format!(
            "DP={:?} MQ={:?} SOR={:?} FS={:?}",
            info_i32(&rec.info, "DP"),
            info_f64(&rec.info, "MQ"),
            info_f64(&rec.info, "SOR"),
            info_f64(&rec.info, "FS")
        ),
    );
    assert_eq!(info_i32(&rec.info, "DP"), Some(1));
    let mq = info_f64(&rec.info, "MQ").unwrap_or(-1.0);
    let sor = info_f64(&rec.info, "SOR").unwrap_or(-1.0);
    assert!((mq - 42.05).abs() < 0.005);
    assert!((sor - 0.6931471805599453).abs() < 1e-9);

    let closed_indel = emitted
        .iter()
        .find(|r| r.position == CLOSED_INDEL && r.reference == "TTC")
        .expect("TTC/T");
    assert_eq!(info_i32(&closed_indel.info, "DP"), Some(1));
    let closed_mq = info_f64(&closed_indel.info, "MQ").unwrap_or(-1.0);
    let closed_sor = info_f64(&closed_indel.info, "SOR").unwrap_or(-1.0);
    assert!((closed_mq - 44.0).abs() < 0.005);
    assert!((closed_sor - 1.6094379124341003).abs() < 1e-3 || (closed_sor - 1.609).abs() < 0.002);

    let closed_specs = parse_intervals_cli_string(&dict, "2:92305500-92305850").expect("closed");
    let closed_walk = traverse_assembly_region_walker(
        &dict,
        &closed_specs,
        &ref_fasta,
        &bam,
        &ReadFilterParams::gatk_standard_hc(),
        &WalkerTraversalConfig::gatk_haplotype_caller_production(100),
    )
    .expect("closed walk");
    let closed_regions = flatten_assembly_regions(&closed_walk);
    let closed_covering = closed_regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= CLOSED_SNP
                && r.end.get() >= CLOSED_SNP
        })
        .expect("closed covering");
    let closed_outcome = HaplotypeCallerEngine::call_region(
        closed_covering,
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("closed call")
    .expect("closed outcome");
    let closed_emitted = try_emit_call_region_variants(
        closed_covering,
        &closed_outcome,
        "NA12878",
        DEFAULT_STAND_EMIT_CONFIDENCE,
    )
    .expect("closed emit");
    let closed_rec = closed_emitted
        .iter()
        .find(|r| r.position == CLOSED_SNP)
        .expect("closed G/T");
    let closed_sample = closed_rec.samples.first().expect("sample");
    assert_eq!(closed_sample.dp.map(|d| d as i32), Some(2));
    assert_eq!(info_i32(&closed_rec.info, "DP"), Some(3));
    let closed_snp_sor = info_f64(&closed_rec.info, "SOR").unwrap_or(-1.0);
    assert!((closed_snp_sor - 0.693).abs() < 0.002);
    assert!(!info_has(&closed_rec.info, "InbreedingCoeff"));
}
