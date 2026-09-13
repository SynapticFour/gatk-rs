//! 6R.178 live: InbreedingCoeff omitted at `2:92305634 G/T` when n<10.
//! Skipped unless `HOLDOUT_6R178=1`.
//!
//! ```text
//! HOLDOUT_6R178=1 cargo test -p gatk-haplotypecaller --test holdout_6r178_inbreeding -- --nocapture --test-threads=1
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
    println!("6R178\t{key}\t{}", value.as_ref());
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
fn holdout_6r178_inbreeding() {
    if std::env::var("HOLDOUT_6R178").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R178=1");
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
        "hc_info_values insert InbreedingCoeff iff n_genotypes >= 10",
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
    assert_eq!(call.genotype.format.pl_as_i32(), vec![90, 6, 0]);
    assert_eq!(call.genotype.format.gq.as_i32(), 6);

    let emitted =
        try_emit_call_region_variants(covering, &outcome, "NA12878", DEFAULT_STAND_EMIT_CONFIDENCE)
            .expect("emit");
    let rec = emitted
        .iter()
        .find(|r| r.position == TARGET)
        .expect("record");
    let sample = rec.samples.first().expect("sample");
    assert_eq!(rec.samples.len(), 1);
    assert_eq!(
        sample.gt.as_ref().map(|g| g.to_string()).as_deref(),
        Some("1/1")
    );
    assert_eq!(sample.ad.as_deref(), Some(&[0u32, 2][..]));
    assert_eq!(sample.dp.map(|d| d as i32), Some(2));
    assert_eq!(sample.gq.map(|g| g as i32), Some(6));
    assert_eq!(sample.pl.as_deref(), Some(&[90u32, 6, 0][..]));
    assert_eq!(info_i32(&rec.info, "DP"), Some(3));
    let sor = info_f64(&rec.info, "SOR").unwrap_or(0.0);
    let fs = info_f64(&rec.info, "FS").unwrap_or(-1.0);
    let mq = info_f64(&rec.info, "MQ").unwrap_or(-1.0);
    assert!(
        (sor - 0.6931471805599453).abs() < 1e-6 || (sor - 0.693).abs() < 0.002,
        "6R.176 SOR stays 0.693"
    );
    assert!(fs < 0.02, "FS stays 0");
    assert!((mq - 41.96).abs() < 0.005, "MQ stays 41.96");
    assert!((rec.quality.unwrap_or(0.0) - 78.32).abs() < 0.02);
    assert!(
        !info_has(&rec.info, "InbreedingCoeff"),
        "n=1 < MIN_SAMPLES=10 must omit InbreedingCoeff"
    );
    kv("info_dp", "3");
    kv("format_dp", "2");
    kv("sor", format!("{sor}"));
    kv("inbreeding_coeff", "OMITTED (n=1 < 10)");
    kv("formula", "UNCHANGED");
}
