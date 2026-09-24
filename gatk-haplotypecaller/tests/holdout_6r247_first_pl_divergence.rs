//! 6R.247 live: `20:29455649 T/TGTTTG` FORMAT PL remains
//! Java `570,0,3517` vs Rust `570,0,3518`. GQ/QUAL/AD/DP unchanged.
//! Proof-only. Skipped unless `HOLDOUT_6R247=1`.
//!
//! ```text
//! HOLDOUT_6R247=1 cargo test -p gatk-haplotypecaller --test holdout_6r247_first_pl_divergence -- --nocapture --test-threads=1
//! ```

use gatk_core::io::vcf::InfoValue;
use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::hc_genotyping_engine::HcGenotypingConfig;
use gatk_haplotypecaller::read_threading_assembler::DEFAULT_NUM_BEST_HAPLOTYPES_PER_GRAPH;
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, traverse_assembly_region_walker,
    try_emit_call_region_variants, AssemblyRegionCallDisposition, CallRegionArgs,
    HaplotypeCallerEngine, ReadFilterParams, WalkerTraversalConfig, DEFAULT_STAND_CALL_CONF,
    DEFAULT_STAND_EMIT_CONFIDENCE,
};
use std::path::{Path, PathBuf};

const INTERVAL: &str = "20:29455000-29456500";
const BAM_REL: &str = "parity/giab/runs/local-pairhmm-diff/HG001.20-29455000-29456500.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const JAVA_VCF_REL: &str = "parity/reports/6r43/chr20_tiny/java.vcf";
const TARGET: u64 = 29_455_649;
const TARGET_REF: &str = "T";
const TARGET_ALT: &str = "TGTTTG";
const CLOSED_GC: u64 = 29_455_314;

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    println!("6R247\t{key}\t{}", value.as_ref());
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

fn vcf_has(path: &Path, pos: u64, r: &str, a: &str) -> bool {
    for line in std::fs::read_to_string(path).unwrap_or_default().lines() {
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let f: Vec<_> = line.split('\t').collect();
        if f.len() < 5 || f[0] != "20" {
            continue;
        }
        if f[1].parse::<u64>().ok() == Some(pos) && f[3] == r && f[4] == a {
            return true;
        }
    }
    false
}

#[test]
fn holdout_6r247_first_pl_divergence() {
    if std::env::var("HOLDOUT_6R247").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R247=1");
        return;
    }
    assert_eq!(DEFAULT_NUM_BEST_HAPLOTYPES_PER_GRAPH, 128);
    assert_eq!(
        HcGenotypingConfig::strict_java().stand_emit_confidence,
        DEFAULT_STAND_CALL_CONF
    );
    let root = repo_root();
    let java_vcf = root.join(JAVA_VCF_REL);
    assert!(vcf_has(&java_vcf, TARGET, TARGET_REF, TARGET_ALT));
    assert!(!vcf_has(&java_vcf, CLOSED_GC, "G", "C"));
    let java_line = std::fs::read_to_string(&java_vcf)
        .unwrap()
        .lines()
        .find(|l| l.starts_with("20\t29455649\t") && l.contains("\tT\tTGTTTG\t"))
        .expect("java line")
        .to_string();
    assert!(java_line.contains("570,0,3517"));
    kv("java_pl", "570,0,3517");

    let ref_fasta = root.join(REF_REL);
    let bam = root.join(BAM_REL);
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
    kv(
        "rust_pl",
        format!(
            "{}",
            sample
                .pl
                .as_ref()
                .map(|p| p
                    .iter()
                    .map(|x| x.to_string())
                    .collect::<Vec<_>>()
                    .join(","))
                .unwrap_or_default()
        ),
    );
    assert_eq!(sample.pl.as_deref(), Some(&[570, 0, 3518][..]));
    assert_eq!(sample.gq, Some(99.0));
    assert_eq!(sample.ad.as_deref(), Some(&[88u32, 22][..]));
    assert_eq!(sample.dp, Some(110));
    assert_eq!(
        sample.gt.as_ref().map(|g| g.alleles.as_slice()),
        Some(&[0, 1][..])
    );
    assert_eq!(info_i32(&emitted.info, "DP"), Some(123));
    let qual = emitted.quality.expect("QUAL");
    assert!((qual - 562.60).abs() < 0.01);
    kv("production_change", "NONE");
    kv("closed", "GT/AD/DP/GQ/QUAL/INFO DP; remaining FORMAT PL ±1");
}
