//! 6R.252 live: Java 4.4.0.0 vs Rust TGTTTG haplotypes at `20:29455649`.
//! Proof-only. Skipped unless `HOLDOUT_6R252=1`.
//!
//! ```text
//! HOLDOUT_6R252=1 cargo test -p gatk-haplotypecaller --test holdout_6r252_java_actual_haplotype_comparison -- --nocapture --test-threads=1
//! ```

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::hc_genotyping_engine::{
    take_colocated_merge_numerics, HcGenotypingConfig,
};
use gatk_haplotypecaller::read_threading_assembler::DEFAULT_NUM_BEST_HAPLOTYPES_PER_GRAPH;
use gatk_haplotypecaller::{
    begin_hap_list_observe, call_disposition, flatten_assembly_regions, take_hap_list_trim_span,
    traverse_assembly_region_walker, try_emit_call_region_variants, AssemblyRegionCallDisposition,
    CallRegionArgs, HaplotypeCallerEngine, HcLikelihoodEngineConfig, PairHmmBackend,
    ReadFilterParams, WalkerTraversalConfig, DEFAULT_STAND_CALL_CONF,
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
const FROZEN_IDX: [usize; 5] = [19, 20, 21, 22, 23];
const RUST_TRIMMED_FNV: [&str; 5] = [
    "55012fcf3b430591",
    "341b2e070ccb5846",
    "6a65e4c02733c2ed",
    "fa07750bf228b3c2",
    "79451c576721a729",
];
const JAVA_TSV: &str = include_str!("forensic_6r252_java_hap_trim.tsv");
const JAVA_TRIMMED_FNV: [&str; 5] = [
    "fcd72c6ce600dd16",
    "c1a4e8204522f645",
    "eb03271fa7548f26",
    "7dc2a8ae5da116e1",
    "343e7c443c1d9732",
];
const JAVA_UNTRIMMED_FNV: [&str; 5] = [
    "a1ed2811a3912892",
    "745492ecfc910685",
    "dae45a854e2d04de",
    "8bdf55500c4f07d1",
    "5cbde2ec4e0a98ea",
];

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    println!("6R252\t{key}\t{}", value.as_ref());
}

fn fnv1a64_hex(bases: &[u8]) -> String {
    let mut h: u64 = 0xcbf29ce484222325;
    for &b in bases {
        h ^= b as u64;
        h = h.wrapping_mul(0x100000001b3);
    }
    format!("{h:016x}")
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

fn java_seq(stage: &str, hash: &str) -> Option<String> {
    for line in JAVA_TSV.lines() {
        if line.contains(&format!("stage={stage}")) && line.contains(&format!("hash={hash}")) {
            return line.split("\tseq=").nth(1).map(str::to_string);
        }
    }
    None
}

#[test]
fn holdout_6r252_java_actual_haplotype_comparison() {
    if std::env::var("HOLDOUT_6R252").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R252=1");
        return;
    }
    assert_eq!(DEFAULT_NUM_BEST_HAPLOTYPES_PER_GRAPH, 128);
    assert_eq!(
        HcGenotypingConfig::strict_java().stand_emit_confidence,
        DEFAULT_STAND_CALL_CONF
    );
    let cfg = HcLikelihoodEngineConfig::gatk_haplotype_caller_production();
    assert_eq!(cfg.resolved_pair_hmm_backend(), PairHmmBackend::NeonF64);

    let root = repo_root();
    assert!(vcf_has(
        &root.join(JAVA_VCF_REL),
        TARGET,
        TARGET_REF,
        TARGET_ALT
    ));
    assert!(!vcf_has(&root.join(JAVA_VCF_REL), CLOSED_GC, "G", "C"));
    assert!(JAVA_TSV.contains("trim_span\t20:29455560-29455728"));

    for i in 0..5 {
        assert!(JAVA_TSV.contains(JAVA_UNTRIMMED_FNV[i]));
        let jseq = java_seq("trimmed", JAVA_TRIMMED_FNV[i]).expect("java trimmed seq");
        assert_eq!(jseq.len(), 174);
    }

    let ref_fasta = root.join(REF_REL);
    let bam = root.join(BAM_REL);
    let dict = SequenceDictionary::from_fasta_path(&ref_fasta).expect("dict");
    let specs = parse_intervals_cli_string(&dict, INTERVAL).expect("interval");
    begin_hap_list_observe();
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
    let trim = take_hap_list_trim_span().expect("trim");
    assert_eq!(trim.trim_start, 29_455_569);
    assert_eq!(trim.trim_end, 29_455_724);
    let snap = take_colocated_merge_numerics()
        .into_iter()
        .find(|s| s.loc == TARGET)
        .expect("snap");
    assert_eq!(snap.n_reads, 123);
    assert_eq!(snap.pool_sizes[2], 5);
    for (i, &idx) in FROZEN_IDX.iter().enumerate() {
        let h = &outcome.assembly.haplotypes[idx];
        assert_eq!(fnv1a64_hex(&h.bases), RUST_TRIMMED_FNV[i]);
        let jseq = java_seq("trimmed", JAVA_TRIMMED_FNV[i]).unwrap();
        assert_eq!(h.bases, jseq.as_bytes()[9..jseq.len() - 4]);
        assert_ne!(RUST_TRIMMED_FNV[i], JAVA_TRIMMED_FNV[i]);
    }
    let recs =
        try_emit_call_region_variants(region, &outcome, "SAMPLE", DEFAULT_STAND_EMIT_CONFIDENCE)
            .unwrap_or_default();
    assert!(!recs.iter().any(|r| r.position == CLOSED_GC
        && r.reference == "G"
        && r.alternate.iter().any(|a| a == "C")));
    kv("java_trim_span", "29455560-29455728");
    kv("rust_trim_span", "29455569-29455724");
    kv("classification", "HAPLOTYPE_TRIM_DIVERGENCE");
    kv("production_change", "NONE");
}
