//! 6R.252: actual GATK 4.4.0.0 Java haplotypes vs Rust at
//! `20:29455649 T/TGTTTG`.
//!
//! Java oracle recovered as official `gatk-package-4.4.0.0-local.jar`
//! (release 4.4.0.0 / source SHA `2dbc025821bc5f686c423ff332a41e6cef892a77`).
//! Dump: `HcFullParityGateDump hap-trim-at-loc` + `eventmap-haps-at-loc`.
//! PRODUCTION CHANGE: NONE.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r252_java_actual_haplotype_comparison -- --nocapture --test-threads=1
//! HOLDOUT_6R252=1 cargo test -p gatk-haplotypecaller --test holdout_6r252_java_actual_haplotype_comparison -- --nocapture --test-threads=1
//! ```

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::event_map::{build_per_haplotype_variation_events, overlapping_events};
use gatk_haplotypecaller::hc_genotyping_engine::take_colocated_merge_numerics;
use gatk_haplotypecaller::read_threading_assembler::DEFAULT_NUM_BEST_HAPLOTYPES_PER_GRAPH;
use gatk_haplotypecaller::{
    begin_hap_list_observe, begin_trim_io_observe, call_disposition, flatten_assembly_regions,
    take_hap_list_trim_span, take_trim_io, traverse_assembly_region_walker,
    AssemblyRegionCallDisposition, CallRegionArgs, HaplotypeCallerEngine, ReadFilterParams,
    WalkerTraversalConfig,
};
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

const JAVA_PIN: &str = "4.4.0.0 / 2dbc025821bc5f686c423ff332a41e6cef892a77";
const INTERVAL: &str = "20:29455000-29456500";
const BAM_REL: &str = "parity/giab/runs/local-pairhmm-diff/HG001.20-29455000-29456500.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const TARGET: u64 = 29_455_649;
const TARGET_REF: &str = "T";
const TARGET_ALT: &str = "TGTTTG";
const FROZEN_IDX: [usize; 5] = [19, 20, 21, 22, 23];
const RUST_TRIMMED_FNV: [&str; 5] = [
    "55012fcf3b430591",
    "341b2e070ccb5846",
    "6a65e4c02733c2ed",
    "fa07750bf228b3c2",
    "79451c576721a729",
];
const RUST_UNTRIMMED_FNV: [&str; 5] = [
    "a1ed2811a3912892",
    "745492ecfc910685",
    "dae45a854e2d04de",
    "8bdf55500c4f07d1",
    "5cbde2ec4e0a98ea",
];
const JAVA_TRIMMED_FNV: [&str; 5] = [
    "fcd72c6ce600dd16",
    "c1a4e8204522f645",
    "eb03271fa7548f26",
    "7dc2a8ae5da116e1",
    "343e7c443c1d9732",
];
const JAVA_TSV: &str = include_str!("forensic_6r252_java_hap_trim.tsv");
const JAVA_EM_TSV: &str = include_str!("forensic_6r252_java_eventmap_haps.tsv");
const JAVA_TRIM_PREFIX: usize = 9;
const JAVA_TRIM_SUFFIX: usize = 4;

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    eprintln!("6R252\t{key}\t{}", value.as_ref());
}

fn fnv1a64_hex(bases: &[u8]) -> String {
    let mut h: u64 = 0xcbf29ce484222325;
    for &b in bases {
        h ^= b as u64;
        h = h.wrapping_mul(0x100000001b3);
    }
    format!("{h:016x}")
}

fn fields(line: &str) -> BTreeMap<String, String> {
    let mut out = BTreeMap::new();
    for p in line.split('\t').skip(2) {
        if let Some((k, v)) = p.split_once('=') {
            out.insert(k.to_string(), v.to_string());
        } else if !p.contains('=') && p.contains(':') {
            out.insert("span".to_string(), p.to_string());
        }
    }
    out
}

#[derive(Clone)]
struct JavaHap {
    stage: String,
    idx: usize,
    hash: String,
    len: usize,
    cigar: String,
    seq: String,
}

fn parse_java_haps(stage: &str) -> Vec<JavaHap> {
    let mut out = Vec::new();
    for line in JAVA_TSV.lines() {
        if !line.starts_with("6R131\thap\t") {
            continue;
        }
        let f = fields(line);
        if f.get("stage").map(String::as_str) != Some(stage) {
            continue;
        }
        out.push(JavaHap {
            stage: stage.to_string(),
            idx: f["idx"].parse().unwrap(),
            hash: f["hash"].clone(),
            len: f["len"].parse().unwrap(),
            cigar: f["cigar"].clone(),
            seq: f.get("seq").cloned().unwrap_or_default(),
        });
    }
    out
}

#[test]
fn forensic_6r252_java_oracle_recovery() {
    kv("java_pin", JAVA_PIN);
    kv("production_change", "NONE");
    kv(
        "executable",
        "parity/build/gatk-oracle/gatk-package-4.4.0.0-local.jar (official GitHub release 4.4.0.0 zip)",
    );
    kv(
        "execution_method",
        "java -cp $GATK_JAR:$class_dir HcFullParityGateDump hap-trim-at-loc REF BAM 20:29455000-29456500 29455649",
    );
    kv("bam", BAM_REL);
    kv("reference", REF_REL);
    kv(
        "gatk_version_observed",
        "The Genome Analysis Toolkit (GATK) v4.4.0.0",
    );
    kv(
        "docker_at_recovery",
        "CLI hung on docker.sock; recovered official 4.4.0.0 JAR instead of pulling an image",
    );
    assert!(JAVA_TSV.contains("6R131\tregion_active\t20:29455560-29455744"));
    assert!(JAVA_TSV.contains("6R131\ttrim_span\t20:29455560-29455728"));
    assert!(JAVA_EM_TSV.contains("ref=T\talt=TGTTTG"));
}

#[test]
fn forensic_6r252_java_actual_haplotype_comparison() {
    kv("java_pin", JAVA_PIN);
    kv("production_change", "NONE");
    kv("target", "20:29455649 T/TGTTTG");
    assert_eq!(DEFAULT_NUM_BEST_HAPLOTYPES_PER_GRAPH, 128);

    let java_untrimmed = parse_java_haps("untrimmed");
    let java_trimmed = parse_java_haps("trimmed");
    kv("java_untrimmed_n", java_untrimmed.len().to_string());
    kv("java_trimmed_n", java_trimmed.len().to_string());
    assert_eq!(java_untrimmed.len(), 128);
    assert_eq!(java_trimmed.len(), 24);

    let mut java_tgtttg_trimmed = Vec::new();
    for line in JAVA_EM_TSV.lines() {
        if line.contains("stage=trimmed")
            && line.contains("alt=TGTTTG")
            && line.contains("start=29455649")
        {
            let f = fields(line);
            java_tgtttg_trimmed.push(f["hash"].clone());
        }
    }
    java_tgtttg_trimmed.sort();
    java_tgtttg_trimmed.dedup();
    kv(
        "java_tgtttg_trimmed_n",
        java_tgtttg_trimmed.len().to_string(),
    );
    assert_eq!(java_tgtttg_trimmed.len(), 5);

    let root = repo_root();
    let ref_fasta = root.join(REF_REL);
    let bam = root.join(BAM_REL);
    if !ref_fasta.is_file() || !bam.is_file() {
        eprintln!("skip: missing chr20_tiny BAM/REF");
        return;
    }
    let dict = SequenceDictionary::from_fasta_path(&ref_fasta).expect("dict");
    let specs = parse_intervals_cli_string(&dict, INTERVAL).expect("interval");
    begin_hap_list_observe();
    begin_trim_io_observe();
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
    let trim = take_hap_list_trim_span().expect("trim span");
    let io = take_trim_io();
    let snap = take_colocated_merge_numerics()
        .into_iter()
        .find(|s| s.loc == TARGET)
        .expect("snap");
    assert_eq!(snap.pool_sizes[2], 5);
    kv(
        "rust_trim_span",
        format!("{}-{}", trim.trim_start, trim.trim_end),
    );
    kv("java_trim_span", "29455560-29455728");
    kv(
        "rust_active",
        format!("{}-{}", region.start.get(), region.end.get()),
    );
    assert_eq!(region.start.get(), 29_455_560);
    assert_eq!(region.end.get(), 29_455_744);
    assert_eq!(trim.trim_start, 29_455_569);
    assert_eq!(trim.trim_end, 29_455_724);

    let haps = &outcome.assembly.haplotypes;
    let hap_events = build_per_haplotype_variation_events(
        haps,
        outcome.assembly.reference_bases(),
        outcome.assembly.padded_reference_start_1based(),
        outcome.assembly.max_mnp_distance(),
        "20",
    );

    let mut n_untrimmed_match = 0usize;
    let mut n_trimmed_fnv_match = 0usize;
    let mut n_interior = 0usize;
    for (i, &idx) in FROZEN_IDX.iter().enumerate() {
        let h = &haps[idx];
        let rust_fnv = fnv1a64_hex(&h.bases);
        assert_eq!(rust_fnv, RUST_TRIMMED_FNV[i]);
        assert_eq!(h.bases.len(), 161);
        let spanning = overlapping_events(hap_events.events_for(idx), TARGET);
        assert!(spanning
            .iter()
            .any(|e| e.ref_allele == TARGET_REF && e.alt_allele == TARGET_ALT));

        let ju = java_untrimmed
            .iter()
            .find(|j| j.hash == RUST_UNTRIMMED_FNV[i])
            .expect("java untrimmed parent");
        n_untrimmed_match += 1;
        kv(
            &format!("H{i}_untrimmed"),
            format!(
                "fnv={} java_idx={} rust_parent_fnv={} len={} MATCH",
                ju.hash, ju.idx, RUST_UNTRIMMED_FNV[i], ju.len
            ),
        );

        let jt = java_trimmed
            .iter()
            .find(|j| j.hash == JAVA_TRIMMED_FNV[i])
            .expect("java trimmed");
        assert_eq!(jt.len, 174);
        assert_eq!(jt.cigar, "90M5I79M");
        assert_eq!(jt.idx, idx);
        if rust_fnv == jt.hash {
            n_trimmed_fnv_match += 1;
        }
        let interior = &jt.seq.as_bytes()[JAVA_TRIM_PREFIX..jt.seq.len() - JAVA_TRIM_SUFFIX];
        let interior_eq = interior == h.bases.as_slice();
        if interior_eq {
            n_interior += 1;
        }
        kv(
            &format!("H{i}_trimmed"),
            format!(
                "rust_fnv={rust_fnv} rust_len=161 rust_cigar=81M5I75M java_fnv={} java_len=174 java_cigar=90M5I79M interior_eq={interior_eq}",
                jt.hash
            ),
        );
        assert!(
            interior_eq,
            "Rust 161-mer must be Java 174-mer without 9 bp prefix and 4 bp suffix"
        );
        assert_ne!(rust_fnv, jt.hash, "trimmed FNVs must differ");

        let parent_kept = io.iter().any(|r| {
            format!("{:016x}", r.input_fnv) == RUST_UNTRIMMED_FNV[i] && r.outcome == "kept"
        });
        assert!(
            parent_kept,
            "Rust trim_io kept parent {}",
            RUST_UNTRIMMED_FNV[i]
        );
    }
    kv("n_untrimmed_fnv_match", n_untrimmed_match.to_string());
    kv("n_trimmed_fnv_match", n_trimmed_fnv_match.to_string());
    kv("n_interior_subsequence", n_interior.to_string());
    assert_eq!(n_untrimmed_match, 5);
    assert_eq!(n_trimmed_fnv_match, 0);
    assert_eq!(n_interior, 5);

    kv(
        "first_divergent_operation",
        "AssemblyRegionTrimmer.trim / trimTo span: Java 20:29455560-29455728 (174 bp, 90M5I79M) vs Rust 20:29455569-29455724 (161 bp, 81M5I75M). Untrimmed Path.getBases() identical.",
    );
    kv("classification", "HAPLOTYPE_TRIM_DIVERGENCE");
    kv(
        "haplotype_input_divergence",
        "CLOSED as sequence-identity at assembleReads; remaining split is the genotyping trim window",
    );
}
