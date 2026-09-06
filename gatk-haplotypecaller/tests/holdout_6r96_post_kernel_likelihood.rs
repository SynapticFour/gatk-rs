//! 6R.96 holdout: first post-kernel likelihood-value boundary (live vs Java dump).
//!
//! Skipped unless `HOLDOUT_6R96=1`. Coordinate-free contract lives in
//! `forensic_6r96_post_kernel_likelihood_pipeline_contract`.
//!
//! Java cells: `6r96_java_seq6_post_kernel.tsv` from GATK 4.4.0.0
//! `AlleleLikelihoods` immediately after `computeLog10Likelihoods`, before
//! `normalizeLikelihoods`, seq=6 (`20:29456294-29456500`).
//!
//! ```text
//! HOLDOUT_6R96=1 cargo test -p gatk-haplotypecaller --test holdout_6r96_post_kernel_likelihood -- --nocapture
//! ```

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::{
    begin_likelihood_pipeline_observe, call_disposition, flatten_assembly_regions,
    take_likelihood_pipeline_cells, take_likelihood_pipeline_snaps,
    traverse_assembly_region_walker, try_emit_call_region_variants, AssemblyRegionCallDisposition,
    CallRegionArgs, HaplotypeCallerEngine, ReadFilterParams, WalkerTraversalConfig,
    DEFAULT_STAND_EMIT_CONFIDENCE,
};
use serde_json::json;
use std::collections::{BTreeSet, HashMap, HashSet};
use std::path::{Path, PathBuf};

const INTERVAL: &str = "20:29455000-29456500";
const BAM_REL: &str = "parity/giab/runs/local-pairhmm-diff/HG001.20-29455000-29456500.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const POS_SNP: u64 = 29_456_344;
const JAVA_DUMP: &str = include_str!("6r96_java_seq6_post_kernel.tsv");

const JAVA_ONLY: [u64; 2] = [0xfa2d2442dde7f8ff, 0x48eb4b18de00d4fd];
const JAVA_LIVE_ONLY_N: usize = 10;
const RUST_ONLY_N: usize = 12;

const READS: &[(&str, u16)] = &[
    ("HISEQ1:11:H8GV6ADXX:2:2216:2203:76921", 147),
    ("HISEQ1:13:H8G92ADXX:1:1111:12251:89078", 83),
    ("HISEQ1:9:H8962ADXX:1:1112:19265:60083", 99),
    ("HWI-D00360:5:H814YADXX:1:1202:11051:34179", 147),
    ("HWI-D00360:5:H814YADXX:1:2207:10890:76583", 147),
    ("HWI-D00360:5:H814YADXX:2:1102:2154:52493", 163),
    ("HWI-D00360:6:H81VLADXX:2:1104:15554:2818", 83),
    ("HWI-D00360:6:H81VLADXX:2:1202:18367:85709", 163),
    ("HWI-D00360:7:H88WKADXX:1:1116:9273:30844", 83),
    ("HWI-D00360:8:H88U0ADXX:1:2108:16806:75328", 163),
    ("HISEQ1:13:H8G92ADXX:1:1205:16330:83279", 163),
    ("HISEQ1:9:H8962ADXX:2:1212:17767:73796", 83),
    ("HWI-D00360:5:H814YADXX:2:2103:4936:45407", 83),
    ("HWI-D00360:6:H81VLADXX:1:1103:1948:22968", 147),
    ("HWI-D00360:6:H81VLADXX:1:1210:4156:72506", 83),
    ("HWI-D00360:7:H88WKADXX:1:2111:4466:65743", 147),
    ("HWI-D00360:7:H88WKADXX:1:2203:20480:101193", 163),
    ("HWI-D00360:7:H88WKADXX:2:1214:6938:52704", 83),
    ("HWI-D00360:8:H88U0ADXX:1:1205:11075:4786", 147),
    ("HWI-D00360:8:H88U0ADXX:1:1213:18559:65935", 163),
    ("HWI-D00360:8:H88U0ADXX:1:2213:15618:11579", 163),
    ("HWI-D00360:8:H88U0ADXX:2:1213:15376:17578", 163),
];

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn parse_kv(s: &str) -> HashMap<&str, &str> {
    let mut d = HashMap::new();
    for part in s.split('\t') {
        if let Some((k, v)) = part.split_once('=') {
            d.insert(k, v);
        }
    }
    d
}

struct JavaDump {
    hap_fnv: Vec<u64>,
    bits: HashMap<(String, u16), Vec<u64>>,
}

fn parse_java_dump(text: &str) -> JavaDump {
    let mut hap_fnv = Vec::new();
    let mut bits = HashMap::new();
    for line in text.lines() {
        let p: Vec<&str> = line.split('\t').collect();
        if p.len() < 4 {
            continue;
        }
        let key = p[2];
        let val = p[3..].join("\t");
        let kv = parse_kv(&val);
        if key.starts_with("hap_") {
            let fnv = u64::from_str_radix(kv.get("fnv").copied().unwrap_or("0"), 16).unwrap_or(0);
            hap_fnv.push(fnv);
        } else if key.starts_with("rowbits_") {
            let q = kv.get("qname").copied().unwrap_or("").to_string();
            let flags: u16 = kv.get("flags").copied().unwrap_or("0").parse().unwrap_or(0);
            let b = kv
                .get("bits")
                .copied()
                .unwrap_or("")
                .split(',')
                .filter(|s| !s.is_empty())
                .map(|s| u64::from_str_radix(s, 16).unwrap_or(0))
                .collect();
            bits.insert((q, flags), b);
        }
    }
    JavaDump { hap_fnv, bits }
}

struct Cmp {
    n: usize,
    eq: usize,
    diff: usize,
    max_abs: f64,
    mean_abs: f64,
    first_qname: String,
    first_flags: u16,
    first_hap: u64,
    first_jv: f64,
    first_rv: f64,
}

fn compare_stage(
    java: &HashMap<(String, u16, u64), f64>,
    rust: &HashMap<(String, u16, u64), f64>,
    reads: &[(&str, u16)],
    common: &BTreeSet<u64>,
) -> Cmp {
    let mut n = 0usize;
    let mut eq = 0usize;
    let mut diff = 0usize;
    let mut sum = 0.0;
    let mut max_abs = 0.0;
    let mut first = None;
    for &(q, flags) in reads.iter().take(22) {
        for &h in common {
            n += 1;
            let jk = (q.to_string(), flags, h);
            match (java.get(&jk), rust.get(&jk)) {
                (Some(&jv), Some(&rv)) => {
                    let d = (jv - rv).abs();
                    sum += d;
                    if d > max_abs {
                        max_abs = d;
                    }
                    if jv.to_bits() == rv.to_bits() {
                        eq += 1;
                    } else {
                        diff += 1;
                        if first.is_none() {
                            first = Some((q.to_string(), flags, h, jv, rv));
                        }
                    }
                }
                _ => {
                    diff += 1;
                    if first.is_none() {
                        first = Some((
                            q.to_string(),
                            flags,
                            h,
                            java.get(&jk).copied().unwrap_or(f64::NAN),
                            rust.get(&jk).copied().unwrap_or(f64::NAN),
                        ));
                    }
                }
            }
        }
    }
    let (fq, ff, fh, fj, fr) = first.unwrap_or_else(|| (String::new(), 0, 0, 0.0, 0.0));
    Cmp {
        n,
        eq,
        diff,
        max_abs,
        mean_abs: if n == 0 { 0.0 } else { sum / n as f64 },
        first_qname: fq,
        first_flags: ff,
        first_hap: fh,
        first_jv: fj,
        first_rv: fr,
    }
}

fn rust_stage_map(
    cells: &[gatk_haplotypecaller::LikelihoodPipelineCell],
    stage: &str,
) -> (HashMap<(String, u16, u64), f64>, usize, usize, HashSet<u64>) {
    let seq = cells
        .iter()
        .filter(|c| c.stage == stage)
        .map(|c| c.seq)
        .min()
        .unwrap_or(0);
    let slice: Vec<_> = cells
        .iter()
        .filter(|c| c.stage == stage && c.seq == seq)
        .collect();
    let mut map = HashMap::new();
    let mut reads: HashSet<(String, u16)> = HashSet::new();
    let mut haps: HashSet<u64> = HashSet::new();
    let mut n_reads_dim = 0usize;
    let mut n_haps_dim = 0usize;
    for c in &slice {
        map.insert((c.qname.clone(), c.flags, c.hap_fnv), c.log10_likelihood);
        reads.insert((c.qname.clone(), c.flags));
        haps.insert(c.hap_fnv);
        n_reads_dim = c.n_reads;
        n_haps_dim = c.n_haps;
    }
    let _ = (reads, seq);
    (map, n_reads_dim, n_haps_dim, haps)
}

fn classify(
    post_diff: usize,
    compact_diff: usize,
    norm_diff: usize,
    refresh_n: usize,
) -> &'static str {
    if post_diff > 0 {
        "POST_KERNEL_LIKELIHOOD_OBJECT"
    } else if compact_diff > 0 {
        "LIKELIHOOD_COMPACTION"
    } else if refresh_n > 0 {
        "LIKELIHOOD_REFRESH"
    } else if norm_diff > 0 {
        "NORMALIZE"
    } else {
        "FILTER_INPUT_CONSTRUCTION"
    }
}

#[test]
fn holdout_6r96_first_post_kernel_likelihood_boundary() {
    if std::env::var("HOLDOUT_6R96").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R96=1");
        return;
    }
    let java = parse_java_dump(JAVA_DUMP);
    assert_eq!(java.hap_fnv.len(), 70);
    assert_eq!(java.bits.len(), 153);

    let root = repo_root();
    let ref_fasta = root.join(REF_REL);
    let bam = root.join(BAM_REL);
    assert!(ref_fasta.is_file() && bam.is_file());
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
    let covering: Vec<_> = regions
        .iter()
        .filter(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= POS_SNP
                && r.end.get() >= POS_SNP
        })
        .collect();
    assert_eq!(covering.len(), 1);

    begin_likelihood_pipeline_observe();
    let outcome = HaplotypeCallerEngine::call_region(
        covering[0],
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("call")
    .expect("Some");
    let cells = take_likelihood_pipeline_cells();
    let snaps = take_likelihood_pipeline_snaps();
    let emitted = try_emit_call_region_variants(
        covering[0],
        &outcome,
        "SAMPLE",
        DEFAULT_STAND_EMIT_CONFIDENCE,
    )
    .unwrap_or_default();
    let vcf = emitted
        .iter()
        .find(|r| {
            r.position == POS_SNP && r.reference == "T" && r.alternate.iter().any(|a| a == "C")
        })
        .expect("canonical T/C");

    for s in &snaps {
        eprintln!(
            "6R.96 rust_snap seq={} stage={} n_reads={} n_haps={} n_ll={}",
            s.seq, s.stage, s.n_reads, s.n_haps, s.n_ll_entries
        );
    }

    let (rust_post, rust_post_r, rust_post_c, rust_post_haps) =
        rust_stage_map(&cells, "post_kernel");
    let (rust_comp, rust_comp_r, rust_comp_c, rust_comp_haps) =
        rust_stage_map(&cells, "compaction");
    let (rust_norm, rust_norm_r, rust_norm_c, _rust_norm_haps) =
        rust_stage_map(&cells, "normalize");
    let refresh_present = cells.iter().any(|c| c.stage == "refresh");

    let java_set: HashSet<u64> = java.hap_fnv.iter().copied().collect();
    let rust_only: HashSet<u64> = rust_post_haps.difference(&java_set).copied().collect();
    let java_only: HashSet<u64> = JAVA_ONLY.into_iter().collect();
    assert!(java_only.iter().all(|h| java_set.contains(h)));
    assert!(rust_only.is_empty());
    assert_eq!(rust_post_haps.len(), 70);
    let common_compare: BTreeSet<u64> = java_set
        .intersection(&rust_post_haps)
        .copied()
        .filter(|h| !java_only.contains(h))
        .collect();

    let mut java_post = HashMap::new();
    for (i, &(q, flags)) in READS.iter().enumerate() {
        if i >= 22 {
            break;
        }
        let bits = java
            .bits
            .get(&(q.to_string(), flags))
            .unwrap_or_else(|| panic!("java post-kernel missing {q} flags={flags}"));
        assert_eq!(bits.len(), 70);
        for (col, &fnv) in java.hap_fnv.iter().enumerate() {
            java_post.insert((q.to_string(), flags, fnv), f64::from_bits(bits[col]));
        }
    }

    let post = compare_stage(&java_post, &rust_post, READS, &common_compare);
    let compact_vs_post = compare_stage(&rust_post, &rust_comp, READS, &common_compare);
    let norm_vs_java = compare_stage(&java_post, &rust_norm, READS, &common_compare);

    let classification = classify(
        post.diff,
        if rust_comp.is_empty() {
            0
        } else {
            compact_vs_post.diff
        },
        if rust_norm.is_empty() {
            0
        } else {
            norm_vs_java.diff
        },
        if refresh_present { 1 } else { 0 },
    );

    let live = READS[0];
    let rust_only_rep = READS[JAVA_LIVE_ONLY_N];
    fn stage_cell(
        map: &HashMap<(String, u16, u64), f64>,
        q: &str,
        flags: u16,
        h: u64,
    ) -> Option<f64> {
        map.get(&(q.to_string(), flags, h)).copied()
    }
    let first_hap = *common_compare.iter().min().unwrap_or(&0);
    eprintln!(
        "6R.96 provenance JAVA_LIVE_ONLY {} flags={} hap={:x} java_post={:?} rust_post={:?} rust_comp={:?} rust_norm={:?}",
        live.0,
        live.1,
        first_hap,
        stage_cell(&java_post, live.0, live.1, first_hap),
        stage_cell(&rust_post, live.0, live.1, first_hap),
        stage_cell(&rust_comp, live.0, live.1, first_hap),
        stage_cell(&rust_norm, live.0, live.1, first_hap)
    );
    eprintln!(
        "6R.96 provenance RUST_ONLY {} flags={} hap={:x} java_post={:?} rust_post={:?} rust_comp={:?} rust_norm={:?}",
        rust_only_rep.0,
        rust_only_rep.1,
        first_hap,
        stage_cell(&java_post, rust_only_rep.0, rust_only_rep.1, first_hap),
        stage_cell(&rust_post, rust_only_rep.0, rust_only_rep.1, first_hap),
        stage_cell(&rust_comp, rust_only_rep.0, rust_only_rep.1, first_hap),
        stage_cell(&rust_norm, rust_only_rep.0, rust_only_rep.1, first_hap)
    );

    eprintln!(
        "{}",
        json!({
            "classification": classification,
            "java_post_kernel_rows": java.bits.len(),
            "java_post_kernel_cols": java.hap_fnv.len(),
            "rust_post_kernel_rows": rust_post_r,
            "rust_post_kernel_cols": rust_post_c,
            "rust_compaction_rows": rust_comp_r,
            "rust_compaction_cols": rust_comp_c,
            "rust_normalize_rows": rust_norm_r,
            "rust_normalize_cols": rust_norm_c,
            "common_haps": common_compare.len(),
            "post_kernel_cells": post.n,
            "post_kernel_equal": post.eq,
            "post_kernel_differing": post.diff,
            "post_kernel_max_abs": post.max_abs,
            "post_kernel_mean_abs": post.mean_abs,
            "post_kernel_first_qname": post.first_qname,
            "post_kernel_first_flags": post.first_flags,
            "post_kernel_first_hap": format!("{:x}", post.first_hap),
            "post_kernel_first_java": post.first_jv,
            "post_kernel_first_rust": post.first_rv,
            "compaction_vs_post_differing": compact_vs_post.diff,
            "refresh_present": refresh_present,
            "vcf_ad": vcf.samples.first().map(|s| s.ad.clone()),
        })
    );
    eprintln!(
        "6R.96 classification={classification} post_kernel equal={} differing={} max_abs={} mean_abs={}",
        post.eq, post.diff, post.max_abs, post.mean_abs
    );

    assert_eq!(java.bits.len(), 153);
    assert_eq!(java.hap_fnv.len(), 70);
    assert_eq!(rust_post_r, 153);
    assert_eq!(rust_post_c, 70);
    assert_eq!(common_compare.len(), 68);
    assert_eq!(post.n, 22 * 68);
    assert_eq!(post.eq, 0);
    assert_eq!(post.diff, 1496);
    assert_eq!(compact_vs_post.diff, 0);
    assert_eq!(classification, "POST_KERNEL_LIKELIHOOD_OBJECT");
    assert!(!refresh_present, "non-P12 span must not refresh PairHMM");
    let _ = (
        JAVA_LIVE_ONLY_N,
        RUST_ONLY_N,
        rust_comp_c,
        rust_norm_c,
        rust_comp_haps,
        rust_comp_r,
        rust_norm_r,
        rust_norm_c,
        norm_vs_java.diff,
    );
}
