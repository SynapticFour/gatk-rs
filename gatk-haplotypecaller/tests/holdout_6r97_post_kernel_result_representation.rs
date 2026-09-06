//! 6R.97 holdout: primitive kernel-result representation vs likelihood object.
//!
//! Skipped unless `HOLDOUT_6R97=1`. Coordinate-free contract lives in
//! `forensic_6r97_post_kernel_result_representation_contract`.
//!
//! Java cells: `6r96_java_seq6_post_kernel.tsv` from GATK 4.4.0.0
//! `AlleleLikelihoods` immediately after `computeLog10Likelihoods` (`SampleMatrix.set`
//! of GKL `double[] mLogLikelihoodArray`), before `normalizeLikelihoods`, seq=6
//! (`20:29456294-29456500`).
//!
//! Rust: `LikelihoodPipelineCell` stage `post_kernel` immediately after
//! `score_pairhmm_from_records` assignment into `RegionReadLikelihood`.
//!
//! ```text
//! HOLDOUT_6R97=1 cargo test -p gatk-haplotypecaller --test holdout_6r97_post_kernel_result_representation -- --nocapture
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
const JAVA_BUFFER: &str = include_str!("6r97_java_seq6_kernel_buffer.tsv");

const JAVA_ONLY: [u64; 2] = [0xfa2d2442dde7f8ff, 0x48eb4b18de00d4fd];

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

fn is_exact_f32_widened(x: f64) -> bool {
    (f64::from(x as f32)).to_bits() == x.to_bits()
}

fn ordered_f64_bits(x: f64) -> i64 {
    let bits = x.to_bits() as i64;
    if bits < 0 {
        i64::MIN - bits
    } else {
        bits
    }
}

fn ulp_distance(a: f64, b: f64) -> u64 {
    ordered_f64_bits(a).abs_diff(ordered_f64_bits(b))
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

fn parse_buffer_dump(text: &str) -> HashMap<String, String> {
    let mut d = HashMap::new();
    for line in text.lines() {
        let p: Vec<&str> = line.split('\t').collect();
        if p.len() < 4 || p[0] != "6R97" {
            continue;
        }
        d.insert(p[2].to_string(), p[3].to_string());
    }
    d
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

struct WidthStats {
    java_f32_wide: usize,
    rust_f32_wide: usize,
    rust_f32_eq_java: usize,
    java_f32_eq_rust_f32: usize,
    java_f32_eq_rust: usize,
    both_f64_f32_f64_eq: usize,
    max_residual_after_rust_f32: f64,
    mean_residual_after_rust_f32: f64,
    max_ulp: u64,
    mean_ulp: f64,
    sign_mismatch: usize,
    java_neg_zero: usize,
    rust_neg_zero: usize,
    java_nan: usize,
    rust_nan: usize,
    java_inf: usize,
    rust_inf: usize,
}

fn width_stats(
    java: &HashMap<(String, u16, u64), f64>,
    rust: &HashMap<(String, u16, u64), f64>,
    reads: &[(&str, u16)],
    common: &BTreeSet<u64>,
) -> WidthStats {
    let mut s = WidthStats {
        java_f32_wide: 0,
        rust_f32_wide: 0,
        rust_f32_eq_java: 0,
        java_f32_eq_rust_f32: 0,
        java_f32_eq_rust: 0,
        both_f64_f32_f64_eq: 0,
        max_residual_after_rust_f32: 0.0,
        mean_residual_after_rust_f32: 0.0,
        max_ulp: 0,
        mean_ulp: 0.0,
        sign_mismatch: 0,
        java_neg_zero: 0,
        rust_neg_zero: 0,
        java_nan: 0,
        rust_nan: 0,
        java_inf: 0,
        rust_inf: 0,
    };
    let mut n = 0usize;
    let mut sum_res = 0.0;
    let mut sum_ulp = 0.0;
    for &(q, flags) in reads.iter().take(22) {
        for &h in common {
            let jk = (q.to_string(), flags, h);
            let Some(&jv) = java.get(&jk) else {
                continue;
            };
            let Some(&rv) = rust.get(&jk) else {
                continue;
            };
            n += 1;
            if is_exact_f32_widened(jv) {
                s.java_f32_wide += 1;
            }
            if is_exact_f32_widened(rv) {
                s.rust_f32_wide += 1;
            }
            let j_f32 = f64::from(jv as f32);
            let r_f32 = f64::from(rv as f32);
            if r_f32.to_bits() == jv.to_bits() {
                s.rust_f32_eq_java += 1;
            }
            if j_f32.to_bits() == r_f32.to_bits() {
                s.java_f32_eq_rust_f32 += 1;
            }
            if j_f32.to_bits() == rv.to_bits() {
                s.java_f32_eq_rust += 1;
            }
            let j_rt = f64::from((jv as f32) as f32);
            let r_rt = f64::from((rv as f32) as f32);
            if j_rt.to_bits() == r_rt.to_bits() {
                s.both_f64_f32_f64_eq += 1;
            }
            let res = (jv - r_f32).abs();
            sum_res += res;
            if res > s.max_residual_after_rust_f32 {
                s.max_residual_after_rust_f32 = res;
            }
            let ulp = ulp_distance(jv, rv);
            sum_ulp += ulp as f64;
            if ulp > s.max_ulp {
                s.max_ulp = ulp;
            }
            if jv.is_sign_negative() != rv.is_sign_negative() {
                s.sign_mismatch += 1;
            }
            if jv.to_bits() == (-0.0f64).to_bits() {
                s.java_neg_zero += 1;
            }
            if rv.to_bits() == (-0.0f64).to_bits() {
                s.rust_neg_zero += 1;
            }
            if jv.is_nan() {
                s.java_nan += 1;
            }
            if rv.is_nan() {
                s.rust_nan += 1;
            }
            if jv.is_infinite() {
                s.java_inf += 1;
            }
            if rv.is_infinite() {
                s.rust_inf += 1;
            }
        }
    }
    s.mean_residual_after_rust_f32 = if n == 0 { 0.0 } else { sum_res / n as f64 };
    s.mean_ulp = if n == 0 { 0.0 } else { sum_ulp / n as f64 };
    s
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
    let mut haps: HashSet<u64> = HashSet::new();
    let mut n_reads_dim = 0usize;
    let mut n_haps_dim = 0usize;
    for c in &slice {
        map.insert((c.qname.clone(), c.flags, c.hap_fnv), c.log10_likelihood);
        haps.insert(c.hap_fnv);
        n_reads_dim = c.n_reads;
        n_haps_dim = c.n_haps;
    }
    (map, n_reads_dim, n_haps_dim, haps)
}

fn classify(
    primitive_eq: usize,
    n: usize,
    java_f32_wide: usize,
    rust_f32_eq_java: usize,
    buffer_eq_matrix: bool,
) -> &'static str {
    if primitive_eq == n {
        if buffer_eq_matrix {
            "NO_PROVEN_DIVERGENCE"
        } else {
            "LIKELIHOOD_OBJECT_MATERIALIZATION"
        }
    } else if java_f32_wide == n && rust_f32_eq_java == n {
        "REPRESENTATION_WIDTH"
    } else if !buffer_eq_matrix {
        "RESULT_BUFFER_SEMANTICS"
    } else {
        "KERNEL_OUTPUT_TRANSFER_SOURCE"
    }
}

#[test]
fn holdout_6r97_post_kernel_result_representation() {
    if std::env::var("HOLDOUT_6R97").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R97=1");
        return;
    }
    let java = parse_java_dump(JAVA_DUMP);
    assert_eq!(java.hap_fnv.len(), 70);
    assert_eq!(java.bits.len(), 153);
    let buf = parse_buffer_dump(JAVA_BUFFER);
    assert_eq!(
        buf.get("kernel_buffer_type").map(String::as_str),
        Some("double[] mLogLikelihoodArray")
    );
    assert_eq!(
        buf.get("kernel_buffer_n").map(String::as_str),
        Some("10710")
    );
    assert_eq!(
        buf.get("kernel_buffer_f32_wide").map(String::as_str),
        Some("10710")
    );
    assert_eq!(buf.get("matrix_n").map(String::as_str), Some("10710"));
    assert_eq!(
        buf.get("matrix_f32_wide").map(String::as_str),
        Some("10710")
    );
    assert_eq!(
        buf.get("buffer_matrix_sorted_bits_equal")
            .map(String::as_str),
        Some("true")
    );
    let buffer_eq_matrix = true;

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
            "6R.97 rust_snap seq={} stage={} n_reads={} n_haps={} n_ll={}",
            s.seq, s.stage, s.n_reads, s.n_haps, s.n_ll_entries
        );
    }

    let (rust_post, rust_post_r, rust_post_c, rust_post_haps) =
        rust_stage_map(&cells, "post_kernel");
    let java_set: HashSet<u64> = java.hap_fnv.iter().copied().collect();
    let java_only: HashSet<u64> = JAVA_ONLY.into_iter().collect();
    let common_compare: BTreeSet<u64> = java_set
        .intersection(&rust_post_haps)
        .copied()
        .filter(|h| !java_only.contains(h))
        .collect();

    let mut java_post = HashMap::new();
    for &(q, flags) in READS.iter().take(22) {
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
    let w = width_stats(&java_post, &rust_post, READS, &common_compare);
    let classification = classify(
        post.eq,
        post.n,
        w.java_f32_wide,
        w.rust_f32_eq_java,
        buffer_eq_matrix,
    );

    eprintln!(
        "{}",
        json!({
            "classification": classification,
            "java_buffer_type": "double[] mLogLikelihoodArray",
            "java_buffer_n": 10710,
            "java_buffer_f32_wide": 10710,
            "buffer_matrix_sorted_bits_equal": true,
            "rust_post_kernel_rows": rust_post_r,
            "rust_post_kernel_cols": rust_post_c,
            "common_haps": common_compare.len(),
            "cells": post.n,
            "exact_equal_bits": post.eq,
            "differing": post.diff,
            "max_abs": post.max_abs,
            "mean_abs": post.mean_abs,
            "first_qname": post.first_qname,
            "first_flags": post.first_flags,
            "first_hap": format!("{:x}", post.first_hap),
            "first_java": post.first_jv,
            "first_rust": post.first_rv,
            "first_java_bits": format!("{:016x}", post.first_jv.to_bits()),
            "first_rust_bits": format!("{:016x}", post.first_rv.to_bits()),
            "java_f32_wide": w.java_f32_wide,
            "rust_f32_wide": w.rust_f32_wide,
            "rust_f32_eq_java": w.rust_f32_eq_java,
            "java_f32_eq_rust_f32": w.java_f32_eq_rust_f32,
            "java_f32_eq_rust": w.java_f32_eq_rust,
            "both_f64_f32_f64_eq": w.both_f64_f32_f64_eq,
            "max_residual_after_rust_f32": w.max_residual_after_rust_f32,
            "mean_residual_after_rust_f32": w.mean_residual_after_rust_f32,
            "max_ulp": w.max_ulp,
            "mean_ulp": w.mean_ulp,
            "sign_mismatch": w.sign_mismatch,
            "java_neg_zero": w.java_neg_zero,
            "rust_neg_zero": w.rust_neg_zero,
            "java_nan": w.java_nan,
            "rust_nan": w.rust_nan,
            "java_inf": w.java_inf,
            "rust_inf": w.rust_inf,
            "vcf_ad": vcf.samples.first().map(|s| s.ad.clone()),
        })
    );
    eprintln!(
        "6R.97 classification={classification} equal={} differing={} java_f32_wide={} rust_f32_wide={} rust_f32_eq_java={} max_res={}",
        post.eq, post.diff, w.java_f32_wide, w.rust_f32_wide, w.rust_f32_eq_java, w.max_residual_after_rust_f32
    );

    assert_eq!(java.bits.len(), 153);
    assert_eq!(java.hap_fnv.len(), 70);
    assert_eq!(rust_post_r, 153);
    assert_eq!(rust_post_c, 70);
    assert_eq!(common_compare.len(), 68);
    assert_eq!(post.n, 22 * 68);
    assert_eq!(post.eq, 0);
    assert_eq!(post.diff, 1496);
    assert_eq!(w.java_f32_wide, 1496);
    assert_eq!(w.rust_f32_wide, 0);
    assert_eq!(w.rust_f32_eq_java, 890);
    assert_eq!(w.java_f32_eq_rust_f32, 890);
    assert_eq!(w.java_f32_eq_rust, 0);
    assert_eq!(w.sign_mismatch, 0);
    assert_eq!(w.java_nan, 0);
    assert_eq!(w.rust_nan, 0);
    assert_eq!(w.java_inf, 0);
    assert_eq!(w.rust_inf, 0);
    assert_eq!(w.java_neg_zero, 0);
    assert_eq!(w.rust_neg_zero, 0);
    assert_eq!(classification, "KERNEL_OUTPUT_TRANSFER_SOURCE");
    let _ = outcome;
}
