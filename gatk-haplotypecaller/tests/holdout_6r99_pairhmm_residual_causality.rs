//! 6R.99 holdout: does the 6R.98 primitive residual move any of 22 reads across −8.0?
//!
//! Skipped unless `HOLDOUT_6R99=1`. Coordinate-free contract lives in
//! `forensic_6r99_pairhmm_residual_causality_contract`.
//!
//! ```text
//! HOLDOUT_6R99=1 cargo test -p gatk-haplotypecaller --test holdout_6r99_pairhmm_residual_causality -- --nocapture
//! ```

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::{
    begin_likelihood_pipeline_observe, begin_poorly_modeled_observe, call_disposition,
    flatten_assembly_regions, take_likelihood_pipeline_cells, take_likelihood_pipeline_snaps,
    take_poorly_modeled_observe, traverse_assembly_region_walker, try_emit_call_region_variants,
    AssemblyRegionCallDisposition, CallRegionArgs, HaplotypeCallerEngine, ReadFilterParams,
    WalkerTraversalConfig, DEFAULT_STAND_EMIT_CONFIDENCE,
};
use serde_json::json;
use std::collections::{BTreeSet, HashMap, HashSet};
use std::path::{Path, PathBuf};

const INTERVAL: &str = "20:29455000-29456500";
const BAM_REL: &str = "parity/giab/runs/local-pairhmm-diff/HG001.20-29455000-29456500.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const POS_SNP: u64 = 29_456_344;
const JAVA_FLOAT: &str = include_str!("6r96_java_seq6_post_kernel.tsv");
const JAVA_DOUBLE: &str = include_str!("6r98_java_seq6_double_post_kernel.tsv");
const JAVA_ONLY: [u64; 2] = [0xfa2d2442dde7f8ff, 0x48eb4b18de00d4fd];
const THRESHOLD: f64 = -8.0;

/// 10 JAVA_LIVE_ONLY (Java KEEP), 12 RUST_ONLY (Java DROP), 2 BOTH_DROP controls.
const READS: &[(&str, u16, bool, usize)] = &[
    ("HISEQ1:11:H8GV6ADXX:2:2216:2203:76921", 147, true, 130),
    ("HISEQ1:13:H8G92ADXX:1:1111:12251:89078", 83, true, 148),
    ("HISEQ1:9:H8962ADXX:1:1112:19265:60083", 99, true, 76),
    ("HWI-D00360:5:H814YADXX:1:1202:11051:34179", 147, true, 96),
    ("HWI-D00360:5:H814YADXX:1:2207:10890:76583", 147, true, 148),
    ("HWI-D00360:5:H814YADXX:2:1102:2154:52493", 163, true, 148),
    ("HWI-D00360:6:H81VLADXX:2:1104:15554:2818", 83, true, 143),
    ("HWI-D00360:6:H81VLADXX:2:1202:18367:85709", 163, true, 134),
    ("HWI-D00360:7:H88WKADXX:1:1116:9273:30844", 83, true, 132),
    ("HWI-D00360:8:H88U0ADXX:1:2108:16806:75328", 163, true, 112),
    ("HISEQ1:13:H8G92ADXX:1:1205:16330:83279", 163, false, 148),
    ("HISEQ1:9:H8962ADXX:2:1212:17767:73796", 83, false, 148),
    ("HWI-D00360:5:H814YADXX:2:2103:4936:45407", 83, false, 84),
    ("HWI-D00360:6:H81VLADXX:1:1103:1948:22968", 147, false, 148),
    ("HWI-D00360:6:H81VLADXX:1:1210:4156:72506", 83, false, 92),
    ("HWI-D00360:7:H88WKADXX:1:2111:4466:65743", 147, false, 144),
    (
        "HWI-D00360:7:H88WKADXX:1:2203:20480:101193",
        163,
        false,
        148,
    ),
    ("HWI-D00360:7:H88WKADXX:2:1214:6938:52704", 83, false, 125),
    ("HWI-D00360:8:H88U0ADXX:1:1205:11075:4786", 147, false, 148),
    ("HWI-D00360:8:H88U0ADXX:1:1213:18559:65935", 163, false, 148),
    ("HWI-D00360:8:H88U0ADXX:1:2213:15618:11579", 163, false, 105),
    ("HWI-D00360:8:H88U0ADXX:2:1213:15376:17578", 163, false, 104),
    ("HISEQ1:11:H8GV6ADXX:2:2105:12137:22761", 163, false, 110),
    ("HISEQ1:13:H8G92ADXX:1:1201:11859:45984", 147, false, 97),
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

fn java_keep(max_ll: f64) -> bool {
    !(max_ll < THRESHOLD)
}

struct JavaDump {
    hap_fnv: Vec<u64>,
    bits: HashMap<(String, u16), Vec<u64>>,
}

fn parse_java_dump(text: &str, prefix: &str) -> JavaDump {
    let mut hap_fnv = Vec::new();
    let mut bits = HashMap::new();
    for line in text.lines() {
        let p: Vec<&str> = line.split('\t').collect();
        if p.len() < 4 || p[0] != prefix {
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

fn cells_from_dump(dump: &JavaDump) -> HashMap<(String, u16, u64), f64> {
    let mut map = HashMap::new();
    for ((q, flags), bits) in &dump.bits {
        for (col, &fnv) in dump.hap_fnv.iter().enumerate() {
            map.insert((q.clone(), *flags, fnv), f64::from_bits(bits[col]));
        }
    }
    map
}

fn rust_post(
    cells: &[gatk_haplotypecaller::LikelihoodPipelineCell],
) -> HashMap<(String, u16, u64), f64> {
    let seq = cells
        .iter()
        .filter(|c| c.stage == "post_kernel")
        .map(|c| c.seq)
        .min()
        .unwrap_or(0);
    let mut map = HashMap::new();
    for c in cells
        .iter()
        .filter(|c| c.stage == "post_kernel" && c.seq == seq)
    {
        map.insert((c.qname.clone(), c.flags, c.hap_fnv), c.log10_likelihood);
    }
    map
}

struct Ranked {
    fnv: u64,
    val: f64,
}

fn ranked(
    row: &HashMap<(String, u16, u64), f64>,
    q: &str,
    flags: u16,
    common: &BTreeSet<u64>,
    hap_order: &[u64],
) -> Vec<Ranked> {
    let mut out = Vec::new();
    for &h in hap_order {
        if !common.contains(&h) {
            continue;
        }
        if let Some(&v) = row.get(&(q.to_string(), flags, h)) {
            out.push(Ranked { fnv: h, val: v });
        }
    }
    out.sort_by(|a, b| {
        b.val
            .partial_cmp(&a.val)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    out
}

/// Java `maximumLikelihoodOverAllAlleles`: first haplotype in list order with strict `>`.
fn max_and_winner(
    row: &HashMap<(String, u16, u64), f64>,
    q: &str,
    flags: u16,
    common: &BTreeSet<u64>,
    hap_order: &[u64],
) -> (f64, u64, f64) {
    let mut best = f64::NEG_INFINITY;
    let mut win = 0u64;
    let mut second = f64::NEG_INFINITY;
    for &h in hap_order {
        if !common.contains(&h) {
            continue;
        }
        let Some(&v) = row.get(&(q.to_string(), flags, h)) else {
            continue;
        };
        if v > best {
            second = best;
            best = v;
            win = h;
        } else if v > second {
            second = v;
        }
    }
    (best, win, best - second)
}

fn classify(
    n_residual_cross: usize,
    n_keep_drop_change: usize,
    n_winner_changes_keep_drop: usize,
) -> &'static str {
    if n_winner_changes_keep_drop > 0 {
        "PAIRHMM_WINNER_SWITCH"
    } else if n_residual_cross > 0 || n_keep_drop_change > 0 {
        "PAIRHMM_RESIDUAL_POTENTIALLY_CAUSAL"
    } else {
        "PAIRHMM_RESIDUAL_NOT_CAUSAL"
    }
}

#[test]
fn holdout_6r99_pairhmm_residual_causality() {
    if std::env::var("HOLDOUT_6R99").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R99=1");
        return;
    }
    let java_f = parse_java_dump(JAVA_FLOAT, "6R96");
    let java_d = parse_java_dump(JAVA_DOUBLE, "6R98");
    assert_eq!(java_f.hap_fnv, java_d.hap_fnv);
    assert_eq!(java_f.hap_fnv.len(), 70);

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
    begin_poorly_modeled_observe();
    let outcome = HaplotypeCallerEngine::call_region(
        covering[0],
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("call")
    .expect("Some");
    let cells = take_likelihood_pipeline_cells();
    let _snaps = take_likelihood_pipeline_snaps();
    let filter_rows = take_poorly_modeled_observe();
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

    let last_pass = filter_rows.iter().map(|r| r.pass).max().unwrap_or(0);
    let rust_filter: HashMap<(String, u16), f64> = filter_rows
        .iter()
        .filter(|r| r.pass == last_pass)
        .map(|r| ((r.qname.clone(), r.flags), r.max_ll))
        .collect();

    let rust = rust_post(&cells);
    let java_only: HashSet<u64> = JAVA_ONLY.into_iter().collect();
    let common: BTreeSet<u64> = java_f
        .hap_fnv
        .iter()
        .copied()
        .filter(|h| !java_only.contains(h))
        .collect();
    assert_eq!(common.len(), 68);

    let jf = cells_from_dump(&java_f);
    let jd = cells_from_dump(&java_d);

    let mut max_cell_abs = 0.0f64;
    let mut n_cell = 0usize;
    let mut sum_abs = 0.0f64;
    for &(q, flags, _, _) in READS.iter().take(22) {
        for &h in &common {
            let k = (q.to_string(), flags, h);
            let Some(&a) = jd.get(&k) else {
                continue;
            };
            let Some(&b) = rust.get(&k) else {
                continue;
            };
            n_cell += 1;
            let d = (a - b).abs();
            sum_abs += d;
            max_cell_abs = max_cell_abs.max(d);
        }
    }
    assert_eq!(n_cell, 1496);

    let mut n_cross = 0usize;
    let mut n_class_change = 0usize;
    let mut n_winner_diff = 0usize;
    let mut n_winner_changes_keep_drop = 0usize;
    let mut n_safe = 0usize;
    let mut min_req = f64::INFINITY;
    let mut rows = Vec::new();

    eprintln!(
        "QNAME\tJava_dec\tRust_post_dec\tJ_win\tR_win\tsame_win\tJ_max68\tR_max68\tdist\tcross\tfilter_max"
    );

    let all70: BTreeSet<u64> = java_f.hap_fnv.iter().copied().collect();
    for (i, &(q, flags, jkeep, qlen)) in READS.iter().enumerate() {
        let jr = ranked(&jd, q, flags, &common, &java_d.hap_fnv);
        let rr = ranked(&rust, q, flags, &common, &java_d.hap_fnv);
        let (jmax, jwin, jmargin) = max_and_winner(&jd, q, flags, &common, &java_d.hap_fnv);
        let (rmax, rwin, rmargin) = max_and_winner(&rust, q, flags, &common, &java_d.hap_fnv);
        let (jmax70, jwin70, _) = max_and_winner(&jf, q, flags, &all70, &java_f.hap_fnv);
        let jdec = java_keep(jmax);
        let rdec = java_keep(rmax);
        assert_eq!(jdec, jkeep, "{q} dump max68 keep matches 6R.93");
        let dist = (jmax - THRESHOLD).abs().min((rmax - THRESHOLD).abs());
        let can_cross = dist <= max_cell_abs;
        if can_cross {
            n_cross += 1;
        }
        if jdec != rdec {
            n_class_change += 1;
        }
        if jwin != rwin {
            n_winner_diff += 1;
            if jdec != rdec {
                n_winner_changes_keep_drop += 1;
            }
        }
        if dist > max_cell_abs * 10.0 {
            n_safe += 1;
        }
        min_req = min_req.min(dist);
        let filter_max = rust_filter.get(&(q.to_string(), flags)).copied();
        eprintln!(
            "{}\t{}\t{}\t{:x}\t{:x}\t{}\t{}\t{}\t{}\t{}\t{:?}",
            q,
            if jdec { "KEEP" } else { "DROP" },
            if rdec { "KEEP" } else { "DROP" },
            jwin,
            rwin,
            jwin == rwin,
            jmax,
            rmax,
            dist,
            can_cross,
            filter_max
        );
        let top3_j: Vec<_> = jr
            .iter()
            .take(3)
            .map(|r| json!({"fnv": format!("{:x}", r.fnv), "v": r.val}))
            .collect();
        let top3_r: Vec<_> = rr
            .iter()
            .take(3)
            .map(|r| json!({"fnv": format!("{:x}", r.fnv), "v": r.val}))
            .collect();
        rows.push(json!({
            "qname": q,
            "flags": flags,
            "qlen": qlen,
            "tag": if i < 10 { "JAVA_LIVE_ONLY" } else if i < 22 { "RUST_ONLY" } else { "BOTH_DROP" },
            "java_keep": jdec,
            "rust_post_keep": rdec,
            "java_max70": jmax70,
            "java_max68": jmax,
            "rust_max68": rmax,
            "java_win70": format!("{:x}", jwin70),
            "java_win68": format!("{:x}", jwin),
            "rust_win68": format!("{:x}", rwin),
            "same_winner": jwin == rwin,
            "java_margin": jmargin,
            "rust_margin": rmargin,
            "dist": dist,
            "can_cross": can_cross,
            "filter_max_ll": filter_max,
            "top3_j_double": top3_j,
            "top3_r_f64": top3_r,
        }));
    }

    let classification = classify(n_cross, n_class_change, n_winner_changes_keep_drop);
    eprintln!(
        "{}",
        json!({
            "classification": classification,
            "n_cells": n_cell,
            "max_abs_j_double_vs_rust": max_cell_abs,
            "mean_abs_j_double_vs_rust": if n_cell == 0 { 0.0 } else { sum_abs / n_cell as f64 },
            "min_required_delta": min_req,
            "n_residual_cross": n_cross,
            "n_keep_drop_change": n_class_change,
            "n_winner_diff": n_winner_diff,
            "n_winner_changes_keep_drop": n_winner_changes_keep_drop,
            "n_safe": n_safe,
            "pairhmm_investigation": if classification == "PAIRHMM_RESIDUAL_NOT_CAUSAL" {
                "NOT JUSTIFIED"
            } else {
                "JUSTIFIED — SEPARATE AUTHORIZATION REQUIRED"
            },
            "vcf_ad": vcf.samples.first().map(|s| s.ad.clone()),
            "rows": rows,
        })
    );

    assert_eq!(n_cell, 1496);
    assert_eq!(n_cross, 0);
    assert_eq!(n_class_change, 0);
    assert_eq!(n_winner_changes_keep_drop, 0);
    assert_eq!(n_safe, 24);
    assert!(min_req > max_cell_abs);
    assert_eq!(classification, "PAIRHMM_RESIDUAL_NOT_CAUSAL");
    let _ = outcome;
}
