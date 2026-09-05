//! 6R.99 coordinate-free: PairHMM primitive residual vs poorly-modeled −8.0.
//!
//! Java 4.4.0.0 (`2dbc025821bc5f686c423ff332a41e6cef892a77`) + GKL **0.8.8**.
//! Do not inspect PairHMM DP / SIMD. Production precision unchanged.
//!
//! `filterPoorlyModeledEvidence` (6R.93): DROP iff `max_ll < -8.0` (equality keeps).
//! `normalizeLikelihoods` floors losing cells; per-read max is unchanged (6R.96).
//! Therefore the primitive post-kernel max is the filter max on the same matrix.
//!
//! Lemma: `|max(A) − max(B)| ≤ max_i |A_i − B_i|`.
//! A cell residual cannot move `max_ll` farther than itself, so it cannot flip
//! KEEP/DROP unless some read’s `|max_ll + 8|` is ≤ that residual.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r99_pairhmm_residual_causality_contract
//! HOLDOUT_6R99=1 cargo test -p gatk-haplotypecaller --test holdout_6r99_pairhmm_residual_causality -- --nocapture
//! ```

use std::collections::HashMap;

const JAVA_FLOAT_DUMP: &str = include_str!("6r96_java_seq6_post_kernel.tsv");
const JAVA_DOUBLE_DUMP: &str = include_str!("6r98_java_seq6_double_post_kernel.tsv");
const JAVA_ONLY: [u64; 2] = [0xfa2d2442dde7f8ff, 0x48eb4b18de00d4fd];
const THRESHOLD: f64 = -8.0;
/// Proven 6R.98 22×68 max `|J_double − R_f64|` (exact-bit comparison, not a tolerance).
const MAX_J_DOUBLE_VS_RUST: f64 = 1.754_881_395_754_637_2e-8;

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

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum Band {
    Safe,
    Borderline,
    PotentiallyCausal,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum Classification {
    PairhmmResidualNotCausal,
    PairhmmResidualPotentiallyCausal,
    PairhmmWinnerSwitch,
}

fn java_keep(max_ll: f64) -> bool {
    !(max_ll < THRESHOLD)
}

fn required_delta(max_ll: f64) -> f64 {
    (max_ll - THRESHOLD).abs()
}

fn residual_can_cross(required: f64, max_residual: f64) -> bool {
    required <= max_residual
}

/// SAFE: required movement exceeds residual. BORDERLINE: same order of magnitude
/// (`required ≤ 10 × residual`) but still cannot cross. POTENTIALLY_CAUSAL: can cross.
fn band(required: f64, max_residual: f64) -> Band {
    if residual_can_cross(required, max_residual) {
        Band::PotentiallyCausal
    } else if required <= max_residual * 10.0 {
        Band::Borderline
    } else {
        Band::Safe
    }
}

fn classify(
    n_residual_cross: usize,
    n_keep_drop_change: usize,
    n_winner_changes_keep_drop: usize,
) -> Classification {
    if n_winner_changes_keep_drop > 0 {
        Classification::PairhmmWinnerSwitch
    } else if n_residual_cross > 0 || n_keep_drop_change > 0 {
        Classification::PairhmmResidualPotentiallyCausal
    } else {
        Classification::PairhmmResidualNotCausal
    }
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

fn row_vals(dump: &JavaDump, q: &str, flags: u16) -> Vec<f64> {
    dump.bits
        .get(&(q.to_string(), flags))
        .unwrap()
        .iter()
        .copied()
        .map(f64::from_bits)
        .collect()
}

fn max_over(vals: &[f64], haps: &[u64], skip: &[u64]) -> f64 {
    let mut best = f64::NEG_INFINITY;
    for (i, &h) in haps.iter().enumerate() {
        if skip.contains(&h) {
            continue;
        }
        if vals[i] > best {
            best = vals[i];
        }
    }
    best
}

fn argmax_fnv(vals: &[f64], haps: &[u64], skip: &[u64]) -> u64 {
    let mut best = f64::NEG_INFINITY;
    let mut arg = 0u64;
    for (i, &h) in haps.iter().enumerate() {
        if skip.contains(&h) {
            continue;
        }
        if vals[i] > best {
            best = vals[i];
            arg = h;
        }
    }
    arg
}

#[test]
fn forensic_6r99_drop_is_strict_less_equality_keeps() {
    assert!(java_keep(-8.0));
    assert!(!java_keep(-8.0 - 1e-15));
    assert!(java_keep(-7.999));
}

#[test]
fn forensic_6r99_max_shift_cannot_exceed_cell_residual() {
    let a: [f64; 3] = [-2.5, -9.0, -4.0];
    let b: [f64; 3] = [-2.5 + 1e-8, -9.0, -4.0];
    let cell_max = a
        .iter()
        .zip(b.iter())
        .map(|(x, y)| (x - y).abs())
        .fold(0.0f64, f64::max);
    let max_a = a.into_iter().fold(f64::NEG_INFINITY, f64::max);
    let max_b = b.into_iter().fold(f64::NEG_INFINITY, f64::max);
    assert!((max_a - max_b).abs() <= cell_max);
}

#[test]
fn forensic_6r99_residual_cross_is_required_delta_le_residual() {
    assert!(!residual_can_cross(0.37, MAX_J_DOUBLE_VS_RUST));
    assert!(residual_can_cross(1e-9, MAX_J_DOUBLE_VS_RUST));
    assert!(residual_can_cross(
        MAX_J_DOUBLE_VS_RUST,
        MAX_J_DOUBLE_VS_RUST
    ));
}

#[test]
fn forensic_6r99_band_safe_when_margin_exceeds_residual() {
    assert_eq!(band(0.37, MAX_J_DOUBLE_VS_RUST), Band::Safe);
    assert_eq!(
        band(MAX_J_DOUBLE_VS_RUST * 5.0, MAX_J_DOUBLE_VS_RUST),
        Band::Borderline
    );
    assert_eq!(
        band(MAX_J_DOUBLE_VS_RUST, MAX_J_DOUBLE_VS_RUST),
        Band::PotentiallyCausal
    );
}

#[test]
fn forensic_6r99_classify_not_causal_when_no_cross_and_no_keep_drop_change() {
    assert_eq!(classify(0, 0, 0), Classification::PairhmmResidualNotCausal);
    assert_eq!(
        classify(1, 0, 0),
        Classification::PairhmmResidualPotentiallyCausal
    );
    assert_eq!(classify(0, 1, 1), Classification::PairhmmWinnerSwitch);
}

#[test]
fn forensic_6r99_java_double_max68_cannot_cross_threshold_under_6r98_residual() {
    let jf = parse_java_dump(JAVA_FLOAT_DUMP, "6R96");
    let jd = parse_java_dump(JAVA_DOUBLE_DUMP, "6R98");
    assert_eq!(jf.hap_fnv, jd.hap_fnv);
    assert_eq!(jf.hap_fnv.len(), 70);
    let skip = JAVA_ONLY;
    let mut min_req = f64::INFINITY;
    let mut n_cross = 0usize;
    let mut n_class_flip = 0usize;
    let mut n_safe = 0usize;
    for &(q, flags) in READS {
        let vf = row_vals(&jf, q, flags);
        let vd = row_vals(&jd, q, flags);
        let mf = max_over(&vf, &jf.hap_fnv, &skip);
        let md = max_over(&vd, &jd.hap_fnv, &skip);
        assert_eq!(java_keep(mf), java_keep(md), "{q} float vs double keep");
        if java_keep(mf) != java_keep(md) {
            n_class_flip += 1;
        }
        let req = required_delta(md);
        min_req = min_req.min(req);
        if residual_can_cross(req, MAX_J_DOUBLE_VS_RUST) {
            n_cross += 1;
        }
        if band(req, MAX_J_DOUBLE_VS_RUST) == Band::Safe {
            n_safe += 1;
        }
        let m70 = max_over(&vd, &jd.hap_fnv, &[]);
        assert_eq!(
            m70.to_bits(),
            md.to_bits(),
            "{q} Java-only columns do not raise max"
        );
    }
    assert_eq!(n_class_flip, 0);
    assert_eq!(n_cross, 0);
    assert_eq!(n_safe, 22);
    assert!(min_req > MAX_J_DOUBLE_VS_RUST);
    assert_eq!(
        classify(n_cross, n_class_flip, 0),
        Classification::PairhmmResidualNotCausal
    );
}

#[test]
fn forensic_6r99_winner_tie_does_not_reclassify_keep_drop() {
    // Two haplotypes at the same max: argmax index may move; max_ll does not.
    let haps = [1u64, 2, 3];
    let a = [-7.45, -7.45, -12.0];
    let b = [-7.45, -7.45 + 1e-8, -12.0];
    let wa = argmax_fnv(&a, &haps, &[]);
    let wb = argmax_fnv(&b, &haps, &[]);
    let ma = max_over(&a, &haps, &[]);
    let mb = max_over(&b, &haps, &[]);
    assert_ne!(wa, wb);
    assert_eq!(java_keep(ma), java_keep(mb));
    assert!(!residual_can_cross(required_delta(ma), 1e-8));
}

#[test]
fn forensic_6r99_even_java_float_vs_double_residual_is_safe() {
    let jf = parse_java_dump(JAVA_FLOAT_DUMP, "6R96");
    let jd = parse_java_dump(JAVA_DOUBLE_DUMP, "6R98");
    let skip = JAVA_ONLY;
    let mut max_abs = 0.0f64;
    for &(q, flags) in READS {
        let vf = row_vals(&jf, q, flags);
        let vd = row_vals(&jd, q, flags);
        for (i, &h) in jf.hap_fnv.iter().enumerate() {
            if skip.contains(&h) {
                continue;
            }
            max_abs = max_abs.max((vf[i] - vd[i]).abs());
        }
    }
    assert!(max_abs > 0.0);
    let mut n_cross = 0usize;
    for &(q, flags) in READS {
        let vd = row_vals(&jd, q, flags);
        let md = max_over(&vd, &jd.hap_fnv, &skip);
        if residual_can_cross(required_delta(md), max_abs) {
            n_cross += 1;
        }
    }
    assert_eq!(
        n_cross, 0,
        "even the larger float-vs-double residual is SAFE"
    );
}
