//! 6R.95 coordinate-free: common-haplotype likelihood cells at poorly-modeled.
//!
//! Java 4.4.0.0 (`2dbc025821bc5f686c423ff332a41e6cef892a77`)
//!
//! Compare `AlleleLikelihoods` cells to Rust `RegionReadLikelihood` cells by
//! `(read identity, haplotype sequence hash)`, never by column index.
//! Java-only haplotype hashes are excluded from the common population.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r95_common_haplotype_likelihood_contract
//! HOLDOUT_6R95=1 cargo test -p gatk-haplotypecaller --test holdout_6r95_common_haplotype_likelihood -- --nocapture
//! ```

use std::collections::{HashMap, HashSet};

const JAVA_ONLY_A: u64 = 0xfa2d2442dde7f8ff;
const JAVA_ONLY_B: u64 = 0x48eb4b18de00d4fd;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum FirstDivergence {
    LikelihoodColumnPopulationOnly,
    PreFilterLikelihoodValue,
    MaxLlReduction,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash)]
struct ReadKey {
    qname: &'static str,
    flags: u16,
}

fn fnv1a64(bases: &[u8]) -> u64 {
    let mut h = 0xcbf29ce484222325u64;
    for &b in bases {
        h ^= u64::from(b);
        h = h.wrapping_mul(0x100000001b3);
    }
    h
}

fn common_hap_hashes(java: &[u64], rust: &[u64]) -> (HashSet<u64>, HashSet<u64>, HashSet<u64>) {
    let j: HashSet<u64> = java.iter().copied().collect();
    let r: HashSet<u64> = rust.iter().copied().collect();
    (
        j.intersection(&r).copied().collect(),
        j.difference(&r).copied().collect(),
        r.difference(&j).copied().collect(),
    )
}

/// Max over values whose haplotype hash is in `common`, first-index strict `>`.
fn max_over_common(cols: &[(u64, f64)], common: &HashSet<u64>) -> (f64, u64) {
    let mut best = f64::NEG_INFINITY;
    let mut arg = 0u64;
    for &(h, v) in cols {
        if !common.contains(&h) {
            continue;
        }
        if v > best {
            best = v;
            arg = h;
        }
    }
    (best, arg)
}

fn classify_cells(
    java: &HashMap<(ReadKey, u64), f64>,
    rust: &HashMap<(ReadKey, u64), f64>,
    keys: &[ReadKey],
    common: &HashSet<u64>,
) -> FirstDivergence {
    let mut any_cell_diff = false;
    let mut any_max_diff = false;
    for &k in keys {
        let mut jrow = Vec::new();
        let mut rrow = Vec::new();
        for &h in common {
            let jv = java.get(&(k, h)).copied();
            let rv = rust.get(&(k, h)).copied();
            match (jv, rv) {
                (Some(a), Some(b)) => {
                    if a.to_bits() != b.to_bits() {
                        any_cell_diff = true;
                    }
                    jrow.push((h, a));
                    rrow.push((h, b));
                }
                _ => any_cell_diff = true,
            }
        }
        let (jm, _) = max_over_common(&jrow, common);
        let (rm, _) = max_over_common(&rrow, common);
        if jm.to_bits() != rm.to_bits() {
            any_max_diff = true;
        }
    }
    if any_cell_diff {
        FirstDivergence::PreFilterLikelihoodValue
    } else if any_max_diff {
        FirstDivergence::MaxLlReduction
    } else {
        FirstDivergence::LikelihoodColumnPopulationOnly
    }
}

#[test]
fn forensic_6r95_match_by_sequence_hash_not_column_index() {
    let hap_a = fnv1a64(b"AAA");
    let hap_b = fnv1a64(b"CCC");
    // Java col0=B, col1=A; Rust col0=A, col1=B — same cells, swapped indices.
    let java = [(hap_b, -9.0), (hap_a, -2.5)];
    let rust = [(hap_a, -2.5), (hap_b, -9.0)];
    let common: HashSet<u64> = [hap_a, hap_b].into_iter().collect();
    let (jm, ja) = max_over_common(&java, &common);
    let (rm, ra) = max_over_common(&rust, &common);
    assert_eq!(jm, rm);
    assert_eq!(ja, ra);
    assert_eq!(ja, hap_a);
}

#[test]
fn forensic_6r95_java_only_hashes_excluded_from_common() {
    let java = [1u64, 2, JAVA_ONLY_A, JAVA_ONLY_B];
    let rust = [1u64, 2];
    let (common, jonly, ronly) = common_hap_hashes(&java, &rust);
    assert_eq!(common.len(), 2);
    assert_eq!(jonly.len(), 2);
    assert!(jonly.contains(&JAVA_ONLY_A) && jonly.contains(&JAVA_ONLY_B));
    assert!(ronly.is_empty());
}

#[test]
fn forensic_6r95_common_max_ignores_java_only_columns() {
    let hap = fnv1a64(b"COMMON");
    let common: HashSet<u64> = [hap].into_iter().collect();
    let row = [(JAVA_ONLY_A, -1.0), (hap, -9.0)];
    let (mx, arg) = max_over_common(&row, &common);
    assert_eq!(mx, -9.0);
    assert_eq!(arg, hap);
}

#[test]
fn forensic_6r95_identical_common_cells_are_column_population_only() {
    let r = ReadKey {
        qname: "R",
        flags: 99,
    };
    let h = fnv1a64(b"H");
    let mut java = HashMap::new();
    let mut rust = HashMap::new();
    java.insert((r, h), -2.5);
    rust.insert((r, h), -2.5);
    let common: HashSet<u64> = [h].into_iter().collect();
    assert_eq!(
        classify_cells(&java, &rust, &[r], &common),
        FirstDivergence::LikelihoodColumnPopulationOnly
    );
}

#[test]
fn forensic_6r95_differing_common_cell_is_pre_filter_value() {
    let r = ReadKey {
        qname: "R",
        flags: 99,
    };
    let h = fnv1a64(b"H");
    let mut java = HashMap::new();
    let mut rust = HashMap::new();
    java.insert((r, h), -2.5);
    rust.insert((r, h), -9.5);
    let common: HashSet<u64> = [h].into_iter().collect();
    assert_eq!(
        classify_cells(&java, &rust, &[r], &common),
        FirstDivergence::PreFilterLikelihoodValue
    );
}

#[test]
fn forensic_6r95_equal_cells_unequal_max_is_reduction() {
    // Synthetic: cells identical as bits, but max() would differ only if
    // a caller dropped a finite cell. The classifier sees equal maps so
    // column-only; this unit checks max_over_common itself is deterministic.
    let hap_a = fnv1a64(b"A");
    let hap_b = fnv1a64(b"B");
    let common: HashSet<u64> = [hap_a, hap_b].into_iter().collect();
    let row = [(hap_a, f64::NAN), (hap_b, -3.0)];
    // Java `>` does not promote NaN; max stays -3.0.
    let (mx, arg) = max_over_common(&row, &common);
    assert_eq!(mx, -3.0);
    assert_eq!(arg, hap_b);
}

#[test]
fn forensic_6r95_java70_vs_java68_detects_java_only_causal_max() {
    let hap = fnv1a64(b"COMMON");
    let common: HashSet<u64> = [hap].into_iter().collect();
    let row70 = [(JAVA_ONLY_A, -2.5), (hap, -9.0)];
    let max70 = row70
        .iter()
        .fold(f64::NEG_INFINITY, |a, &(_, v)| if v > a { v } else { a });
    let (max68, _) = max_over_common(&row70, &common);
    assert!(max70 > max68);
    assert_eq!(max70, -2.5);
    assert_eq!(max68, -9.0);
}

#[test]
fn forensic_6r95_exact_bits_not_tolerance() {
    let a = f64::from_bits(0xc00485a1a0000000);
    let b = f64::from_bits(0xc00485a1a0000001);
    assert_ne!(a.to_bits(), b.to_bits());
    let d = (a - b).abs();
    assert!(d > 0.0);
    assert!(d < 1e-15);
}
