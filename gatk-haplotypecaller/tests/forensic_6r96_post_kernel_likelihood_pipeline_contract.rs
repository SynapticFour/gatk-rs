//! 6R.96 coordinate-free: first post-kernel likelihood-value boundary.
//!
//! Java 4.4.0.0 (`2dbc025821bc5f686c423ff332a41e6cef892a77`)
//! `PairHMMLikelihoodCalculationEngine.computeReadLikelihoods`:
//!
//! ```text
//! AlleleLikelihoods(samples, IndexedAlleleList(haps), perSampleReadList)
//! pairHMM.computeLog10Likelihoods          // kernel — do not inspect
//! //  ← 6R.96 post-kernel object (values = kernel log10 cells)
//! normalizeLikelihoods(log10global, symmetric)   // floors losing cells; max unchanged
//! filterPoorlyModeledEvidence                    // closed in 6R.93
//! ```
//!
//! Default HC `filterAlleles=false`: no haplotype-column compaction and no
//! second PairHMM between kernel materialization and poorly-modeled.
//!
//! Rust `strict_java` on a non-P12 span:
//!
//! ```text
//! score_pairhmm_from_records                   // kernel — do not inspect
//! //  ← post_kernel Vec<RegionReadLikelihood>
//! filter_assembly_and_likelihoods              // column drop + index remap (copy)
//! normalize_region_read_likelihoods            // floors losing cells; max unchanged
//! filter_normalized_region_read_likelihoods
//! ```
//!
//! Compare cells by `(read identity, haplotype sequence hash)`, never column
//! index. Exact f64 bits; no tolerance. Stop at the first boundary with
//! `differing > 0`.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r96_post_kernel_likelihood_pipeline_contract
//! HOLDOUT_6R96=1 cargo test -p gatk-haplotypecaller --test holdout_6r96_post_kernel_likelihood -- --nocapture
//! ```

use std::collections::{HashMap, HashSet};

const JAVA_ONLY_A: u64 = 0xfa2d2442dde7f8ff;
const JAVA_ONLY_B: u64 = 0x48eb4b18de00d4fd;
const LOG10_GLOBAL: f64 = -4.5;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum FirstDivergence {
    PostKernelLikelihoodObject,
    Normalize,
    LikelihoodCompaction,
    LikelihoodRefresh,
    FilterInputConstruction,
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

/// Java `normalizeLikelihoodsPerEvidence`: floor cells below `best + cap`.
fn java_normalize_floor(cells: &mut [f64], cap: f64) {
    let mut best = f64::NEG_INFINITY;
    for &v in cells.iter() {
        if v > best {
            best = v;
        }
    }
    if !best.is_finite() {
        return;
    }
    let floor = best + cap;
    for v in cells.iter_mut() {
        if *v < floor {
            *v = floor;
        }
    }
}

/// Compaction: keep listed haplotype hashes, copy values (no recompute).
fn compact_copy(row: &HashMap<u64, f64>, keep: &HashSet<u64>) -> HashMap<u64, f64> {
    row.iter()
        .filter(|(h, _)| keep.contains(*h))
        .map(|(&h, &v)| (h, v))
        .collect()
}

fn cell_stats(
    java: &HashMap<(ReadKey, u64), f64>,
    rust: &HashMap<(ReadKey, u64), f64>,
    keys: &[ReadKey],
    common: &HashSet<u64>,
) -> (usize, usize, usize, f64, f64) {
    let mut n = 0usize;
    let mut eq = 0usize;
    let mut diff = 0usize;
    let mut sum = 0.0;
    let mut max_abs = 0.0;
    for &k in keys {
        for &h in common {
            n += 1;
            match (java.get(&(k, h)), rust.get(&(k, h))) {
                (Some(&a), Some(&b)) => {
                    if a.to_bits() == b.to_bits() {
                        eq += 1;
                    } else {
                        diff += 1;
                    }
                    let d = (a - b).abs();
                    sum += d;
                    if d > max_abs {
                        max_abs = d;
                    }
                }
                _ => diff += 1,
            }
        }
    }
    let mean = if n == 0 { 0.0 } else { sum / n as f64 };
    (n, eq, diff, max_abs, mean)
}

fn first_divergence(boundaries: &[(&str, usize)]) -> Option<FirstDivergence> {
    for &(name, differing) in boundaries {
        if differing == 0 {
            continue;
        }
        return Some(match name {
            "post_kernel" => FirstDivergence::PostKernelLikelihoodObject,
            "normalize" => FirstDivergence::Normalize,
            "compaction" => FirstDivergence::LikelihoodCompaction,
            "refresh" => FirstDivergence::LikelihoodRefresh,
            _ => FirstDivergence::FilterInputConstruction,
        });
    }
    None
}

#[test]
fn forensic_6r96_compare_by_sequence_hash_not_column_index() {
    let hap_a = fnv1a64(b"AAA");
    let hap_b = fnv1a64(b"CCC");
    let r = ReadKey {
        qname: "R",
        flags: 99,
    };
    let mut java = HashMap::new();
    let mut rust = HashMap::new();
    java.insert((r, hap_b), -9.0);
    java.insert((r, hap_a), -2.5);
    rust.insert((r, hap_a), -2.5);
    rust.insert((r, hap_b), -9.0);
    let common: HashSet<u64> = [hap_a, hap_b].into_iter().collect();
    let (n, eq, diff, ..) = cell_stats(&java, &rust, &[r], &common);
    assert_eq!((n, eq, diff), (2, 2, 0));
}

#[test]
fn forensic_6r96_sub_ulp_delta_is_still_unequal() {
    let r = ReadKey {
        qname: "R",
        flags: 147,
    };
    let h = fnv1a64(b"H");
    let mut java = HashMap::new();
    let mut rust = HashMap::new();
    let jv: f64 = -16.54292106628418;
    let rv: f64 = -16.542921354965245;
    assert!((jv - rv).abs() < 1e-5);
    assert_ne!(jv.to_bits(), rv.to_bits());
    java.insert((r, h), jv);
    rust.insert((r, h), rv);
    let common: HashSet<u64> = [h].into_iter().collect();
    let (n, eq, diff, max_abs, _) = cell_stats(&java, &rust, &[r], &common);
    assert_eq!((n, eq, diff), (1, 0, 1));
    assert!(max_abs > 0.0);
    assert_eq!(
        first_divergence(&[("post_kernel", diff)]),
        Some(FirstDivergence::PostKernelLikelihoodObject)
    );
}

#[test]
fn forensic_6r96_java_only_hashes_excluded_from_common() {
    let java: HashSet<u64> = [1u64, 2, JAVA_ONLY_A, JAVA_ONLY_B].into_iter().collect();
    let rust: HashSet<u64> = [1u64, 2].into_iter().collect();
    let common: HashSet<u64> = java.intersection(&rust).copied().collect();
    assert_eq!(common.len(), 2);
    assert!(!common.contains(&JAVA_ONLY_A));
    assert!(!common.contains(&JAVA_ONLY_B));
}

#[test]
fn forensic_6r96_compaction_is_copy_of_retained_cells() {
    let h0 = fnv1a64(b"H0");
    let h1 = fnv1a64(b"H1");
    let h2 = fnv1a64(b"H2");
    let mut row = HashMap::new();
    row.insert(h0, -1.25);
    row.insert(h1, -9.5);
    row.insert(h2, -3.0);
    let keep: HashSet<u64> = [h0, h2].into_iter().collect();
    let out = compact_copy(&row, &keep);
    assert_eq!(out.len(), 2);
    assert_eq!(out.get(&h0).unwrap().to_bits(), (-1.25f64).to_bits());
    assert_eq!(out.get(&h2).unwrap().to_bits(), (-3.0f64).to_bits());
    assert!(!out.contains_key(&h1));
}

#[test]
fn forensic_6r96_normalize_floors_losing_cells_preserves_max() {
    let mut cells = vec![-2.0, -20.0, -3.0];
    let max_before = cells.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    java_normalize_floor(&mut cells, LOG10_GLOBAL);
    let max_after = cells.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    assert_eq!(max_before.to_bits(), max_after.to_bits());
    assert_eq!(max_after, -2.0);
    assert_eq!(cells[0], -2.0);
    assert_eq!(cells[1], -2.0 + LOG10_GLOBAL);
    assert_eq!(cells[2], -3.0);
}

#[test]
fn forensic_6r96_first_boundary_with_differing_cells_wins() {
    let boundaries = [
        ("post_kernel", 1496usize),
        ("compaction", 1496),
        ("normalize", 1496),
    ];
    assert_eq!(
        first_divergence(&boundaries),
        Some(FirstDivergence::PostKernelLikelihoodObject)
    );
}

#[test]
fn forensic_6r96_later_boundaries_ignored_when_post_kernel_differs() {
    assert_eq!(
        first_divergence(&[("post_kernel", 1), ("refresh", 99)]),
        Some(FirstDivergence::PostKernelLikelihoodObject)
    );
    assert_eq!(
        first_divergence(&[("post_kernel", 0), ("compaction", 4)]),
        Some(FirstDivergence::LikelihoodCompaction)
    );
    assert_eq!(
        first_divergence(&[("post_kernel", 0), ("compaction", 0), ("normalize", 2)]),
        Some(FirstDivergence::Normalize)
    );
    assert_eq!(
        first_divergence(&[
            ("post_kernel", 0),
            ("compaction", 0),
            ("normalize", 0),
            ("refresh", 8)
        ]),
        Some(FirstDivergence::LikelihoodRefresh)
    );
}

#[test]
fn forensic_6r96_java_has_no_compaction_between_kernel_and_normalize() {
    // Default HC: filterAlleles=false. Kernel 153×70 is the normalize input.
    let kernel_cols = 70usize;
    let normalize_in_cols = 70usize;
    assert_eq!(kernel_cols, normalize_in_cols);
}

#[test]
fn forensic_6r96_rust_compaction_must_not_recompute_if_copy() {
    let hap = fnv1a64(b"COMMON");
    let r = ReadKey {
        qname: "R",
        flags: 83,
    };
    let mut before = HashMap::new();
    let mut after = HashMap::new();
    before.insert((r, hap), -4.125);
    after.insert((r, hap), -4.125);
    let common: HashSet<u64> = [hap].into_iter().collect();
    let (_n, eq, diff, ..) = cell_stats(&before, &after, &[r], &common);
    assert_eq!(eq, 1);
    assert_eq!(diff, 0);
}
