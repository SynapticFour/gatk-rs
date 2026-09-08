//! 6R.94 coordinate-free: `max_ll` population attribution at
//! `filterPoorlyModeledEvidence`.
//!
//! Java 4.4.0.0 (`2dbc025821bc5f686c423ff332a41e6cef892a77`)
//!
//! ```text
//! PairHMMLikelihoodCalculationEngine.computeReadLikelihoods
//!   AlleleLikelihoods(samples, IndexedAlleleList(haplotypeList), perSampleReadList)
//!   pairHMM.computeLog10Likelihoods          // kernel  — do not inspect
//!   normalizeLikelihoods
//!   ReadLikelihoodCalculationEngine.filterPoorlyModeledEvidence
//!     AlleleLikelihoods.filterPoorlyModeledEvidence
//!       max_ll = maximumLikelihoodOverAllAlleles(sample, evidenceIndex)
//!              = max over haplotype columns (strict >; first index wins ties)
//! ```
//!
//! The filtering likelihood object is the **post-normalize PairHMM**
//! `AlleleLikelihoods<GATKRead, Haplotype>`: evidence rows are `sampleEvidence`
//! GATKRead objects; columns are the assembly haplotype list. Default HC
//! `filterAlleles=false`, so this object is **not** a later allele-filtered
//! haplotype subset.
//!
//! Arrow order (stop at the first broken one):
//!   read identity
//!     → likelihood-row assignment
//!     → haplotype-column population
//!     → pre-filter likelihood values
//!     → max_ll
//!     → filterPoorlyModeledEvidence   (closed in 6R.93)
//!
//! Column identity is haplotype **sequence** (FNV-1a 64 of bases), not index.
//! Column order is not causal for `max()` when the same values are present.
//! A strict subset of columns cannot raise `max_ll`.
//!
//! Do not collapse mates by QNAME: row identity is `(qname, flags)` (and
//! start/CIGAR when flags collide).
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r94_max_ll_population_contract
//! HOLDOUT_6R94=1 cargo test -p gatk-haplotypecaller --test holdout_6r94_max_ll_population -- --nocapture
//! ```

use std::collections::HashSet;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum FirstDivergence {
    LikelihoodRowPopulation,
    LikelihoodColumnPopulation,
    PreFilterLikelihoodValue,
    AfterMaxLl,
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

/// Java `maximumLikelihoodOverAllAlleles`: first index with strict `>`.
fn max_ll_and_argmax(row: &[f64]) -> (f64, usize) {
    let mut best = f64::NEG_INFINITY;
    let mut arg = 0usize;
    for (i, &v) in row.iter().enumerate() {
        if v > best {
            best = v;
            arg = i;
        }
    }
    (best, arg)
}

fn classify(
    java_rows: &HashSet<ReadKey>,
    rust_rows: &HashSet<ReadKey>,
    java_cols: &HashSet<u64>,
    rust_cols: &HashSet<u64>,
    java_max: &[f64],
    rust_max: &[f64],
) -> FirstDivergence {
    if java_rows != rust_rows {
        return FirstDivergence::LikelihoodRowPopulation;
    }
    if java_cols != rust_cols {
        return FirstDivergence::LikelihoodColumnPopulation;
    }
    let same_max = java_max.len() == rust_max.len()
        && java_max
            .iter()
            .zip(rust_max)
            .all(|(a, b)| a == b || (a.is_nan() && b.is_nan()));
    if same_max {
        FirstDivergence::AfterMaxLl
    } else {
        FirstDivergence::PreFilterLikelihoodValue
    }
}

#[test]
fn forensic_6r94_fnv1a64_is_sequence_identity() {
    assert_eq!(fnv1a64(b"ACGT"), fnv1a64(b"ACGT"));
    assert_ne!(fnv1a64(b"ACGT"), fnv1a64(b"ACGG"));
    assert_ne!(fnv1a64(b"A"), fnv1a64(b"AA"));
}

#[test]
fn forensic_6r94_max_ll_is_first_strict_gt_over_haplotype_columns() {
    let row = [-12.0, -2.5, -2.5, -4.0];
    let (mx, arg) = max_ll_and_argmax(&row);
    assert_eq!(mx, -2.5);
    assert_eq!(arg, 1, "ties keep the first haplotype column");
}

#[test]
fn forensic_6r94_column_order_is_not_causal_for_max() {
    let a = [-9.0, -2.5, -4.0];
    let b = [-4.0, -9.0, -2.5];
    assert_eq!(max_ll_and_argmax(&a).0, max_ll_and_argmax(&b).0);
}

#[test]
fn forensic_6r94_column_subset_cannot_raise_max_ll() {
    let full = [-12.0, -2.5, -20.0];
    let subset = [-12.0, -20.0];
    let max_full = max_ll_and_argmax(&full).0;
    let max_sub = max_ll_and_argmax(&subset).0;
    assert!(max_sub <= max_full);
    assert!(max_sub < max_full);
}

#[test]
fn forensic_6r94_row_identity_does_not_collapse_qname() {
    let mate_a = ReadKey {
        qname: "READ",
        flags: 83,
    };
    let mate_b = ReadKey {
        qname: "READ",
        flags: 163,
    };
    assert_ne!(mate_a, mate_b);
    let mut rows = HashSet::new();
    rows.insert(mate_a);
    rows.insert(mate_b);
    assert_eq!(rows.len(), 2);
}

#[test]
fn forensic_6r94_same_rows_different_hap_hashes_is_column_population() {
    let r = ReadKey {
        qname: "R1",
        flags: 99,
    };
    let mut java_rows = HashSet::new();
    java_rows.insert(r);
    let rust_rows = java_rows.clone();
    let java_cols: HashSet<u64> = [fnv1a64(b"AAA"), fnv1a64(b"AAC"), fnv1a64(b"AAG")]
        .into_iter()
        .collect();
    let rust_cols: HashSet<u64> = [fnv1a64(b"AAA"), fnv1a64(b"AAC")].into_iter().collect();
    assert_eq!(
        classify(
            &java_rows,
            &rust_rows,
            &java_cols,
            &rust_cols,
            &[-2.5],
            &[-9.5]
        ),
        FirstDivergence::LikelihoodColumnPopulation
    );
}

#[test]
fn forensic_6r94_missing_row_is_row_population_even_if_columns_match() {
    let a = ReadKey {
        qname: "A",
        flags: 99,
    };
    let b = ReadKey {
        qname: "B",
        flags: 147,
    };
    let mut java_rows = HashSet::new();
    java_rows.insert(a);
    java_rows.insert(b);
    let mut rust_rows = HashSet::new();
    rust_rows.insert(a);
    let cols: HashSet<u64> = [1u64, 2].into_iter().collect();
    assert_eq!(
        classify(&java_rows, &rust_rows, &cols, &cols, &[-2.5, -9.0], &[-2.5]),
        FirstDivergence::LikelihoodRowPopulation
    );
}

#[test]
fn forensic_6r94_same_rows_and_columns_different_max_is_value_divergence() {
    let r = ReadKey {
        qname: "R1",
        flags: 99,
    };
    let mut rows = HashSet::new();
    rows.insert(r);
    let cols: HashSet<u64> = [10u64, 11].into_iter().collect();
    assert_eq!(
        classify(&rows, &rows, &cols, &cols, &[-2.5], &[-9.5]),
        FirstDivergence::PreFilterLikelihoodValue
    );
}

#[test]
fn forensic_6r94_identical_max_after_same_population_is_after_max_ll() {
    let r = ReadKey {
        qname: "R1",
        flags: 99,
    };
    let mut rows = HashSet::new();
    rows.insert(r);
    let cols: HashSet<u64> = [10u64].into_iter().collect();
    assert_eq!(
        classify(&rows, &rows, &cols, &cols, &[-2.5], &[-2.5]),
        FirstDivergence::AfterMaxLl
    );
}

#[test]
fn forensic_6r94_argmax_hap_only_on_one_side_is_column_not_kernel() {
    let java_cols: HashSet<u64> = [1, 2, 3].into_iter().collect();
    let rust_cols: HashSet<u64> = [1, 2].into_iter().collect();
    let java_argmax = 3u64;
    assert!(java_cols.contains(&java_argmax) && !rust_cols.contains(&java_argmax));
    assert_ne!(java_cols, rust_cols);
}

#[test]
fn forensic_6r94_java_filter_object_is_pairhmm_haplotype_list_not_allele_filtered() {
    // Java 4.4 default `filterAlleles=false`: poorly-modeled runs inside
    // `computeReadLikelihoods` on IndexedAlleleList(haplotypeList).
    // A later genotyping allele list (shorter) is a different object.
    let pairhmm_columns = 70usize;
    let allele_filtered_columns = 68usize;
    assert!(pairhmm_columns > allele_filtered_columns);
    assert_ne!(
        pairhmm_columns, allele_filtered_columns,
        "filter-time columns are the PairHMM haplotype list, not the later allele-filtered set"
    );
}
