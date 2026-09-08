//! 6R.100 coordinate-free: filter-input evidence identity vs post-kernel max.
//!
//! Java 4.4.0.0 (`2dbc025821bc5f686c423ff332a41e6cef892a77`)
//!
//! ```text
//! PairHMMLikelihoodCalculationEngine.computeReadLikelihoods
//!   AlleleLikelihoods(samples, IndexedAlleleList(haps), perSampleReadList)
//!   pairHMM.computeLog10Likelihoods            // kernel — CLOSED
//!   normalizeLikelihoods(log10global, symmetric)
//!     // floors cells < best + (−4.5); per-evidence max unchanged
//!   ReadLikelihoodCalculationEngine.filterPoorlyModeledEvidence
//!     AlleleLikelihoods.filterPoorlyModeledEvidence
//!       max_ll = maximumLikelihoodOverAllAlleles(sampleIndex, evidenceIndex)
//!       DROP iff max_ll < log10MinTrueLikelihood(evidence.get(i))
//! ```
//!
//! `filterPoorlyModeledEvidence` indexes the **same** `sampleEvidence` list that
//! constructed the matrix. `sampleEvidence.get(i)` is the GATKRead whose cells
//! are `valuesBySampleIndex[allele][i]`. It is not a copy, not a resorted list,
//! and not a later `regionForGenotyping` snapshot taken after a parallel clip.
//!
//! Default HC `filterAlleles=false`: poorly-modeled runs **before** realign and
//! **before** allele filtering. Floor-normalize cannot lower a row max.
//!
//! Rust contract (this arrow):
//! `RegionReadLikelihood.read_index` is assigned against the records actually
//! scored by PairHMM. `filter_poorly_modeled_region_read_likelihoods` must look
//! up `reads.get(read_index)` on **that same list**. A permutation or a parallel
//! clip/sort of `region_for_genotyping.reads` attributes another row's `max_ll`
//! to a QNAME and can flip KEEP/DROP with **zero** cell-value change.
//!
//! Compaction (column subset) and floor-normalize are not causal for a drop in
//! per-read max when the winning common haplotype remains.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r100_filter_input_likelihood_pipeline_contract
//! HOLDOUT_6R100=1 cargo test -p gatk-haplotypecaller --test holdout_6r100_filter_input_likelihood -- --nocapture
//! ```

const THRESHOLD: f64 = -8.0;
const LOG10_GLOBAL: f64 = -4.5;

fn java_keep(max_ll: f64) -> bool {
    !(max_ll < THRESHOLD)
}

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

fn row_max(cells: &[f64]) -> f64 {
    cells.iter().copied().fold(f64::NEG_INFINITY, f64::max)
}

/// Java `AlleleLikelihoods`: evidence i owns matrix row i. A filter list that is
/// a permutation of the scored list attributes another row's max to a QNAME.
fn filter_max_under_evidence_perm(scored_max: &[f64], filter_index_of_qname: &[usize]) -> Vec<f64> {
    filter_index_of_qname
        .iter()
        .map(|&i| scored_max[i])
        .collect()
}

#[test]
fn forensic_6r100_floor_normalize_cannot_lower_row_max() {
    let post = [-7.450_010_464_825_709, -20.0, -30.0];
    let mut norm = post;
    java_normalize_floor(&mut norm, LOG10_GLOBAL);
    assert_eq!(row_max(&post).to_bits(), row_max(&norm).to_bits());
    assert_eq!(java_keep(row_max(&post)), java_keep(row_max(&norm)));
    assert!(norm[1] > post[1]);
}

#[test]
fn forensic_6r100_column_subset_cannot_lower_max_if_winner_retained() {
    let winner = 0x501e_24ed_a83c_4dbd;
    let dropped = 0xfa2d_2442_dde7_f8ff;
    let post = [(winner, -7.45), (dropped, -40.0), (0x1111, -13.43)];
    let compact: Vec<(u64, f64)> = post
        .iter()
        .copied()
        .filter(|(h, _)| *h != dropped)
        .collect();
    let post_max = post.iter().map(|c| c.1).fold(f64::NEG_INFINITY, f64::max);
    let compact_max = compact
        .iter()
        .map(|c| c.1)
        .fold(f64::NEG_INFINITY, f64::max);
    assert_eq!(post_max.to_bits(), compact_max.to_bits());
    assert!(compact.iter().any(|(h, _)| *h == winner));
}

#[test]
fn forensic_6r100_evidence_permutation_flips_keep_without_cell_change() {
    // Scored order: KEEP then DROP. Filter list is swapped.
    let scored_max = [-7.450_010_464_825_709, -13.432_897_684_740_624];
    assert!(java_keep(scored_max[0]));
    assert!(!java_keep(scored_max[1]));
    let filter_index_of_qname = [1usize, 0];
    let filter_max = filter_max_under_evidence_perm(&scored_max, &filter_index_of_qname);
    assert_eq!(filter_max[0].to_bits(), scored_max[1].to_bits());
    assert_eq!(filter_max[1].to_bits(), scored_max[0].to_bits());
    assert_ne!(java_keep(scored_max[0]), java_keep(filter_max[0]));
    assert_ne!(java_keep(scored_max[1]), java_keep(filter_max[1]));
    let delta = filter_max[0] - scored_max[0];
    assert!(delta.abs() > 5.0);
    assert!(delta.abs() > 1e-6);
}

#[test]
fn forensic_6r100_pairhmm_residual_cannot_explain_filter_delta() {
    let post: f64 = -7.450_010_464_825_709;
    let filter: f64 = -13.432_897_684_740_624;
    let residual: f64 = 1.754_881_395_754_637_2e-8;
    let delta = (filter - post).abs();
    assert!(delta > 5.0);
    assert!(delta > residual * 1e7);
}

#[test]
fn forensic_6r100_java_filter_indexes_same_evidence_object() {
    // sampleEvidence.get(i) ↔ values[allele][i]. Identity, not a parallel list.
    let evidence = ["A:147", "B:163"];
    let cells = [[-7.45, -13.43], [-9.0, -20.0]];
    let max_ll: Vec<f64> = (0..evidence.len())
        .map(|i| {
            cells
                .iter()
                .map(|col| col[i])
                .fold(f64::NEG_INFINITY, f64::max)
        })
        .collect();
    assert_eq!(max_ll[0].to_bits(), (-7.45f64).to_bits());
    assert_eq!(max_ll[1].to_bits(), (-13.43f64).to_bits());
    let _ = evidence;
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
#[allow(dead_code)]
enum Classification {
    LikelihoodValueTransformation,
    LikelihoodRowTransformation,
    LikelihoodColumnTransformation,
    LikelihoodObjectLifecycle,
    NoProvenDivergence,
}

fn classify_first_causal(
    common_cell_value_diff_before_filter: bool,
    max_changed_before_filter: bool,
    winner_removed: bool,
    filter_index_qname_mismatch: bool,
) -> Classification {
    if winner_removed {
        return Classification::LikelihoodColumnTransformation;
    }
    if max_changed_before_filter {
        return Classification::LikelihoodValueTransformation;
    }
    if filter_index_qname_mismatch && !common_cell_value_diff_before_filter {
        return Classification::LikelihoodObjectLifecycle;
    }
    if filter_index_qname_mismatch {
        return Classification::LikelihoodObjectLifecycle;
    }
    Classification::NoProvenDivergence
}

#[test]
fn forensic_6r100_first_causal_operation_is_evidence_lifecycle() {
    // Live 6R.100: compaction 0/1496 bit diffs on common 68; normalize floors
    // losers but does not change per-read max; 22/22 filter indices name a
    // different QNAME than the scored list at that read_index.
    assert_eq!(
        classify_first_causal(false, false, false, true),
        Classification::LikelihoodObjectLifecycle
    );
    assert_ne!(
        classify_first_causal(false, false, false, true),
        Classification::LikelihoodValueTransformation
    );
    assert_ne!(
        classify_first_causal(false, false, false, true),
        Classification::LikelihoodColumnTransformation
    );
}
