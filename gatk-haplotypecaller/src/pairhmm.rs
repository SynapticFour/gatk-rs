//! PairHMM likelihood scaffold.
//! This module provides a deterministic log10-space PairHMM with:
//! explicit Match/Insertion/Deletion state transitions,
//! base + mapping quality integration via effective quality capping,
//! numerically stable log-sum transitions.

use gatk_common::{GatkError, GatkResult};
use rayon::prelude::*;

const LOG10_NEG_INF: f64 = f64::NEG_INFINITY;

/// PairHMM transition/emission probabilities in linear probability space.
/// # Invariants
/// Probabilities are ∈ (0, 1] where used; converted to log10 at scoring time.
/// # Ownership
/// [`Copy`] parameter bundle.
/// # Mutation
/// Immutable per scoring call unless caller replaces fields.
/// # Biological assumptions
/// Standard gap-open/extend and uniform insertion emission for short-read alignment.
/// # Java equivalence
/// GATK `PairHMM` / `Log10PairHMM` gap and indel emission defaults.
#[derive(Debug, Clone, Copy)]
pub struct PairHmmParams {
    pub gap_open_prob: f64,
    pub gap_extend_prob: f64,
    pub insertion_emission_prob: f64,
}

/// Floating-point tolerance policy when comparing PairHMM likelihoods.
/// # Invariants
/// `abs_epsilon` and `rel_epsilon` are positive comparison thresholds for parity checks.
/// # Ownership
/// [`Copy`] policy snapshot.
/// # Mutation
/// Immutable per comparison.
/// # Biological assumptions
/// None — numeric parity policy only.
/// # Java equivalence
/// Rust-native parity helper for PairHMM float comparisons (algorithm parity gates).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PairHmmFpPolicy {
    pub abs_epsilon: f64,
    pub rel_epsilon: f64,
}

impl Default for PairHmmFpPolicy {
    fn default() -> Self {
        Self {
            abs_epsilon: 1e-10,
            rel_epsilon: 1e-8,
        }
    }
}

impl Default for PairHmmParams {
    fn default() -> Self {
        Self {
            gap_open_prob: 1e-2,
            gap_extend_prob: 1e-1,
            insertion_emission_prob: 0.25,
        }
    }
}

/// One read × haplotype PairHMM scoring input.
/// # Invariants
/// `read_bases.len == read_base_quals.len`; haplotype must be non-empty (validated before score).
/// # Ownership
/// Owns read/haplotype strings and quality vector.
/// # Mutation
/// Immutable per likelihood evaluation unless caller mutates fields.
/// # Biological assumptions
/// Read and haplotype are aligned colinearly in ACGT(N) with Phred qualities and MAPQ.
/// # Java equivalence
/// GATK `PairHMM` input read + haplotype bytes/quals (`ReadLikelihoodCalculationEngine`).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PairHmmInput {
    pub read_bases: String,
    pub read_base_quals: Vec<u8>,
    pub read_mapping_quality: u8,
    pub haplotype_bases: String,
}

impl PairHmmInput {
    fn validate(&self) -> GatkResult<()> {
        if self.read_bases.len() != self.read_base_quals.len() {
            return Err(GatkError::argument(
                "PairHMM input mismatch: read bases length must match base quality length",
            ));
        }
        if self.haplotype_bases.is_empty() {
            return Err(GatkError::argument(
                "PairHMM input invalid: haplotype_bases must be non-empty",
            ));
        }
        Ok(())
    }
}

#[inline]
fn prob_to_log10(p: f64) -> f64 {
    if p <= 0.0 {
        LOG10_NEG_INF
    } else {
        p.log10()
    }
}

#[inline]
fn log10_sum(a: f64, b: f64) -> f64 {
    if a.is_infinite() && a.is_sign_negative() {
        return b;
    }
    if b.is_infinite() && b.is_sign_negative() {
        return a;
    }
    let m = a.max(b);
    m + ((10f64.powf(a - m)) + (10f64.powf(b - m))).log10()
}

#[inline]
fn log10_sum3(a: f64, b: f64, c: f64) -> f64 {
    log10_sum(log10_sum(a, b), c)
}

#[inline]
fn effective_qual(baseq: u8, mapq: u8) -> u8 {
    baseq.min(mapq)
}

/// Precomputed Phred emission log10 tables (match / mismatch) for quals 0..=93.
fn emission_tables() -> &'static ([f64; 94], [f64; 94], f64) {
    use std::sync::OnceLock;
    static TABLES: OnceLock<([f64; 94], [f64; 94], f64)> = OnceLock::new();
    TABLES.get_or_init(|| {
        let mut match_tbl = [0.0; 94];
        let mut mismatch_tbl = [0.0; 94];
        for q in 0..94u8 {
            let error_prob = 10f64.powf(-(q as f64) / 10.0);
            match_tbl[q as usize] = prob_to_log10(1.0 - error_prob);
            mismatch_tbl[q as usize] = prob_to_log10(error_prob / 3.0);
        }
        (match_tbl, mismatch_tbl, prob_to_log10(0.25))
    })
}

#[inline]
fn emission_log10(read_base: u8, hap_base: u8, effective_qual: u8) -> f64 {
    // Treat ambiguous base 'N' as non-informative (uniform over A/C/G/T) to
    // better align with Java PairHMM behavior on artifact-like evidence.
    let (match_tbl, mismatch_tbl, n_emit) = emission_tables();
    if read_base == b'N' || hap_base == b'N' {
        return *n_emit;
    }
    let q = (effective_qual as usize).min(93);
    if read_base == hap_base {
        match_tbl[q]
    } else {
        mismatch_tbl[q]
    }
}

/// Borrowed-byte PairHMM kernel (no `String` / `Vec` clones on the hot path).
pub fn pairhmm_log10_likelihood_slices(
    read_bases: &[u8],
    read_base_quals: &[u8],
    read_mapping_quality: u8,
    haplotype_bases: &[u8],
    params: &PairHmmParams,
) -> GatkResult<f64> {
    if read_bases.len() != read_base_quals.len() {
        return Err(GatkError::argument(
            "PairHMM input mismatch: read bases length must match base quality length",
        ));
    }
    if haplotype_bases.is_empty() {
        return Err(GatkError::argument(
            "PairHMM input invalid: haplotype_bases must be non-empty",
        ));
    }
    if !(0.0..1.0).contains(&params.gap_open_prob)
        || !(0.0..1.0).contains(&params.gap_extend_prob)
        || !(0.0..=1.0).contains(&params.insertion_emission_prob)
    {
        return Err(GatkError::argument(
            "PairHMM params invalid: probabilities must be in (0,1) (insertion emission in [0,1])",
        ));
    }

    let r = read_bases;
    let h = haplotype_bases;
    let rn = r.len();
    let hn = h.len();

    if rn == 0 {
        return Ok(0.0);
    }

    // Flat row-major DP (one allocation per matrix instead of one Vec per row).
    let ncol = hn + 1;
    let cells = (rn + 1) * ncol;
    let mut m = vec![LOG10_NEG_INF; cells];
    let mut ins = vec![LOG10_NEG_INF; cells];
    let mut del = vec![LOG10_NEG_INF; cells];
    m[0] = 0.0;
    // Java LoglessPairHMM allows free leading deletions with an initial
    // normalization factor of 1 / haplotype_length.
    let init_del = -((hn as f64).log10());
    for j in 1..=hn {
        del[j] = init_del;
    }

    let p_go = prob_to_log10(params.gap_open_prob);
    let p_ge = prob_to_log10(params.gap_extend_prob);
    let p_stay = prob_to_log10((1.0 - 2.0 * params.gap_open_prob).max(1e-12));
    let p_ins_emit = prob_to_log10(params.insertion_emission_prob.max(1e-12));

    for i in 1..=rn {
        let q = effective_qual(read_base_quals[i - 1], read_mapping_quality);
        let prev = (i - 1) * ncol;
        let cur = i * ncol;
        for j in 1..=hn {
            let e = emission_log10(r[i - 1], h[j - 1], q);
            let from_match = m[prev + j - 1] + p_stay;
            let from_ins = ins[prev + j - 1] + p_go;
            let from_del = del[prev + j - 1] + p_go;
            m[cur + j] = log10_sum3(from_match, from_ins, from_del) + e;

            let ins_from_m = m[prev + j] + p_go;
            let ins_from_i = ins[prev + j] + p_ge;
            ins[cur + j] = log10_sum(ins_from_m, ins_from_i) + p_ins_emit;

            let del_from_m = m[cur + j - 1] + p_go;
            let del_from_d = del[cur + j - 1] + p_ge;
            del[cur + j] = log10_sum(del_from_m, del_from_d);
        }
    }

    // Java PairHMM sums over terminal Match+Insertion states across the full
    // last read row (free trailing deletions in haplotype).
    let mut terminal = LOG10_NEG_INF;
    let last = rn * ncol;
    for j in 1..=hn {
        terminal = log10_sum(terminal, m[last + j]);
        terminal = log10_sum(terminal, ins[last + j]);
    }
    Ok(terminal)
}

/// Compute log10 P(read | haplotype) with a deterministic scalar PairHMM.
pub fn pairhmm_log10_likelihood(input: &PairHmmInput, params: &PairHmmParams) -> GatkResult<f64> {
    input.validate()?;
    pairhmm_log10_likelihood_slices(
        input.read_bases.as_bytes(),
        &input.read_base_quals,
        input.read_mapping_quality,
        input.haplotype_bases.as_bytes(),
        params,
    )
}

/// Vectorized/batched fast path over multiple haplotypes for the same read evidence (`&[u8]`).
/// Parallelizes across haplotypes and reuses the borrowed-byte scalar kernel (no `String` clones).
pub fn pairhmm_log10_likelihoods_vectorized_slices(
    read_bases: &[u8],
    read_base_quals: &[u8],
    read_mapping_quality: u8,
    haplotype_bases: &[&[u8]],
    params: &PairHmmParams,
) -> GatkResult<Vec<f64>> {
    if read_bases.len() != read_base_quals.len() {
        return Err(GatkError::argument(
            "PairHMM vectorized input mismatch: read bases length must match base quality length",
        ));
    }
    if haplotype_bases.is_empty() {
        return Ok(Vec::new());
    }
    haplotype_bases
        .par_iter()
        .map(|hap| {
            pairhmm_log10_likelihood_slices(
                read_bases,
                read_base_quals,
                read_mapping_quality,
                hap,
                params,
            )
        })
        .collect()
}

/// Vectorized path accepting `String` haplotypes (thin wrapper over slice API).
pub fn pairhmm_log10_likelihoods_vectorized(
    read_bases: &str,
    read_base_quals: &[u8],
    read_mapping_quality: u8,
    haplotype_bases: &[String],
    params: &PairHmmParams,
) -> GatkResult<Vec<f64>> {
    let hap_refs: Vec<&[u8]> = haplotype_bases.iter().map(|h| h.as_bytes()).collect();
    pairhmm_log10_likelihoods_vectorized_slices(
        read_bases.as_bytes(),
        read_base_quals,
        read_mapping_quality,
        &hap_refs,
        params,
    )
}

/// Deterministic floating-point equivalence rule for PairHMM comparisons (Phase 6, step 86).
pub fn pairhmm_fp_eq(a: f64, b: f64, policy: PairHmmFpPolicy) -> bool {
    let diff = (a - b).abs();
    if diff <= policy.abs_epsilon {
        return true;
    }
    let scale = a.abs().max(b.abs()).max(1.0);
    diff <= policy.rel_epsilon * scale
}

#[cfg(test)]
mod tests {
    use super::*;

    fn mk_input(read: &str, bq: u8, mapq: u8, hap: &str) -> PairHmmInput {
        PairHmmInput {
            read_bases: read.to_string(),
            read_base_quals: vec![bq; read.len()],
            read_mapping_quality: mapq,
            haplotype_bases: hap.to_string(),
        }
    }

    #[test]
    fn perfect_match_scores_higher_than_mismatch() {
        let params = PairHmmParams::default();
        let match_ll =
            pairhmm_log10_likelihood(&mk_input("ACGT", 35, 60, "ACGT"), &params).unwrap();
        let mismatch_ll =
            pairhmm_log10_likelihood(&mk_input("ACGT", 35, 60, "AGGT"), &params).unwrap();
        assert!(match_ll > mismatch_ll);
    }

    #[test]
    fn mapping_quality_caps_confidence() {
        let params = PairHmmParams::default();
        let high_mapq =
            pairhmm_log10_likelihood(&mk_input("ACGT", 35, 60, "ACGT"), &params).unwrap();
        let low_mapq = pairhmm_log10_likelihood(&mk_input("ACGT", 35, 8, "ACGT"), &params).unwrap();
        assert!(high_mapq > low_mapq);
    }

    #[test]
    fn long_read_stays_finite() {
        let params = PairHmmParams::default();
        let read = "ACGT".repeat(40);
        let hap = "ACGT".repeat(40);
        let ll = pairhmm_log10_likelihood(&mk_input(&read, 30, 60, &hap), &params).unwrap();
        assert!(ll.is_finite());
    }

    #[test]
    fn vectorized_path_matches_scalar() {
        let params = PairHmmParams::default();
        let read = "ACGTACGT";
        let quals = vec![30; read.len()];
        let haps = vec![
            "ACGTACGT".to_string(),
            "ACGTTCGT".to_string(),
            "ACGTACGA".to_string(),
        ];
        let vec_ll =
            pairhmm_log10_likelihoods_vectorized(read, &quals, 60, &haps, &params).unwrap();
        assert_eq!(vec_ll.len(), haps.len());
        for (idx, hap) in haps.iter().enumerate() {
            let scalar = pairhmm_log10_likelihood(
                &PairHmmInput {
                    read_bases: read.to_string(),
                    read_base_quals: quals.clone(),
                    read_mapping_quality: 60,
                    haplotype_bases: hap.clone(),
                },
                &params,
            )
            .unwrap();
            assert!((scalar - vec_ll[idx]).abs() < 1e-12);
        }
    }

    #[test]
    fn n_base_is_non_informative() {
        let params = PairHmmParams::default();
        let with_n = pairhmm_log10_likelihood(&mk_input("ANAA", 30, 60, "ACAA"), &params).unwrap();
        let strict_mismatch =
            pairhmm_log10_likelihood(&mk_input("AGAA", 30, 60, "ACAA"), &params).unwrap();
        // 'N' should be penalized less harshly than a strict mismatch at the same site.
        assert!(with_n > strict_mismatch);
    }

    #[test]
    fn fp_policy_accepts_tiny_numeric_drift() {
        let p = PairHmmFpPolicy::default();
        assert!(pairhmm_fp_eq(-3.0, -3.0 + 1e-11, p));
        assert!(!pairhmm_fp_eq(-3.0, -2.9, p));
    }
}
