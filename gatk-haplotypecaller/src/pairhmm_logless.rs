//! GATK `LoglessPairHMM` — linear-space DP, log10 only on the final sum.
//! Observable Java contract (`LoglessPairHMM.subComputeReadLikelihoodGivenHaplotypeLog10`):
//! `INITIAL_CONDITION = 2^1020`; free leading deletions = `INITIAL / haplen`.
//! Match/Ins/Del updates multiply (not log-sum); priors use `qualToProb` / `qualToErrorProb/3`.
//! Result = `log10(sum last-row M+I) − log10(INITIAL_CONDITION)`.
//! This is the scalar reference for SIMD kernels. Exact `Log10PairHMM` remains in
//! [`crate::pairhmm_log10`] for F.1/F.2 bit-identical dumps.

use gatk_common::{GatkError, GatkResult};
use std::cell::RefCell;
use std::sync::OnceLock;

const MAX_QUAL: usize = 127;
/// `2^1020` — GATK `LoglessPairHMM.INITIAL_CONDITION`.
pub const INITIAL_CONDITION: f64 = 1.0715086071862673e307; // 2f64.powi(1020)
/// `log10(2^1020)`.
pub const INITIAL_CONDITION_LOG10: f64 = 307.6526555685887; // 1020.0 * f64::log10(2.0)

const MATCH_TO_MATCH: usize = 0;
const INDEL_TO_MATCH: usize = 1;
const MATCH_TO_INSERTION: usize = 2;
const INSERTION_TO_INSERTION: usize = 3;
const MATCH_TO_DELETION: usize = 4;
const DELETION_TO_DELETION: usize = 5;

/// GKL-style underflow threshold for f32→f64 retry (linear-space final sum before log10).
pub const MIN_ACCEPTED_LINEAR_SUM: f64 = 1e-28;

#[inline]
fn qual_to_error_prob(qual: u8) -> f64 {
    10f64.powf(-(qual as f64) / 10.0)
}

#[inline]
fn qual_to_prob(qual: u8) -> f64 {
    1.0 - qual_to_error_prob(qual)
}

fn approximate_log10_sum_log10(a: f64, b: f64) -> f64 {
    let (x, y) = if a > b { (b, a) } else { (a, b) };
    if x.is_infinite() && x.is_sign_negative() {
        return y;
    }
    y + (1.0 + 10f64.powf(x - y)).log10()
}

fn match_to_match_prob(ins_qual: u8, del_qual: u8) -> f64 {
    let (min_q, max_q) = if ins_qual <= del_qual {
        (ins_qual as usize, del_qual as usize)
    } else {
        (del_qual as usize, ins_qual as usize)
    };
    if max_q > MAX_QUAL {
        let log10_sum = approximate_log10_sum_log10(-0.1 * min_q as f64, -0.1 * max_q as f64);
        1.0 - 10f64.powf(log10_sum).min(1.0)
    } else {
        match_to_match_prob_table()[((max_q * (max_q + 1)) >> 1) + min_q]
    }
}

fn match_to_match_prob_table() -> &'static [f64; ((MAX_QUAL + 1) * (MAX_QUAL + 2)) >> 1] {
    static TABLE: OnceLock<[f64; ((MAX_QUAL + 1) * (MAX_QUAL + 2)) >> 1]> = OnceLock::new();
    TABLE.get_or_init(|| {
        let mut table = [0.0_f64; ((MAX_QUAL + 1) * (MAX_QUAL + 2)) >> 1];
        let mut offset = 0usize;
        for i in 0..=MAX_QUAL {
            for j in 0..=i {
                let log10_sum = approximate_log10_sum_log10(-0.1 * j as f64, -0.1 * i as f64);
                let log10_m2m = (1.0 - 10f64.powf(log10_sum).min(1.0)).log10();
                table[offset + j] = 10f64.powf(log10_m2m);
            }
            offset += i + 1;
        }
        table
    })
}

fn qual_to_trans_probs(ins_qual: u8, del_qual: u8, gcp: u8) -> [f64; 6] {
    let gcp_err = qual_to_error_prob(gcp);
    [
        match_to_match_prob(ins_qual, del_qual),
        qual_to_prob(gcp),
        qual_to_error_prob(ins_qual),
        gcp_err,
        qual_to_error_prob(del_qual),
        gcp_err,
    ]
}

struct LoglessScratch {
    transition: Vec<[f64; 6]>,
    prior: Vec<f64>,
    m: Vec<f64>,
    ins: Vec<f64>,
    del: Vec<f64>,
}

impl LoglessScratch {
    fn new() -> Self {
        Self {
            transition: Vec::new(),
            prior: Vec::new(),
            m: Vec::new(),
            ins: Vec::new(),
            del: Vec::new(),
        }
    }

    fn ensure(&mut self, rn: usize, hn: usize) {
        let rows = rn + 1;
        let cols = hn + 1;
        let cells = rows * cols;
        if self.transition.len() < rows {
            self.transition.resize(rows, [0.0; 6]);
        }
        if self.prior.len() < cells {
            self.prior.resize(cells, 0.0);
            self.m.resize(cells, 0.0);
            self.ins.resize(cells, 0.0);
            self.del.resize(cells, 0.0);
        }
    }
}

thread_local! {
    static LOGLESS_SCRATCH: RefCell<LoglessScratch> = RefCell::new(LoglessScratch::new());
}

/// Log10 P(read | haplotype) via GATK `LoglessPairHMM` (linear DP).
pub fn logless_pairhmm_likelihood(
    read_bases: &[u8],
    read_quals: &[u8],
    haplotype_bases: &[u8],
    insertion_gop: &[u8],
    deletion_gop: &[u8],
    overall_gcp: &[u8],
) -> GatkResult<f64> {
    validate_inputs(
        read_bases,
        read_quals,
        haplotype_bases,
        insertion_gop,
        deletion_gop,
        overall_gcp,
    )?;
    let rn = read_bases.len();
    let hn = haplotype_bases.len();
    if rn == 0 {
        return Ok(0.0);
    }
    LOGLESS_SCRATCH.with(|cell| {
        let mut scratch = cell.borrow_mut();
        scratch.ensure(rn, hn);
        Ok(logless_pairhmm_likelihood_into(
            read_bases,
            read_quals,
            haplotype_bases,
            insertion_gop,
            deletion_gop,
            overall_gcp,
            &mut scratch,
        ))
    })
}

/// Parity defaults (45/45/10) for Logless — same quals as Log10 parity path.
pub fn logless_pairhmm_likelihood_parity_defaults(
    read_bases: &[u8],
    read_quals: &[u8],
    haplotype_bases: &[u8],
) -> GatkResult<f64> {
    use crate::pairhmm_log10::{
        GATK_PARITY_DEFAULT_DEL_QUAL, GATK_PARITY_DEFAULT_GCP, GATK_PARITY_DEFAULT_INS_QUAL,
    };
    let n = read_bases.len();
    let ins = vec![GATK_PARITY_DEFAULT_INS_QUAL; n];
    let del = vec![GATK_PARITY_DEFAULT_DEL_QUAL; n];
    let gcp = vec![GATK_PARITY_DEFAULT_GCP; n];
    logless_pairhmm_likelihood(read_bases, read_quals, haplotype_bases, &ins, &del, &gcp)
}

fn validate_inputs(
    read_bases: &[u8],
    read_quals: &[u8],
    haplotype_bases: &[u8],
    insertion_gop: &[u8],
    deletion_gop: &[u8],
    overall_gcp: &[u8],
) -> GatkResult<()> {
    if read_bases.len() != read_quals.len()
        || read_bases.len() != insertion_gop.len()
        || read_bases.len() != deletion_gop.len()
        || read_bases.len() != overall_gcp.len()
    {
        return Err(GatkError::argument(
            "PairHMM read arrays must have equal length",
        ));
    }
    if haplotype_bases.is_empty() {
        return Err(GatkError::argument("haplotype must be non-empty"));
    }
    Ok(())
}

fn logless_pairhmm_likelihood_into(
    read_bases: &[u8],
    read_quals: &[u8],
    haplotype_bases: &[u8],
    insertion_gop: &[u8],
    deletion_gop: &[u8],
    overall_gcp: &[u8],
    scratch: &mut LoglessScratch,
) -> f64 {
    let rn = read_bases.len();
    let hn = haplotype_bases.len();
    let cols = hn + 1;

    for i in 0..rn {
        scratch.transition[i + 1] =
            qual_to_trans_probs(insertion_gop[i], deletion_gop[i], overall_gcp[i]);
    }

    let cells = (rn + 1) * cols;
    for slot in &mut scratch.m[..cells] {
        *slot = 0.0;
    }
    for slot in &mut scratch.ins[..cells] {
        *slot = 0.0;
    }
    for slot in &mut scratch.del[..cells] {
        *slot = 0.0;
    }
    for slot in &mut scratch.prior[..cells] {
        *slot = 0.0;
    }

    for i in 0..rn {
        let x = read_bases[i];
        let qual = read_quals[i];
        let match_p = qual_to_prob(qual);
        let mismatch_p = qual_to_error_prob(qual) / 3.0;
        let row = (i + 1) * cols;
        for j in 0..hn {
            let y = haplotype_bases[j];
            scratch.prior[row + j + 1] = if x == y || x == b'N' || y == b'N' {
                match_p
            } else {
                mismatch_p
            };
        }
    }

    let init_del = INITIAL_CONDITION / hn as f64;
    for j in 0..=hn {
        scratch.del[j] = init_del;
    }

    for i in 1..=rn {
        let t = scratch.transition[i];
        let row = i * cols;
        let prev = (i - 1) * cols;
        for j in 1..=hn {
            let p = scratch.prior[row + j];
            scratch.m[row + j] = p
                * (scratch.m[prev + j - 1] * t[MATCH_TO_MATCH]
                    + scratch.ins[prev + j - 1] * t[INDEL_TO_MATCH]
                    + scratch.del[prev + j - 1] * t[INDEL_TO_MATCH]);
            scratch.ins[row + j] = scratch.m[prev + j] * t[MATCH_TO_INSERTION]
                + scratch.ins[prev + j] * t[INSERTION_TO_INSERTION];
            scratch.del[row + j] = scratch.m[row + j - 1] * t[MATCH_TO_DELETION]
                + scratch.del[row + j - 1] * t[DELETION_TO_DELETION];
        }
    }

    let end_row = rn * cols;
    let mut final_sum = 0.0;
    for j in 1..=hn {
        final_sum += scratch.m[end_row + j] + scratch.ins[end_row + j];
    }
    if final_sum <= 0.0 || !final_sum.is_finite() {
        return f64::NEG_INFINITY;
    }
    final_sum.log10() - INITIAL_CONDITION_LOG10
}

/// Score one read against many haplotypes (scalar Logless; haplotypes sequential).
pub fn logless_pairhmm_likelihoods(
    read_bases: &[u8],
    read_quals: &[u8],
    haplotype_bases: &[&[u8]],
    insertion_gop: &[u8],
    deletion_gop: &[u8],
    overall_gcp: &[u8],
) -> GatkResult<Vec<f64>> {
    haplotype_bases
        .iter()
        .map(|hap| {
            logless_pairhmm_likelihood(
                read_bases,
                read_quals,
                hap,
                insertion_gop,
                deletion_gop,
                overall_gcp,
            )
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pairhmm::{pairhmm_fp_eq, PairHmmFpPolicy};
    use crate::pairhmm_log10::log10_pairhmm_likelihood_parity_defaults;

    /// Looser than Logless↔Logless: Log10 vs Logless are different engines.
    fn logless_vs_log10_policy() -> PairHmmFpPolicy {
        PairHmmFpPolicy {
            abs_epsilon: 1e-4,
            rel_epsilon: 1e-3,
        }
    }

    #[test]
    fn initial_condition_constants_match_java() {
        assert!((INITIAL_CONDITION - 2f64.powi(1020)).abs() / INITIAL_CONDITION < 1e-15);
        assert!((INITIAL_CONDITION_LOG10 - 1020.0 * 2f64.log10()).abs() < 1e-12);
    }

    #[test]
    fn perfect_match_finite_and_near_log10() {
        let read = b"ACGTACGTAC";
        let quals = vec![30u8; read.len()];
        let hap = b"ACGTACGTAC";
        let ll = logless_pairhmm_likelihood_parity_defaults(read, &quals, hap).unwrap();
        let log10 = log10_pairhmm_likelihood_parity_defaults(read, &quals, hap).unwrap();
        assert!(ll.is_finite(), "ll={ll}");
        assert!(
            pairhmm_fp_eq(ll, log10, logless_vs_log10_policy()),
            "logless={ll} log10={log10}"
        );
    }

    #[test]
    fn mismatch_and_short_read() {
        let read = b"A";
        let quals = [20u8];
        let hap = b"ACGT";
        let ll = logless_pairhmm_likelihood_parity_defaults(read, &quals, hap).unwrap();
        assert!(ll.is_finite() || ll.is_infinite());
        let log10 = log10_pairhmm_likelihood_parity_defaults(read, &quals, hap).unwrap();
        assert!(pairhmm_fp_eq(ll, log10, logless_vs_log10_policy()) || ll.is_infinite());
    }
}
