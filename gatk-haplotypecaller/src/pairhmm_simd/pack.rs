//! Portable pairs-in-lanes Logless PairHMM (explicit lanes; used as fallback and
//! for uneven packs). AVX2/NEON specialize full-width packs.
//!
//! Scratch planes are sized once to `max(hap_len)` and reused across haplotypes for
//! cache locality (phenotype: many haps × shared read). Numerics match scalar Logless.

use crate::pairhmm_logless::{
    logless_match_mismatch_prior, logless_pairhmm_likelihood, logless_qual_to_trans_probs,
    INITIAL_CONDITION, INITIAL_CONDITION_LOG10, MIN_ACCEPTED_LINEAR_SUM,
};
use gatk_common::{GatkError, GatkResult};
use std::cell::RefCell;

const MATCH_TO_MATCH: usize = 0;
const INDEL_TO_MATCH: usize = 1;
const MATCH_TO_INSERTION: usize = 2;
const INSERTION_TO_INSERTION: usize = 3;
const MATCH_TO_DELETION: usize = 4;
const DELETION_TO_DELETION: usize = 5;

thread_local! {
    static PACK_F64_SCRATCH: RefCell<F64Scratch> = RefCell::new(F64Scratch::empty());
    static PACK_F32_SCRATCH: RefCell<F32Scratch> = RefCell::new(F32Scratch::empty());
}

struct F64Scratch {
    m: Vec<f64>,
    ins: Vec<f64>,
    del: Vec<f64>,
    prior: Vec<f64>,
}

impl F64Scratch {
    fn empty() -> Self {
        Self {
            m: Vec::new(),
            ins: Vec::new(),
            del: Vec::new(),
            prior: Vec::new(),
        }
    }

    fn ensure_cells(&mut self, cells: usize) {
        if self.m.len() < cells {
            self.m.resize(cells, 0.0);
            self.ins.resize(cells, 0.0);
            self.del.resize(cells, 0.0);
            self.prior.resize(cells, 0.0);
        }
    }

    fn clear_prefix(&mut self, cells: usize) {
        self.m[..cells].fill(0.0);
        self.ins[..cells].fill(0.0);
        self.del[..cells].fill(0.0);
        self.prior[..cells].fill(0.0);
    }
}

/// Score one read against many haplotypes with a portable packed f64 kernel.
/// Haplotypes share one DP scratch sized to the longest hap (cache locality).
pub fn score_haps_logless_packed_f64(
    read_bases: &[u8],
    read_quals: &[u8],
    haplotypes: &[&[u8]],
    insertion_gop: &[u8],
    deletion_gop: &[u8],
    overall_gcp: &[u8],
) -> GatkResult<Vec<f64>> {
    if haplotypes.is_empty() {
        return Ok(Vec::new());
    }
    let rn = read_bases.len();
    if rn == 0 {
        return Ok(vec![0.0; haplotypes.len()]);
    }

    let mut transitions = vec![[0.0f64; 6]; rn + 1];
    for i in 0..rn {
        transitions[i + 1] =
            logless_qual_to_trans_probs(insertion_gop[i], deletion_gop[i], overall_gcp[i]);
    }

    let max_hn = haplotypes.iter().map(|h| h.len()).max().unwrap_or(0);
    let max_cols = max_hn + 1;
    let max_cells = (rn + 1).saturating_mul(max_cols);
    // Match scalar PairHMM fail-closed caps (Peak-RSS on 16 GiB hosts).
    const MAX_PAIRHMM_DIM: usize = 100_000;
    const MAX_PAIRHMM_CELLS: usize = 8_000_000;
    if rn > MAX_PAIRHMM_DIM || max_hn > MAX_PAIRHMM_DIM || max_cells > MAX_PAIRHMM_CELLS {
        return Err(GatkError::algorithm(format!(
            "PairHMM packed-f64 refused oversized DP (read_len={rn}, max_hap_len={max_hn}, cells={max_cells}); \
             inputs must be assembly-region scale, not contig scale"
        )));
    }

    let mut out = Vec::with_capacity(haplotypes.len());
    PACK_F64_SCRATCH.with(|cell| {
        let mut scratch = cell.borrow_mut();
        scratch.ensure_cells(max_cells);
        for &hap in haplotypes {
            out.push(score_one_f64(
                read_bases,
                read_quals,
                hap,
                &transitions,
                &mut scratch,
            ));
        }
    });
    Ok(out)
}

/// Score haplotypes with already-built Logless transitions (avoids rebuild on NEON leftovers).
pub(crate) fn score_haps_logless_packed_f64_with_transitions(
    read_bases: &[u8],
    read_quals: &[u8],
    haplotypes: &[&[u8]],
    transitions: &[[f64; 6]],
) -> GatkResult<Vec<f64>> {
    if haplotypes.is_empty() {
        return Ok(Vec::new());
    }
    let rn = read_bases.len();
    if rn == 0 {
        return Ok(vec![0.0; haplotypes.len()]);
    }
    let max_hn = haplotypes.iter().map(|h| h.len()).max().unwrap_or(0);
    let max_cols = max_hn + 1;
    let max_cells = (rn + 1).saturating_mul(max_cols);
    const MAX_PAIRHMM_DIM: usize = 100_000;
    const MAX_PAIRHMM_CELLS: usize = 8_000_000;
    if rn > MAX_PAIRHMM_DIM || max_hn > MAX_PAIRHMM_DIM || max_cells > MAX_PAIRHMM_CELLS {
        return Err(GatkError::algorithm(format!(
            "PairHMM packed-f64 refused oversized DP (read_len={rn}, max_hap_len={max_hn}, cells={max_cells}); \
             inputs must be assembly-region scale, not contig scale"
        )));
    }
    let mut out = Vec::with_capacity(haplotypes.len());
    PACK_F64_SCRATCH.with(|cell| {
        let mut scratch = cell.borrow_mut();
        scratch.ensure_cells(max_cells);
        for &hap in haplotypes {
            out.push(score_one_f64(
                read_bases,
                read_quals,
                hap,
                transitions,
                &mut scratch,
            ));
        }
    });
    Ok(out)
}

/// Score one hap with prebuilt transition planes (shared across a read's hap pack).
pub(crate) fn score_one_f64(
    read_bases: &[u8],
    read_quals: &[u8],
    hap: &[u8],
    transitions: &[[f64; 6]],
    scratch: &mut F64Scratch,
) -> f64 {
    let rn = read_bases.len();
    let hn = hap.len();
    let cols = hn + 1;
    let cells = (rn + 1) * cols;
    scratch.clear_prefix(cells);
    let m = &mut scratch.m;
    let ins = &mut scratch.ins;
    let del = &mut scratch.del;
    let prior = &mut scratch.prior;

    for i in 0..rn {
        let x = read_bases[i];
        let (match_p, mismatch_p) = logless_match_mismatch_prior(read_quals[i]);
        let row = (i + 1) * cols;
        for j in 0..hn {
            let y = hap[j];
            prior[row + j + 1] = if x == y || x == b'N' || y == b'N' {
                match_p
            } else {
                mismatch_p
            };
        }
    }

    let init_del = INITIAL_CONDITION / hn as f64;
    for j in 0..=hn {
        del[j] = init_del;
    }

    for i in 1..=rn {
        let t = transitions[i];
        let row = i * cols;
        let prev = (i - 1) * cols;
        for j in 1..=hn {
            let p = prior[row + j];
            m[row + j] = p
                * (m[prev + j - 1] * t[MATCH_TO_MATCH]
                    + ins[prev + j - 1] * t[INDEL_TO_MATCH]
                    + del[prev + j - 1] * t[INDEL_TO_MATCH]);
            ins[row + j] =
                m[prev + j] * t[MATCH_TO_INSERTION] + ins[prev + j] * t[INSERTION_TO_INSERTION];
            del[row + j] =
                m[row + j - 1] * t[MATCH_TO_DELETION] + del[row + j - 1] * t[DELETION_TO_DELETION];
        }
    }

    let end_row = rn * cols;
    let mut final_sum = 0.0;
    for j in 1..=hn {
        final_sum += m[end_row + j] + ins[end_row + j];
    }
    if final_sum <= 0.0 || !final_sum.is_finite() {
        return f64::NEG_INFINITY;
    }
    final_sum.log10() - INITIAL_CONDITION_LOG10
}

struct F32Scratch {
    m: Vec<f32>,
    ins: Vec<f32>,
    del: Vec<f32>,
    prior: Vec<f32>,
}

impl F32Scratch {
    fn empty() -> Self {
        Self {
            m: Vec::new(),
            ins: Vec::new(),
            del: Vec::new(),
            prior: Vec::new(),
        }
    }

    fn ensure_cells(&mut self, cells: usize) {
        if self.m.len() < cells {
            self.m.resize(cells, 0.0);
            self.ins.resize(cells, 0.0);
            self.del.resize(cells, 0.0);
            self.prior.resize(cells, 0.0);
        }
    }

    fn clear_prefix(&mut self, cells: usize) {
        self.ensure_cells(cells);
        self.m[..cells].fill(0.0);
        self.ins[..cells].fill(0.0);
        self.del[..cells].fill(0.0);
        self.prior[..cells].fill(0.0);
    }
}

/// f32 packed path with per-haplotype f64 retry when the linear sum underflows.
pub fn score_haps_logless_packed_f32(
    read_bases: &[u8],
    read_quals: &[u8],
    haplotypes: &[&[u8]],
    insertion_gop: &[u8],
    deletion_gop: &[u8],
    overall_gcp: &[u8],
) -> GatkResult<Vec<f64>> {
    if haplotypes.is_empty() {
        return Ok(Vec::new());
    }
    let rn = read_bases.len();
    if rn == 0 {
        return Ok(vec![0.0; haplotypes.len()]);
    }

    let mut transitions_f32 = vec![[0.0f32; 6]; rn + 1];
    for i in 0..rn {
        let t = logless_qual_to_trans_probs(insertion_gop[i], deletion_gop[i], overall_gcp[i]);
        transitions_f32[i + 1] = [
            t[0] as f32,
            t[1] as f32,
            t[2] as f32,
            t[3] as f32,
            t[4] as f32,
            t[5] as f32,
        ];
    }

    let max_hn = haplotypes.iter().map(|h| h.len()).max().unwrap_or(0);
    let max_cells = (rn + 1).saturating_mul(max_hn + 1);
    const MAX_PAIRHMM_DIM: usize = 100_000;
    const MAX_PAIRHMM_CELLS: usize = 8_000_000;
    if rn > MAX_PAIRHMM_DIM || max_hn > MAX_PAIRHMM_DIM || max_cells > MAX_PAIRHMM_CELLS {
        return Err(GatkError::algorithm(format!(
            "PairHMM packed-f32 refused oversized DP (read_len={rn}, max_hap_len={max_hn}, cells={max_cells}); \
             inputs must be assembly-region scale, not contig scale"
        )));
    }
    let mut out = Vec::with_capacity(haplotypes.len());
    let mut err = None;
    PACK_F32_SCRATCH.with(|cell| {
        let mut scratch = cell.borrow_mut();
        scratch.ensure_cells(max_cells);
        for &hap in haplotypes {
            let (ll, linear_sum) =
                score_one_f32(read_bases, read_quals, hap, &transitions_f32, &mut scratch);
            if !linear_sum.is_finite() || (linear_sum as f64) < MIN_ACCEPTED_LINEAR_SUM {
                match logless_pairhmm_likelihood(
                    read_bases,
                    read_quals,
                    hap,
                    insertion_gop,
                    deletion_gop,
                    overall_gcp,
                ) {
                    Ok(v) => out.push(v),
                    Err(e) => {
                        err = Some(e);
                        break;
                    }
                }
            } else {
                out.push(ll);
            }
        }
    });
    if let Some(e) = err {
        return Err(e);
    }
    Ok(out)
}

fn score_one_f32(
    read_bases: &[u8],
    read_quals: &[u8],
    hap: &[u8],
    transitions: &[[f32; 6]],
    scratch: &mut F32Scratch,
) -> (f64, f32) {
    let rn = read_bases.len();
    let hn = hap.len();
    let cols = hn + 1;
    let cells = (rn + 1) * cols;
    scratch.clear_prefix(cells);
    let m = &mut scratch.m;
    let ins = &mut scratch.ins;
    let del = &mut scratch.del;
    let prior = &mut scratch.prior;

    for i in 0..rn {
        let x = read_bases[i];
        let (match_p, mismatch_p) = logless_match_mismatch_prior(read_quals[i]);
        let match_p = match_p as f32;
        let mismatch_p = mismatch_p as f32;
        let row = (i + 1) * cols;
        for j in 0..hn {
            let y = hap[j];
            prior[row + j + 1] = if x == y || x == b'N' || y == b'N' {
                match_p
            } else {
                mismatch_p
            };
        }
    }

    let init_del = (INITIAL_CONDITION / hn as f64) as f32;
    for j in 0..=hn {
        del[j] = init_del;
    }

    for i in 1..=rn {
        let t = transitions[i];
        let row = i * cols;
        let prev = (i - 1) * cols;
        for j in 1..=hn {
            let p = prior[row + j];
            m[row + j] = p
                * (m[prev + j - 1] * t[MATCH_TO_MATCH]
                    + ins[prev + j - 1] * t[INDEL_TO_MATCH]
                    + del[prev + j - 1] * t[INDEL_TO_MATCH]);
            ins[row + j] =
                m[prev + j] * t[MATCH_TO_INSERTION] + ins[prev + j] * t[INSERTION_TO_INSERTION];
            del[row + j] =
                m[row + j - 1] * t[MATCH_TO_DELETION] + del[row + j - 1] * t[DELETION_TO_DELETION];
        }
    }

    let end_row = rn * cols;
    let mut final_sum = 0.0f32;
    for j in 1..=hn {
        final_sum += m[end_row + j] + ins[end_row + j];
    }
    if final_sum <= 0.0 || !final_sum.is_finite() {
        return (f64::NEG_INFINITY, final_sum);
    }
    let ll = (final_sum as f64).log10() - INITIAL_CONDITION_LOG10;
    (ll, final_sum)
}
