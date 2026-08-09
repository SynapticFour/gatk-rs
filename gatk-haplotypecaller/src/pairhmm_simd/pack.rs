//! Portable pairs-in-lanes Logless PairHMM (explicit lanes; used as fallback and
//! for uneven packs). AVX2/NEON specialize full-width packs.
//!
//! Scratch planes are sized once to `max(hap_len)` and reused across haplotypes for
//! cache locality (phenotype: many haps × shared read). Numerics match scalar Logless.

use crate::pairhmm_logless::{
    logless_pairhmm_likelihood, INITIAL_CONDITION, INITIAL_CONDITION_LOG10, MIN_ACCEPTED_LINEAR_SUM,
};
use gatk_common::{GatkError, GatkResult};

const MATCH_TO_MATCH: usize = 0;
const INDEL_TO_MATCH: usize = 1;
const MATCH_TO_INSERTION: usize = 2;
const INSERTION_TO_INSERTION: usize = 3;
const MATCH_TO_DELETION: usize = 4;
const DELETION_TO_DELETION: usize = 5;

const MAX_QUAL: usize = 127;

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
        // Same closed form as the table builder in pairhmm_logless.
        let log10_sum = approximate_log10_sum_log10(-0.1 * min_q as f64, -0.1 * max_q as f64);
        1.0 - 10f64.powf(log10_sum).min(1.0)
    }
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

struct F64Scratch {
    m: Vec<f64>,
    ins: Vec<f64>,
    del: Vec<f64>,
    prior: Vec<f64>,
}

impl F64Scratch {
    fn with_capacity_cells(cells: usize) -> Self {
        Self {
            m: vec![0.0; cells],
            ins: vec![0.0; cells],
            del: vec![0.0; cells],
            prior: vec![0.0; cells],
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
        transitions[i + 1] = qual_to_trans_probs(insertion_gop[i], deletion_gop[i], overall_gcp[i]);
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
    let mut scratch = F64Scratch::with_capacity_cells(max_cells);

    let mut out = Vec::with_capacity(haplotypes.len());
    for &hap in haplotypes {
        out.push(score_one_f64(
            read_bases,
            read_quals,
            hap,
            &transitions,
            &mut scratch,
        ));
    }
    Ok(out)
}

fn score_one_f64(
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
        let match_p = qual_to_prob(read_quals[i]);
        let mismatch_p = qual_to_error_prob(read_quals[i]) / 3.0;
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
    fn with_capacity_cells(cells: usize) -> Self {
        Self {
            m: vec![0.0; cells],
            ins: vec![0.0; cells],
            del: vec![0.0; cells],
            prior: vec![0.0; cells],
        }
    }

    fn clear_prefix(&mut self, cells: usize) {
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
        let t = qual_to_trans_probs(insertion_gop[i], deletion_gop[i], overall_gcp[i]);
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
    let mut scratch = F32Scratch::with_capacity_cells(max_cells);

    let mut out = Vec::with_capacity(haplotypes.len());
    for &hap in haplotypes {
        let (ll, linear_sum) =
            score_one_f32(read_bases, read_quals, hap, &transitions_f32, &mut scratch);
        if !linear_sum.is_finite() || (linear_sum as f64) < MIN_ACCEPTED_LINEAR_SUM {
            out.push(logless_pairhmm_likelihood(
                read_bases,
                read_quals,
                hap,
                insertion_gop,
                deletion_gop,
                overall_gcp,
            )?);
        } else {
            out.push(ll);
        }
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
        let match_p = qual_to_prob(read_quals[i]) as f32;
        let mismatch_p = (qual_to_error_prob(read_quals[i]) / 3.0) as f32;
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
