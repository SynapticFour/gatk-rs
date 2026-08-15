//! aarch64 NEON pairs-in-lanes Logless PairHMM (2×f64).

#![cfg(target_arch = "aarch64")]

use super::pack::{score_haps_logless_packed_f64, score_haps_logless_packed_f64_with_transitions};
use crate::pairhmm_logless::{
    logless_build_transitions, logless_match_mismatch_prior, INITIAL_CONDITION,
    INITIAL_CONDITION_LOG10,
};
use gatk_common::GatkResult;
use std::arch::aarch64::*;
use std::cell::RefCell;

const MATCH_TO_MATCH: usize = 0;
const INDEL_TO_MATCH: usize = 1;
const MATCH_TO_INSERTION: usize = 2;
const INSERTION_TO_INSERTION: usize = 3;
const MATCH_TO_DELETION: usize = 4;
const DELETION_TO_DELETION: usize = 5;

const LANES: usize = 2;

struct NeonScratch {
    m: Vec<f64>,
    ins: Vec<f64>,
    del: Vec<f64>,
    prior: Vec<f64>,
}

impl NeonScratch {
    fn new() -> Self {
        Self {
            m: Vec::new(),
            ins: Vec::new(),
            del: Vec::new(),
            prior: Vec::new(),
        }
    }

    fn ensure(&mut self, cells: usize) {
        let need = cells * LANES;
        if self.m.len() < need {
            self.m.resize(need, 0.0);
            self.ins.resize(need, 0.0);
            self.del.resize(need, 0.0);
            self.prior.resize(need, 0.0);
        }
        // Always clear the active prefix. On grow, `resize` only zeroes the new
        // tail — leaving stale row-0 M/I from a prior smaller pack, which corrupts DP.
        self.m[..need].fill(0.0);
        self.ins[..need].fill(0.0);
        self.del[..need].fill(0.0);
        self.prior[..need].fill(0.0);
    }
}

thread_local! {
    static NEON_SCRATCH: RefCell<NeonScratch> = RefCell::new(NeonScratch::new());
}

/// Drop NEON PairHMM TLS planes (Peak hygiene after a region).
pub fn release_pairhmm_neon_tls_scratch() {
    NEON_SCRATCH.with(|c| {
        let mut s = c.borrow_mut();
        s.m = Vec::new();
        s.ins = Vec::new();
        s.del = Vec::new();
        s.prior = Vec::new();
    });
}

/// Score haplotypes with NEON when available; otherwise portable packed f64.
pub fn score_haps_neon_f64(
    read_bases: &[u8],
    read_quals: &[u8],
    haplotypes: &[&[u8]],
    insertion_gop: &[u8],
    deletion_gop: &[u8],
    overall_gcp: &[u8],
) -> GatkResult<Vec<f64>> {
    // Equal-length look-ahead packs (no full sort). Length-grouping sort was measured
    // slower on dense NA12878 (sort/scatter tax > pack win under SEQUENTIAL=1); keep TLS ensure fix.
    if std::arch::is_aarch64_feature_detected!("neon") {
        // SAFETY: NEON feature detected.
        unsafe {
            return score_haps_neon_f64_unchecked(
                read_bases,
                read_quals,
                haplotypes,
                insertion_gop,
                deletion_gop,
                overall_gcp,
            );
        }
    }
    score_haps_logless_packed_f64(
        read_bases,
        read_quals,
        haplotypes,
        insertion_gop,
        deletion_gop,
        overall_gcp,
    )
}

#[target_feature(enable = "neon")]
unsafe fn score_haps_neon_f64_unchecked(
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

    let transitions = logless_build_transitions(rn, insertion_gop, deletion_gop, overall_gcp);
    let mut out = vec![0.0f64; haplotypes.len()];
    let mut done = vec![false; haplotypes.len()];
    let mut err: Option<gatk_common::GatkError> = None;
    NEON_SCRATCH.with(|cell| {
        let mut scratch = cell.borrow_mut();
        let mut i = 0;
        while i < haplotypes.len() {
            if done[i] {
                i += 1;
                continue;
            }
            // Look ahead for any remaining equal-length partner (no full sort/scatter).
            // Uneven lengths in one NEON pack corrupt lane DP — never pack mismatched lens.
            let partner = ((i + 1)..haplotypes.len())
                .find(|&j| !done[j] && haplotypes[j].len() == haplotypes[i].len());
            if let Some(j) = partner {
                let pack = [haplotypes[i], haplotypes[j]];
                let scores = score_pack2(read_bases, read_quals, &pack, &transitions, &mut scratch);
                out[i] = scores[0];
                out[j] = scores[1];
                done[i] = true;
                done[j] = true;
            } else {
                // Reuse transitions already built for this read — do not rebuild via packed_f64.
                match score_haps_logless_packed_f64_with_transitions(
                    read_bases,
                    read_quals,
                    &haplotypes[i..=i],
                    &transitions,
                ) {
                    Ok(rest) => out[i] = rest[0],
                    Err(e) => {
                        err = Some(e);
                        break;
                    }
                }
                done[i] = true;
            }
            i += 1;
        }
    });
    if let Some(e) = err {
        return Err(e);
    }
    Ok(out)
}

#[target_feature(enable = "neon")]
unsafe fn score_pack2(
    read_bases: &[u8],
    read_quals: &[u8],
    haps: &[&[u8]; 2],
    transitions: &[[f64; 6]],
    scratch: &mut NeonScratch,
) -> [f64; 2] {
    // Equal-length packs only (caller gates). TLS scratch must clear on grow
    // (`NeonScratch::ensure`) — stale row-0 M/I was the equal-length ≠ scalar bug.
    let rn = read_bases.len();
    let hn = haps[0].len();
    debug_assert_eq!(haps[0].len(), haps[1].len());
    let cols = hn + 1;
    let cells = (rn + 1) * cols;
    scratch.ensure(cells);
    let m = &mut scratch.m;
    let ins = &mut scratch.ins;
    let del = &mut scratch.del;
    let prior = &mut scratch.prior;

    let init = INITIAL_CONDITION / hn as f64;
    for j in 0..=hn {
        let base = j * LANES;
        del[base] = init;
        del[base + 1] = init;
    }

    for i in 0..rn {
        let x = read_bases[i];
        let (match_p, mismatch_p) = logless_match_mismatch_prior(read_quals[i]);
        let row = (i + 1) * cols;
        for j in 0..hn {
            let idx = (row + j + 1) * LANES;
            for lane in 0..LANES {
                let y = haps[lane][j];
                prior[idx + lane] = if x == y || x == b'N' || y == b'N' {
                    match_p
                } else {
                    mismatch_p
                };
            }
        }
    }

    for i in 1..=rn {
        let t = transitions[i];
        let tm2m = vdupq_n_f64(t[MATCH_TO_MATCH]);
        let ti2m = vdupq_n_f64(t[INDEL_TO_MATCH]);
        let tm2i = vdupq_n_f64(t[MATCH_TO_INSERTION]);
        let ti2i = vdupq_n_f64(t[INSERTION_TO_INSERTION]);
        let tm2d = vdupq_n_f64(t[MATCH_TO_DELETION]);
        let td2d = vdupq_n_f64(t[DELETION_TO_DELETION]);
        let row = i * cols;
        let prev = (i - 1) * cols;
        for j in 1..=hn {
            let idx = (row + j) * LANES;
            let diag = (prev + j - 1) * LANES;
            let up = (prev + j) * LANES;
            let left = (row + j - 1) * LANES;

            // SAFETY: SoA buffers sized for cells*LANES; NEON enabled; equal-length pack.
            let m_diag = vld1q_f64(m.as_ptr().add(diag));
            let i_diag = vld1q_f64(ins.as_ptr().add(diag));
            let d_diag = vld1q_f64(del.as_ptr().add(diag));
            let m_up = vld1q_f64(m.as_ptr().add(up));
            let i_up = vld1q_f64(ins.as_ptr().add(up));
            let m_left = vld1q_f64(m.as_ptr().add(left));
            let d_left = vld1q_f64(del.as_ptr().add(left));
            let p = vld1q_f64(prior.as_ptr().add(idx));

            let sum = vaddq_f64(
                vaddq_f64(vmulq_f64(m_diag, tm2m), vmulq_f64(i_diag, ti2m)),
                vmulq_f64(d_diag, ti2m),
            );
            let m_new = vmulq_f64(p, sum);
            let i_new = vaddq_f64(vmulq_f64(m_up, tm2i), vmulq_f64(i_up, ti2i));
            // del[i][j] uses match[i][j-1] already stored at previous j (m_left).
            let d_new = vaddq_f64(vmulq_f64(m_left, tm2d), vmulq_f64(d_left, td2d));

            vst1q_f64(m.as_mut_ptr().add(idx), m_new);
            vst1q_f64(ins.as_mut_ptr().add(idx), i_new);
            vst1q_f64(del.as_mut_ptr().add(idx), d_new);
        }
    }

    let end_row = rn * cols;
    let mut sums = [0.0f64; 2];
    for j in 1..=hn {
        let idx = (end_row + j) * LANES;
        let mv = vld1q_f64(m.as_ptr().add(idx));
        let iv = vld1q_f64(ins.as_ptr().add(idx));
        let s = vaddq_f64(mv, iv);
        let mut tmp = [0.0f64; 2];
        vst1q_f64(tmp.as_mut_ptr(), s);
        sums[0] += tmp[0];
        sums[1] += tmp[1];
    }

    let mut out = [0.0f64; 2];
    for lane in 0..LANES {
        if sums[lane] <= 0.0 || !sums[lane].is_finite() {
            out[lane] = f64::NEG_INFINITY;
        } else {
            out[lane] = sums[lane].log10() - INITIAL_CONDITION_LOG10;
        }
    }
    out
}
