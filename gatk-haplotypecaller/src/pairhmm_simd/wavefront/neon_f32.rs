//! NEON f32 row-wavefront: 4-wide Match/Insertion striping; Deletion serial.

#![allow(clippy::too_many_arguments)]

use super::prep::ReadPrep;
use super::rolling_f32::INITIAL_CONDITION_F32;
use super::WavefrontScratch;
use std::arch::aarch64::*;

const LANES: usize = 4;
const MATCH_TO_MATCH: usize = 0;
const INDEL_TO_MATCH: usize = 1;
const MATCH_TO_INSERTION: usize = 2;
const INSERTION_TO_INSERTION: usize = 3;
const MATCH_TO_DELETION: usize = 4;
const DELETION_TO_DELETION: usize = 5;

/// Returns linear-space final sum (caller applies log10 / retry).
///
/// SAFETY: caller must ensure NEON is available (always true on aarch64 targets we build for).
#[target_feature(enable = "neon")]
pub unsafe fn score_one_neon_f32(
    prep: &ReadPrep,
    read_bases: &[u8],
    hap: &[u8],
    scratch: &mut WavefrontScratch,
) -> f32 {
    let rn = read_bases.len();
    let hn = hap.len();
    if rn == 0 {
        return 1.0;
    }
    let cols = hn + 1;
    scratch.ensure_rolling_f32(cols);
    let (m_a, m_b) = scratch.m32[..cols * 2].split_at_mut(cols);
    let (ins_a, ins_b) = scratch.ins32[..cols * 2].split_at_mut(cols);
    let (del_a, del_b) = scratch.del32[..cols * 2].split_at_mut(cols);

    m_a.fill(0.0);
    ins_a.fill(0.0);
    let init_del = INITIAL_CONDITION_F32 / hn as f32;
    for j in 0..=hn {
        del_a[j] = init_del;
    }

    let mut prev_is_a = true;
    for i in 1..=rn {
        let t = prep.transitions_f32[i];
        let x = read_bases[i - 1];
        let (match_p, mismatch_p) = prep.match_mm_f32[i - 1];
        let (m_prev, m_curr, ins_prev, ins_curr, del_prev, del_curr) = if prev_is_a {
            (
                &m_a[..],
                &mut m_b[..],
                &ins_a[..],
                &mut ins_b[..],
                &del_a[..],
                &mut del_b[..],
            )
        } else {
            (
                &m_b[..],
                &mut m_a[..],
                &ins_b[..],
                &mut ins_a[..],
                &del_b[..],
                &mut del_a[..],
            )
        };
        m_curr[0] = 0.0;
        ins_curr[0] = 0.0;
        del_curr[0] = 0.0;
        fill_row_neon_f32(
            hap, x, match_p, mismatch_p, t, m_prev, m_curr, ins_prev, ins_curr, del_prev, del_curr,
            hn,
        );
        prev_is_a = !prev_is_a;
    }

    let (m_end, ins_end) = if prev_is_a {
        (&m_a[..], &ins_a[..])
    } else {
        (&m_b[..], &ins_b[..])
    };
    let mut final_sum = 0.0f32;
    for j in 1..=hn {
        final_sum += m_end[j] + ins_end[j];
    }
    final_sum
}

#[target_feature(enable = "neon")]
unsafe fn fill_row_neon_f32(
    hap: &[u8],
    x: u8,
    match_p: f32,
    mismatch_p: f32,
    t: [f32; 6],
    m_prev: &[f32],
    m_curr: &mut [f32],
    ins_prev: &[f32],
    ins_curr: &mut [f32],
    del_prev: &[f32],
    del_curr: &mut [f32],
    hn: usize,
) {
    let tm2m = vdupq_n_f32(t[MATCH_TO_MATCH]);
    let ti2m = vdupq_n_f32(t[INDEL_TO_MATCH]);
    let tm2i = vdupq_n_f32(t[MATCH_TO_INSERTION]);
    let ti2i = vdupq_n_f32(t[INSERTION_TO_INSERTION]);
    let tm2d = t[MATCH_TO_DELETION];
    let td2d = t[DELETION_TO_DELETION];

    let mut j = 1usize;
    while j + LANES - 1 <= hn {
        let mut prior = [0.0f32; LANES];
        for lane in 0..LANES {
            let y = hap[j - 1 + lane];
            prior[lane] = if x == y || x == b'N' || y == b'N' {
                match_p
            } else {
                mismatch_p
            };
        }
        let p = vld1q_f32(prior.as_ptr());
        let m_diag = vld1q_f32(m_prev.as_ptr().add(j - 1));
        let i_diag = vld1q_f32(ins_prev.as_ptr().add(j - 1));
        let d_diag = vld1q_f32(del_prev.as_ptr().add(j - 1));
        let m_up = vld1q_f32(m_prev.as_ptr().add(j));
        let i_up = vld1q_f32(ins_prev.as_ptr().add(j));

        let from_m = vmulq_f32(m_diag, tm2m);
        let from_i = vmulq_f32(i_diag, ti2m);
        let from_d = vmulq_f32(d_diag, ti2m);
        let sum = vaddq_f32(vaddq_f32(from_m, from_i), from_d);
        let m_new = vmulq_f32(p, sum);
        let i_new = vaddq_f32(vmulq_f32(m_up, tm2i), vmulq_f32(i_up, ti2i));

        vst1q_f32(m_curr.as_mut_ptr().add(j), m_new);
        vst1q_f32(ins_curr.as_mut_ptr().add(j), i_new);

        for lane in 0..LANES {
            let jj = j + lane;
            del_curr[jj] = m_curr[jj - 1] * tm2d + del_curr[jj - 1] * td2d;
        }
        j += LANES;
    }
    while j <= hn {
        let y = hap[j - 1];
        let p = if x == y || x == b'N' || y == b'N' {
            match_p
        } else {
            mismatch_p
        };
        m_curr[j] = p
            * (m_prev[j - 1] * t[MATCH_TO_MATCH]
                + ins_prev[j - 1] * t[INDEL_TO_MATCH]
                + del_prev[j - 1] * t[INDEL_TO_MATCH]);
        ins_curr[j] = m_prev[j] * t[MATCH_TO_INSERTION] + ins_prev[j] * t[INSERTION_TO_INSERTION];
        del_curr[j] = m_curr[j - 1] * tm2d + del_curr[j - 1] * td2d;
        j += 1;
    }
}
