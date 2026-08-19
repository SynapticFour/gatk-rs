//! AVX2 f32 row-wavefront: 8-wide Match/Insertion striping; Deletion serial.

#![allow(clippy::too_many_arguments)]

use super::prep::ReadPrep;
use super::rolling_f32::INITIAL_CONDITION_F32;
use super::WavefrontScratch;
#[cfg(target_arch = "x86")]
use std::arch::x86::*;
#[cfg(target_arch = "x86_64")]
use std::arch::x86_64::*;

const LANES: usize = 8;
const MATCH_TO_MATCH: usize = 0;
const INDEL_TO_MATCH: usize = 1;
const MATCH_TO_INSERTION: usize = 2;
const INSERTION_TO_INSERTION: usize = 3;
const MATCH_TO_DELETION: usize = 4;
const DELETION_TO_DELETION: usize = 5;

/// Returns linear-space final sum (caller applies log10 / retry).
///
/// SAFETY: caller must ensure AVX2 is available.
#[target_feature(enable = "avx2")]
pub unsafe fn score_one_avx2_f32(
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
        fill_row_avx2_f32(
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

#[target_feature(enable = "avx2")]
unsafe fn fill_row_avx2_f32(
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
    let tm2m = _mm256_set1_ps(t[MATCH_TO_MATCH]);
    let ti2m = _mm256_set1_ps(t[INDEL_TO_MATCH]);
    let tm2i = _mm256_set1_ps(t[MATCH_TO_INSERTION]);
    let ti2i = _mm256_set1_ps(t[INSERTION_TO_INSERTION]);
    let tm2d = t[MATCH_TO_DELETION];
    let td2d = t[DELETION_TO_DELETION];
    let match_v = _mm256_set1_ps(match_p);
    let mismatch_v = _mm256_set1_ps(mismatch_p);
    let x_v = _mm256_set1_epi32(x as i32);
    let n_v = _mm256_set1_epi32(b'N' as i32);

    let mut j = 1usize;
    while j + LANES - 1 <= hn {
        // Emission priors for 8 hap bases.
        let mut y_i32 = [0i32; LANES];
        for lane in 0..LANES {
            y_i32[lane] = hap[j - 1 + lane] as i32;
        }
        let y_v = _mm256_loadu_si256(y_i32.as_ptr() as *const __m256i);
        let eq_xy = _mm256_cmpeq_epi32(x_v, y_v);
        let eq_xn = _mm256_cmpeq_epi32(x_v, n_v);
        let eq_yn = _mm256_cmpeq_epi32(y_v, n_v);
        let is_match = _mm256_or_si256(eq_xy, _mm256_or_si256(eq_xn, eq_yn));
        let p = _mm256_blendv_ps(mismatch_v, match_v, _mm256_castsi256_ps(is_match));

        let m_diag = _mm256_loadu_ps(m_prev.as_ptr().add(j - 1));
        let i_diag = _mm256_loadu_ps(ins_prev.as_ptr().add(j - 1));
        let d_diag = _mm256_loadu_ps(del_prev.as_ptr().add(j - 1));
        let m_up = _mm256_loadu_ps(m_prev.as_ptr().add(j));
        let i_up = _mm256_loadu_ps(ins_prev.as_ptr().add(j));

        let from_m = _mm256_mul_ps(m_diag, tm2m);
        let from_i = _mm256_mul_ps(i_diag, ti2m);
        let from_d = _mm256_mul_ps(d_diag, ti2m);
        let sum = _mm256_add_ps(_mm256_add_ps(from_m, from_i), from_d);
        let m_new = _mm256_mul_ps(p, sum);
        let i_new = _mm256_add_ps(_mm256_mul_ps(m_up, tm2i), _mm256_mul_ps(i_up, ti2i));

        _mm256_storeu_ps(m_curr.as_mut_ptr().add(j), m_new);
        _mm256_storeu_ps(ins_curr.as_mut_ptr().add(j), i_new);

        // Deletion is left-dependent — serial across the stripe.
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
