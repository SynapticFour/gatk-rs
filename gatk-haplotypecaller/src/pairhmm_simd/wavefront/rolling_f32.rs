//! Rolling 2-row f32 Logless with f64 retry; dispatches to AVX2/NEON when available.

use super::prep::ReadPrep;
use super::rolling_f64::score_one_rolling_f64;
use super::{select_wavefront_kernel, WavefrontKernel, WavefrontScratch};

/// f32-safe seed (GKL-style `ldexpf(1, 120)`). Final log10 subtracts this scale so
/// scores match f64 Logless up to rounding (`log10(c·I) − log10(I) = log10(c)`).
pub const INITIAL_CONDITION_F32: f32 = f32::from_bits((127 + 120) << 23); // 2^120
pub const INITIAL_CONDITION_F32_LOG10: f64 = 120.0 * std::f64::consts::LOG10_2;
/// Reject f32 linear sums below this (relative to 2^120 scale) and retry f64.
const MIN_ACCEPTED_F32_LINEAR_SUM: f32 = 1e-20;

const MATCH_TO_MATCH: usize = 0;
const INDEL_TO_MATCH: usize = 1;
const MATCH_TO_INSERTION: usize = 2;
const INSERTION_TO_INSERTION: usize = 3;
const MATCH_TO_DELETION: usize = 4;
const DELETION_TO_DELETION: usize = 5;

/// Host-selected f32 kernel with f64 retry on underflow.
pub fn score_one_rolling_f32_with_retry(
    prep: &ReadPrep,
    read_bases: &[u8],
    hap: &[u8],
    scratch: &mut WavefrontScratch,
) -> f64 {
    let sum = match select_wavefront_kernel() {
        WavefrontKernel::Avx2F32 => {
            #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
            {
                // SAFETY: select_wavefront_kernel verified AVX2.
                unsafe { super::avx2_f32::score_one_avx2_f32(prep, read_bases, hap, scratch) }
            }
            #[cfg(not(any(target_arch = "x86", target_arch = "x86_64")))]
            {
                score_one_portable_f32(prep, read_bases, hap, scratch)
            }
        }
        WavefrontKernel::NeonF32 => {
            #[cfg(target_arch = "aarch64")]
            {
                // SAFETY: select_wavefront_kernel verified NEON.
                unsafe { super::neon_f32::score_one_neon_f32(prep, read_bases, hap, scratch) }
            }
            #[cfg(not(target_arch = "aarch64"))]
            {
                score_one_portable_f32(prep, read_bases, hap, scratch)
            }
        }
        WavefrontKernel::PortableF32 => score_one_portable_f32(prep, read_bases, hap, scratch),
    };

    if !sum.is_finite() || sum < MIN_ACCEPTED_F32_LINEAR_SUM {
        return score_one_rolling_f64(prep, read_bases, hap, scratch);
    }
    let s = sum as f64;
    if s <= 0.0 || !s.is_finite() {
        return f64::NEG_INFINITY;
    }
    s.log10() - INITIAL_CONDITION_F32_LOG10
}

/// Portable scalar f32 rolling kernel (returns linear sum, not log10).
pub fn score_one_portable_f32(
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
        let tm2m = t[MATCH_TO_MATCH];
        let ti2m = t[INDEL_TO_MATCH];
        let tm2i = t[MATCH_TO_INSERTION];
        let ti2i = t[INSERTION_TO_INSERTION];
        let tm2d = t[MATCH_TO_DELETION];
        let td2d = t[DELETION_TO_DELETION];
        for j in 1..=hn {
            let y = unsafe { *hap.get_unchecked(j - 1) };
            let p = if x == y || x == b'N' || y == b'N' {
                match_p
            } else {
                mismatch_p
            };
            // SAFETY: j in 1..=hn; buffers length cols = hn+1.
            unsafe {
                *m_curr.get_unchecked_mut(j) = p
                    * (*m_prev.get_unchecked(j - 1) * tm2m
                        + *ins_prev.get_unchecked(j - 1) * ti2m
                        + *del_prev.get_unchecked(j - 1) * ti2m);
                *ins_curr.get_unchecked_mut(j) =
                    *m_prev.get_unchecked(j) * tm2i + *ins_prev.get_unchecked(j) * ti2i;
                *del_curr.get_unchecked_mut(j) =
                    *m_curr.get_unchecked(j - 1) * tm2d + *del_curr.get_unchecked(j - 1) * td2d;
            }
        }
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
