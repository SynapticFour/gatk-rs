//! Rolling 2-row Logless PairHMM in f64 (oracle layout ≡ scalar full-matrix).

use super::prep::ReadPrep;
use super::WavefrontScratch;
use crate::pairhmm_logless::{INITIAL_CONDITION, INITIAL_CONDITION_LOG10};

const MATCH_TO_MATCH: usize = 0;
const INDEL_TO_MATCH: usize = 1;
const MATCH_TO_INSERTION: usize = 2;
const INSERTION_TO_INSERTION: usize = 3;
const MATCH_TO_DELETION: usize = 4;
const DELETION_TO_DELETION: usize = 5;

/// Score one (read, hap) with precomputed [`ReadPrep`] using rolling 2-row f64 DP.
pub fn score_one_rolling_f64(
    prep: &ReadPrep,
    read_bases: &[u8],
    hap: &[u8],
    scratch: &mut WavefrontScratch,
) -> f64 {
    let rn = read_bases.len();
    let hn = hap.len();
    if rn == 0 {
        return 0.0;
    }
    let cols = hn + 1;
    scratch.ensure_rolling_f64(cols);
    let (m_a, m_b) = scratch.m[..cols * 2].split_at_mut(cols);
    let (ins_a, ins_b) = scratch.ins[..cols * 2].split_at_mut(cols);
    let (del_a, del_b) = scratch.del[..cols * 2].split_at_mut(cols);

    m_a.fill(0.0);
    ins_a.fill(0.0);
    let init_del = INITIAL_CONDITION / hn as f64;
    for j in 0..=hn {
        del_a[j] = init_del;
    }

    let mut prev_is_a = true;
    for i in 1..=rn {
        let t = prep.transitions_f64[i];
        let x = read_bases[i - 1];
        let (match_p, mismatch_p) = prep.match_mm_f64[i - 1];
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
        for j in 1..=hn {
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
            ins_curr[j] =
                m_prev[j] * t[MATCH_TO_INSERTION] + ins_prev[j] * t[INSERTION_TO_INSERTION];
            del_curr[j] =
                m_curr[j - 1] * t[MATCH_TO_DELETION] + del_curr[j - 1] * t[DELETION_TO_DELETION];
        }
        prev_is_a = !prev_is_a;
    }

    let (m_end, ins_end) = if prev_is_a {
        (&m_a[..], &ins_a[..])
    } else {
        (&m_b[..], &ins_b[..])
    };
    let mut final_sum = 0.0f64;
    for j in 1..=hn {
        final_sum += m_end[j] + ins_end[j];
    }
    if final_sum <= 0.0 || !final_sum.is_finite() {
        return f64::NEG_INFINITY;
    }
    final_sum.log10() - INITIAL_CONDITION_LOG10
}
