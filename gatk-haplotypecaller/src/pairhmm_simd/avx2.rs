//! x86_64 AVX2 pairs-in-lanes Logless PairHMM (4×f64).

#![cfg(any(target_arch = "x86", target_arch = "x86_64"))]

use super::pack::score_haps_logless_packed_f64;
use crate::pairhmm_logless::{INITIAL_CONDITION, INITIAL_CONDITION_LOG10};
use gatk_common::GatkResult;

#[cfg(target_arch = "x86")]
use std::arch::x86::*;
#[cfg(target_arch = "x86_64")]
use std::arch::x86_64::*;

const MATCH_TO_MATCH: usize = 0;
const INDEL_TO_MATCH: usize = 1;
const MATCH_TO_INSERTION: usize = 2;
const INSERTION_TO_INSERTION: usize = 3;
const MATCH_TO_DELETION: usize = 4;
const DELETION_TO_DELETION: usize = 5;

const LANES: usize = 4;

/// Score haplotypes with AVX2 when available; otherwise portable packed f64.
pub fn score_haps_avx2_f64(
    read_bases: &[u8],
    read_quals: &[u8],
    haplotypes: &[&[u8]],
    insertion_gop: &[u8],
    deletion_gop: &[u8],
    overall_gcp: &[u8],
) -> GatkResult<Vec<f64>> {
    if !is_x86_feature_detected!("avx2") {
        return score_haps_logless_packed_f64(
            read_bases,
            read_quals,
            haplotypes,
            insertion_gop,
            deletion_gop,
            overall_gcp,
        );
    }
    // SAFETY: `is_x86_feature_detected!("avx2")` is true above.
    unsafe {
        score_haps_avx2_f64_unchecked(
            read_bases,
            read_quals,
            haplotypes,
            insertion_gop,
            deletion_gop,
            overall_gcp,
        )
    }
}

#[target_feature(enable = "avx2")]
unsafe fn score_haps_avx2_f64_unchecked(
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

    let transitions = build_transitions(rn, insertion_gop, deletion_gop, overall_gcp);
    let mut out = vec![0.0f64; haplotypes.len()];
    let mut i = 0;
    while i < haplotypes.len() {
        let remaining = haplotypes.len() - i;
        if remaining >= LANES {
            let pack = [
                haplotypes[i],
                haplotypes[i + 1],
                haplotypes[i + 2],
                haplotypes[i + 3],
            ];
            let scores = score_pack4(read_bases, read_quals, &pack, &transitions);
            out[i..i + LANES].copy_from_slice(&scores);
            i += LANES;
        } else {
            let slice = &haplotypes[i..];
            let rest = score_haps_logless_packed_f64(
                read_bases,
                read_quals,
                slice,
                insertion_gop,
                deletion_gop,
                overall_gcp,
            )?;
            out[i..].copy_from_slice(&rest);
            break;
        }
    }
    Ok(out)
}

fn build_transitions(
    rn: usize,
    insertion_gop: &[u8],
    deletion_gop: &[u8],
    overall_gcp: &[u8],
) -> Vec<[f64; 6]> {
    let mut transitions = vec![[0.0f64; 6]; rn + 1];
    for i in 0..rn {
        let ins_q = insertion_gop[i];
        let del_q = deletion_gop[i];
        let gcp = overall_gcp[i];
        let gcp_err = 10f64.powf(-(gcp as f64) / 10.0);
        let ins_err = 10f64.powf(-(ins_q as f64) / 10.0);
        let del_err = 10f64.powf(-(del_q as f64) / 10.0);
        let (min_q, max_q) = if ins_q <= del_q {
            (ins_q as f64, del_q as f64)
        } else {
            (del_q as f64, ins_q as f64)
        };
        let log10_sum = {
            let a = -0.1 * min_q;
            let b = -0.1 * max_q;
            let (x, y) = if a > b { (b, a) } else { (a, b) };
            y + (1.0 + 10f64.powf(x - y)).log10()
        };
        let m2m = 1.0 - 10f64.powf(log10_sum).min(1.0);
        transitions[i + 1] = [m2m, 1.0 - gcp_err, ins_err, gcp_err, del_err, gcp_err];
    }
    transitions
}

/// SAFETY: caller must ensure AVX2 is available (`#[target_feature(enable = "avx2")]`).
#[target_feature(enable = "avx2")]
unsafe fn score_pack4(
    read_bases: &[u8],
    read_quals: &[u8],
    haps: &[&[u8]; 4],
    transitions: &[[f64; 6]],
) -> [f64; 4] {
    let rn = read_bases.len();
    let hn = [haps[0].len(), haps[1].len(), haps[2].len(), haps[3].len()];
    let hn_max = hn.iter().copied().max().unwrap_or(1);
    let cols = hn_max + 1;
    let cells = (rn + 1) * cols;

    // SoA: cell-major, 4 lanes contiguous for SIMD load.
    let mut m = vec![0.0f64; cells * LANES];
    let mut ins = vec![0.0f64; cells * LANES];
    let mut del = vec![0.0f64; cells * LANES];
    let mut prior = vec![0.0f64; cells * LANES];

    for lane in 0..LANES {
        let init = INITIAL_CONDITION / hn[lane] as f64;
        for j in 0..=hn[lane] {
            del[j * LANES + lane] = init;
        }
    }

    for i in 0..rn {
        let x = read_bases[i];
        let match_p = 1.0 - 10f64.powf(-(read_quals[i] as f64) / 10.0);
        let mismatch_p = 10f64.powf(-(read_quals[i] as f64) / 10.0) / 3.0;
        let row = (i + 1) * cols;
        for lane in 0..LANES {
            for j in 0..hn[lane] {
                let y = haps[lane][j];
                let p = if x == y || x == b'N' || y == b'N' {
                    match_p
                } else {
                    mismatch_p
                };
                prior[(row + j + 1) * LANES + lane] = p;
            }
        }
    }

    for i in 1..=rn {
        let t = transitions[i];
        let tm2m = _mm256_set1_pd(t[MATCH_TO_MATCH]);
        let ti2m = _mm256_set1_pd(t[INDEL_TO_MATCH]);
        let tm2i = _mm256_set1_pd(t[MATCH_TO_INSERTION]);
        let ti2i = _mm256_set1_pd(t[INSERTION_TO_INSERTION]);
        let tm2d = _mm256_set1_pd(t[MATCH_TO_DELETION]);
        let td2d = _mm256_set1_pd(t[DELETION_TO_DELETION]);
        let row = i * cols;
        let prev = (i - 1) * cols;
        for j in 1..=hn_max {
            // Active mask: lane participates if j <= hn[lane]
            let active = [
                if j <= hn[0] { 1.0 } else { 0.0 },
                if j <= hn[1] { 1.0 } else { 0.0 },
                if j <= hn[2] { 1.0 } else { 0.0 },
                if j <= hn[3] { 1.0 } else { 0.0 },
            ];
            let mask = _mm256_set_pd(active[3], active[2], active[1], active[0]);

            let idx = (row + j) * LANES;
            let diag = ((prev + j - 1) * LANES) as isize;
            let up = ((prev + j) * LANES) as isize;
            let left = ((row + j - 1) * LANES) as isize;

            // SAFETY: indices are within allocated SoA buffers; AVX2 enabled by target_feature.
            let m_diag = _mm256_loadu_pd(m.as_ptr().offset(diag));
            let i_diag = _mm256_loadu_pd(ins.as_ptr().offset(diag));
            let d_diag = _mm256_loadu_pd(del.as_ptr().offset(diag));
            let m_up = _mm256_loadu_pd(m.as_ptr().offset(up));
            let i_up = _mm256_loadu_pd(ins.as_ptr().offset(up));
            let m_left = _mm256_loadu_pd(m.as_ptr().offset(left));
            let d_left = _mm256_loadu_pd(del.as_ptr().offset(left));
            let p = _mm256_loadu_pd(prior.as_ptr().add(idx));

            let from_m = _mm256_mul_pd(m_diag, tm2m);
            let from_i = _mm256_mul_pd(i_diag, ti2m);
            let from_d = _mm256_mul_pd(d_diag, ti2m);
            let sum = _mm256_add_pd(_mm256_add_pd(from_m, from_i), from_d);
            let mut m_new = _mm256_mul_pd(p, sum);
            m_new = _mm256_mul_pd(m_new, mask);

            let mut i_new = _mm256_add_pd(_mm256_mul_pd(m_up, tm2i), _mm256_mul_pd(i_up, ti2i));
            i_new = _mm256_mul_pd(i_new, mask);

            // Deletion depends on newly written M at (i,j-1) — already in m_left for active left.
            let mut d_new = _mm256_add_pd(_mm256_mul_pd(m_left, tm2d), _mm256_mul_pd(d_left, td2d));
            // For j==1, m_left is row+0 which is 0 — correct.
            // After writing m_new we need d to use m_new for next j; store m first.
            _mm256_storeu_pd(m.as_mut_ptr().add(idx), m_new);
            // Reload m_left equivalent: for this cell's del we needed m[row+j-1], already loaded.
            // But wait — standard DP: del[i][j] uses match[i][j-1] which for same row was just
            // computed at previous j. m_left was loaded before we stored m_new at current j,
            // and points to j-1 — correct.
            d_new = _mm256_mul_pd(d_new, mask);
            _mm256_storeu_pd(ins.as_mut_ptr().add(idx), i_new);
            _mm256_storeu_pd(del.as_mut_ptr().add(idx), d_new);
        }
    }

    let end_row = rn * cols;
    let mut sums = [0.0f64; 4];
    for j in 1..=hn_max {
        let idx = (end_row + j) * LANES;
        let mv = _mm256_loadu_pd(m.as_ptr().add(idx));
        let iv = _mm256_loadu_pd(ins.as_ptr().add(idx));
        let s = _mm256_add_pd(mv, iv);
        let mut tmp = [0.0f64; 4];
        _mm256_storeu_pd(tmp.as_mut_ptr(), s);
        for lane in 0..LANES {
            if j <= hn[lane] {
                sums[lane] += tmp[lane];
            }
        }
    }

    let mut out = [0.0f64; 4];
    for lane in 0..LANES {
        if sums[lane] <= 0.0 || !sums[lane].is_finite() {
            out[lane] = f64::NEG_INFINITY;
        } else {
            out[lane] = sums[lane].log10() - INITIAL_CONDITION_LOG10;
        }
    }
    out
}
