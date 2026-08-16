//! x86_64 AVX2 pairs-in-lanes Logless PairHMM (4×f64).

#![cfg(any(target_arch = "x86", target_arch = "x86_64"))]

use super::pack::{
    score_haps_logless_packed_f64, score_haps_logless_packed_f64_with_transitions,
    score_one_hap_logless_f64_with_transitions,
};
use crate::pairhmm_logless::{
    logless_build_transitions, logless_match_mismatch_prior, INITIAL_CONDITION,
    INITIAL_CONDITION_LOG10,
};
use gatk_common::GatkResult;
use std::cell::RefCell;
use std::collections::HashMap;

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

struct Avx2Scratch {
    m: Vec<f64>,
    ins: Vec<f64>,
    del: Vec<f64>,
    prior: Vec<f64>,
}

impl Avx2Scratch {
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
    static AVX2_SCRATCH: RefCell<Avx2Scratch> = RefCell::new(Avx2Scratch::new());
    static AVX2_BY_LEN: RefCell<HashMap<usize, Vec<usize>>> =
        RefCell::new(HashMap::new());
}

/// Drop AVX2 PairHMM TLS planes (Peak hygiene after a region).
pub fn release_pairhmm_avx2_tls_scratch() {
    AVX2_SCRATCH.with(|c| {
        let mut s = c.borrow_mut();
        s.m = Vec::new();
        s.ins = Vec::new();
        s.del = Vec::new();
        s.prior = Vec::new();
    });
    AVX2_BY_LEN.with(|c| c.borrow_mut().clear());
}

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

    let transitions = logless_build_transitions(rn, insertion_gop, deletion_gop, overall_gcp);
    let mut out = vec![0.0f64; haplotypes.len()];
    let mut err = None;
    AVX2_BY_LEN.with(|by_len_cell| {
        let mut by_len = by_len_cell.borrow_mut();
        for v in by_len.values_mut() {
            v.clear();
        }
        for (i, h) in haplotypes.iter().enumerate() {
            by_len.entry(h.len()).or_default().push(i);
        }
        let lengths: Vec<usize> = by_len
            .iter()
            .filter(|(_, idxs)| !idxs.is_empty())
            .map(|(len, _)| *len)
            .collect();
        AVX2_SCRATCH.with(|cell| {
            let mut scratch = cell.borrow_mut();
            for len in lengths {
                let Some(ordered) = by_len.get_mut(&len) else {
                    continue;
                };
                ordered.sort_by(|&a, &b| haplotypes[a].cmp(haplotypes[b]));
                // Long same-length chains: prefix reuse beats repeated pack4 when haps share prefixes.
                if ordered.len() >= 5 {
                    let mut subset = Vec::with_capacity(ordered.len());
                    for &i in ordered.iter() {
                        subset.push(haplotypes[i]);
                    }
                    match score_haps_logless_packed_f64_with_transitions(
                        read_bases,
                        read_quals,
                        &subset,
                        &transitions,
                    ) {
                        Ok(scores) => {
                            for (k, &i) in ordered.iter().enumerate() {
                                out[i] = scores[k];
                            }
                        }
                        Err(e) => {
                            err = Some(e);
                            return;
                        }
                    }
                    continue;
                }
                let mut chunks = ordered.chunks_exact(LANES);
                for pack_src in chunks.by_ref() {
                    let pack = [
                        haplotypes[pack_src[0]],
                        haplotypes[pack_src[1]],
                        haplotypes[pack_src[2]],
                        haplotypes[pack_src[3]],
                    ];
                    let scores =
                        score_pack4(read_bases, read_quals, &pack, &transitions, &mut scratch);
                    for (k, &idx) in pack_src.iter().enumerate() {
                        out[idx] = scores[k];
                    }
                }
                for &i in chunks.remainder() {
                    match score_one_hap_logless_f64_with_transitions(
                        read_bases,
                        read_quals,
                        haplotypes[i],
                        &transitions,
                    ) {
                        Ok(score) => out[i] = score,
                        Err(e) => {
                            err = Some(e);
                            return;
                        }
                    }
                }
            }
        });
    });
    if let Some(e) = err {
        return Err(e);
    }
    Ok(out)
}

/// SAFETY: caller must ensure AVX2 is available (`#[target_feature(enable = "avx2")]`).
#[target_feature(enable = "avx2")]
unsafe fn score_pack4(
    read_bases: &[u8],
    read_quals: &[u8],
    haps: &[&[u8]; 4],
    transitions: &[[f64; 6]],
    scratch: &mut Avx2Scratch,
) -> [f64; 4] {
    let rn = read_bases.len();
    let hn = [haps[0].len(), haps[1].len(), haps[2].len(), haps[3].len()];
    let hn_max = hn.iter().copied().max().unwrap_or(1);
    let cols = hn_max + 1;
    let cells = (rn + 1) * cols;
    scratch.ensure(cells);
    let m = &mut scratch.m;
    let ins = &mut scratch.ins;
    let del = &mut scratch.del;
    let prior = &mut scratch.prior;

    for lane in 0..LANES {
        let init = INITIAL_CONDITION / hn[lane] as f64;
        for j in 0..=hn[lane] {
            del[j * LANES + lane] = init;
        }
    }

    for i in 0..rn {
        let x = read_bases[i];
        let (match_p, mismatch_p) = logless_match_mismatch_prior(read_quals[i]);
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
