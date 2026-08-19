//! x86_64 AVX2 pairs-in-lanes Logless PairHMM (4×f64).

#![cfg(any(target_arch = "x86", target_arch = "x86_64"))]

use super::pack::{
    mean_consecutive_prefix_frac, score_haps_logless_packed_f64,
    score_haps_logless_packed_f64_with_transitions, score_one_hap_logless_f64_with_transitions,
};
use crate::pairhmm_logless::{
    logless_fill_transitions, logless_match_mismatch_prior, INITIAL_CONDITION,
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
}

impl Avx2Scratch {
    fn new() -> Self {
        Self {
            m: Vec::new(),
            ins: Vec::new(),
            del: Vec::new(),
        }
    }

    fn ensure(&mut self, cells: usize) {
        let need = cells * LANES;
        if self.m.len() < need {
            self.m.resize(need, 0.0);
            self.ins.resize(need, 0.0);
            self.del.resize(need, 0.0);
        }
    }
}

thread_local! {
    static AVX2_SCRATCH: RefCell<Avx2Scratch> = RefCell::new(Avx2Scratch::new());
    static AVX2_TRANSITIONS: RefCell<Vec<[f64; 6]>> = const { RefCell::new(Vec::new()) };
    static AVX2_BY_LEN: RefCell<HashMap<usize, Vec<usize>>> =
        RefCell::new(HashMap::new());
    static AVX2_PACK4: std::cell::Cell<u64> = const { std::cell::Cell::new(0) };
    static AVX2_PREFIX_REUSE: std::cell::Cell<u64> = const { std::cell::Cell::new(0) };
    static AVX2_LEFTOVER: std::cell::Cell<u64> = const { std::cell::Cell::new(0) };
}

/// Keep AVX2 PairHMM TLS high-water (see `run::release_region_tls_scratch`).
pub fn release_pairhmm_avx2_tls_scratch() {}

/// Take-and-reset AVX2 pack occupancy: `(pack4, prefix_reuse_haps, leftover_singles)`.
pub fn take_avx2_pack_stats() -> (u64, u64, u64) {
    (
        AVX2_PACK4.replace(0),
        AVX2_PREFIX_REUSE.replace(0),
        AVX2_LEFTOVER.replace(0),
    )
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

    let mut out = vec![0.0f64; haplotypes.len()];
    let mut err = None;
    AVX2_TRANSITIONS.with(|tcell| {
        let mut transitions = tcell.borrow_mut();
        logless_fill_transitions(
            &mut transitions,
            rn,
            insertion_gop,
            deletion_gop,
            overall_gcp,
        );
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
                    let mut subset = Vec::with_capacity(ordered.len());
                    for &i in ordered.iter() {
                        subset.push(haplotypes[i]);
                    }
                    let use_prefix = ordered.len() >= 5
                        && mean_consecutive_prefix_frac(&subset)
                            >= super::pack::PREFIX_REUSE_OVER_SIMD_FRAC;
                    if use_prefix {
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
                                AVX2_PREFIX_REUSE
                                    .set(AVX2_PREFIX_REUSE.get() + ordered.len() as u64);
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
                        AVX2_PACK4.set(AVX2_PACK4.get() + 1);
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
                        AVX2_LEFTOVER.set(AVX2_LEFTOVER.get() + 1);
                    }
                }
            });
        });
    });
    if let Some(e) = err {
        return Err(e);
    }
    Ok(out)
}

/// SAFETY: caller must ensure AVX2 is available (`#[target_feature(enable = "avx2")]`).
/// Equal-length packs only (`by_len` groups) — no per-j lane masks.
#[target_feature(enable = "avx2")]
unsafe fn score_pack4(
    read_bases: &[u8],
    read_quals: &[u8],
    haps: &[&[u8]; 4],
    transitions: &[[f64; 6]],
    scratch: &mut Avx2Scratch,
) -> [f64; 4] {
    let rn = read_bases.len();
    let hn = haps[0].len();
    debug_assert!(haps.iter().all(|h| h.len() == hn));
    let cols = hn + 1;
    let cells = (rn + 1) * cols;
    scratch.ensure(cells);
    let m = &mut scratch.m;
    let ins = &mut scratch.ins;
    let del = &mut scratch.del;

    m[..cols * LANES].fill(0.0);
    ins[..cols * LANES].fill(0.0);
    let init = INITIAL_CONDITION / hn as f64;
    for j in 0..=hn {
        let base = j * LANES;
        for lane in 0..LANES {
            del[base + lane] = init;
        }
    }

    for i in 1..=rn {
        let t = transitions[i];
        let x = read_bases[i - 1];
        let (match_p, mismatch_p) = logless_match_mismatch_prior(read_quals[i - 1]);
        let tm2m = _mm256_set1_pd(t[MATCH_TO_MATCH]);
        let ti2m = _mm256_set1_pd(t[INDEL_TO_MATCH]);
        let tm2i = _mm256_set1_pd(t[MATCH_TO_INSERTION]);
        let ti2i = _mm256_set1_pd(t[INSERTION_TO_INSERTION]);
        let tm2d = _mm256_set1_pd(t[MATCH_TO_DELETION]);
        let td2d = _mm256_set1_pd(t[DELETION_TO_DELETION]);
        let row = i * cols;
        let prev = (i - 1) * cols;
        let col0 = row * LANES;
        for lane in 0..LANES {
            m[col0 + lane] = 0.0;
            ins[col0 + lane] = 0.0;
            del[col0 + lane] = 0.0;
        }
        for j in 1..=hn {
            let idx = (row + j) * LANES;
            let diag = ((prev + j - 1) * LANES) as isize;
            let up = ((prev + j) * LANES) as isize;
            let left = ((row + j - 1) * LANES) as isize;

            let mut prior_arr = [0.0f64; 4];
            for lane in 0..LANES {
                let y = haps[lane][j - 1];
                prior_arr[lane] = if x == y || x == b'N' || y == b'N' {
                    match_p
                } else {
                    mismatch_p
                };
            }

            // SAFETY: indices within SoA; AVX2 enabled; equal-length pack.
            let m_diag = _mm256_loadu_pd(m.as_ptr().offset(diag));
            let i_diag = _mm256_loadu_pd(ins.as_ptr().offset(diag));
            let d_diag = _mm256_loadu_pd(del.as_ptr().offset(diag));
            let m_up = _mm256_loadu_pd(m.as_ptr().offset(up));
            let i_up = _mm256_loadu_pd(ins.as_ptr().offset(up));
            let m_left = _mm256_loadu_pd(m.as_ptr().offset(left));
            let d_left = _mm256_loadu_pd(del.as_ptr().offset(left));
            let p = _mm256_loadu_pd(prior_arr.as_ptr());

            let from_m = _mm256_mul_pd(m_diag, tm2m);
            let from_i = _mm256_mul_pd(i_diag, ti2m);
            let from_d = _mm256_mul_pd(d_diag, ti2m);
            let sum = _mm256_add_pd(_mm256_add_pd(from_m, from_i), from_d);
            let m_new = _mm256_mul_pd(p, sum);
            let i_new = _mm256_add_pd(_mm256_mul_pd(m_up, tm2i), _mm256_mul_pd(i_up, ti2i));
            let d_new = _mm256_add_pd(_mm256_mul_pd(m_left, tm2d), _mm256_mul_pd(d_left, td2d));

            _mm256_storeu_pd(m.as_mut_ptr().add(idx), m_new);
            _mm256_storeu_pd(ins.as_mut_ptr().add(idx), i_new);
            _mm256_storeu_pd(del.as_mut_ptr().add(idx), d_new);
        }
    }

    let end_row = rn * cols;
    let mut sums = [0.0f64; 4];
    for j in 1..=hn {
        let idx = (end_row + j) * LANES;
        let mv = _mm256_loadu_pd(m.as_ptr().add(idx));
        let iv = _mm256_loadu_pd(ins.as_ptr().add(idx));
        let s = _mm256_add_pd(mv, iv);
        let mut tmp = [0.0f64; 4];
        _mm256_storeu_pd(tmp.as_mut_ptr(), s);
        for lane in 0..LANES {
            sums[lane] += tmp[lane];
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
