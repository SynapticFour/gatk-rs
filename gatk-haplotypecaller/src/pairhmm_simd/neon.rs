//! aarch64 NEON pairs-in-lanes Logless PairHMM (2×f64).

#![cfg(target_arch = "aarch64")]

use super::pack::{
    mean_consecutive_prefix_frac, score_haps_logless_packed_f64,
    score_haps_logless_packed_f64_with_transitions, score_one_hap_logless_f64_with_transitions,
};
use crate::pairhmm_logless::{
    logless_build_transitions, logless_match_mismatch_prior, INITIAL_CONDITION,
    INITIAL_CONDITION_LOG10,
};
use gatk_common::GatkResult;
use std::arch::aarch64::*;
use std::cell::{Cell, RefCell};
use std::collections::HashMap;

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
    }
}

thread_local! {
    static NEON_SCRATCH: RefCell<NeonScratch> = RefCell::new(NeonScratch::new());
    /// Reused equal-length index groups (cleared each call; keeps HashMap + Vec capacity).
    static NEON_BY_LEN: RefCell<HashMap<usize, Vec<usize>>> = RefCell::new(HashMap::new());
    /// TRACE occupancy: pack2 SIMD hits, hapStartIndex prefix-reuse haps, scalar singles.
    static NEON_PACK2: Cell<u64> = const { Cell::new(0) };
    static NEON_PREFIX_REUSE: Cell<u64> = const { Cell::new(0) };
    static NEON_LEFTOVER: Cell<u64> = const { Cell::new(0) };
}

/// Keep NEON PairHMM TLS high-water (see `run::release_region_tls_scratch`).
pub fn release_pairhmm_neon_tls_scratch() {}

/// Take-and-reset NEON pack occupancy: `(pack2, prefix_reuse_haps, leftover_singles)`.
pub fn take_neon_pack_stats() -> (u64, u64, u64) {
    (
        NEON_PACK2.replace(0),
        NEON_PREFIX_REUSE.replace(0),
        NEON_LEFTOVER.replace(0),
    )
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
    // Group equal lengths (AVX2-style `by_len`), not a full hap sort. Sort/scatter was
    // slower under SEQUENTIAL=1; HashMap group + chunks keeps pack density without that tax.
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
    // TLS by_len: clear buckets but keep HashMap + Vec capacity across reads.
    let mut err: Option<gatk_common::GatkError> = None;
    NEON_BY_LEN.with(|by_len_cell| {
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
        NEON_SCRATCH.with(|cell| {
            let mut scratch = cell.borrow_mut();
            for len in lengths {
                let Some(ordered) = by_len.get_mut(&len) else {
                    continue;
                };
                // Sort in place — no idxs.clone(); score-invariant hap order for prefix reuse.
                ordered.sort_by(|&a, &b| haplotypes[a].cmp(haplotypes[b]));
                let mut subset = Vec::with_capacity(ordered.len());
                for &i in ordered.iter() {
                    subset.push(haplotypes[i]);
                }
                // Java Logless hapStartIndex vs GKL SIMD lanes: prefix reuse wins when
                // sorted same-length haps share long prefixes; otherwise pack2.
                let use_prefix = ordered.len() >= 3
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
                            NEON_PREFIX_REUSE.set(NEON_PREFIX_REUSE.get() + ordered.len() as u64);
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
                    let pack = [haplotypes[pack_src[0]], haplotypes[pack_src[1]]];
                    let scores =
                        score_pack2(read_bases, read_quals, &pack, &transitions, &mut scratch);
                    out[pack_src[0]] = scores[0];
                    out[pack_src[1]] = scores[1];
                    NEON_PACK2.set(NEON_PACK2.get() + 1);
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
                    NEON_LEFTOVER.set(NEON_LEFTOVER.get() + 1);
                }
            }
        });
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
    // Equal-length packs only (caller gates). Skip full SoA memset: overwrite interiors,
    // seed row-0 free dels, zero col-0 each read row (stale row-0 M/I corrupts DP).
    let rn = read_bases.len();
    let hn = haps[0].len();
    debug_assert_eq!(haps[0].len(), haps[1].len());
    let cols = hn + 1;
    let cells = (rn + 1) * cols;
    scratch.ensure(cells);
    let m = &mut scratch.m;
    let ins = &mut scratch.ins;
    let del = &mut scratch.del;

    // Row-0 M/I must be 0 (may be stale after grow); free leading dels.
    m[..cols * LANES].fill(0.0);
    ins[..cols * LANES].fill(0.0);
    let init = INITIAL_CONDITION / hn as f64;
    for j in 0..=hn {
        let base = j * LANES;
        del[base] = init;
        del[base + 1] = init;
    }

    for i in 1..=rn {
        let t = transitions[i];
        let x = read_bases[i - 1];
        let (match_p, mismatch_p) = logless_match_mismatch_prior(read_quals[i - 1]);
        let tm2m = vdupq_n_f64(t[MATCH_TO_MATCH]);
        let ti2m = vdupq_n_f64(t[INDEL_TO_MATCH]);
        let tm2i = vdupq_n_f64(t[MATCH_TO_INSERTION]);
        let ti2i = vdupq_n_f64(t[INSERTION_TO_INSERTION]);
        let tm2d = vdupq_n_f64(t[MATCH_TO_DELETION]);
        let td2d = vdupq_n_f64(t[DELETION_TO_DELETION]);
        let row = i * cols;
        let prev = (i - 1) * cols;
        // Col 0 stays 0 for del[j=1] left term.
        let col0 = row * LANES;
        m[col0] = 0.0;
        m[col0 + 1] = 0.0;
        ins[col0] = 0.0;
        ins[col0 + 1] = 0.0;
        del[col0] = 0.0;
        del[col0 + 1] = 0.0;
        for j in 1..=hn {
            let idx = (row + j) * LANES;
            let diag = (prev + j - 1) * LANES;
            let up = (prev + j) * LANES;
            let left = (row + j - 1) * LANES;

            let y0 = haps[0][j - 1];
            let y1 = haps[1][j - 1];
            let p0 = if x == y0 || x == b'N' || y0 == b'N' {
                match_p
            } else {
                mismatch_p
            };
            let p1 = if x == y1 || x == b'N' || y1 == b'N' {
                match_p
            } else {
                mismatch_p
            };
            let prior_arr = [p0, p1];
            let p = vld1q_f64(prior_arr.as_ptr());

            // SAFETY: SoA buffers sized for cells*LANES; NEON enabled; equal-length pack.
            let m_diag = vld1q_f64(m.as_ptr().add(diag));
            let i_diag = vld1q_f64(ins.as_ptr().add(diag));
            let d_diag = vld1q_f64(del.as_ptr().add(diag));
            let m_up = vld1q_f64(m.as_ptr().add(up));
            let i_up = vld1q_f64(ins.as_ptr().add(up));
            let m_left = vld1q_f64(m.as_ptr().add(left));
            let d_left = vld1q_f64(del.as_ptr().add(left));

            let sum = vaddq_f64(
                vaddq_f64(vmulq_f64(m_diag, tm2m), vmulq_f64(i_diag, ti2m)),
                vmulq_f64(d_diag, ti2m),
            );
            let m_new = vmulq_f64(p, sum);
            let i_new = vaddq_f64(vmulq_f64(m_up, tm2i), vmulq_f64(i_up, ti2i));
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
