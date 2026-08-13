//! aarch64 NEON pairs-in-lanes Logless PairHMM (2×f64).

#![cfg(target_arch = "aarch64")]

use super::pack::score_haps_logless_packed_f64;
use crate::pairhmm_logless::{INITIAL_CONDITION, INITIAL_CONDITION_LOG10};
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
        } else {
            self.m[..need].fill(0.0);
            self.ins[..need].fill(0.0);
            self.del[..need].fill(0.0);
            self.prior[..need].fill(0.0);
        }
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
    // NEON is baseline on aarch64 Darwin/Linux for our targets.
    // Do not reorder-by-length for packing: NEON pack2 still diverges from scalar on some
    // equal-length packs (pairhmm_simd_vs_scalar_test); keep adjacent-only packing.
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

    let transitions = build_transitions(rn, insertion_gop, deletion_gop, overall_gcp);
    let mut out = vec![0.0f64; haplotypes.len()];
    let mut err: Option<gatk_common::GatkError> = None;
    NEON_SCRATCH.with(|cell| {
        let mut scratch = cell.borrow_mut();
        let mut i = 0;
        while i < haplotypes.len() {
            let remaining = haplotypes.len() - i;
            // Uneven hap lengths in one NEON pack corrupt lane DP (mask path); fall back.
            if remaining >= LANES && haplotypes[i].len() == haplotypes[i + 1].len() {
                let pack = [haplotypes[i], haplotypes[i + 1]];
                let scores = score_pack2(read_bases, read_quals, &pack, &transitions, &mut scratch);
                out[i..i + LANES].copy_from_slice(&scores);
                i += LANES;
            } else {
                match score_haps_logless_packed_f64(
                    read_bases,
                    read_quals,
                    &haplotypes[i..=i],
                    insertion_gop,
                    deletion_gop,
                    overall_gcp,
                ) {
                    Ok(rest) => out[i] = rest[0],
                    Err(e) => {
                        err = Some(e);
                        break;
                    }
                }
                i += 1;
            }
        }
    });
    if let Some(e) = err {
        return Err(e);
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

#[target_feature(enable = "neon")]
unsafe fn score_pack2(
    read_bases: &[u8],
    read_quals: &[u8],
    haps: &[&[u8]; 2],
    transitions: &[[f64; 6]],
    scratch: &mut NeonScratch,
) -> [f64; 2] {
    let rn = read_bases.len();
    let hn = [haps[0].len(), haps[1].len()];
    let hn_max = hn[0].max(hn[1]);
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
        let tm2m = vdupq_n_f64(t[MATCH_TO_MATCH]);
        let ti2m = vdupq_n_f64(t[INDEL_TO_MATCH]);
        let tm2i = vdupq_n_f64(t[MATCH_TO_INSERTION]);
        let ti2i = vdupq_n_f64(t[INSERTION_TO_INSERTION]);
        let tm2d = vdupq_n_f64(t[MATCH_TO_DELETION]);
        let td2d = vdupq_n_f64(t[DELETION_TO_DELETION]);
        let row = i * cols;
        let prev = (i - 1) * cols;
        for j in 1..=hn_max {
            let active0 = if j <= hn[0] { 1.0 } else { 0.0 };
            let active1 = if j <= hn[1] { 1.0 } else { 0.0 };
            let mask = vsetq_lane_f64(active1, vdupq_n_f64(active0), 1);

            let idx = (row + j) * LANES;
            let diag = (prev + j - 1) * LANES;
            let up = (prev + j) * LANES;
            let left = (row + j - 1) * LANES;

            // SAFETY: SoA buffers sized for cells*LANES; NEON enabled.
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
            let m_new = vmulq_f64(vmulq_f64(p, sum), mask);
            let i_new = vmulq_f64(
                vaddq_f64(vmulq_f64(m_up, tm2i), vmulq_f64(i_up, ti2i)),
                mask,
            );
            let d_new = vmulq_f64(
                vaddq_f64(vmulq_f64(m_left, tm2d), vmulq_f64(d_left, td2d)),
                mask,
            );

            vst1q_f64(m.as_mut_ptr().add(idx), m_new);
            vst1q_f64(ins.as_mut_ptr().add(idx), i_new);
            vst1q_f64(del.as_mut_ptr().add(idx), d_new);
        }
    }

    let end_row = rn * cols;
    let mut sums = [0.0f64; 2];
    for j in 1..=hn_max {
        let idx = (end_row + j) * LANES;
        let mv = vld1q_f64(m.as_ptr().add(idx));
        let iv = vld1q_f64(ins.as_ptr().add(idx));
        let s = vaddq_f64(mv, iv);
        let mut tmp = [0.0f64; 2];
        vst1q_f64(tmp.as_mut_ptr(), s);
        for lane in 0..LANES {
            if j <= hn[lane] {
                sums[lane] += tmp[lane];
            }
        }
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
