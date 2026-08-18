//! Portable pairs-in-lanes Logless PairHMM (explicit lanes; used as fallback and
//! for uneven packs). AVX2/NEON specialize full-width packs.
//!
//! Scratch planes are sized once to `max(hap_len)` and reused across haplotypes for
//! cache locality (phenotype: many haps × shared read). Numerics match scalar Logless.

use crate::pairhmm_logless::{
    logless_match_mismatch_prior, logless_pairhmm_likelihood, logless_qual_to_trans_probs,
    INITIAL_CONDITION, INITIAL_CONDITION_LOG10, MIN_ACCEPTED_LINEAR_SUM,
};
use gatk_common::{GatkError, GatkResult};
use std::cell::RefCell;

const MATCH_TO_MATCH: usize = 0;
const INDEL_TO_MATCH: usize = 1;
const MATCH_TO_INSERTION: usize = 2;
const INSERTION_TO_INSERTION: usize = 3;
const MATCH_TO_DELETION: usize = 4;
const DELETION_TO_DELETION: usize = 5;

thread_local! {
    static PACK_F64_SCRATCH: RefCell<F64Scratch> = RefCell::new(F64Scratch::empty());
    static PACK_F32_SCRATCH: RefCell<F32Scratch> = RefCell::new(F32Scratch::empty());
}

struct F64Scratch {
    m: Vec<f64>,
    ins: Vec<f64>,
    del: Vec<f64>,
    prior: Vec<f64>,
}

impl F64Scratch {
    fn empty() -> Self {
        Self {
            m: Vec::new(),
            ins: Vec::new(),
            del: Vec::new(),
            prior: Vec::new(),
        }
    }

    fn ensure_cells(&mut self, cells: usize) {
        if self.m.len() < cells {
            self.m.resize(cells, 0.0);
            self.ins.resize(cells, 0.0);
            self.del.resize(cells, 0.0);
            self.prior.resize(cells, 0.0);
        }
    }

    /// Rolling 2-row path only needs `2 * cols` in m/ins/del (prior unused).
    fn ensure_rolling_cols(&mut self, cols: usize) {
        let need = cols.saturating_mul(2);
        if self.m.len() < need {
            self.m.resize(need, 0.0);
            self.ins.resize(need, 0.0);
            self.del.resize(need, 0.0);
        }
    }
}

/// Score one read against many haplotypes with a portable packed f64 kernel.
/// Haplotypes share one DP scratch sized to the longest hap (cache locality).
pub fn score_haps_logless_packed_f64(
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

    let mut transitions = vec![[0.0f64; 6]; rn + 1];
    for i in 0..rn {
        transitions[i + 1] =
            logless_qual_to_trans_probs(insertion_gop[i], deletion_gop[i], overall_gcp[i]);
    }

    let max_hn = haplotypes.iter().map(|h| h.len()).max().unwrap_or(0);
    let max_cols = max_hn + 1;
    let max_cells = (rn + 1).saturating_mul(max_cols);
    // Match scalar PairHMM fail-closed caps (Peak-RSS on 16 GiB hosts).
    const MAX_PAIRHMM_DIM: usize = 100_000;
    const MAX_PAIRHMM_CELLS: usize = 8_000_000;
    if rn > MAX_PAIRHMM_DIM || max_hn > MAX_PAIRHMM_DIM || max_cells > MAX_PAIRHMM_CELLS {
        return Err(GatkError::algorithm(format!(
            "PairHMM packed-f64 refused oversized DP (read_len={rn}, max_hap_len={max_hn}, cells={max_cells}); \
             inputs must be assembly-region scale, not contig scale"
        )));
    }

    let mut out = Vec::with_capacity(haplotypes.len());
    PACK_F64_SCRATCH.with(|cell| {
        let mut scratch = cell.borrow_mut();
        scratch.ensure_cells(max_cells);
        let mut prev: Option<&[u8]> = None;
        for &hap in haplotypes {
            let (hap_start, reinit_del) = match prev {
                Some(p) if p.len() == hap.len() => (first_hap_divergence(p, hap), false),
                _ => (0, true),
            };
            out.push(score_one_f64(
                read_bases,
                read_quals,
                hap,
                &transitions,
                hap_start,
                reinit_del,
                &mut scratch,
            ));
            prev = Some(hap);
        }
    });
    Ok(out)
}

/// First index where two haplotypes differ (GATK `PairHMM.findFirstPositionWhereHaplotypesDiffer`).
/// If identical through `min(len)`, returns that min length.
#[inline]
pub(crate) fn first_hap_divergence(a: &[u8], b: &[u8]) -> usize {
    let n = a.len().min(b.len());
    let mut i = 0;
    // Word-at-a-time reject before byte scan (same-length assembly haps share long prefixes).
    while i + 8 <= n {
        let aw = u64::from_ne_bytes(a[i..i + 8].try_into().unwrap());
        let bw = u64::from_ne_bytes(b[i..i + 8].try_into().unwrap());
        if aw != bw {
            break;
        }
        i += 8;
    }
    while i < n {
        if a[i] != b[i] {
            return i;
        }
        i += 1;
    }
    n
}

/// Mean fraction of haplotype bases covered by consecutive `hapStartIndex` prefixes
/// after lexicographic sort. High ⇒ scalar prefix reuse beats 2/4-wide SIMD packs
/// (Java Logless `hapStartIndex`); low ⇒ SIMD packs win (GKL-style lane throughput).
pub(crate) fn mean_consecutive_prefix_frac(haps: &[&[u8]]) -> f64 {
    if haps.len() < 2 {
        return 0.0;
    }
    let mut prefix = 0usize;
    let mut total = 0usize;
    for w in haps.windows(2) {
        total += w[1].len();
        prefix += first_hap_divergence(w[0], w[1]);
    }
    if total == 0 {
        0.0
    } else {
        prefix as f64 / total as f64
    }
}

/// Prefer Java hapStartIndex scalar reuse when consecutive same-length haps share
/// this fraction of bases; otherwise use SIMD packs (GKL lane throughput).
pub(crate) const PREFIX_REUSE_OVER_SIMD_FRAC: f64 = 0.35;

/// Score one hap with prebuilt transition planes (shared across a read's hap pack).
///
/// `hap_start`: 0-based haplotype index to (re)compute from — columns before this are
/// assumed valid from a prior same-length hap with an identical prefix (Java contract).
/// `reinit_del`: when true, refresh free leading-deletion row (`INITIAL / hn`).
///
/// Always uses the full-matrix path so consecutive same-length haps can reuse prefix
/// columns. Leftover singles that never reuse should call [`score_one_f64_rolling`].
pub(crate) fn score_one_f64(
    read_bases: &[u8],
    read_quals: &[u8],
    hap: &[u8],
    transitions: &[[f64; 6]],
    hap_start: usize,
    reinit_del: bool,
    scratch: &mut F64Scratch,
) -> f64 {
    let start = hap_start.min(hap.len());
    score_one_f64_prefix_reuse(
        read_bases,
        read_quals,
        hap,
        transitions,
        start,
        reinit_del,
        scratch,
    )
}

/// Score one haplotype with prebuilt transitions (NEON leftover singles; avoids Vec alloc).
/// Uses rolling 2-row DP — safe because leftovers do not participate in hapStartIndex reuse.
pub(crate) fn score_one_hap_logless_f64_with_transitions(
    read_bases: &[u8],
    read_quals: &[u8],
    hap: &[u8],
    transitions: &[[f64; 6]],
) -> GatkResult<f64> {
    let rn = read_bases.len();
    if rn == 0 {
        return Ok(0.0);
    }
    let hn = hap.len();
    let max_cells = (rn + 1).saturating_mul(hn + 1);
    const MAX_PAIRHMM_DIM: usize = 100_000;
    const MAX_PAIRHMM_CELLS: usize = 8_000_000;
    if rn > MAX_PAIRHMM_DIM || hn > MAX_PAIRHMM_DIM || max_cells > MAX_PAIRHMM_CELLS {
        return Err(GatkError::algorithm(format!(
            "PairHMM packed-f64 refused oversized DP (read_len={rn}, max_hap_len={hn}, cells={max_cells}); \
             inputs must be assembly-region scale, not contig scale"
        )));
    }
    Ok(PACK_F64_SCRATCH.with(|cell| {
        let mut scratch = cell.borrow_mut();
        score_one_f64_rolling(read_bases, read_quals, hap, transitions, &mut scratch)
    }))
}

/// Score haplotypes with already-built Logless transitions (avoids rebuild on NEON leftovers).
///
/// Same-length consecutive haplotypes reuse DP columns before the first divergence
/// (Java `LoglessPairHMM` `hapStartIndex` / `nextHapStartIndex` contract). Measured
/// faster than independent rolling on dense mega packs.
pub(crate) fn score_haps_logless_packed_f64_with_transitions(
    read_bases: &[u8],
    read_quals: &[u8],
    haplotypes: &[&[u8]],
    transitions: &[[f64; 6]],
) -> GatkResult<Vec<f64>> {
    if haplotypes.is_empty() {
        return Ok(Vec::new());
    }
    let rn = read_bases.len();
    if rn == 0 {
        return Ok(vec![0.0; haplotypes.len()]);
    }
    let max_hn = haplotypes.iter().map(|h| h.len()).max().unwrap_or(0);
    let max_cols = max_hn + 1;
    let max_cells = (rn + 1).saturating_mul(max_cols);
    const MAX_PAIRHMM_DIM: usize = 100_000;
    const MAX_PAIRHMM_CELLS: usize = 8_000_000;
    if rn > MAX_PAIRHMM_DIM || max_hn > MAX_PAIRHMM_DIM || max_cells > MAX_PAIRHMM_CELLS {
        return Err(GatkError::algorithm(format!(
            "PairHMM packed-f64 refused oversized DP (read_len={rn}, max_hap_len={max_hn}, cells={max_cells}); \
             inputs must be assembly-region scale, not contig scale"
        )));
    }
    let mut out = Vec::with_capacity(haplotypes.len());
    PACK_F64_SCRATCH.with(|cell| {
        let mut scratch = cell.borrow_mut();
        scratch.ensure_cells(max_cells);
        let mut prev: Option<&[u8]> = None;
        for &hap in haplotypes {
            let (hap_start, reinit_del) = match prev {
                Some(p) if p.len() == hap.len() => (first_hap_divergence(p, hap), false),
                _ => (0, true),
            };
            out.push(score_one_f64(
                read_bases,
                read_quals,
                hap,
                transitions,
                hap_start,
                reinit_del,
                &mut scratch,
            ));
            prev = Some(hap);
        }
    });
    Ok(out)
}

#[inline(always)]
fn logless_emission(x: u8, y: u8, match_p: f64, mismatch_p: f64) -> f64 {
    if x == y || x == b'N' || y == b'N' {
        match_p
    } else {
        mismatch_p
    }
}

/// One Logless DP read-row from column `j0` through `hn`.
///
/// M and I depend only on the previous row, so consecutive haplotype columns are
/// independent — GKL-style SIMD along the hap axis inside Java `hapStartIndex`.
/// D stays scalar (left-cell dependence on the current row).
///
/// SAFETY: `m`/`ins`/`del` cover `(row|prev)+hn`; `hap.len()==hn`; `j0>=1`.
unsafe fn fill_prefix_row(
    hap: &[u8],
    x: u8,
    match_p: f64,
    mismatch_p: f64,
    t: [f64; 6],
    m: &mut [f64],
    ins: &mut [f64],
    del: &mut [f64],
    row: usize,
    prev: usize,
    j0: usize,
    hn: usize,
) {
    #[cfg(target_arch = "aarch64")]
    {
        fill_prefix_row_neon(
            hap, x, match_p, mismatch_p, t, m, ins, del, row, prev, j0, hn,
        );
        return;
    }
    #[cfg(target_arch = "x86_64")]
    {
        fill_prefix_row_sse2(
            hap, x, match_p, mismatch_p, t, m, ins, del, row, prev, j0, hn,
        );
        return;
    }
    #[cfg(not(any(target_arch = "aarch64", target_arch = "x86_64")))]
    fill_prefix_row_scalar(
        hap, x, match_p, mismatch_p, t, m, ins, del, row, prev, j0, hn,
    );
}

unsafe fn fill_prefix_row_scalar(
    hap: &[u8],
    x: u8,
    match_p: f64,
    mismatch_p: f64,
    t: [f64; 6],
    m: &mut [f64],
    ins: &mut [f64],
    del: &mut [f64],
    row: usize,
    prev: usize,
    j0: usize,
    hn: usize,
) {
    for j in j0..=hn {
        let p = logless_emission(x, *hap.get_unchecked(j - 1), match_p, mismatch_p);
        let mj = row + j;
        let diag = prev + j - 1;
        let up = prev + j;
        *m.get_unchecked_mut(mj) = p
            * (*m.get_unchecked(diag) * t[MATCH_TO_MATCH]
                + *ins.get_unchecked(diag) * t[INDEL_TO_MATCH]
                + *del.get_unchecked(diag) * t[INDEL_TO_MATCH]);
        *ins.get_unchecked_mut(mj) = *m.get_unchecked(up) * t[MATCH_TO_INSERTION]
            + *ins.get_unchecked(up) * t[INSERTION_TO_INSERTION];
        *del.get_unchecked_mut(mj) = *m.get_unchecked(mj - 1) * t[MATCH_TO_DELETION]
            + *del.get_unchecked(mj - 1) * t[DELETION_TO_DELETION];
    }
}

/// 4-wide then 2-wide M/I via NEON; D remains a scalar left-to-right sweep.
#[cfg(target_arch = "aarch64")]
unsafe fn fill_prefix_row_neon(
    hap: &[u8],
    x: u8,
    match_p: f64,
    mismatch_p: f64,
    t: [f64; 6],
    m: &mut [f64],
    ins: &mut [f64],
    del: &mut [f64],
    row: usize,
    prev: usize,
    j0: usize,
    hn: usize,
) {
    use std::arch::aarch64::*;
    let tm2m = vdupq_n_f64(t[MATCH_TO_MATCH]);
    let ti2m = vdupq_n_f64(t[INDEL_TO_MATCH]);
    let tm2i = vdupq_n_f64(t[MATCH_TO_INSERTION]);
    let ti2i = vdupq_n_f64(t[INSERTION_TO_INSERTION]);
    let tm2d = t[MATCH_TO_DELETION];
    let td2d = t[DELETION_TO_DELETION];
    let mut j = j0;
    while j + 3 <= hn {
        let p01 = [
            logless_emission(x, *hap.get_unchecked(j - 1), match_p, mismatch_p),
            logless_emission(x, *hap.get_unchecked(j), match_p, mismatch_p),
        ];
        let p23 = [
            logless_emission(x, *hap.get_unchecked(j + 1), match_p, mismatch_p),
            logless_emission(x, *hap.get_unchecked(j + 2), match_p, mismatch_p),
        ];
        let diag = prev + j - 1;
        let up = prev + j;
        let mj = row + j;
        let m_d0 = vld1q_f64(m.as_ptr().add(diag));
        let m_d1 = vld1q_f64(m.as_ptr().add(diag + 2));
        let i_d0 = vld1q_f64(ins.as_ptr().add(diag));
        let i_d1 = vld1q_f64(ins.as_ptr().add(diag + 2));
        let d_d0 = vld1q_f64(del.as_ptr().add(diag));
        let d_d1 = vld1q_f64(del.as_ptr().add(diag + 2));
        let m_u0 = vld1q_f64(m.as_ptr().add(up));
        let m_u1 = vld1q_f64(m.as_ptr().add(up + 2));
        let i_u0 = vld1q_f64(ins.as_ptr().add(up));
        let i_u1 = vld1q_f64(ins.as_ptr().add(up + 2));
        let sum0 = vaddq_f64(
            vaddq_f64(vmulq_f64(m_d0, tm2m), vmulq_f64(i_d0, ti2m)),
            vmulq_f64(d_d0, ti2m),
        );
        let sum1 = vaddq_f64(
            vaddq_f64(vmulq_f64(m_d1, tm2m), vmulq_f64(i_d1, ti2m)),
            vmulq_f64(d_d1, ti2m),
        );
        vst1q_f64(
            m.as_mut_ptr().add(mj),
            vmulq_f64(vld1q_f64(p01.as_ptr()), sum0),
        );
        vst1q_f64(
            m.as_mut_ptr().add(mj + 2),
            vmulq_f64(vld1q_f64(p23.as_ptr()), sum1),
        );
        vst1q_f64(
            ins.as_mut_ptr().add(mj),
            vaddq_f64(vmulq_f64(m_u0, tm2i), vmulq_f64(i_u0, ti2i)),
        );
        vst1q_f64(
            ins.as_mut_ptr().add(mj + 2),
            vaddq_f64(vmulq_f64(m_u1, tm2i), vmulq_f64(i_u1, ti2i)),
        );
        *del.get_unchecked_mut(mj) =
            *m.get_unchecked(mj - 1) * tm2d + *del.get_unchecked(mj - 1) * td2d;
        *del.get_unchecked_mut(mj + 1) =
            *m.get_unchecked(mj) * tm2d + *del.get_unchecked(mj) * td2d;
        *del.get_unchecked_mut(mj + 2) =
            *m.get_unchecked(mj + 1) * tm2d + *del.get_unchecked(mj + 1) * td2d;
        *del.get_unchecked_mut(mj + 3) =
            *m.get_unchecked(mj + 2) * tm2d + *del.get_unchecked(mj + 2) * td2d;
        j += 4;
    }
    while j + 1 <= hn {
        let p0 = logless_emission(x, *hap.get_unchecked(j - 1), match_p, mismatch_p);
        let p1 = logless_emission(x, *hap.get_unchecked(j), match_p, mismatch_p);
        let priors = [p0, p1];
        let p = vld1q_f64(priors.as_ptr());
        let diag = prev + j - 1;
        let up = prev + j;
        let mj = row + j;
        let m_diag = vld1q_f64(m.as_ptr().add(diag));
        let i_diag = vld1q_f64(ins.as_ptr().add(diag));
        let d_diag = vld1q_f64(del.as_ptr().add(diag));
        let m_up = vld1q_f64(m.as_ptr().add(up));
        let i_up = vld1q_f64(ins.as_ptr().add(up));
        let sum = vaddq_f64(
            vaddq_f64(vmulq_f64(m_diag, tm2m), vmulq_f64(i_diag, ti2m)),
            vmulq_f64(d_diag, ti2m),
        );
        vst1q_f64(m.as_mut_ptr().add(mj), vmulq_f64(p, sum));
        vst1q_f64(
            ins.as_mut_ptr().add(mj),
            vaddq_f64(vmulq_f64(m_up, tm2i), vmulq_f64(i_up, ti2i)),
        );
        *del.get_unchecked_mut(mj) =
            *m.get_unchecked(mj - 1) * tm2d + *del.get_unchecked(mj - 1) * td2d;
        *del.get_unchecked_mut(mj + 1) =
            *m.get_unchecked(mj) * tm2d + *del.get_unchecked(mj) * td2d;
        j += 2;
    }
    if j <= hn {
        fill_prefix_row_scalar(
            hap, x, match_p, mismatch_p, t, m, ins, del, row, prev, j, hn,
        );
    }
}

/// 4-wide then 2-wide M/I via SSE2 (x86_64 baseline); D remains scalar.
#[cfg(target_arch = "x86_64")]
unsafe fn fill_prefix_row_sse2(
    hap: &[u8],
    x: u8,
    match_p: f64,
    mismatch_p: f64,
    t: [f64; 6],
    m: &mut [f64],
    ins: &mut [f64],
    del: &mut [f64],
    row: usize,
    prev: usize,
    j0: usize,
    hn: usize,
) {
    use std::arch::x86_64::*;
    let tm2m = _mm_set1_pd(t[MATCH_TO_MATCH]);
    let ti2m = _mm_set1_pd(t[INDEL_TO_MATCH]);
    let tm2i = _mm_set1_pd(t[MATCH_TO_INSERTION]);
    let ti2i = _mm_set1_pd(t[INSERTION_TO_INSERTION]);
    let tm2d = t[MATCH_TO_DELETION];
    let td2d = t[DELETION_TO_DELETION];
    let mut j = j0;
    while j + 3 <= hn {
        let p01 = [
            logless_emission(x, *hap.get_unchecked(j - 1), match_p, mismatch_p),
            logless_emission(x, *hap.get_unchecked(j), match_p, mismatch_p),
        ];
        let p23 = [
            logless_emission(x, *hap.get_unchecked(j + 1), match_p, mismatch_p),
            logless_emission(x, *hap.get_unchecked(j + 2), match_p, mismatch_p),
        ];
        let diag = prev + j - 1;
        let up = prev + j;
        let mj = row + j;
        let m_d0 = _mm_loadu_pd(m.as_ptr().add(diag));
        let m_d1 = _mm_loadu_pd(m.as_ptr().add(diag + 2));
        let i_d0 = _mm_loadu_pd(ins.as_ptr().add(diag));
        let i_d1 = _mm_loadu_pd(ins.as_ptr().add(diag + 2));
        let d_d0 = _mm_loadu_pd(del.as_ptr().add(diag));
        let d_d1 = _mm_loadu_pd(del.as_ptr().add(diag + 2));
        let m_u0 = _mm_loadu_pd(m.as_ptr().add(up));
        let m_u1 = _mm_loadu_pd(m.as_ptr().add(up + 2));
        let i_u0 = _mm_loadu_pd(ins.as_ptr().add(up));
        let i_u1 = _mm_loadu_pd(ins.as_ptr().add(up + 2));
        let sum0 = _mm_add_pd(
            _mm_add_pd(_mm_mul_pd(m_d0, tm2m), _mm_mul_pd(i_d0, ti2m)),
            _mm_mul_pd(d_d0, ti2m),
        );
        let sum1 = _mm_add_pd(
            _mm_add_pd(_mm_mul_pd(m_d1, tm2m), _mm_mul_pd(i_d1, ti2m)),
            _mm_mul_pd(d_d1, ti2m),
        );
        _mm_storeu_pd(
            m.as_mut_ptr().add(mj),
            _mm_mul_pd(_mm_loadu_pd(p01.as_ptr()), sum0),
        );
        _mm_storeu_pd(
            m.as_mut_ptr().add(mj + 2),
            _mm_mul_pd(_mm_loadu_pd(p23.as_ptr()), sum1),
        );
        _mm_storeu_pd(
            ins.as_mut_ptr().add(mj),
            _mm_add_pd(_mm_mul_pd(m_u0, tm2i), _mm_mul_pd(i_u0, ti2i)),
        );
        _mm_storeu_pd(
            ins.as_mut_ptr().add(mj + 2),
            _mm_add_pd(_mm_mul_pd(m_u1, tm2i), _mm_mul_pd(i_u1, ti2i)),
        );
        *del.get_unchecked_mut(mj) =
            *m.get_unchecked(mj - 1) * tm2d + *del.get_unchecked(mj - 1) * td2d;
        *del.get_unchecked_mut(mj + 1) =
            *m.get_unchecked(mj) * tm2d + *del.get_unchecked(mj) * td2d;
        *del.get_unchecked_mut(mj + 2) =
            *m.get_unchecked(mj + 1) * tm2d + *del.get_unchecked(mj + 1) * td2d;
        *del.get_unchecked_mut(mj + 3) =
            *m.get_unchecked(mj + 2) * tm2d + *del.get_unchecked(mj + 2) * td2d;
        j += 4;
    }
    while j + 1 <= hn {
        let p0 = logless_emission(x, *hap.get_unchecked(j - 1), match_p, mismatch_p);
        let p1 = logless_emission(x, *hap.get_unchecked(j), match_p, mismatch_p);
        let priors = [p0, p1];
        let p = _mm_loadu_pd(priors.as_ptr());
        let diag = prev + j - 1;
        let up = prev + j;
        let mj = row + j;
        let m_diag = _mm_loadu_pd(m.as_ptr().add(diag));
        let i_diag = _mm_loadu_pd(ins.as_ptr().add(diag));
        let d_diag = _mm_loadu_pd(del.as_ptr().add(diag));
        let m_up = _mm_loadu_pd(m.as_ptr().add(up));
        let i_up = _mm_loadu_pd(ins.as_ptr().add(up));
        let sum = _mm_add_pd(
            _mm_add_pd(_mm_mul_pd(m_diag, tm2m), _mm_mul_pd(i_diag, ti2m)),
            _mm_mul_pd(d_diag, ti2m),
        );
        _mm_storeu_pd(m.as_mut_ptr().add(mj), _mm_mul_pd(p, sum));
        _mm_storeu_pd(
            ins.as_mut_ptr().add(mj),
            _mm_add_pd(_mm_mul_pd(m_up, tm2i), _mm_mul_pd(i_up, ti2i)),
        );
        *del.get_unchecked_mut(mj) =
            *m.get_unchecked(mj - 1) * tm2d + *del.get_unchecked(mj - 1) * td2d;
        *del.get_unchecked_mut(mj + 1) =
            *m.get_unchecked(mj) * tm2d + *del.get_unchecked(mj) * td2d;
        j += 2;
    }
    if j <= hn {
        fill_prefix_row_scalar(
            hap, x, match_p, mismatch_p, t, m, ins, del, row, prev, j, hn,
        );
    }
}

/// Full-matrix path for Java `hapStartIndex` prefix reuse (columns `< start` kept).
fn score_one_f64_prefix_reuse(
    read_bases: &[u8],
    read_quals: &[u8],
    hap: &[u8],
    transitions: &[[f64; 6]],
    start: usize,
    reinit_del: bool,
    scratch: &mut F64Scratch,
) -> f64 {
    let rn = read_bases.len();
    let hn = hap.len();
    let cols = hn + 1;
    let m = &mut scratch.m;
    let ins = &mut scratch.ins;
    let del = &mut scratch.del;

    if reinit_del {
        // Fresh hap: do NOT memset the full rn×hn planes — every used cell is overwritten.
        // Only seed row-0 free deletions + keep col-0 zeros via per-row writes below.
        let init_del = INITIAL_CONDITION / hn as f64;
        m[..cols].fill(0.0);
        ins[..cols].fill(0.0);
        for j in 0..=hn {
            del[j] = init_del;
        }
    }

    let j0 = start + 1; // 1-based DP column; Java uses hapStartIndex+1
    for i in 1..=rn {
        let t = transitions[i];
        let x = read_bases[i - 1];
        let (match_p, mismatch_p) = logless_match_mismatch_prior(read_quals[i - 1]);
        let row = i * cols;
        let prev = (i - 1) * cols;
        // Col 0 stays 0 (never written in the j-loop); required for del[j=1] left term.
        m[row] = 0.0;
        ins[row] = 0.0;
        del[row] = 0.0;
        // SAFETY: row = i * cols, cols = hn+1, scratch sized to (rn+1)*cols; j in j0..=hn.
        unsafe {
            fill_prefix_row(
                hap, x, match_p, mismatch_p, t, m, ins, del, row, prev, j0, hn,
            );
        }
    }

    let end_row = rn * cols;
    let mut final_sum = 0.0;
    for j in 1..=hn {
        final_sum += m[end_row + j] + ins[end_row + j];
    }
    if final_sum <= 0.0 || !final_sum.is_finite() {
        return f64::NEG_INFINITY;
    }
    final_sum.log10() - INITIAL_CONDITION_LOG10
}

/// Compact DP: only previous + current read rows (priors computed inline).
fn score_one_f64_rolling(
    read_bases: &[u8],
    read_quals: &[u8],
    hap: &[u8],
    transitions: &[[f64; 6]],
    scratch: &mut F64Scratch,
) -> f64 {
    let rn = read_bases.len();
    let hn = hap.len();
    let cols = hn + 1;
    scratch.ensure_rolling_cols(cols);
    // Layout: [0..cols) = prev row, [cols..2*cols) = curr row for m/ins/del.
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
        let t = transitions[i];
        let x = read_bases[i - 1];
        let (match_p, mismatch_p) = logless_match_mismatch_prior(read_quals[i - 1]);
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
    let mut final_sum = 0.0;
    for j in 1..=hn {
        final_sum += m_end[j] + ins_end[j];
    }
    if final_sum <= 0.0 || !final_sum.is_finite() {
        return f64::NEG_INFINITY;
    }
    final_sum.log10() - INITIAL_CONDITION_LOG10
}

struct F32Scratch {
    m: Vec<f32>,
    ins: Vec<f32>,
    del: Vec<f32>,
    prior: Vec<f32>,
}

impl F32Scratch {
    fn empty() -> Self {
        Self {
            m: Vec::new(),
            ins: Vec::new(),
            del: Vec::new(),
            prior: Vec::new(),
        }
    }

    fn ensure_cells(&mut self, cells: usize) {
        if self.m.len() < cells {
            self.m.resize(cells, 0.0);
            self.ins.resize(cells, 0.0);
            self.del.resize(cells, 0.0);
            self.prior.resize(cells, 0.0);
        }
    }

    fn clear_prefix(&mut self, cells: usize) {
        self.ensure_cells(cells);
        self.m[..cells].fill(0.0);
        self.ins[..cells].fill(0.0);
        self.del[..cells].fill(0.0);
        self.prior[..cells].fill(0.0);
    }
}

/// f32 packed path with per-haplotype f64 retry when the linear sum underflows.
pub fn score_haps_logless_packed_f32(
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

    let mut transitions_f32 = vec![[0.0f32; 6]; rn + 1];
    for i in 0..rn {
        let t = logless_qual_to_trans_probs(insertion_gop[i], deletion_gop[i], overall_gcp[i]);
        transitions_f32[i + 1] = [
            t[0] as f32,
            t[1] as f32,
            t[2] as f32,
            t[3] as f32,
            t[4] as f32,
            t[5] as f32,
        ];
    }

    let max_hn = haplotypes.iter().map(|h| h.len()).max().unwrap_or(0);
    let max_cells = (rn + 1).saturating_mul(max_hn + 1);
    const MAX_PAIRHMM_DIM: usize = 100_000;
    const MAX_PAIRHMM_CELLS: usize = 8_000_000;
    if rn > MAX_PAIRHMM_DIM || max_hn > MAX_PAIRHMM_DIM || max_cells > MAX_PAIRHMM_CELLS {
        return Err(GatkError::algorithm(format!(
            "PairHMM packed-f32 refused oversized DP (read_len={rn}, max_hap_len={max_hn}, cells={max_cells}); \
             inputs must be assembly-region scale, not contig scale"
        )));
    }
    let mut out = Vec::with_capacity(haplotypes.len());
    let mut err = None;
    PACK_F32_SCRATCH.with(|cell| {
        let mut scratch = cell.borrow_mut();
        scratch.ensure_cells(max_cells);
        for &hap in haplotypes {
            let (ll, linear_sum) =
                score_one_f32(read_bases, read_quals, hap, &transitions_f32, &mut scratch);
            if !linear_sum.is_finite() || (linear_sum as f64) < MIN_ACCEPTED_LINEAR_SUM {
                match logless_pairhmm_likelihood(
                    read_bases,
                    read_quals,
                    hap,
                    insertion_gop,
                    deletion_gop,
                    overall_gcp,
                ) {
                    Ok(v) => out.push(v),
                    Err(e) => {
                        err = Some(e);
                        break;
                    }
                }
            } else {
                out.push(ll);
            }
        }
    });
    if let Some(e) = err {
        return Err(e);
    }
    Ok(out)
}

fn score_one_f32(
    read_bases: &[u8],
    read_quals: &[u8],
    hap: &[u8],
    transitions: &[[f32; 6]],
    scratch: &mut F32Scratch,
) -> (f64, f32) {
    let rn = read_bases.len();
    let hn = hap.len();
    let cols = hn + 1;
    let cells = (rn + 1) * cols;
    scratch.clear_prefix(cells);
    let m = &mut scratch.m;
    let ins = &mut scratch.ins;
    let del = &mut scratch.del;
    let prior = &mut scratch.prior;

    for i in 0..rn {
        let x = read_bases[i];
        let (match_p, mismatch_p) = logless_match_mismatch_prior(read_quals[i]);
        let match_p = match_p as f32;
        let mismatch_p = mismatch_p as f32;
        let row = (i + 1) * cols;
        for j in 0..hn {
            let y = hap[j];
            prior[row + j + 1] = if x == y || x == b'N' || y == b'N' {
                match_p
            } else {
                mismatch_p
            };
        }
    }

    let init_del = (INITIAL_CONDITION / hn as f64) as f32;
    for j in 0..=hn {
        del[j] = init_del;
    }

    for i in 1..=rn {
        let t = transitions[i];
        let row = i * cols;
        let prev = (i - 1) * cols;
        for j in 1..=hn {
            let p = prior[row + j];
            m[row + j] = p
                * (m[prev + j - 1] * t[MATCH_TO_MATCH]
                    + ins[prev + j - 1] * t[INDEL_TO_MATCH]
                    + del[prev + j - 1] * t[INDEL_TO_MATCH]);
            ins[row + j] =
                m[prev + j] * t[MATCH_TO_INSERTION] + ins[prev + j] * t[INSERTION_TO_INSERTION];
            del[row + j] =
                m[row + j - 1] * t[MATCH_TO_DELETION] + del[row + j - 1] * t[DELETION_TO_DELETION];
        }
    }

    let end_row = rn * cols;
    let mut final_sum = 0.0f32;
    for j in 1..=hn {
        final_sum += m[end_row + j] + ins[end_row + j];
    }
    if final_sum <= 0.0 || !final_sum.is_finite() {
        return (f64::NEG_INFINITY, final_sum);
    }
    let ll = (final_sum as f64).log10() - INITIAL_CONDITION_LOG10;
    (ll, final_sum)
}
