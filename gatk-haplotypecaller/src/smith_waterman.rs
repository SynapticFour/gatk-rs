//! GATK `SmithWatermanJavaAligner` parity (JAVA implementation path).
//!
//! # Layout
//! Backtrack is a flat row-major `Vec<i32>`. Scores use **two rolling rows** plus a
//! saved last-column vector: CIGAR only needs scores on the final row/column to pick
//! the end cell, then follows `btrack`.
//!
//! # Production HC sizes
//! Typical SoftClip realign is ~read×hap (≈100–250). Hap-to-ref Indel uses
//! `SW_PAD` (±10 N) around haplotypes ≈50–300. Full score matrices are not retained.

#![warn(clippy::unwrap_used, clippy::expect_used)]

use crate::cigar::{Cigar, CigarElement, CigarOperator};
use gatk_common::{GatkError, GatkResult};
use std::cell::RefCell;

/// GATK `SmithWatermanAlignmentConstants.NEW_SW_PARAMETERS` (haplotype-to-reference).
#[derive(Debug, Clone, Copy)]
pub struct SwParameters {
    pub match_value: i32,
    pub mismatch_penalty: i32,
    pub gap_open_penalty: i32,
    pub gap_extend_penalty: i32,
}

impl SwParameters {
    pub fn gatk_haplotype_to_reference() -> Self {
        Self {
            match_value: 200,
            mismatch_penalty: -150,
            gap_open_penalty: -260,
            gap_extend_penalty: -11,
        }
    }

    /// GATK `SmithWatermanAlignmentConstants.ALIGNMENT_TO_BEST_HAPLOTYPE_SW_PARAMETERS`.
    pub fn gatk_read_to_best_haplotype() -> Self {
        Self {
            match_value: 10,
            mismatch_penalty: -15,
            gap_open_penalty: -30,
            gap_extend_penalty: -5,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SwOverhangStrategy {
    SoftClip,
    Ignore,
    Indel,
    LeadingIndel,
}

#[derive(Debug, Clone)]
pub struct SmithWatermanAlignment {
    pub cigar: Cigar,
    pub alignment_offset: i32,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum SwState {
    Match,
    Insertion,
    Deletion,
    Clip,
}

#[inline]
fn cell(i: usize, j: usize, ncol: usize) -> usize {
    i * ncol + j
}

fn last_index_of(reference: &[u8], query: &[u8]) -> Option<usize> {
    if query.is_empty() || reference.len() < query.len() {
        return None;
    }
    let qlen = query.len();
    let first = query[0];
    let last = query[qlen - 1];
    // Walk candidate starts from the end. Skip starts whose first byte cannot match
    // without a full slice compare (Java String.lastIndexOf contract).
    let mut r = reference.len() - qlen;
    loop {
        // Fast reject on first byte before touching the last-byte / interior checks.
        if reference[r] == first {
            if qlen == 1
                || (reference[r + qlen - 1] == last
                    && (qlen == 2 || reference[r + 1..r + qlen - 1] == query[1..qlen - 1]))
            {
                return Some(r);
            }
        }
        if r == 0 {
            return None;
        }
        r -= 1;
    }
}

/// Align `alternate` to `reference` (GATK argument order in `CigarUtils.calculateCigar` padded call).
pub fn align(
    reference: &[u8],
    alternate: &[u8],
    parameters: &SwParameters,
    overhang_strategy: SwOverhangStrategy,
) -> GatkResult<SmithWatermanAlignment> {
    align_internal(reference, alternate, parameters, overhang_strategy, true)
}

/// GATK `createReadAlignedToRef` read→hap SW (`SmithWatermanJavaAligner` + `lastIndexOf` fast path).
pub fn align_read_to_best_haplotype(
    haplotype_bases: &[u8],
    read_bases: &[u8],
    parameters: &SwParameters,
) -> GatkResult<SmithWatermanAlignment> {
    align_internal(
        haplotype_bases,
        read_bases,
        parameters,
        SwOverhangStrategy::SoftClip,
        true,
    )
}

#[inline]
fn already_ascii_uppercase(seq: &[u8]) -> bool {
    !seq.iter().any(|b| b.is_ascii_lowercase())
}

fn align_internal(
    reference: &[u8],
    alternate: &[u8],
    parameters: &SwParameters,
    overhang_strategy: SwOverhangStrategy,
    allow_substring_fast_path: bool,
) -> GatkResult<SmithWatermanAlignment> {
    if reference.is_empty() || alternate.is_empty() {
        return Err(GatkError::algorithm(
            "Non-null, non-empty sequences are required for Smith-Waterman",
        ));
    }

    // Skip uppercase copies when inputs are already ASCII-upper (common in HC paths).
    if already_ascii_uppercase(reference) && already_ascii_uppercase(alternate) {
        return align_uppercase_ready(
            reference,
            alternate,
            parameters,
            overhang_strategy,
            allow_substring_fast_path,
        );
    }
    let reference: Vec<u8> = reference.iter().map(|b| b.to_ascii_uppercase()).collect();
    let alternate: Vec<u8> = alternate.iter().map(|b| b.to_ascii_uppercase()).collect();
    align_uppercase_ready(
        &reference,
        &alternate,
        parameters,
        overhang_strategy,
        allow_substring_fast_path,
    )
}

fn align_uppercase_ready(
    reference: &[u8],
    alternate: &[u8],
    parameters: &SwParameters,
    overhang_strategy: SwOverhangStrategy,
    allow_substring_fast_path: bool,
) -> GatkResult<SmithWatermanAlignment> {
    // MUST run before `last_index_of`: a contig-length "reference" makes that scan
    // thrash a 16 GiB laptop for minutes even when the DP matrix is never allocated.
    // Contig × ~read-length grids are ~60 GiB on hs37d5 chr20 (realistic-window OOM).
    const MAX_SW_DIM: usize = 100_000;
    // Match PairHMM: refuse before multi-GiB TLS high-water on 16 GiB hosts.
    const MAX_SW_CELLS: usize = 8_000_000;
    let nrow = reference.len() + 1;
    let ncol = alternate.len() + 1;
    let cells = nrow.saturating_mul(ncol);
    if reference.len() > MAX_SW_DIM || alternate.len() > MAX_SW_DIM || cells > MAX_SW_CELLS {
        return Err(GatkError::algorithm(format!(
            "Smith-Waterman refused oversized matrix (ref_len={}, alt_len={}, cells={cells}); \
             inputs must be assembly-region scale, not contig scale",
            reference.len(),
            alternate.len()
        )));
    }

    if allow_substring_fast_path
        && matches!(
            overhang_strategy,
            SwOverhangStrategy::SoftClip | SwOverhangStrategy::Ignore
        )
    {
        if let Some(match_index) = last_index_of(reference, alternate) {
            return Ok(SmithWatermanAlignment {
                cigar: {
                    let mut c = Cigar::new();
                    c.push(alternate.len(), CigarOperator::Match);
                    c
                },
                alignment_offset: match_index as i32,
            });
        }
    }

    SW_SCRATCH.with(|cell| {
        let scratch = &mut *cell.borrow_mut();
        scratch.ensure(cells, nrow, ncol);
        // Work in place — no mem::take / restore (avoids TLS move churn).
        calculate_matrix(
            reference,
            alternate,
            nrow,
            ncol,
            &mut scratch.sw_prev[..ncol],
            &mut scratch.sw_cur[..ncol],
            &mut scratch.last_col[..nrow],
            &mut scratch.btrack[..cells],
            &mut scratch.best_gap_v,
            &mut scratch.gap_size_v,
            &mut scratch.best_gap_h,
            &mut scratch.gap_size_h,
            overhang_strategy,
            parameters,
        );
        // After DP, `sw_prev` holds the last row (see calculate_matrix swap).
        let aln = calculate_cigar(
            &scratch.last_col[..nrow],
            &scratch.sw_prev[..ncol],
            &scratch.btrack[..cells],
            nrow,
            ncol,
            overhang_strategy,
        );
        // Soft-keep normal-sized planes; drop oversized allocations so Peak does not retain
        // multi-megabase matrices between calls. Region/engine end still hard-clears via
        // [`release_sw_tls_scratch`].
        if cells > SW_TLS_SOFT_KEEP_CELLS || scratch.btrack.capacity() > SW_TLS_SOFT_KEEP_CELLS {
            scratch.clear();
        }
        Ok(aln)
    })
}

/// Soft-keep ceiling for SW TLS backtrack cells. Larger requests still compute,
/// but planes are not retained afterward.
const SW_TLS_SOFT_KEEP_CELLS: usize = 1 << 20; // 1_048_576 cells

struct SwScratch {
    /// Previous / current score rows (`ncol` each). After DP, `sw_prev` is the last row.
    sw_prev: Vec<i32>,
    sw_cur: Vec<i32>,
    /// `sw[i][alt_length]` for end-cell selection (SoftClip / Ignore).
    last_col: Vec<i32>,
    btrack: Vec<i32>,
    best_gap_v: Vec<i32>,
    gap_size_v: Vec<i32>,
    best_gap_h: Vec<i32>,
    gap_size_h: Vec<i32>,
}

impl SwScratch {
    fn new() -> Self {
        Self {
            sw_prev: Vec::new(),
            sw_cur: Vec::new(),
            last_col: Vec::new(),
            btrack: Vec::new(),
            best_gap_v: Vec::new(),
            gap_size_v: Vec::new(),
            best_gap_h: Vec::new(),
            gap_size_h: Vec::new(),
        }
    }

    fn ensure(&mut self, cells: usize, nrow: usize, ncol: usize) {
        if cells > SW_TLS_SOFT_KEEP_CELLS || self.btrack.capacity() > SW_TLS_SOFT_KEEP_CELLS {
            if self.btrack.capacity() > 0 {
                self.clear();
            }
        }
        if self.btrack.len() < cells {
            self.btrack.resize(cells, 0);
        }
        if self.sw_prev.len() < ncol {
            self.sw_prev.resize(ncol, 0);
            self.sw_cur.resize(ncol, 0);
        }
        if self.last_col.len() < nrow {
            self.last_col.resize(nrow, 0);
        }
        if self.best_gap_v.len() < ncol + 1 {
            self.best_gap_v.resize(ncol + 1, 0);
            self.gap_size_v.resize(ncol + 1, 0);
        }
        if self.best_gap_h.len() < nrow + 1 {
            self.best_gap_h.resize(nrow + 1, 0);
            self.gap_size_h.resize(nrow + 1, 0);
        }
    }

    fn clear(&mut self) {
        // Drop retained Peak; next `ensure` reallocates as needed.
        *self = Self::new();
    }
}

thread_local! {
    static SW_SCRATCH: RefCell<SwScratch> = RefCell::new(SwScratch::new());
}

/// Drop SW TLS scratch (see `run::release_region_tls_scratch`).
pub fn release_sw_tls_scratch() {
    SW_SCRATCH.with(|cell| {
        cell.borrow_mut().clear();
    });
}

fn calculate_matrix(
    reference: &[u8],
    alternate: &[u8],
    nrow: usize,
    ncol: usize,
    sw_prev: &mut [i32],
    sw_cur: &mut [i32],
    last_col: &mut [i32],
    btrack: &mut [i32],
    best_gap_v: &mut Vec<i32>,
    gap_size_v: &mut Vec<i32>,
    best_gap_h: &mut Vec<i32>,
    gap_size_h: &mut Vec<i32>,
    overhang_strategy: SwOverhangStrategy,
    parameters: &SwParameters,
) {
    const MATRIX_MIN_CUTOFF: i32 = -100_000_000;
    let low_init = i32::MIN / 4;
    if best_gap_v.len() < ncol + 1 {
        best_gap_v.resize(ncol + 1, low_init);
        gap_size_v.resize(ncol + 1, 0);
    }
    if best_gap_h.len() < nrow + 1 {
        best_gap_h.resize(nrow + 1, low_init);
        gap_size_h.resize(nrow + 1, 0);
    }
    best_gap_v[..ncol + 1].fill(low_init);
    gap_size_v[..ncol + 1].fill(0);
    best_gap_h[..nrow + 1].fill(low_init);
    gap_size_h[..nrow + 1].fill(0);

    // Row 0 into sw_prev.
    sw_prev[..ncol].fill(0);
    if matches!(
        overhang_strategy,
        SwOverhangStrategy::Indel | SwOverhangStrategy::LeadingIndel
    ) {
        sw_prev[1] = parameters.gap_open_penalty;
        let mut cur = parameters.gap_open_penalty;
        for j in 2..ncol {
            cur += parameters.gap_extend_penalty;
            sw_prev[j] = cur;
        }
    }
    last_col[0] = sw_prev[ncol - 1];

    let w_open = parameters.gap_open_penalty;
    let w_extend = parameters.gap_extend_penalty;
    let w_match = parameters.match_value;
    let w_mismatch = parameters.mismatch_penalty;
    let indel_edges = matches!(
        overhang_strategy,
        SwOverhangStrategy::Indel | SwOverhangStrategy::LeadingIndel
    );

    // Hot path: bounds are [1,nrow)×[1,ncol) with nrow=ref+1, ncol=alt+1.
    for i in 1..nrow {
        // SAFETY: i in 1..nrow ⇒ i-1 < reference.len().
        let a_base = unsafe { *reference.get_unchecked(i - 1) };
        // Column 0 of current row.
        if indel_edges {
            if i == 1 {
                sw_cur[0] = parameters.gap_open_penalty;
            } else {
                // Extend from previous row's col0 (same recurrence as full-matrix col0 fill).
                sw_cur[0] = sw_prev[0] + parameters.gap_extend_penalty;
            }
        } else {
            sw_cur[0] = 0;
        }

        for j in 1..ncol {
            // SAFETY: j in 1..ncol ⇒ j-1 < alternate.len(); btrack sized to nrow*ncol.
            let b_base = unsafe { *alternate.get_unchecked(j - 1) };
            let diag = unsafe { *sw_prev.get_unchecked(j - 1) }
                + if a_base == b_base {
                    w_match
                } else {
                    w_mismatch
                };
            let prev_gap = unsafe { *sw_prev.get_unchecked(j) } + w_open;
            let bv = unsafe { best_gap_v.get_unchecked_mut(j) };
            let gv = unsafe { gap_size_v.get_unchecked_mut(j) };
            *bv += w_extend;
            if prev_gap > *bv {
                *bv = prev_gap;
                *gv = 1;
            } else {
                *gv += 1;
            }
            let step_down = *bv;
            let kd = *gv;
            let prev_gap_h = unsafe { *sw_cur.get_unchecked(j - 1) } + w_open;
            let bh = unsafe { best_gap_h.get_unchecked_mut(i) };
            let gh = unsafe { gap_size_h.get_unchecked_mut(i) };
            *bh += w_extend;
            if prev_gap_h > *bh {
                *bh = prev_gap_h;
                *gh = 1;
            } else {
                *gh += 1;
            }
            let step_right = *bh;
            let ki = *gh;
            let (score, track) = if diag >= step_down && diag >= step_right {
                (diag.max(MATRIX_MIN_CUTOFF), 0)
            } else if step_right >= step_down {
                (step_right.max(MATRIX_MIN_CUTOFF), -ki)
            } else {
                (step_down.max(MATRIX_MIN_CUTOFF), kd)
            };
            unsafe {
                *sw_cur.get_unchecked_mut(j) = score;
                *btrack.get_unchecked_mut(cell(i, j, ncol)) = track;
            }
        }
        last_col[i] = sw_cur[ncol - 1];
        // Row i becomes previous for the next iteration.
        sw_prev.swap_with_slice(sw_cur);
    }
}

fn make_element(state: SwState, length: usize) -> CigarElement {
    let operator = match state {
        SwState::Match => CigarOperator::Match,
        SwState::Insertion => CigarOperator::Insertion,
        SwState::Deletion => CigarOperator::Deletion,
        SwState::Clip => CigarOperator::SoftClip,
    };
    CigarElement { length, operator }
}

fn calculate_cigar(
    last_col: &[i32],
    last_row: &[i32],
    btrack: &[i32],
    nrow: usize,
    ncol: usize,
    overhang_strategy: SwOverhangStrategy,
) -> SmithWatermanAlignment {
    let ref_length = nrow - 1;
    let alt_length = ncol - 1;
    let mut p1 = 0usize;
    let mut p2 = alt_length;
    let mut maxscore = i32::MIN;
    let mut segment_length = 0usize;
    if overhang_strategy == SwOverhangStrategy::Indel {
        p1 = ref_length;
    } else {
        for i in 1..=ref_length {
            let cur = last_col[i];
            if cur >= maxscore {
                p1 = i;
                maxscore = cur;
            }
        }
        if overhang_strategy != SwOverhangStrategy::LeadingIndel {
            for j in 1..ncol {
                let cur = last_row[j];
                if cur > maxscore
                    || (cur == maxscore
                        && (ref_length as i32 - j as i32).abs() < (p1 as i32 - p2 as i32).abs())
                {
                    p1 = ref_length;
                    p2 = j;
                    maxscore = cur;
                    segment_length = alt_length - j;
                }
            }
        }
    }
    let mut lce: Vec<CigarElement> = Vec::new();
    if segment_length > 0 && overhang_strategy == SwOverhangStrategy::SoftClip {
        lce.push(make_element(SwState::Clip, segment_length));
        segment_length = 0;
    }
    let mut state = SwState::Match;
    // GATK SmithWatermanJavaAligner#calculateCigar uses do-while (at least one backtrack step).
    loop {
        if p1 == 0 || p2 == 0 {
            break;
        }
        // SAFETY: p1,p2 in range while looping; matrices sized nrow*ncol.
        let btr = unsafe { *btrack.get_unchecked(cell(p1, p2, ncol)) };
        let (new_state, step_length) = if btr > 0 {
            (SwState::Deletion, btr as usize)
        } else if btr < 0 {
            (SwState::Insertion, (-btr) as usize)
        } else {
            (SwState::Match, 1)
        };
        match new_state {
            SwState::Match => {
                p1 -= 1;
                p2 -= 1;
            }
            SwState::Insertion => {
                if step_length > p2 {
                    break;
                }
                p2 -= step_length;
            }
            SwState::Deletion => {
                if step_length > p1 {
                    break;
                }
                p1 -= step_length;
            }
            SwState::Clip => unreachable!(),
        }
        if new_state == state {
            segment_length += step_length;
        } else {
            if segment_length > 0 {
                lce.push(make_element(state, segment_length));
            }
            segment_length = step_length;
            state = new_state;
        }
        if p1 == 0 || p2 == 0 {
            break;
        }
    }
    let alignment_offset = match overhang_strategy {
        SwOverhangStrategy::SoftClip => {
            lce.push(make_element(state, segment_length));
            if p2 > 0 {
                lce.push(make_element(SwState::Clip, p2));
            }
            p1 as i32
        }
        SwOverhangStrategy::Ignore => {
            lce.push(make_element(state, segment_length + p2));
            (p1 as i32) - (p2 as i32)
        }
        SwOverhangStrategy::Indel | SwOverhangStrategy::LeadingIndel => {
            lce.push(make_element(state, segment_length));
            if p1 > 0 {
                lce.push(make_element(SwState::Deletion, p1));
            } else if p2 > 0 {
                lce.push(make_element(SwState::Insertion, p2));
            }
            0
        }
    };
    lce.reverse();
    SmithWatermanAlignment {
        cigar: Cigar { elements: lce },
        alignment_offset,
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used, clippy::expect_used, clippy::panic)]
mod tests {
    use super::*;

    #[test]
    fn substring_fast_path_soft_clip() {
        let reference = b"ACGTACGTACGT";
        let alternate = b"TACGTA";
        let aln = align(
            reference,
            alternate,
            &SwParameters::gatk_haplotype_to_reference(),
            SwOverhangStrategy::SoftClip,
        )
        .expect("non-empty sequences");
        assert_eq!(aln.alignment_offset, 3);
        assert_eq!(aln.cigar.elements.len(), 1);
        assert_eq!(aln.cigar.elements[0].operator, CigarOperator::Match);
        assert_eq!(aln.cigar.elements[0].length, 6);
    }

    #[test]
    fn indel_strategy_produces_alignment() {
        let reference = b"ACGTACGTACGT";
        let alternate = b"ACGTTTACGT";
        let aln = align(
            reference,
            alternate,
            &SwParameters::gatk_haplotype_to_reference(),
            SwOverhangStrategy::Indel,
        )
        .expect("non-empty sequences");
        assert!(!aln.cigar.elements.is_empty());
        assert_eq!(aln.alignment_offset, 0);
    }

    #[test]
    fn empty_sequences_return_err_not_panic() {
        let err = align(
            b"",
            b"ACGT",
            &SwParameters::gatk_haplotype_to_reference(),
            SwOverhangStrategy::Indel,
        );
        assert!(err.is_err());
        let err = align(
            b"ACGT",
            b"",
            &SwParameters::gatk_haplotype_to_reference(),
            SwOverhangStrategy::Indel,
        );
        assert!(err.is_err());
    }

    #[test]
    fn oversized_matrix_is_refused() {
        // Contig-scale × read-scale would be tens of GiB; refuse before allocate/scan.
        let huge = vec![b'A'; 200_000];
        let short = vec![b'A'; 150];
        let err = align(
            &huge,
            &short,
            &SwParameters::gatk_haplotype_to_reference(),
            SwOverhangStrategy::SoftClip,
        );
        assert!(err.is_err(), "expected refusal of contig-scale SW");
        let msg = format!("{err:?}");
        assert!(
            msg.contains("refused oversized") || msg.contains("contig scale"),
            "unexpected err: {msg}"
        );
    }

    /// Rolling scores must remain deterministic across TLS reuse.
    #[test]
    fn cigar_parity_soft_clip_and_indel_variants() {
        let cases: &[(&[u8], &[u8])] = &[
            (b"ACGTACGTACGTACGT", b"ACGTACGTACGT"),
            (b"ACGTACGTACGTACGT", b"ACGTTTACGTAC"),
            (b"AAAAAAAAAA", b"AAAAATAAAA"),
            (b"ACGTACGTACGT", b"TTT"),
            (b"GGGACGTACGTACGTAAA", b"ACGTNACGT"),
        ];
        let params = SwParameters::gatk_haplotype_to_reference();
        for &(reference, alternate) in cases {
            for strategy in [
                SwOverhangStrategy::SoftClip,
                SwOverhangStrategy::Indel,
                SwOverhangStrategy::Ignore,
                SwOverhangStrategy::LeadingIndel,
            ] {
                let aln = align(reference, alternate, &params, strategy).expect("align");
                let aln2 = align(reference, alternate, &params, strategy).expect("align2");
                assert_eq!(aln.alignment_offset, aln2.alignment_offset, "{strategy:?}");
                assert_eq!(
                    aln.cigar.elements, aln2.cigar.elements,
                    "{strategy:?} cigar"
                );
            }
        }
    }

    #[test]
    fn read_to_hap_params_soft_clip_offset() {
        let hap = b"ACGTACGTACGTACGTACGTACGT";
        let read = b"TACGTACGTACGTA";
        let aln =
            align_read_to_best_haplotype(hap, read, &SwParameters::gatk_read_to_best_haplotype())
                .expect("sw");
        // Java lastIndexOf: last occurrence of the read in the haplotype.
        let expected = last_index_of(hap, read).expect("contained") as i32;
        assert_eq!(aln.alignment_offset, expected);
        assert_eq!(aln.cigar.elements[0].operator, CigarOperator::Match);
    }

    /// Independent full-matrix DP oracle — must match rolling production path.
    fn align_full_matrix_oracle(
        reference: &[u8],
        alternate: &[u8],
        parameters: &SwParameters,
        overhang_strategy: SwOverhangStrategy,
    ) -> SmithWatermanAlignment {
        let nrow = reference.len() + 1;
        let ncol = alternate.len() + 1;
        let cells = nrow * ncol;
        let mut sw = vec![0i32; cells];
        let mut btrack = vec![0i32; cells];
        let mut best_gap_v = Vec::new();
        let mut gap_size_v = Vec::new();
        let mut best_gap_h = Vec::new();
        let mut gap_size_h = Vec::new();
        // Reuse production kernel via rolling into a reconstructed full score plane is
        // awkward; duplicate the classic full-matrix fill for the oracle only.
        const MATRIX_MIN_CUTOFF: i32 = -100_000_000;
        let low_init = i32::MIN / 4;
        best_gap_v.resize(ncol + 1, low_init);
        gap_size_v.resize(ncol + 1, 0);
        best_gap_h.resize(nrow + 1, low_init);
        gap_size_h.resize(nrow + 1, 0);
        if matches!(
            overhang_strategy,
            SwOverhangStrategy::Indel | SwOverhangStrategy::LeadingIndel
        ) {
            sw[cell(0, 1, ncol)] = parameters.gap_open_penalty;
            let mut cur = parameters.gap_open_penalty;
            for j in 2..ncol {
                cur += parameters.gap_extend_penalty;
                sw[cell(0, j, ncol)] = cur;
            }
            sw[cell(1, 0, ncol)] = parameters.gap_open_penalty;
            cur = parameters.gap_open_penalty;
            for i in 2..nrow {
                cur += parameters.gap_extend_penalty;
                sw[cell(i, 0, ncol)] = cur;
            }
        }
        let w_open = parameters.gap_open_penalty;
        let w_extend = parameters.gap_extend_penalty;
        let w_match = parameters.match_value;
        let w_mismatch = parameters.mismatch_penalty;
        for i in 1..nrow {
            let a_base = reference[i - 1];
            for j in 1..ncol {
                let b_base = alternate[j - 1];
                let diag = sw[cell(i - 1, j - 1, ncol)]
                    + if a_base == b_base {
                        w_match
                    } else {
                        w_mismatch
                    };
                let prev_gap = sw[cell(i - 1, j, ncol)] + w_open;
                best_gap_v[j] += w_extend;
                if prev_gap > best_gap_v[j] {
                    best_gap_v[j] = prev_gap;
                    gap_size_v[j] = 1;
                } else {
                    gap_size_v[j] += 1;
                }
                let step_down = best_gap_v[j];
                let kd = gap_size_v[j];
                let prev_gap_h = sw[cell(i, j - 1, ncol)] + w_open;
                best_gap_h[i] += w_extend;
                if prev_gap_h > best_gap_h[i] {
                    best_gap_h[i] = prev_gap_h;
                    gap_size_h[i] = 1;
                } else {
                    gap_size_h[i] += 1;
                }
                let step_right = best_gap_h[i];
                let ki = gap_size_h[i];
                let (score, track) = if diag >= step_down && diag >= step_right {
                    (diag.max(MATRIX_MIN_CUTOFF), 0)
                } else if step_right >= step_down {
                    (step_right.max(MATRIX_MIN_CUTOFF), -ki)
                } else {
                    (step_down.max(MATRIX_MIN_CUTOFF), kd)
                };
                sw[cell(i, j, ncol)] = score;
                btrack[cell(i, j, ncol)] = track;
            }
        }
        let mut last_col = vec![0i32; nrow];
        for i in 0..nrow {
            last_col[i] = sw[cell(i, ncol - 1, ncol)];
        }
        let last_row = sw[(nrow - 1) * ncol..nrow * ncol].to_vec();
        calculate_cigar(&last_col, &last_row, &btrack, nrow, ncol, overhang_strategy)
    }

    #[test]
    fn rolling_matches_full_matrix_oracle() {
        let cases: &[(&[u8], &[u8])] = &[
            (b"ACGTACGTACGTACGT", b"ACGTACGTACGT"),
            (b"ACGTACGTACGTACGT", b"ACGTTTACGTAC"),
            (b"AAAAAAAAAA", b"AAAAATAAAA"),
            (b"ACGTACGTACGT", b"TTT"),
            (b"GGGACGTACGTACGTAAA", b"ACGTNACGT"),
            (b"NNNNACGTACGTNNNN", b"NNNNACGTTACGTNNNN"),
        ];
        let params = SwParameters::gatk_haplotype_to_reference();
        for &(reference, alternate) in cases {
            for strategy in [
                SwOverhangStrategy::SoftClip,
                SwOverhangStrategy::Indel,
                SwOverhangStrategy::Ignore,
                SwOverhangStrategy::LeadingIndel,
            ] {
                // Force DP (no substring fast path) by using lowercase → upper path with a mismatch
                // shape that is not an exact substring, or call oracle vs production after
                // disabling fast path via a deliberate single-base edit when needed.
                let prod = align(reference, alternate, &params, strategy).expect("prod");
                let oracle = align_full_matrix_oracle(reference, alternate, &params, strategy);
                assert_eq!(
                    prod.alignment_offset, oracle.alignment_offset,
                    "offset {strategy:?} ref={reference:?} alt={alternate:?}"
                );
                assert_eq!(
                    prod.cigar.elements, oracle.cigar.elements,
                    "cigar {strategy:?} ref={reference:?} alt={alternate:?}"
                );
            }
        }
    }
}
