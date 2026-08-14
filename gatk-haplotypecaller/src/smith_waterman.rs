//! GATK `SmithWatermanJavaAligner` parity (JAVA implementation path).
//! Score / backtrack grids are flat row-major `Vec<i32>` (one allocation each) instead
//! `Vec<Vec<i32>>` (one allocation per row).

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
    for r in (0..=reference.len() - qlen).rev() {
        if reference[r] != first || reference[r + qlen - 1] != last {
            continue;
        }
        if reference[r..r + qlen] == *query {
            return Some(r);
        }
    }
    None
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
        let mut scratch = cell.borrow_mut();
        scratch.ensure(cells, nrow, ncol);
        // Take buffers out so we can mutably borrow planes without aliasing `scratch`.
        let mut sw = std::mem::take(&mut scratch.sw);
        let mut btrack = std::mem::take(&mut scratch.btrack);
        let mut best_gap_v = std::mem::take(&mut scratch.best_gap_v);
        let mut gap_size_v = std::mem::take(&mut scratch.gap_size_v);
        let mut best_gap_h = std::mem::take(&mut scratch.best_gap_h);
        let mut gap_size_h = std::mem::take(&mut scratch.gap_size_h);
        if sw.len() < cells {
            sw.resize(cells, 0);
            btrack.resize(cells, 0);
        } else {
            sw[..cells].fill(0);
            btrack[..cells].fill(0);
        }
        calculate_matrix(
            reference,
            alternate,
            nrow,
            ncol,
            &mut sw[..cells],
            &mut btrack[..cells],
            &mut best_gap_v,
            &mut gap_size_v,
            &mut best_gap_h,
            &mut gap_size_h,
            overhang_strategy,
            parameters,
        );
        let aln = calculate_cigar(
            &sw[..cells],
            &btrack[..cells],
            nrow,
            ncol,
            overhang_strategy,
        );
        // Soft-keep normal-sized planes; drop oversized allocations so Peak does not retain
        // multi-megabase matrices between calls. Region/engine end still hard-clears via
        // [`release_sw_tls_scratch`].
        if cells > SW_TLS_SOFT_KEEP_CELLS
            || sw.capacity() > SW_TLS_SOFT_KEEP_CELLS
            || btrack.capacity() > SW_TLS_SOFT_KEEP_CELLS
        {
            scratch.clear();
        } else {
            scratch.sw = sw;
            scratch.btrack = btrack;
            scratch.best_gap_v = best_gap_v;
            scratch.gap_size_v = gap_size_v;
            scratch.best_gap_h = best_gap_h;
            scratch.gap_size_h = gap_size_h;
        }
        Ok(aln)
    })
}

/// Soft-keep ceiling for SW TLS matrix cells (sw + btrack). Larger requests still compute,
/// but planes are not retained afterward.
const SW_TLS_SOFT_KEEP_CELLS: usize = 1 << 20; // 1_048_576 cells

struct SwScratch {
    sw: Vec<i32>,
    btrack: Vec<i32>,
    best_gap_v: Vec<i32>,
    gap_size_v: Vec<i32>,
    best_gap_h: Vec<i32>,
    gap_size_h: Vec<i32>,
}

impl SwScratch {
    fn new() -> Self {
        Self {
            sw: Vec::new(),
            btrack: Vec::new(),
            best_gap_v: Vec::new(),
            gap_size_v: Vec::new(),
            best_gap_h: Vec::new(),
            gap_size_h: Vec::new(),
        }
    }

    fn ensure(&mut self, cells: usize, nrow: usize, ncol: usize) {
        // Do not grow retained capacity past the soft-keep ceiling without clearing first.
        if cells > SW_TLS_SOFT_KEEP_CELLS
            || self.sw.capacity() > SW_TLS_SOFT_KEEP_CELLS
            || self.btrack.capacity() > SW_TLS_SOFT_KEEP_CELLS
        {
            if !self.sw.is_empty() || self.sw.capacity() > 0 {
                self.clear();
            }
        }
        if self.sw.len() < cells {
            self.sw.resize(cells, 0);
            self.btrack.resize(cells, 0);
        } else {
            self.sw[..cells].fill(0);
            self.btrack[..cells].fill(0);
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
        *self = Self::new();
    }
}

thread_local! {
    static SW_SCRATCH: RefCell<SwScratch> = RefCell::new(SwScratch::new());
}

/// Drop Smith-Waterman TLS arenas (full drop).
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
    sw: &mut [i32],
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
        let prev_base = (i - 1) * ncol;
        let cur_base = i * ncol;
        for j in 1..ncol {
            let b_base = alternate[j - 1];
            let step_diag = sw[prev_base + j - 1]
                + if a_base == b_base {
                    w_match
                } else {
                    w_mismatch
                };
            let prev_gap = sw[prev_base + j] + w_open;
            best_gap_v[j] += w_extend;
            if prev_gap > best_gap_v[j] {
                best_gap_v[j] = prev_gap;
                gap_size_v[j] = 1;
            } else {
                gap_size_v[j] += 1;
            }
            let step_down = best_gap_v[j];
            let kd = gap_size_v[j];
            let prev_gap_h = sw[cur_base + j - 1] + w_open;
            best_gap_h[i] += w_extend;
            if prev_gap_h > best_gap_h[i] {
                best_gap_h[i] = prev_gap_h;
                gap_size_h[i] = 1;
            } else {
                gap_size_h[i] += 1;
            }
            let step_right = best_gap_h[i];
            let ki = gap_size_h[i];
            let cur_idx = cur_base + j;
            if step_diag >= step_down && step_diag >= step_right {
                sw[cur_idx] = step_diag.max(MATRIX_MIN_CUTOFF);
                btrack[cur_idx] = 0;
            } else if step_right >= step_down {
                sw[cur_idx] = step_right.max(MATRIX_MIN_CUTOFF);
                btrack[cur_idx] = -ki;
            } else {
                sw[cur_idx] = step_down.max(MATRIX_MIN_CUTOFF);
                btrack[cur_idx] = kd;
            }
        }
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
    sw: &[i32],
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
            let cur = sw[cell(i, alt_length, ncol)];
            if cur >= maxscore {
                p1 = i;
                maxscore = cur;
            }
        }
        if overhang_strategy != SwOverhangStrategy::LeadingIndel {
            let bottom = ref_length * ncol;
            for j in 1..ncol {
                let cur = sw[bottom + j];
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
        let btr = btrack[cell(p1, p2, ncol)];
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
}
