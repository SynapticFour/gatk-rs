//! GATK `CigarUtils.calculateCigar` + `AlignmentUtils` trim/left-align for haplotype realignment.

use crate::cigar::{length_on_read, length_on_reference, Cigar, CigarElement, CigarOperator};
use crate::cigar_builder::CigarBuilder;
use crate::smith_waterman::{align, SmithWatermanAlignment, SwOverhangStrategy, SwParameters};

pub const SW_PAD: &[u8] = b"NNNNNNNNNN";

/// Optional haplotype CIGAR after Smith–Waterman alignment to reference.
/// # Invariants
/// `None` indicates SW failure (soft-clip in alignment or positive offset per GATK rules).
/// # Ownership
/// Owns optional [`Cigar`]; cheap clone for haplotype assembly paths.
/// # Mutation
/// Immutable result of haplotype CIGAR calculation helpers.
/// # Biological assumptions
/// CIGAR describes gapped alignment of alt haplotype bytes to padded reference span.
/// # Java equivalence
/// GATK `CigarUtils.calculateCigar` result wrapper for haplotype realignment.
#[derive(Debug, Clone)]
pub struct HaplotypeCigarResult {
    pub cigar: Option<Cigar>,
}

fn is_sw_failure(alignment: &SmithWatermanAlignment) -> bool {
    if alignment.alignment_offset > 0 {
        return true;
    }
    alignment
        .cigar
        .elements
        .iter()
        .any(|e| e.operator == CigarOperator::SoftClip)
}

/// GATK `CigarUtils.calculateCigar` for haplotype-to-reference (`NEW_SW_PARAMETERS`).
pub fn calculate_haplotype_cigar(
    ref_seq: &[u8],
    alt_seq: &[u8],
    parameters: &SwParameters,
) -> Option<Cigar> {
    calculate_haplotype_cigar_with_strategy(
        ref_seq,
        alt_seq,
        parameters,
        SwOverhangStrategy::SoftClip,
    )
}

/// Same as [`calculate_haplotype_cigar`] with an explicit SW overhang strategy.
pub fn calculate_haplotype_cigar_with_strategy(
    ref_seq: &[u8],
    alt_seq: &[u8],
    parameters: &SwParameters,
    strategy: SwOverhangStrategy,
) -> Option<Cigar> {
    calculate_haplotype_cigar_sw(ref_seq, alt_seq, parameters, strategy).map(|r| r.cigar)
}

/// Leading ref bases trimmed from SW padding + left-align (GATK `trimmedLeadingDeletions` + left-align).
pub fn alignment_start_from_haplotype_sw(
    alignment: &SmithWatermanAlignment,
    trimmed_leading_deletions: usize,
    left_align_leading_deletions: usize,
) -> usize {
    let leading = trimmed_leading_deletions.saturating_add(left_align_leading_deletions);
    if leading > 0 {
        leading
    } else {
        alignment.alignment_offset.max(0) as usize
    }
}

/// GATK `ReadThreadingAssembler` cigar path: SOFTCLIP first, INDEL retry on ref-span mismatch.
pub fn calculate_haplotype_cigar_for_assembly(
    ref_seq: &[u8],
    alt_seq: &[u8],
    ref_cigar_length: usize,
    parameters: &SwParameters,
) -> Option<Cigar> {
    calculate_haplotype_cigar_for_assembly_with_offset(
        ref_seq,
        alt_seq,
        ref_cigar_length,
        parameters,
    )
    .map(|r| r.cigar)
}

/// Like [`calculate_haplotype_cigar_for_assembly`] but also returns GATK `alignmentStartHapwrtRef`.
pub fn calculate_haplotype_cigar_for_assembly_with_offset(
    ref_seq: &[u8],
    alt_seq: &[u8],
    ref_cigar_length: usize,
    parameters: &SwParameters,
) -> Option<HaplotypeAssemblyCigar> {
    let soft =
        calculate_haplotype_cigar_sw(ref_seq, alt_seq, parameters, SwOverhangStrategy::SoftClip)?;
    if soft.cigar.reference_length() == ref_cigar_length {
        return Some(soft);
    }
    let indel =
        calculate_haplotype_cigar_sw(ref_seq, alt_seq, parameters, SwOverhangStrategy::Indel)?;
    if indel.cigar.elements.iter().any(|e| e.operator.is_indel()) {
        return Some(indel);
    }
    None
}

/// CIGAR + hap offset in padded reference coordinates (GATK `Haplotype.setAlignmentStartHapwrtRef`).
/// # Invariants
/// `alignment_start_hap_wrt_ref` is the haplotype start offset in padded reference coordinates.
/// `cigar` reference length should match assembly expectations when SoftClip/Indel strategy succeeds.
/// # Ownership
/// Owns [`Cigar`] and scalar offset.
/// # Mutation
/// Immutable SW/assembly CIGAR product.
/// # Biological assumptions
/// Aligns assembled alt haplotype to padded reference for EventMap/CIGAR genotyping.
/// # Java equivalence
/// GATK haplotype CIGAR + `alignmentStartHapwrtRef` from `CigarUtils` / assembler path.
#[derive(Debug, Clone)]
pub struct HaplotypeAssemblyCigar {
    pub cigar: Cigar,
    pub alignment_start_hap_wrt_ref: usize,
}

fn calculate_haplotype_cigar_sw(
    ref_seq: &[u8],
    alt_seq: &[u8],
    parameters: &SwParameters,
    strategy: SwOverhangStrategy,
) -> Option<HaplotypeAssemblyCigar> {
    if alt_seq.is_empty() {
        let mut c = Cigar::new();
        c.push(ref_seq.len(), CigarOperator::Deletion);
        return Some(HaplotypeAssemblyCigar {
            cigar: c,
            alignment_start_hap_wrt_ref: 0,
        });
    }
    if ref_seq == alt_seq {
        let mut c = Cigar::new();
        c.push(ref_seq.len(), CigarOperator::Match);
        return Some(HaplotypeAssemblyCigar {
            cigar: c,
            alignment_start_hap_wrt_ref: 0,
        });
    }
    // Equal-length SNP/MNP: EventMap reads mismatches from Match ops — skip padded SW.
    // Length-changing alts still need SoftClip/Indel SW below.
    if ref_seq.len() == alt_seq.len() {
        let mut c = Cigar::new();
        c.push(ref_seq.len(), CigarOperator::Match);
        return Some(HaplotypeAssemblyCigar {
            cigar: c,
            alignment_start_hap_wrt_ref: 0,
        });
    }
    let mut padded_ref = Vec::with_capacity(SW_PAD.len() + ref_seq.len() + SW_PAD.len());
    padded_ref.extend_from_slice(SW_PAD);
    padded_ref.extend_from_slice(ref_seq);
    padded_ref.extend_from_slice(SW_PAD);
    let mut padded_alt = Vec::with_capacity(SW_PAD.len() + alt_seq.len() + SW_PAD.len());
    padded_alt.extend_from_slice(SW_PAD);
    padded_alt.extend_from_slice(alt_seq);
    padded_alt.extend_from_slice(SW_PAD);

    let Ok(alignment) = align(&padded_ref, &padded_alt, parameters, strategy) else {
        return None;
    };
    if is_sw_failure(&alignment) {
        return None;
    }
    let base_start = SW_PAD.len();
    let base_end = padded_alt.len() - SW_PAD.len() - 1;
    let trimmed = trim_cigar_by_bases(&alignment.cigar, base_start, base_end);
    let mut non_standard = trimmed.cigar;
    let trimmed_leading = trimmed.leading_deletions_removed;
    let trimmed_trailing = trimmed.trailing_deletions_removed;
    if trimmed_trailing > 0 {
        non_standard.push(trimmed_trailing, CigarOperator::Deletion);
    }
    let left = left_align_indels(&non_standard, ref_seq, alt_seq, trimmed_leading);
    let total_leading = trimmed_leading + left.leading_deletions_removed;
    let total_trailing = left.trailing_deletions_removed;
    let cigar = if total_leading == 0 && total_trailing == 0 {
        left.cigar
    } else {
        let mut result = Cigar::new();
        if total_leading > 0 {
            result.push(total_leading, CigarOperator::Deletion);
        }
        for e in left.cigar.elements {
            result.push(e.length, e.operator);
        }
        if total_trailing > 0 {
            result.push(total_trailing, CigarOperator::Deletion);
        }
        result
    };
    Some(HaplotypeAssemblyCigar {
        alignment_start_hap_wrt_ref: alignment_start_from_haplotype_sw(
            &alignment,
            trimmed_leading,
            left.leading_deletions_removed,
        ),
        cigar,
    })
}

pub fn trim_cigar_by_reference(
    cigar: &Cigar,
    start: usize,
    end: usize,
) -> crate::cigar_builder::CigarBuilderResult {
    trim_cigar(cigar, start, end, true)
}

fn trim_cigar_by_bases(
    cigar: &Cigar,
    start: usize,
    end: usize,
) -> crate::cigar_builder::CigarBuilderResult {
    trim_cigar(cigar, start, end, false)
}

/// GATK `Haplotype.getConsolidatedPaddedCigar` — append `pad_size` match ops on hap read axis.
pub fn consolidated_padded_cigar(cigar: &Cigar, pad_size: usize) -> Cigar {
    let mut out = cigar.clone();
    if pad_size > 0 {
        out.push(pad_size, CigarOperator::Match);
    }
    out
}

/// GATK `AlignmentUtils.readStartOnReferenceHaplotype`.
pub fn read_start_on_reference_haplotype(
    haplotype_vs_ref_cigar: &Cigar,
    read_start_on_haplotype: usize,
) -> usize {
    if read_start_on_haplotype == 0 {
        return 0;
    }
    let mut ref_bases = 0usize;
    let mut hap_bases = 0usize;
    for e in &haplotype_vs_ref_cigar.elements {
        ref_bases += length_on_reference(e.operator, e.length);
        hap_bases += length_on_read(e.operator, e.length);
        if hap_bases >= read_start_on_haplotype {
            let excess = if e.operator.consumes_reference_bases() {
                hap_bases - read_start_on_haplotype
            } else {
                0
            };
            return ref_bases - excess;
        }
    }
    ref_bases
}

struct CigarPairTransform {
    op13: Option<CigarOperator>,
    advance12: usize,
    advance23: usize,
}

fn cigar_pair_transform(op12: CigarOperator, op23: CigarOperator) -> CigarPairTransform {
    use CigarOperator::{Deletion, Insertion, Match};
    match (op12, op23) {
        (Match, Match) => CigarPairTransform {
            op13: Some(Match),
            advance12: 1,
            advance23: 1,
        },
        (Match, Insertion) => CigarPairTransform {
            op13: Some(Insertion),
            advance12: 1,
            advance23: 1,
        },
        (Match, Deletion) => CigarPairTransform {
            op13: Some(Deletion),
            advance12: 0,
            advance23: 1,
        },
        (Deletion, Match) => CigarPairTransform {
            op13: Some(Deletion),
            advance12: 1,
            advance23: 1,
        },
        (Deletion, Deletion) => CigarPairTransform {
            op13: Some(Deletion),
            advance12: 0,
            advance23: 1,
        },
        (Deletion, Insertion) => CigarPairTransform {
            op13: None,
            advance12: 1,
            advance23: 1,
        },
        (Insertion, Match) => CigarPairTransform {
            op13: Some(Insertion),
            advance12: 1,
            advance23: 0,
        },
        (Insertion, Deletion) => CigarPairTransform {
            op13: Some(Insertion),
            advance12: 1,
            advance23: 0,
        },
        (Insertion, Insertion) => CigarPairTransform {
            op13: Some(Insertion),
            advance12: 1,
            advance23: 0,
        },
        _ => CigarPairTransform {
            op13: None,
            advance12: 1,
            advance23: 1,
        },
    }
}

/// GATK `AlignmentUtils.applyCigarToCigar` (read→hap composed with hap→ref).
pub fn apply_cigar_to_cigar(first_to_second: &Cigar, second_to_third: &Cigar) -> Cigar {
    let mut out = Cigar::new();
    let n12 = first_to_second.elements.len();
    let n23 = second_to_third.elements.len();
    let mut i12 = 0usize;
    let mut i23 = 0usize;
    let mut elt12 = 0usize;
    let mut elt23 = 0usize;
    while i12 < n12 && i23 < n23 {
        let e12 = &first_to_second.elements[i12];
        let e23 = &second_to_third.elements[i23];
        let transform = cigar_pair_transform(e12.operator, e23.operator);
        if let Some(op13) = transform.op13 {
            out.push(1, op13);
        }
        elt12 += transform.advance12;
        elt23 += transform.advance23;
        if elt12 == e12.length {
            i12 += 1;
            elt12 = 0;
        }
        if elt23 == e23.length {
            i23 += 1;
            elt23 = 0;
        }
    }
    out
}

/// GATK `AlignmentUtils.leftAlignIndels` for read→ref CIGAR after hap projection.
pub fn left_align_indels_for_read(
    cigar: &Cigar,
    ref_bases: &[u8],
    read_bases: &[u8],
    read_start_on_reference_haplotype: usize,
) -> crate::cigar_builder::CigarBuilderResult {
    left_align_indels(
        cigar,
        ref_bases,
        read_bases,
        read_start_on_reference_haplotype,
    )
}

/// Trim cigar on hap/read axis between `start` and `end` inclusive (GATK `trimCigarByBases`).
pub fn trim_cigar_by_bases_public(
    cigar: &Cigar,
    start: usize,
    end: usize,
) -> crate::cigar_builder::CigarBuilderResult {
    trim_cigar_by_bases(cigar, start, end)
}

/// GATK `AlignmentUtils.getBasesCoveringRefInterval` (returns `None` if interval starts/ends in deletion).
pub fn get_bases_covering_ref_interval(
    ref_start: usize,
    ref_end: usize,
    bases: &[u8],
    bases_start_on_ref: usize,
    cigar: &Cigar,
) -> Option<Vec<u8>> {
    if ref_start > ref_end {
        return None;
    }
    let read_len: usize = cigar
        .elements
        .iter()
        .map(|e| length_on_read(e.operator, e.length))
        .sum();
    if bases.len() != read_len {
        return None;
    }
    let mut ref_pos = bases_start_on_ref;
    let mut bases_pos = 0usize;
    let mut bases_start_idx = None;
    let mut bases_stop_idx = None;
    let mut done = false;
    for el in &cigar.elements {
        if done {
            break;
        }
        match el.operator {
            CigarOperator::Insertion => {
                bases_pos += el.length;
            }
            CigarOperator::SoftClip => {
                bases_pos += el.length;
            }
            CigarOperator::Match => {
                for _ in 0..el.length {
                    if ref_pos == ref_start {
                        bases_start_idx = Some(bases_pos);
                    }
                    if ref_pos == ref_end {
                        bases_stop_idx = Some(bases_pos);
                        done = true;
                        break;
                    }
                    ref_pos += 1;
                    bases_pos += 1;
                }
            }
            CigarOperator::Deletion => {
                for _ in 0..el.length {
                    if ref_pos == ref_end || ref_pos == ref_start {
                        return None;
                    }
                    ref_pos += 1;
                }
            }
            _ => return None,
        }
    }
    let start_i = bases_start_idx?;
    let stop_i = bases_stop_idx?;
    Some(bases[start_i..=stop_i].to_vec())
}

fn trim_cigar(
    cigar: &Cigar,
    start: usize,
    end: usize,
    by_reference: bool,
) -> crate::cigar_builder::CigarBuilderResult {
    let mut new_elements = CigarBuilder::default_trim();
    let mut element_end = 0usize;
    for elt in &cigar.elements {
        let element_start = element_end;
        element_end = element_start
            + if by_reference {
                length_on_reference(elt.operator, elt.length)
            } else {
                length_on_read(elt.operator, elt.length)
            };
        if element_end < start || (element_end == start && element_start < start) {
            continue;
        }
        if element_start > end && element_end > end + 1 {
            break;
        }
        let overlap = if element_end == element_start {
            elt.length
        } else {
            (end + 1).min(element_end) - start.max(element_start)
        };
        new_elements.add(CigarElement {
            length: overlap,
            operator: elt.operator,
        });
    }
    new_elements.make_and_record()
}

#[derive(Clone, Copy)]
struct IndexRange {
    start: usize,
    end: usize,
}

impl IndexRange {
    fn new(start: usize, end: usize) -> Self {
        Self { start, end }
    }
    fn size(self) -> usize {
        self.end.saturating_sub(self.start)
    }
    fn shift_left(&mut self, n: usize) {
        // Java IndexRange validates from>=0; usize wrap previously produced start≈usize::MAX
        // and panics in next_base_on_left (`start - 1`).
        self.start = self.start.saturating_sub(n);
        self.end = self.end.saturating_sub(n);
        if self.start > self.end {
            self.end = self.start;
        }
    }
    fn shift_start(&mut self, n: usize) {
        self.start = self.start.saturating_add(n);
        if self.start > self.end {
            self.start = self.end;
        }
    }
    fn shift_end_left(&mut self, n: usize) {
        self.end = self.end.saturating_sub(n);
        if self.end < self.start {
            self.end = self.start;
        }
    }
    fn shift_start_left(&mut self, n: usize) {
        self.start = self.start.saturating_sub(n);
    }
}

/// GATK `readStartOnReferenceHaplotype` (legacy name).
pub fn reference_offset_at_hap_read_offset(cigar: &Cigar, hap_read_offset: usize) -> usize {
    read_start_on_reference_haplotype(cigar, hap_read_offset)
}

fn left_align_indels(
    cigar: &Cigar,
    ref_bases: &[u8],
    read_bases: &[u8],
    read_start: usize,
) -> crate::cigar_builder::CigarBuilderResult {
    if !cigar.elements.iter().any(|e| e.operator.is_indel()) {
        return crate::cigar_builder::CigarBuilderResult {
            cigar: cigar.clone(),
            leading_deletions_removed: 0,
            trailing_deletions_removed: 0,
        };
    }
    // INVARIANT: caller checked `any(|e| e.operator.is_indel)` above.
    #[allow(clippy::unwrap_used)]
    let last_indel = cigar
        .elements
        .iter()
        .rposition(|e| e.operator.is_indel())
        .unwrap();
    let necessary_ref: usize = read_start
        + cigar
            .elements
            .iter()
            .take(last_indel + 1)
            .map(|e| length_on_reference(e.operator, e.length))
            .sum::<usize>();
    assert!(
        necessary_ref <= ref_bases.len(),
        "read goes past end of reference"
    );
    let ref_length = cigar.reference_length();
    let mut result_rtl: Vec<CigarElement> = Vec::new();
    let mut ref_indel = IndexRange::new(read_start + ref_length, read_start + ref_length);
    let mut read_indel = IndexRange::new(read_bases.len(), read_bases.len());
    for n in (0..cigar.elements.len()).rev() {
        let element = cigar.elements[n];
        if element.operator.is_indel() {
            ref_indel.shift_start_left(length_on_reference(element.operator, element.length));
            read_indel.shift_start_left(length_on_read(element.operator, element.length));
        } else if ref_indel.size() == 0 && read_indel.size() == 0 {
            result_rtl.push(element);
            ref_indel.shift_left(length_on_reference(element.operator, element.length));
            read_indel.shift_left(length_on_read(element.operator, element.length));
        } else {
            let max_shift = if element.operator.is_alignment() {
                element.length
            } else {
                0
            };
            let (start_shift, end_shift) = normalize_alleles(
                &[ref_bases, read_bases],
                &mut [ref_indel, read_indel],
                max_shift,
                true,
            );
            result_rtl.push(CigarElement {
                length: end_shift,
                operator: CigarOperator::Match,
            });
            let emit_indel =
                n == 0 || start_shift < max_shift as i32 || !element.operator.is_alignment();
            let new_match_left = if start_shift < 0 {
                (-start_shift) as usize
            } else {
                0
            };
            let remaining = if start_shift < 0 {
                element.length
            } else {
                element.length - start_shift as usize
            };
            if emit_indel {
                result_rtl.push(CigarElement {
                    length: ref_indel.size(),
                    operator: CigarOperator::Deletion,
                });
                result_rtl.push(CigarElement {
                    length: read_indel.size(),
                    operator: CigarOperator::Insertion,
                });
                ref_indel.shift_end_left(ref_indel.size());
                read_indel.shift_end_left(read_indel.size());
                ref_indel.shift_left(
                    new_match_left
                        + if element.operator.consumes_reference_bases() {
                            remaining
                        } else {
                            0
                        },
                );
                read_indel.shift_left(
                    new_match_left
                        + if element.operator.consumes_read_bases() {
                            remaining
                        } else {
                            0
                        },
                );
            }
            result_rtl.push(CigarElement {
                length: new_match_left,
                operator: CigarOperator::Match,
            });
            result_rtl.push(CigarElement {
                length: remaining,
                operator: element.operator,
            });
        }
    }
    result_rtl.push(CigarElement {
        length: ref_indel.size(),
        operator: CigarOperator::Deletion,
    });
    result_rtl.push(CigarElement {
        length: read_indel.size(),
        operator: CigarOperator::Insertion,
    });
    result_rtl.reverse();
    let mut builder = CigarBuilder::default_trim();
    for e in result_rtl {
        builder.add(e);
    }
    builder.make_and_record()
}

fn normalize_alleles(
    sequences: &[&[u8]],
    bounds: &mut [IndexRange],
    max_shift: usize,
    trim: bool,
) -> (i32, usize) {
    debug_assert_eq!(sequences.len(), bounds.len());
    // Java `AlignmentUtils.normalizeAlleles`:
    //   Utils.validateArg(maxShift <= bound.getStart(), ...)
    // Cap instead of aborting: callers that overshoot still left-align only as far as
    // start allows, preserving the invariant that `nextBaseOnLeft` never indexes -1.
    let max_shift = max_shift.min(bounds.iter().map(|b| b.start).min().unwrap_or(0));
    let mut start_shift: i32 = 0;
    let mut end_shift = 0usize;
    let mut min_size = bounds.iter().map(|b| b.size()).min().unwrap_or(0);
    while trim && min_size > 0 && last_base_on_right_is_same(sequences, bounds) {
        for b in bounds.iter_mut() {
            b.shift_end_left(1);
        }
        min_size -= 1;
        end_shift += 1;
    }
    while trim && min_size > 0 && first_base_on_left_is_same(sequences, bounds) {
        for b in bounds.iter_mut() {
            b.shift_start(1);
        }
        min_size -= 1;
        start_shift -= 1;
    }
    // Signed compare matches Java `startShift < maxShift` (negative startShift after
    // left-trim must still enter the left-align loop). Casting negative i32 → usize
    // made the old condition always false after trim and could skip left-align.
    while start_shift < max_shift as i32
        && next_base_on_left_is_same(sequences, bounds)
        && last_base_on_right_is_same(sequences, bounds)
    {
        for b in bounds.iter_mut() {
            b.shift_left(1);
        }
        start_shift += 1;
        end_shift += 1;
    }
    (start_shift, end_shift)
}

fn last_base_on_right_is_same(sequences: &[&[u8]], bounds: &[IndexRange]) -> bool {
    // Empty allele (size 0): Java uses end-1 (== start-1), same as next-base-on-left.
    if bounds.iter().any(|b| b.end == 0) {
        return false;
    }
    let last = sequences[0][bounds[0].end - 1];
    sequences
        .iter()
        .zip(bounds)
        .all(|(seq, b)| seq[b.end - 1] == last)
}

fn first_base_on_left_is_same(sequences: &[&[u8]], bounds: &[IndexRange]) -> bool {
    let first = sequences[0][bounds[0].start];
    sequences
        .iter()
        .zip(bounds)
        .all(|(seq, b)| seq[b.start] == first)
}

fn next_base_on_left_is_same(sequences: &[&[u8]], bounds: &[IndexRange]) -> bool {
    // Guard: Java never reaches here when any start==0 because maxShift<=start and
    // startShift < maxShift together stop the loop first. Also reject wrapped starts
    // from historical usize underflow in IndexRange::shift_left.
    if bounds
        .iter()
        .zip(sequences.iter())
        .any(|(b, seq)| b.start == 0 || b.start > seq.len())
    {
        return false;
    }
    let next = sequences[0][bounds[0].start - 1];
    sequences
        .iter()
        .zip(bounds)
        .all(|(seq, b)| seq[b.start - 1] == next)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::smith_waterman::SwParameters;

    fn cigar_str(ref_seq: &str, alt: &str) -> String {
        let p = SwParameters::gatk_haplotype_to_reference();
        calculate_haplotype_cigar(ref_seq.as_bytes(), alt.as_bytes(), &p)
            .map(|c| c.to_gatk_string())
            .unwrap_or_else(|| "null".to_string())
    }

    #[test]
    fn matches_gatk_java_cigar_fixtures() {
        assert_eq!(cigar_str("ACGTTACGT", "ACGTTACGT"), "9M");
        assert_eq!(cigar_str("ACGTTACGT", "ACGTGACGT"), "9M");
        let ins = cigar_str("ACGTT", "ACGTTT");
        assert_eq!(ins, "3M1I2M");
        assert_eq!(cigar_str("ACGTT", "ACGT"), "3M1D1M");
        assert_eq!(cigar_str("ACGTT", "ACGTA"), "5M");
        assert_eq!(cigar_str("TTTACGTTACGT", "ACGTT"), "3D4M4D1M");
    }

    /// Regression: GIAB ci-subset panics when left-align reached `start==0` and
    /// indexed `start - 1` (usize underflow → index usize::MAX).
    #[test]
    fn normalize_alleles_does_not_underflow_at_sequence_start() {
        let ref_seq = b"ACGTACGT";
        let read_seq = b"ACGTACGT";
        // Empty indel ranges at start=0 with a large requested max_shift (Java would
        // reject via validateArg; we cap and must not panic).
        let mut bounds = [IndexRange::new(0, 0), IndexRange::new(0, 0)];
        let (start_shift, end_shift) = normalize_alleles(
            &[ref_seq.as_slice(), read_seq.as_slice()],
            &mut bounds,
            8,
            true,
        );
        assert_eq!(start_shift, 0);
        assert_eq!(end_shift, 0);
        assert_eq!(bounds[0].start, 0);

        // Homopolymer left-align that would walk into start without the signed/capped loop.
        let ref_h = b"AAAAA";
        let read_h = b"AAAA"; // 1D in A-run
        let mut bounds = [IndexRange::new(1, 2), IndexRange::new(1, 1)]; // D of one A
        let _ = normalize_alleles(
            &[ref_h.as_slice(), read_h.as_slice()],
            &mut bounds,
            10,
            true,
        );
        assert!(bounds.iter().all(|b| b.start <= b.end));
    }

    #[test]
    fn left_align_indel_at_read_start_no_panic() {
        let p = SwParameters::gatk_haplotype_to_reference();
        // Insertion near / at the left edge of the alignment block.
        let c = calculate_haplotype_cigar(b"ACGTACGTACGT", b"AACGTACGTACGT", &p);
        assert!(c.is_some());
        let c = calculate_haplotype_cigar(b"TTTTAAAA", b"TTTTAAA", &p);
        assert!(c.is_some());
    }
}
