//! GATK `PileupElement` flank flags for `ReferenceConfidenceModel.isAltBeforeAssembly`.

use crate::read_projection::{cigar_soft_clip_ends, query_index_at_reference_position};
use rust_htslib::bam::record::Cigar;

/// Context at a 0-based reference coordinate for one read (pre-assembly / HC activity).
/// # Invariants
/// Indel flank flags are mutually exclusive with `is_deletion` for a given ref coordinate.
/// `read_base` is uppercase ASCII when the read contributes a match base; `-` for indel contexts.
/// # Ownership
/// [`Copy`] value produced by [`pileup_element_flags_at_ref`]; borrows no BAM data.
/// # Mutation
/// Immutable per pileup lookup; callers use flags read-only.
/// # Biological assumptions
/// One read contributes at most one pileup element per reference position.
/// # Java equivalence
/// GATK `PileupElement` flank flags (`isBeforeDeletionStart`, soft-clip adjacency, etc.).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PileupElementFlags {
    pub read_base: u8,
    pub qual: u8,
    pub is_deletion: bool,
    pub is_before_deletion_start: bool,
    pub is_after_deletion_end: bool,
    pub is_before_insertion: bool,
    pub is_after_insertion: bool,
    pub is_next_to_soft_clip: bool,
}

/// Returns `None` when the read does not contribute a pileup element at `ref_pos0`.
pub fn pileup_element_flags_at_ref(
    alignment_start: i64,
    cigar: &[Cigar],
    seq: &[u8],
    qual: &[u8],
    ref_pos0: i64,
) -> Option<PileupElementFlags> {
    let mut r = alignment_start;
    let mut q: usize = 0;
    let (lead_sc, trail_sc) = cigar_soft_clip_ends_from_slice(cigar);

    for (op_idx, op) in cigar.iter().enumerate() {
        match *op {
            Cigar::SoftClip(n) | Cigar::HardClip(n) => {
                q += n as usize;
            }
            Cigar::Ins(n) => {
                // GATK: base before insertion is at r-1 (last ref base before I).
                if ref_pos0 == r.saturating_sub(1) {
                    return Some(PileupElementFlags {
                        read_base: b'-',
                        qual: 0,
                        is_deletion: false,
                        is_before_deletion_start: false,
                        is_after_deletion_end: false,
                        is_before_insertion: true,
                        is_after_insertion: false,
                        is_next_to_soft_clip: false,
                    });
                }
                if ref_pos0 == r {
                    return Some(PileupElementFlags {
                        read_base: b'-',
                        qual: 0,
                        is_deletion: false,
                        is_before_deletion_start: false,
                        is_after_deletion_end: false,
                        is_before_insertion: false,
                        is_after_insertion: true,
                        is_next_to_soft_clip: false,
                    });
                }
                q += n as usize;
            }
            Cigar::Del(n) => {
                let dr = n as i64;
                if ref_pos0 == r.saturating_sub(1) {
                    return Some(PileupElementFlags {
                        read_base: b'-',
                        qual: 0,
                        is_deletion: false,
                        is_before_deletion_start: true,
                        is_after_deletion_end: false,
                        is_before_insertion: false,
                        is_after_insertion: false,
                        is_next_to_soft_clip: false,
                    });
                }
                if ref_pos0 >= r && ref_pos0 < r + dr {
                    return Some(PileupElementFlags {
                        read_base: b'-',
                        qual: 0,
                        is_deletion: true,
                        is_before_deletion_start: false,
                        is_after_deletion_end: false,
                        is_before_insertion: false,
                        is_after_insertion: false,
                        is_next_to_soft_clip: false,
                    });
                }
                if ref_pos0 == r + dr {
                    return Some(PileupElementFlags {
                        read_base: b'-',
                        qual: 0,
                        is_deletion: false,
                        is_before_deletion_start: false,
                        is_after_deletion_end: true,
                        is_before_insertion: false,
                        is_after_insertion: false,
                        is_next_to_soft_clip: false,
                    });
                }
                r += dr;
            }
            Cigar::Match(n) | Cigar::Equal(n) | Cigar::Diff(n) => {
                let dr = n as i64;
                if ref_pos0 >= r && ref_pos0 < r + dr {
                    let qi = q + (ref_pos0 - r) as usize;
                    if qi >= seq.len() {
                        return None;
                    }
                    let read_base = seq[qi].to_ascii_uppercase();
                    let qv = qual.get(qi).copied().unwrap_or(0);
                    let next_op = cigar.get(op_idx + 1);
                    let prev_op = op_idx.checked_sub(1).and_then(|i| cigar.get(i));
                    // GATK LIBS: only last base of M before D / first base of M after D.
                    let is_before_deletion_start =
                        matches!(next_op, Some(Cigar::Del(_))) && ref_pos0 == r + dr - 1;
                    let is_after_deletion_end =
                        matches!(prev_op, Some(Cigar::Del(_))) && ref_pos0 == r;
                    let is_before_insertion =
                        matches!(next_op, Some(Cigar::Ins(_))) && ref_pos0 == r + dr - 1;
                    let is_after_insertion =
                        matches!(prev_op, Some(Cigar::Ins(_))) && ref_pos0 == r;
                    // GATK `PileupElement.isNextToSoftClip` (at start/end of current cigar + adjacent S).
                    let at_start_of_cigar = ref_pos0 == r;
                    let at_end_of_cigar = ref_pos0 == r + dr - 1;
                    let is_after_soft_clip =
                        at_start_of_cigar && matches!(prev_op, Some(Cigar::SoftClip(_)));
                    let is_before_soft_clip =
                        at_end_of_cigar && matches!(next_op, Some(Cigar::SoftClip(_)));
                    let is_next_to_soft_clip = is_after_soft_clip || is_before_soft_clip;
                    let _ = (lead_sc, trail_sc);
                    return Some(PileupElementFlags {
                        read_base,
                        qual: qv,
                        is_deletion: false,
                        is_before_deletion_start,
                        is_after_deletion_end,
                        is_before_insertion,
                        is_after_insertion,
                        is_next_to_soft_clip,
                    });
                }
                r += dr;
                q += n as usize;
            }
            Cigar::RefSkip(n) | Cigar::Pad(n) => {
                r += n as i64;
            }
        }
    }

    // Fallback: query-index path for odd CIGARs.
    let qi = query_index_at_reference_position(
        alignment_start,
        &cigar_slice_to_string(cigar),
        ref_pos0,
    )?;
    if qi >= seq.len() {
        return None;
    }
    let read_base = seq[qi].to_ascii_uppercase();
    let qv = qual.get(qi).copied().unwrap_or(0);
    let is_next_to_soft_clip = (qi == lead_sc as usize && lead_sc > 0)
        || (qi + 1 >= seq.len().saturating_sub(trail_sc as usize) && trail_sc > 0);
    Some(PileupElementFlags {
        read_base,
        qual: qv,
        is_deletion: false,
        is_before_deletion_start: false,
        is_after_deletion_end: false,
        is_before_insertion: false,
        is_after_insertion: false,
        is_next_to_soft_clip,
    })
}

fn cigar_soft_clip_ends_from_slice(cigar: &[Cigar]) -> (u32, u32) {
    use rust_htslib::bam::record::CigarString;
    cigar_soft_clip_ends(&CigarString::from(cigar.to_vec()))
}

fn cigar_slice_to_string(cigar: &[Cigar]) -> rust_htslib::bam::record::CigarString {
    rust_htslib::bam::record::CigarString::from(cigar.to_vec())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn deletion_span_and_flanks() {
        // 3M 2D 3M at pos 10: ref 10-12 match, 13-14 del, 15-17 match
        let cigar = [Cigar::Match(3), Cigar::Del(2), Cigar::Match(3)];
        let seq = b"AAAGGTTT";
        let qual = [30u8; 8];
        let start = 10i64;
        assert!(
            !pileup_element_flags_at_ref(start, &cigar, seq, &qual, 10)
                .unwrap()
                .is_before_deletion_start
        );
        assert!(
            !pileup_element_flags_at_ref(start, &cigar, seq, &qual, 11)
                .unwrap()
                .is_before_deletion_start
        );
        assert!(
            pileup_element_flags_at_ref(start, &cigar, seq, &qual, 12)
                .unwrap()
                .is_before_deletion_start
        );
        assert!(
            pileup_element_flags_at_ref(start, &cigar, seq, &qual, 13)
                .unwrap()
                .is_deletion
        );
        assert!(
            pileup_element_flags_at_ref(start, &cigar, seq, &qual, 14)
                .unwrap()
                .is_deletion
        );
        assert!(
            pileup_element_flags_at_ref(start, &cigar, seq, &qual, 15)
                .unwrap()
                .is_after_deletion_end
        );
        assert!(
            !pileup_element_flags_at_ref(start, &cigar, seq, &qual, 16)
                .unwrap()
                .is_after_deletion_end
        );
    }

    #[test]
    fn hom_ref_in_match_before_deletion_is_not_before_deletion_start() {
        let cigar = [Cigar::Match(8), Cigar::Del(2), Cigar::Match(3)];
        let seq = b"ACGTACGTGGGTTT";
        let qual = [30u8; 14];
        let start = 92307220i64;
        for pos in 92307220..92307227 {
            let f = pileup_element_flags_at_ref(start, &cigar, seq, &qual, pos).unwrap();
            assert!(!f.is_before_deletion_start, "pos {pos}");
            assert!(!f.is_deletion);
        }
        assert!(
            pileup_element_flags_at_ref(start, &cigar, seq, &qual, 92307227)
                .unwrap()
                .is_before_deletion_start
        );
    }

    #[test]
    fn insertion_before_after_flank_matches_gatk_rcm() {
        let cigar = [Cigar::Match(2), Cigar::Ins(2), Cigar::Match(4)];
        let seq = b"AAAACCCC";
        let qual = [30u8; 8];
        let start = 5i64;
        assert!(
            pileup_element_flags_at_ref(start, &cigar, seq, &qual, 6)
                .unwrap()
                .is_before_insertion
        );
        assert!(
            pileup_element_flags_at_ref(start, &cigar, seq, &qual, 7)
                .unwrap()
                .is_after_insertion
        );
    }
}
