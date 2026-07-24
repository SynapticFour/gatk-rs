//! Genome-wide indel evidence from BAM CIGARs (L7-B1).
//! SNP pileup AD is always 0 for indels; production genotyping and EventMap prune use these
//! helpers so dense GIAB indels are algorithmically recoverable without contig-2 shaping.

use crate::event_map::VariationEvent;
use rust_htslib::bam::{self, record::Cigar, record::CigarString};

fn query_subseq(seq: &[u8], start: usize, len: usize) -> Option<&[u8]> {
    seq.get(start..start.saturating_add(len))
}

/// Read pileup AD for a simple indel from BAM CIGAR I/D.
/// Counts reads whose CIGAR encodes the same insertion/deletion at the VCF anchor. Ref AD counts
/// reads that span the anchor with match ops and no matching indel at that locus.
pub(crate) fn read_indel_allele_depths_from_cigars(
    reads: &[bam::Record],
    event: &VariationEvent,
) -> (i32, i32) {
    let ref_bytes = event.ref_allele.as_bytes();
    let alt_bytes = event.alt_allele.as_bytes();
    if ref_bytes.is_empty() || alt_bytes.is_empty() || ref_bytes == alt_bytes {
        return (0, 0);
    }
    let is_ins = alt_bytes.len() > ref_bytes.len();
    let is_del = ref_bytes.len() > alt_bytes.len();
    if !is_ins && !is_del {
        return (0, 0);
    }
    let shared = ref_bytes
        .iter()
        .zip(alt_bytes.iter())
        .take_while(|(a, b)| a.eq_ignore_ascii_case(b))
        .count();
    if shared == 0 {
        return (0, 0);
    }
    let ins_seq: &[u8] = if is_ins { &alt_bytes[shared..] } else { &[] };
    // Suffix deletions (ALT is a proper prefix of REF) are written right-trimmed in VCF but
    // reads carry left-aligned D after the anchor — match CIGAR at shared=1 (L10 STR).
    let (del_len, shared_for_pos) = if is_del {
        let del_len = ref_bytes.len() - shared;
        let shared_for_pos = if alt_bytes.len() > 1 && ref_bytes.starts_with(alt_bytes) {
            1
        } else {
            shared
        };
        (del_len, shared_for_pos)
    } else {
        (0, shared)
    };
    if is_ins && ins_seq.is_empty() {
        return (0, 0);
    }
    if is_del && del_len == 0 {
        return (0, 0);
    }
    let indel_ref_pos0 = event.start_1based.get() as i64 - 1 + shared_for_pos as i64;
    let mut ref_count = 0i32;
    let mut alt_count = 0i32;
    for rec in reads {
        if rec.is_unmapped() || rec.tid() < 0 {
            continue;
        }
        let cigar = CigarString(rec.cigar().iter().copied().collect());
        let seq = rec.seq().as_bytes();
        let mut ref_pos0 = rec.pos();
        let mut query_pos: usize = 0;
        let mut saw_alt = false;
        let mut spans_anchor = false;
        for op in cigar.iter() {
            match op {
                Cigar::Ins(n) => {
                    let len = *n as usize;
                    if is_ins && ref_pos0 == indel_ref_pos0 && len == ins_seq.len() {
                        if let Some(ins) = query_subseq(&seq, query_pos, len) {
                            if ins.eq_ignore_ascii_case(ins_seq) {
                                saw_alt = true;
                            }
                        }
                    }
                    query_pos += len;
                }
                Cigar::Del(n) => {
                    let len = *n as usize;
                    if is_del && ref_pos0 == indel_ref_pos0 && len == del_len {
                        saw_alt = true;
                    }
                    ref_pos0 += len as i64;
                }
                Cigar::Match(n) | Cigar::Equal(n) | Cigar::Diff(n) => {
                    let len = *n as i64;
                    let start = ref_pos0;
                    let end = ref_pos0 + len;
                    let anchor0 = event.start_1based.get() as i64 - 1;
                    if start <= anchor0 && anchor0 < end {
                        spans_anchor = true;
                    }
                    ref_pos0 = end;
                    query_pos += *n as usize;
                }
                Cigar::SoftClip(n) => {
                    query_pos += *n as usize;
                }
                Cigar::HardClip(_) | Cigar::RefSkip(_) | Cigar::Pad(_) => {}
            }
        }
        // L11: long insertions may be encoded as M-span sequence (not CIGAR I) — plug match.
        if !saw_alt && is_ins && ins_seq.len() >= 10 {
            let anchor0 = event.start_1based.get() as i64 - 1;
            if let Some(qi) = crate::read_projection::query_index_at_reference_position(
                rec.pos(),
                &cigar,
                anchor0,
            ) {
                if seq
                    .get(qi)
                    .is_some_and(|b| b.eq_ignore_ascii_case(&ref_bytes[0]))
                {
                    if let Some(ins) = query_subseq(&seq, qi + 1, ins_seq.len()) {
                        if ins.eq_ignore_ascii_case(ins_seq) {
                            saw_alt = true;
                        }
                    }
                }
            }
        }
        if saw_alt {
            alt_count += 1;
        } else if spans_anchor {
            ref_count += 1;
        }
    }
    (ref_count, alt_count)
}

/// Deterministic haplotype CIGAR for a single left-aligned indel on `ref_bases`.
pub(crate) fn cigar_for_single_indel_event(
    ref_bases: &[u8],
    pad_start_1based: u64,
    event: &VariationEvent,
) -> Option<crate::cigar::Cigar> {
    use crate::cigar::{Cigar, CigarOperator};
    if !event.is_indel() {
        return None;
    }
    let off = event.start_1based.get().saturating_sub(pad_start_1based) as usize;
    let ref_len = event.ref_allele.len();
    let alt_len = event.alt_allele.len();
    if off + ref_len > ref_bases.len() {
        return None;
    }
    let slice = &ref_bases[off..off + ref_len];
    if !slice.eq_ignore_ascii_case(event.ref_allele.as_bytes()) {
        return None;
    }
    let trailing = ref_bases.len() - off - ref_len;
    let mut c = Cigar::new();
    if off > 0 {
        c.push(off, CigarOperator::Match);
    }
    if alt_len > ref_len {
        c.push(ref_len, CigarOperator::Match);
        c.push(alt_len - ref_len, CigarOperator::Insertion);
    } else if ref_len > alt_len {
        if alt_len > 0 {
            c.push(alt_len, CigarOperator::Match);
        }
        c.push(ref_len - alt_len, CigarOperator::Deletion);
    } else {
        return None;
    }
    if trailing > 0 {
        c.push(trailing, CigarOperator::Match);
    }
    Some(c)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::genome_loc::GenomePosition;

    #[test]
    fn cigar_indel_evidence_builds_insertion_cigar() {
        let ev = VariationEvent {
            contig: "20".into(),
            start_1based: GenomePosition::new_1based(9),
            end_1based: GenomePosition::new_1based(9),
            ref_allele: "A".into(),
            alt_allele: "AGG".into(),
        };
        // pad_start=1 → offset 8 lands on 'A'
        let ref_bases = b"NNNNNNNNANNNNNNNN";
        let cig = cigar_for_single_indel_event(ref_bases, 1, &ev).expect("cigar");
        assert!(!cig.elements.is_empty());
    }

    #[test]
    fn indel_ad_empty_reads_is_zero() {
        let ev = VariationEvent {
            contig: "20".into(),
            start_1based: GenomePosition::new_1based(100),
            end_1based: GenomePosition::new_1based(100),
            ref_allele: "AT".into(),
            alt_allele: "A".into(),
        };
        assert_eq!(read_indel_allele_depths_from_cigars(&[], &ev), (0, 0));
    }
}
