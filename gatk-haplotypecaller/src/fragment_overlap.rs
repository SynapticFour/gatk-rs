//! Overlapping paired-fragment base-quality correction.
//! Mirrors `AssemblyBasedCallerUtils.cleanOverlappingReadPairs` and
//! `FragmentUtils.adjustQualsOfOverlappingPairedFragments`.

#![warn(clippy::unwrap_used, clippy::expect_used)]

use crate::read_unclip::{
    alignment_end_1based, cigar_len, consumes_read_bases, consumes_ref_bases,
};
use gatk_common::{GatkError, GatkResult};
use rust_htslib::bam;
use rust_htslib::bam::record::{Aux, Cigar, CigarString};

/// SAM `OC`/`OP` written once before haplotype realignment overwrites CIGAR/POS.
/// RankSum `getElementForRead` retries this alignment when the attached CIGAR
/// does not cover `vc.getStart()` (GATK 4.4 `createReadAlignedToRef` keeps a
/// covering CIGAR on those reads; Rust's later haplotype CIGAR may not).
pub(crate) const RANK_SUM_ORIGINAL_CIGAR_TAG: &[u8] = b"OC";
pub(crate) const RANK_SUM_ORIGINAL_POS_TAG: &[u8] = b"OP";

/// `FragmentUtils.HALF_OF_DEFAULT_PCR_SNV_ERROR_QUAL` (phred(1e-4)/2 = 20).
pub const HALF_OF_DEFAULT_PCR_SNV_ERROR_QUAL: u8 = 20;

const FLAG_PAIRED: u16 = 0x1;
const FLAG_MATE_UNMAPPED: u16 = 0x8;

fn soft_start_1based(rec: &bam::Record) -> i32 {
    let mut soft_start = rec.pos() + 1;
    for c in rec.cigar().iter() {
        if let Cigar::SoftClip(sc_len) = c {
            soft_start -= i64::from(*sc_len);
        } else if !matches!(c, Cigar::HardClip(_)) {
            break;
        }
    }
    soft_start as i32
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum RefCoordIndex {
    NotFound,
    Clipping,
    Index(usize),
}

/// Soft-start from alignment pos0 + leading SoftClip (GATK `getSoftStart`).
fn soft_start_1based_from_cigar(alignment_pos0: i64, cigar: &CigarString) -> i32 {
    let mut soft_start = alignment_pos0 + 1;
    for c in cigar.iter() {
        if let Cigar::SoftClip(sc_len) = c {
            soft_start -= i64::from(*sc_len);
        } else if !matches!(c, Cigar::HardClip(_)) {
            break;
        }
    }
    soft_start as i32
}

/// `ReadUtils.getReadIndexForReferenceCoordinate` on cached CIGAR (SoftClip consumes ref).
///
/// Distinct from AD `query_index_at_reference_position` (BAM-aligned SoftClip = query only).
pub(crate) fn softclip_read_index_for_ref_1based(
    alignment_pos0: i64,
    cigar: &CigarString,
    ref_coord_1based: i32,
) -> RefCoordIndex {
    let alignment_start = soft_start_1based_from_cigar(alignment_pos0, cigar);
    if ref_coord_1based < alignment_start {
        return RefCoordIndex::NotFound;
    }
    let mut last_read = 0usize;
    let mut last_ref = alignment_start;
    let mut first_read;
    let mut first_ref;
    for c in cigar.iter() {
        let op_len = cigar_len(c) as i32;
        let op_consumes_read = consumes_read_bases(c);
        let op_consumes_ref = consumes_ref_bases(c) || matches!(c, Cigar::SoftClip(_));
        first_read = last_read;
        first_ref = last_ref;
        if op_consumes_read {
            last_read += op_len as usize;
        }
        if op_consumes_ref {
            last_ref += op_len;
        }
        if first_ref <= ref_coord_1based && ref_coord_1based < last_ref {
            if matches!(c, Cigar::SoftClip(_) | Cigar::HardClip(_)) {
                return RefCoordIndex::Clipping;
            }
            let read_pos = if op_consumes_read {
                first_read + (ref_coord_1based - first_ref) as usize
            } else {
                first_read
            };
            return RefCoordIndex::Index(read_pos);
        }
    }
    RefCoordIndex::NotFound
}

/// `ReadUtils.getReadIndexForReferenceCoordinate` (1-based ref coordinate).
fn read_index_for_ref_coord_1based(rec: &bam::Record, ref_coord_1based: i32) -> RefCoordIndex {
    softclip_read_index_for_ref_1based(rec.pos(), &rec.cigar(), ref_coord_1based)
}

/// Base at a 1-based reference coordinate, if aligned and not clipped.
pub fn read_base_at_ref_coord_1based(rec: &bam::Record, ref_coord_1based: i32) -> Option<u8> {
    match read_index_for_ref_coord_1based(rec, ref_coord_1based) {
        RefCoordIndex::Index(i) => rec.seq().as_bytes().get(i).copied(),
        _ => None,
    }
}

/// Persist the current POS/CIGAR before `createReadAlignedToRef` overwrites them.
pub(crate) fn preserve_pre_realign_alignment(rec: &mut bam::Record) {
    if rec.aux(RANK_SUM_ORIGINAL_CIGAR_TAG).is_ok() {
        return;
    }
    let oc = rec.cigar().to_string();
    let op = i32::try_from(rec.pos() + 1).unwrap_or(1);
    let _ = rec.push_aux(RANK_SUM_ORIGINAL_CIGAR_TAG, Aux::String(&oc));
    let _ = rec.push_aux(RANK_SUM_ORIGINAL_POS_TAG, Aux::I32(op));
}

fn i32_from_aux(aux: Aux<'_>) -> Option<i32> {
    match aux {
        Aux::I8(v) => Some(i32::from(v)),
        Aux::U8(v) => Some(i32::from(v)),
        Aux::I16(v) => Some(i32::from(v)),
        Aux::U16(v) => Some(i32::from(v)),
        Aux::I32(v) => Some(v),
        Aux::U32(v) => i32::try_from(v).ok(),
        _ => None,
    }
}

fn parse_hts_cigar(s: &str) -> Option<CigarString> {
    let mut ops = Vec::new();
    let mut n = 0u32;
    let mut saw_digit = false;
    for c in s.chars() {
        if let Some(d) = c.to_digit(10) {
            saw_digit = true;
            n = n.saturating_mul(10).saturating_add(d);
            continue;
        }
        if !saw_digit {
            return None;
        }
        let op = match c {
            'M' => Cigar::Match(n),
            'I' => Cigar::Ins(n),
            'D' => Cigar::Del(n),
            'N' => Cigar::RefSkip(n),
            'S' => Cigar::SoftClip(n),
            'H' => Cigar::HardClip(n),
            'P' => Cigar::Pad(n),
            '=' => Cigar::Equal(n),
            'X' => Cigar::Diff(n),
            _ => return None,
        };
        ops.push(op);
        n = 0;
        saw_digit = false;
    }
    if saw_digit || ops.is_empty() {
        return None;
    }
    Some(CigarString::from(ops))
}

/// POS0 + CIGAR from before haplotype realignment, if preserved.
pub(crate) fn pre_realign_alignment(rec: &bam::Record) -> Option<(i64, CigarString)> {
    let op = i32_from_aux(rec.aux(RANK_SUM_ORIGINAL_POS_TAG).ok()?)?;
    let oc = match rec.aux(RANK_SUM_ORIGINAL_CIGAR_TAG).ok()? {
        Aux::String(s) => s.to_string(),
        _ => return None,
    };
    let cigar = parse_hts_cigar(&oc)?;
    Some((i64::from(op) - 1, cigar))
}

fn alignment_end_1based_from_cigar(alignment_pos0: i64, cigar: &CigarString) -> i64 {
    let pos1 = alignment_pos0 + 1;
    let mut ref_len = 0i64;
    for c in cigar.iter() {
        if consumes_ref_bases(c) {
            ref_len += i64::from(cigar_len(c));
        }
    }
    if ref_len > 0 {
        pos1 + ref_len - 1
    } else {
        pos1
    }
}

/// GATK `ReadUtils.getReadBaseQualityAtReferenceCoordinate` on an explicit CIGAR.
pub(crate) fn read_base_quality_at_ref_coord_on_cigar(
    alignment_pos0: i64,
    cigar: &CigarString,
    qual: &[u8],
    ref_coord_1based: i32,
) -> Option<u8> {
    let start = alignment_pos0 + 1;
    let end = alignment_end_1based_from_cigar(alignment_pos0, cigar);
    if i64::from(ref_coord_1based) < start || end < i64::from(ref_coord_1based) {
        return None;
    }
    let alignment_start = soft_start_1based_from_cigar(alignment_pos0, cigar);
    if ref_coord_1based < alignment_start {
        return None;
    }
    let mut last_read = 0usize;
    let mut last_ref = alignment_start;
    for c in cigar.iter() {
        let op_len = cigar_len(c) as i32;
        let op_consumes_read = consumes_read_bases(c);
        let op_consumes_ref = consumes_ref_bases(c) || matches!(c, Cigar::SoftClip(_));
        let first_read = last_read;
        let first_ref = last_ref;
        if op_consumes_read {
            last_read += op_len as usize;
        }
        if op_consumes_ref {
            last_ref += op_len;
        }
        if first_ref <= ref_coord_1based && ref_coord_1based < last_ref {
            if !op_consumes_read {
                return None;
            }
            let idx = first_read + (ref_coord_1based - first_ref) as usize;
            return qual.get(idx).copied();
        }
    }
    None
}

/// GATK `ReadUtils.getReadBaseQualityAtReferenceCoordinate`.
///
/// Empty when `refCoord` is outside unclipped `getStart`/`getEnd`, or when the
/// CIGAR operator at that coordinate does not consume read bases (deletion).
pub fn read_base_quality_at_ref_coord_1based(
    rec: &bam::Record,
    ref_coord_1based: i32,
) -> Option<u8> {
    read_base_quality_at_ref_coord_on_cigar(rec.pos(), &rec.cigar(), rec.qual(), ref_coord_1based)
}

/// Softclip-aware base lookup on cached CIGAR/seq (same contract as
/// [`read_base_at_ref_coord_1based`]).
pub(crate) fn softclip_base_at_ref_1based_cached(
    alignment_pos0: i64,
    cigar: &CigarString,
    seq: &[u8],
    ref_coord_1based: i32,
) -> Option<u8> {
    match softclip_read_index_for_ref_1based(alignment_pos0, cigar, ref_coord_1based) {
        RefCoordIndex::Index(i) => seq.get(i).copied(),
        _ => None,
    }
}

fn mate_start_1based(rec: &bam::Record) -> i32 {
    if rec.flags() & FLAG_MATE_UNMAPPED != 0 {
        0
    } else {
        (rec.mpos() + 1) as i32
    }
}

/// `FragmentCollection.create` for sorted `GATKRead` lists.
/// Returns [`GatkError::read`] if `reads` are not sorted by alignment start
/// (unsorted BAM shards / bad caller input — must not panic).
pub fn overlapping_pairs_indices(reads: &[bam::Record]) -> GatkResult<Vec<(usize, usize)>> {
    let mut pairs = Vec::new();
    let mut pending: std::collections::HashMap<Vec<u8>, usize> = std::collections::HashMap::new();
    let mut last_start = -1i32;
    for (idx, rec) in reads.iter().enumerate() {
        let start_1based = (rec.pos() + 1) as i32;
        if start_1based < last_start {
            return Err(GatkError::read(format!(
                "fragment_overlap: reads must be sorted by alignment start ({start_1based} < {last_start})"
            )));
        }
        last_start = start_1based;
        let end_1based = alignment_end_1based(rec);
        if rec.flags() & FLAG_PAIRED == 0
            || rec.flags() & FLAG_MATE_UNMAPPED != 0
            || mate_start_1based(rec) == 0
            || mate_start_1based(rec) > end_1based
        {
            continue;
        }
        let qname = rec.qname().to_vec();
        if let Some(other) = pending.remove(&qname) {
            pairs.push((other, idx));
        } else {
            pending.insert(qname, idx);
        }
    }
    Ok(pairs)
}

/// `FragmentUtils.adjustQualsOfOverlappingPairedFragments` (SNV path only; no indel qual tags).
pub fn adjust_quals_of_overlapping_pair(
    first: &mut bam::Record,
    second: &mut bam::Record,
    set_conflicting_to_zero: bool,
    half_of_pcr_snv_qual: u8,
) {
    if first.qname() != second.qname() {
        return;
    }
    let first_end = alignment_end_1based(first);
    let second_start = (second.pos() + 1) as i32;
    if first_end < second_start || first.tid() != second.tid() {
        return;
    }
    let (left, right) = if soft_start_1based(first) < soft_start_1based(second) {
        (first, second)
    } else {
        (second, first)
    };
    let offset = match read_index_for_ref_coord_1based(left, (right.pos() + 1) as i32) {
        RefCoordIndex::Index(i) => i,
        _ => return,
    };
    let left_end = match read_index_for_ref_coord_1based(left, alignment_end_1based(left)) {
        RefCoordIndex::Index(i) => i,
        _ => return,
    };
    let right_offset = match read_index_for_ref_coord_1based(right, (right.pos() + 1) as i32) {
        RefCoordIndex::Index(i) => i,
        _ => return,
    };
    let right_end = match read_index_for_ref_coord_1based(right, alignment_end_1based(right)) {
        RefCoordIndex::Index(i) => i,
        _ => return,
    };
    let num_overlap = left_end
        .saturating_sub(offset)
        .min(right_end.saturating_sub(right_offset))
        + 1;
    let mut left_quals = left.qual().to_vec();
    let mut right_quals = right.qual().to_vec();
    let left_bases = left.seq().as_bytes();
    let right_bases = right.seq().as_bytes();
    for i in 0..num_overlap {
        let li = offset + i;
        let ri = right_offset + i;
        if li >= left_bases.len() || ri >= right_bases.len() {
            break;
        }
        if left_bases[li] == right_bases[ri] {
            left_quals[li] = left_quals[li].min(half_of_pcr_snv_qual);
            right_quals[ri] = right_quals[ri].min(half_of_pcr_snv_qual);
        } else if set_conflicting_to_zero {
            left_quals[li] = 0;
            right_quals[ri] = 0;
        }
    }
    set_record_quals(left, &left_quals);
    set_record_quals(right, &right_quals);
}

fn set_record_quals(rec: &mut bam::Record, quals: &[u8]) {
    let cigar_vec: Vec<Cigar> = rec.cigar().iter().copied().collect();
    let cigar = rust_htslib::bam::record::CigarString::from(cigar_vec);
    let seq_bytes = rec.seq().as_bytes();
    let qname = rec.qname().to_vec();
    let flags = rec.flags();
    let tid = rec.tid();
    let pos = rec.pos();
    let mapq = rec.mapq();
    let mpos = rec.mpos();
    let isize = rec.insert_size();
    rec.set(&qname, Some(&cigar), &seq_bytes, quals);
    rec.set_tid(tid);
    rec.set_pos(pos);
    rec.set_mapq(mapq);
    rec.set_flags(flags);
    rec.set_mpos(mpos);
    rec.set_insert_size(isize);
}

/// `AssemblyBasedCallerUtils.cleanOverlappingReadPairs` for one sample's read list.
pub fn clean_overlapping_read_pairs(
    records: &mut [bam::Record],
    set_conflicting_to_zero: bool,
) -> GatkResult<()> {
    let pairs = overlapping_pairs_indices(records)?;
    let half = HALF_OF_DEFAULT_PCR_SNV_ERROR_QUAL;
    for (a, b) in pairs {
        let (i, j) = if a < b { (a, b) } else { (b, a) };
        let (left, right) = records.split_at_mut(j);
        let first = &mut left[i];
        let second = &mut right[0];
        adjust_quals_of_overlapping_pair(first, second, set_conflicting_to_zero, half);
    }
    Ok(())
}

pub fn format_quals(rec: &bam::Record) -> String {
    rec.qual()
        .iter()
        .map(|q| q.to_string())
        .collect::<Vec<_>>()
        .join(",")
}

#[cfg(test)]
#[allow(clippy::unwrap_used, clippy::expect_used, clippy::panic)]
mod tests {
    use super::*;
    use rust_htslib::bam::record::{Cigar, CigarString};

    fn set_pair(
        r: &mut bam::Record,
        qname: &[u8],
        flags: u16,
        pos: i64,
        mpos: i64,
        cigar: &[Cigar],
        quals: &[u8],
    ) {
        let seq = vec![b'A'; quals.len()];
        r.set(qname, Some(&CigarString::from(cigar.to_vec())), &seq, quals);
        r.set_tid(0);
        r.set_mtid(0);
        r.set_pos(pos);
        r.set_mpos(mpos);
        r.set_flags(flags);
    }

    #[test]
    fn agree_overlap_caps_quals_at_twenty() {
        let quals = vec![30u8; 11];
        let mut r1 = bam::Record::new();
        set_pair(&mut r1, b"f", 65, 0, 3, &[Cigar::Match(11)], &quals);
        let mut r2 = bam::Record::new();
        set_pair(&mut r2, b"f", 129, 3, 0, &[Cigar::Match(11)], &quals);
        adjust_quals_of_overlapping_pair(&mut r1, &mut r2, true, 20);
        for i in 0..3 {
            assert_eq!(r1.qual()[i], 30);
        }
        for i in 3..11 {
            assert_eq!(r1.qual()[i], 20);
        }
        for i in 0..8 {
            assert_eq!(r2.qual()[i], 20);
        }
    }
}
