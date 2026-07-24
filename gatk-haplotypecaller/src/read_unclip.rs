//! GATK read soft-clip preparation before assembly.
//! Mirrors `AssemblyBasedCallerUtils.finalizeRegion` soft-clip branch and
//! `ReadClipper` / `CigarUtils.revertSoftClips`.

#![warn(clippy::unwrap_used, clippy::expect_used)]

use rust_htslib::bam;
use rust_htslib::bam::record::{Cigar, CigarString};

/// HC default: `dontUseSoftClippedBases=false`, `overrideSoftclipFragmentCheck=false`.
/// # Invariants
/// [`Self::use_hard_clip_soft_bases`] mirrors GATK soft-clip vs hard-clip decision per read flags.
/// Defaults match HC `AssemblyRegionArgumentCollection` soft-clip settings.
/// # Ownership
/// [`Copy`] policy passed into finalize/unclip helpers; records are mutably borrowed.
/// # Mutation
/// Policy immutable; read CIGAR/bases may be hard-clipped when policy demands.
/// # Biological assumptions
/// Soft-clipped bases are excluded when fragment size is ill-defined or `dontUseSoftClippedBases` is set.
/// # Java equivalence
/// GATK `AssemblyBasedCallerUtils.finalizeRegion` + `ReadClipper` / `dontUseSoftClippedBases`.
#[derive(Debug, Clone, Copy)]
pub struct HcSoftclipPolicy {
    pub dont_use_soft_clipped_bases: bool,
    pub override_softclip_fragment_check: bool,
}

impl HcSoftclipPolicy {
    pub fn haplotype_caller_defaults() -> Self {
        Self {
            dont_use_soft_clipped_bases: false,
            override_softclip_fragment_check: false,
        }
    }

    pub fn use_hard_clip_soft_bases(&self, rec: &bam::Record) -> bool {
        self.dont_use_soft_clipped_bases
            || !(self.override_softclip_fragment_check || has_well_defined_fragment_size(rec))
    }
}

/// `ReadUtils.hasWellDefinedFragmentSize`.
pub fn has_well_defined_fragment_size(rec: &bam::Record) -> bool {
    let flags = rec.flags();
    const PAIRED: u16 = 0x1;
    const UNMAPPED: u16 = 0x4;
    const MATE_UNMAPPED: u16 = 0x8;
    const REVERSE: u16 = 0x10;
    const MATE_REVERSE: u16 = 0x20;
    if rec.insert_size() == 0 {
        return false;
    }
    if flags & PAIRED == 0 {
        return false;
    }
    if flags & UNMAPPED != 0 || flags & MATE_UNMAPPED != 0 {
        return false;
    }
    let read_rev = flags & REVERSE != 0;
    let mate_rev = flags & MATE_REVERSE != 0;
    read_rev != mate_rev
}

pub(crate) fn cigar_len(c: &Cigar) -> u32 {
    match c {
        Cigar::Match(n)
        | Cigar::Ins(n)
        | Cigar::Del(n)
        | Cigar::SoftClip(n)
        | Cigar::HardClip(n)
        | Cigar::Equal(n)
        | Cigar::Diff(n)
        | Cigar::RefSkip(n)
        | Cigar::Pad(n) => *n,
    }
}

pub(crate) fn consumes_read_bases(c: &Cigar) -> bool {
    matches!(
        c,
        Cigar::Match(_) | Cigar::Ins(_) | Cigar::SoftClip(_) | Cigar::Equal(_) | Cigar::Diff(_)
    )
}

pub(crate) fn consumes_ref_bases(c: &Cigar) -> bool {
    matches!(
        c,
        Cigar::Match(_) | Cigar::Del(_) | Cigar::RefSkip(_) | Cigar::Equal(_) | Cigar::Diff(_)
    )
}

fn is_clipping_op(c: &Cigar) -> bool {
    matches!(c, Cigar::SoftClip(_) | Cigar::HardClip(_))
}

fn revert_soft_clips_cigar(cigar: &[Cigar]) -> Vec<Cigar> {
    cigar
        .iter()
        .map(|c| match c {
            Cigar::SoftClip(n) => Cigar::Match(*n),
            Cigar::Match(n) => Cigar::Match(*n),
            Cigar::Ins(n) => Cigar::Ins(*n),
            Cigar::Del(n) => Cigar::Del(*n),
            Cigar::HardClip(n) => Cigar::HardClip(*n),
            Cigar::Equal(n) => Cigar::Equal(*n),
            Cigar::Diff(n) => Cigar::Diff(*n),
            Cigar::RefSkip(n) => Cigar::RefSkip(*n),
            Cigar::Pad(n) => Cigar::Pad(*n),
        })
        .collect()
}

fn cigar_op_key(c: Cigar) -> u8 {
    match c {
        Cigar::Match(_) => 0,
        Cigar::Ins(_) => 1,
        Cigar::Del(_) => 2,
        Cigar::SoftClip(_) => 3,
        Cigar::HardClip(_) => 4,
        Cigar::Equal(_) => 5,
        Cigar::Diff(_) => 6,
        Cigar::RefSkip(_) => 7,
        Cigar::Pad(_) => 8,
    }
}

/// `CigarBuilder.make` — merge consecutive identical operators.
fn normalize_cigar(elems: Vec<Cigar>) -> Vec<Cigar> {
    let mut out = Vec::new();
    for c in elems {
        if c == Cigar::Match(0) {
            continue;
        }
        if let Some(last) = out.last_mut() {
            if cigar_op_key(*last) == cigar_op_key(c) {
                let merged = cigar_len(last) + cigar_len(&c);
                *last = op_with_len(c, merged);
                continue;
            }
        }
        out.push(c);
    }
    out
}

/// GATK `getSoftStart` / `getStart` (1-based ref coords; HTSJDK `POS` is 0-based).
/// GATK `GATKRead.getSoftStart` (1-based).
pub fn gatk_soft_start_1based(rec: &bam::Record) -> i64 {
    let mut soft_start = rec.pos() + 1;
    for c in rec.cigar().iter() {
        if let Cigar::SoftClip(sc_len) = c {
            soft_start -= i64::from(*sc_len);
        } else if !matches!(c, Cigar::HardClip(_)) {
            break;
        }
    }
    soft_start
}

fn replace_record_body(
    template: &bam::Record,
    cigar: &CigarString,
    seq: &[u8],
    qual: &[u8],
) -> bam::Record {
    let mut out = bam::Record::new();
    out.set(template.qname(), Some(cigar), seq, qual);
    out.set_tid(template.tid());
    out.set_pos(template.pos());
    out.set_mapq(template.mapq());
    out.set_flags(template.flags());
    out.set_mpos(template.mpos());
    out.set_insert_size(template.insert_size());
    out
}

/// `CigarUtils.clipCigar` (hard or soft clip operator).
fn clip_cigar(cigar: &[Cigar], start: usize, stop: usize, clip_op: Cigar) -> Vec<Cigar> {
    debug_assert!(matches!(clip_op, Cigar::SoftClip(_) | Cigar::HardClip(_)));
    let clip_left = start == 0;
    let mut new_cigar = Vec::new();
    let mut element_start = 0usize;
    for c in cigar {
        if matches!(c, Cigar::HardClip(_)) {
            new_cigar.push(*c);
            continue;
        }
        let element_end = element_start
            + if consumes_read_bases(c) {
                cigar_len(c) as usize
            } else {
                0
            };
        if element_end <= start || element_start >= stop {
            if consumes_read_bases(c) || (element_start != start && element_start != stop) {
                new_cigar.push(*c);
            }
        } else {
            let unclipped_len = if clip_left {
                element_end.saturating_sub(stop)
            } else {
                start.saturating_sub(element_start)
            };
            let clipped_len = cigar_len(c) as usize - unclipped_len;
            if unclipped_len == 0 {
                if consumes_read_bases(c) {
                    new_cigar.push(clip_op_with_len(clip_op, cigar_len(c)));
                }
            } else if clip_left {
                new_cigar.push(clip_op_with_len(clip_op, clipped_len as u32));
                new_cigar.push(op_with_len(*c, unclipped_len as u32));
            } else {
                new_cigar.push(op_with_len(*c, unclipped_len as u32));
                new_cigar.push(clip_op_with_len(clip_op, clipped_len as u32));
            }
        }
        element_start = element_end;
    }
    new_cigar
}

fn clip_op_with_len(op: Cigar, len: u32) -> Cigar {
    match op {
        Cigar::SoftClip(_) => Cigar::SoftClip(len),
        Cigar::HardClip(_) => Cigar::HardClip(len),
        _ => op,
    }
}

fn op_with_len(c: Cigar, len: u32) -> Cigar {
    match c {
        Cigar::Match(_) => Cigar::Match(len),
        Cigar::Ins(_) => Cigar::Ins(len),
        Cigar::Del(_) => Cigar::Del(len),
        Cigar::SoftClip(_) => Cigar::SoftClip(len),
        Cigar::HardClip(_) => Cigar::HardClip(len),
        Cigar::Equal(_) => Cigar::Equal(len),
        Cigar::Diff(_) => Cigar::Diff(len),
        Cigar::RefSkip(_) => Cigar::RefSkip(len),
        Cigar::Pad(_) => Cigar::Pad(len),
    }
}

/// `CigarUtils.alignmentStartShift` for left hard-clip (read-base index `num_clipped`).
fn alignment_start_shift(cigar: &[Cigar], num_clipped: usize) -> i32 {
    let mut ref_bases = 0i32;
    let mut element_start = 0usize;
    for c in cigar {
        if matches!(c, Cigar::HardClip(_)) {
            continue;
        }
        let element_end = element_start
            + if consumes_read_bases(c) {
                cigar_len(c) as usize
            } else {
                0
            };
        if element_end <= num_clipped {
            if consumes_ref_bases(c) {
                ref_bases += cigar_len(c) as i32;
            }
        } else if element_start < num_clipped {
            let clipped_len = (num_clipped - element_start) as u32;
            if consumes_ref_bases(c) {
                ref_bases += clipped_len as i32;
            }
            break;
        }
        element_start = element_end;
    }
    ref_bases
}

/// `ClippingOp.applyHardClipBases`.
fn apply_hard_clip_bases(rec: &bam::Record, start: usize, stop: usize) -> bam::Record {
    let seq = rec.seq().as_bytes();
    let qual = rec.qual().to_vec();
    let read_len = seq.len();
    if read_len == 0 {
        return rec.clone();
    }
    let clip_len = stop.saturating_sub(start).saturating_add(1);
    let new_length = read_len.saturating_sub(clip_len);
    if new_length == 0 {
        let mut out = bam::Record::new();
        out.set_qname(rec.qname());
        out.set_flags(rec.flags() | 0x4);
        return out;
    }
    let copy_start = if start == 0 { stop + 1 } else { 0 };
    let new_seq = &seq[copy_start..copy_start + new_length];
    let new_qual = &qual[copy_start..copy_start + new_length.min(qual.len())];
    const UNMAPPED: u16 = 0x4;
    let old_cigar: Vec<Cigar> = rec.cigar().iter().copied().collect();
    let cigar = if rec.flags() & UNMAPPED != 0 {
        CigarString::from(vec![])
    } else {
        let elems = clip_cigar(&old_cigar, start, stop + 1, Cigar::HardClip(0));
        CigarString::from(normalize_cigar(elems))
    };
    let mut out = replace_record_body(rec, &cigar, new_seq, new_qual);
    if start == 0 && out.flags() & UNMAPPED == 0 {
        let shift = alignment_start_shift(&old_cigar, stop + 1);
        out.set_pos(rec.pos() + i64::from(shift));
    }
    normalize_record_cigar(&mut out);
    out
}

/// Merge consecutive identical cigar ops (e.g. multiple hard clips from finalize steps).
pub fn normalize_record_cigar(rec: &mut bam::Record) {
    if rec.flags() & 0x4 != 0 {
        return;
    }
    let elems: Vec<Cigar> = rec.cigar().iter().copied().collect();
    let cigar = CigarString::from(normalize_cigar(elems));
    let qname = rec.qname().to_vec();
    let seq = rec.seq().as_bytes();
    let qual = rec.qual().to_vec();
    rec.set(&qname, Some(&cigar), &seq, &qual);
}

/// `ReadClipper.revertSoftClippedBases` + `ClippingOp.applyRevertSoftClippedBases`.
pub fn revert_soft_clipped_bases(rec: &bam::Record) -> bam::Record {
    let cigar_vec: Vec<Cigar> = rec.cigar().iter().copied().collect();
    if cigar_vec.is_empty()
        || (!cigar_vec.first().is_some_and(is_clipping_op)
            && !cigar_vec.last().is_some_and(is_clipping_op))
    {
        return rec.clone();
    }
    // GATK uses `read.getSoftStart` on the original read (before S→M), not the reverted cigar.
    let new_start = gatk_soft_start_1based(rec);
    let new_cigar = CigarString::from(normalize_cigar(revert_soft_clips_cigar(&cigar_vec)));
    let qual = rec.qual().to_vec();
    let seq_bytes = rec.seq().as_bytes();
    let mut out = replace_record_body(rec, &new_cigar, &seq_bytes, &qual);
    if new_start <= 0 {
        if out.flags() & 0x4 == 0 {
            out.set_pos(0);
        }
        let stop = (-new_start) as usize;
        out = apply_hard_clip_bases(&out, 0, stop);
        if out.flags() & 0x4 == 0 {
            out.set_pos(0);
        }
        return out;
    }
    out.set_pos(new_start - 1);
    out
}

/// `ReadClipper.hardClipSoftClippedBases`.
pub fn hard_clip_soft_clipped_bases(rec: &bam::Record) -> bam::Record {
    let seq_len = rec.seq().as_bytes().len();
    if seq_len == 0 {
        return rec.clone();
    }
    let mut read_index = 0usize;
    let mut cut_left = None::<usize>;
    let mut cut_right = None::<usize>;
    let mut right_tail = false;
    for c in rec.cigar().iter() {
        if matches!(c, Cigar::SoftClip(_)) {
            if right_tail {
                cut_right = Some(read_index);
            } else {
                cut_left = Some(read_index + cigar_len(c) as usize - 1);
            }
        } else if !matches!(c, Cigar::HardClip(_)) {
            right_tail = true;
        }
        if consumes_read_bases(c) {
            read_index += cigar_len(c) as usize;
        }
    }
    let mut out = rec.clone();
    if let Some(right) = cut_right {
        out = apply_hard_clip_bases(&out, right, seq_len.saturating_sub(1));
    }
    if let Some(left) = cut_left {
        out = apply_hard_clip_bases(&out, 0, left);
    }
    out
}

pub const ORIGINAL_SOFTCLIP_START_TAG: &[u8] = b"os";
pub const ORIGINAL_SOFTCLIP_END_TAG: &[u8] = b"oe";

/// GATK `ReadClipper.hardClipLowQualEnds` (`minBaseQualityScore - 1` in HC).
pub fn hard_clip_low_qual_ends(rec: &bam::Record, low_qual: u8) -> bam::Record {
    clip_low_qual_ends(rec, low_qual, Cigar::HardClip(0))
}

/// GATK `ReadClipper.softClipLowQualEnds` (`minBaseQualityScore - 1` in HC).
pub fn soft_clip_low_qual_ends(rec: &bam::Record, low_qual: u8) -> bam::Record {
    clip_low_qual_ends(rec, low_qual, Cigar::SoftClip(0))
}

fn clip_low_qual_ends(rec: &bam::Record, low_qual: u8, clip_op: Cigar) -> bam::Record {
    let hard = matches!(clip_op, Cigar::HardClip(_));
    let seq_len = rec.seq().as_bytes().len();
    if seq_len == 0 {
        return rec.clone();
    }
    let quals = rec.qual();
    let mut left = 0isize;
    let mut right = seq_len.saturating_sub(1) as isize;
    while right >= 0 && quals.get(right as usize).copied().unwrap_or(0) <= low_qual {
        right -= 1;
    }
    while (left as usize) < seq_len && quals.get(left as usize).copied().unwrap_or(0) <= low_qual {
        left += 1;
    }
    if left > right {
        return empty_unmapped_read(rec);
    }
    let left = left as usize;
    let right = right as usize;
    if hard {
        let mut out = rec.clone();
        if right < seq_len.saturating_sub(1) {
            out = apply_hard_clip_bases(&out, right.saturating_add(1), seq_len.saturating_sub(1));
        }
        if left > 0 {
            out = apply_hard_clip_bases(&out, 0, left.saturating_sub(1));
        }
        return out;
    }
    let old_cigar: Vec<Cigar> = rec.cigar().iter().copied().collect();
    let mut cigar = old_cigar;
    if right < seq_len.saturating_sub(1) {
        cigar = clip_cigar(&cigar, right.saturating_add(1), seq_len, clip_op);
    }
    if left > 0 {
        cigar = clip_cigar(&cigar, 0, left, clip_op);
    }
    let cigar = CigarString::from(normalize_cigar(cigar));
    let seq = rec.seq().as_bytes();
    let qual = rec.qual();
    let mut out = replace_record_body(rec, &cigar, &seq, qual);
    normalize_record_cigar(&mut out);
    out
}

/// First finalizeRegion soft-clip step (no tail-qual or adaptor clipping).
pub fn apply_hc_softclip_pre_step(
    rec: &bam::Record,
    policy: &HcSoftclipPolicy,
) -> (bam::Record, &'static str, Option<(i32, i32)>) {
    if policy.use_hard_clip_soft_bases(rec) {
        return (hard_clip_soft_clipped_bases(rec), "hard_clip", None);
    }
    let soft_start = (rec.pos() + 1) as i32;
    let soft_end = alignment_end_1based(rec);
    let out = revert_soft_clipped_bases(rec);
    if out.flags() & 0x4 == 0 {
        return (out, "revert", Some((soft_start, soft_end)));
    }
    (out, "revert", None)
}

pub fn alignment_end_1based(rec: &bam::Record) -> i32 {
    let pos1 = (rec.pos() + 1) as i32;
    let mut ref_len = 0i32;
    for c in rec.cigar().iter() {
        if consumes_ref_bases(c) {
            ref_len += cigar_len(c) as i32;
        }
    }
    if ref_len > 0 {
        pos1 + ref_len - 1
    } else {
        pos1
    }
}

/// GATK `ReadUtils.CANNOT_COMPUTE_ADAPTOR_BOUNDARY`.
pub const CANNOT_COMPUTE_ADAPTOR_BOUNDARY: i32 = i32::MIN;

/// GATK `ReadUtils.getAdaptorBoundary` (1-based reference coordinate).
pub fn adaptor_boundary_1based(rec: &bam::Record) -> Option<i64> {
    if !has_well_defined_fragment_size(rec) {
        return None;
    }
    const REVERSE: u16 = 0x10;
    if rec.flags() & REVERSE != 0 {
        // GATK `getMateStart - 1` (0-based SAM `mpos` = 1-based mate start − 1).
        Some(rec.mpos())
    } else {
        let insert = rec.insert_size().unsigned_abs() as i64;
        Some((rec.pos() + 1) + insert)
    }
}

fn is_inside_read_sam(rec: &bam::Record, ref_coord: i64) -> bool {
    let start = rec.pos() + 1;
    let end = i64::from(alignment_end_1based(rec));
    ref_coord >= start && ref_coord <= end
}

/// Returns read index and whether the cigar op at that ref coord consumes read bases.
fn read_index_and_op_for_ref_coord_1based(
    rec: &bam::Record,
    ref_coord_1based: i64,
) -> Option<(usize, bool)> {
    let alignment_start = gatk_soft_start_1based(rec);
    if ref_coord_1based < alignment_start {
        return None;
    }
    let mut last_read = 0usize;
    let mut last_ref = alignment_start;
    for c in rec.cigar().iter() {
        let op_len = i64::from(cigar_len(c));
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
            if matches!(c, Cigar::SoftClip(_) | Cigar::HardClip(_)) {
                return None;
            }
            let read_pos = if op_consumes_read {
                first_read + (ref_coord_1based - first_ref) as usize
            } else {
                first_read
            };
            return Some((read_pos, op_consumes_read));
        }
    }
    None
}

fn read_index_at_ref_for_left_tail(rec: &bam::Record, ref_coord_1based: i64) -> Option<usize> {
    // GATK left-tail hard clip uses `getReadIndexForReferenceCoordinate` then decrements stop only
    // when the coordinate falls in a deletion (non-read-consuming op).
    let (idx, op_consumes_read) = read_index_and_op_for_ref_coord_1based(rec, ref_coord_1based)?;
    Some(if op_consumes_read {
        idx
    } else {
        idx.saturating_sub(1)
    })
}

fn read_index_at_ref_or_none(rec: &bam::Record, ref_coord_1based: i64) -> Option<usize> {
    read_index_and_op_for_ref_coord_1based(rec, ref_coord_1based).map(|(i, _)| i)
}

/// GATK `ReadClipper.hardClipAdaptorSequence`.
pub fn hard_clip_adaptor_sequence(rec: &bam::Record) -> bam::Record {
    let Some(boundary) = adaptor_boundary_1based(rec) else {
        return rec.clone();
    };
    if !is_inside_read_sam(rec, boundary) {
        return rec.clone();
    }
    const REVERSE: u16 = 0x10;
    if rec.flags() & REVERSE != 0 {
        hard_clip_by_ref_left_tail(rec, boundary)
    } else {
        hard_clip_by_ref_right_tail(rec, boundary)
    }
}

fn hard_clip_by_ref_left_tail(rec: &bam::Record, ref_stop_1based: i64) -> bam::Record {
    let Some(stop_idx) = read_index_at_ref_for_left_tail(rec, ref_stop_1based) else {
        return rec.clone();
    };
    apply_hard_clip_bases(rec, 0, stop_idx)
}

fn hard_clip_by_ref_right_tail(rec: &bam::Record, ref_start_1based: i64) -> bam::Record {
    let Some(start_idx) = read_index_at_ref_or_none(rec, ref_start_1based) else {
        return rec.clone();
    };
    let len = rec.seq().as_bytes().len();
    if start_idx >= len {
        return empty_unmapped_read(rec);
    }
    apply_hard_clip_bases(rec, start_idx, len.saturating_sub(1))
}

fn empty_unmapped_read(rec: &bam::Record) -> bam::Record {
    let mut out = bam::Record::new();
    out.set_qname(rec.qname());
    out.set_flags(rec.flags() | 0x4);
    out
}

/// GATK `ReadClipper.hardClipToRegion` on the assembly region padded span (1-based inclusive).
pub fn hard_clip_to_region(
    rec: &bam::Record,
    ref_start_1based: u64,
    ref_stop_1based: u64,
) -> bam::Record {
    let rs = ref_start_1based as i64;
    let re = ref_stop_1based as i64;
    let aln_start = rec.pos() + 1;
    let aln_end = i64::from(alignment_end_1based(rec));
    if aln_start > re || aln_end < rs {
        return empty_unmapped_read(rec);
    }
    if aln_start < rs && aln_end > re {
        let mut r = hard_clip_by_ref_left_tail(rec, rs - 1);
        r = hard_clip_by_ref_right_tail(&r, re + 1);
        return r;
    }
    if aln_start < rs {
        return hard_clip_by_ref_left_tail(rec, rs - 1);
    }
    if aln_end > re {
        return hard_clip_by_ref_right_tail(rec, re + 1);
    }
    rec.clone()
}

pub(crate) fn read_has_positive_cigar_length(rec: &bam::Record) -> bool {
    rec.cigar()
        .iter()
        .any(|c| consumes_read_bases(c) && cigar_len(c) > 0)
}

#[cfg(test)]
#[allow(clippy::unwrap_used, clippy::expect_used, clippy::panic)]
mod read_unclip_tests {
    use super::*;
    use rust_htslib::bam::record::Cigar;

    #[test]
    fn normalize_merges_consecutive_hard_clips() {
        let elems = vec![Cigar::HardClip(1), Cigar::HardClip(30), Cigar::Match(219)];
        let out = normalize_cigar(elems);
        assert_eq!(out, vec![Cigar::HardClip(31), Cigar::Match(219)]);
    }
}
