//! HQ soft-clip statistics for HC activity (`AVERAGE_HQ_SOFTCLIPS_HQ_BASES_THRESHOLD` path).

use crate::read_binding::record_overlaps_closed_interval_1based;
use crate::read_model::{passes_hc_read_filters_with_header, ReadFilterParams};
use crate::read_projection::cigar_soft_clip_ends;
use rust_htslib::bam;
use rust_htslib::bam::record::Cigar;

/// `ReferenceConfidenceModel` HQ soft-clip path: `AlignmentUtils.countHighQualitySoftClips(..., threshold)`
/// with [`RCM_HQ_SOFT_CLIP_QUAL_THRESHOLD`] — base must be **strictly greater** than this (Java `>`).
pub const RCM_HQ_SOFT_CLIP_QUAL_THRESHOLD: u8 = 28;

/// Legacy threshold used by read-level heuristics outside the RCM accumulator (kept for older tests/tools).
pub const GATK_HQ_SOFT_CLIP_BASE_QUAL: u8 = 25;

/// GATK `AlignmentUtils.countHighQualitySoftClips`: total count of soft-clip bases with qual
/// **>** [`RCM_HQ_SOFT_CLIP_QUAL_THRESHOLD`], summed over every `S` CIGAR block.
pub fn count_high_quality_soft_clip_bases_rcm(rec: &bam::Record) -> u32 {
    let cigar = rec.cigar();
    let qual = rec.qual();
    if qual.is_empty() {
        return 0;
    }
    let mut align_pos = 0usize;
    let mut num_hq = 0u32;
    for op in cigar.iter() {
        let n = op.len() as usize;
        match *op {
            Cigar::SoftClip(_) => {
                for i in 0..n {
                    let q = qual.get(align_pos + i).copied().unwrap_or(0);
                    if q > RCM_HQ_SOFT_CLIP_QUAL_THRESHOLD {
                        num_hq += 1;
                    }
                }
                align_pos += n;
            }
            Cigar::HardClip(_) | Cigar::Pad(_) | Cigar::Del(_) | Cigar::RefSkip(_) => {}
            Cigar::Match(_) | Cigar::Ins(_) | Cigar::Diff(_) | Cigar::Equal(_) => {
                align_pos += n;
            }
        }
    }
    num_hq
}

/// Longest contiguous HQ soft-clip run on either end of the read (bases with qual >= [`GATK_HQ_SOFT_CLIP_BASE_QUAL`]).
pub fn max_hq_soft_clip_bases(rec: &bam::Record) -> u32 {
    let qual = rec.qual();
    let seq_len = rec.seq().len();
    if qual.is_empty() || seq_len == 0 {
        return 0;
    }
    let (lead, trail) = cigar_soft_clip_ends(&rec.cigar());
    let lead_hq = hq_prefix_len(qual, lead as usize);
    let trail_hq = hq_suffix_len(qual, trail as usize);
    lead_hq.max(trail_hq) as u32
}

fn hq_prefix_len(qual: &[u8], clip_len: usize) -> usize {
    let n = clip_len.min(qual.len());
    let mut k = 0;
    for &q in &qual[..n] {
        if q >= GATK_HQ_SOFT_CLIP_BASE_QUAL {
            k += 1;
        } else {
            break;
        }
    }
    k
}

fn hq_suffix_len(qual: &[u8], clip_len: usize) -> usize {
    let n = clip_len.min(qual.len());
    if n == 0 {
        return 0;
    }
    let start = qual.len() - n;
    let mut k = 0;
    for &q in qual[start..].iter().rev() {
        if q >= GATK_HQ_SOFT_CLIP_BASE_QUAL {
            k += 1;
        } else {
            break;
        }
    }
    k
}

/// GATK-style running mean of per-read max HQ soft-clip lengths at a 1-based locus.
pub fn hq_soft_clip_running_mean_at_locus(
    records: &[bam::Record],
    header: &bam::HeaderView,
    contig: &str,
    pos1: u64,
    read_filters: &ReadFilterParams,
) -> f64 {
    let mut sum = 0.0;
    let mut n = 0u32;
    for rec in records {
        if !passes_hc_read_filters_with_header(rec, header, read_filters) {
            continue;
        }
        if !record_overlaps_closed_interval_1based(rec, header, contig, pos1, pos1, read_filters) {
            continue;
        }
        let clip = max_hq_soft_clip_bases(rec);
        if clip > 0 {
            sum += clip as f64;
            n += 1;
        }
    }
    if n == 0 {
        0.0
    } else {
        sum / f64::from(n)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::header::HeaderRecord;
    use rust_htslib::bam::{Header, HeaderView, Record};

    fn header() -> HeaderView {
        let mut h = Header::new();
        h.push_record(
            HeaderRecord::new(b"SQ")
                .push_tag(b"SN", &"chr1")
                .push_tag(b"LN", &100),
        );
        HeaderView::from_header(&h)
    }

    #[test]
    fn two_soft_clip_bases_counted_on_leading_clip() {
        use rust_htslib::bam::record::{Cigar, CigarString};
        let hv = header();
        let seq = b"ACGTACGTAC";
        let qual = [30u8, 30, 10, 30, 30, 30, 30, 30, 30, 30];
        let mut rec = Record::new();
        rec.set(
            b"sc",
            Some(&CigarString::from(vec![
                Cigar::SoftClip(2),
                Cigar::Match(8),
            ])),
            seq,
            &qual,
        );
        rec.set_tid(hv.tid(b"chr1").unwrap() as i32);
        rec.set_pos(9);
        assert_eq!(max_hq_soft_clip_bases(&rec), 2);
    }
}
