//! Assembly-path read length gate.
//! Mirrors `HaplotypeCallerEngine.filterNonPassingReads` length check:
//! `AlignmentUtils.unclippedReadLength(read) < READ_LENGTH_FILTER_THRESHOLD`.

use rust_htslib::bam;
use rust_htslib::bam::record::Cigar;

/// `HaplotypeCallerEngine.READ_LENGTH_FILTER_THRESHOLD`.
pub const GATK_READ_LENGTH_FILTER_THRESHOLD: usize = 10;

/// `AlignmentUtils.unclippedReadLength` — read length minus soft-clip bases only.
pub fn unclipped_read_length(rec: &bam::Record) -> usize {
    let read_len = rec.seq().as_bytes().len();
    let mut soft_clipped = 0usize;
    for c in rec.cigar().iter() {
        if let Cigar::SoftClip(n) = c {
            soft_clipped += *n as usize;
        }
    }
    read_len.saturating_sub(soft_clipped)
}

/// Length portion of `filterNonPassingReads` (other filters are separate PRE items).
pub fn passes_read_length_filter(rec: &bam::Record) -> bool {
    unclipped_read_length(rec) >= GATK_READ_LENGTH_FILTER_THRESHOLD
}

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::record::{Cigar, CigarString};

    fn rec_with(cigar: &[Cigar], seq_len: usize) -> bam::Record {
        let seq = vec![b'A'; seq_len];
        let qual = vec![30u8; seq_len];
        let mut r = bam::Record::new();
        r.set(b"t", Some(&CigarString::from(cigar.to_vec())), &seq, &qual);
        r
    }

    #[test]
    fn unclipped_subtracts_soft_clips_only() {
        let r = rec_with(&[Cigar::SoftClip(3), Cigar::Match(7)], 10);
        assert_eq!(unclipped_read_length(&r), 7);
        assert!(!passes_read_length_filter(&r));
    }

    #[test]
    fn ten_match_bases_pass() {
        let r = rec_with(&[Cigar::Match(10)], 10);
        assert_eq!(unclipped_read_length(&r), 10);
        assert!(passes_read_length_filter(&r));
    }
}
