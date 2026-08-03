//! GATK `ReadFilterLibrary.MATE_ON_SAME_CONTIG_OR_NO_MAPPED_MATE`.

use rust_htslib::bam::Record;

/// GATK `MateOnSameContigOrNoMappedMateReadFilter.test`.
#[inline]
pub fn passes_mate_on_same_contig_or_no_mapped_mate(rec: &Record) -> bool {
    if !rec.is_paired() {
        return true;
    }
    if rec.is_mate_unmapped() {
        return true;
    }
    if rec.is_unmapped() {
        return true;
    }
    rec.tid() == rec.mtid()
}

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::record::{Cigar, CigarString};
    use rust_htslib::bam::{HeaderView, Record};
    use std::sync::Arc;

    fn paired_rec(mate_tid: i32, unmapped: bool, mate_unmapped: bool) -> Record {
        let mut r = Record::new();
        // rust-htslib 1.x: Record::set_header takes Arc<HeaderView> (was Rc).
        r.set_header(Arc::new(HeaderView::from_bytes(
            b"@HD\tVN:1.0\n@SQ\tSN:chr1\tLN:1000\n@SQ\tSN:chr2\tLN:1000\n",
        )));
        r.set(
            b"r",
            Some(&CigarString::from(vec![Cigar::Match(10)])),
            b"AAAAAAAAAA",
            &vec![30u8; 10],
        );
        r.set_tid(0);
        r.set_mtid(mate_tid);
        r.set_insert_size(100);
        r.set_paired();
        if unmapped {
            r.set_unmapped();
        }
        if mate_unmapped {
            r.set_mate_unmapped();
        } else {
            r.unset_mate_unmapped();
        }
        if !unmapped {
            r.unset_unmapped();
        }
        r
    }

    #[test]
    fn unpaired_passes() {
        let mut r = paired_rec(0, false, false);
        r.unset_paired();
        assert!(passes_mate_on_same_contig_or_no_mapped_mate(&r));
    }

    #[test]
    fn mate_unmapped_passes() {
        assert!(passes_mate_on_same_contig_or_no_mapped_mate(&paired_rec(
            1, false, true
        )));
    }

    #[test]
    fn mate_different_contig_fails() {
        let r = paired_rec(1, false, false);
        assert!(r.is_paired(), "fixture: paired");
        assert!(!r.is_unmapped(), "fixture: read mapped");
        assert!(
            !r.is_mate_unmapped(),
            "fixture: mate mapped tid={} mtid={}",
            r.tid(),
            r.mtid()
        );
        assert_ne!(r.tid(), r.mtid(), "fixture: different contigs");
        assert!(
            !passes_mate_on_same_contig_or_no_mapped_mate(&r),
            "paired mapped reads on different contigs must fail filter"
        );
    }

    #[test]
    fn unmapped_read_passes() {
        assert!(passes_mate_on_same_contig_or_no_mapped_mate(&paired_rec(
            1, true, false
        )));
    }

    #[test]
    fn mate_same_contig_passes() {
        assert!(passes_mate_on_same_contig_or_no_mapped_mate(&paired_rec(
            0, false, false
        )));
    }
}
