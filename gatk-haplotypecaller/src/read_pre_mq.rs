//! Assembly-path mapping-quality gate.
//! Mirrors `HaplotypeCallerEngine.filterNonPassingReads`:
//! `rec.getMappingQuality < hcArgs.mappingQualityThreshold`.

use gatk_core::MappingQuality;
use rust_htslib::bam;

/// `HaplotypeCallerEngine.DEFAULT_READ_QUALITY_FILTER_THRESHOLD` / `hcArgs.mappingQualityThreshold` default.
pub const GATK_ASSEMBLY_MQ_FILTER_THRESHOLD: u8 = 20;

/// MQ portion of `filterNonPassingReads` (strict `<` in Java → keep when `mapq >= threshold`).
/// SAM MAPQ `255` ([`MappingQuality::Unavailable`]) is treated as passing, matching Java's
/// unsigned comparison (`255 >= threshold` for normal HC thresholds).
#[inline]
pub fn passes_assembly_mq_filter(rec: &bam::Record, threshold: u8) -> bool {
    match MappingQuality::from_sam_mapq(rec.mapq()) {
        MappingQuality::Unavailable => true,
        MappingQuality::Score(q) => q >= threshold,
    }
}

/// Default HC assembly-path threshold (20).
#[allow(dead_code)] // convenience wrapper; call sites pass an explicit threshold today
#[inline]
pub fn passes_assembly_mq_filter_default(rec: &bam::Record) -> bool {
    passes_assembly_mq_filter(rec, GATK_ASSEMBLY_MQ_FILTER_THRESHOLD)
}

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::record::{Cigar, CigarString};

    fn rec_with_mapq(mapq: u8) -> bam::Record {
        let mut r = bam::Record::new();
        r.set(
            b"t",
            Some(&CigarString::from(vec![Cigar::Match(10)])),
            b"AAAAAAAAAA",
            &vec![30u8; 10],
        );
        r.set_mapq(mapq);
        r
    }

    #[test]
    fn threshold_twenty_is_inclusive() {
        assert!(passes_assembly_mq_filter_default(&rec_with_mapq(20)));
        assert!(!passes_assembly_mq_filter_default(&rec_with_mapq(19)));
    }

    #[test]
    fn mq_unavailable_passes_assembly_path() {
        assert!(passes_assembly_mq_filter_default(&rec_with_mapq(255)));
        assert_eq!(
            MappingQuality::from_sam_mapq(255),
            MappingQuality::Unavailable
        );
    }
}
