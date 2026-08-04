//! Shared BAM records for shard/region ownership without deep clones.
//!
//! # Observable contract
//! Same read evidence as owned `bam::Record` lists (qname/pos/flags/bases/cigar).
//! Ownership only — not a downsampler.
//!
//! # Mutation
//! Callers that rewrite CIGAR/bases/MAPQ must use [`record_make_mut`] / [`BamRecordSlot::make_mut`].

use rust_htslib::bam;
use std::sync::Arc;

/// Arc-backed BAM record. [`Clone`] is cheap; mutate via [`record_make_mut`].
pub type SharedBamRecord = Arc<bam::Record>;

#[inline]
pub fn share_record(rec: bam::Record) -> SharedBamRecord {
    Arc::new(rec)
}

/// Unique mutable access (copy-on-write when still shared with the shard).
#[inline]
pub fn record_make_mut(rec: &mut SharedBamRecord) -> &mut bam::Record {
    Arc::make_mut(rec)
}

/// Wrap an owned record list after load / transform / downsample.
#[inline]
pub fn share_records(records: Vec<bam::Record>) -> Vec<SharedBamRecord> {
    records.into_iter().map(share_record).collect()
}

/// Uniform mutable access for owned or Arc-backed BAM records.
pub trait BamRecordSlot {
    fn as_record(&self) -> &bam::Record;
    fn make_mut(&mut self) -> &mut bam::Record;
}

impl BamRecordSlot for bam::Record {
    #[inline]
    fn as_record(&self) -> &bam::Record {
        self
    }
    #[inline]
    fn make_mut(&mut self) -> &mut bam::Record {
        self
    }
}

impl BamRecordSlot for SharedBamRecord {
    #[inline]
    fn as_record(&self) -> &bam::Record {
        self.as_ref()
    }
    #[inline]
    fn make_mut(&mut self) -> &mut bam::Record {
        Arc::make_mut(self)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn share_clone_is_cheap_arc() {
        let mut rec = bam::Record::new();
        rec.set_pos(42);
        let a = share_record(rec);
        let b = a.clone();
        assert!(Arc::ptr_eq(&a, &b));
        assert_eq!(Arc::strong_count(&a), 2);
        assert_eq!(a.pos(), 42);
    }

    #[test]
    fn make_mut_cow_when_shared() {
        let a = share_record(bam::Record::new());
        let mut b = a.clone();
        assert!(Arc::ptr_eq(&a, &b));
        record_make_mut(&mut b).set_pos(7);
        assert!(!Arc::ptr_eq(&a, &b));
        assert_eq!(b.pos(), 7);
        assert_eq!(a.pos(), -1);
    }
}
