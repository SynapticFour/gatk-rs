//! Shared BAM records for shard/region ownership without deep clones.
//!
//! # Observable contract
//! Same read evidence as owned `bam::Record` lists (qname/pos/flags/bases/cigar).
//! Ownership only — not a downsampler.
//!
//! # Mutation
//! Callers that rewrite CIGAR/bases/MAPQ must use [`record_make_mut`] / [`BamRecordSlot::make_mut`].

use rust_htslib::bam;
use std::sync::{Arc, OnceLock};

/// Arc-backed BAM record. [`Clone`] is cheap; mutate via [`record_make_mut`].
pub type SharedBamRecord = Arc<bam::Record>;

/// Process-wide empty placeholder for progressive shard release.
/// Replacing spent slots with this sentinel (instead of `Arc::new(Record::new())` per index)
/// avoids allocating millions of tiny empty records while keeping index maps stable.
static EMPTY_SHARED_BAM: OnceLock<SharedBamRecord> = OnceLock::new();

#[inline]
pub fn share_record(rec: bam::Record) -> SharedBamRecord {
    Arc::new(rec)
}

/// Process-wide empty placeholder reference (no clone).
#[inline]
pub fn empty_shared_record_ref() -> &'static SharedBamRecord {
    EMPTY_SHARED_BAM.get_or_init(|| Arc::new(bam::Record::new()))
}

/// Cheap clone of the shared empty BAM placeholder (see [`EMPTY_SHARED_BAM`]).
#[inline]
pub fn empty_shared_record() -> SharedBamRecord {
    empty_shared_record_ref().clone()
}

/// True when `rec` is the progressive-release empty sentinel (not a real alignment).
#[inline]
pub fn is_empty_shared_record(rec: &SharedBamRecord) -> bool {
    Arc::ptr_eq(rec, empty_shared_record_ref())
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

/// Take unique ownership of shared BAM records when possible (no CoW clone if strong_count==1).
#[inline]
pub fn into_unique_records(reads: Vec<SharedBamRecord>) -> Vec<bam::Record> {
    reads
        .into_iter()
        .map(|arc| match Arc::try_unwrap(arc) {
            Ok(rec) => rec,
            Err(shared) => (*shared).clone(),
        })
        .collect()
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

    #[test]
    fn empty_shared_record_is_one_sentinel() {
        let a = empty_shared_record();
        let b = empty_shared_record();
        assert!(Arc::ptr_eq(&a, &b));
        assert_eq!(a.pos(), -1);
    }
}
