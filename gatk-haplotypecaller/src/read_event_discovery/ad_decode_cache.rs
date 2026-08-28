//! Region-scoped CIGAR/seq decode cache for multi-pass pileup AD.
//!
//! Observable contract unchanged: same AD counts. Avoids repeated `rec.cigar()` /
//! `rec.seq().as_bytes()` allocations when `try_genotype` rescans the same reads.

use crate::shared_bam::SharedBamRecord;
use rust_htslib::bam::record::CigarString;
use std::cell::RefCell;
use std::collections::HashMap;

#[derive(Default)]
pub struct AdDecodeCache {
    /// Keyed by `Arc::as_ptr` of the BAM record.
    entries: HashMap<usize, (CigarString, Vec<u8>)>,
}

impl AdDecodeCache {
    pub fn clear(&mut self) {
        self.entries.clear();
    }

    #[inline]
    pub fn cigar_and_seq<'a>(&'a mut self, rec: &SharedBamRecord) -> (&'a CigarString, &'a [u8]) {
        let key = std::sync::Arc::as_ptr(rec) as usize;
        let entry = self.entries.entry(key).or_insert_with(|| {
            let cigar = CigarString(rec.cigar().iter().copied().collect());
            let seq = rec.seq().as_bytes();
            (cigar, seq)
        });
        (&entry.0, entry.1.as_slice())
    }

    /// Softclip-aware base at 1-based ref (GATK SoftClip-as-ref), reusing cached CIGAR/seq.
    ///
    /// Not interchangeable with AD `query_index_at_reference_position` (aligned SoftClip).
    #[inline]
    pub fn softclip_base_at_ref_1based(
        &mut self,
        rec: &SharedBamRecord,
        ref_coord_1based: i32,
    ) -> Option<u8> {
        let pos0 = rec.pos();
        let (cigar, seq) = self.cigar_and_seq(rec);
        crate::fragment_overlap::softclip_base_at_ref_1based_cached(
            pos0,
            cigar,
            seq,
            ref_coord_1based,
        )
    }
}

thread_local! {
    static AD_DECODE_CACHE: RefCell<AdDecodeCache> = RefCell::new(AdDecodeCache::default());
}

/// Clear the TLS AD decode cache (call once per genotyping region).
pub fn clear_ad_decode_cache() {
    AD_DECODE_CACHE.with(|c| c.borrow_mut().clear());
    crate::read_event_discovery::ad_result_memo::clear_ad_result_memo();
}

/// Borrow TLS cache for a pileup scan.
///
/// Nested callers (genotyping after allele-filter keep) must not panic if an outer
/// scan still holds the `RefCell`. Inner scans use a temporary map: AD counts are
/// unchanged; only CIGAR/seq decode may be repeated.
pub(crate) fn with_ad_decode_cache<R>(f: impl FnOnce(&mut AdDecodeCache) -> R) -> R {
    AD_DECODE_CACHE.with(|c| match c.try_borrow_mut() {
        Ok(mut cache) => f(&mut cache),
        Err(_) => {
            let mut nested = AdDecodeCache::default();
            f(&mut nested)
        }
    })
}
