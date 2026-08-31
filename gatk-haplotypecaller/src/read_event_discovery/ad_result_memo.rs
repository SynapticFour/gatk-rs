//! Memoize pileup AD / softclip scan *results* within a genotyping region.
//!
//! [`super::ad_decode_cache::AdDecodeCache`] only avoids re-decoding CIGAR/seq. Production `try_genotype`
//! still rescans the same `(reads, locus, pad, alleles)` many times. This cache
//! returns identical `(ref,alt)` or softclip counts so observable AD is unchanged.

use std::cell::RefCell;
use std::collections::HashMap;

#[derive(Clone, Copy, Hash, Eq, PartialEq)]
struct AdScanKey {
    reads_ptr: usize,
    reads_len: usize,
    loc: u64,
    pad: u64,
    ref_b: u8,
    alt_b: u8,
    /// 0 = per-read, 1 = QNAME-dedupe.
    mode: u8,
}

#[derive(Clone, Copy, Hash, Eq, PartialEq)]
struct SoftclipKey {
    reads_ptr: usize,
    reads_len: usize,
    loc: u64,
    margin: i32,
    alt_b: u8,
}

#[derive(Clone, Hash, Eq, PartialEq)]
struct IndelAdKey {
    reads_ptr: usize,
    reads_len: usize,
    loc: u64,
    ref_allele: String,
    alt_allele: String,
}

#[derive(Default)]
struct AdResultMemo {
    ad: HashMap<AdScanKey, (i32, i32)>,
    indel: HashMap<IndelAdKey, (i32, i32)>,
    softclip: HashMap<SoftclipKey, (i32, i32)>,
}

impl AdResultMemo {
    fn clear(&mut self) {
        self.ad.clear();
        self.indel.clear();
        self.softclip.clear();
    }
}

thread_local! {
    static AD_RESULT_MEMO: RefCell<AdResultMemo> = RefCell::new(AdResultMemo::default());
}

/// Clear AD result memo (call with [`super::clear_ad_decode_cache`] once per region).
pub fn clear_ad_result_memo() {
    AD_RESULT_MEMO.with(|c| c.borrow_mut().clear());
}

#[inline]
fn reads_key(reads: &[crate::shared_bam::SharedBamRecord]) -> (usize, usize) {
    (reads.as_ptr() as usize, reads.len())
}

/// Look up or compute SNP AD for a read slice.
pub fn memo_snp_ad(
    reads: &[crate::shared_bam::SharedBamRecord],
    loc: u64,
    pad: u64,
    ref_b: u8,
    alt_b: u8,
    dedupe_qname: bool,
    compute: impl FnOnce() -> (i32, i32),
) -> (i32, i32) {
    let (reads_ptr, reads_len) = reads_key(reads);
    let key = AdScanKey {
        reads_ptr,
        reads_len,
        loc,
        pad,
        ref_b,
        alt_b,
        mode: u8::from(dedupe_qname),
    };
    AD_RESULT_MEMO.with(|cell| {
        let memo = cell.borrow_mut();
        if let Some(&v) = memo.ad.get(&key) {
            return v;
        }
        drop(memo);
        let v = compute();
        AD_RESULT_MEMO.with(|cell| {
            cell.borrow_mut().ad.insert(key, v);
        });
        v
    })
}

/// Look up or compute indel CIGAR AD for a read slice (same event, same reads).
pub fn memo_indel_ad(
    reads: &[crate::shared_bam::SharedBamRecord],
    loc: u64,
    ref_allele: &str,
    alt_allele: &str,
    compute: impl FnOnce() -> (i32, i32),
) -> (i32, i32) {
    let (reads_ptr, reads_len) = reads_key(reads);
    let key = IndelAdKey {
        reads_ptr,
        reads_len,
        loc,
        ref_allele: ref_allele.to_owned(),
        alt_allele: alt_allele.to_owned(),
    };
    AD_RESULT_MEMO.with(|cell| {
        let memo = cell.borrow_mut();
        if let Some(&v) = memo.indel.get(&key) {
            return v;
        }
        drop(memo);
        let v = compute();
        AD_RESULT_MEMO.with(|cell| {
            cell.borrow_mut().indel.insert(key, v);
        });
        v
    })
}

/// Look up or compute softclip alt counts `(deduped, fragments)`.
pub fn memo_softclip_alt(
    reads: &[crate::shared_bam::SharedBamRecord],
    loc: u64,
    margin: i32,
    alt_b: u8,
    compute: impl FnOnce() -> (i32, i32),
) -> (i32, i32) {
    let (reads_ptr, reads_len) = reads_key(reads);
    let key = SoftclipKey {
        reads_ptr,
        reads_len,
        loc,
        margin,
        alt_b,
    };
    AD_RESULT_MEMO.with(|cell| {
        let memo = cell.borrow_mut();
        if let Some(&v) = memo.softclip.get(&key) {
            return v;
        }
        drop(memo);
        let v = compute();
        AD_RESULT_MEMO.with(|cell| {
            cell.borrow_mut().softclip.insert(key, v);
        });
        v
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn memo_returns_cached_without_recompute() {
        clear_ad_result_memo();
        let reads: &[crate::shared_bam::SharedBamRecord] = &[];
        let mut calls = 0u32;
        let a = memo_snp_ad(reads, 100, 1, b'A', b'T', false, || {
            calls += 1;
            (3, 5)
        });
        let b = memo_snp_ad(reads, 100, 1, b'A', b'T', false, || {
            calls += 1;
            (99, 99)
        });
        assert_eq!(a, (3, 5));
        assert_eq!(b, (3, 5));
        assert_eq!(calls, 1);
    }

    #[test]
    fn indel_memo_returns_cached_without_recompute() {
        clear_ad_result_memo();
        let reads: &[crate::shared_bam::SharedBamRecord] = &[];
        let mut calls = 0u32;
        let a = memo_indel_ad(reads, 200, "AT", "A", || {
            calls += 1;
            (7, 2)
        });
        let b = memo_indel_ad(reads, 200, "AT", "A", || {
            calls += 1;
            (0, 0)
        });
        assert_eq!(a, (7, 2));
        assert_eq!(b, (7, 2));
        assert_eq!(calls, 1);
    }
}
