//! Region-scoped RT extract cache: avoid rebuilding the same k-mer graph twice.
//!
//! # Invariants
//! - Active only between [`begin_assemble_region`] and [`end_assemble_region`].
//! - Cache key includes kmer + allow_lc/nu + before_remove mode (wrong key ⇒ wrong alleles).
//! - Stores haplotype lists only (graphs dropped after extract) so Peak-RSS stays lean.
//! - `empty_configured_kmers` records production configured k-mers whose before_remove
//!   extract already yielded no alts (RT-first miss) so supplement/merge_rt skip re-probe.
//!
//! # Java equivalence
//! Java pays one graph family per k-mer per region; this cache restores that cost shape
//! when Rust orchestration would otherwise rebuild the same before_remove extract.

use crate::haplotype::Haplotype;
use std::cell::RefCell;
use std::collections::{HashMap, HashSet};
use std::sync::Arc;

/// Extract identity for one RT before/after-remove haplotype pull.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub(crate) struct RtExtractKey {
    pub kmer_size: usize,
    pub allow_low_complexity: bool,
    pub allow_non_unique_ref: bool,
    pub before_remove_paths: bool,
}

thread_local! {
    static ACTIVE: RefCell<bool> = const { RefCell::new(false) };
    static EXTRACTS: RefCell<HashMap<RtExtractKey, Arc<Vec<Haplotype>>>> =
        RefCell::new(HashMap::new());
    static EMPTY_CONFIGURED: RefCell<HashSet<usize>> = RefCell::new(HashSet::new());
}

/// Start a fresh per-region assemble scope (call at `assemble_from_ref_and_reads` entry).
pub(crate) fn begin_assemble_region() {
    ACTIVE.with(|a| *a.borrow_mut() = true);
    EXTRACTS.with(|c| c.borrow_mut().clear());
    EMPTY_CONFIGURED.with(|s| s.borrow_mut().clear());
}

/// Drop cache before PairHMM / next region (Peak-RSS + correctness).
pub(crate) fn end_assemble_region() {
    EXTRACTS.with(|c| c.borrow_mut().clear());
    EMPTY_CONFIGURED.with(|s| s.borrow_mut().clear());
    ACTIVE.with(|a| *a.borrow_mut() = false);
}

fn is_active() -> bool {
    ACTIVE.with(|a| *a.borrow())
}

/// Return a cached extract when present.
pub(crate) fn get_cached(key: &RtExtractKey) -> Option<Vec<Haplotype>> {
    if !is_active() {
        return None;
    }
    EXTRACTS.with(|c| {
        c.borrow().get(key).map(|arc| {
            // CLONE: owned haplotype list for caller merge/dedup.
            Vec::clone(arc)
        })
    })
}

/// Store extract result for the remainder of this region assemble; returns `haps` by move.
pub(crate) fn put_cached(key: RtExtractKey, haps: Vec<Haplotype>) -> Vec<Haplotype> {
    if is_active() {
        let n = haps.len();
        EXTRACTS.with(|c| {
            // CLONE: cache retains a copy; caller keeps the original list.
            c.borrow_mut().insert(key, Arc::new(Vec::clone(&haps)));
        });
        crate::runtime_config::rss_trace_checkpoint(
            "rt_extract_cache_store",
            &format!(
                "kmer={} before_remove={} haps={n}",
                key.kmer_size, key.before_remove_paths
            ),
        );
    }
    haps
}

/// Mark a configured k-mer whose before_remove extract had no usable alts (RT-first miss).
pub(crate) fn mark_configured_kmer_empty(kmer_size: usize) {
    if !is_active() {
        return;
    }
    EMPTY_CONFIGURED.with(|s| {
        s.borrow_mut().insert(kmer_size);
    });
    crate::runtime_config::rss_trace_checkpoint(
        "rt_configured_empty",
        &format!("kmer={kmer_size}"),
    );
}

/// True when RT-first (or equivalent) already proved this configured k-mer empty.
pub(crate) fn configured_kmer_already_empty(kmer_size: usize) -> bool {
    if !is_active() {
        return false;
    }
    EMPTY_CONFIGURED.with(|s| s.borrow().contains(&kmer_size))
}
