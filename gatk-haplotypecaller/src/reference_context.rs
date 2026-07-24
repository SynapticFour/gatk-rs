//! GATK `ReferenceContext` for assembly-region `apply` / `callRegion`.
//! Wraps a 1-based inclusive reference interval and lazily-loaded bases from FASTA, matching
//! `new ReferenceContext(reference, assemblyRegion.getPaddedSpan)` in Java.

use gatk_common::{GatkError, GatkResult};
use gatk_core::reference::{ReferenceWindowCache, SequenceDictionary};
/// Reference bases for one locatable interval (typically the assembly region padded span).
/// # Invariants
/// `window_start`/`window_end` clip `[start, end]` to contig bounds; `bases` length matches inclusive window.
/// Empty `bases` when the clipped window is inverted (`window_start > window_end`).
/// # Ownership
/// Owns contig name and reference byte vector loaded from FASTA cache.
/// # Mutation
/// Immutable after [`Self::from_interval`] construction.
/// # Biological assumptions
/// Uppercase reference bytes match padded assembly region span used in `apply` / `callRegion`.
/// # Java equivalence
/// GATK `ReferenceContext` over assembly region padded span.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ReferenceContext {
    pub contig: String,
    /// Core interval (1-based inclusive).
    pub start: u64,
    pub end: u64,
    /// Query window clipped to contig bounds (1-based inclusive).
    pub window_start: u64,
    pub window_end: u64,
    /// Bases for `[window_start, window_end]` inclusive.
    pub bases: Vec<u8>,
}

impl ReferenceContext {
    pub fn empty() -> Self {
        Self {
            contig: String::new(),
            start: 0,
            end: 0,
            window_start: 0,
            window_end: 0,
            bases: Vec::new(),
        }
    }

    /// Load reference bases for a closed 1-based interval, clipped to `contig` length.
    pub fn from_interval(
        dictionary: &SequenceDictionary,
        ref_cache: &mut ReferenceWindowCache,
        contig: &str,
        start1: u64,
        end1: u64,
    ) -> GatkResult<Self> {
        let contig_len = dictionary
            .contig(contig)
            .map(|c| c.length)
            .ok_or_else(|| GatkError::argument(format!("unknown contig {contig}")))?;
        let window_start = start1.max(1);
        let window_end = end1.min(contig_len);
        let bases = if window_start > window_end {
            Vec::new()
        } else {
            ref_cache
                .get_interval_bytes(dictionary, contig, window_start, window_end)?
                .to_vec()
        };
        Ok(Self {
            contig: contig.to_string(),
            start: start1,
            end: end1,
            window_start,
            window_end,
            bases,
        })
    }

    pub fn len(&self) -> usize {
        self.bases.len()
    }

    /// Uppercase ASCII reference string for parity dumps.
    pub fn bases_ascii(&self) -> String {
        self.bases
            .iter()
            .map(|&b| (b as char).to_ascii_uppercase())
            .collect()
    }

    /// Expected length for a contiguous window (inclusive ends).
    pub fn expected_window_len(&self) -> u64 {
        self.window_end
            .saturating_sub(self.window_start)
            .saturating_add(1)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn padded_span_length_matches_bases() {
        let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../parity/fixtures");
        let ref_fa = root.join("reference.fa");
        let dict = SequenceDictionary::from_fasta_path(&ref_fa).unwrap();
        let mut cache = ReferenceWindowCache::new(&ref_fa, 4);
        let ctx = ReferenceContext::from_interval(&dict, &mut cache, "chr1", 1, 11).unwrap();
        assert_eq!(ctx.len(), 11);
        assert_eq!(ctx.expected_window_len(), 11);
        assert_eq!(ctx.bases_ascii().len(), 11);
    }

    #[test]
    fn extended_span_loads_full_window() {
        let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../parity/fixtures");
        let ref_fa = root.join("reference.fa");
        let dict = SequenceDictionary::from_fasta_path(&ref_fa).unwrap();
        let mut cache = ReferenceWindowCache::new(&ref_fa, 2);
        let ctx = ReferenceContext::from_interval(&dict, &mut cache, "chr1", 1, 32).unwrap();
        assert_eq!(ctx.window_start, 1);
        assert_eq!(ctx.window_end, 32);
        assert_eq!(ctx.len(), 32);
    }
}
