//! Per-locus pileup snapshots attached to assembly regions (GAP-B-07 / B.5.8).

/// One locus pileup retained on an [`crate::assembly_region_iterator::AssemblyRegion`]
/// when pileup tracking is enabled (GATK `AlignmentAndReferenceContext`).
/// # Invariants
/// `depth` is the pileup read count at `pos` on `contig`.
/// `pos` is 1-based reference coordinate matching iterator pileup keys.
/// # Ownership
/// Owns contig name; stored in region-side pileup vectors when tracking is on.
/// # Mutation
/// Immutable snapshot attached to an assembly region after iterator emission.
/// # Biological assumptions
/// Sparse pileup depth summaries support RCM/activity diagnostics without full per-read pileup retention.
/// # Java equivalence
/// GATK `shouldTrackPileupsForAssemblyRegions` / `AlignmentAndReferenceContext` pileup attachment.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RegionPileupLocus {
    pub contig: String,
    pub pos: u64,
    pub depth: usize,
}
