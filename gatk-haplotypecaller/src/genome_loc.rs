//! 1-based inclusive genomic interval (GATK `SimpleInterval` subset).
//! # Invariants
//! Coordinates are **1-based** and **inclusive** on both ends (`[start, end]`).
//! For non-empty spans, `end >= start`; length is `end.get - start.get + 1`
//! ([`GenomeLoc::reference_span_length`]).
//! Contig identity lives outside this type (callers pair `GenomeLoc` with a contig name).

/// 1-based genomic coordinate on a contig (contig name lives outside this type).
/// HC crate newtype for GATK/VCF **1-based locus** coordinates (contig name paired at call sites).
/// Distinct from [`gatk_core::GenomicPosition`], which uses a contig **index** (`u32`) plus position
/// for core I/O and interval types — do not re-export core's type from HC.
/// # Invariants
/// Value is **1-based** (VCF / GATK locus convention). Valid loci are `≥ 1`.
/// Does not encode contig; pair with a contig string at the call site.
/// # Ownership
/// [`Copy`] newtype wrapper around `u64`.
/// # Mutation
/// Immutable newtype; construct via [`Self::new_1based`] or [`Self::try_new_1based`].
/// # Biological assumptions
/// Single reference locus on one contig.
/// # Java equivalence
/// GATK / HTSJDK 1-based locus coordinate (paired with contig elsewhere).
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct GenomePosition(u64);

impl GenomePosition {
    /// Unchecked 1-based constructor for callers that already validated `pos ≥ 1`.
    #[inline]
    pub const fn new_1based(pos: u64) -> Self {
        Self(pos)
    }

    /// Rejects `0` (invalid as a 1-based VCF/GATK locus).
    #[inline]
    pub const fn try_new_1based(pos: u64) -> Option<Self> {
        if pos == 0 {
            None
        } else {
            Some(Self(pos))
        }
    }

    #[inline]
    pub const fn get(self) -> u64 {
        self.0
    }
}

impl From<GenomePosition> for u64 {
    #[inline]
    fn from(pos: GenomePosition) -> Self {
        pos.get()
    }
}

/// 1-based inclusive genomic interval on a contig (contig name is external).
/// # Invariants
/// Coordinates are **1-based** and **inclusive** on both ends (`[start, end]`).
/// For non-empty spans constructed with [`Self::try_new`], `end >= start`.
/// [`Self::new`] assumes the caller already guarantees `end >= start`.
/// # Ownership
/// [`Copy`] value type; contig identity is paired at call sites.
/// # Mutation
/// Immutable via public fields unless the caller replaces coordinates directly.
/// # Biological assumptions
/// Interval refers to reference coordinates on one contig (VCF/GATK convention).
/// # Java equivalence
/// GATK `SimpleInterval` subset (`GenomeLoc` / `Locatable` coordinate contract).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct GenomeLoc {
    pub start: GenomePosition,
    pub end: GenomePosition,
}

impl GenomeLoc {
    /// Construct from raw 1-based endpoints (caller guarantees `end >= start`).
    pub fn new(start_1based: u64, end_1based: u64) -> Self {
        Self {
            start: GenomePosition::new_1based(start_1based),
            end: GenomePosition::new_1based(end_1based),
        }
    }

    /// Construct from typed positions; rejects `end < start`.
    pub fn try_new(start: GenomePosition, end: GenomePosition) -> Option<Self> {
        if end.get() < start.get() {
            None
        } else {
            Some(Self { start, end })
        }
    }

    /// Raw 1-based start for wire / arithmetic edges.
    #[inline]
    pub fn start_1based(self) -> u64 {
        self.start.get()
    }

    /// Raw 1-based inclusive end for wire / arithmetic edges.
    #[inline]
    pub fn end_1based(self) -> u64 {
        self.end.get()
    }

    pub fn contains(&self, other: &GenomeLoc) -> bool {
        other.start.get() >= self.start.get() && other.end.get() <= self.end.get()
    }

    pub fn reference_span_length(&self) -> u64 {
        self.end.get().saturating_sub(self.start.get()) + 1
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn try_new_1based_rejects_zero() {
        assert!(GenomePosition::try_new_1based(0).is_none());
        assert_eq!(GenomePosition::try_new_1based(42).unwrap().get(), 42);
    }

    #[test]
    fn loc_fields_are_genome_positions() {
        let loc = GenomeLoc::new(10, 20);
        assert_eq!(loc.start.get(), 10);
        assert_eq!(loc.end.get(), 20);
        assert_eq!(loc.reference_span_length(), 11);
        assert!(GenomeLoc::try_new(
            GenomePosition::new_1based(20),
            GenomePosition::new_1based(10)
        )
        .is_none());
    }
}
