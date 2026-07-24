//! Strongly typed biological identifiers for HaplotypeCaller.
//! These newtypes exist to make domain mix-ups a compile error:
//! haplotype index vs PL genotype index, 1-based genome locus vs pad-relative
//! offset, MAPQ score vs SAM unavailable (255), etc.
//! Prefer explicit constructors (`new` / `try_new`) over blanket `From` integer
//! conversions. Convert at BAM/VCF wire edges only.

use crate::genome_loc::GenomePosition;

pub use gatk_core::MappingQuality;

/// 1-based reference locus (alias of [`GenomePosition`]).
pub type ReferenceCoordinate = GenomePosition;

/// Index into a haplotype list for one active region.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct HaplotypeIndex(usize);

impl HaplotypeIndex {
    #[inline]
    pub const fn new(index: usize) -> Self {
        Self(index)
    }

    #[inline]
    pub const fn get(self) -> usize {
        self.0
    }
}

/// Index into an allele list at a site (`0` = REF).
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct AlleleIndex(usize);

impl AlleleIndex {
    #[inline]
    pub const fn new(index: usize) -> Self {
        Self(index)
    }

    #[inline]
    pub const fn get(self) -> usize {
        self.0
    }

    /// Reference allele index (always 0).
    #[inline]
    pub const fn reference() -> Self {
        Self(0)
    }
}

/// Index into a genotyping / PairHMM read list for one region.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct ReadIndex(usize);

impl ReadIndex {
    #[inline]
    pub const fn new(index: usize) -> Self {
        Self(index)
    }

    #[inline]
    pub const fn get(self) -> usize {
        self.0
    }
}

/// Index into a multi-sample list.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct SampleIndex(usize);

impl SampleIndex {
    #[inline]
    pub const fn new(index: usize) -> Self {
        Self(index)
    }

    #[inline]
    pub const fn get(self) -> usize {
        self.0
    }
}

/// Biallelic diploid genotype index: `0=REF/REF`, `1=REF/ALT`, `2=ALT/ALT`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct DiploidGenotypeIndex(u8);

impl DiploidGenotypeIndex {
    pub const HOM_REF: Self = Self(0);
    pub const HET: Self = Self(1);
    pub const HOM_ALT: Self = Self(2);

    /// Construct only `{0,1,2}`; rejects other values.
    #[inline]
    pub const fn try_new(index: u8) -> Option<Self> {
        if index <= 2 {
            Some(Self(index))
        } else {
            None
        }
    }

    #[inline]
    pub const fn get(self) -> u8 {
        self.0
    }

    #[inline]
    pub const fn as_usize(self) -> usize {
        self.0 as usize
    }
}

/// 0-based offset into padded reference / haplotype bytes (EventMap axis).
/// Not a genome locus — convert with pad start before VCF emit.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct PadOffset0(usize);

impl PadOffset0 {
    #[inline]
    pub const fn new(offset: usize) -> Self {
        Self(offset)
    }

    #[inline]
    pub const fn get(self) -> usize {
        self.0
    }
}

/// 0-based offset into read bases / qualities.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct ReadCoordinate(usize);

impl ReadCoordinate {
    #[inline]
    pub const fn new(offset: usize) -> Self {
        Self(offset)
    }

    #[inline]
    pub const fn get(self) -> usize {
        self.0
    }
}

/// Per-allele depth (AD element); non-negative.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct AlleleDepth(u32);

impl AlleleDepth {
    #[inline]
    pub const fn new(depth: u32) -> Self {
        Self(depth)
    }

    #[inline]
    pub const fn get(self) -> u32 {
        self.0
    }

    /// VCF FORMAT wire value (HTSJDK uses signed ints).
    #[inline]
    pub fn as_i32(self) -> i32 {
        i32::try_from(self.0).unwrap_or(i32::MAX)
    }

    /// Build from a non-negative VCF/FORMAT integer; negatives clamp to 0.
    #[inline]
    pub fn from_i32_saturating(v: i32) -> Self {
        Self(if v < 0 { 0 } else { v as u32 })
    }
}

/// Total read depth (DP); non-negative.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct ReadDepth(u32);

impl ReadDepth {
    #[inline]
    pub const fn new(depth: u32) -> Self {
        Self(depth)
    }

    #[inline]
    pub const fn get(self) -> u32 {
        self.0
    }

    #[inline]
    pub fn as_i32(self) -> i32 {
        i32::try_from(self.0).unwrap_or(i32::MAX)
    }

    #[inline]
    pub fn from_i32_saturating(v: i32) -> Self {
        Self(if v < 0 { 0 } else { v as u32 })
    }
}

/// Genotype quality (GQ); non-negative Phred-scaled.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct GenotypeQuality(u32);

impl GenotypeQuality {
    #[inline]
    pub const fn new(gq: u32) -> Self {
        Self(gq)
    }

    #[inline]
    pub const fn get(self) -> u32 {
        self.0
    }

    #[inline]
    pub fn as_i32(self) -> i32 {
        i32::try_from(self.0).unwrap_or(i32::MAX)
    }

    #[inline]
    pub fn from_i32_saturating(v: i32) -> Self {
        Self(if v < 0 { 0 } else { v as u32 })
    }
}

/// One PL vector element; non-negative Phred-scaled likelihood.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct PhredLikelihood(u32);

impl PhredLikelihood {
    #[inline]
    pub const fn new(pl: u32) -> Self {
        Self(pl)
    }

    #[inline]
    pub const fn get(self) -> u32 {
        self.0
    }

    #[inline]
    pub fn as_i32(self) -> i32 {
        i32::try_from(self.0).unwrap_or(i32::MAX)
    }

    #[inline]
    pub fn from_i32_saturating(v: i32) -> Self {
        Self(if v < 0 { 0 } else { v as u32 })
    }
}

/// Sample / evaluation ploidy (`≥ 1`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Ploidy(u8);

impl Ploidy {
    /// Diploid sample ploidy (GATK HC default).
    pub const DIPLOID: Self = Self(2);

    #[inline]
    pub const fn try_new(ploidy: u8) -> Option<Self> {
        if ploidy >= 1 {
            Some(Self(ploidy))
        } else {
            None
        }
    }

    /// Unchecked constructor for callers that already validated ploidy ≥ 1.
    #[inline]
    pub const fn new_unchecked(ploidy: u8) -> Self {
        Self(ploidy)
    }

    #[inline]
    pub const fn get(self) -> u8 {
        self.0
    }

    #[inline]
    pub const fn as_u32(self) -> u32 {
        self.0 as u32
    }
}

/// Assembly k-mer size (`≥ 2`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct KmerSize(u16);

impl KmerSize {
    /// Default assembly k-mer size used by [`crate::assembly::AssemblyGraphParams::default`].
    pub const DEFAULT_ASSEMBLY: Self = Self(11);

    #[inline]
    pub const fn try_new(k: u16) -> Option<Self> {
        if k >= 2 {
            Some(Self(k))
        } else {
            None
        }
    }

    /// Unchecked constructor for callers that already validated k ≥ 2.
    #[inline]
    pub const fn new_unchecked(k: u16) -> Self {
        Self(k)
    }

    /// Fallible conversion from wire/`usize` k with a user-facing argument error.
    pub fn try_from_usize(k: usize) -> Result<Self, gatk_common::GatkError> {
        u16::try_from(k)
            .ok()
            .and_then(Self::try_new)
            .ok_or_else(|| {
                gatk_common::GatkError::invalid_argument(
                    "kmer_size",
                    format!("kmer_size must be ≥ 2 and ≤ {}, got {k}", u16::MAX),
                )
            })
    }

    #[inline]
    pub const fn get(self) -> u16 {
        self.0
    }

    #[inline]
    pub const fn as_usize(self) -> usize {
        self.0 as usize
    }
}

impl std::fmt::Display for HaplotypeIndex {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl std::fmt::Display for ReadIndex {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl std::fmt::Display for DiploidGenotypeIndex {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl std::fmt::Display for AlleleDepth {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}
impl std::fmt::Display for ReadDepth {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}
impl std::fmt::Display for GenotypeQuality {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}
impl std::fmt::Display for PhredLikelihood {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}

/// Compile-time guarantee: biological newtypes stay zero-cost over their wire integers.
const _: () = {
    use crate::genome_loc::GenomePosition;
    use core::mem::size_of;

    assert!(size_of::<GenomePosition>() == size_of::<u64>());
    assert!(size_of::<HaplotypeIndex>() == size_of::<usize>());
    assert!(size_of::<AlleleIndex>() == size_of::<usize>());
    assert!(size_of::<SampleIndex>() == size_of::<usize>());
    assert!(size_of::<ReadIndex>() == size_of::<usize>());
    assert!(size_of::<AlleleDepth>() == size_of::<u32>());
    assert!(size_of::<ReadDepth>() == size_of::<u32>());
    assert!(size_of::<GenotypeQuality>() == size_of::<u32>());
    assert!(size_of::<PhredLikelihood>() == size_of::<u32>());
    // MappingQuality / BaseQuality live in gatk-core; asserted there.
};

#[cfg(test)]
mod tests {
    use super::*;
    use crate::genome_loc::GenomePosition;

    #[test]
    fn diploid_genotype_index_rejects_out_of_range() {
        assert_eq!(
            DiploidGenotypeIndex::try_new(0),
            Some(DiploidGenotypeIndex::HOM_REF)
        );
        assert_eq!(
            DiploidGenotypeIndex::try_new(2),
            Some(DiploidGenotypeIndex::HOM_ALT)
        );
        assert_eq!(DiploidGenotypeIndex::try_new(3), None);
    }

    #[test]
    fn mapping_quality_unavailable_is_sam_255() {
        let mq = MappingQuality::from_sam_mapq(255);
        assert_eq!(mq, MappingQuality::Unavailable);
        assert_eq!(mq.score(), None);
        assert_eq!(mq.as_sam_mapq(), 255);
        assert_eq!(MappingQuality::from_sam_mapq(60).score(), Some(60));
    }

    #[test]
    fn ploidy_and_kmer_reject_zero() {
        assert!(Ploidy::try_new(0).is_none());
        assert_eq!(Ploidy::try_new(2).map(|p| p.get()), Some(2));
        assert_eq!(Ploidy::DIPLOID.get(), 2);
        assert!(KmerSize::try_new(0).is_none());
        assert!(KmerSize::try_new(1).is_none());
        assert_eq!(KmerSize::try_new(10).map(|k| k.as_usize()), Some(10));
        assert_eq!(KmerSize::DEFAULT_ASSEMBLY.as_usize(), 11);
        assert!(KmerSize::try_from_usize(1).is_err());
        assert_eq!(KmerSize::try_from_usize(10).unwrap().as_usize(), 10);
    }

    #[test]
    fn allele_depth_saturates_negative_wire_values() {
        assert_eq!(AlleleDepth::from_i32_saturating(-3).get(), 0);
        assert_eq!(AlleleDepth::new(7).as_i32(), 7);
    }

    #[test]
    fn reference_coordinate_aliases_genome_position() {
        let p: ReferenceCoordinate = GenomePosition::new_1based(100);
        assert_eq!(p.get(), 100);
    }

    #[test]
    fn index_domains_are_distinct_types() {
        let h = HaplotypeIndex::new(1);
        let a = AlleleIndex::new(1);
        let r = ReadIndex::new(1);
        assert_eq!(h.get(), a.get());
        assert_eq!(h.get(), r.get());
        // Distinct newtypes: cannot assign across domains without.get/new.
        let _ = (h, a, r);
    }

    #[test]
    fn genome_position_try_new_rejects_zero() {
        assert!(GenomePosition::try_new_1based(0).is_none());
        assert_eq!(
            GenomePosition::try_new_1based(1).map(GenomePosition::get),
            Some(1)
        );
    }
}
