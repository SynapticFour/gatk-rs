//! Core data types for genomic analysis

use indexmap::IndexMap;
use serde::{Deserialize, Serialize};

/// Represents a genomic position on a reference contig.
/// # Invariants
/// `position` is 1-based and must refer to a valid base within the contig when used with a dictionary.
/// `contig` is an opaque numeric index; callers map names via [`crate::reference::SequenceDictionary`].
/// # Ownership
/// `Copy` value type; no heap allocation; safe to pass by value or borrow immutably.
/// # Mutation
/// Fields are public; treat as immutable after construction unless building ad-hoc fixtures.
/// # Biological assumptions
/// Single reference coordinate on one contig; strand is not encoded here.
/// # Java equivalence
/// Approximates `org.broadinstitute.hellbender.utils.SimpleInterval` position facet / `Locatable` (contig + start).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct GenomicPosition {
    pub contig: u32,   // Reference sequence identifier
    pub position: u64, // 1-based position
}

/// Closed genomic interval on a reference contig (1-based, inclusive endpoints).
/// # Invariants
/// `start <= end` when constructed via [`GenomicInterval::new`]; callers must preserve this.
/// Coordinates are 1-based inclusive, matching GATK `-L` interval semantics.
/// # Ownership
/// Owns nothing beyond plain scalars; [`Clone`] for independent copies.
/// # Mutation
/// Public fields allow in-place edits; prefer constructors for validated intervals.
/// # Biological assumptions
/// Interval spans contiguous reference bases on one contig; no strand or padding semantics.
/// # Java equivalence
/// Approximates `org.broadinstitute.hellbender.utils.SimpleInterval`.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct GenomicInterval {
    pub contig: u32,
    pub start: u64, // 1-based, inclusive
    pub end: u64,   // 1-based, inclusive
}

impl GenomicInterval {
    pub fn new(contig: u32, start: u64, end: u64) -> Self {
        Self { contig, start, end }
    }

    pub fn contains(&self, pos: GenomicPosition) -> bool {
        pos.contig == self.contig && self.start <= pos.position && pos.position <= self.end
    }

    pub fn contains_position(&self, pos: GenomicPosition) -> bool {
        self.contains(pos)
    }

    pub fn length(&self) -> u64 {
        self.end - self.start + 1
    }
}

impl GenomicPosition {
    pub fn new(contig: u32, position: u64) -> Self {
        Self { contig, position }
    }
}

/// Base quality score (Phred scale), capped at SAM maximum 93.
/// # Invariants
/// Stored value is in `[0, 93]` after [`BaseQuality::new`].
/// Phred encoding: error probability `10^(-Q/10)`.
/// # Ownership
/// Newtype over `u8`; `Copy` and cheap to pass by value.
/// # Mutation
/// Immutable after construction; create a new value to change quality.
/// # Biological assumptions
/// Per-base sequencing error confidence on the read strand reported in the file.
/// # Java equivalence
/// Same role as Phred qualities in htsjdk `QualityUtils` / SAM `QUAL` field.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
pub struct BaseQuality(u8);

impl BaseQuality {
    pub fn new(quality: u8) -> Self {
        Self(quality.min(93)) // Cap at 93 as per SAM spec
    }

    pub fn value(&self) -> u8 {
        self.0
    }

    pub fn phred_score(&self) -> u8 {
        self.0
    }

    pub fn error_probability(&self) -> f64 {
        10.0f64.powf(-(self.0 as f64) / 10.0)
    }
}

/// SAM/BAM mapping quality, distinguishing unavailable (`255`) from a real score.
/// # Invariants
/// [`MappingQuality::Unavailable`] is exactly SAM MAPQ `255`.
/// [`MappingQuality::Score`] holds a reported MAPQ (`0..=254`); typical aligners use `0..=60`.
/// # Ownership
/// `Copy` enum; cheap to pass by value.
/// # Mutation
/// Immutable after construction; convert at the BAM boundary via [`Self::from_sam_mapq`].
/// # Biological assumptions
/// MAPQ is Phred-scaled mapping confidence; `255` means “unavailable”, not high quality.
/// # Java equivalence
/// htsjdk `SAMRecord.getMappingQuality` / GATK MQ filters (treat `255` as special).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum MappingQuality {
    Score(u8),
    Unavailable,
}

impl MappingQuality {
    pub const UNAVAILABLE_SAM: u8 = 255;

    #[inline]
    pub const fn from_sam_mapq(mapq: u8) -> Self {
        if mapq == Self::UNAVAILABLE_SAM {
            Self::Unavailable
        } else {
            Self::Score(mapq)
        }
    }

    #[inline]
    pub const fn as_sam_mapq(self) -> u8 {
        match self {
            Self::Score(q) => q,
            Self::Unavailable => Self::UNAVAILABLE_SAM,
        }
    }

    #[inline]
    pub const fn score(self) -> Option<u8> {
        match self {
            Self::Score(q) => Some(q),
            Self::Unavailable => None,
        }
    }
}

impl PartialOrd for MappingQuality {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for MappingQuality {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        // SAM numeric order: Unavailable (255) sorts above any reported score.
        self.as_sam_mapq().cmp(&other.as_sam_mapq())
    }
}

/// Compile-time guarantee: quality newtypes stay as small as their wire layout.
const _: () = {
    use core::mem::size_of;

    assert!(size_of::<BaseQuality>() == size_of::<u8>());
    // Score(u8) | Unavailable: Score covers all u8 values, so no niche — two bytes.
    assert!(size_of::<MappingQuality>() == 2);
};

/// DNA base (IUPAC subset: A/C/G/T/N).
/// # Invariants
/// Only enumerated nucleotides; unknown input maps to `None` via [`Base::from_char`].
/// # Ownership
/// `Copy` enum; no allocation.
/// # Mutation
/// Immutable; complement via [`Base::complement`].
/// # Biological assumptions
/// Uppercase canonical bases; `N` marks ambiguous/unknown base calls.
/// # Java equivalence
/// Similar to htsjdk `SequenceUtil` base constants / GATK `BaseUtils` (conceptual).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum Base {
    A,
    C,
    G,
    T,
    N, // Unknown
}

impl Base {
    pub fn from_char(c: char) -> Option<Self> {
        match c.to_ascii_uppercase() {
            'A' => Some(Base::A),
            'C' => Some(Base::C),
            'G' => Some(Base::G),
            'T' => Some(Base::T),
            'N' => Some(Base::N),
            _ => None,
        }
    }

    pub fn to_char(self) -> char {
        match self {
            Base::A => 'A',
            Base::C => 'C',
            Base::G => 'G',
            Base::T => 'T',
            Base::N => 'N',
        }
    }

    pub fn complement(self) -> Self {
        match self {
            Base::A => Base::T,
            Base::C => Base::G,
            Base::G => Base::C,
            Base::T => Base::A,
            Base::N => Base::N,
        }
    }
}

/// Allele sequence at a variant locus (reference or alternate).
/// # Invariants
/// `bases` holds the literal allele sequence; empty alleles are allowed but may be invalid for callers.
/// Equality is sequence-wise over [`Base`] values.
/// # Ownership
/// Owns `Vec<Base>`; clone to share independently; borrow via `bases` slice.
/// # Mutation
/// Public `bases` field is mutable; prefer builders for validated alleles.
/// # Biological assumptions
/// Left-aligned allele representation as stored in VCF `REF`/`ALT`; no normalization enforced here.
/// # Java equivalence
/// Approximates htsjdk `Allele` / GATK `Allele` wrapper around byte sequence.
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct Allele {
    pub bases: Vec<Base>,
}

impl Allele {
    pub fn new(bases: Vec<Base>) -> Self {
        Self { bases }
    }

    pub fn from_string(s: &str) -> Option<Self> {
        let bases: Option<Vec<Base>> = s.chars().map(Base::from_char).collect();
        bases.map(Allele::new)
    }

    pub fn length(&self) -> usize {
        self.bases.len()
    }

    pub fn is_empty(&self) -> bool {
        self.bases.is_empty()
    }
}

impl std::fmt::Display for Allele {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let text: String = self.bases.iter().map(|b| b.to_char()).collect();
        f.write_str(&text)
    }
}

/// High-level variant classification derived from allele lengths.
/// # Invariants
/// Classification is heuristic from ref/alt lengths in [`VariantContext::variant_type`]; not VCF `SVTYPE`.
/// # Ownership
/// `Copy`-less small enum; cheap to clone.
/// # Mutation
/// N/A (enum value).
/// # Biological assumptions
/// SNP = single-base substitution; indels/complex from length differences among alleles.
/// # Java equivalence
/// Similar intent to GATK `VariantContext` type inference (not a direct Java enum).
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub enum VariantType {
    SNP,
    Insertion,
    Deletion,
    Complex, // Both insertion and deletion
    Mixed,   // Multiple alternative alleles
}

/// Sample genotype as allele indices and declared ploidy.
/// # Invariants
/// `alleles` entries index into a variant's allele list (0 = reference).
/// `ploidy` should match `alleles.len` for well-formed VCF genotypes (not enforced).
/// # Ownership
/// Owns `Vec<usize>`; clone to duplicate sample genotype state.
/// # Mutation
/// Public fields; mutable when assembling [`VariantContext`] sample maps.
/// # Biological assumptions
/// Unphased allele indices unless callers encode phasing externally.
/// # Java equivalence
/// Approximates htsjdk `Genotype` / GATK `Genotype` index list.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Genotype {
    pub alleles: Vec<usize>, // Indices into the alleles array
    pub ploidy: u8,
}

impl Genotype {
    pub fn new(alleles: Vec<usize>, ploidy: u8) -> Self {
        Self { alleles, ploidy }
    }

    pub fn is_hom_ref(&self) -> bool {
        self.alleles.iter().all(|&a| a == 0)
    }

    pub fn is_hom_var(&self) -> bool {
        self.alleles.iter().all(|&a| a > 0) && self.alleles.iter().all(|&a| a == self.alleles[0])
    }

    pub fn is_het(&self) -> bool {
        !self.is_hom_ref() && !self.is_hom_var()
    }
}

/// Variant at a genomic locus with alleles, per-sample genotypes, and INFO-like attributes.
/// # Invariants
/// `position` anchors the variant; `reference` is allele index 0 conceptually.
/// `id` is a VCF-style identifier (`"."` when unset).
/// # Ownership
/// Owns allele vectors, genotype map, and attribute map; clone for snapshots.
/// # Mutation
/// Mutate via [`VariantContext::add_genotype`] and [`VariantContext::set_attribute`]; fields also public.
/// # Biological assumptions
/// One variant record per locus; multi-allelic sites use `alternate_alleles`.
/// # Java equivalence
/// Approximates GATK `VariantContext` / htsjdk `VariantContext`.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct VariantContext {
    pub id: String,
    pub position: GenomicPosition,
    pub reference: Allele,
    pub alternate_alleles: Vec<Allele>,
    /// Sample → genotype in insertion order (stable for report/emit iteration).
    pub genotypes: IndexMap<String, Genotype>,
    /// Attribute key → value in insertion order (stable for INFO-like dumps).
    pub attributes: IndexMap<String, String>,
}

impl VariantContext {
    pub fn new(
        position: GenomicPosition,
        reference: Allele,
        alternate_alleles: Vec<Allele>,
    ) -> Self {
        Self {
            id: ".".to_string(),
            position,
            reference,
            alternate_alleles,
            genotypes: IndexMap::new(),
            attributes: IndexMap::new(),
        }
    }

    pub fn variant_type(&self) -> VariantType {
        if self.alternate_alleles.is_empty() {
            return VariantType::SNP; // Default
        }

        let ref_len = self.reference.length();
        let alt_lengths: Vec<usize> = self.alternate_alleles.iter().map(|a| a.length()).collect();

        if alt_lengths.iter().all(|&len| len == ref_len) && ref_len == 1 {
            VariantType::SNP
        } else if alt_lengths.iter().all(|&len| len > ref_len) {
            VariantType::Insertion
        } else if alt_lengths.iter().all(|&len| len < ref_len) {
            VariantType::Deletion
        } else {
            VariantType::Complex
        }
    }

    pub fn add_genotype(&mut self, sample: String, genotype: Genotype) {
        self.genotypes.insert(sample, genotype);
    }

    pub fn set_attribute(&mut self, key: String, value: String) {
        self.attributes.insert(key, value);
    }

    pub fn alternates(&self) -> &[Allele] {
        &self.alternate_alleles
    }

    pub fn contig(&self) -> u32 {
        self.position.contig
    }

    pub fn chromosome(&self) -> String {
        format!("chr{}", self.position.contig)
    }
}

/// Per-read mapping quality and base-quality vector.
/// # Invariants
/// `base_qualities.len` should match read sequence length when attached to [`SequenceRead`] (caller's duty).
/// Mapping quality follows SAM MAPQ semantics `[0, 255]`.
/// # Ownership
/// Owns `Vec<BaseQuality>`; clone to duplicate quality tracks.
/// # Mutation
/// Public fields; typically set once when parsing SAM/BAM/FASTQ.
/// # Biological assumptions
/// Qualities apply to sequenced bases on the stored read orientation.
/// # Java equivalence
/// SAM MAPQ + QUAL fields via htsjdk `SAMRecord`.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ReadQuality {
    pub mapping_quality: MappingQuality,
    pub base_qualities: Vec<BaseQuality>,
}

impl ReadQuality {
    pub fn new(mapping_quality: MappingQuality, base_qualities: Vec<BaseQuality>) -> Self {
        Self {
            mapping_quality,
            base_qualities,
        }
    }

    /// Construct from a raw SAM MAPQ byte and uncapped base qualities.
    pub fn from_sam(mapq: u8, base_qualities: Vec<u8>) -> Self {
        Self {
            mapping_quality: MappingQuality::from_sam_mapq(mapq),
            base_qualities: base_qualities.into_iter().map(BaseQuality::new).collect(),
        }
    }

    pub fn from_vec(base_qualities: Vec<u8>) -> Self {
        Self {
            mapping_quality: MappingQuality::Score(0),
            base_qualities: base_qualities.into_iter().map(BaseQuality::new).collect(),
        }
    }
}

/// In-memory sequence read with alignment anchor and pairing flags.
/// # Invariants
/// `sequence` and `qualities.base_qualities` should have equal length when populated from alignments.
/// `position` is the 1-based leftmost mapped position on `contig`.
/// # Ownership
/// Owns `id`, `sequence`, and nested [`ReadQuality`]; clone for parallel pipelines.
/// # Mutation
/// Public fields; readers/writers may mutate when normalizing records.
/// # Biological assumptions
/// Single fragment; paired-end metadata is flag-only (`is_paired`), not mate coordinates.
/// # Java equivalence
/// Approximates htsjdk `SAMRecord` core fields (name, seq, qual, pos, flags).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SequenceRead {
    pub id: String,
    pub sequence: Vec<Base>,
    pub qualities: ReadQuality,
    pub position: GenomicPosition,
    pub is_reverse_strand: bool,
    pub is_paired: bool,
}

impl SequenceRead {
    pub fn new(
        id: String,
        sequence: Vec<Base>,
        qualities: ReadQuality,
        position: GenomicPosition,
        is_reverse_strand: bool,
        is_paired: bool,
    ) -> Self {
        Self {
            id,
            sequence,
            qualities,
            position,
            is_reverse_strand,
            is_paired,
        }
    }

    pub fn length(&self) -> usize {
        self.sequence.len()
    }

    pub fn len(&self) -> usize {
        self.length()
    }

    pub fn is_empty(&self) -> bool {
        self.sequence.is_empty()
    }

    pub fn base_at(&self, index: usize) -> Option<Base> {
        self.sequence.get(index).copied()
    }

    pub fn quality_at(&self, index: usize) -> Option<BaseQuality> {
        self.qualities.base_qualities.get(index).copied()
    }

    pub fn is_valid_dna(&self) -> bool {
        self.sequence
            .iter()
            .all(|&base| matches!(base, Base::A | Base::C | Base::G | Base::T))
    }

    pub fn name(&self) -> &str {
        &self.id
    }

    pub fn average_quality(&self) -> f64 {
        if self.qualities.base_qualities.is_empty() {
            return 0.0;
        }
        let sum: u64 = self
            .qualities
            .base_qualities
            .iter()
            .map(|q| u64::from(q.value()))
            .sum();
        sum as f64 / self.qualities.base_qualities.len() as f64
    }

    pub fn has_minimum_quality(&self, min_quality: u8) -> bool {
        self.qualities
            .base_qualities
            .iter()
            .all(|q| q.value() >= min_quality)
    }

    pub fn filter_bases_by_quality(&self, min_quality: u8) -> Vec<(u8, u8)> {
        self.sequence
            .iter()
            .zip(self.qualities.base_qualities.iter())
            .filter_map(|(base, quality)| {
                if quality.value() >= min_quality {
                    Some((base.to_char() as u8, quality.value()))
                } else {
                    None
                }
            })
            .collect()
    }
}
