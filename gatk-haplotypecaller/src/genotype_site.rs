//! Genotype-site orchestration surface (L7-B3).
//! Input bundle for region genotyping; production calls [`assign_genotype_likelihoods_for_region`] directly.

use crate::event_map::VariationEvent;
use crate::genome_loc::GenomePosition;
use crate::haplotype::Haplotype;
use crate::region_read_likelihood::RegionReadLikelihood;
use rust_htslib::bam::Record;

/// Inputs for region genotyping (active-window EventMap walk).
/// # Invariants
/// Likelihood indices align with `likelihood_reads` / `haplotypes`.
/// Active window and pad starts are 1-based; `ref_bytes` covers the padded span.
/// # Ownership
/// Lifetime-bound borrows of likelihoods, BAM records, haplotypes, and events.
/// # Mutation
/// View-only bundle; engines must not mutate borrowed slices.
/// # Biological assumptions
/// One assembly region ready for `assignGenotypeLikelihoods`-style genotyping.
/// # Java equivalence
/// Rust-native input bundle for GATK genotyping engine region walk (L7-B3).
pub struct GenotypeSiteRegion<'a> {
    pub likelihoods: &'a [RegionReadLikelihood],
    pub likelihood_reads: &'a [Record],
    pub pileup_reads: &'a [Record],
    pub supplemental_pileup_reads: Option<&'a [Record]>,
    pub haplotypes: &'a [Haplotype],
    pub ref_bytes: &'a [u8],
    pub pad_start_1based: GenomePosition,
    pub full_reference_bases: &'a [u8],
    pub full_reference_pad_1based: GenomePosition,
    pub active_start_1based: GenomePosition,
    pub active_end_1based: GenomePosition,
    pub contig: &'a str,
    pub max_mnp_distance: usize,
    pub stored_events: &'a [VariationEvent],
    pub graph_events: &'a [VariationEvent],
}
