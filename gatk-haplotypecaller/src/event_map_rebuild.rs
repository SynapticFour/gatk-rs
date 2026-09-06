//! Pure EventMap rebuild helpers (L7-B2).
//! Keeps genotyping EventMap construction as data transforms over haplotypes + prior events,
//! rather than only mutating [`crate::assembly_result_set::AssemblyResultSet`] in place.

use crate::event_map::{collect_variation_events, VariationEvent};
use crate::haplotype::Haplotype;

/// Options for [`rebuild_variation_events`].
/// # Invariants
/// `event_map_only` ignores prior/supplement lists (CIGAR EventMap only).
/// When supplements are enabled, `merge_read_supplements` selects supplement vs prior source.
/// # Ownership
/// [`Copy`] option flags.
/// # Mutation
/// Immutable per rebuild call.
/// # Biological assumptions
/// Controls whether read-discovered alleles merge into haplotype EventMap sites.
/// # Java equivalence
/// Rust-native L7-B2 rebuild knobs (algorithm parity with Java EventMap regeneration).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RebuildVariationEventsOpts {
    /// When true, ignore prior/supplement lists (CIGAR/EventMap only).
    pub event_map_only: bool,
    /// Merge `preserved_supplement` instead of `prior_events` when supplements are enabled.
    pub merge_read_supplements: bool,
}

/// Rebuild variation events from haplotype CIGARs, optionally merging prior/supplement alleles.
pub fn rebuild_variation_events(
    haplotypes: &[Haplotype],
    full_ref: &[u8],
    full_pad_1based: u64,
    contig: &str,
    max_mnp_distance: usize,
    prior_events: &[VariationEvent],
    preserved_supplement: &[VariationEvent],
    opts: RebuildVariationEventsOpts,
) -> Vec<VariationEvent> {
    let mut events = collect_variation_events(
        haplotypes,
        full_ref,
        full_pad_1based,
        contig,
        max_mnp_distance,
    );
    if !opts.event_map_only {
        let source = if opts.merge_read_supplements {
            preserved_supplement
        } else {
            prior_events
        };
        for e in source {
            if !events.iter().any(|x| {
                x.start_1based == e.start_1based
                    && x.ref_allele == e.ref_allele
                    && x.alt_allele == e.alt_allele
            }) {
                // CLONE: needed because owned element into collection.
                events.push(e.clone());
            }
        }
    }
    // 6R.57: do not drop SNPs colocated with another haplotype's indel after union.
    events.sort_by_key(|e| e.start_1based);
    events.dedup_by(|a, b| {
        a.start_1based == b.start_1based
            && a.ref_allele == b.ref_allele
            && a.alt_allele == b.alt_allele
    });
    if !opts.event_map_only {
        crate::read_event_discovery::scrub_p12_cluster_phantom_alleles(&mut events);
    }
    events
}
