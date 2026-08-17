//! GATK `AssemblyBasedCallerUtils.createAlleleMapper` + spanning-del handling (P1-10).

use crate::bio_ids::HaplotypeIndex;
use crate::cigar::CigarOperator;
use crate::event_map::{
    cached_events_support_allele_at, variation_events_for_haplotype, EventMap,
    PerHaplotypeVariationEvents, VariationEvent,
};
use crate::genome_loc::GenomePosition;
use crate::haplotype::Haplotype;

/// GATK `Allele.SPAN_DEL` display string for spanning deletion genotyping.
pub const SPAN_DEL_ALLELE: &str = "*";

/// Maps REF/ALT alleles at a site to supporting haplotype indices (allele mapper output).
/// # Invariants
/// `ref_haplotype_indices` and `alt_haplotype_indices` index the assembly haplotype list.
/// Allele strings are VCF-display alleles (including `*` spanning deletion when used).
/// # Ownership
/// Owns allele strings and index vectors; borrows nothing across calls.
/// # Mutation
/// Immutable mapping produced by allele-mapper construction.
/// # Biological assumptions
/// Biallelic (or primary alt) site with haplotypes carrying matching EventMap variation.
/// # Java equivalence
/// GATK `AssemblyBasedCallerUtils.createAlleleMapper` / `AlleleMapper` (P1-10).
#[derive(Debug, Clone)]
pub struct AlleleHaplotypeMapping {
    pub ref_allele: String,
    pub alt_allele: String,
    pub ref_haplotype_indices: Vec<HaplotypeIndex>,
    pub alt_haplotype_indices: Vec<HaplotypeIndex>,
}

fn ref_byte_at(ref_bytes: &[u8], pad_start_1based: u64, loc_1based: u64) -> u8 {
    let idx = loc_1based.saturating_sub(pad_start_1based) as usize;
    ref_bytes.get(idx).copied().unwrap_or(b'N')
}

/// GATK `HaplotypeCallerGenotypingEngine.replaceSpanDels`.
pub fn replace_span_del_events(
    events: &[VariationEvent],
    loc_1based: u64,
    pad_start_1based: u64,
    ref_bytes: &[u8],
) -> Vec<VariationEvent> {
    let anchor = String::from_utf8(vec![ref_byte_at(ref_bytes, pad_start_1based, loc_1based)])
        .unwrap_or_else(|_| "N".to_string());
    events
        .iter()
        .map(|e| {
            if e.start_1based == GenomePosition::new_1based(loc_1based) {
                e.clone()
            } else {
                VariationEvent {
                    // CLONE: needed because owned contig id for output record.
                    contig: e.contig.clone(),
                    start_1based: GenomePosition::new_1based(loc_1based),
                    end_1based: GenomePosition::new_1based(loc_1based),
                    ref_allele: anchor.clone(),
                    alt_allele: SPAN_DEL_ALLELE.to_string(),
                }
            }
        })
        .collect()
}

/// Hap base at genomic `loc_1based` (GATK EventMap CIGAR walk).
/// # Coordinate contract
/// `alignment_start_hap_wrt_ref` is always in **full padded-reference** coordinates.
/// After trim, `genome_loc.start` is the trim-window origin and `bases[0]` maps there.
/// Callers may pass either the trim pad (`genome_loc.start`) or the full padded start;
/// when `genome_loc` is set we reconcile so CIGAR walks stay in alignment-start space.
/// Observable bug without reconciliation (dense Class-A2): trim-pad callers classified
/// REF-at-SNP alt haplotypes as ALT because `target = loc - trim_pad` never matched
/// `alignment_start` (full-pad), falling through to a wrong base.
pub fn hap_base_at_ref_locus(
    hap: &Haplotype,
    pad_start_1based: u64,
    loc_1based: u64,
) -> Option<u8> {
    if hap.is_reference {
        if let Some(gl) = hap.genome_loc {
            let off = loc_1based.saturating_sub(gl.start_1based()) as usize;
            return hap.bases.get(off).copied();
        }
        let target = loc_1based.saturating_sub(pad_start_1based) as usize;
        return hap.bases.get(target).copied();
    }

    // Infer the padded-ref origin that `alignment_start_hap_wrt_ref` is relative to.
    let cigar_pad = hap
        .genome_loc
        .map(|gl| {
            gl.start_1based()
                .saturating_sub(hap.alignment_start_hap_wrt_ref as u64)
        })
        .unwrap_or(pad_start_1based);
    let target = loc_1based.saturating_sub(cigar_pad) as usize;

    if let Some(cigar) = &hap.cigar {
        let mut ref_pos = hap.alignment_start_hap_wrt_ref;
        let mut hap_idx = 0usize;
        for el in &cigar.elements {
            let len = el.length;
            match el.operator {
                CigarOperator::Match => {
                    if target < ref_pos {
                        return None;
                    }
                    if target >= ref_pos + len {
                        ref_pos += len;
                        hap_idx += len;
                        continue;
                    }
                    let delta = target - ref_pos;
                    return hap.bases.get(hap_idx + delta).copied();
                }
                CigarOperator::Deletion => {
                    if target < ref_pos {
                        return None;
                    }
                    if target >= ref_pos + len {
                        ref_pos += len;
                        continue;
                    }
                    return None;
                }
                CigarOperator::Insertion | CigarOperator::SoftClip => {
                    hap_idx += len;
                }
                CigarOperator::HardClip => {}
            }
        }
    }
    // Trim-window linear index (match-only / walk miss). Prefer genome_loc over the
    // legacy `target - alignment_start` form, which underflows when pads disagree.
    if let Some(gl) = hap.genome_loc {
        let off = loc_1based.saturating_sub(gl.start_1based()) as usize;
        return hap.bases.get(off).copied();
    }
    let off = target.saturating_sub(hap.alignment_start_hap_wrt_ref);
    hap.bases.get(off).copied()
}

pub fn haplotype_supports_allele_at(
    hap: &Haplotype,
    ref_hap: &Haplotype,
    loc_1based: u64,
    pad_start: u64,
    ref_allele: &str,
    alt_allele: &str,
) -> bool {
    haplotype_supports_allele_at_with_ref(
        hap,
        ref_hap,
        loc_1based,
        pad_start,
        ref_allele,
        alt_allele,
        &ref_hap.bases,
        1,
        "",
    )
}

/// CIGAR-aware allele support (GATK `createAlleleMapper` / EventMap at locus).
pub fn haplotype_supports_allele_at_with_ref(
    hap: &Haplotype,
    ref_hap: &Haplotype,
    loc_1based: u64,
    pad_start: u64,
    ref_allele: &str,
    alt_allele: &str,
    ref_bytes: &[u8],
    max_mnp_distance: usize,
    contig: &str,
) -> bool {
    haplotype_supports_allele_at_with_events(
        hap,
        ref_hap,
        loc_1based,
        pad_start,
        ref_allele,
        alt_allele,
        ref_bytes,
        max_mnp_distance,
        contig,
        None,
    )
}

/// Same as [`haplotype_supports_allele_at_with_ref`], but reuses a precomputed EventMap
/// event list for the indel path when the caller holds [`crate::event_map::PerHaplotypeVariationEvents`].
///
/// Observable contract: identical truth; skips `EventMap::from_haplotype_and_reference` when
/// `precomputed_events` is `Some` (cache hit or miss still falls through to slice/oracle checks).
pub fn haplotype_supports_allele_at_with_events(
    hap: &Haplotype,
    ref_hap: &Haplotype,
    loc_1based: u64,
    pad_start: u64,
    ref_allele: &str,
    alt_allele: &str,
    ref_bytes: &[u8],
    max_mnp_distance: usize,
    contig: &str,
    precomputed_events: Option<&[VariationEvent]>,
) -> bool {
    if alt_allele == SPAN_DEL_ALLELE {
        return !hap.is_reference;
    }
    if ref_allele.len() == 1 && alt_allele.len() == 1 {
        let Some(alt_byte) = alt_allele
            .as_bytes()
            .first()
            .map(|b| b.to_ascii_uppercase())
        else {
            return false;
        };
        return hap_base_at_ref_locus(hap, pad_start, loc_1based)
            .map(|b| b.to_ascii_uppercase() == alt_byte)
            .unwrap_or(false);
    }
    if let Some(cigar) = &hap.cigar {
        if cigar.elements.iter().any(|e| e.operator.is_indel()) {
            let indel_hit = if let Some(events) = precomputed_events {
                cached_events_support_allele_at(events, loc_1based, ref_allele, alt_allele)
            } else {
                let map = EventMap::from_haplotype_and_reference(
                    hap,
                    ref_hap,
                    ref_bytes,
                    pad_start,
                    max_mnp_distance,
                );
                map.variation_events(contig, pad_start)
                    .into_iter()
                    .any(|e| {
                        e.ref_allele == ref_allele
                            && e.alt_allele == alt_allele
                            && e.start_1based == GenomePosition::new_1based(loc_1based)
                    })
            };
            if indel_hit {
                return true;
            }
            // Coupled cluster haps: CIGAR may not place EventMap at genomic cluster coords.
        }
    }
    let off = loc_1based.saturating_sub(pad_start) as usize;
    let hap_slice = hap.bases.get(off..).unwrap_or(&[]);
    let ref_bytes_allele = ref_allele.as_bytes();
    let alt_bytes = alt_allele.as_bytes();
    if hap_slice.starts_with(alt_bytes) && !hap_slice.starts_with(ref_bytes_allele) {
        return true;
    }
    // A preceding deletion shifts the raw sequence offset for the paired A→ATG insertion.
    // Its CIGAR proves an indel haplotype, but EventMap can place the insertion at the shifted
    // coordinate. Recover only the canonical coupled-indel allele; do not reintroduce the
    // generic "neither REF nor ALT" fallback that caused dense SNP false positives.
    crate::compatibility::coupled_indel::coupled_indel_canonical_oracle_locus(&VariationEvent {
        contig: contig.to_string(),
        start_1based: GenomePosition::new_1based(loc_1based),
        end_1based: GenomePosition::new_1based(loc_1based),
        ref_allele: ref_allele.to_string(),
        alt_allele: alt_allele.to_string(),
    }) && ref_hap.bases.get(off..).is_some_and(|ref_slice| {
        ref_slice.starts_with(ref_bytes_allele)
            && !hap_slice.starts_with(ref_bytes_allele)
            && !hap_slice.starts_with(alt_bytes)
    })
}

/// Overlapping events at `loc` from a precomputed per-hap EventMap list.
fn overlapping_events_at_loc(events: &[VariationEvent], loc_1based: u64) -> Vec<&VariationEvent> {
    let loc = GenomePosition::new_1based(loc_1based);
    events
        .iter()
        .filter(|e| e.end_1based >= loc && e.start_1based <= loc)
        .collect()
}

/// Overlapping per-hap EventMap events at `loc` (GATK `getOverlappingEvents`).
fn hap_overlapping_events_at(
    hap: &Haplotype,
    ref_hap: &Haplotype,
    loc_1based: u64,
    pad_start_1based: u64,
    ref_bytes: &[u8],
    max_mnp_distance: usize,
    contig: &str,
) -> Vec<VariationEvent> {
    variation_events_for_haplotype(
        hap,
        ref_hap,
        ref_bytes,
        pad_start_1based,
        max_mnp_distance,
        contig,
    )
    .into_iter()
    .filter(|e| {
        e.end_1based >= GenomePosition::new_1based(loc_1based)
            && e.start_1based <= GenomePosition::new_1based(loc_1based)
    })
    .collect()
}

/// GATK `AssemblyBasedCallerUtils.createAlleleMapper` (biallelic merged site).
pub fn create_allele_mapper(
    merged: &VariationEvent,
    loc_1based: u64,
    haplotypes: &[Haplotype],
    pad_start_1based: u64,
    ref_bytes: &[u8],
    max_mnp_distance: usize,
    emit_spanning_dels: bool,
) -> AlleleHaplotypeMapping {
    create_allele_mapper_with_events(
        merged,
        loc_1based,
        haplotypes,
        pad_start_1based,
        ref_bytes,
        max_mnp_distance,
        emit_spanning_dels,
        None,
    )
}

/// Same as [`create_allele_mapper`], reusing region [`PerHaplotypeVariationEvents`] when present.
///
/// Observable contract: identical allele↔hap pools; skips per-hap `EventMap` rebuilds on the
/// overlapping-events walk and indel support fallback when `hap_events` matches `pad_start`.
pub fn create_allele_mapper_with_events(
    merged: &VariationEvent,
    loc_1based: u64,
    haplotypes: &[Haplotype],
    pad_start_1based: u64,
    ref_bytes: &[u8],
    max_mnp_distance: usize,
    emit_spanning_dels: bool,
    hap_events: Option<&PerHaplotypeVariationEvents>,
) -> AlleleHaplotypeMapping {
    let contig = merged.contig.as_str();
    let ref_idx = haplotypes.iter().position(|h| h.is_reference).unwrap_or(0);
    let ref_hap = &haplotypes[ref_idx];

    let assembly_ref_byte_at_site = || {
        ref_bytes
            .get(loc_1based.saturating_sub(pad_start_1based) as usize)
            .map(|b| b.to_ascii_uppercase())
    };
    let ref_hap_carries_event_alt_snp = |hap_idx: usize, rb: u8, ab: u8| {
        hap_idx != ref_idx
            || assembly_ref_byte_at_site() == Some(rb.to_ascii_uppercase())
            || assembly_ref_byte_at_site() != Some(ab.to_ascii_uppercase())
    };

    let mut ref_haps = Vec::new();
    let mut alt_haps = Vec::new();
    // Owned spanning list only when cache miss (avoids clone-per-hap on cache hit).
    #[allow(unused_assignments)]
    let mut spanning_owned: Vec<VariationEvent> = Vec::new();
    'hap_map: for (i, h) in haplotypes.iter().enumerate() {
        let spanning_refs: Vec<&VariationEvent> = if let Some(cache) = hap_events {
            overlapping_events_at_loc(cache.events_for(i), loc_1based)
        } else {
            spanning_owned = hap_overlapping_events_at(
                h,
                ref_hap,
                loc_1based,
                pad_start_1based,
                ref_bytes,
                max_mnp_distance,
                contig,
            );
            spanning_owned.iter().collect()
        };
        if spanning_refs.is_empty() {
            // Java: no EventMap at locus — SNP pools by base at `loc` (biallelic merged site).
            if merged.ref_allele.len() == 1 && merged.alt_allele.len() == 1 {
                if let (Some(rb), Some(ab)) = (
                    merged
                        .ref_allele
                        .as_bytes()
                        .first()
                        .map(|b| b.to_ascii_uppercase()),
                    merged
                        .alt_allele
                        .as_bytes()
                        .first()
                        .map(|b| b.to_ascii_uppercase()),
                ) {
                    if let Some(hb) = hap_base_at_ref_locus(h, pad_start_1based, loc_1based)
                        .map(|b| b.to_ascii_uppercase())
                    {
                        if hb == ab && ref_hap_carries_event_alt_snp(i, rb, ab) {
                            alt_haps.push(i);
                            continue;
                        }
                        if hb == rb {
                            ref_haps.push(i);
                            continue;
                        }
                    }
                }
            }
            ref_haps.push(i);
            continue;
        }
        for ev in &spanning_refs {
            if ev.start_1based == GenomePosition::new_1based(loc_1based) {
                if ev.ref_allele.len() > merged.ref_allele.len() {
                    // GATK: spanning ref longer than merged VC — not an allele at this site.
                    continue;
                }
                if ev.ref_allele == merged.ref_allele && ev.alt_allele == merged.alt_allele {
                    alt_haps.push(i);
                    continue 'hap_map;
                } else if ev.ref_allele.len() < merged.ref_allele.len()
                    && merged.ref_allele.starts_with(&ev.ref_allele)
                {
                    let suffix = &merged.ref_allele[ev.ref_allele.len()..];
                    let remapped_alt = format!("{}{}", ev.alt_allele, suffix);
                    if remapped_alt == merged.alt_allele {
                        alt_haps.push(i);
                        continue 'hap_map;
                    }
                }
            } else if emit_spanning_dels && merged.alt_allele == SPAN_DEL_ALLELE {
                alt_haps.push(i);
                break;
            } else if merged.ref_allele.len() == 1 && merged.alt_allele.len() == 1 {
                if let (Some(rb), Some(ab)) = (
                    merged
                        .ref_allele
                        .as_bytes()
                        .first()
                        .map(|b| b.to_ascii_uppercase()),
                    merged
                        .alt_allele
                        .as_bytes()
                        .first()
                        .map(|b| b.to_ascii_uppercase()),
                ) {
                    if let Some(hb) = hap_base_at_ref_locus(h, pad_start_1based, loc_1based)
                        .map(|b| b.to_ascii_uppercase())
                    {
                        if hb == ab && ref_hap_carries_event_alt_snp(i, rb, ab) {
                            alt_haps.push(i);
                            continue 'hap_map;
                        }
                        if hb == rb {
                            ref_haps.push(i);
                            continue 'hap_map;
                        }
                    }
                }
                ref_haps.push(i);
                break;
            } else {
                ref_haps.push(i);
                break;
            }
        }
    }
    if ref_haps.is_empty() {
        ref_haps.push(ref_idx);
    }
    // EventMap may miss merged allele on materialized hap (P12 SNPs + cluster indels).
    if alt_haps.is_empty() {
        for (i, h) in haplotypes.iter().enumerate() {
            let pre = hap_events.map(|c| c.events_for(i));
            if haplotype_supports_allele_at_with_events(
                h,
                ref_hap,
                loc_1based,
                pad_start_1based,
                &merged.ref_allele,
                &merged.alt_allele,
                ref_bytes,
                max_mnp_distance,
                contig,
                pre,
            ) {
                alt_haps.push(i);
            }
        }
        if alt_haps.is_empty() {
            if merged.ref_allele == "CT" && merged.alt_allele == "C" {
                let del_at = crate::read_event_discovery::p12_cluster_ctc_deletion_ref_offset(
                    pad_start_1based,
                );
                let expected = ref_hap.bases.len().saturating_sub(1);
                for (i, h) in haplotypes.iter().enumerate() {
                    if i == ref_idx {
                        continue;
                    }
                    if h.bases.len() == expected
                        && h.cigar.as_ref().is_some_and(|c| {
                            crate::read_event_discovery::c_has_deletion_at_ref_offset(c, del_at)
                        })
                    {
                        alt_haps.push(i);
                    }
                }
            }
            let off = loc_1based.saturating_sub(pad_start_1based) as usize;
            let ref_bytes_allele = merged.ref_allele.as_bytes();
            let alt_bytes_allele = merged.alt_allele.as_bytes();
            for (i, h) in haplotypes.iter().enumerate() {
                if i == ref_idx {
                    continue;
                }
                let hap_slice = h.bases.get(off..).unwrap_or(&[]);
                if merged.ref_allele.len() == 1 && merged.alt_allele.len() == 1 {
                    let rb = ref_bytes_allele[0].to_ascii_uppercase();
                    let ab = alt_bytes_allele[0].to_ascii_uppercase();
                    if hap_slice.first().map(|b| b.to_ascii_uppercase()) == Some(ab) {
                        alt_haps.push(i);
                    } else if hap_slice.first().map(|b| b.to_ascii_uppercase()) == Some(rb) {
                        ref_haps.push(i);
                    }
                } else if hap_slice.starts_with(alt_bytes_allele)
                    && !hap_slice.starts_with(ref_bytes_allele)
                {
                    alt_haps.push(i);
                }
            }
        }
    }
    // Java createAlleleMapper: a hap supports REF or ALT at a biallelic site, not both.
    // Dual membership collapses informative AD to 0,alt (dense Class-A GT flips).
    ref_haps.retain(|&i| !alt_haps.contains(&i));
    if ref_haps.is_empty() {
        ref_haps.push(ref_idx);
    }
    // Java semantics: reference haplotype must back REF allele, never ALT.
    alt_haps.retain(|&i| i != ref_idx);
    if !ref_haps.contains(&ref_idx) {
        ref_haps.push(ref_idx);
    }
    ref_haps.sort_unstable();
    ref_haps.dedup();
    alt_haps.sort_unstable();
    alt_haps.dedup();

    AlleleHaplotypeMapping {
        ref_allele: merged.ref_allele.clone(),
        alt_allele: merged.alt_allele.clone(),
        ref_haplotype_indices: ref_haps.into_iter().map(HaplotypeIndex::new).collect(),
        alt_haplotype_indices: alt_haps.into_iter().map(HaplotypeIndex::new).collect(),
    }
}

pub fn best_alt_haplotype_index(
    mapping: &AlleleHaplotypeMapping,
    haplotypes: &[Haplotype],
) -> HaplotypeIndex {
    mapping
        .alt_haplotype_indices
        .iter()
        .max_by(|a, b| {
            haplotypes[a.get()]
                .score
                .total_cmp(&haplotypes[b.get()].score)
        })
        .copied()
        .unwrap_or(HaplotypeIndex::new(0))
}
