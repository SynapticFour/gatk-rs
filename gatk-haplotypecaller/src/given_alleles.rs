//! GATK `AssemblyBasedCallerUtils.addGivenAlleles` — forced alleles into assembly.
//! **CLI:** `-alleles` is not exposed on `gatk-rs HaplotypeCaller` (Sprint G). Library wiring
//! exists for L2 `c5-force` via internal `GatkConfig` parameters only.
//! See `docs/ARCHITECTURE.md` (T3-5).

use crate::alignment::SwParameters;
use crate::assembly_region_trimmer::TrimVariant;
use crate::assembly_result_set::AssemblyResultSet;
use crate::event_map::{
    collect_variation_events, prefer_indel_over_colocated_snps, VariationEvent,
};
use crate::genome_loc::GenomePosition;
use crate::haplotype::Haplotype;
use crate::haplotype_cigar::calculate_haplotype_cigar_with_strategy;
use crate::read_event_discovery::{
    apply_event_to_ref, events_match, tag_alt_haplotype_from_reference,
};
use crate::smith_waterman::SwOverhangStrategy;
use gatk_common::GatkResult;

/// One forced allele site (`-alleles` / `givenAlleles` VCF input).
/// # Invariants
/// Coordinates are 1-based inclusive; `alt_alleles` are non-empty VCF ALT strings.
/// # Ownership
/// Owns contig and allele strings; merged into assembly/trim event lists.
/// # Mutation
/// Immutable site description after parse/load.
/// # Biological assumptions
/// Forced alleles that must appear in assembly / genotyping even without graph discovery.
/// # Java equivalence
/// GATK `-alleles` / `givenAlleles` input to `AssemblyBasedCallerUtils.addGivenAlleles`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GatkGivenAllele {
    pub contig: String,
    pub start_1based: u64,
    pub end_1based: u64,
    pub ref_allele: String,
    pub alt_alleles: Vec<String>,
}

fn given_to_variation_event(given: &GatkGivenAllele, alt: &str) -> VariationEvent {
    VariationEvent {
        // CLONE: needed because owned contig id for output record.
        contig: given.contig.clone(),
        start_1based: GenomePosition::new_1based(given.start_1based),
        end_1based: GenomePosition::new_1based(given.end_1based.max(given.start_1based)),
        ref_allele: given.ref_allele.clone(),
        alt_allele: alt.to_string(),
    }
}

/// Java: merge `-alleles` into `allVariationEvents` before `trimmer.trim` (not `addGivenAlleles` yet).
pub fn given_alleles_to_trim_variants(
    given: &[GatkGivenAllele],
    contig: &str,
    trim_variants: &mut Vec<TrimVariant>,
) {
    for site in given {
        if site.contig != contig {
            continue;
        }
        for alt in &site.alt_alleles {
            if alt == &site.ref_allele {
                continue;
            }
            let is_indel = site.ref_allele.len() != alt.len();
            let end = site.end_1based.max(site.start_1based);
            if trim_variants
                .iter()
                .any(|t| t.start == site.start_1based && t.end == end)
            {
                continue;
            }
            trim_variants.push(TrimVariant {
                // CLONE: needed because owned contig id for output record.
                contig: site.contig.clone(),
                start: site.start_1based,
                end,
                is_indel,
            });
        }
    }
}

/// Merge forced alleles into the assembly result (Java `addGivenAlleles`).
pub fn merge_given_alleles_into_assembly(
    given: &[GatkGivenAllele],
    assembly: &mut AssemblyResultSet,
) -> GatkResult<()> {
    if given.is_empty() {
        return Ok(());
    }
    let sw = SwParameters::gatk_haplotype_to_reference();
    let pad_start = assembly.padded_reference_start_1based();
    let contig = assembly.contig.clone();
    let ref_bases = assembly.reference_bases_shared();
    let ref_hap = assembly
        .haplotypes
        .iter()
        .find(|h| h.is_reference)
        .cloned()
        // CLONE: needed because fallback owns pileup/value when Option miss.
        .unwrap_or_else(|| Haplotype::new(ref_bases.as_ref().to_vec(), true));
    let apply_bases = assembly.apply_bases_shared();
    let apply_pad = ref_hap
        .genome_loc
        .map(|g| g.start_1based())
        .unwrap_or(pad_start);

    let asm_events = collect_variation_events(
        &assembly.haplotypes,
        &apply_bases,
        apply_pad,
        &contig,
        assembly.max_mnp_distance(),
    );

    let mut seen: std::collections::HashSet<Vec<u8>> = assembly
        .haplotypes
        .iter()
        .map(|h| h.bases.clone())
        .collect();
    let kmer = assembly.kmer_size_for_dump();
    for site in given {
        for alt in &site.alt_alleles {
            if alt == &site.ref_allele {
                continue;
            }
            let event = given_to_variation_event(site, alt);
            if asm_events.iter().any(|e| events_match(e, &event)) {
                continue;
            }
            let Some(alt_bases) = apply_event_to_ref(&apply_bases, &event, apply_pad) else {
                continue;
            };
            // CLONE: needed because owned HashMap/BTree/HashSet key or value.
            if !seen.insert(alt_bases.clone()) {
                continue;
            }
            let cigar = calculate_haplotype_cigar_with_strategy(
                &apply_bases,
                &alt_bases,
                &sw,
                SwOverhangStrategy::Indel,
            );
            let Some(cigar) = cigar else {
                continue;
            };
            let mut h = Haplotype::new(alt_bases, false);
            tag_alt_haplotype_from_reference(&mut h, &ref_hap, kmer);
            h.cigar = Some(cigar);
            h.score = crate::read_event_discovery::SUPPLEMENT_HAPLOTYPE_SCORE;
            assembly.haplotypes.push(h);
        }
    }

    let mut events = collect_variation_events(
        &assembly.haplotypes,
        &apply_bases,
        apply_pad,
        &contig,
        assembly.max_mnp_distance(),
    );
    for site in given {
        for alt in &site.alt_alleles {
            if alt == &site.ref_allele {
                continue;
            }
            let event = given_to_variation_event(site, alt);
            if !events.iter().any(|e| events_match(e, &event)) {
                events.push(event);
            }
        }
    }
    prefer_indel_over_colocated_snps(&mut events);
    events.sort();
    events.dedup();
    assembly.variation_events = events;
    assembly.variation_present =
        assembly.haplotypes.iter().any(|h| !h.is_reference) && assembly.haplotypes.len() > 1;
    Ok(())
}

/// Load forced alleles from a minimal VCF (POS REF ALT per variant line). N3 CLI hook.
pub fn load_given_alleles_from_vcf_path(
    path: &std::path::Path,
) -> GatkResult<Vec<GatkGivenAllele>> {
    let data = std::fs::read_to_string(path)
        .map_err(|e| gatk_common::GatkError::generic(format!("read given VCF: {e}")))?;
    let mut out = Vec::new();
    for line in data.lines() {
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let cols: Vec<_> = line.split('\t').collect();
        if cols.len() < 5 {
            continue;
        }
        let pos: u64 = cols[1]
            .parse()
            .map_err(|_| gatk_common::GatkError::generic("given VCF POS parse error"))?;
        let ref_a = cols[3].to_string();
        let alt = cols[4].split(',').next().unwrap_or("").to_string();
        if alt.is_empty() || alt == ref_a {
            continue;
        }
        out.push(GatkGivenAllele {
            contig: cols[0].to_string(),
            start_1based: pos,
            end_1based: pos,
            ref_allele: ref_a,
            alt_alleles: vec![alt],
        });
    }
    Ok(out)
}

/// Apply `GATK_RS_HC_GIVEN_VCF` when set (N3 harness only — ignored in production release).
pub fn given_alleles_from_env() -> Vec<GatkGivenAllele> {
    crate::parity_harness::env_string("GATK_RS_HC_GIVEN_VCF")
        .and_then(|p| load_given_alleles_from_vcf_path(std::path::Path::new(&p)).ok())
        .unwrap_or_default()
}
