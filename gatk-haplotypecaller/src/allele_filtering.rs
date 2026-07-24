//! GATK `AlleleFilteringHC.filterAlleles` slice — drop weak haplotypes before genotyping.
//! Java ranks haplotypes using per-allele `getAlleleLikelihoodVsInverse` (HC genotyping GL vs
//! inverse allele). Rust uses that score when PairHMM rows exist, else read-LL sum fallback.

use crate::assembly_result_set::AssemblyResultSet;
use crate::event_map::VariationEvent;
use crate::genome_loc::GenomePosition;
use crate::genotyping::emit_genotype_format_fields;
use crate::haplotype::Haplotype;
use crate::hc_allele_mapping::{
    create_allele_mapper, hap_base_at_ref_locus, haplotype_supports_allele_at_with_ref,
};
use crate::hc_genotyping_engine::{
    biallelic_genotype_log10_likelihoods_gatk, marginalize_rows_to_biallelic_alleles,
    region_likelihoods_to_rows,
};
use crate::java_hc_site_semantics::{
    is_cluster_anchor_snp, is_cluster_ctc_del, is_strict_java_production_emit_candidate,
};
use crate::read_event_discovery::is_p12_phase_e_gap_event;
use crate::region_read_likelihood::RegionReadLikelihood;
use gatk_common::{GatkError, GatkResult};

/// GATK `AlleleFilteringHC` cap on non-reference haplotypes retained for genotyping.
pub const MAX_NON_REF_HAPLOTYPES_FOR_GENOTYPING: usize = 12;
/// Tiny log10 tie-break weight so PairHMM LL sum orders equal HC-inverse scores.
const HAPLOTYPE_RANK_LL_SUM_TIEBREAK_SCALE: f64 = 1e-6;

/// Ensure exactly one haplotype is marked reference (trim/filter can drop the flag).
pub fn ensure_reference_haplotype(haplotypes: &mut [Haplotype]) {
    if haplotypes.is_empty() {
        return;
    }
    if haplotypes.iter().any(|h| h.is_reference) {
        return;
    }
    haplotypes[0].is_reference = true;
}

/// Sum per-read log10 likelihoods per haplotype (Java `AlleleFilteringHC` ranking input).
pub fn haplotype_log10_sums_from_region_likelihoods(
    likelihoods: &[RegionReadLikelihood],
    n_haps: usize,
) -> Vec<f64> {
    let mut sums = vec![0.0_f64; n_haps];
    for rl in likelihoods {
        if rl.haplotype_index.get() >= n_haps {
            continue;
        }
        let ll = rl.log10_likelihood;
        if ll.is_finite() && ll > f64::NEG_INFINITY {
            sums[rl.haplotype_index.get()] += ll;
        }
    }
    sums
}

/// Java `AlleleFilteringHC.getAlleleLikelihoodVsInverse`: `min(PL1-PL0, PL2-PL0)` per event.
fn variation_event_hc_inverse_pl(
    event: &VariationEvent,
    likelihoods: &[RegionReadLikelihood],
    haplotypes: &[Haplotype],
    ref_bytes: &[u8],
    pad_start_1based: u64,
    max_mnp: usize,
) -> Option<i32> {
    let mapping = create_allele_mapper(
        event,
        event.start_1based.get(),
        haplotypes,
        pad_start_1based,
        ref_bytes,
        max_mnp,
        true,
    );
    if mapping.alt_haplotype_indices.is_empty() {
        return None;
    }
    let rows = region_likelihoods_to_rows(likelihoods, haplotypes.len());
    if rows.is_empty() {
        return None;
    }
    let marg = marginalize_rows_to_biallelic_alleles(
        &rows,
        &mapping.ref_haplotype_indices,
        &mapping.alt_haplotype_indices,
    );
    let gls = biallelic_genotype_log10_likelihoods_gatk(&marg, 0, 1);
    let format = emit_genotype_format_fields(&gls, &[0, 0]).ok()?;
    if format.pl.len() < 3 {
        return None;
    }
    let d1 = format.pl[1].as_i32().saturating_sub(format.pl[0].as_i32());
    let d2 = format.pl[2].as_i32().saturating_sub(format.pl[0].as_i32());
    Some(d1.min(d2))
}

/// Max per-event HC inverse-PL among events each haplotype supports (higher = stronger).
fn haplotype_hc_inverse_pl_scores(
    assembly: &AssemblyResultSet,
    likelihoods: &[RegionReadLikelihood],
) -> Vec<f64> {
    let (ref_bytes, pad_start) = assembly.event_map_reference();
    let ref_hap = assembly
        .haplotypes
        .iter()
        .find(|h| h.is_reference)
        .unwrap_or(&assembly.haplotypes[0]);
    let max_mnp = assembly.max_mnp_distance();
    let mut out = vec![f64::NEG_INFINITY; assembly.haplotypes.len()];
    if likelihoods.is_empty() || assembly.variation_events.is_empty() {
        return out;
    }
    for (hi, h) in assembly.haplotypes.iter().enumerate() {
        if h.is_reference {
            continue;
        }
        let mut best = f64::NEG_INFINITY;
        for event in &assembly.variation_events {
            if !haplotype_supports_allele_at_with_ref(
                h,
                ref_hap,
                event.start_1based.get(),
                pad_start,
                &event.ref_allele,
                &event.alt_allele,
                ref_bytes,
                max_mnp,
                &event.contig,
            ) {
                continue;
            }
            if let Some(pl) = variation_event_hc_inverse_pl(
                event,
                likelihoods,
                &assembly.haplotypes,
                ref_bytes,
                pad_start,
                max_mnp,
            ) {
                best = best.max(pl as f64);
            }
        }
        out[hi] = best;
    }
    out
}

fn variation_events_in_active_span(
    assembly: &AssemblyResultSet,
    active_start_1based: Option<u64>,
    active_end_1based: Option<u64>,
) -> Vec<&VariationEvent> {
    assembly
        .variation_events
        .iter()
        .filter(|e| match (active_start_1based, active_end_1based) {
            (Some(s), Some(end)) => {
                e.start_1based >= GenomePosition::new_1based(s)
                    && e.start_1based <= GenomePosition::new_1based(end)
            }
            _ => true,
        })
        .collect()
}

fn mark_haplotypes_supporting_variation_events(
    assembly: &AssemblyResultSet,
    keep: &mut [bool],
    ref_idx: usize,
    options: crate::allele_filter_options::AlleleFilterOptions,
) {
    let strict_java_snp_rank_only = options.strict_java_snp_rank_only;
    let events = variation_events_in_active_span(
        assembly,
        options.active_start_1based(),
        options.active_end_1based(),
    );
    if events.is_empty() {
        return;
    }
    let (ref_bytes, pad_start) = assembly.event_map_reference();
    let ref_hap = &assembly.haplotypes[ref_idx];
    let apply_pad = ref_hap
        .genome_loc
        .map(|g| g.start_1based())
        .unwrap_or(pad_start);
    let max_mnp = assembly.max_mnp_distance();
    // R4-2: P12 SNP-rank / emit-band haplotype keep applies only on contig 2.
    // Elsewhere keep haplotypes that support any EventMap allele (hets + indels).
    let p12_scope = assembly.contig == "2" || assembly.contig == "chr2";
    if strict_java_snp_rank_only && p12_scope {
        for e in &events {
            if e.ref_allele.len() != 1 || e.alt_allele.len() != 1 {
                continue;
            }
            if !is_strict_java_production_emit_candidate(e)
                && !is_p12_phase_e_gap_event(e)
                && !is_cluster_anchor_snp(e)
            {
                continue;
            }
            let mut alt_supporters = Vec::new();
            for (i, h) in assembly.haplotypes.iter().enumerate() {
                if i == ref_idx || h.is_reference || keep[i] {
                    continue;
                }
                if haplotype_supports_allele_at_with_ref(
                    h,
                    ref_hap,
                    e.start_1based.get(),
                    pad_start,
                    &e.ref_allele,
                    &e.alt_allele,
                    ref_bytes,
                    max_mnp,
                    &e.contig,
                ) {
                    alt_supporters.push(i);
                }
            }
            let alt_byte = e.alt_allele.as_bytes().first().copied();
            let apply_off = e.start_1based.get().saturating_sub(apply_pad) as usize;
            let exact: Vec<usize> = alt_supporters
                .iter()
                .copied()
                .filter(|&i| {
                    alt_byte.is_some_and(|b| {
                        assembly.haplotypes[i].bases.get(apply_off) == Some(&b)
                            || hap_base_at_ref_locus(
                                &assembly.haplotypes[i],
                                pad_start,
                                e.start_1based.get(),
                            ) == Some(b)
                    })
                })
                .collect();
            if exact.len() == 1 {
                keep[exact[0]] = true;
            } else if alt_supporters.len() == 1 {
                keep[alt_supporters[0]] = true;
            } else if exact.is_empty()
                && (is_strict_java_production_emit_candidate(e)
                    || is_p12_phase_e_gap_event(e)
                    || is_cluster_anchor_snp(e))
            {
                if let Some((i, _)) = assembly.haplotypes.iter().enumerate().find(|(i, h)| {
                    *i != ref_idx
                        && !h.is_reference
                        && alt_byte.is_some_and(|b| h.bases.get(apply_off) == Some(&b))
                }) {
                    keep[i] = true;
                }
            }
        }
        for e in &events {
            if !is_cluster_ctc_del(e) {
                continue;
            }
            if let Some((i, _)) = assembly.haplotypes.iter().enumerate().find(|(i, h)| {
                *i != ref_idx
                    && !h.is_reference
                    && h.bases.len() + 1 == ref_hap.bases.len()
                    && h.cigar.as_ref().is_some_and(|c| {
                        c.elements
                            .iter()
                            .any(|e| e.operator == crate::cigar::CigarOperator::Deletion)
                            && crate::read_event_discovery::c_has_deletion_at_ref_offset(
                                c,
                                crate::read_event_discovery::p12_cluster_ctc_deletion_ref_offset(
                                    apply_pad,
                                ),
                            )
                    })
            }) {
                keep[i] = true;
            }
        }
        for (i, h) in assembly.haplotypes.iter().enumerate() {
            if i == ref_idx || h.is_reference || keep[i] {
                continue;
            }
            if h.cigar
                .as_ref()
                .is_some_and(|c| c.elements.iter().any(|e| e.operator.is_indel()))
            {
                for e in &events {
                    if !e.is_indel()
                        && haplotype_supports_allele_at_with_ref(
                            h,
                            ref_hap,
                            e.start_1based.get(),
                            pad_start,
                            &e.ref_allele,
                            &e.alt_allele,
                            ref_bytes,
                            max_mnp,
                            &e.contig,
                        )
                    {
                        keep[i] = true;
                        break;
                    }
                }
            }
        }
        return;
    }
    for (i, h) in assembly.haplotypes.iter().enumerate() {
        if i == ref_idx || h.is_reference || keep[i] {
            continue;
        }
        for e in &events {
            if haplotype_supports_allele_at_with_ref(
                h,
                ref_hap,
                e.start_1based.get(),
                pad_start,
                &e.ref_allele,
                &e.alt_allele,
                ref_bytes,
                max_mnp,
                &e.contig,
            ) {
                keep[i] = true;
                break;
            }
        }
    }
}

fn select_haplotype_keep_mask(
    assembly: &AssemblyResultSet,
    likelihoods: &[RegionReadLikelihood],
    options: crate::allele_filter_options::AlleleFilterOptions,
) -> GatkResult<Vec<bool>> {
    if assembly.haplotypes.is_empty() {
        return Ok(Vec::new());
    }
    let ref_idx = assembly
        .haplotypes
        .iter()
        .position(|h| h.is_reference)
        .ok_or_else(|| GatkError::algorithm("allele filter: missing reference haplotype"))?;
    let mut keep = vec![false; assembly.haplotypes.len()];
    keep[ref_idx] = true;
    let p12_scope = assembly.contig == "2" || assembly.contig == "chr2";
    // R4-2: retain indel-bearing haplotypes outside contig 2 under strict SNP-rank options.
    if !options.strict_java_snp_rank_only || !p12_scope {
        for (i, h) in assembly.haplotypes.iter().enumerate() {
            if i != ref_idx
                && !h.is_reference
                && h.cigar
                    .as_ref()
                    .is_some_and(|c| c.elements.iter().any(|e| e.operator.is_indel()))
            {
                keep[i] = true;
            }
        }
    }
    mark_haplotypes_supporting_variation_events(assembly, &mut keep, ref_idx, options);
    let sums = haplotype_log10_sums_from_region_likelihoods(likelihoods, assembly.haplotypes.len());
    let hc_scores = haplotype_hc_inverse_pl_scores(assembly, likelihoods);
    let mut scored: Vec<(usize, f64)> = assembly
        .haplotypes
        .iter()
        .enumerate()
        .filter(|(i, h)| *i != ref_idx && !h.is_reference && !keep[*i])
        .map(|(i, _)| {
            let ll_sum = sums.get(i).copied().unwrap_or(f64::NEG_INFINITY);
            let hc = hc_scores.get(i).copied().unwrap_or(f64::NEG_INFINITY);
            let score = if hc.is_finite() {
                hc + ll_sum * HAPLOTYPE_RANK_LL_SUM_TIEBREAK_SCALE
            } else {
                ll_sum
            };
            (i, score)
        })
        .collect();
    // R4-2: score-based top-N fill outside contig 2 even under strict SNP-rank options.
    let allow_score_fill =
        !options.strict_java_snp_rank_only || (assembly.contig != "2" && assembly.contig != "chr2");
    if allow_score_fill {
        scored.sort_by(|a, b| b.1.total_cmp(&a.1));
        let mut non_ref_kept = keep
            .iter()
            .enumerate()
            .filter(|(i, k)| *i != ref_idx && **k)
            .count();
        for (i, _) in scored {
            if non_ref_kept >= MAX_NON_REF_HAPLOTYPES_FOR_GENOTYPING {
                break;
            }
            keep[i] = true;
            non_ref_kept += 1;
        }
    }
    Ok(keep)
}

/// Keep reference + top-scoring non-ref haplotypes (Java filters noisy graph paths).
pub fn filter_haplotypes_for_genotyping(assembly: &mut AssemblyResultSet) -> GatkResult<()> {
    ensure_reference_haplotype(&mut assembly.haplotypes);
    let keep = select_haplotype_keep_mask(
        assembly,
        &[],
        crate::allele_filter_options::AlleleFilterOptions::unrestricted(),
    )?;
    let old = std::mem::take(&mut assembly.haplotypes);
    assembly.haplotypes = old
        .into_iter()
        .enumerate()
        .filter(|(i, _)| keep[*i])
        .map(|(_, h)| h)
        .collect();
    Ok(())
}

/// Remap read likelihood rows after haplotype indices were filtered.
pub fn remap_read_likelihoods_after_hap_filter(
    likelihoods: Vec<RegionReadLikelihood>,
    old_to_new: &[Option<usize>],
) -> Vec<RegionReadLikelihood> {
    likelihoods
        .into_iter()
        .filter_map(|rl| {
            old_to_new.get(rl.haplotype_index.get()).and_then(|o| {
                o.map(|new_hi| RegionReadLikelihood {
                    haplotype_index: crate::bio_ids::HaplotypeIndex::new(new_hi),
                    ..rl
                })
            })
        })
        .collect()
}

/// Filter assembly haplotypes and remap likelihood indices.
pub fn filter_assembly_and_likelihoods(
    assembly: &mut AssemblyResultSet,
    likelihoods: Vec<RegionReadLikelihood>,
    options: crate::allele_filter_options::AlleleFilterOptions,
) -> GatkResult<Vec<RegionReadLikelihood>> {
    ensure_reference_haplotype(&mut assembly.haplotypes);
    let old_len = assembly.haplotypes.len();
    let keep = select_haplotype_keep_mask(assembly, &likelihoods, options)?;
    let mut old_to_new = vec![None; old_len];
    let mut haplotypes = Vec::new();
    let old = std::mem::take(&mut assembly.haplotypes);
    for (i, h) in old.into_iter().enumerate() {
        if keep[i] {
            old_to_new[i] = Some(haplotypes.len());
            haplotypes.push(h);
        }
    }
    assembly.haplotypes = haplotypes;
    Ok(remap_read_likelihoods_after_hap_filter(
        likelihoods,
        &old_to_new,
    ))
}
