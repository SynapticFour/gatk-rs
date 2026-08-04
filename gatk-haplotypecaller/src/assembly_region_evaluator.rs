//! Per-locus HC activity before band-pass smoothing.
//! Production path is free functions [`evaluate_hc_activity_state`] and
//! [`add_locus_for_smoothed_activity`] (no evaluator trait / ZST shell).

use crate::activity_profile::{ActivityProfileState, BandPassActivityProfile};
use crate::activity_scoring::{
    haplotype_caller_activity_profile_state_single_sample, HaplotypeCallerActivityScoringParams,
    PileupObservation,
};
use crate::allele_downsample::apply_contamination_to_pileup;
use crate::gatk_well_rng::Well19937c;
use crate::locus_iterator::{pileup_observation_from_record, LocusPileupState};
use crate::minimal_genotyping::haplotype_caller_activity_profile_state_minimal_genotyping;
use crate::read_header_semantics::ReadHeaderSemantics;
use crate::read_model::{passes_hc_read_filters_with_header, ReadFilterParams};
use crate::shared_bam::SharedBamRecord;
use gatk_common::GatkResult;
use rust_htslib::bam;

/// Shared HC activity path used by the iterator and TSV dump.
pub fn evaluate_hc_activity_state(
    contig: &str,
    pos1: u64,
    pile: &[PileupObservation],
    scoring: &HaplotypeCallerActivityScoringParams,
) -> ActivityProfileState {
    let mut pile = pile.to_vec();
    if scoring.contamination_fraction_to_filter > 0.0 {
        let mut rng = Well19937c::reset_gatk_default();
        apply_contamination_to_pileup(
            &mut pile,
            scoring.contamination_fraction_to_filter,
            &mut rng,
        );
    }
    haplotype_caller_activity_profile_state_single_sample(contig, pos1, &pile, scoring)
}

/// Feed one 1-based locus into the band-pass profile (iterator + activity TSV dumps).
/// When `pileup_state` is `Some`, uses LIBS-style incremental pileup (GAP-B-02); otherwise scans all records.
/// `force_active` mirrors GATK `--alleles` force-calling overlap at this locus.
pub fn add_locus_for_smoothed_activity(
    prof: &mut BandPassActivityProfile,
    records: &[SharedBamRecord],
    header: &bam::HeaderView,
    header_semantics: &ReadHeaderSemantics,
    contig: &str,
    pos1: u64,
    read_filters: &ReadFilterParams,
    ref_base: u8,
    scoring: &HaplotypeCallerActivityScoringParams,
    pileup_state: Option<&mut LocusPileupState>,
    force_active: bool,
) -> GatkResult<()> {
    let st = match pileup_state {
        Some(st) => {
            st.advance_to(records, read_filters, pos1)?;
            hc_activity_after_locus_advance(
                contig,
                pos1,
                st,
                records,
                header_semantics,
                scoring,
                ref_base,
                force_active,
            )?
        }
        None => {
            let pile = pileup_at_locus(records, header, contig, pos1, read_filters, ref_base)?;
            evaluate_hc_activity_state(contig, pos1, &pile, scoring)
        }
    };

    prof.add(st)?;
    Ok(())
}

/// Incremental pileup → activity state (LIBS path). `force_active` sets active_prob=1.0.
pub(crate) fn hc_activity_after_locus_advance(
    contig: &str,
    pos1: u64,
    pileup_state: &mut LocusPileupState,
    records: &[SharedBamRecord],
    header_semantics: &ReadHeaderSemantics,
    scoring: &HaplotypeCallerActivityScoringParams,
    ref_base: u8,
    force_active: bool,
) -> GatkResult<ActivityProfileState> {
    if force_active {
        return Ok(ActivityProfileState {
            contig: std::sync::Arc::from(contig),
            pos: pos1,
            active_prob: 1.0,
            original_active_prob: 1.0,
            evidence: crate::activity_profile::ActivityEvidence::Plain,
        });
    }
    let mut pile = pileup_state.pileup_observations(records, ref_base)?;
    if pile.is_empty() {
        return Ok(ActivityProfileState::new(contig, pos1, 0.0));
    }
    if scoring.contamination_fraction_to_filter > 0.0 {
        let mut rng = Well19937c::reset_gatk_default();
        apply_contamination_to_pileup(
            &mut pile,
            scoring.contamination_fraction_to_filter,
            &mut rng,
        );
    }

    // R4-1: single-sample HC never uses stratified piles — skip the second observation walk.
    if header_semantics.is_single_sample_header() {
        return Ok(haplotype_caller_activity_profile_state_minimal_genotyping(
            contig, pos1, &pile, scoring,
        ));
    }

    let mut stratified = pileup_state.nonempty_stratified_sample_pileups_ordered(
        records,
        header_semantics,
        ref_base,
    )?;
    if scoring.contamination_fraction_to_filter > 0.0 {
        let mut rng = Well19937c::reset_gatk_default();
        for sample_pile in &mut stratified {
            apply_contamination_to_pileup(
                sample_pile,
                scoring.contamination_fraction_to_filter,
                &mut rng,
            );
        }
    }

    Ok(if stratified.len() <= 1 {
        haplotype_caller_activity_profile_state_minimal_genotyping(contig, pos1, &pile, scoring)
    } else {
        {
            let strata: Vec<&[_]> = stratified.iter().map(|p| p.as_slice()).collect();
            crate::hc_joint_is_active::haplotype_caller_joint_multisample_is_active_activity_state(
                contig, pos1, &strata, scoring,
            )
        }
    })
}

fn pileup_at_locus(
    records: &[SharedBamRecord],
    header: &bam::HeaderView,
    contig: &str,
    pos1: u64,
    filters: &ReadFilterParams,
    ref_base: u8,
) -> GatkResult<Vec<PileupObservation>> {
    let ref_pos0 = pos1.saturating_sub(1) as i64;
    let mut obs = Vec::new();

    for rec in records {
        if !passes_hc_read_filters_with_header(rec, header, filters) {
            continue;
        }
        let rn = String::from_utf8_lossy(header.tid2name(rec.tid() as u32)).into_owned();
        if rn != contig {
            continue;
        }
        if let Some(o) = pileup_observation_from_record(rec, ref_pos0, ref_base, false)? {
            obs.push(o);
        }
    }

    Ok(obs)
}
