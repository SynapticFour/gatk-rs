//! PairHMM region scoring used by `call_region` (split from `engine.rs` for N-3).
use super::*;
use crate::assembly_region_finalize::{
    clip_finalized_reads_in_place, finalize_region_reads_for_assembly,
    gatk_min_tail_quality_for_assembly,
};
use crate::likelihood_engine::score_read_against_haplotypes;

pub(super) fn compute_region_read_likelihoods(
    region: &AssemblyRegion,
    haplotypes: &[Haplotype],
    config: &HcLikelihoodEngineConfig,
    apply_normalize: bool,
    pre_finalized: Option<Vec<rust_htslib::bam::Record>>,
) -> GatkResult<Vec<RegionReadLikelihood>> {
    if haplotypes.is_empty() {
        return Ok(Vec::new());
    }
    // A2: consume assemble finalize buffer when present (clip in place — no second owned copy).
    let mut finalized = if let Some(mut pre) = pre_finalized.filter(|p| !p.is_empty()) {
        clip_finalized_reads_in_place(&mut pre, region);
        pre
    } else {
        finalize_region_reads_for_assembly(
            &region.reads,
            region,
            true,
            gatk_min_tail_quality_for_assembly(10),
            false,
        )
    };
    // Java `callRegion` drops these stubs on `regionForGenotyping` before PairHMM.
    finalized.retain(|r| unclipped_read_length(r) >= GATK_MINIMUM_READ_LENGTH_AFTER_TRIMMING);
    let active_span = Some((region.start.get(), region.end.get()));
    // Trim/hard-clip can drop sparse-BAM reads that still overlap the active locus.
    if finalized.is_empty() && !region.reads.is_empty() {
        let out = score_pairhmm_from_records(region.reads.as_slice(), haplotypes, config)?;
        return Ok(post_process_pairhmm_likelihoods(
            out,
            region.reads.as_slice(),
            haplotypes,
            apply_normalize,
            active_span,
        ));
    }
    let out = score_pairhmm_from_records(&finalized, haplotypes, config)?;
    Ok(post_process_pairhmm_likelihoods(
        out,
        &finalized,
        haplotypes,
        apply_normalize,
        active_span,
    ))
}

/// GATK 4.4 `ReadUtils` BI/BD FastQ-33 string → Phred. Absent or non-string → None (Q45).
fn bam_bqsr_indel_quals_phred(rec: &rust_htslib::bam::Record, tag: &[u8]) -> Option<Vec<u8>> {
    match rec.aux(tag) {
        Ok(rust_htslib::bam::record::Aux::String(s)) => {
            Some(s.bytes().map(|b| b.saturating_sub(33)).collect())
        }
        _ => None,
    }
}

/// Score PairHMM from BAM records without `AssemblyRead` / UTF-8 `String` rematerialization.
///
/// # Observable contract
/// Same finalizeRegion evidence and PairHMM inputs as the prior `records_to_assembly_reads` path
/// (BAM seq/qual bytes are ASCII ACGTN — identical to `String::from_utf8_lossy` for valid records).
fn score_pairhmm_from_records<R: std::borrow::Borrow<rust_htslib::bam::Record> + Sync>(
    reads: &[R],
    haplotypes: &[Haplotype],
    config: &HcLikelihoodEngineConfig,
) -> GatkResult<Vec<RegionReadLikelihood>> {
    let _prof = crate::hc_profile::begin(crate::hc_profile::Stage::PairHmm);
    let wall0 = std::time::Instant::now();
    let eligible = pairhmm_eligible_haplotype_indices(haplotypes);
    // L12-A3: zero-copy hap membership for PairHMM (no post-prune `Vec<u8>` rematerialize).
    let hap_refs: Vec<&[u8]> = eligible
        .iter()
        .map(|&hi| haplotypes[hi].bases.as_slice())
        .collect();
    // Parallel across reads when the rayon pool has >1 worker (Java `--native-pair-hmm-threads`).
    // `GATK_RS_HC_SEQUENTIAL` only serializes *regions* for Peak-RSS — PairHMM within a region
    // stays threaded so we can undercut Java wall without stacking mid-size regions.
    // Hap scoring inside each read stays sequential when nested (one parallel axis).
    // Keep dump row-groups contiguous (Java processedReads order) when capturing inputs.
    let parallel = rayon::current_num_threads() > 1
        && reads.len() >= 8
        && crate::runtime_config::pairhmm_input_dump_path().is_none();
    let out = if !parallel {
        let mut out = Vec::with_capacity(reads.len() * eligible.len());
        for (ri, rec) in reads.iter().enumerate() {
            let rec = rec.borrow();
            let bases = rec.seq().as_bytes();
            let ins_tag = bam_bqsr_indel_quals_phred(rec, b"BI");
            let del_tag = bam_bqsr_indel_quals_phred(rec, b"BD");
            let scores = score_read_against_haplotypes(
                config,
                &bases,
                rec.qual(),
                rec.mapq(),
                &hap_refs,
                ins_tag.as_deref(),
                del_tag.as_deref(),
            )?;
            for (score_i, &hi) in eligible.iter().enumerate() {
                out.push(RegionReadLikelihood {
                    read_index: crate::bio_ids::ReadIndex::new(ri),
                    haplotype_index: crate::bio_ids::HaplotypeIndex::new(hi),
                    log10_likelihood: scores[score_i],
                });
            }
        }
        out
    } else {
        // Parallel across reads (Java native PairHMM threads). Each rayon worker has its own
        // PairHMM TLS; collect then flatten in read-index order for stable LL rows.
        use rayon::prelude::*;
        let per_read: Vec<GatkResult<Vec<RegionReadLikelihood>>> = reads
            .par_iter()
            .enumerate()
            .map(|(ri, rec)| {
                let rec = rec.borrow();
                let bases = rec.seq().as_bytes();
                let ins_tag = bam_bqsr_indel_quals_phred(rec, b"BI");
                let del_tag = bam_bqsr_indel_quals_phred(rec, b"BD");
                let scores = score_read_against_haplotypes(
                    config,
                    &bases,
                    rec.qual(),
                    rec.mapq(),
                    &hap_refs,
                    ins_tag.as_deref(),
                    del_tag.as_deref(),
                )?;
                let mut rows = Vec::with_capacity(eligible.len());
                for (score_i, &hi) in eligible.iter().enumerate() {
                    rows.push(RegionReadLikelihood {
                        read_index: crate::bio_ids::ReadIndex::new(ri),
                        haplotype_index: crate::bio_ids::HaplotypeIndex::new(hi),
                        log10_likelihood: scores[score_i],
                    });
                }
                Ok(rows)
            })
            .collect();
        let mut out = Vec::with_capacity(reads.len() * eligible.len());
        for chunk in per_read {
            out.extend(chunk?);
        }
        out
    };

    if crate::hc_profile::enabled() {
        crate::hc_profile::note_pairhmm_region(reads, &hap_refs, wall0.elapsed());
    }
    Ok(out)
}
