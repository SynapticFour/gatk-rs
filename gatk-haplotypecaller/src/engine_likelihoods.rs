//! PairHMM region scoring used by `call_region` (split from `engine.rs` for N-3).
use super::*;
use crate::assembly_region_finalize::{
    clip_finalized_reads_in_place, finalize_region_reads_for_assembly,
    gatk_min_tail_quality_for_assembly,
};
use crate::likelihood_engine::score_read_against_haplotypes;

/// Score PairHMM and return the exact evidence list that owns `read_index`.
///
/// Java 4.4 `AlleleLikelihoods` is constructed from `perSampleReadList` and that same
/// evidence vector is what `filterPoorlyModeledEvidence` indexes (`sampleEvidence.get(i)`
/// ↔ `valuesBySampleIndex[a][i]`). Callers must keep the returned reads as the filter /
/// realign / genotyping evidence list.
///
/// PairHMM membership applies Java `MATE_ON_SAME_CONTIG_OR_NO_MAPPED_MATE` (6R.311):
/// a paired read whose mapped mate is on another contig is omitted. Other
/// filtered reads still receive `addEvidence(..., 0)` cells.
pub(super) fn compute_region_read_likelihoods(
    region: &AssemblyRegion,
    haplotypes: &[Haplotype],
    config: &HcLikelihoodEngineConfig,
    apply_normalize: bool,
    pre_finalized: Option<Vec<rust_htslib::bam::Record>>,
) -> GatkResult<(
    Vec<RegionReadLikelihood>,
    Vec<crate::shared_bam::SharedBamRecord>,
)> {
    if haplotypes.is_empty() {
        return Ok((Vec::new(), Vec::new()));
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
    // Mate contig lives on the pre-filter originals. Clipped copies have no `mtid`.
    finalized.retain(|r| {
        if !super::production_failed_mate_key(r.qname(), r.flags()) {
            return true;
        }
        super::note_production_failed_mate_exclusion(r.qname(), r.flags());
        false
    });
    super::note_forensic_6r307_records("pairhmm_input", &finalized);
    let active_span = Some((region.start.get(), region.end.get()));
    // Trim/hard-clip can drop sparse-BAM reads that still overlap the active locus.
    if finalized.is_empty() && !region.reads.is_empty() {
        let out = score_pairhmm_from_records_java_mate_contig(
            region.reads.as_slice(),
            haplotypes,
            config,
            region.reads.as_slice(),
        )?;
        capture_scored_likelihood_pipeline(&out, region.reads.as_slice(), haplotypes);
        let ll = post_process_pairhmm_likelihoods(
            out,
            region.reads.as_slice(),
            haplotypes,
            apply_normalize,
            active_span,
        );
        return Ok((ll, region.reads.clone()));
    }
    let out = score_pairhmm_from_records_java_mate_contig(
        &finalized,
        haplotypes,
        config,
        region.reads.as_slice(),
    )?;
    capture_scored_likelihood_pipeline(&out, &finalized, haplotypes);
    let ll =
        post_process_pairhmm_likelihoods(out, &finalized, haplotypes, apply_normalize, active_span);
    let scored = finalized
        .into_iter()
        .map(crate::shared_bam::share_record)
        .collect();
    Ok((ll, scored))
}

/// GATK `filterNonPassingReads` mate-contig membership at PairHMM construction.
///
/// Java removes `!MATE_ON_SAME_CONTIG_OR_NO_MAPPED_MATE` reads before
/// `computeReadLikelihoods`. Those reads are dropped from this evidence vector
/// and do not receive a likelihood row. Length, mapping-quality, and read-group
/// removals still receive `addEvidence(..., 0)` cells.
///
/// Mate contig is taken from the original BAM evidence (`mate_originals`), not
/// the clipped PairHMM copy: `hard_clip` rebuilds records without `mtid`.
/// A scored read missing from that list was already dropped by
/// `filterNonPassingReads` and receives `addEvidence(..., 0)` rather than PairHMM.
fn score_pairhmm_from_records_java_mate_contig<
    R: std::borrow::Borrow<rust_htslib::bam::Record> + Sync,
>(
    reads: &[R],
    haplotypes: &[Haplotype],
    config: &HcLikelihoodEngineConfig,
    mate_originals: &[crate::shared_bam::SharedBamRecord],
) -> GatkResult<Vec<RegionReadLikelihood>> {
    use crate::read_pre_mate::passes_mate_on_same_contig_or_no_mapped_mate;
    use std::collections::HashMap;
    let original_pass: HashMap<(Vec<u8>, u16), bool> = mate_originals
        .iter()
        .map(|r| {
            (
                (r.qname().to_vec(), r.flags()),
                passes_mate_on_same_contig_or_no_mapped_mate(r),
            )
        })
        .collect();
    let pass: Vec<bool> = reads
        .iter()
        .map(|r| {
            let rec = r.borrow();
            // Post-`filterNonPassingReads` evidence still has BAM mate contig. Assemble
            // finalize copies lose `mtid`, so identity is QNAME+FLAG. A scored read that
            // is absent from that list was filtered (mate-contig / MAPQ / RG) and must
            // not enter PairHMM (`addEvidence(..., 0)`).
            original_pass
                .get(&(rec.qname().to_vec(), rec.flags()))
                .copied()
                .unwrap_or(false)
        })
        .collect();
    if pass.iter().all(|&ok| ok) {
        return score_pairhmm_from_records(reads, haplotypes, config);
    }
    let passing: Vec<&rust_htslib::bam::Record> = reads
        .iter()
        .zip(pass.iter())
        .filter(|(_, ok)| **ok)
        .map(|(r, _)| r.borrow())
        .collect();
    let orig_idx: Vec<usize> = pass
        .iter()
        .enumerate()
        .filter(|(_, ok)| **ok)
        .map(|(i, _)| i)
        .collect();
    let mut out = if passing.is_empty() {
        Vec::new()
    } else {
        score_pairhmm_from_records(&passing, haplotypes, config)?
    };
    for cell in &mut out {
        let k = cell.read_index.get();
        cell.read_index = crate::bio_ids::ReadIndex::new(orig_idx[k]);
    }
    let eligible = super::pairhmm_eligible_haplotype_indices(haplotypes);
    for (i, ok) in pass.iter().enumerate() {
        if *ok {
            continue;
        }
        let rec = reads[i].borrow();
        if super::forensic_6r306_exclude_zero_row(rec.qname(), rec.flags()) {
            continue;
        }
        for &hi in &eligible {
            out.push(RegionReadLikelihood {
                read_index: crate::bio_ids::ReadIndex::new(i),
                haplotype_index: crate::bio_ids::HaplotypeIndex::new(hi),
                log10_likelihood: 0.0,
            });
        }
    }
    Ok(out)
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
        && crate::runtime_config::pairhmm_input_dump_path().is_none()
        && !crate::likelihood_engine::forensic_6r289_active();
    let out = if !parallel {
        let mut out = Vec::with_capacity(reads.len() * eligible.len());
        for (ri, rec) in reads.iter().enumerate() {
            let rec = rec.borrow();
            let bases = rec.seq().as_bytes();
            let ins_tag = bam_bqsr_indel_quals_phred(rec, b"BI");
            let del_tag = bam_bqsr_indel_quals_phred(rec, b"BD");
            if crate::likelihood_engine::forensic_6r289_matches(rec.qname(), rec.flags()) {
                let slot = eligible.iter().position(|&hi| hi == 0).unwrap_or(0);
                crate::likelihood_engine::forensic_6r289_arm(
                    rec.mapq(),
                    rec.pos() + 1,
                    format!("{}", rec.cigar()),
                    slot,
                );
            }
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

#[cfg(test)]
mod holdout_6r314 {
    use super::super::forensic_6r306_note_pre_filter_reads;
    use super::score_pairhmm_from_records_java_mate_contig;

    fn holdout_6r314_rec(
        qname: &str,
        paired: bool,
        unmapped: bool,
        mate_unmapped: bool,
        mtid: i32,
        mapq: u8,
    ) -> rust_htslib::bam::Record {
        use rust_htslib::bam::record::{Cigar, CigarString};
        use rust_htslib::bam::{HeaderView, Record};
        use std::sync::Arc;
        let mut r = Record::new();
        r.set_header(Arc::new(HeaderView::from_bytes(
            b"@HD\tVN:1.0\n@SQ\tSN:hold-a\tLN:1000\n@SQ\tSN:hold-b\tLN:1000\n",
        )));
        r.set(
            qname.as_bytes(),
            Some(&CigarString::from(vec![Cigar::Match(16)])),
            b"ACGTACGTACGTACGT",
            &vec![30u8; 16],
        );
        r.set_tid(0);
        r.set_pos(100);
        r.set_mtid(mtid);
        r.set_mapq(mapq);
        if paired {
            r.set_paired();
        }
        if unmapped {
            r.set_unmapped();
        } else {
            r.unset_unmapped();
        }
        if mate_unmapped {
            r.set_mate_unmapped();
        } else {
            r.unset_mate_unmapped();
        }
        r
    }

    /// 6R.314 failed-mate checkpoint: mate-fail reads get no likelihood row.
    /// A read absent from the pre-filter originals (length / MAPQ / read-group)
    /// still receives `0.0` cells. Runs in the lib suite; `src/` cannot read
    /// `HOLDOUT_6R314` (`std::env::var` is allowlisted).
    #[test]
    fn holdout_6r314_failed_mate_zero_row() {
        use crate::haplotype::Haplotype;
        use crate::likelihood_engine::HcLikelihoodEngineConfig;
        use crate::shared_bam::share_record;

        let keeper = holdout_6r314_rec("hold314-same", true, false, false, 0, 20);
        let low_mapq = holdout_6r314_rec("hold314-low-mapq", true, false, false, 0, 1);
        let cross = holdout_6r314_rec("hold314-cross", true, false, false, 1, 20);
        let ghost = holdout_6r314_rec("hold314-length-absent", true, false, false, 0, 20);
        let originals = [
            share_record(keeper.clone()),
            share_record(low_mapq.clone()),
            share_record(cross.clone()),
        ];
        forensic_6r306_note_pre_filter_reads(&originals);
        let hap = Haplotype::new(b"ACGTACGTACGTACGT".to_vec(), true);
        let scored = [
            keeper.clone(),
            low_mapq.clone(),
            cross.clone(),
            ghost.clone(),
        ];
        let out = score_pairhmm_from_records_java_mate_contig(
            &scored,
            &[hap],
            &HcLikelihoodEngineConfig::gatk_haplotype_caller_production(),
            &originals,
        )
        .expect("pairhmm");
        let rows = |idx: usize| -> Vec<f64> {
            out.iter()
                .filter(|c| c.read_index.get() == idx)
                .map(|c| c.log10_likelihood)
                .collect()
        };
        let keeper_rows = rows(0);
        let mapq_rows = rows(1);
        let cross_rows = rows(2);
        let ghost_rows = rows(3);
        assert!(!keeper_rows.is_empty(), "same-contig read stays in PairHMM");
        assert!(
            keeper_rows.iter().any(|v| *v != 0.0),
            "same-contig row is a PairHMM score, not an inserted zero"
        );
        assert!(!mapq_rows.is_empty(), "low MAPQ is not a mate exclusion");
        assert!(
            cross_rows.is_empty(),
            "different-contig mate has no zero row"
        );
        assert!(
            !ghost_rows.is_empty(),
            "length/MQ/RG-style miss still has cells"
        );
        assert!(
            ghost_rows.iter().all(|v| *v == 0.0),
            "absent original receives 0.0 cells, got {ghost_rows:?}"
        );
        forensic_6r306_note_pre_filter_reads(&[]);
        println!("6R314\tcheckpoint\tzero_row\tPASS");
    }
}
