//! GATK `AssemblyBasedCallerUtils.realignReadsToTheirBestHaplotype` (P1⁴ REALIGN).

use crate::alignment::SwParameters;
use crate::cigar::{Cigar, CigarElement, CigarOperator};
use crate::cigar_builder::CigarBuilder;
use crate::haplotype::Haplotype;
use crate::haplotype_cigar::{
    apply_cigar_to_cigar, calculate_haplotype_cigar_for_assembly_with_offset,
    consolidated_padded_cigar, left_align_indels_for_read, read_start_on_reference_haplotype,
    trim_cigar_by_bases_public,
};
use crate::read_unclip::hard_clip_soft_clipped_bases;
use crate::region_read_likelihood::RegionReadLikelihood;
use crate::smith_waterman::{align_read_to_best_haplotype, SwParameters as SwParams};
use gatk_common::GatkResult;
use rust_htslib::bam::record::{Cigar as HtsCigar, CigarString};
use rust_htslib::bam::Record;

/// GATK `AlleleLikelihoods.LOG_10_INFORMATIVE_THRESHOLD`.
pub const LOG_10_INFORMATIVE_THRESHOLD: f64 = 0.2;

/// GATK `AssemblyBasedCallerUtils.HAPLOTYPE_ALIGNMENT_TIEBREAKING_PRIORITY`.
pub fn haplotype_alignment_tiebreak_priority(h: &Haplotype) -> f64 {
    let reference_term = if h.is_reference { 1.0 } else { 0.0 };
    let cigar_term = h
        .cigar
        .as_ref()
        .map(|c| 1.0 - c.elements.len() as f64)
        .unwrap_or(0.0);
    reference_term + cigar_term
}

/// Realign each read to its highest-likelihood haplotype projected on the padded reference.
/// GATK order: after `filterAlleles`, before genotyping; `changeEvidence` swaps read objects only.
pub fn realign_reads_to_best_haplotype(
    reads: &mut [Record],
    haplotypes: &[Haplotype],
    read_likelihoods: &[RegionReadLikelihood],
    padded_reference_start_1based: u64,
    hap_to_ref_sw: &SwParameters,
) -> GatkResult<(bool, Vec<usize>)> {
    if reads.is_empty() || haplotypes.is_empty() || read_likelihoods.is_empty() {
        return Ok((false, Vec::new()));
    }
    let ref_idx = haplotypes.iter().position(|h| h.is_reference).unwrap_or(0);
    let ref_hap = &haplotypes[ref_idx];
    let ref_bases = &ref_hap.bases;
    if ref_bases.is_empty() || padded_reference_start_1based == 0 {
        return Ok((false, vec![ref_idx; reads.len()]));
    }

    let priorities: Vec<f64> = haplotypes
        .iter()
        .map(haplotype_alignment_tiebreak_priority)
        .collect();
    let mut changed = false;
    let mut best_per_read = vec![ref_idx; reads.len()];
    for (ri, rec) in reads.iter_mut().enumerate() {
        if rec.is_unmapped() {
            continue;
        }
        let best_hi = best_haplotype_index_for_read(ri, haplotypes, read_likelihoods, &priorities);
        best_per_read[ri] = best_hi;
        if create_read_aligned_to_ref(
            rec,
            &haplotypes[best_hi],
            ref_hap,
            padded_reference_start_1based,
            hap_to_ref_sw,
        )? {
            changed = true;
        }
    }
    Ok((changed, best_per_read))
}

/// GATK `AlleleLikelihoods.searchBestAllele` + `HAPLOTYPE_ALIGNMENT_TIEBREAKING_PRIORITY`.
pub(crate) fn best_haplotype_index_for_read(
    read_index: usize,
    haplotypes: &[Haplotype],
    read_likelihoods: &[RegionReadLikelihood],
    priorities: &[f64],
) -> usize {
    let hap_count = haplotypes.len();
    if hap_count == 0 {
        return 0;
    }
    let ll = |hi: usize| -> f64 {
        read_likelihoods
            .iter()
            .find(|r| r.read_index.get() == read_index && r.haplotype_index.get() == hi)
            .map(|r| r.log10_likelihood)
            .unwrap_or(f64::NEG_INFINITY)
    };

    let mut best_a = 0usize;
    let mut second_a = 0usize;
    let mut best_ll = ll(0);
    let mut second_ll = f64::NEG_INFINITY;
    for a in 1..hap_count {
        let candidate = ll(a);
        if candidate > best_ll {
            second_a = best_a;
            second_ll = best_ll;
            best_a = a;
            best_ll = candidate;
        } else if candidate > second_ll {
            second_a = a;
            second_ll = candidate;
        }
    }

    if best_ll - second_ll < LOG_10_INFORMATIVE_THRESHOLD {
        let mut tie_best = best_a;
        let mut tie_second = second_a;
        let mut tie_best_pri = priorities.get(best_a).copied().unwrap_or(0.0);
        let mut tie_second_pri = priorities.get(second_a).copied().unwrap_or(0.0);
        for a in 0..hap_count {
            let candidate = ll(a);
            if a == best_a || best_ll - candidate > LOG_10_INFORMATIVE_THRESHOLD {
                continue;
            }
            let pri = priorities.get(a).copied().unwrap_or(0.0);
            if pri > tie_best_pri {
                tie_second = tie_best;
                tie_second_pri = tie_best_pri;
                tie_best = a;
                tie_best_pri = pri;
            } else if pri > tie_second_pri {
                tie_second = a;
                tie_second_pri = pri;
            }
        }
        best_a = tie_best;
        second_a = tie_second;
        best_ll = ll(best_a);
        second_ll = if second_a != best_a {
            ll(second_a)
        } else {
            f64::NEG_INFINITY
        };
    }
    let _ = (best_ll, second_ll);
    best_a
}

/// GATK `changeEvidence`: read objects were updated in realign; **PairHMM matrix unchanged**.
pub fn change_evidence_to_best_haplotype(
    likelihoods: Vec<RegionReadLikelihood>,
    _best_hap_per_read: &[usize],
) -> Vec<RegionReadLikelihood> {
    likelihoods
}

fn consolidate_cigar(cigar: &Cigar) -> Cigar {
    let mut builder = CigarBuilder::new(false);
    for e in &cigar.elements {
        builder.add(*e);
    }
    builder.make_and_record().cigar
}

fn edge_clips_from_record(rec: &Record) -> (Vec<CigarElement>, Vec<CigarElement>) {
    let mut leading = Vec::new();
    for hts in rec.cigar().iter() {
        if let Some((len, op)) = clip_from_hts(*hts) {
            leading.push(CigarElement {
                length: len,
                operator: op,
            });
        } else {
            break;
        }
    }
    let mut rev_trailing = Vec::new();
    for hts in rec.cigar().iter().rev() {
        if let Some((len, op)) = clip_from_hts(*hts) {
            rev_trailing.push(CigarElement {
                length: len,
                operator: op,
            });
        } else {
            break;
        }
    }
    rev_trailing.reverse();
    (leading, rev_trailing)
}

/// GATK `AlignmentUtils.appendClippedElementsFromCigarToCigar`.
fn append_clipped_elements_from_original(
    aligned: Cigar,
    leading: &[CigarElement],
    trailing: &[CigarElement],
) -> Cigar {
    let mut out = Cigar::new();
    for e in leading {
        out.push(e.length, e.operator);
    }
    for e in &aligned.elements {
        out.push(e.length, e.operator);
    }
    for e in trailing {
        out.push(e.length, e.operator);
    }
    consolidate_cigar(&out)
}

/// GATK `AlignmentUtils.createReadAlignedToRef` (updates read position + CIGAR on padded ref).
fn create_read_aligned_to_ref(
    rec: &mut Record,
    best_hap: &Haplotype,
    ref_hap: &Haplotype,
    reference_start_1based: u64,
    hap_to_ref_sw: &SwParameters,
) -> GatkResult<bool> {
    const CONSOLIDATED_HAP_PAD: usize = 1000;

    let original_ref_start_1based = rec.pos().max(0) as u64 + 1;
    let read_bases: Vec<u8> = rec.seq().as_bytes().to_vec();
    if read_bases.is_empty() {
        return Ok(false);
    }
    let ref_bases = &ref_hap.bases;
    let (leading_clips, trailing_clips) = edge_clips_from_record(rec);

    let read_minus_soft = hard_clip_soft_clipped_bases(rec);
    let clipped_bases: Vec<u8> = read_minus_soft.seq().as_bytes().to_vec();
    if clipped_bases.is_empty() {
        return Ok(false);
    }

    let read_to_hap_sw = SwParams::gatk_read_to_best_haplotype();
    let read_to_hap =
        match align_read_to_best_haplotype(&best_hap.bases, &clipped_bases, &read_to_hap_sw) {
            Ok(aln) => aln,
            Err(_) => return Ok(false),
        };
    if read_to_hap.alignment_offset < 0 {
        return Ok(false);
    }

    let read_start_hap = read_to_hap.alignment_offset.max(0) as usize;
    let read_to_hap_cigar = consolidate_cigar(&read_to_hap.cigar);

    let (hap_cigar, hap_offset) = if let Some(ref hc) = best_hap.cigar {
        (hc.clone(), best_hap.alignment_start_hap_wrt_ref)
    } else {
        let ref_cigar_len = ref_hap
            .cigar
            .as_ref()
            .map(|c| c.reference_length())
            .unwrap_or_else(|| ref_bases.len());
        let Some(hap_assy) = calculate_haplotype_cigar_for_assembly_with_offset(
            ref_bases,
            &best_hap.bases,
            ref_cigar_len,
            hap_to_ref_sw,
        ) else {
            return Ok(false);
        };
        (hap_assy.cigar, hap_assy.alignment_start_hap_wrt_ref)
    };

    let consolidated = consolidated_padded_cigar(&hap_cigar, CONSOLIDATED_HAP_PAD);
    let read_start_on_ref_hap = read_start_on_reference_haplotype(&consolidated, read_start_hap);
    let read_start_on_ref_1based = reference_start_1based
        .saturating_add(hap_offset as u64)
        .saturating_add(read_start_on_ref_hap as u64);

    // Java trims hap→ref to consolidated padded end, not read end on hap.
    let consolidated_hap_end = consolidated.read_length().saturating_sub(1);
    let hap_to_ref =
        trim_cigar_by_bases_public(&consolidated, read_start_hap, consolidated_hap_end).cigar;
    let read_to_ref = apply_cigar_to_cigar(&read_to_hap_cigar, &hap_to_ref);
    // Java: leftAlignIndels(..., refHaplotype.getBases,..., readStartOnReferenceHaplotype)
    let left = left_align_indels_for_read(
        &read_to_ref,
        ref_bases,
        &clipped_bases,
        read_start_on_ref_hap,
    );
    let read_start_on_ref_1based =
        read_start_on_ref_1based.saturating_add(left.leading_deletions_removed as u64);

    let hap_ref_start_1based = reference_start_1based.saturating_add(hap_offset as u64);
    // GATK `lastIndexOf` can match a read prefix to hap[0] while the read originally aligned
    // before the haplotype reference window — collapsing to hap ref start remaps interior bases.
    if original_ref_start_1based < hap_ref_start_1based
        && read_start_on_ref_1based >= hap_ref_start_1based
        && read_start_hap == 0
    {
        return Ok(false);
    }

    let final_cigar =
        append_clipped_elements_from_original(left.cigar, &leading_clips, &trailing_clips);

    let pos_0based = i64::try_from(read_start_on_ref_1based.saturating_sub(1)).unwrap_or(0);
    let old_pos = rec.pos();
    let old_cigar = format!("{:?}", rec.cigar().iter().collect::<Vec<_>>());
    let hts_cigar = cigar_to_hts(&final_cigar);
    let qname = rec.qname().to_vec();
    let seq = rec.seq().as_bytes().to_vec();
    let qual = rec.qual().to_vec();
    rec.set(&qname, Some(&hts_cigar), &seq, &qual);
    rec.set_pos(pos_0based);
    let new_cigar = format!("{:?}", rec.cigar().iter().collect::<Vec<_>>());
    Ok(old_pos != pos_0based || old_cigar != new_cigar)
}

fn clip_from_hts(hts: HtsCigar) -> Option<(usize, CigarOperator)> {
    match hts {
        HtsCigar::SoftClip(n) => Some((n as usize, CigarOperator::SoftClip)),
        HtsCigar::HardClip(n) => Some((n as usize, CigarOperator::HardClip)),
        _ => None,
    }
}

fn cigar_to_hts(cigar: &Cigar) -> CigarString {
    let v: Vec<HtsCigar> = cigar
        .elements
        .iter()
        .map(|e| match e.operator {
            CigarOperator::Match => HtsCigar::Match(e.length as u32),
            CigarOperator::Insertion => HtsCigar::Ins(e.length as u32),
            CigarOperator::Deletion => HtsCigar::Del(e.length as u32),
            CigarOperator::SoftClip => HtsCigar::SoftClip(e.length as u32),
            CigarOperator::HardClip => HtsCigar::HardClip(e.length as u32),
        })
        .collect();
    CigarString::from(v)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cigar::{Cigar, CigarOperator};
    use crate::genome_loc::GenomeLoc;
    use crate::smith_waterman::SwParameters;
    use rust_htslib::bam::record::Cigar as HtsCigar;
    use rust_htslib::bam::Record;

    fn test_record(seq: &[u8], pos_1based: u64) -> Record {
        let mut rec = Record::new();
        let cigar = CigarString::from(vec![HtsCigar::Match(seq.len() as u32)]);
        let qual: Vec<u8> = vec![30; seq.len()];
        rec.set(b"r1", Some(&cigar), seq, &qual);
        rec.set_pos(i64::try_from(pos_1based.saturating_sub(1)).unwrap_or(0));
        rec
    }

    #[test]
    fn change_evidence_preserves_full_read_hap_matrix() {
        let likelihoods = vec![
            RegionReadLikelihood {
                read_index: crate::bio_ids::ReadIndex::new(0),
                haplotype_index: crate::bio_ids::HaplotypeIndex::new(0),
                log10_likelihood: -1.0,
            },
            RegionReadLikelihood {
                read_index: crate::bio_ids::ReadIndex::new(0),
                haplotype_index: crate::bio_ids::HaplotypeIndex::new(1),
                log10_likelihood: -0.2,
            },
        ];
        let out = change_evidence_to_best_haplotype(likelihoods.clone(), &[1]);
        assert_eq!(out.len(), 2);
        assert_eq!(out, likelihoods);
    }

    #[test]
    fn tiebreak_prefers_reference_within_informative_threshold() {
        let mut ref_hap = Haplotype::new(b"ACGT", true);
        ref_hap.cigar = Some({
            let mut c = Cigar::new();
            c.push(4, CigarOperator::Match);
            c
        });
        let mut alt = Haplotype::new(b"ACGT", false);
        alt.cigar = Some({
            let mut c = Cigar::new();
            c.push(2, CigarOperator::Match);
            c.push(1, CigarOperator::Deletion);
            c.push(1, CigarOperator::Match);
            c
        });
        let haps = vec![ref_hap, alt];
        let priorities: Vec<f64> = haps
            .iter()
            .map(haplotype_alignment_tiebreak_priority)
            .collect();
        assert!(priorities[0] > priorities[1]);
        let likelihoods = vec![
            RegionReadLikelihood {
                read_index: crate::bio_ids::ReadIndex::new(0),
                haplotype_index: crate::bio_ids::HaplotypeIndex::new(0),
                log10_likelihood: -0.5,
            },
            RegionReadLikelihood {
                read_index: crate::bio_ids::ReadIndex::new(0),
                haplotype_index: crate::bio_ids::HaplotypeIndex::new(1),
                log10_likelihood: -0.48,
            },
        ];
        let best = best_haplotype_index_for_read(0, &haps, &likelihoods, &priorities);
        assert_eq!(best, 0, "ref wins tie within 0.2 log10");
    }

    #[test]
    fn create_read_aligned_to_ref_sets_position_on_alt_hap() {
        let sw = SwParameters::gatk_haplotype_to_reference();
        let pad = 1000u64;
        let mut ref_hap = Haplotype::new(b"ACGTACGTACGT", true);
        let mut ref_cigar = Cigar::new();
        ref_cigar.push(ref_hap.bases.len(), CigarOperator::Match);
        ref_hap.cigar = Some(ref_cigar);
        ref_hap.genome_loc = Some(GenomeLoc::new(pad, pad + 11));

        let mut alt_hap = Haplotype::new(b"ACGTACGTACGT", false);
        let mut alt_cigar = Cigar::new();
        alt_cigar.push(alt_hap.bases.len(), CigarOperator::Match);
        alt_hap.cigar = Some(alt_cigar);
        alt_hap.alignment_start_hap_wrt_ref = 0;

        let mut rec = test_record(b"ACGTACGTACGT", pad + 5);
        assert!(create_read_aligned_to_ref(&mut rec, &alt_hap, &ref_hap, pad, &sw,).unwrap());
        assert!(rec.pos() >= i64::try_from(pad.saturating_sub(1)).unwrap_or(0));
    }

    #[test]
    fn create_read_aligned_to_ref_skips_spurious_hap_prefix_match_before_hap_ref_window() {
        let sw = SwParameters::gatk_haplotype_to_reference();
        let pad = 1000u64;
        let mut ref_hap = Haplotype::new(b"ACGTACGTACGTACGTACGT", true);
        let mut ref_cigar = Cigar::new();
        ref_cigar.push(ref_hap.bases.len(), CigarOperator::Match);
        ref_hap.cigar = Some(ref_cigar);
        ref_hap.genome_loc = Some(GenomeLoc::new(pad, pad + 19));

        let mut alt_hap = Haplotype::new(b"ACGTACGTACGTACGTACGT", false);
        let mut alt_cigar = Cigar::new();
        alt_cigar.push(alt_hap.bases.len(), CigarOperator::Match);
        alt_hap.cigar = Some(alt_cigar);
        alt_hap.alignment_start_hap_wrt_ref = 40;

        let mut rec = test_record(b"ACGTACGTACGTACGTACGT", pad + 5);
        let orig_pos = rec.pos();
        assert!(!create_read_aligned_to_ref(&mut rec, &alt_hap, &ref_hap, pad, &sw,).unwrap());
        assert_eq!(rec.pos(), orig_pos);
    }

    #[test]
    fn create_read_aligned_to_ref_preserves_soft_clips() {
        let sw = SwParameters::gatk_haplotype_to_reference();
        let pad = 1000u64;
        let mut ref_hap = Haplotype::new(b"ACGTACGTACGTACGT", true);
        let mut ref_cigar = Cigar::new();
        ref_cigar.push(ref_hap.bases.len(), CigarOperator::Match);
        ref_hap.cigar = Some(ref_cigar);
        ref_hap.genome_loc = Some(GenomeLoc::new(pad, pad + 15));

        let mut alt_hap = Haplotype::new(b"ACGTACGTACGTACGT", false);
        let mut alt_cigar = Cigar::new();
        alt_cigar.push(alt_hap.bases.len(), CigarOperator::Match);
        alt_hap.cigar = Some(alt_cigar);
        alt_hap.alignment_start_hap_wrt_ref = 0;

        let mut rec = Record::new();
        let cigar = CigarString::from(vec![
            HtsCigar::SoftClip(2),
            HtsCigar::Match(8),
            HtsCigar::SoftClip(2),
        ]);
        let seq = b"NNACGTACGTNN";
        let qual: Vec<u8> = vec![30; seq.len()];
        rec.set(b"r1", Some(&cigar), seq, &qual);
        rec.set_pos(i64::try_from(pad + 4).unwrap_or(0));
        assert!(create_read_aligned_to_ref(&mut rec, &alt_hap, &ref_hap, pad, &sw,).unwrap());
        let cigar_view = rec.cigar();
        let out: Vec<_> = cigar_view.into_iter().collect();
        assert!(matches!(out.first(), Some(HtsCigar::SoftClip(2))));
        assert!(matches!(out.last(), Some(HtsCigar::SoftClip(2))));
    }
}
