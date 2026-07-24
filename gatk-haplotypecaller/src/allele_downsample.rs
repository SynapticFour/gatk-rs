//! GATK `AlleleBiasedDownsamplingUtils`.
//! Removes evidence to mitigate contamination bias before genotype likelihoods / `isActive`.

use crate::activity_scoring::PileupObservation;
use crate::gatk_well_rng::Well19937c;
use std::collections::BTreeMap;

/// GATK `AlleleBiasedDownsamplingUtils.targetAlleleCounts`.
#[must_use]
pub fn target_allele_counts(allele_counts: &[usize], num_reads_to_remove: usize) -> Vec<usize> {
    let num_alleles = allele_counts.len();
    if num_alleles == 0 {
        return Vec::new();
    }

    let mut max_score = score_allele_counts(allele_counts);
    let mut allele_counts_of_max = allele_counts.to_vec();

    let num_reads_to_remove_per_allele = num_reads_to_remove / 2;

    for i in 0..num_alleles {
        for j in i..num_alleles {
            let mut new_counts = allele_counts.to_vec();
            if i == j {
                new_counts[i] = new_counts[i].saturating_sub(num_reads_to_remove);
            } else {
                new_counts[i] = new_counts[i].saturating_sub(num_reads_to_remove_per_allele);
                new_counts[j] = new_counts[j].saturating_sub(num_reads_to_remove_per_allele);
            }
            let score = score_allele_counts(&new_counts);
            if score < max_score {
                max_score = score;
                allele_counts_of_max = new_counts;
            }
        }
    }

    allele_counts_of_max
}

fn score_allele_counts(allele_counts: &[usize]) -> i64 {
    if allele_counts.len() < 2 {
        return 0;
    }
    let mut sorted = allele_counts.to_vec();
    sorted.sort_unstable();
    let max_count = sorted[sorted.len() - 1];
    let next_best_count = sorted[sorted.len() - 2];
    let remainder_count: usize = sorted.iter().sum::<usize>() - max_count - next_best_count;

    let term_a = max_count as i64 - next_best_count as i64 + remainder_count as i64;
    let term_b = (next_best_count as i64 + remainder_count as i64).abs();
    term_a.min(term_b)
}

/// Indices into `evidence` lists to **remove** (GATK `selectAlleleBiasedEvidence`).
pub fn select_allele_biased_evidence_indices(
    allele_evidence: &BTreeMap<u8, Vec<usize>>,
    contamination_fraction: f64,
    rng: &mut Well19937c,
) -> Vec<usize> {
    let total: usize = allele_evidence.values().map(|v| v.len()).sum();
    if total == 0 || contamination_fraction <= 0.0 {
        return Vec::new();
    }
    let num_to_remove = (total as f64 * contamination_fraction) as usize;
    if num_to_remove == 0 {
        return Vec::new();
    }

    let alleles: Vec<u8> = allele_evidence.keys().copied().collect();
    let allele_counts: Vec<usize> = alleles
        .iter()
        .map(|a| allele_evidence.get(a).map_or(0, |v| v.len()))
        .collect();
    let target = target_allele_counts(&allele_counts, num_to_remove);

    let mut remove = Vec::new();
    for (i, allele) in alleles.iter().enumerate() {
        let current = allele_counts[i];
        let target_n = target[i];
        if current > target_n {
            let list = allele_evidence.get(allele).expect("allele key");
            remove.extend(downsample_element_indices(list, current - target_n, rng));
        }
    }
    remove
}

fn downsample_element_indices(
    evidence: &[usize],
    num_elements_to_remove: usize,
    rng: &mut Well19937c,
) -> Vec<usize> {
    if num_elements_to_remove == 0 {
        return Vec::new();
    }
    if num_elements_to_remove >= evidence.len() {
        return evidence.to_vec();
    }
    rng.sample_indices_without_replacement(evidence.len(), num_elements_to_remove)
        .into_iter()
        .map(|i| evidence[i])
        .collect()
}

/// Apply contamination downsampling to a single-sample pileup (in-place).
pub fn apply_contamination_to_pileup(
    pile: &mut Vec<PileupObservation>,
    contamination_fraction: f64,
    rng: &mut Well19937c,
) {
    if pile.is_empty() || contamination_fraction <= 0.0 {
        return;
    }
    let mut by_allele: BTreeMap<u8, Vec<usize>> = BTreeMap::new();
    for (idx, obs) in pile.iter().enumerate() {
        let base = obs.read_base.to_ascii_uppercase();
        if base == b'N' {
            continue;
        }
        by_allele.entry(base).or_default().push(idx);
    }
    let mut to_remove =
        select_allele_biased_evidence_indices(&by_allele, contamination_fraction, rng);
    to_remove.sort_unstable();
    to_remove.dedup();
    for &idx in to_remove.iter().rev() {
        if idx < pile.len() {
            pile.remove(idx);
        }
    }
}

/// Dump `targetAlleleCounts` parity row (/ gate `d2c`).
#[cfg(any(feature = "dev-dumps", test))]
pub fn dump_target_allele_counts_tsv(
    counts_csv: &str,
    num_reads_to_remove: usize,
    out: &mut impl std::io::Write,
) -> Result<(), String> {
    let counts: Vec<usize> = counts_csv
        .split(',')
        .map(|s| s.trim().parse().map_err(|_| format!("invalid count: {s}")))
        .collect::<Result<_, _>>()?;
    let target = target_allele_counts(&counts, num_reads_to_remove);
    writeln!(out, "allele_counts\tnum_reads_to_remove\ttarget_counts")
        .map_err(|e| format!("write: {e}"))?;
    writeln!(
        out,
        "{}\t{}\t{}",
        counts_csv,
        num_reads_to_remove,
        target
            .iter()
            .map(|c| c.to_string())
            .collect::<Vec<_>>()
            .join(",")
    )
    .map_err(|e| format!("write: {e}"))?;
    Ok(())
}

/// Dump evidence removal at one locus: `removed_count`, sorted `removed_qnames`.
#[cfg(any(feature = "dev-dumps", test))]
pub fn dump_allele_biased_evidence_locus_tsv(
    reference_fasta: &std::path::Path,
    alignment_path: &std::path::Path,
    contig: &str,
    pos1: u64,
    contamination_fraction: f64,
    out: &mut impl std::io::Write,
    read_filters: &crate::read_model::ReadFilterParams,
) -> Result<(), String> {
    use crate::read_binding::record_overlaps_closed_interval_1based;
    use crate::read_model::passes_hc_read_filters_with_header;
    use crate::read_projection::query_index_at_reference_position;
    use gatk_core::reference::{ReferenceWindowCache, SequenceDictionary};
    use rust_htslib::bam::Read as _;

    let dict = SequenceDictionary::from_fasta_path(reference_fasta).map_err(|e| e.to_string())?;
    let mut ref_cache = ReferenceWindowCache::new(reference_fasta.to_path_buf(), 4);
    let ref_byte = *ref_cache
        .get_interval_bytes(&dict, contig, pos1, pos1)
        .map_err(|e| e.to_string())?
        .first()
        .unwrap_or(&b'N');

    let mut reader =
        rust_htslib::bam::Reader::from_path(alignment_path).map_err(|e| format!("open: {e}"))?;
    let header = reader.header().clone();
    let tid = header
        .tid(contig.as_bytes())
        .ok_or_else(|| format!("missing contig {contig}"))? as i32;
    let mut records = Vec::new();
    for res in reader.records() {
        let rec = res.map_err(|e| format!("read: {e}"))?;
        if rec.tid() == tid {
            records.push(rec);
        }
    }

    let ref_pos0 = pos1.saturating_sub(1) as i64;
    let mut by_allele: BTreeMap<u8, Vec<usize>> = BTreeMap::new();
    for (i, rec) in records.iter().enumerate() {
        if !passes_hc_read_filters_with_header(rec, &header, read_filters) {
            continue;
        }
        if !record_overlaps_closed_interval_1based(rec, &header, contig, pos1, pos1, read_filters) {
            continue;
        }
        let cigar = rec.cigar().to_owned();
        let qi = match query_index_at_reference_position(rec.pos(), &cigar, ref_pos0) {
            Some(q) => q,
            None => continue,
        };
        let seq = rec.seq().as_bytes();
        if qi >= seq.len() {
            continue;
        }
        let read_b = seq[qi].to_ascii_uppercase();
        let _ = ref_byte;
        by_allele.entry(read_b).or_default().push(i);
    }

    let mut rng = Well19937c::reset_gatk_default();
    let removed_idx =
        select_allele_biased_evidence_indices(&by_allele, contamination_fraction, &mut rng);
    let mut names: Vec<String> = removed_idx
        .iter()
        .map(|&i| String::from_utf8_lossy(records[i].qname()).into_owned())
        .collect();
    names.sort();

    writeln!(
        out,
        "contig\tpos\tcontamination_fraction\tremoved_count\tremoved_qnames"
    )
    .map_err(|e| format!("write: {e}"))?;
    writeln!(
        out,
        "{contig}\t{pos1}\t{contamination_fraction}\t{}\t{}",
        names.len(),
        names.join(",")
    )
    .map_err(|e| format!("write: {e}"))?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn target_counts_match_gatk_unit_test_hom_contaminant() {
        let ideal_hom = [0usize, 100, 0, 0];
        let actual = [10usize, 100, 0, 0];
        let target = target_allele_counts(&actual, 10);
        assert_eq!(target, ideal_hom);
    }

    #[test]
    fn target_counts_het_overlapping_unchanged() {
        let actual = [0, 55, 0, 55];
        let target = target_allele_counts(&actual, 10);
        assert_eq!(target, [0, 55, 0, 55]);
    }

    #[test]
    fn select_evidence_removal_count() {
        let mut map: BTreeMap<u8, Vec<usize>> = BTreeMap::new();
        let mut idx = 0usize;
        for _ in 0..10 {
            map.entry(b'A').or_default().push(idx);
            idx += 1;
        }
        for _ in 0..100 {
            map.entry(b'C').or_default().push(idx);
            idx += 1;
        }
        let mut rng = Well19937c::reset_gatk_default();
        let removed = select_allele_biased_evidence_indices(&map, 0.1, &mut rng);
        assert_eq!(removed.len(), 10);
    }
}
