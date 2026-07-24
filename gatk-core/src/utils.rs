//! Utility functions for genomic analysis

use crate::types::*;
use std::collections::{HashMap, HashSet};

/// Reverse complement a DNA sequence
pub fn reverse_complement(sequence: &[Base]) -> Vec<Base> {
    sequence
        .iter()
        .rev()
        .map(|base| base.complement())
        .collect()
}

/// Calculate GC content of a sequence
pub fn gc_content(sequence: &[Base]) -> f64 {
    if sequence.is_empty() {
        return 0.0;
    }

    let gc_count = sequence
        .iter()
        .filter(|&&base| matches!(base, Base::G | Base::C))
        .count();

    gc_count as f64 / sequence.len() as f64
}

/// Check if a position is within any of the given intervals
pub fn position_in_intervals(pos: GenomicPosition, intervals: &[GenomicInterval]) -> bool {
    intervals.iter().any(|interval| interval.contains(pos))
}

/// Merge overlapping intervals
pub fn merge_intervals(intervals: &[GenomicInterval]) -> Vec<GenomicInterval> {
    if intervals.is_empty() {
        return Vec::new();
    }

    let mut sorted_intervals = intervals.to_vec();
    sorted_intervals.sort_by(|a, b| {
        a.contig
            .cmp(&b.contig)
            .then(a.start.cmp(&b.start))
            .then(a.end.cmp(&b.end))
    });

    // Lifetime: `sorted_intervals` is function-local; move intervals into `current`
    // / `merged` instead of cloning while scanning.
    let mut merged = Vec::new();
    let mut iter = sorted_intervals.into_iter();
    let mut current = iter.next().expect("non-empty");

    for interval in iter {
        if interval.contig == current.contig && interval.start <= current.end + 1 {
            // Overlapping or adjacent intervals - merge them
            current.end = current.end.max(interval.end);
        } else {
            // Non-overlapping interval - push current and start new one
            merged.push(current);
            current = interval;
        }
    }
    merged.push(current);

    merged
}

/// Calculate Hamming distance between two sequences
pub fn hamming_distance(seq1: &[Base], seq2: &[Base]) -> Option<usize> {
    if seq1.len() != seq2.len() {
        return None;
    }

    Some(seq1.iter().zip(seq2.iter()).filter(|(a, b)| a != b).count())
}

/// Find all unique bases in a sequence
pub fn unique_bases(sequence: &[Base]) -> HashSet<Base> {
    sequence.iter().cloned().collect()
}

/// Check if a sequence contains only valid DNA bases (A, C, G, T)
pub fn is_valid_dna(sequence: &[Base]) -> bool {
    sequence
        .iter()
        .all(|&base| matches!(base, Base::A | Base::C | Base::G | Base::T))
}

/// Calculate base composition of a sequence
pub fn base_composition(sequence: &[Base]) -> HashMap<Base, usize> {
    let mut composition = HashMap::new();
    for &base in sequence {
        *composition.entry(base).or_insert(0) += 1;
    }
    composition
}

/// Convert Phred quality score to probability
pub fn phred_to_probability(quality: BaseQuality) -> f64 {
    quality.error_probability()
}

/// Convert probability to Phred quality score
pub fn probability_to_phred(prob: f64) -> BaseQuality {
    let quality = if prob <= 0.0 {
        93.0 // Maximum quality
    } else if prob >= 1.0 {
        0.0
    } else {
        (-10.0 * prob.log10()).clamp(0.0, 93.0)
    };
    BaseQuality::new(quality as u8)
}

/// Calculate average base quality
pub fn average_base_quality(qualities: &[BaseQuality]) -> f64 {
    if qualities.is_empty() {
        return 0.0;
    }

    let sum: f64 = qualities.iter().map(|q| q.value() as f64).sum();
    sum / qualities.len() as f64
}

/// Filter bases by minimum quality
pub fn filter_by_quality(
    sequence: &[Base],
    qualities: &[BaseQuality],
    min_quality: u8,
) -> Vec<(Base, BaseQuality)> {
    sequence
        .iter()
        .zip(qualities.iter())
        .filter_map(|(&base, &quality)| {
            if quality.value() >= min_quality {
                Some((base, quality))
            } else {
                None
            }
        })
        .collect()
}

/// Check if two alleles are the same
pub fn alleles_equal(allele1: &Allele, allele2: &Allele) -> bool {
    allele1.bases == allele2.bases
}

/// Create a reference allele (single base)
pub fn ref_allele(base: Base) -> Allele {
    Allele::new(vec![base])
}

/// Create an alternative allele
pub fn alt_allele(bases: Vec<Base>) -> Allele {
    Allele::new(bases)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_reverse_complement() {
        let seq = vec![Base::A, Base::C, Base::G, Base::T];
        let rc = reverse_complement(&seq);
        assert_eq!(rc, vec![Base::A, Base::C, Base::G, Base::T]); // ATGC -> GCAT -> reverse complement

        let seq = vec![Base::A, Base::T, Base::G, Base::C];
        let rc = reverse_complement(&seq);
        assert_eq!(rc, vec![Base::G, Base::C, Base::A, Base::T]);
    }

    #[test]
    fn test_gc_content() {
        let seq = vec![Base::G, Base::C, Base::A, Base::T];
        assert_eq!(gc_content(&seq), 0.5);

        let seq = vec![Base::G, Base::C, Base::G, Base::C];
        assert_eq!(gc_content(&seq), 1.0);

        let seq = vec![Base::A, Base::T, Base::A, Base::T];
        assert_eq!(gc_content(&seq), 0.0);
    }

    #[test]
    fn test_hamming_distance() {
        let seq1 = vec![Base::A, Base::C, Base::G, Base::T];
        let seq2 = vec![Base::A, Base::C, Base::A, Base::T];
        assert_eq!(hamming_distance(&seq1, &seq2), Some(1));

        let seq3 = vec![Base::A, Base::C];
        assert_eq!(hamming_distance(&seq1, &seq3), None);
    }

    #[test]
    fn test_merge_intervals() {
        let intervals = vec![
            GenomicInterval::new(0, 100, 200),
            GenomicInterval::new(0, 150, 250),
            GenomicInterval::new(0, 300, 400),
        ];

        let merged = merge_intervals(&intervals);
        assert_eq!(
            merged,
            vec![
                GenomicInterval::new(0, 100, 250),
                GenomicInterval::new(0, 300, 400),
            ]
        );
    }

    #[test]
    fn test_phred_conversion() {
        let q = BaseQuality::new(20);
        let prob = phred_to_probability(q);
        assert!((prob - 0.01).abs() < 0.001);

        let q_back = probability_to_phred(prob);
        assert_eq!(q_back.value(), 20);
    }
}
