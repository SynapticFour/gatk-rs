//! Genomic interval system for GATK-RS
//! This module provides efficient interval operations and indexing for genomic data,
//! supporting large-scale genomic analyses with optimal performance.

use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};
use std::ops::Range;
use std::sync::Arc;

/// Genomic interval with chromosome, coordinates, strand, and optional labels.
/// # Invariants
/// Coordinates are **0-based inclusive** (`start <= end` expected).
/// Same contig required for overlap/intersection operations.
/// # Ownership
/// Contig / id / name are `Arc<str>` so coordinate-only forks share string storage.
/// # Mutation
/// Public fields; prefer constructors for consistent strand defaults.
/// # Biological assumptions
/// Feature or callable region on reference; strand affects overlap semantics for directed features.
/// # Java equivalence
/// Approximates GATK `SimpleInterval` + strand (BED-like 0-based convention here).
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct GenomicInterval {
    /// Chromosome name (shared)
    pub chromosome: Arc<str>,

    /// Start position (0-based, inclusive)
    pub start: u64,

    /// End position (0-based, inclusive)
    pub end: u64,

    /// Interval strand
    pub strand: Strand,

    /// Optional interval ID (shared when present)
    pub id: Option<Arc<str>>,

    /// Optional interval name (shared when present)
    pub name: Option<Arc<str>>,
}

/// DNA strand orientation for interval features.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
pub enum Strand {
    Forward,
    Reverse,
}

impl Default for Strand {
    fn default() -> Self {
        Self::Forward
    }
}

impl GenomicInterval {
    /// Create a new genomic interval
    pub fn new(chromosome: impl Into<Arc<str>>, start: u64, end: u64) -> Self {
        Self {
            chromosome: chromosome.into(),
            start,
            end,
            strand: Strand::default(),
            id: None,
            name: None,
        }
    }

    /// Create a new genomic interval with strand
    pub fn new_with_strand(
        chromosome: impl Into<Arc<str>>,
        start: u64,
        end: u64,
        strand: Strand,
    ) -> Self {
        Self {
            chromosome: chromosome.into(),
            start,
            end,
            strand,
            id: None,
            name: None,
        }
    }

    /// Create a new genomic interval with ID and name
    pub fn new_with_metadata(
        chromosome: impl Into<Arc<str>>,
        start: u64,
        end: u64,
        strand: Strand,
        id: Option<Arc<str>>,
        name: Option<Arc<str>>,
    ) -> Self {
        Self {
            chromosome: chromosome.into(),
            start,
            end,
            strand,
            id,
            name,
        }
    }

    /// Coordinate/strand fork that bumps `Arc` refs instead of deep-copying strings.
    fn with_span(&self, start: u64, end: u64, strand: Strand) -> Self {
        Self {
            chromosome: Arc::clone(&self.chromosome),
            start,
            end,
            strand,
            id: self.id.as_ref().map(Arc::clone),
            name: self.name.as_ref().map(Arc::clone),
        }
    }

    /// Get the length of the interval
    pub fn length(&self) -> u64 {
        self.end - self.start + 1
    }

    /// Check if interval overlaps with another interval
    pub fn overlaps(&self, other: &GenomicInterval) -> bool {
        self.chromosome == other.chromosome && self.start <= other.end && self.end >= other.start
    }

    /// Check if position is within interval
    pub fn contains(&self, position: u64) -> bool {
        position >= self.start && position <= self.end
    }

    /// Check if interval is valid (start <= end)
    pub fn is_valid(&self) -> bool {
        self.start <= self.end
    }

    /// Half-open range covering this inclusive interval (`start..end+1`).
    pub fn range(&self) -> Range<u64> {
        self.start..self.end.saturating_add(1)
    }

    /// Intersect with another interval
    pub fn intersection(&self, other: &GenomicInterval) -> Option<GenomicInterval> {
        if !self.overlaps(other) {
            return None;
        }

        let start = std::cmp::max(self.start, other.start);
        let end = std::cmp::min(self.end, other.end);

        if start <= end {
            let strand = if self.strand == other.strand {
                self.strand
            } else {
                other.strand
            };
            Some(self.with_span(start, end, strand))
        } else {
            None
        }
    }

    /// Merge with another interval
    pub fn union(&self, other: &GenomicInterval) -> Vec<GenomicInterval> {
        if self.chromosome != other.chromosome {
            return vec![
                self.with_span(self.start, self.end, self.strand),
                other.with_span(other.start, other.end, other.strand),
            ];
        }

        if self.overlaps(other) {
            let start = std::cmp::min(self.start, other.start);
            let end = std::cmp::max(self.end, other.end);
            let strand = if self.strand == other.strand {
                self.strand
            } else {
                other.strand
            };
            vec![self.with_span(start, end, strand)]
        } else {
            vec![
                self.with_span(self.start, self.end, self.strand),
                other.with_span(other.start, other.end, other.strand),
            ]
        }
    }

    /// Get the center position of the interval
    pub fn center(&self) -> u64 {
        (self.start + self.end) / 2
    }

    /// Expand interval by given amount on both sides
    pub fn expand(&self, amount: u64) -> GenomicInterval {
        self.with_span(
            self.start.saturating_sub(amount),
            self.end.saturating_add(amount),
            self.strand,
        )
    }

    /// Shrink interval by given amount on both sides
    pub fn shrink(&self, amount: u64) -> GenomicInterval {
        self.with_span(
            self.start.saturating_add(amount),
            self.end.saturating_sub(amount),
            self.strand,
        )
    }

    /// Pad interval to include given position
    pub fn pad_to_include(&self, position: u64) -> GenomicInterval {
        if self.contains(position) {
            self.with_span(self.start, self.end, self.strand)
        } else if position < self.start {
            self.with_span(position, self.end, self.strand)
        } else {
            self.with_span(self.start, position, self.strand)
        }
    }
}

/// Per-chromosome interval index for overlap and intersection queries.
#[derive(Debug, Clone)]
pub struct IntervalTree {
    intervals: BTreeMap<Arc<str>, Vec<GenomicInterval>>,
}

impl IntervalTree {
    /// Create a new interval tree
    pub fn new() -> Self {
        Self {
            intervals: BTreeMap::new(),
        }
    }

    /// Add intervals to the tree
    pub fn add_intervals(&mut self, intervals: Vec<GenomicInterval>) {
        for interval in intervals {
            let key = Arc::clone(&interval.chromosome);
            self.intervals.entry(key).or_default().push(interval);
        }
    }

    /// Find all intervals overlapping a position
    pub fn find_overlapping(&self, chromosome: &str, position: u64) -> Vec<&GenomicInterval> {
        self.intervals
            .get(chromosome)
            .map(|intervals| {
                intervals
                    .iter()
                    .filter(|interval| interval.contains(position))
                    .collect()
            })
            .unwrap_or_default()
    }

    /// Find intervals intersecting a range
    pub fn find_intersecting(&self, chromosome: &str, range: Range<u64>) -> Vec<&GenomicInterval> {
        self.intervals
            .get(chromosome)
            .map(|intervals| {
                intervals
                    .iter()
                    .filter(|interval| interval.start < range.end && interval.end >= range.start)
                    .collect()
            })
            .unwrap_or_default()
    }

    /// Count total intervals stored
    pub fn count(&self) -> usize {
        self.intervals.values().map(|v| v.len()).sum()
    }

    /// Get all intervals for a chromosome
    pub fn get_chromosome_intervals(&self, chromosome: &str) -> &[GenomicInterval] {
        self.intervals
            .get(chromosome)
            .map(|v| v.as_slice())
            .unwrap_or(&[])
    }

    /// Clear all intervals
    pub fn clear(&mut self) {
        self.intervals.clear();
    }
}

impl Default for IntervalTree {
    fn default() -> Self {
        Self::new()
    }
}

/// Query result bucket for overlapping and fully contained intervals.
/// Borrows intervals from the backing tree — no owned interval copies.
#[derive(Debug, Clone)]
pub struct IntervalQuery<'a> {
    pub overlapping: Vec<&'a GenomicInterval>,
    pub contained: Vec<&'a GenomicInterval>,
    pub count: usize,
}

impl<'a> IntervalQuery<'a> {
    /// Create a new interval query result
    pub fn new() -> Self {
        Self {
            overlapping: Vec::new(),
            contained: Vec::new(),
            count: 0,
        }
    }

    /// Add an overlapping interval
    pub fn add_overlapping(&mut self, interval: &'a GenomicInterval) {
        self.overlapping.push(interval);
        self.count += 1;
    }

    /// Add a contained interval
    pub fn add_contained(&mut self, interval: &'a GenomicInterval) {
        self.contained.push(interval);
        self.count += 1;
    }

    /// Check if any intervals found
    pub fn has_results(&self) -> bool {
        self.count > 0
    }
}

impl<'a> Default for IntervalQuery<'a> {
    fn default() -> Self {
        Self::new()
    }
}

/// Facade over [`IntervalTree`] for position/range/containment queries.
#[derive(Debug, Clone)]
pub struct IntervalQueryEngine {
    tree: IntervalTree,
}

impl IntervalQueryEngine {
    /// Create a new interval query engine
    pub fn new() -> Self {
        Self {
            tree: IntervalTree::new(),
        }
    }

    /// Build interval tree from intervals
    pub fn build_tree(&mut self, intervals: Vec<GenomicInterval>) {
        self.tree.add_intervals(intervals);
    }

    /// Query intervals overlapping a position
    pub fn query_position(&self, chromosome: &str, position: u64) -> IntervalQuery<'_> {
        let overlapping = self.tree.find_overlapping(chromosome, position);
        let count = overlapping.len();
        IntervalQuery {
            overlapping,
            contained: Vec::new(),
            count,
        }
    }

    /// Query intervals intersecting a range
    pub fn query_range(&self, chromosome: &str, range: Range<u64>) -> IntervalQuery<'_> {
        let intersecting = self.tree.find_intersecting(chromosome, range);
        let count = intersecting.len();
        IntervalQuery {
            overlapping: intersecting,
            contained: Vec::new(),
            count,
        }
    }

    /// Query intervals containing a position
    pub fn query_contained(&self, chromosome: &str, position: u64) -> IntervalQuery<'_> {
        let contained = self.tree.find_overlapping(chromosome, position);
        let count = contained.len();
        IntervalQuery {
            overlapping: Vec::new(),
            contained,
            count,
        }
    }

    /// Get statistics for the interval tree
    pub fn get_statistics(&self) -> IntervalStatistics {
        let total_intervals = self.tree.count();
        let total_chromosomes = self.tree.intervals.len();

        let mut total_length = 0u64;
        let mut min_length = u64::MAX;
        let mut max_length = 0u64;

        for intervals in self.tree.intervals.values() {
            for interval in intervals {
                let length = interval.length();
                total_length += length;
                min_length = min_length.min(length);
                max_length = max_length.max(length);
            }
        }

        IntervalStatistics {
            total_intervals,
            total_chromosomes,
            total_length,
            min_length: if total_intervals > 0 { min_length } else { 0 },
            max_length,
            average_length: if total_intervals > 0 {
                total_length / total_intervals as u64
            } else {
                0
            },
        }
    }
}

impl Default for IntervalQueryEngine {
    fn default() -> Self {
        Self::new()
    }
}

/// Aggregate statistics over an indexed interval collection.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct IntervalStatistics {
    pub total_intervals: usize,
    pub total_chromosomes: usize,
    pub total_length: u64,
    pub min_length: u64,
    pub max_length: u64,
    pub average_length: u64,
}

/// Ordered set of unique [`GenomicInterval`] values.
#[derive(Debug, Clone)]
pub struct IntervalSet {
    intervals: BTreeSet<GenomicInterval>,
}

impl IntervalSet {
    /// Create a new interval set
    pub fn new() -> Self {
        Self {
            intervals: BTreeSet::new(),
        }
    }

    /// Add interval to the set
    pub fn insert(&mut self, interval: GenomicInterval) {
        self.intervals.insert(interval);
    }

    /// Check if interval is in the set
    pub fn contains(&self, interval: &GenomicInterval) -> bool {
        self.intervals.contains(interval)
    }

    /// Get all intervals in the set
    pub fn get_intervals(&self) -> Vec<&GenomicInterval> {
        self.intervals.iter().collect()
    }

    /// Clear all intervals
    pub fn clear(&mut self) {
        self.intervals.clear();
    }

    /// Get count of intervals
    pub fn len(&self) -> usize {
        self.intervals.len()
    }

    /// True when the set is empty
    pub fn is_empty(&self) -> bool {
        self.intervals.is_empty()
    }
}

impl Default for IntervalSet {
    fn default() -> Self {
        Self::new()
    }
}

/// Parallel interval overlap query.
/// Results are merged into a [`BTreeMap`] keyed by `(chromosome, position)` so
/// return order is deterministic regardless of Rayon thread count.
pub fn find_overlapping_intervals_parallel<'a>(
    intervals: &'a [GenomicInterval],
    queries: &[(String, u64)],
) -> BTreeMap<(String, u64), Vec<&'a GenomicInterval>> {
    let mut rows: Vec<_> = queries
        .par_iter()
        .map(|(chromosome, position)| {
            let overlapping: Vec<&'a GenomicInterval> = intervals
                .iter()
                .filter(|interval| {
                    interval.chromosome.as_ref() == chromosome.as_str()
                        && interval.contains(*position)
                })
                .collect();
            ((chromosome.to_owned(), *position), overlapping)
        })
        .collect();
    rows.sort_by(|a, b| a.0.cmp(&b.0));
    rows.into_iter().collect()
}

/// Merge overlapping intervals (sorted by start; sequential merge for determinism).
pub fn merge_overlapping_intervals_parallel(intervals: &[GenomicInterval]) -> Vec<GenomicInterval> {
    let mut sorted_intervals = intervals.to_vec();
    sorted_intervals
        .sort_by(|a, b| (&*a.chromosome, a.start, a.end).cmp(&(&*b.chromosome, b.start, b.end)));

    let mut merged: Vec<GenomicInterval> = Vec::new();
    for interval in sorted_intervals {
        let can_merge = merged
            .last()
            .is_some_and(|last| last.chromosome == interval.chromosome && last.overlaps(&interval));
        if can_merge {
            if let Some(last) = merged.pop() {
                let mut union_result = last.union(&interval);
                if union_result.len() == 1 {
                    if let Some(u) = union_result.pop() {
                        merged.push(u);
                    }
                } else {
                    merged.push(last);
                    merged.push(interval);
                }
            }
        } else {
            merged.push(interval);
        }
    }

    merged
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn arc_share_on_expand_keeps_same_contig_ptr() {
        let a = GenomicInterval::new("chr1", 10, 20);
        let b = a.expand(5);
        assert!(Arc::ptr_eq(&a.chromosome, &b.chromosome));
        assert_eq!(b.start, 5);
        assert_eq!(b.end, 25);
    }

    #[test]
    fn query_engine_borrows_without_owned_copies() {
        let mut engine = IntervalQueryEngine::new();
        engine.build_tree(vec![
            GenomicInterval::new("chr1", 0, 100),
            GenomicInterval::new("chr1", 50, 150),
        ]);
        let q = engine.query_position("chr1", 75);
        assert_eq!(q.count, 2);
        assert_eq!(q.overlapping.len(), 2);
    }

    #[test]
    fn merge_overlapping_same_contig() {
        let merged = merge_overlapping_intervals_parallel(&[
            GenomicInterval::new("chr1", 0, 10),
            GenomicInterval::new("chr1", 5, 20),
            GenomicInterval::new("chr1", 30, 40),
        ]);
        assert_eq!(merged.len(), 2);
        assert_eq!((merged[0].start, merged[0].end), (0, 20));
        assert_eq!((merged[1].start, merged[1].end), (30, 40));
    }
}
