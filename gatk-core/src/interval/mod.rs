//! Genomic interval system for GATK-RS
//! This module provides efficient interval operations and indexing for genomic data,
//! supporting large-scale genomic analyses with optimal performance.

use std::collections::{BTreeMap, BTreeSet};
use std::ops::{Range, RangeBounds};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

/// Genomic interval with chromosome, coordinates, strand, and optional labels.
/// # Invariants
/// Coordinates are **0-based inclusive** (`start <= end` expected).
/// Same contig required for overlap/intersection operations.
/// # Ownership
/// Owns chromosome/id/name strings; clone for collections and trees.
/// # Mutation
/// Public fields; prefer constructors for consistent strand defaults.
/// # Biological assumptions
/// Feature or callable region on reference; strand affects overlap semantics for directed features.
/// # Java equivalence
/// Approximates GATK `SimpleInterval` + strand (BED-like 0-based convention here).
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
pub struct GenomicInterval {
    /// Chromosome name
    pub chromosome: String,
    
    /// Start position (0-based, inclusive)
    pub start: u64,
    
    /// End position (0-based, inclusive)
    pub end: u64,
    
    /// Interval strand
    pub strand: Strand,
    
    /// Optional interval ID
    pub id: Option<String>,
    
    /// Optional interval name
    pub name: Option<String>,
}

/// DNA strand orientation for interval features.
/// # Invariants
/// Only forward or reverse enumerated.
/// # Ownership
/// `Copy` enum.
/// # Mutation
/// N/A.
/// # Biological assumptions
/// Forward = reference strand direction; reverse = opposite strand.
/// # Java equivalence
/// Similar to htsjdk `Strand` / GATK `Strand` (conceptual).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
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
    pub fn new(chromosome: String, start: u64, end: u64) -> Self {
        Self {
            chromosome,
            start,
            end,
            strand: Strand::default(),
            id: None,
            name: None,
        }
    }

    /// Create a new genomic interval with strand
    pub fn new_with_strand(
        chromosome: String,
        start: u64,
        end: u64,
        strand: Strand,
    ) -> Self {
        Self {
            chromosome,
            start,
            end,
            strand,
            id: None,
            name: None,
        }
    }

    /// Create a new genomic interval with ID and name
    pub fn new_with_metadata(
        chromosome: String,
        start: u64,
        end: u64,
        strand: Strand,
        id: Option<String>,
        name: Option<String>,
    ) -> Self {
        Self {
            chromosome,
            start,
            end,
            strand,
            id,
            name,
        }
    }

    /// Get the length of the interval
    pub fn length(&self) -> u64 {
        self.end - self.start + 1
    }

    /// Check if interval overlaps with another interval
    pub fn overlaps(&self, other: &GenomicInterval) -> bool {
        self.chromosome == other.chromosome && 
        self.start <= other.end && 
        self.end >= other.start
    }

    /// Check if position is within interval
    pub fn contains(&self, position: u64) -> bool {
        position >= self.start && position <= self.end
    }

    /// Check if interval is valid (start <= end)
    pub fn is_valid(&self) -> bool {
        self.start <= self.end
    }

    /// Get the range of the interval
    pub fn range(&self) -> Range<u64> {
        self.start..=self.end + 1
    }

    /// Intersect with another interval
    pub fn intersection(&self, other: &GenomicInterval) -> Option<GenomicInterval> {
        if !self.overlaps(other) || self.chromosome != other.chromosome {
            return None;
        }

        let start = std::cmp::max(self.start, other.start);
        let end = std::cmp::min(self.end, other.end);

        if start <= end {
            Some(GenomicInterval::new_with_strand(
                self.chromosome.clone(),
                start,
                end,
                if self.strand == other.strand { self.strand } else { other.strand },
                self.id.clone(),
                self.name.clone(),
            ))
        } else {
            None
        }
    }

    /// Merge with another interval
    pub fn union(&self, other: &GenomicInterval) -> Vec<GenomicInterval> {
        if self.chromosome != other.chromosome {
            // CLONE: needed because graph fork needs owned duplicate for speculative path.
            return vec![self.clone(), other.clone()];
        }

        if self.overlaps(other) {
            // Overlapping intervals - merge them
            let start = std::cmp::min(self.start, other.start);
            let end = std::cmp::max(self.end, other.end);
            vec![GenomicInterval::new_with_strand(
                self.chromosome.clone(),
                start,
                end,
                if self.strand == other.strand { self.strand } else { other.strand },
                self.id.clone(),
                self.name.clone(),
            )]
        } else {
            // Non-overlapping intervals - return both
            // CLONE: needed because graph fork needs owned duplicate for speculative path.
            vec![self.clone(), other.clone()]
        }
    }

    /// Get the center position of the interval
    pub fn center(&self) -> u64 {
        (self.start + self.end) / 2
    }

    /// Expand interval by given amount on both sides
    pub fn expand(&self, amount: u64) -> GenomicInterval {
        GenomicInterval::new_with_strand(
            self.chromosome.clone(),
            self.start.saturating_sub(amount),
            self.end.saturating_add(amount),
            self.strand,
            self.id.clone(),
            self.name.clone(),
        )
    }

    /// Shrink interval by given amount on both sides
    pub fn shrink(&self, amount: u64) -> GenomicInterval {
        GenomicInterval::new_with_strand(
            self.chromosome.clone(),
            self.start.saturating_add(amount),
            self.end.saturating_sub(amount),
            self.strand,
            self.id.clone(),
            self.name.clone(),
        )
    }

    /// Pad interval to include given position
    pub fn pad_to_include(&self, position: u64) -> GenomicInterval {
        if self.contains(position) {
            // CLONE: needed because graph fork needs owned duplicate for speculative path.
            self.clone()
        } else if position < self.start {
            GenomicInterval::new_with_strand(
                self.chromosome.clone(),
                position,
                self.end,
                self.strand,
                self.id.clone(),
                self.name.clone(),
            )
        } else if position > self.end {
            GenomicInterval::new_with_strand(
                self.chromosome.clone(),
                self.start,
                position,
                self.strand,
                self.id.clone(),
                self.name.clone(),
            )
        } else {
            // CLONE: needed because graph fork needs owned duplicate for speculative path.
            self.clone()
        }
    }
}

/// Per-chromosome interval index for overlap and intersection queries.
/// # Invariants
/// Intervals grouped by `chromosome` key; query methods use 0-based coordinates.
/// # Ownership
/// Owns nested interval vectors; queries return borrowed `&GenomicInterval`.
/// # Mutation
/// `add_intervals` appends; no automatic merge/dedup.
/// # Biological assumptions
/// Stores callable regions, targets, or feature intervals for random access.
/// # Java equivalence
/// Approximates GATK `GenomeLocTree` / interval roster utilities (Rust-native layout).
#[derive(Debug, Clone)]
pub struct IntervalTree {
    intervals: BTreeMap<String, Vec<GenomicInterval>>,
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
            self.intervals
                // CLONE: needed because owned HashMap entry key.
                .entry(interval.chromosome.clone())
                .or_insert_with(Vec::new)
                .push(interval);
        }
    }

    /// Find all intervals overlapping a position
    pub fn find_overlapping(&self, chromosome: &str, position: u64) -> Vec<&GenomicInterval> {
        self.intervals
            .get(chromosome)
            .map_or(&vec![], |intervals| {
                intervals.iter()
                    .filter(|interval| interval.contains(position))
                    .collect()
            })
    }

    /// Find intervals intersecting a range
    pub fn find_intersecting(&self, chromosome: &str, range: Range<u64>) -> Vec<&GenomicInterval> {
        self.intervals
            .get(chromosome)
            .map_or(&vec![], |intervals| {
                intervals.iter()
                    .filter(|interval| {
                        interval.start <= range.end && interval.end >= range.start
                    })
                    .collect()
            })
    }

    /// Count total intervals stored
    pub fn count(&self) -> usize {
        self.intervals.values().map(|v| v.len()).sum()
    }

    /// Get all intervals for a chromosome
    pub fn get_chromosome_intervals(&self, chromosome: &str) -> &[GenomicInterval] {
        self.intervals.get(chromosome).map(|v| v.as_slice()).unwrap_or(&[])
    }

    /// Clear all intervals
    pub fn clear(&mut self) {
        self.intervals.clear();
    }
}

/// Query result bucket for overlapping and fully contained intervals.
/// # Invariants
/// `count` tracks total entries added via `add_*` helpers.
/// # Ownership
/// Owns cloned [`GenomicInterval`] vectors.
/// # Mutation
/// Mutated while assembling query results.
/// # Biological assumptions
/// Overlap vs containment reflects feature intersection semantics on reference.
/// # Java equivalence
/// None / Rust-native query DTO.
#[derive(Debug, Clone)]
pub struct IntervalQuery {
    pub overlapping: Vec<GenomicInterval>,
    pub contained: Vec<GenomicInterval>,
    pub count: usize,
}

impl IntervalQuery {
    /// Create a new interval query result
    pub fn new() -> Self {
        Self {
            overlapping: Vec::new(),
            contained: Vec::new(),
            count: 0,
        }
    }

    /// Add an overlapping interval
    pub fn add_overlapping(&mut self, interval: GenomicInterval) {
        self.overlapping.push(interval);
        self.count += 1;
    }

    /// Add a contained interval
    pub fn add_contained(&mut self, interval: GenomicInterval) {
        self.contained.push(interval);
        self.count += 1;
    }

    /// Check if any intervals found
    pub fn has_results(&self) -> bool {
        self.count > 0
    }
}

/// Facade over [`IntervalTree`] for position/range/containment queries.
/// # Invariants
/// Wraps a single tree rebuilt via `build_tree`.
/// # Ownership
/// Owns [`IntervalTree`]; returns owned [`IntervalQuery`] results.
/// # Mutation
/// `build_tree` replaces tree contents; queries are read-only.
/// # Biological assumptions
/// Powers target region lookup for walkers and QC.
/// # Java equivalence
/// Approximates GATK interval parsing + lookup utilities (Rust-native).
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
    pub fn query_position(&self, chromosome: &str, position: u64) -> IntervalQuery {
        let overlapping = self.tree.find_overlapping(chromosome, position);
        IntervalQuery {
            overlapping,
            contained: Vec::new(),
            count: overlapping.len(),
        }
    }

    /// Query intervals intersecting a range
    pub fn query_range(&self, chromosome: &str, range: Range<u64>) -> IntervalQuery {
        let intersecting = self.tree.find_intersecting(chromosome, &range);
        IntervalQuery {
            overlapping: intersecting,
            contained: Vec::new(),
            count: intersecting.len(),
        }
    }

    /// Query intervals containing a position
    pub fn query_contained(&self, chromosome: &str, position: u64) -> IntervalQuery {
        let mut contained = Vec::new();
        let count = self.tree.intervals
            .get(chromosome)
            .map(|intervals| intervals.len())
            .unwrap_or(0);

        if let Some(intervals) = self.tree.intervals.get(chromosome) {
            for interval in intervals {
                if interval.contains(position) {
                    // CLONE: needed because owned element into collection.
                    contained.push(interval.clone());
                }
            }
        }

        IntervalQuery {
            overlapping: Vec::new(),
            contained,
            count: contained.len(),
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
            min_length,
            max_length,
            average_length: if total_intervals > 0 { total_length / total_intervals as u64 } else { 0 },
        }
    }
}

/// Aggregate statistics over an indexed interval collection.
/// # Invariants
/// Length stats derived from [`GenomicInterval::length`] (0-based inclusive).
/// # Ownership
/// Plain scalars; serde-friendly clone.
/// # Mutation
/// Immutable snapshot from [`IntervalQueryEngine::get_statistics`].
/// # Biological assumptions
/// Summarizes target panel or callable territory footprint.
/// # Java equivalence
/// None / Rust-native.
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
/// # Invariants
/// Uniqueness via `Ord` on intervals; sorted iteration order.
/// # Ownership
/// Owns intervals in a `BTreeSet`; clone duplicates set.
/// # Mutation
/// Insert/remove mutate the set; set algebra ops copy inputs.
/// # Biological assumptions
/// Deduped interval lists (targets, exons, padded regions).
/// # Java equivalence
/// None / Rust-native.
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
}

/// Parallel interval overlap query.
/// Compute uses `par_iter`; results are merged into a [`BTreeMap`] keyed by
/// `(chromosome, position)` so return order is deterministic regardless
/// Rayon thread count. (This module is currently not wired into `gatk_core::lib`
/// kept correct for future use.)
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
                    interval.chromosome == *chromosome && interval.contains(*position)
                })
                .collect();
            ((chromosome.clone(), *position), overlapping)
        })
        .collect();
    rows.sort_by(|a, b| a.0.cmp(&b.0));
    rows.into_iter().collect()
}

/// Merge overlapping intervals in parallel
pub fn merge_overlapping_intervals_parallel(intervals: &[GenomicInterval]) -> Vec<GenomicInterval> {
    let mut sorted_intervals: Vec<_> = intervals.iter().cloned().collect();
    sorted_intervals.sort_by(|a, b| a.start.cmp(&b.start));
    
    let mut merged = Vec::new();
    for interval in sorted_intervals {
        if merged.is_empty() || !merged.last().unwrap().overlaps(interval) {
            merged.push(interval);
        } else {
            let last = merged.last_mut().unwrap();
            if last.chromosome == interval.chromosome {
                // Same chromosome - merge overlapping intervals
                let union_result = last.union(&interval);
                if union_result.len() == 1 {
                    // Overlapping intervals were merged
                    *last = union_result[0].clone();
                } else {
                    // Non-overlapping intervals
                    merged.push(interval);
                }
            } else {
                // Different chromosome - keep both
                merged.push(interval);
            }
        }
    }
    
    merged
}
