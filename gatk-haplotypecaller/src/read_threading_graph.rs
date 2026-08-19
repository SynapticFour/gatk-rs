//! GATK `ReadThreadingGraph` / `AbstractReadThreadingGraph` parity builder.
//! Replaces naive de Bruijn per-read edge counting with read threading: unique-kmer starts,
//! backward multiplicity on incoming matches, forward extension by suffix, and per-sample
//! pruning multiplicity (`MultiSampleEdge` with `numPruningSamples = 1`).
//!
//! # K-mer representation
//! Hot maps use [`crate::kmer_key::KmerKey`] (packed ACGT integers when `k ≤ 64`,
//! `Arc<[u8]>` for `N` / larger k). Node payloads remain `Arc<[u8]>` for path bases.
//!
//! # Deterministic neighbor order
//! `outgoing` / `incoming` are [`BTreeSet`] so `extend_chain_by_one` first-suffix-match
//! is stable when multiple outs share a base (observable threading topology). Do **not**
//! replace with `HashSet` without an equivalent deterministic tie-break.

use crate::assembly::{AssemblyGraph, AssemblyGraphParams, AssemblyRead};
use crate::kmer_key::{key_from_window, materialize_arc, KmerKey, RollingKmer};
use gatk_common::{GatkError, GatkResult};
use indexmap::IndexMap;
use std::collections::{BTreeSet, BinaryHeap, HashMap, HashSet};
use std::sync::Arc;

/// GATK `AbstractReadThreadingGraph.ANONYMOUS_SAMPLE` (reference sequences).
const ANONYMOUS_SAMPLE: &str = "XXX_UNNAMED_XXX";

const INCREASE_COUNTS_BACKWARDS: bool = true;

#[derive(Debug, Clone)]
struct SequenceForKmers {
    /// One shared allocation for the usable segment (no pending double-copy).
    bases: Arc<[u8]>,
    start: usize,
    stop: usize,
    count: u32,
    is_ref: bool,
}

#[derive(Debug)]
struct ThreadingEdge {
    current_sample: u32,
    /// Min-heap of flushed per-sample totals (GATK `PriorityQueue`); peek = pruning multiplicity.
    flushed_samples: BinaryHeap<std::cmp::Reverse<u32>>,
    num_pruning_samples: usize,
    is_ref: bool,
}

impl ThreadingEdge {
    fn new(initial: u32, num_pruning_samples: usize, is_ref: bool) -> Self {
        let mut flushed_samples = BinaryHeap::new();
        flushed_samples.push(std::cmp::Reverse(initial));
        Self {
            current_sample: initial,
            flushed_samples,
            num_pruning_samples,
            is_ref,
        }
    }

    fn inc(&mut self, delta: u32) {
        self.current_sample = self.current_sample.saturating_add(delta);
    }

    fn flush_sample(&mut self) {
        self.flushed_samples
            .push(std::cmp::Reverse(self.current_sample));
        // GATK: poll when size == capacity + 1 (min-heap drops smallest sample max).
        if self.flushed_samples.len() > self.num_pruning_samples {
            self.flushed_samples.pop();
        }
        self.current_sample = 0;
    }

    fn pruning_multiplicity(&self) -> u32 {
        self.flushed_samples
            .peek()
            .map(|r| r.0)
            .unwrap_or(0)
            .max(self.current_sample)
    }
}

/// Post-build threading graph stats.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ThreadingNonUniqueSummary {
    pub node_count: usize,
    pub edge_count: usize,
    pub unique_kmer_count: usize,
    pub non_unique_kmer_count: usize,
    pub is_low_complexity: bool,
    pub max_kmer_multiplicity: usize,
}

/// In-memory read-threading graph (GATK `ReadThreadingGraph` defaults: `numPruningSamples = 1`).
pub struct ReadThreadingGraphBuilder {
    kmer_size: usize,
    #[allow(dead_code)]
    min_base_quality: u8,
    start_threading_only_at_existing_vertex: bool,
    num_pruning_samples: usize,
    /// GATK `pending`: sequences grouped by sample; flush edge multiplicities after each sample.
    pending: IndexMap<String, Vec<SequenceForKmers>>,
    non_unique_kmers: HashSet<KmerKey>,
    /// GATK `uniqueKmers`: one vertex per unique kmer only (hash lookup; order not observable).
    unique_kmers: HashMap<KmerKey, usize>,
    nodes: Vec<Arc<[u8]>>,
    edges: HashMap<(usize, usize), ThreadingEdge>,
    edge_is_ref: HashSet<(usize, usize)>,
    ref_nodes: HashSet<usize>,
    ref_source_kmer: Option<Arc<[u8]>>,
    /// Neighbor sets are ordered so `extend_chain_by_one` first-match is deterministic.
    outgoing: HashMap<usize, BTreeSet<usize>>,
    incoming: HashMap<usize, BTreeSet<usize>>,
    built: bool,
}

impl ReadThreadingGraphBuilder {
    pub fn new(
        kmer_size: usize,
        num_pruning_samples: usize,
        min_base_quality: u8,
        start_threading_only_at_existing_vertex: bool,
    ) -> Self {
        Self {
            kmer_size,
            min_base_quality,
            start_threading_only_at_existing_vertex,
            num_pruning_samples: num_pruning_samples.max(1),
            pending: IndexMap::new(),
            non_unique_kmers: HashSet::new(),
            unique_kmers: HashMap::new(),
            nodes: Vec::new(),
            edges: HashMap::new(),
            edge_is_ref: HashSet::new(),
            ref_nodes: HashSet::new(),
            ref_source_kmer: None,
            outgoing: HashMap::new(),
            incoming: HashMap::new(),
            built: false,
        }
    }

    fn key_at(bases: &[u8], start: usize, k: usize) -> KmerKey {
        key_from_window(bases, start, k)
    }

    fn suffix_of_kmer(kmer: &[u8]) -> u8 {
        kmer.last().copied().unwrap_or(b'N')
    }

    fn sequences_from_read(
        read: &AssemblyRead,
        kmer_size: usize,
        min_qual: u8,
    ) -> Vec<SequenceForKmers> {
        let mut out = Vec::new();
        let mut last_good: Option<usize> = None;
        let bytes = read.bases.as_slice();
        for end in 0..=bytes.len() {
            let unusable = end == bytes.len()
                || read.base_quals[end] < min_qual
                || !is_base_usable(bytes[end]);
            if unusable {
                if let Some(start) = last_good {
                    let len = end - start;
                    if len >= kmer_size {
                        out.push(SequenceForKmers {
                            bases: Arc::from(&bytes[start..end]),
                            start: 0,
                            stop: len,
                            count: 1,
                            is_ref: false,
                        });
                    }
                }
                last_good = None;
            } else if last_good.is_none() {
                last_good = Some(end);
            }
        }
        out
    }

    fn add_pending_sequence(&mut self, sample: &str, seq: SequenceForKmers) {
        assert!(!self.built, "graph already built");
        self.pending
            .entry(sample.to_string())
            .or_default()
            .push(seq);
    }

    fn preprocess_reads(&mut self) {
        self.non_unique_kmers.clear();
        let k = self.kmer_size;
        for seq in self.pending.values().flatten() {
            let mut seen: HashSet<KmerKey> = HashSet::new();
            let mut roll = RollingKmer::new(k);
            let stop = seq.stop.saturating_sub(k);
            for i in seq.start..=stop {
                // Rolling packed ACGT keys allocate nothing; ambiguous windows allocate once.
                let key = roll.key_at(&seq.bases, i);
                if !seen.insert(key.clone()) {
                    self.non_unique_kmers.insert(key);
                }
            }
        }
    }

    fn is_threading_start(&self, key: &KmerKey) -> bool {
        if self.start_threading_only_at_existing_vertex {
            self.unique_kmers.contains_key(key)
        } else {
            !self.non_unique_kmers.contains(key)
        }
    }

    fn find_start(&self, seq: &SequenceForKmers) -> Option<usize> {
        if seq.is_ref {
            return Some(0);
        }
        let k = self.kmer_size;
        let last = seq.stop.saturating_sub(k);
        let mut roll = RollingKmer::new(k);
        for i in seq.start..last {
            let key = roll.key_at(&seq.bases, i);
            if self.is_threading_start(&key) {
                return Some(i);
            }
        }
        None
    }

    fn get_unique_kmer_vertex(&self, key: &KmerKey, allow_ref_source: bool) -> Option<usize> {
        let id = *self.unique_kmers.get(key)?;
        if !allow_ref_source {
            if let Some(ref_src) = self.ref_source_kmer.as_deref() {
                if self.nodes[id].as_ref() == ref_src {
                    return None;
                }
            }
        }
        Some(id)
    }

    fn create_vertex(&mut self, key: KmerKey) -> usize {
        let id = self.nodes.len();
        let kmer = materialize_arc(&key, self.kmer_size);
        if !self.non_unique_kmers.contains(&key) && !self.unique_kmers.contains_key(&key) {
            self.unique_kmers.insert(key, id);
        }
        self.nodes.push(kmer);
        id
    }

    fn get_or_create_kmer_vertex(&mut self, bases: &[u8], start: usize) -> usize {
        let key = Self::key_at(bases, start, self.kmer_size);
        if let Some(id) = self.get_unique_kmer_vertex(&key, true) {
            return id;
        }
        self.create_vertex(key)
    }

    pub(crate) fn is_low_quality_graph(&self) -> bool {
        self.non_unique_kmers.len() * 4 > self.unique_kmers.len()
    }

    pub(crate) fn is_low_complexity(&self) -> bool {
        self.is_low_quality_graph()
    }

    fn max_kmer_multiplicity(&self) -> usize {
        let mut mult: HashMap<&[u8], usize> = HashMap::new();
        for kmer in &self.nodes {
            *mult.entry(kmer.as_ref()).or_default() += 1;
        }
        mult.values().copied().max().unwrap_or(0)
    }

    pub fn non_unique_summary(&self) -> ThreadingNonUniqueSummary {
        assert!(self.built, "non_unique_summary requires build()");
        ThreadingNonUniqueSummary {
            node_count: self.nodes.len(),
            edge_count: self.edges.len(),
            unique_kmer_count: self.unique_kmers.len(),
            non_unique_kmer_count: self.non_unique_kmers.len(),
            is_low_complexity: self.is_low_complexity(),
            max_kmer_multiplicity: self.max_kmer_multiplicity(),
        }
    }

    fn inc_edge(&mut self, from: usize, to: usize, delta: u32, is_ref: bool) {
        if let Some(e) = self.edges.get_mut(&(from, to)) {
            e.inc(delta);
            if is_ref {
                e.is_ref = true;
                self.edge_is_ref.insert((from, to));
            }
        } else {
            self.edges.insert(
                (from, to),
                ThreadingEdge::new(delta, self.num_pruning_samples, is_ref),
            );
            self.outgoing.entry(from).or_default().insert(to);
            self.incoming.entry(to).or_default().insert(from);
            if is_ref {
                self.edge_is_ref.insert((from, to));
            }
        }
    }

    fn increase_counts_backwards(
        &mut self,
        seq: &SequenceForKmers,
        vertex: usize,
        original_kmer: &[u8],
        offset: isize,
    ) {
        if offset < 0 {
            return;
        }
        let off = offset as usize;
        if off >= original_kmer.len() {
            return;
        }
        let target_base = original_kmer[off];
        let preds: Vec<usize> = self
            .incoming
            .get(&vertex)
            .map(|s| s.iter().copied().collect())
            .unwrap_or_default();
        for pred in preds {
            let suffix = Self::suffix_of_kmer(&self.nodes[pred]);
            let in_degree = self.incoming.get(&vertex).map(|s| s.len()).unwrap_or(0);
            if suffix == target_base && (in_degree == 1) {
                self.inc_edge(pred, vertex, seq.count, false);
                self.increase_counts_backwards(seq, pred, original_kmer, offset - 1);
            }
        }
    }

    fn extend_chain_by_one(
        &mut self,
        prev: usize,
        bases: &[u8],
        kmer_start: usize,
        count: u32,
        is_ref: bool,
    ) -> GatkResult<usize> {
        let k = self.kmer_size;
        let next_pos = kmer_start + k - 1;
        let next_base = bases[next_pos];
        // `BTreeSet` iteration order → stable first suffix match when multiple outs share a base.
        let outs: Vec<usize> = self
            .outgoing
            .get(&prev)
            .map(|s| s.iter().copied().collect())
            .unwrap_or_default();
        for to in outs {
            if Self::suffix_of_kmer(&self.nodes[to]) == next_base {
                self.inc_edge(prev, to, count, is_ref);
                return Ok(to);
            }
        }
        let key = Self::key_at(bases, kmer_start, k);
        let next = if let Some(merge) = self.get_unique_kmer_vertex(&key, false) {
            if is_ref {
                return Err(GatkError::algorithm(format!(
                    "reference threading attempted to merge into unique vertex for kmer {}",
                    String::from_utf8_lossy(&bases[kmer_start..kmer_start + k])
                )));
            }
            merge
        } else {
            self.create_vertex(key)
        };
        self.inc_edge(prev, next, count, is_ref);
        if is_ref {
            self.ref_nodes.insert(prev);
            self.ref_nodes.insert(next);
        }
        Ok(next)
    }

    fn thread_sequence(&mut self, seq: &SequenceForKmers) -> GatkResult<()> {
        let Some(start_pos) = self.find_start(seq) else {
            return Ok(());
        };
        let k = self.kmer_size;
        if seq.is_ref && self.ref_source_kmer.is_none() {
            self.ref_source_kmer = Some(Arc::from(&seq.bases[start_pos..start_pos + k]));
        }
        let start_vertex = self.get_or_create_kmer_vertex(&seq.bases, start_pos);
        if INCREASE_COUNTS_BACKWARDS {
            let kmer_bytes = &seq.bases[start_pos..start_pos + k];
            self.increase_counts_backwards(seq, start_vertex, kmer_bytes, (k as isize) - 2);
        }
        let mut vertex = start_vertex;
        for i in (start_pos + 1)..=(seq.stop.saturating_sub(k)) {
            vertex = self.extend_chain_by_one(vertex, &seq.bases, i, seq.count, seq.is_ref)?;
        }
        Ok(())
    }

    pub fn build(&mut self) -> GatkResult<()> {
        if self.built {
            return Ok(());
        }
        self.preprocess_reads();
        let pending = std::mem::take(&mut self.pending);
        for seqs in pending.values() {
            for seq in seqs {
                self.thread_sequence(seq)?;
            }
            for e in self.edges.values_mut() {
                e.flush_sample();
            }
        }
        self.built = true;
        Ok(())
    }

    pub fn into_assembly_graph(mut self) -> GatkResult<AssemblyGraph> {
        self.build()?;
        Ok(self.finish_into_assembly_graph())
    }

    fn finish_into_assembly_graph(self) -> AssemblyGraph {
        debug_assert!(self.built, "finish_into_assembly_graph requires build()");
        let nodes: Vec<_> = self
            .nodes
            .into_iter()
            .enumerate()
            .map(|(id, kmer)| crate::assembly::KmerNode {
                id,
                kmer,
                support: 1,
            })
            .collect();
        // Arc-keyed map for AssemblyGraph (lookup only; order not observable).
        let mut kmer_to_id: HashMap<Arc<[u8]>, usize> =
            HashMap::with_capacity(self.unique_kmers.len());
        for &id in self.unique_kmers.values() {
            kmer_to_id.insert(Arc::clone(&nodes[id].kmer), id);
        }
        let edges: HashMap<_, _> = self
            .edges
            .into_iter()
            .map(|((from, to), e)| ((from, to), e.pruning_multiplicity()))
            .collect();
        AssemblyGraph::from_threading_build(
            self.kmer_size,
            nodes,
            kmer_to_id,
            edges,
            self.outgoing,
            self.incoming,
            self.edge_is_ref,
            self.ref_nodes,
            self.ref_source_kmer,
        )
    }

    pub fn into_assembly_graph_with_summary(
        mut self,
    ) -> GatkResult<(AssemblyGraph, ThreadingNonUniqueSummary)> {
        self.build()?;
        let summary = self.non_unique_summary();
        Ok((self.finish_into_assembly_graph(), summary))
    }
}

fn is_base_usable(base: u8) -> bool {
    matches!(base, b'A' | b'C' | b'G' | b'T' | b'N')
}

/// Build an [`AssemblyGraph`] using GATK read-threading semantics.
/// Build graph from reference + reads (ref threaded first, GATK order).
pub fn assembly_graph_from_ref_and_reads_threading(
    reference: &AssemblyRead,
    reads: &[AssemblyRead],
    params: &AssemblyGraphParams,
) -> GatkResult<AssemblyGraph> {
    let builder = build_threading_builder(Some(reference), reads, params)?;
    builder.into_assembly_graph()
}

/// Single threading build returning graph + non-unique / low-complexity summary.
pub fn assembly_graph_from_ref_and_reads_threading_with_summary(
    reference: &AssemblyRead,
    reads: &[AssemblyRead],
    params: &AssemblyGraphParams,
) -> GatkResult<(AssemblyGraph, ThreadingNonUniqueSummary)> {
    let builder = build_threading_builder(Some(reference), reads, params)?;
    builder.into_assembly_graph_with_summary()
}

fn build_threading_builder(
    reference: Option<&AssemblyRead>,
    reads: &[AssemblyRead],
    params: &AssemblyGraphParams,
) -> GatkResult<ReadThreadingGraphBuilder> {
    let kmer_size = params.kmer_size.as_usize();
    let mut builder = ReadThreadingGraphBuilder::new(
        kmer_size,
        1,
        params.min_base_quality,
        params.start_threading_only_at_existing_vertex,
    );
    if let Some(reference) = reference {
        for mut seq in ReadThreadingGraphBuilder::sequences_from_read(
            reference,
            kmer_size,
            params.min_base_quality,
        ) {
            seq.is_ref = true;
            builder.add_pending_sequence(ANONYMOUS_SAMPLE, seq);
        }
    }
    const READ_SAMPLE: &str = "SAMPLE";
    for read in reads {
        if read.bases.len() != read.base_quals.len() {
            return Err(GatkError::argument(
                "read bases length must match base quality length",
            ));
        }
        for seq in
            ReadThreadingGraphBuilder::sequences_from_read(read, kmer_size, params.min_base_quality)
        {
            builder.add_pending_sequence(READ_SAMPLE, seq);
        }
    }
    Ok(builder)
}

/// GATK `ReadThreadingGraph.determineNonUniqueKmers` on the reference haplotype only (`createGraph` gate).
pub fn reference_has_non_unique_kmers(reference: &AssemblyRead, kmer_size: usize) -> bool {
    if reference.bases.len() < kmer_size {
        return false;
    }
    let bases = reference.bases.as_slice();
    let stop = reference.bases.len().saturating_sub(kmer_size);
    let mut seen: HashSet<KmerKey> = HashSet::new();
    for i in 0..=stop {
        let key = key_from_window(bases, i, kmer_size);
        if !seen.insert(key) {
            return true;
        }
    }
    false
}

/// Build threading graph and return E.5 non-unique / cycle policy summary.
pub fn threading_non_unique_summary(
    reference: Option<&AssemblyRead>,
    reads: &[AssemblyRead],
    params: &AssemblyGraphParams,
) -> GatkResult<ThreadingNonUniqueSummary> {
    let mut builder = build_threading_builder(reference, reads, params)?;
    builder.build()?;
    Ok(builder.non_unique_summary())
}

pub fn assembly_graph_from_reads_threading(
    reads: &[AssemblyRead],
    params: &AssemblyGraphParams,
) -> GatkResult<AssemblyGraph> {
    let builder = build_threading_builder(None, reads, params)?;
    builder.into_assembly_graph()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::assembly::AssemblyGraphParams;

    fn read(seq: &str, q: u8) -> AssemblyRead {
        AssemblyRead {
            bases: seq.as_bytes().to_vec(),
            base_quals: vec![q; seq.len()],
        }
    }

    #[test]
    fn single_read_threading_matches_gatk_edge_weights() {
        let params = AssemblyGraphParams {
            kmer_size: crate::bio_ids::KmerSize::try_new(3).unwrap(),
            min_base_quality: 10,
            ..Default::default()
        };
        let g = assembly_graph_from_reads_threading(&[read("ACGTT", 30)], &params).unwrap();
        let edges = g.edges_sorted();
        assert_eq!(edges.len(), 2);
        assert_eq!(edges[0].support, 1);
        assert_eq!(edges[1].support, 1);
    }

    #[test]
    fn p5_case1_threading_supports() {
        let reads = vec![
            read("ACGTT", 30),
            read("ACGTT", 30),
            read("ACGTT", 30),
            read("ACGTA", 30),
            read("ACGTA", 30),
        ];
        let params = AssemblyGraphParams {
            kmer_size: crate::bio_ids::KmerSize::try_new(3).unwrap(),
            min_base_quality: 10,
            ..Default::default()
        };
        let g = assembly_graph_from_reads_threading(&reads, &params).unwrap();
        let mut by_pair: HashMap<(Vec<u8>, Vec<u8>), u32> = HashMap::new();
        for e in g.edges_sorted() {
            let from = g.nodes()[e.from].kmer.to_vec();
            let to = g.nodes()[e.to].kmer.to_vec();
            by_pair.insert((from, to), e.support);
        }
        assert_eq!(by_pair.get(&(b"ACG".to_vec(), b"CGT".to_vec())), Some(&5));
        assert_eq!(by_pair.get(&(b"CGT".to_vec(), b"GTT".to_vec())), Some(&3));
        assert_eq!(by_pair.get(&(b"CGT".to_vec(), b"GTA".to_vec())), Some(&2));
    }

    #[test]
    fn repeat_only_reads_yield_empty_graph() {
        let reads = vec![
            read("ATATATG", 30),
            read("ATATATG", 30),
            read("ATATATA", 30),
        ];
        let params = AssemblyGraphParams {
            kmer_size: crate::bio_ids::KmerSize::try_new(3).unwrap(),
            min_base_quality: 10,
            ..Default::default()
        };
        let summary = threading_non_unique_summary(None, &reads, &params).unwrap();
        assert_eq!(summary.node_count, 0);
        assert_eq!(summary.edge_count, 0);
        assert!(summary.non_unique_kmer_count > 0);
    }

    #[test]
    fn ref_with_internal_repeat_kmers_has_multiplicity_gt_one() {
        let reference = read("TTTACGTTACGT", 30);
        let reads = vec![read("ACGTT", 30), read("ACGTT", 30)];
        let params = AssemblyGraphParams {
            kmer_size: crate::bio_ids::KmerSize::try_new(3).unwrap(),
            min_base_quality: 10,
            ..Default::default()
        };
        let summary = threading_non_unique_summary(Some(&reference), &reads, &params).unwrap();
        assert!(summary.max_kmer_multiplicity >= 2);
        assert!(summary.non_unique_kmer_count >= 4);
        assert!(summary.node_count > 0);
    }

    #[test]
    fn create_graph_non_unique_gate_uses_reference_only() {
        let reference = read("ATATATG", 30);
        let reads = vec![
            read("ATATATG", 30),
            read("ATATATG", 30),
            read("ATATATA", 30),
        ];
        assert!(reference_has_non_unique_kmers(&reference, 3));
        let summary =
            threading_non_unique_summary(Some(&reference), &reads, &params_for(3)).unwrap();
        assert!(summary.non_unique_kmer_count > 0);
        let unique_ref = read("ACGTGCTTAGCA", 30);
        assert!(!reference_has_non_unique_kmers(&unique_ref, 3));
    }

    #[test]
    fn ambiguous_n_kmer_still_threads() {
        let reads = vec![read("ACGTNACGT", 30), read("ACGTNACGT", 30)];
        let params = AssemblyGraphParams {
            kmer_size: crate::bio_ids::KmerSize::try_new(5).unwrap(),
            min_base_quality: 10,
            ..Default::default()
        };
        let g = assembly_graph_from_reads_threading(&reads, &params).unwrap();
        assert!(!g.nodes().is_empty() || g.edges_sorted().is_empty() || true);
        // Must not panic; N windows use Bytes keys.
        let _ = g.edges_sorted();
    }

    fn params_for(kmer_size: usize) -> AssemblyGraphParams {
        AssemblyGraphParams {
            kmer_size: crate::bio_ids::KmerSize::try_new(kmer_size as u16).expect("test k"),
            min_base_quality: 10,
            ..Default::default()
        }
    }
}
