//! GATK `ReadThreadingGraph` / `AbstractReadThreadingGraph` parity builder.
//! Replaces naive de Bruijn per-read edge counting with read threading: unique-kmer starts,
//! backward multiplicity on incoming matches, forward extension by suffix, and per-sample
//! pruning multiplicity (`MultiSampleEdge` with `numPruningSamples = 1`).

use crate::assembly::{AssemblyGraph, AssemblyGraphParams, AssemblyRead};
use gatk_common::{GatkError, GatkResult};
use indexmap::IndexMap;
use std::collections::{BTreeMap, BTreeSet, BinaryHeap, HashMap, HashSet};

/// GATK `AbstractReadThreadingGraph.ANONYMOUS_SAMPLE` (reference sequences).
const ANONYMOUS_SAMPLE: &str = "XXX_UNNAMED_XXX";

const INCREASE_COUNTS_BACKWARDS: bool = true;

#[derive(Debug, Clone)]
struct SequenceForKmers {
    #[allow(dead_code)]
    name: String,
    bases: Vec<u8>,
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
    non_unique_kmers: HashSet<String>,
    /// GATK `uniqueKmers`: one vertex per unique kmer only.
    unique_kmers: BTreeMap<String, usize>,
    nodes: Vec<String>,
    edges: HashMap<(usize, usize), ThreadingEdge>,
    edge_is_ref: HashSet<(usize, usize)>,
    ref_nodes: HashSet<usize>,
    ref_source_kmer: Option<String>,
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
            unique_kmers: BTreeMap::new(),
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

    fn kmer_at(bases: &[u8], start: usize, k: usize) -> String {
        // Hot path: assembly bases are ASCII ACGTN; avoid lossy UTF-8 scanning.
        match std::str::from_utf8(&bases[start..start + k]) {
            Ok(s) => s.to_owned(),
            Err(_) => String::from_utf8_lossy(&bases[start..start + k]).into_owned(),
        }
    }

    fn suffix_of_kmer(kmer: &str) -> u8 {
        kmer.as_bytes().last().copied().unwrap_or(b'N')
    }

    fn sequences_from_read(
        read: &AssemblyRead,
        kmer_size: usize,
        min_qual: u8,
    ) -> Vec<SequenceForKmers> {
        let mut out = Vec::new();
        let mut last_good: Option<usize> = None;
        let bytes = read.bases.as_bytes();
        for end in 0..=bytes.len() {
            let unusable = end == bytes.len()
                || read.base_quals[end] < min_qual
                || !is_base_usable(bytes[end]);
            if unusable {
                if let Some(start) = last_good {
                    let len = end - start;
                    if len >= kmer_size {
                        // Store only the usable segment (callers treat start/stop as offsets into `bases`).
                        out.push(SequenceForKmers {
                            name: format!("{start}_{end}"),
                            bases: bytes[start..end].to_vec(),
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

    /// GATK `addSequence` stores a subread-only byte array (see `addRead` / parity `addReadThreadingSequence`).
    fn add_pending_subread(
        &mut self,
        sample: &str,
        name: String,
        bases: Vec<u8>,
        count: u32,
        is_ref: bool,
    ) {
        let len = bases.len();
        if len < self.kmer_size {
            return;
        }
        self.add_pending_sequence(
            sample,
            SequenceForKmers {
                name,
                bases,
                start: 0,
                stop: len,
                count,
                is_ref,
            },
        );
    }

    fn preprocess_reads(&mut self) {
        self.non_unique_kmers.clear();
        let k = self.kmer_size;
        for seq in self.pending.values().flatten() {
            let mut seen = HashSet::new();
            // GATK scans 0..stop-k on the subread byte array; we keep [start, stop) in full read coordinates.
            let stop = seq.stop.saturating_sub(k);
            for i in seq.start..=stop {
                let key = Self::kmer_at(&seq.bases, i, k);
                // Lifetime: first sighting moves into `seen`; a later duplicate moves into
                // `non_unique_kmers`. No clone — both sets only need owned keys they retain.
                if seen.contains(&key) {
                    self.non_unique_kmers.insert(key);
                } else {
                    seen.insert(key);
                }
            }
        }
    }

    fn is_threading_start(&self, kmer: &str) -> bool {
        if self.start_threading_only_at_existing_vertex {
            self.unique_kmers.contains_key(kmer)
        } else {
            !self.non_unique_kmers.contains(kmer)
        }
    }

    fn find_start(&self, seq: &SequenceForKmers) -> Option<usize> {
        // GATK `findStart`: reference always threads from index 0 of the ref sequence array.
        if seq.is_ref {
            return Some(0);
        }
        let k = self.kmer_size;
        let last = seq.stop.saturating_sub(k);
        for i in seq.start..last {
            // GATK: i < stop - kmerSize
            let key = Self::kmer_at(&seq.bases, i, k);
            if self.is_threading_start(&key) {
                return Some(i);
            }
        }
        None
    }

    fn get_unique_kmer_vertex(&self, kmer: &str, allow_ref_source: bool) -> Option<usize> {
        if !allow_ref_source && self.ref_source_kmer.as_deref() == Some(kmer) {
            return None;
        }
        self.unique_kmers.get(kmer).copied()
    }

    fn create_vertex(&mut self, kmer: String) -> usize {
        let id = self.nodes.len();
        // Unique kmers need the string in both `nodes` and `unique_kmers` (HashMap key).
        // Non-unique: move into `nodes` only — avoid the prior always-clone-on-push.
        if !self.non_unique_kmers.contains(&kmer) && !self.unique_kmers.contains_key(&kmer) {
            // CLONE: needed because HashMap key and node storage both own the kmer string.
            self.unique_kmers.insert(kmer.clone(), id);
        }
        self.nodes.push(kmer);
        id
    }

    fn get_or_create_kmer_vertex(&mut self, bases: &[u8], start: usize) -> usize {
        let kmer = Self::kmer_at(bases, start, self.kmer_size);
        if let Some(id) = self.get_unique_kmer_vertex(&kmer, true) {
            return id;
        }
        self.create_vertex(kmer)
    }

    /// GATK `ReadThreadingGraph.isLowQualityGraph` / `isLowComplexity`.
    pub(crate) fn is_low_quality_graph(&self) -> bool {
        self.non_unique_kmers.len() * 4 > self.unique_kmers.len()
    }

    /// Alias retained for E.5 parity naming.
    pub(crate) fn is_low_complexity(&self) -> bool {
        self.is_low_quality_graph()
    }

    fn max_kmer_multiplicity(&self) -> usize {
        let mut mult: HashMap<&str, usize> = HashMap::new();
        for kmer in &self.nodes {
            *mult.entry(kmer.as_str()).or_default() += 1;
        }
        mult.values().copied().max().unwrap_or(0)
    }

    /// Post-`build` threading stats for parity.
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
    ) -> usize {
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
                return to;
            }
        }
        let kmer = Self::kmer_at(bases, kmer_start, k);
        let next = if let Some(merge) = self.get_unique_kmer_vertex(&kmer, false) {
            if is_ref {
                panic!("reference threading attempted to merge into unique vertex for kmer {kmer}");
            }
            merge
        } else {
            self.create_vertex(kmer)
        };
        self.inc_edge(prev, next, count, is_ref);
        if is_ref {
            self.ref_nodes.insert(prev);
            self.ref_nodes.insert(next);
        }
        next
    }

    fn thread_sequence(&mut self, seq: &SequenceForKmers) {
        let Some(start_pos) = self.find_start(seq) else {
            return;
        };
        let k = self.kmer_size;
        if seq.is_ref && self.ref_source_kmer.is_none() {
            self.ref_source_kmer = Some(Self::kmer_at(&seq.bases, start_pos, k));
        }
        let start_vertex = self.get_or_create_kmer_vertex(&seq.bases, start_pos);
        if INCREASE_COUNTS_BACKWARDS {
            let kmer_bytes = &seq.bases[start_pos..start_pos + k];
            self.increase_counts_backwards(seq, start_vertex, kmer_bytes, (k as isize) - 2);
        }
        let mut vertex = start_vertex;
        for i in (start_pos + 1)..=(seq.stop.saturating_sub(k)) {
            vertex = self.extend_chain_by_one(vertex, &seq.bases, i, seq.count, seq.is_ref);
        }
    }

    pub fn build(&mut self) {
        if self.built {
            return;
        }
        self.preprocess_reads();
        // Lifetime: after preprocess, pending is only read for threading; take moves it
        // so thread_sequence can mutably borrow self without cloning the sequence lists.
        let pending = std::mem::take(&mut self.pending);
        for seqs in pending.values() {
            for seq in seqs {
                self.thread_sequence(seq);
            }
            for e in self.edges.values_mut() {
                e.flush_sample();
            }
        }
        self.built = true;
    }

    pub fn into_assembly_graph(mut self) -> AssemblyGraph {
        self.build();
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
        let edges: HashMap<_, _> = self
            .edges
            .iter()
            .map(|(&(from, to), e)| ((from, to), e.pruning_multiplicity()))
            .collect();
        // GATK `buildGraphIfNecessary` keeps all vertices in `vertexSet` (no orphan cleanup here).
        AssemblyGraph::from_threading_build(
            self.kmer_size,
            nodes,
            self.unique_kmers,
            edges,
            self.outgoing,
            self.incoming,
            self.edge_is_ref,
            self.ref_nodes,
            self.ref_source_kmer,
        )
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
    Ok(builder.into_assembly_graph())
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
        for seq in ReadThreadingGraphBuilder::sequences_from_read(
            reference,
            kmer_size,
            params.min_base_quality,
        ) {
            builder.add_pending_subread(
                ANONYMOUS_SAMPLE,
                seq.name,
                seq.bases[seq.start..seq.stop].to_vec(),
                seq.count,
                true,
            );
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
            builder.add_pending_subread(
                READ_SAMPLE,
                seq.name,
                seq.bases[seq.start..seq.stop].to_vec(),
                seq.count,
                false,
            );
        }
    }
    Ok(builder)
}

/// GATK `ReadThreadingGraph.determineNonUniqueKmers` on the reference haplotype only (`createGraph` gate).
pub fn reference_has_non_unique_kmers(reference: &AssemblyRead, kmer_size: usize) -> bool {
    if reference.bases.len() < kmer_size {
        return false;
    }
    let bases = reference.bases.as_bytes();
    let stop = reference.bases.len().saturating_sub(kmer_size);
    let mut seen = HashSet::new();
    for i in 0..=stop {
        let kmer = ReadThreadingGraphBuilder::kmer_at(bases, i, kmer_size);
        if !seen.insert(kmer) {
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
    builder.build();
    Ok(builder.non_unique_summary())
}

pub fn assembly_graph_from_reads_threading(
    reads: &[AssemblyRead],
    params: &AssemblyGraphParams,
) -> GatkResult<AssemblyGraph> {
    let builder = build_threading_builder(None, reads, params)?;
    Ok(builder.into_assembly_graph())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::assembly::AssemblyGraphParams;

    fn read(seq: &str, q: u8) -> AssemblyRead {
        AssemblyRead {
            bases: seq.to_string(),
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
        let mut by_pair: HashMap<(String, String), u32> = HashMap::new();
        for e in g.edges_sorted() {
            let from = g.nodes()[e.from].kmer.clone();
            let to = g.nodes()[e.to].kmer.clone();
            by_pair.insert((from, to), e.support);
        }
        assert_eq!(by_pair.get(&("ACG".into(), "CGT".into())), Some(&5));
        assert_eq!(by_pair.get(&("CGT".into(), "GTT".into())), Some(&3));
        assert_eq!(by_pair.get(&("CGT".into(), "GTA".into())), Some(&2));
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

    fn params_for(kmer_size: usize) -> AssemblyGraphParams {
        AssemblyGraphParams {
            kmer_size: crate::bio_ids::KmerSize::try_new(kmer_size as u16).expect("test k"),
            min_base_quality: 10,
            ..Default::default()
        }
    }
}
