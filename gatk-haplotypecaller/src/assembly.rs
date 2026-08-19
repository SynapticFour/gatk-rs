//! Local assembly graph scaffolding.

use crate::bio_ids::KmerSize;
use gatk_common::{GatkError, GatkResult};
use std::cmp::Reverse;
use std::collections::{BTreeSet, BinaryHeap, HashMap, HashSet};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AssemblyRead {
    /// ASCII ACGTN bases (same bytes as prior `String` path for valid BAM/ref).
    pub bases: Vec<u8>,
    pub base_quals: Vec<u8>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AssemblyGraphParams {
    /// K-mer size (`≥ 2`); type-enforced via [`KmerSize`].
    pub kmer_size: KmerSize,
    pub min_base_quality: u8,
    pub min_edge_weight: u32,
    pub dangling_path_max_nodes: usize,
    pub max_haplotypes: usize,
    pub max_haplotype_bases: usize,
    /// GATK `setThreadingStartOnlyAtExistingVertex(!recoverDanglingBranches)`.
    pub start_threading_only_at_existing_vertex: bool,
}

impl Default for AssemblyGraphParams {
    fn default() -> Self {
        Self {
            kmer_size: KmerSize::DEFAULT_ASSEMBLY,
            min_base_quality: 10,
            min_edge_weight: 2,
            dangling_path_max_nodes: 4,
            max_haplotypes: 128,
            max_haplotype_bases: 512,
            start_threading_only_at_existing_vertex: false,
        }
    }
}

/// GATK `min-pruning` / adaptive pruning knobs (`ReadThreadingAssemblerArgumentCollection`).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AssemblyGraphPruningParams {
    /// GATK `minPruneFactor` — minimum edge multiplicity when adaptive pruning is off.
    pub min_prune_factor: u32,
    pub use_adaptive_pruning: bool,
    /// Natural-log LOD (`pruning-lod-threshold`; default `log10ToLog(1.0)` = ln 10).
    pub pruning_log_odds_threshold: f64,
    /// Natural-log LOD for seeding (`pruning-seeding-lod-threshold`; default `log10ToLog(4.0)`).
    pub pruning_seeding_log_odds_threshold: f64,
    pub initial_error_rate_for_pruning: f64,
    pub max_unpruned_variants: usize,
}

impl AssemblyGraphPruningParams {
    pub fn gatk_haplotype_caller_defaults() -> Self {
        const LN_10: f64 = std::f64::consts::LN_10;
        Self {
            min_prune_factor: 2,
            use_adaptive_pruning: false,
            pruning_log_odds_threshold: LN_10,
            pruning_seeding_log_odds_threshold: 4.0 * LN_10,
            initial_error_rate_for_pruning: 0.001,
            max_unpruned_variants: 100,
        }
    }
}

impl Default for AssemblyGraphPruningParams {
    fn default() -> Self {
        Self::gatk_haplotype_caller_defaults()
    }
}

/// Post-pruning graph metrics for L2 parity **E.3** (log-scale edge-mass summary).
#[derive(Debug, Clone, PartialEq)]
pub struct AssemblyGraphSummary {
    pub node_count: usize,
    pub edge_count: usize,
    pub log10_max_edge_support: f64,
    pub log10_sum_edge_support: f64,
    pub pruning_lod_threshold_ln: f64,
    pub adaptive_pruning: bool,
    pub min_prune_factor: u32,
    pub edges_pruned: u32,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct KmerNode {
    pub id: usize,
    /// Shared kmer bytes — one allocation; maps hold `Arc` clones (Java-style reference sharing).
    pub kmer: std::sync::Arc<[u8]>,
    pub support: u32,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct KmerEdge {
    pub from: usize,
    pub to: usize,
    pub support: u32,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CandidateHaplotype {
    pub sequence: Vec<u8>,
    pub support: u32,
}

/// Alt sequence + CIGAR from a successful ASM-1 dangling tail merge (Java merge junction).
#[derive(Debug, Clone)]
pub struct DanglingMergeHaplotype {
    pub alt_bases: Vec<u8>,
    pub cigar: crate::cigar::Cigar,
    pub alignment_start_hap_wrt_ref: usize,
}

#[derive(Debug, Clone, Default)]
pub struct AssemblyGraph {
    pub kmer_size: usize,
    /// Successful dangling merges that carry variation off the ref-spine (ASM-1 → ASM-8).
    pub(crate) dangling_merge_haps: Vec<DanglingMergeHaplotype>,
    nodes: Vec<KmerNode>,
    /// Lookup only — iteration order is not an observable (dumps sort by k-mer bytes).
    kmer_to_id: HashMap<std::sync::Arc<[u8]>, usize>,
    edges: HashMap<(usize, usize), u32>,
    outgoing: HashMap<usize, BTreeSet<usize>>,
    incoming: HashMap<usize, BTreeSet<usize>>,
    /// Edges on the reference path (`MultiSampleEdge.isRef`).
    pub(crate) ref_edges: HashSet<(usize, usize)>,
    /// Vertices visited when threading the reference.
    pub(crate) ref_nodes: HashSet<usize>,
    /// First reference kmer (`AbstractReadThreadingGraph.refSource`).
    #[allow(dead_code)] // carried from threading graph for future parity dumps
    pub(crate) ref_source_kmer: Option<std::sync::Arc<[u8]>>,
}

impl AssemblyGraph {
    pub fn new(kmer_size: usize) -> GatkResult<Self> {
        let _ = crate::bio_ids::KmerSize::try_from_usize(kmer_size)?;
        Ok(Self {
            kmer_size,
            ..Self::default()
        })
    }

    /// Construct from read-threading build output ([`crate::read_threading_graph`]).
    pub(crate) fn from_threading_build(
        kmer_size: usize,
        nodes: Vec<KmerNode>,
        kmer_to_id: HashMap<std::sync::Arc<[u8]>, usize>,
        edges: HashMap<(usize, usize), u32>,
        outgoing: HashMap<usize, BTreeSet<usize>>,
        incoming: HashMap<usize, BTreeSet<usize>>,
        ref_edges: HashSet<(usize, usize)>,
        ref_nodes: HashSet<usize>,
        ref_source_kmer: Option<std::sync::Arc<[u8]>>,
    ) -> Self {
        Self {
            kmer_size,
            dangling_merge_haps: Vec::new(),
            nodes,
            kmer_to_id,
            edges,
            outgoing,
            incoming,
            ref_edges,
            ref_nodes,
            ref_source_kmer,
        }
    }

    pub(crate) fn ensure_node(&mut self, kmer: &[u8]) -> usize {
        if let Some(id) = self.kmer_to_id.get(kmer) {
            let idx = *id;
            self.nodes[idx].support = self.nodes[idx].support.saturating_add(1);
            return idx;
        }
        let id = self.nodes.len();
        let owned: std::sync::Arc<[u8]> = std::sync::Arc::from(kmer);
        self.nodes.push(KmerNode {
            id,
            kmer: std::sync::Arc::clone(&owned),
            support: 1,
        });
        self.kmer_to_id.insert(owned, id);
        id
    }

    #[allow(dead_code)]
    fn add_edge(&mut self, from: usize, to: usize) {
        self.add_edge_support(from, to, 1);
    }

    pub(crate) fn add_edge_support(&mut self, from: usize, to: usize, support: u32) {
        let e = self.edges.entry((from, to)).or_insert(0);
        *e = e.saturating_add(support);
        self.outgoing.entry(from).or_default().insert(to);
        self.incoming.entry(to).or_default().insert(from);
    }

    pub fn node_count(&self) -> usize {
        self.nodes.len()
    }

    pub fn edge_count(&self) -> usize {
        self.edges.len()
    }

    /// Max out-degree over vertices (Peak bushiness gate).
    pub(crate) fn max_out_degree(&self) -> usize {
        self.outgoing.values().map(|s| s.len()).max().unwrap_or(0)
    }

    pub fn nodes(&self) -> &[KmerNode] {
        &self.nodes
    }

    pub(crate) fn vertex_id_for_kmer(&self, kmer: &[u8]) -> Option<usize> {
        self.kmer_to_id.get(kmer).copied()
    }

    pub fn edges_sorted(&self) -> Vec<KmerEdge> {
        let mut out = self
            .edges
            .iter()
            .map(|((from, to), support)| KmerEdge {
                from: *from,
                to: *to,
                support: *support,
            })
            .collect::<Vec<_>>();
        out.sort_by(|a, b| {
            let (af, at) = (&self.nodes[a.from].kmer, &self.nodes[a.to].kmer);
            let (bf, bt) = (&self.nodes[b.from].kmer, &self.nodes[b.to].kmer);
            af.cmp(bf).then(at.cmp(bt))
        });
        out
    }

    pub(crate) fn source_nodes(&self) -> Vec<usize> {
        (0..self.nodes.len())
            .filter(|&n| {
                self.incoming_count(n) == 0 && self.outgoing.get(&n).is_some_and(|o| !o.is_empty())
            })
            .collect()
    }

    pub(crate) fn outgoing_nodes(&self, from: usize) -> Vec<usize> {
        self.outgoing
            .get(&from)
            .map(|s| s.iter().copied().collect())
            .unwrap_or_default()
    }

    /// Borrow outgoing targets without allocating (cycle-strip / DFS hot path).
    pub(crate) fn outgoing_targets(&self, from: usize) -> impl Iterator<Item = usize> + '_ {
        self.outgoing
            .get(&from)
            .into_iter()
            .flat_map(|s| s.iter().copied())
    }

    pub(crate) fn incoming_nodes(&self, to: usize) -> Vec<usize> {
        self.incoming
            .get(&to)
            .map(|s| s.iter().copied().collect())
            .unwrap_or_default()
    }

    /// Borrow incoming sources without allocating.
    pub(crate) fn incoming_sources(&self, to: usize) -> impl Iterator<Item = usize> + '_ {
        self.incoming
            .get(&to)
            .into_iter()
            .flat_map(|s| s.iter().copied())
    }

    pub(crate) fn incoming_count(&self, node: usize) -> usize {
        self.incoming.get(&node).map(|s| s.len()).unwrap_or(0)
    }

    pub fn edge_is_ref(&self, from: usize, to: usize) -> bool {
        self.ref_edges.contains(&(from, to))
    }

    /// GATK `BaseGraph.isRefSource`.
    pub fn is_ref_source_vertex(&self, v: usize) -> bool {
        if self.node_count() == 1 {
            return true;
        }
        if self.incoming_sources(v).any(|p| self.edge_is_ref(p, v)) {
            return false;
        }
        self.outgoing_targets(v).any(|t| self.edge_is_ref(v, t))
    }

    /// GATK `BaseGraph.isRefSink`.
    pub fn is_ref_sink_vertex(&self, v: usize) -> bool {
        if self.node_count() == 1 {
            return true;
        }
        if self.outgoing_targets(v).any(|t| self.edge_is_ref(v, t)) {
            return false;
        }
        self.incoming_sources(v).any(|p| self.edge_is_ref(p, v))
    }

    pub fn has_ref_threading(&self) -> bool {
        !self.ref_nodes.is_empty()
    }

    pub fn reference_source_vertex(&self) -> Option<usize> {
        (0..self.node_count()).find(|&v| self.is_ref_source_vertex(v))
    }

    pub fn reference_sink_vertex(&self) -> Option<usize> {
        (0..self.node_count()).find(|&v| self.is_ref_sink_vertex(v))
    }

    /// GATK `Path.getBases` on a k-mer graph path.
    pub fn path_bases(&self, start: usize, edges: &[(usize, usize)]) -> Vec<u8> {
        if edges.is_empty() {
            return self.nodes[start].kmer.to_vec();
        }
        let first_from = edges[0].0;
        let mut s = self.nodes[first_from].kmer.to_vec();
        for &(_, to) in edges {
            if let Some(&b) = self.nodes[to].kmer.last() {
                s.push(b);
            }
        }
        s
    }

    pub fn has_cycle(&self) -> bool {
        self.has_cycle_internal()
    }

    /// GATK `BaseGraph.removePathsNotConnectedToRef`.
    pub fn remove_paths_not_connected_to_ref(&mut self) -> GatkResult<()> {
        let source = self
            .reference_source_vertex()
            .ok_or_else(|| GatkError::algorithm("removePathsNotConnectedToRef: no ref source"))?;
        let sink = self
            .reference_sink_vertex()
            .ok_or_else(|| GatkError::algorithm("removePathsNotConnectedToRef: no ref sink"))?;

        let mut from_source = HashSet::new();
        let mut stack = vec![source];
        from_source.insert(source);
        while let Some(v) = stack.pop() {
            for to in self.outgoing_nodes(v) {
                if from_source.insert(to) {
                    stack.push(to);
                }
            }
        }

        let mut from_sink = HashSet::new();
        let mut stack = vec![sink];
        from_sink.insert(sink);
        while let Some(v) = stack.pop() {
            for from in self.incoming_nodes(v) {
                if from_sink.insert(from) {
                    stack.push(from);
                }
            }
        }

        from_source.retain(|v| from_sink.contains(v));
        let remove: HashSet<usize> = (0..self.node_count())
            .filter(|v| !from_source.contains(v))
            .collect();
        self.remove_nodes(&remove);
        self.cleanup_isolated_nodes();
        Ok(())
    }

    pub(crate) fn edge_support(&self, from: usize, to: usize) -> Option<u32> {
        self.edges.get(&(from, to)).copied()
    }

    pub(crate) fn is_source(&self, node: usize) -> bool {
        self.incoming_count(node) == 0
    }

    pub(crate) fn is_sink(&self, node: usize) -> bool {
        self.outgoing
            .get(&node)
            .map(|s| s.is_empty())
            .unwrap_or(true)
    }

    pub(crate) fn remove_edge(&mut self, from: usize, to: usize) {
        self.edges.remove(&(from, to));
        if let Some(out) = self.outgoing.get_mut(&from) {
            out.remove(&to);
        }
        if let Some(inc) = self.incoming.get_mut(&to) {
            inc.remove(&from);
        }
    }

    pub(crate) fn remove_isolated_nodes(&mut self) {
        self.cleanup_isolated_nodes();
    }

    pub(crate) fn kmer_at(&self, id: usize) -> &[u8] {
        &self.nodes[id].kmer
    }

    /// Build one graph per k-mer size; graphs are independent de Bruijn views.
    pub fn from_reads_kmer_sizes(
        reads: &[AssemblyRead],
        kmer_sizes: &[usize],
        min_base_quality: u8,
    ) -> GatkResult<Vec<Self>> {
        let mut out = Vec::with_capacity(kmer_sizes.len());
        for &kmer_size in kmer_sizes {
            let kmer_size = crate::bio_ids::KmerSize::try_from_usize(kmer_size)?;
            let params = AssemblyGraphParams {
                kmer_size,
                min_base_quality,
                min_edge_weight: 1,
                dangling_path_max_nodes: 0,
                max_haplotypes: 128,
                max_haplotype_bases: 512,
                start_threading_only_at_existing_vertex: false,
            };
            out.push(Self::from_reads(reads, &params)?);
        }
        Ok(out)
    }

    pub fn from_reads(reads: &[AssemblyRead], params: &AssemblyGraphParams) -> GatkResult<Self> {
        crate::read_threading_graph::assembly_graph_from_reads_threading(reads, params)
    }

    pub fn prune_low_weight_edges(&mut self, min_edge_weight: u32) {
        let mut p = AssemblyGraphPruningParams::gatk_haplotype_caller_defaults();
        p.min_prune_factor = min_edge_weight;
        p.use_adaptive_pruning = false;
        let _ = self.apply_pruning(&p);
    }

    pub fn total_edge_support(&self) -> u64 {
        self.edges.values().map(|&s| s as u64).sum()
    }

    fn log10_max_edge_support(&self) -> f64 {
        self.edges
            .values()
            .map(|&s| (s.max(1) as f64).log10())
            .fold(f64::NEG_INFINITY, f64::max)
    }

    fn log10_sum_edge_support(&self) -> f64 {
        let sum = self.total_edge_support();
        if sum == 0 {
            f64::NEG_INFINITY
        } else {
            (sum as f64).log10()
        }
    }

    /// GATK `ChainPruner` (fixed or adaptive).
    pub fn apply_pruning(&mut self, params: &AssemblyGraphPruningParams) -> u32 {
        crate::assembly_pruning::apply_gatk_pruning(self, params)
    }

    /// Build graph summary after applying [`Self::apply_pruning`] + isolated-node cleanup.
    pub fn summarize_after_pruning(
        &mut self,
        params: &AssemblyGraphPruningParams,
    ) -> AssemblyGraphSummary {
        let edges_pruned = self.apply_pruning(params);
        self.cleanup_isolated_nodes();
        AssemblyGraphSummary {
            node_count: self.node_count(),
            edge_count: self.edge_count(),
            log10_max_edge_support: self.log10_max_edge_support(),
            log10_sum_edge_support: self.log10_sum_edge_support(),
            pruning_lod_threshold_ln: params.pruning_log_odds_threshold,
            adaptive_pruning: params.use_adaptive_pruning,
            min_prune_factor: params.min_prune_factor,
            edges_pruned,
        }
    }

    /// Remove short dead-end chains with only single-path continuation.
    pub fn remove_dangling_paths(&mut self, max_nodes: usize) {
        if max_nodes == 0 {
            return;
        }
        let mut remove_nodes = HashSet::new();
        for node_id in 0..self.nodes.len() {
            let out_degree = self.outgoing.get(&node_id).map(|x| x.len()).unwrap_or(0);
            if out_degree != 0 {
                continue;
            }
            let mut path = vec![node_id];
            let mut cur = node_id;
            while path.len() < max_nodes {
                // Lifetime: only need to inspect predecessor ids; borrow the set.
                let Some(preds) = self.incoming.get(&cur) else {
                    break;
                };
                if preds.len() != 1 {
                    break;
                }
                let Some(&pred) = preds.iter().next() else {
                    break;
                };
                let pred_out = self.outgoing.get(&pred).map(|x| x.len()).unwrap_or(0);
                if pred_out != 1 {
                    break;
                }
                path.push(pred);
                cur = pred;
            }
            if path.len() <= max_nodes {
                for n in path {
                    remove_nodes.insert(n);
                }
            }
        }
        self.remove_nodes(&remove_nodes);
    }

    pub fn cleanup_isolated_nodes(&mut self) {
        let mut keep = HashSet::new();
        for (from, to) in self.edges.keys() {
            keep.insert(*from);
            keep.insert(*to);
        }
        let remove = (0..self.nodes.len())
            .filter(|n| !keep.contains(n))
            .collect::<HashSet<_>>();
        self.remove_nodes(&remove);
    }

    /// One remapping pass: remove `remove` plus any nodes that would be isolated once those
    /// (and their incident edges) are gone. Used by cycle-strip instead of
    /// `remove_nodes` + `cleanup_isolated_nodes` (two full rebuilds).
    pub(crate) fn remove_nodes_and_isolated(&mut self, mut remove: HashSet<usize>) {
        let mut keep = HashSet::new();
        for (from, to) in self.edges.keys() {
            if !remove.contains(from) && !remove.contains(to) {
                keep.insert(*from);
                keep.insert(*to);
            }
        }
        for id in 0..self.nodes.len() {
            if !keep.contains(&id) {
                remove.insert(id);
            }
        }
        self.remove_nodes(&remove);
    }

    pub(crate) fn remove_nodes(&mut self, remove: &HashSet<usize>) {
        if remove.is_empty() {
            return;
        }
        let old_nodes = std::mem::take(&mut self.nodes);
        let old_edges = std::mem::take(&mut self.edges);
        let old_ref_edges = std::mem::take(&mut self.ref_edges);
        let mut new_nodes = Vec::with_capacity(old_nodes.len().saturating_sub(remove.len()));
        let mut remap = HashMap::with_capacity(new_nodes.capacity());
        self.kmer_to_id.clear();
        for n in old_nodes {
            if !remove.contains(&n.id) {
                let new_id = new_nodes.len();
                remap.insert(n.id, new_id);
                // One Arc bump for the index; move the existing Arc into the node (no second clone loop).
                self.kmer_to_id
                    .insert(std::sync::Arc::clone(&n.kmer), new_id);
                new_nodes.push(KmerNode {
                    id: new_id,
                    kmer: n.kmer,
                    support: n.support,
                });
            }
        }
        self.edges.clear();
        self.ref_edges.clear();
        self.outgoing.clear();
        self.incoming.clear();
        self.edges.reserve(old_edges.len());
        for ((from, to), support) in old_edges {
            if let (Some(&nf), Some(&nt)) = (remap.get(&from), remap.get(&to)) {
                self.edges.insert((nf, nt), support);
                if old_ref_edges.contains(&(from, to)) {
                    self.ref_edges.insert((nf, nt));
                }
                self.outgoing.entry(nf).or_default().insert(nt);
                self.incoming.entry(nt).or_default().insert(nf);
            }
        }
        let new_ref_nodes: HashSet<usize> = self
            .ref_nodes
            .iter()
            .filter_map(|n| remap.get(n).copied())
            .collect();
        self.nodes = new_nodes;
        self.ref_nodes = new_ref_nodes;
    }

    pub fn extract_candidate_haplotypes(
        &self,
        max_haplotypes: usize,
        max_bases: usize,
    ) -> Vec<CandidateHaplotype> {
        if self.nodes.is_empty() || max_haplotypes == 0 || max_bases == 0 {
            return Vec::new();
        }
        let mut starts = (0..self.nodes.len())
            .filter(|id| self.incoming.get(id).map(|x| x.is_empty()).unwrap_or(true))
            .collect::<Vec<_>>();
        if starts.is_empty() {
            starts = (0..self.nodes.len()).collect();
        }
        starts.sort_unstable();

        let mut heap: BinaryHeap<(u32, Reverse<Vec<u8>>)> = BinaryHeap::new();
        for s in starts {
            let mut visited = HashSet::new();
            let mut current = self.nodes[s].kmer.to_vec(); // owned scratch for DFS path growth
            self.dfs_haplotypes(s, &mut current, 0, max_bases, &mut visited, &mut heap);
        }
        let mut out = heap
            .into_iter()
            .map(|(support, Reverse(sequence))| CandidateHaplotype { sequence, support })
            .collect::<Vec<_>>();
        out.sort_by(|a, b| {
            b.support
                .cmp(&a.support)
                .then_with(|| a.sequence.cmp(&b.sequence))
        });
        out.truncate(max_haplotypes);
        out
    }

    /// Basic graph invariants used by phase-5 contracts.
    pub fn validate_basic_invariants(&self) -> GatkResult<()> {
        for n in &self.nodes {
            if n.kmer.len() != self.kmer_size {
                return Err(GatkError::algorithm(format!(
                    "node {} has kmer length {}, expected {}",
                    n.id,
                    n.kmer.len(),
                    self.kmer_size
                )));
            }
        }
        for ((from, to), support) in &self.edges {
            if *support == 0 {
                return Err(GatkError::algorithm(format!(
                    "edge ({from}->{to}) has zero support"
                )));
            }
            if *from >= self.nodes.len() || *to >= self.nodes.len() {
                return Err(GatkError::algorithm(format!(
                    "edge ({from}->{to}) points outside node set"
                )));
            }
            let from_k = &self.nodes[*from].kmer;
            let to_k = &self.nodes[*to].kmer;
            if !from_k[1..].eq(&to_k[..self.kmer_size - 1]) {
                return Err(GatkError::algorithm(format!(
                    "edge ({from}->{to}) violates k-1 overlap invariant"
                )));
            }
        }
        if self.has_cycle_internal() {
            return Err(GatkError::algorithm(
                "graph contains a cycle; expected DAG after cleanup for current wave-B fixtures",
            ));
        }
        Ok(())
    }

    fn has_cycle_internal(&self) -> bool {
        let mut indeg = vec![0usize; self.nodes.len()];
        for (_, to) in self.edges.keys() {
            indeg[*to] += 1;
        }
        let mut stack = indeg
            .iter()
            .enumerate()
            .filter_map(|(i, d)| if *d == 0 { Some(i) } else { None })
            .collect::<Vec<_>>();
        let mut seen = 0usize;
        while let Some(n) = stack.pop() {
            seen += 1;
            if let Some(outs) = self.outgoing.get(&n) {
                for to in outs {
                    indeg[*to] -= 1;
                    if indeg[*to] == 0 {
                        stack.push(*to);
                    }
                }
            }
        }
        seen != self.nodes.len()
    }

    fn dfs_haplotypes(
        &self,
        node: usize,
        current: &mut Vec<u8>,
        support_acc: u32,
        max_bases: usize,
        visited: &mut HashSet<usize>,
        out: &mut BinaryHeap<(u32, Reverse<Vec<u8>>)>,
    ) {
        if current.len() > max_bases || visited.contains(&node) {
            return;
        }
        visited.insert(node);
        // Lifetime: DFS only reads adjacency; iterate borrowed successors.
        let outs = self.outgoing.get(&node);
        if outs.map(|s| s.is_empty()).unwrap_or(true) {
            // CLONE: needed — heap owns completed haplotype sequences.
            out.push((support_acc.max(1), Reverse(current.clone())));
            visited.remove(&node);
            return;
        }
        let Some(outs) = outs else {
            visited.remove(&node);
            return;
        };
        for &to in outs {
            let edge_w = *self.edges.get(&(node, to)).unwrap_or(&0);
            let next_base = self.nodes[to].kmer.last().copied().unwrap_or(b'N');
            current.push(next_base);
            self.dfs_haplotypes(
                to,
                current,
                support_acc.saturating_add(edge_w),
                max_bases,
                visited,
                out,
            );
            current.pop();
        }
        visited.remove(&node);
    }
}

/// Wave-A graph-first assembly API (not on production `call_region` path).
/// Prefer [`crate::read_threading_assembler`] / [`crate::assembly_based_caller`].
/// Exercised by `#[cfg(test)]` unit tests; kept for that contract outside test cfg.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub(crate) struct AssemblyEngine {
    params: AssemblyGraphParams,
}

#[allow(dead_code)] // see struct note — unit-tested under cfg(test)
impl AssemblyEngine {
    /// Create a new assembly engine.
    pub fn new(params: AssemblyGraphParams) -> GatkResult<Self> {
        Ok(Self { params })
    }

    /// Assemble candidate haplotypes from reads (delegates to [`crate::read_threading_assembler`] when reference provided).
    pub fn assemble(&self, reads: &[AssemblyRead]) -> GatkResult<Vec<CandidateHaplotype>> {
        self.assemble_with_optional_ref(None, reads)
    }

    /// Assemble with optional reference sequence (GATK ref-threaded path).
    pub fn assemble_with_optional_ref(
        &self,
        reference: Option<&AssemblyRead>,
        reads: &[AssemblyRead],
    ) -> GatkResult<Vec<CandidateHaplotype>> {
        if let Some(reference) = reference {
            let mut args = crate::read_threading_assembler::ReadThreadingAssemblerArgs::default();
            args.kmer_sizes = vec![self.params.kmer_size.as_usize()];
            args.min_base_quality = self.params.min_base_quality;
            args.min_prune_factor = self.params.min_edge_weight;
            args.num_best_haplotypes_per_graph = self.params.max_haplotypes;
            let result = crate::read_threading_assembler::assemble_from_ref_and_reads(
                reference, reads, &args,
            )?;
            return Ok(result
                .haplotypes
                .into_iter()
                .map(|h| CandidateHaplotype {
                    sequence: h.bases,
                    support: h.score.max(0.0) as u32,
                })
                .collect());
        }
        let mut graph = AssemblyGraph::from_reads(reads, &self.params)?;
        let mut pruning = AssemblyGraphPruningParams::gatk_haplotype_caller_defaults();
        pruning.min_prune_factor = self.params.min_edge_weight;
        graph.apply_pruning(&pruning);
        graph.remove_dangling_paths(self.params.dangling_path_max_nodes);
        graph.cleanup_isolated_nodes();
        Ok(graph.extract_candidate_haplotypes(
            self.params.max_haplotypes,
            self.params.max_haplotype_bases,
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn mk_read(seq: &str, q: u8) -> AssemblyRead {
        AssemblyRead {
            bases: seq.as_bytes().to_vec(),
            base_quals: vec![q; seq.len()],
        }
    }

    #[test]
    fn fixed_pruning_drops_low_support_edges() {
        let reads = vec![
            mk_read("ACGTT", 30),
            mk_read("ACGTT", 30),
            mk_read("ACGTA", 30),
        ];
        let p = AssemblyGraphParams {
            kmer_size: KmerSize::try_new(3).unwrap(),
            min_base_quality: 10,
            ..Default::default()
        };
        let mut g = AssemblyGraph::from_reads(&reads, &p).unwrap();
        let mut prune = AssemblyGraphPruningParams::gatk_haplotype_caller_defaults();
        prune.min_prune_factor = 3;
        let removed = g.apply_pruning(&prune);
        assert!(removed >= 1);
        assert!(g.edges_sorted().iter().all(|e| e.support >= 3));
    }

    #[test]
    fn adaptive_pruning_keeps_branch_after_threading() {
        let reads = vec![
            mk_read("ACGTT", 30),
            mk_read("ACGTT", 30),
            mk_read("ACGTT", 30),
            mk_read("ACGTA", 30),
            mk_read("ACGTA", 30),
        ];
        let p = AssemblyGraphParams {
            kmer_size: KmerSize::try_new(3).unwrap(),
            min_base_quality: 10,
            ..Default::default()
        };
        let mut g = AssemblyGraph::from_reads(&reads, &p).unwrap();
        let mut prune = AssemblyGraphPruningParams::gatk_haplotype_caller_defaults();
        prune.min_prune_factor = 4;
        prune.use_adaptive_pruning = true;
        g.apply_pruning(&prune);
        assert!(g.edge_count() >= 3);
    }

    #[test]
    fn summarize_after_pruning_emits_log10_mass() {
        let reads = vec![
            mk_read("ACGTT", 30),
            mk_read("ACGTT", 30),
            mk_read("ACGTA", 30),
        ];
        let p = AssemblyGraphParams {
            kmer_size: KmerSize::try_new(3).unwrap(),
            min_base_quality: 10,
            ..Default::default()
        };
        let mut g = AssemblyGraph::from_reads(&reads, &p).unwrap();
        let mut prune = AssemblyGraphPruningParams::gatk_haplotype_caller_defaults();
        prune.min_prune_factor = 3;
        let s = g.summarize_after_pruning(&prune);
        assert!(s.log10_max_edge_support.is_finite());
        assert!(s.edges_pruned >= 1);
    }

    #[test]
    fn multi_kmer_graphs_produce_independent_edge_sets() {
        let reads = vec![mk_read("ACGTTACGT", 30), mk_read("ACGTAACGA", 30)];
        let graphs = AssemblyGraph::from_reads_kmer_sizes(&reads, &[3, 5], 10).unwrap();
        assert_eq!(graphs.len(), 2);
        assert_eq!(graphs[0].kmer_size, 3);
        assert_eq!(graphs[1].kmer_size, 5);
        assert!(graphs[0].edge_count() >= 2);
        assert!(graphs[1].edge_count() >= 1);
    }

    #[test]
    fn graph_construction_from_reads_builds_kmer_nodes_and_edges() {
        let p = AssemblyGraphParams {
            kmer_size: KmerSize::try_new(3).unwrap(),
            min_base_quality: 10,
            ..Default::default()
        };
        let g = AssemblyGraph::from_reads(&[mk_read("ACGT", 30), mk_read("ACGA", 30)], &p).unwrap();
        assert!(g.node_count() >= 3);
        assert!(g.edge_count() >= 2);
    }

    #[test]
    fn quality_filter_discards_low_quality_windows() {
        let p = AssemblyGraphParams {
            kmer_size: KmerSize::try_new(3).unwrap(),
            min_base_quality: 20,
            ..Default::default()
        };
        let mut low = mk_read("ACGT", 30);
        low.base_quals[1] = 5;
        let g = AssemblyGraph::from_reads(&[low], &p).unwrap();
        assert_eq!(g.node_count(), 0);
        assert_eq!(g.edge_count(), 0);
    }

    #[test]
    fn pruning_and_dangling_removal_reduce_noise_paths() {
        let p = AssemblyGraphParams {
            kmer_size: KmerSize::try_new(3).unwrap(),
            min_base_quality: 10,
            min_edge_weight: 2,
            dangling_path_max_nodes: 3,
            ..Default::default()
        };
        let reads = vec![
            mk_read("ACGTT", 30),
            mk_read("ACGTT", 30),
            mk_read("ACGTA", 30), // low-support branch
        ];
        let mut g = AssemblyGraph::from_reads(&reads, &p).unwrap();
        let before = g.edge_count();
        g.prune_low_weight_edges(p.min_edge_weight);
        g.remove_dangling_paths(p.dangling_path_max_nodes);
        g.cleanup_isolated_nodes();
        assert!(g.edge_count() < before);
    }

    #[test]
    fn candidate_extraction_is_deterministic_and_sorted() {
        let p = AssemblyGraphParams {
            kmer_size: KmerSize::try_new(3).unwrap(),
            min_base_quality: 10,
            min_edge_weight: 1,
            max_haplotypes: 8,
            max_haplotype_bases: 32,
            ..Default::default()
        };
        let engine = AssemblyEngine::new(p).unwrap();
        let reads = vec![
            mk_read("ACGTT", 30),
            mk_read("ACGTT", 30),
            mk_read("ACGTA", 30),
        ];
        let mut reads_rev = reads.clone();
        reads_rev.reverse();
        let h1 = engine.assemble(&reads).unwrap();
        let h2 = engine.assemble(&reads_rev).unwrap();
        assert_eq!(h1, h2);
        assert!(!h1.is_empty());
        for w in h1.windows(2) {
            assert!(
                w[0].support > w[1].support
                    || (w[0].support == w[1].support && w[0].sequence <= w[1].sequence)
            );
        }
    }

    #[test]
    fn graph_invariants_hold_after_wave_a_pipeline() {
        let p = AssemblyGraphParams {
            kmer_size: KmerSize::try_new(4).unwrap(),
            min_base_quality: 10,
            min_edge_weight: 1,
            dangling_path_max_nodes: 3,
            ..Default::default()
        };
        let reads = vec![
            mk_read("ACGTTAC", 30),
            mk_read("ACGTTAC", 30),
            mk_read("ACGTAAC", 30),
        ];
        let mut g = AssemblyGraph::from_reads(&reads, &p).unwrap();
        g.prune_low_weight_edges(p.min_edge_weight);
        g.remove_dangling_paths(p.dangling_path_max_nodes);
        g.cleanup_isolated_nodes();
        g.validate_basic_invariants().unwrap();
    }

    #[test]
    fn homopolymer_repeat_region_produces_stable_candidates() {
        let p = AssemblyGraphParams {
            kmer_size: KmerSize::try_new(4).unwrap(),
            min_base_quality: 10,
            min_edge_weight: 1,
            dangling_path_max_nodes: 0,
            max_haplotypes: 8,
            ..Default::default()
        };
        let engine = AssemblyEngine::new(p).unwrap();
        let reads = vec![
            mk_read("AAAACAAA", 30),
            mk_read("AAAACAAA", 30),
            mk_read("AAAAGAAA", 30),
            mk_read("AAAAGAAA", 30),
        ];
        let hs = engine.assemble(&reads).unwrap();
        assert!(hs
            .iter()
            .any(|h| h.sequence.windows(5).any(|w| w == b"AAAAC")));
        assert!(hs
            .iter()
            .any(|h| h.sequence.windows(5).any(|w| w == b"AAAAG")));
    }

    #[test]
    fn alt_path_recovery_keeps_supported_minor_branch() {
        let p = AssemblyGraphParams {
            kmer_size: KmerSize::try_new(3).unwrap(),
            min_base_quality: 10,
            min_edge_weight: 2,
            dangling_path_max_nodes: 0,
            max_haplotypes: 16,
            ..Default::default()
        };
        let engine = AssemblyEngine::new(p).unwrap();
        let reads = vec![
            mk_read("ACGTT", 30),
            mk_read("ACGTT", 30),
            mk_read("ACGTT", 30),
            mk_read("ACGTA", 30),
            mk_read("ACGTA", 30),
        ];
        let hs = engine.assemble(&reads).unwrap();
        // Both branches are above min_edge_weight and should survive cleanup.
        assert!(hs.iter().any(|h| h.sequence.ends_with(b"T")));
        assert!(hs.iter().any(|h| h.sequence.ends_with(b"A")));
    }

    #[test]
    fn tie_break_is_lexicographic_for_equal_support() {
        let p = AssemblyGraphParams {
            kmer_size: KmerSize::try_new(3).unwrap(),
            min_base_quality: 10,
            min_edge_weight: 1,
            dangling_path_max_nodes: 0,
            max_haplotypes: 16,
            ..Default::default()
        };
        let engine = AssemblyEngine::new(p).unwrap();
        let reads = vec![
            mk_read("ACGTA", 30),
            mk_read("ACGTA", 30),
            mk_read("ACGTC", 30),
            mk_read("ACGTC", 30),
        ];
        let hs = engine.assemble(&reads).unwrap();
        let eq = hs
            .windows(2)
            .find(|w| w[0].support == w[1].support)
            .expect("equal-support pair");
        assert!(eq[0].sequence <= eq[1].sequence);
    }

    #[test]
    fn high_depth_memory_pressure_contract_keeps_outputs_bounded() {
        let p = AssemblyGraphParams {
            kmer_size: KmerSize::try_new(5).unwrap(),
            min_base_quality: 10,
            min_edge_weight: 2,
            dangling_path_max_nodes: 4,
            max_haplotypes: 32,
            max_haplotype_bases: 120,
            start_threading_only_at_existing_vertex: false,
        };
        let engine = AssemblyEngine::new(p).unwrap();
        let mut reads = Vec::new();
        for i in 0..5_000 {
            let mut seq = String::from("ACGTACGT");
            if i % 7 == 0 {
                seq.push_str("TTAA");
            } else if i % 7 == 1 {
                seq.push_str("TCAA");
            } else {
                seq.push_str("TGAA");
            }
            seq.push_str("AAAAAAAAAAAAAAAAAAAA");
            reads.push(AssemblyRead {
                base_quals: vec![30; seq.len()],
                bases: seq.into_bytes(),
            });
        }
        let hs = engine.assemble(&reads).unwrap();
        assert!(hs.len() <= p.max_haplotypes);
        assert!(hs.iter().all(|h| h.sequence.len() <= p.max_haplotype_bases));
    }
}
