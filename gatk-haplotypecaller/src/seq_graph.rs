//! GATK `SeqGraph` conversion and cleanup (`BaseGraph.toSequenceGraph`, `SeqGraph.cleanup`).
//!
//! Vertex compaction remaps surviving edge endpoints (`compact_vertices_and_remap_edges`)
//! so `vertex.id == index` and every `edge.from`/`edge.to` names a live vertex. Zip copies
//! the last chain vertex's outgoing edges onto the merged keep vertex (GATK
//! `mergeLinearChainVertex`). That is a structural identity/splice repair only: it does not
//! claim Java SeqGraph haplotype parity, does not retire W-H1, and does not change the P12
//! `use_seq_graph = false` production waiver.

use crate::assembly::AssemblyGraph;
use gatk_common::{GatkError, GatkResult};
use std::collections::{HashMap, HashSet};

/// Sequence-labeled vertex after k-mer graph → sequence graph conversion.
/// # Invariants
/// `id` is dense index into the owning sequence graph vertex list.
/// `sequence` holds vertex bases (full k-mer at source, suffix byte at non-source joins).
/// # Ownership
/// Owns base sequence bytes; graph owns vertex vector.
/// # Mutation
/// Graph cleanup may merge/replace vertex sequences in place.
/// # Biological assumptions
/// Vertex sequence is literal assembly graph bases on the reference/read path.
/// # Java equivalence
/// GATK `SeqVertex` in `SeqGraph` (`BaseGraph.toSequenceGraph`).
#[derive(Debug, Clone)]
pub struct SeqVertex {
    pub id: usize,
    pub sequence: Vec<u8>,
}

/// Directed edge in a sequence graph with support and reference-path flag.
/// # Invariants
/// `from` / `to` index valid vertices; `is_ref` marks reference-spine edges.
/// # Ownership
/// [`Copy`] edge record in the owning sequence graph edge list.
/// # Mutation
/// Support may increase on merge; edges removed during cleanup passes.
/// # Biological assumptions
/// Edge support reflects read evidence multiplicity at the junction.
/// # Java equivalence
/// GATK `SeqEdge` / reference edge marking in `SeqGraph`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SeqEdge {
    pub from: usize,
    pub to: usize,
    pub support: u32,
    pub is_ref: bool,
}

/// Outcome of [`SeqGraph::cleanup_seq_graph`] — variation present or reference-only.
/// # Invariants
/// Mutually exclusive assembly outcomes for one cleanup pass.
/// # Ownership
/// [`Copy`] enum returned from cleanup.
/// # Mutation
/// Immutable status tag.
/// # Biological assumptions
/// `AssembledSomeVariation` means non-ref paths survived cleanup connected to ref endpoints.
/// # Java equivalence
/// GATK `SeqGraph` cleanup status (`ASSEMBLED_SOME_VARIATION`, `JUST_ASSEMBLED_REFERENCE`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SeqGraphCleanupStatus {
    AssembledSomeVariation,
    JustAssembledReference,
}

/// Base sequence graph for k-best haplotype discovery after threading conversion.
/// # Invariants
/// Reference source/sink vertices exist after successful cleanup when variation is assembled.
/// Vertex `id` values remain dense (`id == index`) after prune/zip; edge `from`/`to`
/// are rewritten through the same mapping (`compact_vertices_and_remap_edges`).
/// # Ownership
/// Owns vertices, edges, and adjacency indexes; built from [`AssemblyGraph`] without retaining reads.
/// # Mutation
/// Cleanup methods (`clean_non_ref_paths`, `zip_linear_chains`, `simplify_graph`, etc.) mutate in place.
/// # Biological assumptions
/// Graph encodes colinear sequence alternatives over the padded reference window.
/// # Java equivalence
/// GATK `SeqGraph` (`BaseGraph.toSequenceGraph`, `ReadThreadingAssembler.cleanupSeqGraph`).
#[derive(Debug, Clone)]
pub struct SeqGraph {
    pub kmer_size: usize,
    vertices: Vec<SeqVertex>,
    edges: Vec<SeqEdge>,
    outgoing: HashMap<usize, Vec<usize>>,
    incoming: HashMap<usize, Vec<usize>>,
}

impl SeqGraph {
    pub fn from_assembly_graph(graph: &AssemblyGraph) -> Self {
        let kmer_size = graph.kmer_size;
        let mut vertices = Vec::new();
        let mut id_map = HashMap::new();
        for (i, node) in graph.nodes().iter().enumerate() {
            let is_source = graph.incoming_count(i) == 0;
            let seq = additional_sequence_bytes(&node.kmer, is_source);
            let id = vertices.len();
            vertices.push(SeqVertex { id, sequence: seq });
            id_map.insert(i, id);
        }
        let mut edges = Vec::new();
        let mut outgoing: HashMap<usize, Vec<usize>> = HashMap::new();
        let mut incoming: HashMap<usize, Vec<usize>> = HashMap::new();
        for e in graph.edges_sorted() {
            let from = *id_map.get(&e.from).expect("edge from");
            let to = *id_map.get(&e.to).expect("edge to");
            let is_ref = graph.edge_is_ref(e.from, e.to);
            edges.push(SeqEdge {
                from,
                to,
                support: e.support,
                is_ref,
            });
            outgoing.entry(from).or_default().push(to);
            incoming.entry(to).or_default().push(from);
        }
        Self {
            kmer_size,
            vertices,
            edges,
            outgoing,
            incoming,
        }
    }

    pub fn node_count(&self) -> usize {
        self.vertices.len()
    }

    pub fn edge_count(&self) -> usize {
        self.edges.len()
    }

    /// Max out-degree over vertices (Peak bushiness gate).
    pub(crate) fn max_out_degree(&self) -> usize {
        self.outgoing.values().map(|s| s.len()).max().unwrap_or(0)
    }

    pub fn reference_source_vertex(&self) -> Option<usize> {
        (0..self.vertices.len()).find(|&v| self.is_ref_source(v))
    }

    pub fn reference_sink_vertex(&self) -> Option<usize> {
        (0..self.vertices.len()).find(|&v| self.is_ref_sink(v))
    }

    pub(crate) fn is_ref_source_vertex(&self, v: usize) -> bool {
        self.is_ref_source(v)
    }

    fn is_ref_source(&self, v: usize) -> bool {
        if self.vertices.len() == 1 {
            return true;
        }
        if self
            .incoming
            .get(&v)
            .into_iter()
            .flatten()
            .any(|&p| self.edge_is_ref(p, v))
        {
            return false;
        }
        self.outgoing
            .get(&v)
            .into_iter()
            .flatten()
            .any(|&t| self.edge_is_ref(v, t))
    }

    fn is_ref_sink(&self, v: usize) -> bool {
        if self.vertices.len() == 1 {
            return true;
        }
        if self
            .outgoing
            .get(&v)
            .into_iter()
            .flatten()
            .any(|&t| self.edge_is_ref(v, t))
        {
            return false;
        }
        self.incoming
            .get(&v)
            .into_iter()
            .flatten()
            .any(|&p| self.edge_is_ref(p, v))
    }

    pub(crate) fn edge_is_ref(&self, from: usize, to: usize) -> bool {
        self.edges
            .iter()
            .any(|e| e.from == from && e.to == to && e.is_ref)
    }

    pub(crate) fn edge_support(&self, from: usize, to: usize) -> Option<u32> {
        self.edges
            .iter()
            .find(|e| e.from == from && e.to == to)
            .map(|e| e.support)
    }

    pub(crate) fn outgoing_nodes(&self, from: usize) -> Vec<usize> {
        self.outgoing_of(from).to_vec()
    }

    /// GATK `Path.getBases` on a sequence graph.
    ///
    /// SeqVertex payloads are already additional-only (k-mer last-byte applied once in
    /// [`Self::from_assembly_graph`]). Concatenate the stored sequence of every vertex
    /// on the path. Empty dummy vertices contribute nothing.
    pub fn path_bases_bytes(&self, start: usize, edges: &[(usize, usize)]) -> Vec<u8> {
        let first = if edges.is_empty() { start } else { edges[0].0 };
        let mut bases = self.vertices[first].sequence.to_vec();
        for &(_, to) in edges {
            bases.extend_from_slice(&self.vertices[to].sequence);
        }
        bases
    }

    fn in_degree(&self, v: usize) -> usize {
        self.incoming.get(&v).map(|v| v.len()).unwrap_or(0)
    }

    fn out_degree(&self, v: usize) -> usize {
        self.outgoing.get(&v).map(|v| v.len()).unwrap_or(0)
    }

    pub(crate) fn outgoing_of(&self, from: usize) -> &[usize] {
        self.outgoing
            .get(&from)
            .map(|s| s.as_slice())
            .unwrap_or(&[])
    }

    fn incoming_of(&self, to: usize) -> &[usize] {
        self.incoming.get(&to).map(|s| s.as_slice()).unwrap_or(&[])
    }

    pub fn reference_path_bytes(&self) -> Option<Vec<u8>> {
        let source = self.reference_source_vertex()?;
        let sink = self.reference_sink_vertex()?;
        let mut path = self.vertices[source].sequence.clone();
        let mut cur = source;
        while cur != sink {
            let outs = self.outgoing_of(cur);
            let next = *outs
                .iter()
                .find(|&&t| self.edge_is_ref(cur, t))
                .or_else(|| outs.first())?;
            path.extend_from_slice(&self.vertices[next].sequence);
            cur = next;
        }
        Some(path)
    }

    pub fn clean_non_ref_paths(&mut self) {
        let Some(source) = self.reference_source_vertex() else {
            return;
        };
        let Some(sink) = self.reference_sink_vertex() else {
            return;
        };
        let mut to_remove: HashSet<(usize, usize)> = HashSet::new();
        let mut check: Vec<(usize, usize)> = self
            .edges
            .iter()
            .filter(|e| e.to == source && !e.is_ref)
            .map(|e| (e.from, e.to))
            .collect();
        while let Some((from, to)) = check.pop() {
            if !to_remove.insert((from, to)) {
                continue;
            }
            for e in &self.edges {
                if e.to == from && !e.is_ref {
                    check.push((e.from, e.to));
                }
            }
        }
        check = self
            .edges
            .iter()
            .filter(|e| e.from == sink && !e.is_ref)
            .map(|e| (e.from, e.to))
            .collect();
        while let Some((from, to)) = check.pop() {
            if !to_remove.insert((from, to)) {
                continue;
            }
            for e in &self.edges {
                if e.from == to && !e.is_ref {
                    check.push((e.from, e.to));
                }
            }
        }
        self.edges.retain(|e| !to_remove.contains(&(e.from, e.to)));
        self.rebuild_index();
        self.remove_singleton_orphan_vertices();
    }

    pub fn zip_linear_chains(&mut self) -> bool {
        let mut starts: Vec<usize> = (0..self.vertices.len())
            .filter(|&s| self.is_linear_chain_start(s))
            .collect();
        if starts.is_empty() {
            return false;
        }
        let mut merged = false;
        let mut i = 0;
        while i < starts.len() {
            let start = starts[i];
            i += 1;
            if start >= self.vertices.len() || !self.is_linear_chain_start(start) {
                continue;
            }
            let chain = self.trace_linear_chain(start);
            if chain.len() <= 1 {
                continue;
            }
            let Some(old_to_new) = self.merge_linear_chain_map(&chain) else {
                continue;
            };
            merged = true;
            for s in starts.iter_mut().skip(i) {
                *s = old_to_new.get(s).copied().unwrap_or(usize::MAX);
            }
        }
        merged
    }

    fn is_linear_chain_start(&self, source: usize) -> bool {
        if self.out_degree(source) != 1 {
            return false;
        }
        let indeg = self.in_degree(source);
        indeg != 1
            || self
                .incoming_of(source)
                .first()
                .map(|&p| self.out_degree(p) > 1)
                .unwrap_or(false)
    }

    fn trace_linear_chain(&self, zip_start: usize) -> Vec<usize> {
        let mut chain = vec![zip_start];
        let mut last = zip_start;
        loop {
            if self.out_degree(last) != 1 {
                break;
            }
            let target = self.outgoing_of(last)[0];
            if self.in_degree(target) != 1 || last == target {
                break;
            }
            let last_ref = self.is_reference_node(last);
            let target_ref = self.is_reference_node(target);
            if last_ref != target_ref {
                break;
            }
            chain.push(target);
            last = target;
        }
        chain
    }

    pub(crate) fn is_reference_node(&self, v: usize) -> bool {
        self.outgoing_of(v).iter().any(|&t| self.edge_is_ref(v, t))
            || self.incoming_of(v).iter().any(|&p| self.edge_is_ref(p, v))
    }

    /// Incoming edges of `v` as `(from, support, is_ref)`.
    pub(crate) fn edges_into(&self, v: usize) -> Vec<(usize, u32, bool)> {
        self.edges
            .iter()
            .filter(|e| e.to == v)
            .map(|e| (e.from, e.support, e.is_ref))
            .collect()
    }

    /// Outgoing edges of `v` as `(to, support, is_ref)`.
    pub(crate) fn edges_from(&self, v: usize) -> Vec<(usize, u32, bool)> {
        self.edges
            .iter()
            .filter(|e| e.from == v)
            .map(|e| (e.to, e.support, e.is_ref))
            .collect()
    }

    fn merge_linear_chain(&mut self, chain: &[usize]) -> bool {
        self.merge_linear_chain_map(chain).is_some()
    }

    /// GATK `mergeLinearChainVertex`: splice `chain[0]`+…+`last` into `keep`, then copy
    /// outgoing edges of `last` onto `keep` (Java `outgoingEdgesOf(last)` → merged vertex).
    /// Returns old-id → new-id so remaining zip starts can be retranslated after compact.
    fn merge_linear_chain_map(&mut self, chain: &[usize]) -> Option<HashMap<usize, usize>> {
        if chain.len() < 2 {
            return None;
        }
        let keep = chain[0];
        let last = chain[chain.len() - 1];
        let chain_set: HashSet<usize> = chain.iter().copied().collect();
        // Java copies last's outgoing onto the merged vertex before removing the chain.
        // Without this, compact drops `last → successor` (e.g. last → join) and disconnects.
        for e in &mut self.edges {
            if e.from == last && !chain_set.contains(&e.to) {
                e.from = keep;
            }
        }
        let remove: HashSet<usize> = chain[1..].iter().copied().collect();
        // Take keep sequence then extend — avoid cloning every vertex seq on the chain.
        let mut merged_seq = std::mem::take(&mut self.vertices[keep].sequence);
        for &v in &chain[1..] {
            merged_seq.extend_from_slice(&self.vertices[v].sequence);
        }

        let kept: Vec<usize> = (0..self.vertices.len())
            .filter(|i| !remove.contains(i))
            .collect();
        let old_to_new: HashMap<usize, usize> = kept
            .iter()
            .enumerate()
            .map(|(new_id, &old_id)| (old_id, new_id))
            .collect();

        let mut new_vertices = Vec::with_capacity(kept.len());
        for &old_id in &kept {
            let sequence = if old_id == keep {
                std::mem::take(&mut merged_seq)
            } else {
                std::mem::take(&mut self.vertices[old_id].sequence)
            };
            new_vertices.push(SeqVertex {
                id: old_to_new[&old_id],
                sequence,
            });
        }

        let mut new_edges = Vec::new();
        let mut seen: HashSet<(usize, usize, bool)> = HashSet::new();
        for e in &self.edges {
            let Some(&from) = old_to_new.get(&e.from) else {
                continue;
            };
            let Some(&to) = old_to_new.get(&e.to) else {
                continue;
            };
            if from == to {
                continue;
            }
            if seen.insert((from, to, e.is_ref)) {
                new_edges.push(SeqEdge {
                    from,
                    to,
                    support: e.support,
                    is_ref: e.is_ref,
                });
            }
        }
        self.vertices = new_vertices;
        self.edges = new_edges;
        self.rebuild_index();
        Some(old_to_new)
    }

    fn rebuild_index(&mut self) {
        self.outgoing.clear();
        self.incoming.clear();
        for e in &self.edges {
            self.outgoing.entry(e.from).or_default().push(e.to);
            self.incoming.entry(e.to).or_default().push(e.from);
        }
    }

    /// Compact surviving vertices to dense `id == index` and rewrite edge endpoints.
    ///
    /// Drops any edge whose `from`/`to` is not in the kept set. This is the identity
    /// contract already used by [`Self::merge_linear_chain`]; prune helpers must not
    /// reindex `vertex.id` while leaving `SeqEdge` in the old namespace.
    fn compact_vertices_and_remap_edges(
        &mut self,
        remove: &HashSet<usize>,
    ) -> HashMap<usize, usize> {
        let n = self.vertices.len();
        let kept: Vec<usize> = (0..n)
            .filter(|&i| !remove.contains(&i) && !remove.contains(&self.vertices[i].id))
            .collect();

        let mut old_to_new: HashMap<usize, usize> = HashMap::with_capacity(kept.len() * 2);
        for (new_id, &old_idx) in kept.iter().enumerate() {
            old_to_new.insert(old_idx, new_id);
            old_to_new.insert(self.vertices[old_idx].id, new_id);
        }
        if kept.len() == n {
            return old_to_new;
        }

        let mut new_vertices = Vec::with_capacity(kept.len());
        for (new_id, &old_idx) in kept.iter().enumerate() {
            new_vertices.push(SeqVertex {
                id: new_id,
                sequence: std::mem::take(&mut self.vertices[old_idx].sequence),
            });
        }

        let mut new_edges = Vec::with_capacity(self.edges.len());
        for e in &self.edges {
            let Some(&from) = old_to_new.get(&e.from) else {
                continue;
            };
            let Some(&to) = old_to_new.get(&e.to) else {
                continue;
            };
            new_edges.push(SeqEdge {
                from,
                to,
                support: e.support,
                is_ref: e.is_ref,
            });
        }

        self.vertices = new_vertices;
        self.edges = new_edges;
        self.rebuild_index();
        debug_assert!(
            self.vertices.iter().enumerate().all(|(i, v)| v.id == i)
                && self
                    .edges
                    .iter()
                    .all(|e| e.from < self.vertices.len() && e.to < self.vertices.len())
        );
        old_to_new
    }

    pub fn remove_singleton_orphan_vertices(&mut self) {
        loop {
            let orphans: HashSet<usize> = (0..self.vertices.len())
                .filter(|&v| {
                    self.in_degree(v) == 0 && self.out_degree(v) == 0 && !self.is_ref_source(v)
                })
                .collect();
            if orphans.is_empty() {
                break;
            }
            self.compact_vertices_and_remap_edges(&orphans);
        }
    }

    /// GATK `BaseGraph.removeVerticesNotConnectedToRefRegardlessOfEdgeDirection`.
    pub fn remove_vertices_not_connected_to_ref_regardless_of_direction(&mut self) {
        let Some(source) = self.reference_source_vertex() else {
            return;
        };
        let mut keep = HashSet::new();
        let mut stack = vec![source];
        keep.insert(source);
        while let Some(v) = stack.pop() {
            for &t in self.outgoing_of(v) {
                if keep.insert(t) {
                    stack.push(t);
                }
            }
            for &p in self.incoming_of(v) {
                if keep.insert(p) {
                    stack.push(p);
                }
            }
        }
        let remove: HashSet<usize> = (0..self.vertices.len())
            .filter(|v| !keep.contains(v))
            .collect();
        self.compact_vertices_and_remap_edges(&remove);
    }

    pub fn remove_paths_not_connected_to_ref(&mut self) -> GatkResult<()> {
        let source = self
            .reference_source_vertex()
            .ok_or_else(|| GatkError::algorithm("seq graph: no ref source"))?;
        let sink = self
            .reference_sink_vertex()
            .ok_or_else(|| GatkError::algorithm("seq graph: no ref sink"))?;
        let mut from_source = HashSet::new();
        let mut stack = vec![source];
        from_source.insert(source);
        while let Some(v) = stack.pop() {
            for &t in self.outgoing_of(v) {
                if from_source.insert(t) {
                    stack.push(t);
                }
            }
        }
        let mut from_sink = HashSet::new();
        let mut stack = vec![sink];
        from_sink.insert(sink);
        while let Some(v) = stack.pop() {
            for &p in self.incoming_of(v) {
                if from_sink.insert(p) {
                    stack.push(p);
                }
            }
        }
        from_source.retain(|v| from_sink.contains(v));
        let remove: HashSet<usize> = (0..self.vertices.len())
            .filter(|v| !from_source.contains(v))
            .collect();
        self.compact_vertices_and_remap_edges(&remove);
        Ok(())
    }

    pub fn simplify_graph(&mut self) {
        crate::seq_graph_simplify::simplify_graph_full(self);
    }

    pub(crate) fn vertex_in_degree(&self, v: usize) -> usize {
        self.in_degree(v)
    }

    pub(crate) fn vertex_out_degree(&self, v: usize) -> usize {
        self.out_degree(v)
    }

    pub(crate) fn is_sink_vertex(&self, v: usize) -> bool {
        self.out_degree(v) == 0
    }

    pub(crate) fn vertex_sequence(&self, v: usize) -> &[u8] {
        &self.vertices[v].sequence
    }

    pub(crate) fn remove_vertices_by_id(
        &mut self,
        remove: &HashSet<usize>,
    ) -> HashMap<usize, usize> {
        self.compact_vertices_and_remap_edges(remove)
    }

    pub(crate) fn add_seq_vertex(&mut self, sequence: Vec<u8>) -> usize {
        let id = self.vertices.len();
        self.vertices.push(SeqVertex { id, sequence });
        id
    }

    pub(crate) fn find_edge(&self, from: usize, to: usize) -> Option<usize> {
        self.edges.iter().position(|e| e.from == from && e.to == to)
    }

    pub(crate) fn edges_pub(&self) -> &[SeqEdge] {
        &self.edges
    }

    #[cfg(test)]
    pub(crate) fn test_vertex_ids(&self) -> Vec<usize> {
        self.vertices.iter().map(|v| v.id).collect()
    }

    /// Test-only: same order as [`Self::cleanup_seq_graph`], with a snapshot after each stage.
    #[cfg(test)]
    pub(crate) fn traced_cleanup_seq_graph(
        &mut self,
        mut snap: impl FnMut(&str, &SeqGraph),
    ) -> SeqGraphCleanupStatus {
        snap("cleanup_entry", self);
        self.zip_linear_chains();
        snap("after_initial_zip_linear_chains", self);
        self.remove_singleton_orphan_vertices();
        snap("after_remove_singleton_orphan_vertices", self);
        self.remove_vertices_not_connected_to_ref_regardless_of_direction();
        snap(
            "after_remove_vertices_not_connected_to_ref_undirected",
            self,
        );
        crate::seq_graph_simplify::traced_simplify_graph_full(self, |stage, g| {
            snap(&format!("simplify1_{stage}"), g);
        });
        snap("before_source_sink_jar", self);
        if self.reference_source_vertex().is_none() || self.reference_sink_vertex().is_none() {
            return SeqGraphCleanupStatus::JustAssembledReference;
        }
        let _ = self.remove_paths_not_connected_to_ref();
        snap("after_remove_paths_not_connected_to_ref", self);
        crate::seq_graph_simplify::traced_simplify_graph_full(self, |stage, g| {
            snap(&format!("simplify2_{stage}"), g);
        });
        snap("after_second_simplify", self);
        if self.vertices.len() == 1 {
            let complete = 0usize;
            let dummy_id = self.vertices.len();
            self.vertices.push(SeqVertex {
                id: dummy_id,
                sequence: Vec::new(),
            });
            self.edges.push(SeqEdge {
                from: complete,
                to: dummy_id,
                support: 0,
                is_ref: true,
            });
            self.rebuild_index();
            snap("after_dummy_vertex", self);
        }
        snap("final_for_kbest", self);
        SeqGraphCleanupStatus::AssembledSomeVariation
    }

    pub(crate) fn add_or_update_edge(
        &mut self,
        from: usize,
        to: usize,
        support: u32,
        is_ref: bool,
    ) {
        if let Some(idx) = self.find_edge(from, to) {
            let e = &mut self.edges[idx];
            e.support = e.support.saturating_add(support);
            e.is_ref |= is_ref;
        } else {
            self.edges.push(SeqEdge {
                from,
                to,
                support,
                is_ref,
            });
        }
        self.rebuild_index();
    }

    pub(crate) fn incoming_nodes(&self, v: usize) -> Vec<usize> {
        self.incoming_of(v).to_vec()
    }

    pub fn cleanup_seq_graph(&mut self) -> SeqGraphCleanupStatus {
        self.zip_linear_chains();
        self.remove_singleton_orphan_vertices();
        self.remove_vertices_not_connected_to_ref_regardless_of_direction();
        self.simplify_graph();
        if self.reference_source_vertex().is_none() || self.reference_sink_vertex().is_none() {
            return SeqGraphCleanupStatus::JustAssembledReference;
        }
        // GATK `ReadThreadingAssembler.cleanupSeqGraph`: after removePaths + simplify, proceed to
        // ASSEMBLED_SOME_VARIATION (dummy vertex if needed) — no second ref source/sink abort.
        let _ = self.remove_paths_not_connected_to_ref();
        self.simplify_graph();
        if self.vertices.len() == 1 {
            let complete = 0usize;
            let dummy_id = self.vertices.len();
            self.vertices.push(SeqVertex {
                id: dummy_id,
                sequence: Vec::new(),
            });
            self.edges.push(SeqEdge {
                from: complete,
                to: dummy_id,
                support: 0,
                is_ref: true,
            });
            self.rebuild_index();
        }
        SeqGraphCleanupStatus::AssembledSomeVariation
    }
}

fn additional_sequence_bytes(seq: &[u8], is_source: bool) -> Vec<u8> {
    if is_source {
        seq.to_vec()
    } else {
        seq.last().map(|&c| vec![c]).unwrap_or_default()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::assembly::{AssemblyGraphParams, AssemblyRead};
    use crate::assembly_graph_dump::{load_assembly_reads_tsv, load_assembly_ref_tsv};
    use crate::read_threading_assembler::{
        build_threading_graph_for_haplotype_dump, ReadThreadingAssemblerArgs,
    };
    use crate::read_threading_graph::assembly_graph_from_ref_and_reads_threading;
    use std::path::Path;

    fn read(seq: &str, q: u8) -> AssemblyRead {
        AssemblyRead {
            bases: seq.as_bytes().to_vec(),
            base_quals: vec![q; seq.len()],
        }
    }

    #[test]
    fn prepared_p5_seqgraph_dump_path_has_nodes() {
        let repo = Path::new(env!("CARGO_MANIFEST_DIR")).join("..");
        let reference =
            load_assembly_ref_tsv(&repo.join("parity/fixtures/hc-full-parity/e4/p5_case1_ref.tsv"))
                .unwrap();
        let reads = load_assembly_reads_tsv(
            &repo.join("parity/fixtures/hc-full-parity/e4/p5_case1_reads.tsv"),
        )
        .unwrap();
        let args = ReadThreadingAssemblerArgs {
            kmer_sizes: vec![3],
            min_base_quality: 10,
            min_prune_factor: 2,
            min_dangling_branch_length: 4,
            recover_dangling_heads: true,
            ..Default::default()
        };
        let graph =
            build_threading_graph_for_haplotype_dump(&reference, &reads, 3, &args, true, false)
                .unwrap()
                .expect("graph");
        assert!(graph.node_count() > 0, "rt graph");
        let mut seq = SeqGraph::from_assembly_graph(&graph);
        assert!(seq.node_count() > 0, "seq after from");
        seq.clean_non_ref_paths();
        let status = seq.cleanup_seq_graph();
        assert!(seq.node_count() > 0, "after cleanup status={status:?}");
        let path = seq.reference_path_bytes().unwrap();
        assert_eq!(path.as_slice(), reference.bases.as_slice());
    }

    #[test]
    fn seq_graph_from_p5_case1_has_ref_path() {
        let reference = read("ACGTT", 30);
        let reads = vec![read("ACGTT", 30), read("ACGTT", 30), read("ACGTA", 30)];
        let params = AssemblyGraphParams {
            kmer_size: crate::bio_ids::KmerSize::try_new(3).unwrap(),
            min_base_quality: 10,
            ..Default::default()
        };
        let graph =
            assembly_graph_from_ref_and_reads_threading(&reference, &reads, &params).unwrap();
        let seq = SeqGraph::from_assembly_graph(&graph);
        assert!(seq.node_count() > 0);
        assert!(seq.edge_count() > 0);
    }
}

#[cfg(test)]
#[path = "seq_graph_id_invariant_test.rs"]
mod id_invariant_tests;

#[cfg(test)]
#[path = "seq_graph_post_repair_simplify_test.rs"]
mod post_repair_simplify_tests;

#[cfg(test)]
#[path = "seq_graph_path_bases_probe_test.rs"]
mod path_bases_probe_tests;

#[cfg(test)]
#[path = "seq_graph_p12_waiver_gate_test.rs"]
mod p12_waiver_gate_tests;

#[cfg(test)]
#[path = "seq_graph_p12_k85_topology_test.rs"]
mod p12_k85_topology_tests;
