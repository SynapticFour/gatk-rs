//! GATK `JunctionTreeLinkedDeBruijnGraph` parity: graph build + `generateJunctionTrees`.

use crate::assembly::{AssemblyGraph, AssemblyRead, KmerNode};
use gatk_common::{GatkError, GatkResult};
use std::collections::{BTreeMap, BTreeSet, BinaryHeap, HashMap, HashSet};

const INCREASE_COUNTS_BACKWARDS: bool = true;
const SYMBOLIC_END_KMER: &str = "_";

#[derive(Debug, Clone)]
struct SequenceForKmers {
    bases: Vec<u8>,
    start: usize,
    stop: usize,
    count: u32,
    is_ref: bool,
}

#[derive(Debug)]
struct ThreadingEdge {
    total: u32,
    current_sample: u32,
    flushed_samples: BinaryHeap<std::cmp::Reverse<u32>>,
    num_pruning_samples: usize,
    is_ref: bool,
}

impl ThreadingEdge {
    fn new(initial: u32, num_pruning_samples: usize, is_ref: bool) -> Self {
        let mut flushed_samples = BinaryHeap::new();
        flushed_samples.push(std::cmp::Reverse(initial));
        Self {
            total: initial,
            current_sample: initial,
            flushed_samples,
            num_pruning_samples,
            is_ref,
        }
    }

    fn inc(&mut self, delta: u32) {
        self.total = self.total.saturating_add(delta);
        self.current_sample = self.current_sample.saturating_add(delta);
    }

    fn total(&self) -> u32 {
        self.total
    }

    fn flush_sample(&mut self) {
        self.flushed_samples
            .push(std::cmp::Reverse(self.current_sample));
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

/// Junction-tree node (GATK `ThreadingNode`).
#[derive(Debug)]
pub struct ThreadingNode {
    /// Ordered by edge `(from, to)` so branch enumeration / evidence ties are deterministic.
    children: BTreeMap<(usize, usize), ThreadingNode>,
    evidence_count: i32,
    prev_edge: Option<(usize, usize)>,
}

impl ThreadingNode {
    fn new(prev_edge: Option<(usize, usize)>) -> Self {
        Self {
            children: BTreeMap::new(),
            evidence_count: 0,
            prev_edge,
        }
    }

    fn increment_count(&mut self) {
        self.evidence_count += 1;
    }

    fn add_edge(&mut self, edge: (usize, usize)) -> &mut ThreadingNode {
        self.increment_count();
        let child = self
            .children
            .entry(edge)
            .or_insert_with(|| ThreadingNode::new(Some(edge)));
        child.increment_count();
        self.children.get_mut(&edge).expect("just inserted")
    }

    pub fn evidence_count(&self) -> i32 {
        self.evidence_count
    }

    pub fn children(&self) -> &BTreeMap<(usize, usize), ThreadingNode> {
        &self.children
    }

    pub fn has_no_evidence(&self) -> bool {
        self.children.is_empty()
    }

    pub fn is_symbolic_end(&self, symbolic_edge: (usize, usize)) -> bool {
        self.prev_edge == Some(symbolic_edge)
    }
}

/// Junction tree at a fork (GATK `ThreadingTree`).
#[derive(Debug)]
pub struct ThreadingTree {
    pub root: ThreadingNode,
    #[allow(dead_code)]
    tree_vertex: usize,
}

impl ThreadingTree {
    fn new(tree_vertex: usize) -> Self {
        Self {
            root: ThreadingNode::new(None),
            tree_vertex,
        }
    }

    fn get_and_increment_root(&mut self) -> &mut ThreadingNode {
        self.root.increment_count();
        &mut self.root
    }

    pub fn is_empty_tree(&self) -> bool {
        !self.root.children.is_empty()
    }
}

/// Assembly graph plus junction trees (GATK `JunctionTreeLinkedDeBruijnGraph`).
pub struct JunctionTreeLinkedGraph {
    pub graph: AssemblyGraph,
    junction_trees: HashMap<usize, ThreadingTree>,
    pub symbolic_end_vertex: usize,
    pub symbolic_end_edge: (usize, usize),
    reference_source: Option<usize>,
    reference_sink: Option<usize>,
}

impl JunctionTreeLinkedGraph {
    pub fn junction_tree_for_node(&self, vertex: usize) -> Option<&ThreadingTree> {
        self.junction_trees.get(&vertex)
    }

    /// GATK `JunctionTreeLinkedDeBruijnGraph.getReferenceSourceVertex` (from `referencePath`).
    pub fn reference_source(&self) -> Option<usize> {
        self.reference_source
    }

    /// GATK `JunctionTreeLinkedDeBruijnGraph.getReferenceSinkVertex` (from `referencePath`).
    pub fn reference_sink(&self) -> Option<usize> {
        self.reference_sink
    }
}

struct JunctionTreeGraphBuilder {
    kmer_size: usize,
    #[allow(dead_code)]
    min_base_quality: u8,
    pending: Vec<SequenceForKmers>,
    kmer_to_vertex: BTreeMap<Vec<u8>, usize>,
    nodes: Vec<Vec<u8>>,
    edges: HashMap<(usize, usize), ThreadingEdge>,
    edge_is_ref: HashSet<(usize, usize)>,
    ref_nodes: HashSet<usize>,
    ref_source_kmer: Option<Vec<u8>>,
    reference_path: Vec<usize>,
    /// Neighbor sets are ordered so first-match threading is deterministic.
    outgoing: HashMap<usize, BTreeSet<usize>>,
    incoming: HashMap<usize, BTreeSet<usize>>,
    junction_trees: HashMap<usize, ThreadingTree>,
    symbolic_end_vertex: Option<usize>,
    symbolic_end_edge: Option<(usize, usize)>,
    built: bool,
}

#[derive(Clone)]
struct ActiveTreePath {
    tree_vertex: usize,
    edges: Vec<(usize, usize)>,
}

#[derive(Clone)]
struct JtGraphSnapshot {
    kmer_size: usize,
    nodes: Vec<Vec<u8>>,
    kmer_to_vertex: BTreeMap<Vec<u8>, usize>,
    outgoing: HashMap<usize, BTreeSet<usize>>,
    incoming: HashMap<usize, BTreeSet<usize>>,
    reference_sink: Option<usize>,
}

impl JtGraphSnapshot {
    fn from_builder(b: &JunctionTreeGraphBuilder) -> Self {
        Self {
            kmer_size: b.kmer_size,
            nodes: b.nodes.clone(),
            kmer_to_vertex: b.kmer_to_vertex.clone(),
            outgoing: b.outgoing.clone(),
            incoming: b.incoming.clone(),
            reference_sink: b.reference_sink(),
        }
    }

    fn suffix_of_kmer(&self, id: usize) -> u8 {
        JunctionTreeGraphBuilder::suffix_of_kmer(&self.nodes[id])
    }

    fn get_kmer_vertex(&self, bases: &[u8], start: usize) -> Option<usize> {
        let kmer = JunctionTreeGraphBuilder::kmer_at(bases, start, self.kmer_size);
        self.kmer_to_vertex.get(&kmer).copied()
    }

    fn extend_junction_threading_by_one(
        &self,
        prev: usize,
        sequence: &[u8],
        kmer_start: usize,
        helper: Option<&mut JunctionThreadingState<'_>>,
        alter_trees: bool,
    ) -> Option<usize> {
        let k = self.kmer_size;
        let outs: Vec<usize> = self
            .outgoing
            .get(&prev)
            .map(|s| s.iter().copied().collect())
            .unwrap_or_default();
        let next_pos = kmer_start + k - 1;
        let next_base = sequence[next_pos];
        for &to in &outs {
            if self.suffix_of_kmer(to) == next_base {
                if alter_trees {
                    if let Some(h) = helper {
                        h.add_tree_if_necessary(prev, self);
                        if outs.len() != 1 {
                            h.add_edge_to_junction_tree_nodes((prev, to));
                        }
                    }
                }
                return Some(to);
            }
        }
        None
    }
}

struct JunctionThreadingState<'a> {
    active_paths: Vec<ActiveTreePath>,
    trees: &'a mut HashMap<usize, ThreadingTree>,
}

impl JunctionThreadingState<'_> {
    fn clear(&mut self) {
        self.active_paths.clear();
    }

    fn vertex_warrants_junction_tree(&self, vertex: usize, snap: &JtGraphSnapshot) -> bool {
        let out_degree = snap.outgoing.get(&vertex).map(|s| s.len()).unwrap_or(0);
        if out_degree > 1 {
            return true;
        }
        let outs: Vec<usize> = snap
            .outgoing
            .get(&vertex)
            .map(|s| s.iter().copied().collect())
            .unwrap_or_default();
        outs.iter()
            .any(|&to| snap.incoming.get(&to).map(|s| s.len()).unwrap_or(0) > 1)
    }

    fn add_tree_if_necessary(&mut self, prev_vertex: usize, snap: &JtGraphSnapshot) {
        if !self.vertex_warrants_junction_tree(prev_vertex, snap) {
            return;
        }
        self.trees
            .entry(prev_vertex)
            .or_insert_with(|| ThreadingTree::new(prev_vertex))
            .get_and_increment_root();
        self.active_paths.push(ActiveTreePath {
            tree_vertex: prev_vertex,
            edges: Vec::new(),
        });
    }

    fn add_edge_to_junction_tree_nodes(&mut self, edge: (usize, usize)) {
        let mut new_paths = Vec::new();
        for path in self.active_paths.drain(..) {
            if let Some(tree) = self.trees.get_mut(&path.tree_vertex) {
                let mut node = &mut tree.root;
                for &e in &path.edges {
                    node = node.children.get_mut(&e).expect("junction tree path");
                }
                node.add_edge(edge);
                let mut edges = path.edges;
                edges.push(edge);
                new_paths.push(ActiveTreePath {
                    tree_vertex: path.tree_vertex,
                    edges,
                });
            }
        }
        self.active_paths = new_paths;
    }
}

impl JunctionTreeGraphBuilder {
    fn kmer_at(bases: &[u8], start: usize, k: usize) -> Vec<u8> {
        bases[start..start + k].to_vec()
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
        for end in 0..=read.bases.len() {
            let unusable = end == read.bases.len()
                || read.base_quals[end] < min_qual
                || !is_base_usable(read.bases[end]);
            if unusable {
                if let Some(start) = last_good {
                    if end - start >= kmer_size {
                        out.push(SequenceForKmers {
                            bases: read.bases.to_vec(),
                            start,
                            stop: end,
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

    fn find_start(&self, seq: &SequenceForKmers) -> Option<usize> {
        if seq.is_ref {
            return Some(seq.start);
        }
        Some(seq.start)
    }

    fn find_start_for_junction_threading(&self, seq: &SequenceForKmers) -> Option<usize> {
        let k = self.kmer_size;
        let last = seq.stop.saturating_sub(k);
        for i in seq.start..last {
            let key = Self::kmer_at(&seq.bases, i, k);
            if self.kmer_to_vertex.contains_key(&key) {
                return Some(i);
            }
        }
        None
    }

    fn track_kmer(&mut self, kmer: Vec<u8>, id: usize) {
        self.kmer_to_vertex.entry(kmer).or_insert(id);
    }

    fn create_vertex(&mut self, kmer: Vec<u8>) -> usize {
        if let Some(&id) = self.kmer_to_vertex.get(&kmer) {
            return id;
        }
        let id = self.nodes.len();
        // CLONE: needed because owned element into collection.
        self.nodes.push(kmer.clone());
        self.track_kmer(kmer, id);
        id
    }

    fn get_or_create_kmer_vertex(&mut self, bases: &[u8], start: usize) -> usize {
        let kmer = Self::kmer_at(bases, start, self.kmer_size);
        if let Some(&id) = self.kmer_to_vertex.get(&kmer) {
            return id;
        }
        self.create_vertex(kmer)
    }

    fn reference_sink(&self) -> Option<usize> {
        self.reference_path.last().copied()
    }

    fn inc_edge(&mut self, from: usize, to: usize, delta: u32, is_ref: bool) {
        if let Some(e) = self.edges.get_mut(&(from, to)) {
            e.inc(delta);
            if is_ref {
                e.is_ref = true;
                self.edge_is_ref.insert((from, to));
            }
        } else {
            self.edges
                .insert((from, to), ThreadingEdge::new(delta, 1, is_ref));
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
            if suffix == target_base && in_degree == 1 {
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
        let next = self
            .kmer_to_vertex
            .get(&kmer)
            .copied()
            .unwrap_or_else(|| self.create_vertex(kmer));
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
        if seq.is_ref {
            self.reference_path.clear();
            self.reference_path.push(start_vertex);
        }
        let mut vertex = start_vertex;
        for i in (start_pos + 1)..=(seq.stop.saturating_sub(k)) {
            vertex = self.extend_chain_by_one(vertex, &seq.bases, i, seq.count, seq.is_ref);
            if seq.is_ref {
                self.reference_path.push(vertex);
            }
        }
    }

    fn build(&mut self) {
        if self.built {
            return;
        }
        // Lifetime: take pending while threading mutates graph fields; restore afterward
        // because generate_junction_trees still needs the sequences.
        let pending = std::mem::take(&mut self.pending);
        for seq in &pending {
            self.thread_sequence(seq);
        }
        for e in self.edges.values_mut() {
            e.flush_sample();
        }
        self.pending = pending;
        self.built = true;
    }

    fn generate_junction_trees(&mut self) {
        assert!(self.built, "build graph before generateJunctionTrees");
        let snap = JtGraphSnapshot::from_builder(self);
        let symbolic_end_vertex = self.create_vertex(SYMBOLIC_END_KMER.as_bytes().to_vec());
        self.symbolic_end_vertex = Some(symbolic_end_vertex);
        let sink = self
            .reference_sink()
            .expect("reference sink for symbolic edge");
        self.inc_edge(sink, symbolic_end_vertex, 0, true);
        let symbolic_end_edge = (sink, symbolic_end_vertex);
        self.symbolic_end_edge = Some(symbolic_end_edge);

        // Lifetime: final consumer of pending sequences; move instead of clone.
        let pending = std::mem::take(&mut self.pending);
        for seq in pending {
            if seq.is_ref {
                continue;
            }
            let Some(start_pos) = self.find_start_for_junction_threading(&seq) else {
                continue;
            };
            let k = snap.kmer_size;
            let starting_vertex = snap
                .get_kmer_vertex(&seq.bases, start_pos)
                .expect("junction threading start in map");
            let mut helper = JunctionThreadingState {
                active_paths: Vec::new(),
                trees: &mut self.junction_trees,
            };
            let mut last_vertex = starting_vertex;
            let mut has_to_rediscover = false;
            for i in (start_pos + 1)..=(seq.stop.saturating_sub(k)) {
                let vertex = if !has_to_rediscover {
                    snap.extend_junction_threading_by_one(
                        last_vertex,
                        &seq.bases,
                        i,
                        Some(&mut helper),
                        true,
                    )
                } else {
                    snap.get_kmer_vertex(&seq.bases, i)
                };

                let vertex = if vertex.is_none() && !has_to_rediscover {
                    let outs: Vec<usize> = snap
                        .outgoing
                        .get(&last_vertex)
                        .map(|s| s.iter().copied().collect())
                        .unwrap_or_default();
                    if outs.len() == 1 {
                        let first = outs[0];
                        let mut tentative = Some(first);
                        let mut j = i + 1;
                        while j <= seq.stop.saturating_sub(k) && j <= i + k && tentative.is_some() {
                            tentative = snap.extend_junction_threading_by_one(
                                tentative.unwrap(),
                                &seq.bases,
                                j,
                                None,
                                false,
                            );
                            j += 1;
                        }
                        if tentative.is_some() {
                            Some(first)
                        } else {
                            None
                        }
                    } else {
                        None
                    }
                } else {
                    vertex
                };

                if let Some(v) = vertex {
                    last_vertex = v;
                    has_to_rediscover = false;
                } else {
                    helper.clear();
                    has_to_rediscover = true;
                }
            }

            if Some(last_vertex) == snap.reference_sink {
                helper.add_edge_to_junction_tree_nodes(symbolic_end_edge);
            }
        }
    }

    fn finish(mut self) -> JunctionTreeLinkedGraph {
        self.build();
        self.generate_junction_trees();
        let symbolic_end_vertex = self.symbolic_end_vertex.expect("symbolic end");
        let symbolic_end_edge = self.symbolic_end_edge.expect("symbolic edge");
        let junction_trees = std::mem::take(&mut self.junction_trees);
        let nodes: Vec<KmerNode> = self
            .nodes
            .iter()
            .enumerate()
            .map(|(id, kmer)| KmerNode {
                id,
                kmer: std::sync::Arc::from(kmer.as_slice()),
                support: 1,
            })
            .collect();
        let kmer_to_id: HashMap<_, _> = self
            .kmer_to_vertex
            .iter()
            .map(|(k, &id)| (std::sync::Arc::from(k.as_slice()), id))
            .collect();
        let ref_source_kmer = self
            .ref_source_kmer
            .as_ref()
            .map(|k| std::sync::Arc::from(k.as_slice()));
        let mut edges = HashMap::with_capacity(self.edges.len());
        let mut pruning_edges = HashMap::with_capacity(self.edges.len());
        for (&(from, to), e) in &self.edges {
            edges.insert((from, to), e.total());
            pruning_edges.insert((from, to), e.pruning_multiplicity());
        }
        let outgoing: HashMap<_, _> = self
            .outgoing
            .iter()
            .map(|(&k, v)| (k, v.iter().copied().collect()))
            .collect();
        let incoming: HashMap<_, _> = self
            .incoming
            .iter()
            .map(|(&k, v)| (k, v.iter().copied().collect()))
            .collect();
        let mut graph = AssemblyGraph::from_threading_build(
            self.kmer_size,
            nodes,
            kmer_to_id,
            edges,
            pruning_edges,
            outgoing,
            incoming,
            self.edge_is_ref.clone(),
            self.ref_nodes.clone(),
            ref_source_kmer,
        );
        graph.cleanup_isolated_nodes();
        let ref_kmer_path: Vec<Vec<u8>> = self
            .reference_path
            .iter()
            .map(|&id| self.nodes[id].clone())
            .collect();
        let reference_source = ref_kmer_path
            .first()
            .and_then(|k| graph.vertex_id_for_kmer(k));
        let reference_sink = ref_kmer_path
            .last()
            .and_then(|k| graph.vertex_id_for_kmer(k));
        JunctionTreeLinkedGraph {
            graph,
            junction_trees,
            symbolic_end_vertex,
            symbolic_end_edge,
            reference_source,
            reference_sink,
        }
    }
}

fn is_base_usable(base: u8) -> bool {
    matches!(base, b'A' | b'C' | b'G' | b'T' | b'N')
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::assembly_graph_dump::{load_assembly_reads_tsv, load_assembly_ref_tsv};
    use std::path::Path;

    #[test]
    fn p5_jt_snp_het_graph_has_ref_endpoints() {
        let repo = Path::new(env!("CARGO_MANIFEST_DIR")).join("..");
        let reference = load_assembly_ref_tsv(
            &repo.join("parity/fixtures/hc-full-parity/e7-junction/p5_jt_snp_het_ref.tsv"),
        )
        .unwrap();
        let reads = load_assembly_reads_tsv(
            &repo.join("parity/fixtures/hc-full-parity/e7-junction/p5_jt_snp_het_reads.tsv"),
        )
        .unwrap();
        let jt = build_junction_tree_graph_from_ref_and_reads(&reference, &reads, 5, 10).unwrap();
        assert!(jt.graph.node_count() > 0, "graph has nodes");
        assert!(jt.graph.edge_count() > 0, "graph has edges");
        assert!(jt.reference_source().is_some(), "ref source");
        assert!(jt.reference_sink().is_some(), "ref sink");
        let src = jt.reference_source().unwrap();
        let snk = jt.reference_sink().unwrap();
        let mut stack = vec![src];
        let mut seen = std::collections::HashSet::new();
        seen.insert(src);
        let mut reaches_sink = false;
        while let Some(v) = stack.pop() {
            if v == snk {
                reaches_sink = true;
                break;
            }
            for to in jt.graph.outgoing_nodes(v) {
                if seen.insert(to) {
                    stack.push(to);
                }
            }
        }
        assert!(reaches_sink, "ref source must reach ref sink in graph");
        if let Some(t) = jt.junction_tree_for_node(11) {
            for (e, c) in t.root.children() {
                eprintln!("  jt@11 child {e:?} ev={}", c.evidence_count());
            }
        }
        for v in 0..jt.graph.node_count() {
            let outs = jt.graph.outgoing_nodes(v);
            if outs.len() > 1 {
                eprintln!(
                    "fork@{v} kmer={}",
                    String::from_utf8_lossy(jt.graph.kmer_at(v))
                );
                for to in outs {
                    eprintln!(
                        "  ->{to} kmer={} sup={} ref={} in_deg={}",
                        String::from_utf8_lossy(jt.graph.kmer_at(to)),
                        jt.graph.edge_support(v, to).unwrap_or(0),
                        jt.graph.edge_is_ref(v, to),
                        jt.graph.incoming_count(to)
                    );
                }
            }
        }
        for target in [24usize, 29, 12] {
            eprintln!("in@{target} preds={:?}", jt.graph.incoming_nodes(target));
        }
        let paths = crate::junction_kbest::find_junction_best_haplotypes(&jt, 5, 1, true).unwrap();
        for p in &paths {
            eprintln!(
                "  path {} score={}",
                String::from_utf8_lossy(&p.bases(&jt.graph)),
                p.score
            );
        }
        if paths.is_empty() {
            eprintln!(
                "nodes={} edges={} jt_tree_vertices={} src={src} snk={snk} out_sink={}",
                jt.graph.node_count(),
                jt.graph.edge_count(),
                (0..jt.graph.node_count())
                    .filter(|v| jt.junction_tree_for_node(*v).is_some())
                    .count(),
                jt.graph.outgoing_nodes(snk).len()
            );
            for v in 0..jt.graph.node_count() {
                if let Some(t) = jt.junction_tree_for_node(v) {
                    let edges: Vec<_> = t.root.children().keys().collect();
                    eprintln!("  jt@{v} root_children={edges:?}");
                }
            }
        }
        assert!(!paths.is_empty(), "expected haplotypes");
    }
}

/// Build junction-tree-linked graph from reference + reads (GATK parity dump path).
pub fn build_junction_tree_graph_from_ref_and_reads(
    reference: &AssemblyRead,
    reads: &[AssemblyRead],
    kmer_size: usize,
    min_base_quality: u8,
) -> GatkResult<JunctionTreeLinkedGraph> {
    let _ = crate::bio_ids::KmerSize::try_from_usize(kmer_size)?;
    let mut builder = JunctionTreeGraphBuilder {
        kmer_size,
        min_base_quality,
        pending: Vec::new(),
        kmer_to_vertex: BTreeMap::new(),
        nodes: Vec::new(),
        edges: HashMap::new(),
        edge_is_ref: HashSet::new(),
        ref_nodes: HashSet::new(),
        ref_source_kmer: None,
        reference_path: Vec::new(),
        outgoing: HashMap::new(),
        incoming: HashMap::new(),
        junction_trees: HashMap::new(),
        symbolic_end_vertex: None,
        symbolic_end_edge: None,
        built: false,
    };

    for seq in JunctionTreeGraphBuilder::sequences_from_read(reference, kmer_size, min_base_quality)
    {
        let mut s = seq;
        s.is_ref = true;
        builder.pending.push(s);
    }
    for read in reads {
        if read.bases.len() != read.base_quals.len() {
            return Err(GatkError::argument(
                "read bases length must match base quality length",
            ));
        }
        for seq in JunctionTreeGraphBuilder::sequences_from_read(read, kmer_size, min_base_quality)
        {
            builder.pending.push(seq);
        }
    }

    Ok(builder.finish())
}
