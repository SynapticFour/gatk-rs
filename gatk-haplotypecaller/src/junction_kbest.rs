//! GATK `JunctionTreeKBestHaplotypeFinder` / `JTBestHaplotype` parity.

use crate::assembly::AssemblyGraph;
use crate::junction_tree_graph::{JunctionTreeLinkedGraph, ThreadingTree};
use crate::kbest_haplotype::log_penalty;
use gatk_common::{GatkError, GatkResult};
use std::cmp::Ordering;
use std::collections::{BinaryHeap, HashMap, HashSet};

pub const DEFAULT_OUTGOING_JT_EVIDENCE_THRESHOLD_TO_BELIEVE: i32 = 3;
const DEFAULT_MAX_ACCEPTABLE_DECISION_EDGES_WITHOUT_JT_GUIDANCE: i32 = 5;
const DEFAULT_MAX_ACCEPTABLE_REPETITIONS_OF_A_KMER_IN_A_PATH: usize = 1;
const DEFAULT_MAX_PATHS_TO_CONSIDER_WITHOUT_RESULT: usize = 1000;
const DEFAULT_MAX_PATHS_TO_EVER_CONSIDER: usize = 10000;

type Edge = (usize, usize);

/// One scored path through the junction-tree k-best haplotype finder.
/// # Invariants
/// `edges` is an ordered `(from, to)` walk from `start` through the [`AssemblyGraph`].
/// `score` ranks paths (higher is better in the finder heap).
/// # Ownership
/// Owns edge list; borrows graph in [`JunctionKBestPath::bases`].
/// # Mutation
/// Immutable path snapshot after search completes.
/// # Biological assumptions
/// Path sequence is a candidate haplotype over the active assembly graph.
/// # Java equivalence
/// GATK `JunctionTreeKBestHaplotypeFinder` / `JTBestHaplotype` path representation.
#[derive(Debug, Clone)]
pub struct JunctionKBestPath {
    pub start: usize,
    pub edges: Vec<Edge>,
    pub score: f64,
    pub is_reference: bool,
}

impl JunctionKBestPath {
    pub fn bases(&self, graph: &AssemblyGraph) -> String {
        graph.path_bases(self.start, &self.edges)
    }
}

#[derive(Debug, Clone)]
struct JunctionTreeSet {
    visited_trees: HashSet<usize>,
    active_nodes: Vec<ActiveNode>,
}

#[derive(Debug, Clone)]
struct ActiveNode {
    tree_vertex: usize,
    path_edges: Vec<Edge>,
}

impl JunctionTreeSet {
    fn new() -> Self {
        Self {
            visited_trees: HashSet::new(),
            active_nodes: Vec::new(),
        }
    }

    fn clone_from(other: &Self) -> Self {
        Self {
            visited_trees: other.visited_trees.clone(),
            active_nodes: other.active_nodes.clone(),
        }
    }

    fn add_junction_tree(&mut self, tree_vertex: usize, tree: &ThreadingTree) -> bool {
        if self.visited_trees.contains(&tree_vertex) || tree.root.has_no_evidence() {
            return false;
        }
        self.visited_trees.insert(tree_vertex);
        self.active_nodes.push(ActiveNode {
            tree_vertex,
            path_edges: Vec::new(),
        });
        true
    }

    fn traverse_edge_for_all_trees(&mut self, jt: &JunctionTreeLinkedGraph, edge: Edge) {
        self.active_nodes = self
            .active_nodes
            .iter()
            .filter_map(|active| {
                let tree = jt.junction_tree_for_node(active.tree_vertex)?;
                let mut node = &tree.root;
                for &e in &active.path_edges {
                    node = node.children().get(&e)?;
                }
                if !node.children().contains_key(&edge) {
                    return None;
                }
                let mut path_edges = Vec::with_capacity(active.path_edges.len() + 1);
                path_edges.extend_from_slice(&active.path_edges);
                path_edges.push(edge);
                let child = node.children().get(&edge)?;
                if child.has_no_evidence() {
                    return None;
                }
                Some(ActiveNode {
                    tree_vertex: active.tree_vertex,
                    path_edges,
                })
            })
            .collect();
    }

    fn prune_empty_nodes(&mut self, jt: &JunctionTreeLinkedGraph) {
        self.active_nodes
            .retain(|active| total_out_for_branch(jt, active) > 0);
    }

    fn active_tree_branches(&self, jt: &JunctionTreeLinkedGraph) -> Vec<(i32, Vec<(Edge, i32)>)> {
        self.active_nodes
            .iter()
            .filter_map(|active| {
                let tree = jt.junction_tree_for_node(active.tree_vertex)?;
                let mut node = &tree.root;
                for &e in &active.path_edges {
                    node = node.children().get(&e)?;
                }
                let total_out = total_out_for_branch(jt, active);
                let branches: Vec<(Edge, i32)> = node
                    .children()
                    .iter()
                    .map(|(&e, c)| (e, c.evidence_count()))
                    .collect();
                Some((total_out, branches))
            })
            .collect()
    }

    fn has_junction_tree_evidence(&self) -> bool {
        !self.active_nodes.is_empty()
    }

    fn mark_trees_as_visited(&mut self, tree_vertices: &[usize]) {
        for &v in tree_vertices {
            self.visited_trees.insert(v);
        }
    }
}

fn total_out_for_branch(jt: &JunctionTreeLinkedGraph, active: &ActiveNode) -> i32 {
    let tree = match jt.junction_tree_for_node(active.tree_vertex) {
        Some(t) => t,
        None => return 0,
    };
    let mut node = &tree.root;
    for &e in &active.path_edges {
        node = match node.children().get(&e) {
            Some(n) => n,
            None => return 0,
        };
    }
    node.children().values().map(|c| c.evidence_count()).sum()
}

#[derive(Debug, Clone)]
struct JtPath {
    start: usize,
    edges: Vec<Edge>,
    score: f64,
    is_reference: bool,
    jt_manager: JunctionTreeSet,
    decision_edges_since_jt: i32,
}

impl JtPath {
    fn new(start: usize) -> Self {
        Self {
            start,
            edges: Vec::new(),
            score: 0.0,
            is_reference: false,
            jt_manager: JunctionTreeSet::new(),
            decision_edges_since_jt: 0,
        }
    }

    fn last_vertex(&self, _graph: &AssemblyGraph) -> usize {
        self.edges.last().map(|&(_, to)| to).unwrap_or(self.start)
    }

    fn vertices(&self, _graph: &AssemblyGraph) -> Vec<usize> {
        let mut v = vec![self.start];
        for &(_, to) in &self.edges {
            v.push(to);
        }
        v
    }

    fn contains_vertex(&self, graph: &AssemblyGraph, vertex: usize) -> bool {
        self.vertices(graph).contains(&vertex)
    }

    fn has_junction_tree_evidence(&self) -> bool {
        self.jt_manager.has_junction_tree_evidence()
    }

    fn was_last_edge_jt(&self) -> bool {
        self.decision_edges_since_jt == 0
    }

    fn extend_with_edges(
        &self,
        graph: &AssemblyGraph,
        jt: &JunctionTreeLinkedGraph,
        new_edges: &[Edge],
        edge_mult: u32,
        total_out: u32,
        from_jt: bool,
        edge_penalty_override: Option<f64>,
    ) -> Self {
        // Avoid `Vec::clone` on the k-best/junction frontier (OOM amplifier).
        let mut edges = Vec::with_capacity(self.edges.len() + new_edges.len());
        edges.extend_from_slice(&self.edges);
        let mut is_reference = self.is_reference;
        let mut penalty = edge_penalty_override.unwrap_or(0.0);
        if edge_penalty_override.is_none() && !new_edges.is_empty() {
            penalty = log_penalty(edge_mult, total_out);
        }
        for &e in new_edges {
            edges.push(e);
        }
        if let Some(&e) = new_edges.last() {
            is_reference = is_reference && graph.edge_is_ref(e.0, e.1);
        }
        let mut jt_manager = JunctionTreeSet::clone_from(&self.jt_manager);
        if let Some(&last) = new_edges.last() {
            jt_manager.traverse_edge_for_all_trees(jt, last);
        }
        let decision_edges_since_jt = if jt_manager.has_junction_tree_evidence() {
            0
        } else if from_jt {
            0
        } else if edge_penalty_override == Some(0.0) {
            self.decision_edges_since_jt
        } else if new_edges.is_empty() {
            self.decision_edges_since_jt
        } else {
            self.decision_edges_since_jt + 1
        };
        Self {
            start: self.start,
            edges,
            score: self.score + penalty,
            is_reference,
            jt_manager,
            decision_edges_since_jt,
        }
    }

    fn add_junction_tree(&mut self, tree_vertex: usize, tree: &ThreadingTree) {
        if self.jt_manager.add_junction_tree(tree_vertex, tree) {
            self.decision_edges_since_jt = 0;
        }
    }

    fn has_stopping_evidence(
        &mut self,
        jt: &JunctionTreeLinkedGraph,
        weight_threshold: i32,
    ) -> bool {
        self.jt_manager.prune_empty_nodes(jt);
        for (total_out, branches) in self.jt_manager.active_tree_branches(jt) {
            if branches
                .iter()
                .any(|&(edge, _)| edge == jt.symbolic_end_edge)
            {
                return true;
            }
            if total_out >= weight_threshold {
                return false;
            }
        }
        true
    }

    fn applicable_next_edges(
        &mut self,
        graph: &AssemblyGraph,
        jt: &JunctionTreeLinkedGraph,
        fork_vertex: usize,
        chain: &[Edge],
        outgoing: &[usize],
        weight_threshold: i32,
    ) -> Vec<Self> {
        let mut output = Vec::new();
        let mut edges_accounted: HashSet<Edge> = HashSet::new();
        let outgoing_edges: Vec<Edge> = outgoing.iter().map(|&to| (fork_vertex, to)).collect();

        self.jt_manager.prune_empty_nodes(jt);
        let tree_nodes = self.jt_manager.active_tree_branches(jt);
        for (total_out, branches) in tree_nodes {
            // GATK: at high-confidence JT forks only the strongest branch is traversed on the heap;
            // weaker branches are recovered via pivotal edges when `recover_edges` is enabled.
            let branch_list: Vec<(Edge, i32)> = if total_out >= weight_threshold {
                branches
                    .into_iter()
                    .filter(|&(edge, _)| {
                        outgoing_edges.contains(&edge) && edge != jt.symbolic_end_edge
                    })
                    // Evidence first; edge `(from,to)` breaks ties (BTreeMap children order).
                    .max_by(|a, b| a.1.cmp(&b.1).then_with(|| a.0.cmp(&b.0)))
                    .into_iter()
                    .collect()
            } else {
                branches
            };
            for (edge, child_evidence) in branch_list {
                if !outgoing_edges.contains(&edge) {
                    continue;
                }
                if edge == jt.symbolic_end_edge {
                    continue;
                }
                if edges_accounted.contains(&edge) {
                    continue;
                }
                edges_accounted.insert(edge);
                let mut chain_copy = chain.to_vec();
                chain_copy.push(edge);
                output.push(self.extend_with_edges(
                    graph,
                    jt,
                    &chain_copy,
                    child_evidence.max(0) as u32,
                    total_out as u32,
                    true,
                    None,
                ));
            }
            if total_out >= weight_threshold {
                return output;
            }
        }

        let total_outgoing_multiplicity: u32 = outgoing_edges
            .iter()
            .filter_map(|&e| graph.edge_support(e.0, e.1))
            .sum();

        for &edge in &outgoing_edges {
            if edges_accounted.contains(&edge) {
                continue;
            }
            if total_outgoing_multiplicity == 0 {
                continue;
            }
            let mult = graph.edge_support(edge.0, edge.1).unwrap_or(0);
            if mult == 0 {
                continue;
            }
            if graph.edge_is_ref(edge.0, edge.1)
                && mult == 1
                && !(edges_accounted.is_empty() && outgoing_edges.len() < 2)
            {
                continue;
            }
            let mut chain_copy = chain.to_vec();
            chain_copy.push(edge);
            output.push(self.extend_with_edges(
                graph,
                jt,
                &chain_copy,
                mult,
                total_outgoing_multiplicity,
                false,
                None,
            ));
        }
        output
    }
}

struct HeapItem {
    score: f64,
    tie: usize,
    path: JtPath,
}

impl Ord for HeapItem {
    fn cmp(&self, other: &Self) -> Ordering {
        self.score
            .partial_cmp(&other.score)
            .unwrap_or(Ordering::Equal)
            .then_with(|| other.tie.cmp(&self.tie))
    }
}

impl PartialOrd for HeapItem {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl PartialEq for HeapItem {
    fn eq(&self, other: &Self) -> bool {
        self.score.to_bits() == other.score.to_bits() && self.tie == other.tie
    }
}

impl Eq for HeapItem {}

/// GATK `JunctionTreeKBestHaplotypeFinder.findBestHaplotypes`.
pub fn find_junction_best_haplotypes(
    jt: &JunctionTreeLinkedGraph,
    max_haplotypes: usize,
    jt_evidence_threshold: i32,
    recover_edges: bool,
) -> GatkResult<Vec<JunctionKBestPath>> {
    if max_haplotypes == 0 {
        return Ok(Vec::new());
    }
    let graph = &jt.graph;
    let source = jt
        .reference_source()
        .ok_or_else(|| GatkError::algorithm("jt-kbest: no reference source"))?;
    let sink = jt
        .reference_sink()
        .ok_or_else(|| GatkError::algorithm("jt-kbest: no reference sink"))?;
    let sinks: HashSet<usize> = HashSet::from([sink]);
    let _ = sink;

    let mut unvisited_pivotal: Vec<Edge> = if recover_edges {
        create_pivotal_edges_in_topological_order(jt)
    } else {
        Vec::new()
    };

    let mut result: Vec<JtPath> = Vec::new();
    let mut heap: BinaryHeap<HeapItem> = BinaryHeap::new();
    heap.push(HeapItem {
        score: 0.0,
        tie: 0,
        path: JtPath::new(source),
    });

    let mut kmer_chain_cache: HashMap<usize, Vec<Edge>> = HashMap::new();
    let mut tie_counter = 1usize;

    let mut debug_closes = 0usize;
    let mut max_last_vertex = 0usize;
    while result.len() < max_haplotypes && (!heap.is_empty() || !unvisited_pivotal.is_empty()) {
        let max_queue = if result.is_empty() {
            DEFAULT_MAX_PATHS_TO_CONSIDER_WITHOUT_RESULT
        } else {
            DEFAULT_MAX_PATHS_TO_EVER_CONSIDER
        };
        if heap.len() > max_queue {
            break;
        }

        if heap.is_empty() {
            enqueue_next_pivotal_edge(
                graph,
                jt,
                &mut unvisited_pivotal,
                &result,
                &mut heap,
                &mut tie_counter,
            );
            continue;
        }

        let item = heap.pop().expect("non-empty heap");
        let mut path = item.path;
        max_last_vertex = max_last_vertex.max(path.last_vertex(graph));

        if path.decision_edges_since_jt > DEFAULT_MAX_ACCEPTABLE_DECISION_EDGES_WITHOUT_JT_GUIDANCE
        {
            continue;
        }

        let path_last = path.last_vertex(graph);
        let mut vertex_to_extend = path_last;
        let mut outgoing = graph.outgoing_nodes(vertex_to_extend);

        let cached_chain = kmer_chain_cache.get(&path_last).cloned();
        let chain = if let Some(c) = cached_chain {
            vertex_to_extend = c.last().map(|&(_, t)| t).unwrap_or(vertex_to_extend);
            outgoing = graph.outgoing_nodes(vertex_to_extend);
            c
        } else {
            let mut chain = Vec::new();
            while outgoing.len() == 1
                && jt.junction_tree_for_node(vertex_to_extend).is_none()
                && !sinks.contains(&vertex_to_extend)
            {
                let from = vertex_to_extend;
                let to = outgoing[0];
                let edge = (from, to);
                if chain.contains(&edge) {
                    break;
                }
                chain.push(edge);
                vertex_to_extend = to;
                outgoing = graph.outgoing_nodes(vertex_to_extend);
            }
            // CLONE: needed because owned HashMap/BTree/HashSet key or value.
            kmer_chain_cache.insert(path_last, chain.clone());
            chain
        };

        if let Some(tree) = jt.junction_tree_for_node(vertex_to_extend) {
            path.add_junction_tree(vertex_to_extend, tree);
        }

        process_vertex(
            graph,
            jt,
            &sinks,
            &mut path,
            vertex_to_extend,
            &outgoing,
            &chain,
            jt_evidence_threshold,
            recover_edges,
            &mut result,
            &mut unvisited_pivotal,
            &mut heap,
            &mut tie_counter,
            &mut debug_closes,
        );
    }

    #[cfg(test)]
    if result.is_empty() {
        eprintln!(
            "jt-kbest: empty result heap={} pivotal={} debug_closes={debug_closes} max_last={max_last_vertex} sink={sink}",
            heap.len(),
            unvisited_pivotal.len()
        );
    }
    let mut out: Vec<JunctionKBestPath> = result
        .into_iter()
        .map(|p| JunctionKBestPath {
            start: p.start,
            edges: p.edges,
            score: p.score,
            is_reference: p.is_reference,
        })
        .collect();
    out.sort_by(|a, b| {
        b.score
            .partial_cmp(&a.score)
            .unwrap_or(Ordering::Equal)
            .then_with(|| b.bases(graph).cmp(&a.bases(graph)))
    });
    Ok(out)
}

#[allow(clippy::too_many_arguments)]
fn process_vertex(
    graph: &AssemblyGraph,
    jt: &JunctionTreeLinkedGraph,
    sinks: &HashSet<usize>,
    path: &mut JtPath,
    vertex_to_extend: usize,
    outgoing: &[usize],
    chain: &[Edge],
    jt_evidence_threshold: i32,
    _recover_edges: bool,
    result: &mut Vec<JtPath>,
    unvisited_pivotal: &mut Vec<Edge>,
    heap: &mut BinaryHeap<HeapItem>,
    tie_counter: &mut usize,
    debug_closes: &mut usize,
) {
    #[cfg(test)]
    static DEBUG_PROC: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);
    #[cfg(test)]
    {
        let n = DEBUG_PROC.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        if n < 5 {
            eprintln!(
                "jt-kbest: process v={vertex_to_extend} path_last={} out={} chain={} sinks_hit={}",
                path.last_vertex(graph),
                outgoing.len(),
                chain.len(),
                sinks.contains(&vertex_to_extend)
            );
        }
    }
    if sinks.contains(&vertex_to_extend) {
        let stopping = path.has_stopping_evidence(jt, jt_evidence_threshold);
        #[cfg(test)]
        if !stopping {
            eprintln!(
                "jt-kbest: at sink no stopping jt_evidence={} has_jt={}",
                path.has_junction_tree_evidence(),
                path.has_junction_tree_evidence()
            );
        }
        if stopping {
            *debug_closes += 1;
            let closed = if chain.is_empty() {
                path.clone()
            } else {
                path.extend_with_edges(graph, jt, chain, 0, 1, false, Some(0.0))
            };
            match closed.vertices(graph).first() {
                Some(&v) if v == path.start => result.push(closed),
                Some(_) => {
                    #[cfg(test)]
                    eprintln!(
                        "jt-kbest: reject close (non-start vertex), start={}",
                        path.start
                    );
                }
                None => {}
            }
        }
        for e in &path.edges {
            if let Some(pos) = unvisited_pivotal.iter().position(|p| p == e) {
                unvisited_pivotal.remove(pos);
            }
        }
    }

    if outgoing.len() > 1 {
        let jt_paths = path.applicable_next_edges(
            graph,
            jt,
            vertex_to_extend,
            chain,
            outgoing,
            jt_evidence_threshold,
        );
        #[cfg(test)]
        if DEBUG_PROC.load(std::sync::atomic::Ordering::Relaxed) <= 5 {
            eprintln!(
                "jt-kbest: fork jt_paths={} filtered_incoming={}",
                jt_paths.len(),
                jt_paths.len()
            );
        }
        let filtered: Vec<JtPath> = jt_paths
            .into_iter()
            .filter(|p| {
                p.has_junction_tree_evidence()
                    || p.was_last_edge_jt()
                    || p.vertices(graph)
                        .iter()
                        .filter(|&&v| v == p.last_vertex(graph))
                        .count()
                        <= DEFAULT_MAX_ACCEPTABLE_REPETITIONS_OF_A_KMER_IN_A_PATH
            })
            .collect();
        #[cfg(test)]
        if DEBUG_PROC.load(std::sync::atomic::Ordering::Relaxed) <= 5 {
            eprintln!("jt-kbest: after_filter={}", filtered.len());
        }
        for p in filtered {
            push_path(heap, p, tie_counter);
        }
    } else if !outgoing.is_empty() {
        let final_vertex = vertex_to_extend;
        let repeat_count = path
            .vertices(graph)
            .iter()
            .filter(|&&v| v == final_vertex)
            .count();
        if path.has_junction_tree_evidence()
            || repeat_count <= DEFAULT_MAX_ACCEPTABLE_REPETITIONS_OF_A_KMER_IN_A_PATH
        {
            let mut chain_copy = chain.to_vec();
            chain_copy.push((vertex_to_extend, outgoing[0]));
            let extended = path.extend_with_edges(graph, jt, &chain_copy, 0, 1, false, Some(0.0));
            push_path(heap, extended, tie_counter);
        }
    }
}

fn push_path(heap: &mut BinaryHeap<HeapItem>, path: JtPath, tie_counter: &mut usize) {
    let tie = *tie_counter;
    *tie_counter += 1;
    heap.push(HeapItem {
        score: path.score,
        tie,
        path,
    });
}

fn enqueue_next_pivotal_edge(
    graph: &AssemblyGraph,
    jt: &JunctionTreeLinkedGraph,
    unvisited_pivotal: &mut Vec<Edge>,
    result: &[JtPath],
    heap: &mut BinaryHeap<HeapItem>,
    tie_counter: &mut usize,
) {
    let Some(first_edge) = unvisited_pivotal.first().copied() else {
        return;
    };
    unvisited_pivotal.remove(0);
    let pivotal_vertex = first_edge.0;

    let best = result
        .iter()
        .filter(|p| p.contains_vertex(graph, pivotal_vertex))
        .max_by(|a, b| a.score.partial_cmp(&b.score).unwrap_or(Ordering::Equal));
    let Some(best_matching) = best else {
        return;
    };

    let incoming: Vec<Edge> = best_matching
        .edges
        .iter()
        .copied()
        .filter(|&e| e.1 == pivotal_vertex)
        .collect();
    if incoming.is_empty() {
        return;
    }
    let last_incoming = incoming.last().copied().expect("non-empty");
    let split_idx = best_matching
        .edges
        .iter()
        .position(|&e| e == last_incoming)
        .expect("edge in path")
        + 1;
    let mut edges_before: Vec<Edge> = best_matching.edges[..split_idx].to_vec();
    edges_before.push(first_edge);

    let mut path_to_add = JtPath {
        start: best_matching.start,
        edges: edges_before,
        score: best_matching.score,
        is_reference: best_matching.is_reference,
        jt_manager: JunctionTreeSet::new(),
        decision_edges_since_jt: 0,
    };

    let trees_passed: Vec<usize> = path_to_add
        .vertices(graph)
        .into_iter()
        .filter(|v| jt.junction_tree_for_node(*v).is_some())
        .collect();
    path_to_add.jt_manager.mark_trees_as_visited(&trees_passed);
    push_path(heap, path_to_add, tie_counter);
}

fn create_pivotal_edges_in_topological_order(jt: &JunctionTreeLinkedGraph) -> Vec<Edge> {
    let graph = &jt.graph;
    let source = match jt.reference_source() {
        Some(s) => s,
        None => return Vec::new(),
    };
    let mut visited_edges: HashSet<Edge> = HashSet::new();
    let mut output: Vec<Edge> = Vec::new();
    let mut queue: BinaryHeap<std::cmp::Reverse<(usize, Edge)>> = BinaryHeap::new();

    for to in graph.outgoing_nodes(source) {
        queue.push(std::cmp::Reverse((0, (source, to))));
    }

    while let Some(std::cmp::Reverse((score, edge))) = queue.pop() {
        let target = edge.1;
        let outgoing: Vec<usize> = graph.outgoing_nodes(target);
        if outgoing.len() > 1 {
            for &to in &outgoing {
                let e = (target, to);
                let mult = graph.edge_support(e.0, e.1).unwrap_or(0);
                let is_ref_only = graph.edge_is_ref(e.0, e.1) && mult == 1;
                if !is_ref_only && !visited_edges.contains(&e) {
                    output.push(e);
                }
            }
        }
        for to in outgoing {
            let e = (target, to);
            if visited_edges.insert(e) {
                queue.push(std::cmp::Reverse((score + 1, e)));
            }
        }
    }
    output
}
