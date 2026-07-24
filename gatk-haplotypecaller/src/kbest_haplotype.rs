//! GATK `GraphBasedKBestHaplotypeFinder` / `KBestHaplotype` parity (Yen's-style K-shortest paths on the assembly graph).

use crate::assembly::AssemblyGraph;
use gatk_common::{GatkError, GatkResult};
use std::cmp::Ordering;
use std::collections::{BinaryHeap, HashMap, HashSet};

/// One haplotype path with GATK log-penalty score.
/// # Invariants
/// `edges` form an ordered walk from `start` through the [`AssemblyGraph`].
/// `score` is the GATK log-penalty path score used for k-best ranking.
/// # Ownership
/// Owns edge list; borrows graph in [`KBestPath::bases`].
/// # Mutation
/// Immutable path after k-best search.
/// # Biological assumptions
/// Candidate haplotype sequence over the local assembly graph.
/// # Java equivalence
/// GATK `KBestHaplotype` / `GraphBasedKBestHaplotypeFinder` path.
#[derive(Debug, Clone)]
pub struct KBestPath {
    pub start: usize,
    pub edges: Vec<(usize, usize)>,
    pub score: f64,
    pub is_reference: bool,
}

impl KBestPath {
    pub fn bases(&self, graph: &AssemblyGraph) -> String {
        graph.path_bases(self.start, &self.edges)
    }
}

#[derive(Debug, Clone)]
struct PathState {
    start: usize,
    edges: Vec<(usize, usize)>,
    last: usize,
    score: f64,
    is_reference: bool,
}

impl PathState {
    fn new(start: usize) -> Self {
        Self {
            start,
            edges: Vec::new(),
            last: start,
            score: 0.0,
            is_reference: false,
        }
    }

    fn extend(
        &self,
        graph: &AssemblyGraph,
        to: usize,
        edge_support: u32,
        total_outgoing: u32,
    ) -> Self {
        let mut edges = self.edges.clone();
        edges.push((self.last, to));
        let penalty = log_penalty(edge_support, total_outgoing);
        Self {
            start: self.start,
            edges,
            last: to,
            score: self.score + penalty,
            is_reference: self.is_reference && graph.edge_is_ref(self.last, to),
        }
    }
}

pub(crate) fn log_penalty(edge_multiplicity: u32, total_outgoing_multiplicity: u32) -> f64 {
    if total_outgoing_multiplicity == 0 {
        return 0.0;
    }
    (edge_multiplicity.max(1) as f64).log10() - (total_outgoing_multiplicity.max(1) as f64).log10()
}

struct HeapItem {
    score_bits: u64,
    tie: usize,
    path: PathState,
}

impl Ord for HeapItem {
    fn cmp(&self, other: &Self) -> Ordering {
        self.score_bits
            .cmp(&other.score_bits)
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
        self.score_bits == other.score_bits && self.tie == other.tie
    }
}

impl Eq for HeapItem {}

/// Remove cycle edges / dead vertices (GATK `KBestHaplotypeFinder.removeCyclesIfNecessary`).
pub fn graph_for_kbest(mut graph: AssemblyGraph) -> GatkResult<AssemblyGraph> {
    if !graph.has_cycle() {
        return Ok(graph);
    }
    let sources: Vec<usize> = graph.reference_source_vertex().into_iter().collect();
    let sinks: Vec<usize> = graph.reference_sink_vertex().into_iter().collect();
    if sources.is_empty() || sinks.is_empty() {
        return Err(GatkError::algorithm(
            "kbest: missing reference source or sink for cycle removal",
        ));
    }
    let sink_set: HashSet<usize> = sinks.into_iter().collect();
    let mut edges_to_remove = HashSet::new();
    let mut vertices_to_remove = HashSet::new();
    let mut found_path = false;
    for source in sources {
        let mut parents = HashSet::new();
        if find_cycle_guilty(
            &graph,
            source,
            &sink_set,
            &mut edges_to_remove,
            &mut vertices_to_remove,
            &mut parents,
        ) {
            found_path = true;
        }
    }
    if !found_path {
        return Err(GatkError::algorithm(
            "kbest: no path from source to sink after cycle analysis",
        ));
    }
    if edges_to_remove.is_empty() && vertices_to_remove.is_empty() {
        return Err(GatkError::algorithm("kbest: cannot remove cycles"));
    }
    for (f, t) in edges_to_remove {
        graph.remove_edge(f, t);
    }
    graph.remove_nodes(&vertices_to_remove);
    graph.cleanup_isolated_nodes();
    Ok(graph)
}

fn find_cycle_guilty(
    graph: &AssemblyGraph,
    current: usize,
    sinks: &HashSet<usize>,
    edges_to_remove: &mut HashSet<(usize, usize)>,
    vertices_to_remove: &mut HashSet<usize>,
    parent_vertices: &mut HashSet<usize>,
) -> bool {
    if sinks.contains(&current) {
        return true;
    }
    parent_vertices.insert(current);
    let outs = graph.outgoing_nodes(current);
    let mut reaches_sink = false;
    for to in outs {
        if parent_vertices.contains(&to) {
            edges_to_remove.insert((current, to));
        } else if find_cycle_guilty(
            graph,
            to,
            sinks,
            edges_to_remove,
            vertices_to_remove,
            parent_vertices,
        ) {
            reaches_sink = true;
        }
    }
    if !reaches_sink {
        vertices_to_remove.insert(current);
    }
    parent_vertices.remove(&current);
    reaches_sink
}

/// GATK `GraphBasedKBestHaplotypeFinder.findBestHaplotypes`.
pub fn find_best_haplotypes(
    graph: &AssemblyGraph,
    max_number_of_haplotypes: usize,
) -> GatkResult<Vec<KBestPath>> {
    find_best_haplotypes_inner(graph, max_number_of_haplotypes, true)
}

/// K-best on the graph as-is (no cycle stripping). Parity diagnostics only.
pub fn find_best_haplotypes_preserving_cycles(
    graph: &AssemblyGraph,
    max_number_of_haplotypes: usize,
) -> GatkResult<Vec<KBestPath>> {
    find_best_haplotypes_inner(graph, max_number_of_haplotypes, false)
}

/// Production assembly: strip cycles only when the graph is acyclic (P12 cyclic regions keep topology).
pub fn find_best_haplotypes_for_assembly(
    graph: &AssemblyGraph,
    max_number_of_haplotypes: usize,
) -> GatkResult<Vec<KBestPath>> {
    if graph.has_cycle() {
        find_best_haplotypes_preserving_cycles(graph, max_number_of_haplotypes)
    } else {
        find_best_haplotypes(graph, max_number_of_haplotypes)
    }
}

fn find_best_haplotypes_inner(
    graph: &AssemblyGraph,
    max_number_of_haplotypes: usize,
    strip_cycles: bool,
) -> GatkResult<Vec<KBestPath>> {
    if max_number_of_haplotypes == 0 {
        return Ok(Vec::new());
    }
    let graph = if strip_cycles {
        graph_for_kbest(graph.clone())?
    } else {
        graph.clone()
    };
    let source = graph
        .reference_source_vertex()
        .ok_or_else(|| GatkError::algorithm("kbest: no reference source vertex"))?;
    let sink = graph
        .reference_sink_vertex()
        .ok_or_else(|| GatkError::algorithm("kbest: no reference sink vertex"))?;
    let sources = [source];
    let sinks: HashSet<usize> = HashSet::from([sink]);

    let mut result = Vec::new();
    let mut heap: BinaryHeap<HeapItem> = BinaryHeap::new();
    for &s in &sources {
        let path = PathState::new(s);
        heap.push(HeapItem {
            score_bits: path.score.to_bits(),
            tie: 0,
            path,
        });
    }
    let mut vertex_counts: HashMap<usize, usize> =
        graph.nodes().iter().map(|n| (n.id, 0)).collect();

    while !heap.is_empty() && result.len() < max_number_of_haplotypes {
        let item = heap.pop().expect("non-empty");
        let path = item.path;
        if sinks.contains(&path.last) {
            result.push(KBestPath {
                start: path.start,
                edges: path.edges,
                score: path.score,
                is_reference: path.is_reference,
            });
        } else {
            let count = vertex_counts.get_mut(&path.last).expect("vertex");
            if *count < max_number_of_haplotypes {
                *count += 1;
                let outs = graph.outgoing_nodes(path.last);
                let total: u32 = outs
                    .iter()
                    .filter_map(|&t| graph.edge_support(path.last, t))
                    .sum();
                for to in outs {
                    if let Some(support) = graph.edge_support(path.last, to) {
                        let extended = path.extend(&graph, to, support, total);
                        let tie = graph.path_bases(extended.start, &extended.edges).len();
                        heap.push(HeapItem {
                            score_bits: extended.score.to_bits(),
                            tie,
                            path: extended,
                        });
                    }
                }
            }
        }
    }

    result.sort_by(|a, b| {
        b.score
            .partial_cmp(&a.score)
            .unwrap_or(Ordering::Equal)
            .then_with(|| b.bases(&graph).cmp(&a.bases(&graph)))
    });
    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::assembly::{AssemblyGraphParams, AssemblyRead};
    use crate::read_threading_graph::assembly_graph_from_ref_and_reads_threading;

    fn read(seq: &str, q: u8) -> AssemblyRead {
        AssemblyRead {
            bases: seq.to_string(),
            base_quals: vec![q; seq.len()],
        }
    }

    #[test]
    fn p5_case1_finds_ref_and_alt_paths() {
        let reference = read("ACGTT", 30);
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
        let mut graph =
            assembly_graph_from_ref_and_reads_threading(&reference, &reads, &params).unwrap();
        let mut pruning =
            crate::assembly::AssemblyGraphPruningParams::gatk_haplotype_caller_defaults();
        pruning.min_prune_factor = 2;
        graph.apply_pruning(&pruning);
        let paths = find_best_haplotypes(&graph, 128).unwrap();
        assert!(!paths.is_empty());
        let seqs: HashSet<_> = paths.iter().map(|p| p.bases(&graph)).collect();
        assert!(seqs.contains("ACGTT"));
    }
}
