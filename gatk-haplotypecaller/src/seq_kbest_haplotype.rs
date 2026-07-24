//! GATK `GraphBasedKBestHaplotypeFinder` on [`SeqGraph`] (`Path.getBases` stitching).

use crate::kbest_haplotype::{log_penalty, KBestPath};
use crate::seq_graph::SeqGraph;
use gatk_common::{GatkError, GatkResult};
use std::cmp::Ordering;
use std::collections::{BinaryHeap, HashSet};

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

    fn extend(&self, graph: &SeqGraph, to: usize, edge_support: u32, total_outgoing: u32) -> Self {
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

/// K-best paths on a cleaned sequence graph (GATK `GraphBasedKBestHaplotypeFinder`).
pub fn find_best_haplotypes_seq_graph(
    graph: &SeqGraph,
    max_number_of_haplotypes: usize,
) -> GatkResult<Vec<KBestPath>> {
    if max_number_of_haplotypes == 0 {
        return Ok(Vec::new());
    }
    let source = graph
        .reference_source_vertex()
        .ok_or_else(|| GatkError::algorithm("seq kbest: no reference source vertex"))?;
    let sink = graph
        .reference_sink_vertex()
        .ok_or_else(|| GatkError::algorithm("seq kbest: no reference sink vertex"))?;
    let sinks: HashSet<usize> = HashSet::from([sink]);

    let mut result = Vec::new();
    let mut heap: BinaryHeap<HeapItem> = BinaryHeap::new();
    heap.push(HeapItem {
        score_bits: 0.0_f64.to_bits(),
        tie: 0,
        path: PathState::new(source),
    });

    let mut vertex_counts = vec![0usize; graph.node_count()];

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
        } else if vertex_counts[path.last] < max_number_of_haplotypes {
            vertex_counts[path.last] += 1;
            let outs = graph.outgoing_nodes(path.last);
            let total: u32 = outs
                .iter()
                .filter_map(|&t| graph.edge_support(path.last, t))
                .sum();
            for to in outs {
                if let Some(support) = graph.edge_support(path.last, to) {
                    let extended = path.extend(graph, to, support, total);
                    let tie = graph
                        .path_bases_bytes(extended.start, &extended.edges)
                        .len();
                    heap.push(HeapItem {
                        score_bits: extended.score.to_bits(),
                        tie,
                        path: extended,
                    });
                }
            }
        }
    }

    result.sort_by(|a, b| {
        b.score
            .partial_cmp(&a.score)
            .unwrap_or(Ordering::Equal)
            .then_with(|| {
                let ab = graph.path_bases_bytes(a.start, &a.edges);
                let bb = graph.path_bases_bytes(b.start, &b.edges);
                bb.cmp(&ab)
            })
    });
    Ok(result)
}
