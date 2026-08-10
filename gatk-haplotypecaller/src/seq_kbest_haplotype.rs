//! GATK `GraphBasedKBestHaplotypeFinder` on [`SeqGraph`] (`Path.getBases` stitching).

use crate::kbest_haplotype::{log_penalty, KBestPath};
use crate::seq_graph::SeqGraph;
use gatk_common::{GatkError, GatkResult};
use std::cmp::Ordering;
use std::collections::{BinaryHeap, HashSet};

/// Same bound as the read-threading k-best finder — see `kbest_haplotype`.
const MAX_KBEST_HEAP_PATHS: usize = 1_024;
const MAX_KBEST_PATH_EDGES: usize = 4_096;

#[derive(Debug, Clone)]
struct PathState {
    start: usize,
    edges: Vec<(usize, usize)>,
    last: usize,
    score: f64,
    is_reference: bool,
    /// Edge count as a cheap heap tie-break (full bases sort happens at the end).
    edge_count: usize,
}

impl PathState {
    fn new(start: usize) -> Self {
        Self {
            start,
            edges: Vec::new(),
            last: start,
            score: 0.0,
            is_reference: false,
            edge_count: 0,
        }
    }

    fn extend(&self, graph: &SeqGraph, to: usize, edge_support: u32, total_outgoing: u32) -> Self {
        // Same as read-threading k-best: no `Vec::clone` on the frontier path.
        let mut edges = Vec::with_capacity(self.edges.len() + 1);
        edges.extend_from_slice(&self.edges);
        edges.push((self.last, to));
        let penalty = log_penalty(edge_support, total_outgoing);
        Self {
            start: self.start,
            edges,
            last: to,
            score: self.score + penalty,
            is_reference: self.is_reference && graph.edge_is_ref(self.last, to),
            edge_count: self.edge_count + 1,
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
    // Same Peak bound as read-threading k-best (50k → 12k for dense GIAB shards).
    const MAX_KBEST_EXPANSIONS: usize = 12_000;
    let mut expansions = 0usize;

    while !heap.is_empty() && result.len() < max_number_of_haplotypes {
        if crate::runtime_config::hc_rss_abort_triggered() {
            crate::runtime_config::rss_trace_checkpoint(
                "seq_kbest_rss_abort",
                &format!("expansions={expansions} results={}", result.len()),
            );
            break;
        }
        let item = heap.pop().expect("non-empty");
        let path = item.path;
        if sinks.contains(&path.last) {
            result.push(KBestPath {
                start: path.start,
                edges: path.edges,
                score: path.score,
                is_reference: path.is_reference,
            });
            continue;
        }
        if path.edges.len() >= MAX_KBEST_PATH_EDGES
            || expansions >= MAX_KBEST_EXPANSIONS
            || heap.len() >= MAX_KBEST_HEAP_PATHS
        {
            continue;
        }
        if vertex_counts[path.last] < max_number_of_haplotypes {
            vertex_counts[path.last] += 1;
            expansions += 1;
            let outs = graph.outgoing_nodes(path.last);
            let total: u32 = outs
                .iter()
                .filter_map(|&t| graph.edge_support(path.last, t))
                .sum();
            for to in outs {
                if let Some(support) = graph.edge_support(path.last, to) {
                    if heap.len() >= MAX_KBEST_HEAP_PATHS {
                        break;
                    }
                    let extended = path.extend(graph, to, support, total);
                    heap.push(HeapItem {
                        score_bits: extended.score.to_bits(),
                        tie: extended.edge_count,
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
