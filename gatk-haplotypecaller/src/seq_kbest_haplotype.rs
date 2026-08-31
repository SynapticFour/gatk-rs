//! GATK `GraphBasedKBestHaplotypeFinder` on [`SeqGraph`] (`Path.getBases` stitching).

use crate::kbest_haplotype::{cmp_graph_kbest_score, log_penalty, KBestPath};
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
    score: f64,
    tie: usize,
    path: PathState,
}

impl Ord for HeapItem {
    fn cmp(&self, other: &Self) -> Ordering {
        cmp_graph_kbest_score(self.score, other.score).then_with(|| other.tie.cmp(&self.tie))
    }
}

impl PartialOrd for HeapItem {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl PartialEq for HeapItem {
    fn eq(&self, other: &Self) -> bool {
        cmp_graph_kbest_score(self.score, other.score).is_eq() && self.tie == other.tie
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
        score: 0.0,
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
                        score: extended.score,
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::assembly::AssemblyGraph;

    /// One bubble: competing branches with supports 40 vs 5. No genomic coordinates.
    fn diamond_40_vs_5() -> crate::seq_graph::SeqGraph {
        let mut g = AssemblyGraph::new(3).unwrap();
        let src = g.ensure_node(b"AAA");
        let ref_mid = g.ensure_node(b"AAC");
        let alt_mid = g.ensure_node(b"AAG");
        let snk = g.ensure_node(b"ACT");
        g.add_edge_support(src, ref_mid, 40);
        g.add_edge_support(src, alt_mid, 5);
        g.add_edge_support(ref_mid, snk, 40);
        g.add_edge_support(alt_mid, snk, 5);
        g.ref_edges.insert((src, ref_mid));
        g.ref_edges.insert((ref_mid, snk));
        g.ref_nodes.insert(src);
        g.ref_nodes.insert(ref_mid);
        g.ref_nodes.insert(snk);
        g.ref_source_kmer = Some(std::sync::Arc::from(b"AAA".as_slice()));
        crate::seq_graph::SeqGraph::from_assembly_graph(&g)
    }

    fn path_has_kmer(seq: &crate::seq_graph::SeqGraph, path: &KBestPath, kmer: &[u8]) -> bool {
        seq.path_bases_bytes(path.start, &path.edges)
            .windows(kmer.len())
            .any(|w| w == kmer)
    }

    /// 6R.52 regression (coordinate-free).
    ///
    /// Java `GraphBasedKBestHaplotypeFinder.findBestHaplotypes` polls **highest**
    /// `KBestHaplotype.score` first (`comparingDouble(score).reversed()`). Scores are
    /// accumulated `log10(edgeMult / outMult)` (≤ 0). k=1 must return the high-support
    /// branch (`AAC`, 40), not `AAG` (5).
    #[test]
    fn seq_kbest_k1_returns_high_multiplicity_branch() {
        let seq = diamond_40_vs_5();
        let paths = find_best_haplotypes_seq_graph(&seq, 1).expect("kbest");
        assert_eq!(paths.len(), 1);
        let bases = seq.path_bases_bytes(paths[0].start, &paths[0].edges);
        assert!(
            path_has_kmer(&seq, &paths[0], b"AAC"),
            "k=1 must return the 40-support branch first; got {:?}",
            String::from_utf8_lossy(&bases)
        );
        assert!(
            !path_has_kmer(&seq, &paths[0], b"AAG"),
            "k=1 must not return the 5-support branch first; got {:?}",
            String::from_utf8_lossy(&bases)
        );
    }

    #[test]
    fn seq_kbest_k2_high_support_path_is_first() {
        let seq = diamond_40_vs_5();
        let paths = find_best_haplotypes_seq_graph(&seq, 2).expect("kbest");
        assert_eq!(paths.len(), 2);
        assert!(path_has_kmer(&seq, &paths[0], b"AAC"));
        assert!(path_has_kmer(&seq, &paths[1], b"AAG"));
        assert!(
            paths[0].score > paths[1].score,
            "high-support path must have the higher (less negative) score"
        );
    }

    /// Equal-support diamond: scores match; Java PQ then uses reversed `getBases`.
    /// When both paths are collected, Rust's post-sort matches that result order.
    fn diamond_equal_support() -> crate::seq_graph::SeqGraph {
        let mut g = AssemblyGraph::new(3).unwrap();
        let src = g.ensure_node(b"AAA");
        let mid_c = g.ensure_node(b"AAC");
        let mid_g = g.ensure_node(b"AAG");
        let snk = g.ensure_node(b"ACT");
        g.add_edge_support(src, mid_c, 10);
        g.add_edge_support(src, mid_g, 10);
        g.add_edge_support(mid_c, snk, 10);
        g.add_edge_support(mid_g, snk, 10);
        g.ref_edges.insert((src, mid_c));
        g.ref_edges.insert((mid_c, snk));
        g.ref_nodes.insert(src);
        g.ref_nodes.insert(mid_c);
        g.ref_nodes.insert(snk);
        g.ref_source_kmer = Some(std::sync::Arc::from(b"AAA".as_slice()));
        crate::seq_graph::SeqGraph::from_assembly_graph(&g)
    }

    #[test]
    fn seq_kbest_equal_score_k2_sorts_by_reverse_bases() {
        let seq = diamond_equal_support();
        let paths = find_best_haplotypes_seq_graph(&seq, 2).expect("kbest");
        assert_eq!(paths.len(), 2);
        assert!((paths[0].score - paths[1].score).abs() < 1e-12);
        let b0 = seq.path_bases_bytes(paths[0].start, &paths[0].edges);
        let b1 = seq.path_bases_bytes(paths[1].start, &paths[1].edges);
        assert!(
            b0 >= b1,
            "Java thenComparing(getBases, reversed) + Rust post-sort: lexicographically larger first; got {:?} then {:?}",
            String::from_utf8_lossy(&b0),
            String::from_utf8_lossy(&b1)
        );
    }

    /// Heap tie during search is still `edge_count`, not Java `getBases`.
    /// Both diamond arms have 2 edges, so k=1 winner is not pinned to Java's bases order.
    /// This test only records that exactly one equal-score arm is returned (K truncation).
    #[test]
    fn seq_kbest_equal_score_k1_returns_exactly_one_arm() {
        let seq = diamond_equal_support();
        let paths = find_best_haplotypes_seq_graph(&seq, 1).expect("kbest");
        assert_eq!(paths.len(), 1);
        let has_c = path_has_kmer(&seq, &paths[0], b"AAC");
        let has_g = path_has_kmer(&seq, &paths[0], b"AAG");
        assert!(has_c ^ has_g, "k=1 must pick exactly one equal-score arm");
    }

    #[test]
    fn seq_kbest_heap_highest_finite_score_is_greater() {
        let high = HeapItem {
            score: -0.1,
            tie: 99,
            path: PathState::new(0),
        };
        let low = HeapItem {
            score: -0.9,
            tie: 1,
            path: PathState::new(0),
        };
        assert!(
            high > low,
            "BinaryHeap must pop the higher numerical score first"
        );
        assert_eq!(cmp_graph_kbest_score(-0.1, -0.9), Ordering::Greater);
        assert_eq!(cmp_graph_kbest_score(0.0, -1.0), Ordering::Greater);
        let nan = HeapItem {
            score: f64::NAN,
            tie: 0,
            path: PathState::new(0),
        };
        let finite = HeapItem {
            score: -100.0,
            tie: 0,
            path: PathState::new(0),
        };
        assert!(finite > nan, "NaN must not occupy the BinaryHeap head");
    }

    #[test]
    fn seq_kbest_heap_equal_score_prefers_fewer_edges() {
        let few = HeapItem {
            score: -0.5,
            tie: 2,
            path: PathState::new(0),
        };
        let many = HeapItem {
            score: -0.5,
            tie: 3,
            path: PathState::new(0),
        };
        assert_eq!(cmp_graph_kbest_score(-0.5, -0.5), Ordering::Equal);
        assert!(
            few > many,
            "existing SeqGraph heap tie: fewer edges is Greater (polled first)"
        );
    }
}
