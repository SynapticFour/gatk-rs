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

    sort_seq_kbest_paths(graph, &mut result);
    Ok(result)
}

/// Production Peak-RSS caps on SeqGraph k-best (absent from Java
/// `GraphBasedKBestHaplotypeFinder`). Forensic comparison only — production search
/// is unchanged.
pub const SEQ_KBEST_PRODUCTION_MAX_HEAP: usize = 1_024;
pub const SEQ_KBEST_PRODUCTION_MAX_EXPANSIONS: usize = 12_000;
pub const SEQ_KBEST_PRODUCTION_MAX_PATH_EDGES: usize = 4_096;

/// Optional Peak-RSS bounds. `None` means unbounded (Java `PriorityQueue` has no cap).
#[derive(Debug, Clone, Copy)]
pub struct SeqKbestCapPolicy {
    pub max_heap_paths: Option<usize>,
    pub max_expansions: Option<usize>,
    pub max_path_edges: Option<usize>,
}

impl SeqKbestCapPolicy {
    pub fn production() -> Self {
        Self {
            max_heap_paths: Some(SEQ_KBEST_PRODUCTION_MAX_HEAP),
            max_expansions: Some(SEQ_KBEST_PRODUCTION_MAX_EXPANSIONS),
            max_path_edges: Some(SEQ_KBEST_PRODUCTION_MAX_PATH_EDGES),
        }
    }

    pub fn unbounded() -> Self {
        Self {
            max_heap_paths: None,
            max_expansions: None,
            max_path_edges: None,
        }
    }
}

/// One edge's contribution to a k-best path score.
#[derive(Debug, Clone)]
pub struct SeqKbestEdgeTerm {
    pub from: usize,
    pub to: usize,
    pub edge_support: u32,
    pub total_outgoing: u32,
    pub penalty: f64,
    pub is_ref: bool,
}

/// First time a needle sequence was collected as a completed sink path.
#[derive(Debug, Clone)]
pub struct SeqKbestNeedleHit {
    pub sink_ordinal: usize,
    pub score: f64,
    pub n_edges: usize,
    pub rank_after_sort: Option<usize>,
}

/// Forensic trace of one SeqGraph k-best run. Does not alter production search.
#[derive(Debug, Clone)]
pub struct SeqKbestForensicReport {
    pub paths: Vec<KBestPath>,
    pub expansions: usize,
    pub max_heap: usize,
    pub pop_count: usize,
    pub skip_heap_full_at_pop: usize,
    pub skip_heap_full_at_expand: usize,
    pub skip_expansion_cap: usize,
    pub skip_path_edge_cap: usize,
    pub vertex_visit_refused: usize,
    pub heap_remaining: usize,
    pub needles_in_result: Vec<Option<SeqKbestNeedleHit>>,
    pub needles_on_remaining_heap: Vec<bool>,
}

fn contains_bases(hay: &[u8], needle: &[u8]) -> bool {
    needle.is_empty() || hay.windows(needle.len()).any(|w| w == needle)
}

fn sort_seq_kbest_paths(graph: &SeqGraph, result: &mut [KBestPath]) {
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
}

/// Reconstruct GATK `KBestHaplotype.score` as the sum of
/// `log10(edgeMult) - log10(totalOutgoing)` along the path.
pub fn seq_kbest_path_score_terms(graph: &SeqGraph, path: &KBestPath) -> Vec<SeqKbestEdgeTerm> {
    let mut terms = Vec::with_capacity(path.edges.len());
    for &(from, to) in &path.edges {
        let outs = graph.outgoing_nodes(from);
        let total: u32 = outs
            .iter()
            .filter_map(|&t| graph.edge_support(from, t))
            .sum();
        let support = graph.edge_support(from, to).unwrap_or(0);
        terms.push(SeqKbestEdgeTerm {
            from,
            to,
            edge_support: support,
            total_outgoing: total,
            penalty: log_penalty(support, total),
            is_ref: graph.edge_is_ref(from, to),
        });
    }
    terms
}

/// Java `GraphBasedKBestHaplotypeFinder` score order on two finite path scores.
pub fn seq_kbest_score_cmp(lhs: f64, rhs: f64) -> Ordering {
    cmp_graph_kbest_score(lhs, rhs)
}

/// SeqGraph k-best with explicit result/visit limits and Peak-RSS cap policy.
///
/// Production is `max_results == vertex_visit_limit == K` plus [`SeqKbestCapPolicy::production`].
/// Java ties result count and per-vertex visit budget to the same `maxNumberOfHaplotypes`
/// and has no heap/expansion caps.
pub fn find_best_haplotypes_seq_graph_forensic(
    graph: &SeqGraph,
    max_results: usize,
    vertex_visit_limit: usize,
    caps: SeqKbestCapPolicy,
    needles: &[&[u8]],
) -> GatkResult<SeqKbestForensicReport> {
    let mut needles_in_result: Vec<Option<SeqKbestNeedleHit>> = vec![None; needles.len()];
    let mut needles_on_remaining_heap = vec![false; needles.len()];
    if max_results == 0 {
        return Ok(SeqKbestForensicReport {
            paths: Vec::new(),
            expansions: 0,
            max_heap: 0,
            pop_count: 0,
            skip_heap_full_at_pop: 0,
            skip_heap_full_at_expand: 0,
            skip_expansion_cap: 0,
            skip_path_edge_cap: 0,
            vertex_visit_refused: 0,
            heap_remaining: 0,
            needles_in_result,
            needles_on_remaining_heap,
        });
    }
    let source = graph
        .reference_source_vertex()
        .ok_or_else(|| GatkError::algorithm("seq kbest forensic: no reference source vertex"))?;
    let sink = graph
        .reference_sink_vertex()
        .ok_or_else(|| GatkError::algorithm("seq kbest forensic: no reference sink vertex"))?;
    let sinks: HashSet<usize> = HashSet::from([sink]);

    let mut result = Vec::new();
    let mut heap: BinaryHeap<HeapItem> = BinaryHeap::new();
    heap.push(HeapItem {
        score: 0.0,
        tie: 0,
        path: PathState::new(source),
    });

    let mut vertex_counts = vec![0usize; graph.node_count()];
    let mut expansions = 0usize;
    let mut max_heap = heap.len();
    let mut pop_count = 0usize;
    let mut skip_heap_full_at_pop = 0usize;
    let mut skip_heap_full_at_expand = 0usize;
    let mut skip_expansion_cap = 0usize;
    let mut skip_path_edge_cap = 0usize;
    let mut vertex_visit_refused = 0usize;

    while !heap.is_empty() && result.len() < max_results {
        max_heap = max_heap.max(heap.len());
        pop_count += 1;
        let item = heap.pop().expect("non-empty");
        let path = item.path;
        if sinks.contains(&path.last) {
            let kp = KBestPath {
                start: path.start,
                edges: path.edges,
                score: path.score,
                is_reference: path.is_reference,
            };
            let bases = graph.path_bases_bytes(kp.start, &kp.edges);
            for (i, needle) in needles.iter().enumerate() {
                if needles_in_result[i].is_none() && contains_bases(&bases, needle) {
                    needles_in_result[i] = Some(SeqKbestNeedleHit {
                        sink_ordinal: result.len(),
                        score: kp.score,
                        n_edges: kp.edges.len(),
                        rank_after_sort: None,
                    });
                }
            }
            result.push(kp);
            continue;
        }
        let heap_full = caps.max_heap_paths.is_some_and(|m| heap.len() >= m);
        let exp_full = caps.max_expansions.is_some_and(|m| expansions >= m);
        let path_full = caps.max_path_edges.is_some_and(|m| path.edges.len() >= m);
        if path_full {
            skip_path_edge_cap += 1;
            continue;
        }
        if exp_full {
            skip_expansion_cap += 1;
            continue;
        }
        if heap_full {
            skip_heap_full_at_pop += 1;
            continue;
        }
        if vertex_counts[path.last] < vertex_visit_limit {
            vertex_counts[path.last] += 1;
            expansions += 1;
            let outs = graph.outgoing_nodes(path.last);
            let total: u32 = outs
                .iter()
                .filter_map(|&t| graph.edge_support(path.last, t))
                .sum();
            for to in outs {
                if let Some(support) = graph.edge_support(path.last, to) {
                    if caps.max_heap_paths.is_some_and(|m| heap.len() >= m) {
                        skip_heap_full_at_expand += 1;
                        break;
                    }
                    let extended = path.extend(graph, to, support, total);
                    heap.push(HeapItem {
                        score: extended.score,
                        tie: extended.edge_count,
                        path: extended,
                    });
                    max_heap = max_heap.max(heap.len());
                }
            }
        } else {
            vertex_visit_refused += 1;
        }
    }

    let heap_remaining = heap.len();
    for item in heap {
        let bases = graph.path_bases_bytes(item.path.start, &item.path.edges);
        for (i, needle) in needles.iter().enumerate() {
            if contains_bases(&bases, needle) {
                needles_on_remaining_heap[i] = true;
            }
        }
    }

    sort_seq_kbest_paths(graph, &mut result);
    for (i, needle) in needles.iter().enumerate() {
        if let Some(hit) = needles_in_result[i].as_mut() {
            hit.rank_after_sort = result
                .iter()
                .position(|p| contains_bases(&graph.path_bases_bytes(p.start, &p.edges), needle));
        }
    }

    Ok(SeqKbestForensicReport {
        paths: result,
        expansions,
        max_heap,
        pop_count,
        skip_heap_full_at_pop,
        skip_heap_full_at_expand,
        skip_expansion_cap,
        skip_path_edge_cap,
        vertex_visit_refused,
        heap_remaining,
        needles_in_result,
        needles_on_remaining_heap,
    })
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

    #[test]
    fn seq_kbest_score_is_sum_of_log10_multiplicity_ratios() {
        let seq = diamond_40_vs_5();
        let paths = find_best_haplotypes_seq_graph(&seq, 2).expect("kbest");
        for p in &paths {
            let terms = seq_kbest_path_score_terms(&seq, p);
            let sum: f64 = terms.iter().map(|t| t.penalty).sum();
            assert!(
                (sum - p.score).abs() < 1e-12,
                "score must be sum of edge log-penalties"
            );
            assert!(
                terms.iter().all(|t| t.penalty <= 0.0),
                "each log10(mult/out) term is ≤ 0"
            );
            for t in &terms {
                let expect = log_penalty(t.edge_support, t.total_outgoing);
                assert!((t.penalty - expect).abs() < 1e-15);
            }
        }
        let high = &paths[0];
        let terms = seq_kbest_path_score_terms(&seq, high);
        assert_eq!(high.edges.len(), 2);
        let expected = log_penalty(40, 45) + log_penalty(40, 40);
        assert!(
            (high.score - expected).abs() < 1e-12,
            "high-support arm: branch then unique continuation (penalty 0); got {} expected {}",
            high.score,
            expected
        );
        assert_eq!(
            terms
                .iter()
                .filter(|t| t.total_outgoing != t.edge_support)
                .count(),
            1
        );
    }

    #[test]
    fn seq_kbest_forensic_production_caps_match_production_search() {
        let seq = diamond_40_vs_5();
        let prod = find_best_haplotypes_seq_graph(&seq, 2).expect("prod");
        let forensic = find_best_haplotypes_seq_graph_forensic(
            &seq,
            2,
            2,
            SeqKbestCapPolicy::production(),
            &[],
        )
        .expect("forensic");
        assert_eq!(forensic.paths.len(), prod.len());
        for (a, b) in forensic.paths.iter().zip(prod.iter()) {
            assert!((a.score - b.score).abs() < 1e-15);
            assert_eq!(a.edges, b.edges);
        }
        assert_eq!(forensic.skip_heap_full_at_pop, 0);
        assert_eq!(forensic.skip_expansion_cap, 0);
    }

    #[test]
    fn seq_kbest_k1_cutoff_omits_walkable_lower_score_arm() {
        let seq = diamond_40_vs_5();
        let k1 = find_best_haplotypes_seq_graph_forensic(
            &seq,
            1,
            1,
            SeqKbestCapPolicy::production(),
            &[b"AAG"],
        )
        .expect("k1");
        let k2 = find_best_haplotypes_seq_graph_forensic(
            &seq,
            2,
            2,
            SeqKbestCapPolicy::production(),
            &[b"AAG"],
        )
        .expect("k2");
        assert!(k1.needles_in_result[0].is_none());
        assert!(k2.needles_in_result[0].is_some());
        assert_eq!(
            k2.needles_in_result[0].as_ref().unwrap().rank_after_sort,
            Some(1)
        );
    }
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
