//! Forensic compact / lazy SeqGraph k-best (6R.136).
//!
//! Does not replace [`crate::seq_kbest_haplotype::find_best_haplotypes_seq_graph`].
//! Ancestry is a parent-pointer arena; sequences are reconstructed only for completed
//! paths. Eager mode still puts every generated successor on the `BinaryHeap`.
//! Lazy mode keeps sibling tails off the heap until the previous sibling is popped.

use crate::kbest_haplotype::{cmp_graph_kbest_score, log_penalty, KBestPath};
use crate::seq_graph::SeqGraph;
use crate::seq_kbest_haplotype::{sort_seq_kbest_paths, SeqKbestFrontierMem};
use gatk_common::{GatkError, GatkResult};
use std::cmp::Ordering;
use std::collections::BinaryHeap;

#[derive(Debug, Clone, Copy)]
struct CompactNode {
    parent: Option<u32>,
    start: usize,
    last: usize,
    from: usize,
    to: usize,
    score: f64,
    edge_count: usize,
    is_reference: bool,
    next_sib: Option<u32>,
}

struct CompactHeapItem {
    score: f64,
    tie: usize,
    idx: u32,
}

impl Ord for CompactHeapItem {
    fn cmp(&self, other: &Self) -> Ordering {
        cmp_graph_kbest_score(self.score, other.score).then_with(|| other.tie.cmp(&self.tie))
    }
}

impl PartialOrd for CompactHeapItem {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl PartialEq for CompactHeapItem {
    fn eq(&self, other: &Self) -> bool {
        cmp_graph_kbest_score(self.score, other.score).is_eq() && self.tie == other.tie
    }
}

impl Eq for CompactHeapItem {}

/// Forensic metrics for a compact-ancestry k-best run.
#[derive(Debug, Clone)]
pub struct SeqKbestCompactReport {
    pub paths: Vec<KBestPath>,
    pub mode: &'static str,
    pub peak_heap: usize,
    pub peak_live_unpopped: usize,
    pub peak_deferred_sibs: usize,
    pub arena_nodes: usize,
    pub reconstructions: usize,
    pub pops: usize,
    pub expansions: usize,
    pub skip_heap_full: usize,
    /// `live_unpopped + deferred_sibs` peak (logical generated-not-yet-popped).
    pub peak_logical_unpopped: usize,
    /// Arena + heap-vector accounting. Production search is unchanged.
    pub mem: SeqKbestFrontierMem,
}

fn reconstruct_edges(arena: &[CompactNode], mut idx: u32) -> (usize, Vec<(usize, usize)>, bool) {
    let start = arena[idx as usize].start;
    let is_ref = arena[idx as usize].is_reference;
    let mut rev = Vec::new();
    loop {
        let n = &arena[idx as usize];
        if n.edge_count == 0 {
            break;
        }
        rev.push((n.from, n.to));
        match n.parent {
            Some(p) => idx = p,
            None => break,
        }
    }
    rev.reverse();
    (start, rev, is_ref)
}

fn push_seed(arena: &mut Vec<CompactNode>, heap: &mut BinaryHeap<CompactHeapItem>, source: usize) {
    let idx = arena.len() as u32;
    arena.push(CompactNode {
        parent: None,
        start: source,
        last: source,
        from: source,
        to: source,
        score: 0.0,
        edge_count: 0,
        is_reference: false,
        next_sib: None,
    });
    heap.push(CompactHeapItem {
        score: 0.0,
        tie: 0,
        idx,
    });
}

fn make_child(
    arena: &[CompactNode],
    parent_idx: u32,
    graph: &SeqGraph,
    to: usize,
    support: u32,
    total: u32,
) -> CompactNode {
    let p = arena[parent_idx as usize];
    let penalty = log_penalty(support, total);
    CompactNode {
        parent: Some(parent_idx),
        start: p.start,
        last: to,
        from: p.last,
        to,
        score: p.score + penalty,
        edge_count: p.edge_count + 1,
        is_reference: p.is_reference && graph.edge_is_ref(p.last, to),
        next_sib: None,
    }
}

/// Shared-ancestry k-best. `lazy_siblings`: only the best unused sibling sits on the
/// `BinaryHeap`; remaining siblings are a linked list in the arena (still live objects).
/// `max_heap_paths`: `None` is unbounded (Java). Cap applies to `BinaryHeap.len()` only.
pub fn find_best_haplotypes_seq_graph_compact(
    graph: &SeqGraph,
    max_results: usize,
    vertex_visit_limit: usize,
    lazy_siblings: bool,
    max_heap_paths: Option<usize>,
) -> GatkResult<SeqKbestCompactReport> {
    if max_results == 0 {
        return Ok(SeqKbestCompactReport {
            paths: Vec::new(),
            mode: if lazy_siblings {
                "compact_lazy"
            } else {
                "compact_eager"
            },
            peak_heap: 0,
            peak_live_unpopped: 0,
            peak_deferred_sibs: 0,
            arena_nodes: 0,
            reconstructions: 0,
            pops: 0,
            expansions: 0,
            skip_heap_full: 0,
            peak_logical_unpopped: 0,
            mem: SeqKbestFrontierMem::default(),
        });
    }
    let source = graph
        .reference_source_vertex()
        .ok_or_else(|| GatkError::algorithm("compact kbest: no source"))?;
    let sink = graph
        .reference_sink_vertex()
        .ok_or_else(|| GatkError::algorithm("compact kbest: no sink"))?;

    let mut arena: Vec<CompactNode> = Vec::new();
    let mut heap: BinaryHeap<CompactHeapItem> = BinaryHeap::new();
    push_seed(&mut arena, &mut heap, source);

    let mut vertex_counts = vec![0usize; graph.node_count()];
    let mut result = Vec::new();
    let mut expansions = 0usize;
    let mut pops = 0usize;
    let mut reconstructions = 0usize;
    let mut skip_heap_full = 0usize;
    let mut peak_heap = heap.len();
    let mut live_unpopped = 1usize;
    let mut peak_live = 1usize;
    let mut deferred_sibs = 0usize;
    let mut peak_deferred = 0usize;
    let mut peak_logical = 1usize;
    let mut mem = SeqKbestFrontierMem::default();
    mem.rss_before_mib = crate::runtime_config::current_rss_mib();
    sample_compact_rss(&mut mem);
    let mut result_payload_cap = 0usize;
    record_compact_peak(&mut mem, &arena, &heap, live_unpopped, deferred_sibs, 0, 0);

    while !heap.is_empty() && result.len() < max_results {
        if crate::runtime_config::hc_rss_abort_triggered() {
            mem.rss_aborted = true;
            break;
        }
        peak_heap = peak_heap.max(heap.len());
        pops += 1;
        if pops == 1 || pops % 32 == 0 {
            sample_compact_rss(&mut mem);
        }
        let Some(item) = heap.pop() else {
            break;
        };
        live_unpopped = live_unpopped.saturating_sub(1);
        let idx = item.idx;
        let node = arena[idx as usize];

        if lazy_siblings {
            if let Some(sib) = node.next_sib {
                deferred_sibs = deferred_sibs.saturating_sub(1);
                if max_heap_paths.is_some_and(|m| heap.len() >= m) {
                    skip_heap_full += 1;
                } else {
                    let s = arena[sib as usize];
                    heap.push(CompactHeapItem {
                        score: s.score,
                        tie: s.edge_count,
                        idx: sib,
                    });
                    live_unpopped += 1;
                }
            }
        }

        if node.last == sink {
            reconstructions += 1;
            let (start, edges, is_reference) = reconstruct_edges(&arena, idx);
            result.push(KBestPath {
                start,
                edges,
                score: node.score,
                is_reference,
            });
            result_payload_cap = result_payload_cap.saturating_add(
                result
                    .last()
                    .map(|p| p.edges.capacity() * std::mem::size_of::<(usize, usize)>())
                    .unwrap_or(0),
            );
            peak_live = peak_live.max(live_unpopped);
            peak_deferred = peak_deferred.max(deferred_sibs);
            peak_logical = peak_logical.max(live_unpopped + deferred_sibs);
            record_compact_peak(
                &mut mem,
                &arena,
                &heap,
                live_unpopped,
                deferred_sibs,
                result.len(),
                result_payload_cap,
            );
            continue;
        }

        if vertex_counts[node.last] < vertex_visit_limit {
            vertex_counts[node.last] += 1;
            expansions += 1;
            let outs = graph.outgoing_nodes(node.last);
            let total: u32 = outs
                .iter()
                .filter_map(|&t| graph.edge_support(node.last, t))
                .sum();
            let mut kids: Vec<CompactNode> = Vec::new();
            for to in outs {
                if let Some(support) = graph.edge_support(node.last, to) {
                    kids.push(make_child(&arena, idx, graph, to, support, total));
                }
            }
            if kids.is_empty() {
                peak_live = peak_live.max(live_unpopped);
                record_compact_peak(
                    &mut mem,
                    &arena,
                    &heap,
                    live_unpopped,
                    deferred_sibs,
                    result.len(),
                    result_payload_cap,
                );
                continue;
            }
            if lazy_siblings {
                kids.sort_by(|a, b| {
                    cmp_graph_kbest_score(b.score, a.score)
                        .then_with(|| a.edge_count.cmp(&b.edge_count))
                });
                let mut ids = Vec::with_capacity(kids.len());
                for k in kids {
                    ids.push(arena.len() as u32);
                    arena.push(k);
                }
                for i in 0..ids.len().saturating_sub(1) {
                    arena[ids[i] as usize].next_sib = Some(ids[i + 1]);
                }
                deferred_sibs += ids.len().saturating_sub(1);
                if max_heap_paths.is_some_and(|m| heap.len() >= m) {
                    skip_heap_full += 1;
                } else {
                    let head = ids[0];
                    let h = arena[head as usize];
                    heap.push(CompactHeapItem {
                        score: h.score,
                        tie: h.edge_count,
                        idx: head,
                    });
                    live_unpopped += 1;
                }
            } else {
                for k in kids {
                    if max_heap_paths.is_some_and(|m| heap.len() >= m) {
                        skip_heap_full += 1;
                        break;
                    }
                    let cidx = arena.len() as u32;
                    heap.push(CompactHeapItem {
                        score: k.score,
                        tie: k.edge_count,
                        idx: cidx,
                    });
                    arena.push(k);
                    live_unpopped += 1;
                }
            }
        }
        peak_heap = peak_heap.max(heap.len());
        peak_live = peak_live.max(live_unpopped);
        peak_deferred = peak_deferred.max(deferred_sibs);
        peak_logical = peak_logical.max(live_unpopped + deferred_sibs);
        record_compact_peak(
            &mut mem,
            &arena,
            &heap,
            live_unpopped,
            deferred_sibs,
            result.len(),
            result_payload_cap,
        );
    }

    mem.rss_after_mib = crate::runtime_config::current_rss_mib();
    sample_compact_rss(&mut mem);

    sort_seq_kbest_paths(graph, &mut result);
    Ok(SeqKbestCompactReport {
        paths: result,
        mode: if lazy_siblings {
            "compact_lazy"
        } else {
            "compact_eager"
        },
        peak_heap,
        peak_live_unpopped: peak_live,
        peak_deferred_sibs: peak_deferred,
        arena_nodes: arena.len(),
        reconstructions,
        pops,
        expansions,
        skip_heap_full,
        peak_logical_unpopped: peak_logical,
        mem,
    })
}

fn sample_compact_rss(mem: &mut SeqKbestFrontierMem) {
    let Some(rss) = crate::runtime_config::current_rss_mib() else {
        return;
    };
    mem.rss_samples = mem.rss_samples.saturating_add(1);
    mem.rss_peak_mib = Some(match mem.rss_peak_mib {
        Some(p) => p.max(rss),
        None => rss,
    });
}

fn record_compact_peak(
    mem: &mut SeqKbestFrontierMem,
    arena: &Vec<CompactNode>,
    heap: &BinaryHeap<CompactHeapItem>,
    live_unpopped: usize,
    deferred_sibs: usize,
    result_n: usize,
    result_payload_cap: usize,
) {
    let node_b = std::mem::size_of::<CompactNode>();
    let item_b = std::mem::size_of::<CompactHeapItem>();
    let arena_bytes = arena.capacity().saturating_mul(node_b);
    let heap_inline = heap.capacity().saturating_mul(item_b);
    let live_nodes = live_unpopped.saturating_add(deferred_sibs);
    let live_frontier = live_nodes
        .saturating_mul(node_b)
        .saturating_add(heap_inline);
    mem.peak_logical = mem.peak_logical.max(live_nodes);
    mem.peak_heap_len = mem.peak_heap_len.max(heap.len());
    mem.peak_heap_capacity = mem.peak_heap_capacity.max(heap.capacity());
    mem.peak_inline_bytes = mem
        .peak_inline_bytes
        .max(arena_bytes.saturating_add(heap_inline));
    mem.peak_payload_len_bytes = 0;
    mem.peak_payload_cap_bytes = 0;
    mem.peak_frontier_bytes = mem.peak_frontier_bytes.max(live_frontier);
    mem.peak_sum_edge_len = mem.peak_sum_edge_len.max(arena.len());
    mem.peak_sum_edge_cap = mem.peak_sum_edge_cap.max(arena.capacity());
    mem.peak_result_n = mem.peak_result_n.max(result_n);
    mem.peak_result_payload_cap_bytes = mem.peak_result_payload_cap_bytes.max(result_payload_cap);
    mem.peak_frontier_plus_results_bytes = mem.peak_frontier_plus_results_bytes.max(
        arena_bytes
            .saturating_add(heap_inline)
            .saturating_add(result_payload_cap),
    );
}

/// Bytes of one production-style copied edge list of `n_edges` pairs, excluding `Vec` header.
pub fn copied_edge_payload_bytes(n_edges: usize) -> usize {
    n_edges.saturating_mul(std::mem::size_of::<(usize, usize)>())
}

/// `size_of` the forensic compact node (parent pointer + score + vertex; no `Vec`).
pub fn compact_node_size() -> usize {
    std::mem::size_of::<CompactNode>()
}

pub fn compact_heap_item_size() -> usize {
    std::mem::size_of::<CompactHeapItem>()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::assembly::AssemblyGraph;
    use crate::seq_kbest_haplotype::find_best_haplotypes_seq_graph;

    fn diamond_40_vs_5() -> SeqGraph {
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
        SeqGraph::from_assembly_graph(&g)
    }

    #[test]
    fn compact_eager_and_lazy_match_production_diamond() {
        let seq = diamond_40_vs_5();
        let prod = find_best_haplotypes_seq_graph(&seq, 2).expect("prod");
        let eager = find_best_haplotypes_seq_graph_compact(&seq, 2, 2, false, None).expect("eager");
        let lazy = find_best_haplotypes_seq_graph_compact(&seq, 2, 2, true, None).expect("lazy");
        assert_eq!(eager.paths.len(), prod.len());
        assert_eq!(lazy.paths.len(), prod.len());
        for i in 0..prod.len() {
            assert!((eager.paths[i].score - prod[i].score).abs() < 1e-15);
            assert_eq!(eager.paths[i].edges, prod[i].edges);
            assert!((lazy.paths[i].score - prod[i].score).abs() < 1e-15);
            assert_eq!(lazy.paths[i].edges, prod[i].edges);
        }
    }

    fn four_way_bubble() -> SeqGraph {
        let mut g = AssemblyGraph::new(3).unwrap();
        let src = g.ensure_node(b"AAA");
        let m0 = g.ensure_node(b"AAC");
        let m1 = g.ensure_node(b"AAG");
        let m2 = g.ensure_node(b"AAT");
        let m3 = g.ensure_node(b"ACA");
        let snk = g.ensure_node(b"ACT");
        for (mid, sup) in [(m0, 40u32), (m1, 30), (m2, 20), (m3, 10)] {
            g.add_edge_support(src, mid, sup);
            g.add_edge_support(mid, snk, sup);
        }
        g.ref_edges.insert((src, m0));
        g.ref_edges.insert((m0, snk));
        g.ref_nodes.insert(src);
        g.ref_nodes.insert(m0);
        g.ref_nodes.insert(snk);
        g.ref_source_kmer = Some(std::sync::Arc::from(b"AAA".as_slice()));
        SeqGraph::from_assembly_graph(&g)
    }

    #[test]
    fn lazy_peak_logical_can_exceed_peak_heap_entries() {
        let seq = four_way_bubble();
        let lazy = find_best_haplotypes_seq_graph_compact(&seq, 4, 4, true, None).expect("lazy");
        let eager = find_best_haplotypes_seq_graph_compact(&seq, 4, 4, false, None).expect("eager");
        assert_eq!(
            eager.peak_logical_unpopped, eager.peak_heap,
            "copied/eager compact: each heap entry is one logical partial"
        );
        assert!(
            lazy.peak_logical_unpopped >= lazy.peak_heap,
            "lazy logical must count deferred siblings"
        );
        assert!(
            lazy.peak_logical_unpopped > lazy.peak_heap,
            "four-way bubble: deferred siblings make logical > heap heads (logical={} heap={})",
            lazy.peak_logical_unpopped,
            lazy.peak_heap
        );
        assert_eq!(
            lazy.peak_logical_unpopped, eager.peak_logical_unpopped,
            "logical occupancy is independent of heap-object representation"
        );
        assert_ne!(
            lazy.mem.peak_frontier_bytes as usize, lazy.peak_logical_unpopped,
            "byte accounting must not silently equal logical-state count"
        );
    }
}
