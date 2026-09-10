//! GATK `GraphBasedKBestHaplotypeFinder` on [`SeqGraph`] (`Path.getBases` stitching).

use crate::kbest_haplotype::{cmp_graph_kbest_score, log_penalty, KBestPath};
use crate::seq_graph::SeqGraph;
use crate::seq_kbest_resource_policy::{
    SeqKbestResourcePolicy, SeqKbestRunDiagnostics, LEGACY_MAX_HEAP_PATHS,
};
use gatk_common::{GatkError, GatkResult};
use std::cmp::Ordering;
use std::collections::{BinaryHeap, HashSet};

/// Same bound as the read-threading k-best finder — see `kbest_haplotype`.
/// Resource-policy default (`legacy_1024`) uses this value; do not change it here.
const MAX_KBEST_HEAP_PATHS: usize = LEGACY_MAX_HEAP_PATHS;
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

const EDGE_PAIR_BYTES: usize = std::mem::size_of::<(usize, usize)>();

/// Forensic: if `dest` already has `max_per_dest` heap entries, drop the worst
/// when `new_score` is strictly better. Returns whether the new path may be inserted.
fn heap_apply_per_dest_cap(
    heap: &mut BinaryHeap<HeapItem>,
    dest: usize,
    new_score: f64,
    max_per_dest: usize,
) -> bool {
    let n_at_dest = heap.iter().filter(|h| h.path.last == dest).count();
    if n_at_dest < max_per_dest {
        return true;
    }
    let worst = heap
        .iter()
        .filter(|h| h.path.last == dest)
        .map(|h| h.score)
        .min_by(|a, b| a.total_cmp(b));
    let Some(worst) = worst else {
        return true;
    };
    if cmp_graph_kbest_score(new_score, worst) != Ordering::Greater {
        return false;
    }
    let mut items: Vec<HeapItem> = heap.drain().collect();
    if let Some(idx) = items
        .iter()
        .position(|h| h.path.last == dest && h.score == worst)
    {
        items.swap_remove(idx);
    }
    for it in items {
        heap.push(it);
    }
    true
}

/// Instrumented SeqGraph k-best result (paths + resource diagnostics).
#[derive(Debug, Clone)]
pub struct SeqKbestPolicyReport {
    pub paths: Vec<KBestPath>,
    pub diagnostics: SeqKbestRunDiagnostics,
}

#[derive(Clone, Copy)]
struct CopiedFrontierSnap {
    heap_entries: usize,
    frontier_bytes: usize,
    pathstate_bytes: usize,
    payload_cap_bytes: usize,
    payload_len_bytes: usize,
    heap_vec_capacity_bytes: usize,
}

fn snap_copied_frontier(heap: &BinaryHeap<HeapItem>) -> CopiedFrontierSnap {
    let mut payload_cap = 0usize;
    let mut payload_len = 0usize;
    for it in heap.iter() {
        payload_cap =
            payload_cap.saturating_add(it.path.edges.capacity().saturating_mul(EDGE_PAIR_BYTES));
        payload_len =
            payload_len.saturating_add(it.path.edges.len().saturating_mul(EDGE_PAIR_BYTES));
    }
    let n = heap.len();
    let item_b = std::mem::size_of::<HeapItem>();
    CopiedFrontierSnap {
        heap_entries: n,
        frontier_bytes: n.saturating_mul(item_b).saturating_add(payload_cap),
        pathstate_bytes: n.saturating_mul(std::mem::size_of::<PathState>()),
        payload_cap_bytes: payload_cap,
        payload_len_bytes: payload_len,
        heap_vec_capacity_bytes: heap.capacity().saturating_mul(item_b),
    }
}

fn extra_path_frontier_bytes(path: &PathState) -> usize {
    std::mem::size_of::<HeapItem>()
        .saturating_add(path.edges.capacity().saturating_mul(EDGE_PAIR_BYTES))
}

/// Frontier resource check. Path-edge / expansion caps are separate Peak-RSS guards.
fn blocks_expand(policy: SeqKbestResourcePolicy, snap: CopiedFrontierSnap) -> bool {
    match policy {
        SeqKbestResourcePolicy::Legacy1024 => snap.heap_entries >= MAX_KBEST_HEAP_PATHS,
        SeqKbestResourcePolicy::UnboundedDiagnostic => false,
        SeqKbestResourcePolicy::ByteBudget { max_frontier_bytes } => {
            snap.frontier_bytes >= max_frontier_bytes
        }
    }
}

fn blocks_push(
    policy: SeqKbestResourcePolicy,
    snap: CopiedFrontierSnap,
    extra: &PathState,
) -> bool {
    match policy {
        SeqKbestResourcePolicy::Legacy1024 => snap.heap_entries >= MAX_KBEST_HEAP_PATHS,
        SeqKbestResourcePolicy::UnboundedDiagnostic => false,
        SeqKbestResourcePolicy::ByteBudget { max_frontier_bytes } => {
            snap.frontier_bytes
                .saturating_add(extra_path_frontier_bytes(extra))
                >= max_frontier_bytes
        }
    }
}

fn sample_rss_peak(peak: &mut Option<f64>) {
    if let Some(rss) = crate::runtime_config::current_rss_mib() {
        *peak = Some(match *peak {
            Some(p) => p.max(rss),
            None => rss,
        });
    }
}

fn finish_policy_report(mut report: SeqKbestPolicyReport) -> SeqKbestPolicyReport {
    report.diagnostics.completed_paths = report.paths.len();
    report
}

/// K-best paths on a cleaned sequence graph (GATK `GraphBasedKBestHaplotypeFinder`).
///
/// Default resource policy is [`SeqKbestResourcePolicy::Legacy1024`] (unset env).
/// Experimental env: `GATK_RS_EXPERIMENTAL_KBEST_POLICY` — not a product contract.
pub fn find_best_haplotypes_seq_graph(
    graph: &SeqGraph,
    max_number_of_haplotypes: usize,
) -> GatkResult<Vec<KBestPath>> {
    let policy = SeqKbestResourcePolicy::from_runtime()?;
    let report =
        find_best_haplotypes_seq_graph_with_policy(graph, max_number_of_haplotypes, policy)?;
    if crate::runtime_config::kbest_diagnostics_enabled() {
        eprintln!("{}", report.diagnostics.to_log_line());
    }
    Ok(report.paths)
}

/// SeqGraph k-best with an explicit resource policy (6R.140). Production callers
/// should use [`find_best_haplotypes_seq_graph`] so the default stays `legacy_1024`.
pub fn find_best_haplotypes_seq_graph_with_policy(
    graph: &SeqGraph,
    max_number_of_haplotypes: usize,
    policy: SeqKbestResourcePolicy,
) -> GatkResult<SeqKbestPolicyReport> {
    let rss_before = crate::runtime_config::current_rss_mib();
    let mut rss_peak = rss_before;
    if max_number_of_haplotypes == 0 {
        return Ok(finish_policy_report(SeqKbestPolicyReport {
            paths: Vec::new(),
            diagnostics: SeqKbestRunDiagnostics {
                requested_k: 0,
                completed_paths: 0,
                policy: policy.name(),
                resource_limit: policy.resource_limit_display(),
                resource_limit_hit: false,
                peak_logical_frontier: 0,
                peak_heap_entries: 0,
                peak_frontier_bytes: 0,
                peak_pathstate_bytes: 0,
                peak_edge_payload_bytes: 0,
                peak_edge_payload_len_bytes: 0,
                peak_heap_vec_capacity_bytes: 0,
                paths_refused: 0,
                expansions: 0,
                rss_before_mib: rss_before,
                rss_peak_mib: rss_peak,
                rss_after_mib: crate::runtime_config::current_rss_mib(),
                accounting: "copied_pathstate_experimental",
            },
        }));
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
    // Not part of the 6R.140 frontier resource policy.
    const MAX_KBEST_EXPANSIONS: usize = 12_000;
    let mut expansions = 0usize;
    let mut paths_refused = 0usize;
    let mut resource_limit_hit = false;
    let mut peak_logical = heap.len();
    let mut peak_heap = heap.len();
    let mut peak_frontier_bytes = 0usize;
    let mut peak_pathstate_bytes = 0usize;
    let mut peak_payload_cap = 0usize;
    let mut peak_payload_len = 0usize;
    let mut peak_heap_vec_cap = 0usize;

    let mut record = |heap: &BinaryHeap<HeapItem>| {
        let snap = snap_copied_frontier(heap);
        peak_logical = peak_logical.max(snap.heap_entries);
        peak_heap = peak_heap.max(snap.heap_entries);
        peak_frontier_bytes = peak_frontier_bytes.max(snap.frontier_bytes);
        peak_pathstate_bytes = peak_pathstate_bytes.max(snap.pathstate_bytes);
        peak_payload_cap = peak_payload_cap.max(snap.payload_cap_bytes);
        peak_payload_len = peak_payload_len.max(snap.payload_len_bytes);
        peak_heap_vec_cap = peak_heap_vec_cap.max(snap.heap_vec_capacity_bytes);
        snap
    };
    record(&heap);

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
            record(&heap);
            continue;
        }
        let snap = record(&heap);
        if path.edges.len() >= MAX_KBEST_PATH_EDGES || expansions >= MAX_KBEST_EXPANSIONS {
            continue;
        }
        if blocks_expand(policy, snap) {
            resource_limit_hit = true;
            paths_refused += 1;
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
                    let snap = snap_copied_frontier(&heap);
                    let extended = path.extend(graph, to, support, total);
                    if blocks_push(policy, snap, &extended) {
                        resource_limit_hit = true;
                        paths_refused += 1;
                        break;
                    }
                    heap.push(HeapItem {
                        score: extended.score,
                        tie: extended.edge_count,
                        path: extended,
                    });
                    record(&heap);
                }
            }
            if expansions % 32 == 0 {
                sample_rss_peak(&mut rss_peak);
            }
        }
    }

    sort_seq_kbest_paths(graph, &mut result);
    let rss_after = crate::runtime_config::current_rss_mib();
    sample_rss_peak(&mut rss_peak);
    Ok(finish_policy_report(SeqKbestPolicyReport {
        paths: result,
        diagnostics: SeqKbestRunDiagnostics {
            requested_k: max_number_of_haplotypes,
            completed_paths: 0,
            policy: policy.name(),
            resource_limit: policy.resource_limit_display(),
            resource_limit_hit,
            peak_logical_frontier: peak_logical,
            peak_heap_entries: peak_heap,
            peak_frontier_bytes,
            peak_pathstate_bytes,
            peak_edge_payload_bytes: peak_payload_cap,
            peak_edge_payload_len_bytes: peak_payload_len,
            peak_heap_vec_capacity_bytes: peak_heap_vec_cap,
            paths_refused,
            expansions,
            rss_before_mib: rss_before,
            rss_peak_mib: rss_peak,
            rss_after_mib: rss_after,
            accounting: "copied_pathstate_experimental",
        },
    }))
}

/// Production Peak-RSS caps on SeqGraph k-best (absent from Java
/// `GraphBasedKBestHaplotypeFinder`). Forensic comparison only — production search
/// is unchanged.
pub const SEQ_KBEST_PRODUCTION_MAX_HEAP: usize = 1_024;
pub const SEQ_KBEST_PRODUCTION_MAX_EXPANSIONS: usize = 12_000;
pub const SEQ_KBEST_PRODUCTION_MAX_PATH_EDGES: usize = 4_096;

/// Layout of production k-best live objects. Forensic accounting only.
#[derive(Debug, Clone, Copy)]
pub struct SeqKbestProductionLayout {
    pub path_state_bytes: usize,
    pub heap_item_bytes: usize,
    pub edge_pair_bytes: usize,
}

pub fn seq_kbest_production_layout() -> SeqKbestProductionLayout {
    SeqKbestProductionLayout {
        path_state_bytes: std::mem::size_of::<PathState>(),
        heap_item_bytes: std::mem::size_of::<HeapItem>(),
        edge_pair_bytes: std::mem::size_of::<(usize, usize)>(),
    }
}

/// Optional Peak-RSS bounds. `None` means unbounded (Java `PriorityQueue` has no cap).
#[derive(Debug, Clone, Copy)]
pub struct SeqKbestCapPolicy {
    pub max_heap_paths: Option<usize>,
    pub max_expansions: Option<usize>,
    pub max_path_edges: Option<usize>,
    /// Forensic only. If `Some(n)`, keep at most `n` heap entries ending at the
    /// same vertex (retain higher scores). `None` is production insert semantics.
    pub max_heap_paths_per_dest: Option<usize>,
}

impl SeqKbestCapPolicy {
    pub fn production() -> Self {
        Self {
            max_heap_paths: Some(SEQ_KBEST_PRODUCTION_MAX_HEAP),
            max_expansions: Some(SEQ_KBEST_PRODUCTION_MAX_EXPANSIONS),
            max_path_edges: Some(SEQ_KBEST_PRODUCTION_MAX_PATH_EDGES),
            max_heap_paths_per_dest: None,
        }
    }

    pub fn unbounded() -> Self {
        Self {
            max_heap_paths: None,
            max_expansions: None,
            max_path_edges: None,
            max_heap_paths_per_dest: None,
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

/// First successor insertion refused because `max_heap_paths` was already reached.
/// Forensic only — production search is unchanged.
#[derive(Debug, Clone)]
pub struct SeqKbestExpandRefuse {
    pub pop_count: usize,
    pub expansions: usize,
    pub result_len: usize,
    pub heap_len_at_refuse: usize,
    pub parent_score: f64,
    pub parent_n_edges: usize,
    pub parent_last: usize,
    pub parent_start: usize,
    pub parent_edges: Vec<(usize, usize)>,
    pub n_outs: usize,
    pub n_inserted_this_expand: usize,
    pub inserted_tos: Vec<(usize, f64, u32)>,
    pub refused_tos: Vec<(usize, f64, u32)>,
}

/// Compact heap-entry view at the first expand-refuse. Forensic only.
#[derive(Debug, Clone)]
pub struct SeqKbestHeapEntryView {
    pub last: usize,
    pub start: usize,
    pub n_edges: usize,
    pub score: f64,
    pub is_sink: bool,
    pub edges_fnv: u64,
    pub bases_fnv: u64,
    pub prefix8_fnv: u64,
    pub prefix16_fnv: u64,
}

/// Heap + visit counts at the first `max_heap_paths` successor refuse.
#[derive(Debug, Clone)]
pub struct SeqKbestFirstRefuseSnapshot {
    pub vertex_counts: Vec<usize>,
    pub heap: Vec<SeqKbestHeapEntryView>,
    pub inserted_tos: Vec<(usize, f64, u32)>,
}

/// Forensic allocator/RSS accounting for one SeqGraph k-best run.
///
/// Payload bytes use `Vec::capacity` (not `len`) × 16 B per edge pair.
/// Inline bytes use `BinaryHeap` capacity × `size_of::<HeapItem>()` (includes
/// `PathState` / `Vec` headers). Allocator rounding is **not** included.
/// RSS fields are process RSS, not object-precise; treat as deltas.
#[derive(Debug, Clone, Default)]
pub struct SeqKbestFrontierMem {
    pub peak_logical: usize,
    pub peak_heap_len: usize,
    pub peak_heap_capacity: usize,
    pub peak_edge_len: usize,
    pub peak_edge_cap: usize,
    pub peak_sum_edge_len: usize,
    pub peak_sum_edge_cap: usize,
    pub peak_inline_bytes: usize,
    pub peak_payload_len_bytes: usize,
    pub peak_payload_cap_bytes: usize,
    pub peak_frontier_bytes: usize,
    pub peak_result_n: usize,
    pub peak_result_payload_cap_bytes: usize,
    pub peak_frontier_plus_results_bytes: usize,
    pub rss_before_mib: Option<f64>,
    pub rss_peak_mib: Option<f64>,
    pub rss_after_mib: Option<f64>,
    pub rss_samples: u32,
    pub rss_aborted: bool,
}

/// Forensic trace of one SeqGraph k-best run. Does not alter production search.
#[derive(Debug, Clone)]
pub struct SeqKbestForensicReport {
    pub paths: Vec<KBestPath>,
    pub expansions: usize,
    pub max_heap: usize,
    /// Copied-`PathState` frontier accounting. Production search is unchanged.
    pub mem: SeqKbestFrontierMem,
    pub pop_count: usize,
    pub skip_heap_full_at_pop: usize,
    pub skip_heap_full_at_expand: usize,
    pub skip_expansion_cap: usize,
    pub skip_path_edge_cap: usize,
    pub vertex_visit_refused: usize,
    pub heap_remaining: usize,
    pub needles_in_result: Vec<Option<SeqKbestNeedleHit>>,
    pub needles_on_remaining_heap: Vec<bool>,
    pub first_expand_refuse: Option<SeqKbestExpandRefuse>,
    pub heap_first_hit_cap_pop: Option<usize>,
    pub expand_refuses: Vec<SeqKbestExpandRefuse>,
    pub first_refuse_snapshot: Option<SeqKbestFirstRefuseSnapshot>,
}

fn contains_bases(hay: &[u8], needle: &[u8]) -> bool {
    needle.is_empty() || hay.windows(needle.len()).any(|w| w == needle)
}

fn fnv1a64_bytes(bytes: &[u8]) -> u64 {
    let mut h = 0xcbf29ce484222325u64;
    for &b in bytes {
        h ^= u64::from(b);
        h = h.wrapping_mul(0x100000001b3);
    }
    h
}

fn fnv1a64_edges(edges: &[(usize, usize)]) -> u64 {
    let mut h = 0xcbf29ce484222325u64;
    for &(a, b) in edges {
        h ^= a as u64;
        h = h.wrapping_mul(0x100000001b3);
        h ^= b as u64;
        h = h.wrapping_mul(0x100000001b3);
    }
    h
}

pub(crate) fn sort_seq_kbest_paths(graph: &SeqGraph, result: &mut [KBestPath]) {
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

fn sample_frontier_rss(mem: &mut SeqKbestFrontierMem) {
    let Some(rss) = crate::runtime_config::current_rss_mib() else {
        return;
    };
    mem.rss_samples = mem.rss_samples.saturating_add(1);
    mem.rss_peak_mib = Some(match mem.rss_peak_mib {
        Some(p) => p.max(rss),
        None => rss,
    });
}

fn scan_copied_frontier(
    mem: &mut SeqKbestFrontierMem,
    heap: &BinaryHeap<HeapItem>,
    result_n: usize,
    result_payload_cap: usize,
) {
    let mut live_len = 0usize;
    let mut live_cap = 0usize;
    for it in heap.iter() {
        live_len = live_len.saturating_add(it.path.edges.len().saturating_mul(EDGE_PAIR_BYTES));
        live_cap =
            live_cap.saturating_add(it.path.edges.capacity().saturating_mul(EDGE_PAIR_BYTES));
        mem.peak_edge_len = mem.peak_edge_len.max(it.path.edges.len());
        mem.peak_edge_cap = mem.peak_edge_cap.max(it.path.edges.capacity());
    }
    record_copied_peak(mem, heap, live_len, live_cap, result_n, result_payload_cap);
}

fn record_copied_peak(
    mem: &mut SeqKbestFrontierMem,
    heap: &BinaryHeap<HeapItem>,
    live_len: usize,
    live_cap: usize,
    result_n: usize,
    result_payload_cap: usize,
) {
    let inline = heap
        .capacity()
        .saturating_mul(std::mem::size_of::<HeapItem>());
    mem.peak_heap_len = mem.peak_heap_len.max(heap.len());
    mem.peak_heap_capacity = mem.peak_heap_capacity.max(heap.capacity());
    mem.peak_logical = mem.peak_logical.max(heap.len());
    mem.peak_sum_edge_len = mem.peak_sum_edge_len.max(live_len / EDGE_PAIR_BYTES);
    mem.peak_sum_edge_cap = mem.peak_sum_edge_cap.max(live_cap / EDGE_PAIR_BYTES);
    mem.peak_inline_bytes = mem.peak_inline_bytes.max(inline);
    mem.peak_payload_len_bytes = mem.peak_payload_len_bytes.max(live_len);
    mem.peak_payload_cap_bytes = mem.peak_payload_cap_bytes.max(live_cap);
    mem.peak_frontier_bytes = mem.peak_frontier_bytes.max(inline.saturating_add(live_cap));
    mem.peak_result_n = mem.peak_result_n.max(result_n);
    mem.peak_result_payload_cap_bytes = mem.peak_result_payload_cap_bytes.max(result_payload_cap);
    mem.peak_frontier_plus_results_bytes = mem.peak_frontier_plus_results_bytes.max(
        inline
            .saturating_add(live_cap)
            .saturating_add(result_payload_cap),
    );
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
            mem: SeqKbestFrontierMem::default(),
            pop_count: 0,
            skip_heap_full_at_pop: 0,
            skip_heap_full_at_expand: 0,
            skip_expansion_cap: 0,
            skip_path_edge_cap: 0,
            vertex_visit_refused: 0,
            heap_remaining: 0,
            needles_in_result,
            needles_on_remaining_heap,
            first_expand_refuse: None,
            heap_first_hit_cap_pop: None,
            expand_refuses: Vec::new(),
            first_refuse_snapshot: None,
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
    let mut first_expand_refuse: Option<SeqKbestExpandRefuse> = None;
    let mut heap_first_hit_cap_pop: Option<usize> = None;
    let mut expand_refuses: Vec<SeqKbestExpandRefuse> = Vec::new();
    let mut first_refuse_snapshot: Option<SeqKbestFirstRefuseSnapshot> = None;
    let mut mem = SeqKbestFrontierMem::default();
    mem.rss_before_mib = crate::runtime_config::current_rss_mib();
    sample_frontier_rss(&mut mem);
    let mut result_payload_cap = 0usize;
    scan_copied_frontier(&mut mem, &heap, 0, 0);

    while !heap.is_empty() && result.len() < max_results {
        if crate::runtime_config::hc_rss_abort_triggered() {
            mem.rss_aborted = true;
            break;
        }
        max_heap = max_heap.max(heap.len());
        pop_count += 1;
        if pop_count == 1 || pop_count % 32 == 0 {
            sample_frontier_rss(&mut mem);
        }
        let Some(item) = heap.pop() else {
            break;
        };
        let path = item.path;
        if sinks.contains(&path.last) {
            result_payload_cap = result_payload_cap
                .saturating_add(path.edges.capacity().saturating_mul(EDGE_PAIR_BYTES));
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
            scan_copied_frontier(&mut mem, &heap, result.len(), result_payload_cap);
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
            let mut inserted_this = 0usize;
            let mut inserted_tos: Vec<(usize, f64, u32)> = Vec::new();
            for (oi, &to) in outs.iter().enumerate() {
                if let Some(support) = graph.edge_support(path.last, to) {
                    let extended = path.extend(graph, to, support, total);
                    let allowed_by_dest = match caps.max_heap_paths_per_dest {
                        Some(m) => heap_apply_per_dest_cap(&mut heap, to, extended.score, m),
                        None => true,
                    };
                    if !allowed_by_dest {
                        continue;
                    }
                    if caps.max_heap_paths.is_some_and(|m| heap.len() >= m) {
                        skip_heap_full_at_expand += 1;
                        let mut refused_tos = Vec::new();
                        for &t in &outs[oi..] {
                            if let Some(sup) = graph.edge_support(path.last, t) {
                                let ext = path.extend(graph, t, sup, total);
                                refused_tos.push((t, ext.score, sup));
                            }
                        }
                        let rec = SeqKbestExpandRefuse {
                            pop_count,
                            expansions,
                            result_len: result.len(),
                            heap_len_at_refuse: heap.len(),
                            parent_score: path.score,
                            parent_n_edges: path.edges.len(),
                            parent_last: path.last,
                            parent_start: path.start,
                            parent_edges: path.edges.clone(),
                            n_outs: outs.len(),
                            n_inserted_this_expand: inserted_this,
                            inserted_tos: inserted_tos.clone(),
                            refused_tos,
                        };
                        if first_expand_refuse.is_none() {
                            first_expand_refuse = Some(rec.clone());
                            first_refuse_snapshot = Some(SeqKbestFirstRefuseSnapshot {
                                vertex_counts: vertex_counts.clone(),
                                heap: heap
                                    .iter()
                                    .map(|it| SeqKbestHeapEntryView {
                                        last: it.path.last,
                                        start: it.path.start,
                                        n_edges: it.path.edges.len(),
                                        score: it.path.score,
                                        is_sink: sinks.contains(&it.path.last),
                                        edges_fnv: fnv1a64_edges(&it.path.edges),
                                        bases_fnv: fnv1a64_bytes(
                                            &graph.path_bases_bytes(it.path.start, &it.path.edges),
                                        ),
                                        prefix8_fnv: fnv1a64_edges(
                                            &it.path.edges[..it.path.edges.len().min(8)],
                                        ),
                                        prefix16_fnv: fnv1a64_edges(
                                            &it.path.edges[..it.path.edges.len().min(16)],
                                        ),
                                    })
                                    .collect(),
                                inserted_tos: rec.inserted_tos.clone(),
                            });
                        }
                        expand_refuses.push(rec);
                        break;
                    }
                    inserted_tos.push((to, extended.score, support));
                    heap.push(HeapItem {
                        score: extended.score,
                        tie: extended.edge_count,
                        path: extended,
                    });
                    inserted_this += 1;
                    max_heap = max_heap.max(heap.len());
                    if heap_first_hit_cap_pop.is_none()
                        && caps.max_heap_paths.is_some_and(|m| heap.len() == m)
                    {
                        heap_first_hit_cap_pop = Some(pop_count);
                    }
                }
            }
            scan_copied_frontier(&mut mem, &heap, result.len(), result_payload_cap);
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

    mem.rss_after_mib = crate::runtime_config::current_rss_mib();
    sample_frontier_rss(&mut mem);

    Ok(SeqKbestForensicReport {
        paths: result,
        expansions,
        max_heap,
        mem,
        pop_count,
        skip_heap_full_at_pop,
        skip_heap_full_at_expand,
        skip_expansion_cap,
        skip_path_edge_cap,
        vertex_visit_refused,
        heap_remaining,
        needles_in_result,
        needles_on_remaining_heap,
        first_expand_refuse,
        heap_first_hit_cap_pop,
        expand_refuses,
        first_refuse_snapshot,
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

    #[test]
    fn resource_policy_default_matches_explicit_legacy_on_diamond() {
        let seq = diamond_40_vs_5();
        let def = find_best_haplotypes_seq_graph(&seq, 2).expect("default");
        let legacy =
            find_best_haplotypes_seq_graph_with_policy(&seq, 2, SeqKbestResourcePolicy::Legacy1024)
                .expect("legacy");
        assert_eq!(legacy.diagnostics.policy, "legacy_1024");
        assert!(!legacy.diagnostics.resource_limit_hit);
        assert_eq!(def.len(), legacy.paths.len());
        for i in 0..def.len() {
            assert_eq!(def[i].edges, legacy.paths[i].edges);
            assert!((def[i].score - legacy.paths[i].score).abs() < 1e-15);
        }
    }

    #[test]
    fn unbounded_disables_heap_len_refusal_predicate() {
        let snap_full = CopiedFrontierSnap {
            heap_entries: 1024,
            frontier_bytes: 1,
            pathstate_bytes: 1,
            payload_cap_bytes: 1,
            payload_len_bytes: 1,
            heap_vec_capacity_bytes: 1,
        };
        assert!(blocks_expand(SeqKbestResourcePolicy::Legacy1024, snap_full));
        assert!(!blocks_expand(
            SeqKbestResourcePolicy::UnboundedDiagnostic,
            snap_full
        ));
        assert!(blocks_expand(
            SeqKbestResourcePolicy::ByteBudget {
                max_frontier_bytes: 1
            },
            CopiedFrontierSnap {
                heap_entries: 2,
                frontier_bytes: 1,
                pathstate_bytes: 1,
                payload_cap_bytes: 0,
                payload_len_bytes: 0,
                heap_vec_capacity_bytes: 1,
            }
        ));
    }

    #[test]
    fn byte_budget_can_be_selected_and_does_not_equal_logical_count() {
        let seq = diamond_40_vs_5();
        let generous = find_best_haplotypes_seq_graph_with_policy(
            &seq,
            2,
            SeqKbestResourcePolicy::ByteBudget {
                max_frontier_bytes:
                    crate::seq_kbest_resource_policy::DIAGNOSTIC_GENEROUS_BYTE_BUDGET,
            },
        )
        .expect("generous");
        assert_eq!(generous.diagnostics.policy, "byte_budget");
        assert!(!generous.diagnostics.resource_limit_hit);
        assert_eq!(generous.paths.len(), 2);
        assert_ne!(
            generous.diagnostics.peak_frontier_bytes, generous.diagnostics.peak_logical_frontier,
            "byte accounting must not silently equal logical-state count"
        );
        let tiny = find_best_haplotypes_seq_graph_with_policy(
            &seq,
            2,
            SeqKbestResourcePolicy::ByteBudget {
                max_frontier_bytes: 1,
            },
        )
        .expect("tiny");
        assert!(tiny.diagnostics.resource_limit_hit);
        assert!(tiny.diagnostics.paths_refused > 0);
    }
}
