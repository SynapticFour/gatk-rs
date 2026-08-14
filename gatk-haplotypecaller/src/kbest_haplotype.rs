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
    pub fn bases(&self, graph: &AssemblyGraph) -> Vec<u8> {
        graph.path_bases(self.start, &self.edges)
    }
}

/// Cap on in-flight k-best frontier paths. Bushy/cyclic local graphs (NA12878
/// `20:10098169-10098441`) otherwise fill this heap with long `PathState`s and reach
/// multi-GiB Peak-RSS — the realistic-window failure mode.
///
/// Keep well below `haplotype_budget × max_path_edges × 16B` ≈ tens of MiB.
const MAX_KBEST_HEAP_PATHS: usize = 1_024;

/// Cap on edges in a single in-flight path. Cyclic graphs otherwise let `edges` grow
/// without bound → GiB-scale PathStates.
///
/// Must clear padded HC spans (active + ~500 bp each side → often 1–2 kb of k-mers,
/// i.e. ≫256 edges). Cap 256 zeroed k-best on normal GIAB windows (chr21 recall
/// collapse: `paths=0` everywhere, dangling fragments only). Cyclic wander remains
/// bounded by heap / expansion caps and per-vertex visit limits.
const MAX_KBEST_PATH_EDGES: usize = 4_096;

#[derive(Debug, Clone)]
struct PathState {
    start: usize,
    edges: Vec<(usize, usize)>,
    last: usize,
    score: f64,
    is_reference: bool,
    /// Cached `path_bases(..).len()` — avoid rebuilding the haplotype string on every heap push.
    bases_len: usize,
}

impl PathState {
    fn new(graph: &AssemblyGraph, start: usize) -> Self {
        Self {
            start,
            edges: Vec::new(),
            last: start,
            score: 0.0,
            is_reference: false,
            bases_len: graph.nodes()[start].kmer.len(),
        }
    }

    fn extend(
        &self,
        graph: &AssemblyGraph,
        to: usize,
        edge_support: u32,
        total_outgoing: u32,
    ) -> Self {
        // Avoid `Vec::clone` on the hot k-best frontier (OOM amplifier on bushy/cyclic graphs).
        let mut edges = Vec::with_capacity(self.edges.len() + 1);
        edges.extend_from_slice(&self.edges);
        edges.push((self.last, to));
        let penalty = log_penalty(edge_support, total_outgoing);
        // Each edge appends the last base of the destination kmer (`AssemblyGraph::path_bases`).
        let add = usize::from(graph.nodes()[to].kmer.last().is_some());
        Self {
            start: self.start,
            edges,
            last: to,
            score: self.score + penalty,
            is_reference: self.is_reference && graph.edge_is_ref(self.last, to),
            bases_len: self.bases_len + add,
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
///
/// On failure returns the **original** graph in [`Err`] so callers can fall back to a
/// preserving search without `graph.clone()` (Peak-RSS on bushy cyclic loci).
///
/// Observable contract: same cycle-edge / non-reaching-vertex removal as GATK; Rust-native
/// cuts avoid per-edge `Vec` allocs and fold isolated cleanup into one remapping pass.
pub fn graph_for_kbest(mut graph: AssemblyGraph) -> Result<AssemblyGraph, AssemblyGraph> {
    if !graph.has_cycle() {
        return Ok(graph);
    }
    let Some(source) = graph.reference_source_vertex() else {
        return Err(graph);
    };
    let Some(sink) = graph.reference_sink_vertex() else {
        return Err(graph);
    };
    let n = graph.node_count();
    let mut is_sink = vec![false; n];
    if sink < n {
        is_sink[sink] = true;
    }
    let mut edges_to_remove = HashSet::new();
    let mut vertices_to_remove = HashSet::new();
    let mut on_path = vec![false; n];
    let found_path = find_cycle_guilty(
        &graph,
        source,
        &is_sink,
        &mut edges_to_remove,
        &mut vertices_to_remove,
        &mut on_path,
    );
    if !found_path || (edges_to_remove.is_empty() && vertices_to_remove.is_empty()) {
        return Err(graph);
    }
    for (f, t) in edges_to_remove {
        graph.remove_edge(f, t);
    }
    // Fold isolated cleanup into one remapping pass (was remove_nodes + cleanup_isolated).
    graph.remove_nodes_and_isolated(vertices_to_remove);
    Ok(graph)
}

fn find_cycle_guilty(
    graph: &AssemblyGraph,
    current: usize,
    is_sink: &[bool],
    edges_to_remove: &mut HashSet<(usize, usize)>,
    vertices_to_remove: &mut HashSet<usize>,
    on_path: &mut [bool],
) -> bool {
    if current < is_sink.len() && is_sink[current] {
        return true;
    }
    if current < on_path.len() {
        on_path[current] = true;
    }
    let mut reaches_sink = false;
    for to in graph.outgoing_targets(current) {
        if to < on_path.len() && on_path[to] {
            edges_to_remove.insert((current, to));
        } else if find_cycle_guilty(
            graph,
            to,
            is_sink,
            edges_to_remove,
            vertices_to_remove,
            on_path,
        ) {
            reaches_sink = true;
        }
    }
    if !reaches_sink {
        vertices_to_remove.insert(current);
    }
    if current < on_path.len() {
        on_path[current] = false;
    }
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

/// Production assembly k-best: prefer an acyclic graph. Cyclic graphs previously kept
/// full topology (P12) and could grow the k-best frontier without bound on bushy
/// loci (e.g. NA12878 `20:10098169`); try cycle removal first, then fall back to the
/// capped preserving search.
///
/// Takes ownership so cycle stripping does not `graph.clone()` (Peak-RSS on bushy loci).
pub fn find_best_haplotypes_for_assembly(
    graph: AssemblyGraph,
    max_number_of_haplotypes: usize,
) -> GatkResult<(Vec<KBestPath>, AssemblyGraph)> {
    // Observable contract: prefer acyclic k-best; when cycle stripping cannot produce a
    // source→sink DAG, fall back to capped preserving search (pre-Peak behavior). Dropping
    // that fallback and running `find_best_haplotypes_inner` on the cyclic graph left L2
    // `g2-subset-live` at haplotype_count=1 (ref only).
    let (paths, graph) = match graph_for_kbest(graph) {
        Ok(acyclic) => {
            crate::runtime_config::rss_trace_checkpoint(
                "kbest_begin",
                &format!(
                    "nodes={} max_haps={} mode=acyclic",
                    acyclic.nodes().len(),
                    max_number_of_haplotypes
                ),
            );
            let paths = find_best_haplotypes_inner(&acyclic, max_number_of_haplotypes, false)?;
            (paths, acyclic)
        }
        Err(cyclic) => {
            // Peak-RSS: cyclic-preserve k-best on bushy graphs (≤8k nodes still) is the
            // multi-GiB amplifier. Keep preserve for small cyclic L2 fixtures only.
            const MAX_CYCLIC_PRESERVE_NODES: usize = 4_000;
            if cyclic.nodes().len() > MAX_CYCLIC_PRESERVE_NODES {
                crate::runtime_config::rss_trace_checkpoint(
                    "kbest_cyclic_skip_huge",
                    &format!(
                        "nodes={} cap={MAX_CYCLIC_PRESERVE_NODES}",
                        cyclic.nodes().len()
                    ),
                );
                (Vec::new(), cyclic)
            } else {
                crate::runtime_config::rss_trace_checkpoint(
                    "kbest_cyclic_preserve",
                    &format!(
                        "nodes={} max_haps={}",
                        cyclic.nodes().len(),
                        max_number_of_haplotypes
                    ),
                );
                let paths =
                    find_best_haplotypes_preserving_cycles(&cyclic, max_number_of_haplotypes)?;
                (paths, cyclic)
            }
        }
    };
    crate::runtime_config::rss_trace_checkpoint("kbest_done", &format!("paths={}", paths.len()));
    Ok((paths, graph))
}

fn find_best_haplotypes_inner(
    graph: &AssemblyGraph,
    max_number_of_haplotypes: usize,
    strip_cycles: bool,
) -> GatkResult<Vec<KBestPath>> {
    if max_number_of_haplotypes == 0 {
        return Ok(Vec::new());
    }
    // Borrow when possible; only clone for cycle stripping (find_best_haplotypes public API).
    let owned;
    let graph: &AssemblyGraph = if strip_cycles {
        owned = match graph_for_kbest(graph.clone()) {
            Ok(g) => g,
            Err(g) => g, // preserve cyclic topology when strip fails
        };
        &owned
    } else {
        graph
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
        let path = PathState::new(graph, s);
        heap.push(HeapItem {
            score_bits: path.score.to_bits(),
            tie: path.bases_len,
            path,
        });
    }
    let mut vertex_counts: HashMap<usize, usize> =
        graph.nodes().iter().map(|n| (n.id, 0)).collect();

    // Bound total expansions: at the heap cap, a bushy/cyclic graph otherwise spins
    // forever in pop/extend/push, fragmenting the allocator into multi-GiB Peak-RSS.
    // 50k was still enough for ≤8k-node dense GIAB shards to climb to ~15 GiB.
    const MAX_KBEST_EXPANSIONS: usize = 12_000;
    let mut expansions = 0usize;
    while !heap.is_empty() && result.len() < max_number_of_haplotypes {
        if crate::runtime_config::hc_rss_abort_triggered() {
            crate::runtime_config::rss_trace_checkpoint(
                "kbest_rss_abort",
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
        if path.edges.len() >= MAX_KBEST_PATH_EDGES || expansions >= MAX_KBEST_EXPANSIONS {
            continue;
        }
        if expansions > 0 && expansions % 5_000 == 0 {
            crate::runtime_config::rss_trace_checkpoint(
                "kbest_expand",
                &format!(
                    "expansions={expansions} heap={} results={}",
                    heap.len(),
                    result.len()
                ),
            );
        }
        // If the frontier is already saturated, only accept sinks (above); do not expand.
        if heap.len() >= MAX_KBEST_HEAP_PATHS {
            continue;
        }
        let count = vertex_counts.get_mut(&path.last).expect("vertex");
        if *count < max_number_of_haplotypes {
            *count += 1;
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
                        tie: extended.bases_len,
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
            .then_with(|| b.bases(graph).cmp(&a.bases(graph)))
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
            bases: seq.as_bytes().to_vec(),
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
        assert!(seqs.contains(&b"ACGTT".to_vec()));
    }

    #[test]
    fn kbest_bounds_are_finite_and_tight() {
        assert!(MAX_KBEST_HEAP_PATHS <= 4_096);
        assert!(MAX_KBEST_PATH_EDGES <= 8_192);
        // Worst-case edge storage alone stays well under 100 MiB.
        let worst_edge_bytes = MAX_KBEST_HEAP_PATHS
            .saturating_mul(MAX_KBEST_PATH_EDGES)
            .saturating_mul(std::mem::size_of::<(usize, usize)>());
        assert!(worst_edge_bytes < 128 * 1024 * 1024);
    }
}
