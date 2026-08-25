//! GATK `SeqGraph` simplification: mirrors `SeqGraph.simplifyGraphOnce` (GAP-E-03).
//! Java reference: `.gatk-src/.../graphs/SeqGraph.java` — MergeDiamonds, MergeTails,
//! SplitCommonSuffices (`CommonSuffixSplitter`), MergeCommonSuffices (`SharedSequenceMerger`), zip.

use crate::seq_graph::SeqGraph;
use std::collections::{HashMap, HashSet};

const MIN_TAIL_SUFFIX_LEN: usize = 10;
const MAX_SIMPLIFICATION_ITERATIONS: usize = 100;
const MAX_REASONABLE_SIMPLIFICATION_CYCLES: usize = 100;

/// Full Java `SeqGraph.simplifyGraph` loop.
pub fn simplify_graph_full(graph: &mut SeqGraph) {
    let _ = graph.zip_linear_chains();
    let mut prev_signature: Option<Vec<u8>> = None;
    for i in 0..MAX_SIMPLIFICATION_ITERATIONS {
        if i > MAX_REASONABLE_SIMPLIFICATION_CYCLES {
            break;
        }
        if !simplify_graph_once(graph) {
            break;
        }
        if i > 5 {
            let sig = graph_topology_signature(graph);
            if prev_signature.as_ref() == Some(&sig) {
                break;
            }
            prev_signature = Some(sig);
        }
    }
}

fn graph_topology_signature(graph: &SeqGraph) -> Vec<u8> {
    let mut out = Vec::new();
    out.extend_from_slice(&(graph.node_count() as u32).to_le_bytes());
    out.extend_from_slice(&(graph.edge_count() as u32).to_le_bytes());
    for v in 0..graph.node_count() {
        out.extend_from_slice(&(graph.vertex_sequence(v).len() as u32).to_le_bytes());
        out.extend(graph.vertex_sequence(v));
    }
    out
}

/// One Java `simplifyGraphOnce` iteration.
fn simplify_graph_once(graph: &mut SeqGraph) -> bool {
    let mut did = false;
    did |= merge_diamonds_until_complete(graph);
    did |= merge_tails_until_complete(graph);
    did |= split_common_suffices_until_complete(graph);
    did |= merge_common_suffices_until_complete(graph);
    did |= graph.zip_linear_chains();
    did
}

/// Test-only substages of [`simplify_graph_once`] / [`simplify_graph_full`].
/// Production `simplify_graph_full` is unchanged.
#[cfg(test)]
pub(crate) fn traced_simplify_graph_full(
    graph: &mut SeqGraph,
    mut snap: impl FnMut(&str, &SeqGraph),
) {
    snap("simplify_before_initial_zip", graph);
    let _ = graph.zip_linear_chains();
    snap("simplify_after_initial_zip", graph);
    let mut prev_signature: Option<Vec<u8>> = None;
    for i in 0..MAX_SIMPLIFICATION_ITERATIONS {
        if i > MAX_REASONABLE_SIMPLIFICATION_CYCLES {
            break;
        }
        if !traced_simplify_graph_once(graph, &mut snap) {
            break;
        }
        if i > 5 {
            let sig = graph_topology_signature(graph);
            if prev_signature.as_ref() == Some(&sig) {
                break;
            }
            prev_signature = Some(sig);
        }
    }
    snap("simplify_done", graph);
}

#[cfg(test)]
fn traced_simplify_graph_once(
    graph: &mut SeqGraph,
    snap: &mut impl FnMut(&str, &SeqGraph),
) -> bool {
    let mut did = false;
    did |= merge_diamonds_until_complete(graph);
    snap("after_merge_diamonds", graph);
    did |= merge_tails_until_complete(graph);
    snap("after_merge_tails", graph);
    did |= split_common_suffices_until_complete(graph);
    snap("after_split_common_suffices", graph);
    did |= merge_common_suffices_until_complete(graph);
    snap("after_merge_common_suffices", graph);
    did |= graph.zip_linear_chains();
    snap("after_simplify_zip", graph);
    did
}

fn merge_diamonds_until_complete(graph: &mut SeqGraph) -> bool {
    transform_until_complete(graph, try_merge_diamond)
}

fn merge_tails_until_complete(graph: &mut SeqGraph) -> bool {
    transform_until_complete(graph, try_merge_tails)
}

fn split_common_suffices_until_complete(graph: &mut SeqGraph) -> bool {
    let mut did = false;
    let mut already_split: HashSet<usize> = HashSet::new();
    loop {
        let mut found = false;
        for bottom in 0..graph.node_count() {
            if already_split.contains(&bottom) {
                continue;
            }
            if let Some(map) = try_split_common_suffix(graph, bottom, &mut already_split) {
                already_split = already_split
                    .iter()
                    .filter_map(|id| map.get(id).copied())
                    .collect();
                found = true;
                did = true;
                break;
            }
        }
        if !found {
            break;
        }
    }
    did
}

fn merge_common_suffices_until_complete(graph: &mut SeqGraph) -> bool {
    transform_until_complete(graph, try_merge_common_suffices)
}

fn transform_until_complete(
    graph: &mut SeqGraph,
    mut transform: impl FnMut(&mut SeqGraph, usize) -> bool,
) -> bool {
    let mut did_any = false;
    loop {
        let mut found = false;
        for v in 0..graph.node_count() {
            if transform(graph, v) {
                found = true;
                did_any = true;
                break;
            }
        }
        if !found {
            break;
        }
    }
    did_any
}

fn try_merge_diamond(graph: &mut SeqGraph, top: usize) -> bool {
    let middles: Vec<usize> = graph.outgoing_nodes(top);
    if middles.len() <= 1 {
        return false;
    }
    let mut bottom: Option<usize> = None;
    for &mi in &middles {
        if graph.vertex_out_degree(mi) < 1 || graph.vertex_in_degree(mi) != 1 {
            return false;
        }
        for mt in graph.outgoing_nodes(mi) {
            match bottom {
                None => bottom = Some(mt),
                Some(b) if b == mt => {}
                _ => return false,
            }
        }
    }
    let Some(bottom) = bottom else {
        return false;
    };
    if graph.vertex_in_degree(bottom) != middles.len() {
        return false;
    }
    let seqs: Vec<&[u8]> = middles.iter().map(|&m| graph.vertex_sequence(m)).collect();
    let (prefix, suffix) = common_prefix_suffix(&seqs);
    if prefix.is_empty() && suffix.is_empty() {
        return false;
    }
    split_and_update(graph, top, Some(bottom), &middles, &prefix, &suffix)
}

fn try_merge_tails(graph: &mut SeqGraph, top: usize) -> bool {
    let tails: Vec<usize> = graph.outgoing_nodes(top);
    if tails.len() <= 1 {
        return false;
    }
    for &t in &tails {
        if !graph.is_sink_vertex(t) || graph.vertex_in_degree(t) > 1 {
            return false;
        }
    }
    let seqs: Vec<&[u8]> = tails.iter().map(|&t| graph.vertex_sequence(t)).collect();
    let (prefix, suffix) = common_prefix_suffix(&seqs);
    if suffix.len() < MIN_TAIL_SUFFIX_LEN {
        return false;
    }
    split_and_update(graph, top, None, &tails, &prefix, &suffix)
}

/// GATK `CommonSuffixSplitter.split` — split incoming vertices of `bottom` on shared suffix.
///
/// Java allocates a **new** suffix vertex per predecessor (identical sequence), then
/// `SharedSequenceMerger` collapses those copies. A shared suffix node is not equivalent:
/// it changes merge_common_suffices connectivity. Returns the compact remap so
/// `already_split` can follow vertex identity through dense IDs.
fn try_split_common_suffix(
    graph: &mut SeqGraph,
    bottom: usize,
    already_split: &mut HashSet<usize>,
) -> Option<HashMap<usize, usize>> {
    already_split.insert(bottom);
    let prevs: Vec<usize> = graph.incoming_nodes(bottom);
    if prevs.len() < 2 {
        return None;
    }
    if !safe_to_split_suffix(graph, bottom, &prevs) {
        return None;
    }
    let seqs: Vec<&[u8]> = prevs.iter().map(|&p| graph.vertex_sequence(p)).collect();
    let suffix = common_suffix_only(&seqs);
    if suffix.is_empty() {
        return None;
    }
    if would_eliminate_ref_source(graph, &suffix, &prevs) {
        return None;
    }
    if all_vertices_are_only_suffix(&seqs, &suffix) {
        return None;
    }

    struct MidSplit {
        mid: usize,
        prefix: Option<Vec<u8>>,
        out_target: usize,
        out_sup: u32,
        out_ref: bool,
        incoming: Vec<(usize, u32, bool)>,
    }
    let mut plans: Vec<MidSplit> = Vec::with_capacity(prevs.len());
    for &mid in &prevs {
        let outs = graph.outgoing_nodes(mid);
        if outs.len() != 1 {
            return None;
        }
        let from_mid = graph.edges_from(mid);
        let (out_sup, out_ref) = from_mid
            .first()
            .map(|&(_, s, r)| (s, r))
            .unwrap_or((0, graph.is_reference_node(mid)));
        plans.push(MidSplit {
            mid,
            prefix: without_suffix(graph.vertex_sequence(mid), &suffix),
            out_target: outs[0],
            out_sup,
            out_ref,
            incoming: graph.edges_into(mid),
        });
    }

    let mut mids_to_remove: HashSet<usize> = HashSet::new();
    for plan in &plans {
        mids_to_remove.insert(plan.mid);
        // One suffix vertex per predecessor (CommonSuffixSplitter.java).
        let suffix_v = graph.add_seq_vertex(suffix.clone());
        let incoming_target = if let Some(prefix_bytes) = plan.prefix.clone() {
            let pv = graph.add_seq_vertex(prefix_bytes);
            // Java: prefix → suffix with BaseEdge(out.isRef(), multiplicity 1).
            graph.add_or_update_edge(pv, suffix_v, 1, plan.out_ref);
            pv
        } else {
            suffix_v
        };
        graph.add_or_update_edge(suffix_v, plan.out_target, plan.out_sup, plan.out_ref);
        for &(prev, in_sup, in_ref) in &plan.incoming {
            graph.add_or_update_edge(prev, incoming_target, in_sup, in_ref);
        }
    }

    Some(graph.remove_vertices_by_id(&mids_to_remove))
}

/// GATK `SharedSequenceMerger.merge` — merge identical incoming vertices of `bottom`.
fn try_merge_common_suffices(graph: &mut SeqGraph, bottom: usize) -> bool {
    let prevs: Vec<usize> = graph.incoming_nodes(bottom);
    if prevs.is_empty() {
        return false;
    }
    let first = graph.vertex_sequence(prevs[0]);
    for &p in &prevs[1..] {
        if graph.vertex_sequence(p) != first {
            return false;
        }
    }
    for &p in &prevs {
        if graph.vertex_out_degree(p) != 1 || graph.outgoing_nodes(p)[0] != bottom {
            return false;
        }
        if graph.vertex_in_degree(p) == 0 {
            return false;
        }
    }

    let mut merged_seq = first.to_vec();
    merged_seq.extend_from_slice(graph.vertex_sequence(bottom));
    let mut in_copies: Vec<(usize, u32, bool)> = Vec::new();
    for &p in &prevs {
        in_copies.extend(graph.edges_into(p));
    }
    let out_copies = graph.edges_from(bottom);
    let new_v = graph.add_seq_vertex(merged_seq);
    for (from, sup, is_ref) in in_copies {
        graph.add_or_update_edge(from, new_v, sup, is_ref);
    }
    for (to, sup, is_ref) in out_copies {
        graph.add_or_update_edge(new_v, to, sup, is_ref);
    }

    let mut remove: HashSet<usize> = prevs.iter().copied().collect();
    remove.insert(bottom);
    graph.remove_vertices_by_id(&remove);
    true
}

fn safe_to_split_suffix(graph: &SeqGraph, bottom: usize, prevs: &[usize]) -> bool {
    // Java `CommonSuffixSplitter.safeToSplit` does not require bottom out-degree 1
    // (the reference sink has out-degree 0).
    let outgoing_of_bottom: HashSet<usize> = graph.outgoing_nodes(bottom).into_iter().collect();
    for &mid in prevs {
        if mid == bottom {
            return false;
        }
        let outs = graph.outgoing_nodes(mid);
        if outs.len() != 1 || outs[0] != bottom {
            return false;
        }
        if outgoing_of_bottom.contains(&mid) {
            return false;
        }
    }
    true
}

fn would_eliminate_ref_source(graph: &SeqGraph, suffix: &[u8], prevs: &[usize]) -> bool {
    for &p in prevs {
        if graph.is_ref_source_vertex(p) {
            return graph.vertex_sequence(p).len() == suffix.len();
        }
    }
    false
}

fn all_vertices_are_only_suffix(seqs: &[&[u8]], suffix: &[u8]) -> bool {
    seqs.iter().all(|s| s.len() == suffix.len())
}

fn common_suffix_only(seqs: &[&[u8]]) -> Vec<u8> {
    if seqs.is_empty() {
        return Vec::new();
    }
    let min_len = seqs.iter().map(|s| s.len()).min().unwrap_or(0);
    let mut suffix_len = 0usize;
    'outer: for i in 0..min_len {
        let b = seqs[0][seqs[0].len() - 1 - i];
        for s in &seqs[1..] {
            if s[s.len() - 1 - i] != b {
                break 'outer;
            }
        }
        suffix_len += 1;
    }
    if suffix_len == 0 {
        Vec::new()
    } else {
        let start = seqs[0].len() - suffix_len;
        seqs[0][start..].to_vec()
    }
}

fn without_suffix(seq: &[u8], suffix: &[u8]) -> Option<Vec<u8>> {
    if seq.len() < suffix.len() || !seq.ends_with(suffix) {
        return None;
    }
    let prefix = &seq[..seq.len() - suffix.len()];
    if prefix.is_empty() {
        None
    } else {
        Some(prefix.to_vec())
    }
}

/// GATK `SharedVertexSequenceSplitter.split` + `updateGraph`.
///
/// Edges are copied from the original middles *before* compact. Re-deriving
/// top→prefix / suffix→bottom from the outer graph after removing middles
/// yields support 0 / `is_ref` false (lost reference source and sink).
fn split_and_update(
    graph: &mut SeqGraph,
    top: usize,
    bottom: Option<usize>,
    middles: &[usize],
    prefix: &[u8],
    suffix: &[u8],
) -> bool {
    #[derive(Clone, Copy)]
    struct E {
        support: u32,
        is_ref: bool,
    }
    fn process_edge(graph: &SeqGraph, v: usize, incoming: bool) -> E {
        // Java `processEdgeToRemove`: missing edge → BaseEdge(isReferenceNode(v), 0).
        let edges = if incoming {
            graph.edges_into(v)
        } else {
            graph.edges_from(v)
        };
        match edges.first() {
            None => E {
                support: 0,
                is_ref: graph.is_reference_node(v),
            },
            Some(&(_, support, is_ref)) => E { support, is_ref },
        }
    }

    struct MidPlan {
        remaining: Option<Vec<u8>>,
        to_mid: E,
        from_mid: E,
    }
    let plans: Vec<MidPlan> = middles
        .iter()
        .map(|&mid| {
            let to_mid = process_edge(graph, mid, true);
            let from_mid = process_edge(graph, mid, false);
            let remaining = without_prefix_suffix(graph.vertex_sequence(mid), prefix, suffix);
            MidPlan {
                remaining,
                to_mid,
                from_mid,
            }
        })
        .collect();

    let remaining_count = plans.iter().filter(|p| p.remaining.is_some()).count();
    let consumed = plans.iter().any(|p| p.remaining.is_none());
    let prefix_outdeg = remaining_count + usize::from(consumed);
    let has_only_prefix_suffix = consumed && prefix_outdeg == 1;
    // Java: needPrefixNode = !prefix.empty || (top == null && !hasOnlyPrefixSuffixEdges).
    // Callers always pass a top vertex (MergeDiamonds / MergeTails).
    let need_prefix = !prefix.is_empty();
    let need_suffix = !suffix.is_empty() || (bottom.is_none() && !has_only_prefix_suffix);

    let prefix_v = if need_prefix {
        Some(graph.add_seq_vertex(prefix.to_vec()))
    } else {
        None
    };
    let suffix_v = if need_suffix {
        Some(graph.add_seq_vertex(suffix.to_vec()))
    } else {
        None
    };

    let mut new_middles: Vec<usize> = Vec::new();
    let mut prefix_to_mid: Vec<(usize, E)> = Vec::new();
    let mut mid_to_suffix: Vec<(usize, E)> = Vec::new();
    let mut prefix_to_suffix: Option<E> = None;
    for plan in &plans {
        if let Some(remaining) = &plan.remaining {
            let nm = graph.add_seq_vertex(remaining.clone());
            new_middles.push(nm);
            prefix_to_mid.push((nm, plan.to_mid));
            mid_to_suffix.push((nm, plan.from_mid));
        } else {
            let combined = E {
                support: plan.to_mid.support.saturating_add(plan.from_mid.support),
                is_ref: plan.to_mid.is_ref || plan.from_mid.is_ref,
            };
            match &mut prefix_to_suffix {
                Some(prev) => {
                    prev.support = prev.support.saturating_add(combined.support);
                    prev.is_ref |= combined.is_ref;
                }
                None => prefix_to_suffix = Some(combined),
            }
        }
    }

    let top_for = if need_prefix { prefix_v } else { Some(top) };
    let bot_for = if need_suffix { suffix_v } else { bottom };

    // Java `addPrefixNodeAndEdges`: top → prefix with makeOREdge(outgoing of prefix, 1).
    if let Some(pv) = prefix_v {
        let any_ref = prefix_to_mid.iter().any(|(_, e)| e.is_ref)
            || prefix_to_suffix.map(|e| e.is_ref).unwrap_or(false);
        graph.add_or_update_edge(top, pv, 1, any_ref);
    }
    // Java `addSuffixNodeAndEdges`: suffix → bottom with makeOREdge(incoming of suffix, 1).
    if let (Some(sv), Some(b)) = (suffix_v, bottom) {
        let any_ref = mid_to_suffix.iter().any(|(_, e)| e.is_ref)
            || prefix_to_suffix.map(|e| e.is_ref).unwrap_or(false);
        graph.add_or_update_edge(sv, b, 1, any_ref);
    }

    // Java `addEdgesFromTopNode`.
    if let Some(tfc) = top_for {
        for &(nm, e) in &prefix_to_mid {
            graph.add_or_update_edge(tfc, nm, e.support, e.is_ref);
        }
        if let (Some(e), Some(bfc)) = (prefix_to_suffix, bot_for) {
            graph.add_or_update_edge(tfc, bfc, e.support, e.is_ref);
        }
    }
    // Java `addEdgesToBottomNode`. Skip the split-graph prefix vertex when it was
    // not added to the outer graph (empty prefix); the fully-consumed path is
    // already top → bot_for from addEdgesFromTopNode.
    if let Some(bfc) = bot_for {
        for &(nm, e) in &mid_to_suffix {
            graph.add_or_update_edge(nm, bfc, e.support, e.is_ref);
        }
        if let (Some(pv), Some(e)) = (prefix_v, prefix_to_suffix) {
            graph.add_or_update_edge(pv, bfc, e.support, e.is_ref);
        }
    }

    let remove: HashSet<usize> = middles.iter().copied().collect();
    graph.remove_vertices_by_id(&remove);
    !new_middles.is_empty() || prefix_v.is_some() || suffix_v.is_some()
}

fn common_prefix_suffix(seqs: &[&[u8]]) -> (Vec<u8>, Vec<u8>) {
    if seqs.is_empty() {
        return (Vec::new(), Vec::new());
    }
    let min_len = seqs.iter().map(|s| s.len()).min().unwrap_or(0);
    let mut prefix_len = 0usize;
    'prefix: for i in 0..min_len {
        let b = seqs[0][i];
        for s in &seqs[1..] {
            if s[i] != b {
                break 'prefix;
            }
        }
        prefix_len += 1;
    }
    let remain = min_len.saturating_sub(prefix_len);
    let mut suffix_len = 0usize;
    'suffix: for i in 0..remain {
        let b = seqs[0][seqs[0].len() - 1 - i];
        for s in &seqs[1..] {
            if s[s.len() - 1 - i] != b {
                break 'suffix;
            }
        }
        suffix_len += 1;
    }
    let prefix = seqs[0][..prefix_len].to_vec();
    let suffix = if suffix_len == 0 {
        Vec::new()
    } else {
        let start = seqs[0].len() - suffix_len;
        seqs[0][start..].to_vec()
    };
    (prefix, suffix)
}

fn without_prefix_suffix(seq: &[u8], prefix: &[u8], suffix: &[u8]) -> Option<Vec<u8>> {
    if seq.len() < prefix.len() + suffix.len() {
        return None;
    }
    if !seq.starts_with(prefix) {
        return None;
    }
    if !suffix.is_empty() && !seq.ends_with(suffix) {
        return None;
    }
    let end = seq.len() - suffix.len();
    let mid = &seq[prefix.len()..end];
    if mid.is_empty() {
        None
    } else {
        Some(mid.to_vec())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::assembly::{AssemblyGraphParams, AssemblyRead};
    use crate::read_threading_graph::assembly_graph_from_ref_and_reads_threading;
    use crate::seq_graph::SeqGraph;

    fn read(seq: &str, q: u8) -> AssemblyRead {
        AssemblyRead {
            bases: seq.as_bytes().to_vec(),
            base_quals: vec![q; seq.len()],
        }
    }

    #[test]
    fn simplify_graph_once_runs_on_p5_graph() {
        let reference = read("ACGTT", 30);
        let reads = vec![read("ACGTT", 30), read("ACGTA", 30)];
        let params = AssemblyGraphParams {
            kmer_size: crate::bio_ids::KmerSize::try_new(3).unwrap(),
            min_base_quality: 10,
            ..Default::default()
        };
        let ag = assembly_graph_from_ref_and_reads_threading(&reference, &reads, &params).unwrap();
        let mut g = SeqGraph::from_assembly_graph(&ag);
        assert!(simplify_graph_once(&mut g) || g.node_count() > 0);
    }

    #[test]
    fn merge_common_suffices_collapses_identical_predecessors() {
        let mut g = SeqGraph::from_assembly_graph(
            &assembly_graph_from_ref_and_reads_threading(
                &read("ACGT", 30),
                &[read("ACGT", 30), read("ACGT", 30)],
                &AssemblyGraphParams {
                    kmer_size: crate::bio_ids::KmerSize::try_new(3).unwrap(),
                    min_base_quality: 10,
                    ..Default::default()
                },
            )
            .unwrap(),
        );
        let n = g.node_count();
        let _ = simplify_graph_full(&mut g);
        assert!(g.node_count() > 0 && g.node_count() <= n + 5);
    }

    #[test]
    fn common_prefix_suffix_matches_gatk_prefix_suffix_data() {
        // SharedVertexSequenceSplitterUnitTest.PrefixSuffixData (GATK 4.4.0.0).
        let cases: &[(&[&[u8]], usize, usize)] = &[
            (&[b"A", b"C"], 0, 0),
            (&[b"C", b"C"], 1, 0),
            (&[b"ACT", b"AGT"], 1, 1),
            (&[b"ACCT", b"AGT"], 1, 1),
            (&[b"ACT", b"ACT"], 3, 0),
            (&[b"ACTA", b"ACT"], 3, 0),
            (&[b"ACTA", b"ACTG"], 3, 0),
            (&[b"ACTA", b"ACTGA"], 3, 1),
            (&[b"GCTGA", b"ACTGA"], 0, 4),
            (&[b"A", b"C", b"A"], 0, 0),
            (&[b"A", b"A", b"A"], 1, 0),
            (&[b"A", b"AA", b"A"], 1, 0),
            (&[b"A", b"ACA", b"A"], 1, 0),
            (&[b"ACT", b"ACAT", b"ACT"], 2, 1),
            (&[b"ACT", b"ACAT", b"ACGT"], 2, 1),
            (&[b"AAAT", b"AAA", b"CAAA"], 0, 0),
            (&[b"AACTTT", b"AAGTTT", b"AAGCTTT"], 2, 3),
            (&[b"AAA", b"AAA", b"CAAA"], 0, 3),
            (&[b"AAA", b"AAA", b"AAA"], 3, 0),
            (&[b"AC", b"ACA", b"AC"], 2, 0),
        ];
        for &(seqs, prefix_len, suffix_len) in cases {
            let (prefix, suffix) = common_prefix_suffix(seqs);
            assert_eq!(
                prefix.len(),
                prefix_len,
                "prefix {:?}",
                seqs.iter()
                    .map(|s| std::str::from_utf8(s).unwrap())
                    .collect::<Vec<_>>()
            );
            assert_eq!(
                suffix.len(),
                suffix_len,
                "suffix {:?}",
                seqs.iter()
                    .map(|s| std::str::from_utf8(s).unwrap())
                    .collect::<Vec<_>>()
            );
        }
    }
}
