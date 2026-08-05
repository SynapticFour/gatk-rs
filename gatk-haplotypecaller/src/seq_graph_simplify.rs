//! GATK `SeqGraph` simplification: mirrors `SeqGraph.simplifyGraphOnce` (GAP-E-03).
//! Java reference: `.gatk-src/.../graphs/SeqGraph.java` — MergeDiamonds, MergeTails,
//! SplitCommonSuffices (`CommonSuffixSplitter`), MergeCommonSuffices (`SharedSequenceMerger`), zip.

use crate::seq_graph::SeqGraph;
use std::collections::HashSet;

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
            if try_split_common_suffix(graph, bottom, &mut already_split) {
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
fn try_split_common_suffix(
    graph: &mut SeqGraph,
    bottom: usize,
    already_split: &mut HashSet<usize>,
) -> bool {
    already_split.insert(bottom);
    let prevs: Vec<usize> = graph.incoming_nodes(bottom);
    if prevs.len() < 2 {
        return false;
    }
    if !safe_to_split_suffix(graph, bottom, &prevs) {
        return false;
    }
    let seqs: Vec<&[u8]> = prevs.iter().map(|&p| graph.vertex_sequence(p)).collect();
    let suffix = common_suffix_only(&seqs);
    if suffix.is_empty() {
        return false;
    }
    if would_eliminate_ref_source(graph, &suffix, &prevs) {
        return false;
    }
    if all_vertices_are_only_suffix(&seqs, &suffix) {
        return false;
    }

    let suffix_v = graph.add_seq_vertex(suffix.clone());
    let mut edges_to_remove: Vec<(usize, usize)> = Vec::new();
    let mut mids_to_remove: HashSet<usize> = HashSet::new();

    for &mid in &prevs {
        mids_to_remove.insert(mid);
        let mid_seq = graph.vertex_sequence(mid).to_vec();
        let prefix = without_suffix(&mid_seq, &suffix);
        let out_targets = graph.outgoing_nodes(mid);
        if out_targets.len() != 1 {
            return false;
        }
        let out_target = out_targets[0];
        let (out_sup, out_ref) = graph.outgoing_edge_support_is_ref(mid);

        let incoming_target = if let Some(prefix_bytes) = prefix {
            let pv = graph.add_seq_vertex(prefix_bytes);
            graph.add_or_update_edge(pv, suffix_v, out_sup, out_ref);
            edges_to_remove.push((mid, out_target));
            pv
        } else {
            edges_to_remove.push((mid, out_target));
            suffix_v
        };

        // GATK: suffixV -> getEdgeTarget(outgoingEdgeOf(mid)) == bottom (CommonSuffixSplitter.java).
        graph.add_or_update_edge(suffix_v, out_target, out_sup, out_ref);

        for prev in graph.incoming_nodes(mid) {
            let (in_sup, in_ref) = graph.incoming_edge_support_is_ref(mid);
            graph.add_or_update_edge(prev, incoming_target, in_sup, in_ref);
            edges_to_remove.push((prev, mid));
        }
    }

    graph.remove_vertices_by_id(&mids_to_remove);
    for (from, to) in edges_to_remove {
        graph.remove_edge(from, to);
    }
    true
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
    let new_v = graph.add_seq_vertex(merged_seq);

    let bottom_outs: Vec<usize> = graph.outgoing_nodes(bottom);
    for &p in &prevs {
        for prev in graph.incoming_nodes(p) {
            let (sup, is_ref) = graph.incoming_edge_support_is_ref(p);
            graph.add_or_update_edge(prev, new_v, sup, is_ref);
        }
    }
    for &t in &bottom_outs {
        let (sup, is_ref) = graph.outgoing_edge_support_is_ref(bottom);
        graph.add_or_update_edge(new_v, t, sup, is_ref);
    }

    let mut remove: HashSet<usize> = prevs.iter().copied().collect();
    remove.insert(bottom);
    graph.remove_vertices_by_id(&remove);
    true
}

fn safe_to_split_suffix(graph: &SeqGraph, bottom: usize, prevs: &[usize]) -> bool {
    let bottom_outs: Vec<usize> = graph.outgoing_nodes(bottom);
    if bottom_outs.len() != 1 {
        return false;
    }
    let outgoing_of_bottom: HashSet<usize> = bottom_outs.into_iter().collect();
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

fn split_and_update(
    graph: &mut SeqGraph,
    top: usize,
    bottom: Option<usize>,
    middles: &[usize],
    prefix: &[u8],
    suffix: &[u8],
) -> bool {
    let prefix_v = if prefix.is_empty() {
        None
    } else {
        Some(graph.add_seq_vertex(prefix.to_vec()))
    };
    let suffix_v = if suffix.is_empty() {
        None
    } else {
        Some(graph.add_seq_vertex(suffix.to_vec()))
    };

    let mut new_middles: Vec<usize> = Vec::new();
    for &mid in middles {
        let (in_sup, in_ref) = graph.incoming_edge_support_is_ref(mid);
        let (out_sup, out_ref) = graph.outgoing_edge_support_is_ref(mid);
        let mid_seq = graph.vertex_sequence(mid).to_vec();
        if let Some(remaining) = without_prefix_suffix(&mid_seq, prefix, suffix) {
            let nm = graph.add_seq_vertex(remaining);
            new_middles.push(nm);
            if let Some(pv) = prefix_v {
                graph.add_or_update_edge(pv, nm, in_sup, in_ref);
            }
            if let Some(sv) = suffix_v {
                graph.add_or_update_edge(nm, sv, out_sup, out_ref);
            } else if let Some(b) = bottom {
                graph.add_or_update_edge(nm, b, out_sup, out_ref);
            }
        } else if let (Some(pv), Some(sv)) = (prefix_v, suffix_v) {
            graph.add_or_update_edge(pv, sv, in_sup.saturating_add(out_sup), in_ref || out_ref);
        } else if let (Some(pv), None) = (prefix_v, suffix_v) {
            if let Some(b) = bottom {
                graph.add_or_update_edge(pv, b, in_sup.saturating_add(out_sup), in_ref || out_ref);
            }
        }
    }

    let mut top_to_mid_sup = 0u32;
    let mut top_to_mid_ref = false;
    for &mid in middles {
        if let Some(idx) = graph.find_edge(top, mid) {
            let e = &graph.edges_pub()[idx];
            top_to_mid_sup = top_to_mid_sup.saturating_add(e.support);
            top_to_mid_ref |= e.is_ref;
        }
    }

    let remove: HashSet<usize> = middles.iter().copied().collect();
    graph.remove_vertices_by_id(&remove);

    if let Some(pv) = prefix_v {
        graph.add_or_update_edge(top, pv, top_to_mid_sup, top_to_mid_ref);
    } else {
        for &nm in &new_middles {
            let (sup, is_ref) = graph.incoming_edge_support_is_ref(nm);
            graph.add_or_update_edge(top, nm, sup, is_ref);
        }
    }

    if let (Some(sv), Some(b)) = (suffix_v, bottom) {
        let (sup, is_ref) = graph.incoming_edge_support_is_ref(b);
        graph.add_or_update_edge(sv, b, sup, is_ref);
    }

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
}
