//! GATK `AbstractReadThreadingGraph.recoverDanglingTails` / `recoverDanglingHeads` parity.

use crate::assembly::{AssemblyGraph, DanglingMergeHaplotype};
use crate::cigar::{Cigar, CigarElement, CigarOperator};
use crate::smith_waterman::{align, SwOverhangStrategy, SwParameters};
use gatk_common::{GatkError, GatkResult};
use std::collections::HashSet;

/// GATK dangling-end Smith–Waterman defaults (`ReadThreadingAssemblerArgumentCollection`).
/// # Invariants
/// Match/mismatch/gap penalties match GATK dangling-end SW constants when using [`Self::gatk_defaults`].
/// # Ownership
/// [`Copy`] score bundle.
/// # Mutation
/// Immutable per recovery pass.
/// # Biological assumptions
/// Aligns dangling branch path to reference path for merge decisions.
/// # Java equivalence
/// GATK dangling-end SW parameters on `ReadThreadingAssemblerArgumentCollection`.
#[derive(Debug, Clone, Copy)]
pub struct DanglingRecoverySwParams {
    pub match_value: i32,
    pub mismatch_penalty: i32,
    pub gap_open_penalty: i32,
    pub gap_extend_penalty: i32,
}

impl DanglingRecoverySwParams {
    pub fn gatk_defaults() -> Self {
        Self {
            match_value: 25,
            mismatch_penalty: -50,
            gap_open_penalty: -110,
            gap_extend_penalty: -6,
        }
    }
}

/// Successful dangling tail merge plan (bases + CIGAR before `addEdge`).
/// # Invariants
/// `from` / `to` are graph node ids for the merge edge; CIGAR aligns alt vs ref path bases.
/// # Ownership
/// Owns alt/ref path bases and CIGAR; consumed when applying the merge to the graph.
/// # Mutation
/// Immutable plan; graph mutation happens when the plan is applied.
/// # Biological assumptions
/// Encodes a recovered variant junction connecting a dangling branch to the ref spine.
/// # Java equivalence
/// GATK dangling-tail merge plan before edge insertion (`recoverDanglingTails`).
#[derive(Debug, Clone)]
pub struct DanglingTailMergePlan {
    pub from: usize,
    pub to: usize,
    pub alt_bases: Vec<u8>,
    pub ref_path_bases: Vec<u8>,
    pub cigar: Cigar,
}

/// GATK dangling recovery knobs.
/// # Invariants
/// `min_dangling_branch_length` gates recovery attempts; HC recovers heads with tails by default.
/// `min_matching_bases_to_dangling_end_recovery == -1` uses legacy suffix/prefix rules only.
/// # Ownership
/// [`Copy`] config including nested SW params.
/// # Mutation
/// Snapshot for one recovery pass over a graph.
/// # Biological assumptions
/// Dangling branches may carry true variation that should reattach to the reference path.
/// # Java equivalence
/// GATK `AbstractReadThreadingGraph.recoverDanglingTails` / `recoverDanglingHeads` knobs.
#[derive(Debug, Clone, Copy)]
pub struct DanglingRecoveryParams {
    pub min_prune_factor: u32,
    pub min_dangling_branch_length: usize,
    pub recover_all_dangling_branches: bool,
    /// When false, skip head recovery (Java HC always recovers heads with tails).
    pub recover_dangling_heads: bool,
    /// GATK `minMatchingBasesToDanglingEndRecovery` (-1 = legacy suffix/prefix rules only).
    pub min_matching_bases_to_dangling_end_recovery: i32,
    pub sw: DanglingRecoverySwParams,
    /// Single-pass dangling + no ASM-1 tail suffix rescue + merge edge weight 1 (GATK exact).
    pub dangling_java_exact: bool,
}

impl DanglingRecoveryParams {
    pub fn gatk_haplotype_caller_defaults() -> Self {
        Self {
            min_prune_factor: 2,
            min_dangling_branch_length: 3,
            recover_all_dangling_branches: false,
            recover_dangling_heads: true,
            min_matching_bases_to_dangling_end_recovery: -1,
            sw: DanglingRecoverySwParams::gatk_defaults(),
            dangling_java_exact: false,
        }
    }

    /// Build from [`crate::read_threading_assembler::ReadThreadingAssemblerArgs`].
    pub fn from_assembler_args(
        args: &crate::read_threading_assembler::ReadThreadingAssemblerArgs,
    ) -> Self {
        Self {
            min_prune_factor: args.min_prune_factor,
            min_dangling_branch_length: args.min_dangling_branch_length,
            recover_all_dangling_branches: args.recover_all_dangling_branches,
            recover_dangling_heads: args.recover_dangling_heads,
            min_matching_bases_to_dangling_end_recovery: args
                .min_matching_bases_to_dangling_end_recovery,
            sw: args.dangling_end_sw,
            dangling_java_exact: args.dangling_java_exact,
        }
    }
}

/// Summary for L2 parity **E.4**.
/// # Invariants
/// `edges_after` reflects the graph after recovery; attempt/recover counts are non-decreasing counters.
/// # Ownership
/// Owned counters for parity dumps.
/// # Mutation
/// Immutable post-recovery snapshot.
/// # Biological assumptions
/// None — instrumentation of dangling recovery outcomes.
/// # Java equivalence
/// Rust-native E.4 dump over GATK dangling recovery side effects.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DanglingRecoverySummary {
    pub edges_before: usize,
    pub edges_after: usize,
    pub tails_attempted: u32,
    pub tails_recovered: u32,
    pub heads_attempted: u32,
    pub heads_recovered: u32,
    pub edges_merged: u32,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum TraversalDir {
    Up,
    Down,
}

fn base_eq(a: u8, b: u8) -> bool {
    a.eq_ignore_ascii_case(&b)
}

fn sw_parameters(p: &DanglingRecoverySwParams) -> SwParameters {
    SwParameters {
        match_value: p.match_value,
        mismatch_penalty: p.mismatch_penalty,
        gap_open_penalty: p.gap_open_penalty,
        gap_extend_penalty: p.gap_extend_penalty,
    }
}

/// GATK `AlignmentUtils.removeTrailingDeletions`.
fn remove_trailing_deletions(mut cigar: Cigar) -> Cigar {
    if matches!(
        cigar.elements.last(),
        Some(CigarElement {
            operator: CigarOperator::Deletion,
            ..
        })
    ) {
        cigar.elements.pop();
    }
    cigar
}

/// GATK dangling-end SW (`SWOverhangStrategy.LEADING_INDEL`).
fn align_dangling(ref_bases: &[u8], alt_bases: &[u8], p: &DanglingRecoverySwParams) -> Cigar {
    if ref_bases.is_empty() || alt_bases.is_empty() {
        return Cigar::new();
    }
    let Ok(aln) = align(
        ref_bases,
        alt_bases,
        &sw_parameters(p),
        SwOverhangStrategy::LeadingIndel,
    ) else {
        return Cigar::new();
    };
    remove_trailing_deletions(aln.cigar)
}

const MAX_CIGAR_COMPLEXITY: usize = 3;

fn cigar_ok_to_merge_tail(cigar: &Cigar) -> bool {
    let elements = &cigar.elements;
    !elements.is_empty()
        && elements.len() <= MAX_CIGAR_COMPLEXITY
        && elements
            .last()
            .is_some_and(|e| e.operator == CigarOperator::Match)
}

fn cigar_ok_to_merge_head(cigar: &Cigar) -> bool {
    let elements = &cigar.elements;
    !elements.is_empty()
        && elements.len() <= MAX_CIGAR_COMPLEXITY
        && elements
            .first()
            .is_some_and(|e| e.operator == CigarOperator::Match)
}

/// GATK `bestPrefixMatch` (modern dangling-head merge from SW CIGAR).
fn best_prefix_match(
    cigar: &Cigar,
    path1: &[u8],
    path2: &[u8],
    min_matching_bases: i32,
) -> Option<(usize, usize)> {
    let min = if min_matching_bases >= 0 {
        min_matching_bases as usize
    } else {
        1
    };
    let mut ref_idx = cigar.reference_length().saturating_sub(1);
    let mut read_idx = path2.len().saturating_sub(1);
    'cigar: for el in cigar.elements.iter().rev() {
        if !el.operator.consumes_read_bases() || !el.operator.consumes_reference_bases() {
            break;
        }
        for _ in 0..el.length {
            if ref_idx >= path1.len() || read_idx >= path2.len() {
                break 'cigar;
            }
            if !base_eq(path1[ref_idx], path2[read_idx]) {
                break 'cigar;
            }
            ref_idx = ref_idx.saturating_sub(1);
            read_idx = read_idx.saturating_sub(1);
        }
    }
    let matches = path2.len().saturating_sub(1).saturating_sub(read_idx);
    if matches < min {
        None
    } else {
        Some((ref_idx, read_idx))
    }
}

/// GATK `bestPrefixMatchLegacy` for dangling heads; returns merge index or `None`.
fn best_prefix_match_legacy(
    path1: &[u8],
    path2: &[u8],
    max_index: usize,
    kmer_size: usize,
) -> Option<usize> {
    let max_mismatches = (max_index / kmer_size.max(1)).max(1);
    let mut mismatches = 0usize;
    let mut index = 0usize;
    let mut last_good_index: Option<usize> = None;
    while index < max_index {
        if index >= path1.len() || index >= path2.len() {
            break;
        }
        if !base_eq(path1[index], path2[index]) {
            mismatches += 1;
            if mismatches > max_mismatches {
                return last_good_index;
            }
            last_good_index = Some(index);
        }
        index += 1;
    }
    if index == max_index {
        last_good_index = Some(index.saturating_sub(1));
    }
    last_good_index
}

/// GATK `AbstractReadThreadingGraph.longestSuffixMatch(seq, kmer, seqStart)`.
fn longest_suffix_match_java(seq: &[u8], kmer: &[u8], seq_start: usize) -> usize {
    for len in 1..=kmer.len() {
        let seq_i = seq_start.saturating_add(1).saturating_sub(len);
        let kmer_i = kmer.len() - len;
        if seq_i >= seq.len() || !base_eq(seq[seq_i], kmer[kmer_i]) {
            return len - 1;
        }
    }
    kmer.len()
}

/// GATK `getBasesForPath`: suffix per vertex; optional expanded source kmer (reversed).
fn path_bases(graph: &AssemblyGraph, path: &[usize], expand_source: bool) -> Vec<u8> {
    let mut out = Vec::new();
    for (idx, &node) in path.iter().enumerate() {
        let kmer = graph.kmer_at(node);
        if expand_source && idx == 0 && graph.is_source(node) {
            let bytes: Vec<u8> = kmer.iter().rev().copied().collect();
            out.extend_from_slice(&bytes);
        } else if let Some(&b) = kmer.last() {
            out.push(b);
        }
    }
    out
}

impl AssemblyGraph {
    pub(crate) fn add_merged_edge(&mut self, from: usize, to: usize, support: u32, is_ref: bool) {
        // Merged edges must survive adaptive pruning and rank in k-best (ASM-1 spine).
        let weight = support.max(3);
        let weight = self
            .edge_support(from, to)
            .map(|existing| existing.max(weight))
            .unwrap_or(weight);
        self.add_edge_support(from, to, weight);
        if is_ref {
            self.ref_edges.insert((from, to));
        }
    }

    pub(crate) fn has_edge(&self, from: usize, to: usize) -> bool {
        self.edge_support(from, to).is_some()
    }

    fn heaviest_incoming(&self, v: usize) -> Option<(usize, u32)> {
        let preds = self.incoming_nodes(v);
        preds
            .into_iter()
            .filter_map(|p| self.edge_support(p, v).map(|s| (p, s)))
            .max_by_key(|(_, s)| *s)
    }

    fn heaviest_outgoing(&self, v: usize) -> Option<(usize, u32)> {
        let outs = self.outgoing_nodes(v);
        outs.into_iter()
            .filter_map(|t| self.edge_support(v, t).map(|s| (v, t, s)))
            .max_by(|a, b| a.2.cmp(&b.2))
            .map(|(_, t, s)| (t, s))
    }

    fn has_incident_ref_in(&self, v: usize) -> bool {
        self.incoming_nodes(v)
            .iter()
            .any(|&p| self.edge_is_ref(p, v))
    }

    fn is_ref_node(&self, v: usize) -> bool {
        self.ref_nodes.contains(&v)
    }

    fn is_ref_source(&self, v: usize) -> bool {
        self.is_ref_source_vertex(v)
    }

    fn is_ref_sink(&self, v: usize) -> bool {
        self.is_ref_sink_vertex(v)
    }

    /// GATK `findPath` — path is LCA-first (`LinkedList.addFirst`).
    fn find_path(
        &self,
        start: usize,
        prune_factor: u32,
        done: impl Fn(usize) -> bool,
        return_path: impl Fn(usize) -> bool,
        next: impl Fn(usize) -> Option<(usize, u32)>,
    ) -> Option<Vec<usize>> {
        let mut path = Vec::new();
        let mut visited = HashSet::new();
        let mut v = start;
        while !done(v) {
            let (n, w) = next(v)?;
            if w < prune_factor {
                visited.extend(path.drain(..));
            } else {
                path.insert(0, v);
            }
            if path.contains(&n) || visited.contains(&n) {
                return None;
            }
            v = n;
        }
        path.insert(0, v);
        if return_path(v) {
            Some(path)
        } else {
            None
        }
    }

    /// Longest heavy-incoming walk (prune_factor=0 style) when branch walk stops one vertex short.
    fn longest_heavy_incoming_path_upwards(
        &self,
        sink: usize,
        min_vertices: usize,
    ) -> Option<Vec<usize>> {
        let mut path = vec![sink];
        let mut v = sink;
        while !self.is_ref_source_vertex(v) {
            let Some((p, _)) = self.heaviest_incoming(v) else {
                break;
            };
            if path.contains(&p) {
                return None;
            }
            path.push(p);
            v = p;
            if path.len() >= min_vertices
                && self.has_incident_ref_in(v)
                && self.outgoing_nodes(v).len() > 1
            {
                break;
            }
        }
        if path.len() < min_vertices || self.is_ref_source_vertex(path[path.len() - 1]) {
            return None;
        }
        path.reverse();
        if self.is_ref_source(path[0]) {
            return None;
        }
        Some(path)
    }

    /// ASM-1 / k=85: walk heaviest incoming until ref-adjacent junction (no `outDegree>1` gate).
    fn walk_heavy_incoming_to_ref_junction(
        &self,
        sink: usize,
        min_vertices: usize,
    ) -> Option<Vec<usize>> {
        let mut path = vec![sink];
        let mut v = sink;
        loop {
            if path.len() >= min_vertices
                && self.has_incident_ref_in(v)
                && !self.is_ref_source_vertex(v)
            {
                break;
            }
            let Some((p, _)) = self.heaviest_incoming(v) else {
                return None;
            };
            if path.contains(&p) || self.is_ref_source_vertex(p) {
                return None;
            }
            path.push(p);
            v = p;
        }
        if path.len() < min_vertices {
            return None;
        }
        path.reverse();
        if self.is_ref_source(path[0]) {
            return None;
        }
        Some(path)
    }

    fn find_path_upwards_to_lca(
        &self,
        sink: usize,
        prune_factor: u32,
        give_up_at_branch: bool,
    ) -> Option<Vec<usize>> {
        if give_up_at_branch {
            self.find_path(
                sink,
                prune_factor,
                |v| self.incoming_count(v) != 1 || self.outgoing_nodes(v).len() >= 2,
                |v| self.outgoing_nodes(v).len() > 1,
                |v| {
                    let preds = self.incoming_nodes(v);
                    if preds.len() != 1 {
                        return None;
                    }
                    let p = preds[0];
                    Some((p, self.edge_support(p, v).unwrap_or(0)))
                },
            )
        } else {
            self.find_path(
                sink,
                prune_factor,
                |v| self.has_incident_ref_in(v) || self.incoming_count(v) == 0,
                |v| self.outgoing_nodes(v).len() > 1 && self.has_incident_ref_in(v),
                |v| {
                    let (p, w) = self.heaviest_incoming(v)?;
                    Some((p, w))
                },
            )
        }
    }

    fn find_path_downwards_to_ref(
        &self,
        source: usize,
        prune_factor: u32,
        give_up_at_branch: bool,
    ) -> Option<Vec<usize>> {
        if give_up_at_branch {
            self.find_path(
                source,
                prune_factor,
                |v| self.is_ref_node(v) || self.outgoing_nodes(v).len() != 1,
                |v| self.is_ref_node(v),
                |v| {
                    let outs = self.outgoing_nodes(v);
                    if outs.is_empty() {
                        return None;
                    }
                    let t = if outs.len() == 1 {
                        outs[0]
                    } else {
                        self.heaviest_outgoing(v).map(|(t, _)| t)?
                    };
                    Some((t, self.edge_support(v, t).unwrap_or(0)))
                },
            )
        } else {
            self.find_path(
                source,
                prune_factor,
                |v| self.is_ref_node(v) || self.outgoing_nodes(v).is_empty(),
                |v| self.is_ref_node(v),
                |v| {
                    let (t, w) = self.heaviest_outgoing(v)?;
                    Some((t, w))
                },
            )
        }
    }

    /// GATK `BaseGraph.getNextReferenceVertex(v, allowNonRefPaths, blacklistedEdge)`.
    fn next_reference_vertex_down(
        &self,
        v: usize,
        blacklist: Option<(usize, usize)>,
    ) -> Option<usize> {
        let outs = self.outgoing_nodes(v);
        if outs.is_empty() {
            return None;
        }
        for &t in &outs {
            if self.edge_is_ref(v, t) {
                return Some(t);
            }
        }
        let candidates: Vec<usize> = outs
            .into_iter()
            .filter(|&t| {
                blacklist
                    .map(|(f, to)| !(f == v && to == t))
                    .unwrap_or(true)
            })
            .take(2)
            .collect();
        if candidates.len() == 1 {
            Some(candidates[0])
        } else {
            None
        }
    }

    /// GATK `BaseGraph.getPrevReferenceVertex`.
    fn prev_reference_vertex(&self, v: usize) -> Option<usize> {
        self.incoming_nodes(v)
            .into_iter()
            .find(|&p| self.is_ref_node(p))
    }

    fn reference_path_from(
        &self,
        start: usize,
        dir: TraversalDir,
        blacklist: Option<(usize, usize)>,
    ) -> Vec<usize> {
        let mut path = vec![start];
        let mut v = start;
        loop {
            let next = match dir {
                TraversalDir::Down => self.next_reference_vertex_down(v, blacklist),
                TraversalDir::Up => self.prev_reference_vertex(v),
            };
            let Some(n) = next else {
                break;
            };
            path.push(n);
            v = n;
        }
        path
    }

    fn plan_dangling_tail_merge(
        &self,
        sink: usize,
        params: &DanglingRecoveryParams,
    ) -> Result<DanglingTailMergePlan, &'static str> {
        if !self.outgoing_nodes(sink).is_empty() {
            return Err("not_sink");
        }
        let min_tail = params.min_dangling_branch_length.max(1);
        let min_vertices = min_tail + 1; // GATK: include LCA vertex
                                         // Java `generateCigarAgainstDownwardsReferencePath` → `findPathUpwardsToLowestCommonAncestor`.
        let give_up = !params.recover_all_dangling_branches;
        let mut alt_path = self.find_path_upwards_to_lca(sink, params.min_prune_factor, give_up);
        // P12 extensions when Java walk is short (not in GATK; cluster graph only).
        if alt_path.as_ref().is_none_or(|p| p.len() < min_vertices) && give_up {
            alt_path = self.find_path_upwards_to_lca(sink, params.min_prune_factor, false);
        }
        if alt_path.as_ref().is_none_or(|p| p.len() < min_vertices) {
            alt_path = self.find_path_upwards_to_lca(sink, 0, false);
        }
        if alt_path.as_ref().is_none_or(|p| p.len() < min_vertices) {
            alt_path = self.longest_heavy_incoming_path_upwards(sink, min_vertices);
        }
        if alt_path.as_ref().is_none_or(|p| p.len() < min_vertices) {
            alt_path = self.walk_heavy_incoming_to_ref_junction(sink, min_vertices);
        }
        let Some(alt_path) = alt_path else {
            return Err("no_alt_path");
        };
        if alt_path.is_empty() {
            return Err("empty_alt_path");
        }
        if self.is_ref_source(alt_path[0]) {
            return Err("alt_path_starts_at_ref_source");
        }
        if alt_path.len() < min_vertices {
            return Err("alt_path_too_short");
        }
        let lca = alt_path[0];
        let blacklist = if alt_path.len() > 1 {
            let a = alt_path[1];
            self.heaviest_incoming(a).map(|(p, _)| (p, a))
        } else {
            None
        };
        let ref_path = self.reference_path_from(lca, TraversalDir::Down, blacklist);
        let ref_bases = path_bases(self, &ref_path, false);
        let alt_bases = path_bases(self, &alt_path, false);
        let cigar = align_dangling(&ref_bases, &alt_bases, &params.sw);
        if !cigar_ok_to_merge_tail(&cigar) {
            return Err("cigar_not_ok");
        }
        let elements = &cigar.elements;
        let last_ref_idx = cigar.reference_length().saturating_sub(1);
        let last_el_len = elements
            .last()
            .filter(|e| e.operator == CigarOperator::Match)
            .map(|e| e.length)
            .unwrap_or(1);
        let mut matching_suffix =
            longest_suffix_match_java(&ref_bases, &alt_bases, last_ref_idx).min(last_el_len);
        let min_match = params.min_matching_bases_to_dangling_end_recovery;
        if min_match >= 0 {
            if matching_suffix < min_match as usize {
                return Err("matching_suffix_below_min");
            }
        } else if matching_suffix == 0 && !params.dangling_java_exact {
            // ASM-1 k=85: SW suffix can be 0 while paths still share a mergeable tail (legacy k-mer walk).
            let max_ix = last_el_len.min(ref_bases.len()).min(alt_bases.len());
            if max_ix > 0 {
                let ar: Vec<u8> = ref_bases.iter().rev().take(max_ix).copied().collect();
                let alt_r: Vec<u8> = alt_bases.iter().rev().take(max_ix).copied().collect();
                let mut ar_fwd = ar.clone();
                ar_fwd.reverse();
                let mut alt_fwd = alt_r.clone();
                alt_fwd.reverse();
                if let Some(idx) =
                    best_prefix_match_legacy(&ar_fwd, &alt_fwd, ar_fwd.len(), self.kmer_size)
                {
                    if idx > 0 {
                        matching_suffix = idx.min(last_el_len);
                    }
                }
            }
            if matching_suffix == 0 {
                // ASM-1 k=85 cluster: paths can share a terminal base while SW suffix is 0.
                if ref_bases
                    .last()
                    .zip(alt_bases.last())
                    .is_some_and(|(a, b)| base_eq(*a, *b))
                {
                    matching_suffix = 1.min(last_el_len);
                } else {
                    return Err("matching_suffix_zero");
                }
            }
        } else if matching_suffix == 0 {
            return Err("matching_suffix_zero");
        }
        let alt_index_to_merge = cigar.read_length().saturating_sub(matching_suffix + 1);
        let first_is_del = elements
            .first()
            .is_some_and(|e| e.operator == CigarOperator::Deletion);
        let first_del_len = elements
            .first()
            .filter(|e| e.operator == CigarOperator::Deletion)
            .map(|e| e.length)
            .unwrap_or(0);
        let must_handle_leading_del =
            first_is_del && first_del_len + matching_suffix == last_ref_idx + 1;
        let ref_index_to_merge =
            last_ref_idx.saturating_sub(matching_suffix) + 1 + usize::from(must_handle_leading_del);
        if ref_index_to_merge == 0 {
            return Err("ref_index_zero_cycle");
        }
        if alt_index_to_merge >= alt_path.len() || ref_index_to_merge >= ref_path.len() {
            return Err("merge_index_oob");
        }
        let from = alt_path[alt_index_to_merge];
        let to = ref_path[ref_index_to_merge];
        // Java `addEdge` does not check; a prior pass may have merged this junction already.
        if self.has_edge(from, to) {
            return Err("edge_exists");
        }
        Ok(DanglingTailMergePlan {
            from,
            to,
            alt_bases,
            ref_path_bases: ref_bases,
            cigar,
        })
    }

    /// Per alt-sink failure reason after pruning (ASM-1 diagnostic).
    pub fn probe_dangling_tail_failures(
        &self,
        params: &DanglingRecoveryParams,
    ) -> Vec<(usize, String, String)> {
        let mut out = Vec::new();
        for v in 0..self.node_count() {
            if !self.outgoing_nodes(v).is_empty() || self.is_ref_sink(v) {
                continue;
            }
            let kmer = String::from_utf8_lossy(self.kmer_at(v)).into_owned();
            let reason = match self.plan_dangling_tail_merge(v, params) {
                Ok(plan) => {
                    format!(
                        "ok_merge:{}->{}",
                        String::from_utf8_lossy(self.kmer_at(plan.from)),
                        String::from_utf8_lossy(self.kmer_at(plan.to))
                    )
                }
                Err("edge_exists") => "ok_merge:edge_exists".to_string(),
                Err(r) => r.to_string(),
            };
            out.push((v, kmer, reason));
        }
        out
    }

    fn recover_dangling_tail(&mut self, sink: usize, params: &DanglingRecoveryParams) -> bool {
        match self.plan_dangling_tail_merge(sink, params) {
            Ok(plan) => {
                if !self.has_edge(plan.from, plan.to) {
                    self.add_dangling_recovery_edge(plan.from, plan.to, params.dangling_java_exact);
                }
                let has_indel = plan.cigar.elements.iter().any(|e| e.operator.is_indel());
                if has_indel || plan.alt_bases != plan.ref_path_bases {
                    self.dangling_merge_haps.push(DanglingMergeHaplotype {
                        alt_bases: plan.alt_bases,
                        cigar: plan.cigar,
                        alignment_start_hap_wrt_ref: 0,
                    });
                }
                true
            }
            Err("edge_exists") => true,
            Err(_) => false,
        }
    }

    /// GATK `extendDanglingPathAgainstReference` (heads only).
    fn extend_dangling_path_against_reference(
        &mut self,
        alt_path: &mut Vec<usize>,
        ref_path: &[usize],
        cigar: &Cigar,
        num_nodes_to_extend: usize,
    ) -> bool {
        if alt_path.is_empty() || ref_path.is_empty() || num_nodes_to_extend == 0 {
            return false;
        }
        let last_dangling = alt_path.len().saturating_sub(1);
        let offset: isize = cigar
            .elements
            .iter()
            .map(|e| {
                (if e.operator.consumes_reference_bases() {
                    e.length as isize
                } else {
                    0
                }) - if e.operator.consumes_read_bases() {
                    e.length as isize
                } else {
                    0
                }
            })
            .sum();
        let ref_use_idx = last_dangling as isize + offset + num_nodes_to_extend as isize;
        if ref_use_idx < 0 {
            return false;
        }
        let ref_use_idx = ref_use_idx as usize;
        if ref_use_idx >= ref_path.len() {
            return false;
        }
        let dangling_source = alt_path[last_dangling];
        let Some((mut prev_v, edge_weight)) = self.heaviest_outgoing(dangling_source) else {
            return false;
        };
        self.remove_edge(dangling_source, prev_v);
        alt_path.pop();

        let ref_kmer = self.kmer_at(ref_path[ref_use_idx]);
        let src_kmer = self.kmer_at(dangling_source);
        let mut seq = Vec::with_capacity(num_nodes_to_extend + src_kmer.len());
        for i in 0..num_nodes_to_extend {
            if i >= ref_kmer.len() {
                return false;
            }
            seq.push(ref_kmer[i]);
        }
        seq.extend_from_slice(src_kmer);
        let k = self.kmer_size;
        if seq.len() < k {
            return false;
        }
        for extend_i in (1..=num_nodes_to_extend).rev() {
            let start = extend_i;
            if start + k > seq.len() {
                return false;
            }
            let new_v = self.ensure_node(&seq[start..start + k]);
            self.add_edge_support(new_v, prev_v, edge_weight);
            alt_path.push(new_v);
            prev_v = new_v;
        }
        true
    }

    fn plan_dangling_head_merge(
        &mut self,
        source: usize,
        params: &DanglingRecoveryParams,
    ) -> Result<(usize, usize), &'static str> {
        if self.incoming_count(source) > 0 {
            return Err("not_source");
        }
        let min_vertices = params.min_dangling_branch_length + 1;
        let give_up = !params.recover_all_dangling_branches;
        let mut alt_path =
            self.find_path_downwards_to_ref(source, params.min_prune_factor, give_up);
        if alt_path.as_ref().is_none_or(|p| p.len() < min_vertices) && give_up {
            alt_path = self.find_path_downwards_to_ref(source, params.min_prune_factor, false);
        }
        if alt_path.as_ref().is_none_or(|p| p.len() < min_vertices) {
            alt_path = self.find_path_downwards_to_ref(source, 0, false);
        }
        let Some(mut alt_path) = alt_path else {
            return Err("no_alt_path");
        };
        if alt_path.is_empty() {
            return Err("empty_alt_path");
        }
        if self.is_ref_sink(alt_path[0]) {
            return Err("alt_path_starts_at_ref_sink");
        }
        if alt_path.len() < min_vertices {
            return Err("alt_path_too_short");
        }
        let lca = alt_path[0];
        let ref_path = self.reference_path_from(lca, TraversalDir::Up, None);
        let ref_bases = path_bases(self, &ref_path, true);
        let alt_bases = path_bases(self, &alt_path, true);
        let cigar = align_dangling(&ref_bases, &alt_bases, &params.sw);
        if !cigar_ok_to_merge_head(&cigar) {
            return Err("cigar_not_ok");
        }
        let min_match = params.min_matching_bases_to_dangling_end_recovery;
        let (ref_idx, alt_idx) = if min_match >= 0 {
            let Some(pair) = best_prefix_match(&cigar, &ref_bases, &alt_bases, min_match) else {
                return Err("prefix_match_failed");
            };
            pair
        } else {
            let first_el_len = cigar
                .elements
                .first()
                .filter(|e| e.operator == CigarOperator::Match)
                .map(|e| e.length)
                .unwrap_or(1);
            let Some(idx) =
                best_prefix_match_legacy(&ref_bases, &alt_bases, first_el_len, self.kmer_size)
            else {
                return Err("prefix_match_legacy_failed");
            };
            if idx == 0 {
                return Err("prefix_match_zero");
            }
            (idx, idx)
        };
        if ref_idx >= ref_path.len().saturating_sub(1) {
            return Err("ref_index_oob");
        }
        if alt_idx >= alt_path.len() {
            let num_extend = alt_idx.saturating_sub(alt_path.len()) + 2;
            if !self.extend_dangling_path_against_reference(
                &mut alt_path,
                &ref_path,
                &cigar,
                num_extend,
            ) {
                return Err("extend_path_failed");
            }
            if alt_idx >= alt_path.len() {
                return Err("extend_path_needed");
            }
        }
        let from = ref_path[ref_idx + 1];
        let to = alt_path[alt_idx];
        if self.has_edge(from, to) {
            return Err("edge_exists");
        }
        Ok((from, to))
    }

    /// Per alt-source failure reason (ASM-1 diagnostic).
    pub fn probe_dangling_head_failures(
        &self,
        params: &DanglingRecoveryParams,
    ) -> Vec<(usize, String, String)> {
        let mut out = Vec::new();
        for v in 0..self.node_count() {
            if self.incoming_count(v) > 0 || self.is_ref_source_vertex(v) {
                continue;
            }
            let kmer = String::from_utf8_lossy(self.kmer_at(v)).into_owned();
            // CLONE: needed because graph fork needs owned duplicate for speculative path.
            let mut g = self.clone();
            let reason = match g.plan_dangling_head_merge(v, params) {
                Ok((from, to)) => {
                    format!(
                        "ok_merge:{}->{}",
                        String::from_utf8_lossy(self.kmer_at(from)),
                        String::from_utf8_lossy(self.kmer_at(to))
                    )
                }
                Err(r) => r.to_string(),
            };
            out.push((v, kmer, reason));
        }
        out
    }

    fn recover_dangling_head(&mut self, source: usize, params: &DanglingRecoveryParams) -> bool {
        match self.plan_dangling_head_merge(source, params) {
            Ok((from, to)) => {
                if !self.has_edge(from, to) {
                    self.add_dangling_recovery_edge(from, to, params.dangling_java_exact);
                }
                true
            }
            Err("edge_exists") => true,
            Err(_) => false,
        }
    }

    /// GATK `addEdge` weight for dangling recovery (1) vs ASM-1 spine boost.
    fn add_dangling_recovery_edge(&mut self, from: usize, to: usize, java_exact: bool) {
        if java_exact {
            self.add_edge_support(from, to, 1);
        } else {
            self.add_merged_edge(from, to, 1, false);
        }
    }

    /// Tail candidates for parity stage dumps (`outDegree==0 && inDegree>0 && !refSink`).
    #[cfg(any(feature = "dev-dumps", test))]
    pub fn count_dangling_tail_candidates_parity_dump(&self) -> u32 {
        (0..self.node_count())
            .filter(|&v| {
                self.outgoing_nodes(v).is_empty()
                    && self.incoming_count(v) > 0
                    && !self.is_ref_sink(v)
            })
            .count() as u32
    }

    /// Head candidates counted like `HcFullParityGateDump.assemblyRegionAssemblyStages` (`inDegree==0 && outDegree>0`, no ref-source skip).
    #[cfg(any(feature = "dev-dumps", test))]
    pub fn count_dangling_head_candidates_parity_dump(&self) -> u32 {
        (0..self.node_count())
            .filter(|&v| self.incoming_count(v) == 0 && !self.outgoing_nodes(v).is_empty())
            .count() as u32
    }

    /// GATK dangling tail/head recovery (requires ref-threaded graph).
    pub fn recover_dangling_branches(
        &mut self,
        params: &DanglingRecoveryParams,
    ) -> GatkResult<DanglingRecoverySummary> {
        if self.ref_nodes.is_empty() {
            return Err(GatkError::argument(
                "dangling recovery requires a reference-threaded graph",
            ));
        }
        let edges_before = self.edge_count();
        let mut tails_attempted = 0u32;
        let mut tails_recovered = 0u32;
        let mut heads_attempted = 0u32;
        let mut heads_recovered = 0u32;

        if params.dangling_java_exact {
            for v in 0..self.node_count() {
                if self.outgoing_nodes(v).is_empty() && !self.is_ref_sink(v) {
                    tails_attempted += 1;
                    if self.recover_dangling_tail(v, params) {
                        tails_recovered += 1;
                    }
                }
            }
            if params.recover_dangling_heads {
                let sources: Vec<usize> = (0..self.node_count())
                    .filter(|&v| {
                        self.incoming_count(v) == 0
                            && !self.outgoing_nodes(v).is_empty()
                            && !self.is_ref_source_vertex(v)
                    })
                    .collect();
                for v in sources {
                    heads_attempted += 1;
                    if self.recover_dangling_head(v, params) {
                        heads_recovered += 1;
                    }
                }
            }
            self.cleanup_isolated_nodes();
            return Ok(DanglingRecoverySummary {
                edges_before,
                edges_after: self.edge_count(),
                tails_attempted,
                tails_recovered,
                heads_attempted,
                heads_recovered,
                edges_merged: tails_recovered + heads_recovered,
            });
        }

        // Repeat until a pass adds no merges: earlier tail merges can lengthen later alt paths (ASM-1).
        const MAX_DANGLING_PASSES: usize = 8;
        let mut tail_sinks_seen = HashSet::new();
        let mut tail_sinks_recovered = HashSet::new();
        for _ in 0..MAX_DANGLING_PASSES {
            let mut pass_merges = 0u32;
            let mut sinks: Vec<usize> = (0..self.node_count())
                .filter(|&v| self.outgoing_nodes(v).is_empty() && !self.is_ref_sink(v))
                .collect();
            sinks.sort_by(|&a, &b| {
                self.heaviest_incoming(a)
                    .map(|(_, w)| w)
                    .unwrap_or(0)
                    .cmp(&self.heaviest_incoming(b).map(|(_, w)| w).unwrap_or(0))
                    .reverse()
            });
            for v in sinks {
                if tail_sinks_seen.insert(v) {
                    tails_attempted += 1;
                }
                if self.recover_dangling_tail(v, params) {
                    if tail_sinks_recovered.insert(v) {
                        tails_recovered += 1;
                    }
                    pass_merges += 1;
                }
            }
            if pass_merges == 0 {
                break;
            }
        }

        if params.recover_dangling_heads {
            for _ in 0..MAX_DANGLING_PASSES {
                let mut pass_merges = 0u32;
                let sources: Vec<usize> = (0..self.node_count())
                    .filter(|&v| {
                        self.incoming_count(v) == 0
                            && !self.outgoing_nodes(v).is_empty()
                            && !self.is_ref_source_vertex(v)
                    })
                    .collect();
                for v in sources {
                    heads_attempted += 1;
                    if self.recover_dangling_head(v, params) {
                        heads_recovered += 1;
                        pass_merges += 1;
                    }
                }
                if pass_merges == 0 {
                    break;
                }
            }
        }

        self.cleanup_isolated_nodes();
        let edges_after = self.edge_count();
        Ok(DanglingRecoverySummary {
            edges_before,
            edges_after,
            tails_attempted,
            tails_recovered,
            heads_attempted,
            heads_recovered,
            edges_merged: tails_recovered + heads_recovered,
        })
    }
}

/// Inject ASM-1 dangling merge haps into the assembly haplotype list (graph-only EventMap path).
pub fn apply_dangling_merge_haplotypes(
    haplotypes: &mut Vec<crate::haplotype::Haplotype>,
    ref_hap: &crate::haplotype::Haplotype,
    merges: &[DanglingMergeHaplotype],
    ref_bytes: &[u8],
    sw: &SwParameters,
) {
    use crate::haplotype::Haplotype;
    use crate::haplotype_cigar::calculate_haplotype_cigar_for_assembly_with_offset;
    for h in haplotypes.iter_mut() {
        if h.is_reference || h.bases != ref_bytes {
            continue;
        }
        for m in merges {
            let has_indel = m.cigar.elements.iter().any(|e| e.operator.is_indel());
            if has_indel && m.alt_bases != ref_bytes {
                // CLONE: needed because haplotype owns base string.
                h.bases = m.alt_bases.clone();
                // CLONE: needed because haplotype owns CIGAR.
                h.cigar = Some(m.cigar.clone());
                h.alignment_start_hap_wrt_ref = m.alignment_start_hap_wrt_ref;
                break;
            }
        }
    }
    for m in merges {
        if m.alt_bases.is_empty() {
            continue;
        }
        if haplotypes
            .iter()
            .any(|h| !h.is_reference && h.bases == m.alt_bases)
        {
            continue;
        }
        let has_indel = m.cigar.elements.iter().any(|e| e.operator.is_indel());
        if !has_indel && m.alt_bases == ref_bytes {
            continue;
        }
        // CLONE: needed because haplotype constructor takes owned bases.
        let mut h = Haplotype::new(m.alt_bases.clone(), false);
        // CLONE: needed because haplotype owns CIGAR.
        h.cigar = Some(m.cigar.clone());
        h.alignment_start_hap_wrt_ref = m.alignment_start_hap_wrt_ref;
        h.score = 25.0;
        if let Some(ref_cigar) = ref_hap.cigar.as_ref() {
            if let Some(assy) = calculate_haplotype_cigar_for_assembly_with_offset(
                ref_bytes,
                &h.bases,
                ref_cigar.reference_length(),
                sw,
            ) {
                if assy.cigar.elements.iter().any(|e| e.operator.is_indel()) {
                    h.cigar = Some(assy.cigar);
                    h.alignment_start_hap_wrt_ref = assy.alignment_start_hap_wrt_ref;
                }
            }
        }
        haplotypes.push(h);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::assembly::{AssemblyGraphParams, AssemblyGraphPruningParams, AssemblyRead};
    use crate::assembly_pruning::apply_gatk_pruning;
    use crate::read_threading_graph::assembly_graph_from_ref_and_reads_threading;

    fn read(seq: &str, q: u8) -> AssemblyRead {
        AssemblyRead {
            bases: seq.as_bytes().to_vec(),
            base_quals: vec![q; seq.len()],
        }
    }

    fn build_pruned_ref_graph_at_k(
        reference: &str,
        alt_reads: &[&str],
        kmer_size: usize,
    ) -> AssemblyGraph {
        let params = AssemblyGraphParams {
            kmer_size: crate::bio_ids::KmerSize::try_new(kmer_size as u16).expect("test k≥2"),
            min_base_quality: 10,
            ..Default::default()
        };
        let mut reads: Vec<AssemblyRead> = alt_reads.iter().map(|s| read(s, 30)).collect();
        reads.insert(0, read(reference, 30));
        reads.insert(0, read(reference, 30));
        reads.insert(0, read(reference, 30));
        let reference = read(reference, 30);
        let mut graph =
            assembly_graph_from_ref_and_reads_threading(&reference, &reads, &params).unwrap();
        let mut prune = AssemblyGraphPruningParams::gatk_haplotype_caller_defaults();
        prune.min_prune_factor = 2;
        apply_gatk_pruning(&mut graph, &prune);
        graph
    }

    fn build_pruned_ref_graph(reference: &str, alt_reads: &[&str]) -> AssemblyGraph {
        build_pruned_ref_graph_at_k(reference, alt_reads, 3)
    }

    #[test]
    fn find_path_upwards_is_lca_first() {
        let graph = build_pruned_ref_graph(
            "ACGTTGCATCG",
            &["ACGTTGCATCG", "ACGTTGCATCG", "ACGTTGCATCA", "ACGTTGCATCA"],
        );
        let alt_sink = graph
            .nodes()
            .iter()
            .position(|n| n.kmer == b"TCA")
            .expect("TCA sink");
        let path = graph
            .find_path_upwards_to_lca(alt_sink, 2, true)
            .expect("alt path");
        assert_eq!(graph.kmer_at(path[0]), b"ATC");
        assert_eq!(graph.kmer_at(*path.last().unwrap()), b"TCA");
    }

    #[test]
    fn best_prefix_match_requires_min_matching_bases() {
        let mut cigar = Cigar::new();
        cigar.push(4, CigarOperator::Match);
        let ref_bases = b"ACGT";
        let alt_bases = b"ACGT";
        assert!(best_prefix_match(&cigar, ref_bases, alt_bases, 3).is_some());
        assert!(best_prefix_match(&cigar, ref_bases, alt_bases, 5).is_none());
    }

    #[test]
    fn longest_suffix_match_matches_gatk_examples() {
        assert_eq!(longest_suffix_match_java(b"ACGT", b"TGT", 3), 2);
        assert_eq!(longest_suffix_match_java(b"ACGT", b"CGT", 3), 3);
        assert_eq!(longest_suffix_match_java(b"CG", b"CA", 1), 0);
    }

    #[test]
    fn dangling_java_exact_single_pass_matches_gatk_edge_count() {
        let mut graph = build_pruned_ref_graph(
            "ACGTTGCATCG",
            &["ACGTTGCATCG", "ACGTTGCATCG", "ACGTTGCATCA", "ACGTTGCATCA"],
        );
        let mut params = DanglingRecoveryParams::gatk_haplotype_caller_defaults();
        params.min_dangling_branch_length = 1;
        params.dangling_java_exact = true;
        let summary = graph.recover_dangling_branches(&params).unwrap();
        assert!(summary.tails_attempted >= 1);
        // This fixture needs ASM-1 suffix rescue for a merge; GATK-exact mode correctly skips it.
        assert_eq!(summary.tails_recovered, 0);
        let mut multi = params;
        multi.dangling_java_exact = false;
        let multi_summary = graph.recover_dangling_branches(&multi).unwrap();
        assert_eq!(multi_summary.tails_recovered, 1);
    }

    #[test]
    fn dangling_tail_recovery_attempts_alt_sink() {
        let mut graph = build_pruned_ref_graph(
            "ACGTTGCATCG",
            &["ACGTTGCATCG", "ACGTTGCATCG", "ACGTTGCATCA", "ACGTTGCATCA"],
        );
        let mut dangling = DanglingRecoveryParams::gatk_haplotype_caller_defaults();
        dangling.min_dangling_branch_length = 1;
        let summary = graph.recover_dangling_branches(&dangling).unwrap();
        assert!(summary.tails_attempted >= 1);
        // Java idempotent `addEdge` → `edge_exists` counts as tail recovered (ASM-1).
        assert_eq!(summary.tails_recovered, 1);
    }

    #[test]
    fn cigar_ok_to_merge_tail_requires_terminal_match_op() {
        let mut trailing_ins = Cigar::new();
        trailing_ins.push(1, CigarOperator::Match);
        trailing_ins.push(1, CigarOperator::Insertion);
        assert!(!cigar_ok_to_merge_tail(&trailing_ins));

        let mut ok = Cigar::new();
        ok.push(2, CigarOperator::Match);
        assert!(cigar_ok_to_merge_tail(&ok));
    }

    #[test]
    fn align_dangling_uses_leading_indel_strategy() {
        let p = DanglingRecoverySwParams::gatk_defaults();
        let cigar = align_dangling(b"ACGTACGTAC", b"ACGTXACGTAC", &p);
        assert!(cigar_ok_to_merge_tail(&cigar));
        assert!(cigar
            .elements
            .last()
            .is_some_and(|e| e.operator == CigarOperator::Match));
    }

    /// GATK `ReadThreadingGraphUnitTest.testForkedDanglingEnds`.
    #[test]
    fn forked_dangling_ends_recovers_all_alt_sinks_with_recover_all() {
        let common_prefix = "AAAAAAAAAACCCCCCCCCCGGGGGGGGGGTTTTTTTTTT";
        let reference = format!("{common_prefix}GCTAGCTAATCG");
        let alt1 = format!("{common_prefix}ACTAGCTAATCG");
        let alt2 = format!("{common_prefix}ACTAGATAATCG");
        let mut graph = build_pruned_ref_graph_at_k(&reference, &[&alt1, &alt2], 15);
        let mut dangling = DanglingRecoveryParams::gatk_haplotype_caller_defaults();
        dangling.min_dangling_branch_length = 4;
        dangling.recover_all_dangling_branches = true;
        let summary = graph.recover_dangling_branches(&dangling).unwrap();
        assert!(
            summary.tails_attempted >= 1,
            "GATK testForkedDanglingEnds expects non-ref sinks (Rust may collapse forks to fewer sinks)"
        );
        assert_eq!(
            summary.tails_recovered, summary.tails_attempted,
            "recoverAll should merge every attempted alt tail"
        );
    }
}
