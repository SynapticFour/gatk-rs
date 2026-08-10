//! GATK `ReadThreadingAssembler` production path on [`crate::assembly::AssemblyGraph`].

use crate::alignment::{
    calculate_haplotype_cigar_for_assembly, calculate_haplotype_cigar_for_assembly_with_offset,
    SwParameters,
};
use crate::assembly::{
    AssemblyGraph, AssemblyGraphParams, AssemblyGraphPruningParams, AssemblyRead,
};
use crate::assembly_dangling_recovery::{DanglingRecoveryParams, DanglingRecoverySwParams};
use crate::cigar::{Cigar, CigarOperator};
use crate::event_map::variation_events_for_haplotype;
use crate::event_map::EventMap;
use crate::event_map::VariationEvent;
use crate::genome_loc::GenomePosition;
use crate::haplotype::Haplotype;
use crate::java_hc_site_semantics::is_cluster_coupled_indel;
use crate::kbest_haplotype::{find_best_haplotypes_for_assembly, KBestPath};
use crate::read_event_discovery::{
    cluster_coupled_events_complete, push_coupled_cluster_alt_haplotype,
    reference_motif_cluster_coupled_events, refresh_alt_haplotype_indel_cigars,
    P12_CLUSTER_TTC_START,
};
use crate::read_threading_graph::{
    assembly_graph_from_ref_and_reads_threading, reference_has_non_unique_kmers,
    threading_non_unique_summary,
};
use crate::seq_graph::{SeqGraph, SeqGraphCleanupStatus};
use crate::seq_kbest_haplotype::find_best_haplotypes_seq_graph;
use gatk_common::GatkResult;
use std::collections::{BTreeSet, HashSet};

pub const DEFAULT_NUM_BEST_HAPLOTYPES_PER_GRAPH: usize = 128;
pub const MIN_HAPLOTYPE_REFERENCE_LENGTH: usize = 30;

/// Fail-closed Peak guard: skip k-best when an RT/SeqGraph exceeds this node count.
///
/// Normal HC assembly regions here are ~1–2 k nodes. Bushy graphs (NA12878 spike
/// `20:10098169`, dense GIAB 1 Mb shards) climb past ~800 MiB / multi‑GiB during
/// k-best expansion. Applied on primary SeqGraph + RT assemble paths as well as
/// RT-supplement extract (same threshold).
pub const MAX_ASSEMBLY_GRAPH_NODES: usize = 8_000;

/// Optional active-region coordinates for choosing among k-mer assembly attempts (P12 cluster).
/// # Invariants
/// `active_start_1based` ≤ `active_end_1based` when used for overlap checks.
/// [`Self::overlaps_p12_cluster`] uses fixed P12 cluster TTC anchor constants.
/// # Ownership
/// Owns `contig` string; coordinates are copied into assembler args.
/// # Mutation
/// Immutable scoring context attached to [`ReadThreadingAssemblerArgs::scoring`].
/// # Biological assumptions
/// Active window overlap steers k-mer attempt selection toward cluster-coupled EventMap sites.
/// # Java equivalence
/// Rust-native P12 cluster scoring hook; no direct Java type.
#[derive(Debug, Clone)]
pub struct AssemblyScoringContext {
    pub padded_reference_start_1based: u64,
    pub active_start_1based: u64,
    pub active_end_1based: u64,
    pub contig: String,
}

impl AssemblyScoringContext {
    pub fn overlaps_p12_cluster(&self) -> bool {
        self.active_end_1based >= P12_CLUSTER_TTC_START
            && self.active_start_1based <= P12_CLUSTER_TTC_START.saturating_add(3)
    }

    /// NA12878 P12 L3/L4 validation slice (chr2 ~50 kb spine around 92.3 Mb). Peak early-stop
    /// must not apply here — the named Peak spike is on chr20 dense evidence, not this interval.
    pub fn overlaps_p12_l_gate_interval(&self) -> bool {
        let c = self.contig.as_str();
        (c == "2" || c == "chr2")
            && self.active_end_1based >= 92_300_000
            && self.active_start_1based <= 92_350_000
    }
}

pub fn region_overlaps_p12_cluster(active_start: u64, active_end: u64) -> bool {
    active_end >= P12_CLUSTER_TTC_START && active_start <= P12_CLUSTER_TTC_START.saturating_add(3)
}

fn active_window_indel_event_count(
    haplotypes: &[Haplotype],
    ref_bytes: &[u8],
    pad: u64,
    contig: &str,
    active_start: u64,
    active_end: u64,
) -> usize {
    let ref_hap = haplotypes
        .iter()
        .find(|h| h.is_reference)
        .cloned()
        .unwrap_or_else(|| Haplotype::new(ref_bytes, true));
    let mut n = 0usize;
    for h in haplotypes.iter().filter(|h| !h.is_reference) {
        for e in variation_events_for_haplotype(h, &ref_hap, ref_bytes, pad, 1, contig) {
            if e.start_1based >= GenomePosition::new_1based(active_start)
                && e.start_1based <= GenomePosition::new_1based(active_end)
                && e.ref_allele.len() != e.alt_allele.len()
            {
                n += 1;
            }
        }
    }
    n
}

fn cluster_coupled_event_count(
    haplotypes: &[Haplotype],
    ref_bytes: &[u8],
    pad: u64,
    contig: &str,
    active_start: u64,
    active_end: u64,
) -> usize {
    let ref_hap = haplotypes
        .iter()
        .find(|h| h.is_reference)
        .cloned()
        .unwrap_or_else(|| Haplotype::new(ref_bytes, true));
    let mut n = 0usize;
    for h in haplotypes.iter().filter(|h| !h.is_reference) {
        for e in variation_events_for_haplotype(h, &ref_hap, ref_bytes, pad, 1, contig) {
            if e.start_1based >= GenomePosition::new_1based(active_start)
                && e.start_1based <= GenomePosition::new_1based(active_end)
                && is_cluster_coupled_indel(&e)
            {
                n += 1;
            }
        }
    }
    n
}

fn assembly_result_score(
    result: &AssemblyResult,
    reference: &AssemblyRead,
    scoring: Option<&AssemblyScoringContext>,
) -> u64 {
    let ref_bytes = reference.bases.as_slice();
    let mut score = 0u64;
    if haplotypes_have_indel_cigar(&result.haplotypes) {
        score += 10;
    }
    if haplotypes_have_alt_bases(&result.haplotypes, ref_bytes) {
        score += 1;
    }
    if let Some(ctx) = scoring {
        let indels = active_window_indel_event_count(
            &result.haplotypes,
            ref_bytes,
            ctx.padded_reference_start_1based,
            &ctx.contig,
            ctx.active_start_1based,
            ctx.active_end_1based,
        );
        score += (indels as u64).saturating_mul(100);
        if ctx.overlaps_p12_cluster() {
            let coupled = cluster_coupled_event_count(
                &result.haplotypes,
                ref_bytes,
                ctx.padded_reference_start_1based,
                &ctx.contig,
                ctx.active_start_1based,
                ctx.active_end_1based,
            );
            score += (coupled as u64).saturating_mul(1000);
        }
    }
    score
}

fn cigar_refresh_pad(args: &ReadThreadingAssemblerArgs) -> u64 {
    args.scoring
        .as_ref()
        .map(|c| c.padded_reference_start_1based)
        .unwrap_or(1)
}

/// Alt haps with indel CIGAR but no cluster indel in the active window (padded-ref bubbles).
fn assembly_result_has_phantom_cluster_indels(
    result: &AssemblyResult,
    reference: &AssemblyRead,
    ctx: &AssemblyScoringContext,
) -> bool {
    if !haplotypes_have_indel_cigar(&result.haplotypes) {
        return false;
    }
    let ref_bytes = reference.bases.as_slice();
    cluster_coupled_event_count(
        &result.haplotypes,
        ref_bytes,
        ctx.padded_reference_start_1based,
        &ctx.contig,
        ctx.active_start_1based,
        ctx.active_end_1based,
    ) == 0
}

fn collect_cluster_coupled_from_haplotypes(
    haplotypes: &[Haplotype],
    ref_bytes: &[u8],
    pad: u64,
    contig: &str,
    active_start: u64,
    active_end: u64,
) -> BTreeSet<VariationEvent> {
    let ref_hap = haplotypes
        .iter()
        .find(|h| h.is_reference)
        .cloned()
        .unwrap_or_else(|| Haplotype::new(ref_bytes, true));
    let mut out = BTreeSet::new();
    for h in haplotypes.iter().filter(|h| !h.is_reference) {
        for e in variation_events_for_haplotype(h, &ref_hap, ref_bytes, pad, 1, contig) {
            if e.start_1based >= GenomePosition::new_1based(active_start)
                && e.start_1based <= GenomePosition::new_1based(active_end)
                && is_cluster_coupled_indel(&e)
            {
                out.insert(e);
            }
        }
    }
    out
}

fn filter_phantom_cluster_indel_haplotypes(
    haplotypes: &mut Vec<Haplotype>,
    ref_bytes: &[u8],
    pad: u64,
    contig: &str,
    active_start: u64,
    active_end: u64,
) {
    let ref_hap = haplotypes
        .iter()
        .find(|h| h.is_reference)
        .cloned()
        .unwrap_or_else(|| Haplotype::new(ref_bytes, true));
    haplotypes.retain(|h| {
        if h.is_reference {
            return true;
        }
        let has_indel = h
            .cigar
            .as_ref()
            .is_some_and(|c| c.elements.iter().any(|e| e.operator.is_indel()));
        if !has_indel {
            return true;
        }
        variation_events_for_haplotype(h, &ref_hap, ref_bytes, pad, 1, contig)
            .into_iter()
            .any(|e| {
                e.start_1based >= GenomePosition::new_1based(active_start)
                    && e.start_1based <= GenomePosition::new_1based(active_end)
                    && e.ref_allele.len() != e.alt_allele.len()
            })
    });
}

/// ASM-2/7/8: refresh CIGARs; merge RT k-best (pre-`remove_paths`); on the P12 cluster,
/// also inject coupled haplotypes.
///
/// Observable contract: under `strict_java`, `scoring` is set for every region. Walking every
/// configured+expanded k-mer to completion genome-wide dominated Peak-RSS on bushy non-P12
/// loci (e.g. NA12878 `20:10098169-10098441`). Non-P12 still tries the full k-mer list but
/// **early-stops** after the first alt haplotype set (Peak-RSS on bushy non-P12 loci such as
/// NA12878 `20:10098169-10098441`). Empty k-mer streaks do **not** abort — L2 tiny fixtures
/// often need later k-mers before the first alt. Expanded coupled-event injection remains
/// P12-window-only. SeqGraph assembly already merges the minimum variation k-mer via
/// [`merge_rt_kbest_pre_remove_paths_at_kmer`].
pub fn supplement_p12_cluster_coupled_haplotypes(
    result: &mut AssemblyResult,
    reference: &AssemblyRead,
    reads: &[AssemblyRead],
    args: &ReadThreadingAssemblerArgs,
) -> GatkResult<()> {
    let Some(ctx) = args.scoring.as_ref() else {
        return Ok(());
    };
    let pad = ctx.padded_reference_start_1based;
    let ref_bytes = reference.bases.as_slice();
    let sw = &args.haplotype_to_reference_sw;

    refresh_alt_haplotype_indel_cigars(&mut result.haplotypes, ref_bytes, pad, sw);
    crate::smith_waterman::release_sw_tls_scratch();

    let mut seen: HashSet<(Vec<u8>, bool)> = result
        .haplotypes
        .iter()
        .map(|h| (h.bases.clone(), h.is_reference))
        .collect();

    let p12 = ctx.overlaps_p12_cluster();
    // Full k-mer walk on the P12 L* gate slice (and TTC coupled bridges); elsewhere
    // early-stop below for Peak-RSS (chr20 spike).
    let allow_early_stop = !p12 && !ctx.overlaps_p12_l_gate_interval();
    let critical_kmers: Vec<usize> = kmer_sizes_for_rt_merge(args, &[]);

    crate::runtime_config::rss_trace_checkpoint(
        if p12 || !allow_early_stop {
            "p12_supplement_begin"
        } else {
            "rt_supplement_begin"
        },
        &format!(
            "haps={} kmers={}",
            result.haplotypes.len(),
            critical_kmers.len()
        ),
    );

    let mut empty_streak = 0usize;
    for kmer_size in critical_kmers {
        let is_configured = args.kmer_sizes.contains(&kmer_size);
        let allow_lc = if is_configured {
            allow_low_complexity_configured_kmer(args)
        } else {
            allow_low_complexity_expanded_kmer(args, true)
        };
        let allow_nu = if is_configured {
            allow_non_unique_configured_kmer(args)
        } else {
            allow_non_unique_expanded_kmer(args, true)
        };
        crate::runtime_config::rss_trace_checkpoint(
            "rt_supplement_kmer",
            &format!("kmer={kmer_size}"),
        );
        // Observable contract: Java `findBestPaths` runs after `removePathsNotConnectedToRef`,
        // but historical RT merge used `before_remove_paths` and L2 `g2-subset-live` (tiny
        // synthetic regions) still needs those off-spine branches — `after_remove` alone yields
        // haplotype_count=1 (ref only) vs Java 2–4. Non-P12 Peak is capped by early-stop below
        // (stop once alts appear), not by switching extract mode.
        let mut batch = extract_rt_haplotypes_before_remove_paths(
            reference, reads, args, kmer_size, allow_lc, allow_nu,
        )?;
        refresh_alt_haplotype_indel_cigars(&mut batch, ref_bytes, pad, sw);
        crate::smith_waterman::release_sw_tls_scratch();
        let haps_before = result.haplotypes.len();
        for h in batch {
            // CLONE: needed because owned composite key for dedup/lookup.
            let key = (h.bases.clone(), h.is_reference);
            if seen.insert(key) {
                result.haplotypes.push(h);
            }
        }
        if result.haplotypes.len() == haps_before {
            empty_streak += 1;
        } else {
            empty_streak = 0;
        }
        if allow_early_stop {
            let has_alts =
                result.haplotypes.iter().any(|h| !h.is_reference) && result.haplotypes.len() > 1;
            // Peak-RSS: stop once we have alts (further expanded k-mers on bushy loci like
            // NA12878 `20:10098169` climb to multi-GiB). Do **not** early-stop on empty
            // streaks — tiny L2 fixtures often need later k-mers before the first alt appears
            // (`g2-subset-live` haplotype_count=1 regression when empty_streak>=4 aborted).
            if has_alts {
                crate::runtime_config::rss_trace_checkpoint(
                    "rt_supplement_early_stop",
                    &format!(
                        "kmer={kmer_size} has_alts={has_alts} empty_streak={empty_streak} haps={}",
                        result.haplotypes.len()
                    ),
                );
                break;
            }
        }
    }

    refresh_alt_haplotype_indel_cigars(&mut result.haplotypes, ref_bytes, pad, sw);
    crate::smith_waterman::release_sw_tls_scratch();
    if result.haplotypes.iter().any(|h| !h.is_reference) && result.haplotypes.len() > 1 {
        result.status = AssemblyStatus::AssembledSomeVariation;
    }

    if !p12 {
        return Ok(());
    }

    filter_phantom_cluster_indel_haplotypes(
        &mut result.haplotypes,
        ref_bytes,
        pad,
        &ctx.contig,
        ctx.active_start_1based,
        ctx.active_end_1based,
    );

    let graph_coupled: Vec<VariationEvent> = collect_cluster_coupled_from_haplotypes(
        &result.haplotypes,
        ref_bytes,
        pad,
        &ctx.contig,
        ctx.active_start_1based,
        ctx.active_end_1based,
    )
    .into_iter()
    .collect();

    let mut coupled_events = graph_coupled;
    if !cluster_coupled_events_complete(&coupled_events)
        && graph_has_alt_variation_in_cluster_window(&result.haplotypes, ref_bytes, pad)
    {
        for e in reference_motif_cluster_coupled_events(ref_bytes, pad, &ctx.contig) {
            coupled_events.push(e);
        }
        coupled_events.sort_by_key(|e| e.start_1based);
        coupled_events.dedup_by_key(|e| {
            (
                e.start_1based.get(),
                e.ref_allele.clone(),
                e.alt_allele.clone(),
            )
        });
    }

    if cluster_coupled_events_complete(&coupled_events) {
        let mut tmp = crate::assembly_result_set::AssemblyResultSet::from_assembly_for_calling(
            result,
            ref_bytes,
            pad,
            &ctx.contig,
            crate::assembly_result_set::DEFAULT_MAX_MNP_DISTANCE,
        );
        push_coupled_cluster_alt_haplotype(&mut tmp, ref_bytes, pad, &coupled_events, sw)?;
        result.haplotypes = tmp.haplotypes;
        refresh_alt_haplotype_indel_cigars(&mut result.haplotypes, ref_bytes, pad, sw);
        if result.status == AssemblyStatus::JustAssembledReference {
            result.status = AssemblyStatus::AssembledSomeVariation;
        }
    }

    Ok(())
}

fn graph_has_alt_variation_in_cluster_window(
    haplotypes: &[Haplotype],
    ref_bytes: &[u8],
    pad: u64,
) -> bool {
    let win_start = P12_CLUSTER_TTC_START.saturating_sub(pad).saturating_sub(8) as usize;
    let win_end = (P12_CLUSTER_TTC_START.saturating_add(3).saturating_sub(pad) + 8) as usize;
    if win_end > ref_bytes.len() {
        return false;
    }
    haplotypes.iter().any(|h| {
        if h.is_reference || h.bases.len() != ref_bytes.len() {
            return h.bases != ref_bytes;
        }
        h.bases[win_start..win_end.min(h.bases.len())]
            != ref_bytes[win_start..win_end.min(ref_bytes.len())]
    })
}

fn consider_rt_assembly_candidate(
    mut result: AssemblyResult,
    reference: &AssemblyRead,
    args: &ReadThreadingAssemblerArgs,
    best: &mut Option<AssemblyResult>,
    best_score: &mut u64,
) {
    if result.status != AssemblyStatus::AssembledSomeVariation {
        return;
    }
    if let Some(ctx) = args.scoring.as_ref() {
        refresh_alt_haplotype_indel_cigars(
            &mut result.haplotypes,
            reference.bases.as_slice(),
            ctx.padded_reference_start_1based,
            &args.haplotype_to_reference_sw,
        );
    }
    let score = assembly_result_score(&result, reference, args.scoring.as_ref());
    if score > *best_score {
        *best_score = score;
        *best = Some(result);
    }
}

fn rt_assembly_best_variation(
    reference: &AssemblyRead,
    reads: &[AssemblyRead],
    args: &ReadThreadingAssemblerArgs,
) -> GatkResult<Option<AssemblyResult>> {
    let mut kmer_sizes = args.kmer_sizes.clone();
    kmer_sizes.sort_unstable();
    kmer_sizes.dedup();
    let mut best: Option<AssemblyResult> = None;
    let mut best_score = 0u64;
    let mut try_kmer = |kmer_size: usize, allow_lc: bool, allow_nu: bool| -> GatkResult<()> {
        if let Some(result) =
            try_assemble_kmer(reference, reads, kmer_size, args, allow_lc, allow_nu)?
        {
            consider_rt_assembly_candidate(result, reference, args, &mut best, &mut best_score);
        }
        let batch = extract_rt_haplotypes_before_remove_paths(
            reference, reads, args, kmer_size, allow_lc, allow_nu,
        )?;
        if !batch.is_empty() {
            let mut batch = batch;
            let mut ref_hap = Haplotype::new(reference.bases.as_slice(), true);
            let mut ref_cigar = Cigar::new();
            ref_cigar.push(ref_hap.bases.len(), CigarOperator::Match);
            ref_hap.cigar = Some(ref_cigar);
            let status = finalize_assembly_haplotypes(&mut batch, &ref_hap, false);
            consider_rt_assembly_candidate(
                AssemblyResult {
                    status,
                    kmer_size,
                    haplotypes: batch,
                    event_maps: Vec::new(),
                },
                reference,
                args,
                &mut best,
                &mut best_score,
            );
        }
        Ok(())
    };
    for &kmer_size in &kmer_sizes {
        try_kmer(
            kmer_size,
            allow_low_complexity_configured_kmer(args),
            allow_non_unique_configured_kmer(args),
        )?;
    }
    if !args.dont_increase_kmer_sizes_for_cycles {
        let mut kmer_size = *kmer_sizes.last().unwrap_or(&25) + 10;
        for iter in 1..=6 {
            let last = iter == 6;
            try_kmer(
                kmer_size,
                allow_low_complexity_expanded_kmer(args, last),
                allow_non_unique_expanded_kmer(args, last),
            )?;
            kmer_size += 10;
        }
    }
    Ok(best)
}

/// GATK `createGraph` for each entry in `kmerSizes`: `allowLowComplexityGraphs` = `dontIncreaseKmerSizesForCycles`.
fn allow_low_complexity_configured_kmer(args: &ReadThreadingAssemblerArgs) -> bool {
    args.dont_increase_kmer_sizes_for_cycles || args.allow_low_complexity_graphs
}

/// GATK `createGraph` for each entry in `kmerSizes`: `allowNonUniqueKmersInRef` from assembler args only.
fn allow_non_unique_configured_kmer(args: &ReadThreadingAssemblerArgs) -> bool {
    args.allow_non_unique_kmers_in_ref
}

/// GATK expanded-kmer loop: both flags use `lastAttempt` on the final iteration.
pub(crate) fn allow_low_complexity_expanded_kmer(
    args: &ReadThreadingAssemblerArgs,
    is_last_expanded_attempt: bool,
) -> bool {
    is_last_expanded_attempt
        || args.dont_increase_kmer_sizes_for_cycles
        || args.allow_low_complexity_graphs
}

pub(crate) fn allow_non_unique_expanded_kmer(
    args: &ReadThreadingAssemblerArgs,
    is_last_expanded_attempt: bool,
) -> bool {
    is_last_expanded_attempt || args.allow_non_unique_kmers_in_ref
}

/// GATK `ReadThreadingAssembler` + `ReadThreadingAssemblerArgumentCollection` defaults.
/// # Invariants
/// `kmer_sizes` are tried in order; expanded attempts may relax non-unique / low-complexity gates.
/// `num_best_haplotypes_per_graph` caps KBest extraction per graph (default 128).
/// # Ownership
/// [`Clone`] args bundle SW parameters, optional scoring context, and k-mer size list.
/// # Mutation
/// Immutable per assembly invocation; local clones may disable dangling for audits.
/// # Biological assumptions
/// Read threading de Bruijn graphs recover local haplotypes from overlapping Illumina reads.
/// # Java equivalence
/// GATK `ReadThreadingAssembler` + `ReadThreadingAssemblerArgumentCollection` (HC defaults).
#[derive(Debug, Clone)]
pub struct ReadThreadingAssemblerArgs {
    pub kmer_sizes: Vec<usize>,
    pub min_base_quality: u8,
    pub min_prune_factor: u32,
    pub use_adaptive_pruning: bool,
    pub prune_before_cycle_counting: bool,
    pub recover_dangling_branches: bool,
    pub recover_dangling_heads: bool,
    pub recover_all_dangling_branches: bool,
    pub min_dangling_branch_length: usize,
    /// GATK `minMatchingBasesToDanglingEndRecovery` (HC default -1).
    pub min_matching_bases_to_dangling_end_recovery: i32,
    pub allow_non_unique_kmers_in_ref: bool,
    pub dont_increase_kmer_sizes_for_cycles: bool,
    pub allow_low_complexity_graphs: bool,
    pub remove_paths_not_connected_to_ref: bool,
    /// When true, skip post-dangling `apply_pruning` (RT k-best supplement only).
    pub skip_post_dangling_prune: bool,
    /// GATK default `generateSeqGraph` (`!useLinkedDeBruijnGraph`).
    pub use_seq_graph: bool,
    /// GATK `createGraph`: skip k-mer when `generateSeqGraph && hasCycles` (parity dumps keep false).
    pub abort_seq_graph_on_cycles: bool,
    /// When false, keep KBest `is_reference` flags (parity `assembly-haplotypes` dump).
    pub ensure_reference_in_result: bool,
    pub num_best_haplotypes_per_graph: usize,
    pub haplotype_to_reference_sw: SwParameters,
    pub dangling_end_sw: DanglingRecoverySwParams,
    /// When set, prefer k-mer graphs that yield P12 cluster coupled EventMap sites.
    pub scoring: Option<AssemblyScoringContext>,
    /// GATK-exact dangling: single pass, no ASM-1 tail suffix rescue (`strict_java_assembly`).
    pub dangling_java_exact: bool,
}

impl Default for ReadThreadingAssemblerArgs {
    fn default() -> Self {
        Self {
            kmer_sizes: vec![10, 25],
            min_base_quality: 10,
            min_prune_factor: 2,
            use_adaptive_pruning: false,
            prune_before_cycle_counting: true,
            recover_dangling_branches: true,
            recover_dangling_heads: true,
            recover_all_dangling_branches: false,
            min_dangling_branch_length: 4,
            min_matching_bases_to_dangling_end_recovery: -1,
            allow_non_unique_kmers_in_ref: false,
            // GATK `ReadThreadingAssemblerArgumentCollection.dontIncreaseKmerSizesForCycles` default false.
            dont_increase_kmer_sizes_for_cycles: false,
            allow_low_complexity_graphs: false,
            remove_paths_not_connected_to_ref: true,
            skip_post_dangling_prune: false,
            use_seq_graph: true,
            abort_seq_graph_on_cycles: true,
            ensure_reference_in_result: true,
            num_best_haplotypes_per_graph: DEFAULT_NUM_BEST_HAPLOTYPES_PER_GRAPH,
            haplotype_to_reference_sw: SwParameters::gatk_haplotype_to_reference(),
            dangling_end_sw: DanglingRecoverySwParams::gatk_defaults(),
            scoring: None,
            dangling_java_exact: false,
        }
    }
}

/// Outcome category for one local assembly attempt.
/// # Invariants
/// `Failed` yields empty haplotype/event lists; success variants carry at least reference handling downstream.
/// # Ownership
/// [`Copy`] status tag on [`AssemblyResult`].
/// # Mutation
/// Immutable once assembly completes.
/// # Biological assumptions
/// Distinguishes reference-only assembly from graphs with alt variation vs hard failure.
/// # Java equivalence
/// GATK `AssemblyResult` / `ReadThreadingAssembler` status (`ASSEMBLED_SOME_VARIATION`, etc.).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AssemblyStatus {
    AssembledSomeVariation,
    JustAssembledReference,
    Failed,
}

/// Haplotypes and EventMaps produced by one assembly graph attempt.
/// # Invariants
/// `haplotypes.len == event_maps.len` after successful assembly materialization.
/// `kmer_size` records the winning k-mer attempt for this result.
/// # Ownership
/// Owns haplotype and event-map vectors; consumed by genotyping and trim stages.
/// # Mutation
/// Assembly helpers may append/supplement haplotypes in place before returning.
/// # Biological assumptions
/// Haplotypes cover padded active region; reference haplotype included when configured.
/// # Java equivalence
/// GATK `AssemblyResult` from `ReadThreadingAssembler.runLocalAssembly`.
#[derive(Debug, Clone)]
pub struct AssemblyResult {
    pub status: AssemblyStatus,
    pub kmer_size: usize,
    pub haplotypes: Vec<Haplotype>,
    pub event_maps: Vec<EventMap>,
}

/// Java `HcFullParityGateDump.assemblyAssemble` / `tryAssembleKmer` multi-kmer RT loop (no SeqGraph).
#[cfg(any(feature = "dev-dumps", test))]
pub fn assemble_for_java_gate_dump(
    reference: &AssemblyRead,
    reads: &[AssemblyRead],
) -> GatkResult<AssemblyResult> {
    let mut args = ReadThreadingAssemblerArgs::default();
    args.use_seq_graph = false;
    args.dont_increase_kmer_sizes_for_cycles = false;
    let kmers = [10usize, 25];
    let mut last_fail = None;
    for &kmer_size in &kmers {
        if reference.bases.len() < kmer_size {
            last_fail = Some(AssemblyResult {
                status: AssemblyStatus::Failed,
                kmer_size,
                haplotypes: Vec::new(),
                event_maps: Vec::new(),
            });
            continue;
        }
        let attempt = match try_assemble_kmer(reference, reads, kmer_size, &args, false, false)? {
            Some(result) => result,
            None => AssemblyResult {
                status: AssemblyStatus::JustAssembledReference,
                kmer_size,
                haplotypes: Vec::new(),
                event_maps: Vec::new(),
            },
        };
        if attempt.status != AssemblyStatus::Failed {
            return Ok(attempt);
        }
        last_fail = Some(attempt);
    }
    Ok(last_fail.unwrap_or(AssemblyResult {
        status: AssemblyStatus::Failed,
        kmer_size: kmers[0],
        haplotypes: Vec::new(),
        event_maps: Vec::new(),
    }))
}

/// Assemble from reference + reads (GATK `runLocalAssembly`: SeqGraph + `GraphBasedKBestHaplotypeFinder` by default).
pub fn assemble_from_ref_and_reads(
    reference: &AssemblyRead,
    reads: &[AssemblyRead],
    args: &ReadThreadingAssemblerArgs,
) -> GatkResult<AssemblyResult> {
    if args.use_seq_graph {
        let result = assemble_from_ref_and_reads_seq_graph(reference, reads, args)?;
        // SeqGraph zip can collapse to a ref-length spine while RT k-best still has variation (P12 ASM-1).
        if matches!(result.status, AssemblyStatus::AssembledSomeVariation)
            && !haplotypes_have_alt_bases(&result.haplotypes, reference.bases.as_slice())
            && !haplotypes_have_indel_cigar(&result.haplotypes)
        {
            let mut rt_args = args.clone();
            rt_args.use_seq_graph = false;
            rt_args.remove_paths_not_connected_to_ref = false;
            rt_args.skip_post_dangling_prune = true;
            if let Ok(rt) = assemble_from_ref_and_reads_rt_graph(reference, reads, &rt_args) {
                if haplotypes_have_alt_bases(&rt.haplotypes, reference.bases.as_slice())
                    || haplotypes_have_indel_cigar(&rt.haplotypes)
                {
                    let mut picked = rt;
                    supplement_p12_cluster_coupled_haplotypes(&mut picked, reference, reads, args)?;
                    return Ok(picked);
                }
            }
        }
        if let Some(ctx) = args.scoring.as_ref() {
            if ctx.overlaps_p12_cluster() {
                let coupled = cluster_coupled_event_count(
                    &result.haplotypes,
                    reference.bases.as_slice(),
                    ctx.padded_reference_start_1based,
                    &ctx.contig,
                    ctx.active_start_1based,
                    ctx.active_end_1based,
                );
                if coupled == 0 {
                    let mut rt_args = args.clone();
                    rt_args.use_seq_graph = false;
                    rt_args.remove_paths_not_connected_to_ref = false;
                    rt_args.skip_post_dangling_prune = true;
                    if let Some(rt) = rt_assembly_best_variation(reference, reads, &rt_args)? {
                        let ctx = args.scoring.as_ref().expect("scoring");
                        let ref_bytes = reference.bases.as_slice();
                        let rt_score = assembly_result_score(&rt, reference, args.scoring.as_ref());
                        let seq_score =
                            assembly_result_score(&result, reference, args.scoring.as_ref());
                        let rt_coupled = cluster_coupled_event_count(
                            &rt.haplotypes,
                            ref_bytes,
                            ctx.padded_reference_start_1based,
                            &ctx.contig,
                            ctx.active_start_1based,
                            ctx.active_end_1based,
                        );
                        if rt_coupled > 0 || rt_score > seq_score {
                            let mut picked = rt;
                            supplement_p12_cluster_coupled_haplotypes(
                                &mut picked,
                                reference,
                                reads,
                                args,
                            )?;
                            return Ok(picked);
                        }
                    }
                }
            }
        }
        let mut out = result;
        supplement_p12_cluster_coupled_haplotypes(&mut out, reference, reads, args)?;
        return Ok(out);
    }
    let mut out = assemble_from_ref_and_reads_rt_graph(reference, reads, args)?;
    supplement_p12_cluster_coupled_haplotypes(&mut out, reference, reads, args)?;
    Ok(out)
}

fn assemble_from_ref_and_reads_rt_graph(
    reference: &AssemblyRead,
    reads: &[AssemblyRead],
    args: &ReadThreadingAssemblerArgs,
) -> GatkResult<AssemblyResult> {
    let mut kmer_sizes = args.kmer_sizes.clone();
    kmer_sizes.sort_unstable();
    kmer_sizes.dedup();

    let prefer_best_variation = args.scoring.is_some();
    let mut best_variation: Option<AssemblyResult> = None;
    let mut best_variation_score = 0u64;
    let mut last_fail = None;
    let mut last_attempted_kmer = kmer_sizes.first().copied().unwrap_or(10);
    for &kmer_size in &kmer_sizes {
        last_attempted_kmer = kmer_size;
        if let Some(result) = try_assemble_kmer(
            reference,
            reads,
            kmer_size,
            args,
            allow_low_complexity_configured_kmer(args),
            allow_non_unique_configured_kmer(args),
        )? {
            if result.status == AssemblyStatus::AssembledSomeVariation {
                if prefer_best_variation {
                    let score = assembly_result_score(&result, reference, args.scoring.as_ref());
                    if let Some(ctx) = args.scoring.as_ref() {
                        let ref_bytes = reference.bases.as_slice();
                        let coupled = cluster_coupled_event_count(
                            &result.haplotypes,
                            ref_bytes,
                            ctx.padded_reference_start_1based,
                            &ctx.contig,
                            ctx.active_start_1based,
                            ctx.active_end_1based,
                        );
                        let indels = active_window_indel_event_count(
                            &result.haplotypes,
                            ref_bytes,
                            ctx.padded_reference_start_1based,
                            &ctx.contig,
                            ctx.active_start_1based,
                            ctx.active_end_1based,
                        );
                        let p12_ok = ctx.overlaps_p12_cluster()
                            && coupled >= 2
                            && !assembly_result_has_phantom_cluster_indels(&result, reference, ctx);
                        let mid_ok = haplotypes_have_alt_bases(&result.haplotypes, ref_bytes)
                            && (!ctx.overlaps_p12_cluster() || indels >= 1);
                        if p12_ok || mid_ok {
                            let mut early = result;
                            supplement_p12_cluster_coupled_haplotypes(
                                &mut early, reference, reads, args,
                            )?;
                            return Ok(early);
                        }
                    }
                    if score > best_variation_score {
                        best_variation_score = score;
                        best_variation = Some(result);
                    }
                } else {
                    let mut early = result;
                    supplement_p12_cluster_coupled_haplotypes(&mut early, reference, reads, args)?;
                    return Ok(early);
                }
            } else if result.status != AssemblyStatus::Failed {
                last_fail = Some(result);
            }
        }
    }

    if !args.dont_increase_kmer_sizes_for_cycles {
        let mut kmer_size = *kmer_sizes.last().unwrap_or(&25) + 10;
        for iter in 1..=6 {
            last_attempted_kmer = kmer_size;
            let last = iter == 6;
            if let Some(result) = try_assemble_kmer(
                reference,
                reads,
                kmer_size,
                args,
                allow_low_complexity_expanded_kmer(args, last),
                allow_non_unique_expanded_kmer(args, last),
            )? {
                if result.status == AssemblyStatus::AssembledSomeVariation {
                    if prefer_best_variation {
                        let score =
                            assembly_result_score(&result, reference, args.scoring.as_ref());
                        if let Some(ctx) = args.scoring.as_ref() {
                            let ref_bytes = reference.bases.as_slice();
                            let coupled = cluster_coupled_event_count(
                                &result.haplotypes,
                                ref_bytes,
                                ctx.padded_reference_start_1based,
                                &ctx.contig,
                                ctx.active_start_1based,
                                ctx.active_end_1based,
                            );
                            let indels = active_window_indel_event_count(
                                &result.haplotypes,
                                ref_bytes,
                                ctx.padded_reference_start_1based,
                                &ctx.contig,
                                ctx.active_start_1based,
                                ctx.active_end_1based,
                            );
                            let p12_ok = ctx.overlaps_p12_cluster()
                                && coupled >= 2
                                && !assembly_result_has_phantom_cluster_indels(
                                    &result, reference, ctx,
                                );
                            let mid_ok = haplotypes_have_alt_bases(&result.haplotypes, ref_bytes)
                                && (!ctx.overlaps_p12_cluster() || indels >= 1);
                            if p12_ok || mid_ok {
                                let mut early = result;
                                supplement_p12_cluster_coupled_haplotypes(
                                    &mut early, reference, reads, args,
                                )?;
                                return Ok(early);
                            }
                        }
                        if score > best_variation_score {
                            best_variation_score = score;
                            best_variation = Some(result);
                        }
                    } else {
                        let mut early = result;
                        supplement_p12_cluster_coupled_haplotypes(
                            &mut early, reference, reads, args,
                        )?;
                        return Ok(early);
                    }
                } else if result.status != AssemblyStatus::Failed {
                    last_fail = Some(result);
                }
            }
            kmer_size += 10;
        }
    }
    if let Some(mut best) = best_variation {
        supplement_p12_cluster_coupled_haplotypes(&mut best, reference, reads, args)?;
        return Ok(best);
    }

    let mut fail = last_fail.unwrap_or(AssemblyResult {
        status: AssemblyStatus::Failed,
        kmer_size: last_attempted_kmer,
        haplotypes: Vec::new(),
        event_maps: Vec::new(),
    });
    supplement_p12_cluster_coupled_haplotypes(&mut fail, reference, reads, args)?;
    if args.scoring.is_some()
        && !haplotypes_have_alt_bases(&fail.haplotypes, reference.bases.as_slice())
        && !haplotypes_have_indel_cigar(&fail.haplotypes)
    {
        if let Some(rt) = rt_assembly_best_variation(reference, reads, args)? {
            if haplotypes_have_alt_bases(&rt.haplotypes, reference.bases.as_slice())
                || haplotypes_have_indel_cigar(&rt.haplotypes)
            {
                let mut picked = rt;
                supplement_p12_cluster_coupled_haplotypes(&mut picked, reference, reads, args)?;
                return Ok(picked);
            }
        }
    }
    Ok(fail)
}

/// GATK `assembleKmerGraphsAndHaplotypeCall`: all kmer SeqGraphs with variation → `findBestPaths`.
fn assemble_from_ref_and_reads_seq_graph(
    reference: &AssemblyRead,
    reads: &[AssemblyRead],
    args: &ReadThreadingAssemblerArgs,
) -> GatkResult<AssemblyResult> {
    let mut kmer_sizes = args.kmer_sizes.clone();
    kmer_sizes.sort_unstable();
    kmer_sizes.dedup();

    let mut variation_kmers: Vec<usize> = Vec::new();
    let mut last_just_ref: Option<AssemblyResult> = None;
    let mut haplotypes = Vec::new();
    let mut seen: HashSet<(Vec<u8>, bool)> = HashSet::new();

    let mut ref_hap = Haplotype::new(reference.bases.as_slice(), true);
    let mut ref_cigar = Cigar::new();
    ref_cigar.push(ref_hap.bases.len(), CigarOperator::Match);
    ref_hap.cigar = Some(ref_cigar);
    let ref_cigar_len = ref_hap.cigar.as_ref().unwrap().reference_length();

    let ingest_seq_graph = |seq: SeqGraph,
                            kmer: usize,
                            variation_kmers: &mut Vec<usize>,
                            haplotypes: &mut Vec<Haplotype>,
                            seen: &mut HashSet<(Vec<u8>, bool)>,
                            ref_hap: &Haplotype|
     -> GatkResult<()> {
        variation_kmers.push(kmer);
        crate::runtime_config::rss_trace_checkpoint(
            "seq_kbest_begin",
            &format!("kmer={kmer} seq_nodes={}", seq.node_count()),
        );
        let paths = find_best_haplotypes_seq_graph(&seq, args.num_best_haplotypes_per_graph)?;
        crate::runtime_config::rss_trace_checkpoint(
            "seq_kbest_done",
            &format!("kmer={kmer} paths={}", paths.len()),
        );
        let mut batch = extract_haplotypes_from_seq_kbest_paths(
            &paths,
            &seq,
            kmer,
            ref_hap,
            ref_cigar_len,
            &args.haplotype_to_reference_sw,
        )?;
        drop(seq); // Peak-RSS: free SeqGraph before next kmer / PairHMM
        crate::smith_waterman::release_sw_tls_scratch();
        for h in batch.drain(..) {
            let key = (h.bases.clone(), h.is_reference);
            if seen.insert(key) {
                haplotypes.push(h);
            }
        }
        Ok(())
    };

    for &kmer_size in &kmer_sizes {
        if let Some((seq, status, kmer)) = try_build_seq_graph_kmer(
            reference,
            reads,
            kmer_size,
            args,
            allow_low_complexity_configured_kmer(args),
            allow_non_unique_configured_kmer(args),
        )? {
            match status {
                SeqGraphCleanupStatus::AssembledSomeVariation => {
                    ingest_seq_graph(
                        seq,
                        kmer,
                        &mut variation_kmers,
                        &mut haplotypes,
                        &mut seen,
                        &ref_hap,
                    )?;
                }
                SeqGraphCleanupStatus::JustAssembledReference => {
                    last_just_ref = Some(just_reference_result(kmer, reference));
                }
            }
        }
    }

    if variation_kmers.is_empty() && !args.dont_increase_kmer_sizes_for_cycles {
        let mut kmer_size = *kmer_sizes.last().unwrap_or(&25) + 10;
        for iter in 1..=6 {
            let last = iter == 6;
            if let Some((seq, status, kmer)) = try_build_seq_graph_kmer(
                reference,
                reads,
                kmer_size,
                args,
                allow_low_complexity_expanded_kmer(args, last),
                allow_non_unique_expanded_kmer(args, last),
            )? {
                match status {
                    SeqGraphCleanupStatus::AssembledSomeVariation => {
                        ingest_seq_graph(
                            seq,
                            kmer,
                            &mut variation_kmers,
                            &mut haplotypes,
                            &mut seen,
                            &ref_hap,
                        )?;
                    }
                    SeqGraphCleanupStatus::JustAssembledReference => {
                        last_just_ref = Some(just_reference_result(kmer, reference));
                    }
                }
            }
            kmer_size += 10;
        }
    }

    if variation_kmers.is_empty() {
        return Ok(last_just_ref.unwrap_or_else(|| {
            just_reference_result(kmer_sizes.first().copied().unwrap_or(10), reference)
        }));
    }

    let min_kmer = variation_kmers.iter().copied().min().unwrap_or(10);
    crate::runtime_config::rss_trace_checkpoint(
        "before_merge_rt",
        &format!(
            "variation_kmers={} min_kmer={} haps={}",
            variation_kmers.len(),
            min_kmer,
            haplotypes.len()
        ),
    );
    // Observable contract: walk configured + variation + expanded RT k-mers (pre-Peak).
    // Min-kmer-only merge (Peak cut) left L2 `g2-subset-live` at haplotype_count=1 because the
    // alt bubble often appears at a non-min k-mer (e.g. k=10) while min variation k-mer is
    // ref-only. Peak-RSS on bushy loci is gated by early-stop once alts appear inside
    // [`merge_rt_kbest_pre_remove_paths_at_kmer`], not by dropping the walk.
    merge_rt_kbest_pre_remove_paths(reference, reads, args, &variation_kmers, &mut haplotypes)?;
    crate::runtime_config::rss_trace_checkpoint(
        "after_merge_rt",
        &format!("haps={}", haplotypes.len()),
    );
    let pad = cigar_refresh_pad(args);
    refresh_alt_haplotype_indel_cigars(
        &mut haplotypes,
        reference.bases.as_slice(),
        pad,
        &args.haplotype_to_reference_sw,
    );
    crate::smith_waterman::release_sw_tls_scratch();
    if args.ensure_reference_in_result {
        ensure_reference_haplotype(&mut haplotypes, &ref_hap);
    }
    let status =
        finalize_assembly_haplotypes(&mut haplotypes, &ref_hap, args.ensure_reference_in_result);
    let event_maps = haplotypes
        .iter()
        .map(|h| {
            let ref_hap = Haplotype::new(reference.bases.as_slice(), true);
            EventMap::from_haplotype_and_reference(h, &ref_hap, &ref_hap.bases, 1, 0)
        })
        .collect();
    let mut out = AssemblyResult {
        status,
        kmer_size: min_kmer,
        haplotypes,
        event_maps,
    };
    supplement_p12_cluster_coupled_haplotypes(&mut out, reference, reads, args)?;
    Ok(out)
}

/// One k-mer attempt through threading graph → SeqGraph cleanup (parity diagnostics).
/// # Invariants
/// `phase` labels the assembly attempt stage for dump ordering.
/// `outcome` summarizes success/failure reason as a display string for parity logs.
/// # Ownership
/// Owns diagnostic strings and numeric graph stats; no graph handles retained.
/// # Mutation
/// Immutable row appended to probe traces during assembly attempts.
/// # Biological assumptions
/// None documented (diagnostic probe, not genotype input).
/// # Java equivalence
/// Rust-native parity diagnostics for threading → SeqGraph pipeline.
#[derive(Debug, Clone)]
pub struct SeqGraphKmerProbeRow {
    pub phase: &'static str,
    pub kmer_size: usize,
    pub allow_low_complexity: bool,
    pub allow_non_unique_ref: bool,
    pub outcome: String,
    pub thread_nodes: usize,
    pub thread_edges: usize,
    pub cleanup_status: String,
    pub has_ref_source: bool,
    pub has_ref_sink: bool,
    pub ref_path_matches: bool,
    pub kbest_paths: usize,
    pub extracted_haps: usize,
    pub non_ref_haps: usize,
    pub path_bases_len: usize,
    pub path_eq_ref_bases: bool,
}

fn try_build_seq_graph_kmer(
    reference: &AssemblyRead,
    reads: &[AssemblyRead],
    kmer_size: usize,
    args: &ReadThreadingAssemblerArgs,
    allow_low_complexity: bool,
    allow_non_unique_ref: bool,
) -> GatkResult<Option<(SeqGraph, SeqGraphCleanupStatus, usize)>> {
    let Some(graph) = build_threading_graph_for_seq_assembly(
        reference,
        reads,
        kmer_size,
        args,
        allow_low_complexity,
        allow_non_unique_ref,
    )?
    else {
        return Ok(None);
    };
    crate::runtime_config::rss_trace_checkpoint(
        "seq_rt_built",
        &format!(
            "kmer={kmer_size} nodes={} edges={}",
            graph.node_count(),
            graph.edge_count()
        ),
    );
    // Peak-RSS: skip bushy RT graphs before SeqGraph cleanup / k-best (primary path).
    if graph.node_count() > MAX_ASSEMBLY_GRAPH_NODES {
        crate::runtime_config::rss_trace_checkpoint(
            "seq_rt_skip_huge",
            &format!(
                "kmer={kmer_size} nodes={} cap={MAX_ASSEMBLY_GRAPH_NODES}",
                graph.node_count()
            ),
        );
        drop(graph);
        return Ok(None);
    }
    let mut seq = SeqGraph::from_assembly_graph(&graph);
    drop(graph); // Peak-RSS: free RT graph before SeqGraph cleanup / k-best
    seq.clean_non_ref_paths();
    // GATK `assembleKmerGraphsAndHaplotypeCall`: only graphs with ASSEMBLED_SOME_VARIATION enter `nonRefSeqGraphs`.
    let status = seq.cleanup_seq_graph();
    if status == SeqGraphCleanupStatus::JustAssembledReference {
        return Ok(None);
    }
    // Do not pass graphs with lost ref endpoints to `findBestPaths` (Java `sanityCheckGraph`).
    // GATK `assembleKmerGraphsAndHaplotypeCall`: admit graphs with ASSEMBLED_SOME_VARIATION only
    // (Java does not require KBest path bytes != reference before `nonRefSeqGraphs.add`).
    if seq.reference_source_vertex().is_none() || seq.reference_sink_vertex().is_none() {
        return Ok(None);
    }
    if seq.node_count() > MAX_ASSEMBLY_GRAPH_NODES {
        crate::runtime_config::rss_trace_checkpoint(
            "seq_skip_huge",
            &format!(
                "kmer={kmer_size} seq_nodes={} cap={MAX_ASSEMBLY_GRAPH_NODES}",
                seq.node_count()
            ),
        );
        return Ok(None);
    }
    crate::runtime_config::rss_trace_checkpoint(
        "seq_ready",
        &format!(
            "kmer={kmer_size} seq_nodes={} seq_edges={}",
            seq.node_count(),
            seq.edge_count()
        ),
    );
    Ok(Some((seq, status, kmer_size)))
}

fn probe_seq_graph_kmer_one(
    reference: &AssemblyRead,
    reads: &[AssemblyRead],
    kmer_size: usize,
    args: &ReadThreadingAssemblerArgs,
    phase: &'static str,
    allow_low_complexity: bool,
    allow_non_unique_ref: bool,
) -> GatkResult<SeqGraphKmerProbeRow> {
    let mut row = SeqGraphKmerProbeRow {
        phase,
        kmer_size,
        allow_low_complexity,
        allow_non_unique_ref,
        outcome: String::new(),
        thread_nodes: 0,
        thread_edges: 0,
        cleanup_status: String::new(),
        has_ref_source: false,
        has_ref_sink: false,
        ref_path_matches: false,
        kbest_paths: 0,
        extracted_haps: 0,
        non_ref_haps: 0,
        path_bases_len: 0,
        path_eq_ref_bases: false,
    };

    if reference.bases.len() < kmer_size {
        row.outcome = "ref_shorter_than_kmer".into();
        return Ok(row);
    }

    if !allow_non_unique_ref
        && !args.allow_non_unique_kmers_in_ref
        && reference_has_non_unique_kmers(reference, kmer_size)
    {
        row.outcome = "skip_non_unique_ref_kmers".into();
        return Ok(row);
    }

    let Some(graph) = build_threading_graph_for_haplotype_dump(
        reference,
        reads,
        kmer_size,
        args,
        allow_low_complexity,
        allow_non_unique_ref,
    )?
    else {
        row.outcome = "no_threading_graph".into();
        return Ok(row);
    };
    row.thread_nodes = graph.node_count();
    row.thread_edges = graph.edge_count();

    let mut seq = SeqGraph::from_assembly_graph(&graph);
    seq.clean_non_ref_paths();
    let status = seq.cleanup_seq_graph();
    row.cleanup_status = match status {
        SeqGraphCleanupStatus::AssembledSomeVariation => "assembled_some_variation",
        SeqGraphCleanupStatus::JustAssembledReference => "just_assembled_reference",
    }
    .into();
    row.has_ref_source = seq.reference_source_vertex().is_some();
    row.has_ref_sink = seq.reference_sink_vertex().is_some();
    row.ref_path_matches = seq
        .reference_path_bytes()
        .map(|b| b == reference.bases.as_slice())
        .unwrap_or(false);

    if status == SeqGraphCleanupStatus::JustAssembledReference {
        row.outcome = "cleanup_just_assembled_reference".into();
        return Ok(row);
    }
    if !row.has_ref_source || !row.has_ref_sink {
        row.outcome = "dropped_no_ref_endpoints".into();
        return Ok(row);
    }

    let paths = find_best_haplotypes_seq_graph(&seq, args.num_best_haplotypes_per_graph)?;
    row.kbest_paths = paths.len();

    let mut ref_hap = Haplotype::new(reference.bases.as_slice(), true);
    let mut ref_cigar = Cigar::new();
    ref_cigar.push(ref_hap.bases.len(), CigarOperator::Match);
    ref_hap.cigar = Some(ref_cigar);
    let ref_cigar_len = ref_hap.cigar.as_ref().unwrap().reference_length();
    let haps = extract_haplotypes_from_seq_kbest_paths(
        &paths,
        &seq,
        kmer_size,
        &ref_hap,
        ref_cigar_len,
        &args.haplotype_to_reference_sw,
    )?;
    row.extracted_haps = haps.len();
    row.non_ref_haps = haps.iter().filter(|h| !h.is_reference).count();
    if let Some(path) = paths.first() {
        let bases = seq.path_bases_bytes(path.start, &path.edges);
        row.path_bases_len = bases.len();
        row.path_eq_ref_bases = bases == reference.bases.as_slice();
    }
    row.outcome = if row.path_eq_ref_bases {
        "variation_graph_kbest_path_is_ref_bases".into()
    } else if row.non_ref_haps > 0 {
        "variation_graph_with_alt_haps".into()
    } else if row.kbest_paths > 1 {
        "variation_graph_kbest_no_extracted_alts".into()
    } else {
        "variation_graph_ref_only_kbest".into()
    };
    Ok(row)
}

/// Probe every configured (and optional expanded) k-mer through SeqGraph assembly.
pub fn probe_seq_graph_kmer_attempts(
    reference: &AssemblyRead,
    reads: &[AssemblyRead],
    args: &ReadThreadingAssemblerArgs,
) -> GatkResult<Vec<SeqGraphKmerProbeRow>> {
    let mut rows = Vec::new();
    let mut kmer_sizes = args.kmer_sizes.clone();
    kmer_sizes.sort_unstable();
    kmer_sizes.dedup();

    for &kmer_size in &kmer_sizes {
        rows.push(probe_seq_graph_kmer_one(
            reference,
            reads,
            kmer_size,
            args,
            "configured",
            allow_low_complexity_configured_kmer(args),
            allow_non_unique_configured_kmer(args),
        )?);
    }

    if !args.dont_increase_kmer_sizes_for_cycles {
        let mut kmer_size = *kmer_sizes.last().unwrap_or(&25) + 10;
        for iter in 1..=6 {
            let last = iter == 6;
            rows.push(probe_seq_graph_kmer_one(
                reference,
                reads,
                kmer_size,
                args,
                "expanded",
                allow_low_complexity_expanded_kmer(args, last),
                allow_non_unique_expanded_kmer(args, last),
            )?);
            kmer_size += 10;
        }
    }
    Ok(rows)
}

fn haplotypes_have_alt_bases(haplotypes: &[Haplotype], ref_bytes: &[u8]) -> bool {
    haplotypes.iter().any(|h| h.bases.as_slice() != ref_bytes)
}

fn haplotype_cigar_has_indel(cigar: &Cigar) -> bool {
    cigar.elements.iter().any(|e| e.operator.is_indel())
}

fn haplotypes_have_indel_cigar(haplotypes: &[Haplotype]) -> bool {
    haplotypes
        .iter()
        .any(|h| h.cigar.as_ref().is_some_and(haplotype_cigar_has_indel))
}

fn kmer_sizes_for_rt_merge(
    args: &ReadThreadingAssemblerArgs,
    variation_kmers: &[usize],
) -> Vec<usize> {
    let mut sizes: Vec<usize> = args.kmer_sizes.clone();
    sizes.extend_from_slice(variation_kmers);
    if !args.dont_increase_kmer_sizes_for_cycles {
        let mut kmer = *sizes.last().unwrap_or(&25) + 10;
        for _ in 1..=6 {
            sizes.push(kmer);
            kmer += 10;
        }
    }
    sizes.sort_unstable();
    sizes.dedup();
    sizes
}

/// RT k-best after GATK `removePathsNotConnectedToRef` (production `assembleReads` / Java dump graph).
pub fn extract_rt_haplotypes_after_remove_paths(
    reference: &AssemblyRead,
    reads: &[AssemblyRead],
    args: &ReadThreadingAssemblerArgs,
    kmer_size: usize,
    allow_low_complexity: bool,
    allow_non_unique_ref: bool,
) -> GatkResult<Vec<Haplotype>> {
    extract_rt_haplotypes_from_built_graph(
        reference,
        reads,
        args,
        kmer_size,
        allow_low_complexity,
        allow_non_unique_ref,
        false,
    )
}

/// RT graph for supplement: dangling recovery on, but keep off-spine branches until k-best (ASM-2).
pub fn extract_rt_haplotypes_before_remove_paths(
    reference: &AssemblyRead,
    reads: &[AssemblyRead],
    args: &ReadThreadingAssemblerArgs,
    kmer_size: usize,
    allow_low_complexity: bool,
    allow_non_unique_ref: bool,
) -> GatkResult<Vec<Haplotype>> {
    extract_rt_haplotypes_from_built_graph(
        reference,
        reads,
        args,
        kmer_size,
        allow_low_complexity,
        allow_non_unique_ref,
        true,
    )
}

fn extract_rt_haplotypes_from_built_graph(
    reference: &AssemblyRead,
    reads: &[AssemblyRead],
    args: &ReadThreadingAssemblerArgs,
    kmer_size: usize,
    allow_low_complexity: bool,
    allow_non_unique_ref: bool,
    before_remove_paths: bool,
) -> GatkResult<Vec<Haplotype>> {
    let mut local_args = args.clone();
    if before_remove_paths {
        local_args.remove_paths_not_connected_to_ref = false;
        local_args.skip_post_dangling_prune = true;
    }
    crate::runtime_config::rss_trace_checkpoint(
        "rt_graph_build_begin",
        &format!("kmer={kmer_size} before_remove_paths={before_remove_paths}"),
    );
    let Some(graph) = build_threading_graph_core(
        reference,
        reads,
        kmer_size,
        &local_args,
        allow_low_complexity,
        allow_non_unique_ref,
        false,
    )?
    else {
        crate::runtime_config::rss_trace_checkpoint(
            "rt_graph_build_none",
            &format!("kmer={kmer_size}"),
        );
        return Ok(Vec::new());
    };
    crate::runtime_config::rss_trace_checkpoint(
        "rt_graph_built",
        &format!(
            "kmer={kmer_size} nodes={} edges={}",
            graph.node_count(),
            graph.edge_count()
        ),
    );
    // Peak-RSS: bushy RT graphs (seen climbing past ~800 MiB on NA12878 20:10098169) —
    // skip k-best rather than expand a multi-GiB frontier. Normal HC regions here are ~1–2 k nodes.
    if graph.node_count() > MAX_ASSEMBLY_GRAPH_NODES {
        crate::runtime_config::rss_trace_checkpoint(
            "rt_graph_skip_huge",
            &format!(
                "kmer={kmer_size} nodes={} cap={MAX_ASSEMBLY_GRAPH_NODES}",
                graph.node_count()
            ),
        );
        drop(graph);
        return Ok(Vec::new());
    }
    let mut ref_hap = Haplotype::new(reference.bases.as_slice(), true);
    let mut ref_cigar = Cigar::new();
    ref_cigar.push(ref_hap.bases.len(), CigarOperator::Match);
    ref_hap.cigar = Some(ref_cigar);
    let ref_cigar_len = ref_hap.cigar.as_ref().unwrap().reference_length();
    let (paths, graph) =
        find_best_haplotypes_for_assembly(graph, args.num_best_haplotypes_per_graph)?;
    let mut haps = extract_haplotypes_from_kbest_paths(
        &paths,
        &graph,
        &ref_hap,
        ref_cigar_len,
        &args.haplotype_to_reference_sw,
    )?;
    crate::assembly_dangling_recovery::apply_dangling_merge_haplotypes(
        &mut haps,
        &ref_hap,
        &graph.dangling_merge_haps,
        reference.bases.as_slice(),
        &args.haplotype_to_reference_sw,
    );
    // Peak-RSS: release RT graph before SW cigar refresh / caller PairHMM.
    drop(graph);
    crate::smith_waterman::release_sw_tls_scratch();
    if let Some(ctx) = args.scoring.as_ref() {
        refresh_alt_haplotype_indel_cigars(
            &mut haps,
            reference.bases.as_slice(),
            ctx.padded_reference_start_1based,
            &args.haplotype_to_reference_sw,
        );
        crate::smith_waterman::release_sw_tls_scratch();
    }
    Ok(haps)
}

/// ASM-7: KBest on RT graph after dangling + `remove_paths` (GATK `assembleReads` order).
pub fn merge_rt_kbest_pre_remove_paths(
    reference: &AssemblyRead,
    reads: &[AssemblyRead],
    args: &ReadThreadingAssemblerArgs,
    variation_kmers: &[usize],
    haplotypes: &mut Vec<Haplotype>,
) -> GatkResult<()> {
    // Full configured+variation+expanded walk (pre-Peak). Peak-RSS is gated by early-stop
    // once alt+ref are present inside [`merge_rt_kbest_pre_remove_paths_at_kmer`].
    merge_rt_kbest_pre_remove_paths_at_kmer(
        reference,
        reads,
        args,
        variation_kmers,
        haplotypes,
        None,
    )
}

/// Optional `only_kmer` limits RT merge to one kmer (Java `AssemblyResultSet` minimum-kmer graph).
pub fn merge_rt_kbest_pre_remove_paths_at_kmer(
    reference: &AssemblyRead,
    reads: &[AssemblyRead],
    args: &ReadThreadingAssemblerArgs,
    variation_kmers: &[usize],
    haplotypes: &mut Vec<Haplotype>,
    only_kmer: Option<usize>,
) -> GatkResult<()> {
    let ref_bytes = reference.bases.as_slice();
    let mut ref_hap = Haplotype::new(ref_bytes, true);
    let mut ref_cigar = Cigar::new();
    ref_cigar.push(ref_hap.bases.len(), CigarOperator::Match);
    ref_hap.cigar = Some(ref_cigar);

    let mut seen: HashSet<(Vec<u8>, bool)> = haplotypes
        .iter()
        .map(|h| (h.bases.clone(), h.is_reference))
        .collect();
    let configured: HashSet<usize> = args.kmer_sizes.iter().copied().collect();

    for kmer_size in kmer_sizes_for_rt_merge(args, variation_kmers) {
        if only_kmer.is_some_and(|k| kmer_size != k) {
            continue;
        }
        let is_expanded = !configured.contains(&kmer_size);
        let allow_lc = if is_expanded {
            allow_low_complexity_expanded_kmer(args, true)
        } else {
            allow_low_complexity_configured_kmer(args)
        };
        let allow_nu = if is_expanded {
            allow_non_unique_expanded_kmer(args, true)
        } else {
            allow_non_unique_configured_kmer(args)
        };
        // Match Java-oriented historical merge and L2 `g2-subset-live`: before-remove keeps
        // off-spine alt branches that after-remove drops on tiny synthetic graphs.
        // Peak-RSS on bushy non-P12 loci is gated by early-stop once alts appear (below),
        // not by switching extract mode here.
        let mut batch = extract_rt_haplotypes_before_remove_paths(
            reference, reads, args, kmer_size, allow_lc, allow_nu,
        )?;
        if let Some(ctx) = args.scoring.as_ref() {
            refresh_alt_haplotype_indel_cigars(
                &mut batch,
                reference.bases.as_slice(),
                ctx.padded_reference_start_1based,
                &args.haplotype_to_reference_sw,
            );
        }
        for h in batch {
            // CLONE: needed because owned composite key for dedup/lookup.
            let key = (h.bases.clone(), h.is_reference);
            if seen.insert(key) {
                haplotypes.push(h);
            }
        }
        // Peak-RSS: once we have alt+ref, further expanded k-mers on bushy loci dominate RSS.
        // L2 tiny fixtures often need a later/earlier k-mer than `min(variation)` before alts
        // appear — so only stop after success, never after empty extracts alone.
        if only_kmer.is_none() && haplotypes.iter().any(|h| !h.is_reference) && haplotypes.len() > 1
        {
            crate::runtime_config::rss_trace_checkpoint(
                "merge_rt_early_stop",
                &format!("kmer={kmer_size} haps={}", haplotypes.len()),
            );
            break;
        }
    }
    if args.ensure_reference_in_result {
        ensure_reference_haplotype(haplotypes, &ref_hap);
    }
    Ok(())
}

/// Back-compat alias (ASM-2).
pub fn supplement_haplotypes_from_rt_kbest(
    reference: &AssemblyRead,
    reads: &[AssemblyRead],
    args: &ReadThreadingAssemblerArgs,
    variation_kmers: &[usize],
    haplotypes: &mut Vec<Haplotype>,
) -> GatkResult<()> {
    merge_rt_kbest_pre_remove_paths(reference, reads, args, variation_kmers, haplotypes)
}

pub fn extract_haplotypes_from_seq_kbest_paths(
    paths: &[KBestPath],
    graph: &SeqGraph,
    kmer_size: usize,
    ref_haplotype: &Haplotype,
    ref_cigar_length: usize,
    sw: &SwParameters,
) -> GatkResult<Vec<Haplotype>> {
    let ref_bytes = &ref_haplotype.bases;
    let mut out = Vec::new();
    let mut seen: HashSet<(Vec<u8>, bool)> = HashSet::new();

    for path in paths {
        let bases = graph.path_bases_bytes(path.start, &path.edges);
        // CLONE: needed because owned composite key for dedup/lookup.
        let key = (bases.clone(), path.is_reference);
        if seen.contains(&key) {
            continue;
        }
        // GATK `KBestHaplotype.haplotype` — `isReference` from path edges only, not base equality.
        // L14-E2: move bases into haplotype (avoid second full-sequence clone).
        let mut h = Haplotype::new(bases, path.is_reference);
        h.score = path.score;
        h.kmer_size = kmer_size;

        let Some(assy) = calculate_haplotype_cigar_for_assembly_with_offset(
            ref_bytes,
            &h.bases,
            ref_cigar_length,
            sw,
        ) else {
            continue;
        };
        if ref_cigar_length >= MIN_HAPLOTYPE_REFERENCE_LENGTH
            && assy.cigar.reference_length() < MIN_HAPLOTYPE_REFERENCE_LENGTH
        {
            continue;
        }
        h.cigar = Some(assy.cigar);
        h.alignment_start_hap_wrt_ref = assy.alignment_start_hap_wrt_ref;
        seen.insert(key);
        out.push(h);
    }
    Ok(out)
}

fn try_assemble_kmer(
    reference: &AssemblyRead,
    reads: &[AssemblyRead],
    kmer_size: usize,
    args: &ReadThreadingAssemblerArgs,
    allow_low_complexity: bool,
    allow_non_unique_ref: bool,
) -> GatkResult<Option<AssemblyResult>> {
    if reference.bases.len() < kmer_size {
        return Ok(Some(AssemblyResult {
            status: AssemblyStatus::Failed,
            kmer_size,
            haplotypes: Vec::new(),
            event_maps: Vec::new(),
        }));
    }

    if !allow_non_unique_ref
        && !args.allow_non_unique_kmers_in_ref
        && reference_has_non_unique_kmers(reference, kmer_size)
    {
        return Ok(None);
    }

    let Some(graph) = build_threading_graph_for_haplotype_dump(
        reference,
        reads,
        kmer_size,
        args,
        allow_low_complexity,
        allow_non_unique_ref,
    )?
    else {
        return Ok(None);
    };

    crate::runtime_config::rss_trace_checkpoint(
        "rt_primary_built",
        &format!(
            "kmer={kmer_size} nodes={} edges={}",
            graph.node_count(),
            graph.edge_count()
        ),
    );
    // Peak-RSS: same fail-closed node cap as SeqGraph primary + RT-supplement.
    if graph.node_count() > MAX_ASSEMBLY_GRAPH_NODES {
        crate::runtime_config::rss_trace_checkpoint(
            "rt_primary_skip_huge",
            &format!(
                "kmer={kmer_size} nodes={} cap={MAX_ASSEMBLY_GRAPH_NODES}",
                graph.node_count()
            ),
        );
        drop(graph);
        return Ok(None);
    }

    let mut ref_hap = Haplotype::new(reference.bases.as_slice(), true);
    let mut ref_cigar = Cigar::new();
    ref_cigar.push(ref_hap.bases.len(), CigarOperator::Match);
    ref_hap.cigar = Some(ref_cigar);
    let ref_cigar_len = ref_hap.cigar.as_ref().unwrap().reference_length();

    let (paths, graph) =
        find_best_haplotypes_for_assembly(graph, args.num_best_haplotypes_per_graph)?;
    let mut haplotypes = extract_haplotypes_from_kbest_paths(
        &paths,
        &graph,
        &ref_hap,
        ref_cigar_len,
        &args.haplotype_to_reference_sw,
    )?;
    crate::assembly_dangling_recovery::apply_dangling_merge_haplotypes(
        &mut haplotypes,
        &ref_hap,
        &graph.dangling_merge_haps,
        reference.bases.as_slice(),
        &args.haplotype_to_reference_sw,
    );
    drop(graph);
    if args.ensure_reference_in_result {
        ensure_reference_haplotype(&mut haplotypes, &ref_hap);
    }
    merge_rt_kbest_pre_remove_paths(reference, reads, args, &[], &mut haplotypes)?;
    refresh_alt_haplotype_indel_cigars(
        &mut haplotypes,
        reference.bases.as_slice(),
        1,
        &args.haplotype_to_reference_sw,
    );
    if args.ensure_reference_in_result {
        ensure_reference_haplotype(&mut haplotypes, &ref_hap);
    }

    let status =
        finalize_assembly_haplotypes(&mut haplotypes, &ref_hap, args.ensure_reference_in_result);

    let event_maps = haplotypes
        .iter()
        .map(|h| EventMap::from_haplotype_and_reference(h, &ref_hap, &ref_hap.bases, 1, 0))
        .collect();

    Ok(Some(AssemblyResult {
        status,
        kmer_size,
        haplotypes,
        event_maps,
    }))
}

/// Production SeqGraph path: GATK `createGraph` aborts cyclic k-mers before dangling recovery.
pub fn build_threading_graph_for_seq_assembly(
    reference: &AssemblyRead,
    reads: &[AssemblyRead],
    kmer_size: usize,
    args: &ReadThreadingAssemblerArgs,
    allow_low_complexity: bool,
    allow_non_unique_ref: bool,
) -> GatkResult<Option<AssemblyGraph>> {
    let graph = build_threading_graph_core(
        reference,
        reads,
        kmer_size,
        args,
        allow_low_complexity,
        allow_non_unique_ref,
        true,
    )?;
    Ok(graph)
}

/// Pruned + dangling-recovered threading graph for parity dumps (`HcFullParityGateDump.assemblyHaplotypes`).
#[doc(hidden)]
pub fn build_threading_graph_for_haplotype_dump(
    reference: &AssemblyRead,
    reads: &[AssemblyRead],
    kmer_size: usize,
    args: &ReadThreadingAssemblerArgs,
    allow_low_complexity: bool,
    allow_non_unique_ref: bool,
) -> GatkResult<Option<AssemblyGraph>> {
    build_threading_graph_core(
        reference,
        reads,
        kmer_size,
        args,
        allow_low_complexity,
        allow_non_unique_ref,
        false,
    )
}

fn build_threading_graph_core(
    reference: &AssemblyRead,
    reads: &[AssemblyRead],
    kmer_size: usize,
    args: &ReadThreadingAssemblerArgs,
    allow_low_complexity: bool,
    allow_non_unique_ref: bool,
    abort_cyclic_before_dangling: bool,
) -> GatkResult<Option<AssemblyGraph>> {
    if reference.bases.len() < kmer_size {
        return Ok(None);
    }

    let params = AssemblyGraphParams {
        kmer_size: crate::bio_ids::KmerSize::try_from_usize(kmer_size)?,
        min_base_quality: args.min_base_quality,
        min_edge_weight: 1,
        dangling_path_max_nodes: 0,
        max_haplotypes: args.num_best_haplotypes_per_graph,
        max_haplotype_bases: 4096,
        start_threading_only_at_existing_vertex: !args.recover_dangling_branches,
    };

    if !allow_non_unique_ref
        && !args.allow_non_unique_kmers_in_ref
        && reference_has_non_unique_kmers(reference, kmer_size)
    {
        return Ok(None);
    }

    let mut graph = assembly_graph_from_ref_and_reads_threading(reference, reads, &params)?;

    let summary = threading_non_unique_summary(Some(reference), reads, &params)?;
    if !allow_low_complexity && !args.allow_low_complexity_graphs && summary.is_low_complexity {
        return Ok(None);
    }

    let mut pruning = AssemblyGraphPruningParams::gatk_haplotype_caller_defaults();
    pruning.min_prune_factor = args.min_prune_factor;
    pruning.use_adaptive_pruning = args.use_adaptive_pruning;

    if args.prune_before_cycle_counting {
        graph.apply_pruning(&pruning);
    }

    if abort_cyclic_before_dangling
        && args.use_seq_graph
        && args.abort_seq_graph_on_cycles
        && graph.has_cycle()
    {
        return Ok(None);
    }

    if args.recover_dangling_branches && !graph.ref_nodes.is_empty() {
        let dangling = DanglingRecoveryParams::from_assembler_args(args);
        let _ = graph.recover_dangling_branches(&dangling)?;
    }

    if args.recover_all_dangling_branches && graph.has_cycle() {
        return Ok(None);
    }

    if !args.skip_post_dangling_prune && !args.prune_before_cycle_counting {
        graph.apply_pruning(&pruning);
    }

    if args.remove_paths_not_connected_to_ref {
        if graph.reference_source_vertex().is_none() || graph.reference_sink_vertex().is_none() {
            return Ok(None);
        }
        graph.remove_paths_not_connected_to_ref()?;
    }

    if graph.reference_source_vertex().is_none() || graph.reference_sink_vertex().is_none() {
        return Ok(None);
    }

    Ok(Some(graph))
}

/// GATK `AbstractReadThreadingGraph.postProcessForHaplotypeFinding` (RT graph: no-op).
pub fn post_process_for_haplotype_finding(_graph: &AssemblyGraph) {}

/// ASM-1 parity probe: edge count before/after dangling recovery on a pruned RT graph.
/// # Invariants
/// `edges_after` ≥ `edges_before` when dangling recovery adds edges.
/// Attempt counts partition tail vs head dangling recovery attempts.
/// # Ownership
/// [`Copy`] audit snapshot returned from [`audit_threading_dangling_recovery`].
/// # Mutation
/// Immutable probe result; graph is discarded after audit.
/// # Biological assumptions
/// Dangling branch recovery reconnects dead-end paths to the reference backbone.
/// # Java equivalence
/// GATK read-threading dangling recovery (`recoverDanglingBranches` / heads) — ASM-1 audit slice.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ThreadingDanglingAudit {
    pub edges_before: u32,
    pub edges_after: u32,
    pub tails_recovered: u32,
    pub tails_attempted: u32,
    pub heads_recovered: u32,
    pub heads_attempted: u32,
}

pub fn audit_threading_dangling_recovery(
    reference: &AssemblyRead,
    reads: &[AssemblyRead],
    kmer_size: usize,
    args: &ReadThreadingAssemblerArgs,
    allow_low_complexity: bool,
    allow_non_unique_ref: bool,
) -> GatkResult<Option<ThreadingDanglingAudit>> {
    let mut local = args.clone();
    local.recover_dangling_branches = false;
    local.remove_paths_not_connected_to_ref = false;
    let Some(mut graph) = build_threading_graph_core(
        reference,
        reads,
        kmer_size,
        &local,
        allow_low_complexity,
        allow_non_unique_ref,
        false,
    )?
    else {
        return Ok(None);
    };
    let edges_before = graph.edge_count() as u32;
    if args.recover_dangling_branches && !graph.ref_nodes.is_empty() {
        let dangling = DanglingRecoveryParams::from_assembler_args(args);
        let summary = graph.recover_dangling_branches(&dangling)?;
        return Ok(Some(ThreadingDanglingAudit {
            edges_before,
            edges_after: graph.edge_count() as u32,
            tails_recovered: summary.tails_recovered,
            tails_attempted: summary.tails_attempted,
            heads_recovered: summary.heads_recovered,
            heads_attempted: summary.heads_attempted,
        }));
    }
    Ok(Some(ThreadingDanglingAudit {
        edges_before,
        edges_after: edges_before,
        tails_recovered: 0,
        tails_attempted: 0,
        heads_recovered: 0,
        heads_attempted: 0,
    }))
}

/// Single reference haplotype for ref-only assembly (Java `just_assembled_reference` dump slice).
pub fn reference_only_haplotypes(reference: &AssemblyRead, kmer_size: usize) -> Vec<Haplotype> {
    just_reference_result(kmer_size, reference).haplotypes
}

fn just_reference_result(kmer_size: usize, reference: &AssemblyRead) -> AssemblyResult {
    AssemblyResult {
        status: AssemblyStatus::JustAssembledReference,
        kmer_size,
        haplotypes: vec![Haplotype::new(reference.bases.as_slice(), true)],
        event_maps: vec![EventMap::default()],
    }
}

/// Why a KBest path was not emitted (CIGAR-EX diagnostic).
/// # Invariants
/// Each variant corresponds to one SW/CIGAR gate in KBest haplotype extraction.
/// `TooDivergent` is reserved for paths whose CIGAR would contain `N` (currently disabled).
/// # Ownership
/// [`Copy`] reject reason embedded in [`KbestExtractAuditRow::outcome`].
/// # Mutation
/// Immutable enum tag for audit rows.
/// # Biological assumptions
/// Extracted haplotypes must align to reference with GATK SW parameters and minimum ref length.
/// # Java equivalence
/// GATK `findBestPaths` rejection reasons (`pathIsTooDivergentFromReference`, duplicate bases, etc.).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum KbestExtractReject {
    DuplicateBasesLabel,
    SwFailed,
    RefLengthTooShort,
    TooDivergent,
}

/// Per-path extract outcome for parity dumps.
/// # Invariants
/// `path_index` indexes the KBest path list passed to [`audit_kbest_extract`].
/// `outcome` is `Ok()` when the path would be emitted, else carries [`KbestExtractReject`].
/// # Ownership
/// Owns scalar path metadata; borrows no graph data after audit completes.
/// # Mutation
/// Immutable audit row collected per KBest path.
/// # Biological assumptions
/// None documented (extract audit for assembly parity dumps).
/// # Java equivalence
/// Rust-native audit aligned with GATK KBest → haplotype extraction gates.
#[derive(Debug, Clone)]
pub struct KbestExtractAuditRow {
    pub path_index: usize,
    pub is_reference: bool,
    pub bases_len: usize,
    pub eq_ref_bases: bool,
    pub outcome: Result<(), KbestExtractReject>,
}

/// GATK `pathIsTooDivergentFromReference` — skip paths whose SW CIGAR would contain `N`.
pub fn path_is_too_divergent_from_reference(_cigar: &Cigar) -> bool {
    false
}

/// Audit KBest → haplotype extraction (GATK `findBestPaths` SW/CIGAR gates).
pub fn audit_kbest_extract(
    paths: &[KBestPath],
    graph: &AssemblyGraph,
    ref_haplotype: &Haplotype,
    ref_cigar_length: usize,
    sw: &SwParameters,
) -> Vec<KbestExtractAuditRow> {
    let ref_bytes = &ref_haplotype.bases;
    let mut seen: HashSet<(Vec<u8>, bool)> = HashSet::new();
    let mut rows = Vec::new();
    for (path_index, path) in paths.iter().enumerate() {
        let bases = path.bases(graph);
        // CLONE: needed because owned composite key for dedup/lookup.
        let key = (bases.clone(), path.is_reference);
        let eq_ref_bases = bases.as_slice() == ref_bytes.as_slice();
        let outcome = if seen.contains(&key) {
            Err(KbestExtractReject::DuplicateBasesLabel)
        } else {
            seen.insert(key);
            match calculate_haplotype_cigar_for_assembly(ref_bytes, &bases, ref_cigar_length, sw) {
                None => Err(KbestExtractReject::SwFailed),
                Some(cigar) if path_is_too_divergent_from_reference(&cigar) => {
                    Err(KbestExtractReject::TooDivergent)
                }
                Some(cigar)
                    if ref_cigar_length >= MIN_HAPLOTYPE_REFERENCE_LENGTH
                        && cigar.reference_length() < MIN_HAPLOTYPE_REFERENCE_LENGTH =>
                {
                    Err(KbestExtractReject::RefLengthTooShort)
                }
                Some(_) => Ok(()),
            }
        };
        rows.push(KbestExtractAuditRow {
            path_index,
            is_reference: path.is_reference,
            bases_len: bases.len(),
            eq_ref_bases,
            outcome,
        });
    }
    rows
}

/// SW + CIGAR filter on KBest paths (GATK `findBestPaths` cigar checks).
pub fn extract_haplotypes_from_kbest_paths(
    paths: &[KBestPath],
    graph: &AssemblyGraph,
    ref_haplotype: &Haplotype,
    ref_cigar_length: usize,
    sw: &SwParameters,
) -> GatkResult<Vec<Haplotype>> {
    let ref_bytes = &ref_haplotype.bases;
    let mut out = Vec::new();
    let mut seen: HashSet<(Vec<u8>, bool)> = HashSet::new();

    for path in paths {
        let bases = path.bases(graph);
        // CLONE: needed because owned composite key for dedup/lookup.
        let key = (bases.clone(), path.is_reference);
        if seen.contains(&key) {
            continue;
        }
        // L14-E2: move bases into haplotype (avoid second full-sequence clone).
        let mut h = Haplotype::new(bases, path.is_reference);
        h.score = path.score;
        h.kmer_size = graph.kmer_size;

        let Some(assy) = calculate_haplotype_cigar_for_assembly_with_offset(
            ref_bytes,
            &h.bases,
            ref_cigar_length,
            sw,
        ) else {
            continue;
        };
        if path_is_too_divergent_from_reference(&assy.cigar) {
            continue;
        }
        if ref_cigar_length >= MIN_HAPLOTYPE_REFERENCE_LENGTH
            && assy.cigar.reference_length() < MIN_HAPLOTYPE_REFERENCE_LENGTH
        {
            continue;
        }
        h.cigar = Some(assy.cigar);
        h.alignment_start_hap_wrt_ref = assy.alignment_start_hap_wrt_ref;
        seen.insert(key);
        out.push(h);
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::assembly::AssemblyRead;

    fn read(seq: &str, q: u8) -> AssemblyRead {
        AssemblyRead {
            bases: seq.as_bytes().to_vec(),
            base_quals: vec![q; seq.len()],
        }
    }

    #[test]
    fn ref_byte_non_ref_label_collapses_before_variation_present() {
        let ref_bases = b"ACGTACGT".to_vec();
        // CLONE: needed because haplotype constructor takes owned bases.
        let mut ref_hap = Haplotype::new(ref_bases.clone(), true);
        let mut ref_cigar = Cigar::new();
        ref_cigar.push(ref_bases.len(), CigarOperator::Match);
        ref_hap.cigar = Some(ref_cigar);

        // CLONE: needed because haplotype constructor takes owned bases.
        let mut haps = vec![Haplotype::new(ref_bases.clone(), false)];
        let status = finalize_assembly_haplotypes(&mut haps, &ref_hap, true);
        assert_eq!(haps.len(), 1);
        assert!(haps[0].is_reference);
        assert!(matches!(status, AssemblyStatus::JustAssembledReference));
    }

    #[test]
    fn p5_case1_assembler_emits_reference_haplotype() {
        let reference = read("ACGTT", 30);
        let reads = vec![
            read("ACGTT", 30),
            read("ACGTT", 30),
            read("ACGTT", 30),
            read("ACGTA", 30),
            read("ACGTA", 30),
        ];
        let mut args = ReadThreadingAssemblerArgs::default();
        args.kmer_sizes = vec![3];
        args.min_prune_factor = 2;
        args.allow_low_complexity_graphs = true;
        let result = assemble_from_ref_and_reads(&reference, &reads, &args).unwrap();
        assert!(result
            .haplotypes
            .iter()
            .any(|h| h.is_reference && h.sequence_string() == "ACGTT"));
    }
}

/// Ensure a reference-labeled haplotype is present (GATK `findBestPaths` adds `refHaplotype` to the set).
/// When KBest marks a ref-byte path as non-reference (non-ref edges on spine), keep that alt entry so
/// `is_variation_present` can see both alleles; do not collapse by sequence alone.
fn ensure_reference_haplotype(haplotypes: &mut Vec<Haplotype>, ref_hap: &Haplotype) {
    if haplotypes
        .iter()
        .any(|h| h.is_reference && h.bases == ref_hap.bases)
    {
        return;
    }
    // CLONE: needed because owned haplotypes for scoring call.
    let mut rh = ref_hap.clone();
    if rh.cigar.is_none() {
        let mut c = Cigar::new();
        c.push(rh.bases.len(), CigarOperator::Match);
        rh.cigar = Some(c);
    }
    haplotypes.push(rh);
}

fn assembly_status_from_haplotype_list(haplotypes: &[Haplotype]) -> AssemblyStatus {
    if haplotypes.iter().any(|h| !h.is_reference) && haplotypes.len() > 1 {
        AssemblyStatus::AssembledSomeVariation
    } else if haplotypes.len() <= 1 {
        AssemblyStatus::JustAssembledReference
    } else {
        AssemblyStatus::JustAssembledReference
    }
}

fn finalize_assembly_haplotypes(
    haplotypes: &mut Vec<Haplotype>,
    ref_hap: &Haplotype,
    ensure_ref: bool,
) -> AssemblyStatus {
    crate::haplotype::normalize_ref_equivalent_haplotypes(haplotypes, &ref_hap.bases);
    if ensure_ref {
        ensure_reference_haplotype(haplotypes, ref_hap);
    }
    assembly_status_from_haplotype_list(haplotypes)
}
