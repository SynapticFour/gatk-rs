//! Assembly graph edge dumps for L2 parity.

use crate::alignment::{calculate_haplotype_cigar, SwParameters};
use crate::assembly::{
    AssemblyGraph, AssemblyGraphParams, AssemblyGraphPruningParams, AssemblyGraphSummary,
    AssemblyRead,
};
use crate::assembly_dangling_recovery::{DanglingRecoveryParams, DanglingRecoverySummary};
use crate::haplotype::Haplotype;
use crate::junction_kbest::find_junction_best_haplotypes;
use crate::junction_tree_graph::build_junction_tree_graph_from_ref_and_reads;
use crate::kbest_haplotype::find_best_haplotypes;
use crate::read_error_correction::{
    correct_reads_pileup_log_odds, load_aligned_assembly_reads_tsv,
};
use crate::read_threading_assembler::{
    assemble_from_ref_and_reads, build_threading_graph_for_haplotype_dump,
    extract_haplotypes_from_kbest_paths, ReadThreadingAssemblerArgs,
};
use crate::read_threading_graph::{
    assembly_graph_from_ref_and_reads_threading, threading_non_unique_summary,
    ThreadingNonUniqueSummary,
};
use crate::seq_graph::{SeqGraph, SeqGraphCleanupStatus};
use gatk_common::{GatkError, GatkResult};
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::Path;

/// Load a single reference row (`sequence\\tqual`) for ref-first threading.
pub fn load_assembly_ref_tsv(path: &Path) -> GatkResult<AssemblyRead> {
    let rows = load_assembly_reads_tsv(path)?;
    rows.into_iter()
        .next()
        .ok_or_else(|| GatkError::argument(format!("ref tsv {}: no sequence row", path.display())))
}

/// Load `sequence\\tmean_qual` rows (same schema as `parity/fixtures/p5_assembly_case1_reads.tsv`).
pub fn load_assembly_reads_tsv(path: &Path) -> GatkResult<Vec<AssemblyRead>> {
    let f = File::open(path)
        .map_err(|e| GatkError::generic(format!("open {}: {e}", path.display())))?;
    let reader = BufReader::new(f);
    let mut out = Vec::new();
    for line in reader.lines() {
        let line = line.map_err(|e| GatkError::generic(format!("read {}: {e}", path.display())))?;
        let t = line.trim();
        if t.is_empty() || t.starts_with('#') {
            continue;
        }
        let mut parts = t.split_whitespace();
        let bases = parts
            .next()
            .ok_or_else(|| GatkError::argument("reads tsv: missing sequence"))?
            .as_bytes()
            .to_vec();
        let q = parts
            .next()
            .ok_or_else(|| GatkError::argument("reads tsv: missing qual"))?
            .parse::<u8>()
            .map_err(|_| GatkError::argument("reads tsv: invalid qual"))?;
        let n = bases.len();
        out.push(AssemblyRead {
            bases,
            base_quals: vec![q; n],
        });
    }
    Ok(out)
}

/// Write sorted k-mer edge rows for each k: `kmer_size\\tfrom_kmer\\tto_kmer\\tsupport`.
pub fn dump_assembly_graph_multi_kmer_edges_tsv(
    reads_path: &Path,
    kmer_sizes: &[usize],
    min_base_quality: u8,
    out: &mut impl Write,
) -> GatkResult<()> {
    let reads = load_assembly_reads_tsv(reads_path)?;
    let graphs = AssemblyGraph::from_reads_kmer_sizes(&reads, kmer_sizes, min_base_quality)?;
    writeln!(out, "kmer_size\tfrom_kmer\tto_kmer\tsupport")
        .map_err(|e| GatkError::generic(format!("write tsv: {e}")))?;
    for graph in &graphs {
        for e in graph.edges_sorted() {
            let from = String::from_utf8_lossy(&graph.nodes()[e.from].kmer);
            let to = String::from_utf8_lossy(&graph.nodes()[e.to].kmer);
            writeln!(out, "{}\t{from}\t{to}\t{}", graph.kmer_size, e.support)
                .map_err(|e| GatkError::generic(format!("write tsv: {e}")))?;
        }
    }
    Ok(())
}

fn format_summary_f64(v: f64) -> String {
    if v == f64::NEG_INFINITY {
        return "-inf".to_string();
    }
    if v == f64::INFINITY {
        return "inf".to_string();
    }
    format!("{:.8}", v)
}

/// Write post-pruning graph summary (`metric\\tvalue`) for.
pub fn dump_assembly_graph_pruned_summary_tsv(
    reads_path: &Path,
    kmer_size: usize,
    min_base_quality: u8,
    pruning: &AssemblyGraphPruningParams,
    out: &mut impl Write,
) -> GatkResult<()> {
    let reads = load_assembly_reads_tsv(reads_path)?;
    let params = AssemblyGraphParams {
        kmer_size: crate::bio_ids::KmerSize::try_from_usize(kmer_size)?,
        min_base_quality,
        min_edge_weight: 1,
        dangling_path_max_nodes: 0,
        max_haplotypes: 128,
        max_haplotype_bases: 512,
        start_threading_only_at_existing_vertex: false,
    };
    let mut graph = AssemblyGraph::from_reads(&reads, &params)?;
    let summary = graph.summarize_after_pruning(pruning);
    write_assembly_graph_summary_tsv(&summary, out)
}

pub fn write_dangling_recovery_summary_tsv(
    summary: &DanglingRecoverySummary,
    out: &mut impl Write,
) -> GatkResult<()> {
    writeln!(out, "metric\tvalue").map_err(|e| GatkError::generic(format!("write tsv: {e}")))?;
    let rows = [
        ("edges_before", summary.edges_before.to_string()),
        ("edges_after", summary.edges_after.to_string()),
        ("tails_attempted", summary.tails_attempted.to_string()),
        ("tails_recovered", summary.tails_recovered.to_string()),
        ("heads_attempted", summary.heads_attempted.to_string()),
        ("heads_recovered", summary.heads_recovered.to_string()),
        ("edges_merged", summary.edges_merged.to_string()),
    ];
    for (metric, value) in rows {
        writeln!(out, "{metric}\t{value}")
            .map_err(|e| GatkError::generic(format!("write tsv: {e}")))?;
    }
    Ok(())
}

pub fn write_threading_non_unique_summary_tsv(
    summary: &ThreadingNonUniqueSummary,
    out: &mut impl Write,
) -> GatkResult<()> {
    writeln!(out, "metric\tvalue").map_err(|e| GatkError::generic(format!("write tsv: {e}")))?;
    let rows = [
        ("node_count", summary.node_count.to_string()),
        ("edge_count", summary.edge_count.to_string()),
        ("unique_kmer_count", summary.unique_kmer_count.to_string()),
        (
            "non_unique_kmer_count",
            summary.non_unique_kmer_count.to_string(),
        ),
        (
            "is_low_complexity",
            if summary.is_low_complexity {
                "true".to_string()
            } else {
                "false".to_string()
            },
        ),
        (
            "max_kmer_multiplicity",
            summary.max_kmer_multiplicity.to_string(),
        ),
    ];
    for (metric, value) in rows {
        writeln!(out, "{metric}\t{value}")
            .map_err(|e| GatkError::generic(format!("write tsv: {e}")))?;
    }
    Ok(())
}

/// Read threading non-unique kmer / cycle policy summary.
/// Pass `ref_path = None` or use `-` in the CLI for reads-only fixtures.
pub fn dump_assembly_graph_non_unique_summary_tsv(
    ref_path: Option<&Path>,
    reads_path: &Path,
    kmer_size: usize,
    min_base_quality: u8,
    out: &mut impl Write,
) -> GatkResult<()> {
    let reference = ref_path.map(load_assembly_ref_tsv).transpose()?;
    let reads = load_assembly_reads_tsv(reads_path)?;
    let params = AssemblyGraphParams {
        kmer_size: crate::bio_ids::KmerSize::try_from_usize(kmer_size)?,
        min_base_quality,
        min_edge_weight: 1,
        dangling_path_max_nodes: 0,
        max_haplotypes: 128,
        max_haplotype_bases: 512,
        start_threading_only_at_existing_vertex: false,
    };
    let summary = threading_non_unique_summary(reference.as_ref(), &reads, &params)?;
    write_threading_non_unique_summary_tsv(&summary, out)
}

/// Ref-first threading, GATK pruning, then dangling recovery summary.
pub fn dump_assembly_graph_dangling_summary_tsv(
    ref_path: &Path,
    reads_path: &Path,
    kmer_size: usize,
    min_base_quality: u8,
    min_prune_factor: u32,
    min_dangling_branch_length: usize,
    recover_dangling_heads: bool,
    out: &mut impl Write,
) -> GatkResult<()> {
    let reference = load_assembly_ref_tsv(ref_path)?;
    let reads = load_assembly_reads_tsv(reads_path)?;
    let params = AssemblyGraphParams {
        kmer_size: crate::bio_ids::KmerSize::try_from_usize(kmer_size)?,
        min_base_quality,
        min_edge_weight: 1,
        dangling_path_max_nodes: 0,
        max_haplotypes: 128,
        max_haplotype_bases: 512,
        start_threading_only_at_existing_vertex: false,
    };
    let mut graph = assembly_graph_from_ref_and_reads_threading(&reference, &reads, &params)?;
    let mut pruning = AssemblyGraphPruningParams::gatk_haplotype_caller_defaults();
    pruning.min_prune_factor = min_prune_factor;
    graph.apply_pruning(&pruning);
    let mut dangling = DanglingRecoveryParams::gatk_haplotype_caller_defaults();
    dangling.min_dangling_branch_length = min_dangling_branch_length;
    dangling.recover_dangling_heads = recover_dangling_heads;
    // L2 E.4 compares to pinned Java `recoverDanglingTails` (single pass, no ASM-1 rescue).
    dangling.dangling_java_exact = true;
    let summary = graph.recover_dangling_branches(&dangling)?;
    write_dangling_recovery_summary_tsv(&summary, out)
}

pub fn write_assembly_graph_summary_tsv(
    summary: &AssemblyGraphSummary,
    out: &mut impl Write,
) -> GatkResult<()> {
    writeln!(out, "metric\tvalue").map_err(|e| GatkError::generic(format!("write tsv: {e}")))?;
    let rows = [
        ("node_count", summary.node_count.to_string()),
        ("edge_count", summary.edge_count.to_string()),
        (
            "log10_max_edge_support",
            format_summary_f64(summary.log10_max_edge_support),
        ),
        (
            "log10_sum_edge_support",
            format_summary_f64(summary.log10_sum_edge_support),
        ),
        (
            "pruning_lod_threshold_ln",
            format_summary_f64(summary.pruning_lod_threshold_ln),
        ),
        (
            "adaptive_pruning",
            if summary.adaptive_pruning {
                "true".to_string()
            } else {
                "false".to_string()
            },
        ),
        ("min_prune_factor", summary.min_prune_factor.to_string()),
        ("edges_pruned", summary.edges_pruned.to_string()),
    ];
    for (metric, value) in rows {
        writeln!(out, "{metric}\t{value}")
            .map_err(|e| GatkError::generic(format!("write tsv: {e}")))?;
    }
    Ok(())
}

/// Write sorted k-mer edge multiset: `from_kmer\\tto_kmer\\tsupport`.
pub fn dump_assembly_graph_edges_tsv(
    reads_path: &Path,
    kmer_size: usize,
    min_base_quality: u8,
    out: &mut impl Write,
) -> GatkResult<()> {
    let reads = load_assembly_reads_tsv(reads_path)?;
    let params = AssemblyGraphParams {
        kmer_size: crate::bio_ids::KmerSize::try_from_usize(kmer_size)?,
        min_base_quality,
        min_edge_weight: 1,
        dangling_path_max_nodes: 0,
        max_haplotypes: 128,
        max_haplotype_bases: 512,
        start_threading_only_at_existing_vertex: false,
    };
    let graph = AssemblyGraph::from_reads(&reads, &params)?;
    writeln!(out, "from_kmer\tto_kmer\tsupport")
        .map_err(|e| GatkError::generic(format!("write tsv: {e}")))?;
    for e in graph.edges_sorted() {
        let from = String::from_utf8_lossy(&graph.nodes()[e.from].kmer);
        let to = String::from_utf8_lossy(&graph.nodes()[e.to].kmer);
        writeln!(out, "{from}\t{to}\t{}", e.support)
            .map_err(|e| GatkError::generic(format!("write tsv: {e}")))?;
    }
    Ok(())
}

/// Smith–Waterman haplotype-to-reference CIGARs (`CigarUtils.calculateCigar`).
pub fn dump_assembly_haplotype_cigars_tsv(
    ref_path: &Path,
    haplotypes_path: &Path,
    out: &mut impl Write,
) -> GatkResult<()> {
    let reference = load_assembly_ref_tsv(ref_path)?;
    let haplotypes = load_assembly_reads_tsv(haplotypes_path)?;
    let sw = SwParameters::gatk_haplotype_to_reference();
    writeln!(out, "haplotype_idx\tsequence\tcigar")
        .map_err(|e| GatkError::generic(format!("write tsv: {e}")))?;
    for (idx, hap) in haplotypes.iter().enumerate() {
        let cigar =
            calculate_haplotype_cigar(reference.bases.as_slice(), hap.bases.as_slice(), &sw)
                .map(|c| c.to_gatk_string())
                .unwrap_or_default();
        writeln!(
            out,
            "{idx}\t{}\t{cigar}",
            String::from_utf8_lossy(&hap.bases)
        )
        .map_err(|e| GatkError::generic(format!("write tsv: {e}")))?;
    }
    Ok(())
}

fn format_score(v: f64) -> String {
    if v == 0.0 {
        "0".to_string()
    } else {
        format!("{:.8}", v)
    }
}

fn haplotype_dump_args(
    kmer_size: usize,
    min_base_quality: u8,
    min_prune: u32,
    min_dangling_branch_length: usize,
    recover_heads: bool,
    max_haplotypes: usize,
) -> ReadThreadingAssemblerArgs {
    let mut args = ReadThreadingAssemblerArgs::default();
    args.kmer_sizes = vec![kmer_size];
    args.min_base_quality = min_base_quality;
    args.min_prune_factor = min_prune;
    args.min_dangling_branch_length = min_dangling_branch_length;
    args.recover_dangling_heads = recover_heads;
    args.recover_dangling_branches = true;
    args.allow_low_complexity_graphs = true;
    args.num_best_haplotypes_per_graph = max_haplotypes;
    args
}

/// Sort order used by `HcFullParityGateDump.assemblyRegionHaplotypes` (score desc, sequence desc).
pub(crate) fn sort_haplotypes_java_dump_order(haplotypes: &mut [Haplotype]) {
    haplotypes.sort_by(|a, b| {
        b.score
            .partial_cmp(&a.score)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| b.sequence_string().cmp(&a.sequence_string()))
    });
}

pub(crate) fn write_haplotype_rows_with_ref_recovery(
    haplotypes: &[Haplotype],
    ref_bases: &[u8],
    ref_sequence: &[u8],
    out: &mut impl Write,
) -> GatkResult<()> {
    writeln!(out, "rank\tsequence\tscore\tis_reference\tcigar")
        .map_err(|e| GatkError::generic(format!("write tsv: {e}")))?;
    let mut owned: Vec<Haplotype> = haplotypes.to_vec();
    sort_haplotypes_java_dump_order(&mut owned);
    let ranked: Vec<_> = owned.iter().collect();
    let mut rank = 0usize;
    let mut wrote_ref_sequence = false;
    for h in ranked {
        if h.is_reference || h.bases == ref_bases {
            wrote_ref_sequence = true;
        }
        let cigar = h
            .cigar
            .as_ref()
            .map(|c| c.to_gatk_string())
            .unwrap_or_default();
        writeln!(
            out,
            "{rank}\t{seq}\t{score}\t{is_ref}\t{cigar}",
            seq = h.sequence_string(),
            score = format_score(h.score),
            is_ref = h.is_reference,
            cigar = cigar
        )
        .map_err(|e| GatkError::generic(format!("write tsv: {e}")))?;
        rank += 1;
    }
    if !wrote_ref_sequence {
        let sw = SwParameters::gatk_haplotype_to_reference();
        let cigar = calculate_haplotype_cigar(ref_bases, ref_bases, &sw)
            .map(|c| c.to_gatk_string())
            .unwrap_or_default();
        writeln!(
            out,
            "{rank}\t{seq}\t0\ttrue\t{cigar}",
            seq = String::from_utf8_lossy(ref_sequence),
            cigar = cigar
        )
        .map_err(|e| GatkError::generic(format!("write tsv: {e}")))?;
    }
    Ok(())
}

/// Raw `GraphBasedKBestHaplotypeFinder` paths (no SW/CIGAR filter).
pub fn dump_assembly_kbest_paths_tsv(
    ref_path: &Path,
    reads_path: &Path,
    kmer_size: usize,
    min_base_quality: u8,
    min_prune: u32,
    min_dangling_branch_length: usize,
    recover_heads: bool,
    max_haplotypes: usize,
    out: &mut impl Write,
) -> GatkResult<()> {
    let reference = load_assembly_ref_tsv(ref_path)?;
    let reads = load_assembly_reads_tsv(reads_path)?;
    let args = haplotype_dump_args(
        kmer_size,
        min_base_quality,
        min_prune,
        min_dangling_branch_length,
        recover_heads,
        max_haplotypes,
    );
    let Some(graph) = build_threading_graph_for_haplotype_dump(
        &reference, &reads, kmer_size, &args, true, false,
    )?
    else {
        writeln!(out, "rank\tsequence\tscore\tis_reference")
            .map_err(|e| GatkError::generic(format!("write tsv: {e}")))?;
        return Ok(());
    };
    let paths = find_best_haplotypes(&graph, max_haplotypes)?;
    writeln!(out, "rank\tsequence\tscore\tis_reference")
        .map_err(|e| GatkError::generic(format!("write tsv: {e}")))?;
    for (rank, path) in paths.iter().enumerate() {
        writeln!(
            out,
            "{rank}\t{seq}\t{score}\t{is_ref}",
            seq = String::from_utf8_lossy(&path.bases(&graph)),
            score = format_score(path.score),
            is_ref = path.is_reference
        )
        .map_err(|e| GatkError::generic(format!("write tsv: {e}")))?;
    }
    Ok(())
}

/// / SeqGraph: post-`toSequenceGraph` + `cleanupSeqGraph` summary.
pub fn dump_assembly_seqgraph_summary_tsv(
    ref_path: &Path,
    reads_path: &Path,
    kmer_size: usize,
    min_base_quality: u8,
    min_prune: u32,
    min_dangling_branch_length: usize,
    recover_heads: bool,
    out: &mut impl Write,
) -> GatkResult<()> {
    let reference = load_assembly_ref_tsv(ref_path)?;
    let reads = load_assembly_reads_tsv(reads_path)?;
    let args = haplotype_dump_args(
        kmer_size,
        min_base_quality,
        min_prune,
        min_dangling_branch_length,
        recover_heads,
        128,
    );
    writeln!(out, "metric\tvalue").map_err(|e| GatkError::generic(format!("write tsv: {e}")))?;
    let Some(graph) = build_threading_graph_for_haplotype_dump(
        &reference, &reads, kmer_size, &args, true, false,
    )?
    else {
        write_seqgraph_metric(out, "status", "no_graph")?;
        write_seqgraph_metric(out, "node_count", "0")?;
        write_seqgraph_metric(out, "edge_count", "0")?;
        write_seqgraph_metric(out, "ref_path_len", "0")?;
        write_seqgraph_metric(out, "ref_path_sequence", "")?;
        return Ok(());
    };
    let mut seq = SeqGraph::from_assembly_graph(&graph);
    seq.clean_non_ref_paths();
    let status = seq.cleanup_seq_graph();
    let status_s = match status {
        SeqGraphCleanupStatus::AssembledSomeVariation => "assembled_some_variation",
        SeqGraphCleanupStatus::JustAssembledReference => "just_assembled_reference",
    };
    let ref_path = seq
        .reference_path_bytes()
        .map(|b| String::from_utf8_lossy(&b).into_owned())
        .unwrap_or_default();
    write_seqgraph_metric(out, "status", status_s)?;
    write_seqgraph_metric(out, "node_count", &seq.node_count().to_string())?;
    write_seqgraph_metric(out, "edge_count", &seq.edge_count().to_string())?;
    write_seqgraph_metric(out, "ref_path_len", &ref_path.len().to_string())?;
    write_seqgraph_metric(out, "ref_path_sequence", &ref_path)?;
    Ok(())
}

fn write_seqgraph_metric(out: &mut impl Write, metric: &str, value: &str) -> GatkResult<()> {
    writeln!(out, "{metric}\t{value}").map_err(|e| GatkError::generic(format!("write tsv: {e}")))
}

/// `ReadThreadingGraph.isLowQualityGraph` / assembler abort predicate.
pub fn dump_assembly_graph_low_quality_tsv(
    ref_path: &Path,
    reads_path: &Path,
    kmer_size: usize,
    out: &mut impl Write,
) -> GatkResult<()> {
    let reference = load_assembly_ref_tsv(ref_path)?;
    let reads = load_assembly_reads_tsv(reads_path)?;
    let args = ReadThreadingAssemblerArgs::default();
    let params = AssemblyGraphParams {
        kmer_size: crate::bio_ids::KmerSize::try_from_usize(kmer_size)?,
        min_base_quality: args.min_base_quality,
        min_edge_weight: 1,
        dangling_path_max_nodes: 0,
        max_haplotypes: args.num_best_haplotypes_per_graph,
        max_haplotype_bases: 4096,
        start_threading_only_at_existing_vertex: !args.recover_dangling_branches,
    };
    let summary = threading_non_unique_summary(Some(&reference), &reads, &params)?;
    let graph = build_threading_graph_for_haplotype_dump(
        &reference, &reads, kmer_size, &args, false, false,
    )?;
    let would_abort =
        graph.is_none() || (!args.allow_low_complexity_graphs && summary.is_low_complexity);
    writeln!(out, "kmer_size\t{kmer_size}")
        .map_err(|e| GatkError::generic(format!("write tsv: {e}")))?;
    writeln!(out, "is_low_quality_graph\t{}", summary.is_low_complexity)
        .map_err(|e| GatkError::generic(format!("write tsv: {e}")))?;
    writeln!(out, "unique_kmer_count\t{}", summary.unique_kmer_count)
        .map_err(|e| GatkError::generic(format!("write tsv: {e}")))?;
    writeln!(
        out,
        "non_unique_kmer_count\t{}",
        summary.non_unique_kmer_count
    )
    .map_err(|e| GatkError::generic(format!("write tsv: {e}")))?;
    writeln!(out, "would_abort_assembly\t{would_abort}")
        .map_err(|e| GatkError::generic(format!("write tsv: {e}")))?;
    Ok(())
}

/// Full `ReadThreadingAssembler` multi-kmer assemble (`assemble_from_ref_and_reads`).
pub fn dump_assembly_assemble_tsv(
    ref_path: &Path,
    reads_path: &Path,
    out: &mut impl Write,
) -> GatkResult<()> {
    let reference = load_assembly_ref_tsv(ref_path)?;
    let reads = load_assembly_reads_tsv(reads_path)?;
    let result = crate::read_threading_assembler::assemble_for_java_gate_dump(&reference, &reads)?;
    let status_s = match result.status {
        crate::read_threading_assembler::AssemblyStatus::AssembledSomeVariation => {
            "assembled_some_variation"
        }
        crate::read_threading_assembler::AssemblyStatus::JustAssembledReference => {
            "just_assembled_reference"
        }
        crate::read_threading_assembler::AssemblyStatus::Failed => "failed",
    };
    writeln!(out, "status\t{status_s}")
        .map_err(|e| GatkError::generic(format!("write tsv: {e}")))?;
    writeln!(out, "kmer_size\t{}", result.kmer_size)
        .map_err(|e| GatkError::generic(format!("write tsv: {e}")))?;
    let ref_bases = reference.bases.as_slice();
    write_haplotype_rows_with_ref_recovery(&result.haplotypes, ref_bases, &reference.bases, out)
}

/// `PileupReadErrorCorrector` corrected read sequences.
pub fn dump_read_error_correction_tsv(
    reads_path: &Path,
    log_odds_threshold: f64,
    out: &mut impl Write,
) -> GatkResult<()> {
    let mut reads = load_aligned_assembly_reads_tsv(reads_path)?;
    correct_reads_pileup_log_odds(&mut reads, log_odds_threshold)?;
    writeln!(out, "read_index\tsequence\tmean_qual")
        .map_err(|e| GatkError::generic(format!("write tsv: {e}")))?;
    for (i, read) in reads.iter().enumerate() {
        let mean_q = if read.base_quals.is_empty() {
            0
        } else {
            read.base_quals.iter().map(|&q| q as u32).sum::<u32>() / read.base_quals.len() as u32
        };
        writeln!(
            out,
            "{i}\t{seq}\t{mean_q}",
            seq = String::from_utf8_lossy(&read.bases),
            mean_q = mean_q
        )
        .map_err(|e| GatkError::generic(format!("write tsv: {e}")))?;
    }
    Ok(())
}

/// `JunctionTreeKBestHaplotypeFinder` haplotypes vs Java goldens.
pub fn dump_assembly_junction_haplotypes_tsv(
    ref_path: &Path,
    reads_path: &Path,
    kmer_size: usize,
    min_base_quality: u8,
    recover_edges: bool,
    max_haplotypes: usize,
    out: &mut impl Write,
) -> GatkResult<()> {
    let reference = load_assembly_ref_tsv(ref_path)?;
    let reads = load_assembly_reads_tsv(reads_path)?;
    let jt = build_junction_tree_graph_from_ref_and_reads(
        &reference,
        &reads,
        kmer_size,
        min_base_quality,
    )?;
    let graph = &jt.graph;
    if jt.reference_source().is_none() || jt.reference_sink().is_none() {
        writeln!(out, "rank\tsequence\tscore\tis_reference")
            .map_err(|e| GatkError::generic(format!("write tsv: {e}")))?;
        return Ok(());
    }
    let mut paths = find_junction_best_haplotypes(&jt, max_haplotypes, 1, recover_edges)?;
    paths.sort_by(|a, b| {
        b.score
            .partial_cmp(&a.score)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| b.bases(graph).cmp(&a.bases(graph)))
    });
    writeln!(out, "rank\tsequence\tscore\tis_reference")
        .map_err(|e| GatkError::generic(format!("write tsv: {e}")))?;
    for (rank, path) in paths.iter().enumerate() {
        writeln!(
            out,
            "{rank}\t{seq}\t{score}\t{is_ref}",
            seq = String::from_utf8_lossy(&path.bases(graph)),
            score = format_score(path.score),
            is_ref = path.is_reference
        )
        .map_err(|e| GatkError::generic(format!("write tsv: {e}")))?;
    }
    Ok(())
}

fn dump_assembly_haplotypes_with_max(
    ref_path: &Path,
    reads_path: &Path,
    kmer_size: usize,
    min_base_quality: u8,
    min_prune: u32,
    min_dangling_branch_length: usize,
    recover_heads: bool,
    max_haplotypes: usize,
    out: &mut impl Write,
) -> GatkResult<()> {
    let reference = load_assembly_ref_tsv(ref_path)?;
    let reads = load_assembly_reads_tsv(reads_path)?;
    let args = haplotype_dump_args(
        kmer_size,
        min_base_quality,
        min_prune,
        min_dangling_branch_length,
        recover_heads,
        max_haplotypes,
    );
    let ref_bases = reference.bases.as_slice().to_vec();
    let Some(graph) = build_threading_graph_for_haplotype_dump(
        &reference, &reads, kmer_size, &args, true, false,
    )?
    else {
        return write_haplotype_rows_with_ref_recovery(
            // CLONE: needed because haplotype constructor takes owned bases.
            &[Haplotype::new(ref_bases.clone(), true)],
            &ref_bases,
            &reference.bases,
            out,
        );
    };
    // CLONE: needed because haplotype constructor takes owned bases.
    let mut ref_hap = Haplotype::new(ref_bases.clone(), true);
    let mut ref_cigar = crate::cigar::Cigar::new();
    ref_cigar.push(ref_hap.bases.len(), crate::cigar::CigarOperator::Match);
    ref_hap.cigar = Some(ref_cigar);
    let ref_cigar_len = ref_hap.cigar.as_ref().unwrap().reference_length();
    let paths = find_best_haplotypes(&graph, max_haplotypes)?;
    let haplotypes = extract_haplotypes_from_kbest_paths(
        &paths,
        &graph,
        &ref_hap,
        ref_cigar_len,
        &args.haplotype_to_reference_sw,
    )?;
    write_haplotype_rows_with_ref_recovery(&haplotypes, &ref_bases, &reference.bases, out)
}

/// Full assembler haplotype set vs reference (`ReadThreadingAssembler` + KBest).
pub fn dump_assembly_haplotypes_tsv(
    ref_path: &Path,
    reads_path: &Path,
    kmer_size: usize,
    min_base_quality: u8,
    min_prune: u32,
    min_dangling_branch_length: usize,
    recover_heads: bool,
    out: &mut impl Write,
) -> GatkResult<()> {
    dump_assembly_haplotypes_with_max(
        ref_path,
        reads_path,
        kmer_size,
        min_base_quality,
        min_prune,
        min_dangling_branch_length,
        recover_heads,
        128,
        out,
    )
}

/// Haplotype dump via production `ReadThreadingAssembler` (`ensure_reference_in_result`).
pub fn dump_assembly_haplotypes_production_tsv(
    ref_path: &Path,
    reads_path: &Path,
    kmer_size: usize,
    min_base_quality: u8,
    min_prune: u32,
    min_dangling_branch_length: usize,
    recover_heads: bool,
    out: &mut impl Write,
) -> GatkResult<()> {
    let reference = load_assembly_ref_tsv(ref_path)?;
    let reads = load_assembly_reads_tsv(reads_path)?;
    let mut args = haplotype_dump_args(
        kmer_size,
        min_base_quality,
        min_prune,
        min_dangling_branch_length,
        recover_heads,
        128,
    );
    args.ensure_reference_in_result = true;
    let result = assemble_from_ref_and_reads(&reference, &reads, &args)?;
    let ref_bases = reference.bases.as_slice();
    write_haplotype_rows_with_ref_recovery(&result.haplotypes, ref_bases, &reference.bases, out)
}

/// Haplotype dump with explicit `maxNumHaplotypesInPopulation` cap.
pub fn dump_assembly_haplotypes_cap_tsv(
    ref_path: &Path,
    reads_path: &Path,
    kmer_size: usize,
    min_base_quality: u8,
    min_prune: u32,
    min_dangling_branch_length: usize,
    recover_heads: bool,
    max_haplotypes: usize,
    out: &mut impl Write,
) -> GatkResult<()> {
    dump_assembly_haplotypes_with_max(
        ref_path,
        reads_path,
        kmer_size,
        min_base_quality,
        min_prune,
        min_dangling_branch_length,
        recover_heads,
        max_haplotypes,
        out,
    )
}
