//! P12 ASM-1: staged assembly metrics (threading graph → SeqGraph cleanup).

use crate::assembly::AssemblyGraph;
use crate::assembly_based_caller::AssembleReadsArgs;
use crate::assembly_region_finalize::{
    assembly_reads_for_java_materialize_dump, assembly_reads_for_production,
    assembly_reference_read, gatk_min_tail_quality_for_assembly, records_to_assembly_reads,
};

fn parity_min_tail_quality(default_min_base_quality: u8) -> u8 {
    std::env::var("PARITY_ASM_MIN_TAIL_QUALITY")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or_else(|| gatk_min_tail_quality_for_assembly(default_min_base_quality))
}
use crate::assembly::{AssemblyGraphParams, AssemblyGraphPruningParams};
use crate::assembly_dangling_recovery::DanglingRecoveryParams;
use crate::cigar::{Cigar, CigarOperator};
use crate::haplotype::Haplotype;
use crate::kbest_haplotype::{
    find_best_haplotypes, find_best_haplotypes_for_assembly, find_best_haplotypes_preserving_cycles,
};
use crate::read_model::ReadFilterParams;
use crate::read_threading_assembler::{
    allow_low_complexity_expanded_kmer, allow_non_unique_expanded_kmer, ReadThreadingAssemblerArgs,
};
use crate::read_threading_graph::{
    assembly_graph_from_ref_and_reads_threading, reference_has_non_unique_kmers,
    threading_non_unique_summary,
};
use crate::seq_graph::SeqGraph;
use crate::seq_kbest_haplotype::find_best_haplotypes_seq_graph;
use crate::walker_apply::{
    call_disposition, select_region_for_asm_dump, AssemblyRegionCallDisposition,
};
use crate::walker_traversal::{
    flatten_assembly_regions, traverse_assembly_region_walker, WalkerTraversalConfig,
};
use gatk_common::{GatkError, GatkResult};
use gatk_core::reference::{parse_intervals_cli_string, ReferenceWindowCache, SequenceDictionary};
use std::collections::HashSet;
use std::io::Write;
use std::path::Path;

/// One row of staged assembly metrics for L2 parity ASM-1 dumps.
/// # Invariants
/// Count fields are post-stage graph/path summaries (nodes, edges, k-best paths, extracted haps).
/// `graph_kind` / `stage` are static labels for TSV comparison, not runtime enums.
/// # Ownership
/// Immutable snapshot row; `graph_kind` and `stage` borrow static str literals.
/// # Mutation
/// None after construction — diagnostic output only.
/// # Biological assumptions
/// None — parity instrumentation over assembly graphs, not a biological model.
/// # Java equivalence
/// Rust-native dump shape; metrics correspond to GATK assembly pipeline stages (threading → SeqGraph).
#[derive(Debug, Clone)]
pub struct AssemblyStageRow {
    pub graph_kind: &'static str,
    pub stage: &'static str,
    pub nodes: usize,
    pub edges: usize,
    pub ref_spine_vertices: usize,
    pub branch_vertices: usize,
    pub branch_vertices_all: usize,
    pub non_ref_edges_on_spine: usize,
    pub non_ref_edges_all: usize,
    pub kbest_paths: usize,
    pub extracted_haps: usize,
    pub non_ref_haps: usize,
    pub top_path_len: usize,
    pub top_path_eq_ref: bool,
}

fn forward_reach(graph: &AssemblyGraph, start: usize) -> HashSet<usize> {
    let mut keep = HashSet::new();
    let mut stack = vec![start];
    keep.insert(start);
    while let Some(v) = stack.pop() {
        for to in graph.outgoing_nodes(v) {
            if keep.insert(to) {
                stack.push(to);
            }
        }
    }
    keep
}

fn backward_reach(graph: &AssemblyGraph, start: usize) -> HashSet<usize> {
    let mut keep = HashSet::new();
    let mut stack = vec![start];
    keep.insert(start);
    while let Some(v) = stack.pop() {
        for from in graph.incoming_nodes(v) {
            if keep.insert(from) {
                stack.push(from);
            }
        }
    }
    keep
}

fn assembly_ref_spine_vertices(graph: &AssemblyGraph) -> HashSet<usize> {
    let Some(source) = graph.reference_source_vertex() else {
        return HashSet::new();
    };
    let Some(sink) = graph.reference_sink_vertex() else {
        return HashSet::new();
    };
    let from_source = forward_reach(graph, source);
    let from_sink = backward_reach(graph, sink);
    from_source
        .into_iter()
        .filter(|v| from_sink.contains(v))
        .collect()
}

fn assembly_spine_metrics(graph: &AssemblyGraph, spine: &HashSet<usize>) -> (usize, usize, usize) {
    let branch = spine
        .iter()
        .filter(|&&v| graph.outgoing_nodes(v).len() > 1)
        .count();
    let mut non_ref_edges = 0usize;
    for e in graph.edges_sorted() {
        if spine.contains(&e.from) && spine.contains(&e.to) && !graph.edge_is_ref(e.from, e.to) {
            non_ref_edges += 1;
        }
    }
    (spine.len(), branch, non_ref_edges)
}

fn assembly_whole_graph_metrics(graph: &AssemblyGraph) -> (usize, usize) {
    let branch = (0..graph.node_count())
        .filter(|&v| graph.outgoing_nodes(v).len() > 1)
        .count();
    let non_ref_edges = graph
        .edges_sorted()
        .iter()
        .filter(|e| !graph.edge_is_ref(e.from, e.to))
        .count();
    (branch, non_ref_edges)
}

fn record_rt_stage(
    rows: &mut Vec<AssemblyStageRow>,
    stage: &'static str,
    graph: &AssemblyGraph,
    ref_bytes: &[u8],
    args: &ReadThreadingAssemblerArgs,
    _ref_hap: &Haplotype,
    _ref_cigar_len: usize,
) -> GatkResult<()> {
    let spine = assembly_ref_spine_vertices(graph);
    let (spine_n, branch_on_spine, non_ref_on_spine) = assembly_spine_metrics(graph, &spine);
    let (branch_all, non_ref_all) = assembly_whole_graph_metrics(graph);
    let paths = find_best_haplotypes_for_assembly(graph, args.num_best_haplotypes_per_graph)?;
    let (top_len, top_eq) = paths.first().map_or((0, false), |p| {
        let b = graph.path_bases(p.start, &p.edges);
        (b.len(), b.as_slice() == ref_bytes)
    });
    // Java `HcFullParityGateDump.rtStageRow`: extractedHaps/nonRefHaps count KBest paths, not post-SW haps.
    rows.push(AssemblyStageRow {
        graph_kind: "rt",
        stage,
        nodes: graph.node_count(),
        edges: graph.edge_count(),
        ref_spine_vertices: spine_n,
        branch_vertices: branch_on_spine,
        branch_vertices_all: branch_all,
        non_ref_edges_on_spine: non_ref_on_spine,
        non_ref_edges_all: non_ref_all,
        kbest_paths: paths.len(),
        extracted_haps: paths.len(),
        non_ref_haps: paths.iter().filter(|p| !p.is_reference).count(),
        top_path_len: top_len,
        top_path_eq_ref: top_eq,
    });
    Ok(())
}

fn seq_forward_reach(graph: &SeqGraph, start: usize) -> HashSet<usize> {
    let mut keep = HashSet::new();
    let mut stack = vec![start];
    keep.insert(start);
    while let Some(v) = stack.pop() {
        for &t in graph.outgoing_of(v) {
            if keep.insert(t) {
                stack.push(t);
            }
        }
    }
    keep
}

fn seq_backward_reach(graph: &SeqGraph, start: usize) -> HashSet<usize> {
    let mut keep = HashSet::new();
    let mut stack = vec![start];
    keep.insert(start);
    while let Some(v) = stack.pop() {
        for p in graph.incoming_nodes(v) {
            if keep.insert(p) {
                stack.push(p);
            }
        }
    }
    keep
}

fn seq_ref_spine_vertices(graph: &SeqGraph) -> HashSet<usize> {
    let Some(source) = graph.reference_source_vertex() else {
        return HashSet::new();
    };
    let Some(sink) = graph.reference_sink_vertex() else {
        return HashSet::new();
    };
    let from_source = seq_forward_reach(graph, source);
    let from_sink = seq_backward_reach(graph, sink);
    from_source
        .into_iter()
        .filter(|v| from_sink.contains(v))
        .collect()
}

fn seq_spine_metrics(graph: &SeqGraph, spine: &HashSet<usize>) -> (usize, usize, usize) {
    let branch = spine
        .iter()
        .filter(|&&v| graph.outgoing_of(v).len() > 1)
        .count();
    let mut non_ref_edges = 0usize;
    for e in graph.edges_pub() {
        if spine.contains(&e.from) && spine.contains(&e.to) && !e.is_ref {
            non_ref_edges += 1;
        }
    }
    (spine.len(), branch, non_ref_edges)
}

fn seq_whole_graph_metrics(graph: &SeqGraph) -> (usize, usize) {
    let branch = (0..graph.node_count())
        .filter(|&v| graph.outgoing_of(v).len() > 1)
        .count();
    let non_ref_edges = graph.edges_pub().iter().filter(|e| !e.is_ref).count();
    (branch, non_ref_edges)
}

fn record_seq_stage(
    rows: &mut Vec<AssemblyStageRow>,
    stage: &'static str,
    graph: &SeqGraph,
    ref_bytes: &[u8],
    paths: &[crate::kbest_haplotype::KBestPath],
) {
    let spine = seq_ref_spine_vertices(graph);
    let (spine_n, branch, non_ref_e) = seq_spine_metrics(graph, &spine);
    let (branch_all, non_ref_all) = seq_whole_graph_metrics(graph);
    let (top_len, top_eq) = paths.first().map_or((0, false), |p| {
        let b = graph.path_bases_bytes(p.start, &p.edges);
        (b.len(), b == ref_bytes)
    });
    rows.push(AssemblyStageRow {
        graph_kind: "seq",
        stage,
        nodes: graph.node_count(),
        edges: graph.edge_count(),
        ref_spine_vertices: spine_n,
        branch_vertices: branch,
        branch_vertices_all: branch_all,
        non_ref_edges_on_spine: non_ref_e,
        non_ref_edges_all: non_ref_all,
        kbest_paths: paths.len(),
        extracted_haps: paths.len(),
        non_ref_haps: paths.iter().filter(|p| !p.is_reference).count(),
        top_path_len: top_len,
        top_path_eq_ref: top_eq,
    });
}

fn ref_hap_setup(reference: &crate::assembly::AssemblyRead) -> (Haplotype, usize) {
    let mut ref_hap = Haplotype::new(reference.bases.as_slice(), true);
    let mut ref_cigar = Cigar::new();
    ref_cigar.push(ref_hap.bases.len(), CigarOperator::Match);
    ref_hap.cigar = Some(ref_cigar);
    let ref_cigar_len = ref_hap.cigar.as_ref().unwrap().reference_length();
    (ref_hap, ref_cigar_len)
}

/// Graphs captured during k=85 ASM-1 probe for k-best path dumps.
/// # Invariants
/// Optional fields hold graph clones at named pipeline stages (RT pre/post remove-paths; SeqGraph before/after cleanup).
/// # Ownership
/// Owns optional [`AssemblyGraph`] / [`SeqGraph`] snapshots for dump-only use.
/// # Mutation
/// Filled during probe; read afterward for k-best path dumps.
/// # Biological assumptions
/// None — parity instrumentation.
/// # Java equivalence
/// Rust-native ASM-1 dump snapshots of GATK threading → SeqGraph stages.
#[derive(Default)]
pub struct K85GraphSnapshots {
    pub rt_after_dangling_pre_remove_paths: Option<AssemblyGraph>,
    pub rt_after_remove_paths: Option<AssemblyGraph>,
    pub seq_after_to_sequence_graph: Option<SeqGraph>,
    pub seq_after_cleanup: Option<SeqGraph>,
}

/// One ranked k-best path row for assembly parity dumps.
/// # Invariants
/// `rank` is 1-based best-first ordering within `(graph_kind, stage, strip_cycles)`.
/// `eq_ref` compares path sequence to the padded reference bytes used for the dump.
/// # Ownership
/// Owns path `sequence`; borrows static `graph_kind` / `stage` labels.
/// # Mutation
/// None — immutable dump record.
/// # Biological assumptions
/// None — characterizes graph search output for parity gates.
/// # Java equivalence
/// Rust-native dump; aligns with GATK k-best haplotype finder path listings in parity harnesses.
#[derive(Debug, Clone)]
pub struct KbestPathDumpRow {
    pub graph_kind: &'static str,
    pub stage: &'static str,
    pub strip_cycles: bool,
    pub rank: usize,
    pub score: f64,
    pub is_reference: bool,
    pub path_len: usize,
    pub eq_ref: bool,
    pub sequence: String,
}

fn write_rt_kbest_rows(
    rows: &mut Vec<KbestPathDumpRow>,
    stage: &'static str,
    strip_cycles: bool,
    graph: &AssemblyGraph,
    ref_bytes: &[u8],
    max_haps: usize,
) -> GatkResult<()> {
    let paths = if strip_cycles {
        find_best_haplotypes(graph, max_haps)?
    } else {
        find_best_haplotypes_preserving_cycles(graph, max_haps)?
    };
    for (rank, path) in paths.iter().enumerate() {
        let bases = graph.path_bases(path.start, &path.edges);
        rows.push(KbestPathDumpRow {
            graph_kind: "rt",
            stage,
            strip_cycles,
            rank,
            score: path.score,
            is_reference: path.is_reference,
            path_len: bases.len(),
            eq_ref: bases.as_slice() == ref_bytes,
            sequence: String::from_utf8_lossy(&bases).into_owned(),
        });
    }
    Ok(())
}

fn write_seq_kbest_rows(
    rows: &mut Vec<KbestPathDumpRow>,
    stage: &'static str,
    graph: &SeqGraph,
    ref_bytes: &[u8],
    max_haps: usize,
) {
    let paths = find_best_haplotypes_seq_graph(graph, max_haps).unwrap_or_default();
    for (rank, path) in paths.iter().enumerate() {
        let bases = graph.path_bases_bytes(path.start, &path.edges);
        rows.push(KbestPathDumpRow {
            graph_kind: "seq",
            stage,
            strip_cycles: false,
            rank,
            score: path.score,
            is_reference: path.is_reference,
            path_len: bases.len(),
            eq_ref: bases == ref_bytes,
            sequence: String::from_utf8_lossy(&bases).into_owned(),
        });
    }
}

fn kbest_rows_from_snapshots(
    snapshots: &K85GraphSnapshots,
    ref_bytes: &[u8],
    max_haps: usize,
) -> GatkResult<Vec<KbestPathDumpRow>> {
    let mut rows = Vec::new();
    if let Some(g) = &snapshots.rt_after_dangling_pre_remove_paths {
        write_rt_kbest_rows(
            &mut rows,
            "threading_after_dangling_pre_remove_paths",
            true,
            g,
            ref_bytes,
            max_haps,
        )?;
        write_rt_kbest_rows(
            &mut rows,
            "threading_after_dangling_pre_remove_paths",
            false,
            g,
            ref_bytes,
            max_haps,
        )?;
    }
    if let Some(g) = &snapshots.rt_after_remove_paths {
        write_rt_kbest_rows(
            &mut rows,
            "threading_after_remove_paths",
            true,
            g,
            ref_bytes,
            max_haps,
        )?;
    }
    if let Some(g) = &snapshots.seq_after_to_sequence_graph {
        write_seq_kbest_rows(
            &mut rows,
            "seq_after_to_sequence_graph",
            g,
            ref_bytes,
            max_haps,
        );
    }
    if let Some(g) = &snapshots.seq_after_cleanup {
        write_seq_kbest_rows(
            &mut rows,
            "seq_after_cleanup_seq_graph",
            g,
            ref_bytes,
            max_haps,
        );
    }
    Ok(rows)
}

fn probe_k85_stages(
    reference: &crate::assembly::AssemblyRead,
    reads: &[crate::assembly::AssemblyRead],
    args: &ReadThreadingAssemblerArgs,
) -> GatkResult<(
    Vec<AssemblyStageRow>,
    Option<crate::assembly_dangling_recovery::DanglingRecoverySummary>,
    Vec<(usize, String, String)>,
    bool,
    K85GraphSnapshots,
)> {
    const KMER: usize = 85;
    let allow_lc = allow_low_complexity_expanded_kmer(args, true);
    let allow_nu = allow_non_unique_expanded_kmer(args, true);
    let ref_bytes = reference.bases.as_slice();
    let mut rows = Vec::new();
    let mut dangling_meta = None;
    let mut dangling_tail_probes = Vec::new();
    let mut snapshots = K85GraphSnapshots::default();

    if reference.bases.len() < KMER {
        return Ok((rows, dangling_meta, dangling_tail_probes, false, snapshots));
    }
    if !allow_nu
        && !args.allow_non_unique_kmers_in_ref
        && reference_has_non_unique_kmers(reference, KMER)
    {
        return Ok((rows, dangling_meta, dangling_tail_probes, false, snapshots));
    }

    let params = AssemblyGraphParams {
        kmer_size: crate::bio_ids::KmerSize::try_new(KMER as u16).expect("KMER≥2"),
        min_base_quality: args.min_base_quality,
        min_edge_weight: 1,
        dangling_path_max_nodes: 0,
        max_haplotypes: args.num_best_haplotypes_per_graph,
        max_haplotype_bases: 4096,
        start_threading_only_at_existing_vertex: !args.recover_dangling_branches,
    };
    let mut graph = assembly_graph_from_ref_and_reads_threading(reference, reads, &params)?;
    let graph_has_cycles = graph.has_cycle();
    let summary = threading_non_unique_summary(Some(reference), reads, &params)?;
    if !allow_lc && !args.allow_low_complexity_graphs && summary.is_low_complexity {
        return Ok((
            rows,
            dangling_meta,
            dangling_tail_probes,
            graph_has_cycles,
            snapshots,
        ));
    }

    let (ref_hap, ref_cigar_len) = ref_hap_setup(reference);
    record_rt_stage(
        &mut rows,
        "threading_after_build",
        &graph,
        ref_bytes,
        args,
        &ref_hap,
        ref_cigar_len,
    )?;

    let mut pruning = AssemblyGraphPruningParams::gatk_haplotype_caller_defaults();
    pruning.min_prune_factor = args.min_prune_factor;
    pruning.use_adaptive_pruning = args.use_adaptive_pruning;
    if args.prune_before_cycle_counting {
        graph.apply_pruning(&pruning);
        record_rt_stage(
            &mut rows,
            "threading_after_prune_before_dangling",
            &graph,
            ref_bytes,
            args,
            &ref_hap,
            ref_cigar_len,
        )?;
    }

    if args.recover_dangling_branches && !graph.ref_nodes.is_empty() {
        let mut dangling = DanglingRecoveryParams::gatk_haplotype_caller_defaults();
        dangling.min_prune_factor = args.min_prune_factor;
        dangling.min_dangling_branch_length = args.min_dangling_branch_length;
        dangling.recover_dangling_heads = args.recover_dangling_heads;
        dangling.recover_all_dangling_branches = args.recover_all_dangling_branches;
        dangling.sw = args.dangling_end_sw;
        dangling_tail_probes = graph.probe_dangling_tail_failures(&dangling);
        let tails_attempted_parity = graph.count_dangling_tail_candidates_parity_dump();
        let heads_attempted_parity = if args.recover_dangling_heads {
            graph.count_dangling_head_candidates_parity_dump()
        } else {
            0
        };
        let mut summary = graph.recover_dangling_branches(&dangling)?;
        summary.tails_attempted = tails_attempted_parity;
        summary.heads_attempted = heads_attempted_parity;
        dangling_meta = Some(summary);
        record_rt_stage(
            &mut rows,
            "threading_after_dangling",
            &graph,
            ref_bytes,
            args,
            &ref_hap,
            ref_cigar_len,
        )?;
        snapshots.rt_after_dangling_pre_remove_paths = Some(graph.clone());
    }

    if !args.prune_before_cycle_counting {
        graph.apply_pruning(&pruning);
    }

    if args.remove_paths_not_connected_to_ref {
        if graph.reference_source_vertex().is_none() || graph.reference_sink_vertex().is_none() {
            return Ok((
                rows,
                dangling_meta,
                dangling_tail_probes,
                graph_has_cycles,
                snapshots,
            ));
        }
        graph.remove_paths_not_connected_to_ref()?;
    }

    if graph.reference_source_vertex().is_none() || graph.reference_sink_vertex().is_none() {
        return Ok((
            rows,
            dangling_meta,
            dangling_tail_probes,
            graph_has_cycles,
            snapshots,
        ));
    }
    snapshots.rt_after_remove_paths = Some(graph.clone());

    record_rt_stage(
        &mut rows,
        "threading_after_prune_dangling",
        &graph,
        ref_bytes,
        args,
        &ref_hap,
        ref_cigar_len,
    )?;

    let rt = graph;
    let mut seq = SeqGraph::from_assembly_graph(&rt);
    seq.clean_non_ref_paths();
    snapshots.seq_after_to_sequence_graph = Some(seq.clone());
    let mut seq_for_cleanup = seq.clone();
    let _ = seq_for_cleanup.cleanup_seq_graph();
    snapshots.seq_after_cleanup = Some(seq_for_cleanup);
    let paths0 = find_best_haplotypes_seq_graph(&seq, args.num_best_haplotypes_per_graph)
        .unwrap_or_default();
    record_seq_stage(
        &mut rows,
        "after_to_sequence_graph",
        &seq,
        ref_bytes,
        &paths0,
    );

    seq.zip_linear_chains();
    seq.remove_singleton_orphan_vertices();
    seq.remove_vertices_not_connected_to_ref_regardless_of_direction();
    let paths1 = find_best_haplotypes_seq_graph(&seq, args.num_best_haplotypes_per_graph)
        .unwrap_or_default();
    record_seq_stage(
        &mut rows,
        "after_zip_orphans_prune",
        &seq,
        ref_bytes,
        &paths1,
    );

    seq.simplify_graph();
    let paths2 = find_best_haplotypes_seq_graph(&seq, args.num_best_haplotypes_per_graph)
        .unwrap_or_default();
    record_seq_stage(&mut rows, "after_first_simplify", &seq, ref_bytes, &paths2);

    if seq.reference_source_vertex().is_none() || seq.reference_sink_vertex().is_none() {
        record_seq_stage(
            &mut rows,
            "would_abort_just_assembled_reference",
            &seq,
            ref_bytes,
            &[],
        );
        return Ok((
            rows,
            dangling_meta,
            dangling_tail_probes,
            graph_has_cycles,
            snapshots,
        ));
    }

    let _ = seq.remove_paths_not_connected_to_ref();
    seq.simplify_graph();
    let paths3 = find_best_haplotypes_seq_graph(&seq, args.num_best_haplotypes_per_graph)
        .unwrap_or_default();
    record_seq_stage(
        &mut rows,
        "after_remove_paths_final_simplify",
        &seq,
        ref_bytes,
        &paths3,
    );

    Ok((
        rows,
        dangling_meta,
        dangling_tail_probes,
        graph_has_cycles,
        snapshots,
    ))
}

/// First active region — k=85 k-best path rows at RT/Seq checkpoints (ASM-1b).
pub fn dump_assembly_region_kbest_paths_tsv(
    ref_fasta: &Path,
    bam_path: &Path,
    interval_cli: &str,
    padding: u64,
    out: &mut impl Write,
) -> GatkResult<()> {
    let dict = SequenceDictionary::from_fasta_path(ref_fasta)?;
    let specs = parse_intervals_cli_string(&dict, interval_cli)?;
    let filters = ReadFilterParams::gatk_standard_hc();
    let cfg = WalkerTraversalConfig::gatk_haplotype_caller_production(padding);
    let walk = traverse_assembly_region_walker(&dict, &specs, ref_fasta, bam_path, &filters, &cfg)?;
    let regions = flatten_assembly_regions(&walk);
    let mut ref_cache = ReferenceWindowCache::new(ref_fasta.to_path_buf(), 4);

    let region = select_region_for_asm_dump(&regions)
        .ok_or_else(|| GatkError::argument("no assembly region with reads in interval"))?;

    let reference = assembly_reference_read(&dict, &mut ref_cache, region)?;
    let assemble_args = AssembleReadsArgs::default();
    let reads = if std::env::var_os("PARITY_ASM_MATERIALIZE_READS").is_some() {
        assembly_reads_for_java_materialize_dump(&region.reads)
    } else {
        records_to_assembly_reads(&assembly_reads_for_production(
            &region.reads,
            region,
            parity_min_tail_quality(assemble_args.assembler.min_base_quality),
            assemble_args.correct_overlapping_base_qualities,
            false,
        ))
    };

    writeln!(out, "region_contig\t{}", region.contig)
        .map_err(|e| GatkError::generic(format!("write: {e}")))?;
    writeln!(out, "region_start\t{}", region.start.get())
        .map_err(|e| GatkError::generic(format!("write: {e}")))?;
    writeln!(out, "region_end\t{}", region.end.get())
        .map_err(|e| GatkError::generic(format!("write: {e}")))?;
    writeln!(out, "padded_ref_len\t{}", reference.bases.len())
        .map_err(|e| GatkError::generic(format!("write: {e}")))?;

    let (_rows, _dangling, _probes, graph_has_cycles, snapshots) =
        probe_k85_stages(&reference, &reads, &assemble_args.assembler)?;
    if graph_has_cycles {
        writeln!(out, "warn\tgraph_has_cycles\ttrue")
            .map_err(|e| GatkError::generic(format!("write: {e}")))?;
    }

    const MAX_HAPS: usize = 128;
    let path_rows = kbest_rows_from_snapshots(&snapshots, reference.bases.as_slice(), MAX_HAPS)?;
    writeln!(
        out,
        "graph\tstage\tstrip_cycles\trank\tscore\tis_reference\tpath_len\teq_ref\tsequence"
    )
    .map_err(|e| GatkError::generic(format!("write: {e}")))?;
    for r in &path_rows {
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            r.graph_kind,
            r.stage,
            r.strip_cycles,
            r.rank,
            format_kbest_score(r.score),
            r.is_reference,
            r.path_len,
            r.eq_ref,
            r.sequence,
        )
        .map_err(|e| GatkError::generic(format!("write: {e}")))?;
    }
    Ok(())
}

fn format_kbest_score(score: f64) -> String {
    if score == 0.0 {
        "0".to_string()
    } else {
        format!("{score:.8}")
    }
}

/// First active region on interval — k=85 assembly stage dump (ASM-1).
pub fn dump_assembly_region_assembly_stages_tsv(
    ref_fasta: &Path,
    bam_path: &Path,
    interval_cli: &str,
    padding: u64,
    out: &mut impl Write,
) -> GatkResult<()> {
    let dict = SequenceDictionary::from_fasta_path(ref_fasta)?;
    let specs = parse_intervals_cli_string(&dict, interval_cli)?;
    let filters = ReadFilterParams::gatk_standard_hc();
    let cfg = WalkerTraversalConfig::gatk_haplotype_caller_production(padding);
    let walk = traverse_assembly_region_walker(&dict, &specs, ref_fasta, bam_path, &filters, &cfg)?;
    let regions = flatten_assembly_regions(&walk);
    let mut ref_cache = ReferenceWindowCache::new(ref_fasta.to_path_buf(), 4);

    let region = select_region_for_asm_dump(&regions)
        .ok_or_else(|| GatkError::argument("no assembly region with reads in interval"))?;

    let reference = assembly_reference_read(&dict, &mut ref_cache, region)?;
    let assemble_args = AssembleReadsArgs::default();
    let reads = if std::env::var_os("PARITY_ASM_MATERIALIZE_READS").is_some() {
        assembly_reads_for_java_materialize_dump(&region.reads)
    } else {
        records_to_assembly_reads(&assembly_reads_for_production(
            &region.reads,
            region,
            parity_min_tail_quality(assemble_args.assembler.min_base_quality),
            assemble_args.correct_overlapping_base_qualities,
            false,
        ))
    };

    writeln!(out, "region_contig\t{}", region.contig)
        .map_err(|e| GatkError::generic(format!("write: {e}")))?;
    writeln!(out, "region_start\t{}", region.start.get())
        .map_err(|e| GatkError::generic(format!("write: {e}")))?;
    writeln!(out, "region_end\t{}", region.end.get())
        .map_err(|e| GatkError::generic(format!("write: {e}")))?;
    writeln!(out, "padded_ref_len\t{}", reference.bases.len())
        .map_err(|e| GatkError::generic(format!("write: {e}")))?;
    writeln!(out, "read_count\t{}", reads.len())
        .map_err(|e| GatkError::generic(format!("write: {e}")))?;
    if !matches!(
        call_disposition(region),
        AssemblyRegionCallDisposition::ActiveFull
    ) {
        writeln!(out, "warn\tregion_inactive_using_reads_for_asm_diagnostic")
            .map_err(|e| GatkError::generic(format!("write: {e}")))?;
    }

    let (rows, dangling_meta, dangling_tail_probes, graph_has_cycles, _snapshots) =
        probe_k85_stages(&reference, &reads, &assemble_args.assembler)?;
    if graph_has_cycles {
        writeln!(out, "warn\tgraph_has_cycles\ttrue")
            .map_err(|e| GatkError::generic(format!("write: {e}")))?;
    }
    if let Some(d) = dangling_meta {
        // Match `HcFullParityGateDump.assemblyRegionAssemblyStages` (counts only; recovery outcome unknown).
        writeln!(
            out,
            "dangling_recovery\ttails_attempted={}\ttails_recovered=unknown\theads_attempted={}\theads_recovered=unknown\tedges_before={}\tedges_after={}",
            d.tails_attempted,
            d.heads_attempted,
            d.edges_before,
            d.edges_after,
        )
        .map_err(|e| GatkError::generic(format!("write: {e}")))?;
    }
    if std::env::var_os("PARITY_ASM_DANGLING_PROBES").is_some() {
        for (_v, kmer, reason) in &dangling_tail_probes {
            writeln!(out, "dangling_tail_probe\t{kmer}\t{reason}")
                .map_err(|e| GatkError::generic(format!("write: {e}")))?;
        }
    }
    writeln!(
        out,
        "graph\tstage\tnodes\tedges\tref_spine_vertices\tbranch_vertices\tbranch_vertices_all\tnon_ref_edges_on_spine\tnon_ref_edges_all\tkbest_paths\textracted_haps\tnon_ref_haps\ttop_path_len\ttop_path_eq_ref"
    )
    .map_err(|e| GatkError::generic(format!("write: {e}")))?;
    for r in &rows {
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            r.graph_kind,
            r.stage,
            r.nodes,
            r.edges,
            r.ref_spine_vertices,
            r.branch_vertices,
            r.branch_vertices_all,
            r.non_ref_edges_on_spine,
            r.non_ref_edges_all,
            r.kbest_paths,
            r.extracted_haps,
            r.non_ref_haps,
            r.top_path_len,
            r.top_path_eq_ref,
        )
        .map_err(|e| GatkError::generic(format!("write: {e}")))?;
    }
    Ok(())
}
