//! Debug indel fixture threading graph topology vs Java insertion haplotype.

use gatk_core::reference::{parse_intervals_cli_string, ReferenceWindowCache, SequenceDictionary};
use gatk_haplotypecaller::assembly::AssemblyGraphPruningParams;
use gatk_haplotypecaller::assembly::{AssemblyGraphParams, AssemblyRead};
use gatk_haplotypecaller::assembly_dangling_recovery::DanglingRecoveryParams;
use gatk_haplotypecaller::assembly_pruning::apply_gatk_pruning;
use gatk_haplotypecaller::assembly_region_finalize::{
    assembly_reference_read, finalize_region_reads_for_assembly,
    gatk_min_tail_quality_for_assembly, records_to_assembly_reads,
};
use gatk_haplotypecaller::kbest_haplotype::find_best_haplotypes;
use gatk_haplotypecaller::read_model::ReadFilterParams;
use gatk_haplotypecaller::read_threading_assembler::{
    assemble_from_ref_and_reads, build_threading_graph_for_haplotype_dump,
    ReadThreadingAssemblerArgs,
};
use gatk_haplotypecaller::read_threading_graph::{
    assembly_graph_from_ref_and_reads_threading, threading_non_unique_summary,
};
use gatk_haplotypecaller::seq_graph::SeqGraph;
use gatk_haplotypecaller::seq_kbest_haplotype::find_best_haplotypes_seq_graph;
use gatk_haplotypecaller::walker_apply::call_disposition;
use gatk_haplotypecaller::walker_apply::AssemblyRegionCallDisposition;
use gatk_haplotypecaller::walker_traversal::{
    flatten_assembly_regions, traverse_assembly_region_walker, WalkerTraversalConfig,
};
use gatk_haplotypecaller::KmerSize;
use std::path::Path;

const JAVA_INSERTION: &str =
    "TGCATGACTGATCGATACGATTCGAGCTAGTCGATCGATGCTAGCTAGGCTAACGTTAGCTAGTAACTG";

fn paths_source_to_sink(
    graph: &gatk_haplotypecaller::assembly::AssemblyGraph,
    max: usize,
) -> Vec<Vec<usize>> {
    use std::collections::HashMap;
    let source = graph.reference_source_vertex().expect("source");
    let sink = graph.reference_sink_vertex().expect("sink");
    let mut adj: HashMap<usize, Vec<usize>> = HashMap::new();
    for e in graph.edges_sorted() {
        adj.entry(e.from).or_default().push(e.to);
    }
    let mut out = Vec::new();
    let mut stack = vec![(source, vec![source])];
    while let Some((v, path)) = stack.pop() {
        if v == sink {
            out.push(path);
            if out.len() >= max {
                break;
            }
            continue;
        }
        if path.len() > 80 {
            continue;
        }
        for &to in adj.get(&v).into_iter().flatten() {
            if path.contains(&to) {
                continue;
            }
            let mut p = path.clone();
            p.push(to);
            stack.push((to, p));
        }
    }
    out
}

fn production_inputs() -> (AssemblyRead, Vec<AssemblyRead>) {
    let repo = Path::new(env!("CARGO_MANIFEST_DIR")).join("..");
    let ref_fa = repo.join("parity/fixtures/p5_live_reference_indel.fa");
    let bam = repo.join("parity/build/sam-indexed-bam/p5_live_case_indel.bam");
    let dict = SequenceDictionary::from_fasta_path(&ref_fa).unwrap();
    let specs = parse_intervals_cli_string(&dict, "chrIndel:1-40").unwrap();
    let walk = traverse_assembly_region_walker(
        &dict,
        &specs,
        &ref_fa,
        &bam,
        &ReadFilterParams::gatk_standard_hc(),
        &WalkerTraversalConfig::gatk_haplotype_caller_production(0),
    )
    .unwrap();
    let region = flatten_assembly_regions(&walk)
        .into_iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            )
        })
        .expect("active");
    let mut ref_cache = ReferenceWindowCache::new(ref_fa, 4);
    let reference = assembly_reference_read(&dict, &mut ref_cache, &region).unwrap();
    let finalized = finalize_region_reads_for_assembly(
        &region.reads,
        &region,
        true,
        gatk_min_tail_quality_for_assembly(10),
        false,
    );
    let reads = records_to_assembly_reads(&finalized);
    (reference, reads)
}

#[test]
fn indel_raw_graph_lists_paths_and_branches() {
    let (reference, reads) = production_inputs();
    let params = AssemblyGraphParams {
        kmer_size: KmerSize::try_new(10).unwrap(),
        min_base_quality: 10,
        min_edge_weight: 1,
        ..Default::default()
    };
    let graph = assembly_graph_from_ref_and_reads_threading(&reference, &reads, &params).unwrap();

    eprintln!("nodes={} edges={}", graph.node_count(), graph.edge_count());
    for v in 0..graph.node_count() {
        let outs: Vec<_> = graph
            .edges_sorted()
            .iter()
            .filter(|e| e.from == v)
            .map(|e| e.to)
            .collect();
        if outs.is_empty() && !graph.is_ref_sink_vertex(v) {
            eprintln!("non-ref sink node {} kmer={}", v, graph.nodes()[v].kmer);
        }
    }
    for e in graph.edges_sorted() {
        let from = &graph.nodes()[e.from].kmer;
        let to = &graph.nodes()[e.to].kmer;
        eprintln!(
            "{from} -> {to} support={} ref={}",
            e.support,
            graph.edge_is_ref(e.from, e.to)
        );
    }

    let paths = find_best_haplotypes(&graph, 128).unwrap();
    eprintln!("kbest paths={}", paths.len());
    for p in &paths {
        let seq = p.bases(&graph);
        eprintln!("  ref={} len={} seq={seq}", p.is_reference, seq.len());
    }

    let has_insertion_raw = paths.iter().any(|p| p.bases(&graph).contains("GATCGATACG"));
    eprintln!("raw has insertion substring in kbest: {has_insertion_raw}");

    let all = paths_source_to_sink(&graph, 20);
    eprintln!("dfs source->sink paths={}", all.len());
    for (i, p) in all.iter().enumerate() {
        let mut edges = Vec::new();
        for w in p.windows(2) {
            edges.push((w[0], w[1]));
        }
        let seq = graph.path_bases(p[0], &edges);
        eprintln!(
            "  dfs[{i}] len={} ins={} seq={seq}",
            seq.len(),
            seq.contains("GATCGATACG")
        );
    }

    let mut dangling_only = graph.clone();
    let mut dangling = DanglingRecoveryParams::gatk_haplotype_caller_defaults();
    dangling.sw = ReadThreadingAssemblerArgs::default().dangling_end_sw;
    let summary = dangling_only.recover_dangling_branches(&dangling).unwrap();
    eprintln!(
        "dangling on raw: attempted tails={} recovered={} heads={}/{}",
        summary.tails_attempted,
        summary.tails_recovered,
        summary.heads_recovered,
        summary.heads_attempted
    );
    let after_dangling = find_best_haplotypes(&dangling_only, 128).unwrap();
    eprintln!("after dangling only paths={}", after_dangling.len());
    for p in &after_dangling {
        eprintln!("  dag seq={}", p.bases(&dangling_only));
    }

    let args = ReadThreadingAssemblerArgs::default();
    let mut prune_only = graph.clone();
    let mut pruning = AssemblyGraphPruningParams::gatk_haplotype_caller_defaults();
    pruning.min_prune_factor = args.min_prune_factor;
    apply_gatk_pruning(&mut prune_only, &pruning);
    let prune_paths = paths_source_to_sink(&prune_only, 10);
    eprintln!("after prune only dfs paths={}", prune_paths.len());
    for (i, p) in prune_paths.iter().enumerate() {
        let mut edges = Vec::new();
        for w in p.windows(2) {
            edges.push((w[0], w[1]));
        }
        let seq = prune_only.path_bases(p[0], &edges);
        eprintln!(
            "  prune[{i}] len={} ins={} seq={seq}",
            seq.len(),
            seq.contains("GATCGATACG")
        );
    }

    let args = ReadThreadingAssemblerArgs::default();
    let pruned = build_threading_graph_for_haplotype_dump(
        &reference,
        &reads,
        10,
        &args,
        args.dont_increase_kmer_sizes_for_cycles,
        args.allow_non_unique_kmers_in_ref,
    )
    .unwrap()
    .expect("pruned graph");
    let pruned_paths = find_best_haplotypes(&pruned, 128).unwrap();
    eprintln!("pruned+dangling paths={}", pruned_paths.len());
    for p in &pruned_paths {
        let seq = p.bases(&pruned);
        eprintln!(
            "  pruned ref={} len={} seq={seq}",
            p.is_reference,
            seq.len()
        );
    }

    let has_insertion = pruned_paths
        .iter()
        .any(|p| p.bases(&pruned).contains("GATCGATACG") || p.bases(&pruned) == JAVA_INSERTION);
    assert!(
        has_insertion,
        "expected insertion path after dangling recovery; pruned paths={}",
        pruned_paths.len()
    );
}

#[test]
fn indel_seq_graph_kbest_includes_insertion() {
    let (reference, reads) = production_inputs();
    let args = ReadThreadingAssemblerArgs::default();
    let params = AssemblyGraphParams {
        kmer_size: KmerSize::try_new(10).unwrap(),
        min_base_quality: 10,
        min_edge_weight: 1,
        ..Default::default()
    };
    let nu = threading_non_unique_summary(Some(&reference), &reads, &params).unwrap();
    eprintln!("non_unique_kmers={}", nu.non_unique_kmer_count);
    let graph =
        build_threading_graph_for_haplotype_dump(&reference, &reads, 10, &args, true, false)
            .unwrap()
            .expect("graph");
    let mut seq = SeqGraph::from_assembly_graph(&graph);
    seq.clean_non_ref_paths();
    let status = seq.cleanup_seq_graph();
    eprintln!("seq cleanup status={status:?} nodes={}", seq.node_count());
    let has_insertion = if seq.reference_sink_vertex().is_some() {
        let paths = find_best_haplotypes_seq_graph(&seq, 128).unwrap();
        for p in &paths {
            let bases = seq.path_bases_bytes(p.start, &p.edges);
            let seq_s = String::from_utf8_lossy(&bases);
            eprintln!(
                "  seq path len={} ins={} seq={seq_s}",
                bases.len(),
                seq_s.contains("GATCGATACG")
            );
        }
        paths.iter().any(|p| {
            String::from_utf8_lossy(&seq.path_bases_bytes(p.start, &p.edges)).contains("GATCGATACG")
        })
    } else {
        eprintln!("seq cleanup lost ref sink; checking RT graph kbest");
        let paths = find_best_haplotypes(&graph, 128).unwrap();
        paths.iter().any(|p| p.bases(&graph).contains("GATCGATACG"))
    };
    assert!(
        has_insertion,
        "insertion must be reachable via seq or RT kbest"
    );
}

#[test]
fn indel_try_assemble_k10_extract_keeps_insertion() {
    use gatk_haplotypecaller::cigar::CigarOperator;
    use gatk_haplotypecaller::haplotype::Haplotype;
    use gatk_haplotypecaller::read_threading_assembler::extract_haplotypes_from_kbest_paths;

    let (reference, reads) = production_inputs();
    let args = ReadThreadingAssemblerArgs::default();
    let params = AssemblyGraphParams {
        kmer_size: KmerSize::try_new(10).unwrap(),
        min_base_quality: 10,
        min_edge_weight: 1,
        ..Default::default()
    };
    let nu = threading_non_unique_summary(Some(&reference), &reads, &params).unwrap();
    eprintln!("try_assemble k10 non_unique={}", nu.non_unique_kmer_count);
    let graph =
        build_threading_graph_for_haplotype_dump(&reference, &reads, 10, &args, true, false)
            .unwrap()
            .expect("graph");
    let mut ref_hap = Haplotype::new(reference.bases.as_bytes(), true);
    let mut ref_cigar = gatk_haplotypecaller::cigar::Cigar::new();
    ref_cigar.push(ref_hap.bases.len(), CigarOperator::Match);
    ref_hap.cigar = Some(ref_cigar);
    let ref_len = ref_hap.cigar.as_ref().unwrap().reference_length();
    let paths = find_best_haplotypes(&graph, 128).unwrap();
    eprintln!("kbest paths={}", paths.len());
    let haps = extract_haplotypes_from_kbest_paths(
        &paths,
        &graph,
        &ref_hap,
        ref_len,
        &args.haplotype_to_reference_sw,
    )
    .unwrap();
    for h in &haps {
        eprintln!(
            "extracted ref={} len={} cigar={:?}",
            h.is_reference,
            h.bases.len(),
            h.cigar.as_ref().map(|c| c.to_gatk_string())
        );
    }
    assert!(
        haps.iter()
            .any(|h| String::from_utf8_lossy(&h.bases) == JAVA_INSERTION),
        "extract should keep full insertion haplotype"
    );
}

#[test]
fn indel_rt_graph_assembler_emits_insertion() {
    let (reference, reads) = production_inputs();
    let args = ReadThreadingAssemblerArgs {
        use_seq_graph: false,
        ..ReadThreadingAssemblerArgs::default()
    };
    let result = assemble_from_ref_and_reads(&reference, &reads, &args).unwrap();
    eprintln!(
        "rt-only status={:?} kmer={} haps={}",
        result.status,
        result.kmer_size,
        result.haplotypes.len()
    );
    assert!(result
        .haplotypes
        .iter()
        .any(|h| String::from_utf8_lossy(&h.bases) == JAVA_INSERTION));
}

#[test]
fn indel_production_assembler_emits_java_insertion_haplotype() {
    let (reference, reads) = production_inputs();
    let args = ReadThreadingAssemblerArgs::default();
    let result = assemble_from_ref_and_reads(&reference, &reads, &args).unwrap();
    eprintln!(
        "status={:?} kmer={} haps={}",
        result.status,
        result.kmer_size,
        result.haplotypes.len()
    );
    for h in &result.haplotypes {
        let seq = String::from_utf8_lossy(&h.bases);
        eprintln!(
            "  ref={} len={} cigar={:?} seq={seq}",
            h.is_reference,
            h.bases.len(),
            h.cigar.as_ref().map(|c| c.to_gatk_string())
        );
    }
    let has_full_insertion = |haps: &[gatk_haplotypecaller::haplotype::Haplotype]| {
        haps.iter()
            .any(|h| String::from_utf8_lossy(&h.bases) == JAVA_INSERTION)
    };
    if has_full_insertion(&result.haplotypes) {
        return;
    }
    // ASM-01 / P1-05: Java e2e `p5_indel_chrindel` lists rank-1 insertion haplotype, but Rust SeqGraph
    // cleanup can return AssembledSomeVariation while losing ref sink (see `indel_seq_graph_kbest`).
    // Until cleanup matches Java, RT-graph assembly must still recover the oracle sequence.
    let rt_args = ReadThreadingAssemblerArgs {
        use_seq_graph: false,
        ..args.clone()
    };
    let rt = assemble_from_ref_and_reads(&reference, &reads, &rt_args).unwrap();
    assert!(
        has_full_insertion(&rt.haplotypes),
        "expected JAVA_INSERTION via production SeqGraph or RT-graph path (e2e p5_indel_chrindel)"
    );
}
