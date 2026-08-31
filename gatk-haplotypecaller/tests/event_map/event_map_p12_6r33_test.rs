//! 6R.33: Java 4.4 `bestPrefixMatchLegacy` mismatch-cap abort on canonical mid-B.
//! Does not change `path_bases` / `getBasesForPath` encoding (separate parity item).

#[cfg(test)]
mod traces {
    use crate::assembly::{
        AssemblyGraph, AssemblyGraphParams, AssemblyGraphPruningParams, AssemblyRead,
    };
    use crate::assembly_based_caller::AssembleReadsArgs;
    use crate::assembly_dangling_recovery::DanglingRecoveryParams;
    use crate::assembly_region_finalize::{
        assembly_reference_read, create_graph_reference_read, finalize_region_reads_for_assembly,
        gatk_min_tail_quality_for_assembly, records_to_assembly_reads,
    };
    use crate::bio_ids::KmerSize;
    use crate::read_threading_assembler::build_threading_graph_for_seq_assembly;
    use crate::read_threading_graph::assembly_graph_from_ref_and_reads_threading_with_summary;
    use crate::seq_graph::SeqGraph;
    use crate::walker_apply::{call_disposition, AssemblyRegionCallDisposition};
    use gatk_core::reference::{
        parse_intervals_cli_string, ReferenceWindowCache, SequenceDictionary,
    };
    use std::collections::HashSet;
    use std::path::Path;

    const SITE_CA: u64 = 92_317_399;
    const JAVA_ACTIVE: (u64, u64) = (92_317_262, 92_317_491);
    const JAVA_EXTENDED: (u64, u64) = (92_317_162, 92_317_591);
    const K: usize = 25;
    const TARGET: &[u8] = b"TAGAGTTGAAG";
    const JAVA_ALT_HEAD: &[u8] = b"CAAATAAAAGGTAGACAGCAGCATT";
    const PATH_CAP: usize = 32;

    fn fixture_paths() -> Option<(std::path::PathBuf, std::path::PathBuf)> {
        let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("..");
        let ref_path = root.join("parity/realworld/assets/hs37d5.simple.fa");
        let bam = root.join("parity/realworld/na12878_20k_b37/NA12878_20k.b37.bam");
        if !ref_path.is_file() || !bam.is_file() {
            return None;
        }
        Some((ref_path, bam))
    }

    fn graph_params(k: usize) -> AssemblyGraphParams {
        AssemblyGraphParams {
            kmer_size: KmerSize::try_from_usize(k).expect("k"),
            min_base_quality: 10,
            min_edge_weight: 1,
            dangling_path_max_nodes: 0,
            max_haplotypes: 128,
            max_haplotype_bases: 4096,
            start_threading_only_at_existing_vertex: false,
        }
    }

    fn has_motif(graph: &AssemblyGraph, needle: &[u8]) -> bool {
        (0..graph.node_count()).any(|v| graph.kmer_at(v).windows(needle.len()).any(|w| w == needle))
    }

    fn bfs_forward(graph: &AssemblyGraph, start: usize) -> HashSet<usize> {
        let mut seen = HashSet::new();
        let mut stack = vec![start];
        seen.insert(start);
        while let Some(v) = stack.pop() {
            for to in graph.outgoing_nodes(v) {
                if seen.insert(to) {
                    stack.push(to);
                }
            }
        }
        seen
    }

    fn from_src_count(graph: &AssemblyGraph) -> usize {
        graph
            .reference_source_vertex()
            .map(|s| bfs_forward(graph, s).len())
            .unwrap_or(0)
    }

    fn count_st_paths(graph: &SeqGraph) -> usize {
        let Some(src) = graph.reference_source_vertex() else {
            return 0;
        };
        let Some(sink) = graph.reference_sink_vertex() else {
            return 0;
        };
        fn walk(g: &SeqGraph, v: usize, sink: usize, on_path: &mut [bool], found: &mut usize) {
            if *found >= PATH_CAP {
                return;
            }
            if v == sink {
                *found += 1;
                return;
            }
            if v >= on_path.len() || on_path[v] {
                return;
            }
            on_path[v] = true;
            for to in g.outgoing_nodes(v) {
                walk(g, to, sink, on_path, found);
            }
            on_path[v] = false;
        }
        let mut on_path = vec![false; graph.node_count()];
        let mut found = 0usize;
        walk(graph, src, sink, &mut on_path, &mut found);
        found
    }

    fn post_prune_mid_b() -> Option<(
        AssemblyGraph,
        DanglingRecoveryParams,
        crate::read_threading_assembler::ReadThreadingAssemblerArgs,
        AssemblyRead,
        Vec<AssemblyRead>,
    )> {
        let (ref_fasta, bam_path) = fixture_paths()?;
        let dict = SequenceDictionary::from_fasta_path(&ref_fasta).expect("dict");
        let walk_iv = parse_intervals_cli_string(&dict, "2:92317000-92319000").expect("iv");
        let filters = crate::read_model::ReadFilterParams::gatk_standard_hc();
        let cfg =
            crate::walker_traversal::WalkerTraversalConfig::gatk_haplotype_caller_production(100);
        let walk = crate::walker_traversal::traverse_assembly_region_walker(
            &dict, &walk_iv, &ref_fasta, &bam_path, &filters, &cfg,
        )
        .expect("walk");
        let region = crate::walker_traversal::flatten_assembly_regions(&walk)
            .into_iter()
            .find(|r| {
                matches!(
                    call_disposition(r),
                    AssemblyRegionCallDisposition::ActiveFull
                ) && r.start.get() <= SITE_CA
                    && r.end.get() >= SITE_CA
            })
            .expect("ActiveFull mid-B");
        assert_eq!((region.start.get(), region.end.get()), JAVA_ACTIVE);
        assert_eq!(
            (region.extended_start.get(), region.extended_end.get()),
            JAVA_EXTENDED
        );
        let mut ref_cache = ReferenceWindowCache::new(ref_fasta.clone(), 4);
        let padded = assembly_reference_read(&dict, &mut ref_cache, &region).expect("pad");
        let graph_ref = create_graph_reference_read(&padded, &region, &dict);
        let mut assemble_args = AssembleReadsArgs::default();
        assemble_args.strict_java_assembly = true;
        let mut rt_args = assemble_args.assembler.clone();
        rt_args.dangling_java_exact = true;
        let finalized = finalize_region_reads_for_assembly(
            &region.reads,
            &region,
            assemble_args.correct_overlapping_base_qualities,
            gatk_min_tail_quality_for_assembly(rt_args.min_base_quality),
            false,
        );
        let assembly_reads = records_to_assembly_reads(&finalized);
        let (mut graph, _) = assembly_graph_from_ref_and_reads_threading_with_summary(
            &graph_ref,
            &assembly_reads,
            &graph_params(K),
        )
        .expect("raw rt");
        assert_eq!(graph.node_count(), 518);
        let mut pruning = AssemblyGraphPruningParams::gatk_haplotype_caller_defaults();
        pruning.min_prune_factor = rt_args.min_prune_factor;
        graph.apply_pruning(&pruning);
        assert_eq!(graph.node_count(), 516);
        let dangling = DanglingRecoveryParams::from_assembler_args(&rt_args);
        Some((graph, dangling, rt_args, graph_ref, assembly_reads))
    }

    #[test]
    fn six_r33_mid_b_head_rejected_target_removed_no_double_diamond() {
        let Some((mut graph, dangling, rt_args, graph_ref, assembly_reads)) = post_prune_mid_b()
        else {
            eprintln!("Real-data mid-B comparison unavailable");
            return;
        };
        let heads: Vec<usize> = (0..graph.node_count())
            .filter(|&v| {
                graph.incoming_count(v) == 0
                    && !graph.outgoing_nodes(v).is_empty()
                    && !graph.is_ref_source_vertex(v)
            })
            .collect();
        assert_eq!(heads.len(), 1);
        let head = heads[0];
        assert_eq!(graph.kmer_at(head), JAVA_ALT_HEAD);
        assert!(has_motif(&graph, TARGET));
        assert_eq!(from_src_count(&graph), 444);

        let dump = graph.test_dangling_head_decision_dump(head, &dangling);
        eprintln!(
            "6R.33 HEAD dump rust_plan={} rust_idx={:?} java_idx={} first_m={} mm={} cap={}",
            dump.rust_plan,
            dump.rust_idx,
            dump.java_idx_on_rust_seq,
            dump.first_m_len,
            dump.mismatches_in_first_m,
            dump.max_mismatches_legacy
        );
        assert_eq!(dump.first_m_len, 97);
        assert_eq!(dump.max_mismatches_legacy, 3);
        assert_eq!(dump.mismatches_in_first_m, 10);
        assert_eq!(dump.rust_idx, None);
        assert_ne!(dump.rust_idx, Some(35));
        assert_eq!(dump.final_rust, "REJECT");
        assert_eq!(dump.rust_plan, "prefix_match_legacy_failed");
        assert_eq!(dump.java_idx_on_rust_seq, -1);
        assert_eq!(
            dump.rust_alt_bases, dump.java_alt_bases,
            "6R.36: path_bases matches getBasesForPath"
        );

        let from_src_before = from_src_count(&graph);
        let mut post_heads_from_src = 0usize;
        let summary = graph
            .test_java_exact_dangling_tails_then_heads(
                &dangling,
                |_, _, _| {},
                |g, attempted, recovered| {
                    eprintln!("HEADS attempted={attempted} recovered={recovered}");
                    post_heads_from_src = from_src_count(g);
                },
            )
            .expect("dangling");
        assert_eq!(summary.tails_attempted, 0);
        assert_eq!(summary.heads_attempted, 1);
        assert_eq!(summary.heads_recovered, 0);
        assert_eq!(post_heads_from_src, from_src_before);
        assert_eq!(post_heads_from_src, 444);
        assert!(
            has_motif(&graph, TARGET),
            "recovery does not delete the 11-mer; removePaths does"
        );

        graph
            .remove_paths_not_connected_to_ref()
            .expect("removePaths");
        assert!(
            !has_motif(&graph, TARGET),
            "TAGAGTTGAAG absent after removePaths (Java 0.2)"
        );
        assert_eq!(from_src_count(&graph), 444);
        assert!(!has_motif(&graph, JAVA_ALT_HEAD));

        let cleaned = build_threading_graph_for_seq_assembly(
            &graph_ref,
            &assembly_reads,
            K,
            &rt_args,
            false,
            false,
        )
        .expect("rt")
        .expect("k=25 graph");
        assert!(
            !has_motif(&cleaned, TARGET),
            "cleaned RT for SeqGraph must not carry TAGAGTTGAAG"
        );
        let mut seq = SeqGraph::from_assembly_graph(&cleaned);
        seq.clean_non_ref_paths();
        let _ = seq.cleanup_seq_graph();
        let n_paths = count_st_paths(&seq);
        eprintln!(
            "6R.33 SEQ nodes={} edges={} st_paths={} (Java zip is a single diamond, 2 paths)",
            seq.node_count(),
            seq.edge_count(),
            n_paths
        );
        assert!(
            n_paths <= 2,
            "must not retain the 361/371 extra diamond; st_paths={n_paths}"
        );
        assert_eq!(
            n_paths, 2,
            "expected Java-like single diamond (2 s-t paths), got {n_paths}"
        );
    }
}
