//! 6R.34: Java `getBasesForPath` vs Rust `path_bases` encoding (fixture E + downstream).
//! Production `path_bases` is intentionally unchanged (Outcome A).

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
    /// OBSERVED RUST production `path_bases` on the canonical 6R.32 alt path.
    const RUST_ALT: &[u8] =
        b"AGTTTGGACGAGATACTTTCCCTTAGAAGTTGAGATACTCAACTTACGTCTGTAGTCTTTCTTTAAAGACTCT";
    /// SOURCE-DERIVED JAVA extra bases: reverse(head k-mer) without the already-emitted suffix.
    const JAVA_ALT_EXTRA: &[u8] = b"TACGACGACAGATGGAAAATAAAC";

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

    fn mismatch_positions(ref_b: &[u8], alt_b: &[u8], window: usize) -> Vec<usize> {
        let n = window.min(ref_b.len()).min(alt_b.len());
        (0..n).filter(|&i| ref_b[i] != alt_b[i]).collect()
    }

    #[test]
    fn six_r34_canonical_path_encoding_and_downstream_same_reject() {
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

        let dump = graph.test_dangling_head_decision_dump(head, &dangling);

        let head_rev: Vec<u8> = JAVA_ALT_HEAD.iter().copied().rev().collect();
        let java_alt_expected: Vec<u8> = RUST_ALT
            .iter()
            .copied()
            .chain(head_rev[1..].iter().copied())
            .collect();
        assert_eq!(dump.rust_alt_bases, java_alt_expected);
        assert_eq!(dump.rust_alt_bases.len(), 97);
        assert_eq!(&dump.java_alt_bases[RUST_ALT.len()..], JAVA_ALT_EXTRA);
        assert_eq!(dump.java_alt_bases, java_alt_expected);
        assert_eq!(dump.java_alt_bases.len(), 97);
        assert_eq!(dump.rust_ref_bases.len(), 235);
        assert_eq!(dump.java_ref_bases.len(), 235);
        assert_eq!(dump.rust_alt_bases, dump.java_alt_bases);
        assert_eq!(dump.rust_ref_bases, dump.java_ref_bases);
        let last_ref = dump.ref_path_vertices.last().expect("ref path");
        assert!(
            last_ref.is_source,
            "ref-path terminal is the expanded source"
        );
        assert!(
            !dump.alt_path_vertices[0].is_source,
            "canonical alt path[0] is the LCA, not a source"
        );
        assert!(
            dump.alt_path_vertices.last().expect("head").is_source,
            "canonical alt path last vertex is the dangling-head source"
        );

        let rust_mm =
            mismatch_positions(&dump.rust_ref_bases, &dump.rust_alt_bases, dump.first_m_len);
        let java_mm = mismatch_positions(
            &dump.java_ref_bases,
            &dump.java_alt_bases,
            dump.java_seq_first_m_len,
        );
        eprintln!(
            "6R.34 CANONICAL rust_alt={} java_alt={} rust_ref={} java_ref={} rust_cigar={} java_cigar={} rust_first_m={} java_first_m={} rust_mm={:?} java_mm={:?} rust_cap={} java_cap={} rust_idx={:?} java_idx_rust={} java_idx_java={} cigar_ok={} java_cigar_ok={} rust={} java={}",
            dump.rust_alt_bases.len(),
            dump.java_alt_bases.len(),
            dump.rust_ref_bases.len(),
            dump.java_ref_bases.len(),
            dump.rust_cigar,
            dump.java_seq_cigar,
            dump.first_m_len,
            dump.java_seq_first_m_len,
            rust_mm,
            java_mm,
            dump.max_mismatches_legacy,
            (dump.java_seq_first_m_len / K).max(1),
            dump.rust_idx,
            dump.java_idx_on_rust_seq,
            dump.java_idx_on_java_seq,
            dump.cigar_ok,
            dump.java_seq_cigar_ok,
            dump.final_rust,
            dump.final_java_source_derived
        );

        assert_eq!(dump.rust_cigar, "97M");
        assert_eq!(dump.java_seq_cigar, "97M");
        assert_eq!(dump.cigar_ok, "PASS");
        assert!(dump.java_seq_cigar_ok);
        assert_eq!(dump.first_m_len, 97);
        assert_eq!(dump.java_seq_first_m_len, 97);
        assert_eq!(dump.mismatches_in_first_m, 10);
        assert_eq!(dump.java_seq_mismatches_in_first_m, 10);
        assert_eq!(rust_mm, vec![25, 35, 39, 49, 63, 70, 85, 87, 92, 94]);
        assert_eq!(java_mm, vec![25, 35, 39, 49, 63, 70, 85, 87, 92, 94]);
        assert_eq!(dump.max_mismatches_legacy, 3);
        assert_eq!((dump.java_seq_first_m_len / K).max(1), 3);
        assert_eq!(dump.rust_idx, None);
        assert_ne!(dump.rust_idx, Some(35));
        assert_eq!(dump.java_idx_on_rust_seq, -1);
        assert_eq!(dump.java_idx_on_java_seq, -1);
        assert_eq!(dump.final_rust, "REJECT");
        assert_eq!(dump.final_java_source_derived, "REJECT");
        assert_eq!(dump.rust_plan, "prefix_match_legacy_failed");
        assert!(
            dump.merge_from_kmer.is_empty() && dump.merge_to_kmer.is_empty(),
            "6R.33 merge CTTTCTGATGTTTGCATTCAAGTCA → TTTCTGATGTCTGCATTCAACTCAT must not return"
        );

        let from_src_before = from_src_count(&graph);
        let mut post_heads_from_src = 0usize;
        let summary = graph
            .test_java_exact_dangling_tails_then_heads(
                &dangling,
                |_, _, _| {},
                |g, attempted, recovered| {
                    eprintln!("6R.34 HEADS attempted={attempted} recovered={recovered}");
                    post_heads_from_src = from_src_count(g);
                },
            )
            .expect("dangling");
        assert_eq!(summary.tails_attempted, 0);
        assert_eq!(summary.heads_attempted, 1);
        assert_eq!(summary.heads_recovered, 0);
        assert_eq!(post_heads_from_src, from_src_before);
        assert_eq!(post_heads_from_src, 444);

        graph
            .remove_paths_not_connected_to_ref()
            .expect("removePaths");
        assert!(
            !has_motif(&graph, TARGET),
            "TAGAGTTGAAG absent after removePaths"
        );
        assert_eq!(from_src_count(&graph), 444);

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
        assert!(!has_motif(&cleaned, TARGET));
        let mut seq = SeqGraph::from_assembly_graph(&cleaned);
        seq.clean_non_ref_paths();
        let _ = seq.cleanup_seq_graph();
        let n_paths = count_st_paths(&seq);
        assert_eq!(n_paths, 2, "Java-like single diamond, got {n_paths}");
    }
}
