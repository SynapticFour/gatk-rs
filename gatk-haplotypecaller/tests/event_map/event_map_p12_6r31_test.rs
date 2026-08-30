//! 6R.31 TEST-ONLY: isolate which RT cleanup stage drops TAGAGTTGAAG (361T+371G).
//! Does not change production dangling recovery, prune, k, or haplotype suppression.

#[cfg(test)]
mod traces {
    use crate::assembly::{AssemblyGraph, AssemblyGraphParams, AssemblyGraphPruningParams};
    use crate::assembly_based_caller::AssembleReadsArgs;
    use crate::assembly_dangling_recovery::DanglingRecoveryParams;
    use crate::assembly_region_finalize::{
        assembly_reference_read, create_graph_reference_read, finalize_region_reads_for_assembly,
        gatk_min_tail_quality_for_assembly, records_to_assembly_reads,
    };
    use crate::bio_ids::KmerSize;
    use crate::read_threading_graph::assembly_graph_from_ref_and_reads_threading_with_summary;
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
    /// Java 0.1 unique dangling head (in-degree 0) that walks into TAGAGTTGAAG.
    const JAVA_ALT_HEAD: &[u8] = b"CAAATAAAAGGTAGACAGCAGCATT";

    #[derive(Debug, Clone)]
    struct RtCleanupSnapshot {
        stage: &'static str,
        node_count: usize,
        edge_count: usize,
        target_present: bool,
        target_vertices: Vec<usize>,
        from_src_count: usize,
        to_sink_count: usize,
        both_count: usize,
        dangling_heads: Vec<usize>,
        dangling_tails: Vec<usize>,
        alt_head_present: bool,
        alt_head_from_src: bool,
    }

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

    fn ascii(b: &[u8]) -> String {
        String::from_utf8_lossy(b).into_owned()
    }

    fn vertices_with_motif(graph: &AssemblyGraph, needle: &[u8]) -> Vec<usize> {
        (0..graph.node_count())
            .filter(|&v| graph.kmer_at(v).windows(needle.len()).any(|w| w == needle))
            .collect()
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

    fn bfs_backward(graph: &AssemblyGraph, start: usize) -> HashSet<usize> {
        let mut seen = HashSet::new();
        let mut stack = vec![start];
        seen.insert(start);
        while let Some(v) = stack.pop() {
            for from in graph.incoming_nodes(v) {
                if seen.insert(from) {
                    stack.push(from);
                }
            }
        }
        seen
    }

    fn from_src_set(graph: &AssemblyGraph) -> HashSet<usize> {
        match graph.reference_source_vertex() {
            Some(s) => bfs_forward(graph, s),
            None => HashSet::new(),
        }
    }

    fn to_sink_set(graph: &AssemblyGraph) -> HashSet<usize> {
        match graph.reference_sink_vertex() {
            Some(s) => bfs_backward(graph, s),
            None => HashSet::new(),
        }
    }

    fn dangling_heads(graph: &AssemblyGraph) -> Vec<usize> {
        (0..graph.node_count())
            .filter(|&v| {
                graph.incoming_count(v) == 0
                    && !graph.outgoing_nodes(v).is_empty()
                    && !graph.is_ref_source_vertex(v)
            })
            .collect()
    }

    fn dangling_tails(graph: &AssemblyGraph) -> Vec<usize> {
        (0..graph.node_count())
            .filter(|&v| graph.outgoing_nodes(v).is_empty() && !graph.is_ref_sink_vertex(v))
            .collect()
    }

    fn role_of(
        graph: &AssemblyGraph,
        v: usize,
        from_src: &HashSet<usize>,
        to_sink: &HashSet<usize>,
    ) -> String {
        let inn = graph.incoming_count(v);
        let out = graph.outgoing_nodes(v).len();
        let mut parts = Vec::new();
        if inn == 0 && out > 0 && !graph.is_ref_source_vertex(v) {
            parts.push("dangling_head");
        }
        if out == 0 && !graph.is_ref_sink_vertex(v) {
            parts.push("dangling_tail");
        }
        if graph.is_ref_source_vertex(v) {
            parts.push("ref_source");
        }
        if graph.is_ref_sink_vertex(v) {
            parts.push("ref_sink");
        }
        if from_src.contains(&v) && to_sink.contains(&v) {
            parts.push("ref_connected");
        } else if from_src.contains(&v) {
            parts.push("from_src_only");
        } else if to_sink.contains(&v) {
            parts.push("to_sink_only");
        } else {
            parts.push("disconnected");
        }
        if parts.is_empty() {
            parts.push("internal");
        }
        parts.join(",")
    }

    fn snapshot(stage: &'static str, graph: &AssemblyGraph) -> RtCleanupSnapshot {
        let targets = vertices_with_motif(graph, TARGET);
        let from_src = from_src_set(graph);
        let to_sink = to_sink_set(graph);
        let both: HashSet<usize> = from_src.intersection(&to_sink).copied().collect();
        let heads = dangling_heads(graph);
        let tails = dangling_tails(graph);
        let alt_head = vertices_with_motif(graph, JAVA_ALT_HEAD);
        let alt_head_from_src = alt_head.iter().any(|v| from_src.contains(v));
        RtCleanupSnapshot {
            stage,
            node_count: graph.node_count(),
            edge_count: graph.edge_count(),
            target_present: !targets.is_empty(),
            target_vertices: targets,
            from_src_count: from_src.len(),
            to_sink_count: to_sink.len(),
            both_count: both.len(),
            dangling_heads: heads,
            dangling_tails: tails,
            alt_head_present: !alt_head.is_empty(),
            alt_head_from_src,
        }
    }

    fn dump_snapshot(graph: &AssemblyGraph, snap: &RtCleanupSnapshot) {
        eprintln!(
            "SNAP[{}] nodes={} edges={} target={} n_target_verts={} from_src={} to_sink={} both={} heads={} tails={} java_alt_head={} alt_head_from_src={}",
            snap.stage,
            snap.node_count,
            snap.edge_count,
            snap.target_present,
            snap.target_vertices.len(),
            snap.from_src_count,
            snap.to_sink_count,
            snap.both_count,
            snap.dangling_heads.len(),
            snap.dangling_tails.len(),
            snap.alt_head_present,
            snap.alt_head_from_src
        );
        let from_src = from_src_set(graph);
        let to_sink = to_sink_set(graph);
        for &v in &snap.target_vertices {
            let inn = graph.incoming_nodes(v);
            let out = graph.outgoing_nodes(v);
            let w_in: u32 = inn.iter().filter_map(|&p| graph.edge_support(p, v)).sum();
            let w_out: u32 = out.iter().filter_map(|&t| graph.edge_support(v, t)).sum();
            eprintln!(
                "  TARGET_V id={v} kmer={} in={} out={} w_in={w_in} w_out={w_out} role={}",
                ascii(graph.kmer_at(v)),
                inn.len(),
                out.len(),
                role_of(graph, v, &from_src, &to_sink)
            );
            for &p in &inn {
                eprintln!(
                    "    <- {} {} support={} ref={}",
                    p,
                    ascii(graph.kmer_at(p)),
                    graph.edge_support(p, v).unwrap_or(0),
                    graph.edge_is_ref(p, v)
                );
            }
            for &t in &out {
                eprintln!(
                    "    -> {} {} support={} ref={}",
                    t,
                    ascii(graph.kmer_at(t)),
                    graph.edge_support(v, t).unwrap_or(0),
                    graph.edge_is_ref(v, t)
                );
            }
        }
        for &h in &snap.dangling_heads {
            eprintln!(
                "  HEAD id={h} kmer={} out={} role={}",
                ascii(graph.kmer_at(h)),
                graph.outgoing_nodes(h).len(),
                role_of(graph, h, &from_src, &to_sink)
            );
        }
        for &t in &snap.dangling_tails {
            eprintln!(
                "  TAIL id={t} kmer={} in={} role={}",
                ascii(graph.kmer_at(t)),
                graph.incoming_count(t),
                role_of(graph, t, &from_src, &to_sink)
            );
        }
    }

    fn dump_head_probes(graph: &AssemblyGraph, params: &DanglingRecoveryParams) {
        let probes = graph.probe_dangling_head_failures(params);
        eprintln!("HEAD_PROBES n={}", probes.len());
        for (v, kmer, reason) in &probes {
            let on_target = kmer.as_bytes().windows(TARGET.len()).any(|w| w == TARGET)
                || kmer.as_bytes() == JAVA_ALT_HEAD;
            if on_target || probes.len() <= 8 {
                eprintln!(
                    "  HEAD_PROBE v={v} kmer={kmer} reason={reason} target_related={on_target}"
                );
            }
        }
        let tails = graph.probe_dangling_tail_failures(params);
        eprintln!("TAIL_PROBES n={}", tails.len());
        for (v, kmer, reason) in &tails {
            let on_target = kmer.as_bytes().windows(TARGET.len()).any(|w| w == TARGET);
            if on_target || tails.len() <= 8 {
                eprintln!(
                    "  TAIL_PROBE v={v} kmer={kmer} reason={reason} target_related={on_target}"
                );
            }
        }
    }

    #[test]
    fn six_r31_java_cleanup_order_is_tails_heads_remove_paths() {
        let assembler = include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/src/read_threading_assembler.rs"
        ));
        assert!(assembler.contains("recover_dangling_branches"));
        assert!(assembler.contains("remove_paths_not_connected_to_ref"));
        let dangling = include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/src/assembly_dangling_recovery.rs"
        ));
        assert!(
            dangling.contains("recover_dangling_tail"),
            "tails exist as a distinct step"
        );
        assert!(dangling.contains("recover_dangling_head"));
        let java_exact_block = dangling
            .split("if params.dangling_java_exact")
            .nth(1)
            .unwrap_or("");
        let tails_at = java_exact_block.find("recover_dangling_tail").unwrap_or(0);
        let heads_at = java_exact_block.find("recover_dangling_head").unwrap_or(0);
        assert!(
            tails_at < heads_at,
            "Java-exact path must recover tails before heads"
        );
        assert!(!dangling.contains("92317361"));
        assert!(!dangling.contains("TAGAGTTGAAG"));
    }

    #[test]
    fn six_r31_mid_b_rt_cleanup_target_tracking() {
        let Some((ref_fasta, bam_path)) = fixture_paths() else {
            eprintln!("Real-data mid-B comparison unavailable");
            return;
        };
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
        assert_eq!(rt_args.min_prune_factor, 2);
        assert_eq!(rt_args.min_dangling_branch_length, 4);
        assert!(rt_args.recover_dangling_heads);
        assert!(rt_args.remove_paths_not_connected_to_ref);

        let finalized = finalize_region_reads_for_assembly(
            &region.reads,
            &region,
            assemble_args.correct_overlapping_base_qualities,
            gatk_min_tail_quality_for_assembly(rt_args.min_base_quality),
            false,
        );
        let assembly_reads = records_to_assembly_reads(&finalized);
        let params = graph_params(K);
        let (mut graph, _) = assembly_graph_from_ref_and_reads_threading_with_summary(
            &graph_ref,
            &assembly_reads,
            &params,
        )
        .expect("raw rt");

        let raw = snapshot("raw", &graph);
        dump_snapshot(&graph, &raw);
        assert_eq!(raw.node_count, 518, "6R.30 pin: Java/Rust raw = 518");
        assert!(raw.target_present, "raw TAGAGTTGAAG must be present");

        let mut pruning = AssemblyGraphPruningParams::gatk_haplotype_caller_defaults();
        pruning.min_prune_factor = rt_args.min_prune_factor;
        graph.apply_pruning(&pruning);
        let post_prune = snapshot("post_prune", &graph);
        dump_snapshot(&graph, &post_prune);
        assert_eq!(
            post_prune.node_count, 516,
            "6R.30 pin: after chain prune = 516"
        );
        assert!(
            post_prune.target_present,
            "post_prune TAGAGTTGAAG must be present"
        );

        let dangling = DanglingRecoveryParams::from_assembler_args(&rt_args);
        assert!(dangling.dangling_java_exact);
        dump_head_probes(&graph, &dangling);

        let mut post_tails = None;
        let mut post_heads = None;
        let summary = graph
            .test_java_exact_dangling_tails_then_heads(
                &dangling,
                |g, attempted, recovered| {
                    let s = snapshot("post_dangling_tails", g);
                    eprintln!("TAILS attempted={attempted} recovered={recovered}");
                    dump_snapshot(g, &s);
                    post_tails = Some(s);
                },
                |g, attempted, recovered| {
                    let s = snapshot("post_dangling_heads", g);
                    eprintln!("HEADS attempted={attempted} recovered={recovered}");
                    dump_snapshot(g, &s);
                    post_heads = Some(s);
                },
            )
            .expect("traced dangling");
        eprintln!(
            "DANGLING_SUMMARY tails={}/{} heads={}/{} edges {}->{}",
            summary.tails_recovered,
            summary.tails_attempted,
            summary.heads_recovered,
            summary.heads_attempted,
            summary.edges_before,
            summary.edges_after
        );
        let post_tails = post_tails.expect("tails snap");
        let post_heads = post_heads.expect("heads snap");
        let post_dangling_cleanup = snapshot("post_dangling_cleanup", &graph);
        dump_snapshot(&graph, &post_dangling_cleanup);

        graph
            .remove_paths_not_connected_to_ref()
            .expect("removePaths");
        let post_remove = snapshot("post_remove_paths", &graph);
        dump_snapshot(&graph, &post_remove);

        eprintln!(
            "TRANSITION raw={} prune={} tails={} heads={} dangling_cleanup={} remove_paths={}",
            raw.target_present,
            post_prune.target_present,
            post_tails.target_present,
            post_heads.target_present,
            post_dangling_cleanup.target_present,
            post_remove.target_present
        );
        eprintln!(
            "NODES raw={} prune={} tails={} heads={} dangling_cleanup={} remove_paths={} java_0_2=444",
            raw.node_count,
            post_prune.node_count,
            post_tails.node_count,
            post_heads.node_count,
            post_dangling_cleanup.node_count,
            post_remove.node_count
        );

        // 6R.33: head recovery no longer attaches this branch; removePaths drops TAGAGTTGAAG.
        assert!(
            !post_remove.target_present,
            "post_remove_paths TAGAGTTGAAG must be absent (Java 0.2); 6R.33 mismatch-cap reject"
        );
        assert!(
            post_tails.target_present,
            "post_dangling_tails TAGAGTTGAAG present; tails remain a no-op"
        );
        assert!(
            post_heads.target_present,
            "post_dangling_heads does not delete vertices; 11-mer still present until removePaths"
        );
        assert!(
            post_dangling_cleanup.target_present,
            "isolated-node cleanup after dangling must not drop TAGAGTTGAAG"
        );
        assert_eq!(post_prune.from_src_count, 444);
        assert_eq!(post_tails.from_src_count, 444);
        assert_eq!(post_tails.node_count, 516);
        assert_eq!(summary.tails_attempted, 0);
        assert_eq!(summary.heads_attempted, 1);
        assert_eq!(summary.heads_recovered, 0);
        assert_eq!(post_heads.from_src_count, 444);
        assert_eq!(post_heads.node_count, 516);
        assert_eq!(post_remove.from_src_count, 444);
        assert!(
            post_prune.alt_head_present && !post_prune.alt_head_from_src,
            "Java dangling-head k-mer is present after prune but not ref-source-reachable"
        );
        assert!(
            !post_remove.alt_head_present,
            "unattached dangling-head chain is dropped by removePaths"
        );
    }
}
