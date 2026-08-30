//! 6R.37 TEST-ONLY: first remaining Java/Rust cleaned-graph divergence after 6R.33 + 6R.36.
//! Diagnostic only. Does not change production assembly, dangling recovery, SeqGraph, or EventMap.

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
    use crate::read_threading_assembler::{
        build_threading_graph_for_seq_assembly, DEFAULT_NUM_BEST_HAPLOTYPES_PER_GRAPH,
    };
    use crate::read_threading_graph::assembly_graph_from_ref_and_reads_threading_with_summary;
    use crate::seq_graph::SeqGraph;
    use crate::seq_kbest_haplotype::find_best_haplotypes_seq_graph;
    use crate::walker_apply::{call_disposition, AssemblyRegionCallDisposition};
    use gatk_core::reference::{
        parse_intervals_cli_string, ReferenceWindowCache, SequenceDictionary,
    };
    use std::collections::HashSet;
    use std::path::Path;

    const SITE_CT: u64 = 92_317_361;
    const SITE_CG: u64 = 92_317_371;
    const SITE_CA: u64 = 92_317_399;
    const SITE_TC: u64 = 92_317_407;
    const SITE_GC: u64 = 92_317_412;
    const JAVA_ACTIVE: (u64, u64) = (92_317_262, 92_317_491);
    const JAVA_EXTENDED: (u64, u64) = (92_317_162, 92_317_591);
    const K: usize = 25;
    const PATH_CAP: usize = 32;
    const TARGET: &[u8] = b"TAGAGTTGAAG";
    const JAVA_ALT_HEAD: &[u8] = b"CAAATAAAAGGTAGACAGCAGCATT";
    const JAVA_REF_SOURCE: &[u8] = b"GAAACTTTTTCATGATGTATCTACT";
    const JAVA_REF_SINK: &[u8] = b"AAATCCTCTTTTTGTAAAATCTGCA";
    const RUST_ALT_SUFFIX: &[u8] =
        b"AGTTTGGACGAGATACTTTCCCTTAGAAGTTGAGATACTCAACTTACGTCTGTAGTCTTTCTTTAAAGACTCT";
    const JAVA_ALT_EXTRA: &[u8] = b"TACGACGACAGATGGAAAATAAAC";
    const JAVA_0_2_KMER_TXT: &str = include_str!("data/java_6r37_0_2_kmers.txt");
    const JAVA_0_1_MINUS_0_2_TXT: &str = include_str!("data/java_6r37_0_1_minus_0_2_kmers.txt");
    const JAVA_1_1_SEQ_TXT: &str = include_str!("data/java_6r37_1_1_seqs.txt");
    const JAVA_1_4_SEQ_TXT: &str = include_str!("data/java_6r37_1_4_seqs.txt");
    /// FNV-1a 64 of newline-joined sorted k-mers (no trailing newline). OBSERVED JAVA DOT.
    const JAVA_0_1_FNV: u64 = 0xa270_fb53_c0ad_ae7d;
    const JAVA_0_2_FNV: u64 = 0xc32e_7209_0b08_1d90;
    const JAVA_0_1_MINUS_0_2_FNV: u64 = 0x4679_5f74_20e3_fc18;
    const JAVA_1_1_FNV: u64 = 0x4599_49a4_2a59_be40;
    const JAVA_1_4_FNV: u64 = 0x9996_d504_b605_8b1c;

    #[derive(Debug, Clone)]
    struct StageSnap {
        stage: &'static str,
        node_count: usize,
        edge_count: usize,
        n_sources: usize,
        n_sinks: usize,
        from_src: usize,
        to_sink: usize,
        both: usize,
        neither: usize,
        target_present: bool,
        target_vertices: Vec<usize>,
        alt_head_present: bool,
        n_361_371_combined: usize,
        n_361_single_t: usize,
        n_371_single_g: usize,
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

    fn parse_lines(s: &str) -> Vec<String> {
        s.lines()
            .map(str::trim)
            .filter(|l| !l.is_empty())
            .map(str::to_string)
            .collect()
    }

    fn fnv64(data: &[u8]) -> u64 {
        let mut h = 0xcbf2_9ce4_8422_2325u64;
        for &b in data {
            h ^= u64::from(b);
            h = h.wrapping_mul(0x0100_0000_01b3);
        }
        h
    }

    fn kmer_blob(kmers: &[String]) -> Vec<u8> {
        kmers.join("\n").into_bytes()
    }

    fn sorted_graph_kmers(graph: &AssemblyGraph) -> Vec<String> {
        let mut v: Vec<String> = (0..graph.node_count())
            .map(|i| String::from_utf8_lossy(graph.kmer_at(i)).into_owned())
            .collect();
        v.sort();
        v
    }

    fn sorted_seq_payloads(graph: &SeqGraph) -> Vec<String> {
        let mut v: Vec<String> = (0..graph.node_count())
            .map(|i| String::from_utf8_lossy(graph.vertex_sequence(i)).into_owned())
            .collect();
        v.sort();
        v
    }

    fn ascii(b: &[u8]) -> String {
        String::from_utf8_lossy(b).into_owned()
    }

    fn vertices_with_motif(graph: &AssemblyGraph, needle: &[u8]) -> Vec<usize> {
        (0..graph.node_count())
            .filter(|&v| graph.kmer_at(v).windows(needle.len()).any(|w| w == needle))
            .collect()
    }

    fn exact_kmer_ids(graph: &AssemblyGraph, kmers: &[Vec<u8>]) -> Vec<usize> {
        kmers
            .iter()
            .filter_map(|k| graph.vertex_id_for_kmer(k))
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
        graph
            .reference_source_vertex()
            .map(|s| bfs_forward(graph, s))
            .unwrap_or_default()
    }

    fn to_sink_set(graph: &AssemblyGraph) -> HashSet<usize> {
        graph
            .reference_sink_vertex()
            .map(|s| bfs_backward(graph, s))
            .unwrap_or_default()
    }

    fn n_sources(graph: &AssemblyGraph) -> usize {
        (0..graph.node_count())
            .filter(|&v| graph.incoming_count(v) == 0)
            .count()
    }

    fn n_sinks(graph: &AssemblyGraph) -> usize {
        (0..graph.node_count())
            .filter(|&v| graph.outgoing_nodes(v).is_empty())
            .count()
    }

    fn overlapping_kmer_starts(seq_len: usize, off: usize, k: usize) -> Vec<usize> {
        if seq_len < k || off >= seq_len {
            return Vec::new();
        }
        let first = off.saturating_sub(k - 1);
        let last = off.min(seq_len - k);
        if first > last {
            Vec::new()
        } else {
            (first..=last).collect()
        }
    }

    fn windows_overlapping_all(seq_len: usize, offs: &[usize], k: usize) -> Vec<usize> {
        if seq_len < k || offs.is_empty() || offs.iter().any(|&o| o >= seq_len) {
            return Vec::new();
        }
        let first = offs.iter().map(|&o| o.saturating_sub(k - 1)).max().unwrap();
        let last = offs.iter().copied().min().unwrap().min(seq_len - k);
        if first > last {
            Vec::new()
        } else {
            (first..=last).collect()
        }
    }

    fn snp_alt_kmers(ref_bases: &[u8], off: usize, k: usize, alt: u8) -> Vec<Vec<u8>> {
        overlapping_kmer_starts(ref_bases.len(), off, k)
            .into_iter()
            .map(|start| {
                let mut rk = ref_bases[start..start + k].to_vec();
                rk[off - start] = alt;
                rk
            })
            .collect()
    }

    fn combined_alt_kmers(ref_bases: &[u8], k: usize, muts: &[(usize, u8)]) -> Vec<Vec<u8>> {
        let offs: Vec<usize> = muts.iter().map(|(o, _)| *o).collect();
        windows_overlapping_all(ref_bases.len(), &offs, k)
            .into_iter()
            .map(|start| {
                let mut ak = ref_bases[start..start + k].to_vec();
                for &(off, alt) in muts {
                    ak[off - start] = alt;
                }
                ak
            })
            .collect()
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

    fn snapshot(
        stage: &'static str,
        graph: &AssemblyGraph,
        combined_361_371: &[Vec<u8>],
        single_361_t: &[Vec<u8>],
        single_371_g: &[Vec<u8>],
    ) -> StageSnap {
        let from_src = from_src_set(graph);
        let to_sink = to_sink_set(graph);
        let both: HashSet<usize> = from_src.intersection(&to_sink).copied().collect();
        let neither = (0..graph.node_count())
            .filter(|v| !from_src.contains(v) && !to_sink.contains(v))
            .count();
        let targets = vertices_with_motif(graph, TARGET);
        StageSnap {
            stage,
            node_count: graph.node_count(),
            edge_count: graph.edge_count(),
            n_sources: n_sources(graph),
            n_sinks: n_sinks(graph),
            from_src: from_src.len(),
            to_sink: to_sink.len(),
            both: both.len(),
            neither,
            target_present: !targets.is_empty(),
            target_vertices: targets,
            alt_head_present: graph.vertex_id_for_kmer(JAVA_ALT_HEAD).is_some(),
            n_361_371_combined: exact_kmer_ids(graph, combined_361_371).len(),
            n_361_single_t: exact_kmer_ids(graph, single_361_t).len(),
            n_371_single_g: exact_kmer_ids(graph, single_371_g).len(),
        }
    }

    fn dump_stage(graph: &AssemblyGraph, s: &StageSnap) {
        eprintln!(
            "SNAP[{}] nodes={} edges={} src={} sink={} from_src={} to_sink={} both={} neither={} target={} n_target={} head={} 361T+371G={} 361T_only={} 371G_only={}",
            s.stage,
            s.node_count,
            s.edge_count,
            s.n_sources,
            s.n_sinks,
            s.from_src,
            s.to_sink,
            s.both,
            s.neither,
            s.target_present,
            s.target_vertices.len(),
            s.alt_head_present,
            s.n_361_371_combined,
            s.n_361_single_t,
            s.n_371_single_g
        );
        let from_src = from_src_set(graph);
        let to_sink = to_sink_set(graph);
        for &v in &s.target_vertices {
            let inn = graph.incoming_nodes(v);
            let out = graph.outgoing_nodes(v);
            let on_st = from_src.contains(&v) && to_sink.contains(&v);
            eprintln!(
                "  TARGET_V id={v} kmer={} in={} out={} from_src={} to_sink={} s_t={}",
                ascii(graph.kmer_at(v)),
                inn.len(),
                out.len(),
                from_src.contains(&v),
                to_sink.contains(&v),
                on_st
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
        if let Some(h) = graph.vertex_id_for_kmer(JAVA_ALT_HEAD) {
            eprintln!(
                "  HEAD id={h} in={} out={} from_src={} to_sink={}",
                graph.incoming_count(h),
                graph.outgoing_nodes(h).len(),
                from_src.contains(&h),
                to_sink.contains(&h)
            );
        }
    }

    fn set_diff<'a>(a: &'a HashSet<String>, b: &'a HashSet<String>) -> Vec<&'a String> {
        let mut v: Vec<&String> = a.difference(b).collect();
        v.sort();
        v
    }

    #[test]
    fn six_r37_java_dot_fixtures_are_well_formed() {
        let cleaned = parse_lines(JAVA_0_2_KMER_TXT);
        let lost = parse_lines(JAVA_0_1_MINUS_0_2_TXT);
        let zipped = parse_lines(JAVA_1_1_SEQ_TXT);
        let final_sg = parse_lines(JAVA_1_4_SEQ_TXT);
        assert_eq!(cleaned.len(), 444);
        assert_eq!(lost.len(), 72);
        assert_eq!(zipped.len(), 4);
        assert_eq!(final_sg.len(), 4);
        assert_eq!(fnv64(&kmer_blob(&cleaned)), JAVA_0_2_FNV);
        assert_eq!(fnv64(&kmer_blob(&lost)), JAVA_0_1_MINUS_0_2_FNV);
        assert_eq!(fnv64(&kmer_blob(&zipped)), JAVA_1_1_FNV);
        assert_eq!(fnv64(&kmer_blob(&final_sg)), JAVA_1_4_FNV);
        let prune: HashSet<String> = cleaned.iter().chain(lost.iter()).cloned().collect();
        assert_eq!(prune.len(), 516);
        let mut prune_sorted: Vec<String> = prune.into_iter().collect();
        prune_sorted.sort();
        assert_eq!(fnv64(&kmer_blob(&prune_sorted)), JAVA_0_1_FNV);
        assert!(lost.iter().any(|k| k.as_bytes() == JAVA_ALT_HEAD));
        assert_eq!(
            lost.iter()
                .filter(|k| k.as_bytes().windows(TARGET.len()).any(|w| w == TARGET))
                .count(),
            15
        );
        assert!(!cleaned.iter().any(|k| k.as_bytes() == JAVA_ALT_HEAD));
        assert!(!cleaned
            .iter()
            .any(|k| k.as_bytes().windows(TARGET.len()).any(|w| w == TARGET)));
        assert!(cleaned.iter().any(|k| k.as_bytes() == JAVA_REF_SOURCE));
        assert!(cleaned.iter().any(|k| k.as_bytes() == JAVA_REF_SINK));
    }

    #[test]
    fn six_r37_production_sources_have_no_locus_pins() {
        let dangling = include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/src/assembly_dangling_recovery.rs"
        ));
        let assembler = include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/src/read_threading_assembler.rs"
        ));
        let seq = include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/src/seq_graph.rs"));
        for src in [dangling, assembler, seq] {
            assert!(!src.contains("TAGAGTTGAAG"));
            assert!(!src.contains("92317361"));
            assert!(!src.contains("92317371"));
        }
        assert_eq!(DEFAULT_NUM_BEST_HAPLOTYPES_PER_GRAPH, 128);
    }

    #[test]
    fn six_r37_canonical_stage_snapshots_after_6r33_6r36() {
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
        let ref_bases = graph_ref.bases.as_slice();
        assert_eq!(ref_bases.len(), 430);
        let off_ct = (SITE_CT - JAVA_EXTENDED.0) as usize;
        let off_cg = (SITE_CG - JAVA_EXTENDED.0) as usize;
        let off_ca = (SITE_CA - JAVA_EXTENDED.0) as usize;
        let off_tc = (SITE_TC - JAVA_EXTENDED.0) as usize;
        let off_gc = (SITE_GC - JAVA_EXTENDED.0) as usize;
        assert_eq!(ref_bases[off_ct], b'C');
        assert_eq!(ref_bases[off_cg], b'C');
        assert_eq!(ref_bases[off_ca], b'C');
        assert_eq!(ref_bases[off_tc], b'T');
        assert_eq!(ref_bases[off_gc], b'G');

        let combined_361_371 = combined_alt_kmers(ref_bases, K, &[(off_ct, b'T'), (off_cg, b'G')]);
        let single_361_t = snp_alt_kmers(ref_bases, off_ct, K, b'T');
        let single_371_g = snp_alt_kmers(ref_bases, off_cg, K, b'G');
        let oracle_399_a = snp_alt_kmers(ref_bases, off_ca, K, b'A');
        assert!(!combined_361_371.is_empty());

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
        let (mut graph, _) = assembly_graph_from_ref_and_reads_threading_with_summary(
            &graph_ref,
            &assembly_reads,
            &graph_params(K),
        )
        .expect("raw rt");

        let snap = |stage, g: &AssemblyGraph| {
            snapshot(stage, g, &combined_361_371, &single_361_t, &single_371_g)
        };

        let raw = snap("raw", &graph);
        dump_stage(&graph, &raw);
        assert_eq!(raw.node_count, 518);
        assert!(raw.target_present);
        assert!(raw.alt_head_present);
        assert_eq!(raw.target_vertices.len(), 15);

        let mut pruning = AssemblyGraphPruningParams::gatk_haplotype_caller_defaults();
        pruning.min_prune_factor = rt_args.min_prune_factor;
        graph.apply_pruning(&pruning);
        let post_prune = snap("post_prune", &graph);
        dump_stage(&graph, &post_prune);
        assert_eq!(post_prune.node_count, 516);
        assert_eq!(post_prune.edge_count, 516);
        assert_eq!(post_prune.n_sources, 2);
        assert_eq!(post_prune.n_sinks, 1);
        assert_eq!(post_prune.from_src, 444);
        assert_eq!(post_prune.to_sink, 516);
        assert_eq!(post_prune.both, 444);
        assert_eq!(post_prune.neither, 0);
        assert!(post_prune.target_present);
        assert_eq!(post_prune.target_vertices.len(), 15);
        assert!(post_prune.alt_head_present);
        assert_eq!(fnv64(&kmer_blob(&sorted_graph_kmers(&graph))), JAVA_0_1_FNV);

        let from_src_prune = from_src_set(&graph);
        let to_sink_prune = to_sink_set(&graph);
        for &v in &post_prune.target_vertices {
            assert!(
                !from_src_prune.contains(&v),
                "TAGAGTTGAAG vertex {v} must not be reachable from REF source at prune"
            );
            assert!(
                to_sink_prune.contains(&v),
                "TAGAGTTGAAG vertex {v} must be able to reach REF sink at prune"
            );
            assert_eq!(graph.incoming_count(v), 1);
            assert_eq!(graph.outgoing_nodes(v).len(), 1);
        }
        let head = graph.vertex_id_for_kmer(JAVA_ALT_HEAD).expect("alt head");
        assert_eq!(graph.incoming_count(head), 0);
        assert!(!graph.is_ref_source_vertex(head));
        assert!(!from_src_prune.contains(&head));
        assert!(to_sink_prune.contains(&head));
        assert_eq!(
            graph.kmer_at(graph.reference_source_vertex().unwrap()),
            JAVA_REF_SOURCE
        );
        assert_eq!(
            graph.kmer_at(graph.reference_sink_vertex().unwrap()),
            JAVA_REF_SINK
        );

        let dangling = DanglingRecoveryParams::from_assembler_args(&rt_args);
        let dump = graph.test_dangling_head_decision_dump(head, &dangling);
        eprintln!(
            "HEAD_DECISION id={} kmer={} classified={} src_reach={} sink_reach={} path_find={} min_len={} not_sink={} aln={} cigar_ok={} prefix_rust={} plan={} final={} alt_len={} ref_len={} cigar={} firstM={} mm={} cap={} rust_idx={:?} java_idx={}",
            dump.head_id,
            ascii(&dump.head_kmer),
            dump.classified_dangling_head,
            dump.source_reachable,
            dump.sink_reachable,
            dump.path_find,
            dump.min_length,
            dump.not_ref_sink,
            dump.alignment,
            dump.cigar_ok,
            dump.prefix_legacy_rust,
            dump.rust_plan,
            dump.final_rust,
            dump.rust_alt_bases.len(),
            dump.rust_ref_bases.len(),
            dump.rust_cigar,
            dump.first_m_len,
            dump.mismatches_in_first_m,
            dump.max_mismatches_legacy,
            dump.rust_idx,
            dump.java_idx_on_rust_seq
        );
        let mm_pos: Vec<usize> = (0..dump
            .first_m_len
            .min(dump.rust_ref_bases.len())
            .min(dump.rust_alt_bases.len()))
            .filter(|&i| dump.rust_ref_bases[i] != dump.rust_alt_bases[i])
            .collect();
        eprintln!("HEAD_MM_POS {mm_pos:?}");
        eprintln!("HEAD_ALT_BYTES {}", ascii(&dump.rust_alt_bases));
        assert!(dump.classified_dangling_head);
        assert!(!dump.source_reachable);
        assert!(dump.sink_reachable);
        assert_eq!(dump.path_find, "PASS");
        assert_eq!(dump.min_length, "PASS");
        assert_eq!(dump.not_ref_sink, "PASS");
        assert_eq!(dump.alignment, "PASS");
        assert_eq!(dump.cigar_ok, "PASS");
        assert_eq!(dump.rust_alt_bases, dump.java_alt_bases);
        assert_eq!(dump.rust_ref_bases, dump.java_ref_bases);
        assert_eq!(dump.rust_alt_bases.len(), 97);
        assert_eq!(&dump.rust_alt_bases[..73], RUST_ALT_SUFFIX);
        assert_eq!(&dump.rust_alt_bases[73..], JAVA_ALT_EXTRA);
        assert_eq!(dump.rust_ref_bases.len(), 235);
        assert_eq!(dump.first_m_len, 97);
        assert_eq!(dump.mismatches_in_first_m, 10);
        assert_eq!(dump.max_mismatches_legacy, 3);
        assert_eq!(dump.prefix_legacy_rust, "FAIL none");
        assert_eq!(dump.rust_idx, None);
        assert_eq!(dump.java_idx_on_rust_seq, -1);
        assert_eq!(dump.final_rust, "REJECT");
        assert_eq!(dump.rust_plan, "prefix_match_legacy_failed");
        assert!(dump.merge_from_kmer.is_empty());
        assert_eq!(mm_pos.len(), 10);

        let mut post_tails = None;
        let mut post_heads = None;
        let summary = graph
            .test_java_exact_dangling_tails_then_heads(
                &dangling,
                |g, attempted, recovered| {
                    let s = snap("post_dangling_tails", g);
                    eprintln!("TAILS attempted={attempted} recovered={recovered}");
                    dump_stage(g, &s);
                    post_tails = Some(s);
                },
                |g, attempted, recovered| {
                    let s = snap("post_dangling_heads", g);
                    eprintln!("HEADS attempted={attempted} recovered={recovered}");
                    dump_stage(g, &s);
                    post_heads = Some(s);
                },
            )
            .expect("traced dangling");
        let post_tails = post_tails.expect("tails");
        let post_heads = post_heads.expect("heads");
        let pre_remove = snap("pre_remove_paths", &graph);
        dump_stage(&graph, &pre_remove);
        eprintln!(
            "DANGLING tails={}/{} heads={}/{}",
            summary.tails_recovered,
            summary.tails_attempted,
            summary.heads_recovered,
            summary.heads_attempted
        );
        assert_eq!(summary.tails_attempted, 0);
        assert_eq!(summary.tails_recovered, 0);
        assert_eq!(summary.heads_attempted, 1);
        assert_eq!(summary.heads_recovered, 0);
        assert_eq!(post_tails.node_count, 516);
        assert_eq!(post_tails.edge_count, 516);
        assert_eq!(post_tails.from_src, 444);
        assert!(post_tails.target_present);
        assert_eq!(post_heads.node_count, 516);
        assert_eq!(post_heads.edge_count, 516);
        assert_eq!(post_heads.from_src, 444);
        assert!(post_heads.target_present);
        assert_eq!(pre_remove.node_count, 516);
        assert_eq!(pre_remove.from_src, 444);
        assert!(pre_remove.target_present);
        assert_eq!(fnv64(&kmer_blob(&sorted_graph_kmers(&graph))), JAVA_0_1_FNV);

        let prune_kmers: HashSet<String> = sorted_graph_kmers(&graph).into_iter().collect();
        graph
            .remove_paths_not_connected_to_ref()
            .expect("removePaths");
        let post_remove = snap("post_remove_paths", &graph);
        dump_stage(&graph, &post_remove);
        let cleaned_kmers = sorted_graph_kmers(&graph);
        let cleaned_set: HashSet<String> = cleaned_kmers.iter().cloned().collect();
        let java_cleaned: HashSet<String> = parse_lines(JAVA_0_2_KMER_TXT).into_iter().collect();
        let java_lost: HashSet<String> = parse_lines(JAVA_0_1_MINUS_0_2_TXT).into_iter().collect();
        let rust_lost: HashSet<String> = prune_kmers.difference(&cleaned_set).cloned().collect();
        let rust_only = set_diff(&cleaned_set, &java_cleaned);
        let java_only = set_diff(&java_cleaned, &cleaned_set);
        eprintln!(
            "SET_CMP cleaned rust={} java={} intersection={} rust_only={} java_only={} rust_lost={} java_lost={}",
            cleaned_set.len(),
            java_cleaned.len(),
            cleaned_set.intersection(&java_cleaned).count(),
            rust_only.len(),
            java_only.len(),
            rust_lost.len(),
            java_lost.len()
        );
        if !rust_only.is_empty() {
            eprintln!("RUST_ONLY {:?}", rust_only);
        }
        if !java_only.is_empty() {
            eprintln!("JAVA_ONLY {:?}", java_only);
        }
        assert_eq!(post_remove.node_count, 444);
        assert_eq!(post_remove.edge_count, 444);
        assert_eq!(post_remove.n_sources, 1);
        assert_eq!(post_remove.n_sinks, 1);
        assert_eq!(post_remove.from_src, 444);
        assert_eq!(post_remove.to_sink, 444);
        assert_eq!(post_remove.both, 444);
        assert_eq!(post_remove.neither, 0);
        assert!(!post_remove.target_present);
        assert!(!post_remove.alt_head_present);
        assert_eq!(post_remove.n_361_371_combined, 0);
        assert_eq!(post_remove.target_vertices.len(), 0);
        assert_eq!(fnv64(&kmer_blob(&cleaned_kmers)), JAVA_0_2_FNV);
        assert_eq!(cleaned_set, java_cleaned);
        assert_eq!(rust_lost, java_lost);
        assert!(rust_only.is_empty());
        assert!(java_only.is_empty());
        assert_eq!(rust_lost.len(), 72);

        eprintln!(
            "TAGAGTTGAAG raw={} prune={} tails={} heads={} pre_remove={} post_remove={}",
            raw.target_present,
            post_prune.target_present,
            post_tails.target_present,
            post_heads.target_present,
            pre_remove.target_present,
            post_remove.target_present
        );
        assert!(raw.target_present);
        assert!(post_prune.target_present);
        assert!(post_tails.target_present);
        assert!(post_heads.target_present);
        assert!(pre_remove.target_present);
        assert!(!post_remove.target_present);

        assert_eq!(
            post_prune.n_361_371_combined, 4,
            "k=25 windows covering both 361 and 371"
        );
        assert_eq!(
            post_prune.n_361_single_t, 0,
            "no 361T with reference 371C; 361T only exists on the combined 361T+371G chain"
        );
        assert_eq!(
            post_prune.n_371_single_g, 10,
            "k=25 windows covering 371 but not 361 (sites are 10 bp apart)"
        );
        let present_371g: HashSet<String> = single_371_g
            .iter()
            .map(|k| ascii(k))
            .filter(|k| prune_kmers.contains(k))
            .collect();
        let present_combined: HashSet<String> = combined_361_371
            .iter()
            .map(|k| ascii(k))
            .filter(|k| prune_kmers.contains(k))
            .collect();
        assert_eq!(present_371g.len(), 10);
        assert_eq!(present_combined.len(), 4);
        assert!(
            present_371g.is_subset(&java_lost),
            "371G-only windows are on the dangling-head chain Java 0.2 deletes"
        );
        assert!(
            present_combined.is_subset(&java_lost),
            "combined 361T+371G windows are on the dangling-head chain Java 0.2 deletes"
        );
        assert_eq!(post_remove.n_361_371_combined, 0);
        assert_eq!(post_remove.n_371_single_g, 0);
        assert_eq!(post_remove.n_361_single_t, 0);
        let n_oracle_399 = exact_kmer_ids(&graph, &oracle_399_a).len();
        assert!(
            n_oracle_399 > 0,
            "cleaned RT must still encode 92317399 C/A"
        );

        let extra_vs_java = post_remove.node_count as i32 - 444;
        assert_eq!(
            extra_vs_java, 0,
            "historical +35 (479-444) must be gone after 6R.33+6R.36"
        );

        let cleaned_prod = build_threading_graph_for_seq_assembly(
            &graph_ref,
            &assembly_reads,
            K,
            &rt_args,
            false,
            false,
        )
        .expect("rt")
        .expect("k=25 graph");
        assert_eq!(cleaned_prod.node_count(), 444);
        assert_eq!(
            fnv64(&kmer_blob(&sorted_graph_kmers(&cleaned_prod))),
            JAVA_0_2_FNV
        );
        assert!(vertices_with_motif(&cleaned_prod, TARGET).is_empty());

        let mut seq = SeqGraph::from_assembly_graph(&cleaned_prod);
        let seq_from_rt = sorted_seq_payloads(&seq);
        eprintln!(
            "SEQ[from_rt] nodes={} edges={} st_paths={} fnv={:016x}",
            seq.node_count(),
            seq.edge_count(),
            count_st_paths(&seq),
            fnv64(&kmer_blob(&seq_from_rt))
        );
        assert_eq!(seq.node_count(), 444, "SeqGraph input must be cleaned RT");

        let mut zip_payloads = None;
        let status = seq.traced_cleanup_seq_graph(|stage, g| {
            if stage == "after_initial_zip_linear_chains" {
                zip_payloads = Some(sorted_seq_payloads(g));
                eprintln!(
                    "SEQ[zip] nodes={} edges={} st_paths={}",
                    g.node_count(),
                    g.edge_count(),
                    count_st_paths(g)
                );
            }
        });
        eprintln!("CLEANUP_STATUS={status:?}");
        let zip_payloads = zip_payloads.expect("zip snapshot");
        let java_zip: Vec<String> = parse_lines(JAVA_1_1_SEQ_TXT);
        assert_eq!(zip_payloads, java_zip);
        assert_eq!(fnv64(&kmer_blob(&zip_payloads)), JAVA_1_1_FNV);

        let final_payloads = sorted_seq_payloads(&seq);
        let java_final: Vec<String> = parse_lines(JAVA_1_4_SEQ_TXT);
        assert_eq!(seq.node_count(), 4);
        assert_eq!(count_st_paths(&seq), 2);
        assert_eq!(final_payloads, java_final);
        assert_eq!(fnv64(&kmer_blob(&final_payloads)), JAVA_1_4_FNV);
        assert!(!final_payloads.iter().any(|s| s.contains("TAGAGTTGAAG")));

        let paths = find_best_haplotypes_seq_graph(&seq, DEFAULT_NUM_BEST_HAPLOTYPES_PER_GRAPH)
            .expect("kbest");
        let mut uniq = HashSet::new();
        let mut n_361_t = 0usize;
        for (i, p) in paths.iter().enumerate() {
            let b = seq.path_bases_bytes(p.start, &p.edges);
            uniq.insert(b.clone());
            if b.len() > off_ct && b[off_ct] == b'T' {
                n_361_t += 1;
            }
            eprintln!(
                "KBEST[{i}] len={} score={:.6} ref={} 361={} 371={} 399={}",
                b.len(),
                p.score,
                p.is_reference,
                b.get(off_ct).copied().unwrap_or(0) as char,
                b.get(off_cg).copied().unwrap_or(0) as char,
                b.get(off_ca).copied().unwrap_or(0) as char
            );
        }
        assert_eq!(paths.len(), 2);
        assert_eq!(uniq.len(), 2);
        assert_eq!(n_361_t, 0);
        let has_399_a = paths.iter().any(|p| {
            let b = seq.path_bases_bytes(p.start, &p.edges);
            b.get(off_ca) == Some(&b'A')
        });
        let has_407_c = paths.iter().any(|p| {
            let b = seq.path_bases_bytes(p.start, &p.edges);
            b.get(off_tc) == Some(&b'C')
        });
        let has_412_c = paths.iter().any(|p| {
            let b = seq.path_bases_bytes(p.start, &p.edges);
            b.get(off_gc) == Some(&b'C')
        });
        assert!(has_399_a && has_407_c && has_412_c);
    }
}
