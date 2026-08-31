//! 6R.32 TEST-ONLY: Java vs Rust `recoverDanglingHead` predicates for mid-B.
//! Does not change production dangling recovery scoring or acceptance.

#[cfg(test)]
mod traces {
    use crate::assembly::{AssemblyGraph, AssemblyGraphParams, AssemblyGraphPruningParams};
    use crate::assembly_based_caller::AssembleReadsArgs;
    use crate::assembly_dangling_recovery::{
        DanglingHeadDecisionDump, DanglingHeadPathVertexDump, DanglingRecoveryParams,
    };
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
    const JAVA_ALT_HEAD: &[u8] = b"CAAATAAAAGGTAGACAGCAGCATT";

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

    fn dump_path_vertex(label: &str, i: usize, v: &DanglingHeadPathVertexDump) {
        eprintln!(
            "  {label}[{i}] id={} kmer={} in={} out={} is_ref={} is_source={} is_sink={} ref_src={} ref_sink={} next_w={:?} n_out={}",
            v.id,
            ascii(&v.kmer),
            v.in_degree,
            v.out_degree,
            v.is_ref,
            v.is_source,
            v.is_sink,
            v.is_ref_source,
            v.is_ref_sink,
            v.edge_to_next_support,
            v.outgoing.len()
        );
        if v.outgoing.len() > 1 {
            for (t, kmer, w, is_ref) in &v.outgoing {
                eprintln!("    OUT id={t} kmer={} w={w} ref={is_ref}", ascii(kmer));
            }
        }
    }

    fn dump_decision(d: &DanglingHeadDecisionDump) {
        eprintln!("=== 6R.32 HEAD CANDIDATE ===");
        eprintln!("HEAD CANDIDATE");
        eprintln!("  head_id = {}", d.head_id);
        eprintln!("  first_vertex_sequence = {}", ascii(&d.head_kmer));
        eprintln!(
            "  classified_dangling_head = {}",
            d.classified_dangling_head
        );
        eprintln!("  source_reachable = {}", d.source_reachable);
        eprintln!("  sink_reachable = {}", d.sink_reachable);
        eprintln!("  path_length (vertices) = {}", d.alt_path_ids.len());
        eprintln!("  kmer_count = {}", d.alt_path_ids.len());
        eprintln!("  edge_count = {}", d.alt_path_ids.len().saturating_sub(1));
        eprintln!("  branch_weight (min edge support) = {:?}", d.branch_weight);
        eprintln!("  kmer_size = {}", d.kmer_size);
        eprintln!(
            "  min_dangling_branch_length = {} (min_vertices={})",
            d.min_dangling_branch_length, d.min_vertices
        );
        eprintln!("  prune_factor = {}", d.prune_factor);
        eprintln!("  give_up_at_branch = {}", d.give_up_at_branch);
        eprintln!(
            "  min_matching_bases_to_dangling_end_recovery = {}",
            d.min_matching_bases
        );
        eprintln!("  rust_plan = {}", d.rust_plan);
        eprintln!("  complete_candidate_path_ids = {:?}", d.alt_path_ids);
        eprintln!("ALT PATH VERTICES n={}", d.alt_path_vertices.len());
        for (i, v) in d.alt_path_vertices.iter().enumerate() {
            dump_path_vertex("ALT", i, v);
        }
        eprintln!("REFERENCE TARGET PATH ids={:?}", d.ref_path_ids);
        eprintln!("REF PATH VERTICES n={}", d.ref_path_vertices.len());
        for (i, v) in d.ref_path_vertices.iter().enumerate() {
            dump_path_vertex("REF", i, v);
        }
        eprintln!("PATH / SEQUENCE");
        eprintln!(
            "  rust_alt_len={} rust_ref_len={} java_alt_len={} java_ref_len={}",
            d.rust_alt_bases.len(),
            d.rust_ref_bases.len(),
            d.java_alt_bases.len(),
            d.java_ref_bases.len()
        );
        eprintln!("  rust_alt = {}", ascii(&d.rust_alt_bases));
        eprintln!("  rust_ref = {}", ascii(&d.rust_ref_bases));
        eprintln!("  java_alt = {}", ascii(&d.java_alt_bases));
        eprintln!("  java_ref = {}", ascii(&d.java_ref_bases));
        eprintln!(
            "  rust_vs_java_alt_same = {}",
            d.rust_alt_bases == d.java_alt_bases
        );
        eprintln!(
            "  rust_vs_java_ref_same = {}",
            d.rust_ref_bases == d.java_ref_bases
        );
        eprintln!("SW (Rust path_bases)");
        eprintln!("  cigar = {}", d.rust_cigar);
        eprintln!("  alignment_offset = {}", d.rust_alignment_offset);
        eprintln!(
            "  matches={} mismatches={} insertions={} deletions={}",
            d.sw_match, d.sw_mismatch, d.sw_ins, d.sw_del
        );
        eprintln!(
            "  score_source_derived = {} (not an acceptance predicate)",
            d.sw_score_source_derived
        );
        eprintln!(
            "  first_m_len={} matches_in_first_m={} mismatches_in_first_m={} max_mismatches_legacy={}",
            d.first_m_len,
            d.matches_in_first_m,
            d.mismatches_in_first_m,
            d.max_mismatches_legacy
        );
        eprintln!("SW (Java getBasesForPath encoding)");
        eprintln!("  cigar = {}", d.java_seq_cigar);
        eprintln!("  alignment_offset = {}", d.java_seq_alignment_offset);
        eprintln!("  cigar_ok = {}", d.java_seq_cigar_ok);
        eprintln!(
            "  first_m_len={} mismatches_in_first_m={}",
            d.java_seq_first_m_len, d.java_seq_mismatches_in_first_m
        );
        eprintln!("ACCEPTANCE");
        eprintln!("  path_find: {}", d.path_find);
        eprintln!("  min_length: {}", d.min_length);
        eprintln!("  not_ref_sink: {}", d.not_ref_sink);
        eprintln!("  alignment: {}", d.alignment);
        eprintln!("  cigar_ok (head first-M, ≤3): {}", d.cigar_ok);
        eprintln!("  prefix_legacy_rust: {}", d.prefix_legacy_rust);
        eprintln!(
            "  prefix_legacy_java_on_rust_seq: {}",
            d.prefix_legacy_java_on_rust_seq
        );
        eprintln!(
            "  prefix_legacy_java_on_java_seq: {}",
            d.prefix_legacy_java_on_java_seq
        );
        eprintln!("  rust_idx = {:?}", d.rust_idx);
        eprintln!("  java_idx_on_rust_seq = {}", d.java_idx_on_rust_seq);
        eprintln!("  java_idx_on_java_seq = {}", d.java_idx_on_java_seq);
        eprintln!("  final_rust: {}", d.final_rust);
        eprintln!(
            "  final_java_source_derived: {}",
            d.final_java_source_derived
        );
        eprintln!(
            "  merge {} -> {}",
            ascii(&d.merge_from_kmer),
            ascii(&d.merge_to_kmer)
        );
    }

    fn post_prune_mid_b_graph() -> Option<(AssemblyGraph, DanglingRecoveryParams)> {
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
        assert_eq!(rt_args.min_prune_factor, 2);
        assert_eq!(rt_args.min_dangling_branch_length, 4);
        assert!(rt_args.recover_dangling_heads);
        assert_eq!(rt_args.min_matching_bases_to_dangling_end_recovery, -1);

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
        assert_eq!(graph.node_count(), 518);
        let mut pruning = AssemblyGraphPruningParams::gatk_haplotype_caller_defaults();
        pruning.min_prune_factor = rt_args.min_prune_factor;
        graph.apply_pruning(&pruning);
        assert_eq!(graph.node_count(), 516);
        let dangling = DanglingRecoveryParams::from_assembler_args(&rt_args);
        assert!(dangling.dangling_java_exact);
        Some((graph, dangling))
    }

    #[test]
    fn six_r32_mid_b_java_rust_dangling_head_parity() {
        let Some((graph, dangling)) = post_prune_mid_b_graph() else {
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
        eprintln!("DANGLING_HEADS n={} ids={heads:?}", heads.len());
        assert_eq!(heads.len(), 1, "post-prune mid-B has one dangling head");
        let head = heads[0];
        assert_eq!(
            graph.kmer_at(head),
            JAVA_ALT_HEAD,
            "unique dangling head is the Java 0.1 alt-head k-mer"
        );
        let target_on_graph = (0..graph.node_count())
            .any(|v| graph.kmer_at(v).windows(TARGET.len()).any(|w| w == TARGET));
        assert!(target_on_graph, "TAGAGTTGAAG present post-prune");

        let dump = graph.test_dangling_head_decision_dump(head, &dangling);
        dump_decision(&dump);

        let from_src = graph
            .reference_source_vertex()
            .map(|s| bfs_forward(&graph, s))
            .unwrap_or_default();
        assert_eq!(from_src.len(), 444);
        assert!(!from_src.contains(&head));

        assert!(dump.classified_dangling_head);
        assert!(!dump.source_reachable);
        assert!(dump.sink_reachable);
        assert_eq!(dump.path_find, "PASS");
        assert_eq!(dump.min_length, "PASS");
        assert_eq!(dump.not_ref_sink, "PASS");
        assert_eq!(dump.alignment, "PASS");
        assert_eq!(dump.cigar_ok, "PASS");
        assert_eq!(dump.final_rust, "REJECT");
        assert_eq!(
            dump.rust_plan, "prefix_match_legacy_failed",
            "6R.33: mismatch cap abort must reject this head (was ok_merge idx=35)"
        );
        assert_eq!(dump.rust_idx, None);
        assert_ne!(dump.rust_idx, Some(35));
        assert_eq!(dump.first_m_len, 97);
        assert_eq!(dump.max_mismatches_legacy, 3);
        assert_eq!(dump.mismatches_in_first_m, 10);
        assert!(
            dump.mismatches_in_first_m > dump.max_mismatches_legacy,
            "this candidate exceeds Java maxMismatchesInDanglingHead"
        );
        assert_eq!(dump.java_idx_on_rust_seq, -1);
        assert_eq!(dump.java_idx_on_java_seq, -1);
        assert_eq!(dump.final_java_source_derived, "REJECT");
        assert_eq!(
            dump.rust_alt_bases, dump.java_alt_bases,
            "6R.36: path_bases matches getBasesForPath"
        );
        assert_eq!(dump.branch_weight, Some(2));

        eprintln!(
            "VERDICT_HINT rust={} java_source_derived={} seqs_same={} first_m_mm={} rust_idx={:?} java_idx_rust_seq={} java_idx_java_seq={}",
            dump.final_rust,
            dump.final_java_source_derived,
            dump.rust_alt_bases == dump.java_alt_bases && dump.rust_ref_bases == dump.java_ref_bases,
            dump.mismatches_in_first_m,
            dump.rust_idx,
            dump.java_idx_on_rust_seq,
            dump.java_idx_on_java_seq
        );
    }
}
