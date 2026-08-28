//! 6R.22 TEST-ONLY: k=85 dangling-candidate audit for mid-B 92317399.
//! Does not change production dangling recovery, k, gates, EventMap, or W-H1.

#[cfg(test)]
mod traces {
    use crate::assembly::{AssemblyGraph, AssemblyGraphParams, AssemblyGraphPruningParams};
    use crate::assembly_dangling_recovery::DanglingRecoveryParams;
    use crate::assembly_region_finalize::{
        assembly_reference_read, finalize_region_reads_for_assembly,
        gatk_min_tail_quality_for_assembly, padded_reference_loc, records_to_assembly_reads,
    };
    use crate::bio_ids::KmerSize;
    use crate::read_projection::query_index_at_reference_position;
    use crate::read_threading_assembler::ReadThreadingAssemblerArgs;
    use crate::read_threading_graph::assembly_graph_from_ref_and_reads_threading_with_summary;
    use crate::walker_apply::{call_disposition, AssemblyRegionCallDisposition};
    use gatk_core::reference::{
        parse_intervals_cli_string, ReferenceWindowCache, SequenceDictionary,
    };
    use rust_htslib::bam::record::CigarString;
    use std::collections::{HashMap, HashSet};
    use std::path::Path;

    const SITE: u64 = 92_317_399;
    const K: usize = 85;

    fn fixture_paths() -> Option<(std::path::PathBuf, std::path::PathBuf)> {
        let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("..");
        let ref_path = root.join("parity/realworld/assets/hs37d5.simple.fa");
        let bam = root.join("parity/realworld/na12878_20k_b37/NA12878_20k.b37.bam");
        if !ref_path.is_file() || !bam.is_file() {
            return None;
        }
        Some((ref_path, bam))
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

    fn load_mid_b_assembly() -> Option<(
        crate::assembly::AssemblyRead,
        Vec<crate::assembly::AssemblyRead>,
        Vec<rust_htslib::bam::Record>,
        u64,
    )> {
        let (ref_fasta, bam) = fixture_paths()?;
        let dict = SequenceDictionary::from_fasta_path(&ref_fasta).ok()?;
        let specs = parse_intervals_cli_string(&dict, "2:92317000-92319000").ok()?;
        let walk = crate::walker_traversal::traverse_assembly_region_walker(
            &dict,
            &specs,
            &ref_fasta,
            &bam,
            &crate::read_model::ReadFilterParams::gatk_standard_hc(),
            &crate::walker_traversal::WalkerTraversalConfig::gatk_haplotype_caller_production(100),
        )
        .ok()?;
        let regions = crate::walker_traversal::flatten_assembly_regions(&walk);
        let region = regions.iter().find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= SITE
                && r.end.get() >= SITE
        })?;
        let mut ref_cache = ReferenceWindowCache::new(ref_fasta.clone(), 4);
        let reference = assembly_reference_read(&dict, &mut ref_cache, region).ok()?;
        let (pad_start, _) = padded_reference_loc(region, &dict);
        let finalized = finalize_region_reads_for_assembly(
            &region.reads,
            region,
            true,
            gatk_min_tail_quality_for_assembly(10),
            false,
        );
        let assembly_reads = records_to_assembly_reads(&finalized);
        Some((reference, assembly_reads, finalized, pad_start))
    }

    fn production_java_exact_dangling() -> DanglingRecoveryParams {
        let args = ReadThreadingAssemblerArgs {
            dangling_java_exact: true,
            ..Default::default()
        };
        DanglingRecoveryParams::from_assembler_args(&args)
    }

    fn undirected_component(graph: &AssemblyGraph, seeds: &[usize]) -> HashSet<usize> {
        let mut adj: Vec<Vec<usize>> = vec![Vec::new(); graph.node_count()];
        for e in graph.edges_sorted() {
            adj[e.from].push(e.to);
            adj[e.to].push(e.from);
        }
        let mut seen = HashSet::new();
        let mut stack = Vec::new();
        for &s in seeds {
            if seen.insert(s) {
                stack.push(s);
            }
        }
        while let Some(v) = stack.pop() {
            for &w in &adj[v] {
                if seen.insert(w) {
                    stack.push(w);
                }
            }
        }
        seen
    }

    /// Java `findPathUpwardsToLowestCommonAncestor` with `giveUpAtBranch=true` (recoverAll=false).
    fn java_tail_give_up_path(
        graph: &AssemblyGraph,
        sink: usize,
        prune: u32,
    ) -> Option<Vec<usize>> {
        let mut path = Vec::new();
        let mut visited = HashSet::new();
        let mut v = sink;
        let done = |g: &AssemblyGraph, n: usize| {
            g.incoming_count(n) != 1 || g.outgoing_nodes(n).len() >= 2
        };
        while !done(graph, v) {
            let preds = graph.incoming_nodes(v);
            if preds.len() != 1 {
                return None;
            }
            let p = preds[0];
            let w = graph.edge_support(p, v).unwrap_or(0);
            if w < prune {
                visited.extend(path.drain(..));
            } else {
                path.insert(0, v);
            }
            if path.contains(&p) || visited.contains(&p) {
                return None;
            }
            v = p;
        }
        path.insert(0, v);
        if graph.outgoing_nodes(v).len() > 1 {
            Some(path)
        } else {
            None
        }
    }

    /// Java `findPathDownwardsToHighestCommonDescendantOfReference` with `giveUpAtBranch=true`.
    fn java_head_give_up_path(
        graph: &AssemblyGraph,
        source: usize,
        prune: u32,
    ) -> Option<Vec<usize>> {
        let mut path = Vec::new();
        let mut visited = HashSet::new();
        let mut v = source;
        let done = |g: &AssemblyGraph, n: usize| {
            g.ref_nodes.contains(&n) || g.outgoing_nodes(n).len() != 1
        };
        while !done(graph, v) {
            let outs = graph.outgoing_nodes(v);
            if outs.is_empty() {
                return None;
            }
            let t = outs[0];
            let w = graph.edge_support(v, t).unwrap_or(0);
            if w < prune {
                visited.extend(path.drain(..));
            } else {
                path.insert(0, v);
            }
            if path.contains(&t) || visited.contains(&t) {
                return None;
            }
            v = t;
        }
        path.insert(0, v);
        if graph.ref_nodes.contains(&v) {
            Some(path)
        } else {
            None
        }
    }

    fn compact_kmer(k: &[u8]) -> String {
        let s = String::from_utf8_lossy(k);
        if s.len() <= 16 {
            s.into_owned()
        } else {
            format!("{}..{}", &s[..8], &s[s.len() - 8..])
        }
    }

    fn n_read_kmers_in_graph(
        graph: &AssemblyGraph,
        reads: &[crate::assembly::AssemblyRead],
        finalized: &[rust_htslib::bam::Record],
    ) -> usize {
        let mut n = 0usize;
        for (ri, ar) in reads.iter().enumerate() {
            let cigar = CigarString(finalized[ri].cigar().iter().copied().collect());
            let Some(qi) =
                query_index_at_reference_position(finalized[ri].pos(), &cigar, (SITE - 1) as i64)
            else {
                continue;
            };
            for st in overlapping_kmer_starts(ar.bases.len(), qi, K) {
                if graph.vertex_id_for_kmer(&ar.bases[st..st + K]).is_some() {
                    n += 1;
                }
            }
        }
        n
    }

    #[test]
    fn six_r22_k85_dangling_candidate_audit() {
        let Some((reference, assembly_reads, finalized, _pad)) = load_mid_b_assembly() else {
            eprintln!("Real-data mid-B comparison unavailable");
            return;
        };
        assert_eq!(assembly_reads.len(), 2);

        let params = AssemblyGraphParams {
            kmer_size: KmerSize::try_from_usize(K).expect("k"),
            min_base_quality: 10,
            min_edge_weight: 1,
            dangling_path_max_nodes: 0,
            max_haplotypes: 128,
            max_haplotype_bases: 4096,
            start_threading_only_at_existing_vertex: false,
        };
        let (mut graph, summary) = assembly_graph_from_ref_and_reads_threading_with_summary(
            &reference,
            &assembly_reads,
            &params,
        )
        .expect("k=85 threading graph");
        eprintln!(
            "RAW k={K} nodes={} edges={} low_complexity={}",
            graph.node_count(),
            graph.edge_count(),
            summary.is_low_complexity
        );
        assert!(
            graph.node_count() > 0,
            "k=85 graph must build (last-attempt allows non-unique REF)"
        );

        let mut pruning = AssemblyGraphPruningParams::gatk_haplotype_caller_defaults();
        pruning.min_prune_factor = 2;
        let edges_removed = graph.apply_pruning(&pruning);
        eprintln!(
            "PRUNE factor=2 edges_removed={edges_removed} nodes={} edges={}",
            graph.node_count(),
            graph.edge_count()
        );
        assert_eq!(
            edges_removed, 0,
            "6R.21: k=85 prune removes 0 edges on this input"
        );

        let mut seed_ids = Vec::new();
        let mut seed_kmers = Vec::new();
        for (ri, ar) in assembly_reads.iter().enumerate() {
            let cigar = CigarString(finalized[ri].cigar().iter().copied().collect());
            let Some(qi) =
                query_index_at_reference_position(finalized[ri].pos(), &cigar, (SITE - 1) as i64)
            else {
                continue;
            };
            for st in overlapping_kmer_starts(ar.bases.len(), qi, K) {
                let km = &ar.bases[st..st + K];
                if let Some(id) = graph.vertex_id_for_kmer(km) {
                    seed_ids.push(id);
                    if seed_kmers.len() < 2 {
                        seed_kmers.push(compact_kmer(km));
                    }
                }
            }
        }
        seed_ids.sort_unstable();
        seed_ids.dedup();
        eprintln!(
            "READ_SNP_KMERS vertices={} sample={:?}",
            seed_ids.len(),
            seed_kmers
        );
        assert!(
            !seed_ids.is_empty(),
            "ALT-containing read 85-mer exists before dangling"
        );
        let n_before = n_read_kmers_in_graph(&graph, &assembly_reads, &finalized);
        assert!(n_before > 0);

        let component = undirected_component(&graph, &seed_ids);
        let n_ref_in_comp = component
            .iter()
            .filter(|v| graph.ref_nodes.contains(v))
            .count();
        let mut component_edges = 0usize;
        let mut min_mult = u32::MAX;
        let mut max_mult = 0u32;
        for &v in &component {
            for n in graph.outgoing_nodes(v) {
                if component.contains(&n) {
                    component_edges += 1;
                    let w = graph.edge_support(v, n).unwrap_or(0);
                    min_mult = min_mult.min(w);
                    max_mult = max_mult.max(w);
                }
            }
        }
        let sinks: Vec<usize> = component
            .iter()
            .copied()
            .filter(|&v| graph.outgoing_nodes(v).is_empty() && !graph.is_ref_sink_vertex(v))
            .collect();
        let sources: Vec<usize> = component
            .iter()
            .copied()
            .filter(|&v| {
                graph.incoming_count(v) == 0
                    && !graph.outgoing_nodes(v).is_empty()
                    && !graph.is_ref_source_vertex(v)
            })
            .collect();
        let path_like = component_edges + 1 == component.len();
        eprintln!(
            "CA_COMPONENT vertices={} edges={} path_like={path_like} path_bases~{} min_mult={} max_mult={} ref_nodes_in_comp={n_ref_in_comp} tail_sinks={} head_sources={}",
            component.len(),
            component_edges,
            K + component.len().saturating_sub(1),
            if component_edges == 0 { 0 } else { min_mult },
            max_mult,
            sinks.len(),
            sources.len()
        );
        assert_eq!(
            n_ref_in_comp, 0,
            "C/A 85-mers form a component with no REF vertices"
        );
        assert!(
            !sinks.is_empty() || !sources.is_empty(),
            "C/A path must present a dangling tail and/or head candidate"
        );

        let dang = production_java_exact_dangling();
        eprintln!(
            "DANGLING_PARAMS min_len={} prune={} recover_all={} java_exact={} min_match={} heads={}",
            dang.min_dangling_branch_length,
            dang.min_prune_factor,
            dang.recover_all_dangling_branches,
            dang.dangling_java_exact,
            dang.min_matching_bases_to_dangling_end_recovery,
            dang.recover_dangling_heads
        );
        assert_eq!(dang.min_dangling_branch_length, 4);
        assert_eq!(dang.min_prune_factor, 2);
        assert!(!dang.recover_all_dangling_branches);
        assert!(dang.dangling_java_exact);
        assert_eq!(dang.min_matching_bases_to_dangling_end_recovery, -1);

        let tail_probes: HashMap<usize, String> = graph
            .probe_dangling_tail_failures(&dang)
            .into_iter()
            .map(|(v, _, r)| (v, r))
            .collect();
        let mut n_ca_logged = 0usize;
        for &sink in &sinks {
            n_ca_logged += 1;
            let kmer = compact_kmer(graph.kmer_at(sink));
            let in_d = graph.incoming_count(sink);
            let out_d = graph.outgoing_nodes(sink).len();
            let java_path = java_tail_give_up_path(&graph, sink, dang.min_prune_factor);
            let rust_reason = tail_probes
                .get(&sink)
                .cloned()
                .unwrap_or_else(|| "missing_probe".into());
            let min_vertices = dang.min_dangling_branch_length.max(1) + 1;
            let java_len = java_path.as_ref().map(|p| p.len());
            let java_len_ok = java_len.map(|n| n >= min_vertices);
            let hits_ref = java_path
                .as_ref()
                .map(|p| p.iter().any(|n| graph.ref_nodes.contains(n)))
                .unwrap_or(false);
            let sw = rust_reason.starts_with("ok_merge")
                || rust_reason == "cigar_not_ok"
                || rust_reason == "matching_suffix_zero"
                || rust_reason == "matching_suffix_below_min"
                || rust_reason == "ref_index_zero_cycle"
                || rust_reason == "merge_index_oob"
                || rust_reason == "edge_exists";
            eprintln!(
                "CANDIDATE tail sink={sink} kmer={kmer} in={in_d} out={out_d} \
                 java_give_up_path={} java_len={java_len:?} min_vertices={min_vertices} \
                 java_len_ok={java_len_ok:?} hits_ref={hits_ref} sw_reached={sw} \
                 rust_reason={rust_reason} merge={}",
                java_path.is_some(),
                if rust_reason.starts_with("ok_merge") {
                    "ACCEPT"
                } else {
                    "REJECT"
                }
            );
            assert!(
                java_path.is_none(),
                "Java giveUpAtBranch tail path must be null on a REF-free linear chain"
            );
            assert_eq!(
                rust_reason, "no_alt_path",
                "Rust must reject the C/A tail before SW (findPath null)"
            );
            assert!(!sw, "SW must not run when give-up path is null");
        }

        for &src in &sources {
            n_ca_logged += 1;
            let kmer = compact_kmer(graph.kmer_at(src));
            let java_path = java_head_give_up_path(&graph, src, dang.min_prune_factor);
            let rust_reason = graph
                .probe_dangling_head_at(src, &dang)
                .unwrap_or_else(|| "missing_probe".into());
            let min_vertices = dang.min_dangling_branch_length + 1;
            let java_len = java_path.as_ref().map(|p| p.len());
            let sw = rust_reason.starts_with("ok_merge")
                || rust_reason == "cigar_not_ok"
                || rust_reason == "prefix_match_failed"
                || rust_reason == "prefix_match_legacy_failed"
                || rust_reason == "prefix_match_zero"
                || rust_reason == "extend_path_failed";
            eprintln!(
                "CANDIDATE head source={src} kmer={kmer} in=0 out={} \
                 java_give_up_path={} java_len={java_len:?} min_vertices={min_vertices} \
                 hits_ref={} sw_reached={sw} rust_reason={rust_reason} merge={}",
                graph.outgoing_nodes(src).len(),
                java_path.is_some(),
                java_path
                    .as_ref()
                    .map(|p| p.iter().any(|n| graph.ref_nodes.contains(n)))
                    .unwrap_or(false),
                if rust_reason.starts_with("ok_merge") {
                    "ACCEPT"
                } else {
                    "REJECT"
                }
            );
            assert!(
                java_path.is_none(),
                "Java giveUpAtBranch head path must be null unless the walk hits a REF node"
            );
            assert_eq!(
                rust_reason, "no_alt_path",
                "Rust must reject the C/A head before SW (findPath null)"
            );
            assert!(!sw, "SW must not run when give-up path is null");
        }
        assert!(
            n_ca_logged > 0,
            "candidate decision must be recorded for the C/A component"
        );

        let mut after_dang = graph.clone();
        let summary = after_dang
            .recover_dangling_branches(&dang)
            .expect("dangling");
        let n_after_dang = n_read_kmers_in_graph(&after_dang, &assembly_reads, &finalized);
        eprintln!(
            "AFTER_DANGLING tails={}/{} heads={}/{} merge_haps={} read_snp_kmers={n_after_dang} (before={n_before})",
            summary.tails_recovered,
            summary.tails_attempted,
            summary.heads_recovered,
            summary.heads_attempted,
            after_dang.dangling_merge_haps.len()
        );
        assert_eq!(summary.tails_recovered, 0);
        assert_eq!(summary.heads_recovered, 0);
        assert!(after_dang.dangling_merge_haps.is_empty());
        assert_eq!(
            n_after_dang, n_before,
            "rejected dangling must not drop ALT kmers"
        );

        let connected_after = seed_ids.iter().any(|&id| {
            after_dang
                .vertex_id_for_kmer(graph.kmer_at(id))
                .is_some_and(|nid| {
                    after_dang.ref_nodes.contains(&nid)
                        || after_dang
                            .incoming_nodes(nid)
                            .iter()
                            .any(|&p| after_dang.edge_is_ref(p, nid))
                        || after_dang
                            .outgoing_nodes(nid)
                            .iter()
                            .any(|&t| after_dang.edge_is_ref(nid, t))
                })
        });
        eprintln!("REF_CONNECTED_AFTER_DANGLING={connected_after}");

        let mut after_rm = after_dang.clone();
        after_rm
            .remove_paths_not_connected_to_ref()
            .expect("remove_paths");
        let n_after_rm = n_read_kmers_in_graph(&after_rm, &assembly_reads, &finalized);
        eprintln!("AFTER_REMOVE_PATHS read_snp_kmers={n_after_rm}");

        if !connected_after {
            assert_eq!(
                n_after_rm, 0,
                "disconnected ALT must be removed by remove_paths_not_connected_to_ref"
            );
        }
    }
}
