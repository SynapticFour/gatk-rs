//! SeqGraph vertex/edge ID invariant tests (structural repair 6R.2).
//!
//! Intended contracts:
//! - `SeqVertex.id` is a dense index into `vertices` (`seq_graph.rs` type docs).
//! - After every vertex compaction, `id == index` and every `edge.from`/`to` is a live vertex.
//! - Source/sink scan `0..vertices.len()` and look up adjacency by that same integer.
//!
//! W-H1 remains open. These tests do not claim Java SeqGraph parity.

#[cfg(test)]
mod traces {
    use super::super::*;
    use crate::assembly::AssemblyRead;
    use crate::assembly_graph_dump::{load_assembly_reads_tsv, load_assembly_ref_tsv};
    use crate::read_threading_assembler::{
        assemble_from_ref_and_reads, build_threading_graph_for_haplotype_dump,
        build_threading_graph_for_seq_assembly, ReadThreadingAssemblerArgs,
    };
    use crate::seq_kbest_haplotype::find_best_haplotypes_seq_graph;
    use std::collections::{HashMap, HashSet};

    #[derive(Debug, Clone)]
    struct IdBreak {
        stage: &'static str,
        n_vertices: usize,
        n_edges: usize,
        live_ids: Vec<usize>,
        duplicate_ids: Vec<usize>,
        dangling: Vec<(usize, usize, bool)>,
        id_index_mismatch: Vec<(usize, usize)>,
        src: Option<usize>,
        sink: Option<usize>,
        max_live_id: Option<usize>,
        max_endpoint: Option<usize>,
    }

    impl IdBreak {
        fn ok(&self) -> bool {
            self.duplicate_ids.is_empty()
                && self.dangling.is_empty()
                && self.id_index_mismatch.is_empty()
        }

        fn compact(&self) -> String {
            format!(
                "{} nv={} ne={} ids={:?} src={:?} sink={:?} dups={:?} dangling={:?} idx_mismatch={:?} max_id={:?} max_ep={:?}",
                self.stage,
                self.n_vertices,
                self.n_edges,
                self.live_ids,
                self.src,
                self.sink,
                self.duplicate_ids,
                self.dangling,
                self.id_index_mismatch,
                self.max_live_id,
                self.max_endpoint
            )
        }
    }

    fn inspect(stage: &'static str, g: &SeqGraph) -> IdBreak {
        let mut seen = HashSet::new();
        let mut duplicate_ids = Vec::new();
        let mut live_ids = Vec::new();
        let mut id_index_mismatch = Vec::new();
        for (i, v) in g.vertices.iter().enumerate() {
            live_ids.push(v.id);
            if !seen.insert(v.id) {
                duplicate_ids.push(v.id);
            }
            if v.id != i {
                id_index_mismatch.push((i, v.id));
            }
        }
        live_ids.sort_unstable();
        let live: HashSet<usize> = g.vertices.iter().map(|v| v.id).collect();
        let dangling: Vec<(usize, usize, bool)> = g
            .edges
            .iter()
            .filter(|e| !live.contains(&e.from) || !live.contains(&e.to))
            .map(|e| (e.from, e.to, e.is_ref))
            .collect();
        let max_live_id = live_ids.iter().copied().max();
        let max_endpoint = g.edges.iter().flat_map(|e| [e.from, e.to]).max();
        IdBreak {
            stage,
            n_vertices: g.vertices.len(),
            n_edges: g.edges.len(),
            live_ids,
            duplicate_ids,
            dangling,
            id_index_mismatch,
            src: g.reference_source_vertex(),
            sink: g.reference_sink_vertex(),
            max_live_id,
            max_endpoint,
        }
    }

    /// Test-only validator. Does not panic in production.
    fn assert_seq_graph_ids_and_edges_valid(g: &SeqGraph, stage: &str) {
        assert_graph_id_invariants(g, stage);
    }

    /// Test-only validator. Does not panic in production.
    fn assert_graph_id_invariants(g: &SeqGraph, stage: &str) {
        let r = inspect("check", g);
        assert!(
            r.duplicate_ids.is_empty(),
            "{stage}: duplicate vertex IDs {:?}",
            r.duplicate_ids
        );
        assert!(
            r.id_index_mismatch.is_empty(),
            "{stage}: vertex.id != index (implementation indexes vertices[v] with the same v used as ID): {:?}",
            r.id_index_mismatch
        );
        assert!(
            r.dangling.is_empty(),
            "{stage}: edges whose from/to are not live vertex IDs: {:?}",
            r.dangling
        );
        if let Some(s) = r.src {
            assert!(
                s < g.vertices.len(),
                "{stage}: source {s} out of vertex vector"
            );
        }
        if let Some(k) = r.sink {
            assert!(
                k < g.vertices.len(),
                "{stage}: sink {k} out of vertex vector"
            );
        }
        if g.vertices.len() == 2 {
            let dummy = g.vertices.iter().position(|v| v.sequence.is_empty());
            if let Some(d) = dummy {
                assert_eq!(
                    g.vertices[d].id, d,
                    "{stage}: dummy vertex id must equal its index"
                );
            }
        }
    }

    fn insert_dummy_like_cleanup(g: &mut SeqGraph) {
        if g.vertices.len() != 1 {
            return;
        }
        let complete = 0usize;
        let dummy_id = g.vertices.len();
        g.vertices.push(SeqVertex {
            id: dummy_id,
            sequence: Vec::new(),
        });
        g.edges.push(SeqEdge {
            from: complete,
            to: dummy_id,
            support: 0,
            is_ref: true,
        });
        g.rebuild_index();
    }

    const LEFT: &str = "ACGTACGGTTAGCCATAACGGTCCATTGCATAGCTGGAACCT";
    const RIGHT: &str = "GCTTAGGAACCGGTTAACCGATCCTGAACCGGATCCATAGCT";

    fn read(seq: &str, q: u8) -> AssemblyRead {
        AssemblyRead {
            bases: seq.as_bytes().to_vec(),
            base_quals: vec![q; seq.len()],
        }
    }

    fn synthetic_args() -> ReadThreadingAssemblerArgs {
        ReadThreadingAssemblerArgs {
            kmer_sizes: vec![10, 25],
            min_prune_factor: 2,
            allow_low_complexity_graphs: true,
            dont_increase_kmer_sizes_for_cycles: true,
            num_best_haplotypes_per_graph: 32,
            ..Default::default()
        }
    }

    struct Fixture {
        name: &'static str,
        reference: AssemblyRead,
        ref_seq: String,
        alt_seq: String,
    }

    fn fixture(name: &'static str, mid_ref: &str, mid_alt: &str) -> Fixture {
        let ref_seq = format!("{LEFT}{mid_ref}{RIGHT}");
        let alt_seq = format!("{LEFT}{mid_alt}{RIGHT}");
        Fixture {
            name,
            reference: read(&ref_seq, 30),
            ref_seq,
            alt_seq,
        }
    }

    fn reads_for(ref_seq: &str, alt_seq: &str) -> Vec<AssemblyRead> {
        let mut out = Vec::new();
        for _ in 0..4 {
            out.push(read(ref_seq, 30));
            out.push(read(alt_seq, 30));
        }
        out
    }

    fn variation_fixtures() -> Vec<Fixture> {
        vec![
            fixture("control_snp", "A", "C"),
            fixture("control_2bp_del", "TTCA", "TA"),
            fixture("control_2bp_ins", "TA", "TTCA"),
            fixture("control_two_snps", "AAAA", "ACCA"),
            fixture("control_del_plus_snp", "TTCA", "TAC"),
            fixture("holdout_a_TTCA_vs_TATG", "TTCA", "TATG"),
            fixture("holdout_b_GGAC_vs_GCTT", "GGAC", "GCTT"),
        ]
    }

    fn build_seq(
        reference: &AssemblyRead,
        reads: &[AssemblyRead],
        kmer: usize,
    ) -> Option<SeqGraph> {
        let args = synthetic_args();
        let graph = match build_threading_graph_for_seq_assembly(
            reference, reads, kmer, &args, true, false,
        ) {
            Ok(g) => g?,
            Err(_) => return None,
        };
        Some(SeqGraph::from_assembly_graph(&graph))
    }

    /// Walk production cleanup boundaries. Returns (first_break_stage, snaps).
    fn walk_cleanup(mut seq: SeqGraph) -> (Option<&'static str>, Vec<IdBreak>) {
        let mut snaps = Vec::new();
        let mut first_break: Option<&'static str> = None;
        let push = |stage: &'static str,
                    g: &SeqGraph,
                    snaps: &mut Vec<IdBreak>,
                    first: &mut Option<&'static str>| {
            let r = inspect(stage, g);
            if !r.ok() && first.is_none() {
                *first = Some(stage);
            }
            snaps.push(r);
        };

        push("constructed", &seq, &mut snaps, &mut first_break);
        seq.clean_non_ref_paths();
        push("after_clean_non_ref", &seq, &mut snaps, &mut first_break);
        seq.zip_linear_chains();
        push("after_zip", &seq, &mut snaps, &mut first_break);
        seq.remove_singleton_orphan_vertices();
        push("after_orphans", &seq, &mut snaps, &mut first_break);
        seq.remove_vertices_not_connected_to_ref_regardless_of_direction();
        push("after_undirected_prune", &seq, &mut snaps, &mut first_break);
        push("before_first_simplify", &seq, &mut snaps, &mut first_break);
        seq.simplify_graph();
        push("after_first_simplify", &seq, &mut snaps, &mut first_break);
        let first_jar =
            seq.reference_source_vertex().is_none() || seq.reference_sink_vertex().is_none();
        if first_jar {
            push(
                "first_simplify_jar_stop",
                &seq,
                &mut snaps,
                &mut first_break,
            );
            return (first_break, snaps);
        }
        let _ = seq.remove_paths_not_connected_to_ref();
        push("after_remove_paths", &seq, &mut snaps, &mut first_break);
        seq.simplify_graph();
        push("after_second_simplify", &seq, &mut snaps, &mut first_break);
        push("before_dummy", &seq, &mut snaps, &mut first_break);
        insert_dummy_like_cleanup(&mut seq);
        push("after_dummy", &seq, &mut snaps, &mut first_break);
        push("before_try_build", &seq, &mut snaps, &mut first_break);
        (first_break, snaps)
    }

    fn chain_like(vertices: usize, edges: &[(usize, usize)]) -> SeqGraph {
        let mut g = SeqGraph {
            kmer_size: 3,
            vertices: (0..vertices)
                .map(|i| SeqVertex {
                    id: i,
                    sequence: vec![b'A' + i as u8],
                })
                .collect(),
            edges: edges
                .iter()
                .map(|&(from, to)| SeqEdge {
                    from,
                    to,
                    support: 1,
                    is_ref: true,
                })
                .collect(),
            outgoing: HashMap::new(),
            incoming: HashMap::new(),
        };
        g.rebuild_index();
        g
    }

    #[test]
    fn compact_remap_rewrites_surviving_edge_0_to_3_as_0_to_2() {
        // Test A: vertices [0,1,2,3], edges 0->1, 1->2, 0->3; remove 1.
        let mut g = chain_like(4, &[(0, 1), (1, 2), (0, 3)]);
        assert_seq_graph_ids_and_edges_valid(&g, "a_before");
        let mut remove = HashSet::new();
        remove.insert(1);
        g.remove_vertices_by_id(&remove);
        let r = inspect("a_after", &g);
        eprintln!("{}", r.compact());
        assert_seq_graph_ids_and_edges_valid(&g, "a_after");
        assert_eq!(r.live_ids, vec![0, 1, 2]);
        assert!(
            r.dangling.is_empty(),
            "no dangling endpoints: {:?}",
            r.dangling
        );
        assert!(
            !g.edges.iter().any(|e| e.from == 3 || e.to == 3),
            "no edge may reference old vertex 3: {:?}",
            g.edges
        );
        assert!(
            g.edges.iter().any(|e| e.from == 0 && e.to == 2),
            "0->3 must rewrite to 0->2, got {:?}",
            g.edges
        );
        assert!(!g.edges.iter().any(|e| e.from == 0 && e.to == 1));
        assert!(!g.edges.iter().any(|e| e.from == 1 && e.to == 2));
    }

    #[test]
    fn compact_remap_discards_edge_whose_from_was_removed() {
        // Test B
        let mut g = chain_like(3, &[(0, 1), (1, 2)]);
        let mut remove = HashSet::new();
        remove.insert(1);
        g.remove_vertices_by_id(&remove);
        assert_seq_graph_ids_and_edges_valid(&g, "b_after");
        assert!(
            !g.edges.iter().any(|e| e.from == 1 || e.to == 1),
            "removed-from edges must be discarded: {:?}",
            g.edges
        );
        assert!(
            g.edges.is_empty(),
            "0->1 and 1->2 both touch removed 1: {:?}",
            g.edges
        );
    }

    #[test]
    fn compact_remap_is_noop_when_remove_set_is_empty() {
        // Test C
        let mut g = chain_like(3, &[(0, 1), (1, 2)]);
        let before_n = g.vertices.len();
        let before_e: Vec<_> = g.edges.clone();
        g.remove_vertices_by_id(&HashSet::new());
        assert_seq_graph_ids_and_edges_valid(&g, "c_after");
        assert_eq!(g.vertices.len(), before_n);
        assert_eq!(g.edges, before_e);
        assert!(g.vertices.iter().enumerate().all(|(i, v)| v.id == i));
    }

    #[test]
    fn dummy_id_one_collides_with_preexisting_stale_from_one() {
        // Documents dummy insertion on an *already invalid* graph. Production
        // prune must not produce this state after the remap repair.
        let mut g = SeqGraph {
            kmer_size: 10,
            vertices: vec![SeqVertex {
                id: 0,
                sequence: b"ACGT".to_vec(),
            }],
            edges: vec![SeqEdge {
                from: 1,
                to: 2,
                support: 1,
                is_ref: true,
            }],
            outgoing: HashMap::new(),
            incoming: HashMap::new(),
        };
        g.rebuild_index();
        insert_dummy_like_cleanup(&mut g);
        let after = inspect("after_dummy_collision", &g);
        assert!(after.sink.is_none());
    }

    #[test]
    fn merge_linear_chain_does_remap_edge_endpoints() {
        let mut g = SeqGraph {
            kmer_size: 3,
            vertices: vec![
                SeqVertex {
                    id: 0,
                    sequence: b"A".to_vec(),
                },
                SeqVertex {
                    id: 1,
                    sequence: b"C".to_vec(),
                },
                SeqVertex {
                    id: 2,
                    sequence: b"G".to_vec(),
                },
            ],
            edges: vec![
                SeqEdge {
                    from: 0,
                    to: 1,
                    support: 1,
                    is_ref: true,
                },
                SeqEdge {
                    from: 1,
                    to: 2,
                    support: 1,
                    is_ref: true,
                },
            ],
            outgoing: HashMap::new(),
            incoming: HashMap::new(),
        };
        g.rebuild_index();
        assert!(g.merge_linear_chain(&[0, 1, 2]));
        let r = inspect("after_zip_merge", &g);
        eprintln!("{}", r.compact());
        assert!(r.ok(), "zip merge_linear_chain must remap; {}", r.compact());
        assert_eq!(r.n_vertices, 1);
        assert_eq!(r.n_edges, 0);
    }

    fn snap_named<'a>(snaps: &'a [IdBreak], stage: &str) -> &'a IdBreak {
        snaps
            .iter()
            .find(|s| s.stage == stage)
            .unwrap_or_else(|| panic!("missing stage {stage}"))
    }

    fn production_cleanup_row(
        name: &str,
        reference: &AssemblyRead,
        reads: &[AssemblyRead],
        seq: SeqGraph,
    ) -> String {
        let (first, snaps) = walk_cleanup(seq.clone());
        eprintln!("=== {name} first_break={first:?} ===");
        for s in &snaps {
            eprintln!("{}", s.compact());
        }
        let zip = snap_named(&snaps, "after_zip");
        let orphans = snap_named(&snaps, "after_orphans");
        assert!(
            zip.ok(),
            "{name}: zip must keep valid IDs: {}",
            zip.compact()
        );
        assert!(
            orphans.ok(),
            "{name}: orphan removal must remap edges (no dangling): {}",
            orphans.compact()
        );
        assert!(
            !(orphans.n_vertices == 1 && !orphans.dangling.is_empty()),
            "{name}: must not collapse to nv=1 with stale edges: {}",
            orphans.compact()
        );
        assert!(
            first.is_none(),
            "{name}: ID invariant must hold through cleanup; first_break={first:?}"
        );

        let mut cleaned = seq;
        cleaned.clean_non_ref_paths();
        let status = cleaned.cleanup_seq_graph();
        assert_seq_graph_ids_and_edges_valid(&cleaned, &format!("{name}_full_cleanup"));
        let src = cleaned.reference_source_vertex();
        let sink = cleaned.reference_sink_vertex();
        let try_build = status == SeqGraphCleanupStatus::AssembledSomeVariation
            && src.is_some()
            && sink.is_some();
        let kbest = if try_build {
            find_best_haplotypes_seq_graph(&cleaned, 32)
                .map(|p| p.len())
                .unwrap_or(0)
        } else {
            0
        };

        let mut seq_args = synthetic_args();
        seq_args.use_seq_graph = true;
        let seq_asm = assemble_from_ref_and_reads(reference, reads, &seq_args);
        let (seq_status, seq_haps) = match &seq_asm {
            Ok(r) => (format!("{:?}", r.status), r.haplotypes.len()),
            Err(e) => (format!("err:{e}"), 0),
        };

        let mut rt_args = synthetic_args();
        rt_args.use_seq_graph = false;
        let rt_asm = assemble_from_ref_and_reads(reference, reads, &rt_args);
        let (rt_status, rt_haps) = match &rt_asm {
            Ok(r) => (format!("{:?}", r.status), r.haplotypes.len()),
            Err(e) => (format!("err:{e}"), 0),
        };

        let last = snaps.last().expect("snaps");
        format!(
            "{name}\tzip nv={} ne={} src={:?} sink={:?}\torphans nv={} ne={} src={:?} sink={:?} dangling={}\tcleanup nv={} ne={} src={:?} sink={:?} status={status:?}\ttry_build={try_build}\tkbest={kbest}\tseq_graph status={seq_status} haps={seq_haps}\trt status={rt_status} haps={rt_haps}\twalk_last nv={} ne={} src={:?} sink={:?}",
            zip.n_vertices,
            zip.n_edges,
            zip.src,
            zip.sink,
            orphans.n_vertices,
            orphans.n_edges,
            orphans.src,
            orphans.sink,
            orphans.dangling.len(),
            cleaned.vertices.len(),
            cleaned.edges.len(),
            src,
            sink,
            last.n_vertices,
            last.n_edges,
            last.src,
            last.sink,
        )
    }

    #[test]
    fn synthetic_matrix_id_invariants_hold_after_orphan_remap() {
        eprintln!(
            "stage order: constructed, clean_non_ref, zip, orphans, undirected, before/after first simplify, remove_paths, second simplify, dummy, try_build"
        );
        let mut rows = Vec::new();
        for fx in variation_fixtures() {
            let reads = reads_for(&fx.ref_seq, &fx.alt_seq);
            let Some(seq) = build_seq(&fx.reference, &reads, 10) else {
                panic!("{}: expected a threading graph at k=10", fx.name);
            };
            rows.push(production_cleanup_row(fx.name, &fx.reference, &reads, seq));
        }

        let ref_only = fixture("reference_only", "A", "A");
        let ref_reads: Vec<_> = (0..8).map(|_| read(&ref_only.ref_seq, 30)).collect();
        let seq = build_seq(&ref_only.reference, &ref_reads, 10)
            .unwrap_or_else(|| panic!("reference_only graph"));
        rows.push(production_cleanup_row(
            "reference_only",
            &ref_only.reference,
            &ref_reads,
            seq,
        ));

        eprintln!("=== synthetic matrix ===");
        for row in &rows {
            eprintln!("{row}");
        }
    }

    #[test]
    fn reference_only_zip_orphan_dummy_topology() {
        // Test D
        let ref_only = fixture("reference_only", "A", "A");
        let ref_reads: Vec<_> = (0..8).map(|_| read(&ref_only.ref_seq, 30)).collect();
        let mut seq = build_seq(&ref_only.reference, &ref_reads, 10)
            .unwrap_or_else(|| panic!("reference_only graph"));
        seq.clean_non_ref_paths();
        seq.zip_linear_chains();
        let zip = inspect("after_zip", &seq);
        eprintln!("{}", zip.compact());
        assert_eq!(zip.n_vertices, 1);
        assert_seq_graph_ids_and_edges_valid(&seq, "after_zip");
        seq.remove_singleton_orphan_vertices();
        let orphans = inspect("after_orphans", &seq);
        eprintln!("{}", orphans.compact());
        assert_eq!(orphans.n_vertices, 1);
        assert_seq_graph_ids_and_edges_valid(&seq, "after_orphans");
        insert_dummy_like_cleanup(&mut seq);
        let dummy = inspect("after_dummy", &seq);
        eprintln!("{}", dummy.compact());
        assert_eq!(dummy.n_vertices, 2);
        assert_eq!(dummy.n_edges, 1);
        assert!(dummy.src.is_some(), "valid source");
        assert!(dummy.sink.is_some(), "valid sink");
        assert!(seq.edges.iter().any(|e| e.is_ref));
        assert_seq_graph_ids_and_edges_valid(&seq, "after_dummy");
    }

    #[test]
    fn p5_case1_stays_id_valid_through_cleanup_and_keeps_sink() {
        // Test E — same topology/k-best/haplotypes as before the structural repair.
        let repo = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("..");
        let p5_ref =
            load_assembly_ref_tsv(&repo.join("parity/fixtures/hc-full-parity/e4/p5_case1_ref.tsv"))
                .unwrap_or_else(|e| panic!("p5 ref: {e}"));
        let p5_reads = load_assembly_reads_tsv(
            &repo.join("parity/fixtures/hc-full-parity/e4/p5_case1_reads.tsv"),
        )
        .unwrap_or_else(|e| panic!("p5 reads: {e}"));
        let p5_args = ReadThreadingAssemblerArgs {
            kmer_sizes: vec![3],
            min_base_quality: 10,
            min_prune_factor: 2,
            min_dangling_branch_length: 4,
            recover_dangling_heads: true,
            ..Default::default()
        };
        let p5_graph =
            build_threading_graph_for_haplotype_dump(&p5_ref, &p5_reads, 3, &p5_args, true, false)
                .unwrap_or_else(|e| panic!("p5 build: {e}"))
                .unwrap_or_else(|| panic!("p5 graph"));
        let seq = SeqGraph::from_assembly_graph(&p5_graph);
        let (first, snaps) = walk_cleanup(seq.clone());
        for s in &snaps {
            eprintln!("{}", s.compact());
        }
        assert!(
            first.is_none(),
            "p5_case1 must not violate ID invariants; first_break={first:?}"
        );
        let last = snaps.last().expect("snaps");
        assert!(last.src.is_some() && last.sink.is_some());
        assert_eq!(last.n_vertices, 2);
        assert_eq!(last.n_edges, 1);

        let mut cleaned = seq;
        cleaned.clean_non_ref_paths();
        let status = cleaned.cleanup_seq_graph();
        assert_eq!(status, SeqGraphCleanupStatus::AssembledSomeVariation);
        assert_eq!(cleaned.vertices.len(), 2);
        assert_eq!(cleaned.edges.len(), 1);
        assert_seq_graph_ids_and_edges_valid(&cleaned, "p5_cleanup");
        let paths = find_best_haplotypes_seq_graph(&cleaned, 32)
            .unwrap_or_else(|e| panic!("p5 k-best: {e}"));
        assert_eq!(
            paths.len(),
            1,
            "p5_case1 k-best path count changed; STOP for review"
        );
        let ref_path = cleaned
            .reference_path_bytes()
            .unwrap_or_else(|| panic!("p5 ref path"));
        assert_eq!(ref_path.as_slice(), p5_ref.bases.as_slice());
    }

    #[test]
    fn p12_production_seq_graph_waiver_assignments_unchanged() {
        let src = include_str!("assembly_based_caller.rs");
        let pin = "            assembler.use_seq_graph = false;\n            assembler.remove_paths_not_connected_to_ref = false;\n            assembler.skip_post_dangling_prune = true;";
        assert!(
            src.contains(pin),
            "P12 SeqGraph waiver (use_seq_graph=false / remove_paths=false / skip_post_dangling_prune=true) must remain under strict_java + region_overlaps_p12_cluster"
        );
        assert!(src.contains("region_overlaps_p12_cluster("));
        let waiver_idx = src.find(pin).expect("pin");
        let window = &src[waiver_idx.saturating_sub(400)..waiver_idx];
        assert!(
            window.contains("if args.strict_java_assembly"),
            "waiver must remain under strict_java_assembly"
        );
        assert!(
            window.contains("region_overlaps_p12_cluster"),
            "waiver must remain gated by region_overlaps_p12_cluster"
        );
    }
}
