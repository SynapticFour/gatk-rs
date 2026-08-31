//! Test-only first-simplify investigation (6R.3). Production behavior is not modified.
//!
//! Coordinate-free unique-flank synthetics. W-H1 remains open. No Java parity claim.

#[cfg(test)]
mod traces {
    use super::super::*;
    use crate::assembly::AssemblyRead;
    use crate::assembly_graph_dump::{load_assembly_reads_tsv, load_assembly_ref_tsv};
    use crate::read_threading_assembler::{
        build_threading_graph_for_haplotype_dump, build_threading_graph_for_seq_assembly,
        ReadThreadingAssemblerArgs,
    };
    use crate::seq_graph_simplify::traced_simplify_graph_full;
    use std::collections::{HashMap, HashSet, VecDeque};

    const LEFT: &str = "ACGTACGGTTAGCCATAACGGTCCATTGCATAGCTGGAACCT";
    const RIGHT: &str = "GCTTAGGAACCGGTTAACCGATCCTGAACCGGATCCATAGCT";

    #[derive(Debug, Clone)]
    struct Snap {
        stage: String,
        nv: usize,
        ne: usize,
        ids: Vec<usize>,
        dangling: usize,
        src: Option<usize>,
        sink: Option<usize>,
        ref_e: usize,
        alt_e: usize,
        components: usize,
        from_src: usize,
        to_sink: usize,
        both: usize,
        branches: Vec<usize>,
        joins: Vec<usize>,
        seq_lens: Vec<usize>,
        edge_dump: String,
        skeleton: String,
    }

    impl Snap {
        fn identity_ok(&self) -> bool {
            self.dangling == 0 && self.ids.iter().enumerate().all(|(i, id)| *id == i)
        }

        fn src_sink_ok(&self) -> bool {
            self.src.is_some() && self.sink.is_some()
        }

        fn line(&self) -> String {
            format!(
                "{} nv={} ne={} src={:?} sink={:?} ref_e={} alt_e={} dang={} cc={} from_src={} to_sink={} both={} branches={:?} joins={:?} id_ok={}",
                self.stage,
                self.nv,
                self.ne,
                self.src,
                self.sink,
                self.ref_e,
                self.alt_e,
                self.dangling,
                self.components,
                self.from_src,
                self.to_sink,
                self.both,
                self.branches,
                self.joins,
                self.identity_ok()
            )
        }
    }

    fn inspect(stage: impl Into<String>, g: &SeqGraph) -> Snap {
        let ids = g.test_vertex_ids();
        let live: HashSet<usize> = ids.iter().copied().collect();
        let dangling = g
            .edges_pub()
            .iter()
            .filter(|e| !live.contains(&e.from) || !live.contains(&e.to))
            .count();
        let ref_e = g.edges_pub().iter().filter(|e| e.is_ref).count();
        let alt_e = g.edges_pub().len().saturating_sub(ref_e);
        let n = g.node_count();
        let src = g.reference_source_vertex();
        let sink = g.reference_sink_vertex();
        let mut branches = Vec::new();
        let mut joins = Vec::new();
        for v in 0..n {
            if g.outgoing_nodes(v).len() > 1 {
                branches.push(v);
            }
            if g.incoming_nodes(v).len() > 1 {
                joins.push(v);
            }
        }
        let from_src = src.map(|s| directed_forward(g, s).len()).unwrap_or(0);
        let to_sink = sink.map(|k| directed_backward(g, k).len()).unwrap_or(0);
        let both = match (src, sink) {
            (Some(s), Some(k)) => {
                let a = directed_forward(g, s);
                let b = directed_backward(g, k);
                a.intersection(&b).count()
            }
            _ => 0,
        };
        let seq_lens = (0..n).map(|v| g.vertex_sequence(v).len()).collect();
        let edge_dump = g
            .edges_pub()
            .iter()
            .map(|e| format!("{}-{}>{}", e.from, if e.is_ref { "R" } else { "A" }, e.to))
            .collect::<Vec<_>>()
            .join(" ");
        Snap {
            stage: stage.into(),
            nv: n,
            ne: g.edge_count(),
            ids,
            dangling,
            src,
            sink,
            ref_e,
            alt_e,
            components: undirected_components(g),
            from_src,
            to_sink,
            both,
            branches,
            joins,
            seq_lens,
            edge_dump,
            skeleton: skeleton(g, src, sink),
        }
    }

    fn directed_forward(g: &SeqGraph, start: usize) -> HashSet<usize> {
        let mut seen = HashSet::new();
        let mut q = VecDeque::new();
        q.push_back(start);
        seen.insert(start);
        while let Some(v) = q.pop_front() {
            for t in g.outgoing_nodes(v) {
                if seen.insert(t) {
                    q.push_back(t);
                }
            }
        }
        seen
    }

    fn directed_backward(g: &SeqGraph, start: usize) -> HashSet<usize> {
        let mut seen = HashSet::new();
        let mut q = VecDeque::new();
        q.push_back(start);
        seen.insert(start);
        while let Some(v) = q.pop_front() {
            for p in g.incoming_nodes(v) {
                if seen.insert(p) {
                    q.push_back(p);
                }
            }
        }
        seen
    }

    fn undirected_components(g: &SeqGraph) -> usize {
        let n = g.node_count();
        let mut seen = vec![false; n];
        let mut cc = 0usize;
        for s in 0..n {
            if seen[s] {
                continue;
            }
            cc += 1;
            let mut q = vec![s];
            seen[s] = true;
            while let Some(v) = q.pop() {
                for t in g.outgoing_nodes(v) {
                    if !seen[t] {
                        seen[t] = true;
                        q.push(t);
                    }
                }
                for p in g.incoming_nodes(v) {
                    if !seen[p] {
                        seen[p] = true;
                        q.push(p);
                    }
                }
            }
        }
        cc
    }

    fn skeleton(g: &SeqGraph, src: Option<usize>, sink: Option<usize>) -> String {
        let mut parts = Vec::new();
        if let Some(s) = src {
            parts.push(format!("src={s}/len={}", g.vertex_sequence(s).len()));
        }
        if let Some(k) = sink {
            parts.push(format!("sink={k}/len={}", g.vertex_sequence(k).len()));
        }
        for v in 0..g.node_count() {
            let outs = g.outgoing_nodes(v);
            if outs.len() > 1 {
                let desc: Vec<String> = outs
                    .iter()
                    .map(|&t| {
                        format!(
                            "{}:{}:len{}",
                            if g.edge_is_ref(v, t) { "R" } else { "A" },
                            t,
                            g.vertex_sequence(t).len()
                        )
                    })
                    .collect();
                parts.push(format!("branch {v}->[{}]", desc.join(",")));
            }
        }
        for v in 0..g.node_count() {
            let ins = g.incoming_nodes(v);
            if ins.len() > 1 {
                let desc: Vec<String> = ins
                    .iter()
                    .map(|&p| {
                        format!(
                            "{}:{}:len{}",
                            if g.edge_is_ref(p, v) { "R" } else { "A" },
                            p,
                            g.vertex_sequence(p).len()
                        )
                    })
                    .collect();
                parts.push(format!("join {v}<-[{}]", desc.join(",")));
            }
        }
        parts.join(" | ")
    }

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

    fn to_post_orphan_undirected(mut seq: SeqGraph) -> SeqGraph {
        seq.clean_non_ref_paths();
        seq.zip_linear_chains();
        seq.remove_singleton_orphan_vertices();
        seq.remove_vertices_not_connected_to_ref_regardless_of_direction();
        seq
    }

    struct Walk {
        snaps: Vec<Snap>,
        first_src_sink_loss: Option<String>,
        first_identity_loss: Option<String>,
        first_nv_le_2: Option<String>,
    }

    fn walk_first_simplify(mut g: SeqGraph) -> Walk {
        let mut snaps = Vec::new();
        let mut first_src_sink_loss = None;
        let mut first_identity_loss = None;
        let mut first_nv_le_2 = None;
        let mut push = |stage: &str, g: &SeqGraph, snaps: &mut Vec<Snap>| {
            let s = inspect(stage, g);
            if first_identity_loss.is_none() && !s.identity_ok() {
                first_identity_loss = Some(s.stage.clone());
            }
            if first_src_sink_loss.is_none() && !s.src_sink_ok() {
                first_src_sink_loss = Some(s.stage.clone());
            }
            if first_nv_le_2.is_none() && s.nv <= 2 {
                first_nv_le_2 = Some(s.stage.clone());
            }
            snaps.push(s);
        };

        push("entering_first_simplify", &g, &mut snaps);
        traced_simplify_graph_full(&mut g, |stage, graph| {
            push(stage, graph, &mut snaps);
        });
        Walk {
            snaps,
            first_src_sink_loss,
            first_identity_loss,
            first_nv_le_2,
        }
    }

    fn report_fixture(name: &str, seq: SeqGraph) -> Walk {
        let g = to_post_orphan_undirected(seq);
        let enter = inspect("entering_first_simplify", &g);
        assert!(
            enter.identity_ok(),
            "{name}: identity must hold entering first simplify: {}",
            enter.line()
        );
        let mut prod = g.clone();
        let walk = walk_first_simplify(g);
        prod.simplify_graph();
        let traced_done = walk.snaps.last().expect("snaps");
        let prod_done = inspect("prod_simplify_graph", &prod);
        assert_eq!(
            traced_done.nv, prod_done.nv,
            "{name}: traced simplify must match production nv"
        );
        assert_eq!(
            traced_done.ne, prod_done.ne,
            "{name}: traced simplify must match production ne"
        );
        assert_eq!(
            traced_done.src, prod_done.src,
            "{name}: traced simplify must match production src"
        );
        assert_eq!(
            traced_done.sink, prod_done.sink,
            "{name}: traced simplify must match production sink"
        );
        eprintln!(
            "=== {name} first_src_sink_loss={:?} first_identity_loss={:?} first_nv_le_2={:?} ===",
            walk.first_src_sink_loss, walk.first_identity_loss, walk.first_nv_le_2
        );
        for s in &walk.snaps {
            eprintln!("{}", s.line());
            if s.nv <= 12 || s.stage.contains("entering") || s.branches.len() + s.joins.len() > 0 {
                eprintln!("  skeleton: {}", s.skeleton);
                if s.ne <= 48 {
                    eprintln!("  edges: {}", s.edge_dump);
                }
                if s.nv <= 8 {
                    eprintln!("  seq_lens={:?}", s.seq_lens);
                }
            }
        }
        walk
    }

    fn chain_like(vertices: usize, edges: &[(usize, usize, bool)]) -> SeqGraph {
        let mut g = SeqGraph {
            kmer_size: 3,
            vertices: (0..vertices)
                .map(|i| SeqVertex {
                    id: i,
                    sequence: vec![b'A' + (i as u8 % 26)],
                })
                .collect(),
            edges: edges
                .iter()
                .map(|&(from, to, is_ref)| SeqEdge {
                    from,
                    to,
                    support: 1,
                    is_ref,
                })
                .collect(),
            outgoing: HashMap::new(),
            incoming: HashMap::new(),
        };
        g.rebuild_index();
        g
    }

    #[test]
    fn zip_preserves_chain_exit_edge_into_join() {
        // 0 -R> 1 -R> 2(join) <-A- 3. Zipping 0,1 must keep merged → join.
        let mut g = chain_like(4, &[(0, 1, true), (1, 2, true), (3, 2, false)]);
        let before = inspect("before", &g);
        eprintln!("{}", before.line());
        eprintln!("  edges: {}", before.edge_dump);
        assert!(g.zip_linear_chains());
        let after = inspect("after_zip", &g);
        eprintln!("{}", after.line());
        eprintln!("  edges: {}", after.edge_dump);
        assert!(after.identity_ok());
        assert_eq!(after.components, 1, "join must stay connected");
        assert!(after.src_sink_ok());
        assert!(
            after.joins.len() == 1 || after.edge_dump.contains("-A>"),
            "alt incoming to join must survive: {}",
            after.edge_dump
        );
        let src = after.src.expect("src");
        let sink = after.sink.expect("sink");
        assert!(
            g.outgoing_nodes(src).contains(&sink) || after.both == after.nv,
            "merged chain must still reach the join/sink; edges={}",
            after.edge_dump
        );
    }

    #[test]
    fn diamond_with_shared_prefix_attaches_top_to_prefix() {
        // top=0 "AAAA", mids 1="XXT" 2="XXG", bottom=3 "CCCC".
        // Diamond merge should extract prefix XX and reattach top → prefix after compact.
        let mut g = SeqGraph {
            kmer_size: 3,
            vertices: vec![
                SeqVertex {
                    id: 0,
                    sequence: b"AAAA".to_vec(),
                },
                SeqVertex {
                    id: 1,
                    sequence: b"XXT".to_vec(),
                },
                SeqVertex {
                    id: 2,
                    sequence: b"XXG".to_vec(),
                },
                SeqVertex {
                    id: 3,
                    sequence: b"CCCC".to_vec(),
                },
            ],
            edges: vec![
                SeqEdge {
                    from: 0,
                    to: 1,
                    support: 4,
                    is_ref: true,
                },
                SeqEdge {
                    from: 0,
                    to: 2,
                    support: 4,
                    is_ref: false,
                },
                SeqEdge {
                    from: 1,
                    to: 3,
                    support: 4,
                    is_ref: true,
                },
                SeqEdge {
                    from: 2,
                    to: 3,
                    support: 4,
                    is_ref: false,
                },
            ],
            outgoing: HashMap::new(),
            incoming: HashMap::new(),
        };
        g.rebuild_index();
        let before = inspect("diamond_before", &g);
        eprintln!("{}", before.line());
        eprintln!("  edges: {}", before.edge_dump);
        assert!(before.src_sink_ok());
        assert!(before.identity_ok());
        let mut after_diamonds: Option<Snap> = None;
        traced_simplify_graph_full(&mut g, |stage, graph| {
            let s = inspect(stage, graph);
            eprintln!("{}", s.line());
            eprintln!("  skeleton: {}", s.skeleton);
            eprintln!("  edges: {}", s.edge_dump);
            if stage == "after_merge_diamonds" && after_diamonds.is_none() {
                after_diamonds = Some(s);
            }
        });
        let diamonds = after_diamonds.expect("merge_diamonds snap");
        assert_eq!(
            diamonds.from_src, diamonds.nv,
            "after merge_diamonds every vertex must be reachable from source; edges={}",
            diamonds.edge_dump
        );
        let after = inspect("diamond_after", &g);
        eprintln!("diamond_after {}", after.line());
        assert!(after.src_sink_ok(), "diamond must keep source/sink");
        assert!(after.identity_ok());
        assert_eq!(
            after.from_src, after.nv,
            "final diamond graph must stay source-connected; edges={}",
            after.edge_dump
        );
        let src = after.src.expect("src");
        assert!(
            g.outgoing_nodes(src).len() >= 2,
            "source should retain a branch after prefix zip, outs={:?} edges={}",
            g.outgoing_nodes(src),
            after.edge_dump
        );
    }

    #[test]
    fn diamond_empty_prefix_shared_suffix_keeps_ref_path() {
        // Deletion-shaped diamond: ref mid = extra 2bp + alt mid (shared suffix, empty prefix).
        // top=0 "AAAA", ref mid=1 "TTCGTACGTAC", alt mid=2 "CGTACGTAC", bottom=3 "GGGG".
        let mut g = SeqGraph {
            kmer_size: 3,
            vertices: vec![
                SeqVertex {
                    id: 0,
                    sequence: b"AAAA".to_vec(),
                },
                SeqVertex {
                    id: 1,
                    sequence: b"TTCGTACGTAC".to_vec(),
                },
                SeqVertex {
                    id: 2,
                    sequence: b"CGTACGTAC".to_vec(),
                },
                SeqVertex {
                    id: 3,
                    sequence: b"GGGG".to_vec(),
                },
            ],
            edges: vec![
                SeqEdge {
                    from: 0,
                    to: 1,
                    support: 4,
                    is_ref: true,
                },
                SeqEdge {
                    from: 0,
                    to: 2,
                    support: 4,
                    is_ref: false,
                },
                SeqEdge {
                    from: 1,
                    to: 3,
                    support: 4,
                    is_ref: true,
                },
                SeqEdge {
                    from: 2,
                    to: 3,
                    support: 4,
                    is_ref: false,
                },
            ],
            outgoing: HashMap::new(),
            incoming: HashMap::new(),
        };
        g.rebuild_index();
        let mut after_diamonds: Option<Snap> = None;
        traced_simplify_graph_full(&mut g, |stage, graph| {
            if stage == "after_merge_diamonds" && after_diamonds.is_none() {
                after_diamonds = Some(inspect(stage, graph));
            }
        });
        let diamonds = after_diamonds.expect("merge_diamonds");
        eprintln!("{}", diamonds.line());
        eprintln!("  edges: {}", diamonds.edge_dump);
        assert!(
            diamonds.src_sink_ok(),
            "empty-prefix diamond must keep src/sink: {}",
            diamonds.line()
        );
        assert!(
            diamonds.ref_e >= 1,
            "empty-prefix diamond must keep a reference edge: {}",
            diamonds.edge_dump
        );
        assert_eq!(
            diamonds.from_src, diamonds.nv,
            "every vertex reachable from source after diamonds: {}",
            diamonds.edge_dump
        );
        let after = inspect("after", &g);
        assert!(
            after.src_sink_ok(),
            "final graph src/sink: {}",
            after.line()
        );
        assert!(
            after.nv >= 3 && (after.branches.len() + after.joins.len() > 0 || after.alt_e > 0),
            "deletion bubble must survive simplify: {}",
            after.line()
        );
    }

    #[test]
    fn diamond_empty_prefix_remaining_on_alt_keeps_ref_path() {
        // Insertion-shaped diamond: alt mid = extra 2bp + ref mid.
        let mut g = SeqGraph {
            kmer_size: 3,
            vertices: vec![
                SeqVertex {
                    id: 0,
                    sequence: b"AAAA".to_vec(),
                },
                SeqVertex {
                    id: 1,
                    sequence: b"CGTACGTAC".to_vec(),
                },
                SeqVertex {
                    id: 2,
                    sequence: b"TTCGTACGTAC".to_vec(),
                },
                SeqVertex {
                    id: 3,
                    sequence: b"GGGG".to_vec(),
                },
            ],
            edges: vec![
                SeqEdge {
                    from: 0,
                    to: 1,
                    support: 4,
                    is_ref: true,
                },
                SeqEdge {
                    from: 0,
                    to: 2,
                    support: 4,
                    is_ref: false,
                },
                SeqEdge {
                    from: 1,
                    to: 3,
                    support: 4,
                    is_ref: true,
                },
                SeqEdge {
                    from: 2,
                    to: 3,
                    support: 4,
                    is_ref: false,
                },
            ],
            outgoing: HashMap::new(),
            incoming: HashMap::new(),
        };
        g.rebuild_index();
        let mut after_diamonds: Option<Snap> = None;
        traced_simplify_graph_full(&mut g, |stage, graph| {
            if stage == "after_merge_diamonds" && after_diamonds.is_none() {
                after_diamonds = Some(inspect(stage, graph));
            }
        });
        let diamonds = after_diamonds.expect("merge_diamonds");
        eprintln!("{}", diamonds.line());
        eprintln!("  edges: {}", diamonds.edge_dump);
        assert!(
            diamonds.src_sink_ok(),
            "insertion diamond must keep src/sink: {}",
            diamonds.line()
        );
        assert!(
            diamonds.ref_e >= 1,
            "insertion diamond must keep a reference edge: {}",
            diamonds.edge_dump
        );
        let after = inspect("after", &g);
        assert!(
            after.src_sink_ok(),
            "final graph src/sink: {}",
            after.line()
        );
        assert!(
            after.nv >= 3 && (after.branches.len() + after.joins.len() > 0 || after.alt_e > 0),
            "insertion bubble must survive simplify: {}",
            after.line()
        );
    }

    #[test]
    fn common_suffix_split_allowed_at_reference_sink() {
        // Java safeToSplit does not require bottom out-degree 1. Sink has out-degree 0.
        let mut g = SeqGraph {
            kmer_size: 3,
            vertices: vec![
                SeqVertex {
                    id: 0,
                    sequence: b"AAAA".to_vec(),
                },
                SeqVertex {
                    id: 1,
                    sequence: b"XXTAA".to_vec(),
                },
                SeqVertex {
                    id: 2,
                    sequence: b"YYTAA".to_vec(),
                },
                SeqVertex {
                    id: 3,
                    sequence: b"GGGG".to_vec(),
                },
            ],
            edges: vec![
                SeqEdge {
                    from: 0,
                    to: 1,
                    support: 2,
                    is_ref: true,
                },
                SeqEdge {
                    from: 0,
                    to: 2,
                    support: 2,
                    is_ref: false,
                },
                SeqEdge {
                    from: 1,
                    to: 3,
                    support: 2,
                    is_ref: true,
                },
                SeqEdge {
                    from: 2,
                    to: 3,
                    support: 2,
                    is_ref: false,
                },
            ],
            outgoing: HashMap::new(),
            incoming: HashMap::new(),
        };
        g.rebuild_index();
        // Prefix XX vs YY, no shared diamond prefix of length>=1 across first base; suffix TAA.
        // MergeDiamonds: common prefix empty, suffix TAA nonempty → diamonds fire first.
        // Either diamonds or suffix-split must keep src/sink and a branch.
        traced_simplify_graph_full(&mut g, |_stage, _graph| {});
        let after = inspect("after", &g);
        eprintln!("{}", after.line());
        eprintln!("  edges: {}", after.edge_dump);
        assert!(
            after.src_sink_ok(),
            "sink suffix split path: {}",
            after.line()
        );
        assert!(after.identity_ok());
    }

    #[test]
    fn first_simplify_matrix_unique_flank_synthetics() {
        let fixtures = vec![
            fixture("control_snp", "A", "C"),
            fixture("control_2bp_del", "TTCA", "TA"),
            fixture("control_2bp_ins", "TA", "TTCA"),
            fixture("control_two_snps", "AAAA", "ACCA"),
            fixture("control_del_plus_snp", "TTCA", "TAC"),
            fixture("holdout_a_TTCA_vs_TATG", "TTCA", "TATG"),
            fixture("holdout_b_GGAC_vs_GCTT", "GGAC", "GCTT"),
        ];
        let mut rows = Vec::new();
        for fx in fixtures {
            let reads = reads_for(&fx.ref_seq, &fx.alt_seq);
            let seq =
                build_seq(&fx.reference, &reads, 10).unwrap_or_else(|| panic!("{} graph", fx.name));
            let walk = report_fixture(fx.name, seq);
            let enter = walk
                .snaps
                .iter()
                .find(|s| s.stage == "entering_first_simplify")
                .expect("enter");
            rows.push(format!(
                "{}\tenter nv={} ne={} src={:?} sink={:?} ref={} alt={} cc={} both={}\tfirst_src_sink_loss={:?}\tfirst_id_loss={:?}\tfirst_nv_le_2={:?}\tfinal nv={} ne={} src={:?} sink={:?}",
                fx.name,
                enter.nv,
                enter.ne,
                enter.src,
                enter.sink,
                enter.ref_e,
                enter.alt_e,
                enter.components,
                enter.both,
                walk.first_src_sink_loss,
                walk.first_identity_loss,
                walk.first_nv_le_2,
                walk.snaps.last().map(|s| s.nv).unwrap_or(0),
                walk.snaps.last().map(|s| s.ne).unwrap_or(0),
                walk.snaps.last().and_then(|s| s.src),
                walk.snaps.last().and_then(|s| s.sink),
            ));
            assert!(
                walk.first_identity_loss.is_none(),
                "{}: identity must hold through first simplify",
                fx.name
            );
            assert!(
                enter.src_sink_ok() && enter.components == 1,
                "{}: enter simplify connected with src/sink: {}",
                fx.name,
                enter.line()
            );
            let after_zip = walk
                .snaps
                .iter()
                .find(|s| s.stage == "simplify_after_initial_zip")
                .expect("zip snap");
            assert_eq!(
                after_zip.components,
                1,
                "{}: zip must keep the graph connected (cc=1), got {}",
                fx.name,
                after_zip.line()
            );
            assert!(
                after_zip.src_sink_ok(),
                "{}: zip must keep src/sink, got {}",
                fx.name,
                after_zip.line()
            );
            assert!(
                walk.first_src_sink_loss.is_none(),
                "{}: src/sink must survive first simplify (diamond splice + common-suffix), loss={:?} final={}",
                fx.name,
                walk.first_src_sink_loss,
                walk.snaps.last().map(|s| s.line()).unwrap_or_default()
            );
            let last = walk.snaps.last().expect("last");
            assert!(
                last.src_sink_ok(),
                "{}: first simplify must end with src/sink: {}",
                fx.name,
                last.line()
            );
            if fx.name.starts_with("control_2bp_") {
                assert!(
                    last.nv >= 3 && (last.branches.len() + last.joins.len() > 0 || last.alt_e > 0),
                    "{}: indel bubble must not collapse to a single ref path: {}",
                    fx.name,
                    last.line()
                );
            }
        }

        let ref_only = fixture("reference_only", "A", "A");
        let ref_reads: Vec<_> = (0..8).map(|_| read(&ref_only.ref_seq, 30)).collect();
        let seq = build_seq(&ref_only.reference, &ref_reads, 10).expect("ref graph");
        let walk = report_fixture("reference_only", seq);
        rows.push(format!(
            "reference_only\tfirst_src_sink_loss={:?}\tfinal nv={} src={:?} sink={:?}",
            walk.first_src_sink_loss,
            walk.snaps.last().map(|s| s.nv).unwrap_or(0),
            walk.snaps.last().and_then(|s| s.src),
            walk.snaps.last().and_then(|s| s.sink),
        ));

        eprintln!("=== MATRIX ===");
        for r in &rows {
            eprintln!("{r}");
        }
    }

    #[test]
    fn p5_case1_first_simplify_keeps_src_sink() {
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
        let walk = report_fixture("p5_case1", seq);
        assert!(
            walk.first_src_sink_loss.is_none(),
            "p5_case1 must keep src/sink through first simplify; loss={:?}",
            walk.first_src_sink_loss
        );
        let last = walk.snaps.last().expect("last");
        assert!(last.src_sink_ok());
        assert!(last.identity_ok());
    }

    #[test]
    fn p12_production_waiver_not_touched_by_this_investigation() {
        let src = include_str!("assembly_based_caller.rs");
        let pin = "            assembler.use_seq_graph = false;\n            assembler.remove_paths_not_connected_to_ref = false;\n            assembler.skip_post_dangling_prune = true;";
        assert!(src.contains(pin));
    }

    #[test]
    fn unused_chain_like_compiles_for_local_probes() {
        let g = chain_like(3, &[(0, 1, true), (1, 2, true)]);
        assert_eq!(g.node_count(), 3);
    }
}
