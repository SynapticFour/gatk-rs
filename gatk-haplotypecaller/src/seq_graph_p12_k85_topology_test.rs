//! 6R.9 TEST-ONLY: real P12 k=85 SeqGraph first-divergence topology trace.
//! Does not change production cleanup, waiver, or W-H1.

#[cfg(test)]
mod traces {
    use super::super::*;
    use crate::assembly::{AssemblyGraph, AssemblyRead};
    use crate::read_event_discovery::P12_CLUSTER_TTC_START;
    use crate::read_threading_assembler::{
        build_threading_graph_for_haplotype_dump, build_threading_graph_for_seq_assembly,
        extract_rt_haplotypes_before_remove_paths, AssemblyScoringContext,
        ReadThreadingAssemblerArgs,
    };
    use crate::seq_kbest_haplotype::find_best_haplotypes_seq_graph;
    use std::path::Path;

    const LEFT: &str = "ACGTACGGTTAGCCATAACGGTCCATTGCATAGCTGGAACCT";
    const RIGHT: &str = "GCTTAGGAACCGGTTAACCGATCCTGAACCGGATCCATAGCT";
    const REF_WIN: &[u8] = b"CTTTTTCATGATGTAT";
    const ALT_WIN: &[u8] = b"CTTTTATGTGATGTAT";
    const REF_CORE: &[u8] = b"TTCATG";
    const ALT_CORE: &[u8] = b"TATGTG";
    /// Real NA12878 20k P12 window, expressed as offsets from TTC start (N-1).
    const REAL_P12_ACTIVE_START: u64 = P12_CLUSTER_TTC_START - 96;
    const REAL_P12_ACTIVE_END: u64 = P12_CLUSTER_TTC_START + 76;
    const REAL_P12_PAD_START: u64 = P12_CLUSTER_TTC_START - 696;
    const REAL_P12_ATG_START: u64 = P12_CLUSTER_TTC_START + 3;

    fn read(seq: &str, q: u8) -> AssemblyRead {
        AssemblyRead {
            bases: seq.as_bytes().to_vec(),
            base_quals: vec![q; seq.len()],
        }
    }

    fn p12_mode_args(use_seq_graph: bool) -> ReadThreadingAssemblerArgs {
        ReadThreadingAssemblerArgs {
            kmer_sizes: vec![10, 25],
            min_prune_factor: 2,
            allow_low_complexity_graphs: true,
            dont_increase_kmer_sizes_for_cycles: true,
            num_best_haplotypes_per_graph: 32,
            use_seq_graph,
            remove_paths_not_connected_to_ref: false,
            skip_post_dangling_prune: true,
            dangling_java_exact: true,
            scoring: Some(AssemblyScoringContext {
                padded_reference_start_1based: P12_CLUSTER_TTC_START - LEFT.len() as u64,
                active_start_1based: P12_CLUSTER_TTC_START,
                active_end_1based: P12_CLUSTER_TTC_START.saturating_add(3),
                contig: "2".into(),
            }),
            ..Default::default()
        }
    }

    fn real_p12_args() -> ReadThreadingAssemblerArgs {
        ReadThreadingAssemblerArgs {
            use_seq_graph: true,
            remove_paths_not_connected_to_ref: false,
            skip_post_dangling_prune: true,
            dangling_java_exact: true,
            scoring: Some(AssemblyScoringContext {
                padded_reference_start_1based: REAL_P12_PAD_START,
                active_start_1based: REAL_P12_ACTIVE_START,
                active_end_1based: REAL_P12_ACTIVE_END,
                contig: "2".into(),
            }),
            ..Default::default()
        }
    }

    fn id_ok(g: &SeqGraph) -> bool {
        g.test_vertex_ids()
            .iter()
            .enumerate()
            .all(|(i, id)| *id == i)
            && g.edges_pub()
                .iter()
                .all(|e| e.from < g.node_count() && e.to < g.node_count())
    }

    fn undirected_components(g: &SeqGraph) -> usize {
        let n = g.node_count();
        if n == 0 {
            return 0;
        }
        let mut adj = vec![Vec::new(); n];
        for e in g.edges_pub() {
            adj[e.from].push(e.to);
            adj[e.to].push(e.from);
        }
        let mut seen = vec![false; n];
        let mut c = 0usize;
        for i in 0..n {
            if seen[i] {
                continue;
            }
            c += 1;
            let mut stack = vec![i];
            seen[i] = true;
            while let Some(v) = stack.pop() {
                for &w in &adj[v] {
                    if !seen[w] {
                        seen[w] = true;
                        stack.push(w);
                    }
                }
            }
        }
        c
    }

    fn ref_path_ok(g: &SeqGraph) -> bool {
        let (Some(src), Some(sink)) = (g.reference_source_vertex(), g.reference_sink_vertex())
        else {
            return false;
        };
        let mut cur = src;
        let mut guard = 0usize;
        while cur != sink {
            guard += 1;
            if guard > g.node_count() + 2 {
                return false;
            }
            let next = g
                .outgoing_nodes(cur)
                .into_iter()
                .find(|&t| g.edge_is_ref(cur, t));
            let Some(n) = next else {
                return false;
            };
            cur = n;
        }
        true
    }

    fn dangling_vertex_count(g: &SeqGraph) -> usize {
        let src = g.reference_source_vertex();
        let sink = g.reference_sink_vertex();
        (0..g.node_count())
            .filter(|&v| {
                let dangling_head = g.vertex_in_degree(v) == 0 && Some(v) != src;
                let dangling_tail = g.vertex_out_degree(v) == 0 && Some(v) != sink;
                dangling_head || dangling_tail
            })
            .count()
    }

    #[derive(Clone, Copy)]
    struct Topo {
        nv: usize,
        ne: usize,
        src: Option<usize>,
        sink: Option<usize>,
        id_ok: bool,
        dangling: usize,
        components: usize,
        ref_path: bool,
        non_ref: usize,
        branches: usize,
        joins: usize,
        alt_seq: bool,
        alt_kbest: bool,
        kbest_n: usize,
    }

    fn alt_in_seq_graph(g: &SeqGraph) -> bool {
        (0..g.node_count()).any(|v| {
            let s = g.vertex_sequence(v);
            s.windows(ALT_CORE.len()).any(|w| w == ALT_CORE)
                || s.windows(ALT_WIN.len()).any(|w| w == ALT_WIN)
        })
    }

    fn alt_in_kbest(g: &SeqGraph) -> (bool, usize) {
        let paths = find_best_haplotypes_seq_graph(g, 32).unwrap_or_default();
        let n = paths.len();
        let hit = paths.iter().any(|p| {
            let b = g.path_bases_bytes(p.start, &p.edges);
            b.windows(ALT_WIN.len()).any(|w| w == ALT_WIN)
                || b.windows(ALT_CORE.len()).any(|w| w == ALT_CORE)
        });
        (hit, n)
    }

    fn snap_topo(g: &SeqGraph) -> Topo {
        let (alt_kbest, kbest_n) = if g.node_count() <= 2500 && g.edge_count() <= 5000 {
            alt_in_kbest(g)
        } else {
            (false, 0)
        };
        Topo {
            nv: g.node_count(),
            ne: g.edge_count(),
            src: g.reference_source_vertex(),
            sink: g.reference_sink_vertex(),
            id_ok: id_ok(g),
            dangling: dangling_vertex_count(g),
            components: undirected_components(g),
            ref_path: ref_path_ok(g),
            non_ref: g.edges_pub().iter().filter(|e| !e.is_ref).count(),
            branches: (0..g.node_count())
                .filter(|&v| g.vertex_out_degree(v) > 1)
                .count(),
            joins: (0..g.node_count())
                .filter(|&v| g.vertex_in_degree(v) > 1)
                .count(),
            alt_seq: alt_in_seq_graph(g),
            alt_kbest,
            kbest_n,
        }
    }

    fn alt_present(t: &Topo) -> bool {
        t.alt_seq || t.alt_kbest || t.non_ref > 0
    }

    fn print_topo(label: &str, t: &Topo) {
        eprintln!(
            "{label}\tnv={}\tne={}\tsrc={:?}\tsink={:?}\tid_ok={}\tdangling={}\tcc={}\tref_path={}\tnon_ref={}\tbranches={}\tjoins={}\talt_seq={}\talt_kbest={}\tkbest={}",
            t.nv, t.ne, t.src, t.sink, t.id_ok, t.dangling, t.components, t.ref_path,
            t.non_ref, t.branches, t.joins, t.alt_seq, t.alt_kbest, t.kbest_n
        );
    }

    fn dump_compact_edges(g: &SeqGraph, cap: usize) {
        let n = g.edges_pub().len().min(cap);
        for (i, e) in g.edges_pub().iter().take(n).enumerate() {
            let fs = String::from_utf8_lossy(g.vertex_sequence(e.from));
            let ts = String::from_utf8_lossy(g.vertex_sequence(e.to));
            let fs = if fs.len() > 24 {
                format!("{}..len{}", &fs[..12], fs.len())
            } else {
                fs.into_owned()
            };
            let ts = if ts.len() > 24 {
                format!("{}..len{}", &ts[..12.min(ts.len())], ts.len())
            } else {
                ts.into_owned()
            };
            eprintln!(
                "  e{i} {}->{} ref={} sup={} from_seq={fs} to_seq={ts}",
                e.from, e.to, e.is_ref, e.support
            );
        }
        if g.edges_pub().len() > cap {
            eprintln!("  ... {} more edges", g.edges_pub().len() - cap);
        }
    }

    fn dump_non_ref_edges(g: &SeqGraph) {
        let non: Vec<_> = g.edges_pub().iter().filter(|e| !e.is_ref).collect();
        eprintln!("  non_ref_edges={}", non.len());
        for (i, e) in non.iter().take(32).enumerate() {
            eprintln!(
                "  nr{i} {}->{} sup={} from_len={} to_len={}",
                e.from,
                e.to,
                e.support,
                g.vertex_sequence(e.from).len(),
                g.vertex_sequence(e.to).len()
            );
        }
    }

    fn rt_alt_stats(graph: &AssemblyGraph, ref_seq: &[u8]) -> (usize, usize, usize, usize) {
        let alt_kmers = graph
            .nodes()
            .iter()
            .filter(|n| {
                let k = n.kmer.as_ref();
                k.windows(ALT_CORE.len()).any(|w| w == ALT_CORE)
                    || !ref_seq.windows(k.len()).any(|w| w == k)
            })
            .count();
        let non_ref = graph
            .edges_sorted()
            .iter()
            .filter(|e| !graph.edge_is_ref(e.from, e.to))
            .count();
        (graph.node_count(), graph.edge_count(), non_ref, alt_kmers)
    }

    fn load_real_p12() -> Option<(AssemblyRead, Vec<AssemblyRead>, Vec<u8>)> {
        use crate::assembly_region_finalize::{
            assembly_reference_read, finalize_region_reads_for_assembly,
            gatk_min_tail_quality_for_assembly, records_to_assembly_reads,
        };
        use crate::read_model::ReadFilterParams;
        use crate::walker_apply::{call_disposition, AssemblyRegionCallDisposition};
        use crate::walker_traversal::{
            flatten_assembly_regions, traverse_assembly_region_walker, WalkerTraversalConfig,
        };
        use gatk_core::reference::{
            parse_intervals_cli_string, ReferenceWindowCache, SequenceDictionary,
        };

        let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("..");
        let ref_path = root.join("parity/realworld/assets/hs37d5.simple.fa");
        let bam = root.join("parity/realworld/na12878_20k_b37/NA12878_20k.b37.bam");
        if !ref_path.is_file() || !bam.is_file() {
            return None;
        }
        let dict = SequenceDictionary::from_fasta_path(&ref_path).ok()?;
        let interval = format!("2:{REAL_P12_ACTIVE_START}-{REAL_P12_ACTIVE_END}");
        let specs = parse_intervals_cli_string(&dict, &interval).ok()?;
        let walk = traverse_assembly_region_walker(
            &dict,
            &specs,
            &ref_path,
            &bam,
            &ReadFilterParams::gatk_standard_hc(),
            &WalkerTraversalConfig::gatk_haplotype_caller_production(100),
        )
        .ok()?;
        let regions = flatten_assembly_regions(&walk);
        let region = regions.iter().find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= P12_CLUSTER_TTC_START
                && r.end.get() >= REAL_P12_ATG_START
        })?;
        let mut ref_cache = ReferenceWindowCache::new(ref_path.clone(), 4);
        let reference = assembly_reference_read(&dict, &mut ref_cache, region).ok()?;
        let finalized = finalize_region_reads_for_assembly(
            &region.reads,
            region,
            true,
            gatk_min_tail_quality_for_assembly(10),
            false,
        );
        let reads = records_to_assembly_reads(&finalized);
        let bases = reference.bases.clone();
        Some((reference, reads, bases))
    }

    fn trace_seqgraph(label: &str, mut seq: SeqGraph, ref_seq: &[u8]) {
        eprintln!("=== {label} SeqGraph after from_assembly_graph ===");
        let mut first_loss: Option<(String, Topo, String, Topo)> = None;
        let mut prev: Option<(String, Topo)> = None;
        let mut record = |stage: &str, g: &SeqGraph| {
            let t = snap_topo(g);
            print_topo(stage, &t);
            if t.nv <= 24 {
                dump_compact_edges(g, 32);
            } else if t.non_ref > 0 && t.non_ref <= 64 {
                dump_non_ref_edges(g);
            }
            if let Some((ps, pt)) = prev.as_ref() {
                if alt_present(pt) && !alt_present(&t) && first_loss.is_none() {
                    first_loss = Some((ps.clone(), *pt, stage.to_string(), t));
                }
            }
            prev = Some((stage.to_string(), t));
        };

        record("from_assembly_graph", &seq);
        seq.clean_non_ref_paths();
        record("after_clean_non_ref_paths", &seq);
        let status = seq.traced_cleanup_seq_graph(|stage, g| record(stage, g));
        eprintln!("cleanup_status={status:?}");

        if let Some((ps, pt, ns, nt)) = first_loss {
            eprintln!("=== FIRST ALT LOSS ===");
            eprintln!("previous={ps}");
            print_topo("prev", &pt);
            eprintln!("next={ns}");
            print_topo("next", &nt);
        } else if let Some((s, t)) = prev {
            eprintln!(
                "=== NO CLEANUP ALT LOSS (final stage={s} alt={}) ===",
                alt_present(&t)
            );
        }
        let _ = ref_seq;
    }

    #[test]
    fn six_r9_real_p12_k85_topology_trace() {
        let Some((reference, reads, ref_bytes)) = load_real_p12() else {
            eprintln!("Real-data P12 comparison unavailable");
            return;
        };
        let args = real_p12_args();
        eprintln!(
            "=== 6R.9 REAL P12 k=85 topology  ref_len={} nreads={} REF_WIN={} ALT_WIN={} ===",
            ref_bytes.len(),
            reads.len(),
            String::from_utf8_lossy(REF_WIN),
            String::from_utf8_lossy(ALT_WIN)
        );
        eprintln!(
            "REF contains REF_WIN={} ALT_WIN={}",
            ref_bytes.windows(REF_WIN.len()).any(|w| w == REF_WIN),
            ref_bytes.windows(ALT_WIN.len()).any(|w| w == ALT_WIN)
        );

        for kmer in [10usize, 25, 85] {
            let seq_g =
                build_threading_graph_for_seq_assembly(&reference, &reads, kmer, &args, true, true);
            let dump_g = build_threading_graph_for_haplotype_dump(
                &reference, &reads, kmer, &args, true, true,
            );
            match (seq_g, dump_g) {
                (Ok(Some(sg)), Ok(Some(dg))) => {
                    let (sn, se, snr, sa) = rt_alt_stats(&sg, &ref_bytes);
                    let (dn, de, dnr, da) = rt_alt_stats(&dg, &ref_bytes);
                    eprintln!(
                        "RT k={kmer} seq_assembly nv={sn} ne={se} non_ref={snr} alt_kmers={sa} cycle={} src={:?} sink={:?}",
                        sg.has_cycle(),
                        sg.reference_source_vertex(),
                        sg.reference_sink_vertex()
                    );
                    eprintln!(
                        "RT k={kmer} haplotype_dump nv={dn} ne={de} non_ref={dnr} alt_kmers={da} cycle={}",
                        dg.has_cycle()
                    );
                }
                (Ok(None), Ok(Some(dg))) => {
                    let (dn, de, dnr, da) = rt_alt_stats(&dg, &ref_bytes);
                    eprintln!(
                        "RT k={kmer} seq_assembly=None (cycle abort likely); haplotype_dump nv={dn} ne={de} non_ref={dnr} alt_kmers={da} cycle={}",
                        dg.has_cycle()
                    );
                }
                (Ok(None), Ok(None)) => eprintln!("RT k={kmer} both graphs None"),
                (a, b) => eprintln!(
                    "RT k={kmer} unexpected seq_ok={} dump_ok={}",
                    a.is_ok(),
                    b.is_ok()
                ),
            }
        }

        let rt_batch =
            extract_rt_haplotypes_before_remove_paths(&reference, &reads, &args, 85, true, true)
                .unwrap_or_default();
        let rt_alt = rt_batch.iter().any(|h| {
            h.bases.windows(ALT_WIN.len()).any(|w| w == ALT_WIN)
                || h.bases.windows(ALT_CORE.len()).any(|w| w == ALT_CORE)
        });
        eprintln!(
            "RT k=85 extract_before_remove_paths haps={} alt_window={rt_alt}",
            rt_batch.len()
        );

        let Some(graph) =
            build_threading_graph_for_seq_assembly(&reference, &reads, 85, &args, true, true)
                .expect("seq graph build")
        else {
            panic!("expected k=85 seq_assembly graph");
        };
        let (nv, ne, non_ref, alt_kmers) = rt_alt_stats(&graph, &ref_bytes);
        eprintln!(
            "=== RT k=85 BEFORE SeqGraph conversion nv={nv} ne={ne} non_ref={non_ref} alt_kmers={alt_kmers} cycle={} ===",
            graph.has_cycle()
        );
        if nv <= 40 {
            for (i, n) in graph.nodes().iter().enumerate() {
                let k = n.kmer.as_ref();
                let in_ref = ref_bytes.windows(k.len()).any(|w| w == k);
                let alt = k.windows(ALT_CORE.len()).any(|w| w == ALT_CORE);
                if !in_ref || alt || i < 4 || i + 4 >= nv {
                    let shown = if k.len() > 20 {
                        format!("{}..len{}", String::from_utf8_lossy(&k[..12]), k.len())
                    } else {
                        String::from_utf8_lossy(k).into_owned()
                    };
                    eprintln!("  rt_v{i} in_ref={in_ref} alt_core={alt} k={shown}");
                }
            }
            for e in graph.edges_sorted() {
                eprintln!(
                    "  rt_e {}->{} ref={} sup={}",
                    e.from,
                    e.to,
                    graph.edge_is_ref(e.from, e.to),
                    e.support
                );
            }
        } else {
            let mut shown = 0usize;
            for (i, n) in graph.nodes().iter().enumerate() {
                let k = n.kmer.as_ref();
                let in_ref = ref_bytes.windows(k.len()).any(|w| w == k);
                if !in_ref {
                    let shown_k = if k.len() > 24 {
                        format!("{}..len{}", String::from_utf8_lossy(&k[..16]), k.len())
                    } else {
                        String::from_utf8_lossy(k).into_owned()
                    };
                    eprintln!("  alt_kmer v{i} {shown_k}");
                    shown += 1;
                    if shown >= 12 {
                        break;
                    }
                }
            }
        }

        let seq = SeqGraph::from_assembly_graph(&graph);
        trace_seqgraph("REAL_P12_k85", seq, &ref_bytes);
    }

    #[test]
    fn six_r9_synthetic_ttca_tatg_k10_control() {
        let ref_seq = format!("{LEFT}TTCA{RIGHT}");
        let alt_seq = format!("{LEFT}TATG{RIGHT}");
        let reference = read(&ref_seq, 30);
        let mut reads = Vec::new();
        for _ in 0..4 {
            reads.push(read(&ref_seq, 30));
            reads.push(read(&alt_seq, 30));
        }
        let args = p12_mode_args(true);
        let graph =
            build_threading_graph_for_seq_assembly(&reference, &reads, 10, &args, true, true)
                .expect("synth build")
                .expect("synth graph");
        let (nv, ne, non_ref, alt_kmers) = rt_alt_stats(&graph, ref_seq.as_bytes());
        eprintln!(
            "=== 6R.9 SYNTHETIC k=10 RT nv={nv} ne={ne} non_ref={non_ref} alt_kmers={alt_kmers} cycle={} ===",
            graph.has_cycle()
        );
        for e in graph.edges_sorted() {
            eprintln!(
                "  syn_rt_e {}->{} ref={} sup={}",
                e.from,
                e.to,
                graph.edge_is_ref(e.from, e.to),
                e.support
            );
        }
        let seq = SeqGraph::from_assembly_graph(&graph);
        trace_seqgraph("SYNTHETIC_k10", seq, ref_seq.as_bytes());
    }

    fn rt_undirected_cc(g: &AssemblyGraph) -> usize {
        let n = g.node_count();
        if n == 0 {
            return 0;
        }
        let mut adj = vec![Vec::new(); n];
        for e in g.edges_sorted() {
            adj[e.from].push(e.to);
            adj[e.to].push(e.from);
        }
        let mut seen = vec![false; n];
        let mut c = 0usize;
        for i in 0..n {
            if seen[i] {
                continue;
            }
            c += 1;
            let mut stack = vec![i];
            seen[i] = true;
            while let Some(v) = stack.pop() {
                for &w in &adj[v] {
                    if !seen[w] {
                        seen[w] = true;
                        stack.push(w);
                    }
                }
            }
        }
        c
    }

    fn rt_alt_present(g: &AssemblyGraph, ref_seq: &[u8]) -> (usize, bool) {
        let mut n = 0usize;
        let mut tatg = false;
        for node in g.nodes() {
            let k = node.kmer.as_ref();
            let not_ref = !ref_seq.windows(k.len()).any(|w| w == k);
            let core = k.windows(ALT_CORE.len()).any(|w| w == ALT_CORE)
                || k.windows(ALT_WIN.len()).any(|w| w == ALT_WIN);
            if not_ref || core {
                n += 1;
            }
            tatg |= core;
        }
        (n, tatg)
    }

    fn rt_heads_tails(g: &AssemblyGraph) -> (usize, usize) {
        let heads = (0..g.node_count())
            .filter(|&v| {
                g.incoming_count(v) == 0
                    && !g.outgoing_nodes(v).is_empty()
                    && !g.is_ref_source_vertex(v)
            })
            .count();
        let tails = (0..g.node_count())
            .filter(|&v| g.outgoing_nodes(v).is_empty() && !g.is_ref_sink_vertex(v))
            .count();
        (heads, tails)
    }

    fn dump_rt(label: &str, g: &AssemblyGraph, ref_seq: &[u8]) {
        let (alt_n, tatg) = rt_alt_present(g, ref_seq);
        let (heads, tails) = rt_heads_tails(g);
        let non_ref = g
            .edges_sorted()
            .iter()
            .filter(|e| !g.edge_is_ref(e.from, e.to))
            .count();
        eprintln!(
            "{label}\tnv={}\tne={}\tcc={}\tsrc={:?}\tsink={:?}\tnon_ref_e={non_ref}\talt_kmers={alt_n}\ttatg_in_kmer={tatg}\theads={heads}\ttails={tails}",
            g.node_count(),
            g.edge_count(),
            rt_undirected_cc(g),
            g.reference_source_vertex(),
            g.reference_sink_vertex()
        );
    }

    fn dump_island_heads(g: &AssemblyGraph, ref_seq: &[u8]) {
        for v in 0..g.node_count() {
            if g.incoming_count(v) > 0
                || g.is_ref_source_vertex(v)
                || g.outgoing_nodes(v).is_empty()
            {
                continue;
            }
            let k = g.kmer_at(v);
            let in_ref = ref_seq.windows(k.len()).any(|w| w == k);
            let tatg = k.windows(ALT_CORE.len()).any(|w| w == ALT_CORE);
            if in_ref && !tatg {
                continue;
            }
            let shown = if k.len() > 24 {
                format!("{}..len{}", String::from_utf8_lossy(&k[..16]), k.len())
            } else {
                String::from_utf8_lossy(k).into_owned()
            };
            let outs: Vec<_> = g
                .outgoing_nodes(v)
                .into_iter()
                .map(|t| {
                    format!(
                        "{}(ref={},sup={})",
                        t,
                        g.edge_is_ref(v, t),
                        g.edge_support(v, t).unwrap_or(0)
                    )
                })
                .collect();
            eprintln!("  island_head v{v} in_ref={in_ref} tatg={tatg} k={shown} outs={outs:?}");
        }
    }

    fn seqgraph_zip_before_undirected(g: &AssemblyGraph) -> SeqGraph {
        let mut seq = SeqGraph::from_assembly_graph(g);
        seq.clean_non_ref_paths();
        seq.zip_linear_chains();
        seq.remove_singleton_orphan_vertices();
        seq
    }

    fn seq_src_sink_tatg(seq: &SeqGraph) -> bool {
        let paths = find_best_haplotypes_seq_graph(seq, 32).unwrap_or_default();
        paths.iter().any(|p| {
            let b = seq.path_bases_bytes(p.start, &p.edges);
            b.windows(ALT_WIN.len()).any(|w| w == ALT_WIN)
                || b.windows(ALT_CORE.len()).any(|w| w == ALT_CORE)
        })
    }

    fn dump_seq_pre_undirected(label: &str, seq: &SeqGraph) {
        let t = snap_topo(seq);
        eprintln!(
            "{label}\tnv={}\tne={}\tcc={}\tsrc={:?}\tsink={:?}\tnon_ref={}\talt_seq={}\tsrc_sink_tatg={}\tkbest={}",
            t.nv,
            t.ne,
            t.components,
            t.src,
            t.sink,
            t.non_ref,
            t.alt_seq,
            seq_src_sink_tatg(seq),
            t.kbest_n
        );
        if t.nv <= 16 {
            dump_compact_edges(seq, 24);
        }
    }

    #[test]
    fn six_r10_production_dangling_java_exact_pin_unchanged() {
        let abc = include_str!("assembly_based_caller.rs");
        assert!(
            abc.contains("assembler.dangling_java_exact = true;"),
            "strict_java still forces dangling_java_exact"
        );
        assert!(
            abc.contains("assembler.use_seq_graph = false;"),
            "P12 SeqGraph waiver must remain"
        );
        assert!(ReadThreadingAssemblerArgs::default().dangling_java_exact == false);
    }

    #[test]
    fn six_r10_real_p12_k85_dangling_java_exact_ab() {
        use crate::assembly_dangling_recovery::DanglingRecoveryParams;

        let Some((reference, reads, ref_bytes)) = load_real_p12() else {
            eprintln!("Real-data P12 comparison unavailable");
            return;
        };

        eprintln!("=== 6R.10 REAL P12 k=85 dangling_java_exact A/B ===");
        eprintln!(
            "note: recover=true threading uses start_threading_only_at_existing_vertex=false"
        );
        let live_args = real_p12_args();
        let thread_params = crate::assembly::AssemblyGraphParams {
            kmer_size: crate::bio_ids::KmerSize::try_from_usize(85).expect("k=85"),
            min_base_quality: live_args.min_base_quality,
            min_edge_weight: 1,
            dangling_path_max_nodes: 0,
            max_haplotypes: live_args.num_best_haplotypes_per_graph,
            max_haplotype_bases: 4096,
            start_threading_only_at_existing_vertex: false,
        };
        let (mut unpruned, _) =
            crate::read_threading_graph::assembly_graph_from_ref_and_reads_threading_with_summary(
                &reference,
                &reads,
                &thread_params,
            )
            .expect("thread k=85");
        eprintln!("BEFORE_PRUNE cycle={}", unpruned.has_cycle());
        dump_rt("BEFORE_PRUNE recover-true-threading", &unpruned, &ref_bytes);

        let mut pruning =
            crate::assembly::AssemblyGraphPruningParams::gatk_haplotype_caller_defaults();
        pruning.min_prune_factor = live_args.min_prune_factor;
        pruning.use_adaptive_pruning = live_args.use_adaptive_pruning;
        unpruned.apply_pruning(&pruning);
        let pruned = unpruned;
        dump_rt("AFTER_PRUNE_BEFORE_DANGLING", &pruned, &ref_bytes);
        dump_island_heads(&pruned, &ref_bytes);
        eprintln!("AFTER_PRUNE cycle={}", pruned.has_cycle());

        match build_threading_graph_for_seq_assembly(&reference, &reads, 85, &live_args, true, true)
        {
            Ok(Some(live)) => {
                dump_rt("SEQ_ASSEMBLY recover=true (6R.9 path)", &live, &ref_bytes);
            }
            Ok(None) => eprintln!("SEQ_ASSEMBLY recover=true graph=None"),
            Err(e) => eprintln!("SEQ_ASSEMBLY recover=true err={e}"),
        }

        for exact in [true, false] {
            let mut params = DanglingRecoveryParams::from_assembler_args(&live_args);
            params.dangling_java_exact = exact;
            params.recover_dangling_heads = true;
            let probe = pruned.probe_dangling_head_failures(&params);
            eprintln!(
                "=== HEAD PROBE dangling_java_exact={exact} candidates={} ===",
                probe.len()
            );
            let mut reason_counts: std::collections::BTreeMap<&str, usize> =
                std::collections::BTreeMap::new();
            for (v, kmer, reason) in &probe {
                *reason_counts.entry(reason.as_str()).or_insert(0) += 1;
                let tatg = kmer
                    .as_bytes()
                    .windows(ALT_CORE.len())
                    .any(|w| w == ALT_CORE);
                let in_ref = ref_bytes.windows(kmer.len()).any(|w| w == kmer.as_bytes());
                if !tatg && in_ref {
                    continue;
                }
                let shown = if kmer.len() > 24 {
                    format!("{}..len{}", &kmer[..16], kmer.len())
                } else {
                    kmer.clone()
                };
                eprintln!("  head v{v} tatg={tatg} in_ref={in_ref} reason={reason} k={shown}");
            }
            eprintln!("  head_reason_counts={reason_counts:?}");
            let tail_probe = pruned.probe_dangling_tail_failures(&params);
            eprintln!("  tail_candidates={}", tail_probe.len());
            for (v, kmer, reason) in tail_probe.iter().take(8) {
                let tatg = kmer
                    .as_bytes()
                    .windows(ALT_CORE.len())
                    .any(|w| w == ALT_CORE);
                let shown = if kmer.len() > 20 {
                    format!("{}..len{}", &kmer[..16.min(kmer.len())], kmer.len())
                } else {
                    kmer.clone()
                };
                eprintln!("  tail v{v} tatg={tatg} reason={reason} k={shown}");
            }
        }

        for exact in [true, false] {
            let mut g = pruned.clone();
            let mut params = DanglingRecoveryParams::from_assembler_args(&live_args);
            params.dangling_java_exact = exact;
            params.recover_dangling_heads = true;
            let summary = g.recover_dangling_branches(&params).expect("recover");
            eprintln!(
                "=== AFTER_DANGLING exact={exact} tails={}/{} heads={}/{} edges {}→{} ===",
                summary.tails_recovered,
                summary.tails_attempted,
                summary.heads_recovered,
                summary.heads_attempted,
                summary.edges_before,
                summary.edges_after
            );
            dump_rt(&format!("AFTER_DANGLING exact={exact}"), &g, &ref_bytes);
            dump_island_heads(&g, &ref_bytes);

            let seq = seqgraph_zip_before_undirected(&g);
            dump_seq_pre_undirected(
                &format!("AFTER_ZIP_BEFORE_UNDIRECTED_PRUNE exact={exact}"),
                &seq,
            );
            let attached = rt_undirected_cc(&g) == 1;
            let src_sink = seq_src_sink_tatg(&seq);
            eprintln!(
                "VERDICT_SLICE exact={exact} rt_cc1={attached} seq_src_sink_tatg={src_sink} seq_cc={}",
                undirected_components(&seq)
            );
        }
    }
}
