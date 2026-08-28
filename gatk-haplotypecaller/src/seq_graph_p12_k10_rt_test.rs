//! 6R.12 TEST-ONLY: k=10/25 RT (non-SeqGraph) TATG path, cycles, and haplotype stages.
//! Does not change production assembler behavior, the P12 waiver, or W-H1.

#[cfg(test)]
mod traces {
    use crate::assembly::{AssemblyGraph, AssemblyRead};
    use crate::haplotype::Haplotype;
    use crate::kbest_haplotype::{
        find_best_haplotypes_for_assembly, find_best_haplotypes_preserving_cycles, graph_for_kbest,
    };
    use crate::read_event_discovery::P12_CLUSTER_TTC_START;
    use crate::read_threading_assembler::{
        assemble_from_ref_and_reads, build_threading_graph_for_haplotype_dump,
        build_threading_graph_for_seq_assembly, extract_haplotypes_from_kbest_paths,
        extract_rt_haplotypes_before_remove_paths, merge_rt_kbest_pre_remove_paths,
        merge_rt_kbest_pre_remove_paths_at_kmer, supplement_p12_cluster_coupled_haplotypes,
        AssemblyResult, AssemblyScoringContext, AssemblyStatus, ReadThreadingAssemblerArgs,
    };
    use std::collections::HashSet;
    use std::path::Path;

    const ALT_WIN: &[u8] = b"CTTTTATGTGATGTAT";
    const REF_WIN: &[u8] = b"CTTTTTCATGATGTAT";
    const ALT_CORE: &[u8] = b"TATGTG";
    const READ_WIN: &[u8] = b"CTTTTATGTGATGGAT";
    const REAL_P12_ACTIVE_START: u64 = P12_CLUSTER_TTC_START - 96;
    const REAL_P12_ACTIVE_END: u64 = P12_CLUSTER_TTC_START + 76;
    const REAL_P12_ATG_START: u64 = P12_CLUSTER_TTC_START + 3;

    fn p12_waiver_args() -> ReadThreadingAssemblerArgs {
        ReadThreadingAssemblerArgs {
            kmer_sizes: vec![10, 25],
            min_prune_factor: 2,
            allow_low_complexity_graphs: false,
            dont_increase_kmer_sizes_for_cycles: false,
            num_best_haplotypes_per_graph: 32,
            use_seq_graph: false,
            remove_paths_not_connected_to_ref: false,
            skip_post_dangling_prune: true,
            dangling_java_exact: true,
            scoring: Some(AssemblyScoringContext {
                padded_reference_start_1based: P12_CLUSTER_TTC_START - 696,
                active_start_1based: REAL_P12_ACTIVE_START,
                active_end_1based: REAL_P12_ACTIVE_END,
                contig: "2".into(),
            }),
            ..Default::default()
        }
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

    fn rt_cc(g: &AssemblyGraph) -> usize {
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

    fn contains_alt(seq: &[u8]) -> bool {
        seq.windows(ALT_WIN.len()).any(|w| w == ALT_WIN)
            || seq.windows(ALT_CORE.len()).any(|w| w == ALT_CORE)
            || seq.windows(b"CTTTTATGTG".len()).any(|w| w == b"CTTTTATGTG")
    }

    fn kmers(seq: &[u8], k: usize) -> HashSet<Vec<u8>> {
        if seq.len() < k {
            return HashSet::new();
        }
        (0..=seq.len() - k)
            .map(|i| seq[i..i + k].to_vec())
            .collect()
    }

    /// One directed cycle via DFS back-edge, or none.
    fn one_cycle(g: &AssemblyGraph) -> Option<Vec<usize>> {
        let n = g.node_count();
        let mut color = vec![0u8; n];
        let mut parent = vec![None; n];
        fn dfs(
            g: &AssemblyGraph,
            v: usize,
            color: &mut [u8],
            parent: &mut [Option<usize>],
        ) -> Option<Vec<usize>> {
            color[v] = 1;
            for t in g.outgoing_nodes(v) {
                if color[t] == 1 {
                    let mut cyc = vec![t];
                    let mut x = v;
                    cyc.push(x);
                    while x != t {
                        x = parent[x]?;
                        cyc.push(x);
                    }
                    cyc.reverse();
                    return Some(cyc);
                }
                if color[t] == 0 {
                    parent[t] = Some(v);
                    if let Some(c) = dfs(g, t, color, parent) {
                        return Some(c);
                    }
                }
            }
            color[v] = 2;
            None
        }
        for i in 0..n {
            if color[i] == 0 {
                if let Some(c) = dfs(g, i, &mut color, &mut parent) {
                    return Some(c);
                }
            }
        }
        None
    }

    fn dump_stage(stage: &str, haps: &[Haplotype], ref_bytes: &[u8]) {
        let uniq: HashSet<&[u8]> = haps.iter().map(|h| h.bases.as_slice()).collect();
        let tatg = haps.iter().any(|h| contains_alt(&h.bases));
        let alt_win = haps
            .iter()
            .any(|h| h.bases.windows(ALT_WIN.len()).any(|w| w == ALT_WIN));
        let lens: Vec<usize> = haps.iter().map(|h| h.bases.len()).collect();
        eprintln!(
            "STAGE {stage} n={} uniq={} tatg_any={tatg} alt_win={alt_win} ref_len={} hap_lens={lens:?}",
            haps.len(),
            uniq.len(),
            ref_bytes.len()
        );
        for (i, h) in haps.iter().enumerate() {
            if contains_alt(&h.bases) {
                let win = h
                    .bases
                    .windows(ALT_WIN.len())
                    .position(|w| w == ALT_WIN)
                    .map(|_| "ALT_WIN")
                    .unwrap_or("TATG/core");
                eprintln!(
                    "  hap{i} is_ref={} len={} {win} cigar={:?}",
                    h.is_reference,
                    h.bases.len(),
                    h.cigar.as_ref().map(|c| format!("{c:?}"))
                );
            }
        }
    }

    fn dump_graph(label: &str, g: &AssemblyGraph, ref_bytes: &[u8], k: usize) {
        let ref_kmers = kmers(ref_bytes, k);
        let mut alt_kmers = 0usize;
        let mut shared = 0usize;
        let mut non_ref_v = 0usize;
        for n in g.nodes() {
            let km = n.kmer.as_ref();
            let in_ref = ref_kmers.contains(km);
            if !in_ref {
                non_ref_v += 1;
            } else {
                shared += 1;
            }
            if km.windows(ALT_CORE.len()).any(|w| w == ALT_CORE)
                || km.windows(b"CTTTTATGTG".len()).any(|w| w == b"CTTTTATGTG")
            {
                alt_kmers += 1;
            }
        }
        let branches = (0..g.node_count())
            .filter(|&v| g.outgoing_nodes(v).len() > 1)
            .count();
        let joins = (0..g.node_count())
            .filter(|&v| g.incoming_count(v) > 1)
            .count();
        let src = g.reference_source_vertex();
        let sink = g.reference_sink_vertex();
        let dangling = (0..g.node_count())
            .filter(|&v| {
                (g.incoming_count(v) == 0 && Some(v) != src)
                    || (g.outgoing_nodes(v).is_empty() && Some(v) != sink)
            })
            .count();
        eprintln!(
            "{label} nv={} ne={} cc={} cycle={} src={src:?} sink={sink:?} ref_v≈{shared} non_ref_v={non_ref_v} dangling={dangling} branches={branches} joins={joins} alt_kmers={alt_kmers}",
            g.node_count(),
            g.edge_count(),
            rt_cc(g),
            g.has_cycle()
        );
        if let Some(cyc) = one_cycle(g) {
            let tatg_on = cyc
                .iter()
                .any(|&v| g.kmer_at(v).windows(ALT_CORE.len()).any(|w| w == ALT_CORE));
            eprintln!(
                "  cycle_len={} tatg_on_cycle={tatg_on} verts={:?}",
                cyc.len().saturating_sub(1),
                &cyc[..cyc.len().min(16)]
            );
            for &v in cyc.iter().take(8) {
                let km = g.kmer_at(v);
                eprintln!(
                    "    cyc v{v} in_ref={} k={}",
                    ref_kmers.contains(km),
                    String::from_utf8_lossy(km)
                );
            }
        }
        let interesting: Vec<usize> = (0..g.node_count())
            .filter(|&v| {
                let km = g.kmer_at(v);
                km.windows(ALT_WIN.len()).any(|w| w == ALT_WIN)
                    || km.windows(REF_WIN.len()).any(|w| w == REF_WIN)
                    || km.windows(ALT_CORE.len()).any(|w| w == ALT_CORE)
            })
            .collect();
        eprintln!("  cluster_vertices={}", interesting.len());
        for &v in interesting.iter().take(12) {
            let km = g.kmer_at(v);
            let ins: Vec<_> = g
                .incoming_nodes(v)
                .into_iter()
                .map(|p| format!("{p}/s{}", g.edge_support(p, v).unwrap_or(0)))
                .collect();
            let outs: Vec<_> = g
                .outgoing_nodes(v)
                .into_iter()
                .map(|t| {
                    format!(
                        "{t}/s{}/ref={}",
                        g.edge_support(v, t).unwrap_or(0),
                        g.edge_is_ref(v, t)
                    )
                })
                .collect();
            eprintln!(
                "    v{v} in_ref={} k={} in={ins:?} out={outs:?}",
                ref_kmers.contains(km),
                String::from_utf8_lossy(km)
            );
        }
    }

    fn src_sink_alt(g: &AssemblyGraph) -> bool {
        let Some(src) = g.reference_source_vertex() else {
            return false;
        };
        let Some(sink) = g.reference_sink_vertex() else {
            return false;
        };
        // BFS of paths is unbounded on cycles; use k-best preserve instead.
        let paths = find_best_haplotypes_preserving_cycles(g, 32).unwrap_or_default();
        let _ = (src, sink);
        paths.iter().any(|p| contains_alt(&p.bases(g)))
    }

    #[test]
    fn six_r12_production_pins_unchanged() {
        let abc = include_str!("assembly_based_caller.rs");
        assert!(abc.contains("assembler.use_seq_graph = false;"));
        assert!(abc.contains("assembler.dangling_java_exact = true;"));
        assert!(abc.contains("assembler.remove_paths_not_connected_to_ref = false;"));
        assert!(abc.contains("assembler.skip_post_dangling_prune = true;"));
        assert!(ReadThreadingAssemblerArgs::default().abort_seq_graph_on_cycles);
        assert!(ReadThreadingAssemblerArgs::default().use_seq_graph);
    }

    #[test]
    fn six_r12_real_p12_k10_k25_rt_mechanism() {
        let Some((reference, reads, ref_bytes)) = load_real_p12() else {
            eprintln!("Real-data P12 comparison unavailable");
            return;
        };
        let args = p12_waiver_args();
        eprintln!("=== 6R.12 REAL P12 k=10/25 RT mechanism ===");
        eprintln!("reads={} ref_len={}", reads.len(), ref_bytes.len());
        for k in [10usize, 25, 85] {
            eprintln!(
                "ref_non_unique_kmers k={k} {}",
                crate::read_threading_graph::reference_has_non_unique_kmers(&reference, k)
            );
        }

        let mut saw_tatg_after_strip = false;
        let mut saw_tatg_preserve = false;
        let mut saw_tatg_extract = false;

        for k in [10usize, 25] {
            let mut seq_args = args.clone();
            seq_args.use_seq_graph = true;
            let seq_graph = build_threading_graph_for_seq_assembly(
                &reference, &reads, k, &seq_args, true, true,
            )
            .expect("seq build");
            eprintln!(
                "SEQ_ASSEMBLY_GATE use_seq_graph=true k={k} graph={} cycle_on_rt={}",
                if seq_graph.is_some() { "Some" } else { "None" },
                build_threading_graph_for_haplotype_dump(&reference, &reads, k, &args, true, true,)
                    .ok()
                    .flatten()
                    .is_some_and(|g| g.has_cycle())
            );

            let Some(graph) =
                build_threading_graph_for_haplotype_dump(&reference, &reads, k, &args, true, true)
                    .expect("rt build")
            else {
                panic!("RT haplotype-dump graph missing at k={k}");
            };
            dump_graph(&format!("RT_DUMP k={k}"), &graph, &ref_bytes, k);

            let preserve_hit = src_sink_alt(&graph);
            saw_tatg_preserve |= preserve_hit;
            eprintln!("  kbest_preserve_src_sink_alt={preserve_hit}");

            let stripped = graph_for_kbest(graph.clone());
            match &stripped {
                Ok(dag) => {
                    eprintln!(
                        "  cycle_strip=OK nv {}→{} cycle={}",
                        graph.node_count(),
                        dag.node_count(),
                        dag.has_cycle()
                    );
                    let paths =
                        crate::kbest_haplotype::find_best_haplotypes(dag, 32).unwrap_or_default();
                    let hit = paths.iter().any(|p| contains_alt(&p.bases(dag)));
                    saw_tatg_after_strip |= hit;
                    eprintln!("  kbest_on_stripped n={} src_sink_alt={hit}", paths.len());
                }
                Err(_) => {
                    eprintln!("  cycle_strip=FAIL (Rust falls back to cyclic preserve)");
                }
            }

            let g_owned = graph;
            let (paths, g_used) =
                find_best_haplotypes_for_assembly(g_owned, 32).expect("kbest assembly");
            let kbest_hit = paths.iter().any(|p| contains_alt(&p.bases(&g_used)));
            eprintln!(
                "  find_best_haplotypes_for_assembly n={} src_sink_alt={kbest_hit}",
                paths.len()
            );
            let mut ref_hap = Haplotype::new(ref_bytes.as_slice(), true);
            let mut cig = crate::cigar::Cigar::new();
            cig.push(ref_hap.bases.len(), crate::cigar::CigarOperator::Match);
            ref_hap.cigar = Some(cig);
            let ref_cl = ref_hap.cigar.as_ref().unwrap().reference_length();
            let extracted = extract_haplotypes_from_kbest_paths(
                &paths,
                &g_used,
                &ref_hap,
                ref_cl,
                &args.haplotype_to_reference_sw,
            )
            .expect("extract");
            dump_stage(&format!("extract_kbest k={k}"), &extracted, &ref_bytes);
            saw_tatg_extract |= extracted.iter().any(|h| contains_alt(&h.bases));

            let before =
                extract_rt_haplotypes_before_remove_paths(&reference, &reads, &args, k, true, true)
                    .expect("before_remove");
            dump_stage(
                &format!("extract_rt_before_remove k={k}"),
                &before,
                &ref_bytes,
            );
        }

        let assembled = assemble_from_ref_and_reads(&reference, &reads, &args).expect("assemble");
        dump_stage(
            "assemble_from_ref_and_reads waiver",
            &assembled.haplotypes,
            &ref_bytes,
        );
        let final_alt_win = assembled
            .haplotypes
            .iter()
            .any(|h| h.bases.windows(ALT_WIN.len()).any(|w| w == ALT_WIN));
        let final_tatg = assembled.haplotypes.iter().any(|h| contains_alt(&h.bases));
        eprintln!(
            "VERDICT strip_tatg={saw_tatg_after_strip} preserve_tatg={saw_tatg_preserve} extract_tatg={saw_tatg_extract} final_tatg={final_tatg} final_alt_win={final_alt_win} status={:?} kmer_size={}",
            assembled.status,
            assembled.kmer_size
        );

        assert!(
            saw_tatg_preserve,
            "k=10/25 cyclic-preserve k-best must see TATGTG before SeqGraph"
        );
        assert!(
            saw_tatg_extract,
            "extract_haplotypes_from_kbest_paths must keep TATGTG"
        );
    }

    fn hap_flags(haps: &[Haplotype]) -> (bool, bool, bool, usize) {
        let core = haps.iter().any(|h| contains_alt(&h.bases));
        let win = haps
            .iter()
            .any(|h| h.bases.windows(ALT_WIN.len()).any(|w| w == ALT_WIN));
        let read_win = haps
            .iter()
            .any(|h| h.bases.windows(READ_WIN.len()).any(|w| w == READ_WIN));
        (core, win, read_win, haps.len())
    }

    fn dump_alt_kbest_paths(label: &str, g: &crate::assembly::AssemblyGraph) {
        let paths = find_best_haplotypes_preserving_cycles(g, 32).unwrap_or_default();
        let mut n_core = 0usize;
        for p in &paths {
            let b = p.bases(g);
            if !contains_alt(&b) {
                continue;
            }
            n_core += 1;
            if n_core > 3 {
                continue;
            }
            let mut verts = vec![p.start];
            for &(_, t) in &p.edges {
                verts.push(t);
            }
            eprintln!(
                "  {label} kbest_preserve alt_path#{n_core} n_verts={} score={:.4} is_ref={} alt_win={} read_win={} verts={:?}",
                verts.len(),
                p.score,
                p.is_reference,
                b.windows(ALT_WIN.len()).any(|w| w == ALT_WIN),
                b.windows(READ_WIN.len()).any(|w| w == READ_WIN),
                &verts[..verts.len().min(24)]
            );
        }
        eprintln!(
            "  {label} kbest_preserve n_paths={} n_tatg_core={n_core}",
            paths.len()
        );
        match graph_for_kbest(g.clone()) {
            Ok(dag) => {
                let dag_paths =
                    crate::kbest_haplotype::find_best_haplotypes(&dag, 32).unwrap_or_default();
                let n_dag_core = dag_paths
                    .iter()
                    .filter(|p| contains_alt(&p.bases(&dag)))
                    .count();
                eprintln!(
                    "  {label} kbest_after_strip n_paths={} n_tatg_core={n_dag_core} nv {}→{}",
                    dag_paths.len(),
                    g.node_count(),
                    dag.node_count()
                );
                for (i, p) in dag_paths
                    .iter()
                    .filter(|p| contains_alt(&p.bases(&dag)))
                    .take(2)
                    .enumerate()
                {
                    let mut verts = vec![p.start];
                    for &(_, t) in &p.edges {
                        verts.push(t);
                    }
                    let b = p.bases(&dag);
                    eprintln!(
                        "    strip_alt#{i} n_verts={} alt_win={} read_win={} verts={:?}",
                        verts.len(),
                        b.windows(ALT_WIN.len()).any(|w| w == ALT_WIN),
                        b.windows(READ_WIN.len()).any(|w| w == READ_WIN),
                        &verts[..verts.len().min(24)]
                    );
                }
            }
            Err(_) => eprintln!("  {label} cycle_strip=ERR (source→sink DAG not produced)"),
        }
    }

    #[test]
    fn six_r12_tatg_haplotype_provenance() {
        let Some((reference, reads, ref_bytes)) = load_real_p12() else {
            eprintln!("Real-data P12 comparison unavailable");
            return;
        };
        let args = p12_waiver_args();
        eprintln!("=== 6R.12 TATG provenance ===");
        eprintln!(
            "rt_first_skipped_on_p12=true use_seq_graph={} scoring_p12={}",
            args.use_seq_graph,
            args.scoring
                .as_ref()
                .is_some_and(|c| c.overlaps_p12_cluster())
        );

        let mut first_core: Option<(String, usize)> = None;
        let mut first_win: Option<(String, usize)> = None;
        let mut first_read_win: Option<(String, usize)> = None;
        let mut note = |label: &str, k: usize, haps: &[Haplotype]| {
            let (core, win, read_win, n) = hap_flags(haps);
            eprintln!(
                "PROV {label} k={k} n={n} tatg_core={core} alt_win={win} read_win={read_win}"
            );
            if core && first_core.is_none() {
                first_core = Some((label.to_string(), k));
            }
            if win && first_win.is_none() {
                first_win = Some((label.to_string(), k));
            }
            if read_win && first_read_win.is_none() {
                first_read_win = Some((label.to_string(), k));
            }
            for (i, h) in haps.iter().enumerate() {
                if h.bases.windows(ALT_WIN.len()).any(|w| w == ALT_WIN)
                    || h.bases.windows(READ_WIN.len()).any(|w| w == READ_WIN)
                    || contains_alt(&h.bases)
                {
                    eprintln!(
                        "    hap{i} is_ref={} len={} cigar={:?}",
                        h.is_reference,
                        h.bases.len(),
                        h.cigar.as_ref().map(|c| format!("{c:?}"))
                    );
                }
            }
        };

        for k in [10usize, 25, 35, 85] {
            match build_threading_graph_for_haplotype_dump(&reference, &reads, k, &args, true, true)
            {
                Ok(Some(g)) => {
                    dump_graph(&format!("PROV_GRAPH k={k}"), &g, &ref_bytes, k);
                    dump_alt_kbest_paths(&format!("k={k}"), &g);
                }
                Ok(None) => eprintln!("PROV_GRAPH k={k} None even with allow_nu=true"),
                Err(e) => eprintln!("PROV_GRAPH k={k} ERR {e}"),
            }
            match build_threading_graph_for_haplotype_dump(
                &reference, &reads, k, &args, false, false,
            ) {
                Ok(Some(_)) => eprintln!("PROV_GRAPH_configured_flags k={k} Some"),
                Ok(None) => eprintln!(
                    "PROV_GRAPH_configured_flags k={k} None (non-unique/low-complexity/cycle-abort)"
                ),
                Err(e) => eprintln!("PROV_GRAPH_configured_flags k={k} ERR {e}"),
            }
        }

        for k in [10usize, 25, 35, 45, 55, 65, 75, 85] {
            let configured = args.kmer_sizes.contains(&k);
            let last_expanded = k == 85;
            let try_lc = if configured { false } else { last_expanded };
            let try_nu = try_lc;
            let merge_lc = !configured;
            let merge_nu = !configured;
            let try_haps = extract_rt_haplotypes_before_remove_paths(
                &reference, &reads, &args, k, try_lc, try_nu,
            )
            .unwrap_or_else(|_| Vec::new());
            note(
                if configured {
                    "try_assemble_configured"
                } else {
                    "try_assemble_expanded"
                },
                k,
                &try_haps,
            );
            let merge_haps = extract_rt_haplotypes_before_remove_paths(
                &reference, &reads, &args, k, merge_lc, merge_nu,
            )
            .unwrap_or_else(|_| Vec::new());
            note(
                if configured {
                    "merge_rt_configured"
                } else {
                    "merge_rt_expanded"
                },
                k,
                &merge_haps,
            );
        }

        let mut merged: Vec<Haplotype> = Vec::new();
        for k in [10usize, 25, 35, 45, 55, 65, 75, 85] {
            merge_rt_kbest_pre_remove_paths_at_kmer(
                &reference,
                &reads,
                &args,
                &[],
                &mut merged,
                Some(k),
            )
            .expect("merge_rt_at_k");
            note("merge_rt_incremental", k, &merged);
        }

        let k85_only =
            extract_rt_haplotypes_before_remove_paths(&reference, &reads, &args, 85, true, true)
                .unwrap_or_else(|_| Vec::new());
        note("k85_kbest_only_before_merge", 85, &k85_only);
        let mut k85_then_merge = k85_only.clone();
        merge_rt_kbest_pre_remove_paths(&reference, &reads, &args, &[], &mut k85_then_merge)
            .expect("k85_merge");
        note("k85_kbest_plus_merge_rt", 85, &k85_then_merge);

        let mut after_merge = AssemblyResult {
            status: AssemblyStatus::AssembledSomeVariation,
            kmer_size: 85,
            haplotypes: merged.clone(),
            event_maps: Vec::new(),
        };
        supplement_p12_cluster_coupled_haplotypes(&mut after_merge, &reference, &reads, &args)
            .expect("supplement");
        note("after_supplement_on_merge", 0, &after_merge.haplotypes);

        let mut ref_only = AssemblyResult {
            status: AssemblyStatus::JustAssembledReference,
            kmer_size: 10,
            haplotypes: {
                let mut h = Haplotype::new(ref_bytes.as_slice(), true);
                let mut c = crate::cigar::Cigar::new();
                c.push(h.bases.len(), crate::cigar::CigarOperator::Match);
                h.cigar = Some(c);
                vec![h]
            },
            event_maps: Vec::new(),
        };
        let n_before = ref_only.haplotypes.len();
        supplement_p12_cluster_coupled_haplotypes(&mut ref_only, &reference, &reads, &args)
            .expect("supplement_ref");
        note("after_supplement_on_ref_only", 0, &ref_only.haplotypes);
        eprintln!(
            "  supplement_ref_only grew {}→{}",
            n_before,
            ref_only.haplotypes.len()
        );

        let assembled = assemble_from_ref_and_reads(&reference, &reads, &args).expect("assemble");
        note(
            "assemble_from_ref_and_reads",
            assembled.kmer_size,
            &assembled.haplotypes,
        );

        eprintln!(
            "FIRST_CORE={first_core:?} FIRST_ALT_WIN={first_win:?} FIRST_READ_WIN={first_read_win:?} assemble_k={}",
            assembled.kmer_size
        );
        assert!(
            first_core.is_some(),
            "some RT extract/merge/supplement stage must reconstruct TATGTG"
        );
    }

    /// Unique-flank bubble: RT k-best (non-SeqGraph) recovers TATG when k is small enough to share k-mers.
    #[test]
    fn six_r12_synthetic_small_k_rt_recovers_bubble() {
        let left = "ACGTACGGTTAGCCATAACGGTCCATTGCATAGCTGGAACCT";
        let right = "GCTTAGGAACCGGTTAACCGATCCTGAACCGGATCCATAGCT";
        let ref_seq = format!("{left}TTCA{right}");
        let alt_seq = format!("{left}TATG{right}");
        let reference = AssemblyRead {
            bases: ref_seq.as_bytes().to_vec(),
            base_quals: vec![30; ref_seq.len()],
        };
        let alt = AssemblyRead {
            bases: alt_seq.as_bytes().to_vec(),
            base_quals: vec![30; alt_seq.len()],
        };
        let g = crate::read_threading_graph::assembly_graph_from_ref_and_reads_threading(
            &reference,
            &[alt],
            &crate::assembly::AssemblyGraphParams {
                kmer_size: crate::bio_ids::KmerSize::try_from_usize(10).expect("k=10"),
                min_base_quality: 10,
                start_threading_only_at_existing_vertex: false,
                ..Default::default()
            },
        )
        .expect("thread");
        let alt_attached = g.nodes().iter().any(|n| {
            n.kmer.as_ref().windows(4).any(|w| w == b"TATG")
                && !ref_seq
                    .as_bytes()
                    .windows(n.kmer.len())
                    .any(|w| w == n.kmer.as_ref())
        });
        assert!(
            rt_cc(&g) == 1 && alt_attached,
            "synthetic k=10 RT graph must attach a TATG k-mer to REF cc={} alt={alt_attached}",
            rt_cc(&g)
        );
    }
}
