//! 6R.30 TEST-ONLY: SeqGraph / k-best topology at 92317361 C/T and 92317371 C/G.
//! Does not change production assembly, simplify, k-best, EventMap, or haplotype suppression.

#[cfg(test)]
mod traces {
    use crate::assembly::{AssemblyGraph, AssemblyGraphParams, AssemblyGraphPruningParams};
    use crate::assembly_based_caller::{assemble_reads_with_finalized, AssembleReadsArgs};
    use crate::assembly_region_finalize::{
        assembly_reference_read, create_graph_reference_read, finalize_region_reads_for_assembly,
        gatk_min_tail_quality_for_assembly, records_to_assembly_reads,
    };
    use crate::bio_ids::KmerSize;
    use crate::haplotype::Haplotype;
    use crate::read_projection::query_index_at_reference_position;
    use crate::read_threading_assembler::{
        assemble_from_ref_and_reads, build_threading_graph_for_seq_assembly,
        extract_rt_haplotypes_before_remove_paths, AssemblyScoringContext,
        DEFAULT_NUM_BEST_HAPLOTYPES_PER_GRAPH,
    };
    use crate::read_threading_graph::assembly_graph_from_ref_and_reads_threading_with_summary;
    use crate::seq_graph::SeqGraph;
    use crate::seq_kbest_haplotype::find_best_haplotypes_seq_graph;
    use crate::walker_apply::{call_disposition, AssemblyRegionCallDisposition};
    use gatk_core::reference::{
        parse_intervals_cli_string, ReferenceWindowCache, SequenceDictionary,
    };
    use rust_htslib::bam::record::CigarString;
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
    const MOTIF_361_REF: &[u8] = b"AGTCAC";
    const MOTIF_361_ALT: &[u8] = b"AGTCAT";
    const JAVA_ALT_MID: &[u8] = b"AACTCTTTCTGTAC";
    const JAVA_REF_MID: &[u8] = b"CACTCTTTTTGTAG";

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

    fn snp_kmers(ref_bases: &[u8], off: usize, k: usize, alt: u8) -> (Vec<Vec<u8>>, Vec<Vec<u8>>) {
        let mut refs = Vec::new();
        let mut alts = Vec::new();
        for start in overlapping_kmer_starts(ref_bases.len(), off, k) {
            let mut rk = ref_bases[start..start + k].to_vec();
            refs.push(rk.clone());
            rk[off - start] = alt;
            alts.push(rk);
        }
        (refs, alts)
    }

    fn combined_kmers(
        ref_bases: &[u8],
        k: usize,
        muts: &[(usize, u8)],
    ) -> (Vec<Vec<u8>>, Vec<Vec<u8>>) {
        let offs: Vec<usize> = muts.iter().map(|(o, _)| *o).collect();
        let mut refs = Vec::new();
        let mut alts = Vec::new();
        for start in windows_overlapping_all(ref_bases.len(), &offs, k) {
            let rk = ref_bases[start..start + k].to_vec();
            let mut ak = rk.clone();
            for &(off, alt) in muts {
                ak[off - start] = alt;
            }
            refs.push(rk);
            alts.push(ak);
        }
        (refs, alts)
    }

    fn ascii(b: &[u8]) -> String {
        String::from_utf8_lossy(b).into_owned()
    }

    fn kmer_row(graph: &AssemblyGraph, kmer: &[u8]) -> String {
        match graph.vertex_id_for_kmer(kmer) {
            None => "absent".into(),
            Some(id) => {
                let n = &graph.nodes()[id];
                format!("id={id} support={}", n.support)
            }
        }
    }

    fn present_count(graph: &AssemblyGraph, kmers: &[Vec<u8>]) -> usize {
        kmers
            .iter()
            .filter(|k| graph.vertex_id_for_kmer(k).is_some())
            .count()
    }

    fn motif_node_count(graph: &AssemblyGraph, needle: &[u8]) -> usize {
        graph
            .nodes()
            .iter()
            .filter(|n| n.kmer.windows(needle.len()).any(|w| w == needle))
            .count()
    }

    fn seq_motif_vertex_count(graph: &SeqGraph, needle: &[u8]) -> usize {
        (0..graph.node_count())
            .filter(|&v| {
                graph
                    .vertex_sequence(v)
                    .windows(needle.len())
                    .any(|w| w == needle)
            })
            .count()
    }

    fn dump_rt_kmers(
        label: &str,
        graph: &AssemblyGraph,
        ref_bases: &[u8],
        off_ct: usize,
        off_cg: usize,
        off_ca: usize,
        off_tc: usize,
        off_gc: usize,
    ) {
        eprintln!(
            "RT[{label}] nodes={} edges={} max_out={} motif361C={} motif361T={} zip_alt_11={} java_alt_mid={} java_ref_mid={}",
            graph.node_count(),
            graph.edge_count(),
            graph.max_out_degree(),
            motif_node_count(graph, MOTIF_361_REF),
            motif_node_count(graph, MOTIF_361_ALT),
            motif_node_count(graph, b"TAGAGTTGAAG"),
            motif_node_count(graph, JAVA_ALT_MID),
            motif_node_count(graph, JAVA_REF_MID),
        );
        for (name, refs, alts) in [
            (
                "361C/T_single_flip",
                snp_kmers(ref_bases, off_ct, K, b'T').0,
                snp_kmers(ref_bases, off_ct, K, b'T').1,
            ),
            (
                "371C/G_single_flip",
                snp_kmers(ref_bases, off_cg, K, b'G').0,
                snp_kmers(ref_bases, off_cg, K, b'G').1,
            ),
            (
                "361T+371G_combined",
                combined_kmers(ref_bases, K, &[(off_ct, b'T'), (off_cg, b'G')]).0,
                combined_kmers(ref_bases, K, &[(off_ct, b'T'), (off_cg, b'G')]).1,
            ),
            (
                "399C/A_single_flip",
                snp_kmers(ref_bases, off_ca, K, b'A').0,
                snp_kmers(ref_bases, off_ca, K, b'A').1,
            ),
            (
                "399A+407C+412C_combined",
                combined_kmers(
                    ref_bases,
                    K,
                    &[(off_ca, b'A'), (off_tc, b'C'), (off_gc, b'C')],
                )
                .0,
                combined_kmers(
                    ref_bases,
                    K,
                    &[(off_ca, b'A'), (off_tc, b'C'), (off_gc, b'C')],
                )
                .1,
            ),
        ] {
            let n_ref = present_count(graph, &refs);
            let n_alt = present_count(graph, &alts);
            eprintln!(
                "  KMER {name} windows={} ref_present={n_ref} alt_present={n_alt}",
                refs.len()
            );
            if let (Some(rk), Some(ak)) = (refs.first(), alts.first()) {
                eprintln!("    first_ref {} {}", ascii(rk), kmer_row(graph, rk));
                eprintln!("    first_alt {} {}", ascii(ak), kmer_row(graph, ak));
            }
            if name == "361T+371G_combined" {
                for (i, ak) in alts.iter().enumerate() {
                    if graph.vertex_id_for_kmer(ak).is_some() {
                        eprintln!(
                            "    ALT_HIT[{i}] {} agtcat={}",
                            ascii(ak),
                            ak.windows(MOTIF_361_ALT.len()).any(|w| w == MOTIF_361_ALT)
                        );
                    }
                }
            }
        }
    }

    fn count_st_paths(graph: &SeqGraph) -> (usize, bool) {
        let Some(src) = graph.reference_source_vertex() else {
            return (0, false);
        };
        let Some(sink) = graph.reference_sink_vertex() else {
            return (0, false);
        };
        fn walk(
            g: &SeqGraph,
            v: usize,
            sink: usize,
            on_path: &mut [bool],
            cap: usize,
            found: &mut usize,
        ) {
            if *found >= cap {
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
                walk(g, to, sink, on_path, cap, found);
            }
            on_path[v] = false;
        }
        let mut on_path = vec![false; graph.node_count()];
        let mut found = 0usize;
        walk(graph, src, sink, &mut on_path, PATH_CAP, &mut found);
        (found, found >= PATH_CAP)
    }

    fn n_branch(graph: &SeqGraph) -> usize {
        (0..graph.node_count())
            .filter(|&v| graph.outgoing_nodes(v).len() > 1)
            .count()
    }

    fn hap_snps(bases: &[u8], ref_bases: &[u8]) -> String {
        if bases.len() != ref_bases.len() {
            return format!("len={}!={}", bases.len(), ref_bases.len());
        }
        let sites = [
            (SITE_CT, b'T', "361T"),
            (SITE_CG, b'G', "371G"),
            (SITE_CA, b'A', "399A"),
            (SITE_TC, b'C', "407C"),
            (SITE_GC, b'C', "412C"),
        ];
        let mut parts = Vec::new();
        for &(pos, alt, name) in &sites {
            let i = (pos - JAVA_EXTENDED.0) as usize;
            if i >= bases.len() {
                parts.push(format!("{name}=OOB"));
                continue;
            }
            let hit = bases[i] == alt;
            parts.push(format!("{name}={}({})", hit, bases[i] as char));
        }
        parts.join(" ")
    }

    fn n_with_base(haps: &[Haplotype], off: usize, alt: u8) -> usize {
        haps.iter()
            .filter(|h| h.bases.len() > off && h.bases[off] == alt)
            .count()
    }

    fn dump_hap_bases(label: &str, bases: &[u8], ref_bases: &[u8], extra: &str) {
        eprintln!(
            "  {label} len={} {extra} {}",
            bases.len(),
            hap_snps(bases, ref_bases)
        );
    }

    fn dump_haps(label: &str, haps: &[Haplotype], ref_bases: &[u8]) {
        eprintln!("{label} n={}", haps.len());
        for (i, h) in haps.iter().enumerate() {
            dump_hap_bases(
                &format!("{label}[{i}]"),
                &h.bases,
                ref_bases,
                &format!("ref={}", h.is_reference),
            );
        }
        let off_ct = (SITE_CT - JAVA_EXTENDED.0) as usize;
        let off_ca = (SITE_CA - JAVA_EXTENDED.0) as usize;
        eprintln!(
            "{label} n_361T={} n_399A={}",
            n_with_base(haps, off_ct, b'T'),
            n_with_base(haps, off_ca, b'A')
        );
    }

    fn dump_kbest(label: &str, graph: &SeqGraph, ref_bases: &[u8]) -> (usize, usize) {
        let Ok(paths) =
            find_best_haplotypes_seq_graph(graph, DEFAULT_NUM_BEST_HAPLOTYPES_PER_GRAPH)
        else {
            eprintln!("{label} kbest=ERR");
            return (0, 0);
        };
        let mut uniq = HashSet::new();
        let off_ct = (SITE_CT - JAVA_EXTENDED.0) as usize;
        let mut n_361 = 0usize;
        eprintln!(
            "{label} kbest_paths={} cap={DEFAULT_NUM_BEST_HAPLOTYPES_PER_GRAPH}",
            paths.len()
        );
        for (i, p) in paths.iter().enumerate() {
            let b = graph.path_bases_bytes(p.start, &p.edges);
            if b.len() > off_ct && b[off_ct] == b'T' {
                n_361 += 1;
            }
            uniq.insert(b.clone());
            eprintln!(
                "  PATH[{i}] len={} score={:.6} ref_flag={} {}",
                b.len(),
                p.score,
                p.is_reference,
                hap_snps(&b, ref_bases)
            );
        }
        eprintln!("{label} unique_seq={} n_361T_paths={n_361}", uniq.len());
        (uniq.len(), n_361)
    }

    fn dump_seq_snap(label: &str, graph: &SeqGraph, ref_bases: &[u8]) -> (usize, usize, usize) {
        let (n_paths, capped) = count_st_paths(graph);
        eprintln!(
            "SEQ[{label}] nodes={} edges={} branch={} src={:?} sink={:?} st_paths={} capped={capped} motif361C={} motif361T={} java_alt_mid={} java_ref_mid={}",
            graph.node_count(),
            graph.edge_count(),
            n_branch(graph),
            graph.reference_source_vertex(),
            graph.reference_sink_vertex(),
            n_paths,
            seq_motif_vertex_count(graph, MOTIF_361_REF),
            seq_motif_vertex_count(graph, MOTIF_361_ALT),
            seq_motif_vertex_count(graph, JAVA_ALT_MID),
            seq_motif_vertex_count(graph, JAVA_REF_MID),
        );
        if graph.node_count() <= 8 {
            for v in 0..graph.node_count() {
                let seq = graph.vertex_sequence(v);
                if seq.len() <= 48 {
                    eprintln!(
                        "  VERT[{v}] out={} seq={}",
                        graph.outgoing_nodes(v).len(),
                        ascii(seq)
                    );
                } else {
                    eprintln!(
                        "  VERT[{v}] out={} len={} head={}...tail={}",
                        graph.outgoing_nodes(v).len(),
                        seq.len(),
                        ascii(&seq[..12.min(seq.len())]),
                        ascii(&seq[seq.len().saturating_sub(12)..])
                    );
                }
            }
        }
        if graph.reference_source_vertex().is_some() && graph.reference_sink_vertex().is_some() {
            let _ = dump_kbest(&format!("SEQ[{label}]"), graph, ref_bases);
        }
        (graph.node_count(), graph.edge_count(), n_paths)
    }

    fn dump_bam_site(rec: &rust_htslib::bam::Record, pos1: u64, name: &str) {
        let cigar = CigarString(rec.cigar().iter().copied().collect());
        let qname = String::from_utf8_lossy(rec.qname()).into_owned();
        let seq = rec.seq().as_bytes();
        let quals = rec.qual();
        let qi = query_index_at_reference_position(rec.pos(), &cigar, (pos1 - 1) as i64);
        match qi {
            Some(i) if i < seq.len() => {
                let q = if i < quals.len() { quals[i] } else { 255 };
                eprintln!(
                    "BAM[{name}] qname={qname} qi={i} base={} qual={q} usable_q10={}",
                    seq[i] as char,
                    q >= 10
                );
            }
            Some(i) => eprintln!("BAM[{name}] qname={qname} qi={i} OUT_OF_SEQ"),
            None => eprintln!("BAM[{name}] qname={qname} no query base"),
        }
    }

    #[test]
    fn six_r30_java_kbest_cap_is_not_two() {
        assert_eq!(DEFAULT_NUM_BEST_HAPLOTYPES_PER_GRAPH, 128);
        let simplify = include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/src/seq_graph_simplify.rs"
        ));
        assert!(simplify.contains("merge_diamonds_until_complete"));
        assert!(
            !simplify.contains("92317361"),
            "must not special-case 361 in SeqGraph simplify"
        );
        assert!(!simplify.contains("92317371"));
        let ctx = AssemblyScoringContext {
            padded_reference_start_1based: JAVA_EXTENDED.0,
            active_start_1based: JAVA_ACTIVE.0,
            active_end_1based: JAVA_ACTIVE.1,
            contig: "2".into(),
        };
        assert!(
            ctx.overlaps_p12_l_gate_interval(),
            "mid-B sits in the contig-2 L-gate window: RT-first is skipped; supplement does not early-stop"
        );
        assert!(
            !ctx.overlaps_p12_cluster(),
            "mid-B is not the P12 TTC cluster"
        );
    }

    #[test]
    fn six_r30_mid_b_seqgraph_stages_at_361_371() {
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
        let regions = crate::walker_traversal::flatten_assembly_regions(&walk);
        let region = regions
            .iter()
            .find(|r| {
                matches!(
                    call_disposition(r),
                    AssemblyRegionCallDisposition::ActiveFull
                ) && r.start.get() <= SITE_CA
                    && r.end.get() >= SITE_CA
            })
            .expect("ActiveFull mid-B")
            .clone();
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

        let mut assemble_args = AssembleReadsArgs::default();
        assemble_args.strict_java_assembly = true;
        let mut rt_args = assemble_args.assembler.clone();
        rt_args.dangling_java_exact = true;
        assert_eq!(rt_args.min_prune_factor, 2);
        assert!(rt_args.use_seq_graph);
        assert!(!rt_args.allow_non_unique_kmers_in_ref);

        for rec in &region.reads {
            dump_bam_site(rec, SITE_CT, "361");
            dump_bam_site(rec, SITE_CG, "371");
            dump_bam_site(rec, SITE_CA, "399");
        }

        let finalized = finalize_region_reads_for_assembly(
            &region.reads,
            &region,
            assemble_args.correct_overlapping_base_qualities,
            gatk_min_tail_quality_for_assembly(rt_args.min_base_quality),
            false,
        );
        let assembly_reads = records_to_assembly_reads(&finalized);
        eprintln!(
            "INPUT graph_ref_len={} n_bam={} n_assembly_reads={} prune={} k={K}",
            ref_bases.len(),
            region.reads.len(),
            assembly_reads.len(),
            rt_args.min_prune_factor
        );
        for (i, r) in assembly_reads.iter().enumerate() {
            let has_t = r
                .bases
                .windows(MOTIF_361_ALT.len())
                .any(|w| w == MOTIF_361_ALT);
            let has_c = r
                .bases
                .windows(MOTIF_361_REF.len())
                .any(|w| w == MOTIF_361_REF);
            let has_java_alt = r
                .bases
                .windows(JAVA_ALT_MID.len())
                .any(|w| w == JAVA_ALT_MID);
            let has_zip_alt = r.bases.windows(11).any(|w| w == b"TAGAGTTGAAG");
            eprintln!(
                "AREAD[{i}] len={} motif361T={has_t} motif361C={has_c} java_alt_mid={has_java_alt} zip_alt_11={has_zip_alt} seq={}",
                r.bases.len(),
                ascii(&r.bases)
            );
        }

        let params = graph_params(K);
        let (mut raw, summary) = assembly_graph_from_ref_and_reads_threading_with_summary(
            &graph_ref,
            &assembly_reads,
            &params,
        )
        .expect("raw rt");
        eprintln!(
            "RAW_RT low_complexity={} non_unique_ref skipped later if any",
            summary.is_low_complexity
        );
        dump_rt_kmers(
            "raw", &raw, ref_bases, off_ct, off_cg, off_ca, off_tc, off_gc,
        );
        let mut pruning = AssemblyGraphPruningParams::gatk_haplotype_caller_defaults();
        pruning.min_prune_factor = rt_args.min_prune_factor;
        raw.apply_pruning(&pruning);
        dump_rt_kmers(
            "after_prune",
            &raw,
            ref_bases,
            off_ct,
            off_cg,
            off_ca,
            off_tc,
            off_gc,
        );

        let graph = build_threading_graph_for_seq_assembly(
            &graph_ref,
            &assembly_reads,
            K,
            &rt_args,
            false,
            false,
        )
        .expect("rt")
        .expect("k=25 graph");
        dump_rt_kmers(
            "seq_assembly",
            &graph,
            ref_bases,
            off_ct,
            off_cg,
            off_ca,
            off_tc,
            off_gc,
        );
        assert!(
            motif_node_count(&graph, b"TAGAGTTGAAG") == 0,
            "6R.33: cleaned RT must drop TAGAGTTGAAG (Java 0.2)"
        );

        let mut seq = SeqGraph::from_assembly_graph(&graph);
        dump_seq_snap("from_rt", &seq, ref_bases);
        seq.clean_non_ref_paths();
        let mut last_topo = dump_seq_snap("after_clean_non_ref", &seq, ref_bases);

        let status = seq.traced_cleanup_seq_graph(|stage, g| {
            let (n_paths, _) = count_st_paths(g);
            let topo = (g.node_count(), g.edge_count(), n_paths);
            if topo == last_topo && stage != "final_for_kbest" {
                return;
            }
            last_topo = topo;
            dump_seq_snap(stage, g, ref_bases);
        });
        eprintln!("CLEANUP_STATUS={status:?}");
        let (seq_uniq, seq_n361) = dump_kbest("SEQ[final_repeat]", &seq, ref_bases);
        let (final_paths, _) = count_st_paths(&seq);
        assert_eq!(
            final_paths, 2,
            "6R.33: SeqGraph is a single diamond (2 s-t paths), matching Java 1.4.final"
        );
        assert_eq!(seq_uniq, 2);
        assert_eq!(seq_n361, 0);

        let seq_only = assemble_from_ref_and_reads(&graph_ref, &assembly_reads, &rt_args)
            .expect("seq-graph assemble without scoring");
        dump_haps("SEQ_ONLY_ASSEMBLE", &seq_only.haplotypes, ref_bases);

        let before_remove = extract_rt_haplotypes_before_remove_paths(
            &graph_ref,
            &assembly_reads,
            &rt_args,
            K,
            false,
            false,
        )
        .expect("before_remove");
        dump_haps("RT_BEFORE_REMOVE_k25", &before_remove, ref_bases);

        let mut owned = region.clone();
        let mut prod_args = AssembleReadsArgs::default();
        prod_args.strict_java_assembly = true;
        let prod = assemble_reads_with_finalized(&mut owned, &dict, &mut ref_cache, &prod_args)
            .expect("production assemble")
            .assembly;
        dump_haps("PROD_ASSEMBLE", &prod.haplotypes, ref_bases);

        eprintln!(
            "COMPARE seq_kbest_unique={seq_uniq} seq_kbest_361T={seq_n361} seq_only_n={} seq_only_361T={} before_remove_n={} before_remove_361T={} prod_n={} prod_361T={}",
            seq_only.haplotypes.len(),
            n_with_base(&seq_only.haplotypes, off_ct, b'T'),
            before_remove.len(),
            n_with_base(&before_remove, off_ct, b'T'),
            prod.haplotypes.len(),
            n_with_base(&prod.haplotypes, off_ct, b'T'),
        );

        assert_eq!(seq_only.haplotypes.len(), 2);
        assert_eq!(n_with_base(&seq_only.haplotypes, off_ct, b'T'), 0);

        assert!(
            prod.haplotypes.len() >= 2,
            "production retains REF + oracle ALT after 6R.33"
        );
        assert_eq!(
            n_with_base(&prod.haplotypes, off_ct, b'T'),
            0,
            "361T branch is not handed to production haplotypes"
        );
        assert!(
            n_with_base(&prod.haplotypes, off_ca, b'A') >= 1,
            "production must still retain 399A paths"
        );
    }
}
