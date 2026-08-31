//! SeqGraph Path.getBases regression (6R.5 probe / 6R.6 repair).
//! Distinguishes last-byte vs seq[k-1..] vs full stored SeqVertex concat.

#[cfg(test)]
mod traces {
    use super::super::*;
    use crate::assembly::AssemblyRead;
    use crate::assembly_graph_dump::{load_assembly_reads_tsv, load_assembly_ref_tsv};
    use crate::cigar::{Cigar, CigarOperator};
    use crate::haplotype::Haplotype;
    use crate::read_threading_assembler::{
        assemble_from_ref_and_reads, build_threading_graph_for_haplotype_dump,
        build_threading_graph_for_seq_assembly, extract_haplotypes_from_seq_kbest_paths,
        merge_rt_kbest_pre_remove_paths, ReadThreadingAssemblerArgs,
    };
    use crate::seq_kbest_haplotype::find_best_haplotypes_seq_graph;
    use std::collections::{HashMap, HashSet};

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

    fn graph_from(seqs: &[&[u8]], edges: &[(usize, usize, bool)], kmer_size: usize) -> SeqGraph {
        let mut g = SeqGraph {
            kmer_size,
            vertices: seqs
                .iter()
                .enumerate()
                .map(|(i, s)| SeqVertex {
                    id: i,
                    sequence: s.to_vec(),
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

    fn java_seqvertex_path_bases(g: &SeqGraph, start: usize, edges: &[(usize, usize)]) -> Vec<u8> {
        // GATK Path.getBases on SeqGraph: SeqVertex.getAdditionalSequence returns getSequence().
        let first = if edges.is_empty() { start } else { edges[0].0 };
        let mut bases = g.vertices[first].sequence.clone();
        for &(_, to) in edges {
            bases.extend_from_slice(&g.vertices[to].sequence);
        }
        bases
    }

    fn java_seq_kminus1_path_bases(
        g: &SeqGraph,
        start: usize,
        edges: &[(usize, usize)],
        k: usize,
    ) -> Vec<u8> {
        let first = if edges.is_empty() { start } else { edges[0].0 };
        let mut bases = g.vertices[first].sequence.clone();
        for &(_, to) in edges {
            let seq = &g.vertices[to].sequence;
            let skip = k.saturating_sub(1).min(seq.len());
            bases.extend_from_slice(&seq[skip..]);
        }
        bases
    }

    fn rust_additional(seq: &[u8], is_source: bool) -> Vec<u8> {
        additional_sequence_bytes(seq, is_source)
    }

    #[test]
    fn six_r5_minimal_zipped_vertex_path_bytes() {
        // Zipped-style SeqVertices (multi-base payloads already stored, as after zip).
        // last-byte vs seq[k-1..] vs full stored sequence are distinguishable.
        let g = graph_from(
            &[b"ABCDE", b"FGHIJXYZ", b"KLMNOPQR"],
            &[(0, 1, true), (1, 2, true)],
            5,
        );
        let edges = [(0, 1), (1, 2)];
        let rust = g.path_bases_bytes(0, &edges);
        let java_seqvertex = java_seqvertex_path_bases(&g, 0, &edges);
        let java_kminus1 = java_seq_kminus1_path_bases(&g, 0, &edges, 5);
        let last_byte = {
            let mut b = rust_additional(b"ABCDE", true);
            b.extend(rust_additional(b"FGHIJXYZ", false));
            b.extend(rust_additional(b"KLMNOPQR", false));
            b
        };
        eprintln!("=== 6R.5 MINIMAL PROBE ===");
        eprintln!("vertices: 0=ABCDE 1=FGHIJXYZ 2=KLMNOPQR k=5");
        eprintln!("edges: 0-R>1 1-R>2 path order 0,1,2");
        eprintln!(
            "rust additional_sequence_bytes(mid,false)={:?}",
            String::from_utf8_lossy(&rust_additional(b"FGHIJXYZ", false))
        );
        eprintln!(
            "rust additional_sequence_bytes(sink,false)={:?}",
            String::from_utf8_lossy(&rust_additional(b"KLMNOPQR", false))
        );
        eprintln!("rust path_bases_bytes={}", String::from_utf8_lossy(&rust));
        eprintln!(
            "expected last-byte concat={}",
            String::from_utf8_lossy(&last_byte)
        );
        eprintln!(
            "expected Java SeqVertex Path.getBases (full stored)={}",
            String::from_utf8_lossy(&java_seqvertex)
        );
        eprintln!(
            "expected seq[k-1..] concat={}",
            String::from_utf8_lossy(&java_kminus1)
        );

        assert_eq!(java_seqvertex, b"ABCDEFGHIJXYZKLMNOPQR");
        assert_eq!(last_byte, b"ABCDEZR");
        assert_eq!(java_kminus1, b"ABCDEJXYZOPQR");
        assert_eq!(
            rust, java_seqvertex,
            "path_bases_bytes must concat full stored SeqVertex sequences"
        );
        assert_eq!(rust, b"ABCDEFGHIJXYZKLMNOPQR");
        assert_ne!(
            rust, last_byte,
            "full stored concat must not equal last-byte restitch"
        );
        assert_ne!(
            rust, java_kminus1,
            "full stored concat must not equal seq[k-1..] concat"
        );
    }

    fn fixture(
        _name: &'static str,
        mid_ref: &str,
        mid_alt: &str,
    ) -> (String, String, AssemblyRead) {
        let ref_seq = format!("{LEFT}{mid_ref}{RIGHT}");
        let alt_seq = format!("{LEFT}{mid_alt}{RIGHT}");
        (ref_seq.clone(), alt_seq, read(&ref_seq, 30))
    }

    fn reads_for(ref_seq: &str, alt_seq: &str) -> Vec<AssemblyRead> {
        let mut out = Vec::new();
        for _ in 0..4 {
            out.push(read(ref_seq, 30));
            out.push(read(alt_seq, 30));
        }
        out
    }

    fn classify(bases: &[u8], ref_seq: &str, alt_seq: &str) -> &'static str {
        if bases == ref_seq.as_bytes() {
            "REF"
        } else if bases == alt_seq.as_bytes() {
            "ALT"
        } else if bases.len() < ref_seq.len().min(alt_seq.len()) {
            "TRUNCATED"
        } else {
            "OTHER"
        }
    }

    #[test]
    fn six_r5_unique_flank_and_assembler_origin() {
        let cases = vec![
            ("control_snp", "A", "C"),
            ("control_2bp_del", "TTCA", "TA"),
            ("control_2bp_ins", "TA", "TTCA"),
            ("control_two_snps", "AAAA", "ACCA"),
            ("control_del_plus_snp", "TTCA", "TAC"),
            ("holdout_a_TTCA_vs_TATG", "TTCA", "TATG"),
            ("holdout_b_GGAC_vs_GCTT", "GGAC", "GCTT"),
        ];
        let args = synthetic_args();
        eprintln!("=== 6R.5 UNIQUE-FLANK + ASSEMBLER ORIGIN ===");
        for (name, mid_ref, mid_alt) in cases {
            let (ref_seq, alt_seq, reference) = fixture(name, mid_ref, mid_alt);
            let reads = reads_for(&ref_seq, &alt_seq);
            let graph =
                build_threading_graph_for_seq_assembly(&reference, &reads, 10, &args, true, false)
                    .expect("build")
                    .expect("graph");
            let mut seq = SeqGraph::from_assembly_graph(&graph);
            seq.clean_non_ref_paths();
            let _ = seq.cleanup_seq_graph();
            let paths = find_best_haplotypes_seq_graph(&seq, 32).unwrap_or_default();
            eprintln!(
                "--- {name} nv={} ne={} kbest={} ---",
                seq.node_count(),
                seq.edge_count(),
                paths.len()
            );
            for (i, p) in paths.iter().enumerate() {
                let rust = seq.path_bases_bytes(p.start, &p.edges);
                let java = java_seqvertex_path_bases(&seq, p.start, &p.edges);
                let mut vids = vec![p.start];
                for &(_, to) in &p.edges {
                    vids.push(to);
                }
                let vseqs: Vec<String> = vids
                    .iter()
                    .map(|&v| String::from_utf8_lossy(seq.vertex_sequence(v)).into_owned())
                    .collect();
                let edge_flags: Vec<String> = p
                    .edges
                    .iter()
                    .map(|&(f, t)| {
                        format!("{f}-{}>{t}", if seq.edge_is_ref(f, t) { "R" } else { "A" })
                    })
                    .collect();
                eprintln!(
                    "  kbest{i} verts={:?} vseqs={:?} edges={} vcount={} flags={:?} rust_len={} java_full_len={} rust_class={} java_class={} rust={} java_full={}",
                    vids,
                    vseqs,
                    p.edges.len(),
                    vids.len(),
                    edge_flags,
                    rust.len(),
                    java.len(),
                    classify(&rust, &ref_seq, &alt_seq),
                    classify(&java, &ref_seq, &alt_seq),
                    String::from_utf8_lossy(&rust),
                    String::from_utf8_lossy(&java)
                );
                eprintln!(
                    "    rust_eq_java_full={} rust_eq_ref={} rust_eq_alt={} java_eq_ref={} java_eq_alt={} flag_is_ref={}",
                    rust == java,
                    rust == ref_seq.as_bytes(),
                    rust == alt_seq.as_bytes(),
                    java == ref_seq.as_bytes(),
                    java == alt_seq.as_bytes(),
                    p.is_reference
                );
                assert_eq!(
                    rust, java,
                    "{name} kbest{i}: path_bases_bytes must concat full stored SeqVertex sequences"
                );
            }
            if name == "control_snp" {
                let rust_set: HashSet<Vec<u8>> = paths
                    .iter()
                    .map(|p| seq.path_bases_bytes(p.start, &p.edges))
                    .collect();
                assert_eq!(paths.len(), 2, "SNP k-best path count");
                assert_eq!(
                    rust_set.len(),
                    2,
                    "SNP k-best must yield two distinct haplotypes"
                );
                assert!(
                    rust_set.contains(ref_seq.as_bytes()),
                    "SNP k-best must include full REF LEFT+A+RIGHT"
                );
                assert!(
                    rust_set.contains(alt_seq.as_bytes()),
                    "SNP k-best must include full ALT LEFT+C+RIGHT"
                );
                for p in &paths {
                    let rust = seq.path_bases_bytes(p.start, &p.edges);
                    assert_eq!(
                        rust.len(),
                        85,
                        "SNP path must be full biological length, not truncated"
                    );
                    assert_eq!(
                        rust,
                        java_seqvertex_path_bases(&seq, p.start, &p.edges),
                        "SNP path_bases_bytes must match SeqVertex full concat"
                    );
                }
            }

            let kbest_set: HashSet<Vec<u8>> = paths
                .iter()
                .map(|p| seq.path_bases_bytes(p.start, &p.edges))
                .collect();

            let mut ref_hap = Haplotype::new(ref_seq.as_bytes(), true);
            let mut ref_cigar = Cigar::new();
            ref_cigar.push(ref_hap.bases.len(), CigarOperator::Match);
            ref_hap.cigar = Some(ref_cigar);
            let ref_cigar_len = ref_hap.cigar.as_ref().unwrap().reference_length();
            match extract_haplotypes_from_seq_kbest_paths(
                &paths,
                &seq,
                10,
                &ref_hap,
                ref_cigar_len,
                &args.haplotype_to_reference_sw,
            ) {
                Ok(extracted) => {
                    let extracted_set: HashSet<Vec<u8>> =
                        extracted.iter().map(|h| h.bases.clone()).collect();
                    eprintln!(
                        "  seq_extract n={} bytes={:?}",
                        extracted.len(),
                        extracted
                            .iter()
                            .map(|h| format!(
                                "{}(len={},is_ref={})",
                                classify(&h.bases, &ref_seq, &alt_seq),
                                h.bases.len(),
                                h.is_reference
                            ))
                            .collect::<Vec<_>>()
                    );
                    let mut merged = extracted.clone();
                    if let Err(e) = merge_rt_kbest_pre_remove_paths(
                        &reference,
                        &reads,
                        &args,
                        &[10],
                        &mut merged,
                    ) {
                        eprintln!("  merge_rt err={e}");
                    } else {
                        let after: HashSet<Vec<u8>> =
                            merged.iter().map(|h| h.bases.clone()).collect();
                        let added: Vec<_> = after.difference(&extracted_set).cloned().collect();
                        eprintln!(
                            "  after_merge_rt n={} added={:?} all={:?}",
                            merged.len(),
                            added
                                .iter()
                                .map(|b| format!(
                                    "{}(len={})",
                                    classify(b, &ref_seq, &alt_seq),
                                    b.len()
                                ))
                                .collect::<Vec<_>>(),
                            merged
                                .iter()
                                .map(|h| format!(
                                    "{}(len={},is_ref={})",
                                    classify(&h.bases, &ref_seq, &alt_seq),
                                    h.bases.len(),
                                    h.is_reference
                                ))
                                .collect::<Vec<_>>()
                        );
                    }
                }
                Err(e) => eprintln!("  seq_extract err={e}"),
            }

            for (label, kmer_sizes, use_seq) in [
                ("sg_k10", vec![10usize], true),
                ("sg_k25", vec![25], true),
                ("sg_k10_25", vec![10, 25], true),
                ("rt_k10_25", vec![10, 25], false),
            ] {
                let mut a = synthetic_args();
                a.kmer_sizes = kmer_sizes;
                a.use_seq_graph = use_seq;
                match assemble_from_ref_and_reads(&reference, &reads, &a) {
                    Ok(r) => {
                        let mut origin = Vec::new();
                        for h in &r.haplotypes {
                            let from_kbest = kbest_set.contains(&h.bases);
                            origin.push(format!(
                                "{}(len={},kbest={},full_allele={})",
                                classify(&h.bases, &ref_seq, &alt_seq),
                                h.bases.len(),
                                from_kbest,
                                h.bases == ref_seq.as_bytes() || h.bases == alt_seq.as_bytes()
                            ));
                        }
                        eprintln!(
                            "  {label} n={} origin={:?} uniq={}",
                            r.haplotypes.len(),
                            origin,
                            r.haplotypes
                                .iter()
                                .map(|h| h.bases.as_slice())
                                .collect::<HashSet<_>>()
                                .len()
                        );
                    }
                    Err(e) => eprintln!("  {label} err={e}"),
                }
            }
        }

        {
            let (ref_seq, alt_seq, reference) = fixture("kmer_cmp_snp", "A", "C");
            let reads = reads_for(&ref_seq, &alt_seq);
            let mut sets = Vec::new();
            for (label, kmer_sizes) in [
                ("sg_k10", vec![10usize]),
                ("sg_k25", vec![25]),
                ("sg_k10_25", vec![10, 25]),
            ] {
                let mut a = synthetic_args();
                a.kmer_sizes = kmer_sizes;
                a.use_seq_graph = true;
                let r = assemble_from_ref_and_reads(&reference, &reads, &a).expect(label);
                let uniq: HashSet<Vec<u8>> = r.haplotypes.iter().map(|h| h.bases.clone()).collect();
                eprintln!(
                    "=== kmer_cmp {label} n={} uniq={} seqs={:?} ===",
                    r.haplotypes.len(),
                    uniq.len(),
                    uniq.iter()
                        .map(|b| format!("{}(len={})", classify(b, &ref_seq, &alt_seq), b.len()))
                        .collect::<Vec<_>>()
                );
                sets.push((label, uniq));
            }
            assert_eq!(
                sets[0].1, sets[1].1,
                "sg_k10 unique haplotypes must equal sg_k25"
            );
            assert_eq!(
                sets[1].1, sets[2].1,
                "sg_k25 unique haplotypes must equal sg_k10_25"
            );
        }

        let ref_only = fixture("reference_only", "A", "A");
        let ref_reads: Vec<_> = (0..8).map(|_| read(&ref_only.0, 30)).collect();
        if let Ok(Some(graph)) =
            build_threading_graph_for_seq_assembly(&ref_only.2, &ref_reads, 10, &args, true, false)
        {
            let mut seq = SeqGraph::from_assembly_graph(&graph);
            seq.clean_non_ref_paths();
            let _ = seq.cleanup_seq_graph();
            let paths = find_best_haplotypes_seq_graph(&seq, 32).unwrap_or_default();
            eprintln!(
                "--- reference_only nv={} ne={} kbest={} ---",
                seq.node_count(),
                seq.edge_count(),
                paths.len()
            );
            for (i, p) in paths.iter().enumerate() {
                let rust = seq.path_bases_bytes(p.start, &p.edges);
                let java = java_seqvertex_path_bases(&seq, p.start, &p.edges);
                eprintln!(
                    "  kbest{i} rust_len={} java_full_len={} eq={} rust={}",
                    rust.len(),
                    java.len(),
                    rust == java,
                    String::from_utf8_lossy(&rust)
                );
                assert_eq!(rust, java);
                assert_eq!(rust.as_slice(), ref_only.0.as_bytes());
            }
        }

        let repo = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("..");
        let p5_ref =
            load_assembly_ref_tsv(&repo.join("parity/fixtures/hc-full-parity/e4/p5_case1_ref.tsv"))
                .expect("p5 ref");
        let p5_reads = load_assembly_reads_tsv(
            &repo.join("parity/fixtures/hc-full-parity/e4/p5_case1_reads.tsv"),
        )
        .expect("p5 reads");
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
                .expect("p5 build")
                .expect("p5 graph");
        let mut seq = SeqGraph::from_assembly_graph(&p5_graph);
        seq.clean_non_ref_paths();
        let _ = seq.cleanup_seq_graph();
        let paths = find_best_haplotypes_seq_graph(&seq, 32).unwrap_or_default();
        eprintln!(
            "--- p5_case1 nv={} ne={} kbest={} ---",
            seq.node_count(),
            seq.edge_count(),
            paths.len()
        );
        for (i, p) in paths.iter().enumerate() {
            let rust = seq.path_bases_bytes(p.start, &p.edges);
            let java = java_seqvertex_path_bases(&seq, p.start, &p.edges);
            eprintln!(
                "  kbest{i} rust={} java_full={} eq={} verts_nv={}",
                String::from_utf8_lossy(&rust),
                String::from_utf8_lossy(&java),
                rust == java,
                seq.node_count()
            );
            assert_eq!(rust, java);
            assert_eq!(rust.as_slice(), b"ACGTT");
        }
    }

    fn snap_haps(stage: &str, haps: &[Haplotype], ref_seq: &str, alt_seq: &str) {
        let uniq: HashSet<&[u8]> = haps.iter().map(|h| h.bases.as_slice()).collect();
        let classes: Vec<String> = haps
            .iter()
            .map(|h| {
                format!(
                    "{}(len={},is_ref={})",
                    classify(&h.bases, ref_seq, alt_seq),
                    h.bases.len(),
                    h.is_reference
                )
            })
            .collect();
        let mut identical_pairs = 0usize;
        for i in 0..haps.len() {
            for j in (i + 1)..haps.len() {
                if haps[i].bases == haps[j].bases {
                    identical_pairs += 1;
                }
            }
        }
        eprintln!(
            "  {stage} objects={} unique={} identical_pairs={} classes={:?} uniq_lens={:?}",
            haps.len(),
            uniq.len(),
            identical_pairs,
            classes,
            uniq.iter().map(|b| b.len()).collect::<Vec<_>>()
        );
    }

    #[test]
    fn six_r7_haplotype_identity_stages() {
        use crate::haplotype::normalize_ref_equivalent_haplotypes;
        use crate::read_threading_assembler::finalize_assembly_haplotypes;

        let cases = vec![
            ("control_snp", "A", "C"),
            ("control_2bp_del", "TTCA", "TA"),
            ("control_2bp_ins", "TA", "TTCA"),
            ("holdout_a_TTCA_vs_TATG", "TTCA", "TATG"),
            ("holdout_b_GGAC_vs_GCTT", "GGAC", "GCTT"),
        ];
        let args = synthetic_args();
        eprintln!("=== 6R.7 HAPLOTYPE IDENTITY STAGES ===");
        for (name, mid_ref, mid_alt) in cases {
            let (ref_seq, alt_seq, reference) = fixture(name, mid_ref, mid_alt);
            let reads = reads_for(&ref_seq, &alt_seq);
            let graph =
                build_threading_graph_for_seq_assembly(&reference, &reads, 10, &args, true, false)
                    .expect("build")
                    .expect("graph");
            let mut seq = SeqGraph::from_assembly_graph(&graph);
            seq.clean_non_ref_paths();
            let _ = seq.cleanup_seq_graph();
            let paths = find_best_haplotypes_seq_graph(&seq, 32).unwrap_or_default();
            let recon: Vec<Vec<u8>> = paths
                .iter()
                .map(|p| seq.path_bases_bytes(p.start, &p.edges))
                .collect();
            let recon_uniq: HashSet<&[u8]> = recon.iter().map(|b| b.as_slice()).collect();
            eprintln!(
                "--- {name} kbest_paths={} recon_unique={} ---",
                paths.len(),
                recon_uniq.len()
            );
            assert!(
                recon_uniq.contains(ref_seq.as_bytes()),
                "{name}: reconstructed k-best must include REF"
            );
            assert!(
                recon_uniq.contains(alt_seq.as_bytes()),
                "{name}: reconstructed k-best must include ALT"
            );
            assert_eq!(
                recon_uniq.len(),
                2,
                "{name}: k-best reconstruction must be exactly REF+ALT unique sequences"
            );

            let mut ref_hap = Haplotype::new(ref_seq.as_bytes(), true);
            let mut ref_cigar = Cigar::new();
            ref_cigar.push(ref_hap.bases.len(), CigarOperator::Match);
            ref_hap.cigar = Some(ref_cigar);
            let ref_cigar_len = ref_hap.cigar.as_ref().unwrap().reference_length();
            let extracted = extract_haplotypes_from_seq_kbest_paths(
                &paths,
                &seq,
                10,
                &ref_hap,
                ref_cigar_len,
                &args.haplotype_to_reference_sw,
            )
            .expect("extract");
            snap_haps("extract", &extracted, &ref_seq, &alt_seq);
            let extract_uniq: HashSet<&[u8]> =
                extracted.iter().map(|h| h.bases.as_slice()).collect();
            assert!(
                extract_uniq.contains(alt_seq.as_bytes()),
                "{name}: ALT lost at extract"
            );
            assert_eq!(extract_uniq.len(), 2, "{name}: extract unique");

            let mut merged = extracted.clone();
            merge_rt_kbest_pre_remove_paths(&reference, &reads, &args, &[10], &mut merged)
                .expect("merge_rt");
            snap_haps("after_merge_rt", &merged, &ref_seq, &alt_seq);
            let merge_uniq: HashSet<&[u8]> = merged.iter().map(|h| h.bases.as_slice()).collect();
            assert!(
                merge_uniq.contains(alt_seq.as_bytes()),
                "{name}: ALT lost at merge_rt"
            );
            assert_eq!(
                merge_uniq.len(),
                2,
                "{name}: merge_rt must not add a third unique sequence"
            );

            let mut finalized = merged.clone();
            normalize_ref_equivalent_haplotypes(&mut finalized, ref_seq.as_bytes());
            snap_haps("after_normalize", &finalized, &ref_seq, &alt_seq);
            let status = finalize_assembly_haplotypes(&mut finalized, &ref_hap, true);
            snap_haps("after_finalize", &finalized, &ref_seq, &alt_seq);
            let final_uniq: HashSet<&[u8]> = finalized.iter().map(|h| h.bases.as_slice()).collect();
            assert!(
                final_uniq.contains(alt_seq.as_bytes()),
                "{name}: ALT lost at finalize status={status:?}"
            );
            assert!(
                final_uniq.contains(ref_seq.as_bytes()),
                "{name}: REF lost at finalize"
            );
            assert_eq!(final_uniq.len(), 2, "{name}: final unique must be REF+ALT");
            assert_eq!(finalized.len(), 2, "{name}: final objects after normalize");

            let mut sets = Vec::new();
            for (label, kmer_sizes, use_seq) in [
                ("sg_k10", vec![10usize], true),
                ("sg_k25", vec![25], true),
                ("sg_k10_25", vec![10, 25], true),
                ("rt_k10_25", vec![10, 25], false),
            ] {
                let mut a = synthetic_args();
                a.kmer_sizes = kmer_sizes;
                a.use_seq_graph = use_seq;
                let r = assemble_from_ref_and_reads(&reference, &reads, &a).expect(label);
                let uniq: HashSet<Vec<u8>> = r.haplotypes.iter().map(|h| h.bases.clone()).collect();
                eprintln!(
                    "  {label} objects={} unique={} lens={:?}",
                    r.haplotypes.len(),
                    uniq.len(),
                    uniq.iter().map(|b| b.len()).collect::<Vec<_>>()
                );
                assert_eq!(uniq.len(), 2, "{name} {label} unique");
                assert!(
                    uniq.contains(ref_seq.as_bytes()),
                    "{name} {label} missing REF"
                );
                assert!(
                    uniq.contains(alt_seq.as_bytes()),
                    "{name} {label} missing ALT"
                );
                sets.push(uniq);
            }
            assert_eq!(sets[0], sets[1], "{name}: k10 unique bytes == k25");
            assert_eq!(sets[1], sets[2], "{name}: k25 unique bytes == k10+25");
            assert_eq!(
                sets[2], sets[3],
                "{name}: SeqGraph combined unique bytes == RT"
            );
        }

        let (ref_seq, _, reference) = fixture("reference_only", "A", "A");
        let ref_reads: Vec<_> = (0..8).map(|_| read(&ref_seq, 30)).collect();
        let mut a = synthetic_args();
        a.use_seq_graph = true;
        let sg = assemble_from_ref_and_reads(&reference, &ref_reads, &a).expect("sg ref-only");
        a.use_seq_graph = false;
        let rt = assemble_from_ref_and_reads(&reference, &ref_reads, &a).expect("rt ref-only");
        let sg_uniq: HashSet<&[u8]> = sg.haplotypes.iter().map(|h| h.bases.as_slice()).collect();
        let rt_uniq: HashSet<&[u8]> = rt.haplotypes.iter().map(|h| h.bases.as_slice()).collect();
        eprintln!(
            "--- reference_only sg_objects={} sg_unique={} rt_objects={} rt_unique={} ---",
            sg.haplotypes.len(),
            sg_uniq.len(),
            rt.haplotypes.len(),
            rt_uniq.len()
        );
        assert_eq!(sg_uniq.len(), 1);
        assert_eq!(rt_uniq.len(), 1);
        assert!(sg_uniq.contains(ref_seq.as_bytes()));
        assert_eq!(sg_uniq, rt_uniq);
        assert!(
            sg.haplotypes.iter().all(|h| h.is_reference),
            "reference-only must not invent a non-ref haplotype"
        );
    }

    #[test]
    fn six_r7_path_bases_source_does_not_reapply_kmer_last_byte() {
        let src = include_str!("seq_graph.rs");
        let fn_start = src
            .find("pub fn path_bases_bytes")
            .expect("path_bases_bytes");
        let fn_body = &src[fn_start..];
        let fn_end = fn_body
            .find("\n    fn in_degree")
            .expect("end of path_bases_bytes");
        let body = &fn_body[..fn_end];
        assert!(
            !body.contains("additional_sequence_bytes"),
            "path_bases_bytes must not re-apply k-mer last-byte"
        );
        assert!(body.contains("extend_from_slice"));
        assert!(src.contains("let seq = additional_sequence_bytes(&node.kmer, is_source);"));
    }

    #[test]
    fn six_r7_p12_waiver_untouched() {
        let src = include_str!("assembly_based_caller.rs");
        let pin = "            assembler.use_seq_graph = false;\n            assembler.remove_paths_not_connected_to_ref = false;\n            assembler.skip_post_dangling_prune = true;";
        assert!(src.contains(pin));
        assert!(src.contains("strict_java_assembly"));
    }
}
