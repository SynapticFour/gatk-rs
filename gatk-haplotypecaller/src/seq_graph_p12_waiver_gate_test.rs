//! 6R.8 TEST-ONLY P12 SeqGraph waiver behavioral gate.
//!
//! Test-only assembler mode differs from production P12 **only** in `use_seq_graph = true`.
//! Production P12 flags kept on both arms:
//!   `remove_paths_not_connected_to_ref = false`
//!   `skip_post_dangling_prune = true`
//!   `dangling_java_exact = true`
//!   `scoring` overlapping the P12 TTC/ATG cluster (rt_first skips P12).
//!
//! Does not edit `assembly_based_caller` waiver. W-H1 remains OPEN.

#[cfg(test)]
mod traces {
    use super::super::*;
    use crate::assembly::AssemblyRead;
    use crate::assembly_graph_dump::{load_assembly_reads_tsv, load_assembly_ref_tsv};
    use crate::assembly_result_set::DEFAULT_MAX_MNP_DISTANCE;
    use crate::compatibility::coupled_indel::CoupledIndelCluster;
    use crate::event_map::{variation_events_for_haplotype, VariationEvent};
    use crate::haplotype::Haplotype;
    use crate::read_event_discovery::P12_CLUSTER_TTC_START;
    use crate::read_threading_assembler::{
        assemble_from_ref_and_reads, build_threading_graph_for_haplotype_dump,
        build_threading_graph_for_seq_assembly, AssemblyScoringContext, ReadThreadingAssemblerArgs,
    };
    use crate::seq_kbest_haplotype::find_best_haplotypes_seq_graph;
    use std::collections::HashSet;
    use std::path::Path;

    const LEFT: &str = "ACGTACGGTTAGCCATAACGGTCCATTGCATAGCTGGAACCT";
    const RIGHT: &str = "GCTTAGGAACCGGTTAACCGATCCTGAACCGGATCCATAGCT";
    /// Genomic pad so synthetic TTCA starts at P12 TTC locus.
    const COUPLED_PAD: u64 = P12_CLUSTER_TTC_START - LEFT.len() as u64;
    /// Real NA12878 20k P12 window, expressed as offsets from TTC start (N-1).
    const REAL_P12_ACTIVE_START: u64 = P12_CLUSTER_TTC_START - 96;
    const REAL_P12_ACTIVE_END: u64 = P12_CLUSTER_TTC_START + 76;
    const REAL_P12_ATG_START: u64 = P12_CLUSTER_TTC_START + 3;
    const REAL_P12_ALT_WIN_LO: u64 = P12_CLUSTER_TTC_START - 4;
    const REAL_P12_ALT_WIN_HI: u64 = P12_CLUSTER_TTC_START + 11;
    const REAL_P12_EVENT_WIN_LO: u64 = P12_CLUSTER_TTC_START - 24;
    const REAL_P12_EVENT_WIN_HI: u64 = P12_CLUSTER_TTC_START + 26;

    fn read(seq: &str, q: u8) -> AssemblyRead {
        AssemblyRead {
            bases: seq.as_bytes().to_vec(),
            base_quals: vec![q; seq.len()],
        }
    }

    fn p12_scoring() -> AssemblyScoringContext {
        AssemblyScoringContext {
            padded_reference_start_1based: COUPLED_PAD,
            active_start_1based: P12_CLUSTER_TTC_START,
            active_end_1based: P12_CLUSTER_TTC_START.saturating_add(3),
            contig: "2".into(),
        }
    }

    /// Production P12 flags; `use_seq_graph` is the only intended toggle.
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
            scoring: Some(p12_scoring()),
            ..Default::default()
        }
    }

    fn fixture(mid_ref: &str, mid_alt: &str) -> (String, String, AssemblyRead, Vec<AssemblyRead>) {
        let ref_seq = format!("{LEFT}{mid_ref}{RIGHT}");
        let alt_seq = format!("{LEFT}{mid_alt}{RIGHT}");
        let reference = read(&ref_seq, 30);
        let mut reads = Vec::new();
        for _ in 0..4 {
            reads.push(read(&ref_seq, 30));
            reads.push(read(&alt_seq, 30));
        }
        (ref_seq, alt_seq, reference, reads)
    }

    fn cigar_of(h: &Haplotype) -> String {
        h.cigar
            .as_ref()
            .map(|c| c.to_gatk_string())
            .unwrap_or_else(|| "-".into())
    }

    fn events_for(h: &Haplotype, ref_hap: &Haplotype, ref_bytes: &[u8]) -> Vec<VariationEvent> {
        variation_events_for_haplotype(
            h,
            ref_hap,
            ref_bytes,
            COUPLED_PAD,
            DEFAULT_MAX_MNP_DISTANCE,
            "2",
        )
    }

    fn classify_events(events: &[VariationEvent]) -> (bool, bool, bool, bool) {
        let ttc_t = events
            .iter()
            .any(|e| e.ref_allele == "TTC" && e.alt_allele == "T");
        let a_atg = events
            .iter()
            .any(|e| e.ref_allele == "A" && e.alt_allele == "ATG");
        let mnv = events.iter().any(|e| {
            (e.ref_allele == "TTCA" && e.alt_allele == "TATG")
                || (e.ref_allele.len() == e.alt_allele.len() && e.ref_allele.len() >= 3)
        });
        let coupled = events.iter().any(|d| {
            events
                .iter()
                .any(|i| CoupledIndelCluster::try_from_pair(d, i).is_some())
        });
        (ttc_t, a_atg, mnv, coupled)
    }

    fn uniq_bases(haps: &[Haplotype]) -> HashSet<Vec<u8>> {
        haps.iter().map(|h| h.bases.clone()).collect()
    }

    fn dump_haps(label: &str, haps: &[Haplotype], ref_seq: &str, alt_seq: &str) {
        let uniq = uniq_bases(haps);
        eprintln!("  {label} objects={} unique={}", haps.len(), uniq.len());
        for (i, h) in haps.iter().enumerate() {
            let kind = if h.bases == ref_seq.as_bytes() {
                "REF"
            } else if h.bases == alt_seq.as_bytes() {
                "ALT"
            } else {
                "OTHER"
            };
            eprintln!(
                "    hap{i} {kind} len={} is_ref={} cigar={} seq={}",
                h.bases.len(),
                h.is_reference,
                cigar_of(h),
                String::from_utf8_lossy(&h.bases)
            );
        }
    }

    #[test]
    fn six_r8_production_waiver_pin_unchanged() {
        let src = include_str!("assembly_based_caller.rs");
        let pin = "            assembler.use_seq_graph = false;\n            assembler.remove_paths_not_connected_to_ref = false;\n            assembler.skip_post_dangling_prune = true;";
        assert!(src.contains(pin));
        assert!(src.contains("region_overlaps_p12_cluster("));
        assert!(src.contains("strict_java_assembly"));
    }

    #[test]
    fn six_r8_synthetic_waiver_vs_seqgraph() {
        let cases = [
            ("reference_only", "A", "A"),
            ("snp", "A", "C"),
            ("del2", "TTCA", "TA"),
            ("ins2", "TA", "TTCA"),
            ("coupled_ttca_tatg", "TTCA", "TATG"),
        ];
        eprintln!("=== 6R.8 SYNTHETIC WAIVER vs SEQGRAPH ===");
        eprintln!(
            "test-only delta: use_seq_graph only; remove_paths=false skip_post_dangling=true dangling_java_exact=true P12 scoring"
        );
        for (name, mid_ref, mid_alt) in cases {
            let (ref_seq, alt_seq, reference, reads) = if name == "reference_only" {
                let (r, _, reference, _) = fixture(mid_ref, mid_alt);
                let reads: Vec<_> = (0..8).map(|_| read(&r, 30)).collect();
                (r.clone(), r, reference, reads)
            } else {
                fixture(mid_ref, mid_alt)
            };
            let variation = name != "reference_only";
            eprintln!("--- {name} ---");
            for (label, use_sg) in [("waiver_on", false), ("seqgraph_on", true)] {
                let args = p12_mode_args(use_sg);
                let r = assemble_from_ref_and_reads(&reference, &reads, &args).expect(label);
                eprintln!("  {label} status={:?}", r.status);
                dump_haps(label, &r.haplotypes, &ref_seq, &alt_seq);
                let uniq = uniq_bases(&r.haplotypes);
                if variation {
                    assert!(
                        uniq.contains(ref_seq.as_bytes()),
                        "{name} {label} missing REF"
                    );
                    assert!(
                        uniq.contains(alt_seq.as_bytes()),
                        "{name} {label} missing ALT bytes"
                    );
                } else {
                    assert_eq!(
                        uniq.len(),
                        1,
                        "{name} {label} must be one unique REF sequence"
                    );
                    assert!(uniq.contains(ref_seq.as_bytes()));
                    // PathState.is_reference may be false on a REF-byte object; unique bytes matter.
                }
            }
        }

        // Coupled fixture: haplotype bytes vs EventMap coupled vs MNV.
        let (ref_seq, alt_seq, reference, reads) = fixture("TTCA", "TATG");
        eprintln!("=== 6R.8 COUPLED INDEL EVENT/CIGAR TRACE ===");
        eprintln!(
            "haplotype ALT bytes TATG vs REF TTCA are identical for MNV TTCA→TATG and coupled TTC→T + A→ATG"
        );
        let mut ref_hap = Haplotype::new(ref_seq.as_bytes(), true);
        let mut c = crate::cigar::Cigar::new();
        c.push(ref_hap.bases.len(), crate::cigar::CigarOperator::Match);
        ref_hap.cigar = Some(c);

        for (label, use_sg) in [("waiver_on", false), ("seqgraph_on", true)] {
            let args = p12_mode_args(use_sg);
            if use_sg {
                if let Ok(Some(graph)) = build_threading_graph_for_seq_assembly(
                    &reference, &reads, 10, &args, true, false,
                ) {
                    let mut seq = SeqGraph::from_assembly_graph(&graph);
                    seq.clean_non_ref_paths();
                    let status = seq.cleanup_seq_graph();
                    let paths = find_best_haplotypes_seq_graph(&seq, 32).unwrap_or_default();
                    eprintln!(
                        "  seqgraph_on k=10 nv={} ne={} src={:?} sink={:?} cleanup={status:?} kbest={}",
                        seq.node_count(),
                        seq.edge_count(),
                        seq.reference_source_vertex(),
                        seq.reference_sink_vertex(),
                        paths.len()
                    );
                    for (i, p) in paths.iter().enumerate() {
                        let b = seq.path_bases_bytes(p.start, &p.edges);
                        eprintln!(
                            "    kbest{i} len={} eq_ref={} eq_alt={}",
                            b.len(),
                            b == ref_seq.as_bytes(),
                            b == alt_seq.as_bytes()
                        );
                    }
                }
            }
            let r = assemble_from_ref_and_reads(&reference, &reads, &args).expect(label);
            let mut saw_coupled = false;
            let mut saw_mnv = false;
            let mut saw_ttc = false;
            let mut saw_atg = false;
            for h in r
                .haplotypes
                .iter()
                .filter(|h| !h.bases.eq(ref_seq.as_bytes()))
            {
                let ev = events_for(h, &ref_hap, ref_seq.as_bytes());
                let (ttc, atg, mnv, coupled) = classify_events(&ev);
                saw_coupled |= coupled;
                saw_mnv |= mnv;
                saw_ttc |= ttc;
                saw_atg |= atg;
                eprintln!(
                    "  {label} alt cigar={} events={:?} ttc_t={ttc} a_atg={atg} mnv={mnv} coupled_pair={coupled}",
                    cigar_of(h),
                    ev.iter()
                        .map(|e| format!(
                            "{}:{} {}→{}",
                            e.start_1based.get(),
                            e.end_1based.get(),
                            e.ref_allele,
                            e.alt_allele
                        ))
                        .collect::<Vec<_>>()
                );
            }
            eprintln!(
                "  {label} summary coupled_pair={saw_coupled} ttc_t={saw_ttc} a_atg={saw_atg} mnvish={saw_mnv}"
            );
        }
    }

    #[test]
    fn six_r8_p5_case1_unchanged_under_seqgraph() {
        let repo = Path::new(env!("CARGO_MANIFEST_DIR")).join("..");
        let p5_ref =
            load_assembly_ref_tsv(&repo.join("parity/fixtures/hc-full-parity/e4/p5_case1_ref.tsv"))
                .expect("p5 ref");
        let p5_reads = load_assembly_reads_tsv(
            &repo.join("parity/fixtures/hc-full-parity/e4/p5_case1_reads.tsv"),
        )
        .expect("p5 reads");
        let args = ReadThreadingAssemblerArgs {
            kmer_sizes: vec![3],
            min_base_quality: 10,
            min_prune_factor: 2,
            min_dangling_branch_length: 4,
            recover_dangling_heads: true,
            ..Default::default()
        };
        let graph =
            build_threading_graph_for_haplotype_dump(&p5_ref, &p5_reads, 3, &args, true, false)
                .expect("p5 build")
                .expect("p5 graph");
        let mut seq = SeqGraph::from_assembly_graph(&graph);
        seq.clean_non_ref_paths();
        let _ = seq.cleanup_seq_graph();
        let paths = find_best_haplotypes_seq_graph(&seq, 32).unwrap_or_default();
        assert_eq!(paths.len(), 1);
        let rust = seq.path_bases_bytes(paths[0].start, &paths[0].edges);
        assert_eq!(rust.as_slice(), b"ACGTT");
        eprintln!(
            "p5_case1 nv={} ne={} kbest=1 bases=ACGTT",
            seq.node_count(),
            seq.edge_count()
        );
    }

    #[test]
    fn six_r8_real_p12_na12878_waiver_vs_seqgraph() {
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
            eprintln!("Real-data P12 comparison unavailable (missing BAM or hs37d5.simple.fa)");
            return;
        }
        let dict = SequenceDictionary::from_fasta_path(&ref_path).expect("dict");
        let interval = format!("2:{REAL_P12_ACTIVE_START}-{REAL_P12_ACTIVE_END}");
        let specs = parse_intervals_cli_string(&dict, &interval).expect("interval");
        let filters = ReadFilterParams::gatk_standard_hc();
        let cfg = WalkerTraversalConfig::gatk_haplotype_caller_production(100);
        let walk = traverse_assembly_region_walker(&dict, &specs, &ref_path, &bam, &filters, &cfg)
            .expect("walk");
        let regions = flatten_assembly_regions(&walk);
        let region = regions
            .iter()
            .find(|r| {
                matches!(
                    call_disposition(r),
                    AssemblyRegionCallDisposition::ActiveFull
                ) && r.start.get() <= P12_CLUSTER_TTC_START
                    && r.end.get() >= REAL_P12_ATG_START
            })
            .expect("cluster active region");
        let mut ref_cache = ReferenceWindowCache::new(ref_path.clone(), 4);
        let reference = assembly_reference_read(&dict, &mut ref_cache, region).expect("ref");
        let finalized = finalize_region_reads_for_assembly(
            &region.reads,
            region,
            true,
            gatk_min_tail_quality_for_assembly(10),
            false,
        );
        let reads = records_to_assembly_reads(&finalized);
        let (pad, _) = crate::assembly_region_finalize::padded_reference_loc(region, &dict);
        let win_lo = REAL_P12_ALT_WIN_LO;
        let win_hi = REAL_P12_ALT_WIN_HI;
        let off_lo = win_lo.saturating_sub(pad) as usize;
        let off_hi = (win_hi.saturating_sub(pad) as usize + 1).min(reference.bases.len());
        eprintln!(
            "=== 6R.8 REAL P12 NA12878 region={}:{}-{} pad={} ref_len={} nreads={} ===",
            region.contig,
            region.start.get(),
            region.end.get(),
            pad,
            reference.bases.len(),
            reads.len()
        );
        if off_lo < off_hi {
            eprintln!(
                "  REF window {win_lo}-{win_hi} {}",
                String::from_utf8_lossy(&reference.bases[off_lo..off_hi])
            );
        }

        let dump_args = ReadThreadingAssemblerArgs {
            use_seq_graph: true,
            remove_paths_not_connected_to_ref: false,
            skip_post_dangling_prune: true,
            dangling_java_exact: true,
            scoring: Some(AssemblyScoringContext {
                padded_reference_start_1based: pad,
                active_start_1based: region.start.get(),
                active_end_1based: region.end.get(),
                contig: region.contig.clone(),
            }),
            ..Default::default()
        };
        eprintln!("=== 6R.8 REAL P12 SEQGRAPH-ONLY k-best (no inner RT retry) ===");
        for kmer in [10usize, 25, 85] {
            match build_threading_graph_for_seq_assembly(
                &reference, &reads, kmer, &dump_args, true, true,
            ) {
                Ok(Some(graph)) => {
                    let mut seq = SeqGraph::from_assembly_graph(&graph);
                    seq.clean_non_ref_paths();
                    let status = seq.cleanup_seq_graph();
                    let paths = find_best_haplotypes_seq_graph(&seq, 32).unwrap_or_default();
                    eprintln!(
                        "  k={kmer} nv={} ne={} src={:?} sink={:?} cleanup={status:?} kbest={}",
                        seq.node_count(),
                        seq.edge_count(),
                        seq.reference_source_vertex(),
                        seq.reference_sink_vertex(),
                        paths.len()
                    );
                    for (i, p) in paths.iter().enumerate() {
                        let b = seq.path_bases_bytes(p.start, &p.edges);
                        let eq_ref = b.as_slice() == reference.bases.as_slice();
                        let slice = if off_lo < off_hi && off_hi <= b.len() {
                            String::from_utf8_lossy(&b[off_lo..off_hi]).into_owned()
                        } else {
                            format!("len={}", b.len())
                        };
                        eprintln!(
                            "    kbest{i} len={} eq_ref={eq_ref} window={slice}",
                            b.len()
                        );
                    }
                }
                Ok(None) => eprintln!("  k={kmer} graph=None"),
                Err(e) => eprintln!("  k={kmer} graph_err={e}"),
            }
        }

        if let Ok(prod) = crate::assemble_reads(
            region,
            &dict,
            &mut ref_cache,
            &crate::AssembleReadsArgs::default(),
        ) {
            let coupled =
                crate::read_event_discovery::cluster_coupled_events_from_assembly_haplotypes(
                    &prod,
                    "2",
                    region.start.get(),
                    region.end.get(),
                );
            eprintln!(
                "=== 6R.8 PRODUCTION assemble_reads (waiver ON, strict_java) haps={} kmer={} coupled={} events={:?} ===",
                prod.haplotypes.len(),
                prod.kmer_size_for_dump(),
                coupled.len(),
                coupled
                    .iter()
                    .map(|e| format!(
                        "{} {}→{}",
                        e.start_1based.get(),
                        e.ref_allele,
                        e.alt_allele
                    ))
                    .collect::<Vec<_>>()
            );
        }

        let mut rows = Vec::new();
        let mut hap_sets: Vec<HashSet<Vec<u8>>> = Vec::new();
        for (label, use_sg) in [("waiver_on", false), ("seqgraph_on", true)] {
            let args = ReadThreadingAssemblerArgs {
                use_seq_graph: use_sg,
                remove_paths_not_connected_to_ref: false,
                skip_post_dangling_prune: true,
                dangling_java_exact: true,
                scoring: Some(AssemblyScoringContext {
                    padded_reference_start_1based: pad,
                    active_start_1based: region.start.get(),
                    active_end_1based: region.end.get(),
                    contig: region.contig.clone(),
                }),
                ..Default::default()
            };
            let r = assemble_from_ref_and_reads(&reference, &reads, &args).expect(label);
            eprintln!("  {label} status={:?} kmer={}", r.status, r.kmer_size);
            let uniq = uniq_bases(&r.haplotypes);
            hap_sets.push(uniq.clone());
            let mut ref_hap = r
                .haplotypes
                .iter()
                .find(|h| h.is_reference)
                .cloned()
                .unwrap_or_else(|| Haplotype::new(reference.bases.as_slice(), true));
            if ref_hap.cigar.is_none() {
                let mut c = crate::cigar::Cigar::new();
                c.push(ref_hap.bases.len(), crate::cigar::CigarOperator::Match);
                ref_hap.cigar = Some(c);
            }
            let mut ttc = false;
            let mut atg = false;
            let mut coupled = false;
            let mut mnv = false;
            let mut alt_cigars = Vec::new();
            for h in &r.haplotypes {
                let window = if off_lo < off_hi && off_hi <= h.bases.len() {
                    String::from_utf8_lossy(&h.bases[off_lo..off_hi]).into_owned()
                } else {
                    format!("len={}", h.bases.len())
                };
                if h.bases == ref_hap.bases {
                    eprintln!("  {label} REF cigar={} window={window}", cigar_of(h));
                    continue;
                }
                alt_cigars.push(cigar_of(h));
                let ev = variation_events_for_haplotype(
                    h,
                    &ref_hap,
                    reference.bases.as_slice(),
                    pad,
                    DEFAULT_MAX_MNP_DISTANCE,
                    &region.contig,
                );
                let (tt, aa, mn, cp) = classify_events(&ev);
                ttc |= tt;
                atg |= aa;
                mnv |= mn;
                coupled |= cp;
                eprintln!(
                    "  {label} alt len={} cigar={} window={window} events={:?}",
                    h.bases.len(),
                    cigar_of(h),
                    ev.iter()
                        .filter(|e| e.start_1based.get() >= REAL_P12_EVENT_WIN_LO
                            && e.start_1based.get() <= REAL_P12_EVENT_WIN_HI)
                        .map(|e| format!(
                            "{} {}→{}",
                            e.start_1based.get(),
                            e.ref_allele,
                            e.alt_allele
                        ))
                        .collect::<Vec<_>>()
                );
            }
            eprintln!(
                "  {label} objects={} unique={} ttc_t={ttc} a_atg={atg} coupled={coupled} mnvish={mnv} alt_cigars={alt_cigars:?}",
                r.haplotypes.len(),
                uniq.len()
            );
            rows.push((label, ttc, atg, coupled, mnv, uniq.len()));
        }
        assert_eq!(rows.len(), 2);
        assert_eq!(
            hap_sets[0], hap_sets[1],
            "assembler haplotype byte sets must match between waiver ON and SeqGraph ON"
        );
        eprintln!(
            "real P12 compare waiver ttc={} atg={} coupled={} | seqgraph ttc={} atg={} coupled={} hap_bytes_equal=true",
            rows[0].1, rows[0].2, rows[0].3, rows[1].1, rows[1].2, rows[1].3
        );
    }
}
