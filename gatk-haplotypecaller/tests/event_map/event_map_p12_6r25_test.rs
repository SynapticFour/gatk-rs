//! 6R.25 TEST-ONLY: Java `refHaplotype.getBases()` k-mer-gate differential vs padded REF.
//! Does not change production k, unique-kmer gates, dangling, EventMap, genotyping, or W-H1.

#[cfg(test)]
mod traces {
    use crate::alignment::SwParameters;
    use crate::assembly::{AssemblyGraph, AssemblyGraphPruningParams, AssemblyRead};
    use crate::assembly_dangling_recovery::DanglingRecoveryParams;
    use crate::assembly_region_finalize::{
        assembly_reference_read, finalize_region_reads_for_assembly,
        gatk_min_tail_quality_for_assembly, padded_reference_loc, records_to_assembly_reads,
        reference_haplotype_for_assembly_region, GATK_REFERENCE_PADDING_FOR_ASSEMBLY,
    };
    use crate::event_map::collect_variation_events;
    use crate::haplotype::Haplotype;
    use crate::read_projection::query_index_at_reference_position;
    use crate::read_threading_assembler::{
        assemble_from_ref_and_reads, build_threading_graph_for_seq_assembly,
        extract_haplotypes_from_seq_kbest_paths, ReadThreadingAssemblerArgs,
    };
    use crate::read_threading_graph::{
        assembly_graph_from_ref_and_reads_threading_with_summary, reference_has_non_unique_kmers,
    };
    use crate::seq_graph::SeqGraph;
    use crate::seq_kbest_haplotype::find_best_haplotypes_seq_graph;
    use crate::walker_apply::{call_disposition, AssemblyRegionCallDisposition};
    use gatk_core::reference::{
        parse_intervals_cli_string, ReferenceWindowCache, SequenceDictionary,
    };
    use rust_htslib::bam::record::CigarString;
    use std::collections::HashSet;
    use std::path::Path;

    const SITE_CA: u64 = 92_317_399;
    const SITE_TC: u64 = 92_317_407;
    const SITE_GC: u64 = 92_317_412;
    const JAVA_ACTIVE: (u64, u64) = (92_317_262, 92_317_491);
    const JAVA_EXTENDED: (u64, u64) = (92_317_162, 92_317_591);
    const K10: usize = 10;
    const K25: usize = 25;

    fn fixture_paths() -> Option<(std::path::PathBuf, std::path::PathBuf)> {
        let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("..");
        let ref_path = root.join("parity/realworld/assets/hs37d5.simple.fa");
        let bam = root.join("parity/realworld/na12878_20k_b37/NA12878_20k.b37.bam");
        if !ref_path.is_file() || !bam.is_file() {
            return None;
        }
        Some((ref_path, bam))
    }

    fn unique_kmers(bases: &[u8], k: usize) -> bool {
        !reference_has_non_unique_kmers(
            &AssemblyRead {
                bases: bases.to_vec(),
                base_quals: vec![30; bases.len()],
            },
            k,
        )
    }

    fn as_assembly_read(bases: &[u8]) -> AssemblyRead {
        AssemblyRead {
            bases: bases.to_vec(),
            base_quals: vec![30; bases.len()],
        }
    }

    fn prod_args() -> ReadThreadingAssemblerArgs {
        let mut a = ReadThreadingAssemblerArgs::default();
        a.dangling_java_exact = true;
        a
    }

    fn base_at(seq: &[u8], span_start_1based: u64, site: u64) -> Option<u8> {
        let off = site.checked_sub(span_start_1based)? as usize;
        seq.get(off).copied()
    }

    fn hap_has_sub(
        hap: &[u8],
        java_ref: &[u8],
        ext_start: u64,
        site: u64,
        expect_ref: u8,
        expect_alt: u8,
    ) -> bool {
        match (
            base_at(hap, ext_start, site),
            base_at(java_ref, ext_start, site),
        ) {
            (Some(h), Some(r)) => r == expect_ref && h == expect_alt,
            _ => false,
        }
    }

    fn alt_kmers_from_reads(reads: &[AssemblyRead], ref_bases: &[u8], k: usize) -> Vec<Vec<u8>> {
        if ref_bases.len() < k {
            return Vec::new();
        }
        let ref_set: HashSet<&[u8]> = (0..=ref_bases.len() - k)
            .map(|i| &ref_bases[i..i + k])
            .collect();
        let mut out = Vec::new();
        for ar in reads {
            if ar.bases.len() < k {
                continue;
            }
            for i in 0..=ar.bases.len() - k {
                let km = &ar.bases[i..i + k];
                if !ref_set.contains(km) {
                    out.push(km.to_vec());
                }
            }
        }
        out
    }

    fn alt_kmer_on_ref_connected_path(
        graph: &AssemblyGraph,
        alt_kmers: &[Vec<u8>],
    ) -> (usize, usize) {
        let mut present = 0usize;
        let mut touching_ref = 0usize;
        for km in alt_kmers {
            let Some(id) = graph.vertex_id_for_kmer(km) else {
                continue;
            };
            present += 1;
            let neighbor_ref = graph
                .outgoing_nodes(id)
                .iter()
                .any(|&n| graph.ref_nodes.contains(&n))
                || graph
                    .incoming_nodes(id)
                    .iter()
                    .any(|&n| graph.ref_nodes.contains(&n))
                || graph.ref_nodes.contains(&id);
            if neighbor_ref {
                touching_ref += 1;
            }
        }
        (present, touching_ref)
    }

    fn cigar_str(h: &Haplotype) -> String {
        h.cigar
            .as_ref()
            .map(|c| c.to_gatk_string())
            .unwrap_or_else(|| "NA".into())
    }

    fn site_base_char(seq: &[u8], span_start_1based: u64, site: u64) -> char {
        base_at(seq, span_start_1based, site)
            .map(|b| b as char)
            .unwrap_or('?')
    }

    fn hap_carries_java_three_snps(hap: &[u8], java_ref: &[u8], ext_start: u64) -> bool {
        hap_has_sub(hap, java_ref, ext_start, SITE_CA, b'C', b'A')
            && hap_has_sub(hap, java_ref, ext_start, SITE_TC, b'T', b'C')
            && hap_has_sub(hap, java_ref, ext_start, SITE_GC, b'G', b'C')
    }

    fn dump_hap_sites(label: &str, haps: &[Haplotype], java_ref: &[u8], ext_start: u64) {
        eprintln!(
            "REF_SITES {label} C@{}={} T@{}={} G@{}={}",
            SITE_CA,
            site_base_char(java_ref, ext_start, SITE_CA),
            SITE_TC,
            site_base_char(java_ref, ext_start, SITE_TC),
            SITE_GC,
            site_base_char(java_ref, ext_start, SITE_GC)
        );
        for (i, h) in haps.iter().enumerate() {
            let seq = h.bases.as_slice();
            eprintln!(
                "  HAP_SITES[{i}] ref={} 399={} 407={} 412={} three_java_snps={}",
                h.is_reference,
                site_base_char(seq, ext_start, SITE_CA),
                site_base_char(seq, ext_start, SITE_TC),
                site_base_char(seq, ext_start, SITE_GC),
                hap_carries_java_three_snps(seq, java_ref, ext_start)
            );
        }
    }

    fn greedy_unique_k25(len: usize) -> Vec<u8> {
        let mut s: Vec<u8> = (0..K25).map(|i| b"ACGT"[i % 4]).collect();
        let mut seen = HashSet::new();
        seen.insert(s.clone());
        while s.len() < len {
            let mut placed = false;
            for &b in &[b'A', b'C', b'G', b'T'] {
                s.push(b);
                let km = s[s.len() - K25..].to_vec();
                if seen.insert(km) {
                    placed = true;
                    break;
                }
                s.pop();
            }
            assert!(placed, "could not extend unique k=25 sequence");
        }
        s
    }

    /// Control C: uniqueness depends on which reference bytes are presented, not on k itself.
    #[test]
    fn six_r25_synthetic_k25_uniqueness_depends_on_ref_length() {
        // 40 unique bases via mixed prime-stride (no period-k repeat).
        let short = greedy_unique_k25(40);
        assert!(
            unique_kmers(&short, K25),
            "synthetic short ref must be unique at k=25 (got {:?})",
            String::from_utf8_lossy(&short)
        );
        let mut padded = short.clone();
        padded.extend_from_slice(&short[..25]);
        assert!(
            !unique_kmers(&padded, K25),
            "duplicating a 25-mer must make k=25 non-unique"
        );
        eprintln!(
            "CONTROL_C short_len={} padded_len={} k25_short_unique={} k25_padded_unique={}",
            short.len(),
            padded.len(),
            unique_kmers(&short, K25),
            unique_kmers(&padded, K25)
        );
    }

    #[test]
    fn six_r25_java_reference_haplotype_k25() {
        let Some((ref_fasta, bam_path)) = fixture_paths() else {
            eprintln!("Real-data mid-B comparison unavailable");
            return;
        };

        let assemble_src = include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/src/assembly_based_caller.rs"
        ));
        assert!(
            assemble_src.contains(
                "let reference = assembly_reference_read(dictionary, ref_cache, region)?;"
            ),
            "production assembleReads must still pad ±500 for error correction / collapse"
        );
        assert!(
            assemble_src.contains("create_graph_reference_read(&reference, region, dictionary)"),
            "production createGraph must slice getAssemblyRegionReference padding=0"
        );
        assert!(
            assemble_src.contains("assemble_from_ref_and_reads(&graph_ref, &reads, &assembler)?"),
            "production uniqueness + addSequence(ref) must share the extended-span graph_ref"
        );
        let rt_src = include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/src/read_threading_assembler.rs"
        ));
        assert!(
            rt_src.contains("&& reference_has_non_unique_kmers(reference, kmer_size)"),
            "production createGraph uniqueness must still be evaluated on the assembler `reference` argument"
        );
        assert!(
            !ReadThreadingAssemblerArgs::default().allow_non_unique_kmers_in_ref,
            "production must not globally lift allowNonUniqueKmersInRef"
        );

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

        assert_eq!(region.start.get(), JAVA_ACTIVE.0);
        assert_eq!(region.end.get(), JAVA_ACTIVE.1);
        assert_eq!(region.extended_start.get(), JAVA_EXTENDED.0);
        assert_eq!(region.extended_end.get(), JAVA_EXTENDED.1);
        assert_eq!(region.reads.len(), 2);

        let mut ref_cache = ReferenceWindowCache::new(ref_fasta.clone(), 4);
        let padded = assembly_reference_read(&dict, &mut ref_cache, &region).expect("pad");
        let (pad_start, pad_end) = padded_reference_loc(&region, &dict);
        let java_hap = reference_haplotype_for_assembly_region(&padded, &region, &dict);
        let java_ref = java_hap.bases.as_slice();
        let ext_start = region.extended_start.get();
        let ext_end = region.extended_end.get();
        let expected_ext_len = (ext_end - ext_start + 1) as usize;

        eprintln!("PINNED_JAVA GATK 4.4.0.0 SHA=2dbc025821bc5f686c423ff332a41e6cef892a77");
        eprintln!(
            "OBSERVED_JAVA_SOURCE createGraph uniqueness = ReadThreadingGraph.determineNonUniqueKmers(refHaplotype.getBases(), k)"
        );
        eprintln!(
            "OBSERVED_JAVA_SOURCE createReferenceHaplotype bases = AssemblyRegion.getAssemblyRegionReference(reader) padding=0 (padded/extended span), not ±500"
        );
        eprintln!(
            "REGION active={}-{} extended={}-{} pad={pad_start}-{pad_end} extra_pad={}",
            region.start.get(),
            region.end.get(),
            ext_start,
            ext_end,
            GATK_REFERENCE_PADDING_FOR_ASSEMBLY
        );
        eprintln!(
            "REF_LEN java_haplotype={} padded_assembly_read={} expected_extended={expected_ext_len} align_start_wrt_pad={}",
            java_ref.len(),
            padded.bases.len(),
            java_hap.alignment_start_hap_wrt_ref
        );

        assert_eq!(
            java_ref.len(),
            expected_ext_len,
            "Java-equivalent haplotype must be the extended-span slice, not inferred from VCF"
        );
        assert_eq!(
            java_hap.alignment_start_hap_wrt_ref,
            (ext_start - pad_start) as usize
        );

        let java_ref_k10_unique = unique_kmers(java_ref, K10);
        let java_ref_k25_unique = unique_kmers(java_ref, K25);
        let production_ref_k10_unique = unique_kmers(&padded.bases, K10);
        let production_ref_k25_unique = unique_kmers(&padded.bases, K25);
        eprintln!("UNIQUE k=10 java_ref={java_ref_k10_unique} padded={production_ref_k10_unique}");
        eprintln!("UNIQUE k=25 java_ref={java_ref_k25_unique} padded={production_ref_k25_unique}");
        assert!(
            !java_ref_k10_unique,
            "java_ref_k10_unique == false (live Java skipped k=10)"
        );
        assert!(
            java_ref_k25_unique,
            "java_ref_k25_unique == true (live Java accepted k=25)"
        );
        assert!(
            !production_ref_k25_unique,
            "production_ref_k25_unique == false"
        );
        assert!(
            !production_ref_k10_unique,
            "padded REF is also non-unique at k=10"
        );

        let raw_reads: Vec<rust_htslib::bam::Record> =
            region.reads.iter().map(|s| s.as_ref().clone()).collect();
        let finalized = finalize_region_reads_for_assembly(
            &region.reads,
            &region,
            true,
            gatk_min_tail_quality_for_assembly(10),
            false,
        );
        let assembly_reads = records_to_assembly_reads(&finalized);
        assert_eq!(assembly_reads.len(), 2);
        eprintln!(
            "PREPROCESS harness_reads=finalized_hardclip n_raw={} n_finalized={} (Java: well-defined fragment => revertSoftClips then assembler hardClipSoftClips)",
            raw_reads.len(),
            finalized.len()
        );
        for (i, rec) in raw_reads.iter().enumerate() {
            let fin = &finalized[i];
            let raw_cigar = CigarString(rec.cigar().iter().copied().collect());
            let fin_cigar = CigarString(fin.cigar().iter().copied().collect());
            eprintln!(
                "READ[{i}] qname={} flags={} mapq={} aln_start1={} raw_cigar={raw_cigar} raw_seq_len={} fin_cigar={fin_cigar} fin_seq_len={} fin_seq={}",
                String::from_utf8_lossy(rec.qname()),
                rec.flags(),
                rec.mapq(),
                rec.pos() + 1,
                rec.seq().len(),
                fin.seq().len(),
                String::from_utf8_lossy(&fin.seq().as_bytes())
            );
            for site in [SITE_CA, SITE_TC, SITE_GC] {
                let qi =
                    query_index_at_reference_position(fin.pos(), &fin_cigar, (site - 1) as i64);
                match qi {
                    Some(j) if j < fin.seq().len() => {
                        eprintln!(
                            "  fin_site {site} qi={j} base={} qual={}",
                            fin.seq().as_bytes()[j] as char,
                            fin.qual().get(j).copied().unwrap_or(255)
                        );
                    }
                    other => eprintln!("  fin_site {site} qi={other:?}"),
                }
            }
        }

        let args = prod_args();
        let java_ref_read = as_assembly_read(java_ref);
        let java_ref_k25_graph_built = build_threading_graph_for_seq_assembly(
            &java_ref_read,
            &assembly_reads,
            K25,
            &args,
            args.allow_low_complexity_graphs,
            args.allow_non_unique_kmers_in_ref,
        )
        .expect("java k25 build")
        .is_some();
        let production_ref_k25_graph_built = build_threading_graph_for_seq_assembly(
            &padded,
            &assembly_reads,
            K25,
            &args,
            args.allow_low_complexity_graphs,
            args.allow_non_unique_kmers_in_ref,
        )
        .expect("prod k25 build")
        .is_some();
        eprintln!(
            "GRAPH_BUILT k=25 java_ref={java_ref_k25_graph_built} production_padded={production_ref_k25_graph_built}"
        );
        assert!(java_ref_k25_graph_built, "java_ref_k25_graph_built == true");
        assert!(
            !production_ref_k25_graph_built,
            "production_ref_k25_graph_built == false"
        );
        let java_ref_k10_graph_built = build_threading_graph_for_seq_assembly(
            &java_ref_read,
            &assembly_reads,
            K10,
            &args,
            args.allow_low_complexity_graphs,
            args.allow_non_unique_kmers_in_ref,
        )
        .expect("java k10 build")
        .is_some();
        assert!(
            !java_ref_k10_graph_built,
            "k=10 must remain rejected on the Java-equivalent haplotype"
        );

        let alt_kmers = alt_kmers_from_reads(&assembly_reads, java_ref, K25);
        let (raw_g, raw_sum) = assembly_graph_from_ref_and_reads_threading_with_summary(
            &java_ref_read,
            &assembly_reads,
            &crate::assembly::AssemblyGraphParams {
                kmer_size: crate::bio_ids::KmerSize::try_from_usize(K25).expect("k"),
                min_base_quality: args.min_base_quality,
                min_edge_weight: 1,
                dangling_path_max_nodes: 0,
                max_haplotypes: args.num_best_haplotypes_per_graph,
                max_haplotype_bases: 4096,
                start_threading_only_at_existing_vertex: !args.recover_dangling_branches,
            },
        )
        .expect("raw k25");
        let (raw_present, raw_touch) = alt_kmer_on_ref_connected_path(&raw_g, &alt_kmers);
        eprintln!(
            "LIFECYCLE raw nodes={} edges={} low_complexity={} alt_25mers={} in_graph={} touching_ref={}",
            raw_g.node_count(),
            raw_g.edge_count(),
            raw_sum.is_low_complexity,
            alt_kmers.len(),
            raw_present,
            raw_touch
        );

        let mut pruned = raw_g.clone();
        let mut pruning = AssemblyGraphPruningParams::gatk_haplotype_caller_defaults();
        pruning.min_prune_factor = args.min_prune_factor;
        let n_pruned = pruned.apply_pruning(&pruning);
        let (pr_present, pr_touch) = alt_kmer_on_ref_connected_path(&pruned, &alt_kmers);
        eprintln!(
            "LIFECYCLE prune removed={n_pruned} nodes={} edges={} alt_in_graph={pr_present} touching_ref={pr_touch}",
            pruned.node_count(),
            pruned.edge_count()
        );

        let mut dangled = pruned.clone();
        let dangling = DanglingRecoveryParams::from_assembler_args(&args);
        let dang_sum = dangled
            .recover_dangling_branches(&dangling)
            .expect("dangling");
        let (da_present, da_touch) = alt_kmer_on_ref_connected_path(&dangled, &alt_kmers);
        eprintln!(
            "LIFECYCLE dangling tails={}/{} heads={}/{} merge_haps={} nodes={} edges={} alt_in_graph={da_present} touching_ref={da_touch}",
            dang_sum.tails_recovered,
            dang_sum.tails_attempted,
            dang_sum.heads_recovered,
            dang_sum.heads_attempted,
            dangled.dangling_merge_haps.len(),
            dangled.node_count(),
            dangled.edge_count()
        );

        let mut connected = dangled.clone();
        connected
            .remove_paths_not_connected_to_ref()
            .expect("remove_paths");
        let (rm_present, rm_touch) = alt_kmer_on_ref_connected_path(&connected, &alt_kmers);
        eprintln!(
            "LIFECYCLE remove_paths nodes={} edges={} alt_in_graph={rm_present} touching_ref={rm_touch}",
            connected.node_count(),
            connected.edge_count()
        );

        let mut seq = SeqGraph::from_assembly_graph(&connected);
        seq.clean_non_ref_paths();
        let cleanup = seq.cleanup_seq_graph();
        eprintln!(
            "LIFECYCLE seqgraph nodes={} edges={} cleanup={cleanup:?} ref_src={} ref_sink={}",
            seq.node_count(),
            seq.edge_count(),
            seq.reference_source_vertex().is_some(),
            seq.reference_sink_vertex().is_some()
        );

        let paths = find_best_haplotypes_seq_graph(&seq, args.num_best_haplotypes_per_graph)
            .expect("seq kbest");
        let n_alt_paths = paths.iter().filter(|p| !p.is_reference).count();
        eprintln!(
            "LIFECYCLE kbest n_paths={} n_alt_flag={}",
            paths.len(),
            n_alt_paths
        );

        let haps = extract_haplotypes_from_seq_kbest_paths(
            &paths,
            &seq,
            K25,
            &java_hap,
            java_ref.len(),
            &SwParameters::gatk_haplotype_to_reference(),
        )
        .expect("extract");
        let n_ref = haps.iter().filter(|h| h.is_reference).count();
        let alts: Vec<&Haplotype> = haps.iter().filter(|h| !h.is_reference).collect();
        eprintln!(
            "HAPLOTYPES n={} n_ref={n_ref} n_alt={} (from seq-kbest extract)",
            haps.len(),
            alts.len()
        );
        for (i, h) in haps.iter().enumerate() {
            eprintln!(
                "  HAP[{i}] ref={} len={} cigar={} score={} align_start={}",
                h.is_reference,
                h.bases.len(),
                cigar_str(h),
                h.score,
                h.alignment_start_hap_wrt_ref
            );
        }

        dump_hap_sites("seq_kbest", &haps, java_ref, ext_start);
        let seq_kbest_three = haps
            .iter()
            .any(|h| hap_carries_java_three_snps(&h.bases, java_ref, ext_start));
        eprintln!("SEQ_KBEST any_hap_has_three_java_snps={seq_kbest_three}");

        let assembled = assemble_from_ref_and_reads(&java_ref_read, &assembly_reads, &args)
            .expect("assemble java-ref");
        eprintln!(
            "ASSEMBLE_JAVA_REF status={:?} kmer={} n_haps={}",
            assembled.status,
            assembled.kmer_size,
            assembled.haplotypes.len()
        );
        let a_ref = assembled
            .haplotypes
            .iter()
            .filter(|h| h.is_reference)
            .count();
        let a_alts: Vec<&Haplotype> = assembled
            .haplotypes
            .iter()
            .filter(|h| !h.is_reference)
            .collect();
        eprintln!("ASSEMBLE_JAVA_REF n_ref={a_ref} n_alt={}", a_alts.len());
        for (i, h) in assembled.haplotypes.iter().enumerate() {
            eprintln!(
                "  ASSEMBLE_HAP[{i}] ref={} len={} cigar={} score={}",
                h.is_reference,
                h.bases.len(),
                cigar_str(h),
                h.score
            );
        }
        dump_hap_sites(
            "assemble_from_ref_and_reads",
            &assembled.haplotypes,
            java_ref,
            ext_start,
        );
        let assemble_three = assembled
            .haplotypes
            .iter()
            .any(|h| hap_carries_java_three_snps(&h.bases, java_ref, ext_start));
        eprintln!("ASSEMBLE any_hap_has_three_java_snps={assemble_three}");
        let n_430m = assembled
            .haplotypes
            .iter()
            .filter(|h| {
                cigar_str(h) == format!("{}M", java_ref.len()) && h.bases.len() == java_ref.len()
            })
            .count();
        eprintln!(
            "ALT_STRUCTURE n_430M={n_430m}/{} java_len={}",
            assembled.haplotypes.len(),
            java_ref.len()
        );

        let event_haps: &[Haplotype] = if !assembled.haplotypes.is_empty() {
            assembled.haplotypes.as_slice()
        } else {
            haps.as_slice()
        };
        let events = collect_variation_events(event_haps, java_ref, ext_start, &region.contig, 0);
        eprintln!("EVENTMAP n={}", events.len());
        for e in &events {
            eprintln!(
                "  EVENT {} {}/{} indel={}",
                e.start_1based.get(),
                e.ref_allele,
                e.alt_allele,
                e.is_indel()
            );
        }
        let has_ca = events
            .iter()
            .any(|e| e.start_1based.get() == SITE_CA && e.ref_allele == "C" && e.alt_allele == "A");
        let has_tc = events
            .iter()
            .any(|e| e.start_1based.get() == SITE_TC && e.ref_allele == "T" && e.alt_allele == "C");
        let has_gc = events
            .iter()
            .any(|e| e.start_1based.get() == SITE_GC && e.ref_allele == "G" && e.alt_allele == "C");
        eprintln!("EVENTMAP_SITES C/A={has_ca} T/C={has_tc} G/C={has_gc}");
        let extra: Vec<_> = events
            .iter()
            .filter(|e| {
                let s = e.start_1based.get();
                !(s == SITE_CA || s == SITE_TC || s == SITE_GC)
            })
            .collect();
        eprintln!("EVENTMAP_EXTRA n={}", extra.len());
        for e in extra {
            eprintln!(
                "  EXTRA {} {}/{}",
                e.start_1based.get(),
                e.ref_allele,
                e.alt_allele
            );
        }
        assert!(has_ca, "EventMap 92317399 C/A");
        assert!(has_tc, "EventMap 92317407 T/C");
        assert!(has_gc, "EventMap 92317412 G/C");
        assert!(
            assemble_three || seq_kbest_three,
            "some extracted ALT haplotype must carry all three Java substitutions"
        );

        let prod_assembled =
            assemble_from_ref_and_reads(&padded, &assembly_reads, &args).expect("assemble padded");
        let prod_events = collect_variation_events(
            &prod_assembled.haplotypes,
            &padded.bases,
            pad_start,
            &region.contig,
            0,
        );
        let prod_has_ca = prod_events
            .iter()
            .any(|e| e.start_1based.get() == SITE_CA && e.ref_allele == "C" && e.alt_allele == "A");
        eprintln!(
            "CONTROL_A production_padded assemble status={:?} kmer={} n_haps={} event_C/A={prod_has_ca}",
            prod_assembled.status, prod_assembled.kmer_size, prod_assembled.haplotypes.len()
        );
        assert!(
            !prod_has_ca,
            "production padded-REF assembly must not be the source of the Java C/A in this harness"
        );
    }
}
