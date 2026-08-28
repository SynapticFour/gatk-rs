//! 6R.26: production `createGraph` uniqueness + `addSequence("ref")` use the
//! Java-equivalent extended-region haplotype (padding 0), not ±500 padded REF.
//! Does not change EventMap, P12 injection, mapper, k schedule, or uniqueness policy.

#[cfg(test)]
mod traces {
    use crate::assembly::AssemblyRead;
    use crate::assembly_based_caller::{assemble_reads_with_finalized, AssembleReadsArgs};
    use crate::assembly_region_finalize::{
        assembly_reference_read, create_graph_reference_read, padded_reference_loc,
        reference_haplotype_for_assembly_region, GATK_REFERENCE_PADDING_FOR_ASSEMBLY,
    };
    use crate::haplotype::Haplotype;
    use crate::read_event_discovery::P12_CLUSTER_TTC_START;
    use crate::read_threading_assembler::{
        build_threading_graph_for_seq_assembly, ReadThreadingAssemblerArgs,
    };
    use crate::read_threading_graph::reference_has_non_unique_kmers;
    use crate::walker_apply::{call_disposition, AssemblyRegionCallDisposition};
    use gatk_core::reference::{
        parse_intervals_cli_string, ReferenceWindowCache, SequenceDictionary,
    };
    use std::collections::hash_map::DefaultHasher;
    use std::hash::{Hash, Hasher};
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

    fn digest(bytes: &[u8]) -> u64 {
        let mut h = DefaultHasher::new();
        bytes.hash(&mut h);
        h.finish()
    }

    fn cigar_str(h: &Haplotype) -> String {
        h.cigar
            .as_ref()
            .map(|c| c.to_gatk_string())
            .unwrap_or_else(|| "NA".into())
    }

    fn event_has(
        events: &[crate::event_map::VariationEvent],
        start: u64,
        r: &str,
        a: &str,
    ) -> bool {
        events
            .iter()
            .any(|e| e.start_1based.get() == start && e.ref_allele == r && e.alt_allele == a)
    }

    #[test]
    fn six_r26_java_reference_gate_production() {
        assert!(
            !ReadThreadingAssemblerArgs::default().allow_non_unique_kmers_in_ref,
            "allowNonUniqueKmersInRef must remain false"
        );
        assert!(
            !AssembleReadsArgs::default()
                .assembler
                .allow_non_unique_kmers_in_ref,
            "production AssembleReadsArgs must not lift allowNonUniqueKmersInRef"
        );

        let assemble_src = include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/src/assembly_based_caller.rs"
        ));
        assert!(
            assemble_src.contains("assemble_from_ref_and_reads(&graph_ref, &reads, &assembler)?"),
            "uniqueness and addSequence(ref) must share graph_ref"
        );
        assert!(
            assemble_src.contains("create_graph_reference_read(&reference, region, dictionary)"),
            "graph_ref must be the padding-0 assembly-region haplotype"
        );
        assert!(
            !assemble_src.contains("assemble_from_ref_and_reads(&reference, &reads, &assembler)?"),
            "must not seed createGraph from the ±500 padded assembly_reference_read"
        );
        let rt_src = include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/src/read_threading_assembler.rs"
        ));
        assert!(
            rt_src.contains("&& reference_has_non_unique_kmers(reference, kmer_size)"),
            "uniqueness checking must still exist inside createGraph"
        );
        assert!(
            rt_src.contains("allow_non_unique_kmers_in_ref: false"),
            "default allow_non_unique_kmers_in_ref remains false"
        );

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

        assert_eq!(region.start.get(), JAVA_ACTIVE.0);
        assert_eq!(region.end.get(), JAVA_ACTIVE.1);
        assert_eq!(region.extended_start.get(), JAVA_EXTENDED.0);
        assert_eq!(region.extended_end.get(), JAVA_EXTENDED.1);
        assert_eq!(region.reads.len(), 2);
        assert!(
            !crate::read_threading_assembler::region_overlaps_p12_cluster(
                region.start.get(),
                region.end.get(),
            ),
            "canonical mid-B must not overlap the P12 cluster; C/A must not come from P12 inject"
        );
        assert!(
            region.end.get() < P12_CLUSTER_TTC_START
                || region.start.get() > P12_CLUSTER_TTC_START.saturating_add(3),
            "mid-B is not the TTC/ATG cluster"
        );

        let mut ref_cache = ReferenceWindowCache::new(ref_fasta.clone(), 4);
        let padded = assembly_reference_read(&dict, &mut ref_cache, &region).expect("pad");
        let (pad_start, pad_end) = padded_reference_loc(&region, &dict);
        let java_hap = reference_haplotype_for_assembly_region(&padded, &region, &dict);
        let graph_ref = create_graph_reference_read(&padded, &region, &dict);
        let ext_start = region.extended_start.get();
        let ext_end = region.extended_end.get();
        let expected_ext_len = (ext_end - ext_start + 1) as usize;

        assert_eq!(graph_ref.bases.as_slice(), java_hap.bases.as_slice());
        assert_eq!(graph_ref.bases.len(), expected_ext_len);
        assert_eq!(graph_ref.bases.len(), java_hap.bases.len());
        assert_eq!(
            java_hap.alignment_start_hap_wrt_ref,
            (ext_start - pad_start) as usize
        );
        eprintln!("PINNED_JAVA GATK 4.4.0.0 SHA=2dbc025821bc5f686c423ff332a41e6cef892a77");
        eprintln!(
            "IDENTITY uniqueness_and_seed_bytes len={} start={} end={} digest={:016x} pad_offset={}",
            graph_ref.bases.len(),
            ext_start,
            ext_end,
            digest(&graph_ref.bases),
            java_hap.alignment_start_hap_wrt_ref
        );
        eprintln!(
            "PADDED_DIAGNOSTIC len={} start={pad_start} end={pad_end} digest={:016x} extra_pad={}",
            padded.bases.len(),
            digest(&padded.bases),
            GATK_REFERENCE_PADDING_FOR_ASSEMBLY
        );

        let graph_k10 = unique_kmers(&graph_ref.bases, K10);
        let graph_k25 = unique_kmers(&graph_ref.bases, K25);
        let padded_k10 = unique_kmers(&padded.bases, K10);
        let padded_k25 = unique_kmers(&padded.bases, K25);
        eprintln!("UNIQUE k=10 graph_ref={graph_k10} padded={padded_k10}");
        eprintln!("UNIQUE k=25 graph_ref={graph_k25} padded={padded_k25}");
        assert!(
            !graph_k10,
            "Control A: k=10 rejected on Java-equivalent REF"
        );
        assert!(graph_k25, "Control A: k=25 accepted on Java-equivalent REF");
        assert!(
            !padded_k25,
            "Control D: old 1430 bp representation still rejects k=25 (diagnostic, not a second gate)"
        );

        let args = ReadThreadingAssemblerArgs::default();
        assert!(!args.allow_non_unique_kmers_in_ref);

        let mut dup = graph_ref.bases.clone();
        dup.extend_from_slice(&graph_ref.bases[..K25]);
        let dup_read = AssemblyRead {
            bases: dup.clone(),
            base_quals: vec![30; dup.len()],
        };
        assert!(
            !unique_kmers(&dup, K25),
            "Control B: extended REF with a duplicated 25-mer is non-unique"
        );
        let dup_built = build_threading_graph_for_seq_assembly(
            &dup_read,
            &[],
            K25,
            &args,
            args.allow_low_complexity_graphs,
            false,
        )
        .expect("dup k25")
        .is_some();
        assert!(
            !dup_built,
            "Control B: allowNonUniqueKmersInRef=false still rejects non-unique extended REF"
        );

        let mut other_gate_changed = 0usize;
        let mut other_n = 0usize;
        let mut stable_other: Option<crate::assembly_region_iterator::AssemblyRegion> = None;
        for r in &regions {
            if r.start.get() == region.start.get() && r.end.get() == region.end.get() {
                continue;
            }
            if !matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) {
                continue;
            }
            other_n += 1;
            let other_pad = assembly_reference_read(&dict, &mut ref_cache, r).expect("other pad");
            let other_graph = create_graph_reference_read(&other_pad, r, &dict);
            let changed = unique_kmers(&other_graph.bases, K10)
                != unique_kmers(&other_pad.bases, K10)
                || unique_kmers(&other_graph.bases, K25) != unique_kmers(&other_pad.bases, K25);
            eprintln!(
                "CONTROL_C other active={}-{} ext_len={} pad_len={} k10_graph={} k10_pad={} k25_graph={} k25_pad={} uniqueness_changed={changed}",
                r.start.get(),
                r.end.get(),
                other_graph.bases.len(),
                other_pad.bases.len(),
                unique_kmers(&other_graph.bases, K10),
                unique_kmers(&other_pad.bases, K10),
                unique_kmers(&other_graph.bases, K25),
                unique_kmers(&other_pad.bases, K25)
            );
            if changed {
                other_gate_changed += 1;
            } else if stable_other.is_none() {
                stable_other = Some(r.clone());
            }
        }
        eprintln!(
            "CONTROL_C n_other_active_full={other_n} uniqueness_gate_changed={other_gate_changed}"
        );
        if other_gate_changed > 0 {
            eprintln!(
                "CONTROL_C NOTE {other_gate_changed}/{other_n} sibling ActiveFull regions have a padded-vs-extended uniqueness delta (same representation class as mid-B, not a second production gate)"
            );
        }
        if let Some(mut stable) = stable_other {
            let assembled_stable = assemble_reads_with_finalized(
                &mut stable,
                &dict,
                &mut ref_cache,
                &AssembleReadsArgs::default(),
            )
            .expect("Control C production assemble on uniqueness-stable region");
            eprintln!(
                "CONTROL_C stable_region {}-{} n_haps={} kmer={:?} n_events={}",
                stable.start.get(),
                stable.end.get(),
                assembled_stable.assembly.haplotypes.len(),
                assembled_stable.assembly.minimum_kmer_size(),
                assembled_stable.assembly.variation_events.len()
            );
        }

        let mut owned = region.clone();
        let assembled = assemble_reads_with_finalized(
            &mut owned,
            &dict,
            &mut ref_cache,
            &AssembleReadsArgs::default(),
        )
        .expect("production assemble");
        let set = assembled.assembly;
        eprintln!(
            "PRODUCTION n_haps={} kmer={} n_events={} variation={} pad_start={} ref_bases_len={}",
            set.haplotypes.len(),
            set.minimum_kmer_size().unwrap_or(0),
            set.variation_events.len(),
            set.has_variation_for_calling(),
            set.padded_reference_start_1based(),
            set.reference_bases().len()
        );
        assert_eq!(
            set.minimum_kmer_size(),
            Some(K25),
            "production must accept k=25 (not stay on k=10 / expanded k)"
        );
        assert_eq!(
            set.padded_reference_start_1based(),
            ext_start,
            "EventMap loc must be the extended-span origin, not ±500 pad start"
        );
        assert_eq!(
            set.reference_bases(),
            graph_ref.bases.as_slice(),
            "EventMap reference bytes must equal createGraph uniqueness/seed bytes"
        );

        eprintln!(
            "COORD_ORIGIN ext={ext_start}-{ext_end} pad={pad_start}-{pad_end} eventmap_start={} hap_align_wrt_extended=0",
            set.padded_reference_start_1based()
        );
        let em_start = set.padded_reference_start_1based();
        let mut any_three = false;
        for (i, h) in set.haplotypes.iter().enumerate() {
            eprintln!(
                "HAP[{i}] ref={} cigar={} len={} align={} kmer={}",
                h.is_reference,
                cigar_str(h),
                h.bases.len(),
                h.alignment_start_hap_wrt_ref,
                h.kmer_size
            );
            if h.is_reference {
                continue;
            }
            let mut ok = true;
            for (site, expect_ref, expect_alt) in [
                (SITE_CA, b'C', b'A'),
                (SITE_TC, b'T', b'C'),
                (SITE_GC, b'G', b'C'),
            ] {
                let em_off = (site - em_start) as usize;
                let hap_off = em_off.saturating_sub(h.alignment_start_hap_wrt_ref);
                let hap_base = h.bases.get(hap_off).copied();
                let ref_base = set.reference_bases().get(em_off).copied();
                eprintln!(
                    "  COORD {site} hap_off={hap_off} em_off={em_off} align={} hap={:?} em_ref={:?}",
                    h.alignment_start_hap_wrt_ref,
                    hap_base.map(|b| b as char),
                    ref_base.map(|b| b as char)
                );
                if hap_base != Some(expect_alt) || ref_base != Some(expect_ref) {
                    ok = false;
                }
                if hap_base == Some(expect_alt) {
                    assert_eq!(
                        em_start + em_off as u64,
                        site,
                        "extended-span offset must reconstruct genomic coordinate {site}"
                    );
                }
            }
            if ok {
                any_three = true;
            }
        }
        assert!(
            any_three,
            "some production ALT haplotype must carry C/A T/C G/C at the Java sites"
        );

        let events = set.variation_events();
        for e in events {
            eprintln!(
                "EVENTMAP {} {}/{}",
                e.start_1based.get(),
                e.ref_allele,
                e.alt_allele
            );
        }
        assert!(
            event_has(events, SITE_CA, "C", "A"),
            "EventMap must contain 92317399 C/A"
        );
        assert!(
            event_has(events, SITE_TC, "T", "C"),
            "EventMap must contain 92317407 T/C"
        );
        assert!(
            event_has(events, SITE_GC, "G", "C"),
            "EventMap must contain 92317412 G/C"
        );
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

        let wh1 = include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/src/compatibility/mod.rs"
        ));
        assert!(wh1.contains("W-H1"), "W-H1 waiver module remains");
        let discovery = include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/src/read_event_discovery/mod.rs"
        ));
        assert!(
            discovery.contains("fn fix_p12_cluster_coupled_alt_haplotype")
                || discovery.contains("fix_p12_cluster_coupled_alt_haplotype"),
            "fix_p12 remains"
        );
        assert!(
            discovery.contains("ensure_p12_cluster_variation_events_for_active_span"),
            "ensure_p12 remains"
        );
        assert!(
            !assemble_src.contains("92317399"),
            "no P12/mid-B coordinate special case in assembleReads"
        );
        eprintln!(
            "P12_SAFEGUARD W-H1 unchanged (module present); fix_p12/ensure_p12 present; mid-B C/A from ordinary k=25 assembly"
        );
    }
}
