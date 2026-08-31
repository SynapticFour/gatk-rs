//! 6R.29 TEST-ONLY: extra SNP emission differential (92317361 C/T, 92317371 C/G).
//! Does not change production algorithms, EventMap, P12, W-H1, k-mer policy, or the 6R.28 keep-mask.

#[cfg(test)]
mod traces {
    use crate::allele_filter_options::AlleleFilterOptions;
    use crate::allele_filtering::filter_assembly_and_likelihoods;
    use crate::assembly_based_caller::{assemble_reads_with_finalized, AssembleReadsArgs};
    use crate::assembly_region_finalize::assembly_reference_read;
    use crate::assembly_region_trimmer::{
        AssemblyRegionTrimmer, AssemblyRegionTrimmerConfig, TrimVariant,
    };
    use crate::assembly_result_set::AssemblyResultSet;
    use crate::engine::{take_call_region_audit, CallRegionArgs, HaplotypeCallerEngine};
    use crate::event_map::{variation_events_for_haplotype, VariationEvent};
    use crate::genome_loc::GenomePosition;
    use crate::haplotype::Haplotype;
    use crate::hc_allele_mapping::hap_base_at_ref_locus;
    use crate::hc_genotyping_engine::HcGenotypingConfig;
    use crate::java_hc_site_semantics::is_strict_java_production_emit_candidate;
    use crate::read_event_discovery::{
        is_strict_java_production_emit_admits, p12_baseline_emit_oracle_blocks,
    };
    use crate::read_threading_assembler::{AssemblyStatus, DEFAULT_NUM_BEST_HAPLOTYPES_PER_GRAPH};
    use crate::region_vcf_emit::try_emit_call_region_variants_with_config;
    use crate::walker_apply::{call_disposition, AssemblyRegionCallDisposition};
    use gatk_core::reference::{
        parse_intervals_cli_string, ReferenceWindowCache, SequenceDictionary,
    };
    use std::collections::hash_map::DefaultHasher;
    use std::hash::{Hash, Hasher};
    use std::path::Path;

    const SITE_CT: u64 = 92_317_361;
    const SITE_CG: u64 = 92_317_371;
    const SITE_CA: u64 = 92_317_399;
    const SITE_TC: u64 = 92_317_407;
    const SITE_GC: u64 = 92_317_412;
    const JAVA_ACTIVE: (u64, u64) = (92_317_262, 92_317_491);
    const JAVA_EXTENDED: (u64, u64) = (92_317_162, 92_317_591);
    /// OBSERVED JAVA EXECUTION (6R.24–6R.28 Docker): trimmed haplotypes 54M,
    /// span `2:92317379-92317432`.
    const JAVA_TRIM: (u64, u64) = (92_317_379, 92_317_432);
    const SNP_PAD: u64 = 20;

    #[derive(Clone, Copy, Debug, PartialEq, Eq)]
    enum HapKind {
        Ref,
        AltA,
        AltB,
        AltC,
        Other,
    }

    impl HapKind {
        fn as_str(self) -> &'static str {
            match self {
                HapKind::Ref => "REF",
                HapKind::AltA => "ALT-A",
                HapKind::AltB => "ALT-B",
                HapKind::AltC => "ALT-C",
                HapKind::Other => "ALT-OTHER",
            }
        }
    }

    fn fixture_paths() -> Option<(std::path::PathBuf, std::path::PathBuf)> {
        let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("..");
        let ref_path = root.join("parity/realworld/assets/hs37d5.simple.fa");
        let bam = root.join("parity/realworld/na12878_20k_b37/NA12878_20k.b37.bam");
        if !ref_path.is_file() || !bam.is_file() {
            return None;
        }
        Some((ref_path, bam))
    }

    fn digest(bytes: &[u8]) -> u64 {
        let mut h = DefaultHasher::new();
        bytes.hash(&mut h);
        h.finish()
    }

    fn snp_event(pos: u64, r: &str, a: &str) -> VariationEvent {
        VariationEvent {
            contig: "2".into(),
            start_1based: GenomePosition::new_1based(pos),
            end_1based: GenomePosition::new_1based(pos),
            ref_allele: r.into(),
            alt_allele: a.into(),
        }
    }

    fn event_has(events: &[VariationEvent], start: u64, r: &str, a: &str) -> bool {
        events
            .iter()
            .any(|e| e.start_1based.get() == start && e.ref_allele == r && e.alt_allele == a)
    }

    fn hap_has_snp(h: &Haplotype, pad: u64, site: u64, alt: u8) -> bool {
        hap_base_at_ref_locus(h, pad, site) == Some(alt)
    }

    fn classify_hap(h: &Haplotype, pad: u64) -> HapKind {
        if h.is_reference {
            return HapKind::Ref;
        }
        let ct = hap_has_snp(h, pad, SITE_CT, b'T');
        let cg = hap_has_snp(h, pad, SITE_CG, b'G');
        let ca = hap_has_snp(h, pad, SITE_CA, b'A');
        if ca && !ct && !cg {
            HapKind::AltA
        } else if ca && (ct || cg) {
            HapKind::AltB
        } else if !ca && (ct || cg) {
            HapKind::AltC
        } else {
            HapKind::Other
        }
    }

    fn n_alt(a: &AssemblyResultSet) -> usize {
        a.haplotypes.iter().filter(|h| !h.is_reference).count()
    }

    fn hap_ll_sums(
        n_haps: usize,
        rows: &[crate::region_read_likelihood::RegionReadLikelihood],
    ) -> Vec<f64> {
        let mut sums = vec![0.0; n_haps];
        for row in rows {
            let i = row.haplotype_index.get();
            if i < n_haps {
                sums[i] += row.log10_likelihood;
            }
        }
        sums
    }

    fn load_mid_b_region() -> Option<(
        crate::assembly_region_iterator::AssemblyRegion,
        SequenceDictionary,
        std::path::PathBuf,
    )> {
        let (ref_fasta, bam_path) = fixture_paths()?;
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
        Some((region, dict, ref_fasta))
    }

    fn assemble_untrimmed(
        region: &crate::assembly_region_iterator::AssemblyRegion,
        dict: &SequenceDictionary,
        ref_fasta: &std::path::Path,
    ) -> AssemblyResultSet {
        let mut owned = region.clone();
        let mut assemble_args = AssembleReadsArgs::default();
        assemble_args.strict_java_assembly = true;
        let mut ref_cache = ReferenceWindowCache::new(ref_fasta.to_path_buf(), 4);
        let _padded = assembly_reference_read(dict, &mut ref_cache, region).expect("pad");
        assemble_reads_with_finalized(&mut owned, dict, &mut ref_cache, &assemble_args)
            .expect("production assemble")
            .assembly
    }

    fn subset_kinds(src: &AssemblyResultSet, kinds: &[HapKind]) -> AssemblyResultSet {
        let pad = src.padded_reference_start_1based();
        let haps: Vec<Haplotype> = src
            .haplotypes
            .iter()
            .filter(|h| kinds.contains(&classify_hap(h, pad)))
            .cloned()
            .collect();
        assert!(
            haps.iter().any(|h| h.is_reference),
            "subset must retain REF"
        );
        AssemblyResultSet::from_assembly_for_calling_owned(
            AssemblyStatus::AssembledSomeVariation,
            src.kmer_size_for_dump(),
            haps,
            src.reference_bases().to_vec(),
            pad,
            "2",
            src.max_mnp_distance(),
        )
    }

    fn dump_hap_trace(label: &str, assembly: &AssemblyResultSet) {
        let pad = assembly.padded_reference_start_1based();
        let ref_hap = assembly
            .haplotypes
            .iter()
            .find(|h| h.is_reference)
            .expect("REF");
        eprintln!(
            "{label} n_haps={} n_alt={} n_events={} pad={pad} ref_len={}",
            assembly.haplotypes.len(),
            n_alt(assembly),
            assembly.variation_events.len(),
            assembly.reference_bases().len()
        );
        for (i, h) in assembly.haplotypes.iter().enumerate() {
            let kind = classify_hap(h, pad);
            let ct = hap_has_snp(h, pad, SITE_CT, b'T');
            let cg = hap_has_snp(h, pad, SITE_CG, b'G');
            let ca = hap_has_snp(h, pad, SITE_CA, b'A');
            let tc = hap_has_snp(h, pad, SITE_TC, b'C');
            let gc = hap_has_snp(h, pad, SITE_GC, b'C');
            eprintln!(
                "  HAP[{i}] {} ref={} len={} digest={:016x} 361C/T={ct} 371C/G={cg} 399C/A={ca} 407T/C={tc} 412G/C={gc}",
                kind.as_str(),
                h.is_reference,
                h.bases.len(),
                digest(&h.bases)
            );
            let emap = variation_events_for_haplotype(
                h,
                ref_hap,
                assembly.reference_bases(),
                pad,
                assembly.max_mnp_distance(),
                "2",
            );
            for e in &emap {
                if e.start_1based.get() >= SITE_CT && e.start_1based.get() <= SITE_GC {
                    eprintln!(
                        "    EVENT {} {}/{}",
                        e.start_1based.get(),
                        e.ref_allele,
                        e.alt_allele
                    );
                }
            }
        }
        for e in assembly.variation_events() {
            if e.start_1based.get() == SITE_CT
                || e.start_1based.get() == SITE_CG
                || e.start_1based.get() == SITE_CA
                || e.start_1based.get() == SITE_TC
                || e.start_1based.get() == SITE_GC
            {
                eprintln!(
                    "  ASSEMBLY_EVENT {} {}/{} emit_candidate={} emit_admits={} oracle_block={}",
                    e.start_1based.get(),
                    e.ref_allele,
                    e.alt_allele,
                    is_strict_java_production_emit_candidate(e),
                    is_strict_java_production_emit_admits(e),
                    p12_baseline_emit_oracle_blocks(e)
                );
            }
        }
    }

    fn pos_in(span: (u64, u64), pos: u64) -> bool {
        pos >= span.0 && pos <= span.1
    }

    fn java_style_padded_span(
        events: &[VariationEvent],
        region: &crate::assembly_region_iterator::AssemblyRegion,
    ) -> (u64, u64) {
        let overlapping: Vec<_> = events
            .iter()
            .filter(|e| {
                e.start_1based.get() <= region.end.get() && e.end_1based.get() >= region.start.get()
            })
            .collect();
        if overlapping.is_empty() {
            return (region.start.get(), region.end.get());
        }
        let mut min_start = overlapping
            .iter()
            .map(|e| e.start_1based.get())
            .min()
            .unwrap();
        let mut max_end = overlapping
            .iter()
            .map(|e| e.end_1based.get())
            .max()
            .unwrap();
        for e in &overlapping {
            let padding = if e.is_indel() { 75 } else { SNP_PAD };
            min_start = min_start.min(e.start_1based.get().saturating_sub(padding).max(1));
            max_end = max_end.max(e.end_1based.get().saturating_add(padding));
        }
        let padded_start = region.extended_start.get().max(min_start);
        let padded_end = region.extended_end.get().min(max_end);
        (padded_start, padded_end)
    }

    #[test]
    fn six_r29_java_kbest_source_is_not_keep_two() {
        assert_eq!(
            DEFAULT_NUM_BEST_HAPLOTYPES_PER_GRAPH, 128,
            "Rust k-best cap must remain GATK numBestHaplotypesPerGraph default 128"
        );
        let assembler = include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/src/read_threading_assembler.rs"
        ));
        assert!(
            assembler.contains("num_best_haplotypes_per_graph"),
            "findBestPaths must remain k-best capped, not a 2-haplotype rule"
        );
        assert!(
            !assembler.contains("92317361"),
            "must not special-case 361 in the assembler"
        );
        assert!(
            !assembler.contains("92317371"),
            "must not special-case 371 in the assembler"
        );
        let filtering = include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/src/allele_filtering.rs"
        ));
        assert!(
            filtering.contains("legacy_unique_snp_rank_filter_assembly"),
            "6R.28 keep-mask must remain (legacy unique-supporter is test-only)"
        );
    }

    #[test]
    fn six_r29_emit_band_excludes_361_371() {
        let ct = snp_event(SITE_CT, "C", "T");
        let cg = snp_event(SITE_CG, "C", "G");
        let ca = snp_event(SITE_CA, "C", "A");
        let tc = snp_event(SITE_TC, "T", "C");
        let gc = snp_event(SITE_GC, "G", "C");
        assert!(
            !is_strict_java_production_emit_candidate(&ct),
            "361 C/T is before MID_B_DENSE_CLUSTER_START=92317399"
        );
        assert!(!is_strict_java_production_emit_candidate(&cg));
        assert!(!is_strict_java_production_emit_admits(&ct));
        assert!(!is_strict_java_production_emit_admits(&cg));
        assert!(!p12_baseline_emit_oracle_blocks(&ct));
        assert!(!p12_baseline_emit_oracle_blocks(&cg));
        assert!(is_strict_java_production_emit_candidate(&ca));
        assert!(is_strict_java_production_emit_candidate(&tc));
        assert!(is_strict_java_production_emit_candidate(&gc));
        eprintln!(
            "EMIT_BAND 361 C/T candidate={} admits={} | 371 C/G candidate={} admits={} | oracle 399/407/412 all candidate",
            is_strict_java_production_emit_candidate(&ct),
            is_strict_java_production_emit_admits(&ct),
            is_strict_java_production_emit_candidate(&cg),
            is_strict_java_production_emit_admits(&cg)
        );
        eprintln!(
            "SPAN 361 in_active={} in_java_trim={} | 371 in_active={} in_java_trim={}",
            pos_in(JAVA_ACTIVE, SITE_CT),
            pos_in(JAVA_TRIM, SITE_CT),
            pos_in(JAVA_ACTIVE, SITE_CG),
            pos_in(JAVA_TRIM, SITE_CG)
        );
        assert!(pos_in(JAVA_ACTIVE, SITE_CT) && pos_in(JAVA_ACTIVE, SITE_CG));
        assert!(!pos_in(JAVA_TRIM, SITE_CT) && !pos_in(JAVA_TRIM, SITE_CG));
        assert!(
            pos_in(JAVA_TRIM, SITE_CA) && pos_in(JAVA_TRIM, SITE_TC) && pos_in(JAVA_TRIM, SITE_GC)
        );
    }

    #[test]
    fn six_r29_control_a_canonical_mid_b_trace() {
        let Some((region, dict, ref_fasta)) = load_mid_b_region() else {
            eprintln!("Real-data mid-B comparison unavailable");
            return;
        };
        assert_eq!((region.start.get(), region.end.get()), JAVA_ACTIVE);
        assert_eq!(
            (region.extended_start.get(), region.extended_end.get()),
            JAVA_EXTENDED
        );
        assert_eq!(region.reads.len(), 2);

        let untrimmed = assemble_untrimmed(&region, &dict, &ref_fasta);
        dump_hap_trace("6R.29_CONTROL_A_ASSEMBLY", &untrimmed);
        assert_eq!(untrimmed.haplotypes.len(), 2);
        assert_eq!(n_alt(&untrimmed), 1);
        assert!(untrimmed.haplotypes.iter().all(|h| h.bases.len() == 430));
        let pad = untrimmed.padded_reference_start_1based();
        let kinds: Vec<HapKind> = untrimmed
            .haplotypes
            .iter()
            .map(|h| classify_hap(h, pad))
            .collect();
        assert!(kinds.contains(&HapKind::Ref));
        assert!(kinds.contains(&HapKind::AltA));
        assert!(!kinds.contains(&HapKind::AltB));
        assert!(!kinds.contains(&HapKind::AltC));
        assert!(
            !event_has(untrimmed.variation_events(), SITE_CT, "C", "T"),
            "6R.33: 361 C/T is not assembled after Java-faithful prefix-match abort"
        );
        assert!(
            !event_has(untrimmed.variation_events(), SITE_CG, "C", "G"),
            "6R.33: 371 C/G is not assembled after Java-faithful prefix-match abort"
        );
        assert!(event_has(untrimmed.variation_events(), SITE_CA, "C", "A"));
        assert!(event_has(untrimmed.variation_events(), SITE_TC, "T", "C"));
        assert!(event_has(untrimmed.variation_events(), SITE_GC, "G", "C"));

        let args = CallRegionArgs::strict_java();
        let outcome =
            HaplotypeCallerEngine::call_region(&region, &dict, &ref_fasta, &args).expect("call");
        let audit = take_call_region_audit();
        let some = outcome
            .as_ref()
            .expect("call_region must return Some after 6R.28");
        dump_hap_trace("6R.29_CONTROL_A_AFTER_CALL_REGION", &some.assembly);

        let rust_trim = (
            audit.trim_padded_start.unwrap_or(0),
            audit.trim_padded_end.unwrap_or(0),
        );
        eprintln!(
            "TRIM rust_padded={:?}-{:?} java_observed={}..{} 361_in_rust_trim={} 371_in_rust_trim={} 361_in_java_trim={} 371_in_java_trim={}",
            audit.trim_padded_start,
            audit.trim_padded_end,
            JAVA_TRIM.0,
            JAVA_TRIM.1,
            pos_in(rust_trim, SITE_CT),
            pos_in(rust_trim, SITE_CG),
            pos_in(JAVA_TRIM, SITE_CT),
            pos_in(JAVA_TRIM, SITE_CG)
        );

        let ll = hap_ll_sums(some.assembly.haplotypes.len(), &some.read_likelihoods);
        let out_pad = some.assembly.padded_reference_start_1based();
        eprintln!(
            "PAIRHMM n_rows={} n_haps={}",
            some.read_likelihoods.len(),
            some.assembly.haplotypes.len()
        );
        for (i, h) in some.assembly.haplotypes.iter().enumerate() {
            let kind = classify_hap(h, out_pad);
            let ct = hap_has_snp(h, out_pad, SITE_CT, b'T');
            let cg = hap_has_snp(h, out_pad, SITE_CG, b'G');
            eprintln!(
                "  LL_HAP[{i}] {} digest={:016x} ll_sum={:.6} 361C/T={ct} 371C/G={cg} kept=true",
                kind.as_str(),
                digest(&h.bases),
                ll.get(i).copied().unwrap_or(f64::NAN)
            );
        }

        let post_events = some.assembly.variation_events();
        let in_gt = |pos: u64, r: &str, a: &str| {
            some.genotyped_calls.iter().any(|c| {
                c.event.start_1based.get() == pos
                    && c.event.ref_allele == r
                    && c.event.alt_allele == a
            })
        };
        eprintln!("GENOTYPED_CALLS n={}", some.genotyped_calls.len());
        for c in &some.genotyped_calls {
            let e = &c.event;
            if e.start_1based.get() == SITE_CT
                || e.start_1based.get() == SITE_CG
                || e.start_1based.get() == SITE_CA
                || e.start_1based.get() == SITE_TC
                || e.start_1based.get() == SITE_GC
            {
                eprintln!(
                    "  GT {} {}/{} PL={:?} GQ={} AD={:?} DP={} emit_candidate={}",
                    e.start_1based.get(),
                    e.ref_allele,
                    e.alt_allele,
                    c.genotype.format.pl_as_i32(),
                    c.genotype.format.gq.as_i32(),
                    c.genotype.format.ad_as_i32(),
                    c.genotype.format.dp.as_i32(),
                    is_strict_java_production_emit_candidate(e)
                );
            }
        }

        let vcf = try_emit_call_region_variants_with_config(
            &region,
            some,
            "NA12878",
            HcGenotypingConfig::strict_java().stand_emit_confidence,
            &HcGenotypingConfig::strict_java(),
        )
        .expect("vcf emit");
        let vcf_has = |pos: u64, r: &str, a: &str| {
            vcf.iter().any(|rec| {
                rec.position == pos
                    && rec.reference == r
                    && rec.alternate.iter().any(|alt| alt == a)
            })
        };
        eprintln!("VCF n={}", vcf.len());
        for rec in &vcf {
            if rec.position == SITE_CT
                || rec.position == SITE_CG
                || rec.position == SITE_CA
                || rec.position == SITE_TC
                || rec.position == SITE_GC
            {
                eprintln!(
                    "  VCF {} {}/{} QUAL={:?}",
                    rec.position,
                    rec.reference,
                    rec.alternate.join(","),
                    rec.quality
                );
            }
        }

        let row = |pos: u64, r: &str, a: &str, hap_present: bool| {
            let emap = event_has(untrimmed.variation_events(), pos, r, a);
            let trimmed_hap = some
                .assembly
                .haplotypes
                .iter()
                .any(|h| hap_base_at_ref_locus(h, out_pad, pos) == Some(a.as_bytes()[0]));
            let post_emap = event_has(post_events, pos, r, a);
            let gt = in_gt(pos, r, a);
            let emitted = vcf_has(pos, r, a);
            eprintln!(
                "TRACE {pos} {r}/{a} assembly_hap={} eventmap={} trimmed_hap={} post_call_eventmap={} genotyped={} vcf={} emit_candidate={}",
                hap_present,
                emap,
                trimmed_hap,
                post_emap,
                gt,
                emitted,
                is_strict_java_production_emit_candidate(&snp_event(pos, r, a))
            );
            (hap_present, emap, trimmed_hap, post_emap, gt, emitted)
        };

        let ct_hap = untrimmed
            .haplotypes
            .iter()
            .any(|h| hap_has_snp(h, pad, SITE_CT, b'T'));
        let cg_hap = untrimmed
            .haplotypes
            .iter()
            .any(|h| hap_has_snp(h, pad, SITE_CG, b'G'));
        let (ct_asm, ct_em, ct_trim, ct_post, ct_gt, ct_vcf) = row(SITE_CT, "C", "T", ct_hap);
        let (cg_asm, cg_em, cg_trim, cg_post, cg_gt, cg_vcf) = row(SITE_CG, "C", "G", cg_hap);
        assert!(
            !ct_asm && !cg_asm,
            "6R.33: 361/371 extra SNPs are not on assembled haplotypes"
        );
        assert!(
            !ct_em && !cg_em,
            "6R.33: untrimmed EventMap must not contain 361/371"
        );
        assert!(
            !ct_vcf && !cg_vcf,
            "production VCF must not currently emit 361/371 (observe, do not suppress haplotypes)"
        );
        eprintln!(
            "361_STAGE hap={ct_asm} emap={ct_em} trimmed_hap={ct_trim} post_emap={ct_post} gt={ct_gt} vcf={ct_vcf}"
        );
        eprintln!(
            "371_STAGE hap={cg_asm} emap={cg_em} trimmed_hap={cg_trim} post_emap={cg_post} gt={cg_gt} vcf={cg_vcf}"
        );
        assert!(event_has(post_events, SITE_CA, "C", "A"));
        assert!(event_has(post_events, SITE_TC, "T", "C"));
        assert!(event_has(post_events, SITE_GC, "G", "C"));
        assert!(
            n_alt(&some.assembly) >= 1,
            "production must retain assembled extra haplotypes (6R.28 keep-mask)"
        );
        let _ = (ct_trim, ct_post, ct_gt, cg_trim, cg_post, cg_gt);
    }

    #[test]
    fn six_r29_control_b_oracle_three_snp_haplotype_only() {
        let Some((region, dict, ref_fasta)) = load_mid_b_region() else {
            eprintln!("Real-data mid-B comparison unavailable");
            return;
        };
        let untrimmed = assemble_untrimmed(&region, &dict, &ref_fasta);
        let oracle = subset_kinds(&untrimmed, &[HapKind::Ref, HapKind::AltA]);
        dump_hap_trace("6R.29_CONTROL_B_ORACLE_ONLY", &oracle);
        assert_eq!(oracle.haplotypes.len(), 2);
        assert_eq!(n_alt(&oracle), 1);
        assert!(event_has(oracle.variation_events(), SITE_CA, "C", "A"));
        assert!(event_has(oracle.variation_events(), SITE_TC, "T", "C"));
        assert!(event_has(oracle.variation_events(), SITE_GC, "G", "C"));
        assert!(!event_has(oracle.variation_events(), SITE_CT, "C", "T"));
        assert!(!event_has(oracle.variation_events(), SITE_CG, "C", "G"));

        let mut filtered = oracle.clone();
        let _ = filter_assembly_and_likelihoods(
            &mut filtered,
            Vec::new(),
            AlleleFilterOptions::from_strict_java(true, Some(JAVA_ACTIVE.0), Some(JAVA_ACTIVE.1)),
        )
        .expect("filter");
        dump_hap_trace("6R.29_CONTROL_B_AFTER_FILTER", &filtered);
        assert_eq!(
            n_alt(&filtered),
            1,
            "oracle-only ALT-A must survive keep-mask"
        );
        assert!(event_has(filtered.variation_events(), SITE_CA, "C", "A"));
        let ca = snp_event(SITE_CA, "C", "A");
        assert!(
            is_strict_java_production_emit_candidate(&ca),
            "oracle C/A is an emit-band candidate"
        );
    }

    #[test]
    fn six_r29_control_c_extra_snp_haplotypes() {
        let Some((region, dict, ref_fasta)) = load_mid_b_region() else {
            eprintln!("Real-data mid-B comparison unavailable");
            return;
        };
        let untrimmed = assemble_untrimmed(&region, &dict, &ref_fasta);
        dump_hap_trace("6R.29_CONTROL_C_UNTRIMMED", &untrimmed);
        assert_eq!(untrimmed.haplotypes.len(), 2);
        assert_eq!(n_alt(&untrimmed), 1);
        let alt_b = subset_kinds(&untrimmed, &[HapKind::Ref, HapKind::AltB]);
        dump_hap_trace("6R.29_CONTROL_C_ALT_B", &alt_b);
        assert_eq!(n_alt(&alt_b), 0, "6R.33: ALT-B (361/371) is not assembled");
        let alt_c = subset_kinds(&untrimmed, &[HapKind::Ref, HapKind::AltC]);
        dump_hap_trace("6R.29_CONTROL_C_ALT_C", &alt_c);
        assert_eq!(n_alt(&alt_c), 0, "6R.33: ALT-C is not assembled");
        assert!(!is_strict_java_production_emit_candidate(&snp_event(
            SITE_CT, "C", "T"
        )));
        assert!(!is_strict_java_production_emit_candidate(&snp_event(
            SITE_CG, "C", "G"
        )));
        assert!(event_has(untrimmed.variation_events(), SITE_CA, "C", "A"));
    }

    #[test]
    fn six_r29_control_d_emission_span() {
        let Some((region, dict, ref_fasta)) = load_mid_b_region() else {
            eprintln!("Real-data mid-B comparison unavailable");
            return;
        };
        let untrimmed = assemble_untrimmed(&region, &dict, &ref_fasta);
        let trimmer =
            AssemblyRegionTrimmer::new(AssemblyRegionTrimmerConfig::gatk_defaults(), &dict, "2");
        let all_vars: Vec<TrimVariant> = untrimmed
            .variation_events()
            .iter()
            .map(|e| TrimVariant {
                contig: e.contig.clone(),
                start: e.start_1based.get(),
                end: e.end_1based.get(),
                is_indel: e.is_indel(),
            })
            .collect();
        let rust_trim = trimmer.trim(&region, &all_vars, Some(&region.reference));
        let oracle_vars: Vec<TrimVariant> = untrimmed
            .variation_events()
            .iter()
            .filter(|e| {
                let s = e.start_1based.get();
                s == SITE_CA || s == SITE_TC || s == SITE_GC
            })
            .map(|e| TrimVariant {
                contig: e.contig.clone(),
                start: e.start_1based.get(),
                end: e.end_1based.get(),
                is_indel: e.is_indel(),
            })
            .collect();
        let oracle_trim = trimmer.trim(&region, &oracle_vars, Some(&region.reference));
        let java_style_all = java_style_padded_span(untrimmed.variation_events(), &region);
        let oracle_events: Vec<VariationEvent> = untrimmed
            .variation_events()
            .iter()
            .filter(|e| {
                let s = e.start_1based.get();
                s == SITE_CA || s == SITE_TC || s == SITE_GC
            })
            .cloned()
            .collect();
        let java_style_oracle = java_style_padded_span(&oracle_events, &region);

        let rust_span = (
            rust_trim.padded_variant_start.unwrap_or(0),
            rust_trim.padded_variant_end.unwrap_or(0),
        );
        let oracle_span = (
            oracle_trim.padded_variant_start.unwrap_or(0),
            oracle_trim.padded_variant_end.unwrap_or(0),
        );
        eprintln!(
            "CONTROL_D rust_trim={:?}..{:?} oracle_only_trim={:?}..{:?} java_style_all={:?} java_style_oracle={:?} java_observed_trim={:?}",
            rust_trim.padded_variant_start,
            rust_trim.padded_variant_end,
            oracle_trim.padded_variant_start,
            oracle_trim.padded_variant_end,
            java_style_all,
            java_style_oracle,
            JAVA_TRIM
        );
        eprintln!(
            "CONTROL_D 361 in_rust={} in_oracle_trim={} in_java_style_all={} in_java_style_oracle={} in_java_observed={}",
            pos_in(rust_span, SITE_CT),
            pos_in(oracle_span, SITE_CT),
            pos_in(java_style_all, SITE_CT),
            pos_in(java_style_oracle, SITE_CT),
            pos_in(JAVA_TRIM, SITE_CT)
        );
        eprintln!(
            "CONTROL_D 371 in_rust={} in_oracle_trim={} in_java_style_all={} in_java_style_oracle={} in_java_observed={}",
            pos_in(rust_span, SITE_CG),
            pos_in(oracle_span, SITE_CG),
            pos_in(java_style_all, SITE_CG),
            pos_in(java_style_oracle, SITE_CG),
            pos_in(JAVA_TRIM, SITE_CG)
        );
        assert!(
            pos_in(JAVA_ACTIVE, SITE_CT) && pos_in(JAVA_ACTIVE, SITE_CG),
            "361/371 lie inside the active region"
        );
        assert!(
            !pos_in(JAVA_TRIM, SITE_CT) && !pos_in(JAVA_TRIM, SITE_CG),
            "361/371 lie outside Java's observed 54M trim"
        );
        assert!(
            pos_in(oracle_span, SITE_CA),
            "oracle-only trim still covers 399"
        );
        // Diagnostic only: extra SNPs in the EventMap enlarge the Rust pad relative to
        // Java's oracle-only window. Do not treat this as a production contract.
        eprintln!(
            "CONTROL_D span_only_would_exclude_361_from_java_trim=true rust_emission_still_gated_by_emit_band={}",
            !is_strict_java_production_emit_candidate(&snp_event(SITE_CT, "C", "T"))
        );
    }
}
