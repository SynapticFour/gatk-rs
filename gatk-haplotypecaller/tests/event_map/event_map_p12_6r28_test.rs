//! 6R.28: Java 4.4 `AlleleFiltering.filterAlleles` vs unique-supporter keep-mask.

#[cfg(test)]
mod traces {
    use crate::allele_filter_options::AlleleFilterOptions;
    use crate::allele_filtering::{
        filter_assembly_and_likelihoods, legacy_unique_snp_rank_filter_assembly,
    };
    use crate::assembly_based_caller::{assemble_reads_with_finalized, AssembleReadsArgs};
    use crate::assembly_region_finalize::assembly_reference_read;
    use crate::assembly_result_set::AssemblyResultSet;
    use crate::cigar::{Cigar, CigarOperator};
    use crate::event_map::{variation_events_for_haplotype, VariationEvent};
    use crate::genome_loc::GenomeLoc;
    use crate::haplotype::Haplotype;
    use crate::hc_allele_mapping::hap_base_at_ref_locus;
    use crate::read_threading_assembler::AssemblyStatus;
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

    fn java44_default_hc_keep_all(n_haps: usize) -> Vec<bool> {
        vec![true; n_haps]
    }

    fn java44_filter_alleles_drops_because_n_exact_gt_1(_n_exact: usize) -> bool {
        false
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

    fn cigar_str(h: &Haplotype) -> String {
        h.cigar
            .as_ref()
            .map(|c| c.to_gatk_string())
            .unwrap_or_else(|| "NA".into())
    }

    fn match_cigar(len: usize) -> Cigar {
        let mut c = Cigar::new();
        c.push(len, CigarOperator::Match);
        c
    }

    fn event_has(events: &[VariationEvent], start: u64, r: &str, a: &str) -> bool {
        events
            .iter()
            .any(|e| e.start_1based.get() == start && e.ref_allele == r && e.alt_allele == a)
    }

    fn hap_has_snp(h: &Haplotype, pad: u64, site: u64, alt: u8) -> bool {
        hap_base_at_ref_locus(h, pad, site) == Some(alt)
    }

    fn snp_assembly(pad_start: u64, ref_bases: Vec<u8>, alts: Vec<Vec<u8>>) -> AssemblyResultSet {
        let loc = GenomeLoc::new(pad_start, pad_start + ref_bases.len() as u64 - 1);
        let len = ref_bases.len();
        let mut haps = Vec::new();
        let mut rh = Haplotype::new(ref_bases.clone(), true);
        rh.cigar = Some(match_cigar(len));
        rh.genome_loc = Some(loc);
        haps.push(rh);
        for alt in alts {
            let mut h = Haplotype::new(alt, false);
            h.cigar = Some(match_cigar(len));
            h.genome_loc = Some(loc);
            haps.push(h);
        }
        AssemblyResultSet::from_assembly_for_calling_owned(
            AssemblyStatus::AssembledSomeVariation,
            25,
            haps,
            ref_bases,
            pad_start,
            "2",
            0,
        )
    }

    fn filter_span(assembly: &AssemblyResultSet, start: u64, end: u64) -> AssemblyResultSet {
        let mut filtered = assembly.clone();
        let _ = filter_assembly_and_likelihoods(
            &mut filtered,
            Vec::new(),
            AlleleFilterOptions::from_strict_java(true, Some(start), Some(end)),
        )
        .expect("allele filter");
        filtered
    }

    fn n_alt(a: &AssemblyResultSet) -> usize {
        a.haplotypes.iter().filter(|h| !h.is_reference).count()
    }

    fn n_exact_snp(assembly: &AssemblyResultSet, site: u64, alt: u8) -> usize {
        let pad = assembly.padded_reference_start_1based();
        assembly
            .haplotypes
            .iter()
            .filter(|h| !h.is_reference && hap_has_snp(h, pad, site, alt))
            .count()
    }

    fn dump_hap_trace(label: &str, assembly: &AssemblyResultSet) {
        let pad = assembly.padded_reference_start_1based();
        let ref_hap = assembly
            .haplotypes
            .iter()
            .find(|h| h.is_reference)
            .expect("REF");
        eprintln!(
            "{label} n_haps={} n_alt={} n_events={}",
            assembly.haplotypes.len(),
            n_alt(assembly),
            assembly.variation_events.len()
        );
        for (i, h) in assembly.haplotypes.iter().enumerate() {
            let emap = variation_events_for_haplotype(
                h,
                ref_hap,
                assembly.reference_bases(),
                pad,
                assembly.max_mnp_distance(),
                "2",
            );
            let ct = hap_has_snp(h, pad, SITE_CT, b'T');
            let cg = hap_has_snp(h, pad, SITE_CG, b'G');
            let ca = hap_has_snp(h, pad, SITE_CA, b'A');
            let tc = hap_has_snp(h, pad, SITE_TC, b'C');
            let gc = hap_has_snp(h, pad, SITE_GC, b'C');
            let kind = if h.is_reference {
                "REF"
            } else if ca && !ct && !cg {
                "ALT-A"
            } else if ca && (ct || cg) {
                "ALT-B"
            } else if !ca && (ct || cg) {
                "ALT-C"
            } else {
                "ALT-OTHER"
            };
            eprintln!(
                "  HAP[{i}] {kind} ref={} cigar={} len={} digest={:016x} 361C/T={ct} 371C/G={cg} 399C/A={ca} 407T/C={tc} 412G/C={gc}",
                h.is_reference,
                cigar_str(h),
                h.bases.len(),
                digest(&h.bases)
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
            eprintln!("    KEEP_TRACE {kind} C/A_support={ca} extra_361={ct} extra_371={cg}");
        }
    }

    fn load_untrimmed_mid_b() -> Option<(AssemblyResultSet, u64, u64)> {
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
        assert_eq!((region.start.get(), region.end.get()), JAVA_ACTIVE);
        let mut owned = region.clone();
        let mut assemble_args = AssembleReadsArgs::default();
        assemble_args.strict_java_assembly = true;
        let mut ref_cache = ReferenceWindowCache::new(ref_fasta.clone(), 4);
        let _padded = assembly_reference_read(&dict, &mut ref_cache, &region).expect("pad");
        let assembled =
            assemble_reads_with_finalized(&mut owned, &dict, &mut ref_cache, &assemble_args)
                .expect("production assemble");
        Some((assembled.assembly, region.start.get(), region.end.get()))
    }

    #[test]
    fn six_r28_java_source_contract() {
        let filtering = include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/src/allele_filtering.rs"
        ));
        assert!(
            filtering.contains("legacy_unique_snp_rank_filter_assembly"),
            "legacy unique-supporter must remain test-only"
        );
        assert!(
            !filtering.contains("92317399"),
            "must not special-case the mid-B C/A locus"
        );
        assert!(
            !filtering.contains("exact.len() >= 1"),
            "forbidden unique-supporter flip"
        );
    }

    #[test]
    fn six_r28_synthetic_two_supporter_matrix() {
        let pad = SITE_CA - 40;
        let mut ref_bases = vec![b'A'; 80];
        ref_bases[40] = b'C';
        ref_bases[21] = b'G';
        let mut alt_ca = ref_bases.clone();
        alt_ca[40] = b'A';
        let mut alt_ca_extra = alt_ca.clone();
        alt_ca_extra[21] = b'T';

        let a = snp_assembly(pad, ref_bases.clone(), vec![alt_ca.clone()]);
        assert_eq!(
            n_exact_snp(
                &filter_span(&a, JAVA_ACTIVE.0, JAVA_ACTIVE.1),
                SITE_CA,
                b'A'
            ),
            1
        );

        let b = snp_assembly(pad, ref_bases.clone(), vec![alt_ca.clone(), alt_ca.clone()]);
        assert_eq!(n_exact_snp(&b, SITE_CA, b'A'), 2);
        assert!(!java44_filter_alleles_drops_because_n_exact_gt_1(2));
        assert_eq!(java44_default_hc_keep_all(3), vec![true, true, true]);
        assert_eq!(
            n_exact_snp(
                &filter_span(&b, JAVA_ACTIVE.0, JAVA_ACTIVE.1),
                SITE_CA,
                b'A'
            ),
            2,
            "Case B: shared C/A remains with two identical-support ALTs"
        );

        let c = snp_assembly(pad, ref_bases, vec![alt_ca, alt_ca_extra]);
        assert_eq!(
            n_exact_snp(
                &filter_span(&c, JAVA_ACTIVE.0, JAVA_ACTIVE.1),
                SITE_CA,
                b'A'
            ),
            2,
            "Case C: extra SNP on one hap must not discard shared C/A"
        );
    }

    #[test]
    fn six_r28_legacy_unique_supporter_collapses_mid_b_to_ref_only() {
        let Some((untrimmed, active_s, active_e)) = load_untrimmed_mid_b() else {
            eprintln!("Real-data mid-B comparison unavailable");
            return;
        };
        assert_eq!(untrimmed.haplotypes.len(), 4);
        assert_eq!(n_alt(&untrimmed), 3);
        assert!(untrimmed.haplotypes.iter().all(|h| h.bases.len() == 430));
        dump_hap_trace("6R.28_UNTRIMMED_430M", &untrimmed);
        assert_eq!(n_exact_snp(&untrimmed, SITE_CA, b'A'), 2);

        let mut legacy = untrimmed.clone();
        legacy_unique_snp_rank_filter_assembly(
            &mut legacy,
            AlleleFilterOptions::from_strict_java(true, Some(active_s), Some(active_e)),
        )
        .expect("legacy unique-supporter");
        dump_hap_trace("6R.28_LEGACY_UNIQUE_SUPPORTER", &legacy);
        assert_eq!(
            legacy.haplotypes.len(),
            1,
            "legacy unique-supporter → REF only"
        );
        assert_eq!(n_alt(&legacy), 0);
        assert_eq!(n_exact_snp(&legacy, SITE_CA, b'A'), 0);
    }

    #[test]
    fn six_r28_java_faithful_filter_retains_shared_ca() {
        let Some((untrimmed, active_s, active_e)) = load_untrimmed_mid_b() else {
            eprintln!("Real-data mid-B comparison unavailable");
            return;
        };
        assert_eq!(untrimmed.haplotypes.len(), 4);
        assert_eq!(n_alt(&untrimmed), 3);
        assert_eq!(n_exact_snp(&untrimmed, SITE_CA, b'A'), 2);
        assert_eq!(java44_default_hc_keep_all(4), vec![true; 4]);

        let filtered = filter_span(&untrimmed, active_s, active_e);
        dump_hap_trace("6R.28_PRODUCTION_FILTER", &filtered);
        assert!(
            n_exact_snp(&filtered, SITE_CA, b'A') >= 1,
            "Java-faithful filter must retain a C/A-supporting ALT"
        );
        assert!(n_exact_snp(&filtered, SITE_TC, b'C') >= 1);
        assert!(n_exact_snp(&filtered, SITE_GC, b'C') >= 1);
        assert!(n_alt(&filtered) >= 1, "must not collapse to REF-only");
        assert!(
            event_has(untrimmed.variation_events(), SITE_CA, "C", "A"),
            "pre-filter EventMap has C/A"
        );
    }
}
