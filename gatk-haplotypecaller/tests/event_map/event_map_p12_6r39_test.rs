//! 6R.39: Java-parity `trim_modern` max-end (Math.max, not per-event accumulation).
//! Does not change EventMap, PairHMM, genotyping, QUAL, or VCF emission logic.

#[cfg(test)]
mod traces {
    use crate::assembly_based_caller::{assemble_reads_with_finalized, AssembleReadsArgs};
    use crate::assembly_region_finalize::assembly_reference_read;
    use crate::assembly_region_iterator::AssemblyRegion;
    use crate::assembly_region_trimmer::{
        AssemblyRegionTrimmer, AssemblyRegionTrimmerConfig, TrimVariant,
    };
    use crate::engine::{take_call_region_audit, CallRegionArgs, HaplotypeCallerEngine};
    use crate::event_map::variation_events_for_haplotype;
    use crate::feature_context::FeatureContext;
    use crate::genome_loc::GenomePosition;
    use crate::hc_allele_mapping::hap_base_at_ref_locus;
    use crate::hc_genotyping_engine::HcGenotypingConfig;
    use crate::reference_context::ReferenceContext;
    use crate::region_vcf_emit::try_emit_call_region_variants_with_config;
    use crate::walker_apply::{call_disposition, AssemblyRegionCallDisposition};
    use gatk_core::io::vcf::{InfoValue, VcfRecord};
    use gatk_core::reference::{
        parse_intervals_cli_string, ReferenceWindowCache, SequenceDictionary,
    };
    use std::collections::{BTreeSet, HashMap};
    use std::path::Path;

    const SITE_CT: u64 = 92_317_361;
    const SITE_CG: u64 = 92_317_371;
    const SITE_CA: u64 = 92_317_399;
    const SITE_TC: u64 = 92_317_407;
    const SITE_GC: u64 = 92_317_412;
    const JAVA_ACTIVE: (u64, u64) = (92_317_262, 92_317_491);
    const JAVA_EXTENDED: (u64, u64) = (92_317_162, 92_317_591);
    const JAVA_TRIM: (u64, u64) = (92_317_379, 92_317_432);
    const SNP_PAD: u64 = 20;
    const PRE_FIX_ACCUMULATED_END: u64 = 92_317_472;

    fn fixture_paths() -> Option<(std::path::PathBuf, std::path::PathBuf)> {
        let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("..");
        let ref_path = root.join("parity/realworld/assets/hs37d5.simple.fa");
        let bam = root.join("parity/realworld/na12878_20k_b37/NA12878_20k.b37.bam");
        if !ref_path.is_file() || !bam.is_file() {
            return None;
        }
        Some((ref_path, bam))
    }

    fn info_map(rec: &VcfRecord) -> HashMap<String, String> {
        let mut m = HashMap::new();
        for v in &rec.info {
            match v {
                InfoValue::Flag(k) => {
                    m.insert(k.clone(), "true".into());
                }
                InfoValue::Integer(k, xs) => {
                    m.insert(
                        k.clone(),
                        xs.iter()
                            .map(|x| x.to_string())
                            .collect::<Vec<_>>()
                            .join(","),
                    );
                }
                InfoValue::Float(k, xs) => {
                    m.insert(
                        k.clone(),
                        xs.iter()
                            .map(|x| format!("{x:.4}"))
                            .collect::<Vec<_>>()
                            .join(","),
                    );
                }
                InfoValue::String(k, xs) => {
                    m.insert(k.clone(), xs.join(","));
                }
                InfoValue::Character(k, xs) => {
                    m.insert(k.clone(), xs.iter().collect());
                }
            }
        }
        m
    }

    fn close(a: f64, b: f64, eps: f64) -> bool {
        (a - b).abs() <= eps
    }

    fn mid_b_region() -> AssemblyRegion {
        AssemblyRegion {
            contig: "2".into(),
            start: GenomePosition::new_1based(JAVA_ACTIVE.0),
            end: GenomePosition::new_1based(JAVA_ACTIVE.1),
            is_active: true,
            extended_start: GenomePosition::new_1based(JAVA_EXTENDED.0),
            extended_end: GenomePosition::new_1based(JAVA_EXTENDED.1),
            extension: 100,
            reads: Vec::new(),
            read_qnames: Vec::new(),
            reference: ReferenceContext::empty(),
            features: FeatureContext::empty(),
            pileup_loci: Vec::new(),
        }
    }

    fn oracle_snps() -> Vec<TrimVariant> {
        vec![
            TrimVariant {
                contig: "2".into(),
                start: SITE_CA,
                end: SITE_CA,
                is_indel: false,
            },
            TrimVariant {
                contig: "2".into(),
                start: SITE_TC,
                end: SITE_TC,
                is_indel: false,
            },
            TrimVariant {
                contig: "2".into(),
                start: SITE_GC,
                end: SITE_GC,
                is_indel: false,
            },
        ]
    }

    /// Three SNP ends + pad: max(399+20, 407+20, 412+20) = 432, not 412+60.
    #[test]
    fn six_r39_trim_modern_three_snps_pad_is_max_not_sum() {
        let mut dict = SequenceDictionary::new();
        dict.add_contig("2".into(), 243_199_373);
        let trimmer =
            AssemblyRegionTrimmer::new(AssemblyRegionTrimmerConfig::gatk_defaults(), &dict, "2");
        let res = trimmer.trim(&mid_b_region(), &oracle_snps(), None);
        let per_event_end_plus_pad = [SITE_CA + SNP_PAD, SITE_TC + SNP_PAD, SITE_GC + SNP_PAD];
        let java_max = *per_event_end_plus_pad.iter().max().unwrap();
        assert_eq!(java_max, JAVA_TRIM.1);
        assert_eq!(res.padded_variant_start, Some(JAVA_TRIM.0));
        assert_eq!(res.padded_variant_end, Some(java_max));
        assert_ne!(res.padded_variant_end, Some(PRE_FIX_ACCUMULATED_END));
        assert_eq!(
            SITE_GC + SNP_PAD + SNP_PAD + SNP_PAD,
            PRE_FIX_ACCUMULATED_END
        );
    }

    #[test]
    fn six_r39_trim_modern_source_uses_max_not_accumulate() {
        let src = include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/src/assembly_region_trimmer.rs"
        ));
        let modern = src
            .split("fn trim_modern")
            .nth(1)
            .and_then(|s| s.split("fn trim_legacy").next())
            .expect("trim_modern body");
        assert!(
            modern.contains("max_end.max(v.end.saturating_add(padding as u64))"),
            "trim_modern must use Math.max(maxEnd, end+padding)"
        );
        assert!(
            !modern.contains("max_end.saturating_add(padding"),
            "trim_modern must not accumulate padding onto max_end"
        );
    }

    #[test]
    fn six_r39_canonical_mid_b_trim_54m_remeasurement() {
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
        assert_eq!((region.start.get(), region.end.get()), JAVA_ACTIVE);
        assert_eq!(
            (region.extended_start.get(), region.extended_end.get()),
            JAVA_EXTENDED
        );

        let mut owned = region.clone();
        let mut assemble_args = AssembleReadsArgs::default();
        assemble_args.strict_java_assembly = true;
        let mut ref_cache = ReferenceWindowCache::new(ref_fasta.clone(), 4);
        let _padded = assembly_reference_read(&dict, &mut ref_cache, &region).expect("pad");
        let untrimmed =
            assemble_reads_with_finalized(&mut owned, &dict, &mut ref_cache, &assemble_args)
                .expect("assemble")
                .assembly;
        assert_eq!(untrimmed.haplotypes.len(), 2);
        assert!(untrimmed.haplotypes.iter().all(|h| h.bases.len() == 430));

        let args = CallRegionArgs::strict_java();
        let outcome = HaplotypeCallerEngine::call_region(&region, &dict, &ref_fasta, &args)
            .expect("call")
            .expect("call_region Some");
        let audit = take_call_region_audit();
        let rust_trim = (
            audit.trim_padded_start.unwrap_or(0),
            audit.trim_padded_end.unwrap_or(0),
        );
        eprintln!("TRIM rust={:?} java={:?}", rust_trim, JAVA_TRIM);
        assert_eq!(rust_trim, JAVA_TRIM);
        assert!(outcome
            .assembly
            .haplotypes
            .iter()
            .all(|h| h.bases.len() == 54));
        assert_eq!(outcome.assembly.haplotypes.len(), 2);

        let pad = untrimmed.padded_reference_start_1based();
        let ref_hap = untrimmed
            .haplotypes
            .iter()
            .find(|h| h.is_reference)
            .expect("REF");
        let mut keys: BTreeSet<(u64, String, String)> = BTreeSet::new();
        for h in &untrimmed.haplotypes {
            for e in variation_events_for_haplotype(
                h,
                ref_hap,
                untrimmed.reference_bases(),
                pad,
                untrimmed.max_mnp_distance(),
                "2",
            ) {
                keys.insert((
                    e.start_1based.get(),
                    e.ref_allele.clone(),
                    e.alt_allele.clone(),
                ));
            }
        }
        let oracle: BTreeSet<_> = [
            (SITE_CA, "C".into(), "A".into()),
            (SITE_TC, "T".into(), "C".into()),
            (SITE_GC, "G".into(), "C".into()),
        ]
        .into_iter()
        .collect();
        assert_eq!(keys, oracle);
        let post: BTreeSet<_> = outcome
            .assembly
            .variation_events()
            .iter()
            .map(|e| {
                (
                    e.start_1based.get(),
                    e.ref_allele.clone(),
                    e.alt_allele.clone(),
                )
            })
            .collect();
        eprintln!("POST_TRIM_EVENTMAP {post:?}");
        assert_eq!(post, oracle);
        assert!(!post.iter().any(|(p, _, _)| *p == SITE_CT || *p == SITE_CG));

        let out_pad = outcome.assembly.padded_reference_start_1based();
        for h in &outcome.assembly.haplotypes {
            assert!(h
                .genome_loc
                .is_some_and(|g| g.start_1based() == JAVA_TRIM.0 && g.end_1based() == JAVA_TRIM.1));
            if h.is_reference {
                continue;
            }
            assert_eq!(hap_base_at_ref_locus(h, out_pad, SITE_CA), Some(b'A'));
            assert_eq!(hap_base_at_ref_locus(h, out_pad, SITE_TC), Some(b'C'));
            assert_eq!(hap_base_at_ref_locus(h, out_pad, SITE_GC), Some(b'C'));
        }

        assert_eq!(outcome.genotyped_calls.len(), 3);
        for c in &outcome.genotyped_calls {
            eprintln!(
                "  GT {} {}/{} PL={:?} GQ={} AD={:?} DP={}",
                c.event.start_1based.get(),
                c.event.ref_allele,
                c.event.alt_allele,
                c.genotype.format.pl_as_i32(),
                c.genotype.format.gq.as_i32(),
                c.genotype.format.ad_as_i32(),
                c.genotype.format.dp.as_i32(),
            );
            assert_eq!(c.genotype.format.pl_as_i32(), vec![90, 6, 0]);
            assert_eq!(c.genotype.format.gq.as_i32(), 6);
            assert_eq!(c.genotype.format.ad_as_i32(), vec![0, 2]);
            assert_eq!(c.genotype.format.dp.as_i32(), 2);
        }

        let gt_cfg = HcGenotypingConfig::strict_java();
        let vcf = try_emit_call_region_variants_with_config(
            &region,
            &outcome,
            "NA12878",
            gt_cfg.stand_emit_confidence,
            &gt_cfg,
        )
        .expect("vcf");
        assert_eq!(vcf.len(), 3);
        let java_sites = [
            (SITE_CA, "C", "A", 78.32, 25.36),
            (SITE_TC, "T", "C", 78.32, 28.73),
            (SITE_GC, "G", "C", 78.32, 30.97),
        ];
        for (pos, r, a, jqual, jqd) in java_sites {
            let rec = vcf
                .iter()
                .find(|x| x.position == pos)
                .unwrap_or_else(|| panic!("missing {pos}"));
            let s = rec.samples.first().expect("sample");
            let info = info_map(rec);
            let qual = rec.quality.expect("QUAL");
            eprintln!(
                "VCF {}:{} {}/{} QUAL={qual} FILTER={:?} GT={:?} AD={:?} DP={:?} GQ={:?} PL={:?} INFO={:?}",
                rec.chromosome,
                rec.position,
                rec.reference,
                rec.alternate.join(","),
                rec.filter,
                s.gt.as_ref().map(|g| &g.alleles),
                s.ad,
                s.dp,
                s.gq,
                s.pl,
                info
            );
            assert_eq!(rec.chromosome, "2");
            assert_eq!(rec.reference, r);
            assert_eq!(rec.alternate.as_slice(), &[a.to_string()]);
            assert_eq!(
                s.gt.as_ref().map(|g| g.alleles.as_slice()),
                Some([1, 1].as_slice())
            );
            assert_eq!(s.ad.as_deref(), Some([0u32, 2].as_slice()));
            assert_eq!(s.dp, Some(2));
            assert_eq!(s.gq.map(|g| g as i32), Some(6));
            assert_eq!(s.pl.as_deref(), Some([90u32, 6, 0].as_slice()));
            eprintln!(
                "  QUAL_VS_JAVA rust={qual} java={jqual} match={}",
                close(qual, jqual, 0.015)
            );
            eprintln!(
                "  MLEAC={} MLEAF={} QD={} java_QD={jqd}",
                info.get("MLEAC").unwrap_or(&"-".into()),
                info.get("MLEAF").unwrap_or(&"-".into()),
                info.get("QD").unwrap_or(&"-".into()),
            );
        }
    }
}
