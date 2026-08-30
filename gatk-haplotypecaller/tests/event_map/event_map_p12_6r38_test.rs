//! 6R.38 TEST-ONLY: EventMap / genotyping / VCF after graph convergence (6R.37).
//! Diagnostic only. Does not change production EventMap, genotyping, or VCF emission.

#[cfg(test)]
mod traces {
    use crate::assembly_based_caller::{assemble_reads_with_finalized, AssembleReadsArgs};
    use crate::assembly_region_finalize::assembly_reference_read;
    use crate::assembly_region_iterator::AssemblyRegion;
    use crate::assembly_region_trimmer::{
        AssemblyRegionTrimmer, AssemblyRegionTrimmerConfig, TrimVariant,
    };
    use crate::cigar::CigarOperator;
    use crate::engine::{take_call_region_audit, CallRegionArgs, HaplotypeCallerEngine};
    use crate::event_map::{variation_events_for_haplotype, VariationEvent};
    use crate::feature_context::FeatureContext;
    use crate::genome_loc::GenomePosition;
    use crate::haplotype::Haplotype;
    use crate::hc_allele_mapping::hap_base_at_ref_locus;
    use crate::hc_genotyping_engine::HcGenotypingConfig;
    use crate::reference_context::ReferenceContext;
    use crate::region_vcf_emit::try_emit_call_region_variants_with_config;
    use crate::walker_apply::{call_disposition, AssemblyRegionCallDisposition};
    use gatk_core::io::vcf::{InfoValue, VcfRecord};
    use gatk_core::reference::{
        parse_intervals_cli_string, ReferenceWindowCache, SequenceDictionary,
    };
    use std::collections::{BTreeMap, BTreeSet, HashMap};
    use std::path::Path;

    const SITE_CT: u64 = 92_317_361;
    const SITE_CG: u64 = 92_317_371;
    const SITE_CA: u64 = 92_317_399;
    const SITE_TC: u64 = 92_317_407;
    const SITE_GC: u64 = 92_317_412;
    const JAVA_ACTIVE: (u64, u64) = (92_317_262, 92_317_491);
    const JAVA_EXTENDED: (u64, u64) = (92_317_162, 92_317_591);
    /// OBSERVED JAVA EXECUTION (6R.24–6R.30 Docker): trimmed haplotypes 54M.
    const JAVA_TRIM: (u64, u64) = (92_317_379, 92_317_432);
    const JAVA_VCF: &str = include_str!("data/java_6r38_mid_b.vcf");

    #[derive(Clone, Copy)]
    struct JavaSite {
        pos: u64,
        r: &'static str,
        a: &'static str,
        qual: f64,
        ac: i32,
        af: f64,
        an: i32,
        dp: i32,
        excess_het: f64,
        fs: f64,
        mleac: i32,
        mleaf: f64,
        mq: f64,
        qd: f64,
        sor: f64,
        gt: [i32; 2],
        ad: [u32; 2],
        sample_dp: u32,
        gq: i32,
        pl: [u32; 3],
    }

    const JAVA_SITES: [JavaSite; 3] = [
        JavaSite {
            pos: SITE_CA,
            r: "C",
            a: "A",
            qual: 78.32,
            ac: 2,
            af: 1.00,
            an: 2,
            dp: 2,
            excess_het: 0.0,
            fs: 0.0,
            mleac: 1,
            mleaf: 0.5,
            mq: 27.00,
            qd: 25.36,
            sor: 0.693,
            gt: [1, 1],
            ad: [0, 2],
            sample_dp: 2,
            gq: 6,
            pl: [90, 6, 0],
        },
        JavaSite {
            pos: SITE_TC,
            r: "T",
            a: "C",
            qual: 78.32,
            ac: 2,
            af: 1.00,
            an: 2,
            dp: 2,
            excess_het: 0.0,
            fs: 0.0,
            mleac: 1,
            mleaf: 0.5,
            mq: 27.00,
            qd: 28.73,
            sor: 0.693,
            gt: [1, 1],
            ad: [0, 2],
            sample_dp: 2,
            gq: 6,
            pl: [90, 6, 0],
        },
        JavaSite {
            pos: SITE_GC,
            r: "G",
            a: "C",
            qual: 78.32,
            ac: 2,
            af: 1.00,
            an: 2,
            dp: 2,
            excess_het: 0.0,
            fs: 0.0,
            mleac: 1,
            mleaf: 0.5,
            mq: 27.00,
            qd: 30.97,
            sor: 0.693,
            gt: [1, 1],
            ad: [0, 2],
            sample_dp: 2,
            gq: 6,
            pl: [90, 6, 0],
        },
    ];

    fn fixture_paths() -> Option<(std::path::PathBuf, std::path::PathBuf)> {
        let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("..");
        let ref_path = root.join("parity/realworld/assets/hs37d5.simple.fa");
        let bam = root.join("parity/realworld/na12878_20k_b37/NA12878_20k.b37.bam");
        if !ref_path.is_file() || !bam.is_file() {
            return None;
        }
        Some((ref_path, bam))
    }

    fn fnv64(data: &[u8]) -> u64 {
        let mut h = 0xcbf2_9ce4_8422_2325u64;
        for &b in data {
            h ^= u64::from(b);
            h = h.wrapping_mul(0x0100_0000_01b3);
        }
        h
    }

    fn cigar_str(h: &Haplotype) -> String {
        let Some(c) = h.cigar.as_ref() else {
            return "none".into();
        };
        if c.elements.is_empty() {
            return "empty".into();
        }
        c.elements
            .iter()
            .map(|e| {
                let op = match e.operator {
                    CigarOperator::Match => 'M',
                    CigarOperator::Insertion => 'I',
                    CigarOperator::Deletion => 'D',
                    CigarOperator::SoftClip => 'S',
                    CigarOperator::HardClip => 'H',
                };
                format!("{}{op}", e.length)
            })
            .collect()
    }

    fn event_kind(e: &VariationEvent) -> &'static str {
        let rl = e.ref_allele.len();
        let al = e.alt_allele.len();
        if rl == 1 && al == 1 {
            "SNP"
        } else if rl == al {
            "MNP"
        } else if al > rl {
            "INS"
        } else {
            "DEL"
        }
    }

    fn event_key(e: &VariationEvent) -> (u64, String, String) {
        (
            e.start_1based.get(),
            e.ref_allele.clone(),
            e.alt_allele.clone(),
        )
    }

    fn dump_haplotypes(label: &str, haps: &[Haplotype], pad: u64, contig: &str, ref_bytes: &[u8]) {
        let ref_hap = haps.iter().find(|h| h.is_reference);
        eprintln!(
            "{label} n_haps={} pad={pad} ref_len={}",
            haps.len(),
            ref_bytes.len()
        );
        for (i, h) in haps.iter().enumerate() {
            let loc = h.genome_loc.map(|g| (g.start_1based(), g.end_1based()));
            eprintln!(
                "  HAP[{i}] ref={} len={} digest={:016x} cigar={} align0={} loc={:?} 361={} 371={} 399={} 407={} 412={}",
                h.is_reference,
                h.bases.len(),
                fnv64(&h.bases),
                cigar_str(h),
                h.alignment_start_hap_wrt_ref,
                loc,
                hap_base_at_ref_locus(h, pad, SITE_CT)
                    .map(|b| b as char)
                    .unwrap_or('?'),
                hap_base_at_ref_locus(h, pad, SITE_CG)
                    .map(|b| b as char)
                    .unwrap_or('?'),
                hap_base_at_ref_locus(h, pad, SITE_CA)
                    .map(|b| b as char)
                    .unwrap_or('?'),
                hap_base_at_ref_locus(h, pad, SITE_TC)
                    .map(|b| b as char)
                    .unwrap_or('?'),
                hap_base_at_ref_locus(h, pad, SITE_GC)
                    .map(|b| b as char)
                    .unwrap_or('?'),
            );
            if let Some(rh) = ref_hap {
                let ev = variation_events_for_haplotype(h, rh, ref_bytes, pad, 0, contig);
                for e in &ev {
                    eprintln!(
                        "    EVENTMAP hap={i} {} {} {}/{} kind={} len_ref={} len_alt={}",
                        e.start_1based.get(),
                        e.end_1based.get(),
                        e.ref_allele,
                        e.alt_allele,
                        event_kind(e),
                        e.ref_allele.len(),
                        e.alt_allele.len()
                    );
                }
            }
        }
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

    fn info_int_diff(
        info: &HashMap<String, String>,
        pos: u64,
        k: &str,
        exp: i32,
    ) -> Option<String> {
        let got = info.get(k).and_then(|v| v.parse::<i32>().ok());
        if got == Some(exp) {
            None
        } else {
            Some(format!("{pos} {k} rust={got:?} java={exp}"))
        }
    }

    fn info_float_diff(
        info: &HashMap<String, String>,
        pos: u64,
        k: &str,
        exp: f64,
        eps: f64,
    ) -> Option<String> {
        let got = info.get(k).and_then(|v| v.parse::<f64>().ok());
        if got.is_some_and(|g| close(g, exp, eps)) {
            None
        } else {
            Some(format!("{pos} {k} rust={got:?} java={exp}"))
        }
    }

    fn parse_java_vcf_line_count() -> usize {
        JAVA_VCF
            .lines()
            .filter(|l| !l.starts_with('#') && !l.is_empty())
            .count()
    }

    #[test]
    fn six_r38_production_sources_have_no_new_locus_pins() {
        let event_map = include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/src/event_map.rs"));
        let geno = include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/src/hc_genotyping_engine/mod.rs"
        ));
        let emit = include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/src/region_vcf_emit.rs"
        ));
        for src in [event_map, geno, emit] {
            assert!(!src.contains("92317361"));
            assert!(!src.contains("92317371"));
            assert!(!src.contains("TAGAGTTGAAG"));
        }
    }

    #[test]
    fn six_r38_java_vcf_fixture_is_well_formed() {
        assert_eq!(parse_java_vcf_line_count(), 3);
        assert!(JAVA_VCF.contains("92317399\t.\tC\tA\t78.32"));
        assert!(JAVA_VCF.contains("1/1:0,2:2:6:90,6,0"));
    }

    /// 6R.38 found accumulation (`412+20+20+20=472`). 6R.39 matches Java
    /// `Math.max(maxEnd, vc.getEnd() + padding)` → 432.
    #[test]
    fn six_r38_trim_modern_accumulates_snp_padding_per_event() {
        let mut dict = SequenceDictionary::new();
        dict.add_contig("2".into(), 243_199_373);
        let trimmer =
            AssemblyRegionTrimmer::new(AssemblyRegionTrimmerConfig::gatk_defaults(), &dict, "2");
        let region = AssemblyRegion {
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
        };
        let vars = vec![
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
        ];
        let res = trimmer.trim(&region, &vars, None);
        assert!(res.variation_present);
        assert_eq!(res.variant_start, Some(SITE_CA));
        assert_eq!(res.variant_end, Some(SITE_GC));
        assert_eq!(res.padded_variant_start, Some(JAVA_TRIM.0));
        let java_contract_end = SITE_GC + 20;
        assert_eq!(java_contract_end, JAVA_TRIM.1);
        assert_eq!(res.padded_variant_end, Some(java_contract_end));
        assert_ne!(
            res.padded_variant_end,
            Some(SITE_GC + 20 + 20 + 20),
            "6R.39: pad is max(end+pad), not sum"
        );
    }

    #[test]
    fn six_r38_canonical_eventmap_genotyping_vcf() {
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
        assert_eq!(region.reads.len(), 2);

        let mut owned = region.clone();
        let mut assemble_args = AssembleReadsArgs::default();
        assemble_args.strict_java_assembly = true;
        let mut ref_cache = ReferenceWindowCache::new(ref_fasta.clone(), 4);
        let _padded = assembly_reference_read(&dict, &mut ref_cache, &region).expect("pad");
        let untrimmed =
            assemble_reads_with_finalized(&mut owned, &dict, &mut ref_cache, &assemble_args)
                .expect("production assemble")
                .assembly;

        let pad = untrimmed.padded_reference_start_1based();
        dump_haplotypes(
            "UNTRIMMED",
            &untrimmed.haplotypes,
            pad,
            "2",
            untrimmed.reference_bases(),
        );
        assert_eq!(untrimmed.haplotypes.len(), 2, "6R.37: k-best stays 2");
        assert!(untrimmed.haplotypes.iter().any(|h| h.is_reference));
        assert_eq!(
            untrimmed
                .haplotypes
                .iter()
                .filter(|h| !h.is_reference)
                .count(),
            1
        );
        assert!(untrimmed.haplotypes.iter().all(|h| h.bases.len() == 430));
        let n_361 = untrimmed
            .haplotypes
            .iter()
            .filter(|h| hap_base_at_ref_locus(h, pad, SITE_CT) == Some(b'T'))
            .count();
        let n_371 = untrimmed
            .haplotypes
            .iter()
            .filter(|h| hap_base_at_ref_locus(h, pad, SITE_CG) == Some(b'G'))
            .count();
        assert_eq!(n_361, 0);
        assert_eq!(n_371, 0);

        let mut merged: BTreeMap<(u64, String, String), BTreeSet<usize>> = BTreeMap::new();
        let ref_hap = untrimmed
            .haplotypes
            .iter()
            .find(|h| h.is_reference)
            .expect("REF");
        for (i, h) in untrimmed.haplotypes.iter().enumerate() {
            for e in variation_events_for_haplotype(
                h,
                ref_hap,
                untrimmed.reference_bases(),
                pad,
                untrimmed.max_mnp_distance(),
                "2",
            ) {
                merged.entry(event_key(&e)).or_default().insert(i);
            }
        }
        eprintln!("MERGED_EVENTMAP_UNTRIMMED n={}", merged.len());
        for ((pos, r, a), haps) in &merged {
            let java = JAVA_SITES
                .iter()
                .any(|j| j.pos == *pos && j.r == r && j.a == a);
            eprintln!("  pos={pos} {r}/{a} rust_haps={haps:?} java_oracle_snp={java}");
        }
        let rust_keys: BTreeSet<_> = merged.keys().cloned().collect();
        let java_keys: BTreeSet<_> = JAVA_SITES
            .iter()
            .map(|j| (j.pos, j.r.to_string(), j.a.to_string()))
            .collect();
        let rust_only: Vec<_> = rust_keys.difference(&java_keys).cloned().collect();
        let java_only: Vec<_> = java_keys.difference(&rust_keys).cloned().collect();
        eprintln!("EVENTMAP_SET rust_only={rust_only:?} java_only={java_only:?}");
        assert!(
            java_only.is_empty(),
            "Rust untrimmed EventMap missing Java oracle SNPs: {java_only:?}"
        );
        assert!(
            rust_only.is_empty(),
            "Rust untrimmed EventMap has extra alleles vs Java: {rust_only:?}"
        );

        let args = CallRegionArgs::strict_java();
        assert!(args.enable_allele_filtering);
        assert!(!args.enable_read_event_supplement);
        let outcome = HaplotypeCallerEngine::call_region(&region, &dict, &ref_fasta, &args)
            .expect("call")
            .expect("call_region Some");
        let audit = take_call_region_audit();
        let rust_trim = (
            audit.trim_padded_start.unwrap_or(0),
            audit.trim_padded_end.unwrap_or(0),
        );
        eprintln!(
            "TRIM rust={:?} java={:?} match={} n_trim_vars={} overlapping={} n_haps_after_trim={:?} cigars={:?}",
            rust_trim,
            JAVA_TRIM,
            rust_trim == JAVA_TRIM,
            audit.trim_variants.len(),
            audit.n_trim_overlapping,
            audit.n_haps_after_trim,
            audit.hap_cigars
        );
        for v in &audit.trim_variants {
            eprintln!(
                "  TRIM_VAR start={} end={} indel={} overlap={}",
                v.start, v.end, v.is_indel, v.overlaps_active
            );
        }
        // 6R.39: Java Math.max(end+padding) → 92317379–92317432 (54M).
        assert_eq!(
            rust_trim, JAVA_TRIM,
            "trim span must match Java 54M pad after 6R.39"
        );

        let out_pad = outcome.assembly.padded_reference_start_1based();
        dump_haplotypes(
            "AFTER_CALL_REGION",
            &outcome.assembly.haplotypes,
            out_pad,
            "2",
            outcome.assembly.reference_bases(),
        );
        assert_eq!(outcome.assembly.haplotypes.len(), 2);
        assert!(
            outcome
                .assembly
                .haplotypes
                .iter()
                .all(|h| h.bases.len() == 54),
            "6R.39: Java-parity trim is 54M"
        );
        // 361/371 sit left of the trim window. Do not use hap_base_at_ref_locus here:
        // REF haplotypes saturating_sub a pre-window locus to offset 0.
        assert!(outcome
            .assembly
            .haplotypes
            .iter()
            .all(|h| h.genome_loc.is_some_and(|g| g.start_1based() > SITE_CG)));

        let post_keys: BTreeSet<_> = outcome
            .assembly
            .variation_events()
            .iter()
            .map(event_key)
            .collect();
        eprintln!("POST_TRIM_EVENTMAP n={}", post_keys.len());
        for (pos, r, a) in &post_keys {
            eprintln!("  {pos} {r}/{a}");
        }
        assert_eq!(post_keys, java_keys);

        let n_reads = outcome.genotyping_reads.len();
        let n_haps = outcome.assembly.haplotypes.len();
        eprintln!(
            "PAIRHMM n_reads={n_reads} n_haps={n_haps} n_ll_rows={}",
            outcome.read_likelihoods.len()
        );
        let mut matrix = vec![vec![f64::NAN; n_haps]; n_reads.max(1)];
        for row in &outcome.read_likelihoods {
            let r = row.read_index.get();
            let h = row.haplotype_index.get();
            if r < matrix.len() && h < n_haps {
                matrix[r][h] = row.log10_likelihood;
            }
            eprintln!(
                "  LL read={} hap={} log10={:.6}",
                r, h, row.log10_likelihood
            );
        }
        assert_eq!(n_reads, 2);
        assert_eq!(n_haps, 2);
        assert_eq!(outcome.read_likelihoods.len(), 4);
        for row in &matrix {
            assert!(row.iter().all(|x| x.is_finite()));
        }

        eprintln!("GENOTYPED_CALLS n={}", outcome.genotyped_calls.len());
        for c in &outcome.genotyped_calls {
            let e = &c.event;
            eprintln!(
                "  GT {} {}/{} PL={:?} GQ={} AD={:?} DP={} gl={:?}",
                e.start_1based.get(),
                e.ref_allele,
                e.alt_allele,
                c.genotype.format.pl_as_i32(),
                c.genotype.format.gq.as_i32(),
                c.genotype.format.ad_as_i32(),
                c.genotype.format.dp.as_i32(),
                c.genotype.genotype_log10_likelihoods
            );
        }
        assert_eq!(outcome.genotyped_calls.len(), 3);

        let gt_cfg = HcGenotypingConfig::strict_java();
        let vcf = try_emit_call_region_variants_with_config(
            &region,
            &outcome,
            "NA12878",
            gt_cfg.stand_emit_confidence,
            &gt_cfg,
        )
        .expect("vcf");
        eprintln!("VCF n={}", vcf.len());
        for rec in &vcf {
            let s = rec.samples.first();
            let gt = s
                .and_then(|x| x.gt.as_ref())
                .map(|g| {
                    format!(
                        "{}/{}",
                        g.alleles.first().unwrap_or(&-1),
                        g.alleles.get(1).unwrap_or(&-1)
                    )
                })
                .unwrap_or_else(|| ".".into());
            eprintln!(
                "  VCF {}:{} {}/{} QUAL={:?} FILTER={:?} GT={} AD={:?} DP={:?} GQ={:?} PL={:?} INFO={:?}",
                rec.chromosome,
                rec.position,
                rec.reference,
                rec.alternate.join(","),
                rec.quality,
                rec.filter,
                gt,
                s.and_then(|x| x.ad.clone()),
                s.and_then(|x| x.dp),
                s.and_then(|x| x.gq),
                s.and_then(|x| x.pl.clone()),
                info_map(rec)
            );
        }
        assert_eq!(vcf.len(), 3);
        assert!(!vcf
            .iter()
            .any(|r| r.position == SITE_CT || r.position == SITE_CG));

        let mut qual_diffs: Vec<(u64, f64, f64)> = Vec::new();
        let mut info_diffs: Vec<String> = Vec::new();
        for js in &JAVA_SITES {
            let rec = vcf
                .iter()
                .find(|r| r.position == js.pos)
                .unwrap_or_else(|| panic!("missing VCF {}", js.pos));
            assert_eq!(rec.chromosome, "2");
            assert_eq!(rec.reference, js.r);
            assert_eq!(rec.alternate.as_slice(), &[js.a.to_string()]);
            let qual = rec.quality.expect("QUAL");
            if !close(qual, js.qual, 0.015) {
                eprintln!("QUAL_DIFF {} rust={qual} java={}", js.pos, js.qual);
                qual_diffs.push((js.pos, qual, js.qual));
            }
            // FORMAT / alleles must match. QUAL/INFO may still differ (6R.38 remainder;
            // 6R.39 does not change QUAL/MLE/QD). Do not pin pre-fix 78.583.
            let s = rec.samples.first().expect("sample");
            let gt = s.gt.as_ref().expect("GT");
            assert_eq!(gt.alleles.as_slice(), &js.gt);
            assert!(!gt.phased);
            assert_eq!(s.ad.as_deref(), Some(js.ad.as_slice()));
            assert_eq!(s.dp, Some(js.sample_dp));
            assert_eq!(s.gq.map(|g| g as i32), Some(js.gq));
            assert_eq!(s.pl.as_deref(), Some(js.pl.as_slice()));
            let info = info_map(rec);
            info_diffs.extend(info_int_diff(&info, js.pos, "AC", js.ac));
            info_diffs.extend(info_float_diff(&info, js.pos, "AF", js.af, 0.005));
            info_diffs.extend(info_int_diff(&info, js.pos, "AN", js.an));
            info_diffs.extend(info_int_diff(&info, js.pos, "DP", js.dp));
            info_diffs.extend(info_float_diff(
                &info,
                js.pos,
                "ExcessHet",
                js.excess_het,
                0.00015,
            ));
            info_diffs.extend(info_float_diff(&info, js.pos, "FS", js.fs, 0.001));
            info_diffs.extend(info_int_diff(&info, js.pos, "MLEAC", js.mleac));
            info_diffs.extend(info_float_diff(&info, js.pos, "MLEAF", js.mleaf, 0.001));
            info_diffs.extend(info_float_diff(&info, js.pos, "MQ", js.mq, 0.01));
            info_diffs.extend(info_float_diff(&info, js.pos, "QD", js.qd, 0.015));
            info_diffs.extend(info_float_diff(&info, js.pos, "SOR", js.sor, 0.001));
        }
        eprintln!(
            "QUAL_DIFFS n={} {qual_diffs:?} INFO_DIFFS n={} {info_diffs:?}",
            qual_diffs.len(),
            info_diffs.len()
        );
        let _ = (qual_diffs, info_diffs);
    }
}
