//! 6R.40: AF-calculator QUAL / MLEAC / MLEAF (and QD jitter) after 6R.39 trim convergence.

#[cfg(test)]
mod traces {
    use crate::engine::{CallRegionArgs, HaplotypeCallerEngine};
    use crate::hc_genotyping_engine::HcGenotypingConfig;
    use crate::region_vcf_emit::try_emit_call_region_variants_with_config;
    use crate::walker_apply::{call_disposition, AssemblyRegionCallDisposition};
    use gatk_core::io::vcf::{InfoValue, VcfRecord};
    use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
    use std::collections::HashMap;
    use std::path::Path;

    const SITE_CA: u64 = 92_317_399;
    const SITE_TC: u64 = 92_317_407;
    const SITE_GC: u64 = 92_317_412;

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
                _ => {}
            }
        }
        m
    }

    fn close(a: f64, b: f64, eps: f64) -> bool {
        (a - b).abs() <= eps
    }

    #[test]
    fn six_r40_qd_java_jitter_first_draw_from_gatk_seed() {
        // 6R.41: QualByDepth uses Utils.getRandomGenerator().nextGaussian().
        crate::annotator::plugins::qual_by_depth::reset_gatk_qual_by_depth_rng();
        let qd = crate::annotator::plugins::qual_by_depth::qual_by_depth(78.32, 2);
        assert!(78.32 / 2.0 > 35.0);
        assert!((qd - 25.36).abs() < 0.005, "qd={qd}");
    }

    #[test]
    fn six_r40_canonical_mid_b_qual_mleac_after_af_loop() {
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
        let args = CallRegionArgs::strict_java();
        let outcome = HaplotypeCallerEngine::call_region(&region, &dict, &ref_fasta, &args)
            .expect("call")
            .expect("Some");
        let gt_cfg = HcGenotypingConfig::strict_java();
        crate::annotator::plugins::qual_by_depth::reset_gatk_qual_by_depth_rng();
        let vcf = try_emit_call_region_variants_with_config(
            &region,
            &outcome,
            "NA12878",
            gt_cfg.stand_emit_confidence,
            &gt_cfg,
        )
        .expect("vcf");
        assert_eq!(vcf.len(), 3);
        for (pos, r, a, jqd) in [
            (SITE_CA, "C", "A", 25.36),
            (SITE_TC, "T", "C", 28.73),
            (SITE_GC, "G", "C", 30.97),
        ] {
            let rec = vcf
                .iter()
                .find(|x| x.position == pos)
                .unwrap_or_else(|| panic!("missing {pos}"));
            let s = rec.samples.first().expect("sample");
            let info = info_map(rec);
            let qual = rec.quality.expect("QUAL");
            eprintln!(
                "VCF {}:{} {}/{} QUAL={qual:.4} GT={:?} AD={:?} PL={:?} MLEAC={} MLEAF={} QD={} java_QD={jqd}",
                rec.chromosome,
                rec.position,
                rec.reference,
                rec.alternate.join(","),
                s.gt.as_ref().map(|g| &g.alleles),
                s.ad,
                s.pl,
                info.get("MLEAC").unwrap_or(&"-".into()),
                info.get("MLEAF").unwrap_or(&"-".into()),
                info.get("QD").unwrap_or(&"-".into()),
            );
            assert_eq!(rec.reference, r);
            assert_eq!(rec.alternate.as_slice(), &[a.to_string()]);
            assert_eq!(
                s.gt.as_ref().map(|g| g.alleles.as_slice()),
                Some([1, 1].as_slice())
            );
            assert_eq!(s.ad.as_deref(), Some([0u32, 2].as_slice()));
            assert_eq!(s.pl.as_deref(), Some([90u32, 6, 0].as_slice()));
            assert!(close(qual, 78.32, 0.02), "QUAL rust={qual} java=78.32");
            assert_eq!(info.get("MLEAC").map(String::as_str), Some("1"));
            let mleaf: f64 = info.get("MLEAF").unwrap().parse().unwrap();
            assert!(close(mleaf, 0.5, 0.001));
            let qd: f64 = info.get("QD").unwrap().parse().unwrap();
            assert!(close(qd, jqd, 0.015), "QD rust={qd} java={jqd}");
        }
    }
}
