//! 6R.41: GATK 4.4 `QualByDepth.fixTooHighQD` jitter (`java.util.Random` seed 47382911).

#[cfg(test)]
mod traces {
    use crate::annotator::plugins::qual_by_depth::{
        fix_too_high_qd_with_rng, hold_process_qd_rng_for_test, qual_by_depth,
        qual_by_depth_with_rng, reset_gatk_qual_by_depth_rng, MAX_QD_BEFORE_FIXING,
    };
    use crate::engine::{CallRegionArgs, HaplotypeCallerEngine};
    use crate::hc_genotyping_engine::HcGenotypingConfig;
    use crate::read_downsample::GatkJavaRng;
    use crate::region_vcf_emit::try_emit_call_region_variants_with_config;
    use crate::walker_apply::{call_disposition, AssemblyRegionCallDisposition};
    use gatk_core::io::vcf::{InfoValue, VcfRecord};
    use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
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

    fn info_qd(rec: &VcfRecord) -> f64 {
        for v in &rec.info {
            if let InfoValue::Float(k, xs) = v {
                if k == "QD" {
                    return *xs.first().expect("QD");
                }
            }
        }
        panic!("missing QD");
    }

    fn close(a: f64, b: f64, eps: f64) -> bool {
        (a - b).abs() <= eps
    }

    #[test]
    fn six_r41_qd_below_threshold_is_qual_over_depth() {
        let mut rng = GatkJavaRng::reset_gatk_default();
        let qd = qual_by_depth_with_rng(40.0, 2, &mut rng);
        assert!((qd - 20.0).abs() < 1e-12);
        assert!(qd < MAX_QD_BEFORE_FIXING);
        let qd2 = qual_by_depth_with_rng(40.0, 2, &mut rng);
        assert!(
            (qd2 - 20.0).abs() < 1e-12,
            "below-threshold must not consume RNG"
        );
    }

    #[test]
    fn six_r41_qd_equal_to_threshold_applies_jitter() {
        let mut rng = GatkJavaRng::reset_gatk_default();
        // Java: `if (QD < 35) return QD` — 35.0 is not strictly less, so jitter.
        let qd = fix_too_high_qd_with_rng(MAX_QD_BEFORE_FIXING, &mut rng);
        assert_eq!(format!("{qd:.2}"), "25.36");
    }

    #[test]
    fn six_r41_qd_above_threshold_uses_shared_gaussian_stream() {
        let mut rng = GatkJavaRng::reset_gatk_default();
        let a = qual_by_depth_with_rng(78.32, 2, &mut rng);
        let b = qual_by_depth_with_rng(78.32, 2, &mut rng);
        let c = qual_by_depth_with_rng(78.32, 2, &mut rng);
        assert_eq!(format!("{a:.2}"), "25.36");
        assert_eq!(format!("{b:.2}"), "28.73");
        assert_eq!(format!("{c:.2}"), "30.97");
    }

    #[test]
    fn six_r41_qd_thread_local_reset_is_deterministic() {
        let _qd = hold_process_qd_rng_for_test();
        reset_gatk_qual_by_depth_rng();
        let first = qual_by_depth(78.32, 2);
        reset_gatk_qual_by_depth_rng();
        let again = qual_by_depth(78.32, 2);
        assert!((first - again).abs() < 1e-15);
        assert_eq!(format!("{first:.2}"), "25.36");
    }

    #[test]
    fn six_r41_canonical_mid_b_qd_matches_java_seeded_stream() {
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
        let _qd = hold_process_qd_rng_for_test();
        reset_gatk_qual_by_depth_rng();
        let vcf = try_emit_call_region_variants_with_config(
            &region,
            &outcome,
            "NA12878",
            gt_cfg.stand_emit_confidence,
            &gt_cfg,
        )
        .expect("vcf");
        assert_eq!(vcf.len(), 3);
        let expect = [(SITE_CA, 25.36), (SITE_TC, 28.73), (SITE_GC, 30.97)];
        for (pos, jqd) in expect {
            let rec = vcf
                .iter()
                .find(|x| x.position == pos)
                .unwrap_or_else(|| panic!("missing {pos}"));
            let s = rec.samples.first().expect("sample");
            let qual = rec.quality.expect("QUAL");
            let qd = info_qd(rec);
            eprintln!(
                "6R.41 {}:{} QUAL={qual:.4} GT={:?} AD={:?} PL={:?} QD={qd:.4} java_QD={jqd}",
                rec.chromosome,
                rec.position,
                s.gt.as_ref().map(|g| &g.alleles),
                s.ad,
                s.pl,
            );
            assert_eq!(s.pl.as_deref(), Some([90u32, 6, 0].as_slice()));
            assert!(close(qual, 78.32, 0.02));
            assert!(
                close(qd, jqd, 0.005),
                "QD rust={qd} java={jqd} (%.2f of 30+3*nextGaussian)"
            );
        }
    }
}
