//! 6R.48: QualByDepth process-global RNG stream (Java `Utils.getRandomGenerator`).
//!
//! Canonical `p12_snp_cluster` evidence: two ActiveFulls on one interval must continue
//! one Gaussian stream (draws 5–9 on the 716-cluster), not reseed per region.

#[cfg(test)]
mod traces {
    use crate::annotator::plugins::qual_by_depth::{
        hold_process_qd_rng_for_test, qd_process_gaussian_draw_count, reset_gatk_qual_by_depth_rng,
        take_qd_draw_log,
    };
    use crate::engine::{CallRegionArgs, HaplotypeCallerEngine};
    use crate::hc_genotyping_engine::HcGenotypingConfig;
    use crate::region_vcf_emit::try_emit_call_region_variants_with_config;
    use crate::walker_apply::{call_disposition, AssemblyRegionCallDisposition};
    use gatk_core::io::vcf::{InfoValue, VcfRecord};
    use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
    use std::path::Path;

    const INTERVAL: &str = "2:92305500-92305850";
    const JAVA_QD: [(u64, &str); 9] = [
        (92_305_634, "25.36"),
        (92_305_635, "28.73"),
        (92_305_653, "30.97"),
        (92_305_670, "27.24"),
        (92_305_716, "28.20"),
        (92_305_719, "25.00"),
        (92_305_722, "29.56"),
        (92_305_726, "30.62"),
        (92_305_728, "28.17"),
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

    fn info_qd_printed(rec: &VcfRecord) -> String {
        for v in &rec.info {
            if let InfoValue::Float(k, xs) = v {
                if k == "QD" {
                    return format!("{:.2}", xs.first().copied().unwrap_or(0.0));
                }
            }
        }
        panic!("missing QD");
    }

    #[test]
    fn six_r48_snp_cluster_two_activefulls_share_one_qd_stream() {
        let Some((ref_fasta, bam_path)) = fixture_paths() else {
            eprintln!("Real-data snp_cluster comparison unavailable");
            return;
        };
        let dict = SequenceDictionary::from_fasta_path(&ref_fasta).expect("dict");
        let walk_iv = parse_intervals_cli_string(&dict, INTERVAL).expect("iv");
        let filters = crate::read_model::ReadFilterParams::gatk_standard_hc();
        let cfg =
            crate::walker_traversal::WalkerTraversalConfig::gatk_haplotype_caller_production(100);
        let walk = crate::walker_traversal::traverse_assembly_region_walker(
            &dict, &walk_iv, &ref_fasta, &bam_path, &filters, &cfg,
        )
        .expect("walk");
        let mut regions: Vec<_> = crate::walker_traversal::flatten_assembly_regions(&walk)
            .into_iter()
            .filter(|r| {
                matches!(
                    call_disposition(r),
                    AssemblyRegionCallDisposition::ActiveFull
                )
            })
            .collect();
        regions.sort_by(|a, b| a.start.get().cmp(&b.start.get()));
        assert!(
            regions.len() >= 2,
            "expected two ActiveFull clusters, got {}",
            regions.len()
        );

        let args = CallRegionArgs::strict_java();
        let gt_cfg = HcGenotypingConfig::strict_java();
        let _qd = hold_process_qd_rng_for_test();
        reset_gatk_qual_by_depth_rng();
        let _ = take_qd_draw_log();
        let mut all = Vec::new();
        for region in &regions {
            let outcome = HaplotypeCallerEngine::call_region(region, &dict, &ref_fasta, &args)
                .expect("call")
                .expect("Some");
            let vcf = try_emit_call_region_variants_with_config(
                region,
                &outcome,
                "NA12878",
                gt_cfg.stand_emit_confidence,
                &gt_cfg,
            )
            .expect("vcf");
            all.extend(vcf);
        }
        all.sort_by(|a, b| a.position.cmp(&b.position));
        eprintln!(
            "6R.48 snp_cluster sequential emit: {} records, {} gaussian draws",
            all.len(),
            qd_process_gaussian_draw_count()
        );
        for rec in &all {
            eprintln!(
                "  {}:{} QD={}",
                rec.chromosome,
                rec.position,
                info_qd_printed(rec)
            );
        }
        for (pos, jqd) in JAVA_QD {
            let rec = all
                .iter()
                .find(|x| x.position == pos)
                .unwrap_or_else(|| panic!("missing {pos}"));
            let got = info_qd_printed(rec);
            assert_eq!(got, jqd, "pos {pos} QD rust={got} java={jqd}");
        }
        assert_eq!(
            qd_process_gaussian_draw_count(),
            9,
            "Java 716-cluster is draws 5–9 of one process stream; no pre-cluster draws on this interval"
        );
        let log = take_qd_draw_log();
        let first_716 = log
            .iter()
            .find(|t| t.site.as_ref().is_some_and(|(_, p)| *p == 92_305_716))
            .expect("716 in draw log");
        assert_eq!(first_716.gaussian_ordinal, 5);
        assert_eq!(format!("{:.2}", first_716.result_qd), "28.20");
        assert_eq!(first_716.stream_gaussians_before, 4);
    }
}
