/// L12-C: AD/PL reshape stage between PairHMM genotyping and [`GenotypeFinalize`].
/// Owns the Class-A / A2 / A3 family that used to sit inline in
/// `try_genotype_variation_event`. Orchestration still decides *when* to call this;
/// reshape policy lives here.
pub struct SiteReshape;

impl SiteReshape {
    /// Class-A / A2 / A3: when informative AD is skewed vs balanced pileup, restore
    /// pileup AD (keep PairHMM PLs if already het) or sparse-pileup genotype (GT flip).
    /// No-op outside the evidence class (P12 production emit scope excluded).
    pub fn apply_class_a_family(
        mut gt: RegionGenotypeResult,
        event: &VariationEvent,
        read_ref_ad: i32,
        read_alt_ad: i32,
        config: &HcGenotypingConfig,
    ) -> GatkResult<RegionGenotypeResult> {
        // 6R.155 diagnostic only: prove Class-A family is the first post-AD mutation.
        // Default off. Not a product contract (`HOLDOUT_6R155` causality).
        if crate::parity_harness::env_flag_true("GATK_RS_DIAGNOSTIC_SKIP_SITE_RESHAPE")
        {
            return Ok(gt);
        }
        // Dense GIAB hets often have balanced pileup AD but skewed informative AD.
        if crate::read_event_discovery::is_strict_java_p12_production_emit_scope(event)
            || gt.format.ad.len() < 2
            || read_ref_ad < 1
            || read_alt_ad < 1
            || read_ref_ad.saturating_mul(2) < read_alt_ad
        {
            return Ok(gt);
        }
        let info_ref = gt.format.ad[0].as_i32();
        let info_alt = gt.format.ad[1].as_i32();
        let class_a = info_ref == 0 && info_alt >= 2;
        let class_a2 = event.is_snp()
            && read_alt_ad.saturating_mul(2) >= read_ref_ad
            && (biallelic_genotype_index_from_pl(&gt.format.pl).get() == 2
                || (info_alt >= 2 && (info_ref == 0 || info_alt >= info_ref.saturating_mul(3))));
        // Class-A3: pileup supports alt/het but informative AD is heavily REF-skewed
        // (e.g. 20:10037037 pileup ~10,14 vs informative ~20,3).
        let class_a3 = event.is_snp()
            && read_alt_ad.saturating_mul(2) >= read_ref_ad
            && info_ref >= 1
            && info_alt >= 1
            && info_ref >= info_alt.saturating_mul(3);
        // 6R.158: Java 4.4 never pileup-replaces an already-assigned calculator genotype.
        // SiteScore has already assigned GT/PL/AD/GQ/DP. Class-A3 must not mutate it.
        // Do not use pl_gt / GT / AD-ratio / coordinate heuristics here.
        if class_a3 {
            return Ok(gt);
        }
        if !(class_a || class_a2 || class_a3) {
            return Ok(gt);
        }
        let pl_gt = biallelic_genotype_index_from_pl(&gt.format.pl).get();
        // L9 PL: when PairHMM already called het, keep its PLs and only restore
        // pileup AD. SparsePlShape (81,0,36) is reserved for GT flips (PL=1/1).
        if pl_gt == 1 {
            gt = reshape_genotype_allele_depths_keep_pls(gt, read_ref_ad, read_alt_ad);
        } else {
            gt = sparse_snp_genotype_from_read_depths(read_ref_ad, read_alt_ad, config)?;
        }
        Ok(gt)
    }
}

/// 6R.158 Case E: the generic sparse helper is unchanged (not the A3 caller).
#[cfg(test)]
#[test]
fn forensic_6r158_case_e_sparse_helper_unchanged() {
    let Ok(shaped) =
        sparse_snp_genotype_from_read_depths(22, 24, &HcGenotypingConfig::strict_java())
    else {
        assert!(false, "sparse helper should succeed");
        return;
    };
    assert_eq!(best_pl_index(&shaped.format.pl), 1);
    assert_eq!(shaped.format.pl_as_i32(), vec![81, 0, 36]);
    assert_eq!(shaped.format.ad_as_i32(), vec![22, 24]);
    assert_eq!(shaped.format.gq.as_i32(), 36);
    assert_eq!(shaped.format.dp.as_i32(), 46);
    assert_eq!(SparsePlShape::from_pileup_depths(22, 24), SparsePlShape::Het);
}

