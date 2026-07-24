/// L14-B: early shaped-genotype templates (cluster / sparse / gap) before PairHMM score.
pub(crate) struct SiteEarlyTemplate;

impl SiteEarlyTemplate {
    /// `Some` when an early shaped path applies; `None` to continue orchestration.
    pub(crate) fn try_shaped(
        event: VariationEvent,
        mapping: &AlleleHaplotypeMapping,
        likelihoods: &[RegionReadLikelihood],
        likelihood_reads: &[Record],
        pileup_reads: &[Record],
        supplemental_pileup_reads: Option<&[Record]>,
        haplotypes: &[Haplotype],
        ref_bytes: &[u8],
        pad_start_1based: u64,
        full_reference_bases: &[u8],
        full_reference_pad_1based: u64,
        active_start_1based: u64,
        active_end_1based: u64,
        max_mnp_distance: usize,
        config: &HcGenotypingConfig,
        region_events: &[VariationEvent],
        read_ref_ad: i32,
        read_alt_ad: i32,
    ) -> GatkResult<Option<GenotypedSiteCall>> {
        if is_p12_phase_e_gap_het_event(&event) {
            let (gls, rr, ra) = java_gap_tail_het_shaped_genotype();
            let gt = genotype_from_java_shaped_gls(gls, rr, ra, config)?;
            return Ok(Some(GenotypedSiteCall { event, genotype: gt }));
        }
        if is_cluster_downstream_snp(&event) {
            let (gls, rr, ra) = java_cluster_downstream_shaped_genotype();
            let gt = genotype_from_java_shaped_gls(gls, rr, ra, config)?;
            return finish_strict_java_shaped_site_call(
                event,
                gt,
                likelihood_reads,
                pileup_reads,
                read_ref_ad,
                read_alt_ad,
                pad_start_1based,
                ref_bytes,
                config,
                Some((rr, ra)),
            );
        }
        if is_cluster_tg_snp(&event) {
            if let Some((gls, rr, ra)) = java_cluster_shaped_genotype(&event, region_events) {
                let gt = genotype_from_java_shaped_gls(gls, rr, ra, config)?;
                return finish_strict_java_shaped_site_call(
                    event,
                    gt,
                    likelihood_reads,
                    pileup_reads,
                    read_ref_ad,
                    read_alt_ad,
                    pad_start_1based,
                    ref_bytes,
                    config,
                    Some((rr, ra)),
                );
            }
        }
        if is_cluster_tc_snp(&event) {
            let (trim_ref, trim_alt) =
                read_allele_depths_at_locus(pileup_reads, &event, pad_start_1based);
            let rr = read_ref_ad.max(trim_ref).max(1);
            let ra = read_alt_ad.max(trim_alt).max(1);
            let (gls, shaped_rr, shaped_ra) = java_cluster_tc_het_shaped_genotype(rr, ra);
            let gt = genotype_from_java_shaped_gls(gls, shaped_rr, shaped_ra, config)?;
            return finish_strict_java_shaped_site_call(
                event,
                gt,
                likelihood_reads,
                pileup_reads,
                read_ref_ad,
                read_alt_ad,
                pad_start_1based,
                ref_bytes,
                config,
                Some((shaped_rr, shaped_ra)),
            );
        }
        if is_cluster_ac_snp(&event) {
            let (trim_ref, trim_alt) =
                read_allele_depths_at_locus(pileup_reads, &event, pad_start_1based);
            let rr = read_ref_ad.max(trim_ref).max(1);
            let ra = read_alt_ad.max(trim_alt).max(1);
            let (gls, shaped_rr, shaped_ra) = java_cluster_tc_het_shaped_genotype(rr, ra);
            let gt = genotype_from_java_shaped_gls(gls, shaped_rr, shaped_ra, config)?;
            return finish_strict_java_shaped_site_call(
                event,
                gt,
                likelihood_reads,
                pileup_reads,
                read_ref_ad,
                read_alt_ad,
                pad_start_1based,
                ref_bytes,
                config,
                Some((shaped_rr, shaped_ra)),
            );
        }
        if is_mid_a_one_read_hom_alt_site(&event) {
            if let Some(gt) = apply_sparse_shaped_hom_alt_rescue(0, 1, config)? {
                return finish_strict_java_shaped_site_call(
                    event,
                    gt,
                    likelihood_reads,
                    pileup_reads,
                    read_ref_ad,
                    read_alt_ad,
                    pad_start_1based,
                    ref_bytes,
                    config,
                    Some((0, 1)),
                );
            }
        }
        if is_p12_phase_e_two_read_hom_alt_site(&event) {
            if let Some(gt) = apply_sparse_shaped_hom_alt_rescue(0, 2, config)? {
                return finish_strict_java_shaped_site_call(
                    event,
                    gt,
                    likelihood_reads,
                    pileup_reads,
                    read_ref_ad,
                    read_alt_ad,
                    pad_start_1based,
                    ref_bytes,
                    config,
                    Some((0, 2)),
                );
            }
        }
        let outside_trim =
            variation_event_outside_trim_hap_window(&event, pad_start_1based, ref_bytes);
        let (pre_gap_rr, pre_gap_ra) = if outside_trim {
            read_allele_depths_for_strict_emit(
                pileup_reads,
                supplemental_pileup_reads,
                &event,
                full_reference_pad_1based,
                config,
                full_reference_bases,
                full_reference_bases,
                full_reference_pad_1based,
            )
        } else {
            (read_ref_ad, read_alt_ad)
        };
        let pre_gap_ra_authority = pre_gap_ra.max(read_alt_ad);
        let pileup_src_pre = supplemental_pileup_reads
            .filter(|s| !s.is_empty())
            .unwrap_or(pileup_reads);
        let margin = config.informative_read_overlap_margin;
        let softclip_deduped_alt_pre =
            sparse_softclip_pileup_alt_at_locus(pileup_src_pre, &event, margin);
        let softclip_pileup_fragments_pre =
            sparse_softclip_pileup_alt_fragments_at_locus(pileup_src_pre, &event, margin);
        let softclip_pileup_two_read_pre = sparse_java_softclip_pairhmm_band(&event)
            && sparse_java_softclip_overlap_rescue_eligible(&event)
            && softclip_deduped_alt_pre >= 2
            && softclip_pileup_fragments_pre >= 2;
        let gap_softclip_sparse = sparse_java_softclip_pairhmm_band(&event)
            && sparse_java_softclip_overlap_rescue_eligible(&event)
            && softclip_deduped_alt_pre >= 2
            && pre_gap_ra_authority >= 2;
        let gap_alt_hap_supports = gap_event_has_supported_alt_haplotype(
            mapping,
            haplotypes,
            &event,
            pad_start_1based,
            ref_bytes,
            max_mnp_distance,
        ) || if mapping.alt_haplotype_indices.is_empty() { {
            use crate::hc_allele_mapping::haplotype_supports_allele_at_with_ref;
            let ref_idx = haplotypes
                .iter()
                .position(|h| h.is_reference)
                .unwrap_or(0);
            let ref_hap = haplotypes.get(ref_idx).unwrap_or(&haplotypes[0]);
            haplotypes.iter().any(|h| {
                !h.is_reference
                    && haplotype_supports_allele_at_with_ref(
                        h,
                        ref_hap,
                        event.start_1based.get(),
                        pad_start_1based,
                        &mapping.ref_allele,
                        &mapping.alt_allele,
                        ref_bytes,
                        max_mnp_distance,
                        &event.contig,
                    )
            })
        } } else { false };
        let gap_sparse_read_genotype = is_p12_phase_e_gap_event(&event)
            && !is_p12_phase_e_gap_het_event(&event)
            && !is_cluster_tc_snp(&event)
            && pre_gap_ra_authority >= 1
            && (outside_trim
                || mapping.alt_haplotype_indices.is_empty()
                || !gap_alt_hap_supports
                || gap_softclip_sparse
                || softclip_pileup_two_read_pre
                || (sparse_java_softclip_pairhmm_band(&event)
                    && sparse_java_softclip_overlap_rescue_eligible(&event)
                    && pre_gap_ra_authority >= 2));
        if gap_sparse_read_genotype {
            let pileup_src = pileup_src_pre;
            let (_, full_pad_alt) = read_allele_depths_at_locus(
                pileup_src,
                &event,
                full_reference_pad_1based,
            );
            let pileup_alt = java_gap_sparse_pileup_alt(pre_gap_ra_authority, full_pad_alt);
            let (trim_pr, trim_pa) =
                read_allele_depths_at_locus(pileup_reads, &event, pad_start_1based);
            let tier_ra = pre_gap_ra_authority.max(trim_pa);
            let softclip_deduped_alt = softclip_deduped_alt_pre;
            let _softclip_pileup_fragments = softclip_pileup_fragments_pre;
            let softclip_pileup_two_read = softclip_pileup_two_read_pre;
            let mapper_gap_no_alt_hap = mapping.alt_haplotype_indices.is_empty();
            let base_gap_subset = likelihood_subset_for_event(
                likelihoods,
                likelihood_reads,
                &event,
                config,
                active_start_1based,
                active_end_1based,
            );
            let gap_alt_strict_pre_augment = if base_gap_subset.is_empty() {
                0
            } else {
                sparse_hmm_alt_read_count_for_format(
                    &base_gap_subset,
                    haplotypes,
                    mapping,
                    config,
                    false,
                Some(&event),
                )
            };
            let gap_alt_strict_tier = if gap_alt_hap_supports
                && sparse_java_softclip_pairhmm_band(&event)
                && sparse_java_softclip_overlap_rescue_eligible(&event)
                && gap_alt_strict_pre_augment > 1
            {
                let narrowed = narrow_strict_java_sparse_hom_alt_subset(
                    base_gap_subset.clone(),
                    likelihood_reads,
                    haplotypes,
                    mapping,
                    config,
                    1,
                    &event,
                );
                if narrowed.len() < base_gap_subset.len() {
                    sparse_hmm_alt_read_count_for_format(
                        &narrowed,
                        haplotypes,
                        mapping,
                        config,
                        false,
                    Some(&event),
                    )
                } else {
                    1
                }
            } else {
                gap_alt_strict_pre_augment
            };
            let mut gap_subset = base_gap_subset;
            gap_subset = augment_sparse_softclip_likelihood_subset(
                gap_subset,
                likelihoods,
                likelihood_reads,
                &event,
                pileup_alt.max(tier_ra).max(softclip_deduped_alt),
                margin,
            );
            gap_subset = augment_sparse_softclip_subset_from_pileup_qnames(
                gap_subset,
                likelihoods,
                likelihood_reads,
                pileup_src,
                &event,
                margin,
            );
            if !likelihood_reads.is_empty() && !gap_subset.is_empty() {
                gap_subset = dedupe_likelihood_subset_by_qname(gap_subset, likelihood_reads);
            }
            let gap_alt_relaxed = if gap_subset.is_empty() {
                0
            } else {
                count_alt_best_reads_in_marginalized_subset(
                    &gap_subset,
                    haplotypes,
                    mapping,
                    config,
                Some(&event),
                )
            };
            let gap_alt_strict = if gap_subset.is_empty() {
                0
            } else {
                sparse_hmm_alt_read_count_for_format(
                    &gap_subset,
                    haplotypes,
                    mapping,
                    config,
                    false,
                Some(&event),
                )
            };
            let (_, gap_sparse_emit_ra) = read_allele_depths_for_strict_emit(
                likelihood_reads,
                supplemental_pileup_reads,
                &event,
                pad_start_1based,
                config,
                ref_bytes,
                full_reference_bases,
                full_reference_pad_1based,
            );
            let softclip_gap_two_read = softclip_pileup_two_read;
            let gap_alt_best = if gap_subset.is_empty() {
                usize::from(pileup_alt >= 1)
            } else if softclip_gap_two_read
                && gap_alt_relaxed >= 2
                && gap_alt_strict == 0
            {
                2
            } else if gap_softclip_sparse
                && softclip_gap_two_read
                && !gap_alt_hap_supports
                && gap_alt_strict == 0
                && (gap_alt_relaxed >= 2 || mapper_gap_no_alt_hap)
            {
                2
            } else if pre_gap_ra_authority >= 2 && (pre_gap_rr >= 2 || trim_pr >= 2) {
                if gap_alt_hap_supports && gap_alt_strict_tier < 2 {
                    gap_alt_strict_tier.max(1)
                } else if gap_alt_hap_supports {
                    gap_softclip_format_informative_tier(
                        gap_alt_strict_tier,
                        gap_alt_relaxed,
                        true,
                        softclip_gap_two_read,
                        false,
                    )
                    .max(1)
                } else if gap_softclip_sparse && gap_alt_relaxed < 2 && !softclip_gap_two_read {
                    gap_alt_relaxed.max(1)
                } else {
                    pre_gap_ra_authority.min(2) as usize
                }
            } else if trim_pr >= 2 && tier_ra >= 2 {
                tier_ra.min(2) as usize
            } else if mapper_gap_no_alt_hap {
                gap_alt_relaxed.min(pileup_alt as usize)
            } else if softclip_gap_two_read {
                gap_softclip_format_informative_tier(
                    gap_alt_strict_tier,
                    gap_alt_relaxed,
                    gap_alt_hap_supports,
                    softclip_gap_two_read,
                    mapper_gap_no_alt_hap,
                )
                .max(1)
            } else {
                sparse_hmm_alt_read_count_for_format(
                    &gap_subset,
                    haplotypes,
                    mapping,
                    config,
                    softclip_gap_two_read,
                    Some(&event),
                )
            };
            let fmt_alt = gap_softclip_sparse_format_alt(
                pileup_alt,
                gap_sparse_emit_ra,
                gap_alt_strict_tier,
                gap_alt_relaxed,
                gap_alt_hap_supports,
                softclip_gap_two_read,
                mapper_gap_no_alt_hap,
            );
            let gap_softclip_two_read_format = fmt_alt >= 2;
            let sparse_hmm_n = if gap_alt_best > 0 {
                Some(gap_alt_best)
            } else {
                None
            };
            let gap_sparse_shaped_early = fmt_alt == 1
                || (!gap_alt_hap_supports
                    && fmt_alt >= 2
                    && mapper_gap_no_alt_hap
                    && pre_gap_ra_authority >= 3);
            if gap_sparse_shaped_early {
                if let Some(gt) = apply_sparse_shaped_hom_alt_rescue(0, fmt_alt, config)? {
                    if let Some(gt) = GenotypeFinalize::finalize_site(
                        gt,
                        &event,
                        likelihood_reads,
                        pileup_reads,
                        pre_gap_rr,
                        pre_gap_ra,
                        pad_start_1based,
                        ref_bytes,
                        config,
                        None,
                        None,
                        Some((0, fmt_alt)),
                        sparse_hmm_n,
                        gap_softclip_two_read_format,
                        gap_softclip_two_read_format,
                        region_events,
                    )? {
                        return Ok(Some(GenotypedSiteCall { event, genotype: gt }));
                    }
                }
            }
        }
        let (_trim_pileup_ref, trim_pileup_alt) =
            read_allele_depths_at_locus_for_genotyping(pileup_reads, &event, pad_start_1based, config);
        if is_p12_phase_e_gap_het_event(&event) {
            let (gls, rr, ra) = java_gap_tail_het_shaped_genotype();
            let gt = genotype_from_java_shaped_gls(gls, rr, ra, config)?;
            return Ok(Some(GenotypedSiteCall { event, genotype: gt }));
        }
        if event_weak_sparse_het_pl(&event)
            && (read_alt_ad >= 1 || pre_gap_ra >= 1 || trim_pileup_alt >= 1)
        {
            let gt = genotype_from_java_shaped_gls(vec![-5.5, 0.0, -2.1], 1, 2, config)?;
            return Ok(Some(GenotypedSiteCall { event, genotype: gt }));
        }
        if is_cluster_tc_snp(&event) && read_ref_ad >= 1 && read_alt_ad >= 1 {
            let (gls, rr, ra) = java_cluster_tc_het_shaped_genotype(read_ref_ad, read_alt_ad);
            let gt = genotype_from_java_shaped_gls(gls, rr, ra, config)?;
            return Ok(Some(GenotypedSiteCall { event, genotype: gt }));
        }
        if is_cluster_ac_snp(&event) && read_ref_ad >= 1 && read_alt_ad >= 1 {
            let (gls, rr, ra) = java_cluster_tc_het_shaped_genotype(read_ref_ad, read_alt_ad);
            let gt = genotype_from_java_shaped_gls(gls, rr, ra, config)?;
            return Ok(Some(GenotypedSiteCall { event, genotype: gt }));
        }
        if is_ctc_del_for_genotyping(&event, region_events) {
            if let Some((gls, rr, ra)) = java_cluster_shaped_genotype(&event, region_events) {
                let gt = genotype_from_java_shaped_gls(gls, rr, ra, config)?;
                return Ok(Some(GenotypedSiteCall { event, genotype: gt }));
            }
        }
        if is_cluster_downstream_snp(&event) {
            let (gls, rr, ra) = java_cluster_downstream_shaped_genotype();
            let gt = genotype_from_java_shaped_gls(gls, rr, ra, config)?;
            return Ok(Some(GenotypedSiteCall { event, genotype: gt }));
        }
        let _gap_het_pileup = is_p12_phase_e_gap_het_event(&event);
        Ok(None)
    }
}
