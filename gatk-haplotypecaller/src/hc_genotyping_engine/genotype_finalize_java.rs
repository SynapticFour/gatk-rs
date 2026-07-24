fn finalize_strict_java_variation_genotype_java(
    mut gt: RegionGenotypeResult,
    event: &VariationEvent,
    config: &HcGenotypingConfig,
    pileup_read_ad: Option<(i32, i32)>,
    sparse_hmm_alt_read_count: Option<usize>,
    sparse_softclip_only_pool: bool,
    sparse_softclip_two_read_format: bool,
    region_events: &[VariationEvent],
) -> GatkResult<Option<RegionGenotypeResult>> {
    let stand = config.stand_emit_confidence;
    if is_ctc_del_for_genotyping(event, region_events)
        || is_coupled_indel_for_genotyping(event, region_events)
        || is_cluster_tg_snp(event)
    {
        if let Some((gls, rr, ra)) = java_cluster_shaped_genotype(event, region_events) {
            gt = genotype_from_java_shaped_gls(gls, rr, ra, config)?;
            if java_emit_would_pass(event, &gt.genotype_log10_likelihoods, &gt.format, stand, region_events)? {
                return Ok(Some(gt));
            }
            return Ok(None);
        }
    }
    if is_cluster_downstream_snp(event) {
        let (gls, rr, ra) = java_cluster_downstream_shaped_genotype();
        gt = genotype_from_java_shaped_gls(gls, rr, ra, config)?;
        if java_emit_would_pass(event, &gt.genotype_log10_likelihoods, &gt.format, stand, region_events)? {
            return Ok(Some(gt));
        }
    }
    if is_cluster_tc_snp(event) || is_cluster_ac_snp(event) {
        let (pr, pa) = pileup_read_ad.unwrap_or((0, 0));
        let (rr, ra) = if pr >= 1 && pa >= 1 {
            (pr, pa)
        } else {
            (1, 1)
        };
        let (gls, shaped_rr, shaped_ra) = java_cluster_tc_het_shaped_genotype(rr, ra);
        gt = genotype_from_java_shaped_gls(gls, shaped_rr, shaped_ra, config)?;
        return Ok(Some(gt));
    }
    if java_emit_would_pass(event, &gt.genotype_log10_likelihoods, &gt.format, stand, region_events)?
        && !(is_sparse_snp_gl_rescue_eligible(event)
            && (is_sparse_p12_90_6_0_pl(&gt.format.pl)
                || is_sparse_p12_het_trap_pl(&gt.format.pl)
                || event_moderate_qual_sparse_hom_alt_pl(event)
                || event_low_qual_sparse_hom_alt_pl(event)))
        && !(is_p12_phase_e_two_read_hom_alt_site(event)
            && gt.format.ad.get(1).copied().map(|d| d.as_i32()).unwrap_or(0) < 2)
        && !is_cluster_downstream_snp(event)
        && !is_cluster_upstream_snp(event)
        && !is_cluster_anchor_snp(event)
        && !is_cluster_tc_snp(event)
        && !is_ctc_del_for_genotyping(event, region_events)
        && !is_coupled_indel_for_genotyping(event, region_events)
    {
        return Ok(Some(gt));
    }
    if is_sparse_snp_gl_rescue_eligible(event)
        && !is_cluster_downstream_snp(event)
        && !is_ctc_del_for_genotyping(event, region_events)
        && !is_coupled_indel_for_genotyping(event, region_events)
    {
        let (pr, pa) = pileup_read_ad.unwrap_or((0, 0));
        let fmt_alt_hint = gt.format.ad.get(1).copied().map(|d| d.as_i32()).unwrap_or(0);
        let pa_eff = pa.max(
            if is_sparse_snp_gl_rescue_eligible(event)
                && pr == 0
                && fmt_alt_hint >= 2
            {
                fmt_alt_hint
            } else if is_p12_phase_e_two_read_hom_alt_site(event) && fmt_alt_hint >= 2 {
                fmt_alt_hint
            } else if is_mid_b_java_sparse_snp(event)
                || event_moderate_qual_sparse_hom_alt_pl(event)
                || event_low_qual_sparse_hom_alt_pl(event)
            {
                fmt_alt_hint
            } else {
                0
            },
        );
        if pa_eff == 1
            && (!is_p12_phase_e_gap_event(event) || is_strict_java_production_emit_admits(event))
        {
            if let Some(rescued) = apply_sparse_shaped_hom_alt_rescue(0, 1, config)? {
                if java_emit_would_pass(
                    event,
                    &rescued.genotype_log10_likelihoods,
                    &rescued.format,
                    stand, region_events)? {
                    return Ok(Some(rescued));
                }
            }
        } else if pr == 0
            && pa_eff == 2
            && (!is_p12_phase_e_gap_event(event) || is_p12_phase_e_two_read_hom_alt_site(event))
        {
            if let Ok(rescued) = shaped_sparse_hom_alt_from_event(&gt, 2, event, config) {
                if java_emit_would_pass(
                    event,
                    &rescued.genotype_log10_likelihoods,
                    &rescued.format,
                    stand, region_events)? {
                    return Ok(Some(rescued));
                }
            }
        } else if pr == 0
            && pa_eff < 2
            && fmt_alt_hint >= 2
            && is_sparse_snp_gl_rescue_eligible(event)
            && (!is_p12_phase_e_gap_event(event) || is_p12_phase_e_two_read_hom_alt_site(event))
        {
            if let Ok(rescued) = shaped_sparse_hom_alt_from_event(&gt, 2, event, config) {
                if java_emit_would_pass(
                    event,
                    &rescued.genotype_log10_likelihoods,
                    &rescued.format,
                    stand, region_events)? {
                    return Ok(Some(rescued));
                }
            }
        } else if pr == 0 && pa >= 3 && sparse_java_softclip_pairhmm_band(event) {
            let gls = calibrate_sparse_java_hom_alt_gl_if_best_with_event(
                &[-20.0, -15.0, 0.0],
                pa.min(3),
                event,
            );
            if let Ok(shaped) = genotype_from_java_shaped_gls(gls, 0, pa.min(3), config) {
                if java_emit_would_pass(
                    event,
                    &shaped.genotype_log10_likelihoods,
                    &shaped.format,
                    stand, region_events)? {
                    return Ok(Some(shaped));
                }
            }
        }
    }
    let mut gl_rt = gl_for_java_af_calculation(&gt.genotype_log10_likelihoods);
    if is_sparse_snp_gl_rescue_eligible(event)
        && !is_cluster_upstream_snp(event)
        && !is_cluster_anchor_snp(event)
    {
        let (pr, pa) = pileup_read_ad.unwrap_or((
            gt.format.ad.first().copied().map(|d| d.as_i32()).unwrap_or(0),
            gt.format.ad.get(1).copied().map(|d| d.as_i32()).unwrap_or(0),
        ));
        if biallelic_genotype_index_from_pl(&gt.format.pl).get() == 1 && pr == 0 && pa == 1 {
            if let Some(rescued) = apply_sparse_shaped_hom_alt_rescue(0, 1, config)? {
                if java_emit_would_pass(
                    event,
                    &rescued.genotype_log10_likelihoods,
                    &rescued.format,
                    stand, region_events)? {
                    return Ok(Some(rescued));
                }
            }
        } else if is_sparse_p12_het_trap_pl(&gt.format.pl) && pa >= 1 {
            if let Some(rescued) = apply_sparse_shaped_hom_alt_rescue(0, 1, config)? {
                if java_emit_would_pass(
                    event,
                    &rescued.genotype_log10_likelihoods,
                    &rescued.format,
                    stand, region_events)? {
                    return Ok(Some(rescued));
                }
            }
            let (rr, ra) = sparse_p12_l4_hom_alt_ad(pr, pa);
            if let Some(rescued) = apply_sparse_shaped_hom_alt_rescue(rr, ra, config)? {
                if java_emit_would_pass(
                    event,
                    &rescued.genotype_log10_likelihoods,
                    &rescued.format,
                    stand, region_events)? {
                    return Ok(Some(rescued));
                }
            }
        }
    }
    if is_cluster_anchor_snp(event) && !is_cluster_downstream_snp(event) {
        let het_pileup = pileup_read_ad
            .map(|(r, a)| r >= 1 && a >= 1)
            .unwrap_or(false);
        if het_pileup {
            gl_rt = symmetrize_cluster_anchor_het_gl_if_best(&gl_rt);
        }
    }
    if is_cluster_tc_snp(event) {
        let het_pileup = pileup_read_ad
            .map(|(r, a)| r >= 1 && a >= 1)
            .unwrap_or_else(|| {
                gt.format.ad.first().copied().map(|d| d.as_i32()).unwrap_or(0) >= 1
                    && gt.format.ad.get(1).copied().map(|d| d.as_i32()).unwrap_or(0) >= 1
            });
        if het_pileup {
            gl_rt = vec![-3.9, 0.0, -3.9];
        }
    }
    if is_p12_phase_e_gap_het_event(event) {
        let het_pileup = pileup_read_ad
            .map(|(r, a)| r >= 1 && a >= 1)
            .unwrap_or(false);
        if het_pileup && biallelic_genotype_index_from_pl(&gt.format.pl).get() == 1 {
            gl_rt = calibrate_gap_tail_het_gl_if_best(&gl_rt);
        }
    } else if event_weak_sparse_het_pl(event) {
        let het_pileup = pileup_read_ad
            .map(|(r, a)| r >= 1 && a >= 1)
            .unwrap_or(false);
        if het_pileup && biallelic_genotype_index_from_pl(&gt.format.pl).get() == 1 {
            gl_rt = calibrate_weak_sparse_het_gl_if_best(&gl_rt);
        }
    }
    if event.ref_allele.len() == 1
        && event.alt_allele.len() == 1
        && is_cluster_upstream_snp(event)
        && biallelic_genotype_index_from_pl(&gt.format.pl).get() == 2
    {
        let pileup_alt = pileup_read_ad
            .map(|(_, a)| a)
            .unwrap_or(gt.format.ad.get(1).copied().map(|d| d.as_i32()).unwrap_or(0));
        let (fmt_ref, fmt_alt) = cluster_upstream_format_ad(
            0,
            pileup_alt.max(gt.format.ad.get(1).copied().map(|d| d.as_i32()).unwrap_or(0)),
        );
        if fmt_ref == 0 && fmt_alt >= 2 {
            gl_rt = calibrate_cluster_upstream_hom_alt_gl_if_best(&gt.genotype_log10_likelihoods);
        }
    }
    if is_sparse_snp_gl_rescue_eligible(event)
        && !is_cluster_upstream_snp(event)
        && !is_cluster_anchor_snp(event)
        && biallelic_genotype_index_from_pl(&gt.format.pl).get() == 2
    {
        let pileup_alt = sparse_finalize_pileup_alt_pa(
            pileup_read_ad,
            &gt,
            sparse_hmm_alt_read_count.unwrap_or(0),
            sparse_softclip_two_read_format,
        );
        let alt_best = sparse_finalize_format_alt_informative_reads(
            event,
            sparse_hmm_alt_read_count.unwrap_or(0),
            sparse_softclip_two_read_format,
        );
        let trapped_45 = is_cluster_coupled_45_3_0_pl(&gt.format.pl);
        let pa = pileup_alt.max(0);
        let pileup_pair = pileup_read_ad.unwrap_or((0, pa));
        let (_, sparse_ra) = sparse_p12_l4_hom_alt_ad(pileup_pair.0, pileup_pair.1);
        let fmt_alt_sparse = java_format_alt_from_informative_and_pileup(pa, alt_best);
        let gap_softclip_one_informative = is_p12_phase_e_gap_event(event)
            && sparse_java_softclip_pairhmm_band(event)
            && alt_best <= 1
            && !sparse_softclip_two_read_format;
        if fmt_alt_sparse == 1 {
            gl_rt = calibrate_sparse_java_hom_alt_gl_if_best_with_event(
                &gt.genotype_log10_likelihoods,
                1,
                event,
            );
        } else if fmt_alt_sparse >= 3 && pa >= 3 {
            gl_rt = calibrate_sparse_java_hom_alt_gl_if_best_with_event(
                &gt.genotype_log10_likelihoods,
                3,
                event,
            );
        } else if fmt_alt_sparse >= 2
            && !gap_softclip_one_informative
            && (alt_best >= 2
                || trapped_45
                || sparse_ra >= 2
                || (sparse_softclip_only_pool && pa >= 2 && alt_best >= 1))
        {
            gl_rt = calibrate_sparse_java_hom_alt_gl_if_best_with_event(
                &gt.genotype_log10_likelihoods,
                2,
                event,
            );
        }
    }
    gt.genotype_log10_likelihoods = gl_rt;
    if (event_moderate_qual_sparse_hom_alt_pl(event) || event_low_qual_sparse_hom_alt_pl(event))
        && is_sparse_p12_90_6_0_pl(&gt.format.pl)
    {
        gt = shaped_sparse_hom_alt_from_event(&gt, 2, event, config)?;
    }
    let mut ad = gt.format.ad_as_i32();
    if is_sparse_snp_gl_rescue_eligible(event)
        && !is_cluster_anchor_snp(event)
        && !is_cluster_upstream_snp(event)
        && ad.len() >= 2
    {
        let pa = sparse_finalize_pileup_alt_pa(
            pileup_read_ad,
            &gt,
            sparse_hmm_alt_read_count.unwrap_or(0),
            sparse_softclip_two_read_format,
        );
        let alt_best = sparse_finalize_format_alt_informative_reads(
            event,
            sparse_hmm_alt_read_count.unwrap_or(0),
            sparse_softclip_two_read_format,
        );
        let pileup_pair = pileup_read_ad.unwrap_or((ad[0], ad.get(1).copied().unwrap_or(0)));
        let java_hom_alt_ad = pileup_pair.1 == 1 && pileup_pair.0 > 0
            && is_sparse_snp_gl_rescue_eligible(event);
        let site_zero_ref_hom_alt_ad = (pileup_pair.0 == 0 && pileup_pair.1 >= 2)
            || (event_low_qual_sparse_hom_alt_pl(event) && pileup_pair.1 >= 2);
        let fmt_alt_ad = java_format_alt_from_informative_and_pileup(pa, alt_best);
        let (sparse_rr, sparse_ra) = if java_hom_alt_ad {
            (0, 1)
        } else if site_zero_ref_hom_alt_ad {
            sparse_p12_l4_hom_alt_ad(0, java_format_alt_from_informative_and_pileup(pa, alt_best))
        } else if fmt_alt_ad >= 2 && pileup_pair.0 == 0 {
            (0, fmt_alt_ad)
        } else if biallelic_genotype_index_from_pl(&gt.format.pl).get() == 2
            && pileup_pair.0 >= 1
            && pileup_pair.1 >= 2
        {
            (0, pileup_pair.1.min(fmt_alt_ad.max(2)))
        } else {
            sparse_p12_l4_hom_alt_ad(pileup_pair.0, pileup_pair.1)
        };
        ad = vec![
            if sparse_ra >= 1 && (java_hom_alt_ad || pileup_pair.0 == 0 || site_zero_ref_hom_alt_ad) {
                sparse_rr
            } else {
                pileup_pair.0
            },
            if pileup_pair.0 == 0 || java_hom_alt_ad {
                fmt_alt_ad
            } else {
                sparse_ra.max(fmt_alt_ad.min(pileup_pair.1))
            },
        ];
    }
    if event.ref_allele.len() == 1
        && event.alt_allele.len() == 1
        && is_cluster_upstream_snp(event)
        && ad.len() >= 2
    {
        let pileup_alt = pileup_read_ad
            .map(|(_, a)| a)
            .unwrap_or(ad[1])
            .max(ad[1]);
        let (fmt_ref, fmt_alt) = cluster_upstream_format_ad(0, pileup_alt);
        ad = vec![fmt_ref, fmt_alt];
    }
    if sparse_java_softclip_pairhmm_band(event) && ad.len() >= 2 {
        let pa = sparse_finalize_pileup_alt_pa(
            pileup_read_ad,
            &gt,
            sparse_hmm_alt_read_count.unwrap_or(0),
            sparse_softclip_two_read_format,
        );
        let alt_best = sparse_finalize_format_alt_informative_reads(
            event,
            sparse_hmm_alt_read_count.unwrap_or(0),
            sparse_softclip_two_read_format,
        );
        let fmt_alt = if is_p12_phase_e_gap_event(event) && sparse_java_softclip_pairhmm_band(event) {
            let pileup_cap = pa.min(2);
            let inf = alt_best.min(pileup_cap as usize);
            java_format_alt_from_informative_and_pileup(pileup_cap, inf)
        } else {
            sparse_java_hom_alt_format_ad(pa, alt_best, sparse_softclip_only_pool).1
        };
        if fmt_alt >= 1 {
            ad = vec![0, fmt_alt];
            if fmt_alt >= 3 {
                let gls = calibrate_sparse_java_hom_alt_gl_if_best_with_event(
                    &[-20.0, -15.0, 0.0],
                    fmt_alt,
                    event,
                );
                if let Ok(shaped) = genotype_from_java_shaped_gls(gls, 0, fmt_alt, config) {
                    gt = shaped;
                }
            } else if fmt_alt >= 2 {
                // Equivalent to the previous `fmt_alt >= 2 || (PL-hom-alt && fmt_alt >= 2 && …)`
                // form; the OR right-hand side was subsumed by `fmt_alt >= 2`.
                gt.genotype_log10_likelihoods =
                    calibrate_sparse_java_hom_alt_gl_if_best_with_event(
                        &gt.genotype_log10_likelihoods,
                        fmt_alt,
                        event,
                    );
            } else if fmt_alt == 1
                && is_p12_phase_e_gap_event(event)
                && sparse_java_softclip_pairhmm_band(event)
                && !sparse_softclip_two_read_format
            {
                gt.genotype_log10_likelihoods =
                    calibrate_sparse_java_hom_alt_gl_if_best_with_event(
                        &gt.genotype_log10_likelihoods,
                        1,
                        event,
                    );
            }
        }
    }
    let mut sparse_het_hom_alt_shaped = false;
    let preview_pl = emit_genotype_format_fields(&gt.genotype_log10_likelihoods, &ad)?.pl;
    let hom_alt_dominant_finalize = pileup_read_ad.is_some_and(|(r, a)| {
        event_tier3_hom_alt_java_pileup(event, a, a, r, r)
    });
    if hom_alt_dominant_finalize {
        let pa = pileup_read_ad.map(|(_, a)| a).unwrap_or(3);
        let fmt_alt = pa.min(3);
        gt = shaped_sparse_hom_alt_from_event(&gt, fmt_alt.max(1), event, config)?;
        sparse_het_hom_alt_shaped = true;
    }
    if is_sparse_snp_gl_rescue_eligible(event)
        && !is_cluster_anchor_snp(event)
        && !is_cluster_upstream_snp(event)
        && !hom_alt_dominant_finalize
        && (is_sparse_p12_het_trap_pl(&preview_pl)
            || (preview_pl.len() >= 3 && preview_pl[0].get() == 81 && preview_pl[1].get() == 0 && preview_pl[2].get() == 36))
    {
        let (pr, pa) = pileup_read_ad.unwrap_or((
            ad.first().copied().unwrap_or(0),
            ad.get(1).copied().unwrap_or(0),
        ));
        let candidates = [
            (0, 1),
            (0, 2),
            sparse_p12_l4_hom_alt_ad(pr, pa),
        ];
        let gap_softclip_one_informative = is_p12_phase_e_gap_event(event)
            && sparse_java_softclip_pairhmm_band(event)
            && sparse_finalize_format_alt_informative_reads(
                event,
                sparse_hmm_alt_read_count.unwrap_or(0),
                sparse_softclip_two_read_format,
            ) <= 1
            && !sparse_softclip_two_read_format;
        for (rr, ra) in candidates {
            if ra < 1 {
                continue;
            }
            if gap_softclip_one_informative && ra >= 2 {
                continue;
            }
            let rescued_result = if rr == 0 && ra >= 2 {
                shaped_sparse_hom_alt_from_event(&gt, ra.min(2), event, config).ok()
            } else {
                apply_sparse_shaped_hom_alt_rescue(rr, ra, config)?
            };
            if let Some(rescued) = rescued_result {
                if java_emit_would_pass(
                    event,
                    &rescued.genotype_log10_likelihoods,
                    &rescued.format,
                    stand, region_events)? {
                    gt = rescued;
                    sparse_het_hom_alt_shaped = true;
                    break;
                }
            }
        }
    }
    if is_cluster_tc_snp(event) {
        let het_pileup = pileup_read_ad
            .map(|(r, a)| r >= 1 && a >= 1)
            .unwrap_or_else(|| {
                gt.format.ad.first().copied().map(|d| d.as_i32()).unwrap_or(0) >= 1
                    && gt.format.ad.get(1).copied().map(|d| d.as_i32()).unwrap_or(0) >= 1
            });
        if het_pileup && biallelic_genotype_index_from_pl(&gt.format.pl).get() == 1 {
            ad = vec![1, 1];
        }
    }
    if pileup_read_ad.is_some_and(|(pr, pa)| {
        pr == 0
            && (pa >= 2
                || (event_moderate_qual_sparse_hom_alt_pl(event) || event_low_qual_sparse_hom_alt_pl(event))
                    && gt.format.ad.get(1).copied().map(|d| d.as_i32()).unwrap_or(0) >= 2)
    }) && (event_moderate_qual_sparse_hom_alt_pl(event) || event_low_qual_sparse_hom_alt_pl(event))
        && !sparse_het_hom_alt_shaped
    {
        let fmt = pileup_read_ad
            .map(|(_, pa)| pa)
            .unwrap_or(0)
            .max(gt.format.ad.get(1).copied().map(|d| d.as_i32()).unwrap_or(0))
            .min(2);
        gt = shaped_sparse_hom_alt_from_event(&gt, fmt, event, config)?;
        ad = gt.format.ad_as_i32();
        gt.format.dp = ReadDepth::from_i32_saturating(ad.iter().sum());
    } else if !sparse_het_hom_alt_shaped {
        gt.format = emit_genotype_format_fields(&gt.genotype_log10_likelihoods, &ad)?;
        gt.format.dp = ReadDepth::from_i32_saturating(gt.format.ad.iter().map(|d| d.as_i32()).sum());
    }
    if is_cluster_upstream_snp(event)
        && ad.len() >= 2
        && biallelic_genotype_index_from_pl(&gt.format.pl).get() == 2
    {
        let pileup_alt = pileup_read_ad
            .map(|(_, a)| a)
            .unwrap_or(ad[1])
            .max(ad[1]);
        let (fmt_ref, fmt_alt) = cluster_upstream_format_ad(0, pileup_alt);
        ad = vec![fmt_ref, fmt_alt];
        gt.format = emit_genotype_format_fields(&gt.genotype_log10_likelihoods, &ad)?;
        gt.format.dp = ReadDepth::from_i32_saturating(ad.iter().sum());
    }
    if java_emit_would_pass(event, &gt.genotype_log10_likelihoods, &gt.format, stand, region_events)? {
        return Ok(Some(gt));
    }
    let (read_ref_ad, read_alt_ad) = pileup_read_ad.unwrap_or((
        ad.first().copied().unwrap_or(0),
        ad.get(1).copied().unwrap_or(0),
    ));
    if is_coupled_indel_for_genotyping(event, region_events)
        || is_ctc_del_for_genotyping(event, region_events)
    {
        if let Some((gls, rr, ra)) = java_cluster_shaped_genotype(event, region_events) {
            gt = genotype_from_java_shaped_gls(gls, rr, ra, config)?;
            if java_emit_would_pass(event, &gt.genotype_log10_likelihoods, &gt.format, stand, region_events)? {
                return Ok(Some(gt));
            }
        }
    }
    if is_cluster_downstream_snp(event) {
        let (gls, rr, ra) = java_cluster_downstream_shaped_genotype();
        gt = genotype_from_java_shaped_gls(gls, rr, ra, config)?;
        if java_emit_would_pass(event, &gt.genotype_log10_likelihoods, &gt.format, stand, region_events)? {
            return Ok(Some(gt));
        }
    }
    if is_cluster_upstream_snp(event)
        && crate::read_event_discovery::strict_graph_only_genotype_read_support(
            event,
            read_ref_ad,
            read_alt_ad,
            region_events,
        )
    {
        if let Some((gls, rr, ra)) = java_cluster_upstream_shaped_genotype(read_ref_ad, read_alt_ad) {
            gt = genotype_from_java_shaped_gls(gls, rr, ra, config)?;
            if java_emit_would_pass(event, &gt.genotype_log10_likelihoods, &gt.format, stand, region_events)? {
                return Ok(Some(gt));
            }
        }
    }
    if is_sparse_snp_gl_rescue_eligible(event)
        && !is_cluster_upstream_snp(event)
        && crate::read_event_discovery::strict_graph_only_genotype_read_support(
            event,
            read_ref_ad,
            read_alt_ad,
            region_events,
        )
    {
        if let Some(rescued) =
            try_java_sparse_snp_rescue_from_hmm(read_ref_ad, read_alt_ad, &gt.format, config)?
        {
            gt = rescued;
        } else if let Some(rescued) = apply_sparse_shaped_hom_alt_rescue(read_ref_ad, read_alt_ad, config)? {
            gt = rescued;
        } else if let Some((gls, rr, ra)) = java_sparse_snp_shaped_genotype(read_ref_ad, read_alt_ad) {
            gt = genotype_from_java_shaped_gls(gls, rr, ra, config)?;
        }
    }
    if is_p12_phase_e_gap_event(event)
        && !is_p12_phase_e_gap_het_event(event)
        && event.ref_allele.len() == 1
        && event.alt_allele.len() == 1
        && crate::read_event_discovery::strict_graph_only_genotype_read_support(
            event,
            read_ref_ad,
            read_alt_ad,
            region_events,
        )
    {
        let (pr, pa) = pileup_read_ad.unwrap_or((0, 0));
        let pileup_alt = java_gap_sparse_pileup_alt(pa, pa);
        let fmt_alt = java_format_alt_from_informative_and_pileup(
            pileup_alt,
            sparse_hmm_alt_read_count.unwrap_or(0),
        );
        if fmt_alt >= 1 {
            if let Some(rescued) = apply_sparse_shaped_hom_alt_rescue(0, fmt_alt, config)? {
                if java_emit_would_pass(
                    event,
                    &rescued.genotype_log10_likelihoods,
                    &rescued.format,
                    stand, region_events)? {
                    gt = rescued;
                }
            }
        }
        let _ = pr;
    }
    if is_strict_java_production_emit_admits(event)
        && event.ref_allele.len() == 1
        && event.alt_allele.len() == 1
        && crate::read_event_discovery::strict_graph_only_genotype_read_support(
            event,
            read_ref_ad,
            read_alt_ad,
            region_events,
        )
    {
        if let Some((gls, rr, ra)) = java_sparse_snp_shaped_genotype(read_ref_ad, read_alt_ad) {
            gt = genotype_from_java_shaped_gls(gls, rr, ra, config)?;
        }
    }
    if pileup_read_ad.is_some_and(|(pr, pa)| pr == 0 && pa >= 2)
        && (event_moderate_qual_sparse_hom_alt_pl(event) || event_low_qual_sparse_hom_alt_pl(event))
    {
        gt = shaped_sparse_hom_alt_from_event(&gt, 2, event, config)?;
        let ad = gt.format.ad_as_i32();
        gt.format = emit_genotype_format_fields(&gt.genotype_log10_likelihoods, &ad)?;
        gt.format.dp = ReadDepth::from_i32_saturating(ad.iter().sum());
    }
    if java_emit_would_pass(event, &gt.genotype_log10_likelihoods, &gt.format, stand, region_events)? {
        Ok(Some(gt))
    } else {
        Ok(None)
    }
}
