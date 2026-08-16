
/// L11-B: `post_finalize_strict_java_call` is a no-op on production `enable_java_strict`.
#[inline]
fn maybe_post_finalize_strict_java_call(
    call: GenotypedSiteCall,
    pileup_reads: &[SharedBamRecord],
    supplemental_pileup_reads: Option<&[SharedBamRecord]>,
    pad_start_1based: u64,
    ref_bases: &[u8],
    config: &HcGenotypingConfig,
) -> GatkResult<GenotypedSiteCall> {
    if config.enable_java_strict() {
        return Ok(call);
    }
    post_finalize_strict_java_call(
        call,
        pileup_reads,
        supplemental_pileup_reads,
        pad_start_1based,
        ref_bases,
        config,
    )
}

/// GATK `HaplotypeCallerGenotypingEngine.assignGenotypeLikelihoods` (active-window event walk).
pub fn assign_genotype_likelihoods_for_region(
    likelihoods: &[RegionReadLikelihood],
    likelihood_reads: &[SharedBamRecord],
    pileup_reads: &[SharedBamRecord],
    supplemental_pileup_reads: Option<&[SharedBamRecord]>,
    haplotypes: &[Haplotype],
    ref_bytes: &[u8],
    pad_start_1based: u64,
    full_reference_bases: &[u8],
    full_reference_pad_1based: u64,
    active_start_1based: u64,
    active_end_1based: u64,
    contig: &str,
    max_mnp_distance: usize,
    config: &HcGenotypingConfig,
    stored_events: &[VariationEvent],
    _graph_events: &[VariationEvent],
) -> GatkResult<AssignGenotypeLikelihoodsResult> {
    if haplotypes.is_empty() {
        return Err(gatk_common::GatkError::algorithm(
            "assignGenotypeLikelihoods: haplotype list is empty",
        ));
    }
    // Multi-pass pileup AD in try_genotype reuses CIGAR/seq decode across events.
    crate::read_event_discovery::clear_ad_decode_cache();
    let rows = region_likelihoods_to_rows(likelihoods, haplotypes.len());
    let region_summary = if rows.is_empty() {
        sparse_snp_genotype_from_read_depths(0, 0, config)?
    } else {
        genotype_from_read_rows(&rows, haplotypes, config)?
    };
    let emit_spanning = !config.disable_spanning_event_genotyping;
    let hap_events = build_per_haplotype_variation_events(
        haplotypes,
        ref_bytes,
        pad_start_1based,
        max_mnp_distance,
        contig,
    );
    let supplement_events = if config.enable_java_strict() {
        stored_events_with_p12_cluster_anchors(
            stored_events,
            ref_bytes,
            pad_start_1based,
            active_start_1based,
            active_end_1based,
            contig,
            config,
        )
    } else {
        stored_events.to_vec()
    };

    if config.genotype_stored_events_only {
        let mut calls = Vec::new();
        let stored_events = stored_events_with_p12_cluster_anchors(
            stored_events,
            ref_bytes,
            pad_start_1based,
            active_start_1based,
            active_end_1based,
            contig,
            config,
        );
        for event in &stored_events {
            if event.start_1based < GenomePosition::new_1based(active_start_1based) || event.start_1based > GenomePosition::new_1based(active_end_1based) {
                continue;
            }
            if let Some(call) = try_genotype_variation_event(
                // CLONE: needed because try_genotype_variation_event currently takes owned VariationEvent.
                event.clone(),
                likelihoods,
                likelihood_reads,
                pileup_reads,
                supplemental_pileup_reads,
                haplotypes,
                ref_bytes,
                pad_start_1based,
                full_reference_bases,
                full_reference_pad_1based,
                active_start_1based,
                active_end_1based,
                max_mnp_distance,
                config,
                &stored_events,
                Some(&hap_events),
            )? {
                let call = maybe_post_finalize_strict_java_call(
                    call,
                    pileup_reads,
                    supplemental_pileup_reads,
                    pad_start_1based,
                    ref_bytes,
                    config,
                )?;
                if !event_already_called(&calls, &call.event) {
                    calls.push(call);
                }
            }
        }
        return Ok(AssignGenotypeLikelihoodsResult {
            calls,
            region_summary,
        });
    }

    let mut positions = build_event_start_positions_from_cache(&hap_events);
    for e in &supplement_events {
        if e.start_1based >= GenomePosition::new_1based(active_start_1based) && e.start_1based <= GenomePosition::new_1based(active_end_1based) {
            positions.insert(e.start_1based.get());
        }
    }

    let mut calls = Vec::new();
    for loc in positions {
        if loc < active_start_1based || loc > active_end_1based {
            continue;
        }
        let mut raw = variation_events_at_position_from_cache(&hap_events, loc, emit_spanning);
        merge_stored_variation_events_at_position(&mut raw, &supplement_events, loc, emit_spanning);
        let with_span = replace_span_del_events(&raw, loc, pad_start_1based, ref_bytes);
        let mut merged = merged_biallelic_sites_at_position(&with_span, loc);
        merged.retain(|e| e.ref_allele != e.alt_allele && e.alt_allele != SPAN_DEL_ALLELE);
        // Non-empty checked at function entry (K-3).
        let ref_hap = haplotypes
            .iter()
            .find(|h| h.is_reference)
            .or_else(|| haplotypes.first())
            .expect("non-empty haplotypes checked at entry");
        merged.sort_by(|a, b| {
            let sa = haplotypes
                .iter()
                .enumerate()
                .filter(|(hi, h)| {
                    hap_supports_allele_for_sort(
                        h,
                        *hi,
                        &hap_events,
                        ref_hap,
                        loc,
                        pad_start_1based,
                        &a.ref_allele,
                        &a.alt_allele,
                        ref_bytes,
                        max_mnp_distance,
                        contig,
                    )
                })
                .map(|(_, h)| h.score)
                .fold(f64::NEG_INFINITY, f64::max);
            let sb = haplotypes
                .iter()
                .enumerate()
                .filter(|(hi, h)| {
                    hap_supports_allele_for_sort(
                        h,
                        *hi,
                        &hap_events,
                        ref_hap,
                        loc,
                        pad_start_1based,
                        &b.ref_allele,
                        &b.alt_allele,
                        ref_bytes,
                        max_mnp_distance,
                        contig,
                    )
                })
                .map(|(_, h)| h.score)
                .fold(f64::NEG_INFINITY, f64::max);
            sb.total_cmp(&sa)
        });
        if !config.enable_java_strict() && merged.len() > 3 {
            merged.truncate(3);
        }

        for event in merged {
            if let Some(call) = try_genotype_variation_event(
                event,
                likelihoods,
                likelihood_reads,
                pileup_reads,
                supplemental_pileup_reads,
                haplotypes,
                ref_bytes,
                pad_start_1based,
                full_reference_bases,
                full_reference_pad_1based,
                active_start_1based,
                active_end_1based,
                max_mnp_distance,
                config,
                &supplement_events,
                Some(&hap_events),
            )? {
                let call = maybe_post_finalize_strict_java_call(
                    call,
                    pileup_reads,
                    supplemental_pileup_reads,
                    pad_start_1based,
                    ref_bytes,
                    config,
                )?;
                calls.push(call);
            }
        }
        if config.enable_java_strict() {
            for event in &supplement_events {
                if event.start_1based != GenomePosition::new_1based(loc) || !is_cluster_anchor_snp(event) {
                    continue;
                }
                if event_already_called(&calls, event) {
                    continue;
                }
                if let Some(call) = try_genotype_variation_event(
                    // CLONE: needed because try_genotype_variation_event takes owned VariationEvent.
                    event.clone(),
                    likelihoods,
                    likelihood_reads,
                    pileup_reads,
                    supplemental_pileup_reads,
                    haplotypes,
                    ref_bytes,
                    pad_start_1based,
                    full_reference_bases,
                    full_reference_pad_1based,
                    active_start_1based,
                    active_end_1based,
                    max_mnp_distance,
                    config,
                    &supplement_events,
                    Some(&hap_events),
                )? {
                    let call = maybe_post_finalize_strict_java_call(
                        call,
                        pileup_reads,
                        supplemental_pileup_reads,
                        pad_start_1based,
                        ref_bytes,
                        config,
                    )?;
                    calls.push(call);
                }
            }
        }
    }

    if config.enable_java_strict() {
        for event in &supplement_events {
            if event.start_1based < GenomePosition::new_1based(active_start_1based) || event.start_1based > GenomePosition::new_1based(active_end_1based) {
                continue;
            }
            if event_already_called(&calls, event) {
                continue;
            }
            if let Some(call) = try_genotype_variation_event(
                // CLONE: needed because try_genotype_variation_event takes owned VariationEvent.
                event.clone(),
                likelihoods,
                likelihood_reads,
                pileup_reads,
                supplemental_pileup_reads,
                haplotypes,
                ref_bytes,
                pad_start_1based,
                full_reference_bases,
                full_reference_pad_1based,
                active_start_1based,
                active_end_1based,
                max_mnp_distance,
                config,
                &supplement_events,
                Some(&hap_events),
            )? {
                let call = maybe_post_finalize_strict_java_call(
                    call,
                    pileup_reads,
                    supplemental_pileup_reads,
                    pad_start_1based,
                    ref_bytes,
                    config,
                )?;
                calls.push(call);
            }
        }
    }

    if !config.enable_java_strict() {
        for event in stored_events {
            if event.start_1based < GenomePosition::new_1based(active_start_1based) || event.start_1based > GenomePosition::new_1based(active_end_1based) {
                continue;
            }
            if event_already_called(&calls, event) {
                continue;
            }
            if let Some(call) = try_genotype_variation_event(
                // CLONE: needed because try_genotype_variation_event takes owned VariationEvent.
                event.clone(),
                likelihoods,
                likelihood_reads,
                pileup_reads,
                supplemental_pileup_reads,
                haplotypes,
                ref_bytes,
                pad_start_1based,
                full_reference_bases,
                full_reference_pad_1based,
                active_start_1based,
                active_end_1based,
                max_mnp_distance,
                config,
                stored_events,
                Some(&hap_events),
            )? {
                let call = maybe_post_finalize_strict_java_call(
                    call,
                    pileup_reads,
                    supplemental_pileup_reads,
                    pad_start_1based,
                    ref_bytes,
                    config,
                )?;
                calls.push(call);
            }
        }
    }

    Ok(AssignGenotypeLikelihoodsResult {
        calls,
        region_summary,
    })
}

/// Sort-key allele support using the region EventMap cache when possible.
///
/// Observable contract: same truth as [`haplotype_supports_allele_at_with_ref`]; avoid
/// rebuilding EventMaps O(alleles × haps) on the dense position walk.
fn hap_supports_allele_for_sort(
    hap: &Haplotype,
    hap_index: usize,
    hap_events: &crate::event_map::PerHaplotypeVariationEvents,
    ref_hap: &Haplotype,
    loc_1based: u64,
    pad_start: u64,
    ref_allele: &str,
    alt_allele: &str,
    ref_bytes: &[u8],
    max_mnp_distance: usize,
    contig: &str,
) -> bool {
    if alt_allele == SPAN_DEL_ALLELE {
        return !hap.is_reference;
    }
    if ref_allele.len() == 1 && alt_allele.len() == 1 {
        let Some(alt_byte) = alt_allele
            .as_bytes()
            .first()
            .map(|b| b.to_ascii_uppercase())
        else {
            return false;
        };
        return hap_base_at_ref_locus(hap, pad_start, loc_1based)
            .map(|b| b.to_ascii_uppercase() == alt_byte)
            .unwrap_or(false);
    }
    if cached_events_support_allele_at(
        hap_events.events_for(hap_index),
        loc_1based,
        ref_allele,
        alt_allele,
    ) {
        return true;
    }
    crate::hc_allele_mapping::haplotype_supports_allele_at_with_events(
        hap,
        ref_hap,
        loc_1based,
        pad_start,
        ref_allele,
        alt_allele,
        ref_bytes,
        max_mnp_distance,
        contig,
        Some(hap_events.events_for(hap_index)),
    )
}

fn merge_stored_variation_events_at_position(
    from_haps: &mut Vec<VariationEvent>,
    stored: &[VariationEvent],
    loc_1based: u64,
    include_spanning: bool,
) {
    let mut seen: BTreeSet<(u64, String, String)> = from_haps
        .iter()
        .map(|e| (e.start_1based.get(), e.ref_allele.clone(), e.alt_allele.clone()))
        .collect();
    for e in stored {
        let overlaps = e.end_1based >= GenomePosition::new_1based(loc_1based) && e.start_1based <= GenomePosition::new_1based(loc_1based);
        if !overlaps {
            continue;
        }
        if !include_spanning && e.start_1based != GenomePosition::new_1based(loc_1based) {
            continue;
        }
        // CLONE: needed because owned composite key for dedup/lookup.
        let key = (e.start_1based.get(), e.ref_allele.clone(), e.alt_allele.clone());
        if seen.insert(key) {
            // CLONE: needed because owned element into collection.
            from_haps.push(e.clone());
        }
    }
    crate::event_map::prefer_indel_over_colocated_snps(from_haps);
}

/// Biallelic genotyping for an explicit REF/ALT haplotype pair (per variation event).
pub fn genotype_biallelic_indices(
    likelihoods: &[RegionReadLikelihood],
    haplotypes: &[Haplotype],
    ref_idx: usize,
    alt_idx: usize,
    config: &HcGenotypingConfig,
) -> GatkResult<RegionGenotypeResult> {
    let rows = region_likelihoods_to_rows(likelihoods, haplotypes.len());
    let aggregation = aggregate_haplotype_log10_likelihoods(&rows)?;
    let best = best_haplotype_index(&aggregation).unwrap_or(crate::bio_ids::HaplotypeIndex::new(0)).get();
    let gls = biallelic_genotype_log10_likelihoods_gatk(&rows, ref_idx, alt_idx);
    let depths = biallelic_allele_depths_from_rows(&rows, ref_idx, alt_idx);
    let priors = biallelic_diploid_log10_priors(config.priors)?;
    let _posterior = genotype_posteriors_from_log10_likelihoods(&gls, &priors)?;
    let format = emit_genotype_format_fields(&gls, &depths)?;
    Ok(RegionGenotypeResult {
        aggregation,
        best_haplotype_index: best,
        ref_haplotype_index: ref_idx,
        alt_haplotype_index: alt_idx,
        genotype_log10_likelihoods: gls,
        format,
    })
}

