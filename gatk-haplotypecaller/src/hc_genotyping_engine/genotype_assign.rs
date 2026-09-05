
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
    let _prof = crate::hc_profile::begin(crate::hc_profile::Stage::GenotypeAssignment);
    crate::hc_profile::reset_genotype_nested_walls();
    let geno_wall0 = std::time::Instant::now();
    COLOCATED_MERGE_NUMERICS.with(|slot| slot.borrow_mut().clear());
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
    // Java `EventMap.buildEventMapsForHaplotypes(haplotypes, ref, refLoc)` uses
    // `assemblyResult.getFullReferenceWithPadding()` / `getPaddedReferenceLoc()`.
    // `alignmentStartHapwrtRef` is an offset into that padded array. Rebuilding
    // against the trimmed apply window shifts event starts (6R.67: SNP T/C and
    // deletion TG/T at 20:29456344 vanished from the genotyping EventMap; merge
    // never saw both alleles).
    let hap_events = {
        let t0 = crate::hc_profile::enabled().then(std::time::Instant::now);
        let ev = build_per_haplotype_variation_events(
            haplotypes,
            full_reference_bases,
            full_reference_pad_1based,
            max_mnp_distance,
            contig,
        );
        if let Some(t0) = t0 {
            crate::hc_profile::note_event_rebuild_wall(t0.elapsed());
        }
        ev
    };
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
        return finish_assign_genotype_profile(calls, region_summary, geno_wall0);
    }

    let mut positions = build_event_start_positions_from_cache(&hap_events);
    for e in &supplement_events {
        if e.start_1based >= GenomePosition::new_1based(active_start_1based) && e.start_1based <= GenomePosition::new_1based(active_end_1based) {
            positions.insert(e.start_1based.get());
        }
    }

    let mut calls = Vec::new();
    let mut merged_handled_locs: HashSet<u64> = HashSet::new();
    for loc in positions {
        if loc < active_start_1based || loc > active_end_1based {
            continue;
        }
        let mut raw = variation_events_at_position_from_cache(&hap_events, loc, emit_spanning);
        merge_stored_variation_events_at_position(&mut raw, &supplement_events, loc, emit_spanning);
        let with_span_pre = replace_span_del_events(&raw, loc, pad_start_1based, ref_bytes);
        let merged_site_handled = match try_genotype_colocated_snp_indel_merge(
            &with_span_pre,
            loc,
            likelihoods,
            likelihood_reads,
            haplotypes,
            ref_bytes,
            pad_start_1based,
            active_start_1based,
            active_end_1based,
            max_mnp_distance,
            config,
            Some(&hap_events),
        )? {
            ColocatedMergeGenotype::Call(call) => {
                merged_handled_locs.insert(loc);
                calls.push(call);
                true
            }
            ColocatedMergeGenotype::MergedNoEmit => {
                merged_handled_locs.insert(loc);
                true
            }
            ColocatedMergeGenotype::NotApplicable => false,
        };
        if !merged_site_handled {
            crate::event_map::prefer_indel_over_colocated_snps(&mut raw);
        let with_span = replace_span_del_events(&raw, loc, pad_start_1based, ref_bytes);
        let mut merged = merged_biallelic_sites_at_position(&with_span, loc);
        merged.retain(|e| e.ref_allele != e.alt_allele && e.alt_allele != SPAN_DEL_ALLELE);
        // Non-empty checked at function entry (K-3).
        let ref_hap = haplotypes
            .iter()
            .find(|h| h.is_reference)
            .or_else(|| haplotypes.first())
            .expect("non-empty haplotypes checked at entry");
        // Precompute SNP hap bases once per locus — sort otherwise rewalks CIGAR O(A log A · H).
        let hap_snp_bases: Vec<Option<u8>> = haplotypes
            .iter()
            .map(|h| {
                hap_base_at_ref_locus(h, pad_start_1based, loc).map(|b| b.to_ascii_uppercase())
            })
            .collect();
        merged.sort_by(|a, b| {
            let sa = haplotypes
                .iter()
                .enumerate()
                .filter(|(hi, h)| {
                    hap_supports_allele_for_sort_cached(
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
                        &hap_snp_bases,
                    )
                })
                .map(|(_, h)| h.score)
                .fold(f64::NEG_INFINITY, f64::max);
            let sb = haplotypes
                .iter()
                .enumerate()
                .filter(|(hi, h)| {
                    hap_supports_allele_for_sort_cached(
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
                        &hap_snp_bases,
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
        }
        if config.enable_java_strict() {
            for event in &supplement_events {
                if event.start_1based != GenomePosition::new_1based(loc) || !is_cluster_anchor_snp(event) {
                    continue;
                }
                if event_already_called(&calls, event)
                    || merged_handled_locs.contains(&loc)
                {
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
            if event_already_called(&calls, event)
                || merged_handled_locs.contains(&event.start_1based.get())
            {
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
            if event_already_called(&calls, event)
                || merged_handled_locs.contains(&event.start_1based.get())
            {
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

    finish_assign_genotype_profile(calls, region_summary, geno_wall0)
}

fn finish_assign_genotype_profile(
    calls: Vec<GenotypedSiteCall>,
    region_summary: RegionGenotypeResult,
    geno_wall0: std::time::Instant,
) -> GatkResult<AssignGenotypeLikelihoodsResult> {
    record_genotype_profile_samples(&calls, geno_wall0.elapsed());
    Ok(AssignGenotypeLikelihoodsResult {
        calls,
        region_summary,
    })
}

fn record_genotype_profile_samples(
    calls: &[GenotypedSiteCall],
    region_wall: std::time::Duration,
) {
    if !crate::hc_profile::enabled() {
        return;
    }
    let ad_ns = crate::hc_profile::take_ad_wall_ns();
    let event_ns = crate::hc_profile::take_event_rebuild_wall_ns();
    let map_ns = crate::hc_profile::take_allele_map_wall_ns();
    let marg_ns = crate::hc_profile::take_marginalize_wall_ns();
    let enum_ns = crate::hc_profile::take_genotype_enum_wall_ns();
    if calls.is_empty() {
        crate::hc_profile::record_genotype_site(crate::hc_profile::GenotypeSiteSample {
            candidate_alleles: 0,
            genotype_states: 0,
            pl_vector_len: 0,
            samples: 1,
            wall_ns: region_wall.as_nanos() as u64,
            ad_wall_ns: ad_ns,
            event_rebuild_wall_ns: event_ns,
            allele_map_wall_ns: map_ns,
            marginalize_wall_ns: marg_ns,
            genotype_enum_wall_ns: enum_ns,
        });
        return;
    }
    let n = calls.len() as u64;
    let per_site_ns = region_wall.as_nanos() as u64 / n.max(1);
    let per_ad = ad_ns / n.max(1);
    let per_ev = event_ns / n.max(1);
    let per_map = map_ns / n.max(1);
    let per_marg = marg_ns / n.max(1);
    let per_enum = enum_ns / n.max(1);
    for call in calls {
        let alleles = 1u64 + u64::from(!call.event.alt_allele.is_empty());
        // Diploid PL length: typically 3 for biallelic; use actual GL vector.
        let pl_len = call.genotype.genotype_log10_likelihoods.len() as u64;
        let states = pl_len; // PL entries ≡ genotype states for diploid unphased
        crate::hc_profile::record_genotype_site(crate::hc_profile::GenotypeSiteSample {
            candidate_alleles: alleles.max(2),
            genotype_states: states,
            pl_vector_len: pl_len,
            samples: 1,
            wall_ns: per_site_ns,
            ad_wall_ns: per_ad,
            event_rebuild_wall_ns: per_ev,
            allele_map_wall_ns: per_map,
            marginalize_wall_ns: per_marg,
            genotype_enum_wall_ns: per_enum,
        });
    }
}

// NOTE: helper kept next to assign return sites — call from each Ok(...) path below.

/// Sort-key allele support using the region EventMap cache when possible.
///
/// Observable contract: same truth as [`haplotype_supports_allele_at_with_ref`]; avoid
/// rebuilding EventMaps O(alleles × haps) on the dense position walk.
/// When `hap_snp_bases` is populated (one entry per haplotype), SNP checks use it
/// instead of rewalking CIGAR.
fn hap_supports_allele_for_sort_cached(
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
    hap_snp_bases: &[Option<u8>],
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
        if let Some(base) = hap_snp_bases.get(hap_index).copied().flatten() {
            return base == alt_byte;
        }
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
    // Membership-only dedupe — HashSet matches BTreeSet cardinality/content.
    let mut seen: std::collections::HashSet<(u64, String, String)> = from_haps
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
}

/// GATK 4.4 `createAlleleMapper` haplotype pools for a merged SNP+indel VC.
///
/// Empty overlapping EventMap → REF. Matching event at `loc` → that ALT.
/// Spanning-elsewhere with `emitSpanningDels` → SPAN_DEL (`*`) when that allele is
/// in the merged list (6R.85). Unmatched at-loc events → no pool. Leftovers are
/// **not** dumped into REF (6R.84).
fn colocated_merge_allele_pools(
    haplotypes: &[Haplotype],
    loc: u64,
    long_ref: &str,
    alts: &[String],
    emit_spanning: bool,
    hap_events: Option<&crate::event_map::PerHaplotypeVariationEvents>,
    contig: &str,
    end_1based: u64,
    pad_start_1based: u64,
    ref_bytes: &[u8],
    max_mnp_distance: usize,
) -> Vec<Vec<HaplotypeIndex>> {
    let ref_idx = haplotypes.iter().position(|h| h.is_reference).unwrap_or(0);
    if let Some(cache) = hap_events {
        let (mut java_pools, _, _) = java_create_allele_mapper_pools(
            haplotypes.len(),
            cache,
            loc,
            long_ref,
            alts,
            emit_spanning,
        );
        if java_pools[0].is_empty() {
            java_pools[0].push(HaplotypeIndex::new(ref_idx));
        }
        return java_pools;
    }
    let mut allele_pools: Vec<Vec<HaplotypeIndex>> = Vec::with_capacity(1 + alts.len());
    let mut assigned = HashSet::new();
    for alt in alts {
        let ev = VariationEvent {
            contig: contig.to_string(),
            start_1based: GenomePosition::new_1based(loc),
            end_1based: GenomePosition::new_1based(end_1based),
            ref_allele: long_ref.to_string(),
            alt_allele: alt.clone(),
        };
        let mapping = create_allele_mapper_with_events(
            &ev,
            loc,
            haplotypes,
            pad_start_1based,
            ref_bytes,
            max_mnp_distance,
            emit_spanning,
            None,
        );
        let mut pool = Vec::new();
        for hi in mapping.alt_haplotype_indices {
            if assigned.insert(hi.get()) {
                pool.push(hi);
            }
        }
        allele_pools.push(pool);
    }
    let mut ref_pool: Vec<HaplotypeIndex> = Vec::new();
    for i in 0..haplotypes.len() {
        if !assigned.contains(&i) {
            ref_pool.push(HaplotypeIndex::new(i));
        }
    }
    if ref_pool.is_empty() {
        ref_pool.push(HaplotypeIndex::new(ref_idx));
    }
    let mut pools = Vec::with_capacity(1 + alts.len());
    pools.push(ref_pool);
    pools.extend(allele_pools);
    pools
}

/// Outcome of 6R.61 merged genotyping + 6R.62 unused-ALT subset.
///
/// `NotApplicable` falls through to the biallelic prefer-indel walk.
/// `MergedNoEmit` means the colocated merge applied but unused-ALT subset left only REF
/// — do **not** genotype T/C independently (that would recreate the 6R.60 lifecycle).
enum ColocatedMergeGenotype {
    NotApplicable,
    MergedNoEmit,
    Call(GenotypedSiteCall),
}

/// GATK `makeMergedVariantContext` then `calculateGLsForThisEvent` for colocated SNP+indel,
/// then `calculateOutputAlleleSubset` / `AlleleSubsettingUtils.subsetAlleles` (6R.62).
fn try_genotype_colocated_snp_indel_merge(
    events_at_loc: &[VariationEvent],
    loc: u64,
    likelihoods: &[RegionReadLikelihood],
    likelihood_reads: &[SharedBamRecord],
    haplotypes: &[Haplotype],
    ref_bytes: &[u8],
    pad_start_1based: u64,
    active_start_1based: u64,
    active_end_1based: u64,
    max_mnp_distance: usize,
    config: &HcGenotypingConfig,
    hap_events: Option<&crate::event_map::PerHaplotypeVariationEvents>,
) -> GatkResult<ColocatedMergeGenotype> {
    let Some((long_ref, alts)) = merged_alleles_for_genotyping(events_at_loc, loc) else {
        return Ok(ColocatedMergeGenotype::NotApplicable);
    };
    if !is_colocated_snp_indel_merged_site(&long_ref, &alts) {
        return Ok(ColocatedMergeGenotype::NotApplicable);
    }
    let contig = events_at_loc
        .iter()
        .find(|e| e.start_1based.get() == loc)
        .map(|e| e.contig.as_str())
        .unwrap_or("");
    let end_1based = loc.saturating_add(long_ref.len().saturating_sub(1) as u64);
    let primary = VariationEvent {
        contig: contig.to_string(),
        start_1based: GenomePosition::new_1based(loc),
        end_1based: GenomePosition::new_1based(end_1based),
        ref_allele: long_ref.clone(),
        alt_allele: alts[0].clone(),
    };
    let emit_spanning = !config.disable_spanning_event_genotyping;
    let ref_idx = haplotypes.iter().position(|h| h.is_reference).unwrap_or(0);
    let pools = colocated_merge_allele_pools(
        haplotypes,
        loc,
        &long_ref,
        &alts,
        emit_spanning,
        hap_events,
        contig,
        end_1based,
        pad_start_1based,
        ref_bytes,
        max_mnp_distance,
    );
    if pools.iter().skip(1).all(|p| p.is_empty()) {
        return Ok(ColocatedMergeGenotype::NotApplicable);
    }
    let subset = likelihood_subset_for_event(
        likelihoods,
        likelihood_reads,
        &primary,
        config,
        active_start_1based,
        active_end_1based,
    );
    if subset.is_empty() {
        return Ok(ColocatedMergeGenotype::NotApplicable);
    }
    let hap_rows = region_likelihoods_to_rows(subset.as_ref(), haplotypes.len());
    let mut marg: Vec<ReadLikelihoodRow> = hap_rows
        .iter()
        .map(|row| {
            let lls: Vec<f64> = pools.iter().map(|pool| pool_max_log10(pool, row)).collect();
            ReadLikelihoodRow {
                read_index: row.read_index,
                read_id: row.read_id.clone(),
                haplotype_log10_likelihoods: lls,
            }
        })
        .collect();
    apply_java_marginal_normalize_n(&mut marg);
    let n_alleles = 1 + alts.len();
    let gls = diploid_genotype_log10_likelihoods_from_allele_rows(&marg, n_alleles);
    let depths = informative_ad_n_alleles(&marg, n_alleles);
    let merged_format = emit_genotype_format_fields(&gls, &depths)?;
    let assigned_gt = diploid_genotype_alleles_from_pl_index(
        n_alleles,
        best_pl_index(&merged_format.pl),
    );
    // GATK `GenotypingEngine.calculateGenotypes`: QUAL from AFCalculator on the
    // pre-subset merged VC (including SPAN_DEL). `builder.log10PError` is copied
    // through unused-ALT subset and reverse-trim; it is not recomputed from emitted PL.
    let qual_log10_p_error = if alts.iter().any(|a| a == SPAN_DEL_ALLELE) {
        let merged_alleles: Vec<&str> = std::iter::once(long_ref.as_str())
            .chain(alts.iter().map(String::as_str))
            .collect();
        let gl_rt = genotype_log10_likelihoods_after_java_genotype_pl_roundtrip(&gls);
        Some(diploid_af_log10_prob_only_ref_allele_exists(
            &gl_rt,
            &merged_alleles,
            &AfCalculatorConfig::default(),
        )?)
    } else {
        None
    };
    // GATK `calculateOutputAlleleSubset` + `AlleleSubsettingUtils.subsetAlleles` after
    // merged GLs. Reverse-trim only after this subset, and only when allele count
    // changed (`makeAnnotatedCall`). Do not recalculate GLs from reads.
    let unused = crate::allele_subsetting_pl::subset_unused_alts_after_merged_genotyping(
        &alts,
        &assigned_gt,
        &gls,
        &depths,
    )?;
    if unused.alt_alleles.is_empty() {
        return Ok(ColocatedMergeGenotype::MergedNoEmit);
    }
    // GATK `AlleleSubsettingUtils.subsetAlleles` does not write AD (`genotype.hasAD()` is
    // false after `calculateGLsForThisEvent`). The first AD write is
    // `DepthPerAlleleBySample.annotateWithLikelihoods`: identity remarg of remaining
    // `vc.getAlleles()`, then `bestAllelesBreakingTies` + `isInformative` (`> 0.2`).
    // Slicing 4-way informative counts is not that operation.
    let keep_idx = remaining_keep_indices(&alts, &unused.alt_alleles);
    let remarg_ad = remarg_informative_ad(&marg, &keep_idx);
    if let Ok(mut snap) = colocated_merge_numerics_snapshot(
        loc,
        long_ref.clone(),
        &alts,
        &pools,
        &marg,
        gls.clone(),
        depths.clone(),
        merged_format.pl_as_i32(),
        assigned_gt.clone(),
        &unused,
        likelihood_reads,
    ) {
        fill_read_set_audit(
            &mut snap,
            likelihoods,
            likelihood_reads,
            haplotypes,
            &pools,
            &alts,
            &primary,
            config,
            hap_events,
        );
        COLOCATED_MERGE_NUMERICS.with(|slot| slot.borrow_mut().push(snap));
    }
    let format = emit_genotype_format_fields(&unused.log10_gls, &remarg_ad)?;
    let aggregation = aggregate_haplotype_log10_likelihoods(&hap_rows)?;
    let best = best_haplotype_index(&aggregation)
        .unwrap_or(crate::bio_ids::HaplotypeIndex::new(0))
        .get();
    let kept_first = unused.alt_alleles[0].clone();
    let kept_pool_idx = alts
        .iter()
        .position(|a| a == &kept_first)
        .map(|i| i + 1)
        .unwrap_or(1);
    // GATK `reverseTrimAlleles` after unused-ALT subset (not before). Java invokes
    // it iff allele count shrank vs the merged VC. Does not recalculate GLs/AD/GT.
    // Start is unchanged; END follows new REF length.
    let (trim_ref, trim_alts) = if unused.alt_alleles.len() < alts.len() {
        crate::reverse_trim::reverse_trim_alleles(&long_ref, &unused.alt_alleles)
    } else {
        (long_ref.clone(), unused.alt_alleles.clone())
    };
    let mut event = primary;
    event.ref_allele = trim_ref;
    event.alt_allele = trim_alts
        .first()
        .cloned()
        .unwrap_or(kept_first);
    let extra: Vec<String> = trim_alts.into_iter().skip(1).collect();
    event.end_1based = GenomePosition::new_1based(VariationEvent::vcf_end_1based(
        event.start_1based.get(),
        &event.ref_allele,
    ));
    Ok(ColocatedMergeGenotype::Call(GenotypedSiteCall {
        event,
        genotype: RegionGenotypeResult {
            aggregation,
            best_haplotype_index: best,
            ref_haplotype_index: ref_idx,
            alt_haplotype_index: pools
                .get(kept_pool_idx)
                .and_then(|p| p.first())
                .map(|h| h.get())
                .unwrap_or(0),
            genotype_log10_likelihoods: unused.log10_gls,
            format,
        },
        extra_alt_alleles: extra,
        post_merge_unused_alt_subset: unused.alt_alleles.len() < alts.len(),
        qual_log10_p_error,
    }))
}

/// 6R.64 forensic snapshot of colocated SNP+indel merge numerics.
///
/// Does not change production genotyping. Captured from the live merge path and
/// also via [`audit_colocated_snp_indel_merge_numerics`]. Includes the Java
/// annotation-time AD remarginalize (`DepthPerAlleleBySample` on remaining alleles).
#[derive(Debug, Clone)]
pub struct ColocatedMergeNumerics {
    pub loc: u64,
    pub long_ref: String,
    pub alts: Vec<String>,
    pub n_reads: usize,
    pub pool_sizes: Vec<usize>,
    pub merged_pl: Vec<i32>,
    pub merged_gls: Vec<f64>,
    pub merged_ad: Vec<i32>,
    pub assigned_gt: Vec<i32>,
    pub subset_pl: Vec<i32>,
    pub subset_ad_permuted: Vec<i32>,
    pub subset_ad_remarginalized: Vec<i32>,
    pub n_uninformative_3way: usize,
    /// Unique PairHMM read indices in the region matrix (6R.65 A).
    pub n_pairhmm_reads: usize,
    /// Alignment-overlap retainEvidence before QNAME collapse (Java default HC).
    pub n_overlap_before_qname_dedupe: usize,
    pub n_overlap_unique_qnames: usize,
    pub n_qnames_with_multiple_overlapping_reads: usize,
    /// 6-PL / unused-ALT 3-PL if QNAME collapse is skipped (diagnostic counterfactual).
    pub merged_pl_no_qname_dedupe: Vec<i32>,
    pub subset_pl_no_qname_dedupe: Vec<i32>,
    pub merged_ad_no_qname_dedupe: Vec<i32>,
    /// DepthPerAlleleBySample remarg on overlap-retained reads (no QNAME collapse).
    pub subset_ad_remarginalized_no_qname: Vec<i32>,
    /// Unused-ALT permute of 4-way AD on overlap-retained reads (no QNAME collapse).
    pub subset_ad_permuted_no_qname: Vec<i32>,
    /// Keep indices into the merged allele matrix for remaining call alleles (REF + kept ALTs).
    pub remaining_keep_indices: Vec<usize>,
    /// Production merged-allele matrix rows used for AD (6R.101). Match by (QNAME, flags).
    pub ad_row_read_index: Vec<usize>,
    pub ad_row_qname: Vec<String>,
    pub ad_row_flags: Vec<u16>,
    pub ad_row_lls: Vec<Vec<f64>>,
    pub n_haps: usize,
    pub n_haps_with_multiple_events_at_loc: usize,
    /// Distinct EventMap alleles at `loc` (`ref>alt` → hap count).
    pub hap_event_signatures_at_loc: Vec<(String, usize)>,
    pub n_cache_events: usize,
    pub cache_starts_near_loc: Vec<(u64, usize)>,
    pub nearest_event_start_below: Option<u64>,
    pub nearest_event_start_above: Option<u64>,
    pub java_style_pool_sizes: Vec<usize>,
    pub n_haps_unassigned_java: usize,
    pub n_haps_in_multiple_java_alts: usize,
    pub n_haps_rust_ref_but_java_unassigned: usize,
    pub merged_pl_java_style_pools: Vec<i32>,
    pub subset_pl_java_style_pools: Vec<i32>,
    pub n_allele_floor_clips: usize,
    pub merged_pl_no_allele_floor: Vec<i32>,
    pub subset_pl_no_allele_floor: Vec<i32>,
}

thread_local! {
    static COLOCATED_MERGE_NUMERICS: RefCell<Vec<ColocatedMergeNumerics>> =
        const { RefCell::new(Vec::new()) };
}

/// Drain live-merge numerics recorded by the last `assign_genotype_likelihoods_for_region`.
pub fn take_colocated_merge_numerics() -> Vec<ColocatedMergeNumerics> {
    COLOCATED_MERGE_NUMERICS.with(|slot| std::mem::take(&mut *slot.borrow_mut()))
}

fn remaining_keep_indices(alts: &[String], kept_alts: &[String]) -> Vec<usize> {
    let mut keep = vec![0usize];
    for kept in kept_alts {
        if let Some(i) = alts.iter().position(|a| a == kept) {
            keep.push(i + 1);
        }
    }
    keep
}

/// GATK `DepthPerAlleleBySample` identity remarg: vote over remaining allele columns only.
fn remarg_informative_ad(marg: &[ReadLikelihoodRow], keep_idx: &[usize]) -> Vec<i32> {
    let remarg: Vec<ReadLikelihoodRow> = marg
        .iter()
        .map(|row| ReadLikelihoodRow {
            read_index: row.read_index,
            read_id: row.read_id.clone(),
            haplotype_log10_likelihoods: keep_idx
                .iter()
                .filter_map(|&i| row.haplotype_log10_likelihoods.get(i).copied())
                .collect(),
        })
        .collect();
    informative_ad_n_alleles(&remarg, keep_idx.len())
}

fn colocated_merge_numerics_snapshot(
    loc: u64,
    long_ref: String,
    alts: &[String],
    pools: &[Vec<HaplotypeIndex>],
    marg: &[ReadLikelihoodRow],
    gls: Vec<f64>,
    depths: Vec<i32>,
    merged_pl: Vec<i32>,
    assigned_gt: Vec<i32>,
    unused: &crate::allele_subsetting_pl::UnusedAltSubsetResult,
    likelihood_reads: &[SharedBamRecord],
) -> GatkResult<ColocatedMergeNumerics> {
    let subset_format = emit_genotype_format_fields(&unused.log10_gls, &unused.ad)?;
    let keep_idx = remaining_keep_indices(alts, &unused.alt_alleles);
    let remarg_ad = remarg_informative_ad(marg, &keep_idx);
    let n_informative: i32 = depths.iter().sum();
    let mut ad_row_read_index = Vec::with_capacity(marg.len());
    let mut ad_row_qname = Vec::with_capacity(marg.len());
    let mut ad_row_flags = Vec::with_capacity(marg.len());
    let mut ad_row_lls = Vec::with_capacity(marg.len());
    for row in marg {
        ad_row_read_index.push(row.read_index);
        let rec = likelihood_reads.get(row.read_index);
        ad_row_qname.push(
            rec.map(|r| String::from_utf8_lossy(r.qname()).into_owned())
                .unwrap_or_default(),
        );
        ad_row_flags.push(rec.map(|r| r.flags()).unwrap_or(0));
        ad_row_lls.push(row.haplotype_log10_likelihoods.clone());
    }
    Ok(ColocatedMergeNumerics {
        loc,
        long_ref,
        alts: alts.to_vec(),
        n_reads: marg.len(),
        pool_sizes: pools.iter().map(|p| p.len()).collect(),
        merged_pl,
        merged_gls: gls,
        merged_ad: depths,
        assigned_gt,
        subset_pl: subset_format.pl_as_i32(),
        subset_ad_permuted: unused.ad.clone(),
        subset_ad_remarginalized: remarg_ad,
        n_uninformative_3way: marg.len().saturating_sub(n_informative.max(0) as usize),
        n_pairhmm_reads: 0,
        n_overlap_before_qname_dedupe: 0,
        n_overlap_unique_qnames: 0,
        n_qnames_with_multiple_overlapping_reads: 0,
        merged_pl_no_qname_dedupe: Vec::new(),
        subset_pl_no_qname_dedupe: Vec::new(),
        merged_ad_no_qname_dedupe: Vec::new(),
        subset_ad_remarginalized_no_qname: Vec::new(),
        subset_ad_permuted_no_qname: Vec::new(),
        remaining_keep_indices: keep_idx,
        ad_row_read_index,
        ad_row_qname,
        ad_row_flags,
        ad_row_lls,
        n_haps: 0,
        n_haps_with_multiple_events_at_loc: 0,
        hap_event_signatures_at_loc: Vec::new(),
        n_cache_events: 0,
        cache_starts_near_loc: Vec::new(),
        nearest_event_start_below: None,
        nearest_event_start_above: None,
        java_style_pool_sizes: Vec::new(),
        n_haps_unassigned_java: 0,
        n_haps_in_multiple_java_alts: 0,
        n_haps_rust_ref_but_java_unassigned: 0,
        merged_pl_java_style_pools: Vec::new(),
        subset_pl_java_style_pools: Vec::new(),
        n_allele_floor_clips: 0,
        merged_pl_no_allele_floor: Vec::new(),
        subset_pl_no_allele_floor: Vec::new(),
    })
}

fn unique_likelihood_read_indices(ll: &[RegionReadLikelihood]) -> HashSet<usize> {
    ll.iter().map(|e| e.read_index.get()).collect()
}

fn qname_overlap_counts(
    overlapping: &[RegionReadLikelihood],
    reads: &[SharedBamRecord],
) -> (usize, usize) {
    let idxs = unique_likelihood_read_indices(overlapping);
    let mut per_qname: HashMap<Vec<u8>, usize> = HashMap::new();
    for idx in idxs {
        let Some(rec) = reads.get(idx) else {
            continue;
        };
        *per_qname.entry(rec.qname().to_owned()).or_insert(0) += 1;
    }
    let n_qnames = per_qname.len();
    let n_multi = per_qname.values().filter(|&&c| c > 1).count();
    (n_qnames, n_multi)
}

fn colocated_pls_from_likelihood_subset(
    subset: &[RegionReadLikelihood],
    haplotypes: &[Haplotype],
    pools: &[Vec<HaplotypeIndex>],
    alts: &[String],
) -> GatkResult<(Vec<i32>, Vec<i32>, Vec<i32>, Vec<i32>, Vec<i32>)> {
    colocated_pls_from_likelihood_subset_floor(subset, haplotypes, pools, alts, true)
}

fn colocated_pls_from_likelihood_subset_floor(
    subset: &[RegionReadLikelihood],
    haplotypes: &[Haplotype],
    pools: &[Vec<HaplotypeIndex>],
    alts: &[String],
    allele_floor: bool,
) -> GatkResult<(Vec<i32>, Vec<i32>, Vec<i32>, Vec<i32>, Vec<i32>)> {
    let hap_rows = region_likelihoods_to_rows(subset, haplotypes.len());
    let mut marg: Vec<ReadLikelihoodRow> = hap_rows
        .iter()
        .map(|row| {
            let lls: Vec<f64> = pools.iter().map(|pool| pool_max_log10(pool, row)).collect();
            ReadLikelihoodRow {
                read_index: row.read_index,
                read_id: row.read_id.clone(),
                haplotype_log10_likelihoods: lls,
            }
        })
        .collect();
    if allele_floor {
        apply_java_marginal_normalize_n(&mut marg);
    }
    let n_alleles = 1 + alts.len();
    let gls = diploid_genotype_log10_likelihoods_from_allele_rows(&marg, n_alleles);
    let depths = informative_ad_n_alleles(&marg, n_alleles);
    let merged_format = emit_genotype_format_fields(&gls, &depths)?;
    let assigned_gt =
        diploid_genotype_alleles_from_pl_index(n_alleles, best_pl_index(&merged_format.pl));
    let unused = crate::allele_subsetting_pl::subset_unused_alts_after_merged_genotyping(
        alts, &assigned_gt, &gls, &depths,
    )?;
    let subset_format = emit_genotype_format_fields(&unused.log10_gls, &unused.ad)?;
    let keep_idx = remaining_keep_indices(alts, &unused.alt_alleles);
    let remarg_ad = remarg_informative_ad(&marg, &keep_idx);
    Ok((
        merged_format.pl_as_i32(),
        subset_format.pl_as_i32(),
        depths,
        remarg_ad,
        unused.ad.clone(),
    ))
}

/// GATK 4.4 `createAlleleMapper` EventMap walk: unmatched events leave the hap in **no** pool.
/// Dual alt membership is allowed. Leftover unmatched haps are **not** dumped into REF.
fn java_create_allele_mapper_pools(
    n_haps: usize,
    hap_events: &crate::event_map::PerHaplotypeVariationEvents,
    loc: u64,
    long_ref: &str,
    alts: &[String],
    emit_spanning: bool,
) -> (Vec<Vec<HaplotypeIndex>>, usize, usize) {
    let loc_pos = GenomePosition::new_1based(loc);
    let mut pools: Vec<Vec<HaplotypeIndex>> = vec![Vec::new(); 1 + alts.len()];
    let mut n_unassigned = 0usize;
    let mut n_multi = 0usize;
    for i in 0..n_haps {
        let spanning = overlapping_events(hap_events.events_for(i), loc);
        if spanning.is_empty() {
            pools[0].push(HaplotypeIndex::new(i));
            continue;
        }
        let mut in_alt = vec![false; alts.len()];
        let mut in_ref = false;
        for ev in &spanning {
            if ev.start_1based == loc_pos {
                if ev.ref_allele.len() == long_ref.len() {
                    if let Some(ai) = alts.iter().position(|a| a == &ev.alt_allele) {
                        in_alt[ai] = true;
                    }
                } else if ev.ref_allele.len() < long_ref.len() {
                    if let Some(remapped) =
                        remap_alt_onto_longer_ref(&ev.ref_allele, &ev.alt_allele, long_ref)
                    {
                        if let Some(ai) = alts.iter().position(|a| a == &remapped) {
                            in_alt[ai] = true;
                        }
                    }
                }
            } else if emit_spanning {
                if let Some(ai) = alts.iter().position(|a| a == "*") {
                    in_alt[ai] = true;
                }
                break;
            } else {
                in_ref = true;
                break;
            }
        }
        let n_alts_hit = in_alt.iter().filter(|&&b| b).count();
        if n_alts_hit > 1 {
            n_multi += 1;
        }
        let any = in_ref || n_alts_hit > 0;
        if !any {
            n_unassigned += 1;
            continue;
        }
        if in_ref {
            pools[0].push(HaplotypeIndex::new(i));
        }
        for (ai, hit) in in_alt.iter().enumerate() {
            if *hit {
                pools[ai + 1].push(HaplotypeIndex::new(i));
            }
        }
    }
    (pools, n_unassigned, n_multi)
}

fn count_allele_floor_clips(rows: &[ReadLikelihoodRow]) -> usize {
    let mut n = 0usize;
    for row in rows {
        let best = row
            .haplotype_log10_likelihoods
            .iter()
            .copied()
            .filter(|v| v.is_finite())
            .fold(f64::NEG_INFINITY, f64::max);
        if !best.is_finite() {
            continue;
        }
        let floor = best + LOG10_GLOBAL_READ_MISMATCHING_RATE;
        n += row
            .haplotype_log10_likelihoods
            .iter()
            .filter(|v| v.is_finite() && **v < floor)
            .count();
    }
    n
}

fn fill_read_set_audit(
    snap: &mut ColocatedMergeNumerics,
    likelihoods: &[RegionReadLikelihood],
    likelihood_reads: &[SharedBamRecord],
    haplotypes: &[Haplotype],
    pools: &[Vec<HaplotypeIndex>],
    alts: &[String],
    primary: &VariationEvent,
    config: &HcGenotypingConfig,
    hap_events: Option<&crate::event_map::PerHaplotypeVariationEvents>,
) {
    snap.n_pairhmm_reads = unique_likelihood_read_indices(likelihoods).len();
    snap.n_haps = haplotypes.len();
    let loc_u64 = primary.start_1based.get();
    if let Some(cache) = hap_events {
        let loc = primary.start_1based;
        snap.n_haps_with_multiple_events_at_loc = (0..haplotypes.len())
            .filter(|&i| {
                cache
                    .events_for(i)
                    .iter()
                    .filter(|e| e.start_1based == loc)
                    .count()
                    > 1
            })
            .count();
        let mut sigs: HashMap<String, usize> = HashMap::new();
        let mut n_cache_events = 0usize;
        let mut near: HashMap<u64, usize> = HashMap::new();
        for i in 0..haplotypes.len() {
            let evs = cache.events_for(i);
            n_cache_events += evs.len();
            let at_loc: Vec<String> = evs
                .iter()
                .filter(|e| e.start_1based == loc)
                .map(|e| format!("{}/{}", e.ref_allele, e.alt_allele))
                .collect();
            let key = if at_loc.is_empty() {
                "none".to_string()
            } else {
                at_loc.join("+")
            };
            *sigs.entry(key).or_insert(0) += 1;
            for e in evs {
                let s = e.start_1based.get();
                if s.abs_diff(loc_u64) <= 25 {
                    *near.entry(s).or_insert(0) += 1;
                }
            }
        }
        snap.n_cache_events = n_cache_events;
        let mut near_rows: Vec<(u64, usize)> = near.into_iter().collect();
        near_rows.sort_by_key(|(s, _)| *s);
        snap.nearest_event_start_below = near_rows
            .iter()
            .rev()
            .find(|(s, _)| *s < loc_u64)
            .map(|(s, _)| *s);
        snap.nearest_event_start_above = near_rows
            .iter()
            .find(|(s, _)| *s > loc_u64)
            .map(|(s, _)| *s);
        snap.cache_starts_near_loc = near_rows;
        let mut rows: Vec<(String, usize)> = sigs.into_iter().collect();
        rows.sort_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(&b.0)));
        snap.hap_event_signatures_at_loc = rows;
        let emit_spanning = !config.disable_spanning_event_genotyping;
        let (java_pools, n_unassigned, n_multi) = java_create_allele_mapper_pools(
            haplotypes.len(),
            cache,
            loc_u64,
            &snap.long_ref,
            alts,
            emit_spanning,
        );
        snap.java_style_pool_sizes = java_pools.iter().map(|p| p.len()).collect();
        snap.n_haps_unassigned_java = n_unassigned;
        snap.n_haps_in_multiple_java_alts = n_multi;
        let rust_ref: HashSet<usize> = pools
            .first()
            .map(|p| p.iter().map(|h| h.get()).collect())
            .unwrap_or_default();
        let java_any: HashSet<usize> = java_pools
            .iter()
            .flat_map(|p| p.iter().map(|h| h.get()))
            .collect();
        snap.n_haps_rust_ref_but_java_unassigned = rust_ref
            .iter()
            .filter(|i| !java_any.contains(i))
            .count();
    }
    let var_end = primary.end_1based.get().max(
        primary
            .start_1based
            .get()
            .saturating_add(primary.ref_allele.len().saturating_sub(1) as u64),
    );
    let overlap = filter_likelihoods_for_variant(
        likelihoods,
        likelihood_reads,
        primary,
        primary.start_1based.get(),
        var_end,
        config.informative_read_overlap_margin,
        config,
    );
    snap.n_overlap_before_qname_dedupe = unique_likelihood_read_indices(&overlap).len();
    let (n_q, n_multi) = qname_overlap_counts(&overlap, likelihood_reads);
    snap.n_overlap_unique_qnames = n_q;
    snap.n_qnames_with_multiple_overlapping_reads = n_multi;
    if let Ok((mpl, spl, ad, remarg, perm)) =
        colocated_pls_from_likelihood_subset(&overlap, haplotypes, pools, alts)
    {
        snap.merged_pl_no_qname_dedupe = mpl;
        snap.subset_pl_no_qname_dedupe = spl;
        snap.merged_ad_no_qname_dedupe = ad;
        snap.subset_ad_remarginalized_no_qname = remarg;
        snap.subset_ad_permuted_no_qname = perm;
    }
    let hap_rows = region_likelihoods_to_rows(&overlap, haplotypes.len());
    let marg_raw: Vec<ReadLikelihoodRow> = hap_rows
        .iter()
        .map(|row| {
            let lls: Vec<f64> = pools.iter().map(|pool| pool_max_log10(pool, row)).collect();
            ReadLikelihoodRow {
                read_index: row.read_index,
                read_id: row.read_id.clone(),
                haplotype_log10_likelihoods: lls,
            }
        })
        .collect();
    snap.n_allele_floor_clips = count_allele_floor_clips(&marg_raw);
    if let Ok((mpl, spl, _, _, _)) =
        colocated_pls_from_likelihood_subset_floor(&overlap, haplotypes, pools, alts, false)
    {
        snap.merged_pl_no_allele_floor = mpl;
        snap.subset_pl_no_allele_floor = spl;
    }
    if let Some(cache) = hap_events {
        let emit_spanning = !config.disable_spanning_event_genotyping;
        let (java_pools, _, _) = java_create_allele_mapper_pools(
            haplotypes.len(),
            cache,
            loc_u64,
            &snap.long_ref,
            alts,
            emit_spanning,
        );
        if let Ok((mpl, spl, _, _, _)) =
            colocated_pls_from_likelihood_subset(&overlap, haplotypes, &java_pools, alts)
        {
            snap.merged_pl_java_style_pools = mpl;
            snap.subset_pl_java_style_pools = spl;
        }
    }
}

/// Replay [`try_genotype_colocated_snp_indel_merge`] numerics for diagnosis (6R.64).
pub fn audit_colocated_snp_indel_merge_numerics(
    events_at_loc: &[VariationEvent],
    loc: u64,
    likelihoods: &[RegionReadLikelihood],
    likelihood_reads: &[SharedBamRecord],
    haplotypes: &[Haplotype],
    ref_bytes: &[u8],
    pad_start_1based: u64,
    active_start_1based: u64,
    active_end_1based: u64,
    max_mnp_distance: usize,
    config: &HcGenotypingConfig,
    hap_events: Option<&crate::event_map::PerHaplotypeVariationEvents>,
) -> GatkResult<Option<ColocatedMergeNumerics>> {
    let Some((long_ref, alts)) = merged_alleles_for_genotyping(events_at_loc, loc) else {
        return Ok(None);
    };
    if !is_colocated_snp_indel_merged_site(&long_ref, &alts) {
        return Ok(None);
    }
    let contig = events_at_loc
        .iter()
        .find(|e| e.start_1based.get() == loc)
        .map(|e| e.contig.as_str())
        .unwrap_or("");
    let end_1based = loc.saturating_add(long_ref.len().saturating_sub(1) as u64);
    let primary = VariationEvent {
        contig: contig.to_string(),
        start_1based: GenomePosition::new_1based(loc),
        end_1based: GenomePosition::new_1based(end_1based),
        ref_allele: long_ref.clone(),
        alt_allele: alts[0].clone(),
    };
    let emit_spanning = !config.disable_spanning_event_genotyping;
    let pools = colocated_merge_allele_pools(
        haplotypes,
        loc,
        &long_ref,
        &alts,
        emit_spanning,
        hap_events,
        contig,
        end_1based,
        pad_start_1based,
        ref_bytes,
        max_mnp_distance,
    );
    if pools.iter().skip(1).all(|p| p.is_empty()) {
        return Ok(None);
    }
    let subset = likelihood_subset_for_event(
        likelihoods,
        likelihood_reads,
        &primary,
        config,
        active_start_1based,
        active_end_1based,
    );
    if subset.is_empty() {
        return Ok(None);
    }
    let hap_rows = region_likelihoods_to_rows(subset.as_ref(), haplotypes.len());
    let mut marg: Vec<ReadLikelihoodRow> = hap_rows
        .iter()
        .map(|row| {
            let lls: Vec<f64> = pools.iter().map(|pool| pool_max_log10(pool, row)).collect();
            ReadLikelihoodRow {
                read_index: row.read_index,
                read_id: row.read_id.clone(),
                haplotype_log10_likelihoods: lls,
            }
        })
        .collect();
    apply_java_marginal_normalize_n(&mut marg);
    let n_alleles = 1 + alts.len();
    let gls = diploid_genotype_log10_likelihoods_from_allele_rows(&marg, n_alleles);
    let depths = informative_ad_n_alleles(&marg, n_alleles);
    let merged_format = emit_genotype_format_fields(&gls, &depths)?;
    let assigned_gt = diploid_genotype_alleles_from_pl_index(
        n_alleles,
        best_pl_index(&merged_format.pl),
    );
    let unused = crate::allele_subsetting_pl::subset_unused_alts_after_merged_genotyping(
        &alts,
        &assigned_gt,
        &gls,
        &depths,
    )?;
    if unused.alt_alleles.is_empty() {
        return Ok(None);
    }
    let mut snap = colocated_merge_numerics_snapshot(
        loc,
        long_ref,
        &alts,
        &pools,
        &marg,
        gls,
        depths,
        merged_format.pl_as_i32(),
        assigned_gt,
        &unused,
        likelihood_reads,
    )?;
    fill_read_set_audit(
        &mut snap,
        likelihoods,
        likelihood_reads,
        haplotypes,
        &pools,
        &alts,
        &primary,
        config,
        hap_events,
    );
    Ok(Some(snap))
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

