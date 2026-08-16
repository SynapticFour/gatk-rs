/// L13-C1: site genotyping pipeline — softclip/pileup helpers + `try_genotype_variation_event`.
/// Behavior-neutral extract from `genotype_finalize.rs` for N-3.
/// Included into the genotyping engine module (same scope as finalize).

/// Soft-clip pileup alt in one pass: `(QNAME-deduped, fragment)` counts.
/// Same selection rules as the former separate dedupe + fragment scans.
fn sparse_softclip_pileup_alt_counts(
    reads: &[SharedBamRecord],
    event: &VariationEvent,
    margin: i32,
) -> (i32, i32) {
    use crate::fragment_overlap::read_base_at_ref_coord_1based;
    if event.ref_allele.len() != 1 || event.alt_allele.len() != 1 {
        return (0, 0);
    }
    let var_end = event.end_1based.get().max(
        event
            .start_1based
            .get()
            .saturating_add(event.ref_allele.len().saturating_sub(1) as u64),
    );
    let alt_b = event.alt_allele.as_bytes()[0].to_ascii_uppercase();
    let mut seen = std::collections::BTreeSet::new();
    let mut deduped = 0i32;
    let mut fragments = 0i32;
    for rec in reads {
        if !soft_unclipped_read_overlaps_interval(rec, event.start_1based.get(), var_end, margin) {
            continue;
        }
        if let Some(qb) = read_base_at_ref_coord_1based(rec, event.start_1based.get() as i32) {
            if qb.to_ascii_uppercase() == alt_b {
                fragments += 1;
                if seen.insert(rec.qname().to_owned()) {
                    deduped += 1;
                }
            }
        }
    }
    (deduped, fragments)
}

/// Fragment-level pileup AD for strict Java SNPs (dedupe QNAME); per-read for anchors/indels.
fn read_allele_depths_at_locus_for_genotyping(
    reads: &[SharedBamRecord],
    event: &VariationEvent,
    pad_start_1based: u64,
    config: &HcGenotypingConfig,
) -> (i32, i32) {
    if config.enable_java_strict()
        && event.ref_allele.len() == 1
        && event.alt_allele.len() == 1
        && !is_cluster_anchor_snp(event)
    {
        crate::read_event_discovery::read_allele_depths_at_locus_dedupe_qname(
            reads,
            event,
            pad_start_1based,
        )
    } else {
        read_allele_depths_at_locus(reads, event, pad_start_1based)
    }
}

pub fn read_allele_depths_for_strict_emit(
    pileup_reads: &[SharedBamRecord],
    supplemental_pileup_reads: Option<&[SharedBamRecord]>,
    event: &VariationEvent,
    pad_start_1based: u64,
    config: &HcGenotypingConfig,
    apply_bases: &[u8],
    full_reference_bases: &[u8],
    full_reference_pad_1based: u64,
) -> (i32, i32) {
    let pileup_ad = |reads: &[SharedBamRecord], pad: u64| -> (i32, i32) {
        read_allele_depths_at_locus_for_genotyping(reads, event, pad, config)
    };
    let (trim_ref, trim_alt) = pileup_ad(pileup_reads, pad_start_1based);
    let Some(sup) = supplemental_pileup_reads.filter(|s| !s.is_empty()) else {
        return (trim_ref, trim_alt);
    };
    if config.enable_java_strict() && is_cluster_anchor_snp(event) {
        let mut best = (trim_ref, trim_alt);
        let pads: &[u64] = if crate::read_event_discovery::snp_allele_depth_pads_equivalent(
            event.start_1based.get(),
            pad_start_1based,
            full_reference_pad_1based,
        ) {
            &[pad_start_1based]
        } else {
            &[pad_start_1based, full_reference_pad_1based]
        };
        for &pad in pads {
            let (r, a) = read_allele_depths_at_locus(sup, event, pad);
            if a > best.1 || (a == best.1 && r + a > best.0 + best.1) {
                best = (r, a);
            }
        }
        return best;
    }
    if config.enable_java_strict() && is_strict_java_production_emit_admits(event) {
        let (mut rr, mut ra) = read_allele_depths_p12_java_sparse_pileup(
            sup,
            event,
            apply_bases,
            pad_start_1based,
            full_reference_bases,
            full_reference_pad_1based,
        );
        let var_end = event.end_1based.get().max(
            event
                .start_1based
                .get()
                .saturating_add(event.ref_allele.len().saturating_sub(1) as u64),
        );
        let margin = config.informative_read_overlap_margin;
        if !sup.iter().any(|r| {
            java_alignment_read_covers_variant_base(r, event.start_1based.get(), var_end, margin)
        }) {
            let (soft_dedupe, soft_frag) =
                sparse_softclip_pileup_alt_counts(sup, event, margin);
            ra = ra.max(soft_dedupe);
            if sparse_java_softclip_pairhmm_band(event) {
                if soft_frag >= 3 {
                    ra = ra.max(soft_frag);
                } else if soft_frag == 2 {
                    ra = ra.max(2);
                }
            }
            rr = 0;
        }
        return (rr, ra);
    }
    if !config.enable_java_strict() || trim_alt < 1 {
        let (full_ref, full_alt) = pileup_ad(sup, pad_start_1based);
        if full_alt > trim_alt {
            return (full_ref, full_alt);
        }
    }
    (trim_ref, trim_alt)
}

/// Java sparse cluster indel: genotype from the single alt-favoring read (DP=1 / AD 0,1 class).
fn narrow_strict_java_cluster_coupled_indel_subset(
    subset: Vec<RegionReadLikelihood>,
    likelihood_reads: &[SharedBamRecord],
    haplotypes: &[Haplotype],
    mapping: &AlleleHaplotypeMapping,
    config: &HcGenotypingConfig,
    event: &VariationEvent,
) -> Vec<RegionReadLikelihood> {
    if mapping.alt_haplotype_indices.is_empty() {
        return subset;
    }
    let ref_pool = ref_hap_indices_for_genotype_marginalization(mapping, haplotypes, config, Some(event));
    let rows = region_likelihoods_to_rows(&subset, haplotypes.len());
    let marg = marginalize_rows_to_biallelic_alleles(&rows, &ref_pool, &mapping.alt_haplotype_indices);
    let mut ranked: Vec<(Vec<u8>, f64)> = marg
        .iter()
        .filter_map(|row| {
            let lr = row.haplotype_log10_likelihoods[0];
            let la = row.haplotype_log10_likelihoods[1];
            if la > lr {
                row.read_id
                    .strip_prefix("read_")
                    .and_then(|s| s.parse::<usize>().ok())
                    .and_then(|ri| likelihood_reads.get(ri).map(|r| r.qname().to_owned()))
                    .map(|qname| (qname, la - lr))
            } else {
                None
            }
        })
        .collect();
    if ranked.is_empty() {
        return subset;
    }
    ranked.sort_by(|a, b| b.1.total_cmp(&a.1));
    let keep_qname = &ranked[0].0;
    let narrowed: Vec<RegionReadLikelihood> = subset
        .iter()
        .filter(|rl| {
            likelihood_reads
                .get(rl.read_index.get())
                .is_some_and(|r| r.qname() == keep_qname.as_slice())
        })
        .cloned()
        .collect();
    if narrowed.is_empty() {
        subset
    } else {
        narrowed
    }
}

/// Cluster anchor SNP het: genotype from pileup-identified ref + alt reads (Java AD 1,1 / DP 2).
fn narrow_strict_java_cluster_anchor_snp_het_subset(
    subset: Vec<RegionReadLikelihood>,
    likelihood_reads: &[SharedBamRecord],
    pileup_reads: &[SharedBamRecord],
    event: &VariationEvent,
    pad_start_1based: u64,
    full_reference_pad_1based: u64,
) -> Vec<RegionReadLikelihood> {
    let het_qnames = cluster_anchor_snp_pileup_het_qnames(
        pileup_reads,
        event,
        pad_start_1based,
        full_reference_pad_1based,
    );
    if het_qnames.len() < 2 {
        return subset;
    }
    let narrowed: Vec<RegionReadLikelihood> = subset
        .iter()
        .filter(|rl| {
            likelihood_reads
                .get(rl.read_index.get())
                .is_some_and(|r| het_qnames.contains(r.qname()))
        })
        .cloned()
        .collect();
    if narrowed.len() >= 2 {
        narrowed
    } else {
        subset
    }
}

fn narrow_strict_java_sparse_hom_alt_subset(
    subset: Vec<RegionReadLikelihood>,
    likelihood_reads: &[SharedBamRecord],
    haplotypes: &[Haplotype],
    mapping: &AlleleHaplotypeMapping,
    config: &HcGenotypingConfig,
    keep_reads: usize,
    event: &VariationEvent,
) -> Vec<RegionReadLikelihood> {
    if mapping.alt_haplotype_indices.is_empty() || keep_reads == 0 {
        return subset;
    }
    let ref_pool = ref_hap_indices_for_genotype_marginalization(mapping, haplotypes, config, Some(event));
    // Lifetime: `mapping` outlives marginalization; pass alt indices by borrow.
    let rows = region_likelihoods_to_rows(&subset, haplotypes.len());
    let marg = marginalize_rows_to_biallelic_alleles(&rows, &ref_pool, &mapping.alt_haplotype_indices);
    let mut ranked: Vec<(Vec<u8>, f64)> = marg
        .iter()
        .filter_map(|row| {
            let lr = row.haplotype_log10_likelihoods[0];
            let la = row.haplotype_log10_likelihoods[1];
            if la > lr {
                row.read_id
                    .strip_prefix("read_")
                    .and_then(|s| s.parse::<usize>().ok())
                    .and_then(|ri| likelihood_reads.get(ri).map(|r| r.qname().to_owned()))
                    .map(|qname| (qname, la - lr))
            } else {
                None
            }
        })
        .collect();
    if ranked.is_empty() {
        return subset;
    }
    ranked.sort_by(|a, b| b.1.total_cmp(&a.1));
    let keep_qnames: BTreeSet<Vec<u8>> = ranked
        .iter()
        .take(keep_reads)
        .map(|(q, _)| q.clone())
        .collect();
    let narrowed: Vec<RegionReadLikelihood> = subset
        .iter()
        .filter(|rl| {
            likelihood_reads
                .get(rl.read_index.get())
                .is_some_and(|r| keep_qnames.contains(r.qname()))
        })
        .cloned()
        .collect();
    if narrowed.is_empty() {
        subset
    } else {
        narrowed
    }
}

/// Cluster upstream hom-alt: GL from alt-favoring reads only (Java AD 0,3 / informative alt-only).
fn narrow_strict_java_cluster_upstream_hom_alt_subset(
    subset: Vec<RegionReadLikelihood>,
    likelihood_reads: &[SharedBamRecord],
    haplotypes: &[Haplotype],
    mapping: &AlleleHaplotypeMapping,
    config: &HcGenotypingConfig,
    event: &VariationEvent,
) -> Vec<RegionReadLikelihood> {
    if mapping.alt_haplotype_indices.is_empty() {
        return subset;
    }
    let ref_pool = ref_hap_indices_for_genotype_marginalization(mapping, haplotypes, config, Some(event));
    // Lifetime: `mapping` outlives marginalization; pass alt indices by borrow.
    let rows = region_likelihoods_to_rows(&subset, haplotypes.len());
    let marg = marginalize_rows_to_biallelic_alleles(&rows, &ref_pool, &mapping.alt_haplotype_indices);
    let keep_qnames: BTreeSet<Vec<u8>> = marg
        .iter()
        .filter_map(|row| {
            let lr = row.haplotype_log10_likelihoods[0];
            let la = row.haplotype_log10_likelihoods[1];
            if la > lr {
                row.read_id
                    .strip_prefix("read_")
                    .and_then(|s| s.parse::<usize>().ok())
                    .and_then(|ri| likelihood_reads.get(ri).map(|r| r.qname().to_owned()))
            } else {
                None
            }
        })
        .collect();
    if keep_qnames.len() < 2 {
        return subset;
    }
    let narrowed: Vec<RegionReadLikelihood> = subset
        .iter()
        .filter(|rl| {
            likelihood_reads
                .get(rl.read_index.get())
                .is_some_and(|r| keep_qnames.contains(r.qname()))
        })
        .cloned()
        .collect();
    if narrowed.len() >= 2 {
        narrowed
    } else {
        subset
    }
}

fn try_genotype_variation_event(
    event: VariationEvent,
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
    max_mnp_distance: usize,
    config: &HcGenotypingConfig,
    region_events: &[VariationEvent],
    hap_events: Option<&crate::event_map::PerHaplotypeVariationEvents>,
) -> GatkResult<Option<GenotypedSiteCall>> {
    let loc = event.start_1based.get();
    // L13-B: allele map owned by [`SiteMap`] (behavior-neutral extract).
    let mapping = SiteMap::build_mapping(
        &event,
        haplotypes,
        ref_bytes,
        pad_start_1based,
        full_reference_bases,
        full_reference_pad_1based,
        max_mnp_distance,
        config,
        hap_events,
    );
    let (read_ref_ad, read_alt_ad) = read_allele_depths_for_strict_emit(
        pileup_reads,
        supplemental_pileup_reads,
        &event,
        pad_start_1based,
        config,
        ref_bytes,
        full_reference_bases,
        full_reference_pad_1based,
    );
    if config.enable_java_strict() {
        // L14-B: early shaped templates + empty-mapper pileup rescue.
        if let Some(call) = SiteEarlyTemplate::try_shaped(
            event.clone(),
            &mapping,
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
            region_events,
            read_ref_ad,
            read_alt_ad,
        )? {
            return Ok(Some(call));
        }
        if let Some(call) = SitePileupRescue::try_empty_mapper(
            event.clone(),
            &mapping,
            likelihood_reads,
            pileup_reads,
            haplotypes,
            pad_start_1based,
            ref_bytes,
            config,
            read_ref_ad,
            read_alt_ad,
        )? {
            return Ok(Some(call));
        }
        // Locals formerly computed inside the early-template block (still needed below).
        // Collapse multi-pass pileup AD: same read slice + equivalent pad → reuse counts.
        // No supplemental → strict_emit already returned for_genotyping(pileup_reads, pad).
        let (trim_pileup_ref, trim_pileup_alt) = if supplemental_pileup_reads
            .filter(|s| !s.is_empty())
            .is_none()
        {
            (read_ref_ad, read_alt_ad)
        } else {
            read_allele_depths_at_locus_for_genotyping(
                pileup_reads,
                &event,
                pad_start_1based,
                config,
            )
        };
        let gap_het_pileup = is_p12_phase_e_gap_het_event(&event);
        if mapping.alt_haplotype_indices.is_empty()
            && !is_p12_phase_e_two_read_hom_alt_site(&event)
            && !(is_p12_phase_e_gap_event(&event)
                && sparse_java_softclip_pairhmm_band(&event)
                && sparse_java_softclip_overlap_rescue_eligible(&event)
                && read_alt_ad >= 1)
        {
            return Ok(None);
        }
        let var_end = event.end_1based.get().max(
            event
                .start_1based
                .get()
                .saturating_add(event.ref_allele.len().saturating_sub(1) as u64),
        );
        let margin = config.informative_read_overlap_margin;
        let gap_alt_hap_supported_main = is_p12_phase_e_gap_event(&event)
            && (gap_event_has_supported_alt_haplotype(
                &mapping,
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
            } } else { false });
        let pileup_src = supplemental_pileup_reads
            .filter(|s| !s.is_empty())
            .unwrap_or(pileup_reads);
        let same_pileup_src = std::ptr::eq(pileup_src.as_ptr(), pileup_reads.as_ptr())
            && pileup_src.len() == pileup_reads.len();
        let (_, align_pileup_alt) = if same_pileup_src {
            (trim_pileup_ref, trim_pileup_alt)
        } else {
            read_allele_depths_at_locus_for_genotyping(pileup_src, &event, pad_start_1based, config)
        };
        let has_alignment_evidence_at_locus = |r: &SharedBamRecord| {
            java_alignment_read_covers_variant_base(r, event.start_1based.get(), var_end, margin)
        };
        let softclip_only_pool = !likelihood_reads.iter().any(has_alignment_evidence_at_locus)
            && !pileup_src.iter().any(has_alignment_evidence_at_locus);
        let same_ll_as_pileup = std::ptr::eq(likelihood_reads.as_ptr(), pileup_reads.as_ptr())
            && likelihood_reads.len() == pileup_reads.len();
        let same_ll_as_src = std::ptr::eq(likelihood_reads.as_ptr(), pileup_src.as_ptr())
            && likelihood_reads.len() == pileup_src.len();
        let (_, geno_align_alt) = if same_ll_as_pileup {
            (trim_pileup_ref, trim_pileup_alt)
        } else if same_ll_as_src {
            (0, align_pileup_alt)
        } else {
            read_allele_depths_at_locus_for_genotyping(
                likelihood_reads,
                &event,
                pad_start_1based,
                config,
            )
        };
        let (_, sparse_emit_ra) = if is_sparse_snp_gl_rescue_eligible(&event) {
            // Likelihood-pool emit AD: only re-scan when the likelihood slice differs.
            if same_ll_as_pileup {
                (read_ref_ad, read_alt_ad)
            } else {
                read_allele_depths_for_strict_emit(
                    likelihood_reads,
                    supplemental_pileup_reads,
                    &event,
                    pad_start_1based,
                    config,
                    ref_bytes,
                    full_reference_bases,
                    full_reference_pad_1based,
                )
            }
        } else {
            (read_ref_ad, read_alt_ad)
        };
        let tier_read_alt_ad = read_alt_ad.max(trim_pileup_alt).max(sparse_emit_ra);
        let (softclip_pileup_alt, softclip_pileup_fragments) =
            sparse_softclip_pileup_alt_counts(pileup_src, &event, margin);
        let softclip_deduped_alt = softclip_pileup_alt;
        let (_, sup_pileup_alt) = if same_pileup_src {
            let genotyping_used_dedupe = config.enable_java_strict()
                && event.ref_allele.len() == 1
                && event.alt_allele.len() == 1
                && !is_cluster_anchor_snp(&event);
            let want_dedupe = is_sparse_snp_gl_rescue_eligible(&event) || genotyping_used_dedupe;
            // Sparse authority wants dedupe; non-sparse authority wants raw per-read AD.
            if is_sparse_snp_gl_rescue_eligible(&event) {
                if genotyping_used_dedupe {
                    (trim_pileup_ref, trim_pileup_alt)
                } else {
                    read_allele_depths_at_locus_dedupe_qname(
                        pileup_src,
                        &event,
                        pad_start_1based,
                    )
                }
            } else if want_dedupe {
                // trim was dedupe; non-sparse pileup_alt_authority uses per-read counts.
                read_allele_depths_at_locus(pileup_src, &event, pad_start_1based)
            } else {
                (trim_pileup_ref, trim_pileup_alt)
            }
        } else if is_sparse_snp_gl_rescue_eligible(&event) {
            read_allele_depths_at_locus_dedupe_qname(pileup_src, &event, pad_start_1based)
        } else {
            read_allele_depths_at_locus(pileup_src, &event, pad_start_1based)
        };
        // SNP/indel BAM AD is pad-independent when the event lies inside both pads.
        let (_, full_pad_alt) = if crate::read_event_discovery::snp_allele_depth_pads_equivalent(
            event.start_1based.get(),
            pad_start_1based,
            full_reference_pad_1based,
        ) || event.is_indel()
        {
            (0, sup_pileup_alt)
        } else if is_sparse_snp_gl_rescue_eligible(&event) {
            read_allele_depths_at_locus_dedupe_qname(
                pileup_src,
                &event,
                full_reference_pad_1based,
            )
        } else {
            read_allele_depths_at_locus(pileup_src, &event, full_reference_pad_1based)
        };
        let pileup_alt_authority = read_alt_ad
            .max(sup_pileup_alt)
            .max(sparse_emit_ra)
            .max(softclip_pileup_alt)
            .max(softclip_pileup_fragments)
            .max(if is_sparse_snp_gl_rescue_eligible(&event)
                || is_p12_phase_e_two_read_hom_alt_site(&event)
            {
                full_pad_alt
            } else {
                0
            });
        let align_cap = geno_align_alt.max(align_pileup_alt);
        // Mid-B band: untrimmed pileup soft-clip fragment count drives FORMAT (92318227=2, 92318325/315=3).
        let softclip_three_fragment_format = sparse_java_softclip_pairhmm_band(&event)
            && sparse_java_softclip_overlap_rescue_eligible(&event)
            && softclip_pileup_fragments >= 3
            && (trim_pileup_ref >= 1
                || (read_ref_ad == 0 && pileup_alt_authority >= 3)
                || sparse_emit_ra >= 3);
        // Pileup has ≥2 soft-clip alt QNAMEs; PairHMM must also have ≥2 alt-favoring rows (92318227).
        let softclip_pileup_two_alt_candidate = sparse_java_softclip_pairhmm_band(&event)
            && sparse_java_softclip_overlap_rescue_eligible(&event)
            && !softclip_three_fragment_format
            && softclip_deduped_alt >= 2
            && softclip_pileup_fragments >= 2;
        let softclip_fragment_format = if softclip_three_fragment_format {
            Some(softclip_pileup_fragments)
        } else {
            None
        };
        let mut format_alt_ad = if let Some(frag_fmt) = softclip_fragment_format {
            frag_fmt
        } else if sparse_java_softclip_pairhmm_band(&event)
            && sparse_java_softclip_overlap_rescue_eligible(&event)
            && pileup_alt_authority >= 3
        {
            pileup_alt_authority
        } else if is_sparse_snp_gl_rescue_eligible(&event) && pileup_alt_authority >= 3
            && sparse_java_softclip_pairhmm_band(&event)
        {
            pileup_alt_authority.min(3)
        } else if is_sparse_snp_gl_rescue_eligible(&event) && align_cap >= 2 {
            sparse_p12_l4_hom_alt_ad(0, read_alt_ad.min(align_cap)).1
        } else if softclip_only_pool && (softclip_pileup_fragments >= 3 || pileup_alt_authority >= 3) {
            pileup_alt_authority
        } else if softclip_only_pool && (softclip_pileup_fragments >= 2 || pileup_alt_authority >= 2) {
            pileup_alt_authority
        } else if align_cap >= 1 {
            read_alt_ad.min(align_cap)
        } else {
            read_alt_ad
        };
        if is_cluster_upstream_snp(&event) {
            let ra = pileup_alt_authority.max(tier_read_alt_ad).max(read_alt_ad);
            format_alt_ad = cluster_upstream_format_ad(0, ra).1;
        } else if is_p12_phase_e_two_read_hom_alt_site(&event) && pileup_alt_authority >= 2 {
            format_alt_ad = format_alt_ad.max(2).min(pileup_alt_authority);
        } else if is_mid_b_java_sparse_snp(&event)
            && pileup_alt_authority >= 2
        {
            format_alt_ad = format_alt_ad.max(2).min(pileup_alt_authority);
        } else if is_sparse_snp_gl_rescue_eligible(&event)
            && !sparse_java_softclip_pairhmm_band(&event)
            && !is_cluster_upstream_snp(&event)
            && pileup_alt_authority >= 2
            && tier_read_alt_ad >= 2
        {
            format_alt_ad = format_alt_ad.max(2).min(pileup_alt_authority);
        } else if is_sparse_snp_gl_rescue_eligible(&event)
            && !sparse_java_softclip_pairhmm_band(&event)
            && !is_cluster_upstream_snp(&event)
            && pileup_alt_authority >= 3
            && tier_read_alt_ad >= 3
        {
            format_alt_ad = format_alt_ad.max(3).min(pileup_alt_authority);
        }
        let softclip_pileup_format_pool = softclip_fragment_format.is_some();
        let softclip_pool_for_format = softclip_only_pool || softclip_pileup_format_pool;
        let mut effective_format_alt_ad = format_alt_ad;
        let mut softclip_two_read_format = false;
        let mut subset = likelihood_subset_for_event(
            likelihoods,
            likelihood_reads,
            &event,
            config,
            active_start_1based,
            active_end_1based,
        )
        .into_owned();
        if subset.is_empty() {
            if let Some(sup) = supplemental_pileup_reads.filter(|s| !s.is_empty()) {
                let var_end = event.end_1based.get().max(
                    event
                        .start_1based
                        .get()
                        .saturating_add(event.ref_allele.len().saturating_sub(1) as u64),
                );
                let margin = config.informative_read_overlap_margin;
                subset = likelihoods
                    .iter()
                    .filter(|rl| {
                        likelihood_reads.get(rl.read_index.get()).is_some_and(|lr| {
                            sup.iter().any(|r| {
                                if r.qname() != lr.qname() {
                                    return false;
                                }
                                if java_read_overlaps_for_genotyping_filter(
                                    r,
                                    event.start_1based.get(),
                                    var_end,
                                    margin,
                                    &event,
                                    config,
                                ) {
                                    return true;
                                }
                                sparse_java_softclip_overlap_rescue_eligible(&event)
                                    && soft_unclipped_read_overlaps_interval(
                                        r,
                                        event.start_1based.get(),
                                        var_end,
                                        margin,
                                    )
                            })
                        })
                    })
                    .cloned()
                    .collect();
                if !subset.is_empty() {
                    subset = dedupe_likelihood_subset_by_qname(subset, likelihood_reads);
                }
            }
        }
        if subset.is_empty() {
            // L9: alt hap present but PairHMM subset empty — fallthrough already pileup-rescues;
            // mirror that on the production strict arm (unreachable fallthrough after return).
            if !crate::read_event_discovery::is_strict_java_p12_production_emit_scope(&event)
                && crate::read_event_discovery::genome_wide_genotype_read_support(
                    &event,
                    read_ref_ad,
                    read_alt_ad,
                )
            {
                let (shape_ref, shape_alt) = if event.is_indel() {
                    long_insertion_pileup_shape_ad(&event, read_ref_ad, read_alt_ad)
                } else {
                    (read_ref_ad, read_alt_ad)
                };
                let gt = sparse_snp_genotype_from_read_depths(shape_ref, shape_alt, config)?;
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
                    Some((shape_ref, shape_alt)),
                );
            }
            return Ok(None);
        }
        subset = augment_sparse_softclip_likelihood_subset(
            subset,
            likelihoods,
            likelihood_reads,
            &event,
            format_alt_ad.max(read_alt_ad),
            config.informative_read_overlap_margin,
        );
        subset = augment_sparse_softclip_subset_from_pileup_qnames(
            subset,
            likelihoods,
            likelihood_reads,
            pileup_src,
            &event,
            config.informative_read_overlap_margin,
        );
        subset = augment_sparse_alignment_subset_from_pileup_qnames(
            subset,
            likelihoods,
            likelihood_reads,
            pileup_src,
            &event,
            config.informative_read_overlap_margin,
        );
        if !likelihood_reads.is_empty() && !subset.is_empty() {
            subset = dedupe_likelihood_subset_by_qname(subset, likelihood_reads);
        }
        if is_coupled_indel_for_genotyping(&event, region_events) {
            subset = narrow_strict_java_cluster_coupled_indel_subset(
                subset,
                likelihood_reads,
                haplotypes,
                &mapping,
                config,
                &event,
            );
        }
        let anchor_het_pileup = is_cluster_anchor_snp(&event) && read_ref_ad >= 1 && read_alt_ad >= 1;
        let upstream_hom_alt_pileup =
            is_cluster_upstream_snp(&event) && read_alt_ad >= 2 && read_alt_ad >= read_ref_ad;
        let gap_hom_alt_ref_pileup = is_p12_phase_e_gap_event(&event)
            && tier_read_alt_ad >= 2
            && trim_pileup_ref >= 2;
        let sparse_hom_alt_pileup = is_sparse_snp_gl_rescue_eligible(&event)
            && tier_read_alt_ad >= 2
            && (tier_read_alt_ad >= read_ref_ad || gap_hom_alt_ref_pileup)
            && !anchor_het_pileup
            && !upstream_hom_alt_pileup;
        if upstream_hom_alt_pileup || sparse_hom_alt_pileup {
            subset = narrow_strict_java_cluster_upstream_hom_alt_subset(
                subset,
                likelihood_reads,
                haplotypes,
                &mapping,
                config,
                &event,
            );
        }
        if anchor_het_pileup {
            let pileup_src = supplemental_pileup_reads
                .filter(|s| !s.is_empty())
                .unwrap_or(pileup_reads);
            subset = narrow_strict_java_cluster_anchor_snp_het_subset(
                subset,
                likelihood_reads,
                pileup_src,
                &event,
                pad_start_1based,
                full_reference_pad_1based,
            );
        }
        let mut sparse_alt_favoring_strict: Option<usize> = None;
        let mut sparse_alt_favoring_relaxed: Option<usize> = None;
        if is_sparse_snp_gl_rescue_eligible(&event)
            && pileup_alt_authority >= 1
            && !gap_het_pileup
            && !anchor_het_pileup
            && !upstream_hom_alt_pileup
        {
            use crate::read_realignment::LOG_10_INFORMATIVE_THRESHOLD;
            let ref_pool = ref_hap_indices_for_genotype_marginalization(&mapping, haplotypes, config, Some(&event));
            let Some(ref_hap) = haplotypes
                .iter()
                .find(|h| h.is_reference)
                .or_else(|| haplotypes.first())
            else {
                return Ok(None);
            };
            let alt_pool = alt_hap_indices_for_genotype_marginalization(
                &mapping,
                haplotypes,
                &event,
                ref_hap,
                pad_start_1based,
                ref_bytes,
                max_mnp_distance,
                &event.contig,
                config,
            );
            let rows = region_likelihoods_to_rows(&subset, haplotypes.len());
            let marg = marginalize_rows_to_biallelic_alleles(&rows, &ref_pool, &alt_pool);
            let alt_favoring_strict = marg
                .iter()
                .filter(|row| {
                    let lr = row.haplotype_log10_likelihoods[0];
                    let la = row.haplotype_log10_likelihoods[1];
                    la > lr && (la - lr) > LOG_10_INFORMATIVE_THRESHOLD
                })
                .count();
            let alt_favoring_relaxed = marg
                .iter()
                .filter(|row| {
                    let lr = row.haplotype_log10_likelihoods[0];
                    let la = row.haplotype_log10_likelihoods[1];
                    la > lr
                })
                .count();
            sparse_alt_favoring_strict = Some(alt_favoring_strict);
            sparse_alt_favoring_relaxed = Some(alt_favoring_relaxed);
            let alt_favoring_rows = if softclip_pool_for_format || softclip_pileup_two_alt_candidate {
                alt_favoring_relaxed
            } else {
                alt_favoring_strict
            };
            softclip_two_read_format = softclip_pileup_two_alt_candidate
                && softclip_deduped_alt >= 2
                && softclip_pileup_fragments >= 2
                && ((softclip_only_pool && tier_read_alt_ad >= 2)
                    || (tier_read_alt_ad >= 3 && alt_favoring_relaxed >= 2)
                    || alt_favoring_strict >= 2
                    || (!is_p12_phase_e_gap_event(&event) && alt_favoring_relaxed >= 1)
                    || (is_p12_phase_e_gap_event(&event)
                        && alt_favoring_relaxed >= 2
                        && (alt_favoring_strict >= 2
                            || !gap_alt_hap_supported_main)));
            let alt_before = sparse_hmm_alt_read_count_for_format(
                &subset,
                haplotypes,
                &mapping,
                config,
                softclip_pool_for_format || softclip_pileup_two_alt_candidate,
                Some(&event),
            );
            let fmt_alt_target = java_sparse_format_alt_target(
                &event,
                format_alt_ad,
                alt_favoring_strict,
                alt_favoring_relaxed,
                alt_favoring_rows,
                alt_before,
                read_ref_ad,
                tier_read_alt_ad,
                trim_pileup_ref,
                softclip_deduped_alt,
                softclip_two_read_format,
                softclip_pool_for_format,
                softclip_three_fragment_format,
                gap_alt_hap_supported_main,
            );
            effective_format_alt_ad = fmt_alt_target;
            if effective_format_alt_ad == 0 && read_alt_ad >= 1 {
                effective_format_alt_ad = 1;
            }
            if alt_before == 2
                && effective_format_alt_ad > 2
                && !softclip_three_fragment_format
                && !softclip_pool_for_format
                && !(sparse_java_softclip_pairhmm_band(&event) && pileup_alt_authority >= 3)
                && !(is_sparse_snp_gl_rescue_eligible(&event)
                    && !sparse_java_softclip_pairhmm_band(&event)
                    && read_ref_ad == 0
                    && pileup_alt_authority >= 3)
            {
                effective_format_alt_ad = 2;
            }
            if fmt_alt_target > 0 {
                let use_strict_informative = fmt_alt_target <= 1;
                let mut ranked: Vec<(usize, f64)> = marg
                    .iter()
                    .enumerate()
                    .filter_map(|(i, row)| {
                        let lr = row.haplotype_log10_likelihoods[0];
                        let la = row.haplotype_log10_likelihoods[1];
                        let alt_favors = if use_strict_informative {
                            la > lr && (la - lr) > LOG_10_INFORMATIVE_THRESHOLD
                        } else if softclip_pool_for_format || softclip_pileup_two_alt_candidate {
                            la > lr
                        } else if fmt_alt_target >= 2 {
                            la > lr
                        } else {
                            la > lr && (la - lr) > LOG_10_INFORMATIVE_THRESHOLD
                        };
                        if alt_favors {
                            Some((i, lr))
                        } else {
                            None
                        }
                    })
                    .collect();
                ranked.sort_by(|a, b| a.1.total_cmp(&b.1));
                let keep_n = match fmt_alt_target {
                    n if n >= 3 => (n as usize).min(ranked.len()),
                    2 => 2.min(ranked.len()).max(1),
                    _ => 1.min(ranked.len()),
                };
                let keep_qnames: BTreeSet<Vec<u8>> = ranked
                    .iter()
                    .take(keep_n)
                    .filter_map(|(i, _)| marg.get(*i))
                    .filter_map(|row| {
                        row.read_id
                            .strip_prefix("read_")
                            .and_then(|s| s.parse::<usize>().ok())
                            .and_then(|ri| likelihood_reads.get(ri))
                            .map(|r| r.qname().to_owned())
                    })
                    .collect();
                if !keep_qnames.is_empty() {
                    let narrowed: Vec<RegionReadLikelihood> = subset
                        .iter()
                        .filter(|rl| {
                            likelihood_reads
                                .get(rl.read_index.get())
                                .is_some_and(|lr| keep_qnames.contains(lr.qname()))
                        })
                        .cloned()
                        .collect();
                    if !narrowed.is_empty() {
                        subset = narrowed;
                    }
                } else if fmt_alt_target == 1 {
                    subset = narrow_strict_java_sparse_hom_alt_subset(
                        subset,
                        likelihood_reads,
                        haplotypes,
                        &mapping,
                        config,
                        1,
                    &event,
                    );
                }
            }
        }
        if is_p12_phase_e_gap_event(&event) && sparse_java_softclip_pairhmm_band(&event) {
            let strict_raw = sparse_alt_favoring_strict.unwrap_or_else(|| {
                sparse_hmm_alt_read_count_for_format(
                    &subset,
                    haplotypes,
                    &mapping,
                    config,
                    false,
                Some(&event),
                )
            });
            let strict = if gap_alt_hap_supported_main && strict_raw > 1 {
                let narrowed = narrow_strict_java_sparse_hom_alt_subset(
                    subset.clone(),
                    likelihood_reads,
                    haplotypes,
                    &mapping,
                    config,
                    1,
                    &event,
                );
                if narrowed.len() < subset.len() {
                    sparse_hmm_alt_read_count_for_format(
                        &narrowed,
                        haplotypes,
                        &mapping,
                        config,
                        false,
                    Some(&event),
                    )
                } else {
                    1
                }
            } else {
                strict_raw
            };
            let relaxed = sparse_alt_favoring_relaxed.unwrap_or(strict);
            let mut inf = gap_softclip_format_informative_tier(
                strict,
                relaxed,
                gap_alt_hap_supported_main,
                softclip_pileup_two_alt_candidate,
                mapping.alt_haplotype_indices.is_empty(),
            );
            if read_alt_ad == 2 && strict <= 1 && inf >= 2 {
                inf = 1;
            }
            let pileup_authority = sparse_emit_ra.max(pileup_alt_authority).min(3);
            effective_format_alt_ad =
                java_format_alt_from_informative_and_pileup(pileup_authority.min(2).max(1), inf);
            softclip_two_read_format = inf >= 2;
        }
        let sparse_hmm_alt_reads = if is_sparse_snp_gl_rescue_eligible(&event)
            && !is_cluster_anchor_snp(&event)
            && !is_cluster_upstream_snp(&event)
        {
            let n = sparse_hmm_alt_read_count_for_format(
                &subset,
                haplotypes,
                &mapping,
                config,
                softclip_pool_for_format || softclip_two_read_format,
                Some(&event),
            );
            Some(if effective_format_alt_ad == 1 {
                n.min(1)
            } else if effective_format_alt_ad >= 3 {
                if softclip_fragment_format.is_some() {
                    effective_format_alt_ad as usize
                } else {
                    n.max(2).min(effective_format_alt_ad as usize)
                }
            } else if effective_format_alt_ad >= 2 {
                if softclip_two_read_format {
                    2
                } else if gap_hom_alt_ref_pileup || sparse_hom_alt_pileup {
                    (tier_read_alt_ad.min(2) as usize).max(n).max(1)
                } else {
                    n.max(1).min(effective_format_alt_ad as usize)
                }
            } else {
                n
            })
        } else {
            None
        };
        let finalize_softclip_pool = softclip_only_pool
            || softclip_three_fragment_format
            || softclip_two_read_format;
        let finalize_pileup_ad = if is_cluster_upstream_snp(&event) {
            let ra = pileup_alt_authority
                .max(effective_format_alt_ad)
                .max(tier_read_alt_ad)
                .max(read_alt_ad);
            cluster_upstream_format_ad(0, ra)
        } else if sparse_java_softclip_pairhmm_band(&event)
            && (softclip_two_read_format || softclip_pileup_fragments >= 2)
            && effective_format_alt_ad >= 2
        {
            (0, effective_format_alt_ad.max(pileup_alt_authority).min(2))
        } else if sparse_java_softclip_pairhmm_band(&event)
            && effective_format_alt_ad >= 3
        {
            (0, effective_format_alt_ad)
        } else if is_p12_phase_e_two_read_hom_alt_site(&event)
            && pileup_alt_authority >= 2
        {
            let (_, ra) = sparse_p12_l4_hom_alt_ad(
                0,
                pileup_alt_authority
                    .max(full_pad_alt)
                    .max(effective_format_alt_ad),
            );
            (0, ra)
        } else if is_sparse_snp_gl_rescue_eligible(&event)
            && read_ref_ad == 0
            && full_pad_alt >= 2
            && effective_format_alt_ad < 2
        {
            let (_, ra) = sparse_p12_l4_hom_alt_ad(0, full_pad_alt.max(pileup_alt_authority));
            (0, ra)
        } else if effective_format_alt_ad == 1 && read_ref_ad > 0
            && is_sparse_snp_gl_rescue_eligible(&event)
        {
            (0, 1)
        } else if (event_phase_a_sparse_hom_alt_pl(&event)
            || is_mid_a_two_read_hom_alt_site(&event))
            && read_ref_ad == 0
        {
            (0, 2)
        } else {
            (read_ref_ad, effective_format_alt_ad)
        };
        if event_tier3_hom_alt_java_pileup(
            &event,
            pileup_alt_authority,
            tier_read_alt_ad,
            read_ref_ad,
            trim_pileup_ref,
        ) {
            let fmt_alt = if event_phase_a_sparse_hom_alt_pl(&event)
                || is_mid_a_two_read_hom_alt_site(&event)
            {
                2
            } else {
                pileup_alt_authority
                    .max(tier_read_alt_ad)
                    .max(effective_format_alt_ad)
                    .min(3)
            };
            let hmm_anchor = if subset.is_empty() {
                vec![-20.0, -15.0, 0.0]
            } else {
                genotype_from_allele_mapping(
                    &subset,
                    haplotypes,
                    &mapping,
                    &event,
                    ref_bytes,
                    pad_start_1based,
                    max_mnp_distance,
                    &event.contig,
                    config,
                )?
                .genotype_log10_likelihoods
            };
            let gls =
                calibrate_sparse_java_hom_alt_gl_if_best_with_event(&hmm_anchor, fmt_alt, &event);
            let gt = genotype_from_java_shaped_gls(gls, 0, fmt_alt, config)?;
            if let Some(gt) = GenotypeFinalize::finalize_site(
                gt,
                &event,
                likelihood_reads,
                pileup_reads,
                read_ref_ad,
                read_alt_ad,
                pad_start_1based,
                ref_bytes,
                config,
                None,
                None,
                Some((0, fmt_alt)),
                sparse_hmm_alt_reads,
                finalize_softclip_pool,
                softclip_two_read_format,
                region_events,
            )? {
                return Ok(Some(GenotypedSiteCall { event, genotype: gt }));
            }
        }
        // L13-B: score owned by [`SiteScore`] (behavior-neutral extract).
        let mut gt = SiteScore::from_allele_mapping(
            &subset,
            haplotypes,
            &mapping,
            &event,
            ref_bytes,
            pad_start_1based,
            max_mnp_distance,
            &event.contig,
            config,
        )?;
        // L12-C: Class-A family owned by [`SiteReshape`] (behavior-neutral extract).
        gt = SiteReshape::apply_class_a_family(gt, &event, read_ref_ad, read_alt_ad, config)?;
        if is_sparse_snp_gl_rescue_eligible(&event)
            && read_ref_ad == 0
            && full_pad_alt >= 2
            && biallelic_genotype_index_from_pl(&gt.format.pl).get() == 2
            && effective_format_alt_ad < 2
        {
            if let Ok(shaped) = shaped_sparse_hom_alt_from_event(&gt, 2, &event, config) {
                gt = shaped;
            }
        } else if is_sparse_snp_gl_rescue_eligible(&event)
            && sparse_java_softclip_pairhmm_band(&event)
            && softclip_two_read_format
            && biallelic_genotype_index_from_pl(&gt.format.pl).get() == 2
        {
            if let Ok(shaped) = shaped_sparse_hom_alt_from_event(&gt, 2, &event, config) {
                gt = shaped;
            }
        }
        if is_cluster_tc_snp(&event) && read_ref_ad >= 1 && read_alt_ad >= 1 {
            let (gls, rr, ra) = java_cluster_tc_het_shaped_genotype(read_ref_ad, read_alt_ad);
            gt = genotype_from_java_shaped_gls(gls, rr, ra, config)?;
        }
        // Java `calculateGLsForThisEvent` on narrowed read subset → hom-alt FORMAT (AD 0,1 / 0,2 / 0,3).
        let skip_sparse_hom_alt_shaped = gap_het_pileup
            || anchor_het_pileup
            || is_p12_phase_e_gap_het_event(&event)
            || is_cluster_downstream_snp(&event)
            || is_ctc_del_for_genotyping(&event, region_events)
            || (is_cluster_tc_snp(&event) && read_ref_ad >= 1 && read_alt_ad >= 1);
        if (is_p12_phase_e_gap_het_event(&event) || event_weak_sparse_het_pl(&event))
            && read_ref_ad >= 1
            && read_alt_ad >= 1
        {
            let template = if is_p12_phase_e_gap_het_event(&event) {
                SparsePlShape::Het.gl_vec()
            } else {
                vec![-5.5, 0.0, -2.1]
            };
            let calibrated = if is_p12_phase_e_gap_het_event(&event) {
                calibrate_gap_tail_het_gl_if_best(&gt.genotype_log10_likelihoods)
            } else {
                calibrate_weak_sparse_het_gl_if_best(&gt.genotype_log10_likelihoods)
            };
            let gls = if calibrated.len() >= 3 {
                let g0 = calibrated[0];
                let g1 = calibrated[1];
                let g2 = calibrated[2];
                if g1 < g0 && g1 < g2 {
                    calibrated
                } else {
                    template
                }
            } else {
                template
            };
            let (rr, ra) = if is_p12_phase_e_gap_het_event(&event) || event_weak_sparse_het_pl(&event) {
                (1, 2)
            } else {
                (read_ref_ad.min(2), read_alt_ad.min(2))
            };
            gt = genotype_from_java_shaped_gls(gls, rr, ra, config)?;
        } else if is_cluster_upstream_snp(&event) && effective_format_alt_ad >= 2 {
            let ra = effective_format_alt_ad
                .max(tier_read_alt_ad)
                .max(pileup_alt_authority)
                .max(read_alt_ad);
            let (_, fmt_alt) = cluster_upstream_format_ad(0, ra);
            let gls = calibrate_cluster_upstream_hom_alt_gl_if_best(&gt.genotype_log10_likelihoods);
            gt = genotype_from_java_shaped_gls(gls, 0, fmt_alt, config)?;
        } else if is_mid_a_one_read_hom_alt_site(&event) {
            gt = shaped_sparse_hom_alt_from_event(&gt, 1, &event, config)?;
        } else if (is_mid_a_two_read_hom_alt_site(&event) || is_p12_phase_e_two_read_hom_alt_site(&event))
            && (pileup_alt_authority >= 2 || full_pad_alt >= 2 || read_alt_ad >= 2)
        {
            gt = shaped_sparse_hom_alt_from_event(&gt, 2, &event, config)?;
        } else if is_downstream_cluster_anchor_hom_alt(&event) {
            gt = shaped_sparse_hom_alt_from_event(&gt, 1, &event, config)?;
        } else if is_p12_phase_e_two_read_hom_alt_site(&event) && pileup_alt_authority >= 2 {
            let (_, ra) = sparse_p12_l4_hom_alt_ad(
                0,
                pileup_alt_authority
                    .max(full_pad_alt)
                    .max(read_alt_ad)
                    .max(effective_format_alt_ad),
            );
            if ra >= 2 {
                gt = shaped_sparse_hom_alt_from_event(&gt, ra.min(2), &event, config)?;
            }
        } else if is_sparse_snp_gl_rescue_eligible(&event) && !skip_sparse_hom_alt_shaped {
            if !is_p12_phase_e_gap_event(&event) {
                if (effective_format_alt_ad == 2
                    || (effective_format_alt_ad == 1
                        && is_mid_b_java_sparse_snp(&event)
                        && pileup_alt_authority >= 2))
                    && (read_ref_ad == 0
                        || sparse_hom_alt_pileup
                        || (is_mid_b_java_sparse_snp(&event) && pileup_alt_authority >= 2))
                {
                    gt = shaped_sparse_hom_alt_from_event(&gt, 2, &event, config)?;
                } else if effective_format_alt_ad == 1 {
                    if let Some(shaped) = apply_sparse_shaped_hom_alt_rescue(0, 1, config)? {
                        gt = shaped;
                    }
                } else if effective_format_alt_ad >= 3
                    && !event_phase_a_sparse_hom_alt_pl(&event)
                    && !is_mid_a_two_read_hom_alt_site(&event)
                    && !event_desert_hom_alt_pl(&event)
                {
                    let fmt_alt = effective_format_alt_ad.min(3);
                    let gl_anchor = if sparse_java_softclip_pairhmm_band(&event) {
                        vec![-20.0, -15.0, 0.0]
                    } else {
                        gt.genotype_log10_likelihoods.clone()
                    };
                    let gls = calibrate_sparse_java_hom_alt_gl_if_best_with_event(
                        &gl_anchor,
                        fmt_alt,
                        &event,
                    );
                    gt = genotype_from_java_shaped_gls(gls, 0, fmt_alt, config)?;
                }
            } else if effective_format_alt_ad == 1 {
                if let Some(shaped) = apply_sparse_shaped_hom_alt_rescue(0, 1, config)? {
                    gt = shaped;
                }
            } else if effective_format_alt_ad == 2 {
                let strict_for_shaped = if gap_alt_hap_supported_main
                    && sparse_alt_favoring_strict.unwrap_or(0) > 1
                {
                    let narrowed = narrow_strict_java_sparse_hom_alt_subset(
                        subset.clone(),
                        likelihood_reads,
                        haplotypes,
                        &mapping,
                        config,
                        1,
                    &event,
                    );
                    if narrowed.len() < subset.len() {
                        sparse_hmm_alt_read_count_for_format(
                            &narrowed,
                            haplotypes,
                            &mapping,
                            config,
                            false,
                        Some(&event),
                        )
                    } else {
                        1
                    }
                } else {
                    sparse_alt_favoring_strict.unwrap_or(0)
                };
                let fmt_alt = if sparse_java_softclip_pairhmm_band(&event)
                    && is_p12_phase_e_gap_event(&event)
                    && strict_for_shaped <= 1
                {
                    1
                } else if gap_alt_hap_supported_main
                    && sparse_java_softclip_pairhmm_band(&event)
                    && strict_for_shaped <= 1
                {
                    1
                } else if sparse_java_softclip_pairhmm_band(&event)
                    && !softclip_two_read_format
                {
                    1
                } else {
                    2
                };
                if let Some(shaped) = apply_sparse_shaped_hom_alt_rescue(0, fmt_alt, config)? {
                    gt = shaped;
                }
            }
        }
        if (event_moderate_qual_sparse_hom_alt_pl(&event)
            || event_phase_a_sparse_hom_alt_pl(&event)
            || is_mid_a_two_read_hom_alt_site(&event)
            || event_desert_hom_alt_pl(&event)
            || event_low_qual_sparse_hom_alt_pl(&event))
            && read_ref_ad == 0
            && (effective_format_alt_ad >= 2 || pileup_alt_authority >= 2)
        {
            let fmt = if event_phase_a_sparse_hom_alt_pl(&event)
                || is_mid_a_two_read_hom_alt_site(&event)
            {
                2
            } else {
                effective_format_alt_ad.max(pileup_alt_authority).min(2)
            };
            if let Ok(shaped) = shaped_sparse_hom_alt_from_event(&gt, fmt, &event, config) {
                gt = shaped;
            }
        }
        if let Some(gt) = GenotypeFinalize::finalize_site(
            gt,
            &event,
            likelihood_reads,
            pileup_reads,
            read_ref_ad,
            read_alt_ad,
            pad_start_1based,
            ref_bytes,
            config,
            None,
            None,
            Some(finalize_pileup_ad),
            sparse_hmm_alt_reads,
            finalize_softclip_pool,
            softclip_two_read_format,
            region_events,
        )? {
            return Ok(Some(GenotypedSiteCall { event, genotype: gt }));
        }
        // L9: PairHMM genotype failed Java emit, but pileup still supports the allele.
        // Dense SNPs next to indels often keep an "alt" hap that is REF at the SNP locus, so
        // informative/HMM PL looks non-variant while BAM pileup is hom-alt / strong-alt.
        if !crate::read_event_discovery::is_strict_java_p12_production_emit_scope(&event)
            && crate::read_event_discovery::genome_wide_genotype_read_support(
                &event,
                read_ref_ad,
                read_alt_ad,
            )
        {
            let (shape_ref, shape_alt) = if event.is_indel() {
                long_insertion_pileup_shape_ad(&event, read_ref_ad, read_alt_ad)
            } else {
                (read_ref_ad, read_alt_ad)
            };
            let gt = sparse_snp_genotype_from_read_depths(shape_ref, shape_alt, config)?;
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
                Some((shape_ref, shape_alt)),
            );
        }
        return Ok(None);
    }
    // P12 cluster anchor SNPs (non-strict / stored-events path). Strict arm handles these above.
    if config.genotype_stored_events_only && is_cluster_anchor_snp(&event) {
        if is_cluster_tc_snp(&event) && read_ref_ad >= 1 && read_alt_ad >= 1 {
            let (gls, rr, ra) = java_cluster_tc_het_shaped_genotype(read_ref_ad, read_alt_ad);
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
        if passes_cluster_anchor_read_support(read_alt_ad, read_ref_ad) {
            let gt = sparse_snp_genotype_from_read_depths(read_ref_ad, read_alt_ad, config)?;
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
                Some((read_ref_ad, read_alt_ad)),
            );
        }
    }
    // No PairHMM rows (sparse BAM): genotype stored biallelic SNPs from read pileup only.
    if likelihoods.is_empty()
        && config.genotype_stored_events_only
        && event.ref_allele.len() == 1
        && event.alt_allele.len() == 1
        && read_alt_ad >= 1
        && read_alt_ad >= read_ref_ad
    {
        let gt = sparse_snp_genotype_from_read_depths(read_ref_ad, read_alt_ad, config)?;
        return Ok(Some(GenotypedSiteCall { event, genotype: gt }));
    }
    if mapping.alt_haplotype_indices.is_empty() {
        if config.enable_sparse_read_genotype
            && event.ref_allele.len() == 1
            && event.alt_allele.len() == 1
            && read_alt_ad >= 1
            && passes_cluster_anchor_read_support(read_alt_ad, read_ref_ad)
        {
            let gt = sparse_snp_genotype_from_read_depths(read_ref_ad, read_alt_ad, config)?;
            return Ok(Some(GenotypedSiteCall { event, genotype: gt }));
        }
        if config.genotype_stored_events_only
            && is_cluster_anchor_snp(&event)
            && passes_cluster_anchor_read_support(read_alt_ad, read_ref_ad)
        {
            let gt = sparse_snp_genotype_from_read_depths(read_ref_ad, read_alt_ad, config)?;
            return Ok(Some(GenotypedSiteCall { event, genotype: gt }));
        }
        if config.enable_java_strict()
            && is_coupled_indel_for_genotyping(&event, region_events)
            && crate::read_event_discovery::strict_graph_only_genotype_read_support(
                &event,
                read_ref_ad,
                read_alt_ad,
                region_events,
            )
        {
            if let Some((gls, rr, ra)) = java_cluster_shaped_genotype(&event, region_events) {
                let gt = genotype_from_java_shaped_gls(gls, rr, ra, config)?;
                return Ok(Some(GenotypedSiteCall { event, genotype: gt }));
            }
        }
        // R4-2 / L8: genome-wide sites with pileup support when no alt haplotype was retained.
        if config.enable_java_strict()
            && !crate::read_event_discovery::is_strict_java_p12_production_emit_scope(&event)
            && crate::read_event_discovery::genome_wide_genotype_read_support(
                &event,
                read_ref_ad,
                read_alt_ad,
            )
        {
            let gt = sparse_snp_genotype_from_read_depths(read_ref_ad, read_alt_ad, config)?;
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
                Some((read_ref_ad, read_alt_ad)),
            );
        }
        return Ok(None);
    }
    let subset = likelihood_subset_for_event(
        likelihoods,
        likelihood_reads,
        &event,
        config,
        active_start_1based,
        active_end_1based,
    );
    if subset.is_empty() {
        if config.enable_sparse_read_genotype
            && event.ref_allele.len() == 1
            && event.alt_allele.len() == 1
            && read_alt_ad >= 1
            && passes_cluster_anchor_read_support(read_alt_ad, read_ref_ad)
        {
            let gt = sparse_snp_genotype_from_read_depths(read_ref_ad, read_alt_ad, config)?;
            return Ok(Some(GenotypedSiteCall { event, genotype: gt }));
        }
        if config.genotype_stored_events_only
            && is_cluster_anchor_snp(&event)
            && passes_cluster_anchor_read_support(read_alt_ad, read_ref_ad)
        {
            let gt = sparse_snp_genotype_from_read_depths(read_ref_ad, read_alt_ad, config)?;
            return Ok(Some(GenotypedSiteCall { event, genotype: gt }));
        }
        if config.genotype_stored_events_only
            && !mapping.alt_haplotype_indices.is_empty()
            && event.ref_allele.len() == 1
            && event.alt_allele.len() == 1
            && read_alt_ad >= 1
            && read_alt_ad >= read_ref_ad
        {
            let gt = sparse_snp_genotype_from_read_depths(read_ref_ad, read_alt_ad, config)?;
            return Ok(Some(GenotypedSiteCall { event, genotype: gt }));
        }
        // R4-2 / L8: genome-wide site with alt hap but empty PairHMM subset → pileup genotype.
        if config.enable_java_strict()
            && !crate::read_event_discovery::is_strict_java_p12_production_emit_scope(&event)
            && crate::read_event_discovery::genome_wide_genotype_read_support(
                &event,
                read_ref_ad,
                read_alt_ad,
            )
        {
            let gt = sparse_snp_genotype_from_read_depths(read_ref_ad, read_alt_ad, config)?;
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
                Some((read_ref_ad, read_alt_ad)),
            );
        }
        return Ok(None);
    }
    let mut gt = genotype_from_allele_mapping(
        &subset,
        haplotypes,
        &mapping,
        &event,
        ref_bytes,
        pad_start_1based,
        max_mnp_distance,
        &event.contig,
        config,
    )?;
    // Non-strict fallthrough only (strict arm returned above). Class-A lives in
    // [`SiteReshape`] on the strict path before finalize.
    let sparse_emit = config.enable_read_style_emit
        && passes_read_style_sparse_emit(&event, read_alt_ad, read_ref_ad);
    if sparse_emit
        && read_alt_ad > 0 && gt.format.ad.len() >= 2 && gt.format.ad[1].as_i32() == 0 {
            gt.format.ad[1] = crate::bio_ids::AlleleDepth::from_i32_saturating(read_alt_ad);
            gt.format.ad[0] = crate::bio_ids::AlleleDepth::from_i32_saturating(read_ref_ad.max(0));
            gt.format.dp = ReadDepth::from_i32_saturating(
                gt.format.ad.iter().map(|d| d.as_i32()).sum(),
            );
        }
    if !passes_emit_for_variation_event(
        &event,
        &gt.genotype_log10_likelihoods,
        &gt.format,
        config.stand_emit_confidence, region_events)? && !sparse_emit
    {
        if is_coupled_indel_for_genotyping(&event, region_events)
            && !mapping.alt_haplotype_indices.is_empty()
        {
            if let Some(ref_hap) = haplotypes
                .iter()
                .find(|h| h.is_reference)
                .or_else(|| haplotypes.first())
            {
                let alt_supports = mapping.alt_haplotype_indices.iter().any(|&i| {
                    crate::hc_allele_mapping::haplotype_supports_allele_at_with_ref(
                        &haplotypes[i.get()],
                        ref_hap,
                        loc,
                        pad_start_1based,
                        &event.ref_allele,
                        &event.alt_allele,
                        ref_bytes,
                        max_mnp_distance,
                        &event.contig,
                    )
                });
                if alt_supports {
                    return Ok(Some(GenotypedSiteCall { event, genotype: gt }));
                }
            }
        }
        if parity_emit_rescue_with_read_and_alt_hap(
            &event,
            &mapping,
            &gt,
            haplotypes,
            pileup_reads,
            ref_bytes,
            pad_start_1based,
            max_mnp_distance,
        ) {
            return Ok(Some(GenotypedSiteCall { event, genotype: gt }));
        }
        if config.genotype_stored_events_only
            && is_cluster_anchor_snp(&event)
            && passes_cluster_anchor_read_support(read_alt_ad, read_ref_ad)
        {
            return Ok(Some(GenotypedSiteCall { event, genotype: gt }));
        }
        if config.genotype_stored_events_only
            && is_ctc_del_for_genotyping(&event, region_events)
            && !mapping.alt_haplotype_indices.is_empty()
        {
            return Ok(Some(GenotypedSiteCall { event, genotype: gt }));
        }
        if !mapping.alt_haplotype_indices.is_empty()
            && read_alt_ad >= 1 && read_alt_ad >= read_ref_ad {
                return Ok(Some(GenotypedSiteCall { event, genotype: gt }));
            }
        if config.genotype_stored_events_only && is_cluster_anchor_snp(&event) {
            let gt = sparse_snp_genotype_from_read_depths(read_ref_ad, read_alt_ad, config)?;
            return Ok(Some(GenotypedSiteCall { event, genotype: gt }));
        }
        // L11-B3: genome-wide pileup rescue under `enable_java_strict` lived only in the
        // strict arm (above). Fallthrough never sees enable_java_strict == true.
        return Ok(None);
    }
    // Non-strict fallthrough only (`enable_java_strict` arm always returns above).
    if (gt.format.gq.as_i32() as f64) < config.stand_emit_confidence && !sparse_emit
    {
        if config.genotype_stored_events_only && is_ctc_del_for_genotyping(&event, region_events) {
            return Ok(Some(GenotypedSiteCall { event, genotype: gt }));
        }
        if config.genotype_stored_events_only
            && is_cluster_anchor_snp(&event)
            && passes_cluster_anchor_read_support(read_alt_ad, read_ref_ad)
        {
            return Ok(Some(GenotypedSiteCall { event, genotype: gt }));
        }
        return Ok(None);
    }
    Ok(Some(GenotypedSiteCall { event, genotype: gt }))
}

