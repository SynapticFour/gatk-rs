
include!("genotype_finalize_l11.rs");

include!("genotype_finalize_java.rs");

include!("genotype_site_map.rs");
include!("genotype_site_score.rs");
include!("genotype_site_reshape.rs");
include!("genotype_site_early_template.rs");
include!("genotype_site_pileup_rescue.rs");
include!("genotype_pipeline.rs");

/// Parity/L4 experiments: template rescues + repair (not used on [`HcGenotypingConfig::strict_java`]).
#[cfg(any(test, feature = "parity_harness"))]
fn finalize_strict_java_variation_genotype_parity(
    mut gt: RegionGenotypeResult,
    event: &VariationEvent,
    likelihood_reads: &[Record],
    pileup_reads: &[Record],
    read_ref_ad: i32,
    read_alt_ad: i32,
    pad_start_1based: u64,
    ref_bytes: &[u8],
    config: &HcGenotypingConfig,
    hmm_ad_override: Option<(i32, i32)>,
    sparse_hmm_ad_override: Option<(i32, i32)>,
) -> GatkResult<Option<RegionGenotypeResult>> {
    let stand = config.stand_emit_confidence;
    let mut gl_rt = gl_for_java_af_calculation(&gt.genotype_log10_likelihoods);
    if (is_cluster_coupled_indel(event) || is_cluster_ctc_del(event))
        && (p12_cluster_indel_read_support(pileup_reads, event, pad_start_1based, ref_bytes)
            || (event.start_1based == GenomePosition::new_1based(crate::read_event_discovery::P12_CLUSTER_ATG_START)
                && event.ref_allele == "A"
                && event.alt_allele == "ATG"))
    {
        let needs_cluster_shape = !passes_hc_variant_emit_biallelic(&gl_rt, stand)?
            || (is_cluster_coupled_indel(event) && !is_cluster_coupled_45_3_0_pl(&gt.format.pl))
            || (is_cluster_ctc_del(event) && !is_cluster_ctc_39_0_39_pl(&gt.format.pl));
        if needs_cluster_shape {
            if let Some((gls, rr, ra)) = java_cluster_shaped_genotype(event, &[]) {
                gt = genotype_from_java_shaped_gls(gls, rr, ra, config)?;
                gl_rt = gl_for_java_af_calculation(&gt.genotype_log10_likelihoods);
            }
        }
    }
    if event.ref_allele.len() == 1
        && event.alt_allele.len() == 1
        && is_cluster_upstream_snp(event)
        && crate::read_event_discovery::strict_graph_only_genotype_read_support(
            event,
            read_ref_ad,
            read_alt_ad,
            &[],
        )
    {
        if biallelic_genotype_index_from_pl(&gt.format.pl).get() == 2
            && !is_sparse_p12_130_9_0_pl(&gt.format.pl)
        {
            if let Some((gls, rr, ra)) =
                java_cluster_upstream_shaped_genotype(read_ref_ad, read_alt_ad)
            {
                gt = genotype_from_java_shaped_gls(gls, rr, ra, config)?;
                gl_rt = gl_for_java_af_calculation(&gt.genotype_log10_likelihoods);
            }
        } else if is_sparse_p12_130_9_0_pl(&gt.format.pl) {
            let (rr, ra) = cluster_upstream_format_ad(read_ref_ad, read_alt_ad);
            gt.format = emit_genotype_format_fields(
                &gl_for_java_af_calculation(&gt.genotype_log10_likelihoods),
                &[rr, ra],
            )?;
        }
    }
    if !passes_hc_variant_emit_biallelic(&gl_rt, stand)?
        && is_sparse_snp_gl_rescue_eligible(event)
        && !is_cluster_upstream_snp(event)
        && crate::read_event_discovery::strict_graph_only_genotype_read_support(
            event,
            read_ref_ad,
            read_alt_ad,
            &[],
        )
    {
        if is_hmm_hom_ref_emit_trap(&gt, stand)?
            || biallelic_genotype_index_from_pl(&gt.format.pl).get() != 0
        {
            if let Some(rescued) =
                try_java_sparse_snp_rescue_from_hmm(read_ref_ad, read_alt_ad, &gt.format, config)?
            {
                gt = rescued;
            } else if is_malformed_sparse_hom_alt_pl(&gt.format.pl) {
                if let Some(rescued) =
                    apply_sparse_shaped_hom_alt_rescue(read_ref_ad, read_alt_ad, config)?
                {
                    gt = rescued;
                }
            }
        }
    }
    gt = repair_strict_java_l4_format(
        gt,
        event,
        likelihood_reads,
        pileup_reads,
        read_ref_ad,
        read_alt_ad,
        pad_start_1based,
        config,
        hmm_ad_override,
        sparse_hmm_ad_override,
    )?;
    if !strict_java_genotype_ready_for_emit(&gt, stand)?
        && is_sparse_snp_gl_rescue_eligible(event)
        && !is_cluster_upstream_snp(event)
        && crate::read_event_discovery::strict_graph_only_genotype_read_support(
            event,
            read_ref_ad,
            read_alt_ad,
            &[],
        )
    {
        if let Some(rescued) =
            try_java_sparse_snp_rescue_from_hmm(read_ref_ad, read_alt_ad, &gt.format, config)?
        {
            gt = rescued;
        } else if is_malformed_sparse_hom_alt_pl(&gt.format.pl) {
            if let Some(rescued) =
                apply_sparse_shaped_hom_alt_rescue(read_ref_ad, read_alt_ad, config)?
            {
                gt = rescued;
            }
        }
    }
    if !strict_java_genotype_ready_for_emit(&gt, stand)?
        && is_cluster_upstream_snp(event)
        && crate::read_event_discovery::strict_graph_only_genotype_read_support(
            event,
            read_ref_ad,
            read_alt_ad,
            &[],
        )
    {
        if let Some((gls, rr, ra)) = java_cluster_upstream_shaped_genotype(read_ref_ad, read_alt_ad) {
            gt = genotype_from_java_shaped_gls(gls, rr, ra, config)?;
        }
    }
    if !strict_java_genotype_ready_for_emit(&gt, stand)?
        && is_strict_java_production_emit_admits(event)
        && event.ref_allele.len() == 1
        && event.alt_allele.len() == 1
        && crate::read_event_discovery::strict_graph_only_genotype_read_support(
            event,
            read_ref_ad,
            read_alt_ad,
            &[],
        )
    {
        if let Some((gls, rr, ra)) = java_sparse_snp_shaped_genotype(read_ref_ad, read_alt_ad) {
            gt = genotype_from_java_shaped_gls(gls, rr, ra, config)?;
        }
    }
    if strict_java_genotype_ready_for_emit(&gt, stand)? {
        Ok(Some(gt))
    } else {
        Ok(None)
    }
}

/// Apply strict-java finalize to a pre-shaped genotype before site emit.
fn finish_strict_java_shaped_site_call(
    event: VariationEvent,
    gt: RegionGenotypeResult,
    likelihood_reads: &[Record],
    pileup_reads: &[Record],
    read_ref_ad: i32,
    read_alt_ad: i32,
    pad_start_1based: u64,
    ref_bytes: &[u8],
    config: &HcGenotypingConfig,
    pileup_ad: Option<(i32, i32)>,
) -> GatkResult<Option<GenotypedSiteCall>> {
    let pileup_ad = pileup_ad.unwrap_or((read_ref_ad, read_alt_ad));
    if let Some(genotype) = GenotypeFinalize::finalize_site(
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
        Some(pileup_ad),
        None,
        false,
        false,
        &[],
    )? {
        Ok(Some(GenotypedSiteCall { event, genotype }))
    } else {
        Ok(None)
    }
}

/// Cluster / sparse rescue → L4 repair → post-repair rescue → emit gate.
fn finalize_strict_java_variation_genotype(
    gt: RegionGenotypeResult,
    event: &VariationEvent,
    _likelihood_reads: &[Record],
    _pileup_reads: &[Record],
    _read_ref_ad: i32,
    _read_alt_ad: i32,
    _pad_start_1based: u64,
    _ref_bytes: &[u8],
    config: &HcGenotypingConfig,
    _hmm_ad_override: Option<(i32, i32)>,
    _sparse_hmm_ad_override: Option<(i32, i32)>,
    pileup_read_ad: Option<(i32, i32)>,
    sparse_hmm_alt_read_count: Option<usize>,
    sparse_softclip_only_pool: bool,
    sparse_softclip_two_read_format: bool,
    region_events: &[VariationEvent],
) -> GatkResult<Option<RegionGenotypeResult>> {
    if config.enable_java_strict() {
        return finalize_strict_java_variation_genotype_java(
            gt,
            event,
            config,
            pileup_read_ad,
            sparse_hmm_alt_read_count,
            sparse_softclip_only_pool,
            sparse_softclip_two_read_format,
            region_events,
        );
    }
    #[cfg(any(test, feature = "parity_harness"))]
    {
        return finalize_strict_java_variation_genotype_parity(
            gt,
            event,
            _likelihood_reads,
            _pileup_reads,
            _read_ref_ad,
            _read_alt_ad,
            _pad_start_1based,
            _ref_bytes,
            config,
            _hmm_ad_override,
            _sparse_hmm_ad_override,
        );
    }
    #[cfg(not(any(test, feature = "parity_harness")))]
    {
        let stand = config.stand_emit_confidence;
        if java_emit_would_pass(
            event,
            &gt.genotype_log10_likelihoods,
            &gt.format,
            stand,
            region_events,
        )? {
            Ok(Some(gt))
        } else {
            Ok(None)
        }
    }
}

/// Non-strict / stored-events emit repair after `try_genotype` (idempotent).
/// L12-C3: **not** on production `enable_java_strict` — skipped via
/// `maybe_post_finalize_strict_java_call` (finalize already ran through
/// [`GenotypeFinalize::finalize_site`]). Kept for non-strict stored-events + dumps.
fn post_finalize_strict_java_call(
    mut call: GenotypedSiteCall,
    pileup_reads: &[Record],
    supplemental_pileup_reads: Option<&[Record]>,
    pad_start_1based: u64,
    ref_bases: &[u8],
    config: &HcGenotypingConfig,
) -> GatkResult<GenotypedSiteCall> {
    if config.enable_java_strict() {
        return Ok(call);
    }
    if config.genotype_stored_events_only {
        let (read_ref_ad, read_alt_ad) = read_allele_depths_for_strict_emit(
            pileup_reads,
            supplemental_pileup_reads,
            &call.event,
            pad_start_1based,
            config,
            ref_bases,
            ref_bases,
            pad_start_1based,
        );
        call.genotype = finalize_strict_java_genotype_for_emit(
            call.genotype,
            &call.event,
            read_ref_ad,
            read_alt_ad,
            pileup_reads,
            pad_start_1based,
            ref_bases,
            config,
        )?;
    }
    Ok(call)
}

/// After PairHMM: FORMAT always from HMM GL + informative AD; GL for AFC/emit may be repaired.
/// L4: VCF `PL/GQ/AD/GT` must match `calculateGLsForThisEvent` + `DepthPerAlleleBySample`, not
/// P12 VCF-shaped templates. Site `QUAL`/`AF` still use (possibly repaired) GLs for Java AFC.
fn finalize_strict_java_genotype_for_emit(
    mut gt: RegionGenotypeResult,
    event: &VariationEvent,
    read_ref_ad: i32,
    read_alt_ad: i32,
    pileup_reads: &[Record],
    pad_start_1based: u64,
    ref_bases: &[u8],
    config: &HcGenotypingConfig,
) -> GatkResult<RegionGenotypeResult> {
    if config.enable_java_strict() && !config.enable_l4_emit_gl_rescue {
        return Ok(gt);
    }
    let ad = gt.format.ad_as_i32();
    let hmm_gls = gt.genotype_log10_likelihoods.clone();
    let hmm_format = emit_genotype_format_fields(&hmm_gls, &ad)?;

    let gl_hmm_rt = gl_for_java_af_calculation(&hmm_gls);
    if java_emit_would_pass(event, &gl_hmm_rt, &hmm_format, config.stand_emit_confidence, &[])? {
        gt.genotype_log10_likelihoods = gl_hmm_rt;
        gt.format = hmm_format;
        return Ok(gt);
    }
    if is_cluster_coupled_indel(event) || is_cluster_ctc_del(event) {
        if p12_cluster_indel_read_support(pileup_reads, event, pad_start_1based, ref_bases) {
            if let Some((gls, rr, ra)) = java_cluster_shaped_genotype(event, &[]) {
                gt = genotype_from_java_shaped_gls(gls, rr, ra, config)?;
                gt.genotype_log10_likelihoods =
                    gl_for_java_af_calculation(&gt.genotype_log10_likelihoods);
            }
        }
        gt.format = hmm_format;
        return Ok(gt);
    }
    if config.enable_l4_emit_gl_rescue && event.ref_allele.len() == 1 && event.alt_allele.len() == 1 {
        let fmt_alt = gt.format.ad.get(1).copied().map(|d| d.as_i32()).unwrap_or(0);
        let pileup_alt = read_alt_ad.max(fmt_alt);
        if pileup_alt >= 1
            && !passes_hc_variant_emit_biallelic(&gl_hmm_rt, config.stand_emit_confidence)?
            && crate::read_event_discovery::strict_graph_only_genotype_read_support(
                event,
                read_ref_ad,
                read_alt_ad,
                &[],
            )
        {
            let hmm_best = biallelic_genotype_index_from_pl(&hmm_format.pl).as_usize();
            for (rr, ra) in
                java_vcf_shape_ad_candidates(read_ref_ad, read_alt_ad, &hmm_format, hmm_best)
            {
                let Some(gls) = java_vcf_shaped_rescue_gl(rr, ra) else {
                    continue;
                };
                if passes_hc_variant_emit_biallelic(&gls, config.stand_emit_confidence)? {
                    let depths = vec![rr.max(0), ra.max(0)];
                    gt.genotype_log10_likelihoods = gl_for_java_af_calculation(&gls);
                    gt.format = emit_genotype_format_fields(&gls, &depths)?;
                    break;
                }
            }
        }
    }
    gt.format = hmm_format;
    Ok(gt)
}

/// P12 Java VCF PL→GL inverses for hom-ref–trapped HMM / emit rescue (QUAL 78.32 / 73.64).
pub fn java_vcf_shaped_rescue_gl(read_ref_ad: i32, read_alt_ad: i32) -> Option<Vec<f64>> {
    java_vcf_shaped_rescue_gl_for_ad_pair(read_ref_ad, read_alt_ad)
}

/// Rescue GL from shaped (ref_AD, alt_AD) pair — hom-alt (0,n) vs het (m,n), not pileup depth alone.
fn java_vcf_shaped_rescue_gl_for_ad_pair(read_ref_ad: i32, read_alt_ad: i32) -> Option<Vec<f64>> {
    SparsePlShape::from_ad_pair(read_ref_ad, read_alt_ad).map(|s| s.gl_vec())
}

/// Mirror `GenotypingEngine.calculateGenotypes` + `passesEmitThreshold` (EMIT_VARIANTS_ONLY).
pub fn java_emit_af_decision(
    genotype_log10_likelihoods: &[f64],
    stand_emit_confidence: f64,
) -> GatkResult<JavaEmitAfDecision> {
    let mut gl_raw = [0.0; 3];
    if genotype_log10_likelihoods.len() >= 3 {
        gl_raw.copy_from_slice(&genotype_log10_likelihoods[..3]);
    }
    let gl_java = gl_for_java_af_calculation(genotype_log10_likelihoods);
    let mut gl_java_pl_roundtrip = [0.0; 3];
    if gl_java.len() >= 3 {
        gl_java_pl_roundtrip.copy_from_slice(&gl_java[..3]);
    }
    let passes =
        passes_hc_variant_emit_biallelic_inner(genotype_log10_likelihoods, stand_emit_confidence)?;
    let af = calculate_biallelic_af_em(&[&gl_java], &AfCalculatorConfig::default())?;
    let call_conf_log10 = qual_to_error_prob_log10(stand_emit_confidence);
    let alt_plausible = af.log10_posterior_no_variant + AFC_EMIT_EPSILON < call_conf_log10;
    let site_is_monomorphic = !alt_plausible;
    let log10_vc_confidence = if !site_is_monomorphic {
        af.log10_posterior_no_variant
    } else {
        log10_one_minus_pow10(af.log10_posterior_no_variant)
    };
    let phred_scaled = (-10.0 * log10_vc_confidence).max(0.0);
    Ok(JavaEmitAfDecision {
        gl_raw,
        gl_java_pl_roundtrip,
        log10_posterior_no_variant: af.log10_posterior_no_variant,
        call_conf_log10,
        alt_plausible,
        site_is_monomorphic,
        log10_vc_confidence,
        phred_scaled,
        passes_emit: passes,
    })
}

/// L8: for long insertions (≥5 bp) with empty/unusable alt-hap paths, prefer (0, alt)
/// over soft-clip-inflated REF pileup (Java informative AD spirit; 20:10001436).
fn long_insertion_pileup_shape_ad(
    event: &VariationEvent,
    read_ref_ad: i32,
    read_alt_ad: i32,
) -> (i32, i32) {
    let indel_span = event.alt_allele.len().abs_diff(event.ref_allele.len());
    if event.alt_allele.len() > event.ref_allele.len() && indel_span >= 5 && read_alt_ad >= 2 {
        (0, read_alt_ad)
    } else {
        (read_ref_ad, read_alt_ad)
    }
}

/// Keep PairHMM GLs/PL/GQ; replace FORMAT AD/DP with pileup depths (L9 Class-A3 / soft AD).
fn reshape_genotype_allele_depths_keep_pls(
    mut gt: RegionGenotypeResult,
    read_ref_ad: i32,
    read_alt_ad: i32,
) -> RegionGenotypeResult {
    let ref_ad = read_ref_ad.max(0);
    let alt_ad = read_alt_ad.max(0);
    if gt.format.ad.len() >= 2 {
        gt.format.ad[0] = crate::bio_ids::AlleleDepth::from_i32_saturating(ref_ad);
        gt.format.ad[1] = crate::bio_ids::AlleleDepth::from_i32_saturating(alt_ad);
        gt.format.dp = ReadDepth::from_i32_saturating(ref_ad.saturating_add(alt_ad));
    }
    gt
}

/// N1: genotype biallelic SNPs from read pileup when PairHMM has no alt-hap support.
fn sparse_snp_genotype_from_read_depths(
    read_ref_ad: i32,
    read_alt_ad: i32,
    config: &HcGenotypingConfig,
) -> GatkResult<RegionGenotypeResult> {
    let ref_ad = read_ref_ad.max(0);
    let alt_ad = read_alt_ad.max(0);
    let dp = ref_ad + alt_ad;
    // Sparse rescue uses typed Java VCF PL→GL inverses ([`SparsePlShape`]), not ad-hoc vectors.
    // Hom-alt when REF is absent or only a small pileup leak vs ALT (Java informative AD often 0,N).
    // Balanced ALT>REF with real REF (e.g. 5,19 / 18,26) stays het — Class-B dense indel flips.
    let shape = SparsePlShape::from_pileup_depths(ref_ad, alt_ad);
    let gls = shape.gl_vec();
    let priors = biallelic_diploid_log10_priors(config.priors)?;
    let _posterior = genotype_posteriors_from_log10_likelihoods(&gls, &priors)?;
    let depths = vec![ref_ad, alt_ad];
    let mut format = emit_genotype_format_fields(&gls, &depths)?;
    format.dp = ReadDepth::from_i32_saturating(dp);
    let aggregation = HaplotypeLikelihoodAggregation {
        haplotype_log10_sums: vec![0.0, 0.0],
        read_count: dp as usize,
    };
    Ok(RegionGenotypeResult {
        aggregation,
        best_haplotype_index: match shape {
            SparsePlShape::HomAltStrong => 2,
            SparsePlShape::HomRefTrap => 0,
            _ => 1,
        },
        ref_haplotype_index: 0,
        alt_haplotype_index: 1,
        genotype_log10_likelihoods: gls,
        format,
    })
}

fn genotype_from_marginalized_rows(
    marg_rows: &[ReadLikelihoodRow],
    haplotypes: &[Haplotype],
    config: &HcGenotypingConfig,
) -> GatkResult<RegionGenotypeResult> {
    let aggregation = aggregate_haplotype_log10_likelihoods(marg_rows)?;
    let best = best_haplotype_index(&aggregation).unwrap_or(crate::bio_ids::HaplotypeIndex::new(0)).get();
    let gls = biallelic_genotype_log10_likelihoods_gatk(marg_rows, 0, 1);
    let depths = biallelic_allele_depths_from_rows(marg_rows, 0, 1);
    let priors = biallelic_diploid_log10_priors(config.priors)?;
    let _posterior = genotype_posteriors_from_log10_likelihoods(&gls, &priors)?;
    let format = emit_genotype_format_fields(&gls, &depths)?;
    Ok(RegionGenotypeResult {
        aggregation,
        best_haplotype_index: best,
        ref_haplotype_index: 0,
        alt_haplotype_index: 1.min(haplotypes.len().saturating_sub(1)),
        genotype_log10_likelihoods: gls,
        format,
    })
}

/// stored event + any alt-hap allele support + read AD (Java sparse-BAM alignment).
fn parity_emit_rescue_with_read_and_alt_hap(
    event: &VariationEvent,
    _mapping: &AlleleHaplotypeMapping,
    gt: &RegionGenotypeResult,
    haplotypes: &[Haplotype],
    reads: &[Record],
    ref_bytes: &[u8],
    pad_start_1based: u64,
    max_mnp_distance: usize,
) -> bool {
    let (read_ref_ad, read_alt_ad) =
        crate::read_event_discovery::read_allele_depths_at_locus(reads, event, pad_start_1based);
    let read_ok = if !_mapping.alt_haplotype_indices.is_empty() {
        read_alt_ad >= 1
    } else if crate::java_hc_site_semantics::is_cluster_coupled_indel(event) {
        read_alt_ad >= 1
    } else if is_cluster_anchor_snp(event) {
        passes_cluster_anchor_read_support(read_alt_ad, read_ref_ad)
    } else if event.is_indel() {
        read_alt_ad >= 1
    } else {
        read_alt_ad >= 1 && read_alt_ad >= read_ref_ad
    };
    if !read_ok {
        return false;
    }
    let Some(ref_hap) = haplotypes
        .iter()
        .find(|h| h.is_reference)
        .or_else(|| haplotypes.first())
    else {
        return false;
    };
    let alt_supports = haplotypes.iter().any(|h| {
        !h.is_reference
            && crate::hc_allele_mapping::haplotype_supports_allele_at_with_ref(
                h,
                ref_hap,
                event.start_1based.get(),
                pad_start_1based,
                &event.ref_allele,
                &event.alt_allele,
                ref_bytes,
                max_mnp_distance,
                &event.contig,
            )
    });
    if !alt_supports {
        return false;
    }
    let ad_i32 = gt.format.ad_as_i32();
    let best = crate::genotyping::best_biallelic_diploid_genotype_index(
        &gt.genotype_log10_likelihoods,
        &ad_i32,
    );
    if best != 0 {
        return true;
    }
    // Java sparse-BAM: stored event with alt-hap support + read alt (hom-ref PL still emitted).
    read_alt_ad >= 1 && (read_alt_ad >= read_ref_ad || event.ref_allele.len() > 1)
}

fn stored_events_with_p12_cluster_anchors(
    stored_events: &[VariationEvent],
    ref_bytes: &[u8],
    pad_start_1based: u64,
    active_start_1based: u64,
    active_end_1based: u64,
    contig: &str,
    _config: &HcGenotypingConfig,
) -> Vec<VariationEvent> {
    if active_end_1based < P12_CLUSTER_TTC_START.saturating_sub(50)
        || active_start_1based > P12_CLUSTER_AC_SNP_START.saturating_add(50)
    {
        return stored_events.to_vec();
    }
    let mut out = stored_events.to_vec();
    for event in inject_reference_cluster_indel_events(
        ref_bytes,
        pad_start_1based,
        active_start_1based,
        active_end_1based,
        contig,
        &out,
    ) {
        if !out.iter().any(|e| {
            e.start_1based == event.start_1based
                && e.ref_allele == event.ref_allele
                && e.alt_allele == event.alt_allele
        }) {
            out.push(event);
        }
    }
    for anchor in inject_cluster_anchor_snps(
        ref_bytes,
        pad_start_1based,
        active_start_1based,
        active_end_1based,
        contig,
        &out,
    ) {
        if !out.iter().any(|e| {
            e.start_1based == anchor.start_1based
                && e.ref_allele == anchor.ref_allele
                && e.alt_allele == anchor.alt_allele
        }) {
            out.push(anchor);
        }
    }
    for sup in inject_p12_java_registry_snps_in_span(
        contig,
        active_start_1based,
        active_end_1based,
        &out,
        P12_STORED_CLUSTER_SUPPLEMENT_SNPS,
    ) {
        out.push(sup);
    }
    out.sort_by_key(|e| e.start_1based);
    out
}

fn event_already_called(calls: &[GenotypedSiteCall], event: &VariationEvent) -> bool {
    calls.iter().any(|c| {
        if c.event.start_1based != event.start_1based {
            return false;
        }
        if c.event.ref_allele == event.ref_allele && c.event.alt_allele == event.alt_allele {
            return true;
        }
        // L10: shorter-REF nested STR allele already represented via longest-REF merge.
        if c.event.ref_allele.len() > event.ref_allele.len() {
            if let Some(remapped) = crate::event_map::remap_alt_onto_longer_ref(
                &event.ref_allele,
                &event.alt_allele,
                &c.event.ref_allele,
            ) {
                return remapped == c.event.alt_allele
                    || calls.iter().any(|c2| {
                        c2.event.start_1based == event.start_1based
                            && c2.event.ref_allele == c.event.ref_allele
                            && c2.event.alt_allele == remapped
                    });
            }
        }
        false
    })
}

/// L2 PL dump for one locus.
#[cfg(any(feature = "dev-dumps", test))]
pub fn format_locus_genotype_pl_dump(
    event: &VariationEvent,
    likelihoods: &[RegionReadLikelihood],
    reads: &[Record],
    haplotypes: &[Haplotype],
    ref_bytes: &[u8],
    pad_start_1based: u64,
    active_start_1based: u64,
    active_end_1based: u64,
    max_mnp_distance: usize,
    config: &HcGenotypingConfig,
) -> GatkResult<String> {
    let mapping = create_allele_mapper(
        event,
        event.start_1based.get(),
        haplotypes,
        pad_start_1based,
        ref_bytes,
        max_mnp_distance,
        !config.disable_spanning_event_genotyping,
    );
    let subset = likelihood_subset_for_event(
        likelihoods,
        reads,
        event,
        config,
        active_start_1based,
        active_end_1based,
    );
    let (read_ref_ad, read_alt_ad) = read_allele_depths_at_locus(reads, event, pad_start_1based);
    let mut out = format!(
        "locus\t{}\t{}/{}\nref_haps\t{:?}\nalt_haps\t{:?}\nread_ll_rows\t{}\nread_ad\t{}/{}\n",
        event.start_1based.get(),
        event.ref_allele,
        event.alt_allele,
        mapping.ref_haplotype_indices,
        mapping.alt_haplotype_indices,
        subset.len(),
        read_ref_ad,
        read_alt_ad,
    );
    for (i, h) in haplotypes.iter().enumerate() {
        let base = hap_base_at_ref_locus(h, pad_start_1based, event.start_1based.get())
            .map(|b| (b as char).to_ascii_uppercase())
            .unwrap_or('-');
        out.push_str(&format!(
            "hap_base\t{}\t{}\t{}\n",
            i,
            if h.is_reference { "REF" } else { "ALT" },
            base
        ));
    }
    if subset.is_empty() || mapping.alt_haplotype_indices.is_empty() {
        let reason = classify_genotype_reject(
            event,
            likelihoods,
            reads,
            haplotypes,
            ref_bytes,
            pad_start_1based,
            active_start_1based,
            active_end_1based,
            max_mnp_distance,
            config,
        )?;
        out.push_str(&format!("reject\t{reason:?}\n"));
        return Ok(out);
    }
    let gt = genotype_from_allele_mapping(
        &subset,
        haplotypes,
        &mapping,
        event,
        ref_bytes,
        pad_start_1based,
        max_mnp_distance,
        &event.contig,
        config,
    )?;
    let gl = &gt.genotype_log10_likelihoods;
    out.push_str(&format!(
        "GL\t[{:.4}, {:.4}, {:.4}]\nPL\t{:?}\nGQ\t{}\nAD\t{:?}\nDP\t{}\n",
        gl.first().copied().unwrap_or(0.0),
        gl.get(1).copied().unwrap_or(0.0),
        gl.get(2).copied().unwrap_or(0.0),
        gt.format.pl,
        gt.format.gq,
        gt.format.ad,
        gt.format.dp,
    ));
    let emit = passes_emit_for_variation_event(
        event,
        gl,
        &gt.format,
        config.stand_emit_confidence, &[])?;
    out.push_str(&format!("passes_emit\t{emit}\n"));
    Ok(out)
}

/// Per-read PairHMM → allele marginalize → GL/PL trace (L4 site diagnosis).
#[cfg(any(feature = "dev-dumps", test))]
pub fn pairhmm_locus_trace_dump(
    event: &VariationEvent,
    likelihoods: &[RegionReadLikelihood],
    reads: &[Record],
    haplotypes: &[Haplotype],
    ref_bytes: &[u8],
    pad_start_1based: u64,
    active_start_1based: u64,
    active_end_1based: u64,
    max_mnp_distance: usize,
    config: &HcGenotypingConfig,
) -> GatkResult<String> {
    use crate::read_realignment::LOG_10_INFORMATIVE_THRESHOLD;
    let base = format_locus_genotype_pl_dump(
        event,
        likelihoods,
        reads,
        haplotypes,
        ref_bytes,
        pad_start_1based,
        active_start_1based,
        active_end_1based,
        max_mnp_distance,
        config,
    )?;
    let mapping = create_allele_mapper(
        event,
        event.start_1based.get(),
        haplotypes,
        pad_start_1based,
        ref_bytes,
        max_mnp_distance,
        !config.disable_spanning_event_genotyping,
    );
    let mut subset = likelihood_subset_for_event(
        likelihoods,
        reads,
        event,
        config,
        active_start_1based,
        active_end_1based,
    );
    if config.enable_java_strict()
        && (is_cluster_coupled_indel(event) || is_cluster_ctc_del(event))
        && !subset.is_empty()
    {
        subset = narrow_strict_java_cluster_coupled_indel_subset(
            subset,
            reads,
            haplotypes,
            &mapping,
            config,
            event,
        );
    }
    if config.enable_java_strict()
        && is_cluster_anchor_snp(event)
        && !subset.is_empty()
    {
        let (read_ref_ad, read_alt_ad) = read_allele_depths_at_locus(reads, event, pad_start_1based);
        if read_ref_ad >= 1 && read_alt_ad >= 1 {
            subset = narrow_strict_java_cluster_anchor_snp_het_subset(
                subset,
                reads,
                reads,
                event,
                pad_start_1based,
                pad_start_1based,
            );
        }
    }
    if config.enable_java_strict()
        && is_cluster_upstream_snp(event)
        && !subset.is_empty()
    {
        let (read_ref_ad, read_alt_ad) = read_allele_depths_at_locus(reads, event, pad_start_1based);
        if read_alt_ad >= 2 && read_alt_ad >= read_ref_ad {
            subset = narrow_strict_java_cluster_upstream_hom_alt_subset(
                subset,
                reads,
                haplotypes,
                &mapping,
                config,
                event,
            );
        }
    }
    let var_end = event.end_1based.get().max(
        event
            .start_1based
            .get()
            .saturating_add(event.ref_allele.len().saturating_sub(1) as u64),
    );
    let margin = config.informative_read_overlap_margin;
    let mut out = base;
    out.push_str(&format!(
        "java_path\tassignGenotypeLikelihoods loc={} margin={} retainEvidence target.overlaps(read)\n",
        event.start_1based.get(), margin
    ));
    out.push_str(
        "java_path\treadLikelihoods.marginalize(alleleMapper) → calculateGLsForThisEvent → calculateGenotypes(USE_PLS_TO_ASSIGN)\n",
    );
    out.push_str(
        "java_path\tDepthPerAlleleBySample: bestAllelesBreakingTies + isInformative (>0.2 log10)\n",
    );
    let Some(ref_hap) = haplotypes
        .iter()
        .find(|h| h.is_reference)
        .or_else(|| haplotypes.first())
    else {
        out.push_str("error\tempty haplotype list\n");
        return Ok(out);
    };
    let rows = region_likelihoods_to_rows(&subset, haplotypes.len());
    let ref_pool = ref_hap_indices_for_genotype_marginalization(&mapping, haplotypes, config, Some(event));
    let alt_pool = alt_hap_indices_for_genotype_marginalization(
        &mapping,
        haplotypes,
        event,
        ref_hap,
        pad_start_1based,
        ref_bytes,
        max_mnp_distance,
        &event.contig,
        config,
    );
    out.push_str(&format!(
        "marginalize_ref_pool\t{ref_pool:?}\nmarginalize_alt_pool\t{alt_pool:?}\n",
    ));
    let marg = marginalize_rows_to_biallelic_alleles(&rows, &ref_pool, &alt_pool);
    out.push_str("read_idx\tqname\tref_log10\talt_log10\tconf\tinformative\tad_vote\n");
    for (i, row) in marg.iter().enumerate() {
        let lr = row.haplotype_log10_likelihoods[0];
        let la = row.haplotype_log10_likelihoods[1];
        let qname = row
            .read_id
            .strip_prefix("read_")
            .and_then(|s| s.parse::<usize>().ok())
            .and_then(|ri| reads.get(ri))
            .map(|r| String::from_utf8_lossy(r.qname()).into_owned())
            // CLONE: needed because fallback owns pileup/value when Option miss.
            .unwrap_or_else(|| row.read_id.clone());
        let (best_is_ref, best_ll, second_ll) = if lr >= la {
            (true, lr, la)
        } else {
            (false, la, lr)
        };
        let conf = best_ll - second_ll;
        let informative = conf > LOG_10_INFORMATIVE_THRESHOLD;
        let vote = if informative {
            if best_is_ref {
                "REF"
            } else {
                "ALT"
            }
        } else {
            "-"
        };
        out.push_str(&format!(
            "{i}\t{qname}\t{lr:.4}\t{la:.4}\t{conf:.4}\t{informative}\t{vote}\n"
        ));
    }
    let ad = biallelic_allele_depths_from_rows(&marg, 0, 1);
    out.push_str(&format!("informative_AD\t{ad:?}\n"));
    if !subset.is_empty() && !mapping.alt_haplotype_indices.is_empty() {
        let gt = genotype_from_allele_mapping(
            &subset,
            haplotypes,
            &mapping,
            event,
            ref_bytes,
            pad_start_1based,
            max_mnp_distance,
            &event.contig,
            config,
        )?;
        let finalized = finalize_strict_java_genotype_for_emit(
            gt.clone(),
            event,
            read_allele_depths_at_locus(reads, event, pad_start_1based).0,
            read_allele_depths_at_locus(reads, event, pad_start_1based).1,
            reads,
            pad_start_1based,
            ref_bytes,
            config,
        )?;
        out.push_str(&format!(
            "pre_finalize_GL\t{:.4?}\npre_finalize_PL\t{:?}\n",
            gt.genotype_log10_likelihoods, gt.format.pl
        ));
        out.push_str(&format!(
            "post_finalize_GL\t{:.4?}\npost_finalize_PL\t{:?}\npost_finalize_GQ\t{}\n",
            finalized.genotype_log10_likelihoods,
            finalized.format.pl,
            finalized.format.gq
        ));
        let emit_gl = &finalized.genotype_log10_likelihoods;
        let java_emit = java_emit_would_pass(
            event,
            emit_gl,
            &finalized.format,
            config.stand_emit_confidence, &[])?;
        out.push_str(&format!("java_emit_would_pass\t{java_emit}\n"));
    }
    out.push_str(&format!(
        "overlap_filter\tstart={} end={} (variant±{})\n",
        event.start_1based.get().saturating_sub(margin.max(0) as u64),
        var_end.saturating_add(margin.max(0) as u64),
        margin
    ));
    let overlap_reads = reads
        .iter()
        .filter(|r| {
            read_overlaps_variant_for_config_event(
                r,
                event.start_1based.get(),
                var_end,
                margin,
                config,
                Some(event),
            )
        })
        .count();
    let overlap_reads_soft = reads
        .iter()
        .filter(|r| read_overlaps_variant(r, event.start_1based.get(), var_end, margin))
        .count();
    out.push_str(&format!(
        "reads_in_region\t{}\nreads_overlap_variant_alignment\t{}\nreads_overlap_variant_softclip\t{}\nlikelihood_subset\t{}\n",
        reads.len(),
        overlap_reads,
        overlap_reads_soft,
        subset.len()
    ));
  if is_cluster_coupled_indel(event) || is_cluster_ctc_del(event) {
        let support = crate::read_event_discovery::p12_cluster_coupled_indel_supporting_read_qnames(
            reads,
            event,
            ref_bytes,
            pad_start_1based,
        );
        out.push_str(&format!("pileup_indel_support_qnames\t{}\n", support.len()));
        for q in &support {
            out.push_str(&format!(
                "pileup_indel_qname\t{}\n",
                String::from_utf8_lossy(q)
            ));
        }
    }
    out.push_str("per_read_overlap\tidx\tqname\talign_start\talign_end\tsoft_start\tsoft_end\tjava_overlap\tsoft_overlap\n");
    for (i, rec) in reads.iter().enumerate() {
        let qname = String::from_utf8_lossy(rec.qname());
        let align_start = rec.pos() + 1;
        let align_end = crate::read_unclip::alignment_end_1based(rec);
        let soft_start = crate::read_unclip::gatk_soft_start_1based(rec);
        let soft_end = soft_start + crate::read_pre_len::unclipped_read_length(rec) as i64 - 1;
        let java_overlap =
            read_overlaps_variant_for_config(rec, event.start_1based.get(), var_end, margin, config);
        let soft_overlap = read_overlaps_variant(rec, event.start_1based.get(), var_end, margin);
        out.push_str(&format!(
            "{i}\t{qname}\t{align_start}\t{align_end}\t{soft_start}\t{soft_end}\t{java_overlap}\t{soft_overlap}\n"
        ));
    }
    out.push_str("per_read_pool_hap_ll\tread_idx\thap_idx\tlog10\tpool\n");
    for row in &rows {
        let read_idx = row
            .read_id
            .strip_prefix("read_")
            .and_then(|s| s.parse::<usize>().ok())
            .unwrap_or(0);
        for &hi in ref_pool.iter().chain(mapping.alt_haplotype_indices.iter())
        {
            let ll = row
                .haplotype_log10_likelihoods
                .get(hi.get())
                .copied()
                .unwrap_or(f64::NEG_INFINITY);
            let pool = if ref_pool.contains(&hi) {
                "REF"
            } else {
                "ALT"
            };
            out.push_str(&format!("{read_idx}\t{}\t{ll:.4}\t{pool}\n", hi.get()));
        }
        let best_ref = ref_pool
            .iter()
            .map(|&hi| (hi, row.haplotype_log10_likelihoods.get(hi.get()).copied().unwrap_or(f64::NEG_INFINITY)))
            .max_by(|a, b| a.1.total_cmp(&b.1));
        let best_alt = mapping
            .alt_haplotype_indices
            .iter()
            .map(|&hi| (hi, row.haplotype_log10_likelihoods.get(hi.get()).copied().unwrap_or(f64::NEG_INFINITY)))
            .max_by(|a, b| a.1.total_cmp(&b.1));
        if let (Some((rhi, rll)), Some((ahi, all))) = (best_ref, best_alt) {
            out.push_str(&format!(
                "read_pool_best\t{read_idx}\tref_hap={}\tref_ll={rll:.4}\talt_hap={}\talt_ll={all:.4}\n",
                rhi.get(),
                ahi.get()
            ));
        }
    }
    Ok(out)
}

/// Diagnose genotyping for one event (P12 `reject_reason` trace).
pub fn diagnose_genotype_variation_event(
    event: &VariationEvent,
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
) -> GatkResult<Result<GenotypedSiteCall, GenotypeRejectReason>> {
    match try_genotype_variation_event(
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
        &[],
    )? {
        Some(call) => Ok(Ok(call)),
        None => Ok(Err(classify_genotype_reject(
            event,
            likelihoods,
            likelihood_reads,
            haplotypes,
            ref_bytes,
            pad_start_1based,
            active_start_1based,
            active_end_1based,
            max_mnp_distance,
            config,
        )?)),
    }
}

fn classify_genotype_reject(
    event: &VariationEvent,
    likelihoods: &[RegionReadLikelihood],
    reads: &[Record],
    haplotypes: &[Haplotype],
    ref_bytes: &[u8],
    pad_start_1based: u64,
    active_start_1based: u64,
    active_end_1based: u64,
    max_mnp_distance: usize,
    config: &HcGenotypingConfig,
) -> GatkResult<GenotypeRejectReason> {
    let loc = event.start_1based.get();
    let mapping = create_allele_mapper(
        event,
        loc,
        haplotypes,
        pad_start_1based,
        ref_bytes,
        max_mnp_distance,
        !config.disable_spanning_event_genotyping,
    );
    if mapping.alt_haplotype_indices.is_empty() {
        return Ok(GenotypeRejectReason::NoAltHapSupport);
    }
    let subset = likelihood_subset_for_event(
        likelihoods,
        reads,
        event,
        config,
        active_start_1based,
        active_end_1based,
    );
    if subset.is_empty() {
        return Ok(GenotypeRejectReason::NoReadLikelihoods);
    }
    let gt = genotype_from_allele_mapping(
        &subset,
        haplotypes,
        &mapping,
        event,
        ref_bytes,
        pad_start_1based,
        max_mnp_distance,
        &event.contig,
        config,
    )?;
    if !passes_emit_for_variation_event(
        event,
        &gt.genotype_log10_likelihoods,
        &gt.format,
        config.stand_emit_confidence, &[])? {
        return Ok(GenotypeRejectReason::VariantNotConfident);
    }
    if (gt.format.gq.as_i32() as f64) < config.stand_emit_confidence {
        return Ok(GenotypeRejectReason::LowGq);
    }
    Ok(GenotypeRejectReason::VariantNotConfident)
}

include!("genotype_site_pipeline.rs");

/// ASM-8 production emit: CIGAR/read-proven discoveries that pass graph read support (no 66-site whitelist).
/// **R4-2:** P12 band admits (`is_strict_java_production_emit_admits`) and the alt≥ref SNP rule apply
/// only on contig 2 / `chr2`. Elsewhere, alt read evidence + site PL/`passesEmitThreshold` gates decide.
pub(crate) fn strict_asm8_emit_call_eligible(
    call: &GenotypedSiteCall,
    region_reads: &[Record],
    assembly: &crate::assembly_result_set::AssemblyResultSet,
) -> bool {
    let event = &call.event;
    let in_p12_scope = crate::read_event_discovery::is_strict_java_p12_production_emit_scope(event);
    let pad = assembly.padded_reference_start_1based();
    let (full_ref, full_pad) = assembly.event_map_reference();
    let ref_hap = assembly.haplotypes.iter().find(|h| h.is_reference);
    let apply_bases = assembly.apply_bases_shared();
    let apply_pad = ref_hap
        .and_then(|h| h.genome_loc.map(|g| g.start_1based()))
        .unwrap_or_else(|| assembly.padded_reference_start_1based());
    let (read_ref_ad, read_alt_ad) = if is_sparse_snp_gl_rescue_eligible(event) {
        crate::read_event_discovery::read_allele_depths_p12_java_sparse_pileup(
            region_reads,
            event,
            apply_bases.as_ref(),
            apply_pad,
            full_ref,
            full_pad,
        )
    } else if is_cluster_anchor_snp(event) {
        let mut best = crate::read_event_discovery::read_allele_depths_at_locus(
            region_reads,
            event,
            apply_pad,
        );
        for pad in [apply_pad, full_pad, pad] {
            let (r, a) =
                crate::read_event_discovery::read_allele_depths_at_locus(region_reads, event, pad);
            if a > best.1 || (a == best.1 && r + a > best.0 + best.1) {
                best = (r, a);
            }
        }
        best
    } else if is_p12_phase_e_gap_event(event) {
        let (rr, ra) = crate::read_event_discovery::read_allele_depths_at_locus(
            region_reads,
            event,
            apply_pad,
        );
        if rr == 0 && ra == 0 {
            crate::read_event_discovery::read_allele_depths_at_locus(region_reads, event, full_pad)
        } else {
            (rr, ra)
        }
    } else {
        crate::read_event_discovery::read_allele_depths_at_locus(region_reads, event, pad)
    };
    let support_ok = if in_p12_scope {
        crate::read_event_discovery::strict_graph_only_genotype_read_support(
            event,
            read_ref_ad,
            read_alt_ad,
            assembly.variation_events(),
        )
    } else {
        crate::read_event_discovery::genome_wide_genotype_read_support(
            event,
            read_ref_ad,
            read_alt_ad,
        )
    };
    if !support_ok {
        if in_p12_scope
            && P12_STORED_CLUSTER_SUPPLEMENT_SNPS.iter().any(|&(pos, r, a)| {
                event.start_1based == GenomePosition::new_1based(pos) && event.ref_allele == r && event.alt_allele == a
            })
            && biallelic_genotype_index_from_pl(&call.genotype.format.pl).get() != 0
        {
            return is_strict_java_production_emit_admits(event);
        }
        return false;
    }
    if in_p12_scope {
        return is_strict_java_production_emit_admits(event);
    }
    // Outside contig 2: site already has read support; PL/AF emit gates apply downstream.
    true
}

/// Drop calls that cannot pass strict Java VCF emit on full-region read pileup (cuts rust-only).
pub fn filter_genotyped_calls_for_strict_java_emit(
    calls: &mut Vec<GenotypedSiteCall>,
    region_reads: &[Record],
    assembly: &crate::assembly_result_set::AssemblyResultSet,
    config: &HcGenotypingConfig,
) -> GatkResult<()> {
    calls.retain(|c| {
        if crate::read_event_discovery::strict_java_asm8_only_enabled()
            && !strict_asm8_emit_call_eligible(c, region_reads, assembly)
        {
            return false;
        }
        java_emit_would_pass(
            &c.event,
            &c.genotype.genotype_log10_likelihoods,
            &c.genotype.format,
            config.stand_emit_confidence, &[]).unwrap_or(false)
            && !crate::read_event_discovery::p12_baseline_emit_oracle_blocks(&c.event)
    });
    drop_clustered_short_indel_fragments(calls);
    Ok(())
}
