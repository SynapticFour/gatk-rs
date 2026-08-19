/// L13-B: PairHMM → biallelic GL / informative-AD score stage.
/// Owns marginalize + normalize + [`genotype_from_marginalized_rows`] (via the former
/// `genotype_from_allele_mapping` body). Called after [`SiteMap`] and before [`SiteReshape`].
pub(crate) struct SiteScore;

/// Typed ALT haplotype index pool for biallelic marginalize (L13-C4).
#[derive(Debug, Clone)]
pub(crate) struct AltHapSubset(pub Vec<HaplotypeIndex>);

impl AltHapSubset {
    pub(crate) fn as_slice(&self) -> &[HaplotypeIndex] {
        &self.0
    }
}

impl SiteScore {
    /// Score a site from allele↔haplotype mapping and region likelihood rows.
    pub(crate) fn from_allele_mapping(
        likelihoods: &[RegionReadLikelihood],
        haplotypes: &[Haplotype],
        mapping: &AlleleHaplotypeMapping,
        event: &VariationEvent,
        ref_bytes: &[u8],
        pad_start_1based: u64,
        max_mnp_distance: usize,
        contig: &str,
        config: &HcGenotypingConfig,
    ) -> GatkResult<RegionGenotypeResult> {
        let ref_hap = haplotypes
            .iter()
            .find(|h| h.is_reference)
            .or_else(|| haplotypes.first())
            .ok_or_else(|| {
                gatk_common::GatkError::algorithm(
                    "SiteScore::from_allele_mapping: haplotype list is empty",
                )
            })?;
        let profiling = crate::hc_profile::enabled();
        let t_marg = profiling.then(std::time::Instant::now);
        let out = with_region_likelihood_rows(likelihoods, haplotypes.len(), |rows| {
            let ref_pool = ref_hap_indices_for_genotype_marginalization(
                mapping,
                haplotypes,
                config,
                Some(event),
            );
            let alt_pool = AltHapSubset(alt_hap_indices_for_genotype_marginalization(
                mapping,
                haplotypes,
                event,
                ref_hap,
                pad_start_1based,
                ref_bytes,
                max_mnp_distance,
                contig,
                config,
            ));
            let mut marg =
                marginalize_rows_to_biallelic_alleles(rows, &ref_pool, alt_pool.as_slice());
            if config.enable_java_strict() {
                // Java AlleleLikelihoodMatrixMapper normalize applies at all strict sites (not only sparse).
                apply_java_marginal_normalize_gap(&mut marg);
            }
            if let Some(t0) = t_marg {
                crate::hc_profile::note_marginalize_wall(t0.elapsed());
            }
            let t_enum = profiling.then(std::time::Instant::now);
            let out = genotype_from_marginalized_rows(&marg, haplotypes, config);
            if let Some(t0) = t_enum {
                crate::hc_profile::note_genotype_enum_wall(t0.elapsed());
            }
            out
        });
        out
    }
}

/// Compatibility wrapper — production call sites may use [`SiteScore::from_allele_mapping`].
fn genotype_from_allele_mapping(
    likelihoods: &[RegionReadLikelihood],
    haplotypes: &[Haplotype],
    mapping: &AlleleHaplotypeMapping,
    event: &VariationEvent,
    ref_bytes: &[u8],
    pad_start_1based: u64,
    max_mnp_distance: usize,
    contig: &str,
    config: &HcGenotypingConfig,
) -> GatkResult<RegionGenotypeResult> {
    SiteScore::from_allele_mapping(
        likelihoods,
        haplotypes,
        mapping,
        event,
        ref_bytes,
        pad_start_1based,
        max_mnp_distance,
        contig,
        config,
    )
}
