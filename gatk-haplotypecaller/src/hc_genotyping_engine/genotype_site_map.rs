/// L13-B: allele-map stage before PairHMM score / [`SiteReshape`] / finalize.
/// Owns `create_allele_mapper` + L9 full-pad indel retry. Orchestration still decides
/// *when* to map; mapper construction lives here.
pub(crate) struct SiteMap;

impl SiteMap {
    /// Build the allele↔haplotype mapping for a variation event.
    /// When the trim-pad EventMap/linear slice misses indel alt haps (empty mapper),
    /// retries against the full padded reference (L9 SparsePlShape avoidance).
    ///
    /// `hap_events` must be built against `pad_start_1based` (trim pad); full-pad retry
    /// rebuilds EventMaps (coordinates differ).
    pub(crate) fn build_mapping(
        event: &VariationEvent,
        haplotypes: &[Haplotype],
        ref_bytes: &[u8],
        pad_start_1based: u64,
        full_reference_bases: &[u8],
        full_reference_pad_1based: u64,
        max_mnp_distance: usize,
        config: &HcGenotypingConfig,
        hap_events: Option<&crate::event_map::PerHaplotypeVariationEvents>,
    ) -> AlleleHaplotypeMapping {
        let profiling = crate::hc_profile::enabled();
        let t0 = profiling.then(std::time::Instant::now);
        let loc = event.start_1based.get();
        let mut mapping = create_allele_mapper_with_events(
            event,
            loc,
            haplotypes,
            pad_start_1based,
            ref_bytes,
            max_mnp_distance,
            !config.disable_spanning_event_genotyping,
            hap_events,
        );
        if mapping.alt_haplotype_indices.is_empty()
            && event.is_indel()
            && (full_reference_pad_1based != pad_start_1based
                || full_reference_bases.as_ptr() != ref_bytes.as_ptr())
        {
            let mapping_full = create_allele_mapper(
                event,
                loc,
                haplotypes,
                full_reference_pad_1based,
                full_reference_bases,
                max_mnp_distance,
                !config.disable_spanning_event_genotyping,
            );
            if !mapping_full.alt_haplotype_indices.is_empty() {
                if crate::parity_harness::env_flag_set("GATK_RS_INDEL_HAP_TRACE") {
                    eprintln!(
                        "L9-indel-hap full_pad_mapper {}:{}>{} alt_haps={}",
                        event.start_1based.get(),
                        event.ref_allele,
                        event.alt_allele,
                        mapping_full.alt_haplotype_indices.len()
                    );
                }
                mapping = mapping_full;
            }
        }
        if let Some(t0) = t0 {
            crate::hc_profile::note_allele_map_wall(t0.elapsed());
        }
        mapping
    }
}
