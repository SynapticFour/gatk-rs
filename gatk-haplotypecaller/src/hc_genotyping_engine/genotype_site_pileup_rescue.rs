/// L14-B: empty-mapper pileup rescue (genome-wide indel/SNP) before score.
pub(crate) struct SitePileupRescue;

impl SitePileupRescue {
    pub(crate) fn try_empty_mapper(
        event: VariationEvent,
        mapping: &AlleleHaplotypeMapping,
        likelihood_reads: &[Record],
        pileup_reads: &[Record],
        haplotypes: &[Haplotype],
        pad_start_1based: u64,
        ref_bytes: &[u8],
        config: &HcGenotypingConfig,
        read_ref_ad: i32,
        read_alt_ad: i32,
    ) -> GatkResult<Option<GenotypedSiteCall>> {
        // R4-2: genome-wide indels with CIGAR pileup support — genotype from reads even
        // when the allele mapper has no alt haplotype (assembly often misses dense indels).
        if mapping.alt_haplotype_indices.is_empty()
            && event.is_indel()
            && !crate::read_event_discovery::is_strict_java_p12_production_emit_scope(&event)
            && crate::read_event_discovery::genome_wide_genotype_read_support(
                &event,
                read_ref_ad,
                read_alt_ad,
            )
        {
            if crate::parity_harness::env_flag_set("GATK_RS_INDEL_HAP_TRACE") {
                eprintln!(
                    "L9-indel-hap empty_mapper_rescue {}:{}>{} AD={},{} n_haps={}",
                    event.start_1based.get(),
                    event.ref_allele,
                    event.alt_allele,
                    read_ref_ad,
                    read_alt_ad,
                    haplotypes.len()
                );
            }
            let (shape_ref, shape_alt) =
                long_insertion_pileup_shape_ad(&event, read_ref_ad, read_alt_ad);
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
        // L9: same R4-2 pileup rescue for genome-wide SNPs with empty mapper (indels already
        // rescued above). Without this, SNPs near assembly-retained indels die at Ok(None).
        if mapping.alt_haplotype_indices.is_empty()
            && event.is_snp()
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
        Ok(None)
    }
}
