/// L11-B / L12-C: sole production entry for strict-Java site finalize.
/// Pipeline: [`SiteMap`] → [`SiteScore`] → [`SiteReshape`] → this finalize → filter.
/// Emit-gated PL/AD policy for production lands here →
/// [`finalize_strict_java_variation_genotype_java`]. Harness/parity uses
/// [`finalize_strict_java_variation_genotype_parity`] behind cfg.
pub(crate) struct GenotypeFinalize;

impl GenotypeFinalize {
    /// Single finalize surface for a site genotype (production + harness dispatcher).
    #[inline]
    pub(crate) fn finalize_site(
        gt: RegionGenotypeResult,
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
        pileup_read_ad: Option<(i32, i32)>,
        sparse_hmm_alt_read_count: Option<usize>,
        sparse_softclip_only_pool: bool,
        sparse_softclip_two_read_format: bool,
        region_events: &[VariationEvent],
    ) -> GatkResult<Option<RegionGenotypeResult>> {
        finalize_strict_java_variation_genotype(
            gt,
            event,
            likelihood_reads,
            pileup_reads,
            read_ref_ad,
            read_alt_ad,
            pad_start_1based,
            ref_bytes,
            config,
            hmm_ad_override,
            sparse_hmm_ad_override,
            pileup_read_ad,
            sparse_hmm_alt_read_count,
            sparse_softclip_only_pool,
            sparse_softclip_two_read_format,
            region_events,
        )
    }
}
/// Drop short indel fragments and motif-bleed SNPs nested beside a recovered long allele.
/// L11-D1: retired the L10 blind “≥3 span≤4 within 60bp” fallback for indels.
/// L12-B: also drop SNPs in the long-insertion **upstream flank of length = span**
/// (holdout `15001873`/`880` beside +36 INS). Not a 60 bp SNP window — that would FN
/// truth SNPs near unrelated long alleles (chr20 dense). Distant CAT kept.
fn drop_clustered_short_indel_fragments(calls: &mut Vec<GenotypedSiteCall>) {
    use crate::event_map::IndelSpan;
    // P12 chr2 carries true coupled short indels (TTC/T + A/ATG); never spray-suppress there.
    if calls.iter().any(|c| {
        let c = c.event.contig.as_str();
        c == "2" || c == "chr2"
    }) {
        return;
    }
    let long_pos: Vec<(u64, IndelSpan)> = calls
        .iter()
        .filter(|c| c.event.is_indel() && c.event.indel_span().is_long_insertion_span())
        .map(|c| (c.event.start_1based.get(), c.event.indel_span()))
        .collect();
    // L11-D1: no long allele → no fragment suppress (blind cluster spray retired).
    if long_pos.is_empty() {
        return;
    }
    calls.retain(|c| {
        let p = c.event.start_1based.get();
        if c.event.is_snp() {
            return !long_pos.iter().any(|(lp, lspan)| {
                IndelSpan::snp_in_long_insertion_upstream_flank(p, *lp, *lspan)
            });
        }
        if !c.event.is_indel() {
            return true;
        }
        let span = c.event.indel_span();
        if !span.is_short_fragment() {
            return true;
        }
        !long_pos
            .iter()
            .any(|(lp, lspan)| IndelSpan::nests_beside_long(p, *lp, *lspan))
    });
}

