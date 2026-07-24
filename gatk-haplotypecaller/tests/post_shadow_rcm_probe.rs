//! Post-cluster shadow RCM probe (`92305755–92305823`).
//! Run: `P12_REFERENCE=…/hs37d5.simple.fa cargo test -p gatk-haplotypecaller post_shadow_rcm_parity_probe --release -- --ignored --nocapture`

use gatk_core::reference::{parse_intervals_cli_string, ReferenceWindowCache, SequenceDictionary};
use gatk_haplotypecaller::engine::{CallRegionArgs, HaplotypeCallerEngine};
use gatk_haplotypecaller::locus_iterator::LocusPileupState;
use gatk_haplotypecaller::read_model::ReadFilterParams;
use gatk_haplotypecaller::ref_confidence::{
    capped_genotype_likelihoods_by_hom_ref, reference_confidence_loci_for_active_region,
    reference_gq_from_log10_gl, ClusterRcmEvidenceMode, ReferenceConfidenceConfig,
};
use gatk_haplotypecaller::walker_apply::{call_disposition, AssemblyRegionCallDisposition};
use gatk_haplotypecaller::walker_traversal::{
    flatten_assembly_regions, traverse_assembly_region_walker, WalkerTraversalConfig,
};
use rust_htslib::bam::Read;
use std::path::PathBuf;

fn p12_paths() -> Option<(PathBuf, PathBuf)> {
    let ref_fasta = std::env::var("P12_REFERENCE").ok().map(PathBuf::from)?;
    let bam = PathBuf::from("parity/realworld/na12878_20k_b37/NA12878_20k.b37.bam");
    if ref_fasta.is_file() && bam.is_file() {
        Some((ref_fasta, bam))
    } else {
        None
    }
}

#[test]
#[ignore = "P12 BAM: post-shadow GQ stripe parity"]
fn post_shadow_rcm_parity_probe() {
    let Some((ref_fasta, bam)) = p12_paths() else {
        eprintln!("skip: set P12_REFERENCE");
        return;
    };

    let dict = SequenceDictionary::from_fasta_path(&ref_fasta).expect("dict");
    let specs = parse_intervals_cli_string(&dict, "2:92305500-92305830").expect("interval");
    let walk = traverse_assembly_region_walker(
        &dict,
        &specs,
        &ref_fasta,
        &bam,
        &ReadFilterParams::gatk_standard_hc(),
        &WalkerTraversalConfig::gatk_haplotype_caller_production(100),
    )
    .expect("walk");
    let regions = flatten_assembly_regions(&walk);
    let region = regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= 92305755
                && r.end.get() >= 92305823
        })
        .expect("active region spanning post-shadow band");

    let header = rust_htslib::bam::Reader::from_path(&bam)
        .expect("bam")
        .header()
        .clone();
    let outcome = HaplotypeCallerEngine::call_region(
        region,
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("call")
    .expect("outcome");

    let emitted: Vec<u64> = outcome
        .genotyped_calls
        .iter()
        .map(|c| c.event.start_1based.get())
        .collect();
    let first_variant = emitted.iter().copied().min();
    let filters = ReadFilterParams::genotyping_evidence_rcm_pileup();
    let config = ReferenceConfidenceConfig::default();
    let mut ref_cache = ReferenceWindowCache::new(ref_fasta.clone(), 4);
    let loci = reference_confidence_loci_for_active_region(
        region,
        &outcome.genotyping_reads,
        first_variant,
        &emitted,
        &header,
        &config,
        &filters,
        &mut ref_cache,
        &dict,
        ClusterRcmEvidenceMode::ReconcileGenotypingFirst,
    )
    .expect("loci");

    let gt_filters = ReadFilterParams::genotyping_evidence_rcm_pileup();
    let mut gt_after = LocusPileupState::from_genotyping_evidence_records(
        &outcome.genotyping_reads,
        &header,
        &region.contig,
        &gt_filters,
    );
    let mut gt_before = LocusPileupState::from_records_qname_deduped(
        &outcome.genotyping_reads,
        &header,
        &region.contig,
        &gt_filters,
    );
    let ref_bytes = ref_cache
        .get_interval_bytes(&dict, &region.contig, 92305755, 92305779)
        .expect("ref");

    let java_expect: &[(u64, i32)] = &[
        (92305755, 6),
        (92305756, 0),
        (92305759, 6),
        (92305760, 0),
        (92305762, 6),
        (92305779, 0),
    ];

    for &(pos, java_gq) in java_expect {
        let offset = (pos - 92305755) as usize;
        let ref_base = *ref_bytes.get(offset).unwrap_or(&b'N');
        let pile_after = gt_after
            .pileup_at(&outcome.genotyping_reads, &gt_filters, pos, ref_base)
            .expect("after");
        let pile_before = gt_before
            .pileup_at(&outcome.genotyping_reads, &gt_filters, pos, ref_base)
            .expect("before");
        let gl_after = capped_genotype_likelihoods_by_hom_ref(&pile_after, &config);
        let gl_before = capped_genotype_likelihoods_by_hom_ref(&pile_before, &config);
        let locus = loci
            .iter()
            .find(|l| l.position_1based as u64 == pos)
            .unwrap_or_else(|| panic!("missing locus {pos}"));
        eprintln!(
            "pos={pos} java_gq={java_gq} locus_gq={} dp={} after_dp={} alt_after={} gq_after={} alt_before={} gq_before={} gl_after={gl_after:?}",
            locus.gq,
            locus.dp,
            pile_after.len(),
            pile_after.iter().filter(|o| o.is_alt).count(),
            reference_gq_from_log10_gl(&gl_after),
            pile_before.iter().filter(|o| o.is_alt).count(),
            reference_gq_from_log10_gl(&gl_before),
        );
    }
}
