//! Java-aligned genotyping (stored events, changeEvidence, emit threshold).

use gatk_haplotypecaller::bio_ids::{HaplotypeIndex, ReadIndex};
use gatk_haplotypecaller::cigar::{Cigar, CigarOperator};
use gatk_haplotypecaller::engine::RegionReadLikelihood;
use gatk_haplotypecaller::event_map::VariationEvent;
use gatk_haplotypecaller::genome_loc::GenomePosition;
use gatk_haplotypecaller::genotyping::ReadLikelihoodRow;
use gatk_haplotypecaller::haplotype::Haplotype;
use gatk_haplotypecaller::hc_emit_policy::passes_emit_for_variation_event;
use gatk_haplotypecaller::hc_genotyping_engine::{
    assign_genotype_likelihoods_for_region, biallelic_genotype_log10_likelihoods_gatk,
    marginalize_rows_to_biallelic_alleles, HcGenotypingConfig,
};
use gatk_haplotypecaller::read_realignment::change_evidence_to_best_haplotype;

fn ref_hap_with_cigar(bases: &[u8]) -> Haplotype {
    let mut h = Haplotype::new(bases, true);
    let mut c = Cigar::new();
    c.push(bases.len(), CigarOperator::Match);
    h.cigar = Some(c);
    h
}

fn del_alt_hap(bases: &[u8], del_at_ref: usize, del_len: usize) -> Haplotype {
    let mut alt_bases = bases[0..del_at_ref].to_vec();
    alt_bases.extend_from_slice(&bases[del_at_ref + del_len..]);
    let mut h = Haplotype::new(alt_bases, false);
    h.score = 100.0;
    let mut c = Cigar::new();
    c.push(del_at_ref, CigarOperator::Match);
    c.push(del_len, CigarOperator::Deletion);
    c.push(
        h.bases.len().saturating_sub(del_at_ref),
        CigarOperator::Match,
    );
    h.cigar = Some(c);
    h
}

#[test]
fn change_evidence_preserves_full_read_hap_matrix() {
    let likelihoods = vec![
        RegionReadLikelihood {
            read_index: ReadIndex::new(0),
            haplotype_index: HaplotypeIndex::new(0),
            log10_likelihood: -1.0,
        },
        RegionReadLikelihood {
            read_index: ReadIndex::new(0),
            haplotype_index: HaplotypeIndex::new(1),
            log10_likelihood: -0.2,
        },
    ];
    let out = change_evidence_to_best_haplotype(likelihoods.clone(), &[1]);
    assert_eq!(out.len(), 2);
    assert_eq!(out, likelihoods);
}

#[test]
fn biallelic_gl_favors_het_when_reads_support_alt_hap() {
    let rows = vec![
        ReadLikelihoodRow {
            read_index: 0,
            read_id: "r0".into(),
            haplotype_log10_likelihoods: vec![-3.0, -0.1],
        },
        ReadLikelihoodRow {
            read_index: 0,
            read_id: "r1".into(),
            haplotype_log10_likelihoods: vec![-3.0, -0.1],
        },
    ];
    let gls = biallelic_genotype_log10_likelihoods_gatk(&rows, 0, 1);
    assert!(gls[1] > gls[0], "0/1 should beat 0/0: {gls:?}");
    let config = HcGenotypingConfig::strict_java();
    let format = gatk_haplotypecaller::genotyping::emit_genotype_format_fields(
        &gls,
        &gatk_haplotypecaller::hc_genotyping_engine::biallelic_allele_depths_from_rows(&rows, 0, 1),
    )
    .expect("format");
    let event = VariationEvent {
        contig: "chr1".into(),
        start_1based: GenomePosition::new_1based(10),
        end_1based: GenomePosition::new_1based(12),
        ref_allele: "TTC".into(),
        alt_allele: "T".into(),
    };
    assert!(
        passes_emit_for_variation_event(&event, &gls, &format, config.stand_emit_confidence, &[])
            .unwrap(),
        "Java-style AF emit should pass for alt-favored GLs"
    );
}

/// Java hap EventMap walk (`genotype_stored_events_only: false`) + PairHMM GL for TTC/T.
#[test]
fn assign_genotype_hap_walk_deletion_event_emits_call() {
    let pad = 100u64;
    let mut ref_bases = vec![b'N'; 30];
    ref_bases[9] = b'T';
    ref_bases[10] = b'T';
    ref_bases[11] = b'C';
    let ref_hap = ref_hap_with_cigar(&ref_bases);
    let alt = del_alt_hap(&ref_bases, 10, 2);
    let haps = vec![alt, ref_hap];
    let ref_bytes = ref_bases.clone();

    let mut likelihoods = Vec::new();
    for ri in 0..4 {
        likelihoods.push(RegionReadLikelihood {
            read_index: ReadIndex::new(ri),
            haplotype_index: HaplotypeIndex::new(0),
            log10_likelihood: -0.15,
        });
        likelihoods.push(RegionReadLikelihood {
            read_index: ReadIndex::new(ri),
            haplotype_index: HaplotypeIndex::new(1),
            log10_likelihood: -4.0,
        });
    }
    let ttc = VariationEvent {
        contig: "chr1".into(),
        start_1based: GenomePosition::new_1based(pad + 9),
        end_1based: GenomePosition::new_1based(pad + 11),
        ref_allele: "TTC".into(),
        alt_allele: "T".into(),
    };
    let config = HcGenotypingConfig::strict_java();
    let rows = gatk_haplotypecaller::hc_genotyping_engine::region_likelihoods_to_rows(
        &likelihoods,
        haps.len(),
    );
    let marg = marginalize_rows_to_biallelic_alleles(
        &rows,
        &[HaplotypeIndex::new(1)],
        &[HaplotypeIndex::new(0)],
    );
    assert_eq!(marg.len(), 4);
    assert!(marg[0].haplotype_log10_likelihoods[1] > -1.0);

    let result = assign_genotype_likelihoods_for_region(
        &likelihoods,
        &[],
        &[],
        None,
        &haps,
        &ref_bytes,
        pad,
        &ref_bytes,
        pad,
        pad,
        pad + 29,
        "chr1",
        1,
        &config,
        std::slice::from_ref(&ttc),
        &[],
    )
    .expect("assign");
    assert_eq!(result.calls.len(), 1);
    assert_eq!(result.calls[0].event.ref_allele, "TTC");
    assert_eq!(result.calls[0].event.alt_allele, "T");
}
