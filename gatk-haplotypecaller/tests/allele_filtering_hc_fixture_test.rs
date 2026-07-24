//! `AlleleFilteringHC`-style haplotype pruning on synthetic likelihoods.

use gatk_haplotypecaller::allele_filtering::filter_assembly_and_likelihoods;
use gatk_haplotypecaller::assembly_result_set::AssemblyResultSet;
use gatk_haplotypecaller::bio_ids::{HaplotypeIndex, ReadIndex};
use gatk_haplotypecaller::engine::RegionReadLikelihood;
use gatk_haplotypecaller::haplotype::Haplotype;
use gatk_haplotypecaller::read_threading_assembler::AssemblyStatus;

#[test]
fn filter_keeps_ref_plus_top_scored_non_ref_haplotypes() {
    let mut haps = vec![Haplotype::new(b"ACGTACGTACGT", true)];
    for i in 0..14 {
        let mut alt = Haplotype::new(format!("ACGTACGTACG{i}").as_bytes(), false);
        alt.score = if i == 0 { 50.0 } else { 0.001 * (i as f64) };
        haps.push(alt);
    }
    let mut assembly = AssemblyResultSet::from_assembly_result(
        &gatk_haplotypecaller::read_threading_assembler::AssemblyResult {
            status: AssemblyStatus::AssembledSomeVariation,
            kmer_size: 10,
            haplotypes: haps,
            event_maps: vec![],
        },
    );
    let likelihoods: Vec<RegionReadLikelihood> = (0..assembly.haplotypes.len())
        .map(|hi| RegionReadLikelihood {
            read_index: ReadIndex::new(0),
            haplotype_index: HaplotypeIndex::new(hi),
            log10_likelihood: if hi == 0 {
                -0.1
            } else if hi == 1 {
                -0.2
            } else {
                -40.0
            },
        })
        .collect();
    let filtered = filter_assembly_and_likelihoods(
        &mut assembly,
        likelihoods,
        gatk_haplotypecaller::allele_filter_options::AlleleFilterOptions::unrestricted(),
    )
    .expect("filter");
    assert!(assembly.haplotypes.len() <= 13);
    assert!(assembly.haplotypes.iter().any(|h| h.is_reference));
    assert!(filtered.iter().any(|r| r.haplotype_index.get() == 1));
    assert!(!filtered.iter().any(|r| r.haplotype_index.get() >= 13));
}

#[test]
fn filter_ranks_by_read_likelihood_sum_not_haplotype_score() {
    let mut haps = vec![Haplotype::new(b"ACGTACGTACGT", true)];
    let mut high_score_weak_ll = Haplotype::new(b"ACGTACGTACG0", false);
    high_score_weak_ll.score = 100.0;
    let mut strong_ll_low_score = Haplotype::new(b"ACGTACGTACG1", false);
    strong_ll_low_score.score = 0.001;
    haps.push(high_score_weak_ll);
    haps.push(strong_ll_low_score);
    for i in 2..16 {
        let mut alt = Haplotype::new(format!("ACGTACGTACG{i}").as_bytes(), false);
        alt.score = 50.0;
        haps.push(alt);
    }
    let mut assembly = AssemblyResultSet::from_assembly_result(
        &gatk_haplotypecaller::read_threading_assembler::AssemblyResult {
            status: AssemblyStatus::AssembledSomeVariation,
            kmer_size: 10,
            haplotypes: haps,
            event_maps: vec![],
        },
    );
    let likelihoods: Vec<RegionReadLikelihood> = assembly
        .haplotypes
        .iter()
        .enumerate()
        .flat_map(|(hi, _)| {
            (0..3).map(move |ri| RegionReadLikelihood {
                read_index: ReadIndex::new(ri),
                haplotype_index: HaplotypeIndex::new(hi),
                log10_likelihood: if hi == 2 {
                    -0.1
                } else if hi == 1 {
                    -50.0
                } else {
                    -0.5
                },
            })
        })
        .collect();
    filter_assembly_and_likelihoods(
        &mut assembly,
        likelihoods,
        gatk_haplotypecaller::allele_filter_options::AlleleFilterOptions::unrestricted(),
    )
    .expect("filter");
    assert!(
        assembly
            .haplotypes
            .iter()
            .any(|h| h.bases == b"ACGTACGTACG1"),
        "hap with best read-LL sum must survive filter"
    );
}
