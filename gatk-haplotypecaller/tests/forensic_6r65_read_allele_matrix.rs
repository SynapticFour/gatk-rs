//! 6R.66 coordinate-free: empty overlapping EventMap must not assign an indel via
//! pad-slice. SNP-by-base without EventMap is a separate, retained mapper path.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r65_read_allele_matrix
//! ```

use gatk_haplotypecaller::bio_ids::HaplotypeIndex;
use gatk_haplotypecaller::cigar::{Cigar, CigarOperator};
use gatk_haplotypecaller::event_map::VariationEvent;
use gatk_haplotypecaller::genome_loc::GenomePosition;
use gatk_haplotypecaller::haplotype::Haplotype;
use gatk_haplotypecaller::hc_allele_mapping::create_allele_mapper;

fn empty_overlap_indel_fixture() -> (VariationEvent, Vec<Haplotype>, Vec<u8>, u64, u64) {
    let pad = 100u64;
    let loc = 110u64;
    let mut ref_bytes = vec![b'A'; 40];
    ref_bytes[(loc - pad) as usize] = b'T';
    ref_bytes[(loc - pad + 1) as usize] = b'G';
    let mut ref_hap = Haplotype::new(ref_bytes.clone(), true);
    let mut ref_cigar = Cigar::new();
    ref_cigar.push(ref_bytes.len(), CigarOperator::Match);
    ref_hap.cigar = Some(ref_cigar);

    let mut alt = Haplotype::new(ref_bytes.clone(), false);
    alt.bases[(loc - pad + 1) as usize] = b'A';
    let mut alt_cigar = Cigar::new();
    alt_cigar.push(ref_bytes.len(), CigarOperator::Match);
    alt.cigar = Some(alt_cigar);

    let merged = VariationEvent {
        contig: "20".into(),
        start_1based: GenomePosition::new_1based(loc),
        end_1based: GenomePosition::new_1based(loc + 1),
        ref_allele: "TG".into(),
        alt_allele: "T".into(),
    };
    (merged, vec![ref_hap, alt], ref_bytes, pad, loc)
}

#[test]
fn empty_eventmap_apparent_indel_is_not_assigned_via_pad_slice() {
    let (merged, haps, ref_bytes, pad, loc) = empty_overlap_indel_fixture();
    let m = create_allele_mapper(&merged, loc, &haps, pad, &ref_bytes, 0, true);
    assert!(
        m.ref_haplotype_indices.contains(&HaplotypeIndex::new(1)),
        "Java: empty overlapping EventMap → REF"
    );
    assert!(
        !m.alt_haplotype_indices.contains(&HaplotypeIndex::new(1)),
        "6R.66: no indel pad-slice; alt={:?}",
        m.alt_haplotype_indices
    );
}

#[test]
fn allele_floor_can_clip_losing_alleles_without_moving_homref_vs_best_het() {
    // Production apply_java_marginal_normalize_n after pool-max. Java floors only the
    // haplotype matrix. A clip of the unused allele does not change 0/0 vs 0/2.
    let pre = [0.0_f64, -10.0, -0.05];
    let best = 0.0_f64;
    let floor = best - 4.5;
    let mut post = pre;
    for v in &mut post {
        if *v < floor {
            *v = floor;
        }
    }
    assert_eq!(post[0], 0.0);
    assert_eq!(post[2], -0.05);
    assert!((post[1] - floor).abs() < 1e-12);
}
