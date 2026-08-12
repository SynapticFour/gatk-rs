use super::*;
use crate::alignment::SwParameters;
use crate::assembly_region_iterator::AssemblyRegion;
use crate::assembly_result_set::AssemblyResultSet;
use crate::cigar::{Cigar, CigarOperator};
use crate::genome_loc::{GenomeLoc, GenomePosition};
use crate::haplotype::Haplotype;
use crate::read_threading_assembler::{AssemblyResult, AssemblyStatus};

#[test]
fn preserve_untrimmed_indel_haplotypes_reattaches_when_trim_loses_indel_cigar() {
    let sw = SwParameters::gatk_haplotype_to_reference();
    let ref_bases = b"ACGTACGTACGTACGT".to_vec();
    let span = GenomeLoc::new(100, 115);
    let mut ref_cigar = Cigar::new();
    ref_cigar.push(ref_bases.len(), CigarOperator::Match);
    // CLONE: needed because haplotype constructor takes owned bases.
    let mut untrimmed_ref = Haplotype::new(ref_bases.clone(), true);
    // CLONE: needed because haplotype owns CIGAR.
    untrimmed_ref.cigar = Some(ref_cigar.clone());
    untrimmed_ref.genome_loc = Some(span);
    let mut untrimmed_alt = Haplotype::new(b"ACGTACGTTACGTACGT".to_vec(), false);
    let mut indel_cigar = Cigar::new();
    indel_cigar.push(8, CigarOperator::Match);
    indel_cigar.push(1, CigarOperator::Insertion);
    indel_cigar.push(8, CigarOperator::Match);
    untrimmed_alt.cigar = Some(indel_cigar);
    untrimmed_alt.alignment_start_hap_wrt_ref = 0;
    let untrimmed = AssemblyResultSet::from_assembly_for_calling(
        &AssemblyResult {
            status: AssemblyStatus::AssembledSomeVariation,
            kmer_size: 10,
            haplotypes: vec![untrimmed_alt, untrimmed_ref],
            event_maps: Vec::new(),
        },
        ref_bases.as_slice(),
        100,
        "2",
        0,
    );
    let trimmed_region = AssemblyRegion {
        contig: "2".into(),
        start: GenomePosition::new_1based(102),
        end: GenomePosition::new_1based(110),
        is_active: true,
        extended_start: GenomePosition::new_1based(100),
        extended_end: GenomePosition::new_1based(115),
        extension: 0,
        reads: Vec::new(),
        read_qnames: Vec::new(),
        reference: crate::reference_context::ReferenceContext::empty(),
        features: crate::feature_context::FeatureContext::empty(),
        pileup_loci: Vec::new(),
    };
    let mut trimmed_ref = Haplotype::new(ref_bases[2..11].to_vec(), true);
    trimmed_ref.cigar = Some({
        let mut c = Cigar::new();
        c.push(9, CigarOperator::Match);
        c
    });
    trimmed_ref.genome_loc = Some(GenomeLoc::new(100, 115));
    let mut assembly = AssemblyResultSet::from_assembly_for_calling(
        &AssemblyResult {
            status: AssemblyStatus::AssembledSomeVariation,
            kmer_size: 10,
            haplotypes: vec![trimmed_ref],
            event_maps: Vec::new(),
        },
        &ref_bases[2..11],
        100,
        "2",
        0,
    );
    preserve_untrimmed_indel_haplotypes(&untrimmed, &mut assembly, &trimmed_region, &sw);
    assert!(
        assembly
            .haplotypes
            .iter()
            .filter(|h| !h.is_reference)
            .any(|h| {
                h.cigar
                    .as_ref()
                    .is_some_and(|c| c.elements.iter().any(|e| e.operator.is_indel()))
            }),
        "G-6: indel alt hap must survive trim when assembly lost indel CIGAR"
    );
}
