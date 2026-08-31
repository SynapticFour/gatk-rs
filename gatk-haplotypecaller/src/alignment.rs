//! Haplotype-to-reference Smith–Waterman alignment (GATK `CigarUtils.calculateCigar`).

pub use crate::cigar::{Cigar, CigarElement, CigarOperator};
pub use crate::haplotype_cigar::{
    calculate_haplotype_cigar, calculate_haplotype_cigar_for_assembly,
    calculate_haplotype_cigar_for_assembly_with_offset, calculate_haplotype_cigar_java_padded_sw,
    trace_find_best_paths_gates, FindBestPathsGateTrace, HaplotypeAssemblyCigar,
};
pub use crate::smith_waterman::{SwOverhangStrategy, SwParameters};
