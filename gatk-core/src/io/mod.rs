//! I/O module for GATK-RS
//! This module provides efficient parsers and writers for common genomic
//! file formats including FASTA, FASTQ, BAM/SAM, and VCF.

pub mod bam;
pub mod fasta;
pub mod fastq;
pub mod hts_copy;
pub mod sam;
pub mod vcf;

// Re-export commonly used items
pub use bam::*;
pub use fasta::*;
pub use fastq::*;
pub use hts_copy::{
    alignment_header_canonical_hd_sq_rg, copy_alignments_with_htslib,
    copy_alignments_with_htslib_with_reference, count_records_in_region_indexed,
    qnames_in_region_indexed, validate_bam_file,
};
pub use sam::*;
pub use vcf::*;
