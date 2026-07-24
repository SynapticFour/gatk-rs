//! GATK `HaplotypeCallerEngine.filterNonPassingReads` (assembly / genotyping path).

use crate::assembly_region_iterator::AssemblyRegion;
use crate::read_pre_len::passes_read_length_filter;
use crate::read_pre_mate::passes_mate_on_same_contig_or_no_mapped_mate;
use crate::read_pre_mq::{passes_assembly_mq_filter, GATK_ASSEMBLY_MQ_FILTER_THRESHOLD};
use rust_htslib::bam::record::Aux;
use rust_htslib::bam::Record;

/// GATK `HaplotypeCallerEngine.filterNonPassingReads` options.
/// # Invariants
/// MQ filter uses `mq_threshold` (HC default 20); length and mate checks are fixed GATK rules.
/// When `keep_read_group` is set, reads without matching `RG` tag are dropped.
/// # Ownership
/// [`Clone`] config borrowed by [`filter_non_passing_reads`]; region owns mutable read records.
/// # Mutation
/// Config is immutable; filtering mutates the region's read vector in place via `retain`.
/// # Biological assumptions
/// Low-MQ, wrong-contig mates, and truncated reads should not enter local assembly.
/// # Java equivalence
/// GATK `HaplotypeCallerEngine.filterNonPassingReads` + `hcArgs.mappingQualityThreshold` / `keepRG`.
#[derive(Debug, Clone)]
pub struct AssemblyReadFilterConfig {
    /// GATK `hcArgs.mappingQualityThreshold` (default 20).
    pub mq_threshold: u8,
    /// GATK `hcArgs.keepRG` — when set, only reads with this read group are kept.
    pub keep_read_group: Option<String>,
}

impl Default for AssemblyReadFilterConfig {
    fn default() -> Self {
        Self::gatk_defaults()
    }
}

impl AssemblyReadFilterConfig {
    pub fn gatk_defaults() -> Self {
        Self {
            mq_threshold: GATK_ASSEMBLY_MQ_FILTER_THRESHOLD,
            keep_read_group: None,
        }
    }
}

/// Remove reads failing length / MQ / mate / keepRG gates before genotyping.
pub fn filter_non_passing_reads(region: &mut AssemblyRegion, config: &AssemblyReadFilterConfig) {
    region.reads.retain(|rec| passes_assembly_read(rec, config));
}

#[inline]
pub fn passes_assembly_read(rec: &Record, config: &AssemblyReadFilterConfig) -> bool {
    passes_read_length_filter(rec)
        && passes_assembly_mq_filter(rec, config.mq_threshold)
        && passes_mate_on_same_contig_or_no_mapped_mate(rec)
        && passes_keep_read_group(rec, config.keep_read_group.as_deref())
}

fn passes_keep_read_group(rec: &Record, keep: Option<&str>) -> bool {
    let Some(expected) = keep else {
        return true;
    };
    match rec.aux(b"RG") {
        Ok(Aux::String(rg)) => rg == expected,
        _ => false,
    }
}

pub fn filter_non_passing_reads_default(region: &mut AssemblyRegion) {
    filter_non_passing_reads(region, &AssemblyReadFilterConfig::gatk_defaults());
}

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::record::{Aux, Cigar, CigarString};
    use rust_htslib::bam::Record;

    fn rec_with_rg(rg: &str) -> Record {
        let mut r = Record::new();
        r.set(
            b"r",
            Some(&CigarString::from(vec![Cigar::Match(10)])),
            b"AAAAAAAAAA",
            &vec![30u8; 10],
        );
        r.set_mapq(30);
        r.push_aux(b"RG", Aux::String(rg)).unwrap();
        r
    }

    #[test]
    fn keep_rg_filters_other_groups() {
        let cfg = AssemblyReadFilterConfig {
            keep_read_group: Some("RG1".to_string()),
            ..AssemblyReadFilterConfig::gatk_defaults()
        };
        assert!(passes_assembly_read(&rec_with_rg("RG1"), &cfg));
        assert!(!passes_assembly_read(&rec_with_rg("RG2"), &cfg));
    }
}
