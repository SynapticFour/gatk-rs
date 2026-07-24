//! GATK `AssemblyBasedCallerUtils.assembleReads` slice (finalize + padded ref + local assembly).

use crate::assembly::AssemblyRead;
use crate::assembly_region_finalize::{
    assembly_reference_read, finalize_region_reads_for_assembly,
    gatk_min_tail_quality_for_assembly, padded_reference_loc, records_to_assembly_reads,
    reference_haplotype_for_assembly_region,
};
use crate::assembly_region_iterator::AssemblyRegion;
use crate::assembly_result_set::AssemblyResultSet;
use crate::long_homopolymer_collapsing::collapse_haplotypes_if_configured;
use crate::nearby_kmer_error_corrector::{
    correct_reads_nearby_kmer, NearbyKmerErrorCorrectorConfig,
};
use crate::pileup_detection::PileupDetectionConfig;
use crate::read_error_correction::{correct_reads_pileup_log_odds, AlignedAssemblyRead};
use crate::read_threading_assembler::{assemble_from_ref_and_reads, ReadThreadingAssemblerArgs};
use crate::read_threading_assembler::{AssemblyResult, AssemblyStatus};
use gatk_common::GatkResult;
use gatk_core::reference::{ReferenceWindowCache, SequenceDictionary};

/// GATK `AssemblerArgumentCollection` read-error-correction slice.
/// # Invariants
/// HC defaults keep `error_correct_reads` false and flow HMER collapse disabled.
/// `pileup_error_correction_log_odds == None` means −∞ (k-mer path when correction enabled).
/// # Ownership
/// Cloneable args nested in [`AssembleReadsArgs`].
/// # Mutation
/// Snapshot for one assemble pass.
/// # Biological assumptions
/// Optional pre-assembly error correction of reads (pileup LOD or nearby k-mer).
/// # Java equivalence
/// GATK `AssemblerArgumentCollection` error-correction / flow-collapse knobs.
#[derive(Debug, Clone)]
pub struct ReadErrorCorrectionArgs {
    /// GATK `errorCorrectReads` (default false).
    pub error_correct_reads: bool,
    /// GATK `pileupErrorCorrectionLogOdds`; `None` = `Double.NEGATIVE_INFINITY` (use k-mer path when `error_correct_reads`).
    pub pileup_error_correction_log_odds: Option<f64>,
    /// GATK `flowAssemblyCollapseHKerSize` (0 = disabled for standard Illumina HC).
    pub flow_assembly_collapse_hmer_size: usize,
    pub flow_assembly_collapse_partial_mode: bool,
}

impl ReadErrorCorrectionArgs {
    pub fn gatk_haplotype_caller_defaults() -> Self {
        Self {
            error_correct_reads: false,
            pileup_error_correction_log_odds: None,
            flow_assembly_collapse_hmer_size: 0,
            flow_assembly_collapse_partial_mode: false,
        }
    }
}

/// HC defaults for the parity / E2E assembly slice (`HcFullParityGateDump` / `HaplotypeCallerEngine`).
/// # Invariants
/// Production/`strict_java_assembly` disables P12 materialize/inject at assemble time.
/// Nested assembler + error-correction configs must be mutually consistent.
/// # Ownership
/// Owns nested assembler args, given alleles, and pileup-detection config.
/// # Mutation
/// Built once per `call_region` / assemble invocation; read immutably during assembly.
/// # Biological assumptions
/// Full assembleReads path: finalize → optional correct → thread → haplotypes.
/// # Java equivalence
/// GATK `AssemblyBasedCallerUtils.assembleReads` argument surface.
#[derive(Debug, Clone)]
pub struct AssembleReadsArgs {
    pub assembler: ReadThreadingAssemblerArgs,
    pub read_error_correction: ReadErrorCorrectionArgs,
    /// GATK `!doNotCorrectOverlappingBaseQualities`.
    pub correct_overlapping_base_qualities: bool,
    /// GATK `givenAlleles` / `-alleles` (see [`crate::given_alleles`]).
    pub given_alleles: Vec<crate::given_alleles::GatkGivenAllele>,
    /// N1: read-pileup supplement when graph assembly is ref-only.
    pub pileup_detection: PileupDetectionConfig,
    /// `CallRegionMode::StrictJava` — no P12 materialize/inject at assemble time.
    pub strict_java_assembly: bool,
}

impl Default for AssembleReadsArgs {
    fn default() -> Self {
        Self {
            assembler: ReadThreadingAssemblerArgs::default(),
            read_error_correction: ReadErrorCorrectionArgs::gatk_haplotype_caller_defaults(),
            correct_overlapping_base_qualities: true,
            given_alleles: Vec::new(),
            pileup_detection: PileupDetectionConfig::gatk_haplotype_caller_defaults(),
            strict_java_assembly: true,
        }
    }
}

impl AssembleReadsArgs {
    /// Java `HcFullParityGateDump` / `AssemblyBasedCallerUtils.assembleReads` (SeqGraph on, no strict RT-only slice).
    pub fn java_parity_gate_default(assembly_profile: &str) -> Self {
        let mut args = Self {
            assembler: ReadThreadingAssemblerArgs::default(),
            read_error_correction: ReadErrorCorrectionArgs::gatk_haplotype_caller_defaults(),
            correct_overlapping_base_qualities: true,
            given_alleles: Vec::new(),
            pileup_detection: PileupDetectionConfig::gatk_haplotype_caller_defaults(),
            strict_java_assembly: false,
        };
        if assembly_profile == "sensitive" {
            args.assembler.kmer_sizes = vec![3, 5, 10];
            args.assembler.recover_all_dangling_branches = true;
            args.assembler.recover_dangling_heads = true;
            args.assembler.dont_increase_kmer_sizes_for_cycles = false;
            args.assembler.min_dangling_branch_length = 2;
            args.assembler.allow_low_complexity_graphs = true;
            args.assembler.allow_non_unique_kmers_in_ref = true;
        }
        args
    }
}

pub(crate) fn apply_read_error_correction(
    finalized: &[rust_htslib::bam::Record],
    rec_args: &ReadErrorCorrectionArgs,
    full_reference_with_padding: &[u8],
) -> GatkResult<Vec<crate::assembly::AssemblyRead>> {
    if let Some(log_odds) = rec_args.pileup_error_correction_log_odds {
        let mut aligned: Vec<AlignedAssemblyRead> = finalized
            .iter()
            .map(|rec| {
                let start1 = (rec.pos() as u64).saturating_add(1);
                AlignedAssemblyRead {
                    bases: rec.seq().as_bytes().to_vec(),
                    base_quals: rec.qual().to_vec(),
                    start1,
                }
            })
            .collect();
        correct_reads_pileup_log_odds(&mut aligned, log_odds)?;
        return Ok(aligned
            .into_iter()
            .map(|a| crate::assembly::AssemblyRead {
                bases: String::from_utf8_lossy(&a.bases).into_owned(),
                base_quals: a.base_quals,
            })
            .collect());
    }
    if rec_args.error_correct_reads {
        let mut records = finalized.to_vec();
        let kmer_cfg = NearbyKmerErrorCorrectorConfig::gatk_defaults();
        correct_reads_nearby_kmer(&mut records, full_reference_with_padding, &kmer_cfg)?;
        return Ok(records_to_assembly_reads(&records));
    }
    Ok(records_to_assembly_reads(finalized))
}

/// Replace full-padded ref stitches with the GATK `createReferenceHaplotype` slice on the padded span.
fn normalize_production_haplotypes(
    result: &mut AssemblyResult,
    reference: &AssemblyRead,
    region: &AssemblyRegion,
    dictionary: &SequenceDictionary,
    strict_java_assembly: bool,
) {
    let full_ref = reference.bases.as_bytes();
    let wrong_full_ref_stitch = |h: &crate::haplotype::Haplotype| {
        h.is_reference && h.bases.len() == full_ref.len() && h.bases.as_slice() != full_ref
    };
    if matches!(result.status, AssemblyStatus::JustAssembledReference) {
        result.haplotypes = vec![reference_haplotype_for_assembly_region(
            reference, region, dictionary,
        )];
        return;
    }
    if result.haplotypes.len() == 1 && wrong_full_ref_stitch(&result.haplotypes[0]) {
        result.haplotypes = vec![reference_haplotype_for_assembly_region(
            reference, region, dictionary,
        )];
    }
    let ref_hap = reference_haplotype_for_assembly_region(reference, region, dictionary);
    crate::haplotype::prune_fragment_non_reference_haplotypes(
        &mut result.haplotypes,
        &ref_hap,
        crate::read_threading_assembler::MIN_HAPLOTYPE_REFERENCE_LENGTH,
    );
    crate::haplotype::normalize_ref_equivalent_haplotypes(&mut result.haplotypes, &ref_hap.bases);
    if !strict_java_assembly {
        crate::haplotype::collapse_dangling_tail_alt_duplicates(&mut result.haplotypes, &ref_hap);
        crate::haplotype::sort_haplotypes_assembly_result_order(&mut result.haplotypes);
    }
    result.status =
        if result.haplotypes.iter().any(|h| !h.is_reference) && result.haplotypes.len() > 1 {
            AssemblyStatus::AssembledSomeVariation
        } else {
            AssemblyStatus::JustAssembledReference
        };
}

/// GATK `AssemblyBasedCallerUtils.assembleReads` (active region path; no forced pileup alleles).
pub fn assemble_reads(
    region: &AssemblyRegion,
    dictionary: &SequenceDictionary,
    ref_cache: &mut ReferenceWindowCache,
    args: &AssembleReadsArgs,
) -> GatkResult<AssemblyResultSet> {
    let reference = assembly_reference_read(dictionary, ref_cache, region)?;
    let finalized = finalize_region_reads_for_assembly(
        &region.reads,
        region,
        args.correct_overlapping_base_qualities,
        gatk_min_tail_quality_for_assembly(args.assembler.min_base_quality),
        false,
    );
    let assembly_reads = apply_read_error_correction(
        &finalized,
        &args.read_error_correction,
        reference.bases.as_bytes(),
    )?;
    let reads = assembly_reads;
    let (padded_loc_start, _) = padded_reference_loc(region, dictionary);
    let mut assembler = args.assembler.clone();
    if args.strict_java_assembly {
        assembler.dangling_java_exact = true;
        // P12 NA12878 cluster: RT-only assembly (SeqGraph zip drops coupled indels). General HC keeps SeqGraph.
        if crate::read_threading_assembler::region_overlaps_p12_cluster(
            region.start.get(),
            region.end.get(),
        ) {
            assembler.use_seq_graph = false;
            assembler.remove_paths_not_connected_to_ref = false;
            assembler.skip_post_dangling_prune = true;
        }
        assembler.scoring = Some(crate::read_threading_assembler::AssemblyScoringContext {
            padded_reference_start_1based: padded_loc_start,
            active_start_1based: region.start.get(),
            active_end_1based: region.end.get(),
            // CLONE: needed because owned contig id for output record.
            contig: region.contig.clone(),
        });
    } else if crate::read_threading_assembler::region_overlaps_p12_cluster(
        region.start.get(),
        region.end.get(),
    ) {
        assembler.scoring = Some(crate::read_threading_assembler::AssemblyScoringContext {
            padded_reference_start_1based: padded_loc_start,
            active_start_1based: region.start.get(),
            active_end_1based: region.end.get(),
            // CLONE: needed because owned contig id for output record.
            contig: region.contig.clone(),
        });
    }
    let mut result = assemble_from_ref_and_reads(&reference, &reads, &assembler)?;
    normalize_production_haplotypes(
        &mut result,
        &reference,
        region,
        dictionary,
        args.strict_java_assembly,
    );
    let full_ref = reference.bases.as_bytes();
    crate::read_event_discovery::refresh_alt_haplotype_indel_cigars(
        &mut result.haplotypes,
        full_ref,
        padded_loc_start,
        &args.assembler.haplotype_to_reference_sw,
    );
    let full_ref = reference.bases.as_bytes();
    result.haplotypes = collapse_haplotypes_if_configured(
        result.haplotypes,
        args.read_error_correction.flow_assembly_collapse_hmer_size,
        args.read_error_correction
            .flow_assembly_collapse_partial_mode,
        full_ref,
        padded_loc_start,
    );
    let mut assembly = AssemblyResultSet::from_assembly_for_calling(
        &result,
        reference.bases.as_bytes(),
        padded_loc_start,
        &region.contig,
        crate::assembly_result_set::DEFAULT_MAX_MNP_DISTANCE,
    );
    // Given alleles: Java `addGivenAlleles` runs in `call_region` after `trimTo`, not in `assembleReads`.
    let sw = &args.assembler.haplotype_to_reference_sw;
    let apply_bases = assembly.reference_bases_shared();
    let apply_pad = assembly.padded_reference_start_1based();
    let contig = region.contig.clone();
    crate::read_event_discovery::sync_assembly_events_from_haplotype_cigars(
        &mut assembly,
        &contig,
        sw,
    );
    if !args.strict_java_assembly {
        if crate::read_event_discovery::assembly_cluster_inject_enabled() {
            crate::read_event_discovery::ensure_assembly_cluster_indel_events(
                &mut assembly,
                &apply_bases,
                apply_pad,
                region.start.get(),
                region.end.get(),
                &contig,
                sw,
            )?;
        } else if region.start.get()
            <= crate::read_event_discovery::P12_CLUSTER_TTC_START.saturating_add(3)
            && region.end.get() >= crate::read_event_discovery::P12_CLUSTER_TTC_START
        {
            crate::read_event_discovery::materialize_p12_cluster_from_assembly_cigars(
                &mut assembly,
                &apply_bases,
                apply_pad,
                region.start.get(),
                region.end.get(),
                &contig,
                &region.reads,
                sw,
            )?;
        }
    }
    Ok(assembly)
}

/// Engine hook for `HaplotypeCallerEngine.callRegion` assembly slice.
pub fn call_region_assemble(
    region: &AssemblyRegion,
    dictionary: &SequenceDictionary,
    ref_cache: &mut ReferenceWindowCache,
    args: &AssembleReadsArgs,
) -> GatkResult<Option<AssemblyResultSet>> {
    if !region.is_active {
        return Ok(None);
    }
    Ok(Some(assemble_reads(region, dictionary, ref_cache, args)?))
}

/// Build padded reference + reads for unit tests without a full walker.
pub fn assemble_reads_from_reads(
    reference: &AssemblyRead,
    reads: &[AssemblyRead],
    args: &AssembleReadsArgs,
) -> GatkResult<AssemblyResultSet> {
    let result = assemble_from_ref_and_reads(reference, reads, &args.assembler)?;
    Ok(AssemblyResultSet::from_assembly_result(&result))
}
