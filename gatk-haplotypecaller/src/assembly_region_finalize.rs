//! GATK `AssemblyBasedCallerUtils.finalizeRegion` + assembly reference padding (GAP-E-06).

use crate::assembly::AssemblyRead;
use crate::assembly_region_iterator::AssemblyRegion;
use crate::cigar::{Cigar, CigarOperator};
use crate::fragment_overlap::clean_overlapping_read_pairs;
use crate::haplotype::Haplotype;
use crate::read_pre_len::unclipped_read_length;
use crate::read_unclip::{
    apply_hc_softclip_pre_step, hard_clip_adaptor_sequence, hard_clip_low_qual_ends,
    hard_clip_to_region, normalize_record_cigar, read_has_positive_cigar_length,
    soft_clip_low_qual_ends, HcSoftclipPolicy,
};
use crate::shared_bam::SharedBamRecord;
use gatk_common::{GatkError, GatkResult};
use gatk_core::reference::{ReferenceWindowCache, SequenceDictionary};
use rust_htslib::bam;

/// `AssemblyBasedCallerUtils.REFERENCE_PADDING_FOR_ASSEMBLY`.
pub const GATK_REFERENCE_PADDING_FOR_ASSEMBLY: u64 = 500;

/// Padded reference bytes for `assembleReads` / `runLocalAssembly` (±500 bp around extended span).
pub fn assembly_reference_read(
    dictionary: &SequenceDictionary,
    ref_cache: &mut ReferenceWindowCache,
    region: &AssemblyRegion,
) -> GatkResult<AssemblyRead> {
    let contig_len = dictionary
        .contig(&region.contig)
        .map(|c| c.length)
        .ok_or_else(|| GatkError::argument(format!("unknown contig {}", region.contig)))?;
    let pad_start = region
        .extended_start
        .get()
        .saturating_sub(GATK_REFERENCE_PADDING_FOR_ASSEMBLY)
        .max(1);
    let pad_end = (region.extended_end.get() + GATK_REFERENCE_PADDING_FOR_ASSEMBLY).min(contig_len);
    let bases: Vec<u8> = if pad_start > pad_end {
        Vec::new()
    } else {
        ref_cache
            .get_interval_bytes(dictionary, &region.contig, pad_start, pad_end)?
            .to_vec()
    };
    let n = bases.len();
    Ok(AssemblyRead {
        bases,
        base_quals: vec![30; n],
    })
}

/// GATK `minBaseQualityScore - 1` tail clip threshold in `finalizeRegion`.
pub fn gatk_min_tail_quality_for_assembly(min_base_quality_score: u8) -> u8 {
    min_base_quality_score.saturating_sub(1)
}

/// GATK `AssemblyBasedCallerUtils.finalizeRegion` (soft-clip, low-qual tails, adaptor, padded-span clip, overlap).
pub fn finalize_region_reads_for_assembly(
    reads: &[SharedBamRecord],
    region: &AssemblyRegion,
    correct_overlapping_base_qualities: bool,
    min_tail_quality: u8,
    soft_clip_low_quality_ends: bool,
) -> Vec<bam::Record> {
    let policy = HcSoftclipPolicy::haplotype_caller_defaults();
    let ref_start = region.extended_start.get();
    let ref_stop = region.extended_end.get();
    let mut out: Vec<bam::Record> = Vec::new();
    for original in reads {
        let mut read = apply_hc_softclip_pre_step(original, &policy).0;
        if read.tid() < 0 || read.is_unmapped() || read.pos() < 0 {
            continue;
        }
        if unclipped_read_length(&read) == 0 {
            continue;
        }
        read = if soft_clip_low_quality_ends {
            soft_clip_low_qual_ends(&read, min_tail_quality)
        } else {
            hard_clip_low_qual_ends(&read, min_tail_quality)
        };
        if read.is_unmapped() || !read_has_positive_cigar_length(&read) {
            continue;
        }
        let adaptor_clipped = hard_clip_adaptor_sequence(&read);
        if adaptor_clipped.is_unmapped() || !read_has_positive_cigar_length(&adaptor_clipped) {
            continue;
        }
        if !read_overlaps_padded_span(&adaptor_clipped, ref_start, ref_stop) {
            continue;
        }
        read = hard_clip_to_region(&adaptor_clipped, ref_start, ref_stop);
        if read.is_unmapped() || !read_has_positive_cigar_length(&read) {
            continue;
        }
        let aln_start = read.pos() + 1;
        let aln_end = i64::from(alignment_end_1based(&read));
        if aln_start > aln_end || aln_start > ref_stop as i64 || aln_end < ref_start as i64 {
            continue;
        }
        normalize_record_cigar(&mut read);
        out.push(read);
    }
    out.sort_by(|a, b| {
        a.tid()
            .cmp(&b.tid())
            .then_with(|| a.pos().cmp(&b.pos()))
            .then_with(|| a.qname().cmp(b.qname()))
    });
    if correct_overlapping_base_qualities {
        // INVARIANT: `out` was sorted by (tid, pos, qname) immediately above.
        #[allow(clippy::expect_used)]
        clean_overlapping_read_pairs(&mut out, true)
            .expect("reads sorted by alignment start before overlap correction");
    }
    out
}

/// Re-clip already-finalized records to a (possibly trimmed) region's padded span.
///
/// Softclip / adaptor / overlap-quality correction are **not** re-applied — those ran in
/// [`finalize_region_reads_for_assembly`]. Only hard-clip-to-region + overlap/unmapped filters + sort.
pub fn clip_finalized_reads_to_region(
    reads: &[bam::Record],
    region: &AssemblyRegion,
) -> Vec<bam::Record> {
    let mut owned = reads.to_vec();
    clip_finalized_reads_in_place(&mut owned, region);
    owned
}

/// Consume-and-clip path for the assemble `finalizeRegion` buffer (no second full copy when
/// most reads already lie inside the genotyping/padded span).
pub fn clip_finalized_reads_in_place(reads: &mut Vec<bam::Record>, region: &AssemblyRegion) {
    let ref_start = region.extended_start.get();
    let ref_stop = region.extended_end.get();
    reads.retain_mut(|original| {
        if original.tid() < 0 || original.is_unmapped() || original.pos() < 0 {
            return false;
        }
        if !read_overlaps_padded_span(original, ref_start, ref_stop) {
            return false;
        }
        let aln_start = original.pos() + 1;
        let aln_end = i64::from(alignment_end_1based(original));
        let needs_clip = aln_start < ref_start as i64 || aln_end > ref_stop as i64;
        if needs_clip {
            let clipped = hard_clip_to_region(original, ref_start, ref_stop);
            if clipped.is_unmapped() || !read_has_positive_cigar_length(&clipped) {
                return false;
            }
            let aln_start = clipped.pos() + 1;
            let aln_end = i64::from(alignment_end_1based(&clipped));
            if aln_start > aln_end || aln_start > ref_stop as i64 || aln_end < ref_start as i64 {
                return false;
            }
            *original = clipped;
        }
        normalize_record_cigar(original);
        true
    });
    reads.sort_by(|a, b| {
        a.tid()
            .cmp(&b.tid())
            .then_with(|| a.pos().cmp(&b.pos()))
            .then_with(|| a.qname().cmp(b.qname()))
    });
}

fn read_overlaps_padded_span(rec: &bam::Record, ref_start: u64, ref_stop: u64) -> bool {
    if rec.tid() < 0 || rec.is_unmapped() || rec.pos() < 0 {
        return false;
    }
    let aln_start = rec.pos() + 1;
    let aln_end = i64::from(alignment_end_1based(rec));
    aln_start <= ref_stop as i64 && aln_end >= ref_start as i64
}

fn alignment_end_1based(rec: &bam::Record) -> i32 {
    crate::read_unclip::alignment_end_1based(rec)
}

/// GATK `getPaddedReferenceLoc` span (±500 bp around the assembly region padded span).
pub fn padded_reference_loc(
    region: &AssemblyRegion,
    dictionary: &SequenceDictionary,
) -> (u64, u64) {
    let contig_len = dictionary
        .contig(&region.contig)
        .map(|c| c.length)
        .unwrap_or(u64::MAX);
    let left = region
        .extended_start
        .get()
        .saturating_sub(GATK_REFERENCE_PADDING_FOR_ASSEMBLY)
        .max(1);
    let right = (region.extended_end.get() + GATK_REFERENCE_PADDING_FOR_ASSEMBLY).min(contig_len);
    (left, right)
}

/// Reference hap aligned to the assembly region padded span (GATK `createReferenceHaplotype`).
pub fn reference_haplotype_for_assembly_region(
    reference: &AssemblyRead,
    region: &AssemblyRegion,
    dictionary: &SequenceDictionary,
) -> Haplotype {
    let (loc_start, _) = padded_reference_loc(region, dictionary);
    let alignment_start = region.extended_start.get().saturating_sub(loc_start) as usize;
    let span_len = (region.extended_end.get() - region.extended_start.get() + 1) as usize;
    let ref_bytes = reference.bases.as_slice();
    let end = (alignment_start + span_len).min(ref_bytes.len());
    let bases = ref_bytes[alignment_start..end].to_vec();
    // CLONE: needed because haplotype constructor takes owned bases.
    let mut h = Haplotype::new(bases.clone(), true);
    h.alignment_start_hap_wrt_ref = alignment_start;
    let mut cigar = Cigar::new();
    cigar.push(bases.len(), CigarOperator::Match);
    h.cigar = Some(cigar);
    h
}

/// Alias for parity gate dumps (same as production slice).
#[cfg(any(feature = "dev-dumps", test))]
pub fn materialize_reference_haplotype_for_dump(
    reference: &AssemblyRead,
    region: &AssemblyRegion,
    dictionary: &SequenceDictionary,
) -> Haplotype {
    reference_haplotype_for_assembly_region(reference, region, dictionary)
}

/// Production `assembleReads` read set (`finalizeRegion` + optional error correction).
pub fn assembly_reads_for_production(
    records: &[SharedBamRecord],
    region: &AssemblyRegion,
    min_tail_quality: u8,
    correct_overlapping_base_qualities: bool,
    soft_clip_low_quality_ends: bool,
) -> Vec<bam::Record> {
    finalize_region_reads_for_assembly(
        records,
        region,
        correct_overlapping_base_qualities,
        min_tail_quality,
        soft_clip_low_quality_ends,
    )
}

/// Reads for ASM-1 Java `RegionAssemblyMaterial` parity (`HcFullParityGateDump` uses iterator reads, not `finalizeRegion`).
#[cfg(any(feature = "dev-dumps", test))]
pub fn assembly_reads_for_java_materialize_dump(records: &[SharedBamRecord]) -> Vec<AssemblyRead> {
    records_to_assembly_reads(records)
}

pub fn records_to_assembly_reads<R: std::borrow::Borrow<bam::Record>>(
    records: &[R],
) -> Vec<AssemblyRead> {
    records
        .iter()
        .map(|rec| {
            let rec = rec.borrow();
            let bases: Vec<u8> = rec.seq().as_bytes().to_vec();
            AssemblyRead {
                bases,
                base_quals: rec.qual().to_vec(),
            }
        })
        .collect()
}
