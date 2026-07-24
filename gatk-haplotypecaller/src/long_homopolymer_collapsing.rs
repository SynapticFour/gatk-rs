//! GATK `LongHomopolymerHaplotypeCollapsingEngine` (flow-based HC).
//! Standard Illumina `HaplotypeCaller` uses `flowAssemblyCollapseHKerSize <= 0` and skips this engine.
//! When enabled, Java uncollapses long homopolymer haplotypes against the padded reference.

use crate::cigar::CigarOperator;
use crate::genome_loc::GenomeLoc;
use crate::haplotype::Haplotype;
use crate::smith_waterman::{align, SwOverhangStrategy, SwParameters};
use std::collections::BTreeMap;

/// GATK `LongHomopolymerHaplotypeCollapsingEngine.needsCollapsing`.
pub fn needs_collapsing(ref_bases: &[u8], hmer_size_threshold: usize) -> bool {
    if hmer_size_threshold == 0 {
        return false;
    }
    let mut last_base = 0u8;
    let mut base_same_count = 0usize;
    for &base in ref_bases {
        if base == last_base {
            base_same_count += 1;
            if base_same_count >= hmer_size_threshold {
                return true;
            }
        } else {
            last_base = base;
            base_same_count = 0;
        }
    }
    false
}

/// GATK `LongHomopolymerHaplotypeCollapsingEngine.identicalBySequence`.
pub fn identical_by_sequence(haplotypes: &[Haplotype]) -> BTreeMap<usize, Vec<usize>> {
    let mut by_seq: BTreeMap<String, Vec<usize>> = BTreeMap::new();
    for (i, h) in haplotypes.iter().enumerate() {
        by_seq.entry(h.sequence_string()).or_default().push(i);
    }
    let ref_idx = haplotypes
        .iter()
        .position(|h| h.is_reference)
        .expect("reference haplotype");
    // Lifetime: consume `by_seq` so index vectors move into `out` without cloning.
    let mut out: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for (_, indices) in by_seq {
        let key = if indices.contains(&ref_idx) {
            ref_idx
        } else {
            indices[0]
        };
        out.insert(key, indices);
    }
    out
}

fn complement_base(b: u8) -> u8 {
    match b.to_ascii_uppercase() {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        _ => b'N',
    }
}

fn reverse_complement(bases: &mut [u8]) {
    let n = bases.len();
    for i in 0..n / 2 {
        let j = n - 1 - i;
        let (left, right) = (complement_base(bases[j]), complement_base(bases[i]));
        bases[i] = left;
        bases[j] = right;
    }
    if n % 2 == 1 {
        let m = n / 2;
        bases[m] = complement_base(bases[m]);
    }
}

fn collapse_bases(full_bases: &[u8], hmer_size_threshold: usize) -> Vec<u8> {
    let mut collapsed = Vec::with_capacity(full_bases.len());
    let mut last_base = 0u8;
    let mut base_same_count = 0usize;
    let mut first_homopolymer = true;
    for &base in full_bases {
        if base == last_base {
            base_same_count += 1;
            if first_homopolymer || base_same_count < hmer_size_threshold {
                collapsed.push(base);
            }
        } else {
            if last_base != 0 {
                first_homopolymer = false;
            }
            last_base = base;
            base_same_count = 0;
            collapsed.push(base);
        }
    }
    collapsed
}

struct UncollapseResult {
    bases: Vec<u8>,
    offset: usize,
}

fn on_homopolymer(
    bases: &[u8],
    ofs: isize,
    base: u8,
    length: usize,
    hmer_size_threshold: usize,
) -> bool {
    for tick in 0..hmer_size_threshold {
        if same_base(bases, ofs + tick as isize, base, length) {
            return true;
        }
    }
    false
}

fn same_base(bases: &[u8], ofs: isize, base: u8, length: usize) -> bool {
    let mut ofs = ofs;
    let mut length = length;
    if ofs < 0 {
        return false;
    }
    if (ofs as usize).saturating_add(length) > bases.len() {
        return false;
    }
    while length != 0 {
        if bases[ofs as usize] != base {
            return false;
        }
        ofs += 1;
        length -= 1;
    }
    true
}

fn uncollapse_by_ref(
    bases_arg: &[u8],
    ref_arg: &[u8],
    rev: bool,
    hmer_size_threshold: usize,
    partial_mode: bool,
) -> UncollapseResult {
    let mut bases = bases_arg.to_vec();
    let mut reference = ref_arg.to_vec();
    if rev {
        reverse_complement(&mut bases);
        reverse_complement(&mut reference);
    }

    let sw_params = SwParameters::gatk_haplotype_to_reference();
    let Ok(alignment) = align(&reference, &bases, &sw_params, SwOverhangStrategy::Indel) else {
        return UncollapseResult {
            bases: bases_arg.to_vec(),
            offset: 0,
        };
    };

    let mut result_length = bases.len();
    for elem in &alignment.cigar.elements {
        if elem.operator == CigarOperator::Deletion {
            result_length += elem.length;
        }
    }
    let mut result = vec![0u8; result_length];
    let mut bases_ofs = alignment.alignment_offset.max(0) as usize;
    let mut ref_ofs = 0usize;
    let mut result_ofs = 0usize;
    for elem in &alignment.cigar.elements {
        if elem.operator != CigarOperator::Deletion {
            if elem.operator.consumes_read_bases() {
                let len = elem.length.min(bases.len().saturating_sub(bases_ofs));
                result[result_ofs..result_ofs + len]
                    .copy_from_slice(&bases[bases_ofs..bases_ofs + len]);
                bases_ofs += len;
                result_ofs += len;
            }
        } else {
            let fwd_end = (bases_ofs + hmer_size_threshold).min(bases.len());
            let bck_start = bases_ofs.saturating_sub(hmer_size_threshold);
            let fwd_slice = &bases[bases_ofs..fwd_end];
            let bck_slice = &bases[bck_start..bases_ofs];
            if needs_collapsing(fwd_slice, hmer_size_threshold.saturating_sub(1))
                || needs_collapsing(bck_slice, hmer_size_threshold.saturating_sub(1))
            {
                if on_homopolymer(
                    &reference,
                    ref_ofs as isize - hmer_size_threshold as isize,
                    reference[ref_ofs],
                    hmer_size_threshold,
                    hmer_size_threshold,
                ) {
                    let base = reference[ref_ofs];
                    for size in 0..elem.length {
                        if partial_mode && reference.get(ref_ofs + size) != Some(&base) {
                            break;
                        }
                        if result_ofs < result.len() {
                            result[result_ofs] = base;
                            result_ofs += 1;
                        }
                    }
                } else if on_homopolymer(
                    &reference,
                    (ref_ofs + elem.length) as isize,
                    reference[ref_ofs + elem.length.saturating_sub(1)],
                    hmer_size_threshold,
                    hmer_size_threshold,
                ) {
                    let base = reference[ref_ofs + elem.length.saturating_sub(1)];
                    for size in 0..elem.length {
                        let ri = ref_ofs + elem.length - 1 - size;
                        if partial_mode && reference.get(ri) != Some(&base) {
                            break;
                        }
                        if result_ofs < result.len() {
                            result[result_ofs] = base;
                            result_ofs += 1;
                        }
                    }
                }
            }
        }
        if elem.operator.consumes_reference_bases() {
            ref_ofs += elem.length;
        }
    }

    let mut final_result = if result_ofs == result.len() {
        result
    } else {
        result[..result_ofs].to_vec()
    };
    if rev {
        reverse_complement(&mut final_result);
    }
    UncollapseResult {
        bases: final_result,
        offset: alignment.alignment_offset.max(0) as usize,
    }
}

fn uncollapsed_partial_ref(full_ref: &[u8], ref_loc_start: u64, uc_loc: &GenomeLoc) -> Vec<u8> {
    let uc_ofs = uc_loc.start_1based().saturating_sub(ref_loc_start) as usize;
    let size = uc_loc.reference_span_length() as usize;
    full_ref
        .get(uc_ofs..uc_ofs.saturating_add(size))
        .unwrap_or(&[])
        .to_vec()
}

fn uncollapse_single_haplotype(
    h: &Haplotype,
    limit_to_threshold: bool,
    full_ref: &[u8],
    ref_loc_start: u64,
    hmer_size_threshold: usize,
    partial_mode: bool,
    ref_map: &mut BTreeMap<u64, Vec<u8>>,
) -> Haplotype {
    if h.is_reference {
        return h.clone();
    }
    let loc = h
        .genome_loc
        .unwrap_or_else(|| GenomeLoc::new(ref_loc_start, ref_loc_start));
    let ref_bases = ref_map
        .entry(loc.start_1based())
        .or_insert_with(|| uncollapsed_partial_ref(full_ref, ref_loc_start, &loc))
        .clone();

    let fwd = uncollapse_by_ref(
        &h.bases,
        &ref_bases,
        false,
        hmer_size_threshold,
        partial_mode,
    );
    let rev = uncollapse_by_ref(
        &h.bases,
        &ref_bases,
        true,
        hmer_size_threshold,
        partial_mode,
    );
    let picked = if rev.bases.len() > fwd.bases.len() {
        rev
    } else {
        fwd
    };
    let mut bases = picked.bases;
    if limit_to_threshold {
        bases = collapse_bases(&bases, hmer_size_threshold);
    }
    let mut out = Haplotype::new(bases, h.is_reference);
    out.score = h.score;
    out.genome_loc = h.genome_loc;
    out.alignment_start_hap_wrt_ref = picked.offset;
    // CLONE: needed because multi-owner or ownership transfer into new structure.
    out.cigar = h.cigar.clone();
    out
}

/// GATK `uncollapseHmersInHaplotypes`.
pub fn uncollapse_hmers_in_haplotypes(
    haplotypes: &[Haplotype],
    limit_to_hmer_size_threshold: bool,
    full_ref: &[u8],
    ref_loc_start_1based: u64,
    hmer_size_threshold: usize,
    partial_mode: bool,
) -> Vec<Haplotype> {
    let sw_params = SwParameters::gatk_haplotype_to_reference();
    let mut ref_map: BTreeMap<u64, Vec<u8>> = BTreeMap::new();
    let mut ref_bases: Option<&[u8]> = None;
    let mut alignment_start = 0usize;

    for h in haplotypes {
        if h.is_reference {
            ref_bases = Some(&h.bases);
            alignment_start = h.alignment_start_hap_wrt_ref;
        }
    }

    let mut result: Vec<Haplotype> = Vec::with_capacity(haplotypes.len());
    for h in haplotypes {
        let unc = uncollapse_single_haplotype(
            h,
            limit_to_hmer_size_threshold,
            full_ref,
            ref_loc_start_1based,
            hmer_size_threshold,
            partial_mode,
            &mut ref_map,
        );
        result.push(unc);
    }

    if let Some(ref_bytes) = ref_bases {
        for h in &mut result {
            if h.is_reference {
                continue;
            }
            if let Ok(alignment) = align(ref_bytes, &h.bases, &sw_params, SwOverhangStrategy::Indel)
            {
                h.cigar = Some(alignment.cigar);
                h.alignment_start_hap_wrt_ref =
                    alignment.alignment_offset.max(0) as usize + alignment_start;
            }
        }
    }
    result
}

/// Apply flow homopolymer uncollapse when configured (GATK `AssemblyBasedCallerUtils` path).
pub fn collapse_haplotypes_if_configured(
    haplotypes: Vec<Haplotype>,
    flow_assembly_collapse_hmer_size: usize,
    partial_mode: bool,
    full_ref: &[u8],
    ref_loc_start_1based: u64,
) -> Vec<Haplotype> {
    if flow_assembly_collapse_hmer_size == 0 {
        return haplotypes;
    }
    let ref_hap = haplotypes.iter().find(|h| h.is_reference);
    let ref_slice = ref_hap.map(|h| h.bases.as_slice()).unwrap_or(full_ref);
    if !needs_collapsing(ref_slice, flow_assembly_collapse_hmer_size) {
        return haplotypes;
    }
    uncollapse_hmers_in_haplotypes(
        &haplotypes,
        false,
        full_ref,
        ref_loc_start_1based,
        flow_assembly_collapse_hmer_size,
        partial_mode,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn needs_collapsing_matches_java_unit_test() {
        assert!(!needs_collapsing(b"CCAATTGG", 12));
        assert!(!needs_collapsing(b"CCAAAAAAAAAAAATTGG", 12));
        assert!(needs_collapsing(b"CCAAAAAAAAAAAAATTGG", 12));
    }

    #[test]
    fn identical_by_sequence_groups() {
        let haps = vec![
            Haplotype::new(b"AAAA", false),
            Haplotype::new(b"AAAA", true),
            Haplotype::new(b"CCCC", false),
            Haplotype::new(b"AAAA", false),
        ];
        let m = identical_by_sequence(&haps);
        assert_eq!(m.len(), 2);
        assert_eq!(m.get(&1).map(|v| v.len()), Some(3));
        assert_eq!(m.get(&2).map(|v| v.len()), Some(1));
    }
}
