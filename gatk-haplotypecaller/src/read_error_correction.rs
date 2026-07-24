//! GATK `PileupReadErrorCorrector` parity (pileup log-odds correction before assembly).

use crate::pileup_element::pileup_element_flags_at_ref;
use gatk_common::{GatkError, GatkResult};
use std::collections::HashMap;

const GOOD_QUAL: u8 = 30;
const INDEL_SPAN: usize = 15;
const INDEL_MISMATCHES: usize = 3;

/// One read aligned to a common reference interval (1-based `start`).
/// # Invariants
/// `bases` and `base_quals` have equal length when constructed from valid BAM records.
/// `start1` is 1-based alignment start on the shared reference interval.
/// # Ownership
/// Owns base and qual vectors; pileup correction mutates them in place.
/// # Mutation
/// [`correct_reads_pileup_log_odds`] may change bases/quals per locus log-odds threshold.
/// # Biological assumptions
/// Log-odds pileup correction targets isolated base errors before local assembly.
/// # Java equivalence
/// GATK `PileupReadErrorCorrector` aligned-read representation (pileup log-odds path).
#[derive(Debug, Clone)]
pub struct AlignedAssemblyRead {
    pub bases: Vec<u8>,
    pub base_quals: Vec<u8>,
    pub start1: u64,
}

fn qual_to_error_prob(qual: u8) -> f64 {
    10_f64.powf(-(qual as f64) / 10.0)
}

fn qual_to_log_error_prob(qual: u8) -> f64 {
    qual_to_error_prob(qual).ln()
}

fn qual_to_log_prob(qual: u8) -> f64 {
    (1.0 - qual_to_error_prob(qual)).ln()
}

fn fast_bernoulli_entropy(z: f64) -> f64 {
    if z <= 0.0 || z >= 1.0 {
        return 0.0;
    }
    let o = 1.0 - z;
    -z * z.ln() - o * o.ln()
}

fn log_binomial(n: usize, k: usize) -> f64 {
    crate::activity_scoring::log_binomial_coefficient_natural(n as u32, k as u32)
}

/// GATK `Mutect2Engine.logLikelihoodRatio(nRef, altQuals, repeatFactor)`.
pub fn mutect_log_likelihood_ratio_alt_quals(
    n_ref: usize,
    alt_quals: &[u8],
    repeat_factor: usize,
) -> f64 {
    let n_alt = repeat_factor * alt_quals.len();
    if n_alt == 0 {
        return 0.0;
    }
    let n = n_ref + n_alt;
    let digamma = statrs::function::gamma::digamma;
    let f_tilde_ratio = (digamma((n_ref + 1) as f64) - digamma((n_alt + 1) as f64)).exp();
    let mut read_sum = 0.0;
    for &qual in alt_quals {
        let epsilon = qual_to_error_prob(qual);
        let z_bar_alt = (1.0 - epsilon) / (1.0 - epsilon + epsilon * f_tilde_ratio);
        let log_epsilon = qual_to_log_error_prob(qual);
        let log_one_minus_epsilon = qual_to_log_prob(qual);
        read_sum +=
            z_bar_alt * (log_one_minus_epsilon - log_epsilon) + fast_bernoulli_entropy(z_bar_alt);
    }
    let beta_entropy = -(n as f64 + 1.0).ln() - log_binomial(n, n_alt);
    beta_entropy + read_sum * repeat_factor as f64
}

fn skip_pileup_element(flags: &crate::pileup_element::PileupElementFlags) -> bool {
    flags.is_deletion
        || flags.is_before_deletion_start
        || flags.is_after_deletion_end
        || flags.is_before_insertion
        || flags.is_after_insertion
        || flags.is_next_to_soft_clip
}

/// Correct reads in-place (GATK `PileupReadErrorCorrector.correctReads`).
pub fn correct_reads_pileup_log_odds(
    reads: &mut [AlignedAssemblyRead],
    log_odds_threshold: f64,
) -> GatkResult<()> {
    if reads.is_empty() {
        return Ok(());
    }
    let min_start = reads.iter().map(|r| r.start1).min().unwrap();
    let max_end = reads
        .iter()
        .map(|r| r.start1 + r.bases.len() as u64)
        .max()
        .unwrap();

    let mut potential: HashMap<usize, Vec<(usize, u8)>> = HashMap::new();
    for (read_idx, _read) in reads.iter().enumerate() {
        potential.insert(read_idx, Vec::new());
    }

    for pos1 in min_start..max_end {
        let ref_pos0 = pos1.saturating_sub(1) as i64;
        let mut counter = [0u64; 4];
        let mut elements: Vec<(
            usize,
            usize,
            u8,
            u8,
            crate::pileup_element::PileupElementFlags,
        )> = Vec::new();

        for (read_idx, read) in reads.iter().enumerate() {
            let offset = pos1.saturating_sub(read.start1) as usize;
            if offset >= read.bases.len() {
                continue;
            }
            let alignment_start = read.start1.saturating_sub(1) as i64;
            let cigar = [rust_htslib::bam::record::Cigar::Match(
                read.bases.len() as u32
            )];
            let Some(flags) = pileup_element_flags_at_ref(
                alignment_start,
                &cigar,
                &read.bases,
                &read.base_quals,
                ref_pos0,
            ) else {
                continue;
            };
            if skip_pileup_element(&flags) {
                continue;
            }
            let base = flags.read_base;
            match base {
                b'A' => counter[0] += 1,
                b'C' => counter[1] += 1,
                b'G' => counter[2] += 1,
                b'T' => counter[3] += 1,
                _ => {}
            }
            elements.push((read_idx, offset, base, flags.qual, flags));
        }

        let plurality = [
            (b'A', counter[0]),
            (b'C', counter[1]),
            (b'G', counter[2]),
            (b'T', counter[3]),
        ]
        .into_iter()
        .max_by_key(|(_, c)| *c)
        .map(|(b, _)| b);
        let Some(ref_base) = plurality else {
            continue;
        };

        let mut ref_count = 0usize;
        let mut alt_quals = Vec::new();
        for (_, _, base, qual, _) in &elements {
            if *base == ref_base {
                ref_count += 1;
            } else {
                alt_quals.push(*qual);
            }
        }
        if alt_quals.is_empty() {
            continue;
        }
        let log_odds = mutect_log_likelihood_ratio_alt_quals(ref_count, &alt_quals, 1);
        if log_odds >= log_odds_threshold {
            continue;
        }
        for (read_idx, offset, base, _, flags) in elements {
            if base != ref_base && !skip_pileup_element(&flags) {
                potential
                    .get_mut(&read_idx)
                    .expect("read_idx")
                    .push((offset, ref_base));
            }
        }
    }

    for (read_idx, read) in reads.iter_mut().enumerate() {
        let mut edits = potential.remove(&read_idx).unwrap_or_default();
        edits.sort_by_key(|(off, _)| *off);
        let size = edits.len();
        if size == 0 {
            continue;
        }
        let mut first_edit = 0usize;
        for n in 0..size.saturating_sub(INDEL_MISMATCHES) {
            if n + INDEL_MISMATCHES < size
                && edits[n + INDEL_MISMATCHES - 1].0.saturating_sub(edits[n].0) < INDEL_SPAN
            {
                first_edit = n + INDEL_MISMATCHES;
            }
        }
        let mut last_edit = size.saturating_sub(1);
        for n in (INDEL_MISMATCHES.saturating_sub(1)..size).rev() {
            if n >= INDEL_MISMATCHES - 1
                && edits[n]
                    .0
                    .saturating_sub(edits[n - (INDEL_MISMATCHES - 1)].0)
                    < INDEL_SPAN
            {
                last_edit = n.saturating_sub(INDEL_MISMATCHES);
            }
        }
        for n in first_edit..=last_edit.min(size.saturating_sub(1)) {
            let (offset, base) = edits[n];
            if offset < read.bases.len() {
                read.bases[offset] = base;
                read.base_quals[offset] = GOOD_QUAL;
            }
        }
    }
    Ok(())
}

/// Load `sequence qual start1` rows for read-correction parity dumps.
pub fn load_aligned_assembly_reads_tsv(
    path: &std::path::Path,
) -> GatkResult<Vec<AlignedAssemblyRead>> {
    use std::fs::File;
    use std::io::{BufRead, BufReader};
    let f = File::open(path)
        .map_err(|e| GatkError::generic(format!("open {}: {e}", path.display())))?;
    let reader = BufReader::new(f);
    let mut out = Vec::new();
    for line in reader.lines() {
        let line = line.map_err(|e| GatkError::generic(format!("read {}: {e}", path.display())))?;
        let t = line.trim();
        if t.is_empty() || t.starts_with('#') {
            continue;
        }
        let mut parts = t.split_whitespace();
        let bases: Vec<u8> = parts
            .next()
            .ok_or_else(|| GatkError::argument("reads tsv: missing sequence"))?
            .as_bytes()
            .to_vec();
        let q = parts
            .next()
            .ok_or_else(|| GatkError::argument("reads tsv: missing qual"))?
            .parse::<u8>()
            .map_err(|_| GatkError::argument("reads tsv: invalid qual"))?;
        let start1: u64 = parts
            .next()
            .ok_or_else(|| GatkError::argument("reads tsv: missing start1"))?
            .parse()
            .map_err(|_| GatkError::argument("reads tsv: invalid start1"))?;
        let n = bases.len();
        out.push(AlignedAssemblyRead {
            bases,
            base_quals: vec![q; n],
            start1,
        });
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pileup_correction_fixes_singleton_errors() {
        let consensus = b"ACGTACGTACGTACGTACGTACGTACGTACGT";
        let err_pos = 15usize;
        let mut reads = Vec::new();
        for _ in 0..20 {
            reads.push(AlignedAssemblyRead {
                bases: consensus.to_vec(),
                base_quals: vec![30; consensus.len()],
                start1: 1,
            });
        }
        for _ in 0..2 {
            let mut bases = consensus.to_vec();
            bases[err_pos] = b'T';
            reads.push(AlignedAssemblyRead {
                bases,
                base_quals: vec![30; consensus.len()],
                start1: 1,
            });
        }
        correct_reads_pileup_log_odds(&mut reads, 3.0).unwrap();
        for read in &reads {
            assert_eq!(read.bases[err_pos], consensus[err_pos]);
        }
    }
}
