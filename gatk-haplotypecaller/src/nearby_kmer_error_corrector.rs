//! GATK `NearbyKmerErrorCorrector` — k-mer spectrum read correction before assembly.

use gatk_common::GatkResult;
use rust_htslib::bam::Record;
use std::collections::{BTreeMap, HashMap};

/// GATK `AssemblerArgumentCollection` defaults for k-mer read error correction.
pub const GATK_KMER_LENGTH_FOR_READ_ERROR_CORRECTION: usize = 25;
pub const GATK_MIN_OBSERVATIONS_FOR_KMER_TO_BE_SOLID: usize = 20;
/// GATK `HaplotypeCallerEngine.MIN_TAIL_QUALITY_WITH_ERROR_CORRECTION`.
pub const GATK_MIN_TAIL_QUALITY_WITH_ERROR_CORRECTION: u8 = 6;

/// Configuration matching Java `NearbyKmerErrorCorrector` simple constructor.
/// # Invariants
/// `kmer_length` must be positive for k-mer counting; reads shorter than this are skipped.
/// `min_observations_for_kmer_to_be_solid` gates whether a k-mer is trusted for correction.
/// # Ownership
/// [`Clone`] config passed into [`correct_reads_nearby_kmer`]; does not own reads or reference.
/// # Mutation
/// Immutable per correction pass; read bases are mutated in the caller's slice.
/// # Biological assumptions
/// Illumina-like errors cluster in read tails; solid k-mers come from local read pileup agreement.
/// # Java equivalence
/// GATK `NearbyKmerErrorCorrector` simple constructor + `AssemblerArgumentCollection` k-mer defaults.
#[derive(Debug, Clone)]
pub struct NearbyKmerErrorCorrectorConfig {
    pub kmer_length: usize,
    pub min_tail_quality: u8,
    pub min_observations_for_kmer_to_be_solid: usize,
    pub debug: bool,
}

impl NearbyKmerErrorCorrectorConfig {
    pub fn gatk_defaults() -> Self {
        Self {
            kmer_length: GATK_KMER_LENGTH_FOR_READ_ERROR_CORRECTION,
            min_tail_quality: GATK_MIN_TAIL_QUALITY_WITH_ERROR_CORRECTION,
            min_observations_for_kmer_to_be_solid: GATK_MIN_OBSERVATIONS_FOR_KMER_TO_BE_SOLID,
            debug: false,
        }
    }
}

struct KmerCounter {
    kmer_length: usize,
    /// Ordered map so solid-kmer ties resolve by lexicographic sequence.
    counts: BTreeMap<Vec<u8>, u32>,
}

impl KmerCounter {
    fn new(kmer_length: usize) -> Self {
        Self {
            kmer_length,
            counts: BTreeMap::new(),
        }
    }

    fn add_read(&mut self, bases: &[u8]) {
        if bases.len() < self.kmer_length {
            return;
        }
        for offset in 0..=bases.len() - self.kmer_length {
            let kmer = bases[offset..offset + self.kmer_length].to_vec();
            *self.counts.entry(kmer).or_insert(0) += 1;
        }
    }
}

fn max_homopolymer_run(bases: &[u8]) -> usize {
    let mut best = 1usize;
    let mut run = 1usize;
    let mut prev = bases.first().copied();
    for &b in bases.iter().skip(1) {
        if Some(b) == prev {
            run += 1;
            best = best.max(run);
        } else {
            run = 1;
            prev = Some(b);
        }
    }
    best
}

fn hamming_distance(a: &[u8], b: &[u8]) -> usize {
    a.iter().zip(b.iter()).filter(|(x, y)| x != y).count()
}

/// Correct reads when `errorCorrectReads` is enabled (GATK `NearbyKmerErrorCorrector.correctReads`).
pub fn correct_reads_nearby_kmer(
    reads: &mut [Record],
    full_reference_with_padding: &[u8],
    cfg: &NearbyKmerErrorCorrectorConfig,
) -> GatkResult<()> {
    if max_homopolymer_run(full_reference_with_padding) > 12 {
        return Ok(());
    }
    let mut counter = KmerCounter::new(cfg.kmer_length);
    for rec in reads.iter() {
        let bases = rec.seq().as_bytes();
        counter.add_read(&bases);
    }
    let solid_threshold = cfg.min_observations_for_kmer_to_be_solid as u32;

    let mut insolid_to_solid: HashMap<Vec<u8>, Vec<u8>> = HashMap::new();
    for (kmer, &count) in &counter.counts {
        if count >= solid_threshold {
            continue;
        }
        let mut best: Option<(&[u8], u32)> = None;
        for (solid, &sc) in &counter.counts {
            if sc < solid_threshold {
                continue;
            }
            let dist = hamming_distance(kmer, solid);
            if dist > 0 && dist <= 2 {
                let better = match best {
                    None => true,
                    Some((prev, c)) => sc > c || (sc == c && solid.as_slice() < prev),
                };
                if better {
                    best = Some((solid.as_slice(), sc));
                }
            }
        }
        if let Some((solid, _)) = best {
            // CLONE: needed because owned HashMap/BTree/HashSet key or value.
            insolid_to_solid.insert(kmer.clone(), solid.to_vec());
        }
    }

    let mut reads_corrected = 0usize;
    for rec in reads.iter_mut() {
        let mut bases = rec.seq().as_bytes().to_vec();
        let mut quals = rec.qual().to_vec();
        let mut corrections: Vec<(usize, u8)> = Vec::new();
        for offset in 0..bases.len().saturating_sub(cfg.kmer_length) {
            let kmer = bases[offset..offset + cfg.kmer_length].to_vec();
            let Some(solid) = insolid_to_solid.get(&kmer) else {
                continue;
            };
            for (i, (&ob, &sb)) in kmer.iter().zip(solid.iter()).enumerate() {
                if ob != sb {
                    corrections.push((offset + i, sb));
                }
            }
        }
        if corrections.is_empty() {
            continue;
        }
        for (pos, base) in corrections {
            bases[pos] = base;
            quals[pos] = 30;
        }
        let qname = rec.qname().to_vec();
        let cigar = rec.cigar();
        rec.set(&qname, Some(&cigar), &bases, &quals);
        reads_corrected += 1;
    }
    if cfg.debug {
        eprintln!(
            "[NearbyKmerErrorCorrector] corrected {reads_corrected}/{} reads ({} kmer mappings)",
            reads.len(),
            insolid_to_solid.len()
        );
    }
    Ok(())
}
