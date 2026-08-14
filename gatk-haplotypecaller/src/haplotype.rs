//! GATK `Haplotype` representation for assembly output.
//! # Invariants (assembly result lists)
//! Exactly one haplotype should be marked `is_reference` before genotyping filters run.
//! [`sort_haplotypes_assembly_result_order`]: non-reference haplotypes by descending
//! `score` (base-sequence tie-break), **reference last**.
//! When present, `genome_loc` uses 1-based inclusive coordinates ([`GenomeLoc`]).

use crate::cigar::{Cigar, CigarElement};
use crate::cigar_builder::CigarBuilder;
use crate::genome_loc::GenomeLoc;
use crate::haplotype_cigar::{get_bases_covering_ref_interval, trim_cigar_by_reference};
use std::collections::HashMap;
use std::hash::{Hash, Hasher};

/// Assembled haplotype with optional CIGAR vs padded reference.
/// # Invariants
/// Exactly one haplotype should be marked `is_reference` before genotyping filters run.
/// Non-reference haplotypes sort by descending `score` (base-sequence tie-break), reference last.
/// When present, `genome_loc` uses 1-based inclusive coordinates ([`GenomeLoc`]).
/// # Ownership
/// Owns base bytes and optional [`Cigar`]; borrowed by EventMap and PairHMM stages.
/// # Mutation
/// Assembly/finalize paths mutate `cigar`, `genome_loc`, and `alignment_start_hap_wrt_ref`; trim produces new instances.
/// # Biological assumptions
/// Bases are haplotype sequence over the padded active region; CIGAR is vs reference haplotype.
/// # Java equivalence
/// GATK `Haplotype` (assembly result, CIGAR, `alignmentStartHapwrtRef`, trim).
#[derive(Debug, Clone, PartialEq)]
pub struct Haplotype {
    pub bases: Vec<u8>,
    pub is_reference: bool,
    pub score: f64,
    pub kmer_size: usize,
    pub cigar: Option<Cigar>,
    /// GATK `genomeLocation` (required for [`Self::trim`]).
    pub genome_loc: Option<GenomeLoc>,
    /// GATK `alignmentStartHapwrtRef` (offset of haplotype vs reference haplotype).
    pub alignment_start_hap_wrt_ref: usize,
}

impl Haplotype {
    pub fn new(bases: impl Into<Vec<u8>>, is_reference: bool) -> Self {
        Self {
            bases: bases.into(),
            is_reference,
            score: 0.0,
            kmer_size: 0,
            cigar: None,
            genome_loc: None,
            alignment_start_hap_wrt_ref: 0,
        }
    }

    pub fn sequence_string(&self) -> String {
        String::from_utf8_lossy(&self.bases).into_owned()
    }

    /// Tag all haplotypes with the padded reference span (GATK `activeRegionExtendedLocation`).
    pub fn tag_padded_reference_span(haplotypes: &mut [Haplotype], pad_start_1based: u64) {
        let ref_len = haplotypes
            .iter()
            .find(|h| h.is_reference)
            .and_then(|h| h.cigar.as_ref().map(|c| c.reference_length()))
            .or_else(|| {
                haplotypes
                    .iter()
                    .find(|h| h.is_reference)
                    .map(|h| h.bases.len())
            })
            .unwrap_or(0);
        let end = pad_start_1based
            .saturating_add(ref_len as u64)
            .saturating_sub(1);
        let loc = GenomeLoc::new(pad_start_1based, end);
        for h in haplotypes {
            h.genome_loc = Some(loc);
            if h.alignment_start_hap_wrt_ref == 0 && !h.is_reference {
                h.alignment_start_hap_wrt_ref = 0;
            }
        }
    }

    /// GATK `Haplotype.trim(loc, ignoreRefState)`.
    pub fn trim(&self, loc: &GenomeLoc, ignore_ref_state: bool) -> Option<Haplotype> {
        let genome = self.genome_loc?;
        if !genome.contains(loc) {
            return None;
        }
        let cigar = self.cigar.as_ref()?;
        let new_start = loc.start_1based().saturating_sub(genome.start_1based()) as usize;
        let new_stop = new_start + loc.reference_span_length().saturating_sub(1) as usize;
        let new_bases =
            get_bases_covering_ref_interval(new_start, new_stop, &self.bases, 0, cigar)?;
        if new_bases.is_empty() {
            return None;
        }
        let trimmed = trim_cigar_by_reference(cigar, new_start, new_stop);
        let new_cigar = trimmed.cigar;
        let leading_insertion = new_cigar
            .elements
            .first()
            .is_some_and(|e| !e.operator.consumes_reference_bases());
        let trailing_insertion = new_cigar
            .elements
            .last()
            .is_some_and(|e| !e.operator.consumes_reference_bases());
        let first_keep = if leading_insertion { 1 } else { 0 };
        let last_keep_exclusive = new_cigar
            .elements
            .len()
            .saturating_sub(if trailing_insertion { 1 } else { 0 });
        if last_keep_exclusive <= first_keep {
            return None;
        }
        let final_cigar = if leading_insertion || trailing_insertion {
            let mut b = CigarBuilder::new(false);
            for e in &new_cigar.elements[first_keep..last_keep_exclusive] {
                b.add(CigarElement {
                    length: e.length,
                    operator: e.operator,
                });
            }
            b.make_and_record().cigar
        } else {
            new_cigar
        };
        let mut ret = Haplotype::new(new_bases, !ignore_ref_state && self.is_reference);
        ret.cigar = Some(final_cigar);
        ret.genome_loc = Some(*loc);
        ret.score = self.score;
        ret.kmer_size = self.kmer_size;
        ret.alignment_start_hap_wrt_ref = new_start + self.alignment_start_hap_wrt_ref;
        Some(ret)
    }
}

/// Sort key matching GATK `Haplotype.SIZE_AND_BASE_ORDER`.
pub fn haplotype_size_and_base_order(a: &Haplotype, b: &Haplotype) -> std::cmp::Ordering {
    a.bases
        .len()
        .cmp(&b.bases.len())
        .then_with(|| a.bases.cmp(&b.bases))
}

/// Collapse base-identical haplotypes before `variationPresent` / assembly status (GATK ref path).
/// KBest can label a ref-spine path `is_reference=false` while bases match the reference haplotype;
/// Java does not treat that as a distinct alt allele.
pub fn normalize_ref_equivalent_haplotypes(haplotypes: &mut Vec<Haplotype>, ref_bytes: &[u8]) {
    for h in haplotypes.iter_mut() {
        if h.bases.as_slice() == ref_bytes {
            h.is_reference = true;
        }
    }
    // Preserve assembly discovery order (Java `getHaplotypeList`) for parity genotype indexing.
    let mut out: Vec<Haplotype> = Vec::with_capacity(haplotypes.len());
    let mut index_by_bases: std::collections::HashMap<Vec<u8>, usize> =
        std::collections::HashMap::new();
    for h in haplotypes.drain(..) {
        // Lifetime: `h` is owned for this iteration; look up by `&h.bases` so duplicates
        // never allocate a key. Clone bases only when inserting a new HashMap entry that
        // must outlive `h` after `out.push(h)` moves the haplotype.
        let is_ref = h.is_reference;
        let score = h.score;
        if let Some(&idx) = index_by_bases.get(&h.bases) {
            if is_ref {
                out[idx] = h;
            } else if !out[idx].is_reference && score > out[idx].score {
                out[idx] = h;
            }
        } else {
            let idx = out.len();
            // CLONE: needed because owned HashMap/BTree/HashSet key or value.
            index_by_bases.insert(h.bases.clone(), idx);
            out.push(h);
        }
    }
    *haplotypes = out;
}

/// Collapse alt haplotypes that share a prefix and differ only in the last 8 bp (dangling tail recovery).
pub fn collapse_dangling_tail_alt_duplicates(haplotypes: &mut Vec<Haplotype>, ref_hap: &Haplotype) {
    let n = ref_hap.bases.len();
    if n < 8 {
        return;
    }
    let tail_len = 8;
    let mut best_by_prefix: std::collections::HashMap<Vec<u8>, usize> =
        std::collections::HashMap::new();
    let mut keep = vec![true; haplotypes.len()];
    for (i, h) in haplotypes.iter().enumerate() {
        if h.is_reference || h.bases.len() != n {
            continue;
        }
        let prefix = h.bases[..n - tail_len].to_vec();
        if let Some(&prev) = best_by_prefix.get(&prefix) {
            let ref_tail = &ref_hap.bases[n - tail_len..];
            let prev_tail = &haplotypes[prev].bases[n - tail_len..];
            let cur_tail = &h.bases[n - tail_len..];
            let prefer = if prev_tail == ref_tail && cur_tail != ref_tail {
                prev
            } else if cur_tail == ref_tail && prev_tail != ref_tail {
                i
            } else if haplotypes[prev].score >= h.score {
                prev
            } else {
                i
            };
            if prefer == prev {
                keep[i] = false;
            } else {
                keep[prev] = false;
                best_by_prefix.insert(prefix, i);
            }
        } else {
            best_by_prefix.insert(prefix, i);
        }
    }
    let ref_tail = &ref_hap.bases[n - tail_len..];
    let mut out = Vec::with_capacity(haplotypes.len());
    for (i, h) in haplotypes.drain(..).enumerate() {
        if !keep[i] {
            continue;
        }
        // Dangling recovery can emit a ref-spine path with only a tail mismatch (e.g. ATGT vs ACGT).
        if !h.is_reference
            && h.bases.len() == n
            && h.bases[..n - tail_len] == ref_hap.bases[..n - tail_len]
            && h.bases[n - tail_len..] != *ref_tail
        {
            continue;
        }
        out.push(h);
    }
    *haplotypes = out;
}

/// GATK assembly list order: non-reference by descending score, reference last.
/// Invariant: after sort, if any `is_reference` haplotype exists, it is at the end;
/// preceding entries are ordered by `score` desc, then `bases` asc.
pub fn sort_haplotypes_assembly_result_order(haplotypes: &mut Vec<Haplotype>) {
    haplotypes.sort_by(|a, b| match (a.is_reference, b.is_reference) {
        (true, false) => std::cmp::Ordering::Greater,
        (false, true) => std::cmp::Ordering::Less,
        _ => b
            .score
            .partial_cmp(&a.score)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| a.bases.cmp(&b.bases)),
    });
}

/// Drop non-reference k-best fragments (GATK `findBestPaths` / full-length alt paths only).
pub fn prune_fragment_non_reference_haplotypes(
    haplotypes: &mut Vec<Haplotype>,
    ref_hap: &Haplotype,
    min_bases: usize,
) {
    let ref_len = ref_hap.bases.len();
    if ref_len == 0 {
        return;
    }
    let min_alt_len = ref_len.saturating_mul(3).saturating_div(4).max(min_bases);
    haplotypes.retain(|h| h.is_reference || h.bases.len() >= min_alt_len);
}

fn hash_hap_bases(bases: &[u8]) -> u64 {
    let mut hasher = std::collections::hash_map::DefaultHasher::new();
    bases.hash(&mut hasher);
    hasher.finish()
}

/// Sequence-identity set that does not clone haplotype bases into the key.
/// Hash collisions fall back to byte equality against already-accepted haplotypes.
pub(crate) struct HapSeqSet {
    by_hash: HashMap<(u64, bool), Vec<usize>>,
}

impl HapSeqSet {
    pub(crate) fn from_haps(haps: &[Haplotype]) -> Self {
        let mut s = Self {
            by_hash: HashMap::new(),
        };
        for (i, h) in haps.iter().enumerate() {
            s.by_hash
                .entry((hash_hap_bases(&h.bases), h.is_reference))
                .or_default()
                .push(i);
        }
        s
    }

    pub(crate) fn new() -> Self {
        Self {
            by_hash: HashMap::new(),
        }
    }

    /// Push `h` onto `haps` if no existing haplotype has the same bases + ref flag.
    pub(crate) fn insert(&mut self, haps: &mut Vec<Haplotype>, h: Haplotype) -> bool {
        let key = (hash_hap_bases(&h.bases), h.is_reference);
        if let Some(idxs) = self.by_hash.get(&key) {
            if idxs.iter().any(|&i| haps[i].bases == h.bases) {
                return false;
            }
        }
        self.by_hash.entry(key).or_default().push(haps.len());
        haps.push(h);
        true
    }

    pub(crate) fn contains(&self, haps: &[Haplotype], bases: &[u8], is_reference: bool) -> bool {
        let key = (hash_hap_bases(bases), is_reference);
        self.by_hash
            .get(&key)
            .is_some_and(|idxs| idxs.iter().any(|&i| haps[i].bases == bases))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cigar::{Cigar, CigarOperator};

    #[test]
    fn collapse_dangling_tail_dupes_quad4_shape() {
        let ref_bases = b"TGCATGACTGATCGTACGATTCGAGCTAGTCGATCGATGCTAGCTAGGCTAACGTTAGCTAGTAACTGATCGATCGATACGTACGT";
        let ref_hap = Haplotype::new(ref_bases.to_vec(), true);
        let good_alt = b"TGCATGACTGATCGTACGATTCGAGCTAGTCGATCGATGCTAGCTAGGCTAACGTTAGCGAGTAACTGATCGATCGATACGTACGT";
        let typo_alt = b"TGCATGACTGATCGTACGATTCGAGCTAGTCGATCGATGCTAGCTAGGCTAACGTTAGCGAGTAACTGATCGATCGATATGTACGT";
        let mut haps = vec![
            Haplotype::new(ref_bases.to_vec(), true),
            Haplotype::new(typo_alt.to_vec(), false),
            {
                let mut h = Haplotype::new(good_alt.to_vec(), false);
                h.score = 10.0;
                h
            },
            {
                let mut h = Haplotype::new(typo_alt.to_vec(), false);
                h.score = 20.0;
                h
            },
        ];
        collapse_dangling_tail_alt_duplicates(&mut haps, &ref_hap);
        assert_eq!(haps.len(), 2);
        assert!(haps.iter().any(|h| h.bases.as_slice() == good_alt));
        assert!(!haps
            .iter()
            .any(|h| h.bases.as_slice() == typo_alt && !h.is_reference));
    }

    #[test]
    fn prune_drops_short_non_reference_fragments() {
        let ref_bases = vec![b'A'; 68];
        // CLONE: needed because haplotype constructor takes owned bases.
        let ref_hap = Haplotype::new(ref_bases.clone(), true);
        let mut haps = vec![
            Haplotype::new(b"GATACGATTCG".to_vec(), false),
            // CLONE: needed because haplotype constructor takes owned bases.
            Haplotype::new(ref_bases.clone(), true),
            Haplotype::new(vec![b'C'; 68], false),
        ];
        prune_fragment_non_reference_haplotypes(&mut haps, &ref_hap, 30);
        assert_eq!(haps.len(), 2);
        assert!(haps.iter().any(|h| h.is_reference));
    }

    #[test]
    fn ref_equivalent_non_ref_label_becomes_single_reference() {
        let ref_bases = b"ACGTACGT".to_vec();
        let mut haps = vec![
            // CLONE: needed because haplotype constructor takes owned bases.
            Haplotype::new(ref_bases.clone(), false),
            // CLONE: needed because haplotype constructor takes owned bases.
            Haplotype::new(ref_bases.clone(), true),
        ];
        normalize_ref_equivalent_haplotypes(&mut haps, &ref_bases);
        assert_eq!(haps.len(), 1);
        assert!(haps[0].is_reference);
    }

    #[test]
    fn trim_subinterval_keeps_bases_and_cigar() {
        let mut h = Haplotype::new(b"ACGTACGT".to_vec(), false);
        let mut cigar = Cigar::new();
        cigar.push(8, CigarOperator::Match);
        h.cigar = Some(cigar);
        h.genome_loc = Some(GenomeLoc::new(1, 8));
        let sub = GenomeLoc::new(3, 5);
        let t = h.trim(&sub, false).expect("trim");
        assert_eq!(t.bases, b"GTA");
        assert_eq!(t.genome_loc, Some(sub));
    }
}
