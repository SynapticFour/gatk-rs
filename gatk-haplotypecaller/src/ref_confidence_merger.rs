//! GATK `ReferenceConfidenceVariantContextMerger` (H-D02) — diploid gVCF merge slice.

use gatk_common::{GatkError, GatkResult};
use std::collections::HashMap;

fn diploid_genotype_pairs(allele_count: usize) -> Vec<(usize, usize)> {
    let mut pairs = Vec::new();
    for j in 0..allele_count {
        for i in 0..=j {
            pairs.push((i, j));
        }
    }
    pairs
}

pub const NON_REF_ALLELE: &str = "<NON_REF>";
pub const SPAN_DEL_ALLELE: &str = "*";

/// One input gVCF record for merge.
/// # Invariants
/// `alleles` and per-sample `genotypes` must be consistent cardinality before merge remapping.
/// `start` is 1-based VCF position for this input record.
/// # Ownership
/// Owns source label, alleles, and genotype field vectors.
/// # Mutation
/// Immutable merge input; merger produces [`RefConfidenceMergeResult`].
/// # Biological assumptions
/// Diploid gVCF records with optional PL/AD per sample at one genomic start.
/// # Java equivalence
/// GATK `ReferenceConfidenceVariantContextMerger` input VC (H-D02).
#[derive(Debug, Clone)]
pub struct MergeVcInput {
    pub source: String,
    pub start: u64,
    pub alleles: Vec<MergeAllele>,
    pub genotypes: Vec<MergeGenotype>,
}

/// One allele entry in a gVCF merge input (`ReferenceConfidenceVariantContextMerger`).
/// # Invariants
/// Exactly one allele per input VC should have `is_reference == true` when representing a SNP site.
/// `<NON_REF>` and `*` use canonical base strings from merger constants.
/// # Ownership
/// Owns display/base string and reference flag.
/// # Mutation
/// Immutable merge input cell.
/// # Biological assumptions
/// Allele string is VCF display bases including symbolic gVCF alleles.
/// # Java equivalence
/// GATK `ReferenceConfidenceVariantContextMerger` allele representation (H-D02).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MergeAllele {
    pub bases: String,
    pub is_reference: bool,
}

/// Per-sample genotype fields on one input gVCF for reference-confidence merge.
/// # Invariants
/// PL and AD lengths must match merged allele cardinality after remapping (enforced by merger).
/// # Ownership
/// Owns sample name and optional PL/AD vectors.
/// # Mutation
/// Immutable input; merger produces [`MergeGenotypeOut`].
/// # Biological assumptions
/// Diploid gVCF sample with PL/AD optional per input record.
/// # Java equivalence
/// GATK `ReferenceConfidenceVariantContextMerger` genotype inputs.
#[derive(Debug, Clone)]
pub struct MergeGenotype {
    pub sample: String,
    pub pl: Option<Vec<i32>>,
    pub ad: Option<Vec<i32>>,
}

/// Merged site for parity dumps.
/// # Invariants
/// `alleles` lists merged site alleles in emission order; `has_non_ref` flags `<NON_REF>` presence.
/// Output genotypes align PL/AD lengths with merged allele cardinality when present.
/// # Ownership
/// Owns contig, position, allele strings, and merged per-sample genotypes.
/// # Mutation
/// Immutable merge output for dumps and tests.
/// # Biological assumptions
/// Merged gVCF site represents union of overlapping reference-confidence inputs at one start.
/// # Java equivalence
/// GATK `ReferenceConfidenceVariantContextMerger` merged VC output (H-D02).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RefConfidenceMergeResult {
    pub contig: String,
    pub pos: u64,
    pub alleles: Vec<String>,
    pub genotypes: Vec<MergeGenotypeOut>,
    pub has_non_ref: bool,
}

/// Merged per-sample genotype after reference-confidence VC merge.
/// # Invariants
/// PL/AD vectors align with merged site allele list order when present.
/// # Ownership
/// Owns sample name and optional PL/AD outputs.
/// # Mutation
/// Immutable merge output cell.
/// # Biological assumptions
/// Same sample name preserved across merged inputs at one genomic start.
/// # Java equivalence
/// GATK merged gVCF genotype emission from `ReferenceConfidenceVariantContextMerger`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MergeGenotypeOut {
    pub name: String,
    pub pl: Option<Vec<i32>>,
    pub ad: Option<Vec<i32>>,
}

#[derive(Debug, Clone)]
struct RemappedVc {
    vc: MergeVcInput,
    remapped: Vec<MergeAllele>,
    is_spanning: bool,
}

fn allele_key(a: &MergeAllele) -> String {
    if a.bases == NON_REF_ALLELE {
        return NON_REF_ALLELE.to_string();
    }
    if a.bases == SPAN_DEL_ALLELE {
        return SPAN_DEL_ALLELE.to_string();
    }
    a.bases.clone()
}

fn is_non_ref(a: &MergeAllele) -> bool {
    a.bases == NON_REF_ALLELE
}

fn is_span_del(a: &MergeAllele) -> bool {
    a.bases == SPAN_DEL_ALLELE
}

fn replace_with_no_calls_and_dels(vc: &MergeVcInput) -> Vec<MergeAllele> {
    let mut out = Vec::with_capacity(vc.alleles.len());
    out.push(MergeAllele {
        bases: ".".to_string(),
        is_reference: true,
    });
    for a in vc.alleles.iter().skip(1) {
        if is_non_ref(a) {
            // CLONE: needed because owned element into collection.
            out.push(a.clone());
        } else if a.bases.len() < vc.alleles[0].bases.len() {
            out.push(MergeAllele {
                bases: SPAN_DEL_ALLELE.to_string(),
                is_reference: false,
            });
        } else {
            out.push(MergeAllele {
                bases: ".".to_string(),
                is_reference: false,
            });
        }
    }
    out
}

fn extend_allele(allele: &MergeAllele, extra: usize, ref_bases: &[u8]) -> MergeAllele {
    if extra == 0 || is_non_ref(allele) || is_span_del(allele) || allele.bases == "." {
        return allele.clone();
    }
    let mut bytes = allele.bases.as_bytes().to_vec();
    bytes.extend_from_slice(&ref_bases[ref_bases.len() - extra..]);
    MergeAllele {
        bases: String::from_utf8_lossy(&bytes).into_owned(),
        is_reference: false,
    }
}

fn remap_alleles(vc: &MergeVcInput, ref_allele: &MergeAllele) -> GatkResult<Vec<MergeAllele>> {
    let vc_ref = &vc.alleles[0];
    let ref_bases = ref_allele.bases.as_bytes();
    let extra = ref_bases.len().saturating_sub(vc_ref.bases.len());
    if ref_bases.len() < vc_ref.bases.len() {
        return Err(GatkError::argument(
            "remapAlleles: wrong reference selected",
        ));
    }
    let mut out = vec![ref_allele.clone()];
    for a in vc.alleles.iter().skip(1) {
        if is_non_ref(a) || is_span_del(a) || a.bases == "." {
            // CLONE: needed because owned element into collection.
            out.push(a.clone());
        } else {
            out.push(extend_allele(a, extra, ref_bases));
        }
    }
    Ok(out)
}

fn determine_reference_allele(
    vcs: &[MergeVcInput],
    loc_start: u64,
    ref_base: Option<u8>,
) -> GatkResult<Option<MergeAllele>> {
    let mut ref_allele: Option<MergeAllele> = None;
    for vc in vcs {
        if vc.start != loc_start {
            continue;
        }
        // CLONE: needed because MergeAllele is owned by the chosen-reference result.
        let my_ref = vc.alleles[0].clone();
        ref_allele = Some(match ref_allele {
            None => my_ref,
            Some(prev) => {
                if prev.bases.len() < my_ref.bases.len() {
                    my_ref
                } else if my_ref.bases.len() < prev.bases.len() {
                    prev
                } else if prev.bases != my_ref.bases {
                    return Err(GatkError::argument(format!(
                        "inconsistent references at {loc_start}: {} vs {}",
                        prev.bases, my_ref.bases
                    )));
                } else {
                    prev
                }
            }
        });
    }
    Ok(ref_allele.or_else(|| {
        ref_base.map(|b| MergeAllele {
            bases: String::from(b as char),
            is_reference: true,
        })
    }))
}

fn collect_target_alleles(
    pairs: &[RemappedVc],
    ref_allele: &MergeAllele,
    remove_non_ref: bool,
) -> Vec<MergeAllele> {
    let mut ordered: Vec<String> = vec![allele_key(ref_allele)];
    let mut saw_span_del = false;
    let mut saw_non_spanning = false;

    for p in pairs {
        for a in &p.remapped {
            let key = allele_key(a);
            if a.is_reference || is_non_ref(a) || a.bases == "." {
                continue;
            }
            if is_span_del(a) {
                saw_span_del = true;
                continue;
            }
            saw_non_spanning = true;
            if !ordered.contains(&key) {
                ordered.push(key);
            }
        }
        if p.is_spanning
            && p.vc
                .alleles
                .iter()
                .any(|a| !is_non_ref(a) && !a.is_reference)
        {
            saw_span_del = true;
        }
    }

    let mut alleles: Vec<MergeAllele> = ordered
        .iter()
        .map(|k| {
            if k == &allele_key(ref_allele) {
                ref_allele.clone()
            } else if k == SPAN_DEL_ALLELE {
                MergeAllele {
                    bases: SPAN_DEL_ALLELE.to_string(),
                    is_reference: false,
                }
            } else {
                MergeAllele {
                    bases: k.clone(),
                    is_reference: false,
                }
            }
        })
        .collect();

    if saw_span_del && (saw_non_spanning || !remove_non_ref) {
        alleles.push(MergeAllele {
            bases: SPAN_DEL_ALLELE.to_string(),
            is_reference: false,
        });
    }
    if !remove_non_ref {
        alleles.push(MergeAllele {
            bases: NON_REF_ALLELE.to_string(),
            is_reference: false,
        });
    }
    alleles
}

/// GATK `AlleleSubsettingUtils.getIndexesOfRelevantAllelesForGVCF` (diploid, non-somatic).
fn indexes_of_relevant_alleles_for_gvcf(
    remapped: &[MergeAllele],
    target: &[MergeAllele],
) -> GatkResult<Vec<usize>> {
    let non_ref_idx = remapped
        .iter()
        .position(is_non_ref)
        .ok_or_else(|| GatkError::argument("remapped alleles must contain <NON_REF>"))?;
    let mut map = vec![0usize; target.len()];
    for (i, t) in target.iter().enumerate().skip(1) {
        let key = allele_key(t);
        let idx = remapped
            .iter()
            .position(|a| allele_key(a) == key)
            .unwrap_or(non_ref_idx);
        map[i] = idx;
    }
    Ok(map)
}

fn genotype_index_map_diploid(per_sample_map: &[usize], old_allele_count: usize) -> Vec<usize> {
    let new_ac = per_sample_map.len();
    let old_pairs = diploid_genotype_pairs(old_allele_count);
    let new_pairs = diploid_genotype_pairs(new_ac);
    let old_index: HashMap<(usize, usize), usize> = old_pairs
        .iter()
        .enumerate()
        .map(|(idx, &(i, j))| ((i, j), idx))
        .collect();
    new_pairs
        .iter()
        .map(|&(ni, nj)| {
            let oi = per_sample_map[ni];
            let oj = per_sample_map[nj];
            let (oi, oj) = if oi <= oj { (oi, oj) } else { (oj, oi) };
            *old_index.get(&(oi, oj)).expect("genotype index")
        })
        .collect()
}

fn remap_pl(old_pl: &[i32], old_allele_count: usize, per_sample_map: &[usize]) -> Vec<i32> {
    let map = genotype_index_map_diploid(per_sample_map, old_allele_count);
    map.iter().map(|&i| old_pl[i]).collect()
}

/// GATK `AlleleSubsettingUtils.generateAD` / `remapRLengthList`.
fn remap_ad(old_ad: &[i32], per_sample_map: &[usize]) -> Vec<i32> {
    per_sample_map
        .iter()
        .map(|&old_index| {
            if old_index >= old_ad.len() {
                0
            } else {
                old_ad[old_index]
            }
        })
        .collect()
}

fn merge_genotypes(
    vc: &RemappedVc,
    target: &[MergeAllele],
    samples_are_uniquified: bool,
) -> GatkResult<Vec<MergeGenotypeOut>> {
    let mut out = Vec::new();
    let per_sample = indexes_of_relevant_alleles_for_gvcf(&vc.remapped, target)?;
    for g in &vc.vc.genotypes {
        let name = if samples_are_uniquified {
            format!("{}.{}", g.sample, vc.vc.source)
        } else {
            // CLONE: needed because owned sample id for carry/map.
            g.sample.clone()
        };
        let pl =
            g.pl.as_ref()
                .map(|pl| remap_pl(pl, vc.remapped.len(), &per_sample));
        let ad = g.ad.as_ref().map(|ad| remap_ad(ad, &per_sample));
        out.push(MergeGenotypeOut { name, pl, ad });
    }
    Ok(out)
}

/// Production merge for gVCF records at one locus (H-D02).
pub fn merge_reference_confidence(
    contig: &str,
    loc_start: u64,
    vcs: &[MergeVcInput],
    ref_base: Option<u8>,
    remove_non_ref: bool,
    samples_are_uniquified: bool,
) -> GatkResult<Option<RefConfidenceMergeResult>> {
    if vcs.is_empty() {
        return Err(GatkError::argument("merge: empty inputs"));
    }
    let ref_allele = match determine_reference_allele(vcs, loc_start, ref_base)? {
        Some(a) => a,
        None => return Ok(None),
    };

    let mut pairs = Vec::with_capacity(vcs.len());
    for vc in vcs {
        let is_spanning = loc_start != vc.start;
        let remapped = if is_spanning {
            replace_with_no_calls_and_dels(vc)
        } else {
            remap_alleles(vc, &ref_allele)?
        };
        pairs.push(RemappedVc {
            // CLONE: needed because owned variant/locus record in output map.
            vc: vc.clone(),
            remapped,
            is_spanning,
        });
    }

    let target = collect_target_alleles(&pairs, &ref_allele, remove_non_ref);
    let mut genotypes = Vec::new();
    for p in &pairs {
        genotypes.extend(merge_genotypes(p, &target, samples_are_uniquified)?);
    }

    let alleles: Vec<String> = target.iter().map(|a| a.bases.clone()).collect();
    Ok(Some(RefConfidenceMergeResult {
        contig: contig.to_string(),
        pos: loc_start,
        alleles,
        genotypes,
        has_non_ref: target.iter().any(is_non_ref),
    }))
}

// --------------------------------------------------------------------------
// Unit tests — hand-calculated diploid PL remapping (different ALT sets)
// --------------------------------------------------------------------------
// Diploid VCF PL order for allele indices 0..N-1 (htsjdk / GATK):
// for j in 0..N { for i in 0..=j { genotype i/j } }
// So for alleles [REF, A1, A2, <NON_REF>] (N=4) the 10 PLs are:
// 0/0, 0/1, 1/1, 0/2, 1/2, 2/2, 0/3, 1/3, 2/3, 3/3
// Remap contract (GATK `generatePL` + `getIndexesOfRelevantAllelesForGVCF`):
// per_sample_map[new_allele] = old index of same allele, else old <NON_REF> index
// new_PL[k] = old_PL[ old_genotype(map[a], map[b]) ] for new genotype a/b at slot k
// No re-normalization after remap (min need not be recomputed — already 0 in source).

#[cfg(test)]
#[allow(clippy::unwrap_used, clippy::expect_used)]
mod pl_remap_tests {
    use super::*;

    fn allele(bases: &str, is_ref: bool) -> MergeAllele {
        MergeAllele {
            bases: bases.to_string(),
            is_reference: is_ref,
        }
    }

    fn snp_vc(sample: &str, alt: &str, pl: &[i32], ad: &[i32]) -> MergeVcInput {
        MergeVcInput {
            source: sample.to_string(),
            start: 10,
            alleles: vec![
                allele("A", true),
                allele(alt, false),
                allele(NON_REF_ALLELE, false),
            ],
            genotypes: vec![MergeGenotype {
                sample: sample.to_string(),
                pl: Some(pl.to_vec()),
                ad: Some(ad.to_vec()),
            }],
        }
    }

    fn pl_of(merged: &RefConfidenceMergeResult, sample: &str) -> Vec<i32> {
        merged
            .genotypes
            .iter()
            .find(|g| g.name == sample)
            .and_then(|g| g.pl.clone())
            .expect("sample PL")
    }

    fn alts(merged: &RefConfidenceMergeResult) -> Vec<&str> {
        merged.alleles.iter().skip(1).map(|s| s.as_str()).collect()
    }

    /// Historical CombineGVCFs mini failure at chr1:10.
    /// Inputs (Java `endPreviousStates` last→first ⇒ merge order SAMPLE2 then SAMPLE1):
    /// SAMPLE1 alleles `[A,G,<NON_REF>]`, PL `[100,0,100,100,100,100]`
    /// (zero at old 0/1 = REF/G)
    /// SAMPLE2 alleles `[A,T,<NON_REF>]`, PL `[90,0,90,90,90,90]`
    /// (zero at old 0/1 = REF/T)
    /// Target alleles: `[A,T,G,<NON_REF>]` (T discovered before G).
    /// SAMPLE1 map: T→old NON_REF(2), G→1, NON_REF→2
    /// new 0/2 (=REF/G) ← old 0/1 → PL[1]=0
    /// ⇒ `[100,100,100, 0, 100,100,100,100,100,100]`
    /// SAMPLE2 map: T→1, G→old NON_REF(2), NON_REF→2
    /// new 0/1 (=REF/T) ← old 0/1 → PL[1]=0
    /// ⇒ `[90, 0, 90,90,90,90,90,90,90,90]`
    #[test]
    fn pl01_two_samples_different_snp_alts_java_order() {
        let s2 = snp_vc("SAMPLE2", "T", &[90, 0, 90, 90, 90, 90], &[8, 8, 0]);
        let s1 = snp_vc("SAMPLE1", "G", &[100, 0, 100, 100, 100, 100], &[10, 10, 0]);
        let merged = merge_reference_confidence("chr1", 10, &[s2, s1], None, false, false)
            .unwrap()
            .expect("merged");
        assert_eq!(alts(&merged), vec!["T", "G", NON_REF_ALLELE]);
        assert_eq!(
            pl_of(&merged, "SAMPLE1"),
            vec![100, 100, 100, 0, 100, 100, 100, 100, 100, 100]
        );
        assert_eq!(
            pl_of(&merged, "SAMPLE2"),
            vec![90, 0, 90, 90, 90, 90, 90, 90, 90, 90]
        );
    }

    /// Same inputs as pl01 but SAMPLE1 first ⇒ ALT order G,T,<NON_REF>.
    /// SAMPLE1 zero stays at new 0/1; SAMPLE2 zero moves to new 0/2.
    /// (Documents why CombineGVCFs must reverse covering VCs to match Java.)
    #[test]
    fn pl02_opposite_merge_order_swaps_alt_slots() {
        let s1 = snp_vc("SAMPLE1", "G", &[100, 0, 100, 100, 100, 100], &[10, 10, 0]);
        let s2 = snp_vc("SAMPLE2", "T", &[90, 0, 90, 90, 90, 90], &[8, 8, 0]);
        let merged = merge_reference_confidence("chr1", 10, &[s1, s2], None, false, false)
            .unwrap()
            .expect("merged");
        assert_eq!(alts(&merged), vec!["G", "T", NON_REF_ALLELE]);
        // SAMPLE1: G is allele 1 → zero at index 1 (0/1)
        assert_eq!(
            pl_of(&merged, "SAMPLE1"),
            vec![100, 0, 100, 100, 100, 100, 100, 100, 100, 100]
        );
        // SAMPLE2: T is allele 2 → zero at index 3 (0/2)
        assert_eq!(
            pl_of(&merged, "SAMPLE2"),
            vec![90, 90, 90, 0, 90, 90, 90, 90, 90, 90]
        );
    }

    /// Hom-ref-only gVCF + SNP gVCF: missing ALT maps through `<NON_REF>` PL cells.
    /// Hom-ref sample alleles `[A,<NON_REF>]`, PL `[0,30,300]` (0/0,0/1,1/1).
    /// SNP sample alleles `[A,G,<NON_REF>]`, PL `[50,0,50,50,50,50]`.
    /// Target `[A,G,<NON_REF>]`.
    /// Hom-ref map: G→old NON_REF(1), NON_REF→1
    /// 0/0←0 → 0
    /// 0/1←0/1 → 30
    /// 1/1←1/1 → 300
    /// 0/2←0/1 → 30
    /// 1/2←1/1 → 300
    /// 2/2←1/1 → 300
    #[test]
    fn pl03_homref_sample_missing_alt_uses_nonref_pl() {
        let hom = MergeVcInput {
            source: "HOM".into(),
            start: 10,
            alleles: vec![allele("A", true), allele(NON_REF_ALLELE, false)],
            genotypes: vec![MergeGenotype {
                sample: "HOM".into(),
                pl: Some(vec![0, 30, 300]),
                ad: None,
            }],
        };
        let snp = snp_vc("SNP", "G", &[50, 0, 50, 50, 50, 50], &[5, 5, 0]);
        let merged = merge_reference_confidence("chr1", 10, &[hom, snp], None, false, false)
            .unwrap()
            .expect("merged");
        assert_eq!(alts(&merged), vec!["G", NON_REF_ALLELE]);
        assert_eq!(pl_of(&merged, "HOM"), vec![0, 30, 300, 30, 300, 300]);
        assert_eq!(pl_of(&merged, "SNP"), vec![50, 0, 50, 50, 50, 50]);
        // Min PL remains 0 — no spurious re-normalization that would shift other cells.
        assert_eq!(pl_of(&merged, "HOM").iter().copied().min(), Some(0));
        assert_eq!(pl_of(&merged, "SNP").iter().copied().min(), Some(0));
    }

    /// Three samples, three distinct ALTs. Merge order C,B,A ⇒ ALT T,G,C,<NON_REF>.
    /// Each sample's PL zero is at old 0/1 (het REF/own-ALT). After remap:
    /// A (ALT C = new allele 3): zero at new 0/3 → PL index 6
    /// B (ALT G = new allele 2): zero at new 0/2 → PL index 3
    /// C (ALT T = new allele 1): zero at new 0/1 → PL index 1
    #[test]
    fn pl04_three_samples_three_alts_hand_indices() {
        let a = snp_vc("A", "C", &[40, 0, 40, 40, 40, 40], &[4, 4, 0]);
        let b = snp_vc("B", "G", &[50, 0, 50, 50, 50, 50], &[5, 5, 0]);
        let c = snp_vc("C", "T", &[60, 0, 60, 60, 60, 60], &[6, 6, 0]);
        // Java-style last→first discovery order among these three.
        let merged = merge_reference_confidence("chr1", 10, &[c, b, a], None, false, false)
            .unwrap()
            .expect("merged");
        assert_eq!(alts(&merged), vec!["T", "G", "C", NON_REF_ALLELE]);
        let n = 1 + alts(&merged).len(); // 5 alleles → 15 PLs
        assert_eq!(n * (n + 1) / 2, 15);

        let pl_a = pl_of(&merged, "A");
        let pl_b = pl_of(&merged, "B");
        let pl_c = pl_of(&merged, "C");
        assert_eq!(pl_a.len(), 15);
        assert_eq!(pl_a.iter().position(|&x| x == 0), Some(6)); // 0/3
        assert_eq!(pl_b.iter().position(|&x| x == 0), Some(3)); // 0/2
        assert_eq!(pl_c.iter().position(|&x| x == 0), Some(1)); // 0/1
                                                                // Non-zero cells stay at the constant filler from the source gVCF.
        assert!(pl_a.iter().all(|&x| x == 0 || x == 40));
        assert!(pl_b.iter().all(|&x| x == 0 || x == 50));
        assert!(pl_c.iter().all(|&x| x == 0 || x == 60));
    }

    /// AD remaps in lockstep with allele map (R-length), same scenario as pl01.
    /// SAMPLE1 AD `[10,10,0]` on `[A,G,<NON_REF>]` → `[10,0,10,0]` on `[A,T,G,<NON_REF>]`
    /// (T and NON_REF both take old NON_REF depth 0; G keeps 10).
    #[test]
    fn pl05_ad_remaps_with_allele_map() {
        let s2 = snp_vc("SAMPLE2", "T", &[90, 0, 90, 90, 90, 90], &[8, 8, 0]);
        let s1 = snp_vc("SAMPLE1", "G", &[100, 0, 100, 100, 100, 100], &[10, 10, 0]);
        let merged = merge_reference_confidence("chr1", 10, &[s2, s1], None, false, false)
            .unwrap()
            .expect("merged");
        let ad = |sample: &str| {
            merged
                .genotypes
                .iter()
                .find(|g| g.name == sample)
                .and_then(|g| g.ad.clone())
                .expect("AD")
        };
        assert_eq!(ad("SAMPLE1"), vec![10, 0, 10, 0]);
        assert_eq!(ad("SAMPLE2"), vec![8, 8, 0, 0]);
    }

    /// Sample without PL stays `None` after merge (no-call / absent likelihoods).
    #[test]
    fn pl06_missing_pl_stays_absent() {
        let with_pl = snp_vc("HAS", "G", &[10, 0, 10, 10, 10, 10], &[1, 1, 0]);
        let no_pl = MergeVcInput {
            source: "NOPL".into(),
            start: 10,
            alleles: vec![allele("A", true), allele(NON_REF_ALLELE, false)],
            genotypes: vec![MergeGenotype {
                sample: "NOPL".into(),
                pl: None,
                ad: None,
            }],
        };
        let merged = merge_reference_confidence("chr1", 10, &[with_pl, no_pl], None, false, false)
            .unwrap()
            .expect("merged");
        let no = merged.genotypes.iter().find(|g| g.name == "NOPL").unwrap();
        assert!(no.pl.is_none());
        assert!(merged
            .genotypes
            .iter()
            .find(|g| g.name == "HAS")
            .unwrap()
            .pl
            .is_some());
    }
}
