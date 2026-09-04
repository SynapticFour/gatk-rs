//! GATK `GATKVariantContextUtils.reverseTrimAlleles` (4.4.0.0 SHA `2dbc0258`).
//!
//! Called **after** unused-ALT subset, and only when allele count changed versus the
//! merged VC (`HaplotypeCallerGenotypingEngine.makeAnnotatedCall`). Does not recalculate
//! GLs, AD, or GT indices. Start is unchanged (`trimForward=false`).

const SPAN_DEL: &str = "*";

/// Reverse-trim shared trailing bases on REF + ALTs.
///
/// Java `reverseTrimAlleles` → `trimAlleles(vc, false, true)`:
/// - no-op if fewer than two alleles, or any non-span-del allele has length 1;
/// - REF is included in the suffix calculation;
/// - suffix bases common to all non-symbolic / non-`*` alleles are clipped;
/// - if clipping would empty an allele, one base is restored;
/// - start coordinate is unchanged (forward trim is off);
/// - allele order is preserved.
///
/// `AlignmentUtils.normalizeAlleles(..., maxShift=0, trim=true)` also walks a common
/// prefix and a left-align loop. With `maxShift=0` those walks cannot increase
/// `endShift` beyond the suffix loop (suffix exhausts common last bases; left-align
/// then sees either differing last bases or `startShift==0`). Prefix clip is discarded
/// because `trimForward=false`. End-clip + restore-one-base is therefore sufficient.
pub fn reverse_trim_alleles(ref_allele: &str, alt_alleles: &[String]) -> (String, Vec<String>) {
    let mut alleles: Vec<String> = Vec::with_capacity(1 + alt_alleles.len());
    alleles.push(ref_allele.to_string());
    alleles.extend(alt_alleles.iter().cloned());
    if alleles.len() <= 1 {
        return (ref_allele.to_string(), alt_alleles.to_vec());
    }
    if alleles
        .iter()
        .any(|a| a != SPAN_DEL && !is_symbolic(a) && a.len() == 1)
    {
        return (ref_allele.to_string(), alt_alleles.to_vec());
    }
    let clip = reverse_suffix_clip_count(&alleles);
    if clip == 0 {
        return (ref_allele.to_string(), alt_alleles.to_vec());
    }
    let trimmed: Vec<String> = alleles
        .iter()
        .map(|a| {
            if a == SPAN_DEL || is_symbolic(a) {
                a.clone()
            } else {
                a[..a.len() - clip].to_string()
            }
        })
        .collect();
    let new_ref = trimmed[0].clone();
    let new_alts = trimmed[1..].to_vec();
    (new_ref, new_alts)
}

fn is_symbolic(allele: &str) -> bool {
    allele.starts_with('<') && allele.ends_with('>')
}

/// Java `AlignmentUtils.normalizeAlleles` first loop (`trim=true`) + restore-one-base-at-end.
fn reverse_suffix_clip_count(alleles: &[String]) -> usize {
    let seqs: Vec<&[u8]> = alleles
        .iter()
        .filter(|a| *a != SPAN_DEL && !is_symbolic(a))
        .map(|a| a.as_bytes())
        .collect();
    if seqs.len() < 2 {
        return 0;
    }
    let mut remaining: Vec<usize> = seqs.iter().map(|s| s.len()).collect();
    let mut clip = 0usize;
    loop {
        if remaining.iter().any(|&r| r == 0) {
            break;
        }
        let last = seqs[0][remaining[0] - 1];
        if !seqs
            .iter()
            .zip(remaining.iter())
            .all(|(s, &r)| s[r - 1] == last)
        {
            break;
        }
        for r in &mut remaining {
            *r -= 1;
        }
        clip += 1;
    }
    if remaining.iter().any(|&r| r == 0) {
        clip.saturating_sub(1)
    } else {
        clip
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn alts(xs: &[&str]) -> Vec<String> {
        xs.iter().map(|s| (*s).to_string()).collect()
    }

    #[test]
    fn reverse_trim_common_suffix_after_unused_alt_subset() {
        let (r, a) = reverse_trim_alleles("TG", &alts(&["CG"]));
        assert_eq!(r, "T");
        assert_eq!(a, vec!["C".to_string()]);
    }

    #[test]
    fn reverse_trim_no_common_suffix_unchanged() {
        let (r, a) = reverse_trim_alleles("AT", &alts(&["GC"]));
        assert_eq!(r, "AT");
        assert_eq!(a, vec!["GC".to_string()]);
    }

    #[test]
    fn reverse_trim_multiple_base_common_suffix() {
        let (r, a) = reverse_trim_alleles("ACGT", &alts(&["CCGT"]));
        assert_eq!(r, "A");
        assert_eq!(a, vec!["C".to_string()]);
    }

    #[test]
    fn reverse_trim_restores_one_base_when_suffix_would_empty() {
        let (r, a) = reverse_trim_alleles("AAA", &alts(&["AAA"]));
        assert_eq!(r, "A");
        assert_eq!(a, vec!["A".to_string()]);
    }

    #[test]
    fn reverse_trim_skips_when_any_allele_already_length_one() {
        let (r, a) = reverse_trim_alleles("TG", &alts(&["T", "CG"]));
        assert_eq!(r, "TG");
        assert_eq!(a, vec!["T".to_string(), "CG".to_string()]);
    }

    #[test]
    fn reverse_trim_multiallelic_only_suffix_common_to_all() {
        let (r, a) = reverse_trim_alleles("TG", &alts(&["CG", "CGG"]));
        assert_eq!(r, "T");
        assert_eq!(a, vec!["C".to_string(), "CG".to_string()]);
    }

    #[test]
    fn reverse_trim_does_not_left_trim_common_prefix() {
        let (r, a) = reverse_trim_alleles("GAA", &alts(&["GTA"]));
        assert_eq!(r, "GA");
        assert_eq!(a, vec!["GT".to_string()]);
    }

    #[test]
    fn reverse_trim_single_shared_suffix_leaving_one_base() {
        let (r, a) = reverse_trim_alleles("AT", &alts(&["GT"]));
        assert_eq!(r, "A");
        assert_eq!(a, vec!["G".to_string()]);
    }
}
