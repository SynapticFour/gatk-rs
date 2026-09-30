//! GATK `AssemblyRegionTrimmer` (GAP-B-06).
//! Shrinks an assembly region to the span of assembled variation plus genotyping padding
//! (post-assembly in Java `HaplotypeCallerEngine.callRegion`).

use crate::assembly_region_iterator::AssemblyRegion;
use crate::genome_loc::GenomePosition;
use crate::reference_context::ReferenceContext;
use gatk_core::reference::SequenceDictionary;
use std::cmp::{max, min};
use std::io::{BufRead, BufReader};
use std::path::Path;

/// Minimal variant record for trimming (from assembly events or fixture TSV).
/// # Invariants
/// `start` / `end` are 1-based inclusive genomic coordinates on `contig`.
/// `is_indel` selects genotyping padding width in the trimmer.
/// # Ownership
/// Owns contig string; lightweight event proxy for trim.
/// # Mutation
/// Immutable input to [`AssemblyRegionTrimmer::trim`].
/// # Biological assumptions
/// Represents assembled variation span used to shrink the genotyping window.
/// # Java equivalence
/// GATK trim input derived from assembly events / force-calling alleles.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TrimVariant {
    pub contig: String,
    pub start: u64,
    pub end: u64,
    pub is_indel: bool,
    /// Full reference allele, including the shared anchor base.
    pub ref_allele: String,
    /// Full alternate allele, including the shared anchor base.
    pub alt_allele: String,
}

impl TrimVariant {
    pub fn overlaps_active_region(&self, region: &AssemblyRegion) -> bool {
        self.contig == region.contig
            && self.start <= region.end.get()
            && self.end >= region.start.get()
    }
}

/// GATK `AssemblyRegionArgumentCollection` trimming parameters.
/// # Invariants
/// SNP vs indel/STR paddings are non-negative; legacy max-extension `-1` disables the cap.
/// HC defaults: SNP 20, indel/STR 75, legacy mode off.
/// # Ownership
/// Cloneable config nested in [`crate::engine::CallRegionArgs`].
/// # Mutation
/// Snapshot for trimmer construction.
/// # Biological assumptions
/// Extra bases around variation needed for accurate genotyping of SNPs vs indels.
/// # Java equivalence
/// GATK `AssemblyRegionArgumentCollection` genotyping-padding trim knobs.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AssemblyRegionTrimmerConfig {
    pub snp_padding_for_genotyping: u32,
    pub indel_padding_for_genotyping: u32,
    pub str_padding_for_genotyping: u32,
    /// Legacy only; `-1` disables the max-extension cap.
    pub max_extension_into_region_padding_legacy: i32,
    pub enable_legacy_assembly_region_trimming: bool,
}

impl Default for AssemblyRegionTrimmerConfig {
    fn default() -> Self {
        Self::gatk_defaults()
    }
}

impl AssemblyRegionTrimmerConfig {
    pub fn gatk_defaults() -> Self {
        Self {
            snp_padding_for_genotyping: 20,
            indel_padding_for_genotyping: 75,
            str_padding_for_genotyping: 75,
            max_extension_into_region_padding_legacy: 25,
            enable_legacy_assembly_region_trimming: false,
        }
    }
}

/// Result of [`AssemblyRegionTrimmer::trim`].
/// # Invariants
/// When `variation_present` is false, span options are typically `None`.
/// Padded spans include genotyping padding and stay within contig when clipping applies.
/// # Ownership
/// Owned scalar/option bundle.
/// # Mutation
/// Immutable trim outcome.
/// # Biological assumptions
/// Defines the post-assembly genotyping window around called variation.
/// # Java equivalence
/// GATK `AssemblyRegionTrimmer` trim result fields.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AssemblyRegionTrimResult {
    pub variation_present: bool,
    pub variant_start: Option<u64>,
    pub variant_end: Option<u64>,
    pub padded_variant_start: Option<u64>,
    pub padded_variant_end: Option<u64>,
}

/// Shrinks an assembly region to variation span plus genotyping padding.
/// # Invariants
/// Holds [`AssemblyRegionTrimmerConfig`] and contig length resolved at construction.
/// Trim results stay within `[1, contig_len]` when reference context is supplied.
/// # Ownership
/// Owns config and contig length; borrows [`AssemblyRegion`], variants, and optional [`ReferenceContext`] per [`Self::trim`].
/// # Mutation
/// Internal config is fixed after `new`; trimming is a pure function of inputs.
/// # Biological assumptions
/// Variants are left-aligned intervals on the assembly contig; indel vs SNP selects padding width.
/// # Java equivalence
/// GATK 4.4 `AssemblyRegionTrimmer` (post-assembly trim in `HaplotypeCallerEngine.callRegion`).
pub struct AssemblyRegionTrimmer {
    cfg: AssemblyRegionTrimmerConfig,
    contig_len: u64,
}

impl AssemblyRegionTrimmer {
    pub fn new(
        cfg: AssemblyRegionTrimmerConfig,
        dictionary: &SequenceDictionary,
        contig: &str,
    ) -> Self {
        let contig_len = dictionary
            .contig(contig)
            .map(|c| c.length)
            .unwrap_or(u64::MAX);
        Self { cfg, contig_len }
    }

    pub fn trim(
        &self,
        region: &AssemblyRegion,
        variants: &[TrimVariant],
        reference: Option<&ReferenceContext>,
    ) -> AssemblyRegionTrimResult {
        if self.cfg.enable_legacy_assembly_region_trimming {
            self.trim_legacy(region, variants)
        } else {
            self.trim_modern(region, variants, reference)
        }
    }

    /// Apply trim spans to a copy of `region` (active + padded bounds only; reads unchanged).
    pub fn apply_trim(
        region: &AssemblyRegion,
        result: &AssemblyRegionTrimResult,
    ) -> AssemblyRegion {
        if !result.variation_present {
            return region.clone();
        }
        let vs = result.variant_start.expect("variant_start");
        let ve = result.variant_end.expect("variant_end");
        let ps = result.padded_variant_start.expect("padded_variant_start");
        let pe = result.padded_variant_end.expect("padded_variant_end");
        trim_region_bounds(region, vs, ve, ps, pe)
    }

    fn trim_modern(
        &self,
        region: &AssemblyRegion,
        variants: &[TrimVariant],
        reference: Option<&ReferenceContext>,
    ) -> AssemblyRegionTrimResult {
        let in_region: Vec<&TrimVariant> = variants
            .iter()
            .filter(|v| v.overlaps_active_region(region))
            .collect();
        if in_region.is_empty() {
            return Self::no_variation();
        }

        let mut min_start = in_region.iter().map(|v| v.start).min().unwrap();
        let mut max_end = in_region.iter().map(|v| v.end).max().unwrap();
        let variant_start = max(region.start.get(), min_start);
        let variant_end = min(region.end.get(), max_end);

        for v in &in_region {
            let mut padding = if v.is_indel {
                self.cfg.indel_padding_for_genotyping
            } else {
                self.cfg.snp_padding_for_genotyping
            };
            if v.is_indel {
                if let Some(ref_ctx) = reference {
                    if let Some(repeat) = tandem_repeat_at_event(
                        ref_ctx.window_start,
                        ref_ctx.bases.as_slice(),
                        v.start,
                        v.ref_allele.as_bytes(),
                        v.alt_allele.as_bytes(),
                    ) {
                        padding = self
                            .cfg
                            .str_padding_for_genotyping
                            .saturating_add(repeat.padding_bases() as u32);
                    }
                }
            }
            min_start = min_start.min(v.start.saturating_sub(padding as u64).max(1));
            // GATK 4.4 `AssemblyRegionTrimmer.trim` (SHA 2dbc0258):
            // `maxEnd = Math.max(maxEnd, vc.getEnd() + padding)` — pad the farthest
            // event end once, do not accumulate padding per overlapping variant.
            max_end = max_end.max(v.end.saturating_add(padding as u64));
        }

        let padded_start = max(region.extended_start.get(), min_start);
        let padded_end = min(region.extended_end.get(), max_end);

        AssemblyRegionTrimResult {
            variation_present: true,
            variant_start: Some(variant_start),
            variant_end: Some(variant_end),
            padded_variant_start: Some(padded_start),
            padded_variant_end: Some(padded_end),
        }
    }

    fn trim_legacy(
        &self,
        region: &AssemblyRegion,
        variants: &[TrimVariant],
    ) -> AssemblyRegionTrimResult {
        let mut variant_start: Option<u64> = None;
        let mut variant_end: Option<u64> = None;
        let mut found_non_snp = false;

        for v in variants {
            if !v.overlaps_active_region(region) {
                continue;
            }
            found_non_snp |= v.is_indel;
            variant_start = Some(match variant_start {
                None => v.start,
                Some(s) => min(s, v.start),
            });
            variant_end = Some(match variant_end {
                None => v.end,
                Some(e) => max(e, v.end),
            });
        }

        let (Some(vs), Some(ve)) = (variant_start, variant_end) else {
            return Self::no_variation();
        };

        let padding = if found_non_snp {
            self.cfg.indel_padding_for_genotyping
        } else {
            self.cfg.snp_padding_for_genotyping
        };

        let maximum_span = expand_within_contig(
            region.start.get(),
            region.end.get(),
            self.cfg.max_extension_into_region_padding_legacy.max(0) as u64,
            self.contig_len,
        );
        let ideal_span = expand_within_contig(vs, ve, padding as u64, self.contig_len);
        let final_span =
            intersect_closed(maximum_span.0, maximum_span.1, ideal_span.0, ideal_span.1);
        let (fs, fe) = merge_with_contiguous(final_span.0, final_span.1, vs, ve);

        AssemblyRegionTrimResult {
            variation_present: true,
            variant_start: Some(vs),
            variant_end: Some(ve),
            padded_variant_start: Some(fs),
            padded_variant_end: Some(fe),
        }
    }

    fn no_variation() -> AssemblyRegionTrimResult {
        AssemblyRegionTrimResult {
            variation_present: false,
            variant_start: None,
            variant_end: None,
            padded_variant_start: None,
            padded_variant_end: None,
        }
    }
}

/// Identity trim used by iterator scaffold (post-assembly trim uses [`AssemblyRegionTrimmer`]).
pub fn trim_assembly_region(
    region: &AssemblyRegion,
    _cfg: &AssemblyRegionTrimmerConfig,
) -> AssemblyRegion {
    region.clone()
}

fn trim_region_bounds(
    region: &AssemblyRegion,
    variant_start: u64,
    variant_end: u64,
    padded_start: u64,
    padded_end: u64,
) -> AssemblyRegion {
    let new_start = GenomePosition::new_1based(max(
        region.start.get(),
        min(region.end.get(), variant_start),
    ));
    let new_end =
        GenomePosition::new_1based(max(new_start.get(), min(region.end.get(), variant_end)));
    let new_ext_start = GenomePosition::new_1based(max(
        region.extended_start.get(),
        min(region.extended_end.get(), padded_start),
    ));
    let new_ext_end = GenomePosition::new_1based(max(
        new_ext_start.get(),
        min(region.extended_end.get(), padded_end),
    ));
    // A3: `AssemblyRegion::clone` only bumps Arc/refcount on reads + SharedBases; bounds change
    // here does not deep-copy BAM payloads. Unique ownership for realign is via into_unique_records.
    let mut out = region.clone();
    out.start = new_start;
    out.end = new_end;
    out.extended_start = new_ext_start;
    out.extended_end = new_ext_end;
    out
}

fn expand_within_contig(start: u64, end: u64, padding: u64, contig_len: u64) -> (u64, u64) {
    let lo = start.saturating_sub(padding).max(1);
    let hi = end.saturating_add(padding).min(contig_len);
    (lo, hi)
}

fn intersect_closed(a_start: u64, a_end: u64, b_start: u64, b_end: u64) -> (u64, u64) {
    let s = max(a_start, b_start);
    let e = min(a_end, b_end);
    debug_assert!(s <= e, "intersect requires overlap");
    (s, e)
}

/// GATK `mergeWithContiguous` when intervals touch or overlap.
fn merge_with_contiguous(a_start: u64, a_end: u64, b_start: u64, b_end: u64) -> (u64, u64) {
    (min(a_start, b_start), max(a_end, b_end))
}

/// Repeat unit and per-allele unit counts for one indel.
///
/// GATK `TandemRepeat.getNumTandemRepeatUnits` (SHA `2dbc0258`): the reference
/// scan starts at `event_start + 1 - window_start`, so the anchor base is not
/// part of the repeat. The longer of the anchor-stripped ref and alt alleles
/// supplies the unit. Counts are the leading copies of that unit on the
/// reference allele and on the alternate allele, each continued into the
/// reference bases after the anchor. `None` when Java returns null.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TandemRepeatAtEvent {
    pub unit: Vec<u8>,
    /// `[reference count, alternate count]`, Java list order.
    pub counts: Vec<usize>,
}

impl TandemRepeatAtEvent {
    /// Bases added to `str_padding_for_genotyping`: max count times unit length.
    pub fn padding_bases(&self) -> usize {
        self.counts.iter().copied().max().unwrap_or(0) * self.unit.len()
    }
}

/// Java `GATKVariantContextUtils.getNumTandemRepeatUnits` for one alternate allele.
pub fn tandem_repeat_at_event(
    window_start: u64,
    ref_bases: &[u8],
    event_start: u64,
    ref_allele: &[u8],
    alt_allele: &[u8],
) -> Option<TandemRepeatAtEvent> {
    if ref_allele.is_empty() || alt_allele.is_empty() {
        return None;
    }
    if ref_allele.len() == alt_allele.len() {
        return None;
    }
    let start_index = event_start.checked_add(1)?.checked_sub(window_start)? as usize;
    if start_index > ref_bases.len() {
        return None;
    }
    let remaining = &ref_bases[start_index..];
    let ref_after = &ref_allele[1..];
    let alt_after = &alt_allele[1..];
    let long_bases = if alt_after.len() > ref_after.len() {
        alt_after
    } else {
        ref_after
    };
    if long_bases.is_empty() {
        return None;
    }
    let unit_len = find_repeated_substring(long_bases);
    let unit = long_bases[..unit_len].to_vec();
    let repetitions_in_ref = leading_repeat_count(&unit, ref_after);
    let ref_count = leading_repeat_count(&unit, &concat(ref_after, remaining));
    let alt_count = leading_repeat_count(&unit, &concat(alt_after, remaining));
    let ref_count = ref_count.checked_sub(repetitions_in_ref)?;
    let alt_count = alt_count.checked_sub(repetitions_in_ref)?;
    if ref_count == 0 || alt_count == 0 {
        return None;
    }
    Some(TandemRepeatAtEvent {
        unit,
        counts: vec![ref_count, alt_count],
    })
}

fn concat(left: &[u8], right: &[u8]) -> Vec<u8> {
    let mut out = Vec::with_capacity(left.len() + right.len());
    out.extend_from_slice(left);
    out.extend_from_slice(right);
    out
}

/// Java `findRepeatedSubstring`. A string that is not a pure tandem returns its own length.
fn find_repeated_substring(bases: &[u8]) -> usize {
    for rep_len in 1..=bases.len() {
        let unit = &bases[..rep_len];
        let mut start = rep_len;
        let mut all_match = true;
        while start < bases.len() {
            let end = start + rep_len;
            if end > bases.len() || &bases[start..end] != unit {
                all_match = false;
                break;
            }
            start += rep_len;
        }
        if all_match {
            return rep_len;
        }
    }
    bases.len()
}

/// Java `findNumberOfRepetitions(..., leadingRepeats = true)`.
fn leading_repeat_count(unit: &[u8], test: &[u8]) -> usize {
    if unit.is_empty() || test.len() < unit.len() {
        return 0;
    }
    let mut count = 0usize;
    let mut start = 0usize;
    while start + unit.len() <= test.len() {
        if &test[start..start + unit.len()] != unit {
            break;
        }
        count += 1;
        start += unit.len();
    }
    count
}

/// Inclusive-span homopolymer used before 6R.310. Not the trim-padding contract.
/// 6R.292 locks that this returns `None` for a one-base event span.
pub fn longest_str_len_at_variant(
    ref_ctx: &ReferenceContext,
    start: u64,
    end: u64,
) -> Option<usize> {
    let bases = ref_ctx.bases.as_slice();
    if bases.is_empty() {
        return None;
    }
    let win_start = ref_ctx.start;
    let rel_start = start.saturating_sub(win_start) as usize;
    let rel_end = end.saturating_sub(win_start) as usize;
    if rel_start >= bases.len() {
        return None;
    }
    let rel_end = rel_end.min(bases.len().saturating_sub(1));
    let slice = &bases[rel_start..=rel_end];
    if slice.len() < 2 {
        return None;
    }
    let mut best = 1usize;
    let mut run = 1usize;
    for i in 1..slice.len() {
        if slice[i] == slice[i - 1] {
            run += 1;
            best = best.max(run);
        } else {
            run = 1;
        }
    }
    if best >= 3 {
        Some(best)
    } else {
        None
    }
}

/// Load trim variants from gate fixture TSV (`contig`, `start`, `end`, `is_indel`).
pub fn load_trim_variants_tsv(path: &Path) -> Result<Vec<TrimVariant>, String> {
    let f = std::fs::File::open(path).map_err(|e| format!("open {}: {e}", path.display()))?;
    let reader = BufReader::new(f);
    let mut out = Vec::new();
    for line in reader.lines() {
        let line = line.map_err(|e| format!("read {}: {e}", path.display()))?;
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 4 {
            return Err(format!("expected contig\\tstart\\tend\\tis_indel: {line}"));
        }
        let start: u64 = cols[1]
            .parse()
            .map_err(|_| format!("invalid start: {}", cols[1]))?;
        let end: u64 = cols[2]
            .parse()
            .map_err(|_| format!("invalid end: {}", cols[2]))?;
        let is_indel = matches!(cols[3].to_ascii_lowercase().as_str(), "true" | "1" | "yes");
        out.push(TrimVariant {
            contig: cols[0].to_string(),
            start,
            end,
            is_indel,
            ref_allele: String::new(),
            alt_allele: String::new(),
        });
    }
    Ok(out)
}

/// Gate TSV.
#[cfg(any(feature = "dev-dumps", test))]
pub fn dump_assembly_region_trim_tsv(
    region: &AssemblyRegion,
    result: &AssemblyRegionTrimResult,
    trimmed: &AssemblyRegion,
    out: &mut impl std::io::Write,
) -> Result<(), String> {
    writeln!(
        out,
        "contig\torig_start\torig_end\torig_ext_start\torig_ext_end\tvariation_present\tvariant_start\tvariant_end\ttrim_start\ttrim_end\ttrim_ext_start\ttrim_ext_end"
    )
    .map_err(|e| e.to_string())?;
    let (vs, ve, _ps, _pe) = if result.variation_present {
        (
            result
                .variant_start
                .map(|v| v.to_string())
                .unwrap_or_default(),
            result
                .variant_end
                .map(|v| v.to_string())
                .unwrap_or_default(),
            result
                .padded_variant_start
                .map(|v| v.to_string())
                .unwrap_or_default(),
            result
                .padded_variant_end
                .map(|v| v.to_string())
                .unwrap_or_default(),
        )
    } else {
        ("-".into(), "-".into(), "-".into(), "-".into())
    };
    writeln!(
        out,
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        region.contig,
        region.start.get(),
        region.end.get(),
        region.extended_start.get(),
        region.extended_end.get(),
        result.variation_present,
        vs,
        ve,
        trimmed.start.get(),
        trimmed.end.get(),
        trimmed.extended_start.get(),
        trimmed.extended_end.get(),
    )
    .map_err(|e| e.to_string())?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::reference_context::ReferenceContext;

    fn region() -> AssemblyRegion {
        AssemblyRegion {
            contig: "chr1".into(),
            start: GenomePosition::new_1based(5),
            end: GenomePosition::new_1based(15),
            is_active: true,
            extended_start: GenomePosition::new_1based(1),
            extended_end: GenomePosition::new_1based(32),
            extension: 100,
            reads: Vec::new(),
            read_qnames: Vec::new(),
            reference: ReferenceContext::empty(),
            features: crate::feature_context::FeatureContext::empty(),
            pileup_loci: Vec::new(),
        }
    }

    #[test]
    fn no_variants_yields_no_variation() {
        let trimmer = AssemblyRegionTrimmer::new(
            AssemblyRegionTrimmerConfig::gatk_defaults(),
            &SequenceDictionary::default(),
            "chr1",
        );
        let r = region();
        let res = trimmer.trim(&r, &[], None);
        assert!(!res.variation_present);
    }

    #[test]
    fn modern_snp_trim_narrows_active_and_padded() {
        let mut dict = SequenceDictionary::new();
        dict.add_contig("chr1".into(), 32);
        let trimmer =
            AssemblyRegionTrimmer::new(AssemblyRegionTrimmerConfig::gatk_defaults(), &dict, "chr1");
        let r = region();
        let vars = vec![TrimVariant {
            contig: "chr1".into(),
            start: 10,
            end: 10,
            is_indel: false,
            ref_allele: String::new(),
            alt_allele: String::new(),
        }];
        let res = trimmer.trim(&r, &vars, None);
        assert!(res.variation_present);
        assert_eq!(res.variant_start, Some(10));
        assert_eq!(res.variant_end, Some(10));
        assert_eq!(res.padded_variant_start, Some(1));
        assert_eq!(res.padded_variant_end, Some(30));
        let trimmed = AssemblyRegionTrimmer::apply_trim(&r, &res);
        assert_eq!(trimmed.start.get(), 10);
        assert_eq!(trimmed.end.get(), 10);
        assert_eq!(trimmed.extended_start.get(), 1);
        assert_eq!(trimmed.extended_end.get(), 30);
    }

    /// GATK 4.4 `maxEnd = Math.max(maxEnd, vc.getEnd() + padding)`: three SNPs with
    /// pad 20 yield last_end+20, not last_end+20+20+20 (6R.38 Class C / 6R.39 fix).
    #[test]
    fn modern_three_snps_pad_is_max_not_sum() {
        let mut dict = SequenceDictionary::new();
        dict.add_contig("2".into(), 243_199_373);
        let trimmer =
            AssemblyRegionTrimmer::new(AssemblyRegionTrimmerConfig::gatk_defaults(), &dict, "2");
        let r = AssemblyRegion {
            contig: "2".into(),
            start: GenomePosition::new_1based(92_317_262),
            end: GenomePosition::new_1based(92_317_491),
            is_active: true,
            extended_start: GenomePosition::new_1based(92_317_162),
            extended_end: GenomePosition::new_1based(92_317_591),
            extension: 100,
            reads: Vec::new(),
            read_qnames: Vec::new(),
            reference: ReferenceContext::empty(),
            features: crate::feature_context::FeatureContext::empty(),
            pileup_loci: Vec::new(),
        };
        let vars = vec![
            TrimVariant {
                contig: "2".into(),
                start: 92_317_399,
                end: 92_317_399,
                is_indel: false,
                ref_allele: String::new(),
                alt_allele: String::new(),
            },
            TrimVariant {
                contig: "2".into(),
                start: 92_317_407,
                end: 92_317_407,
                is_indel: false,
                ref_allele: String::new(),
                alt_allele: String::new(),
            },
            TrimVariant {
                contig: "2".into(),
                start: 92_317_412,
                end: 92_317_412,
                is_indel: false,
                ref_allele: String::new(),
                alt_allele: String::new(),
            },
        ];
        let res = trimmer.trim(&r, &vars, None);
        assert!(res.variation_present);
        assert_eq!(res.variant_start, Some(92_317_399));
        assert_eq!(res.variant_end, Some(92_317_412));
        assert_eq!(res.padded_variant_start, Some(92_317_379));
        let per_event = [92_317_399 + 20, 92_317_407 + 20, 92_317_412 + 20];
        let java_max = *per_event.iter().max().unwrap();
        assert_eq!(java_max, 92_317_432);
        assert_eq!(res.padded_variant_end, Some(java_max));
        assert_ne!(
            res.padded_variant_end,
            Some(92_317_412 + 20 + 20 + 20),
            "must not accumulate SNP pad per event"
        );
    }

    /// GATK 4.4 `AssemblyRegionTrimmer.trim`: an overlapping deletion event with
    /// indel padding 75 pulls the padded span back across the D. SNP-only events
    /// after the D do not (6R.45/6R.46: that window starts inside the D).
    #[test]
    fn modern_indel_event_extends_padded_span_across_deletion() {
        let mut dict = SequenceDictionary::new();
        dict.add_contig("chr1".into(), 400);
        let trimmer =
            AssemblyRegionTrimmer::new(AssemblyRegionTrimmerConfig::gatk_defaults(), &dict, "chr1");
        let r = AssemblyRegion {
            contig: "chr1".into(),
            start: GenomePosition::new_1based(101),
            end: GenomePosition::new_1based(259),
            is_active: true,
            extended_start: GenomePosition::new_1based(1),
            extended_end: GenomePosition::new_1based(359),
            extension: 100,
            reads: Vec::new(),
            read_qnames: Vec::new(),
            reference: ReferenceContext::empty(),
            features: crate::feature_context::FeatureContext::empty(),
            pileup_loci: Vec::new(),
        };
        let snps = vec![
            TrimVariant {
                contig: "chr1".into(),
                start: 211,
                end: 211,
                is_indel: false,
                ref_allele: String::new(),
                alt_allele: String::new(),
            },
            TrimVariant {
                contig: "chr1".into(),
                start: 212,
                end: 212,
                is_indel: false,
                ref_allele: String::new(),
                alt_allele: String::new(),
            },
            TrimVariant {
                contig: "chr1".into(),
                start: 230,
                end: 230,
                is_indel: false,
                ref_allele: String::new(),
                alt_allele: String::new(),
            },
            TrimVariant {
                contig: "chr1".into(),
                start: 247,
                end: 247,
                is_indel: false,
                ref_allele: String::new(),
                alt_allele: String::new(),
            },
        ];
        let snp_only = trimmer.trim(&r, &snps, None);
        assert_eq!(snp_only.padded_variant_start, Some(191));
        assert_eq!(snp_only.padded_variant_end, Some(267));

        let mut with_del = snps;
        with_del.insert(
            0,
            TrimVariant {
                contig: "chr1".into(),
                start: 28,
                end: 199,
                is_indel: true,
                ref_allele: String::new(),
                alt_allele: String::new(),
            },
        );
        let with = trimmer.trim(&r, &with_del, None);
        assert_eq!(with.padded_variant_start, Some(1));
        assert_eq!(with.padded_variant_end, Some(274));
        assert_eq!(274, 28 + 171 + 75);
    }

    fn ref_at(start: u64, bases: &[u8]) -> ReferenceContext {
        ReferenceContext {
            contig: "chr1".into(),
            start,
            end: start + bases.len() as u64 - 1,
            window_start: start,
            window_end: start + bases.len() as u64 - 1,
            bases: crate::reference_context::SharedBases::from_slice(bases),
        }
    }

    /// Anchor `A` is not part of the following `T` run.
    #[test]
    fn tandem_repeat_excludes_anchor() {
        let rep = tandem_repeat_at_event(100, b"CATTTTT", 101, b"A", b"AT").expect("repeat");
        assert_eq!(rep.unit, b"T");
        assert_eq!(rep.counts, vec![5, 6]);
        assert_ne!(rep.unit, b"A");
    }

    /// The inserted bases extend the repeat past the reference run.
    #[test]
    fn tandem_repeat_counts_alternate_allele() {
        let rep = tandem_repeat_at_event(100, b"CATTT", 101, b"A", b"ATTT").expect("repeat");
        assert_eq!(rep.unit, b"T");
        assert_eq!(rep.counts, vec![3, 6]);
        assert_eq!(rep.padding_bases(), 6);
    }

    #[test]
    fn non_str_indel_keeps_indel_padding() {
        assert!(tandem_repeat_at_event(100, b"CATGC", 101, b"A", b"AG").is_none());
        let mut dict = SequenceDictionary::new();
        dict.add_contig("chr1".into(), 400);
        let trimmer =
            AssemblyRegionTrimmer::new(AssemblyRegionTrimmerConfig::gatk_defaults(), &dict, "chr1");
        let mut r = region();
        r.start = GenomePosition::new_1based(90);
        r.end = GenomePosition::new_1based(120);
        r.extended_start = GenomePosition::new_1based(1);
        r.extended_end = GenomePosition::new_1based(400);
        let vars = vec![TrimVariant {
            contig: "chr1".into(),
            start: 101,
            end: 101,
            is_indel: true,
            ref_allele: "A".into(),
            alt_allele: "AG".into(),
        }];
        let res = trimmer.trim(&r, &vars, Some(&ref_at(100, b"CATGC")));
        assert_eq!(res.padded_variant_start, Some(101 - 75));
        assert_eq!(res.padded_variant_end, Some(101 + 75));
    }

    #[test]
    fn snp_padding_ignores_following_repeat() {
        let mut dict = SequenceDictionary::new();
        dict.add_contig("chr1".into(), 400);
        let trimmer =
            AssemblyRegionTrimmer::new(AssemblyRegionTrimmerConfig::gatk_defaults(), &dict, "chr1");
        let mut r = region();
        r.start = GenomePosition::new_1based(50);
        r.end = GenomePosition::new_1based(200);
        r.extended_start = GenomePosition::new_1based(1);
        r.extended_end = GenomePosition::new_1based(400);
        let vars = vec![TrimVariant {
            contig: "chr1".into(),
            start: 101,
            end: 101,
            is_indel: false,
            ref_allele: "A".into(),
            alt_allele: "G".into(),
        }];
        let res = trimmer.trim(&r, &vars, Some(&ref_at(100, b"CATTTTTTTTT")));
        assert_eq!(res.padded_variant_start, Some(81));
        assert_eq!(res.padded_variant_end, Some(121));
    }
}
