//! GATK `FeatureContext` for assembly-region `apply`.
//! Mirrors `new FeatureContext(featureManager, assemblyRegion.getPaddedSpan)`. HC defaults use
//! no backing resources (`FeatureManager` null); optional VCF paths supply forced / known alleles.

use gatk_common::{GatkError, GatkResult};
use gatk_core::io::vcf::{VcfReader, VcfRecord};
use std::collections::BTreeMap;
use std::path::Path;

/// One locatable feature overlapping a query interval (minimal VCF-shaped record).
/// # Invariants
/// `start` / `end` are 1-based inclusive; `end` derived from REF/ALT lengths when from VCF.
/// `passes_filters` is true for empty/`.`/`PASS` FILTER sets.
/// # Ownership
/// Owns source name, contig, and allele strings.
/// # Mutation
/// Immutable feature row after load.
/// # Biological assumptions
/// Known/forced alleles or resource variants overlapping an assembly region.
/// # Java equivalence
/// GATK `Feature` / VCF-backed locatable for `FeatureContext`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FeatureLocatable {
    pub source: String,
    pub contig: String,
    pub start: u64,
    pub end: u64,
    pub reference: String,
    pub alternates: Vec<String>,
    pub passes_filters: bool,
}

impl FeatureLocatable {
    pub fn from_vcf_record(source: &str, rec: &VcfRecord) -> Self {
        let end = variant_end_1based(rec.position, &rec.reference, &rec.alternate);
        let passes_filters = rec.filter.is_empty()
            || rec
                .filter
                .iter()
                .all(|f| f == "." || f.eq_ignore_ascii_case("PASS"));
        Self {
            source: source.to_string(),
            // CLONE: needed because owned contig id for output record.
            contig: rec.chromosome.clone(),
            start: rec.position,
            end,
            reference: rec.reference.clone(),
            alternates: rec.alternate.clone(),
            passes_filters,
        }
    }

    pub fn overlaps_closed_interval_1based(&self, start1: u64, end1: u64) -> bool {
        self.start <= end1 && self.end >= start1
    }
}

fn variant_end_1based(pos1: u64, reference: &str, alternates: &[String]) -> u64 {
    let mut max_len = reference.len();
    for alt in alternates {
        max_len = max_len.max(alt.len());
    }
    pos1.saturating_add(max_len.saturating_sub(1) as u64)
}

/// Named feature resources (VCF paths) indexed for interval queries.
/// # Invariants
/// Features within a source are sorted by contig/start/end/REF for stable queries.
/// Force-calling overlap uses first-match semantics on the named source.
/// # Ownership
/// Owns named vectors of [`FeatureLocatable`].
/// # Mutation
/// [`Self::load_vcf_source`] inserts/replaces a named source; queries borrow immutably.
/// # Biological assumptions
/// Optional `-alleles` / known-sites resources for force-active and given-allele paths.
/// # Java equivalence
/// GATK `FeatureManager` / resource feature sources for HC.
#[derive(Debug, Clone, Default)]
pub struct FeatureDataSources {
    sources: BTreeMap<String, Vec<FeatureLocatable>>,
}

impl FeatureDataSources {
    pub fn is_empty(&self) -> bool {
        self.sources.is_empty()
    }

    pub fn load_vcf_source(&mut self, name: impl Into<String>, path: &Path) -> GatkResult<()> {
        let name = name.into();
        let mut reader = VcfReader::from_file(path)
            .map_err(|e| GatkError::generic(format!("open VCF {}: {e}", path.display())))?;
        let records = reader
            .read_all_records()
            .map_err(|e| GatkError::generic(format!("read VCF {}: {e}", path.display())))?;
        let mut feats: Vec<FeatureLocatable> = records
            .iter()
            .map(|r| FeatureLocatable::from_vcf_record(&name, r))
            .collect();
        feats.sort_by(|a, b| {
            a.contig
                .cmp(&b.contig)
                .then(a.start.cmp(&b.start))
                .then(a.end.cmp(&b.end))
                .then(a.reference.cmp(&b.reference))
        });
        self.sources.insert(name, feats);
        Ok(())
    }

    pub fn features_overlapping(
        &self,
        contig: &str,
        start1: u64,
        end1: u64,
    ) -> Vec<FeatureLocatable> {
        let mut out = Vec::new();
        for feats in self.sources.values() {
            for f in feats {
                if f.contig == contig && f.overlaps_closed_interval_1based(start1, end1) {
                    // CLONE: needed because owned element into collection.
                    out.push(f.clone());
                }
            }
        }
        out.sort_by(|a, b| {
            a.source
                .cmp(&b.source)
                .then(a.start.cmp(&b.start))
                .then(a.reference.cmp(&b.reference))
        });
        out
    }

    /// `HaplotypeCallerEngine#isActive`: first-match forced alleles from `--alleles` (and friends).
    pub fn force_calling_allele_overlaps_locus(
        &self,
        source_name: &str,
        contig: &str,
        pos1: u64,
        force_call_filtered: bool,
    ) -> bool {
        let Some(feats) = self.sources.get(source_name) else {
            return false;
        };
        feats.iter().any(|f| {
            f.contig == contig
                && f.overlaps_closed_interval_1based(pos1, pos1)
                && (force_call_filtered || f.passes_filters)
        })
    }
}

/// Feature data for one assembly region query interval (typically padded span).
/// # Invariants
/// Query span `start`/`end` is 1-based inclusive (typically padded assembly region).
/// HC default has `has_backing_data_source == false` and empty features.
/// # Ownership
/// Owns contig, span, and overlapping feature list for the region.
/// # Mutation
/// Built once when the assembly region is materialized; treated as immutable in apply.
/// # Biological assumptions
/// Features overlapping the padded span available to `callRegion` / force-calling.
/// # Java equivalence
/// GATK `FeatureContext(featureManager, assemblyRegion.getPaddedSpan)`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FeatureContext {
    pub contig: String,
    pub start: u64,
    pub end: u64,
    pub has_backing_data_source: bool,
    pub features: Vec<FeatureLocatable>,
}

impl FeatureContext {
    pub fn empty() -> Self {
        Self {
            contig: String::new(),
            start: 0,
            end: 0,
            has_backing_data_source: false,
            features: Vec::new(),
        }
    }

    /// HC default: `new FeatureContext(null, interval)` — interval only, no features.
    pub fn without_sources(contig: impl Into<String>, start1: u64, end1: u64) -> Self {
        Self {
            contig: contig.into(),
            start: start1,
            end: end1,
            has_backing_data_source: false,
            features: Vec::new(),
        }
    }

    pub fn from_padded_span(
        contig: &str,
        start1: u64,
        end1: u64,
        sources: Option<&FeatureDataSources>,
    ) -> Self {
        match sources {
            None => Self::without_sources(contig, start1, end1),
            Some(s) if s.is_empty() => Self::without_sources(contig, start1, end1),
            Some(s) => {
                let features = s.features_overlapping(contig, start1, end1);
                Self {
                    contig: contig.to_string(),
                    start: start1,
                    end: end1,
                    has_backing_data_source: true,
                    features,
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    #[test]
    fn without_sources_is_empty() {
        let ctx = FeatureContext::without_sources("chr1", 5, 15);
        assert!(!ctx.has_backing_data_source);
        assert!(ctx.features.is_empty());
    }

    #[test]
    fn force_calling_detects_overlapping_pass_site() {
        let mut s = FeatureDataSources::default();
        let dir = std::env::temp_dir().join(format!("fc_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let vcf = dir.join("a.vcf");
        {
            let mut f = std::fs::File::create(&vcf).unwrap();
            use std::io::Write;
            writeln!(f, "##fileformat=VCFv4.2").unwrap();
            writeln!(f, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO").unwrap();
            writeln!(f, "chr1\t10\t.\tC\tG\t60\tPASS\t.").unwrap();
        }
        s.load_vcf_source("alleles", &vcf).unwrap();
        assert!(s.force_calling_allele_overlaps_locus("alleles", "chr1", 10, false));
        assert!(!s.force_calling_allele_overlaps_locus("alleles", "chr1", 11, false));
        let _ = std::fs::remove_dir_all(dir);
    }

    #[test]
    fn vcf_source_yields_overlapping_feature() {
        let dir = std::env::temp_dir().join(format!("b54_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let vcf = dir.join("alleles.vcf");
        {
            let mut f = std::fs::File::create(&vcf).unwrap();
            writeln!(f, "##fileformat=VCFv4.2").unwrap();
            writeln!(f, "##contig=<ID=chr1,length=32>").unwrap();
            writeln!(f, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO").unwrap();
            writeln!(f, "chr1\t10\t.\tA\tG\t60\tPASS\t.").unwrap();
        }
        let mut sources = FeatureDataSources::default();
        sources.load_vcf_source("alleles", &vcf).unwrap();
        let ctx = FeatureContext::from_padded_span("chr1", 1, 32, Some(&sources));
        assert!(ctx.has_backing_data_source);
        assert_eq!(ctx.features.len(), 1);
        assert_eq!(ctx.features[0].start, 10);
        let _ = std::fs::remove_dir_all(dir);
    }
}
