//! SAM/BAM header semantics (RG/PG) validation for user inputs.

#![warn(clippy::unwrap_used, clippy::expect_used)]

use gatk_common::{GatkError, GatkResult};
use rust_htslib::bam;
use std::collections::{BTreeMap, BTreeSet};

/// Parsed SAM/BAM header semantics for read-group and program validation.
/// # Invariants
/// Each `@RG ID` maps to at most one optional `SM` sample name.
/// `@PG ID` values are unique in the parsed header.
/// # Ownership
/// Owns RG→sample map and PG id set; built once per header text/view.
/// # Mutation
/// Immutable after construction; [`Self::validate_record_links`] returns per-read resolution.
/// # Biological assumptions
/// Sample identity comes from `@RG SM:` tags (multi-sample vs single-sample HC paths).
/// # Java equivalence
/// GATK read header validation / `ReadUtils.alignmentAgreesWithHeader` inputs.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ReadHeaderSemantics {
    read_group_to_sample: BTreeMap<String, Option<String>>,
    program_ids: BTreeSet<String>,
}

/// Resolved `@RG` / `@PG` links for one read record against a parsed header.
/// # Invariants
/// When `read_group_id` is `Some`, it existed in the header map used to resolve `sample_name`.
/// When `program_id` is `Some`, it was present in the header `@PG` set.
/// # Ownership
/// Owns optional tag strings; produced by [`ReadHeaderSemantics::validate_record_links`].
/// # Mutation
/// Immutable per-read snapshot.
/// # Biological assumptions
/// Sample name derives from read group's `SM` when declared.
/// # Java equivalence
/// GATK header-aware read validation (`SAMRecord` RG/PG linkage).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ResolvedReadHeaderSemantics {
    pub read_group_id: Option<String>,
    pub sample_name: Option<String>,
    pub program_id: Option<String>,
}

impl ReadHeaderSemantics {
    pub fn from_sam_header_text(header_text: &str) -> GatkResult<Self> {
        let mut read_group_to_sample = BTreeMap::new();
        let mut program_ids = BTreeSet::new();

        for line in header_text.lines() {
            if line.starts_with("@RG\t") {
                let fields = parse_key_value_fields(line);
                let rg_id = fields.get("ID").cloned().ok_or_else(|| {
                    GatkError::validation("Malformed SAM header: @RG line missing ID")
                })?;
                let sample = fields.get("SM").cloned();
                // Lifetime: detect duplicates with a borrow, then move `rg_id` into the map.
                if read_group_to_sample.contains_key(&rg_id) {
                    return Err(GatkError::validation(format!(
                        "Malformed SAM header: duplicate @RG ID '{rg_id}'"
                    )));
                }
                read_group_to_sample.insert(rg_id, sample);
            } else if line.starts_with("@PG\t") {
                let fields = parse_key_value_fields(line);
                let pg_id = fields.get("ID").cloned().ok_or_else(|| {
                    GatkError::validation("Malformed SAM header: @PG line missing ID")
                })?;
                // Lifetime: `contains` borrows; move `pg_id` into the set on success.
                if program_ids.contains(&pg_id) {
                    return Err(GatkError::validation(format!(
                        "Malformed SAM header: duplicate @PG ID '{pg_id}'"
                    )));
                }
                program_ids.insert(pg_id);
            }
        }

        Ok(Self {
            read_group_to_sample,
            program_ids,
        })
    }

    pub fn from_bam_header_view(header: &bam::HeaderView) -> GatkResult<Self> {
        let text = std::str::from_utf8(header.as_bytes()).map_err(|e| {
            GatkError::validation(format!("Malformed BAM header: non-UTF8 header text ({e})"))
        })?;
        Self::from_sam_header_text(text)
    }

    pub fn sample_for_read_group(&self, rg_id: &str) -> Option<&str> {
        self.read_group_to_sample
            .get(rg_id)
            .and_then(|s| s.as_deref())
    }

    /// Distinct `@RG` SM tags (GATK `samplesList.numberOfSamples` for HC `isActive` shortcut).
    pub fn unique_sample_count(&self) -> usize {
        self.read_group_to_sample
            .values()
            .filter_map(|s| s.as_deref())
            .collect::<BTreeSet<_>>()
            .len()
    }

    /// Single `@RG SM:` when exactly one distinct sample is declared; else `None`.
    pub fn primary_sample_name(&self) -> Option<String> {
        let samples: BTreeSet<_> = self
            .read_group_to_sample
            .values()
            .filter_map(|s| s.as_deref())
            .collect();
        if samples.len() == 1 {
            samples.into_iter().next().map(|s| s.to_string())
        } else {
            None
        }
    }

    /// True when the header declares at most one sample (Java single-sample `isActive` fast path).
    pub fn is_single_sample_header(&self) -> bool {
        self.unique_sample_count() <= 1
    }

    pub fn validate_record_links(
        &self,
        rg_tag: Option<&str>,
        pg_tag: Option<&str>,
    ) -> GatkResult<ResolvedReadHeaderSemantics> {
        let sample_name = match rg_tag {
            Some(rg) => {
                let sample = self.read_group_to_sample.get(rg).ok_or_else(|| {
                    GatkError::validation(format!(
                        "Read header mismatch: record RG '{rg}' not found in header @RG IDs"
                    ))
                })?;
                // CLONE: needed because owned sample id for carry/map.
                sample.clone()
            }
            None => None,
        };

        if let Some(pg) = pg_tag {
            if !self.program_ids.contains(pg) {
                return Err(GatkError::validation(format!(
                    "Read header mismatch: record PG '{pg}' not found in header @PG IDs"
                )));
            }
        }

        Ok(ResolvedReadHeaderSemantics {
            read_group_id: rg_tag.map(ToOwned::to_owned),
            sample_name,
            program_id: pg_tag.map(ToOwned::to_owned),
        })
    }
}

fn parse_key_value_fields(line: &str) -> BTreeMap<String, String> {
    let mut fields = BTreeMap::new();
    for field in line.split('\t').skip(1) {
        if let Some((k, v)) = field.split_once(':') {
            fields.insert(k.to_string(), v.to_string());
        }
    }
    fields
}

#[cfg(test)]
#[allow(clippy::unwrap_used, clippy::expect_used, clippy::panic)]
mod tests {
    use super::*;

    #[test]
    fn parses_rg_pg_and_resolves_sample() {
        let header = "@HD\tVN:1.6\n@RG\tID:rg1\tSM:sampleA\n@PG\tID:pg1\tPN:tool\n";
        let semantics = ReadHeaderSemantics::from_sam_header_text(header).unwrap();
        assert_eq!(semantics.sample_for_read_group("rg1"), Some("sampleA"));
        let resolved = semantics
            .validate_record_links(Some("rg1"), Some("pg1"))
            .unwrap();
        assert_eq!(resolved.sample_name.as_deref(), Some("sampleA"));
    }

    #[test]
    fn rejects_duplicate_rg_id() {
        let header = "@RG\tID:rg1\tSM:s1\n@RG\tID:rg1\tSM:s2\n";
        let err = ReadHeaderSemantics::from_sam_header_text(header).unwrap_err();
        match err {
            GatkError::Validation { message, .. } => assert!(message.contains("duplicate @RG ID")),
            other => panic!("expected Validation error, got {other:?}"),
        }
    }

    #[test]
    fn rejects_unknown_rg_and_pg_links() {
        let header = "@RG\tID:rg1\tSM:s1\n@PG\tID:pg1\n";
        let semantics = ReadHeaderSemantics::from_sam_header_text(header).unwrap();
        let rg_err = semantics
            .validate_record_links(Some("missing-rg"), Some("pg1"))
            .unwrap_err();
        let pg_err = semantics
            .validate_record_links(Some("rg1"), Some("missing-pg"))
            .unwrap_err();
        match rg_err {
            GatkError::Validation { message, .. } => assert!(message.contains("record RG")),
            other => panic!("expected Validation error, got {other:?}"),
        }
        match pg_err {
            GatkError::Validation { message, .. } => assert!(message.contains("record PG")),
            other => panic!("expected Validation error, got {other:?}"),
        }
    }
}
