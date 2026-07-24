//! Parity v1 site INFO annotations (AC/AN/AF/NS/DP).
//! Free-function API — no Java-style plugin engine until a second annotator exists.

use crate::genotyping::{compute_core_variant_annotations, SampleAnnotationInput};
use gatk_common::GatkResult;

/// INFO keys emitted by parity v1 (`i1-core` gate).
pub const PARITY_V1_INFO_KEYS: &[&str] = &["AC", "AN", "AF", "NS", "DP"];

/// FORMAT keys agreed for parity v1 (computed in; listed for manifest completeness).
pub const PARITY_V1_FORMAT_KEYS: &[&str] = &["GT", "GQ", "AD", "DP", "PL"];

/// Input for parity v1 site annotation.
/// # Invariants
/// `alt_allele_count` matches the number of ALT alleles; sample genotypes use VCF allele indices.
/// # Ownership
/// Owns sample annotation inputs for the site.
/// # Mutation
/// Immutable context per annotate call.
/// # Biological assumptions
/// Called genotypes already assigned; annotators aggregate AC/AN/AF/NS/DP-style INFO.
/// # Java equivalence
/// GATK annotation context inputs for the HC default core INFO subset.
#[derive(Debug, Clone, PartialEq)]
pub struct VariantAnnotationContext {
    pub alt_allele_count: usize,
    pub samples: Vec<SampleAnnotationInput>,
}

/// Annotated site after parity v1 core INFO computation.
/// # Invariants
/// Core counts match [`crate::genotyping::CoreVariantAnnotations`]; key lists enumerate emitted INFO/FORMAT.
/// # Ownership
/// Owns annotation values and key name lists.
/// # Mutation
/// Immutable annotate result.
/// # Biological assumptions
/// Site-level summary over called samples at one variant locus.
/// # Java equivalence
/// GATK annotated `VariantContext` INFO fields for the HC default core subset.
#[derive(Debug, Clone, PartialEq)]
pub struct AnnotatedSite {
    pub ac: Vec<i32>,
    pub an: i32,
    pub af: Vec<f64>,
    pub ns: i32,
    pub dp: i32,
    pub info_keys: Vec<String>,
    pub format_keys: Vec<String>,
}

/// Compute parity v1 core INFO (AC/AN/AF/NS/DP) for one site.
pub fn annotate_parity_v1_site(ctx: &VariantAnnotationContext) -> GatkResult<AnnotatedSite> {
    let core = compute_core_variant_annotations(ctx.alt_allele_count, &ctx.samples)?;
    Ok(AnnotatedSite {
        ac: core.ac,
        an: core.an,
        af: core.af,
        ns: core.ns,
        dp: core.dp,
        info_keys: PARITY_V1_INFO_KEYS
            .iter()
            .map(|s| (*s).to_string())
            .collect(),
        format_keys: PARITY_V1_FORMAT_KEYS
            .iter()
            .map(|s| (*s).to_string())
            .collect(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::genotyping::SampleAnnotationInput;

    #[test]
    fn parity_v1_matches_core_annotations() {
        let ctx = VariantAnnotationContext {
            alt_allele_count: 1,
            samples: vec![
                SampleAnnotationInput {
                    genotype_alleles: vec![0, 0],
                    dp: Some(10),
                },
                SampleAnnotationInput {
                    genotype_alleles: vec![0, 1],
                    dp: Some(8),
                },
            ],
        };
        let site = annotate_parity_v1_site(&ctx).expect("annotate");
        assert_eq!(site.ac, vec![1]);
        assert_eq!(site.an, 4);
        assert_eq!(site.ns, 2);
        assert_eq!(site.dp, 18);
        assert!((site.af[0] - 0.25).abs() < 1e-9);
    }
}
