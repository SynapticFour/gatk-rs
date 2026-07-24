use crate::parallel::ParallelConfig;
use polars::prelude::*;

/// Tabular variant row for Polars dataframe export.
/// # Invariants
/// One row per variant; chromosome/position uniquely identify row in typical frames.
/// # Ownership
/// Owns string fields; clone for dataframe building.
/// # Mutation
/// Public fields for ETL staging.
/// # Biological assumptions
/// VCF-like variant summary with depth and allele frequency.
/// # Java equivalence
/// Approximates htsjdk `VariantContext` columns in a dataframe (Rust-native).
#[derive(Debug, Clone)]
pub struct VariantData {
    pub chromosome: String,
    pub position: u64,
    pub id: String,
    pub reference: String,
    pub alternate: String,
    pub quality: f64,
    pub filter: String,
    pub depth: u32,
    pub allele_frequency: f64,
}

/// High-level QC summary computed from variant dataframe analytics.
/// # Invariants
/// `ti_tv_ratio` meaningful only when transition/transversion counts are populated by caller.
/// # Ownership
/// Plain scalars; default zeroed.
/// # Mutation
/// Immutable summary output.
/// # Biological assumptions
/// Ti/Tv and high-quality counts used for SNP QC heuristics.
/// # Java equivalence
/// None / Rust-native analytics.
#[derive(Debug, Clone, Default)]
pub struct GenomicAnalysisSummary {
    pub high_quality_count: usize,
    pub ti_tv_ratio: f64,
}

/// Polars dataframe builder for genomic variant tables.
/// # Invariants
/// Column schema fixed in `create_variant_dataframe` helpers.
/// # Ownership
/// Owns config; returns owned Polars `DataFrame` to caller.
/// # Mutation
/// Methods borrow `&self`; dataframe construction allocates new frames.
/// # Biological assumptions
/// Variant rows follow VCF semantics (CHR/POS/REF/ALT/QUAL/DP/AF).
/// # Java equivalence
/// None / Rust-native (Polars replaces Spark SQL for local analytics).
pub struct PolarsProcessor {
    _config: ParallelConfig,
}

impl PolarsProcessor {
    pub fn new(config: ParallelConfig) -> gatk_common::GatkResult<Self> {
        Ok(Self { _config: config })
    }

    pub fn create_variant_dataframe(
        &self,
        variants: Vec<VariantData>,
    ) -> gatk_common::GatkResult<DataFrame> {
        let chr = Series::new(
            "chromosome",
            variants
                .iter()
                .map(|v| v.chromosome.as_str())
                .collect::<Vec<_>>(),
        );
        let pos = Series::new(
            "position",
            variants.iter().map(|v| v.position).collect::<Vec<_>>(),
        );
        let id = Series::new(
            "id",
            variants.iter().map(|v| v.id.as_str()).collect::<Vec<_>>(),
        );
        let r = Series::new(
            "reference",
            variants
                .iter()
                .map(|v| v.reference.as_str())
                .collect::<Vec<_>>(),
        );
        let a = Series::new(
            "alternate",
            variants
                .iter()
                .map(|v| v.alternate.as_str())
                .collect::<Vec<_>>(),
        );
        let q = Series::new(
            "quality",
            variants.iter().map(|v| v.quality).collect::<Vec<_>>(),
        );
        let filt = Series::new(
            "filter",
            variants
                .iter()
                .map(|v| v.filter.as_str())
                .collect::<Vec<_>>(),
        );
        let dp = Series::new(
            "depth",
            variants.iter().map(|v| v.depth).collect::<Vec<_>>(),
        );
        let af = Series::new(
            "allele_frequency",
            variants
                .iter()
                .map(|v| v.allele_frequency)
                .collect::<Vec<_>>(),
        );
        DataFrame::new(vec![chr, pos, id, r, a, q, filt, dp, af])
            .map_err(|e| gatk_common::GatkError::generic(format!("Failed to build DataFrame: {e}")))
    }

    pub fn filter_variants_by_quality(
        &self,
        df: &DataFrame,
        min_quality: f64,
        min_depth: u32,
    ) -> gatk_common::GatkResult<DataFrame> {
        let out = df
            .clone()
            .lazy()
            .filter(
                col("quality")
                    .gt_eq(lit(min_quality))
                    .and(col("depth").gt_eq(lit(min_depth))),
            )
            .collect()
            .map_err(|e| {
                gatk_common::GatkError::generic(format!("Failed to filter variants: {e}"))
            })?;
        Ok(out)
    }

    pub fn group_by_chromosome(&self, df: &DataFrame) -> gatk_common::GatkResult<DataFrame> {
        df.clone()
            .lazy()
            .group_by([col("chromosome")])
            .agg([col("position").count().alias("count")])
            .collect()
            .map_err(|e| gatk_common::GatkError::generic(format!("Failed to group variants: {e}")))
    }

    pub fn calculate_ti_tv_ratio(&self, df: &DataFrame) -> gatk_common::GatkResult<f64> {
        let ref_col = df.column("reference").map_err(|e| {
            gatk_common::GatkError::generic(format!("Missing reference column: {e}"))
        })?;
        let alt_col = df.column("alternate").map_err(|e| {
            gatk_common::GatkError::generic(format!("Missing alternate column: {e}"))
        })?;

        let refs = ref_col
            .str()
            .map_err(|e| gatk_common::GatkError::generic(format!("{e}")))?;
        let alts = alt_col
            .str()
            .map_err(|e| gatk_common::GatkError::generic(format!("{e}")))?;

        let mut transitions = 0usize;
        let mut transversions = 0usize;
        for (r, a) in refs.into_iter().zip(alts) {
            if let (Some(r), Some(a)) = (r, a) {
                let pair = (r.as_bytes().first().copied(), a.as_bytes().first().copied());
                match pair {
                    (Some(b'A'), Some(b'G'))
                    | (Some(b'G'), Some(b'A'))
                    | (Some(b'C'), Some(b'T'))
                    | (Some(b'T'), Some(b'C')) => transitions += 1,
                    (Some(_), Some(_)) => transversions += 1,
                    _ => {}
                }
            }
        }
        Ok(if transversions == 0 {
            transitions as f64
        } else {
            transitions as f64 / transversions as f64
        })
    }

    pub fn genomic_analysis_pipeline(
        &self,
        variants_df: &DataFrame,
        _reads_df: &DataFrame,
    ) -> gatk_common::GatkResult<GenomicAnalysisSummary> {
        let filtered = self.filter_variants_by_quality(variants_df, 30.0, 10)?;
        let ti_tv = self.calculate_ti_tv_ratio(&filtered)?;
        Ok(GenomicAnalysisSummary {
            high_quality_count: filtered.height(),
            ti_tv_ratio: ti_tv,
        })
    }
}
