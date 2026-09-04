//! Genotyping implementation

use crate::activity_scoring::normalize_from_log10_to_linear_space;
use gatk_common::{GatkError, GatkResult};

/// One read's log10 likelihoods across candidate haplotypes.
/// # Invariants
/// `haplotype_log10_likelihoods.len` is identical across rows in an aggregation batch.
/// # Ownership
/// Owns read id and likelihood vector.
/// # Mutation
/// Immutable input to [`aggregate_haplotype_log10_likelihoods`].
/// # Biological assumptions
/// Per-read haplotype likelihoods already in log10 space (PairHMM or test fixture).
/// # Java equivalence
/// GATK `ReadLikelihoods` row slice (genotyping aggregation primitive).
#[derive(Debug, Clone, PartialEq)]
pub struct ReadLikelihoodRow {
    /// Dense likelihood-matrix row index (production reshape fills this).
    pub read_index: usize,
    /// Diagnostic / fixture id. Production `region_likelihoods_to_rows` leaves this
    /// empty to avoid per-read `format!` in the genotype hot path.
    pub read_id: String,
    pub haplotype_log10_likelihoods: Vec<f64>,
}

impl ReadLikelihoodRow {
    /// Resolve matrix row index: prefer `read_index` when `read_id` is empty.
    #[inline]
    pub fn matrix_read_index(&self) -> Option<usize> {
        if self.read_id.is_empty() {
            return Some(self.read_index);
        }
        self.read_id
            .strip_prefix("read_")
            .and_then(|s| s.parse().ok())
            .or(Some(self.read_index))
    }
}

/// Aggregated log10 likelihood sums over all reads for each haplotype.
/// # Invariants
/// `haplotype_log10_sums.len` equals haplotype width; non-finite likelihoods skipped in the sum.
/// `read_count` is the number of input rows (including rows with all −∞).
/// # Ownership
/// Owns sum vector; returned from aggregation.
/// # Mutation
/// Immutable aggregation result.
/// # Biological assumptions
/// Sums approximate joint support for each haplotype across independent reads (log space).
/// # Java equivalence
/// Rust-native aggregation over GATK-style read×haplotype likelihood matrices.
#[derive(Debug, Clone, PartialEq)]
pub struct HaplotypeLikelihoodAggregation {
    pub haplotype_log10_sums: Vec<f64>,
    pub read_count: usize,
}

/// Aggregate read-level haplotype likelihood vectors into per-haplotype sums.
/// This is the first genotyping core primitive: we keep likelihoods in log10 space
/// and perform deterministic per-haplotype accumulation across reads.
pub fn aggregate_haplotype_log10_likelihoods(
    rows: &[ReadLikelihoodRow],
) -> GatkResult<HaplotypeLikelihoodAggregation> {
    if rows.is_empty() {
        return Ok(HaplotypeLikelihoodAggregation {
            haplotype_log10_sums: Vec::new(),
            read_count: 0,
        });
    }

    let hap_count = rows[0].haplotype_log10_likelihoods.len();
    if hap_count == 0 {
        return Err(GatkError::argument(
            "genotyping aggregation requires at least one haplotype likelihood per read",
        ));
    }

    let mut sums = vec![0.0_f64; hap_count];
    for row in rows {
        if row.haplotype_log10_likelihoods.len() != hap_count {
            return Err(GatkError::argument(format!(
                "genotyping aggregation haplotype width mismatch for read {}: got {} expected {}",
                row.read_id,
                row.haplotype_log10_likelihoods.len(),
                hap_count
            )));
        }
        for (idx, ll) in row.haplotype_log10_likelihoods.iter().copied().enumerate() {
            if !ll.is_finite() || ll == f64::NEG_INFINITY {
                continue;
            }
            sums[idx] += ll;
        }
    }

    Ok(HaplotypeLikelihoodAggregation {
        haplotype_log10_sums: sums,
        read_count: rows.len(),
    })
}

/// Pick the maximum-likelihood haplotype index from aggregated log10 sums.
/// Ties are resolved deterministically by the first (lowest index) haplotype.
pub fn best_haplotype_index(
    aggregation: &HaplotypeLikelihoodAggregation,
) -> Option<crate::bio_ids::HaplotypeIndex> {
    aggregation
        .haplotype_log10_sums
        .iter()
        .enumerate()
        .max_by(|a, b| a.1.total_cmp(b.1))
        .map(|(idx, _)| crate::bio_ids::HaplotypeIndex::new(idx))
}

/// Simple biallelic diploid prior model for genotype states [0/0, 0/1, 1/1].
/// # Invariants
/// `het_prior` and `hom_var_prior` ∈ \[0,1\] and leave positive mass for hom-ref.
/// # Ownership
/// [`Copy`] prior parameters (also nested in [`HcGenotypingConfig`](crate::hc_genotyping_engine::HcGenotypingConfig)).
/// # Mutation
/// Immutable per genotyping call.
/// # Biological assumptions
/// Diploid biallelic site with fixed het / hom-var priors (hom-ref is the remainder).
/// # Java equivalence
/// GATK genotype prior slice used with biallelic diploid PL calculation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BiallelicDiploidPriorModel {
    pub het_prior: f64,
    pub hom_var_prior: f64,
}

impl Default for BiallelicDiploidPriorModel {
    fn default() -> Self {
        Self {
            het_prior: 1e-3,
            hom_var_prior: 5e-4,
        }
    }
}

pub fn biallelic_diploid_log10_priors(model: BiallelicDiploidPriorModel) -> GatkResult<[f64; 3]> {
    if !(0.0..=1.0).contains(&model.het_prior) || !(0.0..=1.0).contains(&model.hom_var_prior) {
        return Err(GatkError::argument(
            "genotype prior model values must be in [0,1]",
        ));
    }
    let hom_ref_prior = 1.0 - model.het_prior - model.hom_var_prior;
    if hom_ref_prior <= 0.0 {
        return Err(GatkError::argument(
            "genotype prior model must leave positive probability for hom-ref",
        ));
    }
    Ok([
        hom_ref_prior.log10(),
        model.het_prior.log10(),
        model.hom_var_prior.log10(),
    ])
}

/// Normalized genotype posterior over diploid (or general) genotype indices.
/// # Invariants
/// `genotype_log10_posteriors` and `genotype_posteriors` have equal length and sum to ~1 in linear space.
/// `most_likely_genotype_index` is argmax of linear posteriors (ties: first index).
/// # Ownership
/// Owns posterior vectors; cheap to clone for emit pipelines.
/// # Mutation
/// Immutable result of [`genotype_posteriors_from_log10_likelihoods`].
/// # Biological assumptions
/// Genotype indices follow VCF diploid ordering for biallelic sites (0/0, 0/1, 1/1).
/// # Java equivalence
/// GATK `GenotypeLikelihoods` posterior normalization / best-genotype selection slice.
#[derive(Debug, Clone, PartialEq)]
pub struct GenotypePosterior {
    pub genotype_log10_posteriors: Vec<f64>,
    pub genotype_posteriors: Vec<f64>,
    pub most_likely_genotype_index: usize,
}

/// Emitted core genotype FORMAT fields for a diploid site.
/// # Invariants
/// `pl` length matches diploid genotype count for the allele set; `gq` ∈ [0, 99] after P7 capping.
/// `dp` equals sum of non-negative `ad` entries.
/// # Ownership
/// Owns PL/AD vectors and GQ/DP scalars.
/// # Mutation
/// Immutable emit snapshot from [`emit_genotype_format_fields`].
/// # Biological assumptions
/// Standard VCF FORMAT PL/GQ/AD/DP for one diploid sample at a site.
/// # Java equivalence
/// GATK / HTSJDK genotype FORMAT emission (`GenotypeBuilder` PL/GQ/AD/DP).
#[derive(Debug, Clone, PartialEq)]
pub struct GenotypeFormatFields {
    pub pl: Vec<crate::bio_ids::PhredLikelihood>,
    pub gq: crate::bio_ids::GenotypeQuality,
    pub ad: Vec<crate::bio_ids::AlleleDepth>,
    pub dp: crate::bio_ids::ReadDepth,
}

impl GenotypeFormatFields {
    /// Build from VCF/Java wire integers (negatives saturate to 0 for AD/DP/GQ/PL).
    pub fn from_wire(pl: Vec<i32>, gq: i32, ad: Vec<i32>, dp: i32) -> Self {
        use crate::bio_ids::{AlleleDepth, GenotypeQuality, PhredLikelihood, ReadDepth};
        Self {
            pl: pl
                .into_iter()
                .map(PhredLikelihood::from_i32_saturating)
                .collect(),
            gq: GenotypeQuality::from_i32_saturating(gq),
            ad: ad
                .into_iter()
                .map(AlleleDepth::from_i32_saturating)
                .collect(),
            dp: ReadDepth::from_i32_saturating(dp),
        }
    }

    /// PL vector as signed ints for VCF emit / Java PL round-trip.
    pub fn pl_as_i32(&self) -> Vec<i32> {
        self.pl.iter().map(|p| p.as_i32()).collect()
    }

    /// AD vector as signed ints for VCF emit.
    pub fn ad_as_i32(&self) -> Vec<i32> {
        self.ad.iter().map(|d| d.as_i32()).collect()
    }
}

/// Left-aligned REF/ALT alleles at a VCF position after normalization.
/// # Invariants
/// `position_1based` is the VCF POS after left-alignment (may shift from input).
/// `alternates` are non-empty for variant sites; REF is non-empty when valid.
/// # Ownership
/// Owns allele strings; produced by normalization helpers.
/// # Mutation
/// Immutable snapshot for emit/sort.
/// # Biological assumptions
/// Alleles are literal base strings on the reference contig (SNP/MNP/indel).
/// # Java equivalence
/// GATK / HTSJDK left-alignment and allele normalization (`VariantContext` allele lists).
#[derive(Debug, Clone, PartialEq)]
pub struct NormalizedVariantAlleles {
    pub position_1based: usize,
    pub reference: String,
    pub alternates: Vec<String>,
}

/// FORMAT phasing tags for one sample genotype (GT, PGT, PID, PS).
/// # Invariants
/// Phased het emits `|` GT with optional PS/PID/PGT; hom and no-call stay unphased.
/// # Ownership
/// Owns GT string and optional phase-set tags.
/// # Mutation
/// Immutable emit snapshot.
/// # Biological assumptions
/// Diploid sample; phase set applies only to heterozygous callable genotypes.
/// # Java equivalence
/// GATK VCF phasing emission (`GenotypeBuilder`, phase-set annotations).
#[derive(Debug, Clone, PartialEq)]
pub struct GenotypePhasingFields {
    pub gt: String,
    pub pgt: Option<String>,
    pub pid: Option<String>,
    pub ps: Option<i32>,
    pub phased: bool,
}

/// Site-level INFO counts: AC, AN, AF, NS, DP across called samples.
/// # Invariants
/// `ac.len` equals alt allele count; `af[i] = ac[i] / an` when `an > 0`.
/// No-call samples do not contribute to AN/NS.
/// # Ownership
/// Owns count vectors; returned from [`compute_core_variant_annotations`].
/// # Mutation
/// Immutable annotation bundle.
/// # Biological assumptions
/// Allele indices follow VCF REF=0, ALT=1.. ordering per sample genotypes.
/// # Java equivalence
/// GATK `ChromosomeCounts` / core INFO annotations (AC, AN, AF, NS, DP).
#[derive(Debug, Clone, PartialEq)]
pub struct CoreVariantAnnotations {
    pub ac: Vec<i32>,
    pub an: i32,
    pub af: Vec<f64>,
    pub ns: i32,
    pub dp: i32,
}

/// Per-sample genotype allele indices and optional depth for site INFO aggregation.
/// # Invariants
/// Allele indices are VCF-style (`-1` = no-call); must not exceed alt count at aggregation time.
/// # Ownership
/// Owns allele index vector; `dp` optional per sample.
/// # Mutation
/// Immutable input to [`compute_core_variant_annotations`].
/// # Biological assumptions
/// Diploid (or general ploidy) sample with REF/ALT indices already assigned.
/// # Java equivalence
/// GATK sample genotype / DP inputs to `VariantContext` INFO calculation.
#[derive(Debug, Clone, PartialEq)]
pub struct SampleAnnotationInput {
    pub genotype_alleles: Vec<i32>,
    pub dp: Option<i32>,
}

/// Minimal variant call record for deterministic sorting and emit scaffolding.
/// # Invariants
/// `position_1based` is 1-based VCF POS; `format_keys` lists emitted FORMAT tags.
/// # Ownership
/// Owns contig, alleles, and format key list.
/// # Mutation
/// Sort helpers reorder vectors in place; fields otherwise treated as immutable snapshots.
/// # Biological assumptions
/// Represents one variant row on a single contig with REF + ALT alleles.
/// # Java equivalence
/// Rust-native emit/sort carrier; mirrors `VariantContext` identity fields for parity dumps.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GenotypedVariantCall {
    pub contig: String,
    pub position_1based: usize,
    pub reference: String,
    pub alternates: Vec<String>,
    pub format_keys: Vec<String>,
}

/// Per-locus reference-confidence metrics for gVCF block construction.
/// # Invariants
/// `gq` and `dp` are non-negative Phred-style / count fields at one reference locus.
/// # Ownership
/// [`Copy`]-friendly scalar bundle in locus vectors.
/// # Mutation
/// Immutable per locus; block builders consume slices.
/// # Biological assumptions
/// Hom-ref confidence at a reference position (no alt allele called).
/// # Java equivalence
/// GATK reference-confidence / gVCF hom-ref site model (`ReferenceConfidenceVariantContext`).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ReferenceConfidenceLocus {
    pub position_1based: usize,
    pub gq: i32,
    pub dp: i32,
}

/// Merged gVCF reference block spanning contiguous hom-ref loci in one GQ band.
/// # Invariants
/// `start_1based <= end_1based`; block min/max DP and min RGQ summarize member loci.
/// # Ownership
/// Owns interval and band summary scalars.
/// # Mutation
/// Immutable output of gVCF block builders.
/// # Biological assumptions
/// Block covers contiguous reference positions with compatible GQ bands for compression.
/// # Java equivalence
/// GATK gVCF block combiner / `GvcfBlock` emission semantics.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GvcfBlock {
    pub start_1based: usize,
    pub end_1based: usize,
    pub gq_band_upper: i32,
    pub min_rgq: i32,
    pub min_dp: i32,
    pub max_dp: i32,
}

/// gVCF block fields as they appear on emitted records (POS, END INFO, DP bounds, GQ band).
/// # Invariants
/// `end_info >= start_1based` (VCF END ≥ POS); used for joint-genotype compatibility checks.
/// # Ownership
/// Scalar field bundle for validation and dumps.
/// # Mutation
/// Immutable record-shaped snapshot.
/// # Biological assumptions
/// Represents one compressed hom-ref block row for joint calling compatibility.
/// # Java equivalence
/// GATK gVCF `END` / GQ-band block records validated by joint genotyping tools.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GvcfBlockRecordFields {
    pub start_1based: usize,
    pub end_info: usize,
    pub min_dp: i32,
    pub max_dp: i32,
    pub gq_band_upper: i32,
    pub min_rgq: i32,
}

/// HC output mode controlling variant vs reference-block emission.
/// # Invariants
/// Each mode maps to a fixed [`LocusEmissionDecision`] for variant vs non-variant loci.
/// # Ownership
/// [`Copy`] enum; no allocation.
/// # Mutation
/// Immutable discriminant selected at engine configuration time.
/// # Biological assumptions
/// VCF emits variants only; GVCF emits blocks; BP-resolution emits per-base hom-ref sites.
/// # Java equivalence
/// GATK `-emit-ref-confidence` / HC emit mode (`VariantCallingMode`, gVCF vs VCF).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EmitMode {
    Vcf,
    Gvcf,
    BpResolution,
}

/// Per-locus emit decision derived from [`EmitMode`] and variant presence.
/// # Invariants
/// Mutually exclusive outcomes for a single locus in the emit planner.
/// # Ownership
/// [`Copy`] enum returned by [`decide_locus_emission`].
/// # Mutation
/// Immutable decision per locus.
/// # Biological assumptions
/// Distinguishes callable variant sites from hom-ref reference output shapes.
/// # Java equivalence
/// GATK HC locus emission branches in `HaplotypeCallerEngine` / gVCF writer paths.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LocusEmissionDecision {
    Skip,
    EmitVariantOnly,
    EmitReferenceBlock,
    EmitReferenceSite,
}

#[inline]
fn expected_diploid_genotype_count(allele_count: usize) -> usize {
    allele_count.saturating_mul(allele_count + 1) / 2
}

fn gq_band_upper(gq: i32, gq_bands: &[i32]) -> i32 {
    for b in gq_bands {
        if gq <= *b {
            return *b;
        }
    }
    gq_bands.last().copied().unwrap_or(100)
}

/// Java `GVCFBlockCombiner.parsePartitions` / `HomRefBlock.withinBounds` — GQ band `[lower, upper)`.
pub fn gvcf_gq_partition(gq: i32, gq_bands: &[i32]) -> (i32, i32) {
    let gq = gq.clamp(0, 99);
    let mut lower = 0;
    for &upper in gq_bands {
        if gq < upper {
            return (lower, upper);
        }
        lower = upper;
    }
    (lower, 100)
}

fn within_gvcf_partition(gq: i32, lower: i32, upper: i32) -> bool {
    let gq = gq.clamp(0, 99);
    gq >= lower && gq < upper
}

/// Merge semantics for reference confidence blocks:
/// loci must be adjacent,
/// GQ band compatible with block `min_rgq` (Java GVCFWriter fringe rule when enabled),
/// RGQ values must not diverge by more than `max_rgq_delta_within_block`.
/// # Invariants
/// `max_rgq_delta_within_block` caps RGQ divergence inside one block.
/// `java_gvcf_band_merge` enables Java fringe merge into min_rgq=0 blocks.
/// # Ownership
/// [`Copy`] policy for block builders.
/// # Mutation
/// Immutable per block-construction pass.
/// # Biological assumptions
/// Compresses contiguous hom-ref loci into gVCF blocks without changing callability.
/// # Java equivalence
/// GATK `GVCFWriter` block-merge / fringe rules (HC emit path).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct GvcfMergeSemantics {
    pub max_rgq_delta_within_block: i32,
    /// When true, merge rising-GQ fringe into blocks whose `min_rgq` is 0 (Java hom-ref span).
    pub java_gvcf_band_merge: bool,
}

impl Default for GvcfMergeSemantics {
    fn default() -> Self {
        Self {
            max_rgq_delta_within_block: 10,
            java_gvcf_band_merge: false,
        }
    }
}

/// Production HC gVCF block merge (Java `GVCFWriter` fringe + P8 base semantics).
pub fn gvcf_merge_semantics_hc_emit() -> GvcfMergeSemantics {
    GvcfMergeSemantics {
        max_rgq_delta_within_block: 10,
        java_gvcf_band_merge: true,
    }
}

fn trim_common_allele_context(position_1based: &mut usize, alleles: &mut [String]) {
    if alleles.is_empty() {
        return;
    }
    loop {
        if alleles.iter().any(|a| a.len() <= 1) {
            break;
        }
        let first = alleles[0].as_bytes()[0];
        if alleles.iter().all(|a| a.as_bytes()[0] == first) {
            for a in alleles.iter_mut() {
                a.remove(0);
            }
            *position_1based += 1;
        } else {
            break;
        }
    }
    loop {
        if alleles.iter().any(|a| a.len() <= 1) {
            break;
        }
        let first = alleles[0].as_bytes()[alleles[0].len() - 1];
        if alleles.iter().all(|a| a.as_bytes()[a.len() - 1] == first) {
            for a in alleles.iter_mut() {
                a.pop();
            }
        } else {
            break;
        }
    }
}

/// Normalize emitted alleles via minimal trimming and repeat-aware left alignment.
/// Inputs are expected in VCF representation (1-based `position_1based`,
/// non-empty REF and ALT alleles). `contig_bases` must contain the full contig sequence.
pub fn normalize_variant_alleles_left_aligned(
    position_1based: usize,
    reference: &str,
    alternates: &[String],
    contig_bases: &str,
) -> GatkResult<NormalizedVariantAlleles> {
    if position_1based == 0 {
        return Err(GatkError::argument(
            "variant position must be 1-based and >= 1",
        ));
    }
    if reference.is_empty() || alternates.is_empty() || alternates.iter().any(|a| a.is_empty()) {
        return Err(GatkError::argument(
            "variant normalization requires non-empty REF and ALT alleles",
        ));
    }
    let contig_len = contig_bases.len();
    if position_1based > contig_len {
        return Err(GatkError::argument(
            "variant normalization position is outside contig sequence",
        ));
    }

    let mut pos = position_1based;
    let mut alleles = std::iter::once(reference.to_string())
        .chain(alternates.iter().cloned())
        .collect::<Vec<_>>();

    trim_common_allele_context(&mut pos, &mut alleles);

    // Left-align in repetitive context for indel-like alleles by rotating one base left.
    let has_length_difference = alleles.iter().skip(1).any(|a| a.len() != alleles[0].len());
    if has_length_difference {
        while pos > 1 {
            let prev_base = contig_bases.as_bytes()[pos - 2] as char;
            let can_shift = alleles.iter().all(|a| {
                if let Some(last) = a.chars().last() {
                    last == prev_base
                } else {
                    false
                }
            });
            if !can_shift {
                break;
            }
            for a in alleles.iter_mut() {
                let mut rotated = String::with_capacity(a.len());
                rotated.push(prev_base);
                rotated.push_str(&a[..a.len() - 1]);
                *a = rotated;
            }
            pos -= 1;
            trim_common_allele_context(&mut pos, &mut alleles);
        }
    }

    Ok(NormalizedVariantAlleles {
        position_1based: pos,
        reference: alleles[0].clone(),
        alternates: alleles[1..].to_vec(),
    })
}

fn allele_to_gt_token(allele: i32) -> String {
    if allele < 0 {
        ".".to_string()
    } else {
        allele.to_string()
    }
}

/// Emit phasing-related genotype fields with conservative parity-minded semantics:
/// set exists.
/// Emit `PGT/PID/PS` only for phased calls.
pub fn emit_genotype_phasing_fields(
    alleles: &[i32],
    phasing_enabled: bool,
    phase_set: Option<i32>,
) -> GatkResult<GenotypePhasingFields> {
    if alleles.is_empty() {
        return Err(GatkError::argument(
            "phasing field emission requires at least one allele",
        ));
    }

    let is_diploid = alleles.len() == 2;
    let has_missing = alleles.iter().any(|a| *a < 0);
    let is_het = is_diploid && alleles[0] != alleles[1];
    let phased = phasing_enabled && phase_set.is_some() && is_het && !has_missing;
    let separator = if phased { "|" } else { "/" };
    let gt = alleles
        .iter()
        .map(|a| allele_to_gt_token(*a))
        .collect::<Vec<_>>()
        .join(separator);

    if phased {
        if let Some(ps) = phase_set {
            let pgt = gt.clone();
            let pid = format!("{ps}_{}_{}", alleles[0], alleles[1]);
            return Ok(GenotypePhasingFields {
                gt,
                pgt: Some(pgt),
                pid: Some(pid),
                ps: Some(ps),
                phased: true,
            });
        }
    }

    Ok(GenotypePhasingFields {
        gt,
        pgt: None,
        pid: None,
        ps: None,
        phased: false,
    })
}

/// Compute HC-core site annotations from sample-level genotype/depth inputs.
/// Semantics:
/// AC: alternate allele counts across called alleles
/// AN: total number of called alleles (non-missing, >= 0)
/// AF: AC / AN for each ALT (0.0 when AN == 0)
/// NS: number of samples with at least one called allele
/// DP: sum of non-negative sample depths provided
pub fn compute_core_variant_annotations(
    alt_allele_count: usize,
    samples: &[SampleAnnotationInput],
) -> GatkResult<CoreVariantAnnotations> {
    if alt_allele_count == 0 {
        return Err(GatkError::argument(
            "core annotations require at least one ALT allele",
        ));
    }

    let mut ac = vec![0_i32; alt_allele_count];
    let mut an = 0_i32;
    let mut ns = 0_i32;
    let mut dp = 0_i32;

    for sample in samples {
        if let Some(sample_dp) = sample.dp {
            if sample_dp < 0 {
                return Err(GatkError::argument(
                    "core annotations do not allow negative sample DP",
                ));
            }
            dp += sample_dp;
        }

        let mut sample_has_called = false;
        for allele in &sample.genotype_alleles {
            if *allele < 0 {
                continue;
            }
            sample_has_called = true;
            an += 1;
            if *allele > 0 {
                let alt_idx = (*allele as usize) - 1;
                if alt_idx >= alt_allele_count {
                    return Err(GatkError::argument(format!(
                        "core annotations allele index {} exceeds ALT allele count {}",
                        allele, alt_allele_count
                    )));
                }
                ac[alt_idx] += 1;
            }
        }
        if sample_has_called {
            ns += 1;
        }
    }

    let af = if an > 0 {
        ac.iter().map(|c| *c as f64 / an as f64).collect::<Vec<_>>()
    } else {
        vec![0.0; alt_allele_count]
    };

    Ok(CoreVariantAnnotations { ac, an, af, ns, dp })
}

/// Canonicalize FORMAT key ordering for deterministic VCF output.
/// Preferred prefix order: `GT`, `GQ`, `DP`, `AD`, `PL`.
/// Remaining keys keep deterministic lexicographic order.
pub fn canonicalize_format_keys(keys: &[String]) -> Vec<String> {
    // Lifetime: `&str` entries borrow from `keys` for this function only; own Strings
    // are allocated once when building the returned Vec.
    let mut seen: std::collections::BTreeSet<&str> = std::collections::BTreeSet::new();
    for k in keys {
        if !k.is_empty() {
            seen.insert(k.as_str());
        }
    }
    let preferred = ["GT", "GQ", "DP", "AD", "PL"];
    let mut out = Vec::new();
    for k in preferred {
        if seen.remove(k) {
            out.push(k.to_string());
        }
    }
    out.extend(seen.into_iter().map(str::to_string));
    out
}

/// Sort genotyped calls deterministically by dictionary order, then position, then allele key.
pub fn sort_genotyped_calls_deterministic(
    calls: &mut [GenotypedVariantCall],
    contig_order: &[String],
) -> GatkResult<()> {
    let rank = contig_order
        .iter()
        .enumerate()
        .map(|(idx, c)| (c.as_str(), idx))
        .collect::<std::collections::BTreeMap<_, _>>();
    for c in calls.iter() {
        if !rank.contains_key(c.contig.as_str()) {
            return Err(GatkError::argument(format!(
                "deterministic sort missing contig '{}' in contig order",
                c.contig
            )));
        }
    }
    calls.sort_by(|a, b| {
        let ra = rank[a.contig.as_str()];
        let rb = rank[b.contig.as_str()];
        ra.cmp(&rb)
            .then(a.position_1based.cmp(&b.position_1based))
            .then(a.reference.cmp(&b.reference))
            .then(a.alternates.join(",").cmp(&b.alternates.join(",")))
    });
    for c in calls {
        c.format_keys = canonicalize_format_keys(&c.format_keys);
    }
    Ok(())
}

/// Build GVCF reference-confidence blocks from per-locus confidence observations.
/// Step-101 contract:
/// adjacent loci in the same GQ band are merged into one block,
/// band changes or coordinate gaps force a new block,
/// each block tracks deterministic bounds (`start`, `end`, `min_dp`, `max_dp`, band upper bound).
pub fn build_gvcf_blocks(
    loci: &[ReferenceConfidenceLocus],
    gq_bands: &[i32],
) -> GatkResult<Vec<GvcfBlock>> {
    build_gvcf_blocks_with_semantics(loci, gq_bands, GvcfMergeSemantics::default())
}

/// Extended Step-102 builder that applies explicit reference-confidence merge semantics.
pub fn build_gvcf_blocks_with_semantics(
    loci: &[ReferenceConfidenceLocus],
    gq_bands: &[i32],
    semantics: GvcfMergeSemantics,
) -> GatkResult<Vec<GvcfBlock>> {
    if gq_bands.is_empty() {
        return Err(GatkError::argument(
            "GVCF block creation requires non-empty GQ bands",
        ));
    }
    if gq_bands.windows(2).any(|w| w[0] > w[1]) {
        return Err(GatkError::argument(
            "GVCF block creation requires non-decreasing GQ bands",
        ));
    }
    if loci.is_empty() {
        return Ok(Vec::new());
    }
    if semantics.max_rgq_delta_within_block < 0 {
        return Err(GatkError::argument(
            "GVCF merge semantics require non-negative max_rgq_delta_within_block",
        ));
    }
    for w in loci.windows(2) {
        if w[0].position_1based >= w[1].position_1based {
            return Err(GatkError::argument(
                "GVCF loci must be strictly increasing by position",
            ));
        }
    }
    if loci.iter().any(|l| l.position_1based == 0 || l.dp < 0) {
        return Err(GatkError::argument(
            "GVCF loci require 1-based positions and non-negative DP",
        ));
    }

    let mut out = Vec::<GvcfBlock>::new();
    let first = &loci[0];
    let mut cur = GvcfBlock {
        start_1based: first.position_1based,
        end_1based: first.position_1based,
        gq_band_upper: gq_band_upper(first.gq, gq_bands),
        min_rgq: first.gq,
        min_dp: first.dp,
        max_dp: first.dp,
    };

    for locus in loci.iter().skip(1) {
        let band = gq_band_upper(locus.gq, gq_bands);
        let contiguous = locus.position_1based == cur.end_1based + 1;
        let min_band = gq_band_upper(cur.min_rgq, gq_bands);
        let same_band = if semantics.java_gvcf_band_merge {
            band == min_band
                || (cur.min_rgq == 0 && locus.gq <= semantics.max_rgq_delta_within_block)
        } else {
            band == cur.gq_band_upper
        };
        let rgq_compatible = (locus.gq - cur.min_rgq).abs() <= semantics.max_rgq_delta_within_block;
        if contiguous && same_band && rgq_compatible {
            cur.end_1based = locus.position_1based;
            cur.min_rgq = cur.min_rgq.min(locus.gq);
            cur.gq_band_upper = gq_band_upper(cur.min_rgq, gq_bands);
            cur.min_dp = cur.min_dp.min(locus.dp);
            cur.max_dp = cur.max_dp.max(locus.dp);
        } else {
            out.push(cur);
            cur = GvcfBlock {
                start_1based: locus.position_1based,
                end_1based: locus.position_1based,
                gq_band_upper: band,
                min_rgq: locus.gq,
                min_dp: locus.dp,
                max_dp: locus.dp,
            };
        }
    }
    out.push(cur);
    Ok(out)
}

/// Java `GVCFBlockCombiner.addHomRefSite` — partition fixed at block start, merge while `withinBounds`.
pub fn build_gvcf_blocks_java_combiner(
    loci: &[ReferenceConfidenceLocus],
    gq_bands: &[i32],
) -> GatkResult<Vec<GvcfBlock>> {
    if gq_bands.is_empty() {
        return Err(GatkError::argument(
            "GVCF block creation requires non-empty GQ bands",
        ));
    }
    if loci.is_empty() {
        return Ok(Vec::new());
    }
    for w in loci.windows(2) {
        if w[0].position_1based >= w[1].position_1based {
            return Err(GatkError::argument(
                "GVCF loci must be strictly increasing by position",
            ));
        }
    }

    let mut out = Vec::<GvcfBlock>::new();
    let first = &loci[0];
    let (part_lo, part_hi) = gvcf_gq_partition(first.gq, gq_bands);
    let mut cur = GvcfBlock {
        start_1based: first.position_1based,
        end_1based: first.position_1based,
        gq_band_upper: part_hi,
        min_rgq: first.gq,
        min_dp: first.dp,
        max_dp: first.dp,
    };
    let mut cur_part_lo = part_lo;

    for locus in loci.iter().skip(1) {
        let contiguous = locus.position_1based == cur.end_1based + 1;
        let same_partition = within_gvcf_partition(locus.gq, cur_part_lo, cur.gq_band_upper);
        if contiguous && same_partition {
            cur.end_1based = locus.position_1based;
            cur.min_rgq = cur.min_rgq.min(locus.gq);
            cur.min_dp = cur.min_dp.min(locus.dp);
            cur.max_dp = cur.max_dp.max(locus.dp);
        } else {
            out.push(cur);
            let (part_lo, part_hi) = gvcf_gq_partition(locus.gq, gq_bands);
            cur_part_lo = part_lo;
            cur = GvcfBlock {
                start_1based: locus.position_1based,
                end_1based: locus.position_1based,
                gq_band_upper: part_hi,
                min_rgq: locus.gq,
                min_dp: locus.dp,
                max_dp: locus.dp,
            };
        }
    }
    out.push(cur);
    Ok(out)
}

/// Production HC gVCF block builder (Java `GVCFBlockCombiner` partition merge).
pub fn build_gvcf_blocks_hc_emit(
    loci: &[ReferenceConfidenceLocus],
    gq_bands: &[i32],
) -> GatkResult<Vec<GvcfBlock>> {
    build_gvcf_blocks_java_combiner(loci, gq_bands)
}

/// Second-pass merge for emitted gVCF blocks (Java `GVCFWriter` / `GvcfBlockCombiner` fringe).
/// Merges adjacent blocks when combined `min_rgq` stays in the same GQ band as the left block
/// and per-block min RGQ values are within `max_rgq_delta` (hom-ref span absorption).
pub fn coalesce_gvcf_blocks_for_emit(
    blocks: Vec<GvcfBlock>,
    gq_bands: &[i32],
    max_rgq_delta: i32,
) -> Vec<GvcfBlock> {
    if blocks.is_empty() {
        return blocks;
    }
    // Lifetime: take ownership of the first block by move; no clone of `GvcfBlock`.
    let mut blocks = blocks.into_iter();
    let Some(first) = blocks.next() else {
        return vec![];
    };
    let mut out = vec![first];
    for block in blocks {
        let Some(cur) = out.last_mut() else {
            continue;
        };
        let contiguous = block.start_1based == cur.end_1based + 1;
        let combined_min = cur.min_rgq.min(block.min_rgq);
        let rgq_ok = (cur.min_rgq - block.min_rgq).abs() <= max_rgq_delta;
        let band_ok = gq_band_upper(combined_min, gq_bands) == gq_band_upper(cur.min_rgq, gq_bands);
        // Only absorb low-GQ fringe into a leading hom-ref (min_rgq==0) span — avoids over-merge in active tails.
        let zero_hom_ref_fringe = cur.min_rgq == 0 && block.min_rgq <= max_rgq_delta;
        if contiguous && zero_hom_ref_fringe && rgq_ok && band_ok {
            cur.end_1based = block.end_1based;
            cur.min_rgq = combined_min;
            cur.gq_band_upper = gq_band_upper(cur.min_rgq, gq_bands);
            cur.min_dp = cur.min_dp.min(block.min_dp);
            cur.max_dp = cur.max_dp.max(block.max_dp);
        } else {
            out.push(block);
        }
    }
    out
}

/// Convert a computed GVCF block into record-facing fields.
/// Step-103 contract:
/// `END` equals block end position (inclusive),
/// single-locus blocks have `END == POS`,
/// emitted bounds/metrics are copied deterministically from block state.
pub fn gvcf_block_to_record_fields(block: &GvcfBlock) -> GatkResult<GvcfBlockRecordFields> {
    if block.start_1based == 0 || block.end_1based == 0 || block.start_1based > block.end_1based {
        return Err(GatkError::argument(
            "gVCF block record emission requires valid 1-based start/end bounds",
        ));
    }
    Ok(GvcfBlockRecordFields {
        start_1based: block.start_1based,
        end_info: block.end_1based,
        min_dp: block.min_dp,
        max_dp: block.max_dp,
        gq_band_upper: block.gq_band_upper,
        min_rgq: block.min_rgq,
    })
}

/// Decide per-locus emission behavior by mode and variant/reference status.
/// Step-105 contract:
/// VCF mode emits only variant loci.
/// GVCF mode emits variants + reference blocks.
/// BP_RESOLUTION-like mode emits variants + per-base reference sites.
pub fn decide_locus_emission(mode: EmitMode, has_variant: bool) -> LocusEmissionDecision {
    match (mode, has_variant) {
        (EmitMode::Vcf, true) => LocusEmissionDecision::EmitVariantOnly,
        (EmitMode::Vcf, false) => LocusEmissionDecision::Skip,
        (EmitMode::Gvcf, true) => LocusEmissionDecision::EmitVariantOnly,
        (EmitMode::Gvcf, false) => LocusEmissionDecision::EmitReferenceBlock,
        (EmitMode::BpResolution, true) => LocusEmissionDecision::EmitVariantOnly,
        (EmitMode::BpResolution, false) => LocusEmissionDecision::EmitReferenceSite,
    }
}

/// Summarize region-level no-variation emission behavior for a given mode.
/// # Invariants
/// Counters sum to modes of [`decide_locus_emission`] over `loci_total` non-variant loci.
/// # Ownership
/// Owned counter bundle from [`summarize_no_variation_region`].
/// # Mutation
/// Immutable summary.
/// # Biological assumptions
/// Region has no variant loci; emission is VCF-skip / gVCF-block / BP-resolution sites.
/// # Java equivalence
/// Documents GATK HC no-variation region emit behavior by mode.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct NoVariationRegionSummary {
    pub loci_total: usize,
    pub reference_blocks_emitted: usize,
    pub reference_sites_emitted: usize,
    pub variant_sites_emitted: usize,
}

/// Result of validating gVCF block records for joint-genotyping compatibility.
/// # Invariants
/// When `compatible` is true, records are sorted, non-overlapping, and satisfy END ≥ POS.
/// # Ownership
/// Lightweight summary returned from [`validate_joint_compatibility_gvcf_records`].
/// # Mutation
/// Immutable validation outcome.
/// # Biological assumptions
/// Input records represent one sample's gVCF blocks on one contig.
/// # Java equivalence
/// GATK joint-genotyping gVCF compatibility checks (GenomicsDB / CombineGVCFs expectations).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct JointCompatibilitySummary {
    pub records_total: usize,
    pub compatible: bool,
}

/// Evaluate expected behavior on a region with no variant loci.
/// Step-106 contract:
/// VCF mode emits nothing for purely non-variant regions.
/// GVCF mode emits block-oriented reference output.
/// BP-resolution mode emits one reference site per locus.
pub fn summarize_no_variation_region(
    mode: EmitMode,
    locus_count: usize,
) -> NoVariationRegionSummary {
    let mut summary = NoVariationRegionSummary {
        loci_total: locus_count,
        reference_blocks_emitted: 0,
        reference_sites_emitted: 0,
        variant_sites_emitted: 0,
    };
    for _ in 0..locus_count {
        match decide_locus_emission(mode, false) {
            LocusEmissionDecision::Skip => {}
            LocusEmissionDecision::EmitVariantOnly => summary.variant_sites_emitted += 1,
            LocusEmissionDecision::EmitReferenceBlock => summary.reference_blocks_emitted += 1,
            LocusEmissionDecision::EmitReferenceSite => summary.reference_sites_emitted += 1,
        }
    }
    summary
}

/// Validate a sequence of gVCF-like records for downstream joint-genotyping compatibility.
/// Step-107 sanity checks:
/// each record must satisfy `END >= POS`,
/// records must be sorted by start position,
/// records must not overlap in genomic span.
pub fn validate_joint_compatibility_gvcf_records(
    records: &[GvcfBlockRecordFields],
) -> GatkResult<JointCompatibilitySummary> {
    if records.is_empty() {
        return Ok(JointCompatibilitySummary {
            records_total: 0,
            compatible: true,
        });
    }

    let mut prev_start = 0usize;
    let mut prev_end = 0usize;
    for (idx, rec) in records.iter().enumerate() {
        if rec.start_1based == 0 {
            return Err(GatkError::argument(format!(
                "joint-compat record[{idx}] has invalid 1-based start=0"
            )));
        }
        if rec.end_info < rec.start_1based {
            return Err(GatkError::argument(format!(
                "joint-compat record[{idx}] has END < POS ({} < {})",
                rec.end_info, rec.start_1based
            )));
        }
        if idx > 0 {
            if rec.start_1based < prev_start {
                return Err(GatkError::argument(format!(
                    "joint-compat records are not sorted by POS at index {idx}"
                )));
            }
            if rec.start_1based <= prev_end {
                return Err(GatkError::argument(format!(
                    "joint-compat records overlap at index {idx}: start={} previous_end={}",
                    rec.start_1based, prev_end
                )));
            }
        }
        prev_start = rec.start_1based;
        prev_end = rec.end_info;
    }

    Ok(JointCompatibilitySummary {
        records_total: records.len(),
        compatible: true,
    })
}

/// Combine genotype likelihoods and priors in log10 space and normalize to posterior probabilities.
pub fn genotype_posteriors_from_log10_likelihoods(
    genotype_log10_likelihoods: &[f64],
    genotype_log10_priors: &[f64],
) -> GatkResult<GenotypePosterior> {
    if genotype_log10_likelihoods.is_empty() {
        return Err(GatkError::argument(
            "genotype posterior requires at least one genotype likelihood",
        ));
    }
    if genotype_log10_likelihoods.len() != genotype_log10_priors.len() {
        return Err(GatkError::argument(
            "genotype posterior requires matching likelihood and prior vector lengths",
        ));
    }

    let mut log10_post = Vec::with_capacity(genotype_log10_likelihoods.len());
    for (idx, (ll, prior)) in genotype_log10_likelihoods
        .iter()
        .copied()
        .zip(genotype_log10_priors.iter().copied())
        .enumerate()
    {
        if !ll.is_finite() || !prior.is_finite() {
            return Err(GatkError::argument(format!(
                "genotype posterior requires finite likelihood/prior values (idx={idx})"
            )));
        }
        log10_post.push(ll + prior);
    }

    // Lifetime: keep owned `log10_post` for the return value; normalize borrows it.
    let post = normalize_from_log10_to_linear_space(&log10_post);
    let most_likely_genotype_index = post
        .iter()
        .enumerate()
        .max_by(|a, b| a.1.total_cmp(b.1))
        .map(|(i, _)| i)
        .unwrap_or(0);
    Ok(GenotypePosterior {
        genotype_log10_posteriors: log10_post,
        genotype_posteriors: post,
        most_likely_genotype_index,
    })
}

/// P7 / HTSJDK GQ: hom-ref-best uses PL gap; het/alt-best uses log10 error from max-LL genotype.
pub fn gq_phred_p7_compatible(
    genotype_log10_likelihoods: &[f64],
    pls: &[i32],
    allele_depths: &[i32],
) -> i32 {
    let raw = if genotype_log10_likelihoods.len() == 3 {
        let best = best_biallelic_diploid_genotype_index(genotype_log10_likelihoods, allele_depths);
        if best == 1 {
            // HTSJDK/GATK FORMAT GQ is always capped at 99.
            if pls.len() >= 3 && pls[1] == 0 && pls[0] >= 39 && pls[2] >= 39 {
                pls[0].min(pls[2])
            } else if pls.len() >= 3 && pls[0] == pls[2] && pls[0] >= 39 {
                gq_phred_from_pl(pls)
            } else {
                gq_phred_from_biallelic_log10(genotype_log10_likelihoods, best)
            }
        } else if best == 0 || best == 2 {
            gq_phred_from_pl(pls)
        } else {
            gq_phred_from_biallelic_log10(genotype_log10_likelihoods, best)
        }
    } else {
        gq_phred_from_pl(pls)
    };
    raw.clamp(0, 99)
}

/// FORMAT GQ from normalized PLs: gap between best and second-best genotype (HTSJDK / P7 contract).
pub fn gq_phred_from_pl(pl: &[i32]) -> i32 {
    if pl.is_empty() {
        return 0;
    }
    let min_pl = pl.iter().copied().min().unwrap_or(0);
    if pl.iter().filter(|p| **p == min_pl).count() > 1 {
        return 0;
    }
    let second = pl
        .iter()
        .copied()
        .filter(|p| *p > min_pl)
        .min()
        .unwrap_or(min_pl);
    (second - min_pl).min(99)
}

/// GATK `GenotypeLikelihoods.getGQLog10FromLikelihoods` → phred GQ for biallelic diploid sites.
pub fn gq_phred_from_biallelic_log10(genotype_log10_likelihoods: &[f64], best_idx: usize) -> i32 {
    if genotype_log10_likelihoods.len() != 3 {
        return 99;
    }
    let log10_p_error = match best_idx {
        0 => crate::activity_scoring::log10_sum_log10(&[
            genotype_log10_likelihoods[1],
            genotype_log10_likelihoods[2],
        ]),
        1 => crate::activity_scoring::log10_sum_log10(&[
            genotype_log10_likelihoods[0],
            genotype_log10_likelihoods[2],
        ]),
        _ => crate::activity_scoring::log10_sum_log10(&[
            genotype_log10_likelihoods[0],
            genotype_log10_likelihoods[1],
        ]),
    }
    .min(0.0);
    (-10.0 * log10_p_error).round().clamp(0.0, 99.0) as i32
}

/// GATK `GenotypeAssignmentMethod.USE_PLS_TO_ASSIGN`: genotype with minimum PL.
pub fn biallelic_genotype_index_from_pl(
    pl: &[crate::bio_ids::PhredLikelihood],
) -> crate::bio_ids::DiploidGenotypeIndex {
    let idx = pl
        .iter()
        .enumerate()
        .min_by_key(|(_, p)| p.get())
        .map(|(i, _)| i)
        .unwrap_or(0);
    crate::bio_ids::DiploidGenotypeIndex::try_new(idx as u8)
        .unwrap_or(crate::bio_ids::DiploidGenotypeIndex::HOM_REF)
}

/// Minimum-PL genotype index for any allele count (not capped at biallelic `{0,1,2}`).
///
/// [`biallelic_genotype_index_from_pl`] uses [`crate::bio_ids::DiploidGenotypeIndex`], which
/// rejects indices `> 2` and falls back to hom-ref. Merged SNP+indel sites have 6 PLs; the
/// best state can be `1/2` (index 4).
pub fn best_pl_index(pl: &[crate::bio_ids::PhredLikelihood]) -> usize {
    pl.iter()
        .enumerate()
        .min_by_key(|(_, p)| p.get())
        .map(|(i, _)| i)
        .unwrap_or(0)
}

/// Diploid GT allele indices for PL index `best` with `n_alleles` (VCF PL ordering).
pub fn diploid_genotype_alleles_from_pl_index(n_alleles: usize, best: usize) -> Vec<i32> {
    let mut k = 0usize;
    for j in 0..n_alleles {
        for i in 0..=j {
            if k == best {
                return vec![i as i32, j as i32];
            }
            k += 1;
        }
    }
    vec![0, 0]
}

/// Wire-integer overload for tests / dump helpers that still hold raw PL vectors.
pub fn biallelic_genotype_index_from_pl_i32(pl: &[i32]) -> crate::bio_ids::DiploidGenotypeIndex {
    let owned: Vec<_> = pl
        .iter()
        .copied()
        .map(crate::bio_ids::PhredLikelihood::from_i32_saturating)
        .collect();
    biallelic_genotype_index_from_pl(&owned)
}

/// Best diploid genotype index; near-tie het/hom-alt resolves toward hom-alt when AD supports alt.
pub fn best_biallelic_diploid_genotype_index(
    genotype_log10_likelihoods: &[f64],
    allele_depths: &[i32],
) -> usize {
    let mut best = 0usize;
    let mut best_ll = f64::NEG_INFINITY;
    for (i, &gl) in genotype_log10_likelihoods.iter().enumerate() {
        if gl > best_ll {
            best_ll = gl;
            best = i;
        }
    }
    if genotype_log10_likelihoods.len() >= 3 && allele_depths.len() >= 2 {
        let ref_d = allele_depths[0].max(0);
        let alt_d = allele_depths[1].max(0);
        let het_gl = genotype_log10_likelihoods[1];
        let hom_alt_gl = genotype_log10_likelihoods[2];
        if alt_d > ref_d + 1 && (het_gl - hom_alt_gl).abs() < 1e-3 {
            return 2;
        }
    }
    best
}

/// Emit PL/GQ/AD/DP semantics from diploid genotype likelihoods and allele depths.
/// `genotype_log10_likelihoods`: log10 likelihoods in GT ordering (e.g. [0/0,0/1,1/1])
/// `allele_depths`: per-allele depths in VCF allele order ([REF, ALT1,...])
pub fn emit_genotype_format_fields(
    genotype_log10_likelihoods: &[f64],
    allele_depths: &[i32],
) -> GatkResult<GenotypeFormatFields> {
    if genotype_log10_likelihoods.is_empty() {
        return Err(GatkError::argument(
            "PL/GQ emission requires at least one genotype likelihood",
        ));
    }
    if allele_depths.is_empty() {
        return Err(GatkError::argument(
            "AD/DP emission requires at least one allele depth",
        ));
    }
    if genotype_log10_likelihoods.iter().any(|v| !v.is_finite()) {
        return Err(GatkError::argument(
            "PL/GQ emission requires finite genotype likelihoods",
        ));
    }
    if allele_depths.iter().any(|d| *d < 0) {
        return Err(GatkError::argument(
            "AD/DP emission does not allow negative allele depths",
        ));
    }
    let expected_gl_count = expected_diploid_genotype_count(allele_depths.len());
    if genotype_log10_likelihoods.len() != expected_gl_count {
        return Err(GatkError::argument(format!(
            "PL/GQ emission genotype likelihood count mismatch: got {} expected {} for {} alleles",
            genotype_log10_likelihoods.len(),
            expected_gl_count,
            allele_depths.len()
        )));
    }

    // PL: phred-scaled normalized likelihoods, min value anchored to 0.
    let max_ll = genotype_log10_likelihoods
        .iter()
        .copied()
        .fold(f64::NEG_INFINITY, f64::max);
    let mut pls = genotype_log10_likelihoods
        .iter()
        .map(|ll| ((-10.0 * (ll - max_ll)).round() as i32).max(0))
        .collect::<Vec<_>>();
    let min_pl = pls.iter().copied().min().unwrap_or(0);
    for pl in &mut pls {
        *pl -= min_pl;
    }

    let gq = gq_phred_p7_compatible(genotype_log10_likelihoods, &pls, allele_depths);

    let dp = allele_depths.iter().sum::<i32>();
    Ok(GenotypeFormatFields::from_wire(
        pls,
        gq,
        allele_depths.to_vec(),
        dp,
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::gvcf_writer::GATK_HC_DEFAULT_GQB;

    #[test]
    fn aggregation_sums_log10_likelihoods_per_haplotype() {
        let rows = vec![
            ReadLikelihoodRow {
                read_index: 0,
                read_id: "r1".to_string(),
                haplotype_log10_likelihoods: vec![-0.1, -1.2, -2.3],
            },
            ReadLikelihoodRow {
                read_index: 0,
                read_id: "r2".to_string(),
                haplotype_log10_likelihoods: vec![-0.3, -0.8, -2.6],
            },
            ReadLikelihoodRow {
                read_index: 0,
                read_id: "r3".to_string(),
                haplotype_log10_likelihoods: vec![-0.2, -1.1, -2.0],
            },
        ];
        let agg = aggregate_haplotype_log10_likelihoods(&rows).expect("aggregation");
        assert_eq!(agg.read_count, 3);
        let expected = [-0.6_f64, -3.1_f64, -6.9_f64];
        for (got, exp) in agg.haplotype_log10_sums.iter().zip(expected.iter()) {
            assert!((got - exp).abs() <= 1e-12, "got={got} expected={exp}");
        }
        assert_eq!(
            best_haplotype_index(&agg),
            Some(crate::bio_ids::HaplotypeIndex::new(0))
        );
    }

    #[test]
    fn aggregation_rejects_haplotype_width_mismatch() {
        let rows = vec![
            ReadLikelihoodRow {
                read_index: 0,
                read_id: "r1".to_string(),
                haplotype_log10_likelihoods: vec![-0.1, -1.2],
            },
            ReadLikelihoodRow {
                read_index: 0,
                read_id: "r2".to_string(),
                haplotype_log10_likelihoods: vec![-0.2],
            },
        ];
        let err = aggregate_haplotype_log10_likelihoods(&rows).expect_err("width mismatch");
        assert!(err.to_string().contains("haplotype width mismatch"));
    }

    #[test]
    fn aggregation_skips_non_finite_likelihoods() {
        let rows = vec![ReadLikelihoodRow {
            read_index: 0,
            read_id: "r1".to_string(),
            haplotype_log10_likelihoods: vec![f64::NEG_INFINITY],
        }];
        let agg = aggregate_haplotype_log10_likelihoods(&rows).expect("skip -inf");
        assert_eq!(agg.haplotype_log10_sums, vec![0.0]);
        assert_eq!(agg.read_count, 1);
    }

    #[test]
    fn biallelic_priors_are_valid_and_finite() {
        let priors = biallelic_diploid_log10_priors(BiallelicDiploidPriorModel::default())
            .expect("default priors");
        assert_eq!(priors.len(), 3);
        assert!(priors.iter().all(|v| v.is_finite()));
    }

    #[test]
    fn genotype_posteriors_combine_likelihoods_and_priors() {
        let priors = biallelic_diploid_log10_priors(BiallelicDiploidPriorModel {
            het_prior: 1e-2,
            hom_var_prior: 1e-3,
        })
        .expect("priors");
        let likelihoods = vec![-0.2, -1.5, -2.2];
        let post =
            genotype_posteriors_from_log10_likelihoods(&likelihoods, &priors).expect("posterior");
        let sum: f64 = post.genotype_posteriors.iter().sum();
        assert!((sum - 1.0).abs() <= 1e-12, "sum={sum}");
        assert_eq!(post.most_likely_genotype_index, 0);
    }

    #[test]
    fn format_field_emission_pl_gq_ad_dp_basic_semantics() {
        let fields = emit_genotype_format_fields(&[-0.2, -1.5, -2.2], &[12, 5]).expect("fields");
        assert_eq!(fields.ad_as_i32(), vec![12, 5]);
        assert_eq!(fields.dp.get(), 17);
        assert_eq!(fields.pl_as_i32(), vec![0, 13, 20]);
        assert_eq!(fields.gq.get(), 13, "gq={}", fields.gq.as_i32());
    }

    #[test]
    fn format_field_emission_synthetic_het_gq_matches_p7() {
        let fields =
            emit_genotype_format_fields(&[-2.0, -0.1, -2.0], &[5, 12]).expect("synthetic het");
        assert_eq!(fields.pl_as_i32(), vec![19, 0, 19]);
        assert_eq!(fields.gq.get(), 17, "gq={}", fields.gq.as_i32());
    }

    #[test]
    fn het_strong_pl_gap_gq_is_capped_at_99() {
        // PLS gap path previously returned uncapped min(PL0,PL2); HTSJDK caps at 99.
        let gq = gq_phred_p7_compatible(&[-5.0, 0.0, -5.0], &[500, 0, 480], &[20, 25]);
        assert_eq!(gq, 99, "gq={gq}");
    }

    #[test]
    fn format_field_emission_rejects_negative_depths() {
        let err = emit_genotype_format_fields(&[-0.2, -1.5, -2.2], &[12, -1]).expect_err("depth");
        assert!(err.to_string().contains("negative allele depths"));
    }

    #[test]
    fn biallelic_genotype_fields_cardinality_and_shape() {
        let fields = emit_genotype_format_fields(&[-0.2, -1.5, -2.2], &[10, 4]).expect("biallelic");
        assert_eq!(fields.pl.len(), 3);
        assert_eq!(fields.ad.len(), 2);
        assert_eq!(fields.dp.get(), 14);
    }

    #[test]
    fn multiallelic_genotype_fields_cardinality_and_shape() {
        // Diploid tri-allelic site has 6 genotype likelihood entries.
        let gl = [-0.3, -1.1, -2.0, -1.6, -2.4, -3.2];
        let fields = emit_genotype_format_fields(&gl, &[18, 6, 3]).expect("multiallelic");
        assert_eq!(fields.pl.len(), 6);
        assert_eq!(fields.ad.len(), 3);
        assert_eq!(fields.dp.get(), 27);
        assert_eq!(fields.pl.get(0).map(|p| p.as_i32()).unwrap_or(0), 0);
        assert!(fields.gq.as_i32() >= 0 && fields.gq.get() <= 99);
    }

    #[test]
    fn multiallelic_rejects_gl_cardinality_mismatch() {
        let err =
            emit_genotype_format_fields(&[-0.3, -1.1, -2.0], &[18, 6, 3]).expect_err("mismatch");
        assert!(err
            .to_string()
            .contains("genotype likelihood count mismatch"));
    }

    #[test]
    fn normalize_variant_trims_shared_prefix_suffix() {
        let out = normalize_variant_alleles_left_aligned(
            100,
            "AAT",
            &[String::from("AGT")],
            &"N".repeat(500),
        )
        .expect("normalize");
        assert_eq!(out.position_1based, 101);
        assert_eq!(out.reference, "A");
        assert_eq!(out.alternates, vec![String::from("G")]);
    }

    #[test]
    fn normalize_variant_left_aligns_deletion_in_homopolymer() {
        // Deleting one A in AAAAAA should left-align to the earliest position.
        let out = normalize_variant_alleles_left_aligned(4, "AA", &[String::from("A")], "AAAAAA")
            .expect("normalize");
        assert_eq!(out.position_1based, 1);
        assert_eq!(out.reference, "AA");
        assert_eq!(out.alternates, vec![String::from("A")]);
    }

    #[test]
    fn phasing_fields_emit_for_diploid_het_with_phase_set() {
        let fields = emit_genotype_phasing_fields(&[0, 1], true, Some(1234)).expect("phasing");
        assert_eq!(fields.gt, "0|1");
        assert_eq!(fields.pgt, Some("0|1".to_string()));
        assert_eq!(fields.ps, Some(1234));
        assert!(fields.pid.is_some());
        assert!(fields.phased);
    }

    #[test]
    fn phasing_fields_remain_unphased_without_phase_set_or_for_hom() {
        let no_ps = emit_genotype_phasing_fields(&[0, 1], true, None).expect("no ps");
        assert_eq!(no_ps.gt, "0/1");
        assert!(!no_ps.phased);
        assert!(no_ps.pgt.is_none());

        let hom = emit_genotype_phasing_fields(&[1, 1], true, Some(1234)).expect("hom");
        assert_eq!(hom.gt, "1/1");
        assert!(!hom.phased);
        assert!(hom.pid.is_none());
    }

    #[test]
    fn core_annotations_biallelic_ac_an_af_ns_dp() {
        let ann = compute_core_variant_annotations(
            1,
            &[
                SampleAnnotationInput {
                    genotype_alleles: vec![0, 1],
                    dp: Some(12),
                },
                SampleAnnotationInput {
                    genotype_alleles: vec![1, 1],
                    dp: Some(20),
                },
                SampleAnnotationInput {
                    genotype_alleles: vec![-1, -1],
                    dp: Some(0),
                },
            ],
        )
        .expect("annotations");
        assert_eq!(ann.ac, vec![3]);
        assert_eq!(ann.an, 4);
        assert_eq!(ann.ns, 2);
        assert_eq!(ann.dp, 32);
        assert!((ann.af[0] - 0.75).abs() <= 1e-12);
    }

    #[test]
    fn core_annotations_multiallelic_counts_match_indices() {
        let ann = compute_core_variant_annotations(
            2,
            &[
                SampleAnnotationInput {
                    genotype_alleles: vec![1, 2],
                    dp: Some(18),
                },
                SampleAnnotationInput {
                    genotype_alleles: vec![0, 2],
                    dp: Some(10),
                },
            ],
        )
        .expect("annotations");
        assert_eq!(ann.ac, vec![1, 2]);
        assert_eq!(ann.an, 4);
        assert_eq!(ann.ns, 2);
        assert_eq!(ann.dp, 28);
        assert!((ann.af[0] - 0.25).abs() <= 1e-12);
        assert!((ann.af[1] - 0.5).abs() <= 1e-12);
    }

    #[test]
    fn core_annotations_reject_invalid_alt_index() {
        let err = compute_core_variant_annotations(
            1,
            &[SampleAnnotationInput {
                genotype_alleles: vec![0, 2],
                dp: Some(8),
            }],
        )
        .expect_err("invalid alt index");
        assert!(err.to_string().contains("exceeds ALT allele count"));
    }

    #[test]
    fn no_call_phasing_fields_emit_dot_gt_without_phase_tags() {
        let fields = emit_genotype_phasing_fields(&[-1, -1], true, Some(4321)).expect("no-call");
        assert_eq!(fields.gt, "./.");
        assert!(!fields.phased);
        assert!(fields.pgt.is_none());
        assert!(fields.pid.is_none());
        assert!(fields.ps.is_none());
    }

    #[test]
    fn no_call_annotations_have_zero_an_af_and_ns() {
        let ann = compute_core_variant_annotations(
            1,
            &[
                SampleAnnotationInput {
                    genotype_alleles: vec![-1, -1],
                    dp: Some(0),
                },
                SampleAnnotationInput {
                    genotype_alleles: vec![-1, -1],
                    dp: None,
                },
            ],
        )
        .expect("no-call annotations");
        assert_eq!(ann.ac, vec![0]);
        assert_eq!(ann.an, 0);
        assert_eq!(ann.ns, 0);
        assert_eq!(ann.dp, 0);
        assert_eq!(ann.af, vec![0.0]);
    }

    #[test]
    fn low_confidence_locus_has_low_gq_when_likelihoods_nearly_equal() {
        let fields =
            emit_genotype_format_fields(&[-1.000, -1.001, -1.002], &[9, 9]).expect("low-conf");
        assert_eq!(fields.pl_as_i32(), vec![0, 0, 0]);
        assert!(
            fields.gq.as_i32() < 10,
            "near-tie GLs should yield low GQ, got {}",
            fields.gq.as_i32()
        );
        assert_eq!(fields.dp.get(), 18);
    }

    #[test]
    fn posterior_selects_expected_genotype_for_low_confidence_case() {
        let priors = biallelic_diploid_log10_priors(BiallelicDiploidPriorModel {
            het_prior: 0.2,
            hom_var_prior: 0.1,
        })
        .expect("priors");
        let post = genotype_posteriors_from_log10_likelihoods(&[-1.000, -1.001, -1.002], &priors)
            .expect("posterior");
        let sum: f64 = post.genotype_posteriors.iter().sum();
        assert!((sum - 1.0).abs() <= 1e-12, "sum={sum}");
        // Priors slightly favor hom-ref in this setup.
        assert_eq!(post.most_likely_genotype_index, 0);
    }

    #[test]
    fn multi_sample_annotation_is_invariant_to_sample_order() {
        let samples_a = vec![
            SampleAnnotationInput {
                genotype_alleles: vec![0, 1],
                dp: Some(11),
            },
            SampleAnnotationInput {
                genotype_alleles: vec![1, 2],
                dp: Some(17),
            },
            SampleAnnotationInput {
                genotype_alleles: vec![0, 2],
                dp: Some(9),
            },
            SampleAnnotationInput {
                genotype_alleles: vec![-1, -1],
                dp: Some(0),
            },
        ];
        let mut samples_b = samples_a.clone();
        samples_b.reverse();

        let a = compute_core_variant_annotations(2, &samples_a).expect("annotations A");
        let b = compute_core_variant_annotations(2, &samples_b).expect("annotations B");
        assert_eq!(a, b);
    }

    #[test]
    fn multi_sample_regression_mixed_called_and_nocall_counts_stable() {
        let ann = compute_core_variant_annotations(
            2,
            &[
                SampleAnnotationInput {
                    genotype_alleles: vec![0, 0],
                    dp: Some(14),
                },
                SampleAnnotationInput {
                    genotype_alleles: vec![0, 1],
                    dp: Some(12),
                },
                SampleAnnotationInput {
                    genotype_alleles: vec![1, 2],
                    dp: Some(19),
                },
                SampleAnnotationInput {
                    genotype_alleles: vec![-1, -1],
                    dp: None,
                },
            ],
        )
        .expect("mixed multi-sample");
        // Called alleles: (0,0),(0,1),(1,2) => AN=6, AC=[2,1]
        assert_eq!(ann.an, 6);
        assert_eq!(ann.ac, vec![2, 1]);
        assert_eq!(ann.ns, 3);
        assert_eq!(ann.dp, 45);
        assert!((ann.af[0] - (2.0 / 6.0)).abs() <= 1e-12);
        assert!((ann.af[1] - (1.0 / 6.0)).abs() <= 1e-12);
    }

    #[test]
    fn deterministic_call_sort_orders_by_contig_then_pos_then_alleles() {
        let mut calls = vec![
            GenotypedVariantCall {
                contig: "chr2".to_string(),
                position_1based: 10,
                reference: "A".to_string(),
                alternates: vec!["G".to_string()],
                format_keys: vec!["AD".to_string(), "GT".to_string(), "PL".to_string()],
            },
            GenotypedVariantCall {
                contig: "chr1".to_string(),
                position_1based: 20,
                reference: "C".to_string(),
                alternates: vec!["T".to_string()],
                format_keys: vec!["PL".to_string(), "GT".to_string(), "DP".to_string()],
            },
            GenotypedVariantCall {
                contig: "chr1".to_string(),
                position_1based: 10,
                reference: "A".to_string(),
                alternates: vec!["C".to_string()],
                format_keys: vec!["GQ".to_string(), "GT".to_string(), "AD".to_string()],
            },
        ];
        sort_genotyped_calls_deterministic(&mut calls, &["chr1".to_string(), "chr2".to_string()])
            .expect("sort");
        assert_eq!(calls[0].contig, "chr1");
        assert_eq!(calls[0].position_1based, 10);
        assert_eq!(calls[1].contig, "chr1");
        assert_eq!(calls[1].position_1based, 20);
        assert_eq!(calls[2].contig, "chr2");
        assert_eq!(
            calls[2].format_keys,
            vec!["GT".to_string(), "AD".to_string(), "PL".to_string()]
        );
    }

    #[test]
    fn deterministic_sort_rejects_unknown_contig() {
        let mut calls = vec![GenotypedVariantCall {
            contig: "chrX".to_string(),
            position_1based: 1,
            reference: "A".to_string(),
            alternates: vec!["G".to_string()],
            format_keys: vec!["GT".to_string()],
        }];
        let err = sort_genotyped_calls_deterministic(&mut calls, &["chr1".to_string()])
            .expect_err("contig");
        assert!(err.to_string().contains("missing contig"));
    }

    #[test]
    fn java_combiner_merges_within_start_partition_only() {
        let loci = vec![
            ReferenceConfidenceLocus {
                position_1based: 100,
                gq: 0,
                dp: 0,
            },
            ReferenceConfidenceLocus {
                position_1based: 101,
                gq: 0,
                dp: 0,
            },
            ReferenceConfidenceLocus {
                position_1based: 102,
                gq: 12,
                dp: 2,
            },
        ];
        let blocks = build_gvcf_blocks_java_combiner(&loci, GATK_HC_DEFAULT_GQB).expect("blocks");
        assert_eq!(blocks.len(), 2);
        assert_eq!(blocks[0].start_1based, 100);
        assert_eq!(blocks[0].end_1based, 101);
        assert_eq!(blocks[0].min_rgq, 0);
        assert_eq!(blocks[1].start_1based, 102);
        assert_eq!(blocks[1].gq_band_upper, 13);
    }

    #[test]
    fn gvcf_hc_emit_semantics_absorb_zero_min_rgq_fringe() {
        let loci: Vec<_> = (1..=100)
            .map(|p| ReferenceConfidenceLocus {
                position_1based: p,
                gq: if p <= 80 { 0 } else { 3 },
                dp: 0,
            })
            .collect();
        let strict = build_gvcf_blocks(&loci, GATK_HC_DEFAULT_GQB).expect("blocks");
        assert!(strict.len() > 1);
        let merged = build_gvcf_blocks_with_semantics(
            &loci,
            GATK_HC_DEFAULT_GQB,
            gvcf_merge_semantics_hc_emit(),
        )
        .expect("blocks");
        assert_eq!(merged.len(), 1);
        assert_eq!(merged[0].end_1based, 100);
    }

    #[test]
    fn gvcf_blocks_merge_adjacent_positions_with_same_band() {
        let loci = vec![
            ReferenceConfidenceLocus {
                position_1based: 100,
                gq: 8,
                dp: 12,
            },
            ReferenceConfidenceLocus {
                position_1based: 101,
                gq: 9,
                dp: 10,
            },
            ReferenceConfidenceLocus {
                position_1based: 102,
                gq: 22,
                dp: 9,
            },
            ReferenceConfidenceLocus {
                position_1based: 103,
                gq: 25,
                dp: 15,
            },
        ];
        let blocks = build_gvcf_blocks(&loci, &[9, 19, 29, 99]).expect("blocks");
        assert_eq!(blocks.len(), 2);
        assert_eq!(
            blocks[0],
            GvcfBlock {
                start_1based: 100,
                end_1based: 101,
                gq_band_upper: 9,
                min_rgq: 8,
                min_dp: 10,
                max_dp: 12,
            }
        );
        assert_eq!(
            blocks[1],
            GvcfBlock {
                start_1based: 102,
                end_1based: 103,
                gq_band_upper: 29,
                min_rgq: 22,
                min_dp: 9,
                max_dp: 15,
            }
        );
    }

    #[test]
    fn gvcf_blocks_split_on_coordinate_gap_even_same_band() {
        let loci = vec![
            ReferenceConfidenceLocus {
                position_1based: 200,
                gq: 5,
                dp: 7,
            },
            ReferenceConfidenceLocus {
                position_1based: 202,
                gq: 6,
                dp: 8,
            },
        ];
        let blocks = build_gvcf_blocks(&loci, &[9, 99]).expect("blocks");
        assert_eq!(blocks.len(), 2);
        assert_eq!(blocks[0].start_1based, 200);
        assert_eq!(blocks[1].start_1based, 202);
    }

    #[test]
    fn gvcf_blocks_split_when_rgq_delta_exceeds_merge_semantics() {
        let loci = vec![
            ReferenceConfidenceLocus {
                position_1based: 300,
                gq: 5,
                dp: 10,
            },
            ReferenceConfidenceLocus {
                position_1based: 301,
                gq: 18, // same band upper with [19], but rgq jump is large
                dp: 9,
            },
        ];
        let blocks = build_gvcf_blocks_with_semantics(
            &loci,
            &[19, 99],
            GvcfMergeSemantics {
                max_rgq_delta_within_block: 4,
                java_gvcf_band_merge: false,
            },
        )
        .expect("blocks");
        assert_eq!(blocks.len(), 2);
    }

    #[test]
    fn gvcf_record_fields_end_matches_single_locus_boundary() {
        let block = GvcfBlock {
            start_1based: 500,
            end_1based: 500,
            gq_band_upper: 19,
            min_rgq: 11,
            min_dp: 8,
            max_dp: 8,
        };
        let rec = gvcf_block_to_record_fields(&block).expect("record");
        assert_eq!(rec.start_1based, 500);
        assert_eq!(rec.end_info, 500);
    }

    #[test]
    fn gvcf_record_fields_end_matches_multi_locus_boundary() {
        let block = GvcfBlock {
            start_1based: 700,
            end_1based: 729,
            gq_band_upper: 29,
            min_rgq: 20,
            min_dp: 5,
            max_dp: 14,
        };
        let rec = gvcf_block_to_record_fields(&block).expect("record");
        assert_eq!(rec.start_1based, 700);
        assert_eq!(rec.end_info, 729);
        assert_eq!(rec.min_dp, 5);
        assert_eq!(rec.max_dp, 14);
    }

    #[test]
    fn emit_mode_vcf_skips_reference_loci() {
        assert_eq!(
            decide_locus_emission(EmitMode::Vcf, false),
            LocusEmissionDecision::Skip
        );
        assert_eq!(
            decide_locus_emission(EmitMode::Vcf, true),
            LocusEmissionDecision::EmitVariantOnly
        );
    }

    #[test]
    fn emit_mode_gvcf_emits_reference_blocks_for_non_variant_loci() {
        assert_eq!(
            decide_locus_emission(EmitMode::Gvcf, false),
            LocusEmissionDecision::EmitReferenceBlock
        );
        assert_eq!(
            decide_locus_emission(EmitMode::Gvcf, true),
            LocusEmissionDecision::EmitVariantOnly
        );
    }

    #[test]
    fn emit_mode_bp_resolution_emits_reference_sites_per_base() {
        assert_eq!(
            decide_locus_emission(EmitMode::BpResolution, false),
            LocusEmissionDecision::EmitReferenceSite
        );
        assert_eq!(
            decide_locus_emission(EmitMode::BpResolution, true),
            LocusEmissionDecision::EmitVariantOnly
        );
    }

    #[test]
    fn no_variation_region_vcf_mode_emits_nothing() {
        let s = summarize_no_variation_region(EmitMode::Vcf, 25);
        assert_eq!(s.loci_total, 25);
        assert_eq!(s.reference_blocks_emitted, 0);
        assert_eq!(s.reference_sites_emitted, 0);
        assert_eq!(s.variant_sites_emitted, 0);
    }

    #[test]
    fn no_variation_region_gvcf_mode_prefers_block_style_emission() {
        let s = summarize_no_variation_region(EmitMode::Gvcf, 25);
        assert_eq!(s.loci_total, 25);
        assert!(s.reference_blocks_emitted > 0);
        assert_eq!(s.reference_sites_emitted, 0);
        assert_eq!(s.variant_sites_emitted, 0);
    }

    #[test]
    fn no_variation_region_bp_resolution_emits_per_base_reference_sites() {
        let s = summarize_no_variation_region(EmitMode::BpResolution, 25);
        assert_eq!(s.loci_total, 25);
        assert_eq!(s.reference_blocks_emitted, 0);
        assert_eq!(s.reference_sites_emitted, 25);
        assert_eq!(s.variant_sites_emitted, 0);
    }

    #[test]
    fn joint_compat_accepts_sorted_non_overlapping_records() {
        let records = vec![
            GvcfBlockRecordFields {
                start_1based: 100,
                end_info: 120,
                min_dp: 8,
                max_dp: 14,
                gq_band_upper: 19,
                min_rgq: 11,
            },
            GvcfBlockRecordFields {
                start_1based: 121,
                end_info: 130,
                min_dp: 9,
                max_dp: 16,
                gq_band_upper: 29,
                min_rgq: 20,
            },
        ];
        let summary = validate_joint_compatibility_gvcf_records(&records).expect("compat");
        assert_eq!(summary.records_total, 2);
        assert!(summary.compatible);
    }

    #[test]
    fn joint_compat_rejects_overlapping_records() {
        let records = vec![
            GvcfBlockRecordFields {
                start_1based: 200,
                end_info: 220,
                min_dp: 7,
                max_dp: 11,
                gq_band_upper: 19,
                min_rgq: 10,
            },
            GvcfBlockRecordFields {
                start_1based: 220,
                end_info: 230,
                min_dp: 8,
                max_dp: 12,
                gq_band_upper: 19,
                min_rgq: 12,
            },
        ];
        let err = validate_joint_compatibility_gvcf_records(&records).expect_err("overlap");
        assert!(err.to_string().contains("overlap"));
    }

    #[test]
    fn joint_compat_rejects_end_before_pos() {
        let records = vec![GvcfBlockRecordFields {
            start_1based: 300,
            end_info: 299,
            min_dp: 5,
            max_dp: 5,
            gq_band_upper: 9,
            min_rgq: 7,
        }];
        let err = validate_joint_compatibility_gvcf_records(&records).expect_err("end-before-pos");
        assert!(err.to_string().contains("END < POS"));
    }
}
