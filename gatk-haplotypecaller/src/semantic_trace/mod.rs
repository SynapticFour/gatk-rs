//! Optional structured semantic checkpoints for HaplotypeCaller (observe-only).
//! Enable with `GATK_RS_SEMANTIC_TRACE=<ndjson-path>`. When unset, every emit helper
//! is a single atomic load and returns immediately — no algorithm branching.
//! Schema: `gatk_rs.hc.semantic_trace/v1` (see [`schema`]).
//! Compare Java vs Rust with `scripts/parity/compare_semantic_trace.py`.

mod schema;
mod sink;

pub use schema::{
    RegionKey, SemanticStage, SemanticTraceEvent, SCHEMA_ID, TRACE_IMPL_JAVA, TRACE_IMPL_RUST,
};
pub use sink::{is_enabled, try_init_from_runtime, TraceSinkHandle};

use serde_json::{json, Value};

use crate::assembly_region_iterator::AssemblyRegion;
use crate::assembly_result_set::AssemblyResultSet;
use crate::haplotype::Haplotype;
use crate::hc_genotyping_engine::GenotypedSiteCall;
use crate::region_read_likelihood::RegionReadLikelihood;
use gatk_core::io::vcf::VcfRecord;

/// Emit one checkpoint when tracing is enabled.
#[inline]
pub fn checkpoint(stage: SemanticStage, region: Option<&RegionKey>, payload: Value) {
    if !is_enabled() {
        return;
    }
    sink::emit(stage, region.cloned(), payload);
}

/// Active-region boundary after the iterator materializes a region.
pub fn emit_active_region(region: &AssemblyRegion) {
    if !is_enabled() {
        return;
    }
    let key = RegionKey::from_assembly_region(region);
    checkpoint(
        SemanticStage::ActiveRegion,
        Some(&key),
        json!({
            "is_active": region.is_active,
            "extended_start": region.extended_start.get(),
            "extended_end": region.extended_end.get(),
            "read_count": region.reads.len(),
            "pileup_locus_count": region.pileup_loci.len(),
        }),
    );
}

/// Activity-profile cut summary (one event per popped profile region).
pub fn emit_activity_profile_cut(
    contig: &str,
    start: u64,
    end: u64,
    is_active: bool,
    padded_start: u64,
    padded_end: u64,
    extension: u32,
) {
    if !is_enabled() {
        return;
    }
    let key = RegionKey {
        contig: contig.to_string(),
        start,
        end,
        is_active: Some(is_active),
    };
    checkpoint(
        SemanticStage::ActivityProfile,
        Some(&key),
        json!({
            "padded_start": padded_start,
            "padded_end": padded_end,
            "extension": extension,
            "span_bp": end.saturating_sub(start).saturating_add(1),
        }),
    );
}

/// Assembly graph / result-set metrics after `call_region_assemble`.
pub fn emit_assembly_graph(region: &AssemblyRegion, assembly: &AssemblyResultSet) {
    if !is_enabled() {
        return;
    }
    let key = RegionKey::from_assembly_region(region);
    let n_ref = assembly
        .haplotypes
        .iter()
        .filter(|h| h.is_reference)
        .count();
    let n_alt = assembly.haplotypes.len().saturating_sub(n_ref);
    let kmer_sizes: Vec<usize> = {
        let mut ks: Vec<usize> = assembly
            .haplotypes
            .iter()
            .map(|h| h.kmer_size)
            .filter(|&k| k > 0)
            .collect();
        ks.sort_unstable();
        ks.dedup();
        ks
    };
    checkpoint(
        SemanticStage::AssemblyGraph,
        Some(&key),
        json!({
            "haplotype_count": assembly.haplotypes.len(),
            "ref_haplotype_count": n_ref,
            "alt_haplotype_count": n_alt,
            "variation_present": assembly.variation_present,
            "variation_event_count": assembly.variation_events().len(),
            "padded_reference_start": assembly.padded_reference_start_1based(),
            "reference_len": assembly.reference_bases().len(),
            "kmer_sizes": kmer_sizes,
        }),
    );
}

/// Reference haplotype path digest.
pub fn emit_reference_path(region: &AssemblyRegion, assembly: &AssemblyResultSet) {
    if !is_enabled() {
        return;
    }
    let key = RegionKey::from_assembly_region(region);
    let ref_hap = assembly.haplotypes.iter().find(|h| h.is_reference);
    let (digest, len, cigar) = match ref_hap {
        Some(h) => (
            bases_digest(&h.bases),
            h.bases.len(),
            h.cigar.as_ref().map(|c| c.to_gatk_string()),
        ),
        None => (
            bases_digest(assembly.reference_bases()),
            assembly.reference_bases().len(),
            None,
        ),
    };
    checkpoint(
        SemanticStage::ReferencePath,
        Some(&key),
        json!({
            "bases_digest": digest,
            "bases_len": len,
            "cigar": cigar,
            "pad_start": assembly.padded_reference_start_1based(),
        }),
    );
}

/// Candidate haplotype list (digests + scores; not full sequences).
pub fn emit_candidate_haplotypes(region: &AssemblyRegion, haplotypes: &[Haplotype]) {
    if !is_enabled() {
        return;
    }
    let key = RegionKey::from_assembly_region(region);
    let mut haps: Vec<Value> = haplotypes
        .iter()
        .enumerate()
        .map(|(i, h)| {
            json!({
                "index": i,
                "is_reference": h.is_reference,
                "bases_len": h.bases.len(),
                "bases_digest": bases_digest(&h.bases),
                "score": round6(h.score),
                "kmer_size": h.kmer_size,
                "cigar": h.cigar.as_ref().map(|c| c.to_gatk_string()),
            })
        })
        .collect();
    haps.sort_by(|a, b| {
        let da = a["bases_digest"].as_str().unwrap_or("");
        let db = b["bases_digest"].as_str().unwrap_or("");
        da.cmp(db)
    });
    checkpoint(
        SemanticStage::CandidateHaplotypes,
        Some(&key),
        json!({
            "count": haplotypes.len(),
            "haplotypes": haps,
        }),
    );
}

/// PairHMM / read×haplotype likelihood matrix summary.
pub fn emit_read_likelihoods(
    region: &AssemblyRegion,
    likelihoods: &[RegionReadLikelihood],
    n_haplotypes: usize,
) {
    if !is_enabled() {
        return;
    }
    let key = RegionKey::from_assembly_region(region);
    let n_reads = likelihoods
        .iter()
        .map(|r| r.read_index.get())
        .collect::<std::collections::BTreeSet<_>>()
        .len();
    let mut sample: Vec<Value> = likelihoods
        .iter()
        .take(64)
        .map(|r| {
            json!({
                "read": r.read_index.get(),
                "hap": r.haplotype_index.get(),
                "ll": round6(r.log10_likelihood),
            })
        })
        .collect();
    // Stable order for compare.
    sample.sort_by(|a, b| {
        let ra = a["read"].as_u64().unwrap_or(0);
        let rb = b["read"].as_u64().unwrap_or(0);
        ra.cmp(&rb).then_with(|| {
            a["hap"]
                .as_u64()
                .unwrap_or(0)
                .cmp(&b["hap"].as_u64().unwrap_or(0))
        })
    });
    let sum_ll: f64 = likelihoods.iter().map(|r| r.log10_likelihood).sum();
    checkpoint(
        SemanticStage::ReadLikelihoods,
        Some(&key),
        json!({
            "matrix_cells": likelihoods.len(),
            "n_reads": n_reads,
            "n_haplotypes": n_haplotypes,
            "sum_ll": round6(sum_ll),
            "sample_cells": sample,
        }),
    );
}

/// Genotype likelihood / site call summary after assignGenotypeLikelihoods.
pub fn emit_genotype_likelihoods(region: &AssemblyRegion, calls: &[GenotypedSiteCall]) {
    if !is_enabled() {
        return;
    }
    let key = RegionKey::from_assembly_region(region);
    let mut sites: Vec<Value> = calls
        .iter()
        .map(|c| {
            json!({
                "pos": c.event.start_1based.get(),
                "ref": c.event.ref_allele,
                "alt": c.event.alt_allele,
                "pl": c.genotype.format.pl_as_i32(),
                "gq": c.genotype.format.gq.as_i32(),
                "ad": c.genotype.format.ad_as_i32(),
                "dp": c.genotype.format.dp.as_i32(),
                "gl": c
                    .genotype
                    .genotype_log10_likelihoods
                    .iter()
                    .map(|x| round6(*x))
                    .collect::<Vec<_>>(),
            })
        })
        .collect();
    sites.sort_by(|a, b| {
        a["pos"]
            .as_u64()
            .unwrap_or(0)
            .cmp(&b["pos"].as_u64().unwrap_or(0))
            .then_with(|| {
                a["ref"]
                    .as_str()
                    .unwrap_or("")
                    .cmp(b["ref"].as_str().unwrap_or(""))
            })
            .then_with(|| {
                a["alt"]
                    .as_str()
                    .unwrap_or("")
                    .cmp(b["alt"].as_str().unwrap_or(""))
            })
    });
    checkpoint(
        SemanticStage::GenotypeLikelihoods,
        Some(&key),
        json!({
            "site_count": calls.len(),
            "sites": sites,
        }),
    );
}

/// VCF records emitted for one active region (or inactive RCM batch).
pub fn emit_vcf_emission(region: Option<&AssemblyRegion>, records: &[VcfRecord]) {
    if !is_enabled() {
        return;
    }
    let key = region.map(RegionKey::from_assembly_region);
    let mut sites: Vec<Value> = records
        .iter()
        .map(|r| {
            json!({
                "chrom": r.chromosome,
                "pos": r.position,
                "ref": r.reference,
                "alt": r.alternate,
                "qual": r.quality.map(round6),
                "filter": r.filter,
            })
        })
        .collect();
    sites.sort_by(|a, b| {
        a["chrom"]
            .as_str()
            .unwrap_or("")
            .cmp(b["chrom"].as_str().unwrap_or(""))
            .then_with(|| {
                a["pos"]
                    .as_u64()
                    .unwrap_or(0)
                    .cmp(&b["pos"].as_u64().unwrap_or(0))
            })
            .then_with(|| {
                a["ref"]
                    .as_str()
                    .unwrap_or("")
                    .cmp(b["ref"].as_str().unwrap_or(""))
            })
    });
    checkpoint(
        SemanticStage::VcfEmission,
        key.as_ref(),
        json!({
            "record_count": records.len(),
            "sites": sites,
        }),
    );
}

/// Inactive reference-confidence model path (no assembly).
pub fn emit_inactive_rcm(region: &AssemblyRegion, locus_count: usize) {
    if !is_enabled() {
        return;
    }
    let key = RegionKey::from_assembly_region(region);
    checkpoint(
        SemanticStage::InactiveRcm,
        Some(&key),
        json!({
            "locus_count": locus_count,
        }),
    );
}

/// Convenience: emit post-assemble checkpoints (graph + ref path + candidates).
pub fn emit_post_assemble(region: &AssemblyRegion, assembly: &AssemblyResultSet) {
    if !is_enabled() {
        return;
    }
    emit_assembly_graph(region, assembly);
    emit_reference_path(region, assembly);
    emit_candidate_haplotypes(region, &assembly.haplotypes);
}

fn round6(x: f64) -> f64 {
    if !x.is_finite() {
        return x;
    }
    (x * 1_000_000.0).round() / 1_000_000.0
}

/// Stable FNV-1a digest of haplotype / reference bases (not cryptographic).
pub fn bases_digest(bases: &[u8]) -> String {
    let mut h: u64 = 14695981039346656037;
    for &b in bases {
        h ^= u64::from(b);
        h = h.wrapping_mul(1099511628211);
    }
    format!("{:016x}:{}", h, bases.len())
}

impl RegionKey {
    pub fn from_assembly_region(region: &AssemblyRegion) -> Self {
        Self {
            // CLONE: needed because owned contig id for output record.
            contig: region.contig.clone(),
            start: region.start.get(),
            end: region.end.get(),
            is_active: Some(region.is_active),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::genome_loc::GenomePosition;
    use crate::runtime_config::RuntimeConfig;
    use std::sync::Mutex;
    use tempfile::NamedTempFile;

    static ENV_LOCK: Mutex<()> = Mutex::new(());

    fn sample_region() -> AssemblyRegion {
        AssemblyRegion {
            contig: "20".into(),
            start: GenomePosition::new_1based(100),
            end: GenomePosition::new_1based(200),
            is_active: true,
            extended_start: GenomePosition::new_1based(50),
            extended_end: GenomePosition::new_1based(250),
            extension: 50,
            reads: Vec::new(),
            read_qnames: Vec::new(),
            reference: crate::reference_context::ReferenceContext::empty(),
            features: crate::feature_context::FeatureContext::empty(),
            pileup_loci: Vec::new(),
        }
    }

    #[test]
    fn disabled_by_default_is_noop() {
        let _g = ENV_LOCK.lock().unwrap();
        std::env::remove_var("GATK_RS_SEMANTIC_TRACE");
        sink::reset_for_tests();
        try_init_from_runtime(&RuntimeConfig::from_env());
        assert!(!is_enabled());
        emit_active_region(&sample_region());
    }

    #[test]
    fn writes_ndjson_when_enabled() {
        let _g = ENV_LOCK.lock().unwrap();
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path().to_string_lossy().into_owned();
        std::env::set_var("GATK_RS_SEMANTIC_TRACE", &path);
        sink::reset_for_tests();
        try_init_from_runtime(&RuntimeConfig::from_env());
        assert!(is_enabled());

        emit_active_region(&sample_region());
        emit_activity_profile_cut("20", 100, 200, true, 50, 250, 100);
        sink::flush_for_tests();

        let text = std::fs::read_to_string(&path).unwrap();
        let lines: Vec<&str> = text.lines().filter(|l| !l.is_empty()).collect();
        assert_eq!(lines.len(), 2);
        let ev: SemanticTraceEvent = serde_json::from_str(lines[0]).unwrap();
        assert_eq!(ev.schema, SCHEMA_ID);
        assert_eq!(ev.stage, SemanticStage::ActiveRegion);
        assert_eq!(ev.impl_name, TRACE_IMPL_RUST);

        std::env::remove_var("GATK_RS_SEMANTIC_TRACE");
        sink::reset_for_tests();
    }

    #[test]
    fn bases_digest_stable() {
        assert_eq!(bases_digest(b"ACGT"), bases_digest(b"ACGT"));
        assert_ne!(bases_digest(b"ACGT"), bases_digest(b"ACGG"));
    }
}
