//! Production assemble shortcut: RT before_remove on configured k-mers before SeqGraph.

use crate::assembly::AssemblyRead;
use crate::cigar::{Cigar, CigarOperator};
use crate::event_map::EventMap;
use crate::haplotype::Haplotype;
use crate::read_event_discovery::refresh_alt_haplotype_indel_cigars;
use crate::read_threading_assembler::{
    finalize_assembly_haplotypes, haplotypes_have_alt_bases, just_reference_result,
    merge_rt_kbest_pre_remove_paths_at_kmer, AssemblyResult, AssemblyStatus,
    ReadThreadingAssemblerArgs,
};
use gatk_common::GatkResult;

/// Non-P12 production shortcut: RT before_remove on configured k-mers (largest first).
///
/// # Invariants
/// Returns `Some` only with alt+ref haplotypes; `None` falls through to SeqGraph.
/// Never runs on P12 / L-gate or tiny L2 k-mer sets (`k < 10`).
/// # Java equivalence
/// Same before_remove RT extract as
/// [`crate::read_threading_assembler::supplement_p12_cluster_coupled_haplotypes`], ordered
/// to avoid a redundant SeqGraph prelude when that extract already yields alts.
pub(crate) fn try_rt_configured_alts_before_seq_graph(
    reference: &AssemblyRead,
    reads: &[AssemblyRead],
    args: &ReadThreadingAssemblerArgs,
) -> GatkResult<Option<AssemblyResult>> {
    // Production HC sets `scoring` under `strict_java_assembly`. Without it, SeqGraph +
    // merge_rt is the sole alt source (supplement is a no-op) — do not short-circuit.
    let Some(ctx) = args.scoring.as_ref() else {
        return Ok(None);
    };
    if ctx.overlaps_p12_cluster() || ctx.overlaps_p12_l_gate_interval() {
        return Ok(None);
    }
    let mut kmers: Vec<usize> = args.kmer_sizes.iter().copied().collect();
    kmers.sort_unstable();
    kmers.dedup();
    let production_kmer_set = kmers.iter().any(|&k| k >= 25) && !kmers.iter().any(|&k| k < 10);
    if !production_kmer_set {
        return Ok(None);
    }
    // Dense ci-subset: early-stop is almost always at k=25 after a useless k=10 probe.
    kmers.sort_unstable_by(|a, b| b.cmp(a));

    crate::runtime_config::rss_trace_checkpoint(
        "rt_first_configured_begin",
        &format!("kmers={kmers:?}"),
    );

    let ref_bytes = reference.bases.as_slice();
    let pad = ctx.padded_reference_start_1based;
    let sw = &args.haplotype_to_reference_sw;
    let mut hit_kmer = kmers.first().copied().unwrap_or(25);
    let mut haplotypes = just_reference_result(hit_kmer, reference).haplotypes;

    for &kmer_size in &kmers {
        if crate::runtime_config::hc_rss_abort_triggered() {
            break;
        }
        crate::runtime_config::rss_trace_checkpoint(
            "rt_first_configured_kmer",
            &format!("kmer={kmer_size}"),
        );
        // Reuse merge_rt extract/dedup (owned-key clones stay in that helper's ratchet count).
        merge_rt_kbest_pre_remove_paths_at_kmer(
            reference,
            reads,
            args,
            &[],
            &mut haplotypes,
            Some(kmer_size),
        )?;
        if !haplotypes.iter().any(|h| !h.is_reference) || haplotypes.len() <= 1 {
            continue;
        }
        refresh_alt_haplotype_indel_cigars(&mut haplotypes, ref_bytes, pad, sw);
        crate::smith_waterman::release_sw_tls_scratch();
        let mut ref_hap = Haplotype::new(ref_bytes, true);
        let mut ref_cigar = Cigar::new();
        ref_cigar.push(ref_hap.bases.len(), CigarOperator::Match);
        ref_hap.cigar = Some(ref_cigar);
        let status = finalize_assembly_haplotypes(
            &mut haplotypes,
            &ref_hap,
            args.ensure_reference_in_result,
        );
        // `normalize_ref_equivalent_haplotypes` may collapse KBest `is_reference=false`
        // paths that match ref bases — do not skip SeqGraph on a phantom hit (L2
        // `g2-subset-live` p11: haplotype_count=1 vs Java 2).
        if !matches!(status, AssemblyStatus::AssembledSomeVariation)
            || !haplotypes_have_alt_bases(&haplotypes, ref_bytes)
        {
            haplotypes = just_reference_result(kmer_size, reference).haplotypes;
            continue;
        }
        hit_kmer = kmer_size;
        let event_maps = haplotypes
            .iter()
            .map(|h| {
                let rh = Haplotype::new(ref_bytes, true);
                EventMap::from_haplotype_and_reference(h, &rh, &rh.bases, 1, 0)
            })
            .collect();
        crate::runtime_config::rss_trace_checkpoint(
            "rt_first_configured_hit",
            &format!("kmer={kmer_size} haps={}", haplotypes.len()),
        );
        return Ok(Some(AssemblyResult {
            status,
            kmer_size: hit_kmer,
            haplotypes,
            event_maps,
        }));
    }

    crate::runtime_config::rss_trace_checkpoint(
        "rt_first_configured_miss",
        &format!("haps={}", haplotypes.len()),
    );
    Ok(None)
}
