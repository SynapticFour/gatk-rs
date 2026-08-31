//! Diagnostic RT before_remove extract on configured k-mers (retired Peak-RSS skip).
//!
//! Production assemble no longer uses a successful extract to bypass SeqGraph
//! (Java 4.4 `assembleKmerGraphsAndHaplotypeCall`). Kept so holdouts can measure
//! whether RT configured alts *would* have succeeded.

use crate::assembly::AssemblyRead;
use crate::cigar::{Cigar, CigarOperator};
use crate::haplotype::Haplotype;
use crate::read_threading_assembler::{
    finalize_assembly_haplotypes, haplotypes_have_alt_bases, just_reference_result,
    merge_rt_kbest_pre_remove_paths_at_kmer, AssemblyResult, AssemblyStatus,
    ReadThreadingAssemblerArgs,
};
use gatk_common::GatkResult;

/// Diagnostic: RT before_remove on configured k-mers (largest first).
///
/// # Invariants
/// Returns `Some` only with alt+ref haplotypes. Never used to skip SeqGraph in production.
/// Never runs on P12 / L-gate or tiny L2 k-mer sets (`k < 10`).
/// On miss, marks each configured k-mer empty in [`crate::rt_region_cache`].
/// # Java equivalence
/// Same before_remove RT extract as
/// [`crate::read_threading_assembler::supplement_p12_cluster_coupled_haplotypes`].
/// Java 4.4 has no “skip SeqGraph when this extract already has alts” branch.
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
            crate::rt_region_cache::mark_configured_kmer_empty(kmer_size);
            continue;
        }
        // Defer indel CIGAR refresh to assemble_reads_with_finalized.
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
            crate::rt_region_cache::mark_configured_kmer_empty(kmer_size);
            haplotypes = just_reference_result(kmer_size, reference).haplotypes;
            continue;
        }
        hit_kmer = kmer_size;
        let event_maps = Vec::new(); // production rebuilds EventMap from CIGARs later
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::assembly::AssemblyRead;
    use crate::read_threading_assembler::{AssemblyScoringContext, ReadThreadingAssemblerArgs};

    fn dummy_read(len: usize) -> AssemblyRead {
        AssemblyRead {
            bases: vec![b'A'; len],
            base_quals: vec![30; len],
        }
    }

    /// 6R.54 characterization (coordinate-free).
    ///
    /// Java 4.4 `assembleKmerGraphsAndHaplotypeCall` always runs SeqGraph
    /// `findBestPaths` (ReadThreadingAssembler.java ~213–228). Rust's RT-first
    /// shortcut is not a Java helper. Without `scoring`, it must not fire.
    #[test]
    fn rt_first_does_not_run_without_scoring_seqgraph_is_the_java_path() {
        let args = ReadThreadingAssemblerArgs::default();
        assert!(args.scoring.is_none());
        let reference = dummy_read(80);
        let reads = [dummy_read(80)];
        let out = try_rt_configured_alts_before_seq_graph(&reference, &reads, &args).unwrap();
        assert!(
            out.is_none(),
            "Java has no RT-before-SeqGraph skip; scoring=None must fall through"
        );
    }

    /// Java 4.4 `generateSeqGraph = !useLinkedDeBruijnGraph` with default `false`.
    #[test]
    fn default_use_seq_graph_matches_java_generate_seq_graph() {
        assert!(
            ReadThreadingAssemblerArgs::default().use_seq_graph,
            "Java default HC always enters assembleKmerGraphsAndHaplotypeCall"
        );
    }

    /// Documents when the retired Peak-RSS helper is eligible to *extract* RT alts.
    /// Production no longer skips SeqGraph on that extract (6R.56).
    #[test]
    fn rt_first_eligible_only_on_production_kmer_set_with_non_p12_scoring() {
        let mut args = ReadThreadingAssemblerArgs::default();
        args.scoring = Some(AssemblyScoringContext {
            padded_reference_start_1based: 1000,
            active_start_1based: 1100,
            active_end_1based: 1200,
            contig: "synth".into(),
        });
        let ctx = args.scoring.as_ref().unwrap();
        assert!(!ctx.overlaps_p12_cluster());
        assert!(!ctx.overlaps_p12_l_gate_interval());
        let production_kmer_set =
            args.kmer_sizes.iter().any(|&k| k >= 25) && !args.kmer_sizes.iter().any(|&k| k < 10);
        assert!(
            production_kmer_set,
            "default k=[10,25] is the production set the retired helper could probe"
        );
    }

    fn greedy_unique_k25(len: usize) -> Vec<u8> {
        const K: usize = 25;
        let mut s: Vec<u8> = (0..K).map(|i| b"ACGT"[i % 4]).collect();
        let mut seen = std::collections::HashSet::new();
        seen.insert(s.clone());
        while s.len() < len {
            let mut placed = false;
            for &b in &[b'A', b'C', b'G', b'T'] {
                s.push(b);
                let km = s[s.len() - K..].to_vec();
                if seen.insert(km) {
                    placed = true;
                    break;
                }
                s.pop();
            }
            assert!(placed, "could not extend unique k=25 sequence");
        }
        s
    }

    fn assembly_read(bases: &[u8]) -> AssemblyRead {
        AssemblyRead {
            bases: bases.to_vec(),
            base_quals: vec![30; bases.len()],
        }
    }

    /// 6R.56: Java 4.4 always runs SeqGraph `findBestPaths` on default HC.
    /// RT configured alts may succeed; that success must not terminate assemble.
    #[test]
    fn rt_configured_alt_success_does_not_skip_seq_graph() {
        let mut ref_bases = greedy_unique_k25(80);
        let snp = 40usize;
        let alt_base = if ref_bases[snp] == b'A' { b'C' } else { b'A' };
        let mut alt_bases = ref_bases.clone();
        alt_bases[snp] = alt_base;
        let snp2 = 55usize;
        let alt2 = if ref_bases[snp2] == b'T' { b'G' } else { b'T' };
        alt_bases[snp2] = alt2;

        let reference = assembly_read(&ref_bases);
        let mut reads = Vec::new();
        for _ in 0..6 {
            reads.push(assembly_read(&alt_bases));
        }
        for _ in 0..2 {
            reads.push(assembly_read(&ref_bases));
        }

        let mut args = ReadThreadingAssemblerArgs::default();
        args.dangling_java_exact = true;
        args.scoring = Some(AssemblyScoringContext {
            padded_reference_start_1based: 1,
            active_start_1based: 10,
            active_end_1based: 70,
            contig: "synth".into(),
        });

        let rt = try_rt_configured_alts_before_seq_graph(&reference, &reads, &args).unwrap();
        assert!(
            rt.is_some(),
            "RT configured extract must succeed so the test exercises the retired skip"
        );
        let rt_has_snp1 = rt
            .as_ref()
            .unwrap()
            .haplotypes
            .iter()
            .any(|h| h.bases.get(snp) == Some(&alt_base));

        crate::read_threading_assembler::seq_graph_assemble_probe::reset();
        let assembled =
            crate::read_threading_assembler::assemble_from_ref_and_reads(&reference, &reads, &args)
                .unwrap();
        assert!(
            crate::read_threading_assembler::seq_graph_assemble_probe::get() >= 1,
            "Java 4.4 assembleKmerGraphsAndHaplotypeCall must still run SeqGraph"
        );
        assert!(
            assembled.haplotypes.len() >= 2,
            "assembled haplotypes must include ref+alt"
        );
        let has_snp1 = assembled
            .haplotypes
            .iter()
            .any(|h| h.bases.get(snp) == Some(&alt_base));
        let has_snp2 = assembled
            .haplotypes
            .iter()
            .any(|h| h.bases.get(snp2) == Some(&alt2));
        assert!(has_snp1, "first SNP must survive SeqGraph assemble");
        assert!(
            has_snp2,
            "second SNP must survive SeqGraph even if RT extract already had an alt (rt_snp1={rt_has_snp1})"
        );
    }
}
