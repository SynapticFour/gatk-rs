//! ASM-1 cluster graph: dangling merge + k-best indel CIGAR (no full P12 trace).
//! Run: `P12_REFERENCE=…/hs37d5.simple.fa cargo test -p gatk-haplotypecaller p12_cluster_asm1 --release`

use gatk_core::reference::ReferenceWindowCache;
use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::assembly_dangling_recovery::DanglingRecoveryParams;
use gatk_haplotypecaller::assembly_region_finalize::{
    assembly_reference_read, finalize_region_reads_for_assembly,
    gatk_min_tail_quality_for_assembly, records_to_assembly_reads,
};
use gatk_haplotypecaller::event_map::variation_events_for_haplotype;
use gatk_haplotypecaller::read_event_discovery::is_p12_cluster_coupled_indel;
use gatk_haplotypecaller::read_threading_assembler::{
    build_threading_graph_for_haplotype_dump, extract_rt_haplotypes_before_remove_paths,
    AssemblyScoringContext,
};
use gatk_haplotypecaller::{
    assemble_reads, audit_threading_dangling_recovery, call_disposition, flatten_assembly_regions,
    traverse_assembly_region_walker, AssembleReadsArgs, AssemblyRegionCallDisposition,
    ReadFilterParams, ReadThreadingAssemblerArgs, WalkerTraversalConfig,
};
use std::path::Path;

fn p12_paths() -> Option<(std::path::PathBuf, std::path::PathBuf)> {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("..");
    let ref_path = std::env::var("P12_REFERENCE")
        .map(std::path::PathBuf::from)
        .unwrap_or_else(|_| root.join("parity/realworld/assets/hs37d5.simple.fa"));
    let ref_path = if ref_path.is_absolute() {
        ref_path
    } else {
        root.join(ref_path)
    };
    let bam = root.join("parity/realworld/na12878_20k_b37/NA12878_20k.b37.bam");
    if !ref_path.is_file() || !bam.is_file() {
        return None;
    }
    Some((ref_path, bam))
}

#[test]
fn p12_cluster_asm1_dangling_and_kbest() {
    let Some((ref_fasta, bam)) = p12_paths() else {
        eprintln!("skip: set P12_REFERENCE");
        return;
    };
    let dict = SequenceDictionary::from_fasta_path(&ref_fasta).expect("dict");
    let specs = parse_intervals_cli_string(&dict, "2:92307228-92307400").expect("interval");
    let filters = ReadFilterParams::gatk_standard_hc();
    let cfg = WalkerTraversalConfig::gatk_haplotype_caller_production(100);
    let walk = traverse_assembly_region_walker(&dict, &specs, &ref_fasta, &bam, &filters, &cfg)
        .expect("walk");
    let regions = flatten_assembly_regions(&walk);
    let region = regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= 92307324
                && r.end.get() >= 92307327
        })
        .expect("cluster active region");

    let mut ref_cache = ReferenceWindowCache::new(ref_fasta.clone(), 4);
    let reference = assembly_reference_read(&dict, &mut ref_cache, region).expect("ref");
    let finalized = finalize_region_reads_for_assembly(
        &region.reads,
        region,
        true,
        gatk_min_tail_quality_for_assembly(10),
        false,
    );
    let reads = records_to_assembly_reads(&finalized);
    let args = ReadThreadingAssemblerArgs::default();
    for &kmer in &[25usize, 85] {
        let audit = audit_threading_dangling_recovery(&reference, &reads, kmer, &args, true, true)
            .expect("audit")
            .expect("graph");
        eprintln!(
            "asm1_dangling\tk={kmer}\tbefore={} after={} tails={}/{}",
            audit.edges_before, audit.edges_after, audit.tails_recovered, audit.tails_attempted
        );
        assert!(
            audit.edges_after >= audit.edges_before,
            "k={kmer}: dangling should not shrink edge set"
        );
    }
    let audit25 = audit_threading_dangling_recovery(&reference, &reads, 25, &args, true, true)
        .expect("audit")
        .expect("graph");
    assert!(
        audit25.tails_recovered >= 1,
        "ASM-1 exit: expect ≥1 tail merge on cluster graph at k=25"
    );
    let audit85 = audit_threading_dangling_recovery(&reference, &reads, 85, &args, true, true)
        .expect("audit")
        .expect("graph");
    assert!(
        audit85.tails_recovered >= 1,
        "ASM-1 k=85: expect dangling tail merge on cluster graph (got {}/{} recovered)",
        audit85.tails_recovered,
        audit85.tails_attempted
    );
    let dangling = DanglingRecoveryParams::from_assembler_args(&args);
    let mut probe_args = args.clone();
    probe_args.recover_dangling_branches = false;
    probe_args.remove_paths_not_connected_to_ref = false;
    let graph =
        build_threading_graph_for_haplotype_dump(&reference, &reads, 85, &probe_args, true, true)
            .expect("probe graph")
            .expect("k=85 graph");
    let tail_probes: Vec<_> = graph.probe_dangling_tail_failures(&dangling);
    assert!(
        tail_probes.iter().any(|(_, _, r)| r != "no_alt_path"),
        "ASM-1 k=85: alt path must be found (got {:?})",
        tail_probes
    );
    eprintln!(
        "asm1_k85_probe\ttails={}/{} probe={:?}",
        audit85.tails_recovered, audit85.tails_attempted, tail_probes
    );

    let assembly = assemble_reads(region, &dict, &mut ref_cache, &AssembleReadsArgs::default())
        .expect("assemble");
    let indel_haps: Vec<_> = assembly
        .haplotypes
        .iter()
        .filter(|h| {
            !h.is_reference
                && h.cigar
                    .as_ref()
                    .is_some_and(|c| c.elements.iter().any(|e| e.operator.is_indel()))
        })
        .collect();
    eprintln!(
        "asm1_assemble\thaps={} indel_alt_haps={}",
        assembly.haplotypes.len(),
        indel_haps.len()
    );
    for h in &indel_haps {
        eprintln!(
            "asm1_indel_hap\tcigar={}",
            h.cigar
                .as_ref()
                .map(|c| c.to_gatk_string())
                .unwrap_or_default()
        );
    }
    let coupled =
        gatk_haplotypecaller::read_event_discovery::cluster_coupled_events_from_assembly_haplotypes(
            &assembly,
            "2",
            region.start.get(),
            region.end.get(),
        );
    eprintln!(
        "asm1_cluster_coupled\tcount={} events={:?}",
        coupled.len(),
        coupled
    );
    assert!(
        coupled.len() >= 2,
        "ASM-8: assembly must emit coupled cluster indels from graph path (got {coupled:?})"
    );

    let pad = assembly.padded_reference_start_1based();
    let ref_hap = assembly
        .haplotypes
        .iter()
        .find(|h| h.is_reference)
        .expect("ref");
    let full_ref = assembly.reference_bases();
    for &kmer in &[25usize, 85] {
        let batch =
            extract_rt_haplotypes_before_remove_paths(&reference, &reads, &args, kmer, true, true)
                .expect("rt batch");
        let mut coupled_n = 0usize;
        for (i, h) in batch.iter().filter(|h| !h.is_reference).enumerate() {
            for e in variation_events_for_haplotype(h, ref_hap, full_ref, pad, 1, "2") {
                if is_p12_cluster_coupled_indel(&e) {
                    coupled_n += 1;
                    eprintln!(
                        "rt_k{kmer}_hap{i}\t{} {}/{} cigar={}",
                        e.start_1based.get(),
                        e.ref_allele,
                        e.alt_allele,
                        h.cigar
                            .as_ref()
                            .map(|c| c.to_gatk_string())
                            .unwrap_or_default()
                    );
                }
            }
        }
        eprintln!("rt_k{kmer}_coupled_hits\t{coupled_n}");
    }

    let mut asm_scored = args.clone();
    asm_scored.scoring = Some(AssemblyScoringContext {
        padded_reference_start_1based: pad,
        active_start_1based: region.start.get(),
        active_end_1based: region.end.get(),
        contig: "2".into(),
    });
    let mut ref_cache2 = ReferenceWindowCache::new(ref_fasta.clone(), 4);
    let assembly_scored = assemble_reads(
        region,
        &dict,
        &mut ref_cache2,
        &AssembleReadsArgs {
            assembler: asm_scored,
            ..AssembleReadsArgs::default()
        },
    )
    .expect("scored assemble");
    let coupled2 =
        gatk_haplotypecaller::read_event_discovery::cluster_coupled_events_from_assembly_haplotypes(
            &assembly_scored,
            "2",
            region.start.get(),
            region.end.get(),
        );
    eprintln!(
        "asm1_scored_coupled\tcount={} kmer={}",
        coupled2.len(),
        assembly_scored.kmer_size_for_dump()
    );

    let read_coupled =
        gatk_haplotypecaller::read_event_discovery::discover_p12_cluster_coupled_events_from_reads(
            &region.reads,
            full_ref,
            pad,
            region.start.get(),
            region.end.get(),
            "2",
        );
    eprintln!(
        "asm1_read_coupled\tcount={} events={:?}",
        read_coupled.len(),
        read_coupled
    );
}
