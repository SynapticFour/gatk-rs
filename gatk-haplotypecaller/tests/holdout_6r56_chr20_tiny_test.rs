//! 6R.56 holdout: after restoring Java SeqGraph control flow, `20:29455902 G/A`
//! must survive assemble / EventMap. Not a production locus pin.
//!
//! Skipped unless `HOLDOUT_6R56=1`.
//!
//! ```text
//! HOLDOUT_6R56=1 cargo test -p gatk-haplotypecaller --test holdout_6r56_chr20_tiny_test -- --nocapture
//! ```

use gatk_core::reference::{parse_intervals_cli_string, ReferenceWindowCache, SequenceDictionary};
use gatk_haplotypecaller::assembly_region_finalize::{
    assembly_reference_read, create_graph_reference_read, records_to_assembly_reads,
};
use gatk_haplotypecaller::hc_allele_mapping::hap_base_at_ref_locus;
use gatk_haplotypecaller::hc_genotyping_engine::DEFAULT_STAND_EMIT_CONFIDENCE;
use gatk_haplotypecaller::read_threading_assembler::{
    AssemblyScoringContext, DEFAULT_NUM_BEST_HAPLOTYPES_PER_GRAPH,
};
use gatk_haplotypecaller::{
    assemble_from_ref_and_reads, assemble_reads_with_finalized, call_disposition,
    diagnostic_rt_first_skip_seq_graph_kmer, flatten_assembly_regions,
    traverse_assembly_region_walker, try_emit_call_region_variants, AssemblyRegionCallDisposition,
    CallRegionArgs, HaplotypeCallerEngine, ReadFilterParams, WalkerTraversalConfig,
};
use serde_json::json;
use std::path::{Path, PathBuf};

const INTERVAL: &str = "20:29455000-29456500";
const BAM_REL: &str = "parity/giab/runs/local-pairhmm-diff/HG001.20-29455000-29456500.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const ACTIVE: (u64, u64) = (29_455_745, 29_455_993);
const TARGET: u64 = 29_455_902;
const TARGET_ALT: u8 = b'A';
const K: usize = DEFAULT_NUM_BEST_HAPLOTYPES_PER_GRAPH;

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

#[test]
fn holdout_6r56_seqgraph_retains_a_that_rt_shortcut_dropped() {
    if std::env::var("HOLDOUT_6R56").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R56=1");
        return;
    }
    let root = repo_root();
    let ref_fasta = root.join(REF_REL);
    let bam = root.join(BAM_REL);
    assert!(ref_fasta.is_file() && bam.is_file());

    let dict = SequenceDictionary::from_fasta_path(&ref_fasta).expect("dict");
    let specs = parse_intervals_cli_string(&dict, INTERVAL).expect("interval");
    let walk = traverse_assembly_region_walker(
        &dict,
        &specs,
        &ref_fasta,
        &bam,
        &ReadFilterParams::gatk_standard_hc(),
        &WalkerTraversalConfig::gatk_haplotype_caller_production(100),
    )
    .expect("walk");
    let regions = flatten_assembly_regions(&walk);
    let region = regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() == ACTIVE.0
                && r.end.get() == ACTIVE.1
        })
        .expect("ActiveFull");

    let args = CallRegionArgs::strict_java();
    let mut ref_cache = ReferenceWindowCache::new(ref_fasta.clone(), 4);
    let mut owned = region.clone();
    let assembled =
        assemble_reads_with_finalized(&mut owned, &dict, &mut ref_cache, &args.assemble)
            .expect("assemble");
    let pad = assembled.assembly.padded_reference_start_1based();
    let padded_ref = assembly_reference_read(&dict, &mut ref_cache, region).expect("pad ref");
    let graph_ref = create_graph_reference_read(&padded_ref, region, &dict);
    let graph_reads = records_to_assembly_reads(&assembled.finalized_reads);

    let mut assembler = args.assemble.assembler.clone();
    assembler.dangling_java_exact = true;
    assembler.scoring = Some(AssemblyScoringContext {
        padded_reference_start_1based: region.extended_start.get(),
        active_start_1based: region.start.get(),
        active_end_1based: region.end.get(),
        contig: region.contig.clone(),
    });

    let would_rt = diagnostic_rt_first_skip_seq_graph_kmer(&graph_ref, &graph_reads, &assembler)
        .expect("diagnostic");
    let prod = assemble_from_ref_and_reads(&graph_ref, &graph_reads, &assembler).expect("prod");
    let has_a = prod
        .haplotypes
        .iter()
        .any(|h| hap_base_at_ref_locus(h, pad, TARGET) == Some(TARGET_ALT));
    let event_ga = assembled
        .assembly
        .variation_events()
        .iter()
        .any(|e| e.start_1based.get() == TARGET && e.ref_allele == "G" && e.alt_allele == "A");

    let call = HaplotypeCallerEngine::call_region(region, &dict, &ref_fasta, &args).expect("call");
    let outcome = call.expect("Some");
    let vcf =
        try_emit_call_region_variants(region, &outcome, "NA12878", DEFAULT_STAND_EMIT_CONFIDENCE)
            .expect("emit");
    let vcf_ga = vcf.iter().any(|r| {
        r.position == TARGET && r.reference == "G" && r.alternate.iter().any(|a| a == "A")
    });

    let doc = json!({
        "k": K,
        "rt_configured_alts_would_succeed_k": would_rt,
        "production_n_haps": prod.haplotypes.len(),
        "haplotypes_have_A": has_a,
        "assemble_reads_event_G_A": event_ga,
        "vcf_has_G_A": vcf_ga,
        "first_loss_moved_past_rt_shortcut": has_a && event_ga,
    });
    println!("{}", serde_json::to_string_pretty(&doc).unwrap());

    assert_eq!(K, 128, "do not raise K");
    assert_eq!(would_rt, Some(25), "RT still finds alts at k=25");
    assert!(
        prod.haplotypes.len() > 2,
        "SeqGraph must run; RT 2-hap set is not the assemble result"
    );
    assert!(
        has_a,
        "SeqGraph-only A must survive into assembled haplotypes"
    );
    assert!(event_ga, "EventMap must see 20:29455902 G/A");
}
