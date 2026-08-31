//! 6R.55 forensic: Java-vs-Rust control-flow contract around RT-configured
//! alternatives vs SeqGraph. No production algorithm change.
//!
//! Proven locus (evidence only): ActiveFull `20:29455745–29455993`,
//! `20:29455902 G/A`.
//!
//! Skipped unless `HOLDOUT_6R55=1`.
//!
//! ```text
//! HOLDOUT_6R55=1 cargo test -p gatk-haplotypecaller --test holdout_6r55_chr20_tiny_test -- --nocapture
//! ```

use gatk_core::reference::{parse_intervals_cli_string, ReferenceWindowCache, SequenceDictionary};
use gatk_haplotypecaller::assembly_region_finalize::{
    assembly_reference_read, create_graph_reference_read, records_to_assembly_reads,
};
use gatk_haplotypecaller::hc_allele_mapping::hap_base_at_ref_locus;
use gatk_haplotypecaller::hc_genotyping_engine::DEFAULT_STAND_EMIT_CONFIDENCE;
use gatk_haplotypecaller::read_threading_assembler::{
    region_overlaps_p12_cluster, AssemblyScoringContext, DEFAULT_NUM_BEST_HAPLOTYPES_PER_GRAPH,
};
use gatk_haplotypecaller::{
    assemble_from_ref_and_reads, assemble_reads_with_finalized, call_disposition,
    diagnostic_rt_first_skip_seq_graph_kmer, flatten_assembly_regions,
    traverse_assembly_region_walker, try_emit_call_region_variants, AssemblyRegion,
    AssemblyRegionCallDisposition, CallRegionArgs, Haplotype, HaplotypeCallerEngine,
    ReadFilterParams, WalkerTraversalConfig,
};
use serde_json::{json, Value};
use std::path::{Path, PathBuf};

const INTERVAL: &str = "20:29455000-29456500";
const BAM_REL: &str = "parity/giab/runs/local-pairhmm-diff/HG001.20-29455000-29456500.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const ACTIVE: (u64, u64) = (29_455_745, 29_455_993);
const TARGET: u64 = 29_455_902;
const TARGET_ALT: u8 = b'A';
const K: usize = DEFAULT_NUM_BEST_HAPLOTYPES_PER_GRAPH;
const L_GATE_START: u64 = 92_300_000;
const L_GATE_END: u64 = 92_350_000;

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn cigar_str(h: &Haplotype) -> String {
    h.cigar
        .as_ref()
        .map(|c| c.to_gatk_string())
        .unwrap_or_else(|| "-".to_string())
}

fn hap_has_alt_a(haps: &[Haplotype], pad: u64) -> bool {
    haps.iter()
        .any(|h| hap_base_at_ref_locus(h, pad, TARGET) == Some(TARGET_ALT))
}

fn event_has_ga(assembly: &gatk_haplotypecaller::AssemblyResultSet) -> bool {
    assembly
        .variation_events()
        .iter()
        .any(|e| e.start_1based.get() == TARGET && e.ref_allele == "G" && e.alt_allele == "A")
}

fn production_assembler(
    region: &AssemblyRegion,
    args: &CallRegionArgs,
) -> gatk_haplotypecaller::ReadThreadingAssemblerArgs {
    let mut assembler = args.assemble.assembler.clone();
    assembler.dangling_java_exact = true;
    assembler.scoring = Some(AssemblyScoringContext {
        padded_reference_start_1based: region.extended_start.get(),
        active_start_1based: region.start.get(),
        active_end_1based: region.end.get(),
        contig: region.contig.clone(),
    });
    if region_overlaps_p12_cluster(region.start.get(), region.end.get()) {
        assembler.use_seq_graph = false;
        assembler.remove_paths_not_connected_to_ref = false;
        assembler.skip_post_dangling_prune = true;
    }
    assembler
}

/// A = SeqGraph reached (Java-equivalent for this helper).
/// B = Rust-only RT-first skip.
/// C = RT-first not involved; another Rust SeqGraph bypass (P12 `use_seq_graph=false`).
fn classify_control_flow(region: &AssemblyRegion, skip_k: Option<usize>) -> &'static str {
    if region_overlaps_p12_cluster(region.start.get(), region.end.get()) {
        "C"
    } else if skip_k.is_some() {
        "B"
    } else {
        "A"
    }
}

fn graph_inputs(
    region: &AssemblyRegion,
    dict: &SequenceDictionary,
    ref_fasta: &Path,
    args: &CallRegionArgs,
) -> (
    gatk_haplotypecaller::AssemblyRead,
    Vec<gatk_haplotypecaller::AssemblyRead>,
    u64,
) {
    let mut ref_cache = ReferenceWindowCache::new(ref_fasta.to_path_buf(), 4);
    let mut owned = region.clone();
    let assembled = assemble_reads_with_finalized(&mut owned, dict, &mut ref_cache, &args.assemble)
        .expect("assemble");
    let pad = assembled.assembly.padded_reference_start_1based();
    let padded_ref = assembly_reference_read(dict, &mut ref_cache, region).expect("pad ref");
    let graph_ref = create_graph_reference_read(&padded_ref, region, dict);
    let graph_reads = records_to_assembly_reads(&assembled.finalized_reads);
    (graph_ref, graph_reads, pad)
}

fn classify_region(
    region: &AssemblyRegion,
    dict: &SequenceDictionary,
    ref_fasta: &Path,
    args: &CallRegionArgs,
) -> Value {
    let (graph_ref, graph_reads, _) = graph_inputs(region, dict, ref_fasta, args);
    let assembler = production_assembler(region, args);
    let skip_k = diagnostic_rt_first_skip_seq_graph_kmer(&graph_ref, &graph_reads, &assembler)
        .expect("rt_first diagnostic");
    let contig = region.contig.as_str();
    let l_gate = (contig == "2" || contig == "chr2")
        && region.end.get() >= L_GATE_START
        && region.start.get() <= L_GATE_END;
    json!({
        "contig": region.contig,
        "active": [region.start.get(), region.end.get()],
        "p12_cluster": region_overlaps_p12_cluster(region.start.get(), region.end.get()),
        "p12_l_gate": l_gate,
        "use_seq_graph": assembler.use_seq_graph,
        "rt_first_skip_k": skip_k,
        "class": classify_control_flow(region, skip_k),
    })
}

fn walk_active_fulls(
    dict: &SequenceDictionary,
    ref_fasta: &Path,
    bam: &Path,
    interval: &str,
) -> Vec<AssemblyRegion> {
    let specs = parse_intervals_cli_string(dict, interval).expect("interval");
    let walk = traverse_assembly_region_walker(
        dict,
        &specs,
        ref_fasta,
        bam,
        &ReadFilterParams::gatk_standard_hc(),
        &WalkerTraversalConfig::gatk_haplotype_caller_production(100),
    )
    .expect("walk");
    flatten_assembly_regions(&walk)
        .into_iter()
        .filter(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            )
        })
        .collect()
}

#[test]
fn holdout_6r55_rt_first_control_flow_contract() {
    if std::env::var("HOLDOUT_6R55").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R55=1");
        return;
    }
    let root = repo_root();
    let ref_fasta = root.join(REF_REL);
    let bam = root.join(BAM_REL);
    assert!(ref_fasta.is_file(), "missing {}", ref_fasta.display());
    assert!(bam.is_file(), "missing {}", bam.display());

    let dict = SequenceDictionary::from_fasta_path(&ref_fasta).expect("dict");
    let args = CallRegionArgs::strict_java();
    let regions = walk_active_fulls(&dict, &ref_fasta, &bam, INTERVAL);
    let region = regions
        .iter()
        .find(|r| r.start.get() == ACTIVE.0 && r.end.get() == ACTIVE.1)
        .expect("ActiveFull 29455745-29455993");

    let (graph_ref, graph_reads, pad) = graph_inputs(region, &dict, &ref_fasta, &args);
    let assembler = production_assembler(region, &args);
    let skip_k = diagnostic_rt_first_skip_seq_graph_kmer(&graph_ref, &graph_reads, &assembler)
        .expect("rt_first diagnostic");

    let prod = assemble_from_ref_and_reads(&graph_ref, &graph_reads, &assembler).expect("prod");
    let mut seq_args = assembler.clone();
    seq_args.scoring = None;
    let seq_only =
        assemble_from_ref_and_reads(&graph_ref, &graph_reads, &seq_args).expect("seq-only");

    let mut counterfactual_assemble = args.assemble.clone();
    counterfactual_assemble.strict_java_assembly = false;
    counterfactual_assemble.assembler.dangling_java_exact = true;
    let mut owned = region.clone();
    let mut ref_cache = ReferenceWindowCache::new(ref_fasta.clone(), 4);
    let cf =
        assemble_reads_with_finalized(&mut owned, &dict, &mut ref_cache, &counterfactual_assemble)
            .expect("counterfactual assemble");
    let cf_pad = cf.assembly.padded_reference_start_1based();

    let call = HaplotypeCallerEngine::call_region(region, &dict, &ref_fasta, &args).expect("call");
    let outcome = call.expect("production call_region Some");
    let vcf =
        try_emit_call_region_variants(region, &outcome, "NA12878", DEFAULT_STAND_EMIT_CONFIDENCE)
            .expect("emit");
    let vcf_has_ga = vcf.iter().any(|r| {
        r.position == TARGET && r.reference == "G" && r.alternate.iter().any(|a| a == "A")
    });

    let locus = json!({
        "active": [region.start.get(), region.end.get()],
        "k_production": K,
        "production": {
            "rt_configured_alts_would_succeed_k": skip_k,
            "seq_graph_not_replaced_by_rt": prod.haplotypes.len() > 2,
            "n_haps": prod.haplotypes.len(),
            "cigars": prod.haplotypes.iter().take(8).map(cigar_str).collect::<Vec<_>>(),
            "has_alt_a": hap_has_alt_a(&prod.haplotypes, pad),
            "vcf_has_20_29455902_G_A": vcf_has_ga,
        },
        "counterfactual_scoring_off": {
            "n_haps": seq_only.haplotypes.len(),
            "has_alt_a": hap_has_alt_a(&seq_only.haplotypes, pad),
            "assemble_reads_event_has_G_A": event_has_ga(&cf.assembly),
            "assemble_reads_has_alt_a": hap_has_alt_a(&cf.assembly.haplotypes, cf_pad),
            "note": "scoring disabled only on AssembleReadsArgs in this test; no production flag",
        },
    });

    let chr20_tiny: Vec<Value> = regions
        .iter()
        .map(|r| classify_region(r, &dict, &ref_fasta, &args))
        .collect();

    let mut other = Vec::new();
    let mid_b_bam = root.join("parity/realworld/na12878_20k_b37/NA12878_20k.b37.bam");
    if mid_b_bam.is_file() {
        let mid = walk_active_fulls(&dict, &ref_fasta, &mid_b_bam, "2:92317250-92317520");
        for r in mid.iter().take(4) {
            let mut row = classify_region(r, &dict, &ref_fasta, &args);
            row["id"] = json!("ctrl_mid_b");
            other.push(row);
        }
    }
    for (id, interval, bam_rel, cap) in [
        (
            "chr20_w47",
            "20:47000000-47050000",
            "parity/giab/runs/hang-repro/w47_50k.bam",
            4usize,
        ),
        (
            "chr21_w10",
            "21:10000000-10035000",
            "parity/giab/runs/hang-repro/w10_50k.bam",
            4usize,
        ),
    ] {
        let other_bam = root.join(bam_rel);
        if !other_bam.is_file() {
            other.push(json!({"id": id, "class": "D", "reason": "bam missing"}));
            continue;
        }
        let regs = walk_active_fulls(&dict, &ref_fasta, &other_bam, interval);
        for r in regs.iter().take(cap) {
            let mut row = classify_region(r, &dict, &ref_fasta, &args);
            row["id"] = json!(id);
            other.push(row);
        }
    }

    let doc = json!({
        "locus_20_29455902": locus,
        "generalization_chr20_tiny": chr20_tiny,
        "generalization_other": other,
    });
    println!("{}", serde_json::to_string_pretty(&doc).unwrap());

    assert_eq!(K, 128, "do not raise K");
    assert_eq!(
        skip_k,
        Some(25),
        "RT configured extract still succeeds at k=25 (diagnostic)"
    );
    assert!(
        hap_has_alt_a(&prod.haplotypes, pad),
        "production SeqGraph path must retain A (6R.56)"
    );
    assert!(
        prod.haplotypes.len() > 2,
        "production must not return the 2-hap RT shortcut set"
    );
    assert!(
        hap_has_alt_a(&seq_only.haplotypes, pad),
        "SeqGraph path must retain A"
    );
    assert!(
        event_has_ga(&cf.assembly),
        "scoring-off EventMap must include G/A"
    );
    // VCF emit is reported; do not require QUAL/AD parity in this diagnostic.
    eprintln!("production VCF has 20:29455902 G/A: {vcf_has_ga}");
}
