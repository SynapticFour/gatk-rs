//! 6R.218: production Java `createGraph` cycle abort on SeqGraph-path RT extract.
//!
//! Frozen Java 4.4.0.0 SHA `2dbc025821bc5f686c423ff332a41e6cef892a77`.
//! ONE production change: `extract_rt_haplotypes_from_built_graph` passes
//! `abort_cyclic_before_dangling = args.use_seq_graph && args.abort_seq_graph_on_cycles`.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r218_java_seqgraph_cycle_gate -- --nocapture --test-threads=1
//! ```

use gatk_core::io::vcf::InfoValue;
use gatk_core::reference::{parse_intervals_cli_string, ReferenceWindowCache, SequenceDictionary};
use gatk_haplotypecaller::assembly::{AssemblyGraphParams, AssemblyGraphPruningParams};
use gatk_haplotypecaller::assembly_region_finalize::{
    assembly_reference_read, create_graph_reference_read, records_to_assembly_reads,
};
use gatk_haplotypecaller::bio_ids::KmerSize;
use gatk_haplotypecaller::read_threading_assembler::{
    assemble_from_ref_and_reads, build_threading_graph_for_haplotype_dump,
    build_threading_graph_for_seq_assembly, extract_haplotypes_from_seq_kbest_paths,
    extract_rt_haplotypes_before_remove_paths,
};
use gatk_haplotypecaller::read_threading_graph::assembly_graph_from_ref_and_reads_threading_with_summary;
use gatk_haplotypecaller::seq_graph::SeqGraph;
use gatk_haplotypecaller::seq_kbest_haplotype::find_best_haplotypes_seq_graph;
use gatk_haplotypecaller::{
    assemble_reads_with_finalized, call_disposition, flatten_assembly_regions,
    traverse_assembly_region_walker, try_emit_call_region_variants, AssemblyRead,
    AssemblyRegionCallDisposition, CallRegionArgs, Cigar, Haplotype, HaplotypeCallerEngine,
    ReadFilterParams, ReadThreadingAssemblerArgs, WalkerTraversalConfig,
    DEFAULT_STAND_EMIT_CONFIDENCE,
};
use std::collections::BTreeSet;
use std::path::{Path, PathBuf};

const JAVA_PIN: &str = "4.4.0.0 / 2dbc025821bc5f686c423ff332a41e6cef892a77";
const INTERVAL: &str = "20:29455000-29456500";
const BAM_REL: &str = "parity/giab/runs/local-pairhmm-diff/HG001.20-29455000-29456500.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const P11_REF_REL: &str = "parity/fixtures/p5_live_reference.fa";
const P11_BAM_REL: &str = "parity/build/sam-indexed-bam/p11_java_positive.bam";
const INDEL4_REF_REL: &str = "parity/fixtures/g2_subset_live_indel4.fa";
const INDEL4_BAM_REL: &str = "parity/build/sam-indexed-bam/g2_subset_live_indel4.bam";
const P12_BAM_REL: &str = "parity/realworld/na12878_20k_b37/NA12878_20k.b37.bam";
const TARGET: u64 = 29_455_015;
const SIB: u64 = 29_455_019;
const JAVA_ACTIVE: &str = "20:29455000-29455145";
const JAVA_UNTRIMMED_N: usize = 78;
const P12_TTC: u64 = 92_307_324;
const TG_INTERVAL: &str = "2:92307200-92307550";

const JAVA_UNTRIMMED_HASHES: &[&str] = &[
    "72369f70c97564fb",
    "9b4130f578a5c3b1",
    "0760164d9ef6dd05",
    "a411db2e5639776a",
    "c3ca8558d183ddd0",
    "ad1727667b42fe9e",
    "c599fb57c27f5a17",
    "31f12d6e2e7365a1",
    "b88de6b98d11d0d7",
    "ed01a60a7f2d5970",
    "b55bcced2fefd423",
    "bfe3c2fb2c2a240d",
    "3990d39c635fde1a",
    "9177506879ad0a9c",
    "810d9fc18079ce30",
    "44d25f40711fe1fa",
    "54afee1dd89e373b",
    "649c345e483e5443",
    "53cb3d3d6311786a",
    "a2f24e4d4a6e8d9e",
    "b9c46fd1377b9041",
    "cdef5e7202584ed0",
    "d1be531684070ff1",
    "eeaca65161357903",
    "5c88c58d266f75d4",
    "cc3f9cadeaeeb4c7",
    "833d4f233780f9a8",
    "9825a1930286fe3e",
    "4afa9304857dc4b1",
    "4685cd1bf0436e72",
    "2cf89addc6b88d5c",
    "8bea11740b6e4204",
    "3bfc7896af17a76f",
    "a14f1b575e96b1fe",
    "e114fb78453ccc5b",
    "3bfbd484c5d0ed45",
    "8aca083d6ff2a5cd",
    "40e4247f4b1a91ec",
    "71e588c6222f0401",
    "8e95b473d340af19",
    "cceb12a3a0c7b29e",
    "20bc439d334fc9ff",
    "de388cc436a18b40",
    "1c39b8c66923c1d2",
    "3fc7de0954270452",
    "c48d5100bcad6748",
    "ae7714c59dafc75e",
    "188da0dd30b021c1",
    "793815145ccc35bb",
    "391503e0ef3650da",
    "dfb7ae3851aed80d",
    "15f1bf3597e17af4",
    "0994a4aaad72cff3",
    "b10c816794fadac2",
    "97cf358727e398f9",
    "284791f3257fae56",
    "93b1d9eda04c9065",
    "985a04a85c946f9a",
    "6ed6c60d101efa23",
    "a7c57619bbaef828",
    "7f563f731e72eaf7",
    "23afacd2bee260a7",
    "edc349af27dd02d0",
    "ac3b46130edb998e",
    "b7ac973e561ef465",
    "8016ead01172b854",
    "9a86c939b608f767",
    "8dea0ce0dc47c9fc",
    "6083fb7b49feebc5",
    "fe10cb506100d9cb",
    "add30f91e861475a",
    "3a3a91f45eb372c8",
    "fdfe44dbde6f93d1",
    "1dc3e207f92421d8",
    "e6a394771b59fa81",
    "2a2b6e4d256f5695",
    "550fe41efa3bd69c",
    "7a4a7948ede22768",
];

const PRE_FIX_EXTRAS: &[&str] = &[
    "70825a603033cd15",
    "20734c56283366c5",
    "e0c1092ca29b8eb6",
    "00bda5d67ed6ce26",
];

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    eprintln!("6R218\t{key}\t{}", value.as_ref());
}

fn fnv1a64_hex(bases: &[u8]) -> String {
    let mut h: u64 = 0xcbf29ce484222325;
    for &b in bases {
        h ^= b as u64;
        h = h.wrapping_mul(0x100000001b3);
    }
    format!("{h:016x}")
}

fn hap_hashes(haps: &[Haplotype]) -> BTreeSet<String> {
    haps.iter().map(|h| fnv1a64_hex(&h.bases)).collect()
}

fn cyclic_before_dangling(
    reference: &AssemblyRead,
    reads: &[AssemblyRead],
    assembler: &ReadThreadingAssemblerArgs,
    kmer: usize,
) -> Option<bool> {
    let kmer_size = KmerSize::try_from_usize(kmer).ok()?;
    let params = AssemblyGraphParams {
        kmer_size,
        min_base_quality: assembler.min_base_quality,
        min_edge_weight: 1,
        dangling_path_max_nodes: 0,
        max_haplotypes: assembler.num_best_haplotypes_per_graph,
        max_haplotype_bases: 4096,
        start_threading_only_at_existing_vertex: !assembler.recover_dangling_branches,
    };
    let (mut graph, summary) =
        assembly_graph_from_ref_and_reads_threading_with_summary(reference, reads, &params).ok()?;
    if !assembler.allow_low_complexity_graphs && summary.is_low_complexity {
        return None;
    }
    let mut pruning = AssemblyGraphPruningParams::gatk_haplotype_caller_defaults();
    pruning.min_prune_factor = assembler.min_prune_factor;
    pruning.use_adaptive_pruning = assembler.use_adaptive_pruning;
    if assembler.prune_before_cycle_counting {
        graph.apply_pruning(&pruning);
    }
    Some(graph.has_cycle())
}

fn seq_kbest_haps(
    graph_ref: &AssemblyRead,
    graph_reads: &[AssemblyRead],
    assembler: &ReadThreadingAssemblerArgs,
    kmer: usize,
) -> Vec<Haplotype> {
    let Some(graph) = build_threading_graph_for_seq_assembly(
        graph_ref,
        graph_reads,
        kmer,
        assembler,
        false,
        false,
    )
    .ok()
    .flatten() else {
        return Vec::new();
    };
    let mut seq = SeqGraph::from_assembly_graph(&graph);
    seq.clean_non_ref_paths();
    let _ = seq.cleanup_seq_graph();
    let Ok(paths) = find_best_haplotypes_seq_graph(&seq, assembler.num_best_haplotypes_per_graph)
    else {
        return Vec::new();
    };
    let mut ref_hap = Haplotype::new(graph_ref.bases.as_slice(), true);
    let mut ref_cigar = Cigar::new();
    ref_cigar.push(
        ref_hap.bases.len(),
        gatk_haplotypecaller::CigarOperator::Match,
    );
    ref_hap.cigar = Some(ref_cigar);
    let ref_cigar_len = ref_hap.cigar.as_ref().unwrap().reference_length();
    extract_haplotypes_from_seq_kbest_paths(
        &paths,
        &seq,
        kmer,
        &ref_hap,
        ref_cigar_len,
        &assembler.haplotype_to_reference_sw,
    )
    .unwrap_or_default()
}

fn info_f64(info: &[InfoValue], key: &str) -> Option<f64> {
    for v in info {
        if let InfoValue::Float(k, xs) = v {
            if k == key {
                return xs.first().copied();
            }
        }
    }
    None
}

#[test]
fn forensic_6r218_production_gate_pin() {
    kv("java_pin", JAVA_PIN);
    kv(
        "production_rule",
        "extract_rt_haplotypes_from_built_graph: abort_cyclic_before_dangling = use_seq_graph && abort_seq_graph_on_cycles",
    );
    kv("classification", "GRAPH_STATE_DIVERGENCE");
    kv(
        "production_change",
        "ONE — SeqGraph-path RT extract cycle abort",
    );
    let args = ReadThreadingAssemblerArgs::default();
    assert!(args.use_seq_graph);
    assert!(args.abort_seq_graph_on_cycles);
    assert_eq!(JAVA_UNTRIMMED_HASHES.len(), JAVA_UNTRIMMED_N);
}

#[test]
fn forensic_6r218_live_java_seqgraph_cycle_gate() {
    let root = repo_root();
    let ref_fasta = root.join(REF_REL);
    let bam = root.join(BAM_REL);
    if !ref_fasta.is_file() || !bam.is_file() {
        eprintln!("skip: missing chr20_tiny BAM/REF");
        return;
    }
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
            ) && r.start.get() <= TARGET
                && r.end.get() >= TARGET
        })
        .expect("covering");
    assert_eq!(
        format!(
            "{}:{}-{}",
            region.contig,
            region.start.get(),
            region.end.get()
        ),
        JAVA_ACTIVE
    );

    let mut owned = region.clone();
    let mut ref_cache = ReferenceWindowCache::new(ref_fasta.clone(), 4);
    let args = CallRegionArgs::strict_java();
    let assembled =
        assemble_reads_with_finalized(&mut owned, &dict, &mut ref_cache, &args.assemble)
            .expect("assemble");
    let padded_ref = assembly_reference_read(&dict, &mut ref_cache, region).expect("pad");
    let graph_ref = create_graph_reference_read(&padded_ref, region, &dict);
    let graph_reads = records_to_assembly_reads(&assembled.finalized_reads);
    let assembler = args.assemble.assembler.clone();
    assert!(assembler.use_seq_graph && assembler.abort_seq_graph_on_cycles);

    let cyclic10 = cyclic_before_dangling(&graph_ref, &graph_reads, &assembler, 10);
    let seq_g10 = build_threading_graph_for_seq_assembly(
        &graph_ref,
        &graph_reads,
        10,
        &assembler,
        false,
        false,
    )
    .expect("seq k10");
    let dump_g10 = build_threading_graph_for_haplotype_dump(
        &graph_ref,
        &graph_reads,
        10,
        &assembler,
        false,
        false,
    )
    .expect("dump k10");
    kv("k10_cyclic_before_dangling", format!("{cyclic10:?}"));
    kv(
        "java_createGraph_k10",
        "null — generateSeqGraph && hasCycles before dangling",
    );
    kv(
        "rust_seq_createGraph_k10",
        if seq_g10.is_none() { "None" } else { "Some" },
    );
    kv(
        "dump_graph_k10",
        if dump_g10.is_some() {
            "Some — dump/fallback still recovers dangling"
        } else {
            "None"
        },
    );
    assert_eq!(cyclic10, Some(true));
    assert!(seq_g10.is_none());
    assert!(dump_g10.is_some());

    let seq25 = seq_kbest_haps(&graph_ref, &graph_reads, &assembler, 25);
    let seq10 = seq_kbest_haps(&graph_ref, &graph_reads, &assembler, 10);
    let java_untrimmed: BTreeSet<&str> = JAVA_UNTRIMMED_HASHES.iter().copied().collect();
    let seq25_h = hap_hashes(&seq25);
    let seq25_refs: BTreeSet<&str> = seq25_h.iter().map(String::as_str).collect();
    kv("seq_kbest_k10_n", seq10.len().to_string());
    kv("seq_kbest_k25_n", seq25.len().to_string());
    assert_eq!(seq10.len(), 0);
    assert_eq!(seq25.len(), JAVA_UNTRIMMED_N);
    assert_eq!(java_untrimmed.difference(&seq25_refs).count(), 0);
    assert_eq!(seq25_refs.difference(&java_untrimmed).count(), 0);

    let rt10_seqgraph = extract_rt_haplotypes_before_remove_paths(
        &graph_ref,
        &graph_reads,
        &assembler,
        10,
        false,
        false,
    )
    .unwrap_or_default();
    let mut fallback_args = assembler.clone();
    fallback_args.use_seq_graph = false;
    let rt10_fallback = extract_rt_haplotypes_before_remove_paths(
        &graph_ref,
        &graph_reads,
        &fallback_args,
        10,
        false,
        false,
    )
    .unwrap_or_default();
    kv("rt_k10_seqgraph_extract_n", rt10_seqgraph.len().to_string());
    kv("rt_k10_fallback_extract_n", rt10_fallback.len().to_string());
    kv(
        "gate_evidence",
        "seqgraph extract n=0 ⇒ createGraph None ⇒ dangling/k-best not entered; fallback n>0 ⇒ use_seq_graph=false unchanged",
    );
    assert!(
        rt10_seqgraph.is_empty(),
        "SeqGraph-path RT extract must abort cyclic k=10 before dangling/k-best"
    );
    assert!(
        !rt10_fallback.is_empty(),
        "use_seq_graph=false must still extract the cyclic graph"
    );
    assert_eq!(rt10_fallback.len(), 54);
    let fallback_h = hap_hashes(&rt10_fallback);
    for pin in PRE_FIX_EXTRAS {
        assert!(
            fallback_h.contains(*pin),
            "pre-fix extra {pin} still exists on use_seq_graph=false extract; SeqGraph-path abort did not hash-delete it"
        );
        assert!(
            !rt10_seqgraph.iter().any(|h| fnv1a64_hex(&h.bases) == *pin),
            "pre-fix extra {pin} must not be extracted on the SeqGraph path"
        );
    }

    let prod = assemble_from_ref_and_reads(&graph_ref, &graph_reads, &assembler).expect("prod");
    let prod_h = hap_hashes(&prod.haplotypes);
    let prod_refs: BTreeSet<&str> = prod_h.iter().map(String::as_str).collect();
    let rust_only: Vec<&str> = prod_refs.difference(&java_untrimmed).copied().collect();
    kv("assemble_n", prod.haplotypes.len().to_string());
    kv("rust_only_n", rust_only.len().to_string());
    kv(
        "java_only_n",
        java_untrimmed.difference(&prod_refs).count().to_string(),
    );
    assert_eq!(prod.haplotypes.len(), JAVA_UNTRIMMED_N);
    assert!(
        rust_only.is_empty(),
        "44 extras must vanish by graph reject, not hash delete"
    );
    assert_eq!(java_untrimmed.difference(&prod_refs).count(), 0);
    for pin in PRE_FIX_EXTRAS {
        assert!(
            !prod_h.contains(*pin),
            "pre-fix extra {pin} must be absent because cyclic graph was rejected"
        );
    }

    let outcome = HaplotypeCallerEngine::call_region(
        region,
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("call")
    .expect("outcome");
    kv(
        "call_region_hap_n",
        outcome.assembly.haplotypes.len().to_string(),
    );
    let tgt = outcome
        .genotyped_calls
        .iter()
        .find(|c| {
            c.event.start_1based.get() == TARGET
                && c.event.ref_allele == "G"
                && c.event.alt_allele == "T"
        })
        .expect("G/T");
    let pl = tgt.genotype.format.pl_as_i32();
    let gq = tgt.genotype.format.gq;
    kv("rust_pl", format!("{pl:?}"));
    kv("java_pl", "69,0,2140");
    kv("rust_gq", format!("{gq:?}"));
    kv("java_gq", "69");
    let emitted =
        try_emit_call_region_variants(region, &outcome, "SAMPLE", DEFAULT_STAND_EMIT_CONFIDENCE)
            .unwrap_or_default();
    let rec = emitted
        .iter()
        .find(|r| r.position == TARGET && r.reference == "G")
        .expect("emit G/T");
    let sib = emitted
        .iter()
        .find(|r| r.position == SIB && r.reference == "G")
        .expect("sib");
    kv("rust_qual", format!("{:?}", rec.quality));
    kv("java_qual", "61.64");
    kv("rust_qd", format!("{:?}", info_f64(&rec.info, "QD")));
    kv("java_qd", "0.96");
    kv("sib_pl", format!("{:?}", sib.samples[0].pl.as_deref()));
    kv("java_sib_pl", "645,0,1147");
    kv("target_pl", "MATCH Java 69,0,2140");
    assert_eq!(pl, vec![69, 0, 2140]);
    assert_eq!(tgt.genotype.format.gq.as_i32(), 69);
    assert_eq!(outcome.assembly.haplotypes.len(), 30);
    assert_eq!(rec.samples[0].pl.as_deref(), Some([69, 0, 2140].as_slice()));
    assert!((rec.quality.expect("QUAL") - 61.64).abs() < 0.01);
    let qd = info_f64(&rec.info, "QD").expect("QD");
    assert!((qd - 0.96).abs() < 0.02);
    assert_eq!(
        sib.samples[0].pl.as_deref(),
        Some([645, 0, 1147].as_slice())
    );
    kv("kbest_policy", "legacy_1024 / K=128 unchanged");
}

#[test]
fn forensic_6r218_l2_p11_acyclic_alt_survives() {
    let root = repo_root();
    let ref_fasta = root.join(P11_REF_REL);
    let bam = root.join(P11_BAM_REL);
    if !ref_fasta.is_file() || !bam.is_file() {
        eprintln!("skip: missing p11");
        return;
    }
    let dict = SequenceDictionary::from_fasta_path(&ref_fasta).expect("dict");
    let specs = parse_intervals_cli_string(&dict, "chrLive:1-63").expect("interval");
    let walk = traverse_assembly_region_walker(
        &dict,
        &specs,
        &ref_fasta,
        &bam,
        &ReadFilterParams::gatk_standard_hc(),
        &WalkerTraversalConfig::gatk_haplotype_caller_production(0),
    )
    .expect("walk");
    let region = flatten_assembly_regions(&walk)
        .into_iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            )
        })
        .expect("p11");
    let mut owned = region.clone();
    let mut ref_cache = ReferenceWindowCache::new(ref_fasta.clone(), 4);
    let args = CallRegionArgs::strict_java();
    let assembled =
        assemble_reads_with_finalized(&mut owned, &dict, &mut ref_cache, &args.assemble)
            .expect("assemble");
    let padded_ref = assembly_reference_read(&dict, &mut ref_cache, &region).expect("pad");
    let graph_ref = create_graph_reference_read(&padded_ref, &region, &dict);
    let graph_reads = records_to_assembly_reads(&assembled.finalized_reads);
    let assembler = args.assemble.assembler.clone();
    let cyclic10 = cyclic_before_dangling(&graph_ref, &graph_reads, &assembler, 10);
    let seq10 = seq_kbest_haps(&graph_ref, &graph_reads, &assembler, 10);
    let rt10 = extract_rt_haplotypes_before_remove_paths(
        &graph_ref,
        &graph_reads,
        &assembler,
        10,
        false,
        false,
    )
    .unwrap_or_default();
    let prod = assemble_from_ref_and_reads(&graph_ref, &graph_reads, &assembler).expect("prod");
    kv("p11_k10_cyclic", format!("{cyclic10:?}"));
    kv("p11_seq10_n", seq10.len().to_string());
    kv("p11_rt10_n", rt10.len().to_string());
    kv("p11_assemble_n", prod.haplotypes.len().to_string());
    assert_eq!(cyclic10, Some(false));
    assert_eq!(seq10.len(), 2);
    assert!(seq10.iter().any(|h| !h.is_reference));
    assert!(rt10.iter().any(|h| !h.is_reference));
    assert!(prod.haplotypes.iter().any(|h| !h.is_reference));
}

#[test]
fn forensic_6r218_l2_indel4_and_p12_ttc() {
    let root = repo_root();
    let args = CallRegionArgs::strict_java();
    let indel4_ref = root.join(INDEL4_REF_REL);
    let indel4_bam = root.join(INDEL4_BAM_REL);
    if indel4_ref.is_file() && indel4_bam.is_file() {
        let dict = SequenceDictionary::from_fasta_path(&indel4_ref).expect("dict");
        let specs = parse_intervals_cli_string(&dict, "chrIndel4:1-68").expect("indel4");
        let walk = traverse_assembly_region_walker(
            &dict,
            &specs,
            &indel4_ref,
            &indel4_bam,
            &ReadFilterParams::gatk_standard_hc(),
            &WalkerTraversalConfig::gatk_haplotype_caller_production(0),
        )
        .expect("walk");
        let region = flatten_assembly_regions(&walk)
            .into_iter()
            .find(|r| {
                matches!(
                    call_disposition(r),
                    AssemblyRegionCallDisposition::ActiveFull
                )
            })
            .expect("indel4");
        let mut owned = region.clone();
        let mut ref_cache = ReferenceWindowCache::new(indel4_ref.clone(), 4);
        let assembled =
            assemble_reads_with_finalized(&mut owned, &dict, &mut ref_cache, &args.assemble)
                .expect("assemble");
        let padded_ref = assembly_reference_read(&dict, &mut ref_cache, &region).expect("pad");
        let graph_ref = create_graph_reference_read(&padded_ref, &region, &dict);
        let graph_reads = records_to_assembly_reads(&assembled.finalized_reads);
        let assembler = args.assemble.assembler.clone();
        let cyclic10 = cyclic_before_dangling(&graph_ref, &graph_reads, &assembler, 10);
        let seq10 = seq_kbest_haps(&graph_ref, &graph_reads, &assembler, 10);
        let prod = assemble_from_ref_and_reads(&graph_ref, &graph_reads, &assembler).expect("prod");
        kv("indel4_k10_cyclic", format!("{cyclic10:?}"));
        kv("indel4_seq10_n", seq10.len().to_string());
        kv("indel4_assemble_n", prod.haplotypes.len().to_string());
        assert_eq!(cyclic10, Some(false));
        assert_eq!(seq10.len(), 3);
        assert_eq!(prod.haplotypes.len(), 3);
    } else {
        kv("indel4", "skip missing BAM");
    }

    let ref_fasta = root.join(REF_REL);
    let p12_bam = root.join(P12_BAM_REL);
    if !ref_fasta.is_file() || !p12_bam.is_file() {
        eprintln!("skip: missing P12");
        return;
    }
    let dict = SequenceDictionary::from_fasta_path(&ref_fasta).expect("dict");
    let specs = parse_intervals_cli_string(&dict, TG_INTERVAL).expect("tg");
    let walk = traverse_assembly_region_walker(
        &dict,
        &specs,
        &ref_fasta,
        &p12_bam,
        &ReadFilterParams::gatk_standard_hc(),
        &WalkerTraversalConfig::gatk_haplotype_caller_production(100),
    )
    .expect("walk");
    let region = flatten_assembly_regions(&walk)
        .into_iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= P12_TTC
                && r.end.get() >= P12_TTC
        })
        .expect("ttc");
    let mut owned = region.clone();
    let mut ref_cache = ReferenceWindowCache::new(ref_fasta.clone(), 4);
    let assembled =
        assemble_reads_with_finalized(&mut owned, &dict, &mut ref_cache, &args.assemble)
            .expect("assemble");
    let padded_ref = assembly_reference_read(&dict, &mut ref_cache, &region).expect("pad");
    let graph_ref = create_graph_reference_read(&padded_ref, &region, &dict);
    let graph_reads = records_to_assembly_reads(&assembled.finalized_reads);
    let assembler = args.assemble.assembler.clone();
    let cyclic10 = cyclic_before_dangling(&graph_ref, &graph_reads, &assembler, 10);
    let seq25 = seq_kbest_haps(&graph_ref, &graph_reads, &assembler, 25);
    let prod = assemble_from_ref_and_reads(&graph_ref, &graph_reads, &assembler).expect("prod");
    kv("p12_ttc_k10_cyclic", format!("{cyclic10:?}"));
    kv("p12_ttc_seq25_n", seq25.len().to_string());
    kv("p12_ttc_assemble_n", prod.haplotypes.len().to_string());
    assert_eq!(cyclic10, Some(true));
    assert_eq!(seq25.len(), 6);
    assert_eq!(prod.haplotypes.len(), 6);
}
