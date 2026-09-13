//! 6R.162: trace extra after_assemble hap `3a53b2a941bbdc43` to the first assembly split.
//! Production k-best `legacy_1024`. Skipped unless `HOLDOUT_6R162=1`.
//!
//! Production change: NONE. Does not chase PairHMM / AD / `45,3,0`.
//!
//! ```text
//! HOLDOUT_6R162=1 cargo test -p gatk-haplotypecaller --test holdout_6r162_extra_deletion_haplotype -- --nocapture --test-threads=1
//! ```

use gatk_core::reference::{parse_intervals_cli_string, ReferenceWindowCache, SequenceDictionary};
use gatk_haplotypecaller::assembly_based_caller::{
    assemble_reads_with_finalized, AssembleReadsArgs,
};
use gatk_haplotypecaller::assembly_region_finalize::{
    create_graph_reference_read, records_to_assembly_reads,
};
use gatk_haplotypecaller::read_threading_assembler::{
    audit_threading_dangling_recovery, build_threading_graph_for_seq_assembly,
    extract_rt_haplotypes_after_remove_paths, extract_rt_haplotypes_before_remove_paths,
    probe_seq_graph_kmer_attempts, ReadThreadingAssemblerArgs,
};
use gatk_haplotypecaller::seq_graph::SeqGraph;
use gatk_haplotypecaller::seq_kbest_haplotype::find_best_haplotypes_seq_graph;
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, traverse_assembly_region_walker,
    AssemblyRegionCallDisposition, ReadFilterParams, WalkerTraversalConfig,
};

use std::path::{Path, PathBuf};

const INTERVAL: &str = "2:92316200-92316580";
const BAM_REL: &str = "parity/realworld/na12878_20k_b37/NA12878_20k.b37.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const TARGET: u64 = 92_316_347;
const EXTRA_FNV: &str = "3a53b2a941bbdc43";
const JAVA_REF: &str = "9b4092f3b20a5de7";
const JAVA_ALT: &str = "8c9a6fc8f302f4a3";
/// Frozen dangling path_bases (6R.163). Used only to prove the sequence is absent from k-best.
const EXTRA_SEQ: &str = "GTCAGCAGCATTCTCAGAAAGTTCTTTGTGATGATTGCATTCAAGTCACAGAATTGAACATTCCCTTTCACAGAGCAGGTTTGAAACACTCTTTTTGTAGTGTCTGTAAGTGAACATTTGGATTGCTTTCAGGCCTAAGGTGAAAAAGGAAATATCTTCCCATAAAAACTAGAC";

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    println!("6R162\t{key}\t{}", value.as_ref());
}

fn fnv1a64_hex(bases: &[u8]) -> String {
    let mut h: u64 = 0xcbf29ce484222325;
    for &b in bases {
        h ^= u64::from(b);
        h = h.wrapping_mul(0x100000001b3);
    }
    format!("{h:016x}")
}

#[test]
fn holdout_6r162_extra_deletion_haplotype() {
    if std::env::var("HOLDOUT_6R162").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R162=1");
        return;
    }
    let root = repo_root();
    let ref_fasta = root.join(REF_REL);
    let bam = root.join(BAM_REL);
    if !ref_fasta.is_file() || !bam.is_file() {
        eprintln!("skip: missing P12 ref/BAM");
        return;
    }

    kv(
        "java_pin",
        "4.4.0.0 / 2dbc025821bc5f686c423ff332a41e6cef892a77",
    );
    kv("variant", format!("2:{TARGET} G/A"));
    kv("extra_fnv", EXTRA_FNV);
    kv(
        "kbest_policy",
        "legacy_1024 (unset GATK_RS_EXPERIMENTAL_KBEST_POLICY)",
    );

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
    let covering = regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= TARGET
                && r.end.get() >= TARGET
        })
        .expect("ActiveFull covering 2:92316347");
    kv(
        "active_full",
        format!(
            "{}:{}-{}\textended={}-{}\tn_reads={}",
            covering.contig,
            covering.start.get(),
            covering.end.get(),
            covering.extended_start.get(),
            covering.extended_end.get(),
            covering.reads.len()
        ),
    );

    let mut region = covering.clone();
    let mut ref_cache = ReferenceWindowCache::new(ref_fasta.clone(), 4);
    let mut args = AssembleReadsArgs::default();
    args.strict_java_assembly = true;
    let assembled =
        assemble_reads_with_finalized(&mut region, &dict, &mut ref_cache, &args).expect("assemble");
    kv(
        "assemble_n",
        assembled.assembly.haplotypes.len().to_string(),
    );
    kv(
        "finalized_reads",
        assembled.finalized_reads.len().to_string(),
    );

    let mut extra_bases: Option<Vec<u8>> = None;
    for (i, h) in assembled.assembly.haplotypes.iter().enumerate() {
        let hash = fnv1a64_hex(&h.bases);
        let cigar = h
            .cigar
            .as_ref()
            .map(|c| c.to_gatk_string())
            .unwrap_or_else(|| ".".into());
        kv(
            "assemble_hap",
            format!(
                "idx={i}\thash={hash}\tisRef={}\tlen={}\tscore={}\tkmer={}\talign0={}\tcigar={cigar}",
                h.is_reference,
                h.bases.len(),
                h.score,
                h.kmer_size,
                h.alignment_start_hap_wrt_ref
            ),
        );
        if hash == EXTRA_FNV {
            extra_bases = Some(h.bases.clone());
            kv("extra_seq", String::from_utf8_lossy(&h.bases).into_owned());
        }
    }
    let hashes: Vec<String> = assembled
        .assembly
        .haplotypes
        .iter()
        .map(|h| fnv1a64_hex(&h.bases))
        .collect();
    kv(
        "shared_java_ref_present",
        hashes.iter().any(|h| h == JAVA_REF).to_string(),
    );
    kv(
        "shared_java_alt_present",
        hashes.iter().any(|h| h == JAVA_ALT).to_string(),
    );
    kv(
        "extra_present",
        hashes.iter().any(|h| h == EXTRA_FNV).to_string(),
    );
    assert_eq!(
        assembled.assembly.haplotypes.len(),
        2,
        "6R.164: Java-exact assemble n=2"
    );
    assert!(
        hashes.iter().any(|h| h == JAVA_REF) && hashes.iter().any(|h| h == JAVA_ALT),
        "shared Java REF/ALT haplotypes must remain"
    );
    assert!(
        extra_bases.is_none(),
        "6R.164: dangling path_bases must not be a standalone haplotype"
    );
    let extra_bases = EXTRA_SEQ.as_bytes().to_vec();

    for rec in &assembled.finalized_reads {
        let qname = String::from_utf8_lossy(rec.qname()).into_owned();
        kv(
            "finalize_read",
            format!(
                "qname={qname}\tflags={}\tstart={}\tcigar={}\tlen={}",
                rec.flags(),
                rec.pos() + 1,
                rec.cigar().to_string(),
                rec.seq_len()
            ),
        );
    }

    let graph_ref = create_graph_reference_read(
        &{
            gatk_haplotypecaller::assembly_region_finalize::assembly_reference_read(
                &dict,
                &mut ReferenceWindowCache::new(ref_fasta.clone(), 4),
                covering,
            )
            .expect("ref")
        },
        covering,
        &dict,
    );
    kv(
        "graph_ref",
        format!(
            "len={}\thash={}",
            graph_ref.bases.len(),
            fnv1a64_hex(&graph_ref.bases)
        ),
    );
    let reads = records_to_assembly_reads(&assembled.finalized_reads);
    kv("assembly_reads", reads.len().to_string());
    for (i, r) in reads.iter().enumerate() {
        kv(
            "asm_read",
            format!(
                "idx={i}\tlen={}\thash={}\tminq={}",
                r.bases.len(),
                fnv1a64_hex(&r.bases),
                r.base_quals.iter().copied().min().unwrap_or(0)
            ),
        );
    }

    let assembler = args.assembler.clone();
    kv(
        "assembler",
        format!(
            "kmers={:?}\tuse_seq_graph={}\tmin_prune={}\trecover_dangling={}\tnum_best={}",
            assembler.kmer_sizes,
            assembler.use_seq_graph,
            assembler.min_prune_factor,
            assembler.recover_dangling_branches,
            assembler.num_best_haplotypes_per_graph
        ),
    );

    let probes = probe_seq_graph_kmer_attempts(&graph_ref, &reads, &assembler).expect("probe");
    for p in &probes {
        kv(
            "seq_probe",
            format!(
                "phase={}\tkmer={}\toutcome={}\trt_nodes={}\trt_edges={}\tcleanup={}\tkbest_n={}\textracted={}\tnon_ref={}",
                p.phase,
                p.kmer_size,
                p.outcome,
                p.thread_nodes,
                p.thread_edges,
                p.cleanup_status,
                p.kbest_paths,
                p.extracted_haps,
                p.non_ref_haps
            ),
        );
    }

    dump_kmer_graph(
        10,
        &graph_ref,
        &reads,
        &assembler,
        &extra_bases,
        false,
        false,
    );
    dump_kmer_graph(
        25,
        &graph_ref,
        &reads,
        &assembler,
        &extra_bases,
        false,
        false,
    );
    dump_kmer_graph(
        35,
        &graph_ref,
        &reads,
        &assembler,
        &extra_bases,
        false,
        false,
    );

    if let Ok(Some(audit)) =
        audit_threading_dangling_recovery(&graph_ref, &reads, 35, &assembler, false, false)
    {
        kv(
            "dangling_audit_k35",
            format!(
                "edges_before={}\tedges_after={}\ttails_rec={}/{}\theads_rec={}/{}",
                audit.edges_before,
                audit.edges_after,
                audit.tails_recovered,
                audit.tails_attempted,
                audit.heads_recovered,
                audit.heads_attempted
            ),
        );
    }

    for &(kmer, allow_nu, label) in &[
        (10usize, false, "configured"),
        (25, false, "configured"),
        (35, false, "seq_expanded_iter1"),
        (35, true, "supplement_expanded"),
    ] {
        dump_rt_extract(
            kmer,
            &graph_ref,
            &reads,
            &assembler,
            &extra_bases,
            false,
            allow_nu,
            label,
        );
    }
}

fn dump_kmer_graph(
    kmer: usize,
    graph_ref: &gatk_haplotypecaller::assembly::AssemblyRead,
    reads: &[gatk_haplotypecaller::assembly::AssemblyRead],
    assembler: &ReadThreadingAssemblerArgs,
    extra_bases: &[u8],
    allow_low_complexity: bool,
    allow_non_unique_ref: bool,
) {
    let extra_hash = fnv1a64_hex(extra_bases);
    kv("kmer_begin", kmer.to_string());
    let rt = match build_threading_graph_for_seq_assembly(
        graph_ref,
        reads,
        kmer,
        assembler,
        allow_low_complexity,
        allow_non_unique_ref,
    ) {
        Ok(Some(g)) => g,
        Ok(None) => {
            kv("rt", format!("kmer={kmer}\tstatus=none"));
            return;
        }
        Err(e) => {
            kv("rt", format!("kmer={kmer}\tstatus=err\t{e}"));
            return;
        }
    };
    kv(
        "rt",
        format!(
            "kmer={kmer}\tnodes={}\tedges={}",
            rt.node_count(),
            rt.edge_count()
        ),
    );
    let non_ref_edges = rt
        .edges_sorted()
        .into_iter()
        .filter(|e| !rt.edge_is_ref(e.from, e.to))
        .count();
    kv("rt_nonref_edges", format!("kmer={kmer}\tn={non_ref_edges}"));

    let mut seq = SeqGraph::from_assembly_graph(&rt);
    seq.clean_non_ref_paths();
    let status = seq.cleanup_seq_graph();
    kv(
        "seq",
        format!(
            "kmer={kmer}\tstatus={status:?}\tnodes={}\tedges={}",
            seq.node_count(),
            seq.edge_count()
        ),
    );
    for (i, v) in seq.vertices().iter().enumerate() {
        kv(
            "seq_vtx",
            format!(
                "kmer={kmer}\tid={i}\thash={}\tlen={}\tis_source={}\tis_sink={}\tseq={}",
                fnv1a64_hex(&v.sequence),
                v.sequence.len(),
                seq.reference_source_vertex() == Some(i),
                seq.reference_sink_vertex() == Some(i),
                String::from_utf8_lossy(&v.sequence)
            ),
        );
    }
    for e in seq.edges() {
        kv(
            "seq_edge",
            format!(
                "kmer={kmer}\tfrom={}\tto={}\ttotal_mult={}\tis_ref={}",
                e.from, e.to, e.support, e.is_ref
            ),
        );
    }

    let paths = find_best_haplotypes_seq_graph(&seq, assembler.num_best_haplotypes_per_graph)
        .unwrap_or_default();
    kv("kbest_n", format!("kmer={kmer}\tn={}", paths.len()));
    let mut extra_rank = None;
    for (rank, p) in paths.iter().enumerate() {
        let bases = seq.path_bases_bytes(p.start, &p.edges);
        let hash = fnv1a64_hex(&bases);
        kv(
            "kbest",
            format!(
                "kmer={kmer}\trank={rank}\thash={hash}\tscore={}\tisRef={}\tlen={}\tedges={}",
                p.score,
                p.is_reference,
                bases.len(),
                p.edges.len()
            ),
        );
        if hash == extra_hash {
            extra_rank = Some(rank);
            kv(
                "extra_path",
                format!(
                    "kmer={kmer}\trank={rank}\tstart={}\tedges={:?}",
                    p.start, p.edges
                ),
            );
            for (from, to) in &p.edges {
                let edge = seq.edges().iter().find(|e| e.from == *from && e.to == *to);
                kv(
                    "extra_path_edge",
                    format!(
                        "kmer={kmer}\tfrom={from}\tto={to}\ttotal_mult={}\tis_ref={}",
                        edge.map(|e| e.support).unwrap_or(0),
                        edge.map(|e| e.is_ref).unwrap_or(false)
                    ),
                );
            }
        }
    }
    kv(
        "extra_in_kbest",
        format!("kmer={kmer}\tpresent={}", extra_rank.is_some()),
    );
}

fn dump_rt_extract(
    kmer: usize,
    graph_ref: &gatk_haplotypecaller::assembly::AssemblyRead,
    reads: &[gatk_haplotypecaller::assembly::AssemblyRead],
    assembler: &ReadThreadingAssemblerArgs,
    extra_bases: &[u8],
    allow_low_complexity: bool,
    allow_non_unique_ref: bool,
    label: &str,
) {
    let extra_hash = fnv1a64_hex(extra_bases);
    let before = extract_rt_haplotypes_before_remove_paths(
        graph_ref,
        reads,
        assembler,
        kmer,
        allow_low_complexity,
        allow_non_unique_ref,
    )
    .unwrap_or_default();
    let after = extract_rt_haplotypes_after_remove_paths(
        graph_ref,
        reads,
        assembler,
        kmer,
        allow_low_complexity,
        allow_non_unique_ref,
    )
    .unwrap_or_default();
    kv(
        "rt_extract",
        format!(
            "label={label}\tkmer={kmer}\tallow_nu={allow_non_unique_ref}\tbefore_n={}\tafter_n={}\tbefore_has_extra={}\tafter_has_extra={}",
            before.len(),
            after.len(),
            before.iter().any(|h| fnv1a64_hex(&h.bases) == extra_hash),
            after.iter().any(|h| fnv1a64_hex(&h.bases) == extra_hash)
        ),
    );
    for (stage, batch) in [("before_remove", &before), ("after_remove", &after)] {
        for (i, h) in batch.iter().enumerate() {
            let hash = fnv1a64_hex(&h.bases);
            let cigar = h
                .cigar
                .as_ref()
                .map(|c| c.to_gatk_string())
                .unwrap_or_else(|| ".".into());
            kv(
                "rt_hap",
                format!(
                    "label={label}\tstage={stage}\tkmer={kmer}\tidx={i}\thash={hash}\tisRef={}\tlen={}\tscore={}\tkmer_sz={}\tcigar={cigar}",
                    h.is_reference,
                    h.bases.len(),
                    h.score,
                    h.kmer_size
                ),
            );
        }
    }
}
