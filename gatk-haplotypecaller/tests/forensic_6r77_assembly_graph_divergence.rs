//! 6R.77: isolate assembly/graph origin of J0/J1 vs R0.
//!
//! Forensic only. No production change in this file.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r77_assembly_graph_divergence
//! HOLDOUT_6R77=1 cargo test -p gatk-haplotypecaller --test forensic_6r77_assembly_graph_divergence live_ -- --nocapture
//! ```

use gatk_haplotypecaller::{
    AssemblyGraphParams, AssemblyRead, Haplotype, KmerSize, ReadThreadingAssemblerArgs,
};

const JAVA_ONLY_J0: &[u8] = b"CATGGAGCCTGACCTTATTTGAAGTAGGGCATTTGCAGATGTATTTAAGATATTTGAGGCTGGGCACAGTGGCTCACGTCTGTAATCCCAGCACTTTGAAAGGCCGAGGCAGGTGGATTCACCTGAGGTCAGGAGTTTGAGACCAGCCTGTCCCACATGGTGAAAAGCCCGTATCTACCAAAAATACAAACGTTAGCTGTGTGTGGTGGTGGCGGCACCTGTAATCCCAGCTACTCGAGAGCCAGAG";
const JAVA_ONLY_J1: &[u8] = b"CATGGAGCCTGACTTTATTTGAAGTAGGGCATTTGCAGATGTATTTAAGATATTTGAGGCCGGGCACAGTGGCTCATGTCTGTAATCCCAGCACTTTAAAAGGCTGAGGCAGGTGTATTCACCTGAGGTCAGGAGTTCGAGACCAGCCTGGCCAACATGGTGAAAAGCCCGTATCTACCAAAAATACAAAAGTTAGCTGGGTGTGGTGGTGGCGGCACCTGTAATCCCAGCTACTCGAGAGCCAGAG";
const RUST_ONLY_R0: &[u8] = b"CATGGAGCCTGACCTTATTTGAAGTAGGGCATTTGCAGATGTATTTAAGATATTTGAGGCTGGGCACAGTGGCTCACGTCTGTAATCCCAGCACTTTGAAAGGCCGAGGCAGGTGGATTCACCTGAGGTCAGGAGTTTGAGACCAGCCTGTCCCACATGGTGAAAAGCCCGTATCTACCAAAAATACAAACGTTAGCTGTGTGTGGTGGTGGCGGCACCTGTAATCGCAGCTACTCGAGAGCCAGAG";

/// Synthetic 25-mers that differ only at the last base (C vs G), matching the
/// J0/R0 Hamming-1 at the divergent site. Graph construction must be able to
/// distinguish them; this is the coordinate-free stand-in for 20:29456494.
const K25_C: &[u8] = b"ACCTGTAATCCCAGCTACTCGAGAG";
const K25_G: &[u8] = b"ACCTGTAATCGCAGCTACTCGAGAG";

fn contains_slice(h: &Haplotype, needle: &[u8]) -> bool {
    h.bases.windows(needle.len()).any(|w| w == needle)
}

fn kmer_in_graph(graph: &gatk_haplotypecaller::AssemblyGraph, kmer: &[u8]) -> bool {
    graph.nodes().iter().any(|n| n.kmer.as_ref() == kmer)
}

fn kmer_support(graph: &gatk_haplotypecaller::AssemblyGraph, kmer: &[u8]) -> u32 {
    graph
        .nodes()
        .iter()
        .find(|n| n.kmer.as_ref() == kmer)
        .map(|n| n.support)
        .unwrap_or(0)
}

/// C- and G-25-mers that differ at one base are distinct graph vertices when both
/// sequences are threaded. Proves the J0/R0 SNP is a graph-input k-mer pair, not a
/// trim artifact.
#[test]
fn forensic_6r77_c_and_g_25mers_are_distinct_graph_vertices() {
    assert_eq!(K25_C.len(), 25);
    assert_eq!(K25_G.len(), 25);
    assert_eq!(
        K25_C
            .iter()
            .zip(K25_G.iter())
            .filter(|(a, b)| a != b)
            .count(),
        1
    );
    assert_eq!(K25_C[10], b'C');
    assert_eq!(K25_G[10], b'G');

    let mut ref_bases = vec![b'T'; 40];
    ref_bases[8..33].copy_from_slice(K25_G);
    let reference = AssemblyRead {
        base_quals: vec![30; ref_bases.len()],
        bases: ref_bases,
    };
    let mut alt_bases = reference.bases.clone();
    alt_bases[8..33].copy_from_slice(K25_C);
    let alt = AssemblyRead {
        base_quals: vec![30; alt_bases.len()],
        bases: alt_bases,
    };
    let mut params = AssemblyGraphParams::default();
    params.kmer_size = KmerSize::try_new(25).unwrap();
    params.min_base_quality = 10;
    params.min_edge_weight = 1;
    params.start_threading_only_at_existing_vertex = false;
    let graph = gatk_haplotypecaller::assembly_graph_from_ref_and_reads_threading(
        &reference,
        &[alt],
        &params,
    )
    .expect("thread");
    assert!(
        kmer_in_graph(&graph, K25_G),
        "G-25-mer (R0/ref allele) must exist as a vertex"
    );
    assert!(
        kmer_in_graph(&graph, K25_C),
        "C-25-mer (J0 allele) must exist as a vertex when a read carries it"
    );
    assert_ne!(
        kmer_support(&graph, K25_C),
        0,
        "C-25-mer must have threading support"
    );
}

/// Java `ReadThreadingAssembler.numHaplotypesInPopulation` / `--kmer-size` KBest cap is 128.
#[test]
fn forensic_6r77_default_kbest_cap_is_128() {
    let args = ReadThreadingAssemblerArgs::default();
    assert_eq!(args.num_best_haplotypes_per_graph, 128);
    assert_eq!(args.kmer_sizes, vec![10, 25]);
}

#[test]
fn forensic_6r77_j0_r0_discriminated_by_25mer_at_offset_226() {
    assert_eq!(JAVA_ONLY_J0.len(), 247);
    assert_eq!(RUST_ONLY_R0.len(), 247);
    let diffs: Vec<usize> = JAVA_ONLY_J0
        .iter()
        .zip(RUST_ONLY_R0)
        .enumerate()
        .filter(|(_, (a, b))| a != b)
        .map(|(i, _)| i)
        .collect();
    assert_eq!(diffs, vec![226]);
    let start = 226usize.saturating_sub(10);
    assert_eq!(&JAVA_ONLY_J0[start..start + 25], K25_C);
    assert_eq!(&RUST_ONLY_R0[start..start + 25], K25_G);
}

/// Live region: assembly inputs, k-mer evidence, per-k haplotype sets.
#[test]
fn live_assembly_graph_divergence() {
    if std::env::var("HOLDOUT_6R77").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R77=1");
        return;
    }
    use gatk_core::reference::{
        parse_intervals_cli_string, ReferenceWindowCache, SequenceDictionary,
    };
    use gatk_haplotypecaller::assembly_region_finalize::{
        create_graph_reference_read, records_to_assembly_reads,
    };
    use gatk_haplotypecaller::{
        assemble_from_ref_and_reads, assemble_reads_with_finalized, call_disposition,
        flatten_assembly_regions, probe_seq_graph_kmer_attempts, query_index_at_reference_position,
        reference_has_non_unique_kmers, traverse_assembly_region_walker,
        AssemblyRegionCallDisposition, CallRegionArgs, ReadFilterParams, WalkerTraversalConfig,
    };
    use rust_htslib::bam::record::CigarString;
    use std::collections::BTreeMap;
    use std::path::Path;

    const INTERVAL: &str = "20:29455000-29456500";
    const BAM_REL: &str = "parity/giab/runs/local-pairhmm-diff/HG001.20-29455000-29456500.bam";
    const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
    const POS_SNP: u64 = 29_456_344;
    const POS_DIV: u64 = 29_456_494;

    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("..");
    let ref_fasta = root.join(REF_REL);
    let bam = root.join(BAM_REL);
    if !ref_fasta.is_file() || !bam.is_file() {
        eprintln!("skip: live BAM/ref missing");
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
    let covering: Vec<_> = regions
        .iter()
        .filter(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= POS_SNP
                && r.end.get() >= POS_SNP
        })
        .collect();
    assert_eq!(covering.len(), 1);
    let region = covering[0];
    let args = CallRegionArgs::strict_java();
    let mut ref_cache = ReferenceWindowCache::new(ref_fasta.clone(), 4);
    let mut owned = region.clone();
    let assembled =
        assemble_reads_with_finalized(&mut owned, &dict, &mut ref_cache, &args.assemble)
            .expect("assemble");

    eprintln!(
        "REGION active={}..{} extended={}..{} n_region_reads={} n_finalized={}",
        region.start.get(),
        region.end.get(),
        region.extended_start.get(),
        region.extended_end.get(),
        region.reads.len(),
        assembled.finalized_reads.len()
    );
    eprintln!(
        "UNTRIMMED n={} unique_k={:?}",
        assembled.assembly.haplotypes.len(),
        {
            let mut ks: Vec<usize> = assembled
                .assembly
                .haplotypes
                .iter()
                .map(|h| h.kmer_size)
                .collect();
            ks.sort_unstable();
            ks.dedup();
            ks
        }
    );

    let mut pile: BTreeMap<char, Vec<(String, u8, u8, bool)>> = BTreeMap::new();
    let pos0 = (POS_DIV - 1) as i64;
    for rec in &assembled.finalized_reads {
        if rec.is_unmapped() {
            continue;
        }
        let cigar = CigarString(rec.cigar().iter().copied().collect());
        let Some(qi) = query_index_at_reference_position(rec.pos(), &cigar, pos0) else {
            continue;
        };
        let seq = rec.seq();
        if qi >= seq.len() {
            continue;
        }
        let b = (seq[qi] as char).to_ascii_uppercase();
        let q = rec.qual().get(qi).copied().unwrap_or(0);
        pile.entry(b).or_default().push((
            String::from_utf8_lossy(rec.qname()).into_owned(),
            q,
            rec.mapq(),
            rec.is_reverse(),
        ));
    }
    for (b, rows) in &pile {
        let n_rev = rows.iter().filter(|r| r.3).count();
        let mean_mq = rows.iter().map(|r| r.2 as f64).sum::<f64>() / rows.len() as f64;
        eprintln!(
            "PILEUP {POS_DIV} {b} n={} mean_bq={:.1} mean_mq={:.1} n_rev={n_rev} n_fwd={} qnames={:?}",
            rows.len(),
            rows.iter().map(|r| r.1 as f64).sum::<f64>() / rows.len() as f64,
            mean_mq,
            rows.len() - n_rev,
            rows.iter().map(|r| r.0.as_str()).collect::<Vec<_>>()
        );
    }
    let n_c = pile.get(&'C').map(|v| v.len()).unwrap_or(0);
    let n_g = pile.get(&'G').map(|v| v.len()).unwrap_or(0);
    eprintln!(
        "PILEUP_SUMMARY C={n_c} G={n_g} other_alleles={}",
        pile.len()
    );

    let padded_ref = gatk_haplotypecaller::assembly_region_finalize::assembly_reference_read(
        &dict,
        &mut ref_cache,
        region,
    )
    .expect("pad ref");
    let graph_ref = create_graph_reference_read(&padded_ref, region, &dict);
    let graph_reads = records_to_assembly_reads(&assembled.finalized_reads);
    eprintln!(
        "GRAPH_REF len={} start={} k10_nonunique={} k25_nonunique={}",
        graph_ref.bases.len(),
        region.extended_start.get(),
        reference_has_non_unique_kmers(&graph_ref, 10),
        reference_has_non_unique_kmers(&graph_ref, 25),
    );

    let off = (POS_DIV - region.extended_start.get()) as usize;
    eprintln!(
        "REF_AT_DIV off={off} base={}",
        graph_ref.bases.get(off).copied().unwrap_or(b'?') as char
    );
    if off >= 10 && off + 15 < graph_ref.bases.len() {
        let gmer = &graph_ref.bases[off - 10..off + 15];
        let mut alt_g = gmer.to_vec();
        alt_g[10] = b'G';
        let mut params = AssemblyGraphParams::default();
        params.kmer_size = KmerSize::try_new(25).unwrap();
        params.min_base_quality = 10;
        params.min_edge_weight = 1;
        params.start_threading_only_at_existing_vertex = false;
        let raw = gatk_haplotypecaller::assembly_graph_from_ref_and_reads_threading(
            &graph_ref,
            &graph_reads,
            &params,
        )
        .expect("k25 graph");
        eprintln!(
            "REF_25MER (fasta allele at div)={} G_ALLELE_25MER={}",
            String::from_utf8_lossy(gmer),
            String::from_utf8_lossy(&alt_g)
        );
        eprintln!(
            "K25_RAW_GRAPH nodes={} edges={} has_fasta_C_kmer={} has_G_kmer={} supp_C={} supp_G={}",
            raw.node_count(),
            raw.edges_sorted().len(),
            kmer_in_graph(&raw, gmer),
            kmer_in_graph(&raw, &alt_g),
            kmer_support(&raw, gmer),
            kmer_support(&raw, &alt_g),
        );
    }

    let mut assembler = args.assemble.assembler.clone();
    assembler.dangling_java_exact = true;
    let probes =
        probe_seq_graph_kmer_attempts(&graph_ref, &graph_reads, &assembler).expect("probe");
    for p in &probes {
        if p.kmer_size == 10 || p.kmer_size == 25 {
            eprintln!(
                "PROBE k={} phase={} outcome={} nodes={} edges={} cleanup={} kbest={} extracted={} nonref={}",
                p.kmer_size,
                p.phase,
                p.outcome,
                p.thread_nodes,
                p.thread_edges,
                p.cleanup_status,
                p.kbest_paths,
                p.extracted_haps,
                p.non_ref_haps
            );
        }
    }

    let mut k128_j0 = 0usize;
    let mut k128_j1 = 0usize;
    let mut k128_r0 = 0usize;
    let mut k512_j0 = 0usize;
    let mut k512_j1 = 0usize;
    let mut k512_r0 = 0usize;
    for (k, khaps) in [(25usize, 128usize), (25, 512), (35, 128)] {
        let mut a = assembler.clone();
        a.kmer_sizes = vec![k];
        a.dont_increase_kmer_sizes_for_cycles = true;
        a.use_seq_graph = true;
        a.num_best_haplotypes_per_graph = khaps;
        match assemble_from_ref_and_reads(&graph_ref, &graph_reads, &a) {
            Ok(res) => {
                let mut dump_hits = |label: &str, needle: &[u8]| {
                    for (hi, h) in res.haplotypes.iter().enumerate() {
                        if contains_slice(h, needle) {
                            eprintln!(
                                "K{k}_KBEST{khaps} {label} hap[{hi}] len={} k={} score={:.6} ref={} cigar={}",
                                h.bases.len(),
                                h.kmer_size,
                                h.score,
                                h.is_reference,
                                h.cigar
                                    .as_ref()
                                    .map(|c| c.to_gatk_string())
                                    .unwrap_or_default()
                            );
                        }
                    }
                };
                dump_hits("J0", JAVA_ONLY_J0);
                dump_hits("J1", JAVA_ONLY_J1);
                dump_hits("R0", RUST_ONLY_R0);
                let nj0 = res
                    .haplotypes
                    .iter()
                    .filter(|h| contains_slice(h, JAVA_ONLY_J0))
                    .count();
                let nj1 = res
                    .haplotypes
                    .iter()
                    .filter(|h| contains_slice(h, JAVA_ONLY_J1))
                    .count();
                let nr0 = res
                    .haplotypes
                    .iter()
                    .filter(|h| contains_slice(h, RUST_ONLY_R0))
                    .count();
                eprintln!(
                    "K{k}_KBEST{khaps} status={:?} n={} n_j0={nj0} n_j1={nj1} n_r0={nr0}",
                    res.status,
                    res.haplotypes.len()
                );
                if k == 25 && khaps == 128 {
                    k128_j0 = nj0;
                    k128_j1 = nj1;
                    k128_r0 = nr0;
                    let off_canon = (POS_SNP - region.extended_start.get()) as usize;
                    let mut combos: BTreeMap<(char, char), usize> = BTreeMap::new();
                    let mut best_j0 = (usize::MAX, 0usize, Vec::new());
                    for (hi, h) in res.haplotypes.iter().enumerate() {
                        if h.bases.len() > off && h.bases.len() > off_canon {
                            let aa = h.bases[off_canon] as char;
                            let bb = h.bases[off] as char;
                            *combos.entry((aa, bb)).or_default() += 1;
                        }
                        if h.bases.len() >= JAVA_ONLY_J0.len() {
                            for w in h.bases.windows(JAVA_ONLY_J0.len()) {
                                let diffs: Vec<usize> = w
                                    .iter()
                                    .zip(JAVA_ONLY_J0)
                                    .enumerate()
                                    .filter(|(_, (x, y))| x != y)
                                    .map(|(i, _)| i)
                                    .collect();
                                if diffs.len() < best_j0.0 {
                                    best_j0 = (diffs.len(), hi, diffs);
                                }
                            }
                        }
                    }
                    eprintln!("K25 combos (base@{POS_SNP}, base@{POS_DIV}) = {combos:?}");
                    eprintln!(
                        "K25 nearest J0 ham={} hap[{}] diff_idx={:?}",
                        best_j0.0, best_j0.1, best_j0.2
                    );
                }
                if k == 25 && khaps == 512 {
                    k512_j0 = nj0;
                    k512_j1 = nj1;
                    k512_r0 = nr0;
                }
            }
            Err(e) => eprintln!("K{k}_KBEST{khaps} err={e}"),
        }
    }
    assert_eq!(k128_j0, 1, "6R.80: K=128 retains J0 (Java rank 124)");
    let _ = k128_j1;
    assert_eq!(k128_r0, 0, "6R.80: R0 is outside production K=128");
    assert!(k512_j0 >= 1, "J0 exists as a rust k=25 graph path");
    let _ = k512_j1;
    let _ = k512_r0;

    let n_j0 = assembled
        .assembly
        .haplotypes
        .iter()
        .filter(|h| contains_slice(h, JAVA_ONLY_J0))
        .count();
    let n_j1 = assembled
        .assembly
        .haplotypes
        .iter()
        .filter(|h| contains_slice(h, JAVA_ONLY_J1))
        .count();
    let n_r0 = assembled
        .assembly
        .haplotypes
        .iter()
        .filter(|h| contains_slice(h, RUST_ONLY_R0))
        .count();
    eprintln!("PRODUCTION_UNTRIMMED n_j0={n_j0} n_j1={n_j1} n_r0={n_r0}");
    assert_eq!(n_j0, 1, "6R.80: J0 retained in rust assembly");
    let _ = n_j1;
    assert_eq!(n_r0, 0, "6R.80: R0 displaced from production K=128");
}
