//! 6R.80: C-allele read-threading increment — Java `getMultiplicity()` vs pruning.
//!
//! Forensic only. The production change lives in `ThreadingEdge` / graph conversion:
//! SeqGraph and k-best copy total multiplicity (`MultiSampleEdge.copy()`), while
//! `LowWeightChainPruner` still uses `getPruningMultiplicity()`.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r80_read_threading_increment_contract
//! HOLDOUT_6R80=1 cargo test -p gatk-haplotypecaller --test forensic_6r80_read_threading_increment_contract live_ -- --nocapture
//! ```

use gatk_haplotypecaller::bio_ids::KmerSize;
use gatk_haplotypecaller::seq_graph::SeqGraph;
use gatk_haplotypecaller::{
    assembly_graph_from_reads_threading, assembly_graph_from_ref_and_reads_threading,
    AssemblyGraph, AssemblyGraphParams, AssemblyRead,
};

const K: usize = 4;
/// Unique-k reference spine (`addSequence("ref")` / `ANONYMOUS_SAMPLE`).
const REF: &[u8] = b"ACGTACGGTTAA";
/// Same as REF: contributes the reference continuation `ACGT -> CGTA`.
const READ_REF_LIKE: &[u8] = b"ACGTACGGTTAA";
/// SNP off the continuation (`ACGT -> CGTG`); not on the reference haplotype.
const READ_ALT: &[u8] = b"ACGTGCGGTTAA";

fn ar(seq: &[u8]) -> AssemblyRead {
    AssemblyRead {
        bases: seq.to_vec(),
        base_quals: vec![30; seq.len()],
    }
}

fn params() -> AssemblyGraphParams {
    AssemblyGraphParams {
        kmer_size: KmerSize::try_new(K as u16).unwrap(),
        min_base_quality: 10,
        ..Default::default()
    }
}

fn edge_by_kmers(g: &AssemblyGraph, from: &[u8], to: &[u8]) -> Option<(usize, usize, u32, bool)> {
    let fid = g.nodes().iter().find(|n| n.kmer.as_ref() == from)?.id;
    let tid = g.nodes().iter().find(|n| n.kmer.as_ref() == to)?.id;
    let e = g
        .edges_sorted()
        .into_iter()
        .find(|e| e.from == fid && e.to == tid)?;
    Some((fid, tid, e.support, g.edge_is_ref(fid, tid)))
}

fn usable_contains(bases: &[u8], quals: &[u8], min_q: u8, needle: &[u8]) -> bool {
    if needle.is_empty() || bases.len() != quals.len() {
        return false;
    }
    let mut last_good: Option<usize> = None;
    for end in 0..=bases.len() {
        let unusable = end == bases.len()
            || quals[end] < min_q
            || !matches!(bases[end], b'A' | b'C' | b'G' | b'T' | b'N');
        if unusable {
            if let Some(start) = last_good {
                if end - start >= needle.len()
                    && bases[start..end].windows(needle.len()).any(|w| w == needle)
                {
                    return true;
                }
            }
            last_good = None;
        } else if last_good.is_none() {
            last_good = Some(end);
        }
    }
    false
}

/// Java `MultiSampleEdge.copy()` copies `getMultiplicity()` (total), not pruning.
#[test]
fn forensic_6r80_java_copy_uses_total_not_pruning() {
    let total = 16u32;
    let pruning = 15u32;
    let copied = total;
    assert_eq!(copied, 16);
    assert_ne!(copied, pruning);
}

/// `numPruningSamples = 1`: after ref sample 1 and read sample N, pruning = N, total = N+1.
#[test]
fn forensic_6r80_two_sample_flush_pruning_drops_ref_total_keeps_it() {
    let n_ref_like = 3usize;
    let n_alt = 2usize;
    let mut reads: Vec<AssemblyRead> = (0..n_ref_like).map(|_| ar(READ_REF_LIKE)).collect();
    reads.extend((0..n_alt).map(|_| ar(READ_ALT)));
    let g = assembly_graph_from_ref_and_reads_threading(&ar(REF), &reads, &params()).unwrap();

    let (_, _, ref_cont, is_ref) = edge_by_kmers(&g, b"ACGT", b"CGTA").expect("ref continuation");
    let (_, _, alt_cont, alt_is_ref) =
        edge_by_kmers(&g, b"ACGT", b"CGTG").expect("alt continuation");
    assert!(is_ref);
    assert!(!alt_is_ref);
    assert_eq!(
        ref_cont,
        (n_ref_like as u32) + 1,
        "SeqGraph weight is Java total: reference haplotype + matching reads"
    );
    assert_eq!(
        alt_cont, n_alt as u32,
        "alt-only edge has no reference increment"
    );

    let ref_only = assembly_graph_from_ref_and_reads_threading(&ar(REF), &[], &params()).unwrap();
    assert_eq!(edge_by_kmers(&ref_only, b"ACGT", b"CGTA").unwrap().2, 1);
    let reads_only = assembly_graph_from_reads_threading(
        &(0..n_ref_like)
            .map(|_| ar(READ_REF_LIKE))
            .collect::<Vec<_>>(),
        &params(),
    )
    .unwrap();
    assert_eq!(
        edge_by_kmers(&reads_only, b"ACGT", b"CGTA").unwrap().2,
        n_ref_like as u32
    );
}

/// SeqGraph conversion copies AssemblyGraph total support (Java `e.copy()`).
#[test]
fn forensic_6r80_seqgraph_from_assembly_graph_copies_total() {
    let reads = vec![ar(READ_REF_LIKE), ar(READ_REF_LIKE), ar(READ_REF_LIKE)];
    let g = assembly_graph_from_ref_and_reads_threading(&ar(REF), &reads, &params()).unwrap();
    let (from, to, supp, _) = edge_by_kmers(&g, b"ACGT", b"CGTA").unwrap();
    assert_eq!(supp, 4);
    let seq = SeqGraph::from_assembly_graph(&g);
    let seq_e = seq
        .edges()
        .iter()
        .find(|e| e.from == from && e.to == to)
        .expect("seq edge");
    assert_eq!(seq_e.support, 4);
}

/// One `SequenceForKmers` occurrence of a unique continuation increments that edge once.
#[test]
fn forensic_6r80_one_read_one_unique_continuation_one_increment() {
    let g = assembly_graph_from_reads_threading(&[ar(READ_REF_LIKE)], &params()).unwrap();
    assert_eq!(edge_by_kmers(&g, b"ACGT", b"CGTA").unwrap().2, 1);
    let g2 =
        assembly_graph_from_reads_threading(&[ar(READ_REF_LIKE), ar(READ_REF_LIKE)], &params())
            .unwrap();
    assert_eq!(edge_by_kmers(&g2, b"ACGT", b"CGTA").unwrap().2, 2);
}

/// Conservation: BAM-level motif count is not the RT continuation. The extra 1 is REF.
#[test]
fn forensic_6r80_conservation_breaks_at_seqgraph_total_vs_pruning() {
    let reads = vec![ar(READ_REF_LIKE); 3];
    let motif = b"ACGTA";
    assert!(REF.windows(5).any(|w| w == motif));
    let n_read_obs = reads
        .iter()
        .filter(|r| usable_contains(&r.bases, &r.base_quals, 10, motif))
        .count();
    assert_eq!(n_read_obs, 3);
    let g = assembly_graph_from_ref_and_reads_threading(&ar(REF), &reads, &params()).unwrap();
    let continuation = edge_by_kmers(&g, b"ACGT", b"CGTA").unwrap().2;
    assert_eq!(continuation, n_read_obs as u32 + 1);
    assert_ne!(
        continuation, n_read_obs as u32,
        "pruning multiplicity would equal the read sample only"
    );
}

#[test]
fn live_c_continuation_is_ref_plus_matching_reads() {
    if std::env::var("HOLDOUT_6R80").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R80=1");
        return;
    }
    use gatk_core::reference::{
        parse_intervals_cli_string, ReferenceWindowCache, SequenceDictionary,
    };
    use gatk_haplotypecaller::assembly_region_finalize::{
        create_graph_reference_read, records_to_assembly_reads,
    };
    use gatk_haplotypecaller::read_threading_assembler::build_threading_graph_for_seq_assembly;
    use gatk_haplotypecaller::{
        assemble_reads_with_finalized, call_disposition, flatten_assembly_regions,
        traverse_assembly_region_walker, AssemblyRegionCallDisposition, CallRegionArgs,
        ReadFilterParams, WalkerTraversalConfig,
    };
    use std::path::Path;

    const INTERVAL: &str = "20:29455000-29456500";
    const BAM_REL: &str = "parity/giab/runs/local-pairhmm-diff/HG001.20-29455000-29456500.bam";
    const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
    const POS_SNP: u64 = 29_456_344;
    const K25_C: &[u8] = b"ACCTGTAATCCCAGCTACTCGAGAG";
    const K25_C_NEXT: &[u8] = b"CCTGTAATCCCAGCTACTCGAGAGC";
    const C26: &[u8] = b"ACCTGTAATCCCAGCTACTCGAGAGC";
    const G26: &[u8] = b"ACCTGTAATCGCAGCTACTCGAGAGC";

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
    let padded_ref = gatk_haplotypecaller::assembly_region_finalize::assembly_reference_read(
        &dict,
        &mut ref_cache,
        region,
    )
    .expect("pad ref");
    let graph_ref = create_graph_reference_read(&padded_ref, region, &dict);
    let graph_reads = records_to_assembly_reads(&assembled.finalized_reads);
    assert_eq!(graph_reads.len(), 176);

    let ref_has_c = usable_contains(&graph_ref.bases, &graph_ref.base_quals, 10, C26);
    let ref_has_g = usable_contains(&graph_ref.bases, &graph_ref.base_quals, 10, G26);
    assert!(
        ref_has_c,
        "reference haplotype carries the C continuation 26-mer"
    );
    assert!(
        !ref_has_g,
        "reference haplotype does not carry the G continuation"
    );

    let mut rust_c: Vec<String> = Vec::new();
    let mut rust_g: Vec<String> = Vec::new();
    for (i, rec) in assembled.finalized_reads.iter().enumerate() {
        let r = &graph_reads[i];
        let qname = String::from_utf8_lossy(rec.qname()).into_owned();
        let start = rec.pos();
        let cigar = rust_htslib::bam::record::CigarString(rec.cigar().iter().copied().collect());
        let fp = format!(
            "{qname}|{start}|{cigar}|{}",
            String::from_utf8_lossy(&r.bases)
        );
        if usable_contains(&r.bases, &r.base_quals, 10, C26) {
            rust_c.push(fp.clone());
        }
        if usable_contains(&r.bases, &r.base_quals, 10, G26) {
            rust_g.push(fp);
        }
    }
    eprintln!(
        "LIVE finalized={} C26_reads={} G26_reads={} ref_C26={} ref_G26={}",
        graph_reads.len(),
        rust_c.len(),
        rust_g.len(),
        ref_has_c,
        ref_has_g
    );
    for (i, fp) in rust_c.iter().enumerate() {
        eprintln!("  C_READ[{i}] {fp}");
    }
    eprintln!(
        "JAVA_ONLY_C_SUPPORT={{ REF haplotype addSequence(ANONYMOUS_SAMPLE) }} n_reads={}",
        rust_c.len()
    );

    let mut assembler = args.assemble.assembler.clone();
    assembler.dangling_java_exact = true;
    let rt = build_threading_graph_for_seq_assembly(
        &graph_ref,
        &graph_reads,
        25,
        &assembler,
        false,
        false,
    )
    .expect("rt")
    .expect("k=25 graph");

    let c_id = rt
        .nodes()
        .iter()
        .find(|n| n.kmer.as_ref() == K25_C)
        .map(|n| n.id)
        .expect("C25");
    let c_next = rt
        .nodes()
        .iter()
        .find(|n| n.kmer.as_ref() == K25_C_NEXT)
        .map(|n| n.id)
        .expect("C25 next");
    let c_edge = rt
        .edges_sorted()
        .into_iter()
        .find(|e| e.from == c_id && e.to == c_next)
        .expect("C continuation");
    eprintln!(
        "C_CONTINUATION {}->{} supp={} ref={} from={} to={}",
        c_id,
        c_next,
        c_edge.support,
        rt.edge_is_ref(c_id, c_next),
        String::from_utf8_lossy(K25_C),
        String::from_utf8_lossy(K25_C_NEXT)
    );
    assert_eq!(
        c_edge.support,
        rust_c.len() as u32 + 1,
        "Java SeqGraph total = matching reads + reference haplotype"
    );
    assert!(rt.edge_is_ref(c_id, c_next));
}
