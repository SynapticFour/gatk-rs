//! Phase 5 live profile: Java EventMap hap signatures vs Rust assembly candidates.
//!
//! Invoked by `scripts/parity/run_p5_live_java_rust_diff.sh` with:
//! - `P5_LIVE_SAM` — fixture SAM path
//! - `P5_LIVE_JAVA_OUT` — Java HC `--debug-assembly` stdout/stderr
//! - `P5_LIVE_JAVA_VCF` — Java HC VCF for the same window
//!
//! Contract (`scripts/parity/README.md`):
//! - Non-empty Java EventMap hap signatures → ≥1 Rust candidate must overlap.
//! - No EventMap signatures → non-applicable only when the Java VCF has no variant records.

use gatk_haplotypecaller::{
    AssemblyGraph, AssemblyGraphParams, AssemblyGraphPruningParams, AssemblyRead, KmerSize,
};
use std::collections::BTreeSet;
use std::env;
use std::fs;
use std::path::Path;

/// Haplotypes under `=== Best Haplotypes ===` whose `Events` map is non-empty.
fn parse_eventmap_hap_signatures(java_out: &str) -> Vec<String> {
    let mut sigs = Vec::new();
    let mut in_best = false;
    let mut pending_seq: Option<String> = None;

    for line in java_out.lines() {
        let Some(msg) = eventmap_message(line) else {
            continue;
        };
        if msg.contains("=== Best Haplotypes ===") {
            in_best = true;
            pending_seq = None;
            continue;
        }
        if !in_best {
            continue;
        }
        // Next assembly-region banner ends the current Best Haplotypes block.
        if msg.contains("=== ") && !msg.contains("Best Haplotypes") {
            in_best = false;
            pending_seq = None;
            continue;
        }
        if is_dna_sequence(&msg) {
            pending_seq = Some(msg.to_ascii_uppercase());
            continue;
        }
        if let Some(rest) = msg.strip_prefix(">> Events = ") {
            if let Some(seq) = pending_seq.take() {
                if event_map_nonempty(rest) {
                    sigs.push(seq);
                }
            }
        }
    }
    sigs
}

fn eventmap_message(line: &str) -> Option<&str> {
    // `16:59:09.016 INFO  EventMap - TGCATG...`
    let idx = line.find("EventMap - ")?;
    Some(line[idx + "EventMap - ".len()..].trim())
}

fn is_dna_sequence(s: &str) -> bool {
    !s.is_empty()
        && s.bytes().all(|b| {
            matches!(
                b,
                b'A' | b'C' | b'G' | b'T' | b'N' | b'a' | b'c' | b'g' | b't' | b'n'
            )
        })
}

fn event_map_nonempty(events_body: &str) -> bool {
    let body = events_body
        .strip_prefix("EventMap{")
        .and_then(|s| s.strip_suffix('}'))
        .unwrap_or(events_body)
        .trim();
    !body.is_empty()
}

fn vcf_has_variant_records(vcf: &str) -> bool {
    vcf.lines()
        .any(|l| !l.starts_with('#') && !l.trim().is_empty())
}

fn sequences_overlap(a: &str, b: &str) -> bool {
    !a.is_empty() && !b.is_empty() && (a.contains(b) || b.contains(a))
}

fn load_assembly_reads_from_sam(path: &Path) -> Vec<AssemblyRead> {
    fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("read {}: {e}", path.display()))
        .lines()
        .filter(|l| !l.starts_with('@') && !l.trim().is_empty())
        .filter_map(|l| {
            let cols: Vec<&str> = l.split('\t').collect();
            if cols.len() < 11 {
                return None;
            }
            let bases = cols[9];
            let quals = cols[10];
            if bases == "*" || quals == "*" {
                return None;
            }
            let base_quals: Vec<u8> = quals.bytes().map(|b| b.saturating_sub(33)).collect();
            Some(AssemblyRead {
                bases: bases.to_ascii_uppercase(),
                base_quals,
            })
        })
        .collect()
}

fn rust_candidate_sequences(sam: &Path) -> BTreeSet<String> {
    let reads = load_assembly_reads_from_sam(sam);
    assert!(
        !reads.is_empty(),
        "SAM {} produced no assembly reads",
        sam.display()
    );
    // Match Java HC default primary k=10 path used by the live debug run.
    let params = AssemblyGraphParams {
        kmer_size: KmerSize::try_new(10).expect("k=10"),
        ..AssemblyGraphParams::default()
    };
    let mut graph = AssemblyGraph::from_reads(&reads, &params).expect("build graph");
    let mut pruning = AssemblyGraphPruningParams::gatk_haplotype_caller_defaults();
    pruning.min_prune_factor = params.min_edge_weight;
    graph.apply_pruning(&pruning);
    graph.remove_dangling_paths(params.dangling_path_max_nodes);
    graph.cleanup_isolated_nodes();
    graph
        .extract_candidate_haplotypes(params.max_haplotypes, params.max_haplotype_bases)
        .into_iter()
        .map(|h| h.sequence.to_ascii_uppercase())
        .collect()
}

#[test]
fn live_java_eventmap_haplotype_signatures_cover_rust_candidates() {
    let Ok(sam) = env::var("P5_LIVE_SAM") else {
        // Not invoked by the live harness — skip in plain `cargo test`.
        return;
    };
    let java_out_path = env::var("P5_LIVE_JAVA_OUT").expect("P5_LIVE_JAVA_OUT");
    let java_vcf_path = env::var("P5_LIVE_JAVA_VCF").expect("P5_LIVE_JAVA_VCF");

    let java_out =
        fs::read_to_string(&java_out_path).unwrap_or_else(|e| panic!("read {java_out_path}: {e}"));
    let java_vcf =
        fs::read_to_string(&java_vcf_path).unwrap_or_else(|e| panic!("read {java_vcf_path}: {e}"));

    let signatures = parse_eventmap_hap_signatures(&java_out);
    if signatures.is_empty() {
        assert!(
            !vcf_has_variant_records(&java_vcf),
            "Java VCF has variant records but EventMap emitted no hap signatures \
             (out={java_out_path}, vcf={java_vcf_path})"
        );
        return;
    }

    let rust = rust_candidate_sequences(Path::new(&sam));
    let overlap = signatures
        .iter()
        .any(|j| rust.iter().any(|r| sequences_overlap(j, r)));
    assert!(
        overlap,
        "no Rust candidate overlaps Java EventMap hap signatures\n\
         java_sigs={signatures:?}\n\
         rust_cands={rust:?}\n\
         sam={sam}"
    );
}

#[test]
fn parses_nonempty_eventmap_as_signature() {
    let log = "\
INFO  EventMap - === Best Haplotypes ===\n\
INFO  EventMap - ACGTACGTAC\n\
INFO  EventMap - > Cigar = 4M1X5M\n\
INFO  EventMap - >> Events = EventMap{pos=5}\n\
";
    assert_eq!(
        parse_eventmap_hap_signatures(log),
        vec!["ACGTACGTAC".to_string()]
    );
}

#[test]
fn empty_eventmap_is_not_a_signature() {
    let log = "\
INFO  EventMap - === Best Haplotypes ===\n\
INFO  EventMap - TGCATGACTGATCGTACGATTCGAGCTAGTCGATCGATGCTAGCTAGGCTAACGTTAGCTAGT\n\
INFO  EventMap - > Cigar = 63M\n\
INFO  EventMap - >> Events = EventMap{}\n\
";
    assert!(parse_eventmap_hap_signatures(log).is_empty());
}
