//! `ReadShard` construction vs GATK `AssemblyRegionWalker` shard semantics.

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::{make_read_shards_default_padding, HaplotypeCallerEngine};
use std::path::PathBuf;

fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../parity/fixtures")
        .join(name)
}

#[test]
fn engine_read_shards_fixture_chr1() {
    let ref_fa = fixture("reference.fa");
    let dict = SequenceDictionary::from_fasta_path(&ref_fa).expect("dict");
    let specs = parse_intervals_cli_string(&dict, "chr1:5-15").expect("parse");
    let engine = HaplotypeCallerEngine::prepare_traversal_default(&dict, specs).expect("engine");
    let shards = engine.read_shards_default_padding(&dict).expect("shards");
    assert_eq!(shards.len(), 1);
    assert_eq!(shards[0].contig, "chr1");
    // len 32, pad 100 → full contig
    assert_eq!(shards[0].padded_spans, vec![(1, 32)]);
}

#[test]
fn standalone_api_matches_engine() {
    let ref_fa = fixture("reference.fa");
    let dict = SequenceDictionary::from_fasta_path(&ref_fa).expect("dict");
    let specs = parse_intervals_cli_string(&dict, "chr1:10-10").expect("parse");
    let a = make_read_shards_default_padding(&dict, &specs).expect("a");
    let engine = HaplotypeCallerEngine::prepare_traversal_default(&dict, specs).expect("engine");
    let b = engine.read_shards_default_padding(&dict).expect("b");
    assert_eq!(a, b);
}
