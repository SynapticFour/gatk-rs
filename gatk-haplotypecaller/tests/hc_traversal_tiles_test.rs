//! Engine produces tiles from `-L` string matching CLI storage.

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::HaplotypeCallerEngine;
use std::path::PathBuf;

fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../parity/fixtures")
        .join(name)
}

#[test]
fn engine_tiles_subinterval_of_chr1() {
    let ref_fa = fixture("reference.fa");
    let dict = SequenceDictionary::from_fasta_path(&ref_fa).expect("dict");
    let specs = parse_intervals_cli_string(&dict, "chr1:1-20").expect("parse");
    let engine = HaplotypeCallerEngine::prepare_traversal_default(&dict, specs).expect("engine");
    assert_eq!(engine.interval_specs.len(), 1);
    assert!(engine.tile_count() >= 1);
    assert_eq!(engine.tiles[0].contig, "chr1");
    assert_eq!(engine.tiles[0].start, 1);
}
