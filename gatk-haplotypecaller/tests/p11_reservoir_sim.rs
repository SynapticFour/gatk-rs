//! One-off parity check: Java reservoir on alt0..alt299 in BAM iteration order.

use gatk_haplotypecaller::read_downsample::{
    apply_positional_downsampler, GatkJavaRng, PositionalDownsamplerConfig,
};
use rust_htslib::bam::{self, Read};
use std::path::Path;

const JAVA_KEPT: &[&str] = &[
    "alt107", "alt115", "alt12", "alt123", "alt131", "alt133", "alt134", "alt137", "alt143",
    "alt147", "alt152", "alt153", "alt160", "alt161", "alt171", "alt180", "alt184", "alt19",
    "alt193", "alt194", "alt208", "alt211", "alt213", "alt22", "alt220", "alt223", "alt231",
    "alt232", "alt233", "alt254", "alt256", "alt257", "alt27", "alt274", "alt277", "alt28",
    "alt280", "alt29", "alt291", "alt31", "alt33", "alt35", "alt53", "alt58", "alt59", "alt63",
    "alt70", "alt91", "alt92", "alt98",
];

#[test]
fn p11_indexed_bam_reservoir_matches_java_dump() {
    let bam = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../parity/build/sam-indexed-bam/p11_java_positive.bam");
    if !bam.is_file() {
        eprintln!("skip: run samtools sort on p11 first");
        return;
    }
    let mut reader = bam::IndexedReader::from_path(&bam).expect("open");
    let header = reader.header().clone();
    let tid = header.tid(b"chrLive").unwrap() as i32;
    reader.fetch((tid, 0, 63)).expect("fetch");
    let mut recs = Vec::new();
    for res in reader.records() {
        recs.push(res.expect("rec"));
    }
    assert_eq!(recs.len(), 300, "expected 300 reads in fixture");

    let mut rng = GatkJavaRng::reset_gatk_default();
    apply_positional_downsampler(
        &mut recs,
        Some(&header),
        &PositionalDownsamplerConfig::gatk_haplotype_caller_defaults(),
        &mut rng,
    )
    .expect("positional downsample");
    assert_eq!(recs.len(), 50);
    let mut kept: Vec<String> = recs
        .iter()
        .map(|r| String::from_utf8_lossy(r.qname()).into_owned())
        .collect();
    kept.sort();
    let mut expected: Vec<String> = JAVA_KEPT.iter().map(|s| (*s).to_string()).collect();
    expected.sort();
    assert_eq!(
        kept, expected,
        "reservoir selection mismatch vs Java b5-reads dump"
    );
}
