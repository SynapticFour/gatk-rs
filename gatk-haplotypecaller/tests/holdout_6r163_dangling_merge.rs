//! 6R.163: live check that frozen dangling path_bases stay a substring of SeqGraph ALT
//! and are not a standalone haplotype after 6R.164.
//! Skipped unless `HOLDOUT_6R163=1`. Production change: NONE.
//!
//! ```text
//! HOLDOUT_6R163=1 cargo test -p gatk-haplotypecaller --test holdout_6r163_dangling_merge -- --nocapture --test-threads=1
//! ```

use gatk_core::reference::{parse_intervals_cli_string, ReferenceWindowCache, SequenceDictionary};
use gatk_haplotypecaller::assembly_based_caller::{
    assemble_reads_with_finalized, AssembleReadsArgs,
};
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
const JAVA_ALT: &str = "8c9a6fc8f302f4a3";
const EXTRA_SEQ: &str = "GTCAGCAGCATTCTCAGAAAGTTCTTTGTGATGATTGCATTCAAGTCACAGAATTGAACATTCCCTTTCACAGAGCAGGTTTGAAACACTCTTTTTGTAGTGTCTGTAAGTGAACATTTGGATTGCTTTCAGGCCTAAGGTGAAAAAGGAAATATCTTCCCATAAAAACTAGAC";

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    println!("6R163\t{key}\t{}", value.as_ref());
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
fn holdout_6r163_dangling_merge() {
    if std::env::var("HOLDOUT_6R163").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R163=1");
        return;
    }
    let root = repo_root();
    let ref_fasta = root.join(REF_REL);
    let bam = root.join(BAM_REL);
    if !ref_fasta.is_file() || !bam.is_file() {
        eprintln!("skip: missing P12 ref/BAM");
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
    let covering = regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= TARGET
                && r.end.get() >= TARGET
        })
        .expect("ActiveFull");
    let mut region = covering.clone();
    let mut ref_cache = ReferenceWindowCache::new(ref_fasta.clone(), 4);
    let mut args = AssembleReadsArgs::default();
    args.strict_java_assembly = true;
    let assembled =
        assemble_reads_with_finalized(&mut region, &dict, &mut ref_cache, &args).expect("assemble");

    let extra_present = assembled
        .assembly
        .haplotypes
        .iter()
        .any(|h| fnv1a64_hex(&h.bases) == EXTRA_FNV);
    let alt = assembled
        .assembly
        .haplotypes
        .iter()
        .find(|h| fnv1a64_hex(&h.bases) == JAVA_ALT)
        .expect("alt hap");
    let extra_s = EXTRA_SEQ;
    let alt_s = String::from_utf8_lossy(&alt.bases);
    let idx = alt_s.find(extra_s);
    kv(
        "assemble_n",
        assembled.assembly.haplotypes.len().to_string(),
    );
    kv("extra_present", extra_present.to_string());
    kv(
        "alt_contains_extra_seq",
        format!(
            "yes={}\toffset={}\textra_len={}\talt_len={}",
            idx.is_some(),
            idx.unwrap_or(usize::MAX),
            extra_s.len(),
            alt.bases.len()
        ),
    );
    kv(
        "java_reachability",
        "B_graph_path_not_source_sink_substring_of_selected_alt",
    );
    assert_eq!(assembled.assembly.haplotypes.len(), 2);
    assert!(
        !extra_present,
        "6R.164: dangling path_bases must not be a standalone haplotype"
    );
    assert!(
        idx.is_some(),
        "frozen dangling path_bases remain a substring of k-best ALT"
    );
    assert_eq!(idx, Some(195));
}
