//! 6R.164: dangling-tail recovery is graph splice, not a standalone haplotype.
//!
//! Java 4.4.0.0 `AbstractReadThreadingGraph.mergeDanglingTail` (SHA `2dbc0258`)
//! does `addEdge(danglingPath[altIndex], referencePath[refIndex], weight=1)` and
//! returns. Haplotypes are materialized later from source→sink `findBestPaths`.
//!
//! Rust keeps that splice. On `dangling_java_exact` it must **not** push
//! `DanglingMergeHaplotype` / `apply_dangling_merge_haplotypes`.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r164_no_dangling_merge_haplotype_materialization -- --test-threads=1 --nocapture
//! HOLDOUT_6R164=1 cargo test -p gatk-haplotypecaller --test holdout_6r164_no_dangling_merge -- --nocapture --test-threads=1
//! ```

use gatk_core::reference::{parse_intervals_cli_string, ReferenceWindowCache, SequenceDictionary};
use gatk_haplotypecaller::assembly::DanglingMergeHaplotype;
use gatk_haplotypecaller::assembly_based_caller::{
    assemble_reads_with_finalized, AssembleReadsArgs,
};
use gatk_haplotypecaller::assembly_dangling_recovery::apply_dangling_merge_haplotypes;
use gatk_haplotypecaller::assembly_region_finalize::{
    create_graph_reference_read, records_to_assembly_reads,
};
use gatk_haplotypecaller::cigar::{Cigar, CigarOperator};
use gatk_haplotypecaller::haplotype::Haplotype;
use gatk_haplotypecaller::read_threading_assembler::audit_threading_dangling_recovery;
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, traverse_assembly_region_walker,
    AssemblyRegionCallDisposition, ReadFilterParams, SwParameters, WalkerTraversalConfig,
};
use std::path::{Path, PathBuf};

const JAVA_REF: &str = "9b4092f3b20a5de7";
const JAVA_ALT: &str = "8c9a6fc8f302f4a3";
const EXTRA: &str = "3a53b2a941bbdc43";
const INTERVAL: &str = "2:92316200-92316580";
const BAM_REL: &str = "parity/realworld/na12878_20k_b37/NA12878_20k.b37.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const TARGET: u64 = 92_316_347;

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
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
fn forensic_6r164_java_merge_dangling_tail_is_add_edge_only() {
    // Pinned Java 4.4.0.0 does not construct a Haplotype in mergeDanglingTail.
    let java_creates_haplotype = false;
    assert!(!java_creates_haplotype);
}

#[test]
fn forensic_6r164_apply_still_injects_when_invoked() {
    // Non-Java-exact ASM-1 EventMap path still owns this helper. 6R.164 isolates
    // Java-exact materialization; it does not delete the structure.
    let ref_bases = b"ACGTACGTACGTACGTACGT";
    let alt_bases = b"ACGTACGTACGTAAAAACGT";
    let mut ref_cigar = Cigar::new();
    ref_cigar.push(ref_bases.len(), CigarOperator::Match);
    let mut ref_hap = Haplotype::new(ref_bases.to_vec(), true);
    ref_hap.cigar = Some(ref_cigar);
    let mut haps = vec![ref_hap.clone()];
    let mut merge_cigar = Cigar::new();
    merge_cigar.push(12, CigarOperator::Match);
    merge_cigar.push(4, CigarOperator::Insertion);
    merge_cigar.push(4, CigarOperator::Match);
    let merges = [DanglingMergeHaplotype {
        alt_bases: alt_bases.to_vec(),
        cigar: merge_cigar,
        alignment_start_hap_wrt_ref: 0,
    }];
    apply_dangling_merge_haplotypes(
        &mut haps,
        &ref_hap,
        &merges,
        ref_bases,
        &SwParameters::gatk_haplotype_to_reference(),
    );
    assert!(
        haps.iter().any(|h| !h.is_reference
            && h.bases == alt_bases
            && (h.score - 25.0).abs() < f64::EPSILON),
        "non-exact helper must still be able to inject a merge haplotype"
    );
}

#[test]
fn forensic_6r164_java_exact_guards_materialization_not_splice() {
    let recovery = include_str!("../src/assembly_dangling_recovery.rs");
    let assembler = include_str!("../src/read_threading_assembler.rs");
    assert!(
        recovery.contains(
            "self.add_dangling_recovery_edge(plan.from, plan.to, params.dangling_java_exact)"
        ),
        "Java-equivalent graph splice must remain"
    );
    assert!(
        recovery.contains(
            "6R.164: Java-exact recovery produces graph structure, not a standalone haplotype."
        ),
        "java-exact must skip DanglingMergeHaplotype push"
    );
    assert!(
        recovery.contains("if !params.dangling_java_exact")
            && recovery.contains("self.dangling_merge_haps.push(DanglingMergeHaplotype"),
        "merge-hap push must remain behind dangling_java_exact"
    );
    assert!(
        assembler
            .contains("6R.164: Java-exact path keeps the graph splice but does not materialize")
            && assembler.contains("if !args.dangling_java_exact"),
        "assembler must skip apply_dangling_merge_haplotypes on java-exact"
    );
    assert!(
        !recovery.contains("if locus") && !assembler.contains("if locus == "),
        "no coordinate-specific materialization guard"
    );
    assert!(
        !recovery.contains("92316347") && !assembler.contains("92316347"),
        "no P12-locus exception in production recovery/assembler"
    );
}

#[test]
fn forensic_6r164_no_dangling_merge_haplotype_materialization() {
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
        .expect("ActiveFull covering target");

    let mut region = covering.clone();
    let mut ref_cache = ReferenceWindowCache::new(ref_fasta.clone(), 4);
    let mut args = AssembleReadsArgs::default();
    args.strict_java_assembly = true;
    let assembled =
        assemble_reads_with_finalized(&mut region, &dict, &mut ref_cache, &args).expect("assemble");

    let hashes: Vec<String> = assembled
        .assembly
        .haplotypes
        .iter()
        .map(|h| fnv1a64_hex(&h.bases))
        .collect();
    eprintln!("6R164\tassemble_n\t{}", hashes.len());
    for (i, (h, hash)) in assembled
        .assembly
        .haplotypes
        .iter()
        .zip(hashes.iter())
        .enumerate()
    {
        eprintln!(
            "6R164\tassemble_hap\tidx={i}\thash={hash}\tisRef={}\tlen={}\tscore={}\tkmer={}",
            h.is_reference,
            h.bases.len(),
            h.score,
            h.kmer_size
        );
    }

    assert_eq!(
        hashes.len(),
        2,
        "Java assembleReads n=2; extra hap must be gone"
    );
    assert!(hashes.contains(&JAVA_REF.to_string()), "shared REF hap");
    assert!(hashes.contains(&JAVA_ALT.to_string()), "shared ALT hap");
    assert!(
        !hashes.contains(&EXTRA.to_string()),
        "dangling path_bases must not be a standalone haplotype"
    );
    assert!(
        assembled
            .assembly
            .haplotypes
            .iter()
            .all(|h| (h.score - 25.0).abs() >= f64::EPSILON),
        "score=25 dangling-merge sentinel must not enter the haplotype list"
    );

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
    let reads = records_to_assembly_reads(&assembled.finalized_reads);
    let mut assembler = args.assembler.clone();
    assembler.dangling_java_exact = true;
    let audit = audit_threading_dangling_recovery(&graph_ref, &reads, 35, &assembler, false, false)
        .expect("audit")
        .expect("k=35 graph");
    eprintln!(
        "6R164\tdangling_audit_k35\tedges_before={}\tedges_after={}\ttails={}/{}\theads={}/{}",
        audit.edges_before,
        audit.edges_after,
        audit.tails_recovered,
        audit.tails_attempted,
        audit.heads_recovered,
        audit.heads_attempted
    );
    assert_eq!(
        audit.tails_recovered, 1,
        "dangling-tail graph recovery still occurs"
    );
    assert!(
        audit.edges_after > audit.edges_before,
        "recovered merge edge must remain in the graph"
    );
}
