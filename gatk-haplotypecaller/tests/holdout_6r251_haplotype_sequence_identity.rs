//! 6R.251 live: five TGTTTG haplotype sequences at `20:29455649 T/TGTTTG`.
//! Proof-only. Skipped unless `HOLDOUT_6R251=1`.
//!
//! ```text
//! HOLDOUT_6R251=1 cargo test -p gatk-haplotypecaller --test holdout_6r251_haplotype_sequence_identity -- --nocapture --test-threads=1
//! ```

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::event_map::{
    build_per_haplotype_variation_events, overlapping_events, VariationEvent,
};
use gatk_haplotypecaller::hc_genotyping_engine::{
    take_colocated_merge_numerics, HcGenotypingConfig,
};
use gatk_haplotypecaller::read_threading_assembler::DEFAULT_NUM_BEST_HAPLOTYPES_PER_GRAPH;
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, traverse_assembly_region_walker,
    try_emit_call_region_variants, AssemblyRegionCallDisposition, CallRegionArgs, CigarOperator,
    HaplotypeCallerEngine, HcLikelihoodEngineConfig, PairHmmBackend, ReadFilterParams,
    WalkerTraversalConfig, DEFAULT_STAND_CALL_CONF, DEFAULT_STAND_EMIT_CONFIDENCE,
};
use std::path::{Path, PathBuf};

const INTERVAL: &str = "20:29455000-29456500";
const BAM_REL: &str = "parity/giab/runs/local-pairhmm-diff/HG001.20-29455000-29456500.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const JAVA_VCF_REL: &str = "parity/reports/6r43/chr20_tiny/java.vcf";
const TARGET: u64 = 29_455_649;
const TARGET_REF: &str = "T";
const TARGET_ALT: &str = "TGTTTG";
const CLOSED_GC: u64 = 29_455_314;
const FROZEN_IDX: [usize; 5] = [19, 20, 21, 22, 23];
const FROZEN_FNV: [&str; 5] = [
    "55012fcf3b430591",
    "341b2e070ccb5846",
    "6a65e4c02733c2ed",
    "fa07750bf228b3c2",
    "79451c576721a729",
];

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    println!("6R251\t{key}\t{}", value.as_ref());
}

fn fnv1a64_hex(bases: &[u8]) -> String {
    let mut h: u64 = 0xcbf29ce484222325;
    for &b in bases {
        h ^= b as u64;
        h = h.wrapping_mul(0x100000001b3);
    }
    format!("{h:016x}")
}

fn vcf_has(path: &Path, pos: u64, r: &str, a: &str) -> bool {
    for line in std::fs::read_to_string(path).unwrap_or_default().lines() {
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let f: Vec<_> = line.split('\t').collect();
        if f.len() < 5 || f[0] != "20" {
            continue;
        }
        if f[1].parse::<u64>().ok() == Some(pos) && f[3] == r && f[4] == a {
            return true;
        }
    }
    false
}

fn cigar_insertion_at(h: &gatk_haplotypecaller::Haplotype, hap_pos: usize) -> Option<String> {
    let cigar = h.cigar.as_ref()?;
    let mut pos = 0usize;
    for el in &cigar.elements {
        match el.operator {
            CigarOperator::Insertion => {
                if pos == hap_pos {
                    let end = (pos + el.length).min(h.bases.len());
                    return Some(String::from_utf8_lossy(&h.bases[pos..end]).into_owned());
                }
                pos += el.length;
            }
            CigarOperator::Match | CigarOperator::SoftClip => {
                pos += el.length;
            }
            CigarOperator::Deletion => {}
            _ => {
                if el.operator.consumes_read_bases() {
                    pos += el.length;
                }
            }
        }
    }
    None
}

fn apply_events(window: &[u8], window_start: u64, events: &[VariationEvent]) -> Option<Vec<u8>> {
    let mut sorted: Vec<&VariationEvent> = events.iter().collect();
    sorted.sort_by_key(|e| e.start_1based.get());
    let mut out = window.to_vec();
    for ev in sorted.into_iter().rev() {
        let start = ev.start_1based.get();
        if start < window_start {
            return None;
        }
        let off = (start - window_start) as usize;
        let rb = ev.ref_allele.as_bytes();
        let ab = ev.alt_allele.as_bytes();
        if off + rb.len() > out.len() || &out[off..off + rb.len()] != rb {
            return None;
        }
        out.splice(off..off + rb.len(), ab.iter().copied());
    }
    Some(out)
}

#[test]
fn holdout_6r251_haplotype_sequence_identity() {
    if std::env::var("HOLDOUT_6R251").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R251=1");
        return;
    }
    assert_eq!(DEFAULT_NUM_BEST_HAPLOTYPES_PER_GRAPH, 128);
    assert_eq!(
        HcGenotypingConfig::strict_java().stand_emit_confidence,
        DEFAULT_STAND_CALL_CONF
    );
    let cfg = HcLikelihoodEngineConfig::gatk_haplotype_caller_production();
    assert_eq!(cfg.resolved_pair_hmm_backend(), PairHmmBackend::NeonF64);

    let root = repo_root();
    assert!(vcf_has(
        &root.join(JAVA_VCF_REL),
        TARGET,
        TARGET_REF,
        TARGET_ALT
    ));
    assert!(!vcf_has(&root.join(JAVA_VCF_REL), CLOSED_GC, "G", "C"));

    let ref_fasta = root.join(REF_REL);
    let bam = root.join(BAM_REL);
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
        .expect("ActiveFull");
    let outcome = HaplotypeCallerEngine::call_region(
        region,
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("call")
    .expect("Some");
    let snap = take_colocated_merge_numerics()
        .into_iter()
        .find(|s| s.loc == TARGET)
        .expect("snap");
    assert_eq!(snap.n_reads, 123);
    assert_eq!(snap.pool_sizes[2], 5);

    let haps = &outcome.assembly.haplotypes;
    let ref_hap = haps.iter().find(|h| h.is_reference).expect("ref");
    let trim_start = ref_hap.genome_loc.expect("loc").start_1based();
    let hap_events = build_per_haplotype_variation_events(
        haps,
        outcome.assembly.reference_bases(),
        outcome.assembly.padded_reference_start_1based(),
        outcome.assembly.max_mnp_distance(),
        "20",
    );
    let mut seqs = Vec::new();
    for (i, &idx) in FROZEN_IDX.iter().enumerate() {
        let h = &haps[idx];
        assert_eq!(fnv1a64_hex(&h.bases), FROZEN_FNV[i]);
        assert_eq!(h.bases.len(), 161);
        assert_eq!(cigar_insertion_at(h, 81).as_deref(), Some("GTTTG"));
        let spanning = overlapping_events(hap_events.events_for(idx), TARGET);
        assert!(spanning
            .iter()
            .any(|e| e.ref_allele == TARGET_REF && e.alt_allele == TARGET_ALT));
        let window_events: Vec<VariationEvent> = hap_events
            .events_for(idx)
            .iter()
            .filter(|e| {
                e.start_1based.get() >= trim_start
                    && e.start_1based.get() <= ref_hap.genome_loc.expect("loc").end_1based()
            })
            .cloned()
            .collect();
        let recon = apply_events(&ref_hap.bases, trim_start, &window_events).expect("recon");
        assert_eq!(recon, h.bases);
        seqs.push(h.bases.clone());
    }
    let unique: std::collections::BTreeSet<_> = seqs.iter().cloned().collect();
    assert_eq!(unique.len(), 5);

    let recs =
        try_emit_call_region_variants(region, &outcome, "SAMPLE", DEFAULT_STAND_EMIT_CONFIDENCE)
            .unwrap_or_default();
    assert!(!recs.iter().any(|r| r.position == CLOSED_GC
        && r.reference == "G"
        && r.alternate.iter().any(|a| a == "C")));
    kv("n_reads", "123");
    kv("n_unique_sequences", "5");
    kv("insertion", "GTTTG");
    kv("eventmap_recon", "EXACT");
    kv("classification", "HAPLOTYPE_SEQUENCE_INTERNALLY_CONSISTENT");
    kv("production_change", "NONE");
}
