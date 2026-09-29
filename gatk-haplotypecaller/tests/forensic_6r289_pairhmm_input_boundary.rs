//! 6R.289: PairHMM inputs for one read and haplotype 0.
//! Stops at the input boundary. No kernel arithmetic.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r289_pairhmm_input_boundary -- --nocapture --test-threads=1
//! ```

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::hc_genotyping_engine::HcGenotypingConfig;
use gatk_haplotypecaller::likelihood_engine::{
    set_forensic_6r289_substitute, set_forensic_6r289_target, take_forensic_6r289_snaps,
    Forensic6r289Field, Forensic6r289Snap,
};
use gatk_haplotypecaller::read_threading_assembler::DEFAULT_NUM_BEST_HAPLOTYPES_PER_GRAPH;
use gatk_haplotypecaller::{
    begin_hap_list_observe, call_disposition, flatten_assembly_regions, take_hap_list_trim_span,
    traverse_assembly_region_walker, try_emit_call_region_variants, AssemblyRegionCallDisposition,
    CallRegionArgs, HaplotypeCallerEngine, ReadFilterParams, WalkerTraversalConfig,
    DEFAULT_STAND_CALL_CONF, DEFAULT_STAND_EMIT_CONFIDENCE,
};
use std::path::{Path, PathBuf};

const INTERVAL: &str = "20:29455000-29456500";
const BAM_REL: &str = "parity/giab/runs/local-pairhmm-diff/HG001.20-29455000-29456500.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const JAVA_TSV_REL: &str = "gatk-haplotypecaller/tests/6r289_java_pairhmm_inputs.tsv";
const TARGET: u64 = 29_455_649;
const TARGET_REF: &str = "T";
const TARGET_ALT: &str = "TGTTTG";
const QNAME: &str = "HISEQ1:11:H8GV6ADXX:1:2116:18670:99941";
const FLAGS: u16 = 99;
const JAVA_H0: f64 = f64::from_bits(0xc00d9faa00000000);
const RUST_H0: f64 = f64::from_bits(0xc004aa3284709780);
const RUST_CLIP: (u64, u64) = (29_455_569, 29_455_724);

struct JavaBoundary {
    n_haps: usize,
    hap0: Vec<u8>,
    hap0_align: i64,
    hap0_is_ref: bool,
    orig_bases: Vec<u8>,
    orig_start: i64,
    orig_cigar: String,
    proc_bases: Vec<u8>,
    proc_bq: Vec<u8>,
    proc_iq: Vec<u8>,
    proc_dq: Vec<u8>,
    proc_gcp: Vec<u8>,
    proc_start: i64,
    proc_cigar: String,
    proc_mapq: u8,
    proc_len: usize,
    constant_gcp: u8,
}

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    println!("6R289\t{key}\t{}", value.as_ref());
}

fn fmt_f64(x: f64) -> String {
    format!("{x:.17} bits=0x{:016x}", x.to_bits())
}

fn fnv_hex(bytes: &[u8]) -> String {
    let mut h: u64 = 0xcbf29ce484222325;
    for &b in bytes {
        h ^= b as u64;
        h = h.wrapping_mul(0x100000001b3);
    }
    format!("{h:016x}")
}

fn parse_csv(s: &str) -> Vec<u8> {
    if s.is_empty() || s == "." {
        return Vec::new();
    }
    s.split(',').map(|p| p.parse::<u8>().unwrap()).collect()
}

fn kv_after<'a>(parts: &'a [String], key: &str) -> &'a str {
    parts
        .iter()
        .find_map(|p| p.strip_prefix(&format!("{key}=")))
        .unwrap_or("")
}

fn load_java(path: &Path) -> JavaBoundary {
    let text =
        std::fs::read_to_string(path).unwrap_or_else(|e| panic!("read {}: {e}", path.display()));
    let mut n_haps = 0usize;
    let mut hap0 = Vec::new();
    let mut hap0_align = 0i64;
    let mut hap0_is_ref = false;
    let mut constant_gcp = 0u8;
    let mut orig: Option<Vec<String>> = None;
    let mut proc_: Option<Vec<String>> = None;
    for line in text.lines() {
        let f: Vec<String> = line.split('\t').map(|s| s.to_string()).collect();
        if f.first().map(String::as_str) != Some("6R289") || f.len() < 3 {
            continue;
        }
        match f[1].as_str() {
            "n_haps" => n_haps = f[2].parse().unwrap_or(0),
            "hap0_bases" => hap0 = f[2].as_bytes().to_vec(),
            "hap0_align_start" => hap0_align = f[2].parse().unwrap_or(0),
            "hap0_is_ref" => hap0_is_ref = f[2] == "true",
            "constant_gcp" => constant_gcp = f[2].parse().unwrap_or(0),
            "orig" => orig = Some(f[2..].to_vec()),
            "proc" => proc_ = Some(f[2..].to_vec()),
            _ => {}
        }
    }
    let orig = orig.expect("java orig read");
    let proc_ = proc_.expect("java proc read");
    JavaBoundary {
        n_haps,
        hap0,
        hap0_align,
        hap0_is_ref,
        orig_bases: kv_after(&orig, "bases").as_bytes().to_vec(),
        orig_start: kv_after(&orig, "start").parse().unwrap_or(0),
        orig_cigar: kv_after(&orig, "cigar").to_string(),
        proc_bases: kv_after(&proc_, "bases").as_bytes().to_vec(),
        proc_bq: parse_csv(kv_after(&proc_, "bq")),
        proc_iq: parse_csv(kv_after(&proc_, "iq")),
        proc_dq: parse_csv(kv_after(&proc_, "dq")),
        proc_gcp: parse_csv(kv_after(&proc_, "gcp")),
        proc_start: kv_after(&proc_, "start").parse().unwrap_or(0),
        proc_cigar: kv_after(&proc_, "cigar").to_string(),
        proc_mapq: kv_after(&proc_, "mapq").parse().unwrap_or(0),
        proc_len: kv_after(&proc_, "len").parse().unwrap_or(0),
        constant_gcp,
    }
}

fn first_byte_diff(java: &[u8], rust: &[u8]) -> Option<(usize, u8, u8)> {
    let n = java.len().min(rust.len());
    for i in 0..n {
        if java[i] != rust[i] {
            return Some((i, java[i], rust[i]));
        }
    }
    if java.len() != rust.len() {
        Some((
            n,
            java.get(n).copied().unwrap_or(0),
            rust.get(n).copied().unwrap_or(0),
        ))
    } else {
        None
    }
}

fn classification(field: Option<Forensic6r289Field>) -> &'static str {
    match field {
        None => "PAIRHMM_INPUTS_IDENTICAL",
        Some(Forensic6r289Field::ReadBases) => "PAIRHMM_INPUT_READ_BASES",
        Some(Forensic6r289Field::BaseQuals) => "PAIRHMM_INPUT_BASE_QUALITIES",
        Some(Forensic6r289Field::InsQuals) => "PAIRHMM_INPUT_INSERTION_QUALITIES",
        Some(Forensic6r289Field::DelQuals) => "PAIRHMM_INPUT_DELETION_QUALITIES",
        Some(Forensic6r289Field::Gcp) => "PAIRHMM_INPUT_GAP_CONTINUATION",
        Some(Forensic6r289Field::Hap0) => "PAIRHMM_INPUT_HAPLOTYPE_0",
    }
}

struct FlagGuard;
impl Drop for FlagGuard {
    fn drop(&mut self) {
        set_forensic_6r289_target(None, 0);
        set_forensic_6r289_substitute(None, None);
    }
}

fn run_region(
    region: &gatk_haplotypecaller::AssemblyRegion,
    dict: &SequenceDictionary,
    ref_fasta: &Path,
) -> gatk_haplotypecaller::CallRegionOutcome {
    HaplotypeCallerEngine::call_region(region, dict, ref_fasta, &CallRegionArgs::strict_java())
        .expect("call")
        .expect("Some")
}

#[test]
fn forensic_6r289_pairhmm_input_boundary() {
    let _guard = FlagGuard;
    assert_eq!(DEFAULT_NUM_BEST_HAPLOTYPES_PER_GRAPH, 128);
    assert_eq!(
        HcGenotypingConfig::strict_java().stand_emit_confidence,
        DEFAULT_STAND_CALL_CONF
    );
    let root = repo_root();
    let java = load_java(&root.join(JAVA_TSV_REL));
    let ref_fasta = root.join(REF_REL);
    let bam = root.join(BAM_REL);
    let dict = SequenceDictionary::from_fasta_path(&ref_fasta).expect("dict");
    let specs = parse_intervals_cli_string(&dict, INTERVAL).expect("interval");
    begin_hap_list_observe();
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
    let raw = region
        .reads
        .iter()
        .find(|r| r.flags() == FLAGS && std::str::from_utf8(r.qname()).ok() == Some(QNAME))
        .expect("raw bam read");
    let raw_bases = raw.seq().as_bytes();
    set_forensic_6r289_target(Some(QNAME), FLAGS);
    set_forensic_6r289_substitute(None, None);
    let _outcome = run_region(region, &dict, &ref_fasta);
    let snaps = take_forensic_6r289_snaps();
    let trim = take_hap_list_trim_span().expect("trim");
    assert_eq!((trim.trim_start, trim.trim_end), RUST_CLIP);
    let snap = snaps.last().expect("rust pairhmm input snap");
    kv("qname", QNAME);
    kv("flags", FLAGS.to_string());
    kv("n_rust_kernels", snaps.len().to_string());
    kv("raw_bam_len", raw_bases.len().to_string());
    kv("raw_bam_fnv", fnv_hex(&raw_bases));
    kv(
        "java_orig_vs_raw",
        if java.orig_bases == raw_bases {
            "IDENTICAL".into()
        } else {
            format!(
                "len_java={} len_raw={}",
                java.orig_bases.len(),
                raw_bases.len()
            )
        },
    );
    report_pair(&java, snap);
    let field = first_field(&java, snap);
    kv("classification", classification(field));
    kv("baseline_h0", fmt_f64(snap.h0_score));
    assert_eq!(snap.h0_score.to_bits(), RUST_H0.to_bits());
    if let Some(field) = field {
        let bytes = match field {
            Forensic6r289Field::ReadBases => java.proc_bases.clone(),
            Forensic6r289Field::BaseQuals => java.proc_bq.clone(),
            Forensic6r289Field::InsQuals => java.proc_iq.clone(),
            Forensic6r289Field::DelQuals => java.proc_dq.clone(),
            Forensic6r289Field::Gcp => java.proc_gcp.clone(),
            Forensic6r289Field::Hap0 => java.hap0.clone(),
        };
        let same_len = match field {
            Forensic6r289Field::Hap0 => true,
            Forensic6r289Field::ReadBases => bytes.len() == snap.base_quals.len(),
            _ => bytes.len() == snap.read_bases.len(),
        };
        kv("first_field", classification(Some(field)));
        if !same_len {
            kv(
                "counterfactual",
                "NOT_RUN length mismatch on the first field",
            );
        } else {
            set_forensic_6r289_substitute(Some(field), Some(bytes));
            let outcome2 = run_region(region, &dict, &ref_fasta);
            let snaps2 = take_forensic_6r289_snaps();
            let snap2 = snaps2.last().expect("cf snap");
            let recs = try_emit_call_region_variants(
                region,
                &outcome2,
                "SAMPLE",
                DEFAULT_STAND_EMIT_CONFIDENCE,
            )
            .unwrap_or_default();
            let pl = recs.iter().find(|r| {
                r.position == TARGET
                    && r.reference == TARGET_REF
                    && r.alternate.iter().any(|a| a == TARGET_ALT)
            });
            kv("counterfactual_h0", fmt_f64(snap2.h0_score));
            kv("java_h0", fmt_f64(JAVA_H0));
            kv(
                "h0_movement_toward_java",
                fmt_f64((snap.h0_score - JAVA_H0).abs() - (snap2.h0_score - JAVA_H0).abs()),
            );
            if let Some(rec) = pl {
                kv(
                    "counterfactual_integer_PL",
                    rec.samples[0]
                        .pl
                        .as_ref()
                        .map(|p| {
                            p.iter()
                                .map(|v| v.to_string())
                                .collect::<Vec<_>>()
                                .join(",")
                        })
                        .unwrap_or_else(|| "NONE".into()),
                );
            }
        }
    } else {
        kv("counterfactual", "NOT_RUN inputs identical");
    }
    kv("production_change", "NONE");
}

fn report_pair(java: &JavaBoundary, snap: &Forensic6r289Snap) {
    kv("rust_read_len", snap.read_bases.len().to_string());
    kv("java_proc_len", java.proc_len.to_string());
    kv("rust_hap_index", snap.hap0_index.to_string());
    kv("rust_hap_len", snap.hap0.len().to_string());
    kv("java_hap_len", java.hap0.len().to_string());
    kv("java_hap0_is_ref", java.hap0_is_ref.to_string());
    kv("java_hap0_align_start", java.hap0_align.to_string());
    kv("java_n_haps", java.n_haps.to_string());
    kv("rust_n_haps", snap.n_haps.to_string());
    kv("rust_mapq", snap.mapq.to_string());
    kv("java_proc_mapq", java.proc_mapq.to_string());
    kv("rust_pos", snap.pos_1based.to_string());
    kv("java_proc_start", java.proc_start.to_string());
    kv("java_orig_start", java.orig_start.to_string());
    kv("rust_cigar", &snap.cigar);
    kv("java_proc_cigar", &java.proc_cigar);
    kv("java_orig_cigar", &java.orig_cigar);
    kv("java_constant_gcp", java.constant_gcp.to_string());
    kv("rust_read_fnv", fnv_hex(&snap.read_bases));
    kv("java_proc_read_fnv", fnv_hex(&java.proc_bases));
    kv("java_orig_read_fnv", fnv_hex(&java.orig_bases));
    kv("rust_hap0_fnv", fnv_hex(&snap.hap0));
    kv("java_hap0_fnv", fnv_hex(&java.hap0));
    for (name, j, r) in [
        (
            "read_bases",
            java.proc_bases.as_slice(),
            snap.read_bases.as_slice(),
        ),
        (
            "base_qualities",
            java.proc_bq.as_slice(),
            snap.base_quals.as_slice(),
        ),
        (
            "insertion_qualities",
            java.proc_iq.as_slice(),
            snap.ins_quals.as_slice(),
        ),
        (
            "deletion_qualities",
            java.proc_dq.as_slice(),
            snap.del_quals.as_slice(),
        ),
        (
            "gap_continuation",
            java.proc_gcp.as_slice(),
            snap.gcp.as_slice(),
        ),
        ("haplotype_0", java.hap0.as_slice(), snap.hap0.as_slice()),
    ] {
        if j == r {
            kv(
                name,
                format!("IDENTICAL len={} fnv={}", j.len(), fnv_hex(j)),
            );
        } else if let Some((i, jv, rv)) = first_byte_diff(j, r) {
            kv(
                name,
                format!(
                    "DIFF len_java={} len_rust={} first_index={i} java={jv} rust={rv} java_fnv={} rust_fnv={}",
                    j.len(),
                    r.len(),
                    fnv_hex(j),
                    fnv_hex(r)
                ),
            );
        }
    }
    kv(
        "rust_read_bases",
        String::from_utf8_lossy(&snap.read_bases).to_string(),
    );
    kv(
        "java_proc_bases",
        String::from_utf8_lossy(&java.proc_bases).to_string(),
    );
    kv("rust_hap0", String::from_utf8_lossy(&snap.hap0).to_string());
    kv("java_hap0", String::from_utf8_lossy(&java.hap0).to_string());
}

fn first_field(java: &JavaBoundary, snap: &Forensic6r289Snap) -> Option<Forensic6r289Field> {
    if java.proc_bases.as_slice() != snap.read_bases.as_slice() {
        return Some(Forensic6r289Field::ReadBases);
    }
    if java.proc_bq.as_slice() != snap.base_quals.as_slice() {
        return Some(Forensic6r289Field::BaseQuals);
    }
    if java.proc_iq.as_slice() != snap.ins_quals.as_slice() {
        return Some(Forensic6r289Field::InsQuals);
    }
    if java.proc_dq.as_slice() != snap.del_quals.as_slice() {
        return Some(Forensic6r289Field::DelQuals);
    }
    if java.proc_gcp.as_slice() != snap.gcp.as_slice() {
        return Some(Forensic6r289Field::Gcp);
    }
    if java.hap0.as_slice() != snap.hap0.as_slice() {
        return Some(Forensic6r289Field::Hap0);
    }
    None
}
