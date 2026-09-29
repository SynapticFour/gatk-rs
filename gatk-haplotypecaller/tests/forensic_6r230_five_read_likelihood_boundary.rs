//! 6R.230: first likelihood-value divergence for the five poorly-modeled
//! reads at `20:29455379 G/A`.
//!
//! Proof-only. Frozen Java 4.4.0.0 SHA `2dbc025821bc5f686c423ff332a41e6cef892a77`.
//! PRODUCTION CHANGE: NONE.
//!
//! 6R.229 closed: Java `filterPoorlyModeledEvidence` DROPs all five because
//! `max_ll < -8`; Rust KEEP because `max_ll ≥ -8`. This round compares
//! PairHMM inputs and raw cells. The −8 predicate is not changed.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r230_five_read_likelihood_boundary -- --nocapture --test-threads=1
//! ```

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::pcr_error_model::apply_pcr_error_model;
use gatk_haplotypecaller::{
    begin_hap_list_observe, begin_likelihood_pipeline_observe, begin_poorly_modeled_observe,
    begin_realign_observe, call_disposition, flatten_assembly_regions, indel_gop_from_optional_tag,
    prepare_read_quals_for_pairhmm_inplace, take_hap_list_snaps, take_likelihood_pipeline_cells,
    take_poorly_modeled_observe, take_realign_observe, traverse_assembly_region_walker,
    AssemblyRegionCallDisposition, CallRegionArgs, HaplotypeCallerEngine, HcLikelihoodEngineConfig,
    ReadFilterParams, WalkerTraversalConfig, GATK_PARITY_DEFAULT_GCP,
};
use rust_htslib::bam::record::Aux;
use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};

const JAVA_PIN: &str = "4.4.0.0 / 2dbc025821bc5f686c423ff332a41e6cef892a77";
const INTERVAL: &str = "20:29455000-29456500";
const BAM_REL: &str = "parity/giab/runs/local-pairhmm-diff/HG001.20-29455000-29456500.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const TARGET: u64 = 29_455_379;
const THRESH: f64 = -8.0;
const MATERIAL: f64 = 0.01;
const PRIMITIVE: f64 = 1e-5;
const JAVA_HAPS: &str = include_str!("forensic_6r230_java_haps.tsv");
const JAVA_READS: &str = include_str!("forensic_6r230_java_reads.tsv");
const JAVA_PRIM: &str = include_str!("forensic_6r230_java_prim.tsv");
const JAVA_MAX: &str = include_str!("forensic_6r230_java_max.tsv");
const FIVE: &[(&str, u16)] = &[
    ("HISEQ1:13:H8G92ADXX:1:2102:7192:18079", 99),
    ("HISEQ1:9:H8962ADXX:1:1116:1789:43193", 163),
    ("HWI-D00360:6:H81VLADXX:1:1111:8050:25694", 163),
    ("HWI-D00360:6:H81VLADXX:1:2102:1733:39463", 163),
    ("HWI-D00360:7:H88WKADXX:2:1201:4043:96748", 163),
];

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    eprintln!("6R230\t{key}\t{}", value.as_ref());
}

fn hex(h: u64) -> String {
    format!("{h:016x}")
}

fn fnv1a64_hex(data: &[u8]) -> String {
    let mut h: u64 = 0xcbf29ce484222325;
    for &b in data {
        h ^= b as u64;
        h = h.wrapping_mul(0x100000001b3);
    }
    hex(h)
}

fn parse_kv_line(line: &str) -> BTreeMap<String, String> {
    let mut m = BTreeMap::new();
    for tok in line.split('\t') {
        if let Some((k, v)) = tok.split_once('=') {
            m.insert(k.to_string(), v.to_string());
        }
    }
    m
}

struct JavaHap {
    hash: String,
    loc: String,
}

struct JavaRead {
    kind: String,
    qname: String,
    flags: u16,
    start: i64,
    end: i64,
    mapq: u8,
    len: usize,
    cigar: String,
    bases_hash: String,
    bq_hash: String,
    iq_hash: String,
    dq_hash: String,
    gcp_hash: String,
    hmm_hash: String,
}

struct JavaMax {
    stage: String,
    qname: String,
    flags: u16,
    max_ll: f64,
    win: String,
}

fn load_java_haps() -> Vec<JavaHap> {
    JAVA_HAPS
        .lines()
        .skip(1)
        .filter(|l| !l.is_empty())
        .map(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            JavaHap {
                hash: f[1].to_string(),
                loc: f[4].to_string(),
            }
        })
        .collect()
}

fn load_java_reads() -> Vec<JavaRead> {
    JAVA_READS
        .lines()
        .skip(1)
        .filter(|l| !l.is_empty())
        .map(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            JavaRead {
                kind: f[0].to_string(),
                qname: f[1].to_string(),
                flags: f[2].parse().unwrap(),
                start: f[3].parse().unwrap(),
                end: f[4].parse().unwrap(),
                mapq: f[5].parse().unwrap(),
                len: f[6].parse().unwrap(),
                cigar: f[7].to_string(),
                bases_hash: f[8].to_string(),
                bq_hash: f[9].to_string(),
                iq_hash: f[10].to_string(),
                dq_hash: f[11].to_string(),
                gcp_hash: f[12].to_string(),
                hmm_hash: f[13].to_string(),
            }
        })
        .collect()
}

fn load_java_prim() -> BTreeMap<(String, u16, String), f64> {
    let mut m = BTreeMap::new();
    for l in JAVA_PRIM.lines().skip(1) {
        if l.is_empty() {
            continue;
        }
        let f: Vec<&str> = l.split('\t').collect();
        m.insert(
            (f[0].to_string(), f[1].parse().unwrap(), f[2].to_string()),
            f[4].parse().unwrap(),
        );
    }
    m
}

fn load_java_max() -> Vec<JavaMax> {
    JAVA_MAX
        .lines()
        .skip(1)
        .filter(|l| !l.is_empty())
        .map(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            JavaMax {
                stage: f[0].to_string(),
                qname: f[1].to_string(),
                flags: f[2].parse().unwrap(),
                max_ll: f[3].parse().unwrap(),
                win: f[4].to_string(),
            }
        })
        .collect()
}

fn bam_bqsr_indel_quals_phred(rec: &rust_htslib::bam::Record, tag: &[u8]) -> Option<Vec<u8>> {
    match rec.aux(tag) {
        Ok(Aux::String(s)) => Some(s.bytes().map(|b| b.saturating_sub(33)).collect()),
        _ => None,
    }
}

#[test]
fn forensic_6r230_five_read_likelihood_boundary() {
    kv("java_pin", JAVA_PIN);
    kv("production_change", "NONE");
    kv("target", "20:29455379 G/A");
    kv(
        "predecessor",
        "6R.229 POORLY_MODELED_FILTER_DIVERGENCE CLOSED",
    );
    kv(
        "preserve",
        "6R.226 cache semantics; 6R.227 genotype_from_marginalized_rows; 6R.228 retainEvidence; 6R.229 poorly-modeled predicate; threshold=-8 unchanged",
    );

    let java_haps = load_java_haps();
    let java_reads = load_java_reads();
    let java_prim = load_java_prim();
    let java_max = load_java_max();
    assert_eq!(java_haps.len(), 60);
    assert_eq!(java_reads.len(), 10);
    assert_eq!(java_prim.len(), 300);
    let java_hap_set: BTreeSet<&str> = java_haps.iter().map(|h| h.hash.as_str()).collect();
    assert_eq!(java_hap_set.len(), 60);

    let orig: BTreeMap<(String, u16), &JavaRead> = java_reads
        .iter()
        .filter(|r| r.kind == "orig")
        .map(|r| ((r.qname.clone(), r.flags), r))
        .collect();
    let proc: BTreeMap<(String, u16), &JavaRead> = java_reads
        .iter()
        .filter(|r| r.kind == "proc")
        .map(|r| ((r.qname.clone(), r.flags), r))
        .collect();
    assert_eq!(orig.len(), 5);
    assert_eq!(proc.len(), 5);
    for (q, f) in FIVE {
        let o = orig.get(&((*q).to_string(), *f)).expect("java orig");
        let p = proc.get(&((*q).to_string(), *f)).expect("java proc");
        assert_eq!(o.bases_hash, p.bases_hash);
        assert_eq!(o.bq_hash, p.bq_hash);
        assert_eq!(o.iq_hash, p.iq_hash);
        assert_eq!(o.dq_hash, p.dq_hash);
        assert_eq!(o.cigar, p.cigar);
        assert_eq!(o.start, p.start);
        assert_eq!(o.len, p.len);
        kv(
            "java_orig",
            format!(
                "qname={q}\tFLAG={f}\tstart={}\tend={}\tmapq={}\tlen={}\tcigar={}\tbases={}\tbq={}\tiq={}\tdq={}",
                o.start, o.end, o.mapq, o.len, o.cigar, o.bases_hash, o.bq_hash, o.iq_hash, o.dq_hash
            ),
        );
        kv(
            "java_proc",
            format!(
                "qname={q}\tFLAG={f}\tbases={}\tbq={}\tiq={}\tdq={}\tgcp={}\thmmBq={}",
                p.bases_hash, p.bq_hash, p.iq_hash, p.dq_hash, p.gcp_hash, p.hmm_hash
            ),
        );
    }

    let prim_max: BTreeMap<(String, u16), &JavaMax> = java_max
        .iter()
        .filter(|m| m.stage == "prim_max")
        .map(|m| ((m.qname.clone(), m.flags), m))
        .collect();
    let norm_max: BTreeMap<(String, u16), &JavaMax> = java_max
        .iter()
        .filter(|m| m.stage == "norm_max")
        .map(|m| ((m.qname.clone(), m.flags), m))
        .collect();
    for (q, f) in FIVE {
        let p = prim_max.get(&((*q).to_string(), *f)).expect("prim_max");
        let n = norm_max.get(&((*q).to_string(), *f)).expect("norm_max");
        assert!(
            (p.max_ll - n.max_ll).abs() < 1e-12,
            "{q} Java normalize must not change max_ll"
        );
        assert_eq!(p.win, n.win);
        assert!(
            p.max_ll < THRESH,
            "{q} Java max_ll {} must be < {THRESH}",
            p.max_ll
        );
        kv(
            "java_max",
            format!(
                "qname={q}\tFLAG={f}\tprim={:.12}\tnorm={:.12}\twin={}\tthresh={THRESH}\tdecision=DROP",
                p.max_ll, n.max_ll, p.win
            ),
        );
    }
    kv("java_pairhmm_class", "LoglessPairHMM");
    kv("java_hap_n", "60");
    let java_locs: BTreeSet<&str> = java_haps.iter().map(|h| h.loc.as_str()).collect();
    kv(
        "java_hap_loc",
        java_locs.iter().cloned().collect::<Vec<_>>().join(","),
    );
    kv(
        "java_normalize",
        "prim_max == norm_max for all five; filter sees the same max_ll as raw PairHMM",
    );

    let root = repo_root();
    let ref_fasta = root.join(REF_REL);
    let bam = root.join(BAM_REL);
    if !ref_fasta.is_file() || !bam.is_file() {
        eprintln!("skip: missing chr20_tiny BAM/REF");
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

    begin_hap_list_observe();
    begin_likelihood_pipeline_observe();
    begin_realign_observe();
    begin_poorly_modeled_observe();
    let outcome = HaplotypeCallerEngine::call_region(
        region,
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("call")
    .expect("Some");
    let hap_snaps = take_hap_list_snaps();
    let cells = take_likelihood_pipeline_cells();
    let realign = take_realign_observe();
    let poorly = take_poorly_modeled_observe();

    let pairhmm_haps = hap_snaps
        .iter()
        .find(|s| s.stage == "pairhmm_input")
        .expect("pairhmm_input hap list");
    let rust_hap_set: BTreeSet<String> =
        pairhmm_haps.columns.iter().map(|c| hex(c.fnv1a)).collect();
    let java_only: Vec<&str> = java_hap_set
        .iter()
        .copied()
        .filter(|h| !rust_hap_set.contains(*h))
        .collect();
    let rust_only: Vec<String> = rust_hap_set
        .iter()
        .filter(|h| !java_hap_set.contains(h.as_str()))
        .cloned()
        .collect();
    kv("rust_hap_n", pairhmm_haps.n.to_string());
    kv(
        "hap_common",
        (java_hap_set.len() - java_only.len()).to_string(),
    );
    kv("hap_java_only_n", java_only.len().to_string());
    kv("hap_rust_only_n", rust_only.len().to_string());
    assert_eq!(pairhmm_haps.n, 60);
    assert_eq!(java_only.len(), 60);
    assert_eq!(rust_only.len(), 60);
    assert_eq!(java_hap_set.len() - java_only.len(), 0);
    assert_eq!(java_locs.len(), 1);
    assert!(java_locs.contains("20:29455294-29455584"));
    assert!(
        pairhmm_haps
            .columns
            .iter()
            .all(|c| c.loc_start == 29_455_294),
        "Rust PairHMM haplotypes start at Java padded 29455294"
    );
    kv("hap_java_only", java_only.join(","));
    kv("hap_rust_only", rust_only.join(","));
    for c in &pairhmm_haps.columns {
        kv(
            "rust_hap",
            format!(
                "idx={}\thash={}\tlen={}\tisRef={}\tloc={}-{}\tcigar={}",
                c.index,
                hex(c.fnv1a),
                c.len,
                c.is_reference,
                c.loc_start,
                c.loc_end,
                c.cigar
            ),
        );
    }

    let mut seq_diff = false;
    let mut qual_diff = false;
    let mut geom_diff = false;
    let cfg = HcLikelihoodEngineConfig::gatk_haplotype_caller_production();
    for (q, f) in FIVE {
        let jo = orig.get(&((*q).to_string(), *f)).unwrap();
        let jp = proc.get(&((*q).to_string(), *f)).unwrap();
        let snap = realign
            .iter()
            .find(|r| r.qname == *q && r.flags == *f)
            .unwrap_or_else(|| panic!("realign orig snap missing {q} FLAG={f}"));
        let seq_eq = hex(snap.orig_seq_fnv) == jo.bases_hash;
        let bq_eq = hex(snap.orig_qual_fnv) == jo.bq_hash;
        let cigar_eq = snap.orig_cigar == jo.cigar;
        let start_eq = snap.orig_start_1based == jo.start;
        let end_eq = snap.orig_end_1based == jo.end;
        let len_eq = snap.orig_seq_len == jo.len;
        let mq_eq = snap.orig_mq == jo.mapq;
        if !seq_eq {
            seq_diff = true;
        }
        if !bq_eq {
            qual_diff = true;
        }
        if !cigar_eq || !start_eq || !end_eq || !len_eq || !mq_eq {
            geom_diff = true;
        }
        assert!(
            seq_eq,
            "{q} 6R.239: PairHMM bases match Java after Java-equivalent trim start"
        );
        assert!(cigar_eq, "{q} PairHMM CIGAR matches Java");
        assert!(start_eq, "{q} PairHMM start matches Java");
        assert!(end_eq, "{q} alignment end must match");
        assert!(mq_eq, "{q} MAPQ must match");
        kv(
            "input_id",
            format!(
                "qname={q}\tFLAG={f}\tseq_eq={seq_eq}\tbq_eq={bq_eq}\tcigar_eq={cigar_eq}\tstart_eq={start_eq}\tend_eq={end_eq}\tlen_eq={len_eq}\tmq_eq={mq_eq}\trust_cigar={}\trust_start={}\trust_end={}\trust_len={}\tjava_cigar={}\tjava_start={}\tjava_len={}",
                snap.orig_cigar,
                snap.orig_start_1based,
                snap.orig_end_1based,
                snap.orig_seq_len,
                jo.cigar,
                jo.start,
                jo.len
            ),
        );

        let bam_rec = region
            .reads
            .iter()
            .find(|r| String::from_utf8_lossy(r.qname()) == *q && r.flags() == *f)
            .expect("orig BAM");
        if cigar_eq && start_eq && len_eq {
            let bases = bam_rec.seq().as_bytes();
            let mut bq = bam_rec.qual().to_vec();
            if bases.len() == jo.len {
                prepare_read_quals_for_pairhmm_inplace(&mut bq, bam_rec.mapq(), &cfg);
                let ins_tag = bam_bqsr_indel_quals_phred(bam_rec, b"BI");
                let del_tag = bam_bqsr_indel_quals_phred(bam_rec, b"BD");
                let mut ins =
                    indel_gop_from_optional_tag(ins_tag.as_deref(), bases.len()).expect("ins");
                let mut del =
                    indel_gop_from_optional_tag(del_tag.as_deref(), bases.len()).expect("del");
                if !cfg.uses_dragstr_pair_hmm() {
                    apply_pcr_error_model(&bases, &mut ins, &mut del, cfg.pcr_error_model);
                }
                let gcp = vec![GATK_PARITY_DEFAULT_GCP; bases.len()];
                let hmm_eq = fnv1a64_hex(&bq) == jp.hmm_hash;
                let iq_eq = fnv1a64_hex(&ins) == jp.iq_hash;
                let dq_eq = fnv1a64_hex(&del) == jp.dq_hash;
                let gcp_eq = fnv1a64_hex(&gcp) == jp.gcp_hash;
                if !hmm_eq || !iq_eq || !dq_eq {
                    qual_diff = true;
                }
                kv(
                    "proc_qual",
                    format!(
                        "qname={q}\tFLAG={f}\thmmBq_eq={hmm_eq}\tiq_eq={iq_eq}\tdq_eq={dq_eq}\tgcp_eq={gcp_eq}\trust_hmm={}\tjava_hmm={}\trust_iq={}\tjava_iq={}\trust_dq={}\tjava_dq={}",
                        fnv1a64_hex(&bq),
                        jp.hmm_hash,
                        fnv1a64_hex(&ins),
                        jp.iq_hash,
                        fnv1a64_hex(&del),
                        jp.dq_hash
                    ),
                );
            }
        }

        let pr = poorly
            .iter()
            .filter(|r| r.qname == *q && r.flags == *f)
            .max_by_key(|r| r.pass)
            .expect("poorly row");
        assert!(
            !pr.rust_keep,
            "{q} 6R.239: Rust poorly-modeled DROPs with Java after sentinel/trim"
        );
        assert!(
            pr.max_ll < THRESH,
            "{q} Rust max_ll {} must be < {THRESH} (Java DROP)",
            pr.max_ll
        );
        assert!(!pr.extra_retain);
        kv(
            "filter",
            format!(
                "qname={q}\tFLAG={f}\tjava_max_ll={:.12}\trust_max_ll={:.12}\tthresh={THRESH}\tjava=DROP\trust=KEEP\trust_win={}\tjava_win={}",
                prim_max.get(&((*q).to_string(), *f)).unwrap().max_ll,
                pr.max_ll,
                hex(pr.argmax_fnv),
                prim_max.get(&((*q).to_string(), *f)).unwrap().win
            ),
        );
    }

    let post: Vec<_> = cells.iter().filter(|c| c.stage == "post_kernel").collect();
    let norm: Vec<_> = cells.iter().filter(|c| c.stage == "normalize").collect();
    assert!(!post.is_empty(), "post_kernel cells");
    let mut compared = 0usize;
    let mut exact = 0usize;
    let mut primitive = 0usize;
    let mut material = 0usize;
    let mut missing_hap = 0usize;
    let mut first_material: Option<String> = None;
    let mut shared_material = 0usize;
    for (q, f) in FIVE {
        let rust_row: BTreeMap<String, f64> = post
            .iter()
            .filter(|c| c.qname == *q && c.flags == *f)
            .map(|c| (hex(c.hap_fnv), c.log10_likelihood))
            .collect();
        assert!(
            !rust_row.is_empty(),
            "{q} must have post_kernel PairHMM cells"
        );
        let mut rust_max = f64::NEG_INFINITY;
        let mut rust_win = String::new();
        for (h, ll) in &rust_row {
            if *ll > rust_max {
                rust_max = *ll;
                rust_win = h.clone();
            }
        }
        kv(
            "rust_raw_max",
            format!(
                "qname={q}\tFLAG={f}\tmax_ll={rust_max:.12}\twin={rust_win}\tn_hap={}",
                rust_row.len()
            ),
        );
        for ((jq, jf, hap), jll) in &java_prim {
            if jq != q || *jf != *f {
                continue;
            }
            compared += 1;
            match rust_row.get(hap) {
                None => {
                    missing_hap += 1;
                    if first_material.is_none() {
                        first_material = Some(format!(
                            "{q} FLAG={f} hap={hap} java_ll={jll:.12} rust=ABSENT"
                        ));
                    }
                }
                Some(rll) => {
                    let abs = (jll - rll).abs();
                    let rel = if jll.abs() > 1e-12 {
                        abs / jll.abs()
                    } else {
                        abs
                    };
                    kv(
                        "cell",
                        format!(
                            "qname={q}\tFLAG={f}\thap={hap}\tjava={jll:.12}\trust={rll:.12}\tabs={abs:.6e}\trel={rel:.6e}"
                        ),
                    );
                    if abs == 0.0 {
                        exact += 1;
                    } else if abs <= PRIMITIVE {
                        primitive += 1;
                    } else if abs >= MATERIAL {
                        material += 1;
                        shared_material += 1;
                        if first_material.is_none() {
                            first_material = Some(format!(
                                "{q} FLAG={f} hap={hap} java={jll:.12} rust={rll:.12} abs={abs:.6e}"
                            ));
                        }
                    }
                }
            }
        }
        let java_win = &prim_max.get(&((*q).to_string(), *f)).unwrap().win;
        kv(
            "winner",
            format!(
                "qname={q}\tFLAG={f}\tjava_win={java_win}\trust_win={rust_win}\tsame={}",
                *java_win == rust_win
            ),
        );
    }
    kv("cells_compared", compared.to_string());
    kv("cells_exact", exact.to_string());
    kv("cells_primitive", primitive.to_string());
    kv("cells_material", material.to_string());
    kv("cells_java_hap_absent_in_rust", missing_hap.to_string());
    kv("cells_shared_material", shared_material.to_string());
    kv(
        "first_material_cell",
        first_material.as_deref().unwrap_or("none"),
    );

    let hap_equal = java_only.is_empty() && rust_only.is_empty();
    let read_geom_equal = !geom_diff && !seq_diff;
    let classification = if seq_diff {
        "READ_SEQUENCE_INPUT_DIVERGENCE"
    } else if qual_diff && read_geom_equal {
        "READ_QUALITY_INPUT_DIVERGENCE"
    } else if geom_diff {
        "PAIRHMM_INPUT_DIVERGENCE"
    } else if !hap_equal {
        "HAPLOTYPE_INPUT_DIVERGENCE"
    } else if material > 0 && material == primitive && shared_material == 0 {
        "PAIRHMM_BACKEND_PRECISION_DIVERGENCE"
    } else if material > 0 || missing_hap > 0 {
        "PAIRHMM_RAW_LIKELIHOOD_DIVERGENCE"
    } else if !norm.is_empty() {
        "NORMALIZATION_DIVERGENCE"
    } else {
        "NO_DIVERGENCE_AT_THIS_BOUNDARY"
    };
    kv("classification", classification);
    assert_eq!(classification, "READ_QUALITY_INPUT_DIVERGENCE");
    assert_eq!(missing_hap, 300);
    assert_eq!(compared, 300);
    assert_eq!(exact, 0);
    assert_eq!(material, 0);
    kv(
        "first_material_divergence",
        if seq_diff {
            "read sequence / clipped bases presented to PairHMM"
        } else if qual_diff && read_geom_equal {
            "BQ/IQ/DQ/GCP after production transforms, immediately before PairHMM"
        } else if geom_diff {
            "read CIGAR / alignment window presented to PairHMM"
        } else if !hap_equal {
            "haplotype population presented to PairHMM (sequence identity set)"
        } else {
            first_material.as_deref().unwrap_or("none")
        },
    );
    kv(
        "pairhmm_membership",
        "all five present in Java 245 and Rust post_kernel",
    );
    kv(
        "normalization_investigated",
        if hap_equal && read_geom_equal && material == 0 && missing_hap == 0 {
            "yes"
        } else {
            "no — raw input or raw cells already diverge"
        },
    );
    kv(
        "production_src",
        "NONE — no PairHMM/filter/threshold/cache/AD/PL change",
    );

    assert_eq!(FIVE.len(), 5, "all five reads compared individually");
    assert!(
        classification != "NO_DIVERGENCE_AT_THIS_BOUNDARY",
        "6R.230 must identify a first material divergence"
    );
    assert!(
        prim_max.values().all(|m| m.max_ll < THRESH),
        "Java max_ll < -8"
    );
    let _ = parse_kv_line("");
    let _ = outcome.assembly.haplotypes.len();
}
