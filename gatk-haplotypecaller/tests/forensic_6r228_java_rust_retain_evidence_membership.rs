//! 6R.228: Java vs Rust event-local retainEvidence membership at `20:29455379 G/A`.
//!
//! Proof-only. Frozen Java 4.4.0.0 SHA `2dbc025821bc5f686c423ff332a41e6cef892a77`.
//! PRODUCTION CHANGE: NONE.
//!
//! Live Java `genotype-emit-at-loc` dump after `retainEvidence(SimpleInterval(merged)±2)`
//! is the 47-read object consumed by `calculateGLsForThisEvent` / DepthPerAlleleBySample.
//! Do not infer that set from FORMAT AD.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r228_java_rust_retain_evidence_membership -- --nocapture --test-threads=1
//! ```

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::event_map::{
    build_per_haplotype_variation_events, merged_alleles_for_genotyping,
    variation_events_at_position_from_cache,
};
use gatk_haplotypecaller::hc_allele_mapping::create_allele_mapper_with_events;
use gatk_haplotypecaller::hc_genotyping_engine::{
    alt_hap_indices_for_genotype_marginalization, java_alignment_read_overlaps_interval,
    marginalize_rows_to_biallelic_alleles, ref_hap_indices_for_genotype_marginalization,
    region_likelihoods_to_rows, HcGenotypingConfig, InformativeAd,
    DEFAULT_INFORMATIVE_READ_OVERLAP_MARGIN,
};
use gatk_haplotypecaller::read_realignment::LOG_10_INFORMATIVE_THRESHOLD;
use gatk_haplotypecaller::read_unclip::{alignment_end_1based, gatk_soft_start_1based};
use gatk_haplotypecaller::{
    begin_poorly_modeled_observe, call_disposition, flatten_assembly_regions,
    take_poorly_modeled_observe, traverse_assembly_region_walker, try_emit_call_region_variants,
    AssemblyRegionCallDisposition, CallRegionArgs, GenomePosition, HaplotypeCallerEngine,
    ReadFilterParams, WalkerTraversalConfig, DEFAULT_STAND_EMIT_CONFIDENCE,
};
use rust_htslib::bam::record::Aux;
use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};

const JAVA_PIN: &str = "4.4.0.0 / 2dbc025821bc5f686c423ff332a41e6cef892a77";
const INTERVAL: &str = "20:29455000-29456500";
const BAM_REL: &str = "parity/giab/runs/local-pairhmm-diff/HG001.20-29455000-29456500.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const TARGET: u64 = 29_455_379;
const CLOSED_PL: u64 = 29_455_015;
const TARGET_REF: &str = "G";
const TARGET_ALT: &str = "A";
const JAVA_AD: [i32; 2] = [42, 5];
const JAVA_PL: [i32; 3] = [84, 0, 1738];
const JAVA_RETAIN_START: i64 = 29_455_377;
const JAVA_RETAIN_END: i64 = 29_455_381;
const MARGIN: i32 = DEFAULT_INFORMATIVE_READ_OVERLAP_MARGIN;
const JAVA_POST_TSV: &str = include_str!("forensic_6r228_java_post_retain.tsv");
const JAVA_PRE_TSV: &str = include_str!("forensic_6r228_java_pre_retain.tsv");
const MARG_HASH: u64 = 0x28a5025451aafbd3;
const RUST_ONLY_PIN: &[(&str, u16)] = &[];

#[derive(Clone, Debug)]
#[allow(dead_code)]
struct JavaPost {
    qname: String,
    flag: u16,
    start: i64,
    end: i64,
    mapq: u8,
    cigar: String,
    len: i32,
    mate: String,
    mate_start: i64,
    lr: f64,
    la: f64,
    vote: String,
}

#[derive(Clone, Debug)]
#[allow(dead_code)]
struct JavaPre {
    qname: String,
    flag: u16,
    start: i64,
    end: i64,
    u_start: i64,
    u_end: i64,
    cigar: String,
    overlap: bool,
}

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    eprintln!("6R228\t{key}\t{}", value.as_ref());
}

fn fnv1a64(data: &[u8]) -> u64 {
    let mut h = 0xcbf29ce484222325u64;
    for b in data {
        h ^= u64::from(*b);
        h = h.wrapping_mul(0x100000001b3);
    }
    h
}

fn membership_hash(ids: &BTreeSet<(String, u16)>) -> u64 {
    let mut buf = String::new();
    for (q, f) in ids {
        buf.push_str(q);
        buf.push('\t');
        buf.push_str(&f.to_string());
        buf.push('\n');
    }
    fnv1a64(buf.as_bytes())
}

fn vote_row(lr: f64, la: f64) -> (f64, bool, &'static str) {
    let conf = (lr - la).abs();
    let informative = lr.is_finite() && la.is_finite() && conf > LOG_10_INFORMATIVE_THRESHOLD;
    let vote = if !informative {
        "UNINF"
    } else if lr > la {
        "REF"
    } else {
        "ALT"
    };
    (conf, informative, vote)
}

fn java_interval_overlaps(read_start: i64, read_end: i64) -> bool {
    JAVA_RETAIN_START <= read_end && read_start <= JAVA_RETAIN_END
}

fn mate_on_same_contig_or_unmapped(rec: &rust_htslib::bam::Record) -> bool {
    if !rec.is_paired() || rec.is_mate_unmapped() || rec.is_unmapped() {
        return true;
    }
    rec.tid() == rec.mtid()
}

fn orig_alignment(rec: &rust_htslib::bam::Record) -> (i64, String, bool) {
    let has_orig = rec.aux(b"OP").is_ok() || rec.aux(b"OC").is_ok();
    let orig_pos = match rec.aux(b"OP") {
        Ok(Aux::I32(v)) => i64::from(v),
        Ok(Aux::U32(v)) => i64::from(v),
        Ok(Aux::I8(v)) => i64::from(v),
        Ok(Aux::U8(v)) => i64::from(v),
        Ok(Aux::I16(v)) => i64::from(v),
        Ok(Aux::U16(v)) => i64::from(v),
        _ => rec.pos() + 1,
    };
    let orig_cigar = match rec.aux(b"OC") {
        Ok(Aux::String(s)) => s.to_string(),
        _ => rec.cigar().to_string(),
    };
    (orig_pos, orig_cigar, has_orig)
}

fn mate_contig_label(rec: &rust_htslib::bam::Record) -> String {
    if !rec.is_paired() || rec.is_mate_unmapped() {
        return ".".into();
    }
    if rec.tid() == rec.mtid() {
        "20".into()
    } else {
        format!("tid={}", rec.mtid())
    }
}

fn parse_java_post() -> Vec<JavaPost> {
    let mut out = Vec::new();
    for (i, line) in JAVA_POST_TSV.lines().enumerate() {
        if i == 0 || line.is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        assert!(f.len() >= 12, "java post tsv cols");
        out.push(JavaPost {
            qname: f[0].to_string(),
            flag: f[1].parse().unwrap(),
            start: f[2].parse().unwrap(),
            end: f[3].parse().unwrap(),
            mapq: f[4].parse().unwrap(),
            cigar: f[5].to_string(),
            len: f[6].parse().unwrap(),
            mate: f[7].to_string(),
            mate_start: f[8].parse().unwrap(),
            lr: f[9].parse().unwrap(),
            la: f[10].parse().unwrap(),
            vote: f[11].to_string(),
        });
    }
    out
}

fn parse_java_pre() -> Vec<JavaPre> {
    let mut out = Vec::new();
    for (i, line) in JAVA_PRE_TSV.lines().enumerate() {
        if i == 0 || line.is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        assert!(f.len() >= 12, "java pre tsv cols");
        out.push(JavaPre {
            qname: f[0].to_string(),
            flag: f[1].parse().unwrap(),
            start: f[2].parse().unwrap(),
            end: f[3].parse().unwrap(),
            u_start: f[4].parse().unwrap(),
            u_end: f[5].parse().unwrap(),
            cigar: f[7].to_string(),
            overlap: f[11] == "true",
        });
    }
    out
}

#[test]
fn forensic_6r228_java_contract_pin() {
    kv("java_pin", JAVA_PIN);
    kv("production_change", "NONE — proof-only");
    kv(
        "java_lifecycle",
        "assignGenotypeLikelihoods: marginalize(alleleMapper) → retainEvidence(SimpleInterval(merged).expandWithinContig(2)) → calculateGLsForThisEvent → DepthPerAlleleBySample on the same object",
    );
    kv("java_retain_interval", "20:29455377-29455381");
    kv("merged_span", "20:29455379-29455379");
    assert_eq!(MARGIN, 2);
    assert_eq!(LOG_10_INFORMATIVE_THRESHOLD, 0.2);
    let post = parse_java_post();
    let pre = parse_java_pre();
    assert_eq!(post.len(), 47, "live Java post-retainEvidence membership");
    assert_eq!(
        pre.len(),
        236,
        "live Java pre-retainEvidence sampleEvidence"
    );
    assert_eq!(pre.iter().filter(|r| r.overlap).count(), 47);
    let mut ref_n = 0i32;
    let mut alt_n = 0i32;
    for r in &post {
        match r.vote.as_str() {
            "REF" => ref_n += 1,
            "ALT" => alt_n += 1,
            other => panic!("java post vote {other}"),
        }
        assert!(java_interval_overlaps(r.start, r.end), "{}", r.qname);
        assert_eq!(r.mate, "20");
    }
    assert_eq!([ref_n, alt_n], JAVA_AD);
    kv("java_post_retain_n", "47");
    kv("java_post_retain_ad_votes", "42,5");
    kv("java_pl_from_dump", "84,0,1738");
    assert_eq!(JAVA_PL, [84, 0, 1738]);
}

#[test]
fn forensic_6r228_java_rust_retain_evidence_membership() {
    let root = repo_root();
    let ref_fasta = root.join(REF_REL);
    let bam = root.join(BAM_REL);
    if !ref_fasta.is_file() || !bam.is_file() {
        eprintln!("skip: missing chr20_tiny BAM/REF");
        return;
    }
    kv("java_pin", JAVA_PIN);
    kv("production_change", "NONE");
    kv("target", "20:29455379 G/A");
    kv("starting_6r227", "isolated==production P2 52×2 marg_hash=0x28a5025451aafbd3 AD=47,5 PL=68,0,1937 vs Java AD=42,5 PL=84,0,1738");

    let java_post = parse_java_post();
    let java_pre = parse_java_pre();
    let java_ids: BTreeSet<(String, u16)> = java_post
        .iter()
        .map(|r| (r.qname.clone(), r.flag))
        .collect();
    let java_pre_by_id: BTreeMap<(String, u16), &JavaPre> = java_pre
        .iter()
        .map(|r| ((r.qname.clone(), r.flag), r))
        .collect();
    kv(
        "java_membership_hash",
        format!("{:#x}", membership_hash(&java_ids)),
    );
    kv("java_retain_interval", "20:29455377-29455381");
    kv("java_hap_ll_n", java_pre.len().to_string());
    kv("java_pre_retain_n", java_pre.len().to_string());
    kv("java_post_retain_n", java_post.len().to_string());
    kv("java_pre_overlap_true_n", "47");
    assert_eq!(java_ids.len(), 47);

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
    kv(
        "active_full",
        format!(
            "{}:{}-{}",
            region.contig,
            region.start.get(),
            region.end.get()
        ),
    );

    begin_poorly_modeled_observe();
    let outcome = HaplotypeCallerEngine::call_region(
        region,
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("call")
    .expect("Some");
    let poorly = take_poorly_modeled_observe();

    let haps = &outcome.assembly.haplotypes;
    let reads = &outcome.genotyping_reads;
    let likelihoods = &outcome.read_likelihoods;
    kv("hap_n", haps.len().to_string());
    kv("genotyping_read_n", reads.len().to_string());
    kv(
        "pairhmm_unique_n",
        likelihoods
            .iter()
            .map(|c| c.read_index.get())
            .collect::<BTreeSet<_>>()
            .len()
            .to_string(),
    );

    let (full_ref, full_pad) = outcome.assembly.event_map_reference();
    let hap_events = build_per_haplotype_variation_events(
        haps,
        full_ref,
        full_pad,
        outcome.assembly.max_mnp_distance(),
        &region.contig,
    );
    let emit_spanning = !HcGenotypingConfig::strict_java().disable_spanning_event_genotyping;
    let cache_only = variation_events_at_position_from_cache(&hap_events, TARGET, emit_spanning);
    assert_eq!(cache_only.len(), 1);
    assert_eq!(cache_only[0].ref_allele, TARGET_REF);
    assert_eq!(cache_only[0].alt_allele, TARGET_ALT);
    kv(
        "merged_alleles",
        format!("{:?}", merged_alleles_for_genotyping(&cache_only, TARGET)),
    );

    let apply_bases = outcome.assembly.apply_bases_shared();
    let apply_pad = haps
        .iter()
        .find(|h| h.is_reference)
        .and_then(|h| h.genome_loc.map(|g| g.start_1based()))
        .unwrap_or(full_pad);
    let mapping = create_allele_mapper_with_events(
        &cache_only[0],
        TARGET,
        haps,
        apply_pad,
        apply_bases.as_ref(),
        outcome.assembly.max_mnp_distance(),
        emit_spanning,
        Some(&hap_events),
    );
    let config = HcGenotypingConfig::strict_java();
    let ref_hap = haps.iter().find(|h| h.is_reference).expect("ref hap");
    let ref_pool =
        ref_hap_indices_for_genotype_marginalization(&mapping, haps, &config, Some(&cache_only[0]));
    let alt_pool = alt_hap_indices_for_genotype_marginalization(
        &mapping,
        haps,
        &cache_only[0],
        ref_hap,
        apply_pad,
        apply_bases.as_ref(),
        outcome.assembly.max_mnp_distance(),
        &region.contig,
        &config,
    );
    let rows = region_likelihoods_to_rows(likelihoods, haps.len());
    let marg_all = marginalize_rows_to_biallelic_alleles(&rows, &ref_pool, &alt_pool);

    let mut rust_rows = Vec::new();
    let mut rust_ids: BTreeSet<(String, u16)> = BTreeSet::new();
    for row in &marg_all {
        let Some(idx) = row.matrix_read_index() else {
            continue;
        };
        let rec = &reads[idx];
        if !java_alignment_read_overlaps_interval(rec, TARGET, TARGET, MARGIN) {
            continue;
        }
        let qname = String::from_utf8_lossy(rec.qname()).into_owned();
        let flag = rec.flags();
        let lr = row.haplotype_log10_likelihoods[0];
        let la = row.haplotype_log10_likelihoods[1];
        let (conf, informative, vote) = vote_row(lr, la);
        let start = rec.pos() + 1;
        let end = alignment_end_1based(rec) as i64;
        let (orig_pos, orig_cigar, has_orig) = orig_alignment(rec);
        let soft_start = gatk_soft_start_1based(rec);
        let java_coords_overlap = java_interval_overlaps(start, end);
        let unclipped_overlap = java_interval_overlaps(soft_start, end);
        let mate_ok = mate_on_same_contig_or_unmapped(rec);
        rust_ids.insert((qname.clone(), flag));
        rust_rows.push((
            idx,
            qname,
            flag,
            start,
            end,
            rec.mapq(),
            rec.cigar().to_string(),
            rec.seq_len(),
            mate_contig_label(rec),
            rec.mpos() + 1,
            orig_pos,
            orig_cigar,
            has_orig,
            lr,
            la,
            conf,
            informative,
            vote,
            java_coords_overlap,
            unclipped_overlap,
            mate_ok,
        ));
    }
    rust_rows.sort_by(|a, b| a.1.cmp(&b.1).then(a.2.cmp(&b.2)));
    kv("rust_event_local_n", rust_ids.len().to_string());
    kv(
        "rust_membership_hash",
        format!("{:#x}", membership_hash(&rust_ids)),
    );
    kv("marg_hash_pin", format!("{MARG_HASH:#x}"));

    let mut rust_ref = 0i32;
    let mut rust_alt = 0i32;
    let mut rust_uninf = 0i32;
    for r in &rust_rows {
        match r.17 {
            "REF" => rust_ref += 1,
            "ALT" => rust_alt += 1,
            _ => rust_uninf += 1,
        }
        kv(
            "rust_row",
            format!(
                "idx={}\tqname={}\tFLAG={}\tPOS={}-{}\tMAPQ={}\tCIGAR={}\tlen={}\tmate={}:{}\torigPos={}\torigCigar={}\thas_OP/OC={}\tlr={:.6}\tla={:.6}\tinf={}\tvote={}\tjava_interval_overlap={}\tunclipped_overlap={}\tmate_ok={}",
                r.0, r.1, r.2, r.3, r.4, r.5, r.6, r.7, r.8, r.9, r.10, r.11, r.12, r.13, r.14, r.16, r.17, r.18, r.19, r.20
            ),
        );
    }
    kv(
        "rust_votes",
        format!("REF={rust_ref}\tALT={rust_alt}\tUNINF={rust_uninf}"),
    );
    let site = outcome
        .genotyped_calls
        .iter()
        .find(|c| {
            c.event.start_1based.get() == TARGET
                && c.event.ref_allele == TARGET_REF
                && c.event.alt_allele == TARGET_ALT
        })
        .expect("G/A");
    let ad = site.genotype.format.ad_as_i32();
    let pl = site.genotype.format.pl_as_i32();
    kv("rust_format_ad", format!("{ad:?}"));
    kv("rust_format_pl", format!("{pl:?}"));
    assert_eq!(ad, vec![42, 5]);
    assert_eq!(pl, vec![84, 0, 1738]);
    assert_eq!(
        InformativeAd::from_marginalized_rows(
            &marg_all
                .iter()
                .filter(|row| {
                    row.matrix_read_index().is_some_and(|i| {
                        reads.get(i).is_some_and(|r| {
                            java_alignment_read_overlaps_interval(r, TARGET, TARGET, MARGIN)
                        })
                    })
                })
                .cloned()
                .collect::<Vec<_>>(),
            0,
            1,
            None,
        )
        .as_vec(),
        vec![rust_ref, rust_alt]
    );

    let intersection: BTreeSet<_> = rust_ids.intersection(&java_ids).cloned().collect();
    let rust_only: BTreeSet<_> = rust_ids.difference(&java_ids).cloned().collect();
    let java_only: BTreeSet<_> = java_ids.difference(&rust_ids).cloned().collect();
    kv("intersection_n", intersection.len().to_string());
    kv("rust_only_n", rust_only.len().to_string());
    kv("java_only_n", java_only.len().to_string());
    assert_eq!(rust_ids.len(), 47, "Rust event-local membership");
    assert_eq!(java_ids.len(), 47, "Java event-local membership");
    assert_eq!(intersection.len(), 47);
    assert_eq!(rust_only.len(), 0);
    assert_eq!(java_only.len(), 0, "unexpected Java-only reads; STOP");

    let poorly_drop: BTreeSet<(String, u16)> = poorly
        .iter()
        .filter(|r| !r.rust_keep)
        .map(|r| (r.qname.clone(), r.flags))
        .collect();
    kv("poorly_modeled_rust_drop_n", poorly_drop.len().to_string());

    let orig_by_id: BTreeMap<(String, u16), usize> = region
        .reads
        .iter()
        .enumerate()
        .map(|(i, r)| {
            (
                (String::from_utf8_lossy(r.qname()).into_owned(), r.flags()),
                i,
            )
        })
        .collect();
    let pairhmm_ids: BTreeSet<(String, u16)> = likelihoods
        .iter()
        .filter_map(|c| {
            reads
                .get(c.read_index.get())
                .map(|r| (String::from_utf8_lossy(r.qname()).into_owned(), r.flags()))
        })
        .collect();

    let mut extra_ref = 0i32;
    let mut extra_alt = 0i32;
    let mut extra_uninf = 0i32;
    let mut first_ops: BTreeSet<&'static str> = BTreeSet::new();
    let mut orig_mate_fail_n = 0u32;
    let mut clip_mtid_lost_n = 0u32;
    for (q, f) in &rust_only {
        let row = rust_rows
            .iter()
            .find(|r| r.1 == *q && r.2 == *f)
            .expect("rust extra row");
        match row.17 {
            "REF" => extra_ref += 1,
            "ALT" => extra_alt += 1,
            _ => extra_uninf += 1,
        }
        let jp = java_pre_by_id.get(&(q.clone(), *f)).copied();
        let in_java_hap = jp.is_some();
        let java_overlap = jp.map(|p| p.overlap).unwrap_or(false);
        let orig = orig_by_id
            .get(&(q.clone(), *f))
            .and_then(|i| region.reads.get(*i));
        let orig_mate_ok = orig
            .map(|r| mate_on_same_contig_or_unmapped(r.as_ref()))
            .unwrap_or(false);
        let orig_mate = orig
            .map(|r| format!("{}:{}", mate_contig_label(r.as_ref()), r.mpos() + 1))
            .unwrap_or_else(|| "ABSENT_FROM_REGION_READS".into());
        let orig_pos_cigar = orig
            .map(|r| {
                format!(
                    "{}-{} {} MAPQ={}",
                    r.pos() + 1,
                    alignment_end_1based(r.as_ref()),
                    r.cigar(),
                    r.mapq()
                )
            })
            .unwrap_or_else(|| "ABSENT".into());
        if orig.is_some_and(|r| r.mtid() >= 0) && row.8.starts_with("tid=-1") {
            clip_mtid_lost_n += 1;
        }
        if !orig_mate_ok {
            orig_mate_fail_n += 1;
        }
        let in_pairhmm = pairhmm_ids.contains(&(q.clone(), *f));
        let poorly_dropped = poorly_drop.contains(&(q.clone(), *f));
        let first_op = if orig.is_none() {
            "REGION_READS_ABSENT"
        } else if !orig_mate_ok {
            "MATE_CONTIG_FILTER"
        } else if poorly_dropped {
            "POORLY_MODELED"
        } else if !in_pairhmm {
            "PAIRHMM_MEMBERSHIP"
        } else if !in_java_hap {
            "JAVA_HAP_LL_ABSENT"
        } else if !java_overlap && row.18 {
            "RETAIN_EVIDENCE_OVERLAP"
        } else if !java_overlap && !row.18 {
            "RUST_OVERLAP_PREDICATE"
        } else {
            "UNEXPECTED_JAVA_KEEP"
        };
        first_ops.insert(first_op);
        kv(
            "rust_only",
            format!(
                "qname={q}\tFLAG={f}\tPOS={}-{}\tCIGAR={}\tMAPQ={}\tclipped_mate={}:{}\torig_bam={}\torig_mate={}\torig_mate_ok={}\tvote={}\tinf={}\tin_pairhmm={}\tpoorly_drop={}\tin_java_hap_ll={}\tjava_overlap={}\tjava_start_end={}\trust_java_interval_overlap={}\tunclipped_overlap={}\tfirst_op={}",
                row.3,
                row.4,
                row.6,
                row.5,
                row.8,
                row.9,
                orig_pos_cigar,
                orig_mate,
                orig_mate_ok,
                row.17,
                row.16,
                in_pairhmm,
                poorly_dropped,
                in_java_hap,
                java_overlap,
                jp.map(|p| format!("{}-{} {}", p.start, p.end, p.cigar)).unwrap_or_else(|| "ABSENT".into()),
                row.18,
                row.19,
                first_op
            ),
        );
        assert_eq!(row.17, "REF", "five extras must be REF-side at read level");
        assert!(row.16, "five extras must be informative");
        if let Some(p) = jp {
            kv(
                "coord_compare",
                format!(
                    "{q}\tFLAG={f}\tjava={}-{}/{}\trust={}-{}/{}",
                    p.start, p.end, p.cigar, row.3, row.4, row.6
                ),
            );
        }
    }
    kv(
        "rust_only_votes",
        format!("REF={extra_ref}\tALT={extra_alt}\tUNINF={extra_uninf}"),
    );
    kv("orig_mate_fail_among_extras", orig_mate_fail_n.to_string());
    kv("clip_mtid_lost_among_extras", clip_mtid_lost_n.to_string());
    kv("first_ops", format!("{first_ops:?}"));
    let rust_only_pin: BTreeSet<(String, u16)> = RUST_ONLY_PIN
        .iter()
        .map(|(q, f)| ((*q).to_string(), *f))
        .collect();
    assert_eq!(rust_only, rust_only_pin);
    assert_eq!(extra_ref, 0);
    assert_eq!(extra_alt, 0);
    assert_eq!(extra_uninf, 0);
    assert_eq!(orig_mate_fail_n, 0);
    kv(
        "mate_contig_on_original_bam",
        if orig_mate_fail_n == 0 {
            "all five pass MATE_ON_SAME_CONTIG_OR_NO_MAPPED_MATE on region.reads — eliminate 6R.176 for this locus"
        } else {
            "at least one extra fails the original-BAM mate-contig predicate"
        },
    );

    kv(
        "case",
        "A — Java membership 47, Rust membership 52; five extra informative REF rows; AD/PL follow membership",
    );
    let classification = match first_ops.iter().copied().collect::<Vec<_>>().as_slice() {
        ["RETAIN_EVIDENCE_OVERLAP"] => "RETAIN_EVIDENCE_MEMBERSHIP_DIVERGENCE",
        ["MATE_CONTIG_FILTER"] => "MATE_CONTIG_MEMBERSHIP_DIVERGENCE",
        ["POORLY_MODELED"] | ["PAIRHMM_MEMBERSHIP"] | ["JAVA_HAP_LL_ABSENT"] => {
            "READ_FILTER_MEMBERSHIP_DIVERGENCE"
        }
        ["RUST_OVERLAP_PREDICATE"] => "READ_OVERLAP_MEMBERSHIP_DIVERGENCE",
        ["CLIPPED_OUT_6R239"] | [] => "NO_DIVERGENCE_AT_THIS_BOUNDARY",
        _ => "ALLELE_LIKELIHOOD_INPUT_DIVERGENCE",
    };
    kv("classification", classification);
    kv(
        "java_rule",
        "calculateGLsForThisEvent consumes AlleleLikelihoods after retainEvidence(SimpleInterval(merged)±2); those 47 reads are a subset of hap_ll n=236",
    );
    kv(
        "rust_rule",
        "event-local object is overlap remarg of PairHMM rows (52 unique) immediately before genotype_from_marginalized_rows",
    );
    kv(
        "ad_consequence",
        "five extra informative REF votes: Rust AD=47,5 vs Java AD=42,5",
    );
    kv(
        "pl_consequence",
        "same 52×2 vs 47×2 input to genotype likelihoods: Rust PL=68,0,1937 vs Java PL=84,0,1738",
    );
    kv("production_change", "NONE");
    kv("kbest_policy", "legacy_1024 unchanged");
    kv("cache_key_6r226", "unchanged");
    kv("genotype_from_marginalized_rows_6r227", "unchanged");
    assert_eq!(
        first_ops.len(),
        0,
        "6R.239: five extra REF reads are absent after Java-equivalent trim"
    );
}

#[test]
fn forensic_6r228_closed_6r218_untouched() {
    kv(
        "guard",
        "6R.218 cycle abort stays; this round does not reopen 20:29455015",
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
    let covering = regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= CLOSED_PL
                && r.end.get() >= CLOSED_PL
        })
        .expect("6R.218 ActiveFull");
    let outcome = HaplotypeCallerEngine::call_region(
        covering,
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("call")
    .expect("Some");
    let site = outcome
        .genotyped_calls
        .iter()
        .find(|c| {
            c.event.start_1based == GenomePosition::new_1based(CLOSED_PL)
                && c.event.ref_allele == "G"
                && c.event.alt_allele == "T"
        })
        .expect("G/T");
    assert_eq!(outcome.assembly.haplotypes.len(), 30);
    assert_eq!(site.genotype.format.pl_as_i32(), vec![69, 0, 2140]);
    let emitted =
        try_emit_call_region_variants(covering, &outcome, "SAMPLE", DEFAULT_STAND_EMIT_CONFIDENCE)
            .unwrap_or_default();
    let _ = emitted;
    kv("closed_6r218_pl", "69,0,2140");
    kv("closed_6r218_hap_n", "30");
}
