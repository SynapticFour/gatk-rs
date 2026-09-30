//! 6R.229: first Java hap_ll membership drop for the five Rust-only reads
//! at `20:29455379 G/A`.
//!
//! Proof-only. Frozen Java 4.4.0.0 SHA `2dbc025821bc5f686c423ff332a41e6cef892a77`.
//! PRODUCTION CHANGE: NONE.
//!
//! Live Java `hap-ll-membership-at-loc` tracks the frozen five through
//! region → trim → filterNonPassingReads → PairHMM → normalize →
//! `filterPoorlyModeledEvidence` → retainEvidence. All five are present
//! through normalize (bitmap `11111`) and absent after poorly-modeled
//! (`00000`). Mate-contig / clip / stub / filterNonPassingReads KEEP.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r229_java_hap_ll_membership_pipeline -- --nocapture --test-threads=1
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
    region_likelihoods_to_rows, HcGenotypingConfig, DEFAULT_INFORMATIVE_READ_OVERLAP_MARGIN,
};
use gatk_haplotypecaller::read_realignment::LOG_10_INFORMATIVE_THRESHOLD;
use gatk_haplotypecaller::read_unclip::alignment_end_1based;
use gatk_haplotypecaller::{
    begin_poorly_modeled_observe, call_disposition, flatten_assembly_regions,
    take_poorly_modeled_observe, traverse_assembly_region_walker, AssemblyRegionCallDisposition,
    CallRegionArgs, GenomePosition, HaplotypeCallerEngine, PoorlyModeledObserveRow,
    ReadFilterParams, WalkerTraversalConfig,
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
const MARGIN: i32 = DEFAULT_INFORMATIVE_READ_OVERLAP_MARGIN;
const JAVA_STAGES: &str = include_str!("forensic_6r229_java_stages.tsv");
const JAVA_POORLY: &str = include_str!("forensic_6r229_java_poorly.tsv");
const FIVE: &[(&str, u16)] = &[
    ("HISEQ1:13:H8G92ADXX:1:2102:7192:18079", 99),
    ("HISEQ1:9:H8962ADXX:1:1116:1789:43193", 163),
    ("HWI-D00360:6:H81VLADXX:1:1111:8050:25694", 163),
    ("HWI-D00360:6:H81VLADXX:1:2102:1733:39463", 163),
    ("HWI-D00360:7:H88WKADXX:2:1201:4043:96748", 163),
];
const JAVA_POORLY_PIN: &[(u16, usize, f64, f64)] = &[
    (99, 148, -12.324832, -8.0),
    (163, 113, -8.732083, -8.0),
    (163, 148, -8.682712, -8.0),
    (163, 148, -12.536733, -8.0),
    (163, 148, -15.402292, -8.0),
];

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    eprintln!("6R229\t{key}\t{}", value.as_ref());
}

fn parse_stages() -> BTreeMap<String, String> {
    let mut out = BTreeMap::new();
    for line in JAVA_STAGES.lines().skip(1) {
        if line.is_empty() {
            continue;
        }
        let (k, v) = line.split_once('\t').expect("stage tsv");
        out.insert(k.to_string(), v.to_string());
    }
    out
}

fn parse_poorly() -> Vec<(String, u16, usize, f64, f64, &'static str)> {
    let mut out = Vec::new();
    for line in JAVA_POORLY.lines().skip(1) {
        if line.is_empty() {
            continue;
        }
        let mut qname = String::new();
        let mut flags = 0u16;
        let mut qlen = 0usize;
        let mut max_ll = 0.0f64;
        let mut thresh = 0.0f64;
        let mut keep = "DROP";
        for field in line.split('\t') {
            if let Some(v) = field.strip_prefix("qname=") {
                qname = v.to_string();
            } else if let Some(v) = field.strip_prefix("flags=") {
                flags = v.parse().unwrap();
            } else if let Some(v) = field.strip_prefix("qlen=") {
                qlen = v.parse().unwrap();
            } else if let Some(v) = field.strip_prefix("max_ll=") {
                max_ll = v.parse().unwrap();
            } else if let Some(v) = field.strip_prefix("thresh=") {
                thresh = v.parse().unwrap();
            } else if let Some(v) = field.strip_prefix("keep=") {
                keep = if v == "KEEP" { "KEEP" } else { "DROP" };
            }
        }
        out.push((qname, flags, qlen, max_ll, thresh, keep));
    }
    out
}

fn vote_row(lr: f64, la: f64) -> &'static str {
    let conf = (lr - la).abs();
    if !(lr.is_finite() && la.is_finite() && conf > LOG_10_INFORMATIVE_THRESHOLD) {
        "UNINF"
    } else if lr > la {
        "REF"
    } else {
        "ALT"
    }
}

fn mate_on_same_contig_or_unmapped(rec: &rust_htslib::bam::Record) -> bool {
    if !rec.is_paired() || rec.is_mate_unmapped() || rec.is_unmapped() {
        return true;
    }
    rec.tid() == rec.mtid()
}

fn mate_contig_label(rec: &rust_htslib::bam::Record) -> String {
    if rec.mtid() < 0 {
        "unmapped".into()
    } else {
        format!("mtid={}", rec.mtid())
    }
}

fn orig_pos_cigar(rec: &rust_htslib::bam::Record) -> (i64, String) {
    let pos = match rec.aux(b"OP") {
        Ok(Aux::I32(v)) => i64::from(v),
        Ok(Aux::U32(v)) => i64::from(v),
        _ => rec.pos() + 1,
    };
    let cigar = match rec.aux(b"OC") {
        Ok(Aux::String(s)) => s.to_string(),
        _ => rec.cigar().to_string(),
    };
    (pos, cigar)
}

fn last_poorly<'a>(
    rows: &'a [PoorlyModeledObserveRow],
    qname: &str,
    flags: u16,
) -> Option<&'a PoorlyModeledObserveRow> {
    rows.iter()
        .filter(|r| r.qname == qname && r.flags == flags)
        .max_by_key(|r| r.pass)
}

#[test]
fn forensic_6r229_java_hap_ll_membership_pipeline() {
    kv("java_pin", JAVA_PIN);
    kv("target", format!("20:{TARGET} {TARGET_REF}/{TARGET_ALT}"));
    kv("production_change", "NONE");
    kv(
        "classification_candidates_eliminated",
        "region/clip/stub/filterNonPassingReads/mate/PairHMM-input/hap_ll-construction",
    );

    let stages = parse_stages();
    let java_poorly = parse_poorly();
    assert_eq!(stages.get("iterator_region_reads_five").unwrap(), "11111");
    assert_eq!(stages.get("after_assemble_finalize_five").unwrap(), "11111");
    assert_eq!(stages.get("after_trim_five").unwrap(), "11111");
    assert_eq!(stages.get("after_stub_five").unwrap(), "11111");
    assert_eq!(
        stages.get("after_filter_non_passing_five").unwrap(),
        "11111"
    );
    assert_eq!(stages.get("pairhmm_input_five").unwrap(), "11111");
    assert_eq!(stages.get("initial_hap_ll_five").unwrap(), "11111");
    assert_eq!(stages.get("after_pairhmm_kernel_five").unwrap(), "11111");
    assert_eq!(stages.get("after_normalize_five").unwrap(), "11111");
    assert_eq!(stages.get("after_poorly_modeled_five").unwrap(), "00000");
    assert_eq!(
        stages.get("before_calculateGLs_hap_ll_five").unwrap(),
        "00000"
    );
    assert_eq!(stages.get("after_marginalize_five").unwrap(), "00000");
    assert_eq!(stages.get("after_retainEvidence_five").unwrap(), "00000");
    assert_eq!(stages.get("after_normalize_count").unwrap(), "245");
    assert_eq!(stages.get("after_poorly_modeled_count").unwrap(), "236");
    assert_eq!(stages.get("after_retainEvidence_count").unwrap(), "47");
    assert_eq!(stages.get("mq_threshold").unwrap(), "20");
    assert_eq!(stages.get("keep_rg").unwrap(), "null");
    kv("java_variant_span", stages.get("variant_span").unwrap());
    kv("java_variant_padded", stages.get("variant_padded").unwrap());
    kv(
        "java_lifecycle",
        "BAM/iterator → finalizeRegion → trim(hardClipToRegion padded) → stub≥10 → filterNonPassingReads → PairHMM AlleleLikelihoods → normalizeLikelihoods → filterPoorlyModeledEvidence → realign/changeEvidence → marginalize → retainEvidence(±2) → calculateGLsForThisEvent",
    );

    assert_eq!(java_poorly.len(), 5);
    for (i, row) in java_poorly.iter().enumerate() {
        assert_eq!(row.0, FIVE[i].0);
        assert_eq!(row.1, FIVE[i].1);
        assert_eq!(row.2, JAVA_POORLY_PIN[i].1);
        assert!((row.3 - JAVA_POORLY_PIN[i].2).abs() < 1e-5);
        assert!((row.4 - JAVA_POORLY_PIN[i].3).abs() < 1e-9);
        assert_eq!(row.5, "DROP");
        assert!(
            row.3 < row.4,
            "{} max_ll {} must be < thresh {}",
            row.0,
            row.3,
            row.4
        );
        kv(
            "java_poorly",
            format!(
                "qname={}\tFLAG={}\tqlen={}\tmax_ll={:.6}\tthresh={:.1}\tkeep=DROP",
                row.0, row.1, row.2, row.3, row.4
            ),
        );
    }

    kv(
        "java_matrix",
        "stage                         count   five\n\
         iterator_region_reads           319    11111\n\
         after_assemble_finalize         319    11111\n\
         after_trim                      251    11111\n\
         after_stub                      246    11111\n\
         after_filter_non_passing        245    11111\n\
         pairhmm_input                   245    11111\n\
         initial_hap_ll                  245    11111\n\
         after_pairhmm_kernel            245    11111\n\
         after_normalize                 245    11111\n\
         after_poorly_modeled            236    00000  FIRST JAVA DROP\n\
         before_calculateGLs_hap_ll      236    00000\n\
         after_marginalize               236    00000\n\
         after_retainEvidence             47    00000",
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
    let orig_ids: BTreeSet<(String, u16)> = region
        .reads
        .iter()
        .map(|r| (String::from_utf8_lossy(r.qname()).into_owned(), r.flags()))
        .collect();
    let pairhmm_ids: BTreeSet<(String, u16)> = likelihoods
        .iter()
        .filter_map(|c| {
            reads
                .get(c.read_index.get())
                .map(|r| (String::from_utf8_lossy(r.qname()).into_owned(), r.flags()))
        })
        .collect();
    kv("rust_region_reads_n", region.reads.len().to_string());
    kv("rust_genotyping_reads_n", reads.len().to_string());
    kv("rust_pairhmm_unique_n", pairhmm_ids.len().to_string());

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
    let mut rust_overlap: BTreeSet<(String, u16)> = BTreeSet::new();
    let mut rust_votes = BTreeMap::<(String, u16), &'static str>::new();
    for row in &marg_all {
        let Some(idx) = row.matrix_read_index() else {
            continue;
        };
        let rec = &reads[idx];
        if !java_alignment_read_overlaps_interval(rec, TARGET, TARGET, MARGIN) {
            continue;
        }
        let id = (
            String::from_utf8_lossy(rec.qname()).into_owned(),
            rec.flags(),
        );
        rust_overlap.insert(id.clone());
        rust_votes.insert(
            id,
            vote_row(
                row.haplotype_log10_likelihoods[0],
                row.haplotype_log10_likelihoods[1],
            ),
        );
    }
    kv("rust_event_local_n", rust_overlap.len().to_string());
    assert_eq!(rust_overlap.len(), 47);

    kv(
        "five_header",
        "QNAME\tFLAG\torigPOS\torigEnd\torigCIGAR\tmate\tMAPQ\tqlen\tclippedPOS\tclippedCIGAR\trust_region\trust_pairhmm\trust_poorly_keep\tjava_hap_ll\tjava_poorly",
    );
    let mut last_common = "after_normalize";
    let mut first_java = BTreeSet::new();
    for (i, (q, f)) in FIVE.iter().enumerate() {
        let orig = region
            .reads
            .iter()
            .find(|r| String::from_utf8_lossy(r.qname()) == *q && r.flags() == *f)
            .unwrap_or_else(|| panic!("five read missing from region.reads: {q} FLAG={f}"));
        let geno = reads
            .iter()
            .find(|r| String::from_utf8_lossy(r.qname()) == *q && r.flags() == *f);
        let orig_start = orig.pos() + 1;
        let orig_end = alignment_end_1based(orig.as_ref()) as i64;
        let (op, oc) = orig_pos_cigar(orig.as_ref());
        let mate_ok = mate_on_same_contig_or_unmapped(orig.as_ref());
        let in_region = orig_ids.contains(&((*q).to_string(), *f));
        let in_pairhmm = pairhmm_ids.contains(&((*q).to_string(), *f));
        let pr = last_poorly(&poorly, q, *f).expect("poorly observe row");
        let java_keep = java_poorly[i].5;
        let clip = geno
            .map(|r| format!("{}-{} {}", r.pos() + 1, alignment_end_1based(r), r.cigar()))
            .unwrap_or_else(|| "ABSENT".into());
        kv(
            "five",
            format!(
                "qname={q}\tFLAG={f}\torigPOS={orig_start}-{orig_end}\torigCIGAR={}\tOP/OC={op}/{oc}\tmate_tid={}\tmatePOS={}\tMAPQ={}\tqlen={}\tclipped={clip}\tregion={in_region}\tpairhmm={in_pairhmm}\trust_java_equiv_keep={}\trust_keep={}\textra_retain={}\trust_max_ll={:.6}\trust_thresh={:.1}\tjava_hap_ll=0\tjava_poorly={java_keep}\tmatePred={}",
                orig.cigar(),
                mate_contig_label(orig.as_ref()),
                orig.mpos() + 1,
                orig.mapq(),
                orig.seq_len(),
                pr.java_equiv_keep,
                pr.rust_keep,
                pr.extra_retain,
                pr.max_ll,
                pr.threshold,
                if mate_ok { "PASS" } else { "FAIL" }
            ),
        );
        assert!(in_region, "{q} must exist in initial region reads");
        assert!(mate_ok, "{q} original-BAM mate predicate must PASS");
        assert!(
            !in_pairhmm,
            "{q} 6R.239: extra five are clipped out of PairHMM with Java trim start"
        );
        assert!(
            rust_votes.get(&((*q).to_string(), *f)).is_none(),
            "{q} must not contribute a remarg vote"
        );
        continue;
        first_java.insert("filterPoorlyModeledEvidence");
        last_common = "after_normalize / initial hap_ll";
        kv(
            "first_drop",
            format!("{q} FLAG={f} | after_normalize YES → filterPoorlyModeledEvidence NO"),
        );
    }
    kv("last_common_membership", last_common);
    kv(
        "first_java_stage",
        "filterPoorlyModeledEvidence (245→236, five 11111→00000)",
    );
    kv(
        "exact_java_predicate",
        "ReadLikelihoodCalculationEngine.filterPoorlyModeledEvidence / log10MinTrueLikelihood: thresh=min(2, ceil(qlen*0.02))*-4; KEEP iff !(max_ll < thresh). All five max_ll < -8 → DROP.",
    );
    kv(
        "rust_counterpart",
        "filter_poorly_modeled_region_read_likelihoods uses the same thresh; extra_retain=false; Rust java_equiv_keep=true because max_ll is ≥ -8 (≈ -2.5…-7.2) while Java max_ll is < -8.",
    );
    assert_eq!(first_java.len(), 0);

    let site = outcome
        .genotyped_calls
        .iter()
        .find(|c| {
            c.event.start_1based == GenomePosition::new_1based(TARGET)
                && c.event.ref_allele == TARGET_REF
                && c.event.alt_allele == TARGET_ALT
        })
        .expect("target site");
    let ad = site.genotype.format.ad_as_i32();
    let pl = site.genotype.format.pl_as_i32();
    kv("rust_ad", format!("{ad:?}"));
    kv("rust_pl", format!("{pl:?}"));
    kv("java_ad", format!("{JAVA_AD:?}"));
    kv("java_pl", format!("{JAVA_PL:?}"));
    assert_eq!(ad, vec![42, 5]);
    assert_eq!(pl, vec![84, 0, 1738]);
    kv(
        "ad_pl_consequence",
        "five extra informative REF votes remain in Rust hap_ll / overlap remarg: AD 47,5 PL 68,0,1937 vs Java AD 42,5 PL 84,0,1738",
    );
    kv("outcome", "C — survive region/clip/filterNonPassingReads/PairHMM/initial hap_ll; Java filterPoorlyModeledEvidence DROPs all five");
    kv("classification", "POORLY_MODELED_FILTER_DIVERGENCE");
    kv("cache_key_6r226", "unchanged");
    kv("genotype_from_marginalized_rows_6r227", "unchanged");
    kv("retainEvidence_6r228", "unchanged");
}

#[test]
fn forensic_6r229_closed_6r218_untouched() {
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
        .expect("closed 6R.218 site");
    assert_eq!(site.genotype.format.pl_as_i32(), vec![69, 0, 2140]);
    kv("closed_6r218_pl", "69,0,2140");
}
