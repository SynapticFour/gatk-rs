//! 6R.157: Class-A3 calculator-boundary on the frozen 20:29456196 A/T site
//! (live check that `pl_gt` is the min-PL index, not a validity flag).
//! Diagnostic-only `unbounded_diagnostic`. Skipped unless `HOLDOUT_6R157=1`.
//!
//! ```text
//! HOLDOUT_6R157=1 GATK_RS_EXPERIMENTAL_KBEST_POLICY=unbounded_diagnostic \
//!   cargo test -p gatk-haplotypecaller --test holdout_6r157_class_a3_boundary -- --nocapture --test-threads=1
//! ```

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::event_map::{build_per_haplotype_variation_events, VariationEvent};
use gatk_haplotypecaller::genotyping::{
    best_pl_index, emit_genotype_format_fields, ReadLikelihoodRow,
};
use gatk_haplotypecaller::hc_allele_mapping::create_allele_mapper_with_events;
use gatk_haplotypecaller::hc_genotyping_engine::java_alignment_read_overlaps_interval;
use gatk_haplotypecaller::query_index_at_reference_position;
use gatk_haplotypecaller::{
    biallelic_genotype_log10_likelihoods_gatk, call_disposition, flatten_assembly_regions,
    marginalize_rows_to_biallelic_alleles, region_likelihoods_to_rows,
    traverse_assembly_region_walker, AssemblyRegionCallDisposition, CallRegionArgs, GenomePosition,
    HaplotypeCallerEngine, InformativeAd, ReadFilterParams, SparsePlShape, WalkerTraversalConfig,
    DEFAULT_INFORMATIVE_READ_OVERLAP_MARGIN,
};
use rust_htslib::bam::record::CigarString;
use rust_htslib::bam::Record;
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::path::{Path, PathBuf};

const INTERVAL: &str = "20:29455000-29456500";
const BAM_REL: &str = "parity/giab/runs/local-pairhmm-diff/HG001.20-29455000-29456500.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const TARGET: u64 = 29_456_196;
const JAVA_ACTIVE_START: u64 = 29_455_995;
const JAVA_ACTIVE_END: u64 = 29_456_293;
const JAVA_PAD_START: u64 = 29_455_895;
const JAVA_PAD_END: u64 = 29_456_393;
const JAVA_FLOAT_DUMP_REL: &str = "parity/giab/runs/local-pairhmm-diff/6r143_java.txt";
const MERGED_REF: &str = "A";
const MERGED_ALT: &str = "T";
const MARGIN: i32 = DEFAULT_INFORMATIVE_READ_OVERLAP_MARGIN;
const MATE_QNAME: &str = "HWI-D00360:7:H88WKADXX:2:2107:6787:30989";

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    println!("6R157\t{key}\t{}", value.as_ref());
}

struct EnvGuard {
    key: &'static str,
    prior: Option<String>,
}

impl EnvGuard {
    fn set(key: &'static str, value: &str) -> Self {
        let prior = std::env::var(key).ok();
        std::env::set_var(key, value);
        Self { key, prior }
    }
}

impl Drop for EnvGuard {
    fn drop(&mut self) {
        match &self.prior {
            Some(v) => std::env::set_var(self.key, v),
            None => std::env::remove_var(self.key),
        }
    }
}

type ReadKey = (String, u16);

fn parse_kv_fields(parts: &[&str]) -> BTreeMap<String, String> {
    let mut m = BTreeMap::new();
    for tok in parts {
        if let Some((k, v)) = tok.split_once('=') {
            m.insert(k.to_string(), v.to_string());
        }
    }
    m
}

fn load_java_stage_reads(path: &Path, stage: &str) -> BTreeSet<ReadKey> {
    let text = std::fs::read_to_string(path).expect("java dump");
    let mut out = BTreeSet::new();
    for line in text.lines() {
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 4 || parts[0] != "6R143" || parts[1] != stage {
            continue;
        }
        let qname = parts[2].to_string();
        let kv = parse_kv_fields(&parts[3..]);
        let flags: u16 = kv.get("flags").and_then(|s| s.parse().ok()).unwrap_or(0);
        out.insert((qname, flags));
    }
    out
}

fn rec_key(rec: &Record) -> ReadKey {
    (
        String::from_utf8_lossy(rec.qname()).into_owned(),
        rec.flags(),
    )
}

fn snp_base_at(rec: &Record) -> Option<u8> {
    let loc0 = (TARGET as i64).saturating_sub(1);
    if rec.is_unmapped() {
        return None;
    }
    let cigar = CigarString(rec.cigar().iter().copied().collect());
    let qi = query_index_at_reference_position(rec.pos(), &cigar, loc0)?;
    let seq = rec.seq();
    if qi >= seq.len() {
        return None;
    }
    Some(seq.as_bytes()[qi].to_ascii_uppercase())
}

fn snp_pileup_ad<'a, I>(reads: I, dedupe_qname: bool) -> (i32, i32)
where
    I: IntoIterator<Item = &'a Record>,
{
    let ref_b = MERGED_REF.as_bytes()[0].to_ascii_uppercase();
    let alt_b = MERGED_ALT.as_bytes()[0].to_ascii_uppercase();
    let mut seen = BTreeSet::new();
    let mut rr = 0i32;
    let mut ra = 0i32;
    for rec in reads {
        if rec.is_unmapped() {
            continue;
        }
        if dedupe_qname && !seen.insert(rec.qname().to_vec()) {
            continue;
        }
        match snp_base_at(rec) {
            Some(b) if b == alt_b => ra += 1,
            Some(b) if b == ref_b => rr += 1,
            _ => {}
        }
    }
    (rr, ra)
}

fn fmt_state(gt: usize, pl: &[i32], ad: &[i32], gq: i32, dp: i32) -> String {
    format!(
        "GT={}/{} PL={} AD={} GQ={gq} DP={dp} alleles=A,T",
        gt / 2,
        gt % 2 + gt / 2,
        pl.iter()
            .map(|x| x.to_string())
            .collect::<Vec<_>>()
            .join(","),
        ad.iter()
            .map(|x| x.to_string())
            .collect::<Vec<_>>()
            .join(","),
    )
}

#[test]
fn holdout_6r157_class_a3_boundary() {
    if std::env::var("HOLDOUT_6R157").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R157=1");
        return;
    }
    assert_ne!(
        std::env::var("GATK_RS_DIAGNOSTIC_SKIP_SITE_RESHAPE")
            .ok()
            .as_deref(),
        Some("1"),
        "6R.157 must not alter the default production result via the diagnostic skip"
    );
    let _policy = EnvGuard::set("GATK_RS_EXPERIMENTAL_KBEST_POLICY", "unbounded_diagnostic");
    kv(
        "parity_mode",
        "unbounded_diagnostic\tnote=diagnostic_not_production_default",
    );
    kv("variant", format!("20:{TARGET} {MERGED_REF}/{MERGED_ALT}"));
    kv("alleles", "A,T");
    kv("ploidy", "2");
    kv("frozen_pre_a3", "GT=0/0 PL=0,5,1174 AD=37,4 GQ=5 DP=41");

    let root = repo_root();
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
    let covering = regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= TARGET
                && r.end.get() >= TARGET
        })
        .expect("ActiveFull")
        .clone();
    let mut java_bounds = covering.clone();
    java_bounds.start = GenomePosition::new_1based(JAVA_ACTIVE_START);
    java_bounds.end = GenomePosition::new_1based(JAVA_ACTIVE_END);
    java_bounds.extended_start = GenomePosition::new_1based(JAVA_PAD_START);
    java_bounds.extended_end = GenomePosition::new_1based(JAVA_PAD_END);

    let args = CallRegionArgs::strict_java();
    let outcome = HaplotypeCallerEngine::call_region(&java_bounds, &dict, &ref_fasta, &args)
        .expect("call")
        .expect("outcome");

    let haps = &outcome.assembly.haplotypes;
    assert_eq!(haps.len(), 25);
    let (full_ref, full_pad) = outcome.assembly.event_map_reference();
    let apply_pad = outcome
        .assembly
        .haplotypes
        .iter()
        .find(|h| h.is_reference)
        .and_then(|h| h.genome_loc.as_ref().map(|g| g.start.get()))
        .unwrap_or(full_pad);
    let hap_cache = build_per_haplotype_variation_events(
        haps,
        full_ref,
        full_pad,
        outcome.assembly.max_mnp_distance(),
        "20",
    );
    let merged = VariationEvent::from_alleles("20", TARGET, MERGED_REF, MERGED_ALT);
    let mapping = create_allele_mapper_with_events(
        &merged,
        TARGET,
        haps,
        apply_pad,
        outcome.assembly.reference_bases(),
        outcome.assembly.max_mnp_distance(),
        true,
        Some(&hap_cache),
    );

    let rust_rows = region_likelihoods_to_rows(&outcome.read_likelihoods, haps.len());
    let rust_marg = marginalize_rows_to_biallelic_alleles(
        &rust_rows,
        &mapping.ref_haplotype_indices,
        &mapping.alt_haplotype_indices,
    );
    let read_by_idx: Vec<ReadKey> = outcome
        .genotyping_reads
        .iter()
        .map(|r| rec_key(r.as_ref()))
        .collect();
    let rec_by_key: HashMap<ReadKey, &Record> = outcome
        .genotyping_reads
        .iter()
        .map(|r| (rec_key(r.as_ref()), r.as_ref()))
        .collect();

    let java_keep = load_java_stage_reads(&root.join(JAVA_FLOAT_DUMP_REL), "stored");
    assert_eq!(java_keep.len(), 201);
    let mut retained = BTreeSet::new();
    for (rk, rec) in &rec_by_key {
        if java_keep.contains(rk)
            && java_alignment_read_overlaps_interval(rec, TARGET, TARGET, MARGIN)
        {
            retained.insert(rk.clone());
        }
    }
    kv("retained_reads", format!("{}", retained.len()));
    let mate99: ReadKey = (MATE_QNAME.to_string(), 99);
    let mate147: ReadKey = (MATE_QNAME.to_string(), 147);

    let mut first_rows = Vec::new();
    for row in &rust_marg {
        let Some(rk) = read_by_idx.get(row.read_index) else {
            continue;
        };
        if !retained.contains(rk) {
            continue;
        }
        first_rows.push(ReadLikelihoodRow {
            read_index: row.read_index,
            read_id: String::new(),
            haplotype_log10_likelihoods: row.haplotype_log10_likelihoods.clone(),
        });
    }
    let first_ad = InformativeAd::from_marginalized_rows(&first_rows, 0, 1, None);
    let first_gls = biallelic_genotype_log10_likelihoods_gatk(&first_rows, 0, 1);
    let first_fmt =
        emit_genotype_format_fields(&first_gls, &first_ad.as_vec()).expect("first write");
    let first_pl = first_fmt.pl_as_i32();
    let first_gt = best_pl_index(&first_fmt.pl);

    let n_region = covering.reads.len();
    let n_secondary = covering.reads.iter().filter(|r| r.is_secondary()).count();
    let n_supp = covering
        .reads
        .iter()
        .filter(|r| r.is_supplementary())
        .count();
    let n_dup = covering.reads.iter().filter(|r| r.is_duplicate()).count();
    let mut seen = BTreeSet::new();
    let mut qname_claimed_without_base = 0usize;
    let mut covering_after_first_wins = 0usize;
    for rec in covering.reads.iter() {
        if rec.is_unmapped() {
            continue;
        }
        let qn = rec.qname().to_vec();
        if !seen.insert(qn) {
            continue;
        }
        match snp_base_at(rec) {
            Some(_) => covering_after_first_wins += 1,
            None => qname_claimed_without_base += 1,
        }
    }
    let region_pileup = snp_pileup_ad(covering.reads.iter().map(|r| r.as_ref()), true);
    kv(
        "pileup_source",
        "QNAME-first-wins region.reads (untrimmed active-region; not realigned)",
    );
    kv(
        "pileup_population",
        format!(
            "n_region={n_region} n_geno={} n_retain={} secondary={n_secondary} supplementary={n_supp} duplicate={n_dup} first_wins_claimed_without_locus_base={qname_claimed_without_base} first_wins_with_base={covering_after_first_wins}",
            outcome.genotyping_reads.len(),
            retained.len()
        ),
    );
    kv(
        "pileup_depths",
        format!("{},{}", region_pileup.0, region_pileup.1),
    );

    let (pr, pa) = region_pileup;
    let info_ref = first_ad.ref_depth;
    let info_alt = first_ad.alt_depth;
    let preds = [
        ("not_p12_chr2", true),
        ("ad_len>=2", first_fmt.ad.len() >= 2),
        ("pileup_ref>=1", pr >= 1),
        ("pileup_alt>=1", pa >= 1),
        ("pileup_ref*2>=pileup_alt", pr.saturating_mul(2) >= pa),
        ("is_snp_A_T", true),
        ("pileup_alt*2>=pileup_ref", pa.saturating_mul(2) >= pr),
        ("info_ref>=1", info_ref >= 1),
        ("info_alt>=1", info_alt >= 1),
        (
            "info_ref>=3*info_alt",
            info_ref >= info_alt.saturating_mul(3),
        ),
        ("class_a", info_ref == 0 && info_alt >= 2),
        ("pl_gt==2", first_gt == 2),
        ("pl_gt==1_het", first_gt == 1),
        ("pl_gt==0_hom_ref", first_gt == 0),
        ("sample_count_checked", false),
        ("ploidy_checked", false),
        ("mapq_in_pileup_fn", false),
        ("baseq_in_pileup_fn", false),
        ("likelihood_derived_pileup", false),
    ];
    kv(
        "a3_predicates",
        preds
            .iter()
            .map(|(k, v)| format!("{k}={v}"))
            .collect::<Vec<_>>()
            .join(" "),
    );
    let class_a3 = first_fmt.ad.len() >= 2
        && pr >= 1
        && pa >= 1
        && pr.saturating_mul(2) >= pa
        && pa.saturating_mul(2) >= pr
        && info_ref >= 1
        && info_alt >= 1
        && info_ref >= info_alt.saturating_mul(3);
    kv("class_a3_fires", format!("{class_a3}"));
    kv("sparse_path", format!("{}", class_a3 && first_gt != 1));

    let cf_a = emit_genotype_format_fields(&first_gls, &[pr, pa]).expect("cf A");
    let het = SparsePlShape::from_pileup_depths(pr, pa);
    let cf_b = emit_genotype_format_fields(&het.gl_vec(), &first_ad.as_vec()).expect("cf B");
    let cf_c = first_fmt.clone();
    let cf_d = first_fmt.clone();
    kv(
        "counterfactual_A",
        fmt_state(
            best_pl_index(&cf_a.pl),
            &cf_a.pl_as_i32(),
            &cf_a.ad_as_i32(),
            cf_a.gq.as_i32(),
            cf_a.dp.as_i32(),
        ),
    );
    kv(
        "counterfactual_B",
        fmt_state(
            best_pl_index(&cf_b.pl),
            &cf_b.pl_as_i32(),
            &cf_b.ad_as_i32(),
            cf_b.gq.as_i32(),
            cf_b.dp.as_i32(),
        ),
    );
    kv(
        "counterfactual_C",
        fmt_state(
            best_pl_index(&cf_c.pl),
            &cf_c.pl_as_i32(),
            &cf_c.ad_as_i32(),
            cf_c.gq.as_i32(),
            cf_c.dp.as_i32(),
        ),
    );
    kv(
        "counterfactual_D",
        fmt_state(
            best_pl_index(&cf_d.pl),
            &cf_d.pl_as_i32(),
            &cf_d.ad_as_i32(),
            cf_d.gq.as_i32(),
            cf_d.dp.as_i32(),
        ),
    );

    let prod = outcome.genotyped_calls.iter().find(|c| {
        c.event.start_1based.get() == TARGET
            && c.event.ref_allele == MERGED_REF
            && c.event.alt_allele == MERGED_ALT
    });
    match prod {
        Some(c) => kv(
            "production_genotyped_calls",
            fmt_state(
                best_pl_index(&c.genotype.format.pl),
                &c.genotype.format.pl_as_i32(),
                &c.genotype.format.ad_as_i32(),
                c.genotype.format.gq.as_i32(),
                c.genotype.format.dp.as_i32(),
            ),
        ),
        None => kv("production_genotyped_calls", "absent"),
    }
    kv(
        "gq_dependency",
        "sparse path: emit_genotype_format_fields on SparsePlShape::Het GLs; GQ=PL gap 36",
    );
    kv(
        "dp_dependency",
        "sparse path: DP = pileup_ref+pileup_alt = 46; not an independent A3 formula",
    );
    kv("classification", "CANDIDATE_B_JAVA_EQUIVALENT_AT_A3_CALLER");
    kv("production_change", "YES_via_6R158");
    kv(
        "pl_gt_meaning",
        "first min-PL index 0=0/0 1=0/1 2=1/1; not a validity flag",
    );

    assert_eq!(retained.len(), 41);
    assert!(retained.contains(&mate99) && retained.contains(&mate147));
    assert_eq!(first_ad.as_vec(), vec![37, 4]);
    assert_eq!(first_pl, vec![0, 5, 1174]);
    assert_eq!(first_gt, 0);
    assert_eq!(first_fmt.gq.as_i32(), 5);
    assert_eq!(first_fmt.dp.as_i32(), 41);
    assert_eq!(region_pileup, (22, 24));
    assert!(class_a3);
    assert_eq!(het, SparsePlShape::Het);
    if let Some(c) = prod {
        assert_ne!(c.genotype.format.pl_as_i32(), vec![81, 0, 36]);
        assert_ne!(c.genotype.format.ad_as_i32(), vec![22, 24]);
        assert_eq!(best_pl_index(&c.genotype.format.pl), 0);
        assert_eq!(c.genotype.format.pl_as_i32(), vec![0, 5, 1174]);
        assert_eq!(c.genotype.format.ad_as_i32(), vec![37, 4]);
    }
}
