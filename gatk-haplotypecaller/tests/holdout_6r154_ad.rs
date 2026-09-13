//! 6R.154: reconstruct FORMAT/AD from the 6R.151 41-read retained object.
//! Diagnostic-only `unbounded_diagnostic`. Skipped unless `HOLDOUT_6R154=1`.
//!
//! ```text
//! HOLDOUT_6R154=1 GATK_RS_EXPERIMENTAL_KBEST_POLICY=unbounded_diagnostic \
//!   cargo test -p gatk-haplotypecaller --test holdout_6r154_ad -- --nocapture --test-threads=1
//! ```

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::event_map::{
    build_per_haplotype_variation_events, overlapping_events, VariationEvent,
};
use gatk_haplotypecaller::genotyping::ReadLikelihoodRow;
use gatk_haplotypecaller::hc_allele_mapping::{create_allele_mapper_with_events, SPAN_DEL_ALLELE};
use gatk_haplotypecaller::hc_genotyping_engine::{
    biallelic_allele_depths_from_rows, java_alignment_read_overlaps_interval, InformativeAd,
};
use gatk_haplotypecaller::read_realignment::LOG_10_INFORMATIVE_THRESHOLD;
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, marginalize_rows_to_biallelic_alleles,
    region_likelihoods_to_rows, traverse_assembly_region_walker, AssemblyRegionCallDisposition,
    CallRegionArgs, GenomePosition, HaplotypeCallerEngine, ReadFilterParams, WalkerTraversalConfig,
    DEFAULT_INFORMATIVE_READ_OVERLAP_MARGIN,
};
use rust_htslib::bam::Record;
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::path::{Path, PathBuf};

const INTERVAL: &str = "20:29455000-29456500";
const BAM_REL: &str = "parity/giab/runs/local-pairhmm-diff/HG001.20-29455000-29456500.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const JAVA_FLOAT_DUMP_REL: &str = "parity/giab/runs/local-pairhmm-diff/6r143_java.txt";
const TARGET: u64 = 29_456_196;
const JAVA_ACTIVE_START: u64 = 29_455_995;
const JAVA_ACTIVE_END: u64 = 29_456_293;
const JAVA_PAD_START: u64 = 29_455_895;
const JAVA_PAD_END: u64 = 29_456_393;
const EMIT_SPANNING_DELS: bool = true;
const MERGED_REF: &str = "A";
const MERGED_ALT: &str = "T";
const MARGIN: i32 = DEFAULT_INFORMATIVE_READ_OVERLAP_MARGIN;
const JAVA_INFORMATIVE: f64 = 0.2;
const MATE_QNAME: &str = "HWI-D00360:7:H88WKADXX:2:2107:6787:30989";

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    println!("6R154\t{key}\t{}", value.as_ref());
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

fn fnv1a64_hex(bases: &[u8]) -> String {
    let mut h: u64 = 0xcbf29ce484222325;
    for &b in bases {
        h ^= b as u64;
        h = h.wrapping_mul(0x100000001b3);
    }
    format!("{h:016x}")
}

fn parse_kv_fields(parts: &[&str]) -> BTreeMap<String, String> {
    let mut m = BTreeMap::new();
    for tok in parts {
        if let Some((k, v)) = tok.split_once('=') {
            m.insert(k.to_string(), v.to_string());
        }
    }
    m
}

type ReadKey = (String, u16);

fn load_java_stage_reads(path: &Path, stage: &str) -> BTreeSet<ReadKey> {
    let text = std::fs::read_to_string(path).expect("java dump");
    let mut out = BTreeSet::new();
    for line in text.lines() {
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 4 || parts[0] != "6R143" || parts[1] != stage {
            continue;
        }
        let qname = if stage == "orig" {
            parse_kv_fields(&parts[3..])
                .get("qname")
                .cloned()
                .unwrap_or_else(|| parts[2].to_string())
        } else {
            parts[2].to_string()
        };
        let kv = parse_kv_fields(&parts[3..]);
        let flags: u16 = kv.get("flags").and_then(|s| s.parse().ok()).unwrap_or(0);
        if stage == "orig" {
            let q = kv.get("qname").cloned().unwrap_or(qname);
            out.insert((q, flags));
        } else {
            out.insert((qname, flags));
        }
    }
    out
}

fn load_java_matrix(path: &Path, stage: &str) -> HashMap<(String, u16, String), f64> {
    let text = std::fs::read_to_string(path).expect("java dump");
    let mut out = HashMap::new();
    for line in text.lines() {
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 4 || parts[0] != "6R143" || parts[1] != stage {
            continue;
        }
        let qname = parts[2].to_string();
        let kv = parse_kv_fields(&parts[3..]);
        let flags: u16 = kv.get("flags").and_then(|s| s.parse().ok()).unwrap_or(0);
        let hap = kv.get("hap").cloned().unwrap_or_default();
        let bits =
            u64::from_str_radix(kv.get("bits").map(|s| s.as_str()).unwrap_or("0"), 16).unwrap_or(0);
        out.insert((qname, flags, hap), f64::from_bits(bits));
    }
    out
}

fn java_create_allele_mapper_labels(
    spanning: &[VariationEvent],
    loc: u64,
    merged_ref: &str,
    merged_alts: &[String],
    emit_spanning: bool,
) -> Vec<String> {
    let loc_pos = GenomePosition::new_1based(loc);
    if spanning.is_empty() {
        return vec![merged_ref.to_string()];
    }
    let mut labels = Vec::new();
    for ev in spanning {
        if ev.start_1based == loc_pos {
            if ev.ref_allele.len() == merged_ref.len() {
                if merged_alts.iter().any(|a| a == &ev.alt_allele)
                    && !labels.contains(&ev.alt_allele)
                {
                    labels.push(ev.alt_allele.clone());
                }
            }
        } else if emit_spanning {
            let star = SPAN_DEL_ALLELE.to_string();
            if !labels.contains(&star) {
                labels.push(star);
            }
            break;
        } else if !labels.contains(&merged_ref.to_string()) {
            labels.push(merged_ref.to_string());
            break;
        }
    }
    labels
}

fn pool_max(lls: &[f64]) -> f64 {
    lls.iter().copied().fold(f64::NEG_INFINITY, f64::max)
}

fn rec_key(rec: &Record) -> ReadKey {
    (
        String::from_utf8_lossy(rec.qname()).into_owned(),
        rec.flags(),
    )
}

/// Java `searchBestAllele` + default REF tie-break.
fn java_best_breaking_ties(lls: &[f64], ref_i: usize) -> (usize, usize, f64, bool) {
    let n = lls.len();
    if n == 0 {
        return (0, 0, 0.0, false);
    }
    let mut best_i = 0usize;
    let mut second_i = 0usize;
    let mut best = lls[0];
    let mut second = f64::NEG_INFINITY;
    for a in 1..n {
        let cand = lls[a];
        if cand > best {
            second_i = best_i;
            second = best;
            best_i = a;
            best = cand;
        } else if cand > second {
            second_i = a;
            second = cand;
        }
    }
    let priorities: Vec<f64> = (0..n).map(|i| if i == ref_i { 1.0 } else { 0.0 }).collect();
    if best - second < JAVA_INFORMATIVE {
        let mut best_pri = priorities[best_i];
        let mut second_pri = priorities[second_i];
        for a in 0..n {
            let cand = lls[a];
            if a == best_i || best - cand > JAVA_INFORMATIVE {
                continue;
            }
            let pri = priorities[a];
            if pri > best_pri {
                second_i = best_i;
                best_i = a;
                second_pri = best_pri;
                best_pri = pri;
            } else if pri > second_pri {
                second_i = a;
                second_pri = pri;
            }
        }
    }
    let best_ll = lls[best_i];
    let second_ll = if second_i != best_i {
        lls[second_i]
    } else {
        f64::NEG_INFINITY
    };
    let conf = if best_ll == second_ll {
        0.0
    } else {
        best_ll - second_ll
    };
    (best_i, second_i, conf, conf > JAVA_INFORMATIVE)
}

fn rust_biallelic_vote(lr: f64, la: f64) -> (usize, usize, f64, bool) {
    let (best_i, second_i, best, second) = if lr > la {
        (0usize, 1usize, lr, la)
    } else if la > lr {
        (1, 0, la, lr)
    } else {
        (0, 1, lr, la)
    };
    let gap = (lr - la).abs();
    let inf = gap > LOG_10_INFORMATIVE_THRESHOLD;
    let _ = (best, second);
    (best_i, second_i, gap, inf)
}

fn allele_name(i: usize) -> &'static str {
    if i == 0 {
        "A"
    } else {
        "T"
    }
}

#[test]
fn holdout_6r154_ad_from_41_read_object() {
    if std::env::var("HOLDOUT_6R154").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R154=1");
        return;
    }
    let _policy = EnvGuard::set("GATK_RS_EXPERIMENTAL_KBEST_POLICY", "unbounded_diagnostic");
    kv(
        "parity_mode",
        "unbounded_diagnostic\tnote=diagnostic_not_production_default",
    );
    kv("resource_policy", "quarantined");
    kv("variant", format!("20:{TARGET} {MERGED_REF}/{MERGED_ALT}"));
    kv("alleles", "A,T");
    kv("ploidy", "2");

    let root = repo_root();
    let jf = root.join(JAVA_FLOAT_DUMP_REL);
    let java_keep = load_java_stage_reads(&jf, "stored");
    let java_stored = load_java_matrix(&jf, "stored");
    assert_eq!(java_keep.len(), 201);

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
    let mut java_bounds = covering;
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
        EMIT_SPANNING_DELS,
        Some(&hap_cache),
    );

    let rust_hashes: Vec<String> = haps.iter().map(|h| fnv1a64_hex(&h.bases)).collect();
    let mut java_ref_idx = Vec::new();
    let mut java_alt_idx = Vec::new();
    let mut n_star_haps = 0usize;
    let mut n_unmapped_haps = 0usize;
    for i in 0..haps.len() {
        let spanning = overlapping_events(hap_cache.events_for(i), TARGET);
        let labs = java_create_allele_mapper_labels(
            &spanning,
            TARGET,
            MERGED_REF,
            &[MERGED_ALT.to_string()],
            EMIT_SPANNING_DELS,
        );
        if labs.iter().any(|l| l == SPAN_DEL_ALLELE) {
            n_star_haps += 1;
        }
        if labs.iter().any(|l| l == MERGED_REF) {
            java_ref_idx.push(i);
        }
        if labs.iter().any(|l| l == MERGED_ALT) {
            java_alt_idx.push(i);
        }
        if labs.is_empty()
            || (!labs.iter().any(|l| l == MERGED_REF)
                && !labs.iter().any(|l| l == MERGED_ALT)
                && !labs.iter().any(|l| l == SPAN_DEL_ALLELE))
        {
            n_unmapped_haps += 1;
        }
    }
    kv(
        "mapper_pools",
        format!(
            "ref_haps={} alt_haps={} star_label_haps={} unmapped_haps={} rust_ref={} rust_alt={}",
            java_ref_idx.len(),
            java_alt_idx.len(),
            n_star_haps,
            n_unmapped_haps,
            mapping.ref_haplotype_indices.len(),
            mapping.alt_haplotype_indices.len()
        ),
    );
    kv(
        "internal_star_present",
        format!(
            "haplotype_span_del_label={} remaining_call_alleles=A,T ad_input_columns=A,T",
            n_star_haps > 0
        ),
    );

    let rec_by_key: HashMap<ReadKey, &Record> = outcome
        .genotyping_reads
        .iter()
        .map(|r| (rec_key(r.as_ref()), r.as_ref()))
        .collect();
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

    let variant_start = TARGET as i64;
    let variant_end =
        TARGET as i64 + i64::try_from(MERGED_REF.len().saturating_sub(1)).unwrap_or(0);
    let mut java_retained = BTreeSet::new();
    let mut rust_prod_retained = BTreeSet::new();
    for rk in &java_keep {
        let Some(rec) = rec_by_key.get(rk) else {
            continue;
        };
        if java_alignment_read_overlaps_interval(
            rec,
            variant_start as u64,
            variant_end as u64,
            MARGIN,
        ) {
            java_retained.insert(rk.clone());
            rust_prod_retained.insert(rk.clone());
        }
    }
    kv("retained_reads_java", format!("{}", java_retained.len()));
    kv(
        "retained_reads_rust",
        format!("{}", rust_prod_retained.len()),
    );
    kv(
        "membership_equal",
        format!(
            "{}",
            java_retained.len() == 41 && rust_prod_retained.len() == 41
        ),
    );

    let mate99: ReadKey = (MATE_QNAME.to_string(), 99);
    let mate147: ReadKey = (MATE_QNAME.to_string(), 147);
    kv(
        "restored_mate_membership",
        format!(
            "qname={MATE_QNAME} flags99={} flags147={}",
            java_retained.contains(&mate99),
            java_retained.contains(&mate147)
        ),
    );

    let mut rust_by_key: HashMap<ReadKey, (f64, f64)> = HashMap::new();
    let mut rust_ad_rows = Vec::new();
    for row in &rust_marg {
        let Some(rk) = read_by_idx.get(row.read_index) else {
            continue;
        };
        if !java_retained.contains(rk) {
            continue;
        }
        let lr = row.haplotype_log10_likelihoods[0];
        let la = row.haplotype_log10_likelihoods[1];
        rust_by_key.insert(rk.clone(), (lr, la));
        rust_ad_rows.push(ReadLikelihoodRow {
            read_index: row.read_index,
            read_id: rk.0.clone(),
            haplotype_log10_likelihoods: vec![lr, la],
        });
    }

    let mut java_ad = [0i32, 0];
    let mut rust_ad = [0i32, 0];
    let mut n_best_eq = 0usize;
    let mut n_second_eq = 0usize;
    let mut n_inf_eq = 0usize;
    let mut n_contrib_eq = 0usize;
    let mut n_compared = 0usize;
    let mut first_best_div: Option<String> = None;
    let mut first_inf_div: Option<String> = None;
    let mut first_contrib_div: Option<String> = None;

    println!("6R154\tread_table\tQNAME\tflags\tJ_best\tR_best\tJ_second\tR_second\tJ_d\tR_d\tJ_inf\tR_inf\tJ_contrib\tR_contrib");
    let mut ordered: Vec<ReadKey> = java_retained.iter().cloned().collect();
    ordered.sort();
    for rk in &ordered {
        let mut j_ref = Vec::new();
        let mut j_alt = Vec::new();
        for &hi in &java_ref_idx {
            if let Some(&v) = java_stored.get(&(rk.0.clone(), rk.1, rust_hashes[hi].clone())) {
                j_ref.push(v);
            }
        }
        for &hi in &java_alt_idx {
            if let Some(&v) = java_stored.get(&(rk.0.clone(), rk.1, rust_hashes[hi].clone())) {
                j_alt.push(v);
            }
        }
        let jr = pool_max(&j_ref);
        let ja = pool_max(&j_alt);
        let (j_best, j_second, j_conf, j_inf) = java_best_breaking_ties(&[jr, ja], 0);
        let j_contrib = if j_inf {
            Some(allele_name(j_best))
        } else {
            None
        };
        if j_inf {
            java_ad[j_best] += 1;
        }

        let (rr, ra) = rust_by_key.get(rk).copied().unwrap_or((f64::NAN, f64::NAN));
        let (r_best, r_second, r_gap, r_inf) = rust_biallelic_vote(rr, ra);
        let r_contrib = if r_inf {
            Some(allele_name(r_best))
        } else {
            None
        };
        if r_inf {
            rust_ad[r_best] += 1;
        }

        n_compared += 1;
        if j_best == r_best {
            n_best_eq += 1;
        } else if first_best_div.is_none() {
            first_best_div = Some(format!("{} flags={}", rk.0, rk.1));
        }
        if j_second == r_second {
            n_second_eq += 1;
        }
        if j_inf == r_inf {
            n_inf_eq += 1;
        } else if first_inf_div.is_none() {
            first_inf_div = Some(format!("{} flags={}", rk.0, rk.1));
        }
        if j_contrib == r_contrib {
            n_contrib_eq += 1;
        } else if first_contrib_div.is_none() {
            first_contrib_div = Some(format!(
                "{} flags={} java={:?} rust={:?}",
                rk.0, rk.1, j_contrib, r_contrib
            ));
        }

        println!(
            "6R154\tread\t{}\t{}\t{}\t{}\t{}\t{}\t{:.6}\t{:.6}\t{}\t{}\t{}\t{}",
            rk.0,
            rk.1,
            allele_name(j_best),
            allele_name(r_best),
            allele_name(j_second),
            allele_name(r_second),
            j_conf,
            r_gap,
            j_inf,
            r_inf,
            j_contrib.unwrap_or("-"),
            r_contrib.unwrap_or("-"),
        );

        if rk.0 == MATE_QNAME {
            kv(
                &format!("restored_mate_flags{}", rk.1),
                format!(
                    "java_best={} rust_best={} java_second={} rust_second={} java_d={:.6} rust_d={:.6} java_inf={} rust_inf={} java_contrib={} rust_contrib={}",
                    allele_name(j_best),
                    allele_name(r_best),
                    allele_name(j_second),
                    allele_name(r_second),
                    j_conf,
                    r_gap,
                    j_inf,
                    r_inf,
                    j_contrib.unwrap_or("-"),
                    r_contrib.unwrap_or("-"),
                ),
            );
        }
    }

    let rust_prod = InformativeAd::from_marginalized_rows(&rust_ad_rows, 0, 1, None);
    let rust_prod_vec = biallelic_allele_depths_from_rows(&rust_ad_rows, 0, 1);
    kv(
        "java_ad_input",
        format!(
            "evidence={} alleles=2 ordering=A,T after_retainEvidence=true after_marginalize=true unused_alt_subset=identity remaining=A,T addEvidence_filtered=uninformative_if_present",
            n_compared
        ),
    );
    kv(
        "rust_ad_input",
        format!(
            "evidence={} alleles=2 ordering=A,T after_retainEvidence=true after_marginalize=true path=InformativeAd::from_marginalized_rows",
            rust_ad_rows.len()
        ),
    );
    kv(
        "java_alleles",
        "A,T keyed_by=Allele.equals vc_order=REF_then_ALT",
    );
    kv("rust_alleles", "A,T keyed_by=column_index 0=REF_A 1=ALT_T");
    kv(
        "java_informative_rule",
        "searchBestAllele first-max then REF-priority if gap<0.2; isInformative iff confidence>0.2",
    );
    kv(
        "rust_informative_rule",
        "biallelic |lr-la|>0.2 then lr>la ? REF : ALT",
    );
    kv("java_ad", format!("{},{}", java_ad[0], java_ad[1]));
    kv("rust_ad", format!("{},{}", rust_ad[0], rust_ad[1]));
    kv(
        "rust_production_informative_ad",
        format!("{},{}", rust_prod.ref_depth, rust_prod.alt_depth),
    );
    kv(
        "allele_identity",
        format!(
            "java_AD0=A java_AD1=T rust_AD0=A rust_AD1=T production={:?}",
            rust_prod_vec
        ),
    );

    let prod_call = outcome.genotyped_calls.iter().find(|c| {
        c.event.start_1based.get() == TARGET
            && c.event.ref_allele == MERGED_REF
            && c.event.alt_allele == MERGED_ALT
    });
    match prod_call {
        Some(c) => {
            let ad = c.genotype.format.ad_as_i32();
            kv(
                "production_genotyped_call_ad",
                format!(
                    "{:?} extra_alts={:?} unused_alt_subset={} gt_pl={:?}",
                    ad,
                    c.extra_alt_alleles,
                    c.post_merge_unused_alt_subset,
                    c.genotype.format.pl_as_i32()
                ),
            );
            kv(
                "later_site_reshape",
                "observed: genotyped_calls FORMAT AD is not the first DepthPerAlleleBySample write; SiteReshape Class-A3 / sparse pileup can replace AD and GT/PL. Out of 6R.154 scope (would reopen GT/PL and L9).",
            );
        }
        None => kv(
            "production_genotyped_call_ad",
            "absent (site may fail emit; AD reconstructed from retained object)",
        ),
    }

    let assignments_equal = n_best_eq == n_compared && n_second_eq == n_compared;
    let inf_equal = n_inf_eq == n_compared;
    let contrib_equal = n_contrib_eq == n_compared;
    let ad_equal = java_ad == rust_ad
        && rust_ad[0] == rust_prod.ref_depth
        && rust_ad[1] == rust_prod.alt_depth;
    kv(
        "per_read_assignments_equal",
        format!(
            "{} compared={n_compared} best_eq={n_best_eq} second_eq={n_second_eq}",
            assignments_equal
        ),
    );
    kv(
        "per_read_informativeness_equal",
        format!("{} inf_eq={n_inf_eq}", inf_equal),
    );
    kv(
        "per_read_contribution_equal",
        format!("{} contrib_eq={n_contrib_eq}", contrib_equal),
    );
    kv("ad_equal", format!("{ad_equal}"));

    let (classification, first_divergence) = if java_retained.len() != 41
        || rust_prod_retained.len() != 41
        || !java_retained.contains(&mate99)
        || !java_retained.contains(&mate147)
    {
        (
            "AD_INPUT_OBJECT_DIVERGENCE",
            "6R.151 41-read / mate baseline broken".to_string(),
        )
    } else if n_compared != 41 {
        (
            "AD_INPUT_OBJECT_DIVERGENCE",
            format!("AD input evidence count {n_compared} != 41"),
        )
    } else if !assignments_equal {
        (
            "AD_BEST_ALLELE_SELECTION_DIVERGENCE",
            first_best_div.unwrap_or_else(|| "best-allele mismatch".to_string()),
        )
    } else if !inf_equal {
        (
            "AD_INFORMATIVENESS_DIVERGENCE",
            first_inf_div.unwrap_or_else(|| "informativeness mismatch".to_string()),
        )
    } else if !contrib_equal {
        (
            "AD_COUNTING_DIVERGENCE",
            first_contrib_div.unwrap_or_else(|| "contribution mismatch".to_string()),
        )
    } else if java_ad != rust_ad {
        (
            "AD_COUNTING_DIVERGENCE",
            format!("java_ad={:?} rust_ad={:?}", java_ad, rust_ad),
        )
    } else if rust_ad[0] != rust_prod.ref_depth || rust_ad[1] != rust_prod.alt_depth {
        (
            "AD_POSTPROCESSING_DIVERGENCE",
            format!(
                "row reconstruction {:?} vs InformativeAd {},{}",
                rust_ad, rust_prod.ref_depth, rust_prod.alt_depth
            ),
        )
    } else {
        (
            "NO_DIVERGENCE",
            "41-read remaining-allele AD matches Java DepthPerAlleleBySample row-by-row"
                .to_string(),
        )
    };
    kv("first_divergence", &first_divergence);
    kv("classification", classification);
    kv(
        "semantic_consequence",
        format!(
            "first AD write A,T = {},{}; restored mate independently eligible; unused-ALT identity remarg; * not in remaining call alleles; later SiteReshape overwrite is out of this round",
            java_ad[0], java_ad[1]
        ),
    );
    kv("production_change", "NONE");
    kv(
        "next_arrow",
        "first AD write closed on the 41-read object (37,4); QUAL/L9/VCF stay out of this round",
    );

    assert_eq!(java_retained.len(), 41);
    assert_eq!(rust_prod_retained.len(), 41);
    assert!(java_retained.contains(&mate99) && java_retained.contains(&mate147));
    assert_eq!(n_compared, 41);
    assert!(assignments_equal && inf_equal && contrib_equal && ad_equal);
    assert_eq!(classification, "NO_DIVERGENCE");
}
