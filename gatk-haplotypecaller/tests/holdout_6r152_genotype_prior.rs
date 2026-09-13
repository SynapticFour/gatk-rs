//! 6R.152: genotype-prior relative order after the 6R.151 41-read evidence object.
//! Diagnostic-only `unbounded_diagnostic`. Skipped unless `HOLDOUT_6R152=1`.
//!
//! ```text
//! HOLDOUT_6R152=1 GATK_RS_EXPERIMENTAL_KBEST_POLICY=unbounded_diagnostic \
//!   cargo test -p gatk-haplotypecaller --test holdout_6r152_genotype_prior -- --nocapture --test-threads=1
//! ```

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::activity_scoring::{
    DEFAULT_INDEL_HETEROZYGOSITY, DEFAULT_SNP_HETEROZYGOSITY,
};
use gatk_haplotypecaller::event_map::{
    build_per_haplotype_variation_events, overlapping_events, VariationEvent,
};
use gatk_haplotypecaller::genotyping::{
    best_pl_index, biallelic_diploid_log10_priors, diploid_genotype_alleles_from_pl_index,
    emit_genotype_format_fields, genotype_posteriors_from_log10_likelihoods,
    BiallelicDiploidPriorModel, ReadLikelihoodRow,
};
use gatk_haplotypecaller::hc_allele_mapping::{create_allele_mapper_with_events, SPAN_DEL_ALLELE};
use gatk_haplotypecaller::hc_genotyping_engine::java_alignment_read_overlaps_interval;
use gatk_haplotypecaller::read_unclip::{alignment_end_1based, gatk_soft_start_1based};
use gatk_haplotypecaller::{
    biallelic_genotype_log10_likelihoods_gatk, call_disposition, flatten_assembly_regions,
    marginalize_rows_to_biallelic_alleles, region_likelihoods_to_rows,
    traverse_assembly_region_walker, AssemblyRegionCallDisposition, CallRegionArgs, GenomePosition,
    HaplotypeCallerEngine, ReadFilterParams, WalkerTraversalConfig,
    DEFAULT_INFORMATIVE_READ_OVERLAP_MARGIN, DEFAULT_STAND_EMIT_CONFIDENCE,
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

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    println!("6R152\t{key}\t{}", value.as_ref());
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

fn load_java_orig_coords(path: &Path) -> HashMap<ReadKey, (i64, i64, i64, i64)> {
    let text = std::fs::read_to_string(path).expect("java dump");
    let mut out = HashMap::new();
    for line in text.lines() {
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 4 || parts[0] != "6R143" || parts[1] != "orig" {
            continue;
        }
        let kv = parse_kv_fields(&parts[3..]);
        let qname = kv.get("qname").cloned().unwrap_or_default();
        let flags: u16 = kv.get("flags").and_then(|s| s.parse().ok()).unwrap_or(0);
        let start = kv.get("start").and_then(|s| s.parse().ok()).unwrap_or(0);
        let end = kv.get("end").and_then(|s| s.parse().ok()).unwrap_or(0);
        let u_start = kv.get("uStart").and_then(|s| s.parse().ok()).unwrap_or(0);
        let u_end = kv.get("uEnd").and_then(|s| s.parse().ok()).unwrap_or(0);
        out.insert((qname, flags), (start, end, u_start, u_end));
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

fn load_java_hap_hashes(path: &Path) -> Vec<String> {
    let text = std::fs::read_to_string(path).expect("java dump");
    let mut out = Vec::new();
    for line in text.lines() {
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 4 || parts[0] != "6R143" || parts[1] != "hap" {
            continue;
        }
        out.push(parts[3].to_string());
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

fn java_expand_within_contig(start: i64, end: i64, padding: i32, contig_len: i64) -> (i64, i64) {
    let p = i64::from(padding.max(0));
    ((start - p).max(1), (end + p).min(contig_len))
}

fn java_assuming_hw_snp_log10_priors(snp_heterozygosity: f64) -> [f64; 3] {
    let snp_het = snp_heterozygosity.log10();
    let log10_snp_norm = 3.0_f64.log10();
    [
        0.0,
        snp_het - log10_snp_norm,
        snp_het * 2.0 - log10_snp_norm,
    ]
}

fn rec_key(rec: &Record) -> ReadKey {
    (
        String::from_utf8_lossy(rec.qname()).into_owned(),
        rec.flags(),
    )
}

fn rec_align(rec: &Record) -> (i64, i64) {
    (
        (rec.pos() + 1).max(1),
        i64::from(alignment_end_1based(rec).max(1)),
    )
}

fn rec_soft(rec: &Record) -> (i64, i64) {
    (
        gatk_soft_start_1based(rec).max(1),
        i64::from(alignment_end_1based(rec).max(1)),
    )
}

fn collapse_qname(
    members: &BTreeSet<ReadKey>,
    max_ll: &HashMap<ReadKey, f64>,
    index: &HashMap<ReadKey, usize>,
) -> BTreeSet<ReadKey> {
    let mut best: HashMap<String, ReadKey> = HashMap::new();
    for k in members {
        let ll = max_ll.get(k).copied().unwrap_or(f64::NEG_INFINITY);
        let idx = index.get(k).copied().unwrap_or(usize::MAX);
        best.entry(k.0.clone())
            .and_modify(|cur| {
                let cur_ll = max_ll.get(cur).copied().unwrap_or(f64::NEG_INFINITY);
                let cur_idx = index.get(cur).copied().unwrap_or(usize::MAX);
                if ll > cur_ll || (ll == cur_ll && idx < cur_idx) {
                    *cur = k.clone();
                }
            })
            .or_insert_with(|| k.clone());
    }
    best.into_values().collect()
}

fn gls_for_keys(
    keys: &BTreeSet<ReadKey>,
    rust_marg: &[ReadLikelihoodRow],
    read_by_idx: &[ReadKey],
    java_keep: &BTreeSet<ReadKey>,
    java_stored: &HashMap<(String, u16, String), f64>,
    java_ref_idx: &[usize],
    java_alt_idx: &[usize],
    rust_hashes: &[String],
) -> (Vec<f64>, Vec<f64>) {
    let mut j_rows = Vec::new();
    let mut r_rows = Vec::new();
    for row in rust_marg {
        let Some(rk) = read_by_idx.get(row.read_index) else {
            continue;
        };
        if !java_keep.contains(rk) || !keys.contains(rk) {
            continue;
        }
        let mut j_ref = Vec::new();
        let mut j_alt = Vec::new();
        for &hi in java_ref_idx {
            if let Some(&v) = java_stored.get(&(rk.0.clone(), rk.1, rust_hashes[hi].clone())) {
                j_ref.push(v);
            }
        }
        for &hi in java_alt_idx {
            if let Some(&v) = java_stored.get(&(rk.0.clone(), rk.1, rust_hashes[hi].clone())) {
                j_alt.push(v);
            }
        }
        let jr = pool_max(&j_ref);
        let ja = pool_max(&j_alt);
        if jr.is_finite() || ja.is_finite() {
            j_rows.push(ReadLikelihoodRow {
                read_index: row.read_index,
                read_id: String::new(),
                haplotype_log10_likelihoods: vec![jr, ja],
            });
        }
        r_rows.push(ReadLikelihoodRow {
            read_index: row.read_index,
            read_id: String::new(),
            haplotype_log10_likelihoods: row.haplotype_log10_likelihoods.clone(),
        });
    }
    let j = if j_rows.is_empty() {
        vec![0.0, 0.0, 0.0]
    } else {
        biallelic_genotype_log10_likelihoods_gatk(&j_rows, 0, 1)
    };
    let r = if r_rows.is_empty() {
        vec![0.0, 0.0, 0.0]
    } else {
        biallelic_genotype_log10_likelihoods_gatk(&r_rows, 0, 1)
    };
    (j, r)
}

fn fmt_gl(gl: &[f64]) -> String {
    format!("0/0={:.8} 0/1={:.8} 1/1={:.8}", gl[0], gl[1], gl[2])
}

fn fmt_arr3(v: &[f64; 3]) -> String {
    format!("0/0={:.8} 0/1={:.8} 1/1={:.8}", v[0], v[1], v[2])
}

fn winner_idx(gl: &[f64]) -> usize {
    gl.iter()
        .enumerate()
        .max_by(|a, b| a.1.total_cmp(b.1))
        .map(|(i, _)| i)
        .unwrap_or(0)
}

fn winner_margin(gl: &[f64]) -> f64 {
    let w = winner_idx(gl);
    let best = gl[w];
    gl.iter()
        .enumerate()
        .filter(|(i, _)| *i != w)
        .map(|(_, v)| (best - *v).abs())
        .fold(f64::INFINITY, f64::min)
}

#[test]
fn holdout_6r152_genotype_prior_relative_order() {
    if std::env::var("HOLDOUT_6R152").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R152=1");
        return;
    }
    let _policy = EnvGuard::set("GATK_RS_EXPERIMENTAL_KBEST_POLICY", "unbounded_diagnostic");
    kv(
        "parity_mode",
        "unbounded_diagnostic\tnote=diagnostic_not_production_default",
    );
    kv("resource_policy", "quarantined");
    kv("variant", format!("20:{TARGET} {MERGED_REF}/{MERGED_ALT}"));

    let root = repo_root();
    let jf = root.join(JAVA_FLOAT_DUMP_REL);
    let java_keep = load_java_stage_reads(&jf, "stored");
    let java_orig_coords = load_java_orig_coords(&jf);
    let java_stored = load_java_matrix(&jf, "stored");
    let java_hap_hashes = load_java_hap_hashes(&jf);
    assert_eq!(java_keep.len(), 201);
    assert_eq!(java_hap_hashes.len(), 25);

    let ref_fasta = root.join(REF_REL);
    let bam = root.join(BAM_REL);
    let dict = SequenceDictionary::from_fasta_path(&ref_fasta).expect("dict");
    let contig_len = dict
        .contig("20")
        .map(|c| c.length as i64)
        .unwrap_or(i64::MAX);
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

    let variant_start = TARGET as i64;
    let variant_end =
        TARGET as i64 + i64::try_from(MERGED_REF.len().saturating_sub(1)).unwrap_or(0);
    let (exp_start, exp_end) =
        java_expand_within_contig(variant_start, variant_end, MARGIN, contig_len);
    kv(
        "retain_interval_java",
        format!(
            "contig=20 variant_start={variant_start} variant_end={variant_end} margin={MARGIN} expanded_start={exp_start} expanded_end={exp_end} contig_len={contig_len} source=SimpleInterval(mergedVC).expandWithinContig then target.overlaps(read)"
        ),
    );
    kv(
        "retain_interval_rust",
        format!(
            "contig=20 variant_start={variant_start} variant_end={variant_end} margin={MARGIN} expanded_start={exp_start} expanded_end={exp_end} predicate=java_alignment_read_overlaps_interval"
        ),
    );

    let rust_hashes: Vec<String> = haps.iter().map(|h| fnv1a64_hex(&h.bases)).collect();
    let mut java_ref_idx = Vec::new();
    let mut java_alt_idx = Vec::new();
    for i in 0..haps.len() {
        let spanning = overlapping_events(hap_cache.events_for(i), TARGET);
        let labs = java_create_allele_mapper_labels(
            &spanning,
            TARGET,
            MERGED_REF,
            &[MERGED_ALT.to_string()],
            EMIT_SPANNING_DELS,
        );
        if labs.iter().any(|l| l == MERGED_REF) {
            java_ref_idx.push(i);
        }
        if labs.iter().any(|l| l == MERGED_ALT) {
            java_alt_idx.push(i);
        }
    }

    let rec_by_key: HashMap<ReadKey, &Record> = outcome
        .genotyping_reads
        .iter()
        .map(|r| (rec_key(r.as_ref()), r.as_ref()))
        .collect();
    let index_by_key: HashMap<ReadKey, usize> = outcome
        .genotyping_reads
        .iter()
        .enumerate()
        .map(|(i, r)| (rec_key(r.as_ref()), i))
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

    let mut max_ll: HashMap<ReadKey, f64> = HashMap::new();
    for row in &rust_rows {
        let Some(rk) = read_by_idx.get(row.read_index) else {
            continue;
        };
        let m = row
            .haplotype_log10_likelihoods
            .iter()
            .copied()
            .fold(f64::NEG_INFINITY, f64::max);
        max_ll
            .entry(rk.clone())
            .and_modify(|e| {
                if m > *e {
                    *e = m;
                }
            })
            .or_insert(m);
    }

    let mut java_retained = BTreeSet::new();
    let mut rust_pred_retained = BTreeSet::new();
    let mut orig_aln_retained = BTreeSet::new();
    let mut orig_soft_retained = BTreeSet::new();
    let mut n_missing_post = 0usize;
    for rk in &java_keep {
        if let Some(&(s, e, us, ue)) = java_orig_coords.get(rk) {
            if exp_start <= e && s <= exp_end {
                orig_aln_retained.insert(rk.clone());
            }
            if exp_start <= ue && us <= exp_end {
                orig_soft_retained.insert(rk.clone());
            }
        }
        let Some(rec) = rec_by_key.get(rk) else {
            n_missing_post += 1;
            continue;
        };
        let java_ov = java_alignment_read_overlaps_interval(
            rec,
            variant_start as u64,
            variant_end as u64,
            MARGIN,
        );
        let rust_ov = java_ov;
        if java_ov {
            java_retained.insert(rk.clone());
        }
        if rust_ov {
            rust_pred_retained.insert(rk.clone());
        }
    }
    let helper_collapsed = collapse_qname(&rust_pred_retained, &max_ll, &index_by_key);
    // 6R.151 production: Java-strict retainEvidence is overlap-only.
    let rust_prod_retained = rust_pred_retained.clone();

    kv(
        "retain_reads_java",
        format!(
            "input=201 retained={} orig_aln={} orig_soft={} missing_post_realign={n_missing_post}",
            java_retained.len(),
            orig_aln_retained.len(),
            orig_soft_retained.len()
        ),
    );
    kv(
        "retain_reads_rust",
        format!(
            "input=201 predicate={} helper_qname_collapse={} production={}",
            rust_pred_retained.len(),
            helper_collapsed.len(),
            rust_prod_retained.len()
        ),
    );

    let common: BTreeSet<_> = java_retained
        .intersection(&rust_pred_retained)
        .cloned()
        .collect();
    let java_only: BTreeSet<_> = java_retained
        .difference(&rust_pred_retained)
        .cloned()
        .collect();
    let rust_only: BTreeSet<_> = rust_pred_retained
        .difference(&java_retained)
        .cloned()
        .collect();
    let prod_java_only: BTreeSet<_> = java_retained
        .difference(&rust_prod_retained)
        .cloned()
        .collect();
    let prod_rust_only: BTreeSet<_> = rust_prod_retained
        .difference(&java_retained)
        .cloned()
        .collect();
    kv(
        "retain_membership_predicate",
        format!(
            "COMMON={} JAVA_ONLY={} RUST_ONLY={}",
            common.len(),
            java_only.len(),
            rust_only.len()
        ),
    );
    kv(
        "retain_membership_production",
        format!(
            "COMMON={} JAVA_ONLY={} RUST_ONLY={} (production = overlap, no QNAME collapse)",
            java_retained.intersection(&rust_prod_retained).count(),
            prod_java_only.len(),
            prod_rust_only.len()
        ),
    );
    let helper_java_only: BTreeSet<_> = java_retained
        .difference(&helper_collapsed)
        .cloned()
        .collect();
    kv(
        "qname_dedupe_helper",
        format!(
            "present_as_helper=true applied_on_java_strict_path=false overlap={} collapsed={} JAVA_ONLY_if_applied={}",
            rust_pred_retained.len(),
            helper_collapsed.len(),
            helper_java_only.len()
        ),
    );
    const DROPPED_QNAME: &str = "HWI-D00360:7:H88WKADXX:2:2107:6787:30989";
    let mate_first: ReadKey = (DROPPED_QNAME.to_string(), 99);
    let mate_keys: Vec<_> = java_retained
        .iter()
        .filter(|(q, _)| q == DROPPED_QNAME)
        .cloned()
        .collect();
    kv(
        "paired_read_semantics",
        format!(
            "qname={DROPPED_QNAME} java_mates={} flags={:?} flags99_in_java={} flags99_in_helper_collapse={} flags99_in_production={}",
            mate_keys.len(),
            mate_keys.iter().map(|k| k.1).collect::<Vec<_>>(),
            java_retained.contains(&mate_first),
            helper_collapsed.contains(&mate_first),
            rust_prod_retained.contains(&mate_first)
        ),
    );

    let mut n_asymm = 0usize;
    for rk in prod_java_only.iter().chain(prod_rust_only.iter()) {
        n_asymm += 1;
        let (js, je) = rec_by_key.get(rk).map(|r| rec_align(r)).unwrap_or((0, 0));
        let (ss, se) = rec_by_key.get(rk).map(|r| rec_soft(r)).unwrap_or((0, 0));
        let (os, oe, _, _) = java_orig_coords.get(rk).copied().unwrap_or((0, 0, 0, 0));
        let java_dec = java_retained.contains(rk);
        let rust_pred = rust_pred_retained.contains(rk);
        let rust_prod = rust_prod_retained.contains(rk);
        let reason = if rust_pred && !helper_collapsed.contains(rk) && rust_prod {
            "helper would drop overlapping mate; production keep matches Java"
        } else if rust_pred && !rust_prod {
            "Rust QNAME collapse dropped overlapping mate; Java retainEvidence keeps both"
        } else if java_dec != rust_pred {
            "alignment overlap predicate mismatch"
        } else {
            "asymmetric"
        };
        kv(
            "asymmetric_read",
            format!(
                "QNAME={} flags={} post_realign={js}-{je} orig={os}-{oe} soft={ss}-{se} java_overlap={java_dec} rust_predicate={rust_pred} rust_production={rust_prod} reason={reason}",
                rk.0, rk.1
            ),
        );
        if n_asymm >= 32 {
            break;
        }
    }

    let predicate_equal = java_only.is_empty() && rust_only.is_empty();
    let production_equal = prod_java_only.is_empty() && prod_rust_only.is_empty();
    kv(
        "retain_membership_equal",
        format!("predicate={predicate_equal} production={production_equal}"),
    );

    // Retained likelihood object on the Java-equivalent (no QNAME collapse) set.
    let mut n_pairs = 0usize;
    let mut exact_equal = 0usize;
    let mut max_abs = 0.0_f64;
    let mut sum_abs = 0.0_f64;
    for row in &rust_marg {
        let Some(rk) = read_by_idx.get(row.read_index) else {
            continue;
        };
        if !java_retained.contains(rk) {
            continue;
        }
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
        let rr = row.haplotype_log10_likelihoods[0];
        let ra = row.haplotype_log10_likelihoods[1];
        for (j, r) in [(jr, rr), (ja, ra)] {
            if !j.is_finite() {
                continue;
            }
            n_pairs += 1;
            let d = (j - r).abs();
            max_abs = max_abs.max(d);
            sum_abs += d;
            if j == r {
                exact_equal += 1;
            }
        }
    }
    let mean_abs = if n_pairs == 0 {
        0.0
    } else {
        sum_abs / n_pairs as f64
    };
    let cells_equal = max_abs < 1e-4;
    kv(
        "retained_likelihood_object",
        format!(
            "alleles=A,T pairs={n_pairs} exact_equal={exact_equal} max_abs_delta={max_abs:.6e} mean_abs_delta={mean_abs:.6e} ordering=REF_then_ALT"
        ),
    );
    kv("retained_likelihood_object_equal", format!("{cells_equal}"));
    kv("max_abs_delta", format!("{max_abs:.6e}"));

    let (java_gl_retain, rust_gl_retain) = gls_for_keys(
        &java_retained,
        &rust_marg,
        &read_by_idx,
        &java_keep,
        &java_stored,
        &java_ref_idx,
        &java_alt_idx,
        &rust_hashes,
    );
    let (java_gl_prod, rust_gl_prod) = gls_for_keys(
        &rust_prod_retained,
        &rust_marg,
        &read_by_idx,
        &java_keep,
        &java_stored,
        &java_ref_idx,
        &java_alt_idx,
        &rust_hashes,
    );
    kv("raw_gl_retained_java", fmt_gl(&java_gl_retain));
    kv("raw_gl_retained_rust", fmt_gl(&rust_gl_retain));
    kv("raw_gl_production_java", fmt_gl(&java_gl_prod));
    kv("raw_gl_production_rust", fmt_gl(&rust_gl_prod));
    let (java_gl_helper, rust_gl_helper) = gls_for_keys(
        &helper_collapsed,
        &rust_marg,
        &read_by_idx,
        &java_keep,
        &java_stored,
        &java_ref_idx,
        &java_alt_idx,
        &rust_hashes,
    );
    kv("raw_gl_helper_collapse_java_ll", fmt_gl(&java_gl_helper));
    kv("raw_gl_helper_collapse_rust_ll", fmt_gl(&rust_gl_helper));
    let gl_retain_delta = java_gl_retain
        .iter()
        .zip(rust_gl_retain.iter())
        .map(|(a, b)| (a - b).abs())
        .fold(0.0_f64, f64::max);
    kv(
        "raw_gl_equal_within_residual",
        format!("{}", gl_retain_delta < 1e-3),
    );

    let sample_name = outcome
        .genotyping_reads
        .first()
        .and_then(|r| match r.as_ref().aux(b"RG") {
            Ok(rust_htslib::bam::record::Aux::String(s)) => Some(s.to_string()),
            _ => None,
        })
        .unwrap_or_else(|| "UNKNOWN_RG".to_string());
    kv(
        "calculate_genotypes_input",
        format!(
            "sample_count=1 sample={sample_name} ploidy=2 alleles=A,T genotyping_model=IndependentSampleGenotypesModel assignment=USE_PLS_TO_ASSIGN gt_at_entry=NO_CALL evidence_count={} representation=log10_GL_as_PL_on_NO_CALL",
            java_retained.len()
        ),
    );
    kv("ploidy", "2");
    kv("sample_count", "1");
    kv("genotyping_model", "IndependentSampleGenotypesModel");
    kv(
        "stand_call_confidence",
        format!(
            "java=30 rust_emit_default={DEFAULT_STAND_EMIT_CONFIDENCE} role=passesEmitThreshold_only_not_GT_PL"
        ),
    );
    kv(
        "java_prior_path",
        "calculateGLsForThisEvent (no prior) → resolveGenotypePriorCalculator (dragstrParams==null → assumingHW(log10(snpHet), log10(indelHet))) → calculateGenotypes → makeGenotypeCall USE_PLS_TO_ASSIGN (gpc unused); USE_POSTERIOR_PROBABILITIES would ebeAdd(getLog10Priors, GLs)",
    );
    kv(
        "rust_prior_path",
        "biallelic_diploid_log10_priors(BiallelicDiploidPriorModel::default) remainder-mass simplex; genotype_posteriors_from_log10_likelihoods = GL+prior then linear-normalize; production discards _posterior and emits from raw GLs",
    );
    kv(
        "prior_parameters_java",
        format!(
            "snpHeterozygosity={DEFAULT_SNP_HETEROZYGOSITY} source=HomoSapiensConstants consumed=assumingHW_snpHet; indelHeterozygosity={DEFAULT_INDEL_HETEROZYGOSITY} unused_on_SNP_A_T; ploidy=2 allele_count=2 sample_count_unused dragstrParams=null"
        ),
    );
    kv(
        "prior_parameters_rust",
        "het_prior=1e-3 hom_var_prior=5e-4 hom_ref=1-het-hom_var=0.9985 source=BiallelicDiploidPriorModel::default consumed=biallelic_diploid_log10_priors",
    );

    let prior_java = java_assuming_hw_snp_log10_priors(DEFAULT_SNP_HETEROZYGOSITY);
    let prior_rust =
        biallelic_diploid_log10_priors(BiallelicDiploidPriorModel::default()).expect("rust priors");
    kv("prior_vector_java", fmt_arr3(&prior_java));
    kv(
        "prior_vector_rust",
        format!(
            "0/0={:.8} 0/1={:.8} 1/1={:.8} indel_het={DEFAULT_INDEL_HETEROZYGOSITY}",
            prior_rust[0], prior_rust[1], prior_rust[2]
        ),
    );
    let prior_max = prior_java
        .iter()
        .zip(prior_rust.iter())
        .map(|(a, b)| (a - b).abs())
        .fold(0.0_f64, f64::max);
    let prior_equal = prior_max == 0.0;
    kv(
        "prior_equal",
        format!("exact={prior_equal} max_abs_delta={prior_max:.6e}"),
    );
    let java_rel = [
        prior_java[0] - prior_java[1],
        prior_java[0] - prior_java[2],
        prior_java[1] - prior_java[2],
    ];
    let rust_rel = [
        prior_rust[0] - prior_rust[1],
        prior_rust[0] - prior_rust[2],
        prior_rust[1] - prior_rust[2],
    ];
    let rel_max = java_rel
        .iter()
        .zip(rust_rel.iter())
        .map(|(a, b)| (a - b).abs())
        .fold(0.0_f64, f64::max);
    let relative_prior_ratios_equal = rel_max < 1e-6;
    kv(
        "relative_prior_log_diffs",
        format!(
            "java_00_vs_01={:.8} java_00_vs_11={:.8} java_01_vs_11={:.8} rust_00_vs_01={:.8} rust_00_vs_11={:.8} rust_01_vs_11={:.8} max_abs_delta={rel_max:.6e}",
            java_rel[0], java_rel[1], java_rel[2], rust_rel[0], rust_rel[1], rust_rel[2]
        ),
    );
    kv(
        "prior_normalization_equivalent",
        format!("{}", relative_prior_ratios_equal),
    );
    kv(
        "relative_prior_order_equal",
        format!(
            "{}",
            java_rel
                .iter()
                .zip(rust_rel.iter())
                .all(|(j, r)| j.signum() == r.signum() || (*j).abs() < 1e-12 && (*r).abs() < 1e-12)
        ),
    );
    kv(
        "relative_prior_ratios_equal",
        format!("{relative_prior_ratios_equal}"),
    );

    let post_java: [f64; 3] = [
        java_gl_retain[0] + prior_java[0],
        java_gl_retain[1] + prior_java[1],
        java_gl_retain[2] + prior_java[2],
    ];
    let rust_post_obj =
        genotype_posteriors_from_log10_likelihoods(&rust_gl_retain, &prior_rust).expect("post");
    kv("posterior_gl_java", fmt_arr3(&post_java));
    kv(
        "posterior_gl_rust",
        fmt_gl(&rust_post_obj.genotype_log10_posteriors),
    );
    let post_max = post_java
        .iter()
        .zip(rust_post_obj.genotype_log10_posteriors.iter())
        .map(|(a, b)| (a - b).abs())
        .fold(0.0_f64, f64::max);
    kv(
        "posterior_equal",
        format!("max_abs_delta={post_max:.6e} java_op=GL_plus_prior rust_op=GL_plus_prior (discarded in production USE_PLS)"),
    );

    let cf_a = genotype_posteriors_from_log10_likelihoods(&rust_gl_retain, &prior_java)
        .expect("cfA rust GL + java prior");
    let cf_b_post: [f64; 3] = [
        java_gl_retain[0] + prior_rust[0],
        java_gl_retain[1] + prior_rust[1],
        java_gl_retain[2] + prior_rust[2],
    ];
    kv(
        "counterfactual_A_rust_gl_java_prior",
        format!(
            "winner={} {}",
            winner_idx(&cf_a.genotype_log10_posteriors),
            fmt_gl(&cf_a.genotype_log10_posteriors)
        ),
    );
    kv(
        "counterfactual_B_java_gl_rust_prior",
        format!("winner={} {}", winner_idx(&cf_b_post), fmt_arr3(&cf_b_post)),
    );

    let java_gt_from_gl = winner_idx(&java_gl_retain);
    let rust_gt_from_gl = winner_idx(&rust_gl_retain);
    let java_gt_from_post = winner_idx(&post_java);
    let rust_gt_from_post = rust_post_obj.most_likely_genotype_index;
    let java_pl = emit_genotype_format_fields(&java_gl_retain, &[0, 0])
        .expect("java pl")
        .pl_as_i32();
    let rust_pl = emit_genotype_format_fields(&rust_gl_retain, &[0, 0])
        .expect("rust pl")
        .pl_as_i32();
    let java_gt_alleles = diploid_genotype_alleles_from_pl_index(
        2,
        best_pl_index(
            &java_pl
                .iter()
                .copied()
                .map(gatk_haplotypecaller::bio_ids::PhredLikelihood::from_i32_saturating)
                .collect::<Vec<_>>(),
        ),
    );
    let rust_gt_alleles = diploid_genotype_alleles_from_pl_index(
        2,
        best_pl_index(
            &rust_pl
                .iter()
                .copied()
                .map(gatk_haplotypecaller::bio_ids::PhredLikelihood::from_i32_saturating)
                .collect::<Vec<_>>(),
        ),
    );
    kv(
        "GT_from_USE_PLS",
        format!(
            "java={}/{} rust={}/{} gl_winner_java={java_gt_from_gl} gl_winner_rust={rust_gt_from_gl} post_winner_java={java_gt_from_post} post_winner_rust={rust_gt_from_post} margin_java={:.6e} margin_rust={:.6e}",
            java_gt_alleles[0],
            java_gt_alleles[1],
            rust_gt_alleles[0],
            rust_gt_alleles[1],
            winner_margin(&java_gl_retain),
            winner_margin(&rust_gl_retain)
        ),
    );
    kv(
        "PL_from_raw_GL",
        format!(
            "java={},{},{} rust={},{},{}",
            java_pl[0], java_pl[1], java_pl[2], rust_pl[0], rust_pl[1], rust_pl[2]
        ),
    );
    let gt_equal = java_gt_alleles == rust_gt_alleles;
    let pl_equal = java_pl == rust_pl;
    kv(
        "GT_equal",
        format!("{gt_equal} (USE_PLS_TO_ASSIGN on retained GLs; priors unused)"),
    );
    kv(
        "PL_equal",
        format!("{pl_equal} (from raw GLs, not posteriors)"),
    );
    kv(
        "comparison_table",
        format!(
            "GT\trawJ\trawR\tdRaw\tpriorJ\tpriorR\tdPrior\tpostJ\tpostR\tdPost | 0/0\t{:.12}\t{:.12}\t{:.6e}\t{:.12}\t{:.12}\t{:.6e}\t{:.12}\t{:.12}\t{:.6e} | 0/1\t{:.12}\t{:.12}\t{:.6e}\t{:.12}\t{:.12}\t{:.6e}\t{:.12}\t{:.12}\t{:.6e} | 1/1\t{:.12}\t{:.12}\t{:.6e}\t{:.12}\t{:.12}\t{:.6e}\t{:.12}\t{:.12}\t{:.6e}",
            java_gl_retain[0],
            rust_gl_retain[0],
            (java_gl_retain[0] - rust_gl_retain[0]).abs(),
            prior_java[0],
            prior_rust[0],
            (prior_java[0] - prior_rust[0]).abs(),
            post_java[0],
            rust_post_obj.genotype_log10_posteriors[0],
            (post_java[0] - rust_post_obj.genotype_log10_posteriors[0]).abs(),
            java_gl_retain[1],
            rust_gl_retain[1],
            (java_gl_retain[1] - rust_gl_retain[1]).abs(),
            prior_java[1],
            prior_rust[1],
            (prior_java[1] - prior_rust[1]).abs(),
            post_java[1],
            rust_post_obj.genotype_log10_posteriors[1],
            (post_java[1] - rust_post_obj.genotype_log10_posteriors[1]).abs(),
            java_gl_retain[2],
            rust_gl_retain[2],
            (java_gl_retain[2] - rust_gl_retain[2]).abs(),
            prior_java[2],
            prior_rust[2],
            (prior_java[2] - prior_rust[2]).abs(),
            post_java[2],
            rust_post_obj.genotype_log10_posteriors[2],
            (post_java[2] - rust_post_obj.genotype_log10_posteriors[2]).abs(),
        ),
    );

    let collapse_changes_winner = winner_idx(&java_gl_retain) != winner_idx(&java_gl_helper);
    let collapse_margin = winner_margin(&java_gl_retain);
    let post_ranking_equal = java_gt_from_post == rust_gt_from_post;
    kv(
        "posterior_ranking_equal",
        format!(
            "{post_ranking_equal} java_winner={java_gt_from_post} rust_winner={rust_gt_from_post} collapse_changes_gl_winner={collapse_changes_winner} gl_margin={collapse_margin:.6e}"
        ),
    );

    let (classification, first_divergence) = if !predicate_equal || !production_equal {
        (
            "RETAIN_EVIDENCE_MEMBERSHIP_DIVERGENCE",
            "6R.151 41-read baseline broken".to_string(),
        )
    } else if !cells_equal || gl_retain_delta >= 1e-3 {
        (
            "EXPECTED_PRECISION_NON_CAUSAL",
            format!("raw GL residual out of band max_abs_delta={max_abs:.6e}"),
        )
    } else if relative_prior_ratios_equal {
        (
            "PRIOR_NORMALIZATION_ONLY",
            "pairwise log10 prior diffs match; vectors differ by a common offset".to_string(),
        )
    } else {
        (
            "GENOTYPE_PRIOR_RELATIVE_DIVERGENCE",
            format!(
                "assumingHW vs remainder-mass pairwise log-diff max_abs_delta={rel_max:.6e}; default USE_PLS_TO_ASSIGN does not apply gpc to GT/PL"
            ),
        )
    };

    kv("classification", classification);
    kv("first_divergence", &first_divergence);
    kv("production_change", "NONE");
    kv(
        "next_arrow",
        "default USE_PLS GT/PL already match on the 41-read object; AD/QUAL/L9/VCF remain out of scope",
    );

    assert_eq!(haps.len(), 25);
    assert_eq!(java_keep.len(), 201);
    assert_eq!(java_retained.len(), 41);
    assert_eq!(rust_prod_retained.len(), 41);
    assert!(
        predicate_equal && production_equal,
        "6R.151 41-read retainEvidence baseline must hold"
    );
    assert!(
        java_retained.contains(&mate_first) && rust_prod_retained.contains(&mate_first),
        "both overlapping mates must remain"
    );
    assert!(
        max_abs < 1e-4,
        "retained allele-likelihood residual must stay in the Java-float / Rust-f64 band"
    );
    assert!(gt_equal && pl_equal);
    assert!(post_ranking_equal);
    assert!(!relative_prior_ratios_equal);
    assert_eq!(winner_idx(&cf_a.genotype_log10_posteriors), 0);
    assert_eq!(winner_idx(&cf_b_post), 0);
    assert_eq!(classification, "GENOTYPE_PRIOR_RELATIVE_DIVERGENCE");
}
