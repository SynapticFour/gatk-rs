//! 6R.149: genotype-likelihood entry after the 6R.148-equivalent 201×25 KEEP object.
//! Diagnostic-only `unbounded_diagnostic`. Production k-best policy unchanged.
//! Skipped unless `HOLDOUT_6R149=1`.
//!
//! ```text
//! HOLDOUT_6R149=1 GATK_RS_EXPERIMENTAL_KBEST_POLICY=unbounded_diagnostic \
//!   cargo test -p gatk-haplotypecaller --test holdout_6r149_genotype_entry -- --nocapture --test-threads=1
//! ```

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::event_map::{
    build_per_haplotype_variation_events, overlapping_events, VariationEvent,
};
use gatk_haplotypecaller::hc_allele_mapping::{
    create_allele_mapper_with_events, hap_base_at_ref_locus, SPAN_DEL_ALLELE,
};
use gatk_haplotypecaller::{
    biallelic_genotype_log10_likelihoods_gatk, call_disposition, flatten_assembly_regions,
    marginalize_rows_to_biallelic_alleles, region_likelihoods_to_rows,
    traverse_assembly_region_walker, AssemblyRegionCallDisposition, CallRegionArgs, GenomePosition,
    HaplotypeCallerEngine, ReadFilterParams, WalkerTraversalConfig,
};
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

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    println!("6R149\t{key}\t{}", value.as_ref());
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
        let qname = parts[2].to_string();
        let kv = parse_kv_fields(&parts[3..]);
        let flags: u16 = kv.get("flags").and_then(|s| s.parse().ok()).unwrap_or(0);
        out.insert((qname, flags));
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

/// GATK 4.4 `AssemblyBasedCallerUtils.createAlleleMapper` assignment for one haplotype.
/// Empty overlapping EventMap → REF. Matching event at `loc` whose alt is in `merged_alts`
/// → that ALT. `start < loc` with `emitSpanningDels` → `*`. Unmatched at-loc events leave
/// the hap unmapped (not dumped into REF). Dual membership is allowed.
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
            } else if ev.ref_allele.len() < merged_ref.len()
                && merged_ref.starts_with(&ev.ref_allele)
            {
                let suffix = &merged_ref[ev.ref_allele.len()..];
                let remapped = format!("{}{}", ev.alt_allele, suffix);
                if merged_alts.iter().any(|a| a == &remapped) && !labels.contains(&remapped) {
                    labels.push(remapped);
                }
            }
        } else if emit_spanning {
            let star = SPAN_DEL_ALLELE.to_string();
            if !labels.contains(&star) {
                labels.push(star);
            }
            break;
        } else {
            if !labels.contains(&merged_ref.to_string()) {
                labels.push(merged_ref.to_string());
            }
            break;
        }
    }
    labels
}

fn rust_labels(
    i: usize,
    ref_set: &BTreeSet<usize>,
    alt_set: &BTreeSet<usize>,
    merged_ref: &str,
    merged_alt: &str,
) -> Vec<String> {
    let mut v = Vec::new();
    if ref_set.contains(&i) {
        v.push(merged_ref.to_string());
    }
    if alt_set.contains(&i) {
        v.push(merged_alt.to_string());
    }
    v
}

fn canonical_label(mut labels: Vec<String>) -> String {
    if labels.is_empty() {
        return "unmapped".to_string();
    }
    labels.sort();
    labels.join(",")
}

fn pool_max(lls: &[f64]) -> f64 {
    lls.iter().copied().fold(f64::NEG_INFINITY, f64::max)
}

fn java_reason(spanning: &[VariationEvent], _loc: u64, assigned: &str) -> String {
    if spanning.is_empty() {
        return format!("empty_overlapping_EventMap → {assigned} (Java REF)");
    }
    let desc: Vec<String> = spanning
        .iter()
        .map(|e| format!("{}:{}/{}", e.start_1based.get(), e.ref_allele, e.alt_allele))
        .collect();
    format!("spanning=[{}] → {assigned}", desc.join(";"))
}

fn rust_reason(
    spanning: &[VariationEvent],
    hap_base: Option<u8>,
    assigned: &str,
    loc: u64,
) -> String {
    if spanning.is_empty() {
        return format!(
            "empty_span hap_base={} assigned={assigned} (Rust SNP-base pool then REF fallback)",
            hap_base.map(|b| b as char).unwrap_or('?')
        );
    }
    let prior = spanning.iter().any(|e| e.start_1based.get() < loc);
    format!(
        "nonempty_span prior_start={prior} hap_base={} assigned={assigned}",
        hap_base.map(|b| b as char).unwrap_or('?')
    )
}

#[test]
fn holdout_6r149_genotype_entry() {
    if std::env::var("HOLDOUT_6R149").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R149=1");
        return;
    }
    let _policy = EnvGuard::set("GATK_RS_EXPERIMENTAL_KBEST_POLICY", "unbounded_diagnostic");
    kv(
        "parity_mode",
        "unbounded_diagnostic\tnote=diagnostic_not_production_default",
    );
    kv("resource_policy", "quarantined");

    let root = repo_root();
    let jf = root.join(JAVA_FLOAT_DUMP_REL);
    let java_keep = load_java_stage_reads(&jf, "stored");
    let java_stored = load_java_matrix(&jf, "stored");
    let java_hap_hashes = load_java_hap_hashes(&jf);
    assert_eq!(java_keep.len(), 201);
    assert_eq!(java_hap_hashes.len(), 25);

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
    let (full_ref, full_pad) = outcome.assembly.event_map_reference();
    let apply_pad = outcome
        .assembly
        .haplotypes
        .iter()
        .find(|h| h.is_reference)
        .and_then(|h| h.genome_loc.as_ref().map(|g| g.start.get()))
        .unwrap_or(full_pad);
    kv(
        "pads",
        format!(
            "event_map_pad={full_pad}\tapply_pad={apply_pad}\thap_n={}\tll_cells={}\tgeno_reads={}",
            haps.len(),
            outcome.read_likelihoods.len(),
            outcome.genotyping_reads.len()
        ),
    );
    assert_eq!(haps.len(), 25);

    let hap_cache = build_per_haplotype_variation_events(
        haps,
        full_ref,
        full_pad,
        outcome.assembly.max_mnp_distance(),
        "20",
    );

    let mut unique_at_loc: BTreeSet<(u64, String, String)> = BTreeSet::new();
    let mut unique_overlap: BTreeSet<(u64, String, String)> = BTreeSet::new();
    for i in 0..haps.len() {
        for ev in overlapping_events(hap_cache.events_for(i), TARGET) {
            unique_overlap.insert((
                ev.start_1based.get(),
                ev.ref_allele.clone(),
                ev.alt_allele.clone(),
            ));
            if ev.start_1based.get() == TARGET {
                unique_at_loc.insert((
                    ev.start_1based.get(),
                    ev.ref_allele.clone(),
                    ev.alt_allele.clone(),
                ));
            }
        }
    }
    kv(
        "events_at_loc",
        unique_at_loc
            .iter()
            .map(|(s, r, a)| format!("{s}:{r}/{a}"))
            .collect::<Vec<_>>()
            .join(","),
    );
    kv(
        "overlapping_events_union",
        unique_overlap
            .iter()
            .map(|(s, r, a)| format!("{s}:{r}/{a}"))
            .collect::<Vec<_>>()
            .join(";"),
    );

    let merged_alts: Vec<String> = unique_at_loc
        .iter()
        .filter(|(_, r, a)| r == MERGED_REF && a != MERGED_REF)
        .map(|(_, _, a)| a.clone())
        .collect();
    let merged_alts = if merged_alts.is_empty() {
        vec![MERGED_ALT.to_string()]
    } else {
        merged_alts
    };
    let merged_alt = merged_alts
        .first()
        .cloned()
        .unwrap_or_else(|| MERGED_ALT.to_string());
    let merged = VariationEvent::from_alleles("20", TARGET, MERGED_REF, &merged_alt);

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
    let ref_set: BTreeSet<usize> = mapping
        .ref_haplotype_indices
        .iter()
        .map(|i| i.get())
        .collect();
    let alt_set: BTreeSet<usize> = mapping
        .alt_haplotype_indices
        .iter()
        .map(|i| i.get())
        .collect();

    kv(
        "genotype_variant",
        format!("20:{TARGET} {MERGED_REF}/{}", merged_alts.join(",")),
    );
    kv(
        "java_alleles",
        format!("{MERGED_REF},{}", merged_alts.join(",")),
    );
    kv(
        "rust_alleles",
        format!("{},{}", mapping.ref_allele, mapping.alt_allele),
    );
    kv(
        "mapper_pool_counts",
        format!(
            "rust_ref={}\trust_alt={}\tjava_merged_alts={}",
            ref_set.len(),
            alt_set.len(),
            merged_alts.len()
        ),
    );

    let rust_hashes: Vec<String> = haps.iter().map(|h| fnv1a64_hex(&h.bases)).collect();
    let java_hash_set: BTreeSet<String> = java_hap_hashes.iter().cloned().collect();
    let rust_hash_set: BTreeSet<String> = rust_hashes.iter().cloned().collect();
    kv(
        "haplotype_identity",
        format!(
            "java=25 rust={} COMMON={} JAVA_ONLY={} RUST_ONLY={}",
            rust_hashes.len(),
            java_hash_set.intersection(&rust_hash_set).count(),
            java_hash_set.difference(&rust_hash_set).count(),
            rust_hash_set.difference(&java_hash_set).count(),
        ),
    );
    assert_eq!(
        java_hash_set.symmetric_difference(&rust_hash_set).count(),
        0
    );

    let mut n_map_eq = 0usize;
    let mut first_mismatch: Option<String> = None;
    let mut n_empty_span = 0usize;
    let mut n_empty_span_java_ref_rust_alt = 0usize;
    let mut n_prior_span = 0usize;
    let mut n_java_star = 0usize;
    let mut n_java_unmapped = 0usize;
    for (i, h) in haps.iter().enumerate() {
        let spanning = overlapping_events(hap_cache.events_for(i), TARGET);
        if spanning.is_empty() {
            n_empty_span += 1;
        }
        if spanning.iter().any(|e| e.start_1based.get() < TARGET) {
            n_prior_span += 1;
        }
        let hap_base = hap_base_at_ref_locus(h, apply_pad, TARGET);
        let j_labels = java_create_allele_mapper_labels(
            &spanning,
            TARGET,
            MERGED_REF,
            &merged_alts,
            EMIT_SPANNING_DELS,
        );
        let r_labels = rust_labels(i, &ref_set, &alt_set, MERGED_REF, &merged_alt);
        let j_c = canonical_label(j_labels.clone());
        let r_c = canonical_label(r_labels.clone());
        if j_c.contains('*') {
            n_java_star += 1;
        }
        if j_c == "unmapped" {
            n_java_unmapped += 1;
        }
        if spanning.is_empty() && j_c == MERGED_REF && r_c == merged_alt {
            n_empty_span_java_ref_rust_alt += 1;
        }
        let equal = j_c == r_c;
        if equal {
            n_map_eq += 1;
        } else if first_mismatch.is_none() {
            first_mismatch = Some(format!(
                "idx={i}\thash={}\tisRef={}\thap_base={}\tjava={j_c}\trust={r_c}\tjava_reason={}\trust_reason={}",
                rust_hashes[i],
                h.is_reference,
                hap_base.map(|b| b as char).unwrap_or('?'),
                java_reason(&spanning, TARGET, &j_c),
                rust_reason(&spanning, hap_base, &r_c, TARGET),
            ));
        }
        kv(
            "hap_map",
            format!(
                "idx={i}\thash={}\tisRef={}\thap_base={}\tspanning_n={}\tjava={j_c}\trust={r_c}\tequal={equal}",
                rust_hashes[i],
                h.is_reference,
                hap_base.map(|b| b as char).unwrap_or('?'),
                spanning.len(),
            ),
        );
    }
    kv(
        "haplotype_to_allele_equal",
        format!(
            "{n_map_eq}/25\tempty_span={n_empty_span}\tempty_span_java_ref_rust_alt={n_empty_span_java_ref_rust_alt}\tprior_span={n_prior_span}\tjava_star={n_java_star}\tjava_unmapped={n_java_unmapped}"
        ),
    );
    if let Some(m) = &first_mismatch {
        kv("first_mapping_mismatch", m);
    }

    let rust_reads: BTreeSet<ReadKey> = outcome
        .genotyping_reads
        .iter()
        .map(|r| (String::from_utf8_lossy(r.qname()).into_owned(), r.flags()))
        .collect();
    kv(
        "genotype_input_reads",
        format!(
            "java=201 rust={} COMMON={} JAVA_ONLY={} RUST_ONLY={}",
            rust_reads.len(),
            java_keep.intersection(&rust_reads).count(),
            java_keep.difference(&rust_reads).count(),
            rust_reads.difference(&java_keep).count(),
        ),
    );

    let sample_name = outcome
        .genotyping_reads
        .first()
        .and_then(|r| match r.aux(b"RG") {
            Ok(rust_htslib::bam::record::Aux::String(s)) => Some(s.to_string()),
            _ => None,
        })
        .unwrap_or_else(|| "UNKNOWN_RG".to_string());
    kv(
        "sample_ploidy",
        format!(
            "sample_count=1 sample_rg={sample_name} ploidy=2 (Java default HC; GT at calculateGLsForThisEvent is NO_CALL)"
        ),
    );

    let mut classification = if n_map_eq == 25 {
        "NO_DIVERGENCE"
    } else {
        "HAPLOTYPE_ALLELE_MAPPING_DIVERGENCE"
    };
    let mut first_divergence = if n_map_eq == 25 {
        "NONE at haplotype→allele mapping".to_string()
    } else {
        first_mismatch
            .clone()
            .unwrap_or_else(|| "mapping mismatch".to_string())
    };

    // Allele-likelihood / marginalize only as measurement after mapping is stated.
    let rust_rows = region_likelihoods_to_rows(&outcome.read_likelihoods, haps.len());
    let rust_marg = marginalize_rows_to_biallelic_alleles(
        &rust_rows,
        &mapping.ref_haplotype_indices,
        &mapping.alt_haplotype_indices,
    );
    kv(
        "allele_likelihood_dimensions_rust",
        format!("{}x2 (biallelic REF/ALT max-pool)", rust_marg.len()),
    );

    let mut java_ref_idx = Vec::new();
    let mut java_alt_idx = Vec::new();
    let mut java_star_idx = Vec::new();
    for i in 0..haps.len() {
        let spanning = overlapping_events(hap_cache.events_for(i), TARGET);
        let labs = java_create_allele_mapper_labels(
            &spanning,
            TARGET,
            MERGED_REF,
            &merged_alts,
            EMIT_SPANNING_DELS,
        );
        if labs.iter().any(|l| l == MERGED_REF) {
            java_ref_idx.push(i);
        }
        if labs.iter().any(|l| l == &merged_alt) {
            java_alt_idx.push(i);
        }
        if labs.iter().any(|l| l == SPAN_DEL_ALLELE) {
            java_star_idx.push(i);
        }
    }
    kv(
        "java_mapper_pools",
        format!(
            "REF={} ALT={} STAR={}",
            java_ref_idx.len(),
            java_alt_idx.len(),
            java_star_idx.len()
        ),
    );

    let n_java_alleles =
        1 + usize::from(!java_alt_idx.is_empty()) + usize::from(!java_star_idx.is_empty());
    kv(
        "allele_likelihood_dimensions_java",
        format!(
            "201x{n_java_alleles} (LinkedHashMap order REF then non-symbolic ALTs then * if added)"
        ),
    );

    let mut exact_equal = 0usize;
    let mut winner_equal = 0usize;
    let mut n_winner_flip = 0usize;
    let mut n_flip_near_tie = 0usize;
    let mut n_flip_material = 0usize;
    let mut n_pairs = 0usize;
    let mut max_abs = 0.0_f64;
    let mut sum_abs = 0.0_f64;
    let read_by_idx: Vec<ReadKey> = outcome
        .genotyping_reads
        .iter()
        .map(|r| (String::from_utf8_lossy(r.qname()).into_owned(), r.flags()))
        .collect();
    for row in &rust_marg {
        let Some(rk) = read_by_idx.get(row.read_index) else {
            continue;
        };
        if !java_keep.contains(rk) {
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
        let j_win = if jr > ja {
            "REF"
        } else if ja > jr {
            "ALT"
        } else {
            "TIE"
        };
        let r_win = if rr > ra {
            "REF"
        } else if ra > rr {
            "ALT"
        } else {
            "TIE"
        };
        if j_win == r_win {
            winner_equal += 1;
        } else {
            n_winner_flip += 1;
            let j_margin = (jr - ja).abs();
            if j_margin < 1e-4 {
                n_flip_near_tie += 1;
            } else {
                n_flip_material += 1;
            }
        }
    }
    let mean_abs = if n_pairs == 0 {
        0.0
    } else {
        sum_abs / n_pairs as f64
    };
    kv(
        "allele_likelihood_compare",
        format!(
            "pairs={n_pairs} exact_equal={exact_equal} max_abs_delta={max_abs:.6e} mean_abs_delta={mean_abs:.6e} winner_equal={winner_equal}/{} winner_flip={n_winner_flip} flip_near_tie_lt_1e-4={n_flip_near_tie} flip_material={n_flip_material}",
            rust_marg.len()
        ),
    );
    kv(
        "likelihood_object_vs_read_list",
        format!(
            "pairhmm_keep_rows={}\tgenotyping_reads={}\tjava_keep=201\tnote=48 extra genotyping_reads are 6R.146 poorly-modeled DROP; not in AlleleLikelihoods matrix",
            rust_marg.len(),
            outcome.genotyping_reads.len()
        ),
    );

    let java_raw_gl = {
        // Reconstruct Java max-marginalize (empty → −Inf) then GATK calculator on finite rows.
        let mut j_rows = Vec::new();
        for row in &rust_rows {
            let Some(rk) = read_by_idx.get(row.read_index) else {
                continue;
            };
            if !java_keep.contains(rk) {
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
            if !jr.is_finite() && !ja.is_finite() {
                continue;
            }
            j_rows.push(gatk_haplotypecaller::ReadLikelihoodRow {
                read_index: row.read_index,
                read_id: String::new(),
                haplotype_log10_likelihoods: vec![jr, ja],
            });
        }
        if j_rows.is_empty() {
            vec![0.0, 0.0, 0.0]
        } else {
            biallelic_genotype_log10_likelihoods_gatk(&j_rows, 0, 1)
        }
    };
    let rust_raw_gl = if rust_marg.is_empty() {
        vec![0.0, 0.0, 0.0]
    } else {
        biallelic_genotype_log10_likelihoods_gatk(&rust_marg, 0, 1)
    };
    kv(
        "raw_gl_java",
        format!(
            "0/0={:.8} 0/1={:.8} 1/1={:.8}",
            java_raw_gl[0], java_raw_gl[1], java_raw_gl[2]
        ),
    );
    kv(
        "raw_gl_rust",
        format!(
            "0/0={:.8} 0/1={:.8} 1/1={:.8}",
            rust_raw_gl[0], rust_raw_gl[1], rust_raw_gl[2]
        ),
    );
    let gl_max_d = java_raw_gl
        .iter()
        .zip(rust_raw_gl.iter())
        .map(|(a, b)| (a - b).abs())
        .fold(0.0_f64, f64::max);
    kv("raw_gl_max_abs_delta", format!("{gl_max_d:.6e}"));

    let target_calls: Vec<_> = outcome
        .genotyped_calls
        .iter()
        .filter(|c| c.event.start_1based.get() == TARGET)
        .collect();
    kv(
        "production_calls_at_target",
        if target_calls.is_empty() {
            "NONE".to_string()
        } else {
            target_calls
                .iter()
                .map(|c| {
                    format!(
                        "{}/{} GT_from_PL_len={} extra_alts={}",
                        c.event.ref_allele,
                        c.event.alt_allele,
                        c.genotype.format.pl.len(),
                        c.extra_alt_alleles.join(",")
                    )
                })
                .collect::<Vec<_>>()
                .join(" ; ")
        },
    );

    // Do not reclassify to RAW_GL if mapping already diverged.
    if n_map_eq == 25 {
        if max_abs > 1e-4 {
            classification = "ALLELE_LIKELIHOOD_DIVERGENCE";
            first_divergence = format!("allele likelihood max_abs_delta={max_abs:.6e}");
        } else if gl_max_d > 1e-3 {
            classification = "RAW_GL_DIVERGENCE";
            first_divergence = format!("raw GL max_abs_delta={gl_max_d:.6e}");
        } else if max_abs > 0.0 || gl_max_d > 0.0 {
            classification = "EXPECTED_PRECISION_NON_CAUSAL";
            first_divergence = "NONE (Java-float / Rust-f64 residual)".to_string();
        }
    }

    kv("classification", classification);
    kv("first_divergence", &first_divergence);
    kv("production_change", "NONE");
    kv(
        "next_arrow",
        if classification == "HAPLOTYPE_ALLELE_MAPPING_DIVERGENCE" {
            "implement Java empty-EventMap→REF / no SNP-base pooling on this evidence class (no production change this round)"
        } else if classification == "NO_DIVERGENCE"
            || classification == "EXPECTED_PRECISION_NON_CAUSAL"
        {
            "retainEvidence overlap + genotype priors / calculateGenotypes (not AD/PL/QUAL/L9/VCF)"
        } else {
            "stop at first_divergence; no production change this round"
        },
    );

    assert_eq!(haps.len(), 25);
    assert_eq!(java_keep.len(), 201);
    assert_eq!(n_map_eq, 25, "haplotype→allele mapping must be 25/25");
    assert_eq!(n_empty_span_java_ref_rust_alt, 0);
    assert_eq!(n_java_star, 0);
    assert_eq!(
        n_flip_material, 0,
        "no material allele-likelihood winner flips"
    );
    assert!(
        max_abs < 1e-4,
        "allele-likelihood residual must stay in the Java-float / Rust-f64 band"
    );
    assert!(
        gl_max_d < 1e-3,
        "raw GL residual must stay in the established precision band"
    );
}
