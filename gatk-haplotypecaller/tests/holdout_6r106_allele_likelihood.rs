//! 6R.106 live dump: allele mapper + read×allele likelihoods at HOLDOUT_6R53
//! `20:29455388 C/T`. Skipped unless `HOLDOUT_6R106=1`.
//!
//! Production change: NONE. Substitution / haplotype restriction is test-only.
//!
//! ```text
//! HOLDOUT_6R106=1 cargo test -p gatk-haplotypecaller --test holdout_6r106_allele_likelihood -- --nocapture
//! ```

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::bio_ids::HaplotypeIndex;
use gatk_haplotypecaller::event_map::{
    build_per_haplotype_variation_events, overlapping_events, VariationEvent,
};
use gatk_haplotypecaller::genome_loc::GenomePosition;
use gatk_haplotypecaller::genotyping::{
    best_pl_index, diploid_genotype_alleles_from_pl_index, emit_genotype_format_fields,
    ReadLikelihoodRow,
};
use gatk_haplotypecaller::hc_allele_mapping::{
    create_allele_mapper_with_events, hap_base_at_ref_locus, AlleleHaplotypeMapping,
};
use gatk_haplotypecaller::hc_genotyping_engine::{
    biallelic_genotype_log10_likelihoods_gatk, java_alignment_read_overlaps_interval,
    marginalize_rows_to_biallelic_alleles, region_likelihoods_to_rows, InformativeAd,
    SparsePlShape, DEFAULT_INFORMATIVE_READ_OVERLAP_MARGIN, DEFAULT_STAND_EMIT_CONFIDENCE,
};
use gatk_haplotypecaller::RegionReadLikelihood;
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, traverse_assembly_region_walker,
    try_emit_call_region_variants, AssemblyRegionCallDisposition, CallRegionArgs,
    HaplotypeCallerEngine, ReadFilterParams, WalkerTraversalConfig,
};
use serde_json::{json, Value};
use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};

const INTERVAL: &str = "20:29455000-29456500";
const BAM_REL: &str = "parity/giab/runs/local-pairhmm-diff/HG001.20-29455000-29456500.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const JAVA_VCF_REL: &str = "parity/reports/6r43/chr20_tiny/java.vcf";
const RUST_VCF_REL: &str = "parity/reports/6r43/chr20_tiny/rust.vcf";
const JAVA_HASH_REL: &str = "parity/reports/6r106/java_hap_hashes.txt";
const JAVA_READ_LL_REL: &str = "parity/reports/6r106/java_read_allele_ll.tsv";
const COVERING: (u64, u64) = (29_455_300, 29_455_559);
const TARGET: u64 = 29_455_388;
const LOG10_GLOBAL_READ_MISMATCHING_RATE: f64 = -4.5;

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn fnv1a64_hex(bases: &[u8]) -> String {
    let mut h: u64 = 0xcbf29ce484222325;
    for &b in bases {
        h ^= b as u64;
        h = h.wrapping_mul(0x100000001b3);
    }
    format!("{h:016x}")
}

fn vcf_record(path: &Path, pos: u64) -> Option<Value> {
    for line in std::fs::read_to_string(path).unwrap_or_default().lines() {
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let f: Vec<_> = line.split('\t').collect();
        if f.len() < 10 || f[0] != "20" {
            continue;
        }
        if f[1].parse::<u64>().ok()? != pos {
            continue;
        }
        let fmt: Vec<_> = f[8].split(':').collect();
        let samp: Vec<_> = f[9].split(':').collect();
        let get = |k: &str| {
            fmt.iter()
                .position(|x| *x == k)
                .and_then(|i| samp.get(i))
                .unwrap_or(&".")
                .to_string()
        };
        return Some(json!({
            "alleles": format!("{}/{}", f[3], f[4]),
            "qual": f[5],
            "gt": get("GT"),
            "ad": get("AD"),
            "pl": get("PL"),
        }));
    }
    None
}

/// GATK 4.4 `createAlleleMapper` EventMap walk for a biallelic SNP `[C,T]`.
/// Empty overlapping events → REF. Spanning del (`start < loc`) with emitSpanningDels
/// → not in C or T (Java `*` pool). Unmatched overlapping events at loc → unmapped.
fn java_snp_mapper_pools(
    n_haps: usize,
    hap_events: &gatk_haplotypecaller::event_map::PerHaplotypeVariationEvents,
    loc: u64,
    ref_allele: &str,
    alt_allele: &str,
    emit_spanning: bool,
) -> (
    Vec<HaplotypeIndex>,
    Vec<HaplotypeIndex>,
    Vec<HaplotypeIndex>,
) {
    let loc_pos = GenomePosition::new_1based(loc);
    let mut ref_p = Vec::new();
    let mut alt_p = Vec::new();
    let mut other = Vec::new();
    for i in 0..n_haps {
        let spanning = overlapping_events(hap_events.events_for(i), loc);
        if spanning.is_empty() {
            ref_p.push(HaplotypeIndex::new(i));
            continue;
        }
        let mut in_alt = false;
        let mut in_ref = false;
        let mut in_star = false;
        for ev in &spanning {
            if ev.start_1based == loc_pos {
                if ev.ref_allele == ref_allele && ev.alt_allele == alt_allele {
                    in_alt = true;
                }
            } else if emit_spanning {
                in_star = true;
                break;
            } else {
                in_ref = true;
                break;
            }
        }
        if in_alt {
            alt_p.push(HaplotypeIndex::new(i));
        } else if in_ref {
            ref_p.push(HaplotypeIndex::new(i));
        } else {
            other.push(HaplotypeIndex::new(i));
            let _ = in_star;
        }
    }
    (ref_p, alt_p, other)
}

fn production_mapped_allele(i: usize, mapping: &AlleleHaplotypeMapping) -> &'static str {
    let hi = HaplotypeIndex::new(i);
    let in_ref = mapping.ref_haplotype_indices.contains(&hi);
    let in_alt = mapping.alt_haplotype_indices.contains(&hi);
    match (in_ref, in_alt) {
        (true, false) => "C",
        (false, true) => "T",
        (true, true) => "C,T",
        (false, false) => "unmapped",
    }
}

fn java_mapped_allele(
    i: usize,
    ref_p: &[HaplotypeIndex],
    alt_p: &[HaplotypeIndex],
    other: &[HaplotypeIndex],
) -> &'static str {
    let hi = HaplotypeIndex::new(i);
    let in_ref = ref_p.contains(&hi);
    let in_alt = alt_p.contains(&hi);
    let in_other = other.contains(&hi);
    match (in_ref, in_alt, in_other) {
        (true, false, _) => "C",
        (false, true, _) => "T",
        (true, true, _) => "C,T",
        (false, false, true) => "other",
        _ => "unmapped",
    }
}

/// Java `AlleleLikelihoods` mismapping floor after biallelic marginalize.
fn apply_mismapping_floor(marg: &mut [ReadLikelihoodRow]) {
    for row in marg {
        let lr = row.haplotype_log10_likelihoods[0];
        let la = row.haplotype_log10_likelihoods[1];
        if !lr.is_finite() || !la.is_finite() {
            continue;
        }
        let best = lr.max(la);
        let floor = best + LOG10_GLOBAL_READ_MISMATCHING_RATE;
        if lr < floor {
            row.haplotype_log10_likelihoods[0] = floor;
        }
        if la < floor {
            row.haplotype_log10_likelihoods[1] = floor;
        }
    }
}

fn gl_pl_from_pools(
    likelihoods: &[gatk_haplotypecaller::RegionReadLikelihood],
    n_haps: usize,
    ref_p: &[HaplotypeIndex],
    alt_p: &[HaplotypeIndex],
) -> (Vec<f64>, Vec<i32>, Vec<i32>) {
    let rows = region_likelihoods_to_rows(likelihoods, n_haps);
    let mut marg = marginalize_rows_to_biallelic_alleles(&rows, ref_p, alt_p);
    apply_mismapping_floor(&mut marg);
    let gls = biallelic_genotype_log10_likelihoods_gatk(&marg, 0, 1);
    let fmt = emit_genotype_format_fields(&gls, &[0, 0]).expect("fmt");
    let gt = diploid_genotype_alleles_from_pl_index(2, best_pl_index(&fmt.pl));
    (gls, fmt.pl_as_i32(), gt)
}

fn restrict_pools(pool: &[HaplotypeIndex], keep: &BTreeSet<usize>) -> Vec<HaplotypeIndex> {
    pool.iter()
        .copied()
        .filter(|hi| keep.contains(&hi.get()))
        .collect()
}

fn load_java_mapper(path: &Path) -> BTreeMap<String, String> {
    let mut out = BTreeMap::new();
    if !path.is_file() {
        return out;
    }
    for line in std::fs::read_to_string(path).unwrap_or_default().lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let mut parts = line.split('\t');
        if let (Some(hash), Some(mapped)) = (parts.next(), parts.next()) {
            out.insert(hash.to_string(), mapped.to_string());
        }
    }
    out
}

fn load_java_hashes(path: &Path) -> BTreeSet<String> {
    if !path.is_file() {
        return BTreeSet::new();
    }
    std::fs::read_to_string(path)
        .unwrap_or_default()
        .lines()
        .map(str::trim)
        .filter(|l| !l.is_empty() && !l.starts_with('#'))
        .map(|s| s.to_string())
        .collect()
}

fn load_java_read_ll(path: &Path) -> BTreeMap<(String, u16), (f64, f64)> {
    let mut out = BTreeMap::new();
    if !path.is_file() {
        return out;
    }
    for (i, line) in std::fs::read_to_string(path)
        .unwrap_or_default()
        .lines()
        .enumerate()
    {
        if i == 0 {
            continue;
        }
        let f: Vec<_> = line.split('\t').collect();
        if f.len() < 5 {
            continue;
        }
        let qname = f[0].to_string();
        let flags: u16 = f[1].parse().unwrap_or(0);
        let c: f64 = f[3].parse().unwrap_or(f64::NAN);
        let t: f64 = f[4].parse().unwrap_or(f64::NAN);
        out.insert((qname, flags), (c, t));
    }
    out
}

fn filter_overlapping_likelihoods(
    likelihoods: &[RegionReadLikelihood],
    reads: &[gatk_haplotypecaller::shared_bam::SharedBamRecord],
) -> Vec<RegionReadLikelihood> {
    likelihoods
        .iter()
        .filter(|rl| {
            reads.get(rl.read_index.get()).is_some_and(|r| {
                java_alignment_read_overlaps_interval(
                    r.as_ref(),
                    TARGET,
                    TARGET,
                    DEFAULT_INFORMATIVE_READ_OVERLAP_MARGIN,
                )
            })
        })
        .cloned()
        .collect()
}

fn hap_event_summary(
    hap_events: &gatk_haplotypecaller::event_map::PerHaplotypeVariationEvents,
    i: usize,
) -> String {
    overlapping_events(hap_events.events_for(i), TARGET)
        .iter()
        .map(|e| format!("{}:{}/{}", e.start_1based.get(), e.ref_allele, e.alt_allele))
        .collect::<Vec<_>>()
        .join(";")
}

#[test]
fn holdout_6r106_allele_likelihood_dump() {
    if std::env::var("HOLDOUT_6R106").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R106=1");
        return;
    }
    let root = repo_root();
    let ref_fasta = root.join(REF_REL);
    let bam = root.join(BAM_REL);
    let java_vcf = root.join(JAVA_VCF_REL);
    let rust_vcf = root.join(RUST_VCF_REL);
    assert!(ref_fasta.is_file(), "missing {}", ref_fasta.display());
    assert!(bam.is_file(), "missing {}", bam.display());

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
            ) && r.start.get() == COVERING.0
                && r.end.get() == COVERING.1
        })
        .expect("covering ActiveFull");
    let args = CallRegionArgs::strict_java();
    let outcome = HaplotypeCallerEngine::call_region(covering, &dict, &ref_fasta, &args)
        .expect("call")
        .expect("ActiveFull outcome");

    let pad = outcome.assembly.padded_reference_start_1based();
    let ref_bytes = outcome.assembly.reference_bases();
    let haps = &outcome.assembly.haplotypes;
    let max_mnp = outcome.assembly.max_mnp_distance();
    let hap_events = build_per_haplotype_variation_events(haps, ref_bytes, pad, max_mnp, "20");

    let event = VariationEvent::from_alleles("20", TARGET, "C", "T");
    let calls: Vec<_> = outcome
        .genotyped_calls
        .iter()
        .filter(|c| {
            c.event.start_1based.get() == TARGET
                && c.event.ref_allele == "C"
                && c.event.alt_allele == "T"
        })
        .collect();

    let mapping = create_allele_mapper_with_events(
        &event,
        TARGET,
        haps,
        pad,
        ref_bytes,
        max_mnp,
        true,
        Some(&hap_events),
    );
    let (java_ref, java_alt, java_other) =
        java_snp_mapper_pools(haps.len(), &hap_events, TARGET, "C", "T", true);

    let mut hashes = Vec::new();
    let mut prod_by_hash = BTreeMap::new();
    let mut java_by_hash = BTreeMap::new();
    let mut n_t_base = 0usize;
    let mut n_c_base = 0usize;
    let mut n_eventmap_ct = 0usize;
    let mut per_hap = Vec::new();
    for (i, h) in haps.iter().enumerate() {
        let hash = fnv1a64_hex(&h.bases);
        let base = hap_base_at_ref_locus(h, pad, TARGET)
            .map(|b| (b as char).to_string())
            .unwrap_or_else(|| ".".to_string());
        if base == "T" {
            n_t_base += 1;
        }
        if base == "C" {
            n_c_base += 1;
        }
        let spanning = hap_event_summary(&hap_events, i);
        let has_ct = overlapping_events(hap_events.events_for(i), TARGET)
            .iter()
            .any(|e| e.start_1based.get() == TARGET && e.ref_allele == "C" && e.alt_allele == "T");
        if has_ct {
            n_eventmap_ct += 1;
        }
        let prod = production_mapped_allele(i, &mapping);
        let java = java_mapped_allele(i, &java_ref, &java_alt, &java_other);
        hashes.push(hash.clone());
        prod_by_hash.insert(hash.clone(), prod.to_string());
        java_by_hash.insert(hash.clone(), java.to_string());
        per_hap.push(json!({
            "idx": i,
            "hash": hash,
            "len": h.bases.len(),
            "is_reference": h.is_reference,
            "base_at_target": base,
            "eventmap_ct": has_ct,
            "spanning": spanning,
            "production_mapped": prod,
            "java_walk_mapped": java,
            "mapper_mismatch": prod != java,
        }));
    }

    let n_mapper_mismatch = per_hap
        .iter()
        .filter(|h| h["mapper_mismatch"].as_bool() == Some(true))
        .count();

    let ll_all = &outcome.read_likelihoods;
    let ll_overlap = filter_overlapping_likelihoods(ll_all, &outcome.genotyping_reads);
    let n_reads_all = region_likelihoods_to_rows(ll_all, haps.len()).len();
    let ll = &ll_overlap;
    let n_reads = region_likelihoods_to_rows(ll, haps.len()).len();
    let (prod_gl, prod_pl, prod_gt) = gl_pl_from_pools(
        ll,
        haps.len(),
        &mapping.ref_haplotype_indices,
        &mapping.alt_haplotype_indices,
    );
    let (java_gl, java_pl, java_gt) = gl_pl_from_pools(ll, haps.len(), &java_ref, &java_alt);

    let live_pl = calls.first().map(|c| c.genotype.format.pl_as_i32());
    let live_gt = calls
        .first()
        .map(|c| diploid_genotype_alleles_from_pl_index(2, best_pl_index(&c.genotype.format.pl)));
    let live_gl = calls
        .first()
        .map(|c| c.genotype.genotype_log10_likelihoods.clone());

    let java_hashes = load_java_hashes(&root.join(JAVA_HASH_REL));
    let rust_hash_set: BTreeSet<String> = hashes.iter().cloned().collect();
    let common: BTreeSet<_> = rust_hash_set.intersection(&java_hashes).cloned().collect();
    let rust_only: BTreeSet<_> = rust_hash_set.difference(&java_hashes).cloned().collect();
    let java_only: BTreeSet<_> = java_hashes.difference(&rust_hash_set).cloned().collect();

    let java_mapper_path = root.join("parity/reports/6r106/java_hap_mapper.tsv");
    let java_mapped_by_hash = load_java_mapper(&java_mapper_path);
    let mut common_map_mismatch = 0usize;
    let mut common_map_rows = Vec::new();
    for hash in &common {
        let rust_prod = prod_by_hash.get(hash).map(String::as_str).unwrap_or("?");
        let rust_walk = java_by_hash.get(hash).map(String::as_str).unwrap_or("?");
        let java_live = java_mapped_by_hash
            .get(hash)
            .map(String::as_str)
            .unwrap_or("?");
        let rust_norm = if rust_prod == "C" { "C*" } else { rust_prod };
        let mismatch = java_live != "?" && rust_norm != java_live && rust_prod != java_live;
        if mismatch {
            common_map_mismatch += 1;
        }
        common_map_rows.push(json!({
            "hash": hash,
            "java_mapped": java_live,
            "rust_production": rust_prod,
            "rust_java_walk": rust_walk,
            "mismatch": mismatch,
        }));
    }

    let keep_common: BTreeSet<usize> = haps
        .iter()
        .enumerate()
        .filter(|(_, h)| java_hashes.contains(&fnv1a64_hex(&h.bases)))
        .map(|(i, _)| i)
        .collect();
    let keep_java_walk_equiv: BTreeSet<usize> = (0..haps.len())
        .filter(|&i| {
            production_mapped_allele(i, &mapping)
                == java_mapped_allele(i, &java_ref, &java_alt, &java_other)
        })
        .collect();

    let (common_gl, common_pl, common_gt) = if java_hashes.is_empty() {
        (vec![], vec![], vec![])
    } else {
        gl_pl_from_pools(
            ll,
            haps.len(),
            &restrict_pools(&mapping.ref_haplotype_indices, &keep_common),
            &restrict_pools(&mapping.alt_haplotype_indices, &keep_common),
        )
    };
    let (equiv_gl, equiv_pl, equiv_gt) = gl_pl_from_pools(
        ll,
        haps.len(),
        &restrict_pools(&mapping.ref_haplotype_indices, &keep_java_walk_equiv),
        &restrict_pools(&mapping.alt_haplotype_indices, &keep_java_walk_equiv),
    );

    let rows = region_likelihoods_to_rows(ll, haps.len());
    let mut prod_marg = marginalize_rows_to_biallelic_alleles(
        &rows,
        &mapping.ref_haplotype_indices,
        &mapping.alt_haplotype_indices,
    );
    apply_mismapping_floor(&mut prod_marg);
    let mut java_marg = marginalize_rows_to_biallelic_alleles(&rows, &java_ref, &java_alt);
    apply_mismapping_floor(&mut java_marg);
    let mut n_equal_c = 0usize;
    let mut n_equal_t = 0usize;
    let mut max_abs_c = 0.0f64;
    let mut max_abs_t = 0.0f64;
    let mut sum_abs_c = 0.0f64;
    let mut sum_abs_t = 0.0f64;
    let mut first_diff: Option<Value> = None;
    let n_cells = prod_marg.len().min(java_marg.len());
    for i in 0..n_cells {
        let pc = prod_marg[i].haplotype_log10_likelihoods[0];
        let pt = prod_marg[i].haplotype_log10_likelihoods[1];
        let jc = java_marg[i].haplotype_log10_likelihoods[0];
        let jt = java_marg[i].haplotype_log10_likelihoods[1];
        let dc = (pc - jc).abs();
        let dt = (pt - jt).abs();
        if dc == 0.0 {
            n_equal_c += 1;
        }
        if dt == 0.0 {
            n_equal_t += 1;
        }
        max_abs_c = max_abs_c.max(dc);
        max_abs_t = max_abs_t.max(dt);
        sum_abs_c += dc;
        sum_abs_t += dt;
        if first_diff.is_none() && (dc > 0.0 || dt > 0.0) {
            let rec = outcome
                .genotyping_reads
                .get(prod_marg[i].read_index)
                .map(|r| {
                    json!({
                        "qname": String::from_utf8_lossy(r.qname()),
                        "flags": r.flags(),
                    })
                });
            first_diff = Some(json!({
                "row": i,
                "read": rec,
                "production_C": pc,
                "production_T": pt,
                "java_walk_C": jc,
                "java_walk_T": jt,
            }));
        }
    }

    let java_live_ll = load_java_read_ll(&root.join(JAVA_READ_LL_REL));
    let mut rust_keys = BTreeSet::new();
    let mut n_eq_c_java = 0usize;
    let mut n_eq_t_java = 0usize;
    let mut max_c_java = 0.0f64;
    let mut max_t_java = 0.0f64;
    let mut sum_c_java = 0.0f64;
    let mut sum_t_java = 0.0f64;
    let mut n_common_reads = 0usize;
    let mut first_java_diff: Option<Value> = None;
    let mut n_rust_t_gt_c = 0usize;
    for row in &prod_marg {
        let Some(rec) = outcome.genotyping_reads.get(row.read_index) else {
            continue;
        };
        let qname = String::from_utf8_lossy(rec.qname()).into_owned();
        let flags = rec.flags();
        rust_keys.insert((qname.clone(), flags));
        let rc = row.haplotype_log10_likelihoods[0];
        let rt = row.haplotype_log10_likelihoods[1];
        if rt > rc {
            n_rust_t_gt_c += 1;
        }
        if let Some((jc, jt)) = java_live_ll.get(&(qname.clone(), flags)) {
            n_common_reads += 1;
            let dc = (rc - jc).abs();
            let dt = (rt - jt).abs();
            if dc == 0.0 {
                n_eq_c_java += 1;
            }
            if dt == 0.0 {
                n_eq_t_java += 1;
            }
            max_c_java = max_c_java.max(dc);
            max_t_java = max_t_java.max(dt);
            sum_c_java += dc;
            sum_t_java += dt;
            if first_java_diff.is_none() && (dc > 0.0 || dt > 0.0) {
                first_java_diff = Some(json!({
                    "qname": qname,
                    "flags": flags,
                    "rust_C": rc,
                    "rust_T": rt,
                    "java_C": jc,
                    "java_T": jt,
                    "delta_C": dc,
                    "delta_T": dt,
                }));
            }
        }
    }
    let java_keys: BTreeSet<_> = java_live_ll.keys().cloned().collect();
    let rust_only_reads = rust_keys.difference(&java_keys).count();
    let java_only_reads = java_keys.difference(&rust_keys).count();

    let rust_t_hashes: Vec<String> = (0..haps.len())
        .filter(|&i| production_mapped_allele(i, &mapping) == "T")
        .map(|i| fnv1a64_hex(&haps[i].bases))
        .collect();
    let java_t_hashes: Vec<String> = java_mapped_by_hash
        .iter()
        .filter(|(_, m)| *m == "T")
        .map(|(h, _)| h.clone())
        .collect();

    let emitted =
        try_emit_call_region_variants(covering, &outcome, "SAMPLE", DEFAULT_STAND_EMIT_CONFIDENCE)
            .unwrap_or_default();
    let rust_emitted = emitted.iter().any(|r| r.position == TARGET);

    let mapper_differs = n_mapper_mismatch > 0
        || mapping.ref_haplotype_indices.len() != java_ref.len()
        || mapping.alt_haplotype_indices.len() != java_alt.len();
    let java_walk_homref = java_gt.as_slice() == [0, 0];
    let production_het =
        prod_gt.as_slice() == [0, 1] || live_gt.as_ref().is_some_and(|g| g.as_slice() == [0, 1]);
    let info = InformativeAd::from_marginalized_rows(&prod_marg, 0, 1, None);
    let info_ref = info.ref_depth;
    let info_alt = info.alt_depth;
    let emitted_ad = calls.first().map(|c| c.genotype.format.ad_as_i32());

    let doc = json!({
        "holdout": "20:29455388 C/T",
        "entering_alleles": format!("{}/{}", event.ref_allele, event.alt_allele),
        "hap_count": haps.len(),
        "read_x_hap_rows_all": n_reads_all,
        "read_x_hap_rows_overlap": n_reads,
        "n_pairhmm_cells_overlap": ll.len(),
        "hap_bases_at_target": {"C": n_c_base, "T": n_t_base},
        "eventmap_ct_haps": n_eventmap_ct,
        "production_mapper": {
            "C": mapping.ref_haplotype_indices.len(),
            "T": mapping.alt_haplotype_indices.len(),
        },
        "java_walk_mapper": {
            "C": java_ref.len(),
            "T": java_alt.len(),
            "other_unmapped_star": java_other.len(),
        },
        "n_haps_mapper_assignment_differs": n_mapper_mismatch,
        "live_genotyped": {
            "pl": live_pl,
            "log10_gl": live_gl,
            "gt": live_gt,
        },
        "reconstructed_production_mapper_gl": {
            "pl": prod_pl,
            "log10_gl": prod_gl,
            "gt": prod_gt,
        },
        "substitution_java_walk_mapper_gl": {
            "pl": java_pl,
            "log10_gl": java_gl,
            "gt": java_gt,
        },
        "haplotype_count_experiments": {
            "java_hash_file_present": !java_hashes.is_empty(),
            "java_hap_count": java_hashes.len(),
            "common_hashes": common.len(),
            "rust_only_hashes": rust_only.len(),
            "java_only_hashes": java_only.len(),
            "common_mapper_mismatches": common_map_mismatch,
            "common_mapper_rows": common_map_rows,
            "rust_restricted_to_java_common": {
                "n_haps_kept": keep_common.len(),
                "pl": common_pl,
                "log10_gl": common_gl,
                "gt": common_gt,
            },
            "rust_restricted_to_java_walk_equivalent_assignment": {
                "n_haps_kept": keep_java_walk_equiv.len(),
                "pl": equiv_pl,
                "log10_gl": equiv_gl,
                "gt": equiv_gt,
            },
        },
        "read_allele_matrix_vs_java_live": {
            "java_reads": java_live_ll.len(),
            "rust_overlap_reads": rust_keys.len(),
            "common_reads": n_common_reads,
            "java_only_reads": java_only_reads,
            "rust_only_reads": rust_only_reads,
            "exact_equal_C": n_eq_c_java,
            "exact_equal_T": n_eq_t_java,
            "max_abs_delta_C": max_c_java,
            "max_abs_delta_T": max_t_java,
            "mean_abs_delta_C": if n_common_reads == 0 { 0.0 } else { sum_c_java / n_common_reads as f64 },
            "mean_abs_delta_T": if n_common_reads == 0 { 0.0 } else { sum_t_java / n_common_reads as f64 },
            "rust_T_gt_C": n_rust_t_gt_c,
            "first_differing_cell": first_java_diff,
        },
        "t_hap_hashes": {
            "rust_production": rust_t_hashes,
            "java_live": java_t_hashes,
        },
        "read_allele_matrix_production_vs_java_walk": {
            "n_reads": n_cells,
            "exact_equal_C": n_equal_c,
            "exact_equal_T": n_equal_t,
            "max_abs_delta_C": max_abs_c,
            "max_abs_delta_T": max_abs_t,
            "mean_abs_delta_C": if n_cells == 0 { 0.0 } else { sum_abs_c / n_cells as f64 },
            "mean_abs_delta_T": if n_cells == 0 { 0.0 } else { sum_abs_t / n_cells as f64 },
            "first_differing_cell": first_diff,
        },
        "final_vcf": {
            "java": vcf_record(&java_vcf, TARGET),
            "rust": vcf_record(&rust_vcf, TARGET),
            "rust_emitted_live": rust_emitted,
        },
        "causality": {
            "mapper_differs": mapper_differs,
            "java_walk_on_rust_haps_is_homref": java_walk_homref,
            "production_is_het": production_het,
            "java_walk_moves_toward_java_gls": java_walk_homref && production_het,
            "pairhmm_overlap_pl_equals_java": prod_pl.as_slice() == [0, 6, 1780],
            "live_pl_equals_sparse_het": live_pl.as_ref().is_some_and(|p| p.as_slice() == SparsePlShape::Het.pl()),
            "pairhmm_informative_ad": [info_ref, info_alt],
            "live_emitted_ad": emitted_ad,
            "sparse_het_pl": SparsePlShape::Het.pl(),
        },
        "per_haplotype": per_hap,
    });
    let out_path = root.join("parity/reports/6r106/rust_holdout.json");
    std::fs::create_dir_all(out_path.parent().unwrap()).ok();
    std::fs::write(&out_path, serde_json::to_string_pretty(&doc).unwrap()).ok();
    let mut summary = doc.clone();
    summary.as_object_mut().unwrap().remove("per_haplotype");
    if let Some(h) = summary.get_mut("haplotype_count_experiments") {
        if let Some(o) = h.as_object_mut() {
            o.remove("common_mapper_rows");
        }
    }
    println!("{}", serde_json::to_string_pretty(&summary).unwrap());

    assert_eq!(event.ref_allele, "C");
    assert_eq!(event.alt_allele, "T");
    assert!(
        calls.is_empty(),
        "6R.107: L9 must not overwrite calculator GLs; C/T must not be genotyped"
    );
    assert!(
        !rust_emitted,
        "6R.107: C/T must not be emitted after preserving calculator PL 0,6,1780"
    );
    assert_eq!(
        prod_pl.as_slice(),
        [0, 6, 1780],
        "Java retainEvidence overlap + max-marginalize + GL calculator must reproduce Java PL"
    );
    assert_eq!(
        java_pl.as_slice(),
        [0, 6, 1780],
        "Java EventMap-walk mapper substitution is not causal on the overlap set"
    );
    assert_eq!(n_common_reads, 49);
    assert_eq!(java_only_reads, 0);
    assert_eq!(rust_only_reads, 0);
    assert_eq!(mapping.alt_haplotype_indices.len(), 7);
    assert_eq!(java_alt.len(), 7);
    assert_eq!(n_mapper_mismatch, 3);
}
