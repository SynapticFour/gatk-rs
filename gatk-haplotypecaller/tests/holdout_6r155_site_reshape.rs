//! 6R.155: first mutation after the 6R.154 first AD write (37,4).
//! Diagnostic-only `unbounded_diagnostic`. Skipped unless `HOLDOUT_6R155=1`.
//!
//! ```text
//! HOLDOUT_6R155=1 GATK_RS_EXPERIMENTAL_KBEST_POLICY=unbounded_diagnostic \
//!   cargo test -p gatk-haplotypecaller --test holdout_6r155_site_reshape -- --nocapture --test-threads=1
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
    println!("6R155\t{key}\t{}", value.as_ref());
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

fn snp_pileup_ad<'a, I>(reads: I, dedupe_qname: bool) -> (i32, i32)
where
    I: IntoIterator<Item = &'a Record>,
{
    let loc0 = (TARGET as i64).saturating_sub(1);
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
        let cigar = CigarString(rec.cigar().iter().copied().collect());
        let Some(qi) = query_index_at_reference_position(rec.pos(), &cigar, loc0) else {
            continue;
        };
        let seq = rec.seq();
        if qi >= seq.len() {
            continue;
        }
        let qb = seq.as_bytes()[qi];
        match qb.to_ascii_uppercase() {
            b if b == alt_b => ra += 1,
            b if b == ref_b => rr += 1,
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
fn holdout_6r155_first_mutation_after_ad_write() {
    if std::env::var("HOLDOUT_6R155").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R155=1");
        return;
    }
    let _policy = EnvGuard::set("GATK_RS_EXPERIMENTAL_KBEST_POLICY", "unbounded_diagnostic");
    kv(
        "parity_mode",
        "unbounded_diagnostic\tnote=diagnostic_not_production_default",
    );
    kv("variant", format!("20:{TARGET} {MERGED_REF}/{MERGED_ALT}"));
    kv("alleles", "A,T");
    kv("ploidy", "2");
    kv("frozen_first_ad", "java=37,4 rust=37,4");

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
    kv(
        "restored_mates",
        format!(
            "flags99={} flags147={}",
            retained.contains(&mate99),
            retained.contains(&mate147)
        ),
    );

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
        emit_genotype_format_fields(&first_gls, &first_ad.as_vec()).expect("first write fields");
    let first_pl = first_fmt.pl_as_i32();
    let first_gt = best_pl_index(&first_fmt.pl);
    kv(
        "stage_first_ad_write",
        fmt_state(
            first_gt,
            &first_pl,
            &first_ad.as_vec(),
            first_fmt.gq.as_i32(),
            first_fmt.dp.as_i32(),
        ),
    );
    kv(
        "java_lifecycle",
        "annotateWithLikelihoods → reverseTrim(no-op A/T) → phaseCalls(off) → emit threshold (QUAL/VCF out of round)",
    );
    kv(
        "rust_lifecycle",
        "SiteScore::from_allele_mapping (InformativeAd) → SiteReshape::apply_class_a_family → GenotypeFinalize::finalize_site → genotyped_calls",
    );

    let likelihood_pileup_per_read = snp_pileup_ad(
        retained.iter().filter_map(|k| rec_by_key.get(k)).copied(),
        false,
    );
    let likelihood_pileup_dedupe = snp_pileup_ad(
        retained.iter().filter_map(|k| rec_by_key.get(k)).copied(),
        true,
    );
    let geno_pileup_dedupe =
        snp_pileup_ad(outcome.genotyping_reads.iter().map(|r| r.as_ref()), true);
    let region_pileup_dedupe = snp_pileup_ad(covering.reads.iter().map(|r| r.as_ref()), true);
    kv(
        "pileup_evidence",
        format!(
            "retainEvidence_per_read={:?} retainEvidence_qname={:?} genotyping_reads_qname={:?} walker_region_reads_qname={:?} n_region={} n_geno={} n_retain={}",
            likelihood_pileup_per_read,
            likelihood_pileup_dedupe,
            geno_pileup_dedupe,
            region_pileup_dedupe,
            covering.reads.len(),
            outcome.genotyping_reads.len(),
            retained.len()
        ),
    );

    let (pr, pa) = if region_pileup_dedupe == (22, 24) {
        region_pileup_dedupe
    } else if geno_pileup_dedupe == (22, 24) {
        geno_pileup_dedupe
    } else {
        region_pileup_dedupe
    };
    let info_ref = first_ad.ref_depth;
    let info_alt = first_ad.alt_depth;
    let early_exit = pr < 1 || pa < 1 || pr.saturating_mul(2) < pa || first_fmt.ad.len() < 2;
    let class_a = info_ref == 0 && info_alt >= 2;
    let class_a2 = first_gt == 2
        || (info_alt >= 2 && (info_ref == 0 || info_alt >= info_ref.saturating_mul(3)));
    let class_a2 = !early_exit && pa.saturating_mul(2) >= pr && class_a2;
    let class_a3 = !early_exit
        && pa.saturating_mul(2) >= pr
        && info_ref >= 1
        && info_alt >= 1
        && info_ref >= info_alt.saturating_mul(3);
    kv(
        "site_reshape_predicates",
        format!(
            "early_exit={early_exit} class_a={class_a} class_a2={class_a2} class_a3={class_a3} pileup={pr},{pa} info={info_ref},{info_alt} pl_gt={first_gt}"
        ),
    );
    kv(
        "site_reshape_invoked",
        format!("{}", class_a || class_a2 || class_a3),
    );
    kv("class_a3_invoked", format!("{class_a3}"));
    let sparse_path = class_a3 && first_gt != 1;
    kv("sparse_path_invoked", format!("{sparse_path}"));

    let shape = SparsePlShape::from_pileup_depths(pr, pa);
    let predicted = if sparse_path {
        emit_genotype_format_fields(&shape.gl_vec(), &[pr, pa]).expect("predicted sparse")
    } else if (class_a || class_a2 || class_a3) && first_gt == 1 {
        emit_genotype_format_fields(&first_gls, &[pr, pa]).expect("predicted ad-only")
    } else {
        first_fmt.clone()
    };
    kv(
        "stage_site_reshape_entry",
        fmt_state(
            first_gt,
            &first_pl,
            &first_ad.as_vec(),
            first_fmt.gq.as_i32(),
            first_fmt.dp.as_i32(),
        ),
    );
    kv(
        "stage_site_reshape_output_predicted",
        fmt_state(
            best_pl_index(&predicted.pl),
            &predicted.pl_as_i32(),
            &predicted.ad_as_i32(),
            predicted.gq.as_i32(),
            predicted.dp.as_i32(),
        ),
    );

    let prod = outcome.genotyped_calls.iter().find(|c| {
        c.event.start_1based.get() == TARGET
            && c.event.ref_allele == MERGED_REF
            && c.event.alt_allele == MERGED_ALT
    });
    match prod {
        Some(c) => {
            kv(
                "stage_final_genotyped_calls",
                fmt_state(
                    best_pl_index(&c.genotype.format.pl),
                    &c.genotype.format.pl_as_i32(),
                    &c.genotype.format.ad_as_i32(),
                    c.genotype.format.gq.as_i32(),
                    c.genotype.format.dp.as_i32(),
                ),
            );
        }
        None => kv("stage_final_genotyped_calls", "absent"),
    }

    let skip_outcome = {
        let _skip = EnvGuard::set("GATK_RS_DIAGNOSTIC_SKIP_SITE_RESHAPE", "1");
        HaplotypeCallerEngine::call_region(&java_bounds, &dict, &ref_fasta, &args)
            .expect("skip call")
            .expect("skip outcome")
    };
    let skip_call = skip_outcome.genotyped_calls.iter().find(|c| {
        c.event.start_1based.get() == TARGET
            && c.event.ref_allele == MERGED_REF
            && c.event.alt_allele == MERGED_ALT
    });
    match skip_call {
        Some(c) => kv(
            "counterfactual_skip_site_reshape",
            fmt_state(
                best_pl_index(&c.genotype.format.pl),
                &c.genotype.format.pl_as_i32(),
                &c.genotype.format.ad_as_i32(),
                c.genotype.format.gq.as_i32(),
                c.genotype.format.dp.as_i32(),
            ),
        ),
        None => kv(
            "counterfactual_skip_site_reshape",
            "absent_from_genotyped_calls (emit gate after restored hom-ref; QUAL/VCF out of round). reshape-step fields restored to first AD write 0/0 0,5,1174 37,4",
        ),
    }

    let prod_ad = prod.map(|c| c.genotype.format.ad_as_i32());
    let prod_pl = prod.map(|c| c.genotype.format.pl_as_i32());
    let predicted_matches_prod = prod_ad.as_ref() == Some(&predicted.ad_as_i32())
        && prod_pl.as_ref() == Some(&predicted.pl_as_i32());
    kv(
        "predicted_reshape_matches_genotyped_calls",
        format!("{predicted_matches_prod}"),
    );
    kv("ad_before", format!("{},{}", info_ref, info_alt));
    kv("ad_after", format!("{pr},{pa}"));
    kv("gt_before", format!("{first_gt} (0/0)"));
    kv(
        "gt_after",
        format!(
            "{} (from SparsePlShape::{shape:?})",
            best_pl_index(&predicted.pl)
        ),
    );
    kv(
        "pl_before",
        first_pl
            .iter()
            .map(|x| x.to_string())
            .collect::<Vec<_>>()
            .join(","),
    );
    kv(
        "pl_after",
        predicted
            .pl_as_i32()
            .iter()
            .map(|x| x.to_string())
            .collect::<Vec<_>>()
            .join(","),
    );

    let first_div = if first_ad.as_vec() != vec![37, 4] || first_pl != vec![0, 5, 1174] {
        "first AD write oracle broken"
    } else if !class_a3 {
        "SiteReshape Class-A3 did not fire; later operation not this family"
    } else {
        "SiteReshape::apply_class_a_family Class-A3 → sparse_snp_genotype_from_read_depths(pileup)"
    };
    kv("first_divergent_operation", first_div);
    kv(
        "java_state_before",
        "GT=0/0 PL=0,5,1174 AD=37,4 GQ=5 (annotation state)",
    );
    kv(
        "rust_state_before",
        fmt_state(
            first_gt,
            &first_pl,
            &first_ad.as_vec(),
            first_fmt.gq.as_i32(),
            first_fmt.dp.as_i32(),
        ),
    );
    kv(
        "java_state_after",
        "unchanged copy through reverseTrim/phase (no pileup reshape)",
    );
    kv(
        "rust_state_after",
        fmt_state(
            best_pl_index(&predicted.pl),
            &predicted.pl_as_i32(),
            &predicted.ad_as_i32(),
            predicted.gq.as_i32(),
            predicted.dp.as_i32(),
        ),
    );
    kv("java_analogue", "none");
    let classification = if class_a3 && !predicted_matches_prod {
        "A3_CALCULATOR_GENOTYPE_PRESERVATION"
    } else if class_a3 && sparse_path && predicted_matches_prod {
        "RUST_ONLY_MULTI_FIELD_MUTATION"
    } else {
        "NO_DIVERGENCE"
    };
    kv("classification", classification);
    kv(
        "semantic_consequence",
        "6R.158: Class-A3 no longer replaces calculator GT/PL/AD; Java keeps 0/0 0,5,1174 37,4",
    );
    kv("production_change", "YES_via_6R158");
    kv(
        "next_arrow",
        "If genotyped_calls is absent, the next first divergence is the hom-ref emit gate vs Java emit 0/0.",
    );

    assert_eq!(retained.len(), 41);
    assert!(retained.contains(&mate99) && retained.contains(&mate147));
    assert_eq!(first_ad.as_vec(), vec![37, 4]);
    assert_eq!(first_pl, vec![0, 5, 1174]);
    assert_eq!(first_gt, 0);
    assert_eq!(first_fmt.gq.as_i32(), 5);
    assert_eq!(first_fmt.dp.as_i32(), 41);
    assert_eq!(region_pileup_dedupe, (22, 24));
    assert!(class_a3);
    assert!(sparse_path);
    assert!(!SparsePlShape::pileup_is_hom_alt_strong(pr, pa));
    assert_eq!(predicted.ad_as_i32(), vec![22, 24]);
    assert_eq!(predicted.pl_as_i32(), vec![81, 0, 36]);
    assert_eq!(best_pl_index(&predicted.pl), 1);
    assert_eq!(predicted.gq.as_i32(), 36);
    assert_eq!(predicted.dp.as_i32(), 46);
    assert_eq!(classification, "A3_CALCULATOR_GENOTYPE_PRESERVATION");
    assert!(!predicted_matches_prod);
    assert!(skip_call.is_none());
    if let Some(ad) = prod_ad {
        assert_eq!(ad, first_ad.as_vec());
        assert_ne!(ad, predicted.ad_as_i32());
    }
    if let Some(pl) = prod_pl {
        assert_eq!(pl, first_pl);
        assert_ne!(pl, predicted.pl_as_i32());
    }
}
