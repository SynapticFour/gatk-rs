//! 6R.160: lifecycle dump for the first genuine VCF divergence `2:92316347 G/A`.
//! Production `call_region` (no diagnostic k-best). Skipped unless `HOLDOUT_6R160=1`.
//!
//! Does not patch production. Does not assert a future FORMAT contract.
//!
//! ```text
//! HOLDOUT_6R160=1 cargo test -p gatk-haplotypecaller --test holdout_6r160_next_vcf_divergence -- --nocapture --test-threads=1
//! ```

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::event_map::{build_per_haplotype_variation_events, VariationEvent};
use gatk_haplotypecaller::hc_allele_mapping::create_allele_mapper_with_events;
use gatk_haplotypecaller::hc_genotyping_engine::java_alignment_read_overlaps_interval;
use gatk_haplotypecaller::query_index_at_reference_position;
use gatk_haplotypecaller::{
    biallelic_genotype_log10_likelihoods_gatk, call_disposition, flatten_assembly_regions,
    java_emit_af_decision, marginalize_rows_to_biallelic_alleles, region_likelihoods_to_rows,
    traverse_assembly_region_walker, AssemblyRegionCallDisposition, CallRegionArgs, GenomePosition,
    HaplotypeCallerEngine, InformativeAd, ReadFilterParams, WalkerTraversalConfig,
    DEFAULT_INFORMATIVE_READ_OVERLAP_MARGIN, DEFAULT_STAND_EMIT_CONFIDENCE,
};
use rust_htslib::bam::record::CigarString;
use rust_htslib::bam::Record;
use std::collections::BTreeSet;
use std::path::{Path, PathBuf};

const INTERVAL: &str = "2:92316200-92316580";
const BAM_REL: &str = "parity/realworld/na12878_20k_b37/NA12878_20k.b37.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const TARGET: u64 = 92_316_347;
const MERGED_REF: &str = "G";
const MERGED_ALT: &str = "A";
const JAVA_STAND_CALL_CONF: f64 = 30.0;

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    println!("6R160\t{key}\t{}", value.as_ref());
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

#[test]
fn holdout_6r160_first_vcf_divergence_lifecycle() {
    if std::env::var("HOLDOUT_6R160").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R160=1");
        return;
    }
    assert_ne!(
        std::env::var("GATK_RS_EXPERIMENTAL_KBEST_POLICY")
            .ok()
            .as_deref(),
        Some("unbounded_diagnostic"),
        "6R.160 VCF corpus is production run_haplotype_caller / legacy_1024"
    );

    let root = repo_root();
    let ref_fasta = root.join(REF_REL);
    let bam = root.join(BAM_REL);
    if !ref_fasta.is_file() || !bam.is_file() {
        eprintln!("skip: missing P12 ref/BAM");
        return;
    }

    kv(
        "java_pin",
        "4.4.0.0 / 2dbc025821bc5f686c423ff332a41e6cef892a77",
    );
    kv(
        "java_cmd",
        "broadinstitute/gatk:4.4.0.0 HaplotypeCaller default --native-pair-hmm-threads 1 -L 2:92316200-92316580 (frozen 6R.43 java.vcf, not regenerated)",
    );
    kv(
        "rust_cmd",
        "HOLDOUT_6R43=1 cargo test holdout_6r43_test; run_haplotype_caller GatkConfig intervals=2:92316200-92316580 ReadFilterParams::gatk_standard_hc WalkerTraversalConfig::gatk_haplotype_caller_production(100) CallRegionArgs::strict_java k-best=legacy_1024",
    );
    kv("variant", format!("2:{TARGET} {MERGED_REF}/{MERGED_ALT}"));
    kv("java_vcf", "GT=1/1 AD=0,3 DP=3 GQ=9 PL=135,9,0 QUAL=121.84");
    kv("rust_vcf", "GT=1/1 AD=0,1 DP=1 GQ=3 PL=45,3,0 QUAL=35.48");

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
        .expect("ActiveFull covering 2:92316347");
    kv(
        "active_full",
        format!(
            "{}:{}-{}\textended={}-{}\tn_reads={}",
            covering.contig,
            covering.start.get(),
            covering.end.get(),
            covering.extended_start.get(),
            covering.extended_end.get(),
            covering.reads.len()
        ),
    );

    let args = CallRegionArgs::strict_java();
    let outcome = HaplotypeCallerEngine::call_region(covering, &dict, &ref_fasta, &args)
        .expect("call")
        .expect("outcome");

    let event = VariationEvent::from_alleles("2", TARGET, MERGED_REF, MERGED_ALT);
    let in_eventmap = outcome.assembly.variation_events().iter().any(|e| {
        e.start_1based == GenomePosition::new_1based(TARGET)
            && e.ref_allele == MERGED_REF
            && e.alt_allele == MERGED_ALT
    });
    kv(
        "eventmap",
        format!(
            "has_G_A={in_eventmap} n_events={}",
            outcome.assembly.variation_events().len()
        ),
    );
    assert!(
        in_eventmap,
        "EventMap must contain 2:92316347 G/A (allele set already 17/17)"
    );

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
        "haplotype_population",
        format!(
            "n_haps={} n_pairhmm_rows={} pad={apply_pad} geno_reads={}",
            haps.len(),
            outcome.read_likelihoods.len(),
            outcome.genotyping_reads.len()
        ),
    );

    let hap_cache = build_per_haplotype_variation_events(
        haps,
        full_ref,
        full_pad,
        outcome.assembly.max_mnp_distance(),
        "2",
    );
    let mapping = create_allele_mapper_with_events(
        &event,
        TARGET,
        haps,
        apply_pad,
        outcome.assembly.reference_bases(),
        outcome.assembly.max_mnp_distance(),
        true,
        Some(&hap_cache),
    );
    kv(
        "allele_mapper",
        format!(
            "n_ref_haps={} n_alt_haps={}",
            mapping.ref_haplotype_indices.len(),
            mapping.alt_haplotype_indices.len()
        ),
    );
    assert!(
        !mapping.alt_haplotype_indices.is_empty(),
        "alt haplotype pool must be non-empty if EventMap has G/A"
    );

    let pileup_region = snp_pileup_ad(covering.reads.iter().map(|r| r.as_ref()), false);
    let pileup_region_qname = snp_pileup_ad(covering.reads.iter().map(|r| r.as_ref()), true);
    let pileup_geno = snp_pileup_ad(outcome.genotyping_reads.iter().map(|r| r.as_ref()), false);
    let pileup_geno_qname =
        snp_pileup_ad(outcome.genotyping_reads.iter().map(|r| r.as_ref()), true);
    kv(
        "pileup_AD",
        format!(
            "region={},{}\tregion_qname={},{}\tgeno={},{}\tgeno_qname={},{}",
            pileup_region.0,
            pileup_region.1,
            pileup_region_qname.0,
            pileup_region_qname.1,
            pileup_geno.0,
            pileup_geno.1,
            pileup_geno_qname.0,
            pileup_geno_qname.1
        ),
    );

    let rust_rows = region_likelihoods_to_rows(&outcome.read_likelihoods, haps.len());
    let rust_marg = marginalize_rows_to_biallelic_alleles(
        &rust_rows,
        &mapping.ref_haplotype_indices,
        &mapping.alt_haplotype_indices,
    );
    let full_ad = InformativeAd::from_marginalized_rows(&rust_marg, 0, 1, None);
    kv(
        "informative_AD_full_pairhmm",
        format!("{},{}", full_ad.ref_depth, full_ad.alt_depth),
    );

    let mut retained_rows = Vec::new();
    for row in &rust_marg {
        let Some(rec) = outcome.genotyping_reads.get(row.read_index) else {
            continue;
        };
        if java_alignment_read_overlaps_interval(
            rec.as_ref(),
            TARGET,
            TARGET,
            DEFAULT_INFORMATIVE_READ_OVERLAP_MARGIN,
        ) {
            retained_rows.push(row.clone());
        }
    }
    let retain_ad = InformativeAd::from_marginalized_rows(&retained_rows, 0, 1, None);
    kv(
        "informative_AD_retainEvidence_overlap",
        format!(
            "n_rows={}\tAD={},{}",
            retained_rows.len(),
            retain_ad.ref_depth,
            retain_ad.alt_depth
        ),
    );
    let retain_gls = biallelic_genotype_log10_likelihoods_gatk(&retained_rows, 0, 1);
    kv("calculateGLs_retainEvidence", format!("{retain_gls:?}"));

    let call = outcome
        .genotyped_calls
        .iter()
        .find(|c| {
            c.event.start_1based == GenomePosition::new_1based(TARGET)
                && c.event.ref_allele == MERGED_REF
                && c.event.alt_allele == MERGED_ALT
        })
        .expect("site is VCF-emitted; genotyped_calls must contain it");
    let fmt = &call.genotype.format;
    kv(
        "genotyped_calls_FORMAT",
        format!(
            "PL={:?} AD={:?} GQ={} DP={}",
            fmt.pl_as_i32(),
            fmt.ad_as_i32(),
            fmt.gq.as_i32(),
            fmt.dp.as_i32()
        ),
    );
    kv(
        "calculator_GLs",
        format!("{:?}", call.genotype.genotype_log10_likelihoods),
    );

    let af10 = java_emit_af_decision(
        &call.genotype.genotype_log10_likelihoods,
        DEFAULT_STAND_EMIT_CONFIDENCE,
    )
    .expect("af10");
    let af30 = java_emit_af_decision(
        &call.genotype.genotype_log10_likelihoods,
        JAVA_STAND_CALL_CONF,
    )
    .expect("af30");
    kv(
        "af_threshold",
        format!(
            "rust_stand_emit={}\tjava_stand_call={}\taf10_mono={}\taf10_pass={}\taf30_mono={}\taf30_pass={}",
            DEFAULT_STAND_EMIT_CONFIDENCE,
            JAVA_STAND_CALL_CONF,
            af10.site_is_monomorphic,
            af10.passes_emit,
            af30.site_is_monomorphic,
            af30.passes_emit
        ),
    );

    let gl_best = retain_gls
        .iter()
        .enumerate()
        .max_by(|a, b| a.1.partial_cmp(b.1).unwrap())
        .map(|(i, _)| i)
        .unwrap_or(0);
    kv(
        "pairhmm_GL_best_index",
        format!("{gl_best}\t0=0/0 1=0/1 2=1/1"),
    );
    kv(
        "first_divergent_arrow",
        "Java DepthPerAlleleBySample AD=0,3 (PL 135,9,0 1/1) ≠ Rust InformativeAd after remarg AD=2,5 (PairHMM GLs het-best). retainEvidence overlap did not change AD. Sparse 45,3,0 / AD 0,1 overwrite is downstream — stop here.",
    );
    kv("af_10_vs_30_causal", "NO");
    kv("production_change", "NONE");
}
