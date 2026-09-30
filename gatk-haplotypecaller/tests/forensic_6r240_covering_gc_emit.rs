//! 6R.240: first Java/Rust divergence at covering `20:29455314 G>C`.
//!
//! Frozen Java 4.4.0.0 SHA `2dbc025821bc5f686c423ff332a41e6cef892a77`.
//! Java covering VCF omits the locus. This diagnosis forces GATK3
//! `stand_emit=10` so GLs remain visible; 6R.241 wires production
//! strict-Java `standardConfidenceForCalling=30`.
//! 6R.239 sentinel arithmetic is unchanged. K=128 is unchanged.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r240_covering_gc_emit -- --nocapture --test-threads=1
//! ```

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::cigar::CigarOperator;
use gatk_haplotypecaller::emit_gates::{
    java_emit_would_pass, passes_hc_variant_emit_biallelic, passes_java_emit_not_hom_ref,
};
use gatk_haplotypecaller::event_map::VariationEvent;
use gatk_haplotypecaller::genotyping::best_pl_index;
use gatk_haplotypecaller::haplotype::Haplotype;
use gatk_haplotypecaller::hc_genotyping_engine::{
    java_emit_af_decision, DEFAULT_STAND_EMIT_CONFIDENCE,
};
use gatk_haplotypecaller::read_threading_assembler::DEFAULT_NUM_BEST_HAPLOTYPES_PER_GRAPH;
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, traverse_assembly_region_walker,
    try_emit_call_region_variants, AssemblyRegionCallDisposition, CallRegionArgs,
    HaplotypeCallerEngine, ReadFilterParams, WalkerTraversalConfig, DEFAULT_STAND_CALL_CONF,
};
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

const JAVA_PIN: &str = "4.4.0.0 / 2dbc025821bc5f686c423ff332a41e6cef892a77";
const INTERVAL: &str = "20:29455000-29456500";
const BAM_REL: &str = "parity/giab/runs/local-pairhmm-diff/HG001.20-29455000-29456500.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const JAVA_VCF_REL: &str = "parity/reports/6r43/chr20_tiny/java.vcf";
const TARGET: u64 = 29_455_314;
const TARGET_REF: &str = "G";
const TARGET_ALT: &str = "C";
const COVERING: (u64, u64) = (29_455_300, 29_455_559);
const JAVA_CARRIER_UNTRIM: &str = "c7acc50dfb9f9ecc";
const JAVA_CARRIER_TRIM: &str = "c8074b3fd2634582";
const JAVA_TSV: &str = include_str!("forensic_6r232_java_eventmap.tsv");

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    eprintln!("6R240\t{key}\t{}", value.as_ref());
}

fn fnv1a64_hex(data: &[u8]) -> String {
    let mut h: u64 = 0xcbf29ce484222325;
    for &b in data {
        h ^= b as u64;
        h = h.wrapping_mul(0x100000001b3);
    }
    format!("{h:016x}")
}

fn parse_kv_fields(line: &str) -> BTreeMap<String, String> {
    let mut out = BTreeMap::new();
    for p in line.split('\t').skip(2) {
        if let Some((k, v)) = p.split_once('=') {
            out.insert(k.to_string(), v.to_string());
        }
    }
    out
}

fn hap_base_at(h: &Haplotype, pos: u64) -> Option<u8> {
    let gl = h.genome_loc?;
    if pos < gl.start_1based() || pos > gl.end_1based() {
        return None;
    }
    let cigar = h.cigar.as_ref()?;
    let mut ref_pos = gl.start_1based();
    let mut hap_pos = 0usize;
    for el in &cigar.elements {
        match el.operator {
            CigarOperator::Match => {
                for _ in 0..el.length {
                    if ref_pos == pos {
                        return h.bases.get(hap_pos).copied();
                    }
                    ref_pos += 1;
                    hap_pos += 1;
                }
            }
            CigarOperator::Insertion | CigarOperator::SoftClip => {
                hap_pos += el.length;
            }
            CigarOperator::Deletion => {
                for _ in 0..el.length {
                    if ref_pos == pos {
                        return None;
                    }
                    ref_pos += 1;
                }
            }
            CigarOperator::HardClip => {}
        }
    }
    None
}

fn cigar_str(h: &Haplotype) -> String {
    match &h.cigar {
        Some(c) => c
            .elements
            .iter()
            .map(|el| {
                let op = match el.operator {
                    CigarOperator::Match => 'M',
                    CigarOperator::Insertion => 'I',
                    CigarOperator::Deletion => 'D',
                    CigarOperator::SoftClip => 'S',
                    CigarOperator::HardClip => 'H',
                };
                format!("{}{op}", el.length)
            })
            .collect(),
        None => ".".to_string(),
    }
}

fn vcf_has_gc(path: &Path) -> bool {
    for line in std::fs::read_to_string(path).unwrap_or_default().lines() {
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let f: Vec<_> = line.split('\t').collect();
        if f.len() < 5 || f[0] != "20" {
            continue;
        }
        if f[1].parse::<u64>().ok() == Some(TARGET) && f[3] == TARGET_REF && f[4] == TARGET_ALT {
            return true;
        }
    }
    false
}

fn event_is_gc(e: &VariationEvent) -> bool {
    e.start_1based.get() == TARGET && e.ref_allele == TARGET_REF && e.alt_allele == TARGET_ALT
}

#[test]
fn forensic_6r240_covering_gc_first_divergence() {
    kv("java_pin", JAVA_PIN);
    kv("gkl", "0.8.8");
    kv(
        "pairhmm",
        "native float; --native-pair-hmm-use-double-precision=false",
    );
    kv("production_change", "NONE");
    kv("target", "20:29455314 G>C");
    kv(
        "predecessor",
        "6R.239 INDEX_SENTINEL_SEMANTICS_DIVERGENCE CLOSED",
    );
    assert_eq!(DEFAULT_NUM_BEST_HAPLOTYPES_PER_GRAPH, 128, "do not raise K");
    assert_eq!(DEFAULT_STAND_EMIT_CONFIDENCE, 10.0);
    assert_eq!(DEFAULT_STAND_CALL_CONF, 30.0);

    let mut java_untrimmed_event = false;
    let mut java_trimmed_event = false;
    let mut java_untrimmed_hap_n = String::new();
    let mut java_region = String::new();
    let mut java_variant_span = String::new();
    let mut java_trim_span = String::new();
    for line in JAVA_TSV.lines() {
        let f = parse_kv_fields(line);
        if line.contains("\tinput\t") && f.get("stage").map(String::as_str) == Some("untrimmed") {
            java_untrimmed_hap_n = f.get("hap_n").cloned().unwrap_or_default();
            java_region = format!(
                "{} padded={}",
                f.get("region").map(String::as_str).unwrap_or("?"),
                f.get("region_padded").map(String::as_str).unwrap_or("?")
            );
        }
        if line.contains("\tevent\t")
            && f.get("start").map(String::as_str) == Some("29455314")
            && f.get("ref").map(String::as_str) == Some(TARGET_REF)
            && f.get("alt").map(String::as_str) == Some(TARGET_ALT)
        {
            match f.get("stage").map(String::as_str) {
                Some("untrimmed") => {
                    java_untrimmed_event = true;
                    assert_eq!(f.get("hash").map(String::as_str), Some(JAVA_CARRIER_UNTRIM));
                }
                Some("trimmed") => {
                    java_trimmed_event = true;
                    assert_eq!(f.get("hash").map(String::as_str), Some(JAVA_CARRIER_TRIM));
                }
                _ => {}
            }
        }
        if line.contains("\tvariant_span\t") {
            java_variant_span = line.split('\t').last().unwrap_or("").to_string();
        }
        if line.contains("\ttrim_span\t") {
            java_trim_span = line.split('\t').last().unwrap_or("").to_string();
        }
    }
    assert!(java_untrimmed_event, "Java untrimmed EventMap has G>C");
    assert!(java_trimmed_event, "Java trimmed EventMap has G>C");
    kv("java_untrimmed_hap_n", &java_untrimmed_hap_n);
    kv("java_active_region", &java_region);
    kv("java_variant_span", &java_variant_span);
    kv("java_trim_span", &java_trim_span);
    kv("java_eventmap_gc", "present untrimmed+trimmed");

    let root = repo_root();
    let java_vcf = root.join(JAVA_VCF_REL);
    if !java_vcf.is_file() {
        eprintln!("skip: missing pinned Java covering VCF");
        return;
    }
    let java_vcf_has = vcf_has_gc(&java_vcf);
    kv(
        "java_covering_vcf",
        if java_vcf_has {
            "HAS 20:29455314 G>C"
        } else {
            "OMIT 20:29455314 G>C"
        },
    );
    assert!(!java_vcf_has, "pinned Java covering VCF must omit G>C");

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
            ) && r.start.get() == COVERING.0
                && r.end.get() == COVERING.1
        })
        .expect("covering ActiveFull");
    kv(
        "rust_active_region",
        format!(
            "{}:{}-{} extended={}-{} reads={}",
            region.contig,
            region.start.get(),
            region.end.get(),
            region.extended_start.get(),
            region.extended_end.get(),
            region.reads.len()
        ),
    );
    assert!(
        region.start.get() <= TARGET && region.end.get() >= TARGET,
        "locus inside Rust active region"
    );

    let mut args10 = CallRegionArgs::strict_java();
    args10.genotyping.stand_emit_confidence = DEFAULT_STAND_EMIT_CONFIDENCE;
    let outcome = HaplotypeCallerEngine::call_region(region, &dict, &ref_fasta, &args10)
        .expect("call")
        .expect("Some");

    let pad = outcome.assembly.padded_reference_start_1based();
    let (full_ref, full_pad) = outcome.assembly.event_map_reference();
    let ref_off = TARGET.saturating_sub(full_pad) as usize;
    let ref_base = full_ref.get(ref_off).copied().unwrap_or(b'?');
    kv(
        "assembly_interval",
        format!(
            "pad_start={pad} event_map_pad={full_pad} ref_len={} ref[{TARGET}]={}",
            full_ref.len(),
            ref_base as char
        ),
    );
    assert_eq!(ref_base, b'G');

    let haps = &outcome.assembly.haplotypes;
    let mut rust_carriers = 0usize;
    let mut rust_has_untrim_hash = false;
    let mut rust_has_trim_hash = false;
    for (i, h) in haps.iter().enumerate() {
        let hash = fnv1a64_hex(&h.bases);
        if hash == JAVA_CARRIER_UNTRIM {
            rust_has_untrim_hash = true;
        }
        if hash == JAVA_CARRIER_TRIM {
            rust_has_trim_hash = true;
        }
        if hap_base_at(h, TARGET) == Some(b'C') {
            rust_carriers += 1;
            let gl = h
                .genome_loc
                .map(|g| format!("{}-{}", g.start_1based(), g.end_1based()));
            kv(
                "rust_carrier",
                format!(
                    "idx={i} hash={hash} ref={} len={} loc={} cigar={} score={:.8} kmer={} C_at_target=1",
                    h.is_reference,
                    h.bases.len(),
                    gl.as_deref().unwrap_or("."),
                    cigar_str(h),
                    h.score,
                    h.kmer_size
                ),
            );
        }
    }
    kv(
        "haplotype_population",
        format!(
            "n={} C_carriers={rust_carriers} untrim_hash={} trim_hash={}",
            haps.len(),
            rust_has_untrim_hash,
            rust_has_trim_hash
        ),
    );
    assert!(
        rust_carriers > 0,
        "Rust haplotypes must carry G>C after 6R.239"
    );
    let carrier = haps
        .iter()
        .find(|h| hap_base_at(h, TARGET) == Some(b'C'))
        .expect("C carrier");
    assert!(
        (carrier.score + 2.64786702).abs() < 1e-8,
        "C carrier must be canonical 6R.239 path score -2.64786702, got {}",
        carrier.score
    );
    kv(
        "carrier_identity",
        "trimmed call_region haplotypes; G>C present; score matches Java untrimmed carrier; FNV may differ after trim",
    );

    let events: Vec<_> = outcome
        .assembly
        .variation_events()
        .iter()
        .filter(|e| e.start_1based.get() == TARGET)
        .cloned()
        .collect();
    kv(
        "rust_eventmap_at_target",
        if events.is_empty() {
            "ABSENT".to_string()
        } else {
            events
                .iter()
                .map(|e| {
                    format!(
                        "{}:{} {}>{}",
                        e.contig,
                        e.start_1based.get(),
                        e.ref_allele,
                        e.alt_allele
                    )
                })
                .collect::<Vec<_>>()
                .join(",")
        },
    );
    assert!(
        events.iter().any(event_is_gc),
        "Rust EventMap must contain 20:29455314 G>C"
    );

    let calls: Vec<_> = outcome
        .genotyped_calls
        .iter()
        .filter(|c| event_is_gc(&c.event))
        .collect();
    kv("genotyped_gc_n", calls.len().to_string());
    assert_eq!(calls.len(), 1, "Rust must genotype G>C as one site call");
    let call = calls[0];
    let ad = call.genotype.format.ad_as_i32();
    let pl = call.genotype.format.pl_as_i32();
    let dp = call.genotype.format.dp.as_i32();
    let gq = call.genotype.format.gq.as_i32();
    let gt_idx = best_pl_index(&call.genotype.format.pl);
    let gls = &call.genotype.genotype_log10_likelihoods;
    let qual = call
        .qual_log10_p_error
        .map(|pe| (-10.0 * pe.min(0.0)).max(0.0));
    kv(
        "rust_genotype",
        format!(
            "gt_idx={gt_idx} ad={ad:?} pl={pl:?} dp={dp} gq={gq} extra_alt={} post_merge={} qual_log10={:?} qual_phred={qual:?}",
            call.extra_alt_alleles.len(),
            call.post_merge_unused_alt_subset,
            call.qual_log10_p_error
        ),
    );

    let af10 = java_emit_af_decision(gls, DEFAULT_STAND_EMIT_CONFIDENCE).expect("af10");
    let af30 = java_emit_af_decision(gls, DEFAULT_STAND_CALL_CONF).expect("af30");
    kv(
        "af_stand10",
        format!(
            "phred={:.4} monomorphic={} alt_plausible={} passes_emit={}",
            af10.phred_scaled, af10.site_is_monomorphic, af10.alt_plausible, af10.passes_emit
        ),
    );
    kv(
        "af_stand30",
        format!(
            "phred={:.4} monomorphic={} alt_plausible={} passes_emit={}",
            af30.phred_scaled, af30.site_is_monomorphic, af30.alt_plausible, af30.passes_emit
        ),
    );
    let emit10 = java_emit_would_pass(
        &call.event,
        gls,
        &call.genotype.format,
        DEFAULT_STAND_EMIT_CONFIDENCE,
        &[],
    )
    .expect("emit10");
    let emit30 = java_emit_would_pass(
        &call.event,
        gls,
        &call.genotype.format,
        DEFAULT_STAND_CALL_CONF,
        &[],
    )
    .expect("emit30");
    kv(
        "java_emit_would_pass",
        format!(
            "stand10={emit10} stand30={emit30} not_hom_ref={}",
            passes_java_emit_not_hom_ref(gls, &call.genotype.format)
        ),
    );
    kv(
        "passes_hc_variant_emit",
        format!(
            "stand10={} stand30={}",
            passes_hc_variant_emit_biallelic(gls, DEFAULT_STAND_EMIT_CONFIDENCE).unwrap(),
            passes_hc_variant_emit_biallelic(gls, DEFAULT_STAND_CALL_CONF).unwrap()
        ),
    );

    let recs10 =
        try_emit_call_region_variants(region, &outcome, "SAMPLE", DEFAULT_STAND_EMIT_CONFIDENCE)
            .unwrap_or_default();
    let recs30 = try_emit_call_region_variants(region, &outcome, "SAMPLE", DEFAULT_STAND_CALL_CONF)
        .unwrap_or_default();
    let emit_live_10 = recs10.iter().any(|r| {
        r.position == TARGET
            && r.reference == TARGET_REF
            && r.alternate.iter().any(|a| a == TARGET_ALT)
    });
    let emit_live_30 = recs30.iter().any(|r| {
        r.position == TARGET
            && r.reference == TARGET_REF
            && r.alternate.iter().any(|a| a == TARGET_ALT)
    });
    kv(
        "rust_vcf_emit",
        format!(
            "stand10={emit_live_10} stand30={emit_live_30} n10={} n30={}",
            recs10.len(),
            recs30.len()
        ),
    );
    assert!(emit_live_10, "live Rust emits G>C at stand_emit=10");

    assert!(
        !call.post_merge_unused_alt_subset,
        "G>C is not unused-ALT subset"
    );
    assert!(
        call.extra_alt_alleles.is_empty(),
        "G>C is a biallelic allele set"
    );

    let first = if rust_carriers == 0 {
        "HAPLOTYPE_CONSTRUCTION_DIVERGENCE"
    } else if !events.iter().any(event_is_gc) {
        "EVENTMAP_DIVERGENCE"
    } else if calls.is_empty() {
        "GENOTYPE_DIVERGENCE"
    } else if emit10
        && !emit30
        && emit_live_10
        && !emit_live_30
        && !af30.passes_emit
        && af10.passes_emit
    {
        "EMISSION_PREDICATE_DIVERGENCE"
    } else if emit10 && emit30 && emit_live_10 {
        "EMISSION_PREDICATE_DIVERGENCE"
    } else {
        "GENOTYPE_DIVERGENCE"
    };
    kv("first_divergent_operation", first);
    kv(
        "classification",
        if first == "EMISSION_PREDICATE_DIVERGENCE" {
            "EMISSION_PREDICATE_DIVERGENCE — Java EventMap/haplotypes already contain G>C; Java covering VCF omit is calculateGenotypes / stand-call-conf=30 vs Rust stand_emit=10"
        } else {
            first
        },
    );
    assert_eq!(
        first, "EMISSION_PREDICATE_DIVERGENCE",
        "first arrow after 6R.239 must be emit threshold, not assembly/EventMap"
    );
    assert!(af10.passes_emit);
    assert!(!af30.passes_emit);
    assert!(emit10);
    assert!(!emit30);
    assert!(!emit_live_30);
}

/// Coordinate-free: live covering G>C PL fails Java `stand-call-conf=30`
/// and passes Rust `stand_emit=10`.
#[test]
fn forensic_6r240_covering_gc_pl_fails_java_stand_call_conf_30() {
    let event = VariationEvent::from_alleles("20", TARGET, TARGET_REF, TARGET_ALT);
    let pl = [21, 0, 1461];
    let gl: Vec<f64> = pl.iter().map(|&p| (p as f64) / -10.0).collect();
    let fmt =
        gatk_haplotypecaller::genotyping::emit_genotype_format_fields(&gl, &[35, 3]).expect("fmt");
    assert!(passes_java_emit_not_hom_ref(&gl, &fmt));
    let af10 = java_emit_af_decision(&gl, DEFAULT_STAND_EMIT_CONFIDENCE).expect("af10");
    let af30 = java_emit_af_decision(&gl, DEFAULT_STAND_CALL_CONF).expect("af30");
    assert!(af10.passes_emit);
    assert!(!af10.site_is_monomorphic);
    assert!((af10.phred_scaled - 13.6327).abs() < 0.01);
    assert!(af30.site_is_monomorphic);
    assert!(!af30.passes_emit);
    assert!(java_emit_would_pass(&event, &gl, &fmt, DEFAULT_STAND_EMIT_CONFIDENCE, &[]).unwrap());
    assert!(!java_emit_would_pass(&event, &gl, &fmt, DEFAULT_STAND_CALL_CONF, &[]).unwrap());
}
