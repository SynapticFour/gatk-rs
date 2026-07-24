//! Map each Java-only P12 variant to Rust `call_region` / emit outcome.
//! Run: `P12_REFERENCE=…/hs37d5.simple.fa cargo test -p gatk-haplotypecaller p12_java_site_trace --release -- --nocapture`

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::genome_loc::GenomePosition;
use gatk_haplotypecaller::{
    call_disposition, diagnose_genotype_variation_event, flatten_assembly_regions,
    hc_emit_policy::explain_strict_java_emit_gates,
    read_event_discovery::{
        p12_java_event_registry_enabled, read_allele_depths_at_locus, strict_java_asm8_only_enabled,
    },
    traverse_assembly_region_walker, try_emit_call_region_variants, AssemblyRegionCallDisposition,
    CallRegionArgs, GenotypeRejectReason, HaplotypeCallerEngine, ReadFilterParams,
    WalkerTraversalConfig,
};
use std::collections::BTreeSet;
use std::io::{BufRead, BufReader};
use std::path::Path;

const STAND_EMIT: f64 = 10.0;

fn p12_paths() -> Option<(std::path::PathBuf, std::path::PathBuf)> {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("..");
    let ref_path = std::env::var("P12_REFERENCE")
        .map(std::path::PathBuf::from)
        .unwrap_or_else(|_| root.join("parity/realworld/assets/hs37d5.simple.fa"));
    let ref_path = if ref_path.is_absolute() {
        ref_path
    } else {
        root.join(ref_path)
    };
    let bam = root.join("parity/realworld/na12878_20k_b37/NA12878_20k.b37.bam");
    if !ref_path.is_file() {
        eprintln!("skip: P12_REFERENCE not found: {}", ref_path.display());
        return None;
    }
    if !bam.is_file() {
        eprintln!("skip: P12 BAM not found: {}", bam.display());
        return None;
    }
    Some((ref_path, bam))
}

fn load_java_only() -> Vec<(u64, String, String)> {
    let path = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../parity/fixtures/p12-java-production-emit/p12_production_emit_sites.tsv");
    let file = std::fs::File::open(path).expect("java_only tsv");
    let mut out = Vec::new();
    for line in BufReader::new(file).lines().skip(1) {
        let line = line.expect("line");
        let cols: Vec<_> = line.split('\t').collect();
        if cols.len() >= 4 {
            out.push((
                cols[1].parse().expect("pos"),
                cols[2].into(),
                cols[3].into(),
            ));
        }
    }
    out
}

fn trace_limit() -> usize {
    std::env::var("P12_TRACE_LIMIT")
        .ok()
        .and_then(|v| v.parse::<usize>().ok())
        .unwrap_or(0)
}

fn apply_trace_limit(mut sites: Vec<(u64, String, String)>) -> Vec<(u64, String, String)> {
    let limit = trace_limit();
    if limit > 0 && sites.len() > limit {
        sites.truncate(limit);
    }
    sites
}

/// Outcome trace. Run: `env -u P12_PHASE_E cargo test p12_java_site_trace --release -- --ignored --nocapture`
#[test]
#[ignore = "L3 trace (~10+ min); default graph-only production or P12_PHASE_E=1"]
fn p12_java_site_trace() {
    if std::env::var("P12_PHASE_E").is_err() && !strict_java_asm8_only_enabled() {
        eprintln!(
            "skip p12_java_site_trace: graph-only off (set P12_PHASE_E=1 or unset GATK_RS_P12_ENSURE_BRIDGES)"
        );
        return;
    }
    let Some((ref_fasta, bam)) = p12_paths() else {
        eprintln!("skip: set P12_REFERENCE");
        return;
    };
    let java_only = apply_trace_limit(load_java_only());
    let interval = "2:92300000-92350000";
    let dict = SequenceDictionary::from_fasta_path(&ref_fasta).expect("dict");
    let specs = parse_intervals_cli_string(&dict, interval).expect("interval");
    let filters = ReadFilterParams::gatk_standard_hc();
    let cfg = WalkerTraversalConfig::gatk_haplotype_caller_production(100);
    let walk = traverse_assembly_region_walker(&dict, &specs, &ref_fasta, &bam, &filters, &cfg)
        .expect("walk");
    let regions = flatten_assembly_regions(&walk);
    let args = CallRegionArgs::strict_java();

    let mut emitted_keys = BTreeSet::new();
    for region in &regions {
        if !matches!(
            call_disposition(region),
            AssemblyRegionCallDisposition::ActiveFull
        ) {
            continue;
        }
        let Some(outcome) =
            HaplotypeCallerEngine::call_region(region, &dict, &ref_fasta, &args).expect("call")
        else {
            continue;
        };
        for rec in
            try_emit_call_region_variants(region, &outcome, "SAMPLE", STAND_EMIT).expect("emit")
        {
            let alt = rec.alternate.first().cloned().unwrap_or_default();
            emitted_keys.insert((rec.position, rec.reference.clone(), alt));
        }
    }

    let mut inactive = 0usize;
    let mut no_call = 0usize;
    let mut active_no_emit = 0usize;
    let mut shared_emit = 0usize;
    let mut bucket_no_event = 0usize;
    let mut bucket_genotyped = 0usize;
    let mut bucket_no_alt_hap = 0usize;
    let mut bucket_no_reads = 0usize;
    let mut bucket_not_confident = 0usize;
    let mut bucket_low_gq = 0usize;
    let mut bucket_other = 0usize;
    let gt_cfg = args.genotyping.clone();

    for (pos, ref_a, alt_a) in &java_only {
        let covering: Vec<_> = regions
            .iter()
            .filter(|r| r.start.get() <= *pos && r.end.get() >= *pos)
            .collect();
        let active = covering.iter().find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            )
        });
        if active.is_none() {
            inactive += 1;
            eprintln!("JAVA_ONLY\t{pos}\t{ref_a}/{alt_a}\tinactive_region");
            continue;
        }
        let region = active.unwrap();
        let key = (*pos, ref_a.clone(), alt_a.clone());
        if emitted_keys.contains(&key) {
            shared_emit += 1;
            eprintln!("JAVA_ONLY\t{pos}\t{ref_a}/{alt_a}\temitted");
            continue;
        }
        let outcome =
            HaplotypeCallerEngine::call_region(region, &dict, &ref_fasta, &args).expect("call");
        let Some(outcome) = outcome else {
            no_call += 1;
            eprintln!("JAVA_ONLY\t{pos}\t{ref_a}/{alt_a}\tcall_none");
            continue;
        };
        let mut event_hit = false;
        let mut matching_event = None;
        for e in outcome.assembly.variation_events() {
            if e.start_1based == GenomePosition::new_1based(*pos)
                && &e.ref_allele == ref_a
                && &e.alt_allele == alt_a
            {
                event_hit = true;
                matching_event = Some(e.clone());
            }
        }
        let mut matching_call = None;
        for c in &outcome.genotyped_calls {
            if c.event.start_1based == GenomePosition::new_1based(*pos)
                && c.event.ref_allele == *ref_a
                && c.event.alt_allele == *alt_a
            {
                matching_call = Some(c);
            }
        }
        let call_hit = matching_call.is_some();
        active_no_emit += 1;
        let reject = if !event_hit {
            bucket_no_event += 1;
            "no_event"
        } else if call_hit {
            bucket_genotyped += 1;
            "genotyped_not_emitted"
        } else if let Some(ref ev) = matching_event {
            let ref_hap = outcome
                .assembly
                .haplotypes
                .iter()
                .find(|h| h.is_reference)
                .expect("ref hap");
            let ref_bytes = ref_hap.bases.clone();
            let pad = ref_hap
                .genome_loc
                .as_ref()
                .map(|g| g.start_1based())
                .unwrap_or_else(|| outcome.assembly.padded_reference_start_1based());
            let max_mnp = outcome.assembly.max_mnp_distance();
            match diagnose_genotype_variation_event(
                ev,
                &outcome.read_likelihoods,
                &outcome.genotyping_reads,
                &outcome.genotyping_reads,
                Some(region.reads.as_slice()),
                &outcome.assembly.haplotypes,
                &ref_bytes,
                pad,
                outcome.assembly.reference_bases(),
                outcome.assembly.padded_reference_start_1based(),
                region.start.get(),
                region.end.get(),
                max_mnp,
                &gt_cfg,
            )
            .expect("diagnose")
            {
                Ok(_) => {
                    bucket_genotyped += 1;
                    "genotyped_not_emitted"
                }
                Err(GenotypeRejectReason::NoAltHapSupport) => {
                    bucket_no_alt_hap += 1;
                    "no_alt_hap_support"
                }
                Err(GenotypeRejectReason::NoReadLikelihoods) => {
                    bucket_no_reads += 1;
                    "no_read_likelihoods"
                }
                Err(GenotypeRejectReason::VariantNotConfident) => {
                    bucket_not_confident += 1;
                    "variant_not_confident"
                }
                Err(GenotypeRejectReason::LowGq) => {
                    bucket_low_gq += 1;
                    "low_gq"
                }
            }
        } else {
            bucket_other += 1;
            "other"
        };
        eprintln!(
            "JAVA_ONLY\t{pos}\t{ref_a}/{alt_a}\tactive:{}-{} event={event_hit} genotyped={call_hit} reject={reject} calls={}",
            region.start.get(),
            region.end.get(),
            outcome.genotyped_calls.len()
        );
        if reject == "genotyped_not_emitted" {
            let ref_hap = outcome
                .assembly
                .haplotypes
                .iter()
                .find(|h| h.is_reference)
                .expect("ref hap");
            let pad = ref_hap
                .genome_loc
                .as_ref()
                .map(|g| g.start_1based())
                .unwrap_or_else(|| outcome.assembly.padded_reference_start_1based());
            let (read_ref_ad, read_alt_ad) =
                read_allele_depths_at_locus(&region.reads, matching_event.as_ref().unwrap(), pad);
            if let Some(call) = matching_call {
                let gates = explain_strict_java_emit_gates(
                    &call.event,
                    &call.genotype.genotype_log10_likelihoods,
                    &call.genotype.format,
                    gt_cfg.stand_emit_confidence,
                    gt_cfg.genotype_stored_events_only,
                    read_ref_ad,
                    read_alt_ad,
                    &[],
                )
                .expect("emit gates");
                eprintln!("JAVA_ONLY\t{pos}\temit_gates\t{gates}");
            } else if let Some(ref ev) = matching_event {
                if let Ok(Ok(call)) = diagnose_genotype_variation_event(
                    ev,
                    &outcome.read_likelihoods,
                    &outcome.genotyping_reads,
                    &outcome.genotyping_reads,
                    Some(region.reads.as_slice()),
                    &outcome.assembly.haplotypes,
                    &ref_hap.bases,
                    pad,
                    outcome.assembly.reference_bases(),
                    outcome.assembly.padded_reference_start_1based(),
                    region.start.get(),
                    region.end.get(),
                    outcome.assembly.max_mnp_distance(),
                    &gt_cfg,
                ) {
                    let gates = explain_strict_java_emit_gates(
                        &call.event,
                        &call.genotype.genotype_log10_likelihoods,
                        &call.genotype.format,
                        gt_cfg.stand_emit_confidence,
                        gt_cfg.genotype_stored_events_only,
                        read_ref_ad,
                        read_alt_ad,
                        &[],
                    )
                    .expect("emit gates");
                    eprintln!(
                        "JAVA_ONLY\t{pos}\temit_gates\t{gates}\tnote=diagnose_ok_not_in_calls"
                    );
                }
            }
        }
    }

    eprintln!("=== java site trace ===");
    eprintln!("java_only\t{}", java_only.len());
    eprintln!("inactive\t{inactive}");
    eprintln!("call_none\t{no_call}");
    eprintln!("active_no_emit\t{active_no_emit}");
    eprintln!("would_match_emit\t{shared_emit}");
    eprintln!("rust_emitted_total\t{}", emitted_keys.len());
    eprintln!("bucket_no_event\t{bucket_no_event}");
    eprintln!("bucket_genotyped_not_emitted\t{bucket_genotyped}");
    eprintln!("bucket_no_alt_hap_support\t{bucket_no_alt_hap}");
    eprintln!("bucket_no_read_likelihoods\t{bucket_no_reads}");
    eprintln!("bucket_variant_not_confident\t{bucket_not_confident}");
    eprintln!("bucket_low_gq\t{bucket_low_gq}");
    eprintln!("bucket_other\t{bucket_other}");
    eprintln!("parity_mode\tstrict_java");
    eprintln!(
        "p12_event_registry\t{}",
        if p12_java_event_registry_enabled() {
            "on"
        } else {
            "off"
        }
    );
    eprintln!(
        "genotype_stored_events_only\t{}",
        gt_cfg.genotype_stored_events_only
    );
    eprintln!(
        "asm8_only\t{}",
        if gatk_haplotypecaller::read_event_discovery::strict_java_asm8_only_enabled() {
            "on"
        } else {
            "off"
        }
    );
    eprintln!(
        "p12_ensure_bridges\t{}",
        if gatk_haplotypecaller::read_event_discovery::strict_java_p12_ensure_bridges_enabled() {
            "on"
        } else {
            "off"
        }
    );
    assert_eq!(
        bucket_not_confident, 0,
        "strict_java variant_not_confident bucket {bucket_not_confident}"
    );
    let short_trace = trace_limit() > 0;
    if short_trace {
        eprintln!(
            "short12_gate\twould_match_emit\t{shared_emit}/{}\tbucket_no_event\t{bucket_no_event}\tbucket_no_alt_hap\t{bucket_no_alt_hap}",
            java_only.len()
        );
        assert_eq!(bucket_no_event, 0, "short12: no_event {bucket_no_event}");
        assert!(
            shared_emit >= 10,
            "short12: would_match_emit {shared_emit}/{} (target ≥10/12)",
            java_only.len()
        );
        return;
    }
    if p12_java_event_registry_enabled() {
        assert_eq!(
            bucket_no_event, 0,
            "registry on: no_event {bucket_no_event}"
        );
        assert_eq!(
            bucket_genotyped, 0,
            "registry on: genotyped_not_emitted {bucket_genotyped}"
        );
        assert_eq!(no_call, 0, "registry on: call_none {no_call}");
        assert_eq!(
            shared_emit, 66,
            "registry on M5: shared_emit {shared_emit}/66"
        );
        assert_eq!(
            emitted_keys.len(),
            66,
            "registry on M5: rust_emitted_total {}",
            emitted_keys.len()
        );
    } else {
        assert_eq!(
            bucket_no_event, 0,
            "graph-only ASM-8: no_event {bucket_no_event}"
        );
        assert_eq!(
            bucket_genotyped, 0,
            "graph-only: genotyped_not_emitted {bucket_genotyped}"
        );
        assert_eq!(no_call, 0, "graph-only: call_none {no_call}");
        assert_eq!(
            shared_emit,
            java_only.len(),
            "graph-only M5: would_match_emit {shared_emit}/{}",
            java_only.len()
        );
        let rust_only_extra = emitted_keys.len().saturating_sub(shared_emit);
        eprintln!(
            "full_gate\twould_match_emit\t{shared_emit}/{}\ttrack_rust_only_extra_vcf\t{rust_only_extra}",
            java_only.len()
        );
        assert_eq!(
            emitted_keys.len(),
            java_only.len(),
            "rust_emitted_total {} != java_only {}",
            emitted_keys.len(),
            java_only.len()
        );
        assert_eq!(rust_only_extra, 0, "rust_only_extra_vcf {rust_only_extra}");
    }
}
