//! 6R.159: hom-ref emission boundary on the frozen 20:29456196 A/T site.
//! Diagnostic-only `unbounded_diagnostic`. Skipped unless `HOLDOUT_6R159=1`.
//!
//! ```text
//! HOLDOUT_6R159=1 GATK_RS_EXPERIMENTAL_KBEST_POLICY=unbounded_diagnostic \
//!   cargo test -p gatk-haplotypecaller --test holdout_6r159_homref_emission -- --nocapture --test-threads=1
//! ```

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::emit_gates::{
    java_emit_would_pass, passes_hc_variant_emit_biallelic, passes_java_emit_not_hom_ref,
};
use gatk_haplotypecaller::event_map::VariationEvent;
use gatk_haplotypecaller::genotyping::{
    best_biallelic_diploid_genotype_index, emit_genotype_format_fields,
};
use gatk_haplotypecaller::hc_genotyping_engine::{
    java_emit_af_decision, l9_may_overwrite_pairhmm_gls_after_emit_fail,
    DEFAULT_STAND_EMIT_CONFIDENCE,
};
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, traverse_assembly_region_walker,
    AssemblyRegionCallDisposition, CallRegionArgs, GenomePosition, HaplotypeCallerEngine,
    ReadFilterParams, SparsePlShape, WalkerTraversalConfig,
};
use std::path::{Path, PathBuf};

const INTERVAL: &str = "20:29455000-29456500";
const BAM_REL: &str = "parity/giab/runs/local-pairhmm-diff/HG001.20-29455000-29456500.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const TARGET: u64 = 29_456_196;
const JAVA_ACTIVE_START: u64 = 29_455_995;
const JAVA_ACTIVE_END: u64 = 29_456_293;
const JAVA_PAD_START: u64 = 29_455_895;
const JAVA_PAD_END: u64 = 29_456_393;
const MERGED_REF: &str = "A";
const MERGED_ALT: &str = "T";
const JAVA_STAND_CALL_CONF: f64 = 30.0;
const FROZEN_GL: [f64; 3] = [0.0, -0.5, -117.4];
const FROZEN_AD: [i32; 2] = [37, 4];
const A3_PILEUP: [i32; 2] = [22, 24];

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    println!("6R159\t{key}\t{}", value.as_ref());
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

#[test]
fn holdout_6r159_homref_emission_boundary() {
    if std::env::var("HOLDOUT_6R159").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R159=1");
        return;
    }
    assert_ne!(
        std::env::var("GATK_RS_DIAGNOSTIC_SKIP_SITE_RESHAPE")
            .ok()
            .as_deref(),
        Some("1"),
        "6R.159 must not alter production via the diagnostic skip"
    );
    let _policy = EnvGuard::set("GATK_RS_EXPERIMENTAL_KBEST_POLICY", "unbounded_diagnostic");
    kv(
        "parity_mode",
        "unbounded_diagnostic\tnote=diagnostic_not_production_default",
    );
    kv("variant", format!("20:{TARGET} {MERGED_REF}/{MERGED_ALT}"));
    kv(
        "frozen_pre_emit",
        "GT=0/0 PL=0,5,1174 AD=37,4 GQ=5 DP=41 GLs=0.0,-0.5,-117.4",
    );

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

    let prod = outcome.genotyped_calls.iter().find(|c| {
        c.event.start_1based.get() == TARGET
            && c.event.ref_allele == MERGED_REF
            && c.event.alt_allele == MERGED_ALT
    });
    kv(
        "production_genotyped_calls",
        if prod.is_some() { "present" } else { "absent" },
    );
    kv(
        "genotyped_calls_semantics",
        "sites that survived GenotypeFinalize::finalize_site (Java returnCalls / call!=null), not all calculator results",
    );

    let event = VariationEvent::from_alleles("20", TARGET, MERGED_REF, MERGED_ALT);
    let fmt = emit_genotype_format_fields(&FROZEN_GL, &FROZEN_AD).expect("fmt");
    let best = best_biallelic_diploid_genotype_index(&FROZEN_GL, &FROZEN_AD);
    let not_hom_ref = passes_java_emit_not_hom_ref(&FROZEN_GL, &fmt);
    let af10 = java_emit_af_decision(&FROZEN_GL, DEFAULT_STAND_EMIT_CONFIDENCE).expect("af10");
    let af30 = java_emit_af_decision(&FROZEN_GL, JAVA_STAND_CALL_CONF).expect("af30");
    let rust_emit =
        java_emit_would_pass(&event, &FROZEN_GL, &fmt, DEFAULT_STAND_EMIT_CONFIDENCE, &[]).unwrap();
    let java_emit =
        java_emit_would_pass(&event, &FROZEN_GL, &fmt, JAVA_STAND_CALL_CONF, &[]).unwrap();
    let af_emit_10 =
        passes_hc_variant_emit_biallelic(&FROZEN_GL, DEFAULT_STAND_EMIT_CONFIDENCE).unwrap();
    let l9 = l9_may_overwrite_pairhmm_gls_after_emit_fail(&event, A3_PILEUP[0], A3_PILEUP[1], true);
    let first_fail = if !not_hom_ref {
        "passes_java_emit_not_hom_ref"
    } else if !af_emit_10 {
        "passes_hc_variant_emit_biallelic"
    } else {
        "none"
    };

    kv("calculator_best_gt_index", format!("{best}"));
    kv("passes_java_emit_not_hom_ref", format!("{not_hom_ref}"));
    kv(
        "gq_gate_in_java_emit_would_pass",
        "false\tnote=GQ is not consulted by java_emit_would_pass",
    );
    kv(
        "af10",
        format!(
            "mono={} plausible={} phred={:.4} passes_emit={}",
            af10.site_is_monomorphic, af10.alt_plausible, af10.phred_scaled, af10.passes_emit
        ),
    );
    kv(
        "af30",
        format!(
            "mono={} plausible={} phred={:.4} passes_emit={}",
            af30.site_is_monomorphic, af30.alt_plausible, af30.phred_scaled, af30.passes_emit
        ),
    );
    kv("java_emit_would_pass_stand10", format!("{rust_emit}"));
    kv("java_emit_would_pass_stand30", format!("{java_emit}"));
    kv("first_failing_predicate", first_fail);
    kv(
        "l9_overwrite",
        format!(
            "{l9} pileup={},{} shape={:?}",
            A3_PILEUP[0],
            A3_PILEUP[1],
            SparsePlShape::from_pileup_depths(A3_PILEUP[0], A3_PILEUP[1])
        ),
    );
    kv("classification", "HOMREF_GT_VETO_BEFORE_AF_EMIT");
    kv("production_change", "NONE");

    assert!(
        prod.is_none(),
        "canonical site must not enter genotyped_calls"
    );
    assert_eq!(best, 0);
    assert!(!not_hom_ref);
    assert_eq!(first_fail, "passes_java_emit_not_hom_ref");
    assert!(!rust_emit);
    assert!(!java_emit);
    assert!(af30.site_is_monomorphic);
    assert!(!af30.passes_emit);
    assert!(!l9);
    assert_eq!(fmt.gq.as_i32(), 5);
    assert_eq!(fmt.pl_as_i32(), vec![0, 5, 1174]);
}
