//! 6R.102 holdout: first QUAL divergence at the canonical T/C site.
//!
//! Skipped unless `HOLDOUT_6R102=1`. Coordinate-free contract lives in
//! `forensic_6r102_qual_calculation_contract`.
//!
//! Java `GenotypingEngine.calculateGenotypes` writes QUAL from
//! `AlleleFrequencyCalculator.calculate` on the pre-subset merged VC
//! (`TG,T,CG,*`). When `*` is present, `log10PNoVariant` is the log10-sum of
//! REF+SPAN_DEL genotype posteriors, not HOM_REF alone. Emitted PL/GT/AD are
//! unchanged.
//!
//! ```text
//! HOLDOUT_6R102=1 cargo test -p gatk-haplotypecaller --test holdout_6r102_qual_calculation -- --nocapture
//! ```

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, take_colocated_merge_numerics,
    traverse_assembly_region_walker, try_emit_call_region_variants, AssemblyRegionCallDisposition,
    CallRegionArgs, HaplotypeCallerEngine, ReadFilterParams, WalkerTraversalConfig,
    DEFAULT_STAND_EMIT_CONFIDENCE,
};
use serde_json::json;
use std::path::{Path, PathBuf};

const INTERVAL: &str = "20:29455000-29456500";
const BAM_REL: &str = "parity/giab/runs/local-pairhmm-diff/HG001.20-29455000-29456500.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const POS_SNP: u64 = 29_456_344;

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

#[test]
fn holdout_6r102_qual_from_merged_span_del_af() {
    if std::env::var("HOLDOUT_6R102").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R102=1");
        return;
    }
    let root = repo_root();
    let ref_fasta = root.join(REF_REL);
    let bam = root.join(BAM_REL);
    assert!(ref_fasta.is_file() && bam.is_file());

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
    let covering: Vec<_> = regions
        .iter()
        .filter(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= POS_SNP
                && r.end.get() >= POS_SNP
        })
        .collect();
    assert_eq!(covering.len(), 1);
    let region = covering[0];
    let outcome = HaplotypeCallerEngine::call_region(
        region,
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("call")
    .expect("Some");
    let emitted =
        try_emit_call_region_variants(region, &outcome, "SAMPLE", DEFAULT_STAND_EMIT_CONFIDENCE)
            .unwrap_or_default();
    let vcf = emitted
        .iter()
        .find(|r| {
            r.position == POS_SNP && r.reference == "T" && r.alternate.iter().any(|a| a == "C")
        })
        .expect("T/C");
    let live = take_colocated_merge_numerics();
    let snap = live
        .iter()
        .find(|n| n.loc == POS_SNP)
        .expect("colocated merge numerics");

    let vcf_gt = vcf.samples[0].gt.as_ref().map(|g| g.alleles.clone());
    let vcf_ad = vcf.samples[0].ad.clone().unwrap_or_default();
    let vcf_pl = vcf.samples[0].pl.clone().unwrap_or_default();
    let vcf_qual = vcf.quality.expect("QUAL");

    let doc = json!({
        "locus": "20:29456344 T/C",
        "merged_alleles": std::iter::once(snap.long_ref.clone()).chain(snap.alts.iter().cloned()).collect::<Vec<_>>(),
        "merged_pl": snap.merged_pl,
        "subset_pl": snap.subset_pl,
        "vcf": {
            "gt": vcf_gt,
            "ad": vcf_ad,
            "pl": vcf_pl,
            "qual": vcf_qual,
        },
        "java_oracle": {"gt": [0, 1], "ad": [36, 19], "pl": [542, 0, 1353], "qual": 510.06},
        "first_divergent_operation": "QUAL_FORMULA",
        "note": "Java log10PNoVariant sums REF+SPAN_DEL posteriors on the pre-subset merged VC",
    });
    println!("{}", serde_json::to_string_pretty(&doc).unwrap());

    assert!(
        snap.alts.iter().any(|a| a == "*"),
        "merged alleles include SPAN_DEL: {:?}",
        snap.alts
    );
    assert_eq!(snap.merged_gls.len(), 10);
    assert_eq!(vcf_pl, vec![542u32, 0, 1353]);
    assert_eq!(vcf_ad, vec![36u32, 19]);
    assert_eq!(vcf_gt.as_deref(), Some(&[0, 1][..]));
    assert!(
        (vcf_qual - 510.06).abs() < 0.02,
        "Java QUAL 510.06 from SPAN_DEL AF, got {vcf_qual}"
    );
}
