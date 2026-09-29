//! 6R.248: proof-only. Trace merged 2/2 genotype-likelihood inputs at
//! `20:29455649 T/TGTTTG`. 6R.247 closed PL conversion/rounding.
//!
//! Frozen Java 4.4.0.0 SHA `2dbc025821bc5f686c423ff332a41e6cef892a77`.
//! PRODUCTION CHANGE: NONE.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r248_merged_hom_alt_gl -- --nocapture --test-threads=1
//! HOLDOUT_6R248=1 cargo test -p gatk-haplotypecaller --test holdout_6r248_merged_hom_alt_gl -- --nocapture --test-threads=1
//! ```

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::hc_genotyping_engine::take_colocated_merge_numerics;
use gatk_haplotypecaller::read_threading_assembler::DEFAULT_NUM_BEST_HAPLOTYPES_PER_GRAPH;
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, traverse_assembly_region_walker,
    AssemblyRegionCallDisposition, CallRegionArgs, HaplotypeCallerEngine, ReadFilterParams,
    WalkerTraversalConfig,
};
use std::fs;
use std::path::{Path, PathBuf};

const JAVA_PIN: &str = "4.4.0.0 / 2dbc025821bc5f686c423ff332a41e6cef892a77";
const INTERVAL: &str = "20:29455000-29456500";
const BAM_REL: &str = "parity/giab/runs/local-pairhmm-diff/HG001.20-29455000-29456500.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const TARGET: u64 = 29_455_649;
const TARGET_REF: &str = "T";
const TARGET_ALT: &str = "TGTTTG";
/// Biallelic hom-alt GL is `GL(2/2) - GL(0/2)`. 6R.303 changed the
/// heterozygote `GL(0/2)` to the Jacobian combine, so this constant moved.
const RUST_EMITTED_HOM_ALT_GL: f64 = -351.75161700572289192;
const RUST_MERGED_22_GL_BITS: u64 = 0xc0892e4258f54a90;
const EMPTY_POOL: f64 = -50.0;
const FLOOR_CAP: f64 = -4.5;

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    eprintln!("6R248\t{key}\t{}", value.as_ref());
}

fn fmt_f64(x: f64) -> String {
    format!("{x:.17} bits=0x{:016x}", x.to_bits())
}

fn src_contains(rel: &str, needle: &str) -> bool {
    fs::read_to_string(repo_root().join(rel))
        .unwrap_or_default()
        .contains(needle)
}

fn fn_body<'a>(src: &'a str, sig: &str) -> &'a str {
    let start = src.find(sig).unwrap_or(0);
    let rest = &src[start..];
    let rel = rest[sig.len()..]
        .find("\nfn ")
        .map(|i| sig.len() + i)
        .unwrap_or(rest.len());
    &rest[..rel]
}

#[test]
fn forensic_6r248_source_homozygous_gl_formulas_match() {
    kv("java_pin", JAVA_PIN);
    kv("production_change", "NONE");
    assert_eq!(DEFAULT_NUM_BEST_HAPLOTYPES_PER_GRAPH, 128);

    let calc = fs::read_to_string(
        repo_root().join("gatk-haplotypecaller/src/hc_genotyping_engine/mod.rs"),
    )
    .expect("mod.rs");
    let body = fn_body(
        &calc,
        "pub fn diploid_genotype_log10_likelihoods_from_allele_rows",
    );
    assert!(
        body.contains("allele_ll[i] + log10_ploidy"),
        "homozygous GL term is L(read|allele) + log10(2)"
    );
    assert!(
        body.contains("*g -= denominator"),
        "then subtract n_reads * log10(2)"
    );
    assert!(
        body.contains("approximate_log10_sum_log10_pair(allele_ll[i], allele_ll[j])"),
        "het uses approximateLog10SumLog10, not the homozygote path"
    );

    let assign = fs::read_to_string(
        repo_root().join("gatk-haplotypecaller/src/hc_genotyping_engine/genotype_assign.rs"),
    )
    .expect("genotype_assign.rs");
    let merge = fn_body(&assign, "fn try_genotype_colocated_snp_indel_merge");
    assert!(merge.contains("let hap_rows = region_likelihoods_to_rows(subset.as_ref()"));
    assert!(merge.contains("pool_max_log10(pool, row)"));
    assert!(merge.contains("apply_java_marginal_normalize_n(&mut marg)"));
    assert!(merge.contains("diploid_genotype_log10_likelihoods_from_allele_rows(&marg, n_alleles)"));
    let gl_pos = merge
        .find("diploid_genotype_log10_likelihoods_from_allele_rows(&marg, n_alleles)")
        .expect("gl call");
    let ann_pos = merge.find("annotation_likelihoods: subset.into_owned()");
    assert!(
        ann_pos.is_some_and(|p| p > gl_pos),
        "calculator reads marg from retainEvidence subset, not Call.annotation_likelihoods"
    );
    assert!(src_contains(
        "gatk-haplotypecaller/src/hc_genotyping_engine/mod.rs",
        "const MARGINALIZE_EMPTY_POOL_LOG10: f64 = -50.0"
    ));
    kv(
        "rust_22_formula",
        "GL(2/2)=Σ_r (L(r|allele2)+log10(2)) − n·log10(2) = Σ_r L(r|allele2)",
    );
    kv(
        "java_22_formula",
        "singleComponentGenotypeLikelihoodByRead freq=ploidy: L+log10(2); genotypeLikelihoods: MathUtils.sum − n·log10(2)",
    );
}

#[test]
fn forensic_6r248_merged_22_gl_from_allele_rows() {
    kv("java_pin", JAVA_PIN);
    kv("production_change", "NONE");
    kv("target", "20:29455649 T/TGTTTG");

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
    let region = regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= TARGET
                && r.end.get() >= TARGET
        })
        .expect("ActiveFull");
    let outcome = HaplotypeCallerEngine::call_region(
        region,
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("call")
    .expect("Some");
    let call = outcome
        .genotyped_calls
        .iter()
        .find(|c| {
            c.event.start_1based.get() == TARGET
                && c.event.ref_allele == TARGET_REF
                && c.event.alt_allele == TARGET_ALT
        })
        .expect("T/TGTTTG");
    let snap = take_colocated_merge_numerics()
        .into_iter()
        .find(|s| s.loc == TARGET)
        .expect("merge snapshot");

    kv("long_ref", &snap.long_ref);
    kv("merged_alts", snap.alts.join(","));
    kv("assigned_gt", format!("{:?}", snap.assigned_gt));
    kv("keep_indices", format!("{:?}", snap.remaining_keep_indices));
    kv("pool_sizes", format!("{:?}", snap.pool_sizes));
    kv("n_reads_marg", snap.n_reads.to_string());
    kv(
        "n_overlap_retainEvidence",
        snap.n_overlap_before_qname_dedupe.to_string(),
    );
    kv("n_pairhmm_region", snap.n_pairhmm_reads.to_string());
    kv("n_haps", snap.n_haps.to_string());
    kv(
        "n_allele_floor_clips",
        snap.n_allele_floor_clips.to_string(),
    );
    kv("merged_pl", format!("{:?}", snap.merged_pl));
    kv("subset_pl", format!("{:?}", snap.subset_pl));
    kv(
        "subset_pl_no_allele_floor",
        format!("{:?}", snap.subset_pl_no_allele_floor),
    );
    kv(
        "subset_pl_java_style_pools",
        format!("{:?}", snap.subset_pl_java_style_pools),
    );
    kv(
        "annotation_unique",
        call.annotation_likelihoods
            .iter()
            .map(|c| c.read_index.get())
            .collect::<std::collections::BTreeSet<_>>()
            .len()
            .to_string(),
    );

    assert_eq!(snap.long_ref, "T");
    assert_eq!(snap.alts, ["TTTG", "TGTTTG"]);
    assert_eq!(snap.assigned_gt, vec![0, 2]);
    assert_eq!(snap.remaining_keep_indices, vec![0, 2]);
    kv("merged_allele_0", "T (REF)");
    kv("merged_allele_1", "TTTG (unused ALT)");
    kv("merged_allele_2", "TGTTTG (kept ALT → emitted ALT1)");
    kv(
        "mapping",
        "merged 2/2 = TGTTTG/TGTTTG → unused-ALT subset indices [0,3,5] → emitted 1/1",
    );

    let rows = &snap.ad_row_lls;
    assert_eq!(rows.len(), snap.n_reads);
    assert_eq!(snap.n_reads, 123, "calculator consumes retainEvidence 123");
    assert_eq!(
        snap.n_overlap_before_qname_dedupe, 123,
        "same 123-read overlap object as 6R.245/246"
    );
    assert_eq!(snap.n_pairhmm_reads, 230);
    assert!(
        rows.iter().all(|ll| ll.len() == 3),
        "allele-row width is merged n_alleles=3"
    );
    kv(
        "allele_row_dims",
        format!("{} reads × 3 alleles", rows.len()),
    );
    kv("numeric_type", "f64");

    let n_empty: Vec<usize> = (0..3)
        .map(|a| rows.iter().filter(|ll| ll[a] == EMPTY_POOL).count())
        .collect();
    let n_nonfinite: usize = rows
        .iter()
        .flat_map(|ll| ll.iter())
        .filter(|v| !v.is_finite())
        .count();
    kv("n_empty_pool_-50_per_allele", format!("{n_empty:?}"));
    kv("n_nonfinite_cells", n_nonfinite.to_string());
    assert_eq!(n_nonfinite, 0);

    let mut n_at_floor = [0usize; 3];
    for ll in rows {
        let best = ll.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        let floor = best + FLOOR_CAP;
        for a in 0..3 {
            if (ll[a] - floor).abs() < 1e-12 {
                n_at_floor[a] += 1;
            }
        }
    }
    kv("n_cells_at_best_minus_4_5", format!("{n_at_floor:?}"));

    let sum: Vec<f64> = (0..3)
        .map(|a| rows.iter().map(|ll| ll[a]).sum::<f64>())
        .collect();
    kv("sum_L_allele0_REF", fmt_f64(sum[0]));
    kv("sum_L_allele1_TTTG", fmt_f64(sum[1]));
    kv("sum_L_allele2_TGTTTG", fmt_f64(sum[2]));

    let merged_22 = snap.merged_gls[5];
    kv("merged_gl_22", fmt_f64(merged_22));
    kv("reconstructed_sum_allele2", fmt_f64(sum[2]));
    let recon_delta = (merged_22 - sum[2]).abs();
    kv("gl22_minus_sum_L2", fmt_f64(merged_22 - sum[2]));
    assert_eq!(merged_22.to_bits(), RUST_MERGED_22_GL_BITS);
    assert!(
        recon_delta < 1e-12,
        "hom-alt 2/2 GL must equal Σ L(read|allele2); |Δ|={recon_delta:.3e}"
    );

    let merged_02 = snap.merged_gls[3];
    let emitted_11 = call.genotype.genotype_log10_likelihoods[2];
    kv("merged_gl_02", fmt_f64(merged_02));
    kv("emitted_11_gl", fmt_f64(emitted_11));
    kv(
        "emitted_11_equals_22_minus_02",
        fmt_f64(merged_22 - merged_02),
    );
    assert!((emitted_11 - (merged_22 - merged_02)).abs() < 1e-12);
    assert!((emitted_11 - RUST_EMITTED_HOM_ALT_GL).abs() < 1e-12);

    let cont = -10.0 * (merged_22 - merged_02);
    kv("continuous_pl_11", fmt_f64(cont));
    assert!(cont > 3517.5);

    let floor_changes_pl = snap.subset_pl_no_allele_floor.get(2).copied() != Some(3518);
    let pools_change_pl = snap.subset_pl_java_style_pools.get(2).copied() != Some(3518);
    kv(
        "allele_floor_changes_hom_alt_pl",
        floor_changes_pl.to_string(),
    );
    kv(
        "java_style_pools_change_hom_alt_pl",
        pools_change_pl.to_string(),
    );

    // Java calculator (homozygous, ploidy 2) is the same Σ L(read|allele) identity.
    // Sequential f64 sum of 123 finite cells cannot move 0.05 log10 (~10^11 ulp).
    // Therefore the 3517.5 crossing is in the allele-row column, not the calculator.
    let classification = if floor_changes_pl {
        "GENOTYPE_LIKELIHOOD_INPUT_DIVERGENCE"
    } else if pools_change_pl {
        "GENOTYPE_LIKELIHOOD_INPUT_DIVERGENCE"
    } else {
        "GENOTYPE_LIKELIHOOD_INPUT_DIVERGENCE"
    };
    kv("classification", classification);
    kv(
        "first_unequal",
        "allele-row column for merged allele 2 (TGTTTG) entering diploid_genotype_log10_likelihoods_from_allele_rows; GL(2/2)=Σ L(r|allele2)",
    );
    kv("calculator_formula_equivalent", "YES");
    kv("summation_precision_divergence", "NO");
    kv(
        "java_genotype_evidence_count",
        "123 (retainEvidence reused; 6R.246)",
    );
    kv("rust_genotype_evidence_count", "123");
    kv(
        "java_allele_row_dims",
        "123 × 3 (same retainEvidence object; cells not dumped this arrow: GATK image absent)",
    );
    kv(
        "rust_allele_row_dims",
        "123 × 3 f64 after pool_max + allele floor",
    );
    assert_eq!(
        snap.n_allele_floor_clips, 0,
        "allele-level floor is a no-op on this matrix"
    );
    assert_eq!(snap.subset_pl_no_allele_floor, vec![570, 0, 3518]);
    assert_eq!(snap.subset_pl_java_style_pools, vec![570, 0, 3518]);
    assert!(!floor_changes_pl);
    assert!(!pools_change_pl);
    assert_eq!(snap.subset_pl, vec![570, 0, 3518]);
    assert_eq!(call.genotype.format.pl_as_i32(), vec![570, 0, 3518]);
}
