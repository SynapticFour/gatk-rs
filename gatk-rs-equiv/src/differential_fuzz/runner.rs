//! Campaign loop: generate → compare Java vs Rust → shrink → fixture → issue.

use super::github;
use super::scenario::{scenario_from_bytes, scenario_to_seed_bytes, Scenario};
use super::shrink::shrink_scenario;
use super::synth::materialize_scenario;
use crate::compare::compare_callsets_with_ad_tol;
use crate::hc::{self, JavaGatk};
use crate::resources;
use crate::types::DirectCompare;
use anyhow::{bail, Context, Result};
use chrono::Utc;
use serde::Serialize;
use std::fs;
use std::path::{Path, PathBuf};

fn bytes_to_hex(bytes: &[u8]) -> String {
    bytes.iter().map(|b| format!("{b:02x}")).collect()
}

#[derive(Debug, Clone)]
pub struct DiffFuzzArgs {
    pub iterations: u32,
    pub out_dir: PathBuf,
    pub rust_binary: PathBuf,
    pub java_gatk_jar: Option<PathBuf>,
    pub java_gatk_bin: Option<PathBuf>,
    pub fixture_root: PathBuf,
    pub format_ad_tol: u32,
    pub open_github_issue: bool,
    pub shrink_steps: usize,
    pub seed_override: Option<u64>,
    pub min_free_gb: u64,
}

#[derive(Debug, Clone)]
pub struct EvalBins {
    pub rust_binary: PathBuf,
    pub java: JavaGatk,
    /// Max |AD_i − AD_i'| (and DP) allowed before FORMAT counts as divergence.
    pub format_ad_tol: u32,
}

#[derive(Debug, Clone, Serialize)]
pub struct Divergence {
    pub kind: String,
    pub summary: String,
    pub direct: DirectCompare,
    pub java_only_sites: u64,
    pub rust_only_sites: u64,
    pub format_mismatch_same_gt: u64,
    pub gt_mismatch: u64,
}

pub fn evaluate_scenario(
    scenario: &Scenario,
    bins: &EvalBins,
    work: &Path,
) -> Result<Option<Divergence>> {
    fs::create_dir_all(work)?;
    let mat = materialize_scenario(scenario, work)?;
    let java_vcf = work.join("java.vcf");
    let rust_vcf = work.join("rust.vcf");

    hc::run_java_hc(
        &bins.java,
        &mat.reference,
        &mat.bam,
        &java_vcf,
        Some(&mat.interval),
        2,
    )?;
    hc::run_rust_hc(
        &bins.rust_binary,
        &mat.reference,
        &mat.bam,
        &rust_vcf,
        Some(&mat.interval),
        2,
    )?;

    let direct = compare_callsets_with_ad_tol(&java_vcf, &rust_vcf, bins.format_ad_tol)?;
    Ok(classify_divergence(&direct))
}

fn classify_divergence(direct: &DirectCompare) -> Option<Divergence> {
    let allele_div = direct.java_only + direct.rust_only + direct.allele_match_gt_mismatch;
    let format_div = direct.format_mismatch_same_gt;
    if allele_div == 0 && format_div == 0 {
        return None;
    }
    let kind = if allele_div > 0 {
        "allele_or_gt".to_string()
    } else {
        "format_fields".to_string()
    };
    let summary = format!(
        "java_only={} rust_only={} gt_mismatch={} format_mismatch_same_gt={} identical={}",
        direct.java_only,
        direct.rust_only,
        direct.allele_match_gt_mismatch,
        direct.format_mismatch_same_gt,
        direct.identical_sites
    );
    Some(Divergence {
        kind,
        summary,
        direct: direct.clone(),
        java_only_sites: direct.java_only,
        rust_only_sites: direct.rust_only,
        format_mismatch_same_gt: direct.format_mismatch_same_gt,
        gt_mismatch: direct.allele_match_gt_mismatch,
    })
}

fn write_fixture(
    fixture_dir: &Path,
    scenario: &Scenario,
    div: &Divergence,
    work: &Path,
) -> Result<()> {
    fs::create_dir_all(fixture_dir)?;
    fs::copy(work.join("reference.fa"), fixture_dir.join("reference.fa"))?;
    let _ = fs::copy(
        work.join("reference.fa.fai"),
        fixture_dir.join("reference.fa.fai"),
    );
    fs::copy(work.join("reads.bam"), fixture_dir.join("reads.bam"))?;
    let _ = fs::copy(
        work.join("reads.bam.bai"),
        fixture_dir.join("reads.bam.bai"),
    );
    fs::copy(work.join("java.vcf"), fixture_dir.join("java.vcf"))?;
    fs::copy(work.join("rust.vcf"), fixture_dir.join("rust.vcf"))?;
    fs::write(
        fixture_dir.join("scenario.json"),
        serde_json::to_string_pretty(scenario)? + "\n",
    )?;
    fs::write(
        fixture_dir.join("seed_bytes.hex"),
        bytes_to_hex(&scenario_to_seed_bytes(scenario)) + "\n",
    )?;
    fs::write(
        fixture_dir.join("diverge.json"),
        serde_json::to_string_pretty(div)? + "\n",
    )?;
    fs::write(
        fixture_dir.join("README.md"),
        format!(
            "# HC differential regression fixture\n\n\
             - **Kind:** {}\n\
             - **Summary:** {}\n\
             - **Seed:** {}\n\
             - **Interval:** synth1:1-{}\n\n\
             Reproduce:\n\n\
             ```bash\n\
             cargo run -p gatk-rs-equiv -- differential-fuzz --replay-fixture {}\n\
             ```\n",
            div.kind,
            div.summary,
            scenario.seed,
            scenario.ref_len,
            fixture_dir.display()
        ),
    )?;
    Ok(())
}

pub fn run_campaign(args: DiffFuzzArgs) -> Result<i32> {
    resources::require_free_gb(&args.out_dir, args.min_free_gb)?;
    resources::apply_process_env(2);
    fs::create_dir_all(&args.out_dir)?;
    fs::create_dir_all(&args.fixture_root)?;

    let java = JavaGatk::resolve(args.java_gatk_jar.clone(), args.java_gatk_bin.clone())?;
    let bins = EvalBins {
        rust_binary: args.rust_binary.clone(),
        java,
        format_ad_tol: args.format_ad_tol,
    };

    let mut divergences = 0u32;
    let mut rng_seed = args.seed_override.unwrap_or_else(|| {
        use std::time::{SystemTime, UNIX_EPOCH};
        SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .map(|d| d.as_nanos() as u64)
            .unwrap_or(1)
    });

    for iter in 0..args.iterations {
        let mut bytes = rng_seed.to_le_bytes().to_vec();
        bytes.extend_from_slice(&(iter as u64).to_le_bytes());
        bytes.extend_from_slice(&((rng_seed >> 17) as u32).to_le_bytes());
        rng_seed = rng_seed.wrapping_mul(6364136223846793005).wrapping_add(1);

        let scenario = scenario_from_bytes(&bytes);
        let work = args.out_dir.join(format!("iter_{iter:04}"));
        eprintln!(
            "[diff-fuzz] iter={iter} seed={} reads={} reflen={}",
            scenario.seed, scenario.n_reads, scenario.ref_len
        );

        let Some(div0) = evaluate_scenario(&scenario, &bins, &work)? else {
            continue;
        };
        eprintln!("[diff-fuzz] DIVERGENCE: {} — {}", div0.kind, div0.summary);

        let shrink_root = work.join("shrink");
        let (mini, div) = shrink_scenario(scenario, &bins, &shrink_root, args.shrink_steps)?;
        let final_work = work.join("minimized");
        let _ = evaluate_scenario(&mini, &bins, &final_work)?;

        let stamp = Utc::now().format("%Y%m%dT%H%M%SZ");
        let hash = format!("{:08x}", (mini.seed as u32) ^ (mini.n_reads << 16));
        let fixture_dir = args.fixture_root.join(format!("{stamp}_{hash}"));
        write_fixture(&fixture_dir, &mini, &div, &final_work)?;
        eprintln!("[diff-fuzz] fixture → {}", fixture_dir.display());

        if args.open_github_issue {
            match github::open_parity_issue(&fixture_dir, &mini, &div) {
                Ok(url) => eprintln!("[diff-fuzz] GitHub issue: {url}"),
                Err(e) => eprintln!("[diff-fuzz] WARNING: gh issue failed: {e:#}"),
            }
        }

        divergences += 1;
    }

    eprintln!(
        "[diff-fuzz] done: iterations={} divergences={}",
        args.iterations, divergences
    );
    Ok(if divergences > 0 { 1 } else { 0 })
}

/// Replay a saved fixture directory (re-run both callers, re-compare).
pub fn replay_fixture(
    fixture_dir: &Path,
    rust_binary: &Path,
    java_jar: Option<PathBuf>,
    java_bin: Option<PathBuf>,
    format_ad_tol: u32,
) -> Result<i32> {
    let scenario: Scenario = serde_json::from_str(
        &fs::read_to_string(fixture_dir.join("scenario.json"))
            .with_context(|| format!("missing scenario.json in {}", fixture_dir.display()))?,
    )?;
    let java = JavaGatk::resolve(java_jar, java_bin)?;
    let bins = EvalBins {
        rust_binary: rust_binary.to_path_buf(),
        java,
        format_ad_tol,
    };
    let work = fixture_dir.join("replay");
    match evaluate_scenario(&scenario, &bins, &work)? {
        Some(div) => {
            eprintln!(
                "[diff-fuzz] REPLAY still diverges: {} — {}",
                div.kind, div.summary
            );
            Ok(1)
        }
        None => {
            eprintln!("[diff-fuzz] REPLAY clean (no divergence)");
            Ok(0)
        }
    }
}

pub fn run_from_cli(args: crate::cli::DiffFuzzCliArgs) -> Result<i32> {
    if let Some(fx) = &args.replay_fixture {
        return replay_fixture(
            fx,
            &args.rust_binary,
            args.java_gatk_jar.clone(),
            args.java_gatk_bin.clone(),
            args.format_ad_tol,
        );
    }
    if !args.rust_binary.exists() {
        bail!("--rust-binary not found: {}", args.rust_binary.display());
    }
    run_campaign(DiffFuzzArgs {
        iterations: args.iterations,
        out_dir: args.out,
        rust_binary: args.rust_binary,
        java_gatk_jar: args.java_gatk_jar,
        java_gatk_bin: args.java_gatk_bin,
        fixture_root: args.fixture_root,
        format_ad_tol: args.format_ad_tol,
        open_github_issue: args.open_github_issue,
        shrink_steps: args.shrink_steps,
        seed_override: args.seed,
        min_free_gb: args.min_free_gb,
    })
}
