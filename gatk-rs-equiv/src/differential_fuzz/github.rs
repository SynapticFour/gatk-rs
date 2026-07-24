//! Open a GitHub issue for a minimized parity divergence (via `gh` CLI).

use super::runner::Divergence;
use super::scenario::Scenario;
use anyhow::{bail, Context, Result};
use std::path::Path;
use std::process::Command;

pub const LABEL: &str = "parity-divergence";

fn ensure_label() -> Result<()> {
    let check = Command::new("gh")
        .args([
            "label", "list", "--search", LABEL, "--json", "name", "-q", ".[].name",
        ])
        .output()
        .context("gh label list (is `gh` installed and authenticated?)")?;
    if !check.status.success() {
        bail!(
            "gh label list failed: {}",
            String::from_utf8_lossy(&check.stderr)
        );
    }
    let names = String::from_utf8_lossy(&check.stdout);
    if names.lines().any(|l| l.trim() == LABEL) {
        return Ok(());
    }
    let create = Command::new("gh")
        .args([
            "label",
            "create",
            LABEL,
            "--description",
            "Java GATK4 vs gatk-rs HaplotypeCaller divergence (differential fuzzer)",
            "--color",
            "D93F0B",
        ])
        .output()
        .context("gh label create")?;
    // Ignore "already exists" races.
    if !create.status.success() {
        let err = String::from_utf8_lossy(&create.stderr);
        if !err.to_ascii_lowercase().contains("already exists") {
            bail!("gh label create failed: {err}");
        }
    }
    Ok(())
}

/// Create an issue with the fixture path and divergence summary. Returns the issue URL.
pub fn open_parity_issue(
    fixture_dir: &Path,
    scenario: &Scenario,
    div: &Divergence,
) -> Result<String> {
    ensure_label()?;

    let title = format!(
        "parity-divergence: {} (seed={}, reads={}, reflen={})",
        div.kind, scenario.seed, scenario.n_reads, scenario.ref_len
    );

    let body = format!(
        "Automated report from `gatk-rs-equiv differential-fuzz` / `fuzz/run_hc_differential.sh`.\n\n\
         ## Summary\n\n\
         - **Kind:** `{kind}`\n\
         - **Details:** {summary}\n\
         - **Seed:** `{seed}`\n\
         - **Reads / ref_len / read_len:** {n_reads} / {ref_len} / {read_len}\n\
         - **java_only / rust_only / gt_mismatch / format_mismatch:** \
           {java_only} / {rust_only} / {gt_mismatch} / {format_mismatch}\n\n\
         ## Minimal fixture\n\n\
         Path (in-repo):\n\n\
         ```\n{fixture}\n```\n\n\
         Contents: `scenario.json`, `reference.fa`, `reads.bam`, `java.vcf`, `rust.vcf`, `diverge.json`.\n\n\
         ## Reproduce\n\n\
         ```bash\n\
         cargo run -p gatk-rs-equiv -- differential-fuzz \\\n\
           --replay-fixture {fixture} \\\n\
           --rust-binary ./target/release/gatk-rs\n\
         ```\n\n\
         Label: `{label}`\n",
        kind = div.kind,
        summary = div.summary,
        seed = scenario.seed,
        n_reads = scenario.n_reads,
        ref_len = scenario.ref_len,
        read_len = scenario.read_len,
        java_only = div.java_only_sites,
        rust_only = div.rust_only_sites,
        gt_mismatch = div.gt_mismatch,
        format_mismatch = div.format_mismatch_same_gt,
        fixture = fixture_dir.display(),
        label = LABEL,
    );

    let output = Command::new("gh")
        .args([
            "issue", "create", "--title", &title, "--body", &body, "--label", LABEL,
        ])
        .output()
        .context("gh issue create")?;
    if !output.status.success() {
        bail!(
            "gh issue create failed: {}",
            String::from_utf8_lossy(&output.stderr)
        );
    }
    let url = String::from_utf8_lossy(&output.stdout).trim().to_string();
    if url.is_empty() {
        bail!("gh issue create succeeded but returned empty URL");
    }
    Ok(url)
}
