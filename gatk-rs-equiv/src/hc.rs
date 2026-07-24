//! Invoke Java GATK4 and gatk-rs HaplotypeCaller on identical inputs.

use anyhow::{bail, Context, Result};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

#[derive(Debug, Clone)]
pub enum JavaGatk {
    Jar(PathBuf),
    Bin(PathBuf),
}

impl JavaGatk {
    pub fn resolve(jar: Option<PathBuf>, bin: Option<PathBuf>) -> Result<Self> {
        if let Some(j) = jar {
            if !j.is_file() {
                bail!("--java-gatk-jar not found: {}", j.display());
            }
            return Ok(Self::Jar(j));
        }
        if let Some(b) = bin {
            if !b.exists() {
                bail!("--java-gatk-bin not found: {}", b.display());
            }
            return Ok(Self::Bin(b));
        }
        // Fall back to `gatk` on PATH.
        if which("gatk").is_some() {
            return Ok(Self::Bin(PathBuf::from("gatk")));
        }
        bail!(
            "Provide --java-gatk-jar or --java-gatk-bin (or install `gatk` on PATH). \
             Pin: docs/GATK_PINNED.env (GATK 4.4.0.0)."
        );
    }
}

fn which(name: &str) -> Option<PathBuf> {
    std::env::var_os("PATH").and_then(|paths| {
        for dir in std::env::split_paths(&paths) {
            let candidate = dir.join(name);
            if candidate.is_file() {
                return Some(candidate);
            }
        }
        None
    })
}

pub fn run_java_hc(
    java: &JavaGatk,
    reference: &Path,
    bam: &Path,
    out_vcf: &Path,
    interval: Option<&str>,
    threads: u32,
) -> Result<()> {
    if let Some(parent) = out_vcf.parent() {
        std::fs::create_dir_all(parent)?;
    }
    let mut cmd = match java {
        JavaGatk::Jar(jar) => {
            let mut c = Command::new("java");
            // Explicit heap even if JAVA_TOOL_OPTIONS is set — keeps GATK from grabbing the host.
            c.arg(format!("-Xmx{}", crate::resources::DEFAULT_JAVA_XMX))
                .arg("-jar")
                .arg(jar);
            c
        }
        JavaGatk::Bin(bin) => Command::new(bin),
    };
    cmd.arg("HaplotypeCaller")
        .arg("-R")
        .arg(reference)
        .arg("-I")
        .arg(bam)
        .arg("-O")
        .arg(out_vcf)
        .arg("--verbosity")
        .arg("ERROR")
        .arg("--native-pair-hmm-threads")
        .arg(threads.to_string());
    if let Some(iv) = interval {
        cmd.arg("-L").arg(iv);
    }
    cmd.env("RAYON_NUM_THREADS", threads.to_string());
    eprintln!("[gatk-rs-equiv] Java HC: {:?}", cmd);
    let status = cmd
        .stdout(Stdio::inherit())
        .stderr(Stdio::inherit())
        .status()
        .with_context(|| "failed to spawn Java GATK HaplotypeCaller")?;
    if !status.success() {
        bail!("Java HaplotypeCaller failed with {status}");
    }
    Ok(())
}

pub fn run_rust_hc(
    rust_bin: &Path,
    reference: &Path,
    bam: &Path,
    out_vcf: &Path,
    interval: Option<&str>,
    threads: u32,
) -> Result<()> {
    if !rust_bin.exists() {
        bail!("--rust-binary not found: {}", rust_bin.display());
    }
    if let Some(parent) = out_vcf.parent() {
        std::fs::create_dir_all(parent)?;
    }
    let mut cmd = Command::new(rust_bin);
    cmd.arg("HaplotypeCaller")
        .arg("-R")
        .arg(reference)
        .arg("-I")
        .arg(bam)
        .arg("-O")
        .arg(out_vcf);
    if let Some(iv) = interval {
        cmd.arg("-L").arg(iv);
    }
    cmd.env("RAYON_NUM_THREADS", threads.to_string());
    eprintln!("[gatk-rs-equiv] Rust HC: {:?}", cmd);
    let status = cmd
        .stdout(Stdio::inherit())
        .stderr(Stdio::inherit())
        .status()
        .with_context(|| "failed to spawn gatk-rs HaplotypeCaller")?;
    if !status.success() {
        bail!("gatk-rs HaplotypeCaller failed with {status}");
    }
    Ok(())
}
