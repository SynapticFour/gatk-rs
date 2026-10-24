//! Build-time Git HEAD for the `gatk-rs` version line.
//!
//! Looks only at this package directory and its parent (the workspace root).
//! A missing Git directory, a missing `git` binary, or a failed `rev-parse`
//! becomes the literal `unknown`. Nothing from the environment, the worktree
//! dirty state, or the machine is embedded.

use std::fs;
use std::path::{Component, Path, PathBuf};
use std::process::{Command, Stdio};

pub fn git_head_sha(manifest_dir: &Path) -> String {
    let Some(dot_git) = find_dot_git(manifest_dir) else {
        return unknown();
    };
    let Some(git_dir) = resolve_git_dir(&dot_git) else {
        return unknown();
    };
    match rev_parse_head(&git_dir) {
        Some(sha) if is_full_sha(&sha) => sha,
        _ => unknown(),
    }
}

/// Ask Cargo to rebuild when HEAD or the branch ref it names changes.
pub fn emit_git_rerun_directives(manifest_dir: &Path) {
    let Some(dot_git) = find_dot_git(manifest_dir) else {
        for candidate in dot_git_candidates(manifest_dir) {
            println!("cargo:rerun-if-changed={}", candidate.display());
        }
        return;
    };
    if dot_git.is_file() {
        println!("cargo:rerun-if-changed={}", dot_git.display());
    }
    let Some(git_dir) = resolve_git_dir(&dot_git) else {
        println!("cargo:rerun-if-changed={}", dot_git.display());
        return;
    };
    let head_path = git_dir.join("HEAD");
    println!("cargo:rerun-if-changed={}", head_path.display());
    let Ok(head) = fs::read_to_string(&head_path) else {
        return;
    };
    let Some(rel) = head.trim().strip_prefix("ref: ") else {
        return;
    };
    let rel = rel.trim();
    if !is_safe_git_ref(rel) {
        return;
    }
    emit_ref(&git_dir, rel);
    if let Ok(common) = fs::read_to_string(git_dir.join("commondir")) {
        let common = common.trim();
        if !common.is_empty() && !common.contains(['\n', '\r']) {
            emit_ref(&git_dir.join(common), rel);
        }
    }
}

fn unknown() -> String {
    "unknown".to_string()
}

fn dot_git_candidates(manifest_dir: &Path) -> Vec<PathBuf> {
    let mut paths = vec![manifest_dir.join(".git")];
    if let Some(parent) = manifest_dir.parent() {
        paths.push(parent.join(".git"));
    }
    paths
}

fn find_dot_git(manifest_dir: &Path) -> Option<PathBuf> {
    dot_git_candidates(manifest_dir)
        .into_iter()
        .find(|path| path.exists())
}

fn resolve_git_dir(dot_git: &Path) -> Option<PathBuf> {
    if dot_git.is_dir() {
        return Some(dot_git.to_path_buf());
    }
    if !dot_git.is_file() {
        return None;
    }
    let text = fs::read_to_string(dot_git).ok()?;
    let rest = text.trim().strip_prefix("gitdir:")?.trim();
    if rest.is_empty() || rest.contains(['\n', '\r']) {
        return None;
    }
    let pointed = PathBuf::from(rest);
    let resolved = if pointed.is_absolute() {
        pointed
    } else {
        dot_git.parent()?.join(pointed)
    };
    resolved.is_dir().then_some(resolved)
}

fn rev_parse_head(git_dir: &Path) -> Option<String> {
    let output = Command::new("git")
        .arg("--git-dir")
        .arg(git_dir)
        .args(["rev-parse", "HEAD"])
        .stdin(Stdio::null())
        .stderr(Stdio::null())
        .output()
        .ok()?;
    if !output.status.success() {
        return None;
    }
    let text = String::from_utf8(output.stdout).ok()?;
    Some(text.trim().to_string())
}

fn is_full_sha(sha: &str) -> bool {
    (sha.len() == 40 || sha.len() == 64)
        && sha.bytes().all(|b| matches!(b, b'0'..=b'9' | b'a'..=b'f'))
}

fn is_safe_git_ref(rel: &str) -> bool {
    let path = Path::new(rel);
    !rel.is_empty()
        && path
            .components()
            .all(|component| matches!(component, Component::Normal(_)))
}

fn emit_ref(git_dir: &Path, rel: &str) {
    let loose = git_dir.join(rel);
    if loose.is_file() {
        println!("cargo:rerun-if-changed={}", loose.display());
    } else {
        println!(
            "cargo:rerun-if-changed={}",
            git_dir.join("packed-refs").display()
        );
    }
}
