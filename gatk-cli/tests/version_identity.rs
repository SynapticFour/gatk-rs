//! `--version` and `-V` share one build identity line.

use assert_cmd::Command;
use std::path::Path;
use std::process::Command as StdCommand;

#[allow(dead_code)]
#[path = "../build/git_identity.rs"]
mod git_identity;

fn identity_line(stdout: &[u8]) -> String {
    let text = String::from_utf8(stdout.to_vec()).unwrap();
    assert_eq!(
        text.lines().count(),
        1,
        "version output must be one line:\n{text}"
    );
    text.trim_end_matches(['\r', '\n']).to_string()
}

fn sha_from_line(line: &str) -> &str {
    let prefix = format!("gatk-rs {} (git ", env!("CARGO_PKG_VERSION"));
    assert!(line.starts_with(&prefix), "{line}");
    let rest = line.strip_prefix(&prefix).unwrap();
    let sha = rest.strip_suffix(')').unwrap_or_else(|| panic!("{line}"));
    assert!(
        sha == "unknown" || is_full_sha(sha),
        "identity suffix must be a full sha or unknown: {sha}"
    );
    sha
}

fn is_full_sha(sha: &str) -> bool {
    (sha.len() == 40 || sha.len() == 64)
        && sha.bytes().all(|b| matches!(b, b'0'..=b'9' | b'a'..=b'f'))
}

#[test]
fn version_flags_print_the_same_identity_line() {
    let long = Command::cargo_bin("gatk-rs")
        .unwrap()
        .arg("--version")
        .assert()
        .success();
    let short = Command::cargo_bin("gatk-rs")
        .unwrap()
        .arg("-V")
        .assert()
        .success();
    assert_eq!(
        long.get_output().stdout,
        short.get_output().stdout,
        "--version and -V must be byte-identical"
    );
    let line = identity_line(&long.get_output().stdout);
    let sha = sha_from_line(&line);
    let head = StdCommand::new("git")
        .args(["rev-parse", "HEAD"])
        .current_dir(Path::new(env!("CARGO_MANIFEST_DIR")).join(".."))
        .output()
        .unwrap();
    assert!(head.status.success());
    let expected = String::from_utf8(head.stdout).unwrap();
    assert_eq!(sha, expected.trim());
}

#[test]
fn help_keeps_disclaimer_and_guide_keeps_pinned_oracle() {
    let help = Command::cargo_bin("gatk-rs")
        .unwrap()
        .arg("--help")
        .assert()
        .success();
    let text = String::from_utf8(help.get_output().stdout.clone()).unwrap();
    assert!(
        text.contains("not affiliated with, endorsed by, or supported by the Broad Institute"),
        "disclaimer must stay on --help"
    );
    let guide = std::fs::read_to_string(
        Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("..")
            .join("docs/USER_GUIDE.md"),
    )
    .unwrap();
    assert!(
        guide.contains("GATK **4.4.0.0**"),
        "pinned Java oracle version must stay in the user guide"
    );
}

#[test]
fn missing_git_metadata_is_unknown() {
    let dir = tempfile::tempdir().unwrap();
    assert_eq!(git_identity::git_head_sha(dir.path()), "unknown");

    let broken = tempfile::tempdir().unwrap();
    std::fs::create_dir(broken.path().join(".git")).unwrap();
    assert_eq!(git_identity::git_head_sha(broken.path()), "unknown");
}
