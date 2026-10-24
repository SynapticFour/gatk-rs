#[path = "build/git_identity.rs"]
mod git_identity;

use std::path::PathBuf;

fn main() {
    let manifest_dir = match std::env::var_os("CARGO_MANIFEST_DIR") {
        Some(dir) => PathBuf::from(dir),
        None => {
            println!("cargo:rustc-env=GATK_RS_GIT_SHA=unknown");
            return;
        }
    };
    let sha = git_identity::git_head_sha(&manifest_dir);
    println!("cargo:rustc-env=GATK_RS_GIT_SHA={sha}");
    git_identity::emit_git_rerun_directives(&manifest_dir);
}
