# Contributing

Thanks for helping with gatk-rs. Keep claims honest: only assert what
[`docs/CLAIM_MATRIX.md`](docs/CLAIM_MATRIX.md) lists as **Yes**.

By participating, you agree to the [`CODE_OF_CONDUCT.md`](CODE_OF_CONDUCT.md).
Report security issues privately per [`SECURITY.md`](SECURITY.md) — do not open
public issues for vulnerabilities.

## Setup

First step after clone — enable the versioned Git hooks (`.githooks/`):

```bash
./scripts/dev/install-hooks.sh
```

This sets `git config core.hooksPath .githooks` so every contributor shares the same
pre-commit checks (fmt, scoped clippy, unwrap/clone ratchets, doc links, >1 MiB guard).

Optional additional hooks via the Python `pre-commit` framework:

```bash
pip install pre-commit   # if needed
pre-commit install
```

## Build

```bash
cargo build -p gatk-cli
# lean hosts (e.g. 16GB):
CARGO_BUILD_JOBS=1 cargo build -p gatk-cli --release
```

Pinned Java oracle for differential work: [`docs/GATK_PINNED.env`](docs/GATK_PINNED.env).

## Test

```bash
cargo test --workspace
# optional harness / L2 dump surface:
cargo test -p gatk-haplotypecaller --features parity_harness
```

Malformed-input and Rayon determinism gates run in CI (`quality.yml` / `ci.yml`).

## Equivalence tooling

| Entry | Notes |
|-------|--------|
| [`gatk-rs-equiv/`](gatk-rs-equiv/) | `cargo run -p gatk-rs-equiv -- --help` |
| [`scripts/parity/`](scripts/parity/) | L2 / P12 / GIAB scripts used by CI |
| [`fuzz/run_hc_differential.sh`](fuzz/run_hc_differential.sh) | Differential fuzz wrapper |
| [`tools/equivalence/README.md`](tools/equivalence/README.md) | Index |

Do not widen P12 bands or market genome-wide equivalence without updating the claim matrix.

## Pull requests

1. Branch from `main`.
2. Keep changes focused; prefer algorithm fixes with a regression test or parity fixture.
3. Run `cargo test --workspace` (and relevant parity scripts if you touch HC emit/genotyping).
4. Update [`docs/CLAIM_MATRIX.md`](docs/CLAIM_MATRIX.md) if you change what the project asserts.
5. Open a PR with a short summary of *why* and how you verified it.

## Release process

Before **every** version tag (`v0.1.0`, `v0.2.0`, …) or GitHub Release:

1. Complete [`docs/RELEASE_CHECKLIST.md`](docs/RELEASE_CHECKLIST.md) against the tag SHA (dashboard freshness, claim-matrix audit, open gates, hygiene).
2. Record evidence paths / CI run URLs in the checklist copy attached to the release PR or notes.
3. Only then: `git tag` / publish the release.
4. Public announcement, profile pins, and external links are **optional and separate** — allowed only when the checklist’s “Public announcement” verdict is Yes.

Do not skip the checklist for “small” tags: the same process applies every time.

Architecture overview: [`docs/ARCHITECTURE.md`](docs/ARCHITECTURE.md).  
Independence / trademark / third-party data: [`NOTICE.md`](NOTICE.md).
