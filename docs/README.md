# Documentation map

Read these first. Everything else is evidence, history, or a specialist note.

| Document | Role |
|----------|------|
| [`../README.md`](../README.md) | What the project is, what it claims, how to build it |
| [`USER_GUIDE.md`](USER_GUIDE.md) | CLI, inputs, outputs, configuration, failures |
| [`ARCHITECTURE.md`](ARCHITECTURE.md) | Crates, pipeline stages, determinism, extension points |
| [`PARITY.md`](PARITY.md) | How validation is structured and what current results mean |
| [`CLAIM_MATRIX.md`](CLAIM_MATRIX.md) | The only product-claim authority |
| [`PILOT_GUIDE.md`](PILOT_GUIDE.md) | Run Java GATK 4.4 and gatk-rs on your own BAMs and compare |
| [`../CONTRIBUTING.md`](../CONTRIBUTING.md) | Hooks, tests, pull requests |
| [`RELEASE_CHECKLIST.md`](RELEASE_CHECKLIST.md) | Gate before any version tag |

## Kept on purpose, not part of the first path

| Path | What it is |
|------|------------|
| [`parity/`](parity/) | 6R/7R investigation notes and k-best policy. Not acceptance by themselves. Index: [`parity/README.md`](parity/README.md). |
| [`PARITY_MILESTONE_6R.md`](PARITY_MILESTONE_6R.md) | Short 2026-08-30 mid-B milestone. The site log lives in [`PARITY.md`](PARITY.md). |
| [`adr/`](adr/) | Scope boundary and removal of a generic toolkit crate. |
| [`../parity/`](../parity/) | Harness, fixtures, acceptance-class vocabulary. |
| [`perf/`](perf/) | Performance notes. A public speed or memory claim still requires the dedicated host in [`ci/PERF_BENCHMARK_HOST.md`](ci/PERF_BENCHMARK_HOST.md). |
| [`EQUIVALENCE_DASHBOARD.md`](EQUIVALENCE_DASHBOARD.md) | Nightly markdown page. Its F1 table is unpopulated, so this file is not a measured result. |
| [`parity-site/`](parity-site/) | GitHub Pages source. [`parity-site/data/history.json`](parity-site/data/history.json) contains historical runs. Those rows are not the claim matrix’s signed HG001 ΔF1 record. Product claims stay in [`CLAIM_MATRIX.md`](CLAIM_MATRIX.md). |

Historical L6–L14 narratives that are **not** reproducible as Yes rows from
`main` are called out in [`CLAIM_MATRIX.md`](CLAIM_MATRIX.md). Some of that
prose exists only on branch `pre-cleanup-archive`.
