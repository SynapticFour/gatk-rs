# gatk-rs

gatk-rs is a Rust caller for germline short variants on an interval you already sequence: HaplotypeCaller, then CombineGVCFs, GenotypeGVCFs, and hard filters.

For someone who already runs GATK HaplotypeCaller on a defined interval and wants a Rust build they can diff against Java GATK 4.4.

**Status: No public status is recorded.** The claims registry has no gatk-rs stage. The workspace version in `Cargo.toml` is 0.1.0. That number is not a status.

[![CI](https://github.com/SynapticFour/gatk-rs/actions/workflows/ci.yml/badge.svg)](https://github.com/SynapticFour/gatk-rs/actions/workflows/ci.yml)
[![License: Apache-2.0](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](LICENSE)

Not affiliated with, endorsed by, or supported by the Broad Institute. "GATK" is a trademark of the Broad Institute. See [NOTICE.md](NOTICE.md).

Maintained by [Synaptic Four](https://synapticfour.com) ([contact@synapticfour.com](mailto:contact@synapticfour.com)).

## Why does this exist

The comparison baseline is pinned GATK 4.4.0.0 ([docs/GATK_PINNED.env](docs/GATK_PINNED.env)). The goal is algorithm parity on agreed gates, not a second copy of the GATK source tree and not a clinical drop-in.

It is not a replacement for BQSR, VQSR, Mutect2, or a whole-genome joint-calling service. Scope: [docs/adr/0001-scope-boundary.md](docs/adr/0001-scope-boundary.md).

## What it does

- The default HaplotypeCaller command runs the `assembly-region-v1` path (`CallRegionArgs::strict_java()`).
- `CombineGVCFs`, `GenotypeGVCFs`, and `VariantFiltration` are real subcommands.
- Product claims are only the **Yes** rows in [docs/CLAIM_MATRIX.md](docs/CLAIM_MATRIX.md).

## Limits

- The source names that default path a production default (`gatk-haplotypecaller/src/run.rs`). The name is the code path. This repository is not clinically validated and it is not a drop-in replacement for Broad Institute GATK.
- Whole-codebase HaplotypeCaller parity is not established. BQSR, VQSR, Mutect2, gCNV, SV, and Funcotator are not in this tree.
- HG001 GIAB ci-subset F1 at `f70abf86` is a historical record pending revalidation. It is not a current-main **Yes** row. Full-autosome GIAB is not complete.
- A CI badge is a workflow link, not those claims.
- Performance numbers are not published here. Do not quote a laptop Criterion run or the 32 bp smoke memory profile. The dedicated-host contract is [docs/perf/FAIR_HC_COMPARISON.md](docs/perf/FAIR_HC_COMPARISON.md).

The full gate table is in [docs/FOR-EVALUATORS.md](docs/FOR-EVALUATORS.md).

## Quick start

Prerequisites: Rust, a C toolchain, and samtools. The workspace `rust-version` is 1.75. Details are in [docs/GETTING-STARTED.md](docs/GETTING-STARTED.md).

`ref.fa` and `reads.bam` are yours. They are not in a clean checkout.

```bash
cargo build -p gatk-cli --release --locked
samtools faidx ref.fa
samtools index reads.bam
./target/release/gatk-rs HaplotypeCaller \
  -R ref.fa -I reads.bam -O out.vcf \
  -L chr20:10000000-10050000
```

With those files indexed, this writes one plain-text VCF. It does not prove GATK equivalence for that interval.

## Build and run

The commands and the input table are in [docs/GETTING-STARTED.md](docs/GETTING-STARTED.md). The same build is the block under Quick start. Windows is not a CI target.

## Validated Scope

Signed product claims are only the **Yes** rows in [docs/CLAIM_MATRIX.md](docs/CLAIM_MATRIX.md). Those rows cover the canonical mid-B HaplotypeCaller path on `2:92317000-92319000`, P12 L3/L4/L5 on `chr2:92300000-92350000`, the L2 synthetic record of 223/223 strict on 2026-07-24, the synthetic joint-genotype ladder at the recommended cap of 100 samples on a 10 kb interval, and hard-filter FILTER parity on the synthetic SNP boundary set. Waivers named on those rows stay part of the claim. The boundary is [docs/adr/0001-scope-boundary.md](docs/adr/0001-scope-boundary.md).

## Documentation

Every file under `docs/` is listed in [docs/README.md](docs/README.md). The groups are a tutorial, how-to guides, reference, and explanation. Notes under `docs/parity/` and `docs/perf/` are investigation records. They are not extra **Yes** rows.

### Tutorial

- [Getting started](docs/GETTING-STARTED.md)
- [User guide](docs/USER_GUIDE.md)

### How-to guides

- [Pilot guide](docs/PILOT_GUIDE.md) — diff your own Java and Rust VCFs
- [Release checklist](docs/RELEASE_CHECKLIST.md)

### Reference

- [Claim matrix](docs/CLAIM_MATRIX.md)
- [For evaluators](docs/FOR-EVALUATORS.md)
- [Parity](docs/PARITY.md)
- [Supported flags](docs/SUPPORTED_FEATURE_MATRIX.md)

### Explanation

- [Architecture](docs/ARCHITECTURE.md)
- [Scope boundary](docs/adr/0001-scope-boundary.md)

## License

This repository is licensed under the Apache License, Version 2.0 ([LICENSE](LICENSE)). The workspace package metadata records the same identifier, `Apache-2.0`. [NOTICE.md](NOTICE.md) states that this is an independent research reimplementation, that "GATK" is a trademark of the Broad Institute, and that GIAB truth sets, the pinned GATK jar, hap.py, and the HTSlib crates keep their own terms. Apache-2.0 on this tree does not relicense those assets. There is no gatk-rs licence claim in the registry.

## Cite

Use [CITATION.cff](CITATION.cff). Prefer that file over a hand-copied line.

## Website

- Open tools: [https://synapticfour.com/en/offene-werkzeuge/](https://synapticfour.com/en/offene-werkzeuge/)
- Technical profile: [https://synapticfour.com/en/technical/code/](https://synapticfour.com/en/technical/code/)

The same pages exist in German (`/de/offene-werkzeuge/`, `/de/technik/code/`) and in French (`/fr/offene-werkzeuge/`, `/fr/technique/code/`).

The map of repositories: [organisation profile](https://github.com/SynapticFour/.github/blob/main/profile/README.md).

Contact: [contact@synapticfour.com](mailto:contact@synapticfour.com)

Changelog: [CHANGELOG.md](CHANGELOG.md).
