# gatk-rs

[![CI](https://github.com/SynapticFour/gatk-rs/actions/workflows/ci.yml/badge.svg)](https://github.com/SynapticFour/gatk-rs/actions/workflows/ci.yml)
[![License: Apache-2.0](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](LICENSE)

**Version 0.1.0 (Alpha).** Independent research reimplementation of a narrow
GATK 4 germline path. Not affiliated with, endorsed by, or supported by the
Broad Institute. “GATK” is a trademark of the Broad Institute. See
[`NOTICE.md`](NOTICE.md).

Maintained by [Synaptic Four](https://synapticfour.com)
([contact@synapticfour.com](mailto:contact@synapticfour.com)).

## Why does this exist

`gatk-rs` is a Rust workspace and `gatk-rs` command-line tool for germline
short variants:

**HaplotypeCaller → CombineGVCFs → GenotypeGVCFs → VariantFiltration** (hard filters).

The comparison baseline is pinned GATK **4.4.0.0**
([`docs/GATK_PINNED.env`](docs/GATK_PINNED.env)). The goal is algorithm parity
on agreed gates, not a second copy of the GATK source tree and not a clinical
drop-in.

## Who it is for

Someone who already runs GATK HaplotypeCaller on a defined interval and wants
a native Rust caller they can build, diff against Java GATK 4.4, and extend
without treating an unscoped “equivalent to GATK” claim as true.

It is not a replacement for BQSR, VQSR, Mutect2, or a WGS joint-calling
service. Scope: [`docs/adr/0001-scope-boundary.md`](docs/adr/0001-scope-boundary.md).

## What works today

The default `HaplotypeCaller` command runs the production
`assembly-region-v1` path (`CallRegionArgs::strict_java()`): activity
regions, local reassembly, PairHMM, genotyping, and VCF or gVCF emission.
`CombineGVCFs`, `GenotypeGVCFs`, and `VariantFiltration` are real subcommands.

Product claims are only the **Yes** rows in
[`docs/CLAIM_MATRIX.md`](docs/CLAIM_MATRIX.md). In short, those rows cover:

- the canonical mid-B HaplotypeCaller path on one ActiveFull region,
  `2:92317000-92319000` (three oracle SNPs);
- P12 L3/L4/L5 on `chr2:92300000-92350000` (66 Java sites);
- L2 synthetic gates, recorded **223/223 strict** on **2026-07-24** (that date
  is the matrix evidence; this README does not re-certify a later commit);
- a synthetic joint-genotype ladder, recommended **≤ 100** samples on a 10 kb
  interval;
- hard-filter FILTER parity on a synthetic SNP boundary set;
- an HG001 GIAB **ci-subset** F1 comparison (not full autosomes).

A CI badge is a workflow link, not those claims.

## What is not GATK equivalence

| In this tree | Not in this tree |
|--------------|------------------|
| Germline HC → CombineGVCFs → GenotypeGVCFs → hard `VariantFiltration` | BQSR, VQSR, Mutect2, gCNV/SV, Funcotator |
| Scoped gates above | Genome-wide or full-autosome HaplotypeCaller equivalence |
| One-sample HC, then Combine | Java-style joint HC that merges multiple `-I` BAMs into one read set |
| Hard filters (`--preset snp` / `indel`) | VQSR |
| Omitted `--pair-hmm` selects `FASTEST_AVAILABLE` (host SIMD when the build can use it, otherwise Logless scalar). Log10 scalar is an explicit `--pair-hmm LOG10_PAIRHMM` | A signed claim that the omitted flag matches GATK’s log10 PairHMM, or a signed SIMD speed claim |

Also not claimed: bitwise-identical QUAL/FORMAT on every site, bamout,
DRAGSTR, DRAGEN mode, allele-specific `AS_*` annotations, or
`--assembly-region-out`. Multi-sample Combine loads gVCFs in memory and has
no GenomicsDB path.

Waivers that limit an otherwise **Yes** row (W-H1, W-H3, W-L7-FORMAT, and the
W-J* tables) are named in the claim matrix. Do not drop them when citing a
gate.

## Compared with GATK HaplotypeCaller

On the signed scopes, the caller is exercised against GATK 4.4.0.0
(source SHA `2dbc025821bc5f686c423ff332a41e6cef892a77`, image
`us.gcr.io/broad-gatk/gatk:4.4.0.0`). Outside those scopes, agreement is
unproven. The CLI accepts familiar `-R` / `-I` / `-O` / `-L` flags.
`--java-options` is accepted and does not start a JVM.

Whole-codebase HaplotypeCaller parity is **not** established. The 6R/7R
notes under [`docs/parity/`](docs/parity/) are investigation records, not
extra Yes rows.

## Prerequisites

- Rust **stable** with rustfmt and clippy (`rust-toolchain.toml`). The
  workspace `rust-version` is 1.75.
- A C toolchain. BAM/VCF/FASTA I/O uses `rust-htslib` / HTSlib, compiled as
  a crate dependency. CI builds on `ubuntu-latest` and `macos-latest`.
  Windows is not a CI target (`hts-sys` OpenSSL headers).
- [samtools](https://www.htslib.org/) to index the BAM and FASTA you pass
  in, and for parity scripts that stage fixtures. It is not the caller.
- Java GATK 4.4.0.0 only when you diff against the oracle. The `gatk-rs`
  binary does not embed it.

## Downloadable builds

No product release is published from this repository yet. The tag workflow
that publishes one is limited to a single native binary, Linux x86_64 GNU
(`x86_64-unknown-linux-gnu`), plus a `SHA256SUMS` file for that executable.
The tag must be `v` plus the workspace version. The uploaded binary's
`gatk-rs --version` line must report that version and the tagged commit.
macOS and Windows binaries are not offered.

This project is an alpha research reimplementation. It is not a drop-in
replacement for Broad Institute GATK and it is not clinically validated.

A version tag can still start the existing Docker publish job. That job is
unchanged. Its multi-architecture image is not a tested product release.

## Build and run

```bash
cargo build -p gatk-cli --release --locked
samtools faidx ref.fa
samtools index reads.bam
./target/release/gatk-rs HaplotypeCaller \
  -R ref.fa -I reads.bam -O out.vcf \
  -L chr20:10000000-10050000
```

`ref.fa` and `reads.bam` are paths you supply. They are not in a clean
checkout. The interval `chr20:10000000-10050000` is an illustration, not the
tracked 32 bp smoke fixture under `parity/fixtures/`. This command has not
been run as part of the documentation. With those files indexed, it writes
one plain-text VCF. Paths ending in `.vcf.gz`, `.bgz`, `.bgzip`, or `.bcf`
are rejected. It does not prove GATK equivalence for that interval. Defaults,
gVCF, Combine, genotype, and hard filters:
[`docs/USER_GUIDE.md`](docs/USER_GUIDE.md).

## Validated Scope

Signed product claims are only the **Yes** rows in
[`docs/CLAIM_MATRIX.md`](docs/CLAIM_MATRIX.md). Those rows cover the canonical
mid-B HaplotypeCaller path on `2:92317000-92319000`, P12 L3/L4/L5 on
`chr2:92300000-92350000`, the L2 synthetic record of 223/223 strict on
2026-07-24, the synthetic joint-genotype ladder at the recommended cap of
100 samples on a 10 kb interval, hard-filter FILTER parity on the synthetic
SNP boundary set, and the HG001 GIAB ci-subset F1 comparison. Waivers named
on those rows stay part of the claim.

That record is not genome-wide or full-autosome HaplotypeCaller equivalence,
and it does not include BQSR, VQSR, Mutect2, gCNV/SV, or Funcotator. The
boundary is [`docs/adr/0001-scope-boundary.md`](docs/adr/0001-scope-boundary.md).

## License

This repository is licensed under the Apache License, Version 2.0
([`LICENSE`](LICENSE)). The workspace package metadata records the same
identifier, `Apache-2.0`. [`NOTICE.md`](NOTICE.md) states that this is an
independent research reimplementation, that “GATK” is a trademark of the
Broad Institute, and that GIAB truth sets, the pinned GATK jar, hap.py, and
the HTSlib crates keep their own terms. Apache-2.0 on this tree does not
relicense those assets.

## Inputs

| File | Requirement |
|------|-------------|
| Reference FASTA (`-R`) | Contig names must match the BAM. A sibling `.fai` is how contig lengths are loaded without scanning a whole genome. This loader does not read a Picard `.dict`. |
| Reads (`-I`) | BAM opened with an index (`.bai`). Repeating `-I` does **not** implement Java joint-read HaplotypeCaller. |
| Intervals (`-L`) | Repeatable `contig:start-end`. Omit `-L` only when you intend the whole reference. |

## Where to go next

| Question | Document |
|----------|----------|
| Install, flags, outputs, failures | [`docs/USER_GUIDE.md`](docs/USER_GUIDE.md) |
| Diff your own Java and Rust VCFs | [`docs/PILOT_GUIDE.md`](docs/PILOT_GUIDE.md) |
| Pipeline stages and crates | [`docs/ARCHITECTURE.md`](docs/ARCHITECTURE.md) |
| What a gate does and does not prove | [`docs/PARITY.md`](docs/PARITY.md) |
| Yes / No claims | [`docs/CLAIM_MATRIX.md`](docs/CLAIM_MATRIX.md) |
| Doc map | [`docs/README.md`](docs/README.md) |
| Change the code | [`CONTRIBUTING.md`](CONTRIBUTING.md) |
| Performance numbers | Not published. Dedicated-host contract: [`docs/perf/FAIR_HC_COMPARISON.md`](docs/perf/FAIR_HC_COMPARISON.md). Do not quote laptop Criterion or the 32 bp smoke RSS profile. |
| Dashboards | [`docs/EQUIVALENCE_DASHBOARD.md`](docs/EQUIVALENCE_DASHBOARD.md) is pending: its F1 table is empty, so that file is not a result. [`docs/parity-site/data/history.json`](docs/parity-site/data/history.json) does contain historical runs. Those rows are not the claim matrix’s signed HG001 ΔF1 record. Claims: [`docs/CLAIM_MATRIX.md`](docs/CLAIM_MATRIX.md). |
| Security, license | [`SECURITY.md`](SECURITY.md), [`LICENSE`](LICENSE), [`NOTICE.md`](NOTICE.md) |

Changelog: [`CHANGELOG.md`](CHANGELOG.md).
