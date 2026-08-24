# NOTICE

## Independent community project (not Broad GATK)

gatk-rs (also referred to as GATK-RS in historical docs) is an **independent,
community-driven reimplementation**. It is **not** affiliated with, endorsed by,
or supported by the Broad Institute or the official GATK project.

“GATK” is a trademark of the Broad Institute. This project’s name, CLI branding,
and flag naming will be revisited if requested by the trademark holder.

CLI familiarity flags such as `--java-options` and GATK-style short options
(`-R`, `-I`, `-O`, `-L`) exist for **user interoperability** with existing
pipelines. The `gatk-rs` binary is native Rust and does **not** launch the Broad
GATK JVM. Presence of those flags is **not** an endorsement or affiliation claim.

## License of this repository

This software is licensed under the **Apache License, Version 2.0**. See
[`LICENSE`](LICENSE).

Apache-2.0 is compatible with typical community redistribution of *original*
Rust source in this tree. It does **not** by itself clear third-party data or
upstream Java code that you download separately for parity testing.

## What this tree does *not* ship as Broad GATK source

This repository aims to reimplement algorithms against a **pinned GATK 4.4**
oracle for differential testing. Contributors must **not** copy Broad GATK /
HTSJDK Java source into this tree. Observing Java behavior (running the
published GATK jar, reading public docs) for parity is distinct from
redistributing Broad code.

If any file were ever found to contain Broad-copyrighted source verbatim,
it should be removed and replaced with an independent implementation.

## Third-party data and tools (separate terms)

Downstream users and CI scripts may download external assets. Those assets are
**not** relicensed by this NOTICE or by Apache-2.0 covering gatk-rs source.

| Asset / tool | Typical use here | Notes |
|--------------|------------------|--------|
| **GIAB / NIST Genome in a Bottle** truth VCFs & high-confidence BEDs | Scientific F1 / equivalence (`gatk-rs-equiv`, GIAB harnesses) | Subject to NIST/GIAB distribution and citation terms. Download from official GIAB channels; do not assume Apache-2.0 covers the truth sets. |
| **Pinned GATK 4.4 jar** (Broad) | Java oracle for differential parity | Licensed under Broad/GATK terms; obtained by the user or CI, not redistributed as part of this source tree’s Apache grant. |
| **hap.py / RTG Tools** | Optional equivalence engines | Their own licenses apply when you install them. |
| **htslib / rust-htslib** and other Rust crates | BAM/VCF/FASTA I/O | Covered by each crate’s license (see `Cargo.lock` / `cargo deny`). |
| **Fixture BAMs / FASTAs under `parity/fixtures/`** | Tiny synthetic or derived test inputs | Intended for automated tests; if a fixture embeds restricted human data, treat upstream terms as controlling and prefer synthetic replacements. |

## Parity Java in this tree

Oracle dump programs under `scripts/parity/java/` are **original Synaptic Four
code**. They compile **against** a user- or CI-supplied pinned GATK 4.4 jar
(`import org.broadinstitute.hellbender…` is interoperability with the published
oracle). They do **not** live in the `org.broadinstitute.*` package namespace
and they are **not** Broad/Hellbender source.

Do not add files under `org/broadinstitute/` in this repository. If a future
parity dump needs a package-private GATK API, prefer a documented reflection
shim or a public GATK entry point — do not re-home Synaptic Four sources into
Broad packages.

## Citation / marketing

Do not describe gatk-rs as “official GATK”, “Broad GATK for Rust”, or a
genome-wide clinical drop-in unless [`docs/CLAIM_MATRIX.md`](docs/CLAIM_MATRIX.md)
explicitly asserts that scope (it currently does not).

See also the disclaimer at the top of [`README.md`](README.md) and CLI `--help`.
