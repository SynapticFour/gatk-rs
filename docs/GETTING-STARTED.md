# Getting started

How to build `gatk-rs` and run one interval. Agreement with GATK is in [CLAIM_MATRIX.md](CLAIM_MATRIX.md), not on this page.

## Prerequisites

- Rust stable with rustfmt and clippy (`rust-toolchain.toml`). The workspace `rust-version` is 1.75.
- A C toolchain. BAM, VCF, and FASTA I/O uses `rust-htslib` / HTSlib, compiled as a crate dependency. CI builds on `ubuntu-latest` and `macos-latest`. Windows is not a CI target (`hts-sys` OpenSSL headers).
- [samtools](https://www.htslib.org/) to index the BAM and FASTA you pass in, and for parity scripts that stage fixtures. It is not the caller.
- Java GATK 4.4.0.0 only when you diff against the oracle. The `gatk-rs` binary does not embed it.

## Build and run

```bash
cargo build -p gatk-cli --release --locked
samtools faidx ref.fa
samtools index reads.bam
./target/release/gatk-rs HaplotypeCaller \
  -R ref.fa -I reads.bam -O out.vcf \
  -L chr20:10000000-10050000
```

`ref.fa` and `reads.bam` are paths you supply. They are not in a clean checkout. The interval `chr20:10000000-10050000` is an illustration, not the tracked 32 bp smoke fixture under `parity/fixtures/`. With those files indexed, the command writes one plain-text VCF. Paths ending in `.vcf.gz`, `.bgz`, `.bgzip`, or `.bcf` are rejected. It does not prove GATK equivalence for that interval.

Defaults, gVCF, Combine, genotype, and hard filters: [USER_GUIDE.md](USER_GUIDE.md).

## Inputs

| File | Requirement |
|------|-------------|
| Reference FASTA (`-R`) | Contig names must match the BAM. A sibling `.fai` is how contig lengths are loaded without scanning a whole genome. This loader does not read a Picard `.dict`. |
| Reads (`-I`) | BAM opened with an index (`.bai`). Repeating `-I` does not implement Java joint-read HaplotypeCaller. |
| Intervals (`-L`) | Repeatable `contig:start-end`. Omit `-L` only when you intend the whole reference. |

## Downloadable builds

No product release is published from this repository yet. A version tag triggers a validated Linux x86_64 GNU binary release (`x86_64-unknown-linux-gnu`). The release assets are that executable and `SHA256SUMS`. The tag must be `v` plus the workspace version. The uploaded binary's `gatk-rs --version` line must report that version and the tagged commit. macOS and Windows binaries are not offered.

The only git tag found on this clone, besides that rule, was `giab-ref-v1`, which peels to `60ddb683` and is not the default branch.

Docker image publication and GitHub Pages deployment are not part of this tag workflow.
