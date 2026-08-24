# gatk-rs

[![CI](https://github.com/SynapticFour/gatk-rs/actions/workflows/ci.yml/badge.svg)](https://github.com/SynapticFour/gatk-rs/actions/workflows/ci.yml)
[![License: Apache-2.0](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](LICENSE)

Built by **[Synaptic Four](https://synapticfour.com)**.

gatk-rs is an independent, community-driven reimplementation and is not
affiliated with, endorsed by, or supported by the Broad Institute.
"GATK" is a trademark of the Broad Institute; this project's name and
branding will be revisited if requested. Parity tests **call** a pinned
GATK 4.4 jar as an oracle; they do not ship Broad source (see [`NOTICE.md`](NOTICE.md)).

> **Maturity: Alpha** — validated on limited genomic regions and fixtures, not as a
> genome-wide clinical drop-in. Authoritative claims and non-claims:
> [`docs/CLAIM_MATRIX.md`](docs/CLAIM_MATRIX.md). Trademark and third-party data notes:
> [`NOTICE.md`](NOTICE.md).

### One screen: what this proves / does not prove

| Proves (on `main`) | Does **not** prove |
|--------------------|--------------------|
| Germline spine vs pinned **GATK 4.4**: HC → CombineGVCFs → GenotypeGVCFs → hard `VariantFiltration` | Full GATK4 toolkit (BQSR, VQSR, Mutect2, gCNV/SV, Funcotator, …) |
| **L2** synthetic gates (223/223) + **P12** L3/L4/L5 on `chr2:92300000–92350000` | Genome-wide / full-autosome HC equivalence |
| Synthetic joint-genotype cohort ladder (**≤100** samples on a tiny interval) | WGS × large-N / GenomicsDB-class joint calling |
| Scoped algorithm parity with honest waivers (W-H1 / W-H3 / W-L7-FORMAT) | Clinical drop-in, bitwise-identical QUAL/FORMAT everywhere, or a product launch |
| Equivalence **harness** green on GIAB **smoke** (hosted CI, RTG F1 Δ=0 on three ~50 kb windows) | A **signed** GIAB `ci-subset` / full-autosome F1 claim — still **unsigned** |

Authority: [`docs/CLAIM_MATRIX.md`](docs/CLAIM_MATRIX.md). Historical L6–L14 narratives live on `pre-cleanup-archive` only — not unqualified **Yes** rows here.

## Why does this exist?

Nobody asked for this. Here it is anyway: a Rust port of the GATK4
HaplotypeCaller.

I run a company built on the claim that AI tools can be used to real
customer benefit — not just for demos, but for genuinely hard, unglamorous
engineering work. At some point that claim needed a stress test. So I went
looking for a project ambitious enough to actually find the limits: where AI
tooling holds up, where it breaks, and how it behaves once a task stops
being a toy problem. Bioinformatics software — a mature, gnarly, decades-old
Java codebase with real scientific stakes — seemed like a fair fight. "Someone
should really port GATK4 to Rust" is a sentence I apparently said out loud
once, and this project is what happened next.

Every line of code in this repository was written by AI. I didn't write any
of it myself. What I *did* do, for months, was something closer to herding
than engineering — reviewing output, redirecting agents that were confidently
heading the wrong way, and reining things in when they got creative in ways
nobody asked for. If this project taught me one thing about my own job, it's
that my actual title should probably be **Agent Wrangler**. Most days, that's
exactly what it felt like.

So — is this now a scientifically proven, drop-in equivalent to GATK4? I
can't claim that with any authority, and I'd be suspicious of anyone who
could after a project like this. What I *can* say: on the **signed scopes**
in [`CLAIM_MATRIX`](docs/CLAIM_MATRIX.md) we have reproducible gate evidence —
and just as importantly, that file lists what we do **not** claim (including
unsigned GIAB `ci-subset` / full-autosome runs). That's not proof of
genome-wide equivalence. It is the most honest proof I can offer that
AI-assisted engineering can produce something real, not just something that
looks real until someone with domain expertise looks closely.

If you know genomics and something here is wrong, that's not a
disappointment — that's the actual experiment working. Open an issue.

### Live equivalence dashboard

**[Equivalence dashboard (GitHub Pages)](https://synapticfour.github.io/gatk-rs/)**
(source: [`docs/parity-site/`](docs/parity-site/)) — Chart.js view of hap.py
metrics **only after** a signed publish lands in `history.json`.
Until then the UI shows “No published runs yet.” Treat the site as
**instrumentation**, not a product claim — GIAB `ci-subset` is still
**unsigned** in [`CLAIM_MATRIX`](docs/CLAIM_MATRIX.md).

## Validated Scope

gatk-rs targets the **germline short-variant** workflow against pinned GATK
**4.4** — a focused experiment, not a toolkit clone:

**HaplotypeCaller → CombineGVCFs → GenotypeGVCFs → VariantFiltration** (hard filters).

Claims and non-claims for what that path has actually proven:
[`docs/CLAIM_MATRIX.md`](docs/CLAIM_MATRIX.md). Structure:
[`docs/ARCHITECTURE.md`](docs/ARCHITECTURE.md).

**External pilots:** run the same spine on your own BAMs and compare Java↔Rust
yourself — [`docs/PILOT_GUIDE.md`](docs/PILOT_GUIDE.md) +
[`scripts/pilot/compare_callsets.py`](scripts/pilot/compare_callsets.py).

**Architecture decision (read this first for scope):**
[**ADR 0001 — Scope boundary**](docs/adr/0001-scope-boundary.md)
(why BQSR, VQSR, Mutect2, gCNV/SV, and Funcotator are intentionally out of
scope). Related: [ADR 0002](docs/adr/0002-remove-gatk-tools.md) (no generic
`gatk-tools` crate — use **samtools** / **bcftools** for sort/index/view).

## What this is

A native Rust workspace focused on **HaplotypeCaller** plus the post-call
helpers above. Generic BAM/VCF utilities belong in **samtools** / **bcftools**
([ADR 0002](docs/adr/0002-remove-gatk-tools.md)); leftover utility subcommands
stay callable for parity harnesses but are hidden from default `--help`.

### VariantFiltration vs VQSR

`gatk-rs variant-filtration` implements GATK-compatible **hard-filtering**
(`--filter-expression` / `--filter-name`, plus `--preset snp|indel` for the
official Best Practices tables). **This does not replace VQSR algorithmically** —
it is the recommended pragmatic fallback for smaller cohorts where VQSR is not
cleanly trainable. That matches official GATK guidance (VQSR is recommended only
once the cohort is large enough to support model training). See
[`docs/ARCHITECTURE.md`](docs/ARCHITECTURE.md) and [`docs/CLAIM_MATRIX.md`](docs/CLAIM_MATRIX.md).

## Build and run

```bash
cargo build -p gatk-cli --release
# `--threads N` / `--nt N` / `-t N` sizes the Rayon pool (default: CPU count).
# Active regions are processed in parallel inside one process; VCF rows are
# sorted/deduped before write (byte-identical across thread counts).
./target/release/gatk-rs HaplotypeCaller -R ref.fa -I reads.bam -O out.vcf \
  -L chr2:92300000-92350000 --threads 8
# Official SNP hard filters (GATK Best Practices table):
./target/release/gatk-rs variant-filtration -V snps.vcf -O snps.filtered.vcf --preset snp
# Or explicit Java-compatible pairs:
./target/release/gatk-rs variant-filtration -V snps.vcf -O snps.filtered.vcf \
  --filter-expression "QD < 2.0" --filter-name QD2 \
  --filter-expression "FS > 60.0" --filter-name FS60
./target/release/gatk-rs --help
```

See [`CONTRIBUTING.md`](CONTRIBUTING.md) for tests and contribution workflow.

## Performance (fair comparison)

End-to-end HaplotypeCaller timings are published only from the dedicated quiet
host ([`docs/ci/PERF_BENCHMARK_HOST.md`](docs/ci/PERF_BENCHMARK_HOST.md)), never
from GitHub-hosted runners or the genomewide correctness VM.

| Contract | Detail |
|----------|--------|
| Workflow | [`.github/workflows/benchmark.yml`](.github/workflows/benchmark.yml) |
| Harness | [`scripts/perf/run_fair_hc_comparison.sh`](scripts/perf/run_fair_hc_comparison.sh) |
| Configs | gatk-rs `LOGLESS_HMM` (scalar), gatk-rs `AVX`/SIMD, Java **`FASTEST_AVAILABLE`**, Java `LOGLESS_CACHING` |
| Stats | ≥5 repeats → **median ± sample stdev** (not best-of-N) |
| Regions | small / medium / large nested GIAB windows |
| Metrics | wall, user, sys, Peak-RSS; RAPL energy when `perf` allows |
| Primary baseline | Java **`FASTEST_AVAILABLE`** (native AVX verified before timing) |
| Raw report | [`docs/perf/FAIR_HC_COMPARISON.md`](docs/perf/FAIR_HC_COMPARISON.md) · JSON · host [`HOST_SPECS.md`](docs/perf/HOST_SPECS.md) |
| Dashboard | Performance tab on the [public site](docs/parity-site/) (`data/perf_history.json`) |

**Production PairHMM default remains `LOG10_PAIRHMM`** until a signed GIAB run
clears SIMD/Logless. SIMD code path: `--pair-hmm AVX` (unit gate
`pairhmm_simd_vs_scalar_test`).

Headline HC speedups appear in `FAIR_HC_COMPARISON.md` after the first successful
dedicated-host run. Until then, treat any laptop Criterion numbers (e.g. local
NEON microbench in [`docs/perf/PAIRHMM_SPEEDUP.md`](docs/perf/PAIRHMM_SPEEDUP.md))
as **dev-host only** — not interchangeable with the fair suite.

## Memory profile (Peak-RSS)

Two labeled Peak-RSS profiles via
[`scripts/perf/run_hc_memory_profile.sh`](scripts/perf/run_hc_memory_profile.sh)
vs pinned Java GATK **4.4.0.0** (`-Xms1g -Xmx4g`). Full tables, commands, and
raw logs: [`docs/perf/HC_MEMORY_PROFILE.md`](docs/perf/HC_MEMORY_PROFILE.md).

### A. Trivial smoke — reproducibility only

Checked-in fixture `parity/fixtures/`, interval `chr1:1-32` (32 bp). Dominated
by JVM/runtime fixed cost — **not** a public “X% less memory” claim.

| Engine | Peak RSS (run `20260724T181512Z`) |
|--------|-------------------------------------|
| gatk-rs (release) | **9.52 MiB (9744 KiB)** |
| Java GATK 4.4.0.0 | **437.49 MiB (447988 KiB)** |

### B. Realistic GIAB-dense window — public-claim basis

Multi-Mb NA12878 window on the known-dense chr20 locus
(default `20:10000000-12000000`, 2 Mb). **Only this profile** may back a public
memory claim, and **only** when measured on the dedicated
`gatk-rs-benchmark` host ([`docs/ci/PERF_BENCHMARK_HOST.md`](docs/ci/PERF_BENCHMARK_HOST.md))
with [`docs/perf/HOST_SPECS.md`](docs/perf/HOST_SPECS.md) populated.

| Engine | Peak RSS |
|--------|----------|
| gatk-rs / Java GATK 4.4 | *Pending dedicated-host run* — see `HC_MEMORY_PROFILE.md` |

Until that host run lands, do **not** advertise a genome-wide memory savings %.

## Equivalence

| Path | Purpose |
|------|---------|
| [`gatk-rs-equiv/`](gatk-rs-equiv/) | GIAB / hap.py / vcfeval + differential fuzz |
| [`scripts/parity/`](scripts/parity/) | L2 / P12 / GIAB harness scripts |
| [`tools/equivalence/README.md`](tools/equivalence/README.md) | Index of the proof surface |
| [`docs/CLAIM_MATRIX.md`](docs/CLAIM_MATRIX.md) | What those tests do and do not claim |

## Changelog

See [`CHANGELOG.md`](CHANGELOG.md).

## License

Apache License 2.0 — see [`LICENSE`](LICENSE). Trademark, GIAB, and oracle-jar notes: [`NOTICE.md`](NOTICE.md).

## Security and conduct

- Vulnerability reporting: [`SECURITY.md`](SECURITY.md)
- Code of conduct: [`CODE_OF_CONDUCT.md`](CODE_OF_CONDUCT.md)
