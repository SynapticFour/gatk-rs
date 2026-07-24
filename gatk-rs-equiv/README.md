# gatk-rs-equiv

Standalone CLI that checks **output equivalence** of **gatk-rs HaplotypeCaller** against **real GATK4 (Java)** using community-standard evaluation tools:

1. **Preferred:** [Illumina hap.py](https://github.com/Illumina/hap.py)
2. **Fallback:** [RTG Tools](https://github.com/RealTimeGenomics/rtg-tools) `vcfeval`

This crate is intentionaly **independent** of the internal L2–L14 parity history. You only need:

- a GATK4 JAR (or `gatk` launcher) at the pinned release,
- a built `gatk-rs` binary,
- reference / BAM / truth VCF / confident BED,
- hap.py **or** RTG on `PATH` (or use `Dockerfile.equiv`).

Authoritative product claims for the wider project still live in
[`docs/CLAIM_MATRIX.md`](../docs/CLAIM_MATRIX.md).
This tool measures **what you ask it to measure** on the inputs you provide.

---

## Limits (read this)

| In scope | Out of scope |
|----------|----------------|
| SNP/INDEL precision, recall, F1 vs a truth VCF (via hap.py / vcfeval) | Runtime / wall-clock parity |
| **Rust−Java F1 delta** as the gate metric | Memory / RSS parity |
| Stratified F1 when you pass stratification BEDs | Genome-wide claims without providing genome-wide inputs |
| Direct site table: exact POS+REF+ALT+GT match vs FORMAT drift | Full QUAL/INFO/annotation bitwise identity |
| Exit code usable as a CI gate | Replacing clinical validation |

**Not affiliated with the Broad Institute.** GATK is a trademark of the Broad; this is an independent community checker.

hap.py / vcfeval perform sophisticated haplotype-aware matching. This tool does **not** invent its own truth-matching heuristics for the F1 numbers — it shells out to those engines and parses their summaries.

---

## Build

From the workspace root:

```bash
cargo build -p gatk-rs-equiv --release
# binary: target/release/gatk-rs-equiv
```

Also build the caller under test:

```bash
cargo build -p gatk-cli --release
# binary: target/release/gatk-rs
```

Pinned GATK version: see [`GATK_PINNED_SHA`](../GATK_PINNED_SHA) and
[`docs/GATK_PINNED.env`](../docs/GATK_PINNED.env)
(`4.4.0.0` / SHA `2dbc0258…`).

---

## Usage

### `run`

```bash
gatk-rs-equiv run \
  --java-gatk-jar /path/to/gatk-package-4.4.0.0-local.jar \
  --rust-binary ./target/release/gatk-rs \
  --reference hs37d5.fa \
  --bam sample.bam \
  --truth-vcf HG001_benchmark.vcf.gz \
  --confident-regions HG001_benchmark.bed \
  --interval 20:10000000-10050000 \
  --out /tmp/equiv_out \
  --f1-delta-threshold 0.02 \
  --stratification-bed low_complexity=/path/to/LowComplexity.bed \
  --stratification-bed segdup=/path/to/SegmentalDuplications.bed
```

What it does:

1. Runs **Java** `HaplotypeCaller` and **gatk-rs** `HaplotypeCaller` on the same `-R/-I/-O/-L`.
2. Runs hap.py (or RTG vcfeval) for **each** callset against the same truth + confident BED.
3. Computes **ΔF1 = rust_f1 − java_f1** per class/stratum.
4. Builds a direct Java↔Rust site comparison (exact GT match vs FORMAT-only drift).
5. Writes `results.json`, `REPORT.md`, `report.json`.
6. Exits **0** only if every configured |ΔF1| ≤ `--f1-delta-threshold` (default `0.02`).

### `report`

Re-render Markdown/JSON from an existing results directory (optional threshold override):

```bash
gatk-rs-equiv report --results-dir /tmp/equiv_out --f1-delta-threshold 0.02
```

### `differential-fuzz`

Truth-set F1 deltas do not prove **identity** on pathological inputs. This command:

1. Builds synthetic FASTA/BAM scenarios (read length, coverage, indels, soft-clips, mate overlap, MAPQ, BQ).
2. Runs Java GATK4 and gatk-rs HaplotypeCaller with the same parameters.
3. On allele/GT/FORMAT divergence (AD/DP within `--format-ad-tol` allowed), **shrinks** the scenario.
4. Writes a minimal fixture under `gatk-haplotypecaller/tests/fixtures/regressions/`.
5. Optionally opens a GitHub issue with label `parity-divergence` (`--open-github-issue`, needs `gh`).

Preferred entrypoint (M4-safe defaults): [`fuzz/run_hc_differential.sh`](../fuzz/run_hc_differential.sh).
LibFuzzer generative surface: `cargo +nightly fuzz run hc_differential` (scenario decode only).

```bash
./fuzz/run_hc_differential.sh --iterations 8
# or:
gatk-rs-equiv differential-fuzz \
  --java-gatk-jar /path/to/gatk-package-4.4.0.0-local.jar \
  --rust-binary ./target/release/gatk-rs \
  --iterations 8 \
  --open-github-issue
```

Requires `samtools` on `PATH`. Exit `1` if any divergence was found (fixtures still written).

### Exit codes

| Code | Meaning |
|------|---------|
| `0` | Gate passed (all \|ΔF1\| ≤ threshold), or differential-fuzz found no divergences |
| `1` | Gate failed (some \|ΔF1\| above threshold), or differential-fuzz found ≥1 divergence |
| `2` | Tool / config / engine error |

---

## Docker (`Dockerfile.equiv`)

Builds an image with:

- GATK **4.4.0.0** (matches `GATK_PINNED.env`),
- hap.py + RTG (fallback),
- `gatk-rs` and `gatk-rs-equiv` release binaries.

```bash
docker build -f Dockerfile.equiv -t gatk-rs-equiv:local .
docker run --rm -v "$PWD/data:/data" gatk-rs-equiv:local \
  gatk-rs-equiv run \
    --java-gatk-bin gatk \
    --rust-binary gatk-rs \
    --reference /data/ref.fa \
    --bam /data/sample.bam \
    --truth-vcf /data/truth.vcf.gz \
    --confident-regions /data/truth.bed \
    --interval 20:10000000-10050000 \
    --out /data/equiv_out
```

---

## Engine selection

| `--engine` | Behavior |
|------------|----------|
| `auto` (default) | Prefer `hap.py` / `HAPPY_BIN`; else `rtg` / `RTG_BIN` |
| `happy` | Require hap.py |
| `vcfeval` | Require RTG `vcfeval` |

---

## Output layout

```
out/
  java.vcf
  rust.vcf
  eval/java.*          # hap.py / vcfeval artifacts
  eval/rust.*
  manifest.json
  results.json
  report.json
  REPORT.md
```

---

## Example GIAB-style slice

Any interval works; a common small dense window used elsewhere in this repo:

- BAM: NA12878 GIAB 30× downsample slice  
- Interval: `20:10000000-10050000`  
- Truth: HG001 GRCh37 v4.2.1 VCF + high-confidence BED  
- Optional strata: GIAB LowComplexity / SegmentalDuplications BEDs  

Provide your own paths — this tool does not download datasets.
