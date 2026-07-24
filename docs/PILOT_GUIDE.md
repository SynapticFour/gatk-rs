# Pilot guide — verify gatk-rs on your own data

This guide is for an **external team** that already runs a Java GATK germline
pipeline and wants to check gatk-rs **on their own BAMs**, instead of trusting
the public dashboard alone.

The trust lever is: **same inputs → parallel Java + Rust → your comparison**.

Authoritative claims / non-claims: [`CLAIM_MATRIX.md`](CLAIM_MATRIX.md).  
Pinned Java GATK: [`GATK_PINNED.env`](GATK_PINNED.env) (**4.4.0.0**).

> gatk-rs is an independent community project, **not** affiliated with the Broad
> Institute. Maturity is **Alpha** — start with a small interval, not a clinical
> cutover.

---

## 0. What you will run

Validated spine (hard filters, not VQSR):

```text
HaplotypeCaller (GVCF)
  → CombineGVCFs
  → GenotypeGVCFs
  → VariantFiltration
```

For a first pilot, restrict to **one interval** (e.g. 50–500 kb) and
**≤ a few samples**. Large WGS cohorts: see Combine limits in
[`CLAIM_MATRIX.md`](CLAIM_MATRIX.md) (in-memory Combine; recommended ≤100
samples on the synthetic scale gate; no GenomicsDBImport).

---

## 1. Build gatk-rs

Requirements: Rust stable (edition as in repo root `rust-toolchain` / `Cargo.toml`),
usual C toolchain. From a clone:

```bash
git clone https://github.com/SynapticFour/gatk-rs.git
cd gatk-rs
cargo build -p gatk-cli --release --locked
# binary:
./target/release/gatk-rs --version
```

Optional (truth-based HC equivalence CLI used by our suite — not required for
the pilot script below):

```bash
cargo build -p gatk-rs-equiv --release --locked
```

---

## 2. Match your Java pin

Use **GATK 4.4.0.0** for apples-to-apples comparisons (Docker pin below matches
CI). Other 4.x lines will diverge for reasons unrelated to gatk-rs.

```bash
# Example: Broad public image
docker pull us.gcr.io/broad-gatk/gatk:4.4.0.0
```

Record the exact Java command lines you use today (ERC mode, stand-call-conf,
intervals, filtration expressions). **Copy those settings** into the Rust
commands — do not “improve” one side.

---

## 3. Parallel pipeline (your BAMs)

Assume:

| Variable | Meaning |
|----------|---------|
| `REF` | Indexed FASTA (+ `.fai`, `.dict`) |
| `BAM` / `BAMs` | Indexed BAM(s) you already trust in Java |
| `INTERVAL` | e.g. `chr20:10000000-10050000` or a BED via `-L` |
| `OUT` | Fresh work directory |

### 3a. Per-sample HaplotypeCaller → gVCF

**Java (your launcher / Docker):**

```bash
gatk --java-options "-Xmx4g" HaplotypeCaller \
  -R "$REF" -I sample.bam -O "$OUT/java/${SAMPLE}.g.vcf.gz" \
  -L "$INTERVAL" \
  -ERC GVCF \
  --standard-min-confidence-threshold-for-calling 30
```

**gatk-rs:**

```bash
./target/release/gatk-rs HaplotypeCaller \
  -R "$REF" -I sample.bam -O "$OUT/rust/${SAMPLE}.g.vcf" \
  -L "$INTERVAL" \
  --emit-ref-confidence GVCF \
  --stand-call-conf 30 \
  --threads 8
```

Notes:

- Index BAMs (`samtools index`) before either caller.
- Default PairHMM in gatk-rs is `LOG10_PAIRHMM`. For fair speed experiments you
  may try `--pair-hmm AVX`, but keep both sides aligned for an equivalence pilot.
- Multi-sample HC by passing multiple `-I` is **not** claimed as Java joint-read
  merge — call **one BAM → one gVCF** per sample, then Combine (see CLAIM_MATRIX).

### 3b. CombineGVCFs

```bash
# Java
gatk CombineGVCFs -R "$REF" \
  -V s1.g.vcf.gz -V s2.g.vcf.gz \
  -O "$OUT/java/combined.g.vcf.gz"

# Rust
./target/release/gatk-rs CombineGVCFs -R "$REF" \
  -V s1.g.vcf -V s2.g.vcf \
  -O "$OUT/rust/combined.g.vcf"
```

### 3c. GenotypeGVCFs

```bash
# Java
gatk GenotypeGVCFs -R "$REF" \
  -V "$OUT/java/combined.g.vcf.gz" \
  -O "$OUT/java/genotyped.vcf.gz" \
  --standard-min-confidence-threshold-for-calling 30

# Rust
./target/release/gatk-rs GenotypeGVCFs -R "$REF" \
  -V "$OUT/rust/combined.g.vcf" \
  -O "$OUT/rust/genotyped.vcf" \
  --stand-call-conf 30
```

### 3d. VariantFiltration (hard filters)

Prefer the same expressions you use in production, or the official SNP/INDEL
presets:

```bash
# Rust — official SNP hard-filter table
./target/release/gatk-rs VariantFiltration \
  -V "$OUT/rust/genotyped.vcf" \
  -O "$OUT/rust/filtered.vcf" \
  --preset snp

# Java — same expressions / --filter-name pairs you already use
gatk VariantFiltration \
  -V "$OUT/java/genotyped.vcf.gz" \
  -O "$OUT/java/filtered.vcf.gz" \
  --filter-expression "QD < 2.0" --filter-name QD2 \
  # … mirror your pipeline …
```

`variant-filtration` is **not** VQSR. If your production path is VQSR-only,
compare at the **GenotypeGVCFs** outputs instead, or add a hard-filter arm on
both sides for the pilot.

---

## 4. Compare (standalone tool)

Use [`scripts/pilot/compare_callsets.py`](../scripts/pilot/compare_callsets.py).
It does **not** require our CI runners, Docker helpers, or GIAB staging.

### 4a. No truth set (pure Rust-vs-Java)

```bash
python3 scripts/pilot/compare_callsets.py \
  --java "$OUT/java/genotyped.vcf.gz" \
  --rust "$OUT/rust/genotyped.vcf" \
  --out "$OUT/compare_gt"
```

Hard failures: missing sites, REF/ALT set disagreement, GT disagreement
(allele-identity), FILTER disagreement.  
Soft (reported, non-fatal by default): modest QUAL drift; AD/DP/GQ/PL within
tolerances aligned with waiver **W-L7-FORMAT**.

### 4b. With your truth (hap.py preferred, RTG vcfeval fallback)

Install [hap.py](https://github.com/Illumina/hap.py) or
[RTG Tools](https://github.com/RealTimeGenomics/rtg-tools) on `PATH`
(or set `HAPPY_BIN` / `RTG_BIN`).

```bash
python3 scripts/pilot/compare_callsets.py \
  --java "$OUT/java/genotyped.vcf.gz" \
  --rust "$OUT/rust/genotyped.vcf" \
  --reference "$REF" \
  --truth /path/to/truth.vcf.gz \
  --confident /path/to/confident.bed \
  --engine auto \
  --f1-delta-threshold 0.02 \
  --out "$OUT/compare_truth"
```

This scores **each** callset against **your** truth, then reports
**ΔF1 = Rust − Java** (same metric family as the internal suite). Exit `1` if
hard site diffs appear or `|ΔF1|` exceeds the threshold.

Artifacts: `$OUT/compare_*/REPORT.md` and `summary.json`.

### 4c. Optional: full HC runner from this repo

If you prefer one binary that also *launches* both callers against a truth set:

```bash
./target/release/gatk-rs-equiv run \
  --java-gatk-jar /path/to/gatk-package-4.4.0.0-local.jar \
  --rust-binary ./target/release/gatk-rs \
  --reference "$REF" --bam sample.bam \
  --truth-vcf truth.vcf.gz --confident-regions truth.bed \
  --interval "$INTERVAL" --out "$OUT/equiv"
```

See [`gatk-rs-equiv/README.md`](../gatk-rs-equiv/README.md). For “we already
have both VCFs”, prefer `compare_callsets.py`.

---

## 5. Expected deviations vs real bugs

### Usually **benign** (do not file as equivalence bugs by default)

| Observation | Why |
|-------------|-----|
| Small AD / DP / GQ differences | Soft FORMAT fields; waiver **W-L7-FORMAT** (dense soft AD/DP band) |
| PL vector not bitwise identical | Permanent soft-PL residual (**W-L7-FORMAT** / L9–L14 policy) |
| QUAL within tens of phred of Java | Not claimed genome-wide bitwise; comparator default `--qual-tol 50` |
| ALT order differs but same alleles + same called bases | Java Combine last→first vs discovery order; compare uses allele-identity GT |
| Header / INFO key ordering, annotation extras (e.g. ExcessHet) | Not part of the hard pilot gate |
| Runtime / Peak-RSS differences | Perf claims are separate ([`docs/ci/PERF_BENCHMARK_HOST.md`](ci/PERF_BENCHMARK_HOST.md)) |

### **Report these** (likely real deviation)

| Observation | Action |
|-------------|--------|
| Site present in only one engine (after matching `-L` / chrom naming) | File issue |
| Different REF/ALT set or different called genotype (allele-identity) | File issue |
| FILTER disagree on the same hard-filter expressions | File issue |
| `\|ΔF1\|` vs your truth above your agreed threshold (e.g. 0.02) with matched pin/interval | File issue |
| Crash / empty output / non-zero exit on inputs Java accepts | File as bug (may use bug template) |

Open with the dedicated template:

**[Equivalence deviation issue template](https://github.com/SynapticFour/gatk-rs/issues/new?template=equivalence_deviation.md)**

(`.github/ISSUE_TEMPLATE/equivalence_deviation.md`)

Include: gatk-rs commit, Java 4.4 pin confirmation, interval, sample IDs,
commands, and `REPORT.md` / a few loci (redact PHI).

---

## 6. Suggested pilot checklist

1. Build release `gatk-rs` on the machine that will run the pilot.  
2. Pick **one** interval already green in your Java pipeline.  
3. Run HC→…→Filtration (or stop at Genotype) on **both** engines with identical knobs.  
4. Run `compare_callsets.py` without truth; triage hard vs soft.  
5. If you have a truth VCF for that interval, re-run with `--truth`.  
6. Only escalate hard failures / large ΔF1 via the equivalence template.  
7. Expand interval / sample count only after the small window is clean.

---

## 7. Related links

| Doc / tool | Role |
|------------|------|
| [`CLAIM_MATRIX.md`](CLAIM_MATRIX.md) | What is / is not asserted |
| [`ARCHITECTURE.md`](ARCHITECTURE.md) | Scope spine + ADRs |
| [`scripts/pilot/compare_callsets.py`](../scripts/pilot/compare_callsets.py) | Standalone pilot compare |
| [`gatk-rs-equiv/`](../gatk-rs-equiv/) | Optional HC+truth orchestrator |
| [Equivalence dashboard](https://gatk-rs.github.io/gatk-rs/) | Public charts — **supplement**, not a substitute for your run |
| [Equivalence deviation template](https://github.com/SynapticFour/gatk-rs/issues/new?template=equivalence_deviation.md) | Report real diffs |
