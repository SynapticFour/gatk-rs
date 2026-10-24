# User guide

How to build `gatk-rs` and run the germline short-variant spine. Claims about
agreement with GATK are in [`CLAIM_MATRIX.md`](CLAIM_MATRIX.md), not in this
page. Comparing a Java GATK 4.4 run with a Rust run on your own BAMs is
[`PILOT_GUIDE.md`](PILOT_GUIDE.md).

The supported interface is the `gatk-rs` binary (`gatk-cli`). The
`gatk-haplotypecaller` crate is the engine that binary calls. Dump and oracle
modules are not part of a default release build’s behavior; see
[`ARCHITECTURE.md`](ARCHITECTURE.md).

## Build

From a clone, with Rust stable (see the repository `rust-toolchain.toml`) and
a C toolchain:

```bash
cargo build -p gatk-cli --release --locked
./target/release/gatk-rs --version
./target/release/gatk-rs --help
```

Workspace version is **0.1.0**. On a memory-tight machine,
`CARGO_BUILD_JOBS=1` matches what CI sets.

Optional allocator, off unless you pass the feature:

```bash
cargo build -p gatk-cli --release --locked --features jemalloc
```

`--help` lists the workflow tools. Experimental SAM/BAM utilities
(`PrintReads`, `Validate`, `CountReadsInRegion`, and similar) stay callable
and are hidden from that help. `gatk-rs --list` shows them. For ordinary
file operations use samtools or bcftools
([ADR 0002](adr/0002-remove-gatk-tools.md)).

## HaplotypeCaller

```bash
samtools faidx ref.fa
samtools index sample.bam

./target/release/gatk-rs HaplotypeCaller \
  -R ref.fa \
  -I sample.bam \
  -O sample.vcf \
  -L chr20:10000000-10050000
```

| Flag | Meaning | Default |
|------|---------|---------|
| `-R`, `--reference` | Reference FASTA | required |
| `-I`, `--input` | Indexed BAM. HaplotypeCaller opens it with HTSlib’s indexed reader. Repeatable. A SAM path is not this input. | required |
| `-O`, `--output` | Plain-text VCF path. `.vcf`, `.g.vcf`, and other names are written as text. Names ending in `.vcf.gz`, `.bgz`, `.bgzip`, or `.bcf` are rejected. | required |
| `-L`, `--intervals` | Interval; repeatable | whole reference if omitted |
| `--output-mode` (`VCF`, `GVCF`, `BP_RESOLUTION`) | Variant records vs reference confidence | `VCF` |
| `--emit-ref-confidence` (`NONE`, `GVCF`, `BP_RESOLUTION`) | GATK `-ERC` equivalent. There is no single-letter `-ERC` flag. `GVCF` and `BP_RESOLUTION` replace `--output-mode`. | unset |
| `--stand-call-conf` | Min confidence to call | `30.0` |
| `--stand-emit-conf` | Min confidence to emit | `10.0` |
| `--min-mapping-quality` | Mapping-quality floor | `20` |
| `--min-base-quality-score` | Base-quality floor | `10` |
| `--pair-hmm` | Omitted selects `FASTEST_AVAILABLE`: a supported host SIMD kernel (AVX2, NEON, or packed f64) when one is available, otherwise Logless scalar. Also accepts `LOG10_PAIRHMM`, `LOGLESS_HMM`, `AVX` / `SIMD`, `SIMD_F32`, and `WAVEFRONT`. `WAVEFRONT` is accepted by the parser. It is not the omitted-flag value. | omitted (`FASTEST_AVAILABLE`) |
| `--alleles` / `-A` | VCF of alleles to force (claim-matrix T3-5; not a full GATK `-alleles` sign-off) | unset |
| `--threads` / `-t` / `--nt` | Rayon worker count. When the flag is set, including `--threads 0`, it overrides `RAYON_NUM_THREADS`. `--threads 0` uses the logical CPU count. | omitted: `RAYON_NUM_THREADS` if that value is a positive integer, otherwise the logical CPU count |
| `--java-options` | Accepted on the top-level command. Does not start a JVM and does not size the Rust process. `GatkConfig::with_java_options` parses the string. If it contains `-Xmx`, validation rejects a memory token that does not start with a digit and end in `M` or `G` before calling starts. Other tokens are stored and not applied. | unset |

gVCF for one sample:

```bash
./target/release/gatk-rs HaplotypeCaller \
  -R ref.fa -I sample.bam -O sample.g.vcf \
  -L chr20:10000000-10050000 \
  --emit-ref-confidence GVCF \
  --stand-call-conf 30 \
  --threads 8
```

Repeating `-I` does not reproduce Java’s merged multi-sample HaplotypeCaller.
Call one BAM to one gVCF, then combine.

## Reference and reads

- Contig names in the FASTA, the BAM, and `-L` must refer to the same
  assembly. The dictionary loader tolerates an optional `chr` prefix when it
  looks a name up.
- If `ref.fa.fai` is present, contig names and lengths come from that index.
  Without it, the FASTA is scanned. A Picard `.dict` is not the file this
  loader reads.
- HaplotypeCaller opens the BAM with HTSlib’s indexed reader. The project’s
  fixtures and scripts use a `.bai`. A missing index fails at open or fetch,
  before variants are written.
- Reads should already be the BAMs you would pass to GATK HaplotypeCaller.
  This project does not recalibrate base qualities (no BQSR).

There is no example genome in the git checkout that is large enough to
demonstrate a realistic call. The tracked files under
[`../parity/fixtures/`](../parity/fixtures/) are tiny smoke fixtures (tens of
bases), not a recommended biological run. Real NA12878 assets used by the
P12 scripts live under `parity/realworld/assets/` and are gitignored.

## Cohort steps

Combine per-sample gVCFs, then joint-genotype, then hard-filter. The signed
cohort evidence is a synthetic ladder (recommended ≤ 100 samples on 10 kb /
400 SNPs), not WGS and not GenomicsDB. See the claim matrix before using
this on a large cohort: Combine holds the gVCFs in memory.

```bash
./target/release/gatk-rs CombineGVCFs \
  -R ref.fa \
  -V sample1.g.vcf -V sample2.g.vcf \
  -O combined.g.vcf \
  -L chr20:10000000-10050000

./target/release/gatk-rs GenotypeGVCFs \
  -R ref.fa \
  -V combined.g.vcf \
  -O genotyped.vcf \
  -L chr20:10000000-10050000 \
  --stand-call-conf 30

./target/release/gatk-rs VariantFiltration \
  -V genotyped.vcf \
  -O filtered.vcf \
  --preset snp
```

`--preset indel` selects the indel table. `--preset` may be combined with
explicit pairs. Each `--filter-expression` needs a matching
`--filter-name`. Prefer one annotation per expression.

SNP preset expressions (from `gatk_core::variant_filtration::GATK_HARD_FILTER_SNP`):

- `QD < 2.0`, `QUAL < 30.0`, `SOR > 3.0`, `FS > 60.0`, `MQ < 40.0`,
  `MQRankSum < -12.5`, `ReadPosRankSum < -8.0`

Indel preset: `QD < 2.0`, `QUAL < 30.0`, `FS > 200.0`,
`ReadPosRankSum < -20.0`, `SOR > 10.0`. `InbreedingCoeff` is not in that
table. VariantFiltration writes FILTER tags. It is not VQSR. `-R` is
accepted and unused for these hard filters.

## Outputs

| Mode | File |
|------|------|
| Default | VCF of called variants |
| `--emit-ref-confidence GVCF` or `--output-mode GVCF` | gVCF, including reference-confidence blocks, for CombineGVCFs |
| `BP_RESOLUTION` | Reference confidence at base resolution |
| VariantFiltration | VCF with FILTER annotations. Records are not deleted by a failing expression. |

HaplotypeCaller, CombineGVCFs, GenotypeGVCFs, and VariantFiltration write that
VCF as plain text. BGZF, BCF, and an index are not implemented. A path whose
file name ends in `.vcf.gz`, `.bgz`, `.bgzip`, or `.bcf` is rejected before
the file is created or truncated.

Do not assume every INFO or FORMAT field matches GATK 4.4 outside the scopes
in the claim matrix. On the signed mid-B region the recorded Java-equivalent
fields are described in [`PARITY_MILESTONE_6R.md`](PARITY_MILESTONE_6R.md).
Elsewhere, a QUAL, PL, or annotation difference is an open comparison until
a gate says otherwise. Header text and tool versions differ from Java GATK
on purpose.

Active regions are processed with Rayon. Before the VCF is written, region
batches are sorted by genomic position and deduplicated. CI runs
`gatk-cli` test `rayon_thread_determinism_test`, which requires the VCF to
be byte-identical across thread counts 1, 2, 4, 8, and 16 after volatile
header lines are normalized. That is a repository gate on the fixture that
test uses. It is not a promise about every BAM.

## Reproducibility

Record the git SHA (`gatk-rs --version` prints the crate version, not the
commit), the exact command, the reference build, and the BAM. For a Java
comparison, use GATK **4.4.0.0** as pinned in [`GATK_PINNED.env`](GATK_PINNED.env).
Other 4.x releases are a different baseline.

Omitting `--pair-hmm` selects `FASTEST_AVAILABLE`, not the log10 scalar
kernel. For an equivalence comparison that intends that kernel, pass
`--pair-hmm LOG10_PAIRHMM`. `LOGLESS_HMM`, `AVX`, and `WAVEFRONT` are other
explicit choices. None of them is a signed speed or accuracy result.

QualByDepth on the mid-B contract follows GATK’s `Random` seed `47382911`.
Do not expect QD to match if the set of sites drawn from that stream differs.

## Resource use

There is no published genome-wide memory or runtime number. The fair timing
contract is [`perf/FAIR_HC_COMPARISON.md`](perf/FAIR_HC_COMPARISON.md) and has
no dedicated-host result yet.

Practical controls that exist in the code:

| Control | Effect |
|---------|--------|
| `--threads N` | Size of the Rayon pool |
| `GATK_RS_HC_SEQUENTIAL=1` | One assembly region at a time |
| `--features jemalloc` at build time | Optional allocator. Off by default. |
| `GATK_RS_HC_RSS_ABORT_MIB=<MiB>` | If set, k-best returns a memory error once process RSS reaches that threshold. Unset means off. Not a tuning recommendation. |
| `-L` a bounded interval | The main way to keep a first run small |

`GATK_RS_EXPERIMENTAL_KBEST_POLICY` changes SeqGraph frontier behavior.
Leave it unset for the production default (`unbounded_diagnostic`).
Do not set `legacy_1024` or `byte_budget` to “match GATK” unless you are
doing the investigation described in
[`parity/KBEST_RESOURCE_POLICY.md`](parity/KBEST_RESOURCE_POLICY.md).

## Failures

| Symptom | What to check |
|---------|----------------|
| Indexed BAM open/fetch error | `samtools index` the file you passed to `-I`. The path must be the BAM, not the index. |
| Unknown contig or empty interval | FASTA contig names vs `-L` and the BAM header. |
| Header-only VCF | `GATK_RS_HC_SCAFFOLD_OUTPUT=1` forces a header-only scaffold. Unset it for real calls. |
| Very large RSS | Restrict `-L`. A missing `.fai` makes reference setup scan the FASTA. `GATK_RS_HC_SEQUENTIAL=1` serializes regions. |
| Output that Java would not emit, or the reverse | Not a crash. Compare with [`PILOT_GUIDE.md`](PILOT_GUIDE.md) and read the claim matrix before calling it a bug. |
| `Validate` / `PrintReads` missing from `--help` | Hidden utilities. `gatk-rs --list`. Prefer samtools for new work. |

A non-zero exit from the binary is the failure signal. Logs use the `tracing`
subscriber; `RUST_LOG=error` is what the parity scripts set.

## What this page does not demonstrate

No command here is a genome-wide validation. The smallest useful workflow is
one indexed BAM, one interval, one VCF. Equivalence of that VCF to GATK is
a separate, scoped question: [`PARITY.md`](PARITY.md).
