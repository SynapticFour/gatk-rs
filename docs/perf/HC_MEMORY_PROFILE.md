# HaplotypeCaller memory profile (reproducible)

**Generated (UTC):** `20260724T051610Z`  
**Host:** `Darwin 25.5.0 arm64`  
**Runner script:** [`scripts/perf/run_hc_memory_profile.sh`](../../scripts/perf/run_hc_memory_profile.sh)  
**Raw run directory:** `docs/perf/runs/20260724T051610Z/`

> **Scope warning:** This profile uses the checked-in p4 smoke fixture
> (`parity/fixtures/sample.bam` + `reference.fa`, interval `chr1:1-32`).
> Absolute Peak-RSS is dominated by runtime/JVM fixed costs on such a tiny
> window. **Do not** advertise these numbers as genome-wide “X% less memory”
> without re-measuring on a realistic GIAB shard.

## Peak-RSS (side by side)

| Engine | Peak RSS | Wall time |
|--------|----------|-----------|
| **gatk-rs** (Rust release) | **9.27 MiB (9488 KiB)** | 0.29 s |
| **Java GATK 4.4.0.0** | **451.92 MiB (462764 KiB)** | 2.644 s |

| Java / Rust Peak-RSS | 48.77× |
| Rust as fraction of Java Peak-RSS | 2.1% |
| Absolute delta (Java − Rust) | 442.65 MiB |


## Exact commands

### Rust

```bash
cargo build -p gatk-cli --release --locked
# rustc: rustc 1.88.0 (6b00bc388 2025-06-23)
# cargo: cargo 1.88.0 (873a06493 2025-05-10)
# git: ecb97e5
<release-bin>/gatk-rs HaplotypeCaller \
  -R <repo>/parity/fixtures/reference.fa \
  -I <repo>/parity/fixtures/sample.bam \
  -O /tmp/rust.hc.vcf \
  -L chr1:1-32
```

Time capture (this host): local `time` log under gitignored `docs/perf/runs/<timestamp>/` (macOS `/usr/bin/time -l` or GNU `/usr/bin/time -v`).

### Java GATK 4.4

- Pin: `GATK_PINNED_SHA=2dbc025821bc5f686c423ff332a41e6cef892a77` (`docs/GATK_PINNED.env`)
- Image: `us.gcr.io/broad-gatk/gatk:4.4.0.0`
- JVM options (pipeline-realistic): `-Xms1g -Xmx4g`

```bash
# Re-run via the harness (preferred):
./scripts/perf/run_hc_memory_profile.sh
# Exact docker/java cmdline is written under gitignored docs/perf/runs/<timestamp>/ by the harness.
```

Time capture: gitignored `docs/perf/runs/<timestamp>/java.time.txt` written by the harness.  
When Docker is used, Peak-RSS is sampled from `/proc/*/status` **VmHWM**
for `java`/`gatk` **inside** the Linux container (the Broad 4.4 image has no
GNU `/usr/bin/time`). Host `time docker …` is never used for RSS.

## Optional deeper profiling

- **macOS Instruments:**  
  `xcrun xctrace record --template 'Allocations' --output /tmp/rust.allocations.trace --launch -- <release-bin>/gatk-rs HaplotypeCaller -R <repo>/parity/fixtures/reference.fa -I <repo>/parity/fixtures/sample.bam -O /tmp/rust.hc.vcf -L chr1:1-32`
- **Linux heaptrack (Docker):**  
  mount the release binary and fixture into a heaptrack image; keep Peak-RSS
  from this script as the primary comparable number.

## Re-run

```bash
./scripts/perf/run_hc_memory_profile.sh
# optional overrides:
#   JAVA_XMX=4g JAVA_XMS=1g HC_MEM_INTERVAL=chr1:1-32 ./scripts/perf/run_hc_memory_profile.sh
```
