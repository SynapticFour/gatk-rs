# Dedicated performance-benchmark host

Timing numbers for PairHMM / HaplotypeCaller must come from a **quiet,
non-burstable, single-tenant-ish** machine — **not** from:

| Wrong place | Why |
|-------------|-----|
| GitHub-hosted `ubuntu-latest` | Shared/virtualized CPUs, noisy neighbors, no fixed core allocation |
| Self-hosted **`gatk-rs-genomewide`** VM ([Prompt H](SELF_HOSTED_RUNNER_SETUP.md)) | Sized for correctness (RAM + disk), often busy with GIAB; different SKU; not tuned for low-jitter timing |

This document provisions a **second** repository-scoped runner with label
**`gatk-rs-benchmark`**, started/stopped by
[`.github/workflows/benchmark.yml`](../../.github/workflows/benchmark.yml).

> **This VM costs money** while powered on (plus disk while it exists). Keep it
> **stopped** between weekly / manual runs. See [Cost model](#cost-model).

---

## Recommended machine (committed choice)

**Primary (preferred for AVX-512 / GKL-class PairHMM):**

| Spec | Value |
|------|--------|
| Provider SKU | **AWS `c7i.2xlarge`** (fallback: **`c6i.2xlarge`**) |
| Arch | **x86_64** (`linux/amd64`) |
| vCPU | **8 dedicated** Intel (non-burstable; **not** T-family / shared) |
| SIMD | **AVX2 + AVX-512** on c7i/c6i |
| RAM | 16 GiB (enough for HC shard benches; not for full WGS GIAB) |
| Disk | **≥100 GB gp3** (fixtures + `target/` + Criterion; keep separate from genomewide `/data`) |
| OS | Ubuntu 22.04/24.04 LTS |
| Purchase option | **On-demand** for scheduled benches (Spot allowed for ad-hoc only — interrupt risk) |

**Explicitly forbidden for this role:** AWS **T3/T4g/…** (credit throttling), any
**shared/burstable** vCPU product, GitHub-hosted runners, the
`gatk-rs-genomewide` correctness VM.

**Optional EU alternative (AVX2 only):** Hetzner Cloud **CCX23/CCX33** (dedicated
AMD). Use only if you accept **no AVX-512** (Rust NEON/AVX2 paths still run;
Java GKL AVX-512 path will not). Prefer AWS c7i when publishing SIMD speedups
that claim AVX-512.

---

## Separation from the genomewide runner

| | Correctness VM | Benchmark VM |
|--|----------------|--------------|
| Label | `gatk-rs-genomewide` | **`gatk-rs-benchmark`** |
| Workflow | `genomewide-validation.yml` | **`benchmark.yml`** |
| Secrets | `HCLOUD_SERVER_ID` / `AWS_INSTANCE_ID` | **`HCLOUD_PERF_SERVER_ID` / `AWS_PERF_INSTANCE_ID`** |
| Tuning | none required | governor `performance`, pinned cores, HT off |
| Typical SKU | CCX43 / r6i (RAM) | **c7i.2xlarge** (compute + AVX-512) |

Never register both roles on one machine. Never reuse the genomewide volume as
the sole bench disk (I/O contention).

---

## Security (same bar as Prompt H)

1. **Repository-scoped** runner only (not org-wide).
2. Dedicated OS user `actions` — not root.
3. Labels: `self-hosted`, `linux`, `x64`, **`gatk-rs-benchmark`**.
4. Cloud tokens only in GitHub Actions secrets.
5. SSH key-only; no unrelated services.

---

## One-time provisioning

### 1. Create the VM

**AWS (console / CLI example):**

```bash
# Illustrative — pick AMI, subnet, SG in your account
aws ec2 run-instances \
  --instance-type c7i.2xlarge \
  --count 1 \
  --tag-specifications 'ResourceType=instance,Tags=[{Key=Role,Value=gatk-rs-benchmark}]' \
  ...
```

Confirm in the instance details / `lscpu`:

```text
Flags: ... avx2 ... avx512f ...
```

**Hetzner (AVX2-only fallback):**

```bash
hcloud server create \
  --name gatk-rs-benchmark \
  --type ccx33 \
  --image ubuntu-24.04 \
  --ssh-key "$YOUR_SSH_KEY_NAME" \
  --label role=gatk-rs-benchmark
```

### 2. OS packages + Docker + Rust

```bash
sudo apt-get update
sudo apt-get install -y --no-install-recommends \
  curl ca-certificates git build-essential pkg-config \
  zlib1g-dev libbz2-dev liblzma-dev \
  samtools bcftools time openjdk-17-jre-headless \
  docker.io jq linux-tools-common linux-tools-generic \
  cpufrequtils util-linux

sudo systemctl enable --now docker
sudo adduser --disabled-password --gecos 'GitHub Actions' actions
sudo usermod -aG docker actions
sudo -u actions bash -lc 'curl --proto "=https" --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y'

# Pre-pull pinned GATK (docs/GATK_PINNED.env)
sudo -u actions docker pull --platform linux/amd64 us.gcr.io/broad-gatk/gatk:4.4.0.0
```

### 3. Host tuning (document + apply)

Run as root after every boot (or install as a oneshot systemd unit):

```bash
sudo ./scripts/perf/tune_host.sh apply
./scripts/perf/tune_host.sh status   # as any user — prints governor / SMT / isolcpus
```

What `tune_host.sh` does / documents:

| Knob | Target |
|------|--------|
| CPU governor | **`performance`** on all online CPUs |
| Turbo | leave on unless you need stricter thermal stability (document if off) |
| Hyperthreading / SMT | **prefer off** for published benches (`/sys/devices/system/cpu/smt/control`) |
| Process pinning | benches wrap with `taskset -c <cores>` (see suite script) |

Persist SMT-off across reboot via BIOS/UEFI when possible; sysfs toggles may
reset on some clouds.

### 4. Verify Java GATK native AVX PairHMM

Pinned image: `us.gcr.io/broad-gatk/gatk:4.4.0.0`
(`GATK_PINNED_SHA=2dbc0258…` in [`docs/GATK_PINNED.env`](../GATK_PINNED.env)).

```bash
./scripts/perf/verify_java_avx_pairhmm.sh
```

This runs a tiny HC with `--pair-hmm-implementation FASTEST_AVAILABLE` and
**fails** unless the log shows a native vector PairHMM path (e.g. `AVX` /
`VectorLoglessPairHMM`) and **not** a silent Java `Log10PairHMM` /
`LOGLESS_HMM` software fallback.

### 5. Register the runner

GitHub → **Settings → Actions → Runners → New self-hosted runner → Linux / x64**

```bash
sudo -u actions -i
mkdir -p ~/actions-runner && cd ~/actions-runner
# download + config.sh from the UI, then:
./config.sh --url https://github.com/OWNER/gatk-rs --token XXX \
  --labels self-hosted,linux,x64,gatk-rs-benchmark \
  --name gatk-rs-benchmark-$(hostname -s)
sudo ./svc.sh install
sudo ./svc.sh start
```

### 6. GitHub secrets / variables (perf VM only)

| Name | Type | Purpose |
|------|------|---------|
| `CLOUD_PROVIDER_PERF` | Variable | `aws` (preferred) or `hetzner` |
| `AWS_ACCESS_KEY_ID` / `AWS_SECRET_ACCESS_KEY` / `AWS_REGION` | Secret | Shared IAM OK if scoped to Start/Stop |
| **`AWS_PERF_INSTANCE_ID`** | Secret | **`i-…` of the c7i/c6i bench instance** (≠ genomewide id) |
| `HCLOUD_TOKEN` | Secret | Hetzner (if used) |
| **`HCLOUD_PERF_SERVER_ID`** | Secret | Hetzner bench server id (≠ `HCLOUD_SERVER_ID`) |
| `GH_RUNNER_WAIT_TOKEN` | Secret | Optional: wait until label `gatk-rs-benchmark` is online |

Reuse of genomewide `AWS_INSTANCE_ID` / `HCLOUD_SERVER_ID` is a **misconfiguration**.

---

## Running benches

### Automated

Actions → **Performance benchmark (dedicated host)** → Run workflow  
or wait for the weekly schedule (see `benchmark.yml`).

The workflow:

1. Starts the **perf** VM (`cloud_vmctl.sh` with perf secrets).
2. Runs on `[self-hosted, gatk-rs-benchmark]`.
3. Captures host specs → `docs/perf/HOST_SPECS.md` + JSON.
4. Re-verifies Java AVX PairHMM for `FASTEST_AVAILABLE`.
5. Runs the **fair HC comparison**
   ([`scripts/perf/run_fair_hc_comparison.sh`](../../scripts/perf/run_fair_hc_comparison.sh)):
   4 configs × 3 region sizes × ≥5 repeats → median ± stdev
   (wall / user / sys / Peak-RSS / optional RAPL energy).
6. Updates the public dashboard **Performance** tab
   (`docs/parity-site/data/perf_history.json`).
7. Uploads artifacts; **always** stops the VM.

Primary speedup baseline is always Java **`FASTEST_AVAILABLE`**.
`LOGLESS_CACHING` is measured as a secondary software reference only.

### Manual (SSH)

```bash
sudo ./scripts/perf/tune_host.sh apply
PERF_REPEATS=5 PERF_THREADS=1 ./scripts/perf/run_fair_hc_comparison.sh
# or full suite (fair + Peak-RSS smoke+realistic + optional microbenches):
PERF_SKIP_MEMORY=0 ./scripts/perf/run_dedicated_benchmark_suite.sh
# Peak-RSS only (smoke + 2 Mb GIAB-dense window — public memory claims):
./scripts/perf/capture_host_specs.sh
HC_MEM_PROFILES=smoke,realistic HC_MEM_THREADS=1 ./scripts/perf/run_hc_memory_profile.sh
```

Public “X% less memory” claims must cite the **realistic** profile from
[`docs/perf/HC_MEMORY_PROFILE.md`](../perf/HC_MEMORY_PROFILE.md) measured on this
host (never the trivial `chr1:1-32` smoke table).

---

## Publishing numbers

Every public timing / Peak-RSS claim must point at:

1. [`docs/perf/HOST_SPECS.md`](../perf/HOST_SPECS.md) (or the stamped copy under `docs/perf/runs/…`)
2. The measurement script + raw logs under `docs/perf/`
3. GATK pin (`docs/GATK_PINNED.env`) and confirmed PairHMM implementation line from the Java log

MacBooks / GitHub-hosted / genomewide VMs may still produce **dev** profiles
(e.g. local M4 NEON); label those hosts clearly and do **not** mix them into
“production SIMD vs Java GKL” marketing tables.

---

## Cost model (order of magnitude)

| Item | Notes |
|------|--------|
| AWS `c7i.2xlarge` on-demand | ~$0.30–$0.45 / h (region-dependent) |
| 100 GB gp3 | ~$8 / mo while volume exists |
| Weekly suite wall | typically **15–45 min** compute if caches warm |

**Ballpark:** a few USD/EUR per weekly run + ~$8/mo disk — far cheaper than the
genomewide correctness VM. Still: **stop after every run**.

---

## Ops cheat sheet

| Goal | Action |
|------|--------|
| Run suite | Actions → **Performance benchmark (dedicated host)** |
| Skip this week | Disable schedule / cancel run; keep VM stopped |
| Emergency stop | `CLOUD_PROVIDER=aws AWS_INSTANCE_ID=$AWS_PERF_INSTANCE_ID ./scripts/ci/cloud_vmctl.sh stop` |
| Confirm AVX in Java | `./scripts/perf/verify_java_avx_pairhmm.sh` |
