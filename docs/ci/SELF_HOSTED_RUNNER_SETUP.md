# Self-hosted runner for full-autosome GIAB validation

GitHub-hosted `ubuntu-latest` runners (~14 GB disk, ~6 hour soft limits in practice)
cannot hold **full 30× WGS BAMs + reference + Java/Rust outputs** for
`GIAB_MODE=autosomes`. Those runs belong on a **dedicated x86_64 cloud VM**
registered as a **repository-scoped** Actions runner with label
`gatk-rs-genomewide`.

Workflow: [`.github/workflows/genomewide-validation.yml`](../../.github/workflows/genomewide-validation.yml)  
VM control: [`scripts/ci/cloud_vmctl.sh`](../../scripts/ci/cloud_vmctl.sh)

> **Performance timings are a second machine.** Do **not** publish PairHMM /
> HC wall-clock or Peak-RSS marketing numbers from this correctness VM (or from
> GitHub-hosted runners). Provision and operate the quiet host separately:
> [`PERF_BENCHMARK_HOST.md`](PERF_BENCHMARK_HOST.md) +
> [`.github/workflows/benchmark.yml`](../../.github/workflows/benchmark.yml)
> (label **`gatk-rs-benchmark`**, secrets `AWS_PERF_INSTANCE_ID` /
> `HCLOUD_PERF_SERVER_ID`).

> **This VM costs real money.** GitHub-hosted minutes for public repos are free;
> this machine is billed by your cloud provider every hour it is **running**,
> plus **storage** while it exists (even when stopped). See [Cost model](#cost-model).

---

## Recommended machine

| Spec | Minimum | Notes |
|------|---------|--------|
| Arch | **x86_64** (`linux/amd64`) | Matches Broad GATK Docker (`linux/amd64`) |
| RAM | **64 GB** | Java HC + Rust HC + hap.py/RTG concurrently |
| vCPU | 8–16 dedicated | Hetzner **CCX43** (16 vCPU / 64 GB) or AWS **r6i.2xlarge** / Spot |
| Disk | **≥500 GB** SSD | Local NVMe **or** attached volume for BAM cache + `target/` + run dirs |
| OS | Ubuntu 22.04/24.04 LTS | |

Example SKUs:

- **Hetzner Cloud:** CCX43 + **500 GB Volume** mounted at `/data` (local disk on CCX43 is only ~360 GB).
- **AWS EC2 Spot:** `r6i.2xlarge` (or `r5.2xlarge`) + 500 GB **gp3** EBS. Prefer Spot for compute; keep the EBS volume.

---

## Security (non-negotiable)

1. **Repository runner only** — when registering, choose *this repo*, **not** the organization. Org-wide runners can pick up jobs from *any* repo that targets `self-hosted`, which is a lateral-movement risk.
2. **Dedicated OS user** (e.g. `actions`) — do **not** run the runner as `root`.
3. **Least privilege sudo** — the `actions` user should not have passwordless root for arbitrary commands. Install packages as root once during provisioning.
4. **Labels** — register with labels `self-hosted`, `linux`, `x64`, **`gatk-rs-genomewide`**. The workflow requires `gatk-rs-genomewide` so random self-hosted boxes cannot steal the job.
5. **Secrets** — cloud API tokens live in **GitHub Actions secrets**, never in the repo. Rotate if a workflow log ever prints them (the control script redacts tokens).
6. **Network** — prefer SSH key-only login; disable password auth; restrict SSH (`ufw`/security group) to your IP or a bastion.
7. **Ephemeral workdirs** — keep large BAMs under `/data` owned by `actions`; wipe run directories after artifact upload if disk pressure grows.
8. **Do not** install unrelated services on this VM. Treat it as a single-purpose build farm node.

---

## One-time VM provisioning

### 1. Create the VM + disk

**Hetzner (console or `hcloud`):**

```bash
# Example — adjust location/image names to your project
hcloud server create \
  --name gatk-rs-genomewide \
  --type ccx43 \
  --image ubuntu-24.04 \
  --ssh-key "$YOUR_SSH_KEY_NAME" \
  --label role=gatk-rs-genomewide

hcloud volume create --name gatk-rs-data --size 500 --server gatk-rs-genomewide --format ext4
```

**AWS:** launch Spot/on-demand instance with 500 GB gp3, tag `Role=gatk-rs-genomewide`.

Mount the data volume (example `/dev/sdb` → `/data`):

```bash
sudo mkdir -p /data
# if not already formatted by the cloud helper:
# sudo mkfs.ext4 -L gatkdata /dev/disk/by-id/...
echo 'LABEL=gatkdata /data ext4 defaults,nofail 0 2' | sudo tee -a /etc/fstab
sudo mount -a
sudo mkdir -p /data/{actions-runner,giab,cargo-target,work}
```

### 2. Create the `actions` user

```bash
sudo adduser --disabled-password --gecos 'GitHub Actions' actions
sudo usermod -aG docker actions   # after Docker is installed
sudo chown -R actions:actions /data
```

### 3. Install dependencies (as root)

```bash
sudo apt-get update
sudo apt-get install -y --no-install-recommends \
  curl ca-certificates git build-essential pkg-config \
  zlib1g-dev libbz2-dev liblzma-dev \
  samtools bcftools time openjdk-17-jre-headless \
  docker.io jq unzip

# Docker
sudo systemctl enable --now docker

# Rust (install for user actions)
sudo -u actions bash -lc 'curl --proto "=https" --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y'
```

Optional but useful: pre-pull heavy images while the machine is up for setup:

```bash
sudo -u actions docker pull --platform linux/amd64 us.gcr.io/broad-gatk/gatk:4.4.0.0
```

### 4. Register the GitHub Actions runner (repo-scoped)

In the GitHub UI:

**Settings → Actions → Runners → New self-hosted runner → Linux / x64**

Copy the token shown (short-lived). On the VM as `actions`:

```bash
sudo -u actions -i
cd /data/actions-runner
curl -fsSL -o actions-runner-linux-x64.tar.gz \
  -L https://github.com/actions/runner/releases/download/v2.321.0/actions-runner-linux-x64-2.321.0.tar.gz
# Use the current version from the GitHub "New runner" page if newer.
tar xzf actions-runner-linux-x64.tar.gz

./config.sh \
  --url https://github.com/<OWNER>/gatk-rs \
  --token <REGISTRATION_TOKEN> \
  --name gatk-rs-genomewide-1 \
  --labels gatk-rs-genomewide \
  --work /data/work \
  --unattended
```

Confirm when prompted that the runner is for **this repository only**.

Install as a systemd service so it comes back when the VM is **started** (power-on):

```bash
# still in /data/actions-runner — must be run with sudo
sudo ./svc.sh install actions
sudo ./svc.sh start
sudo systemctl status actions.runner.*.service
```

Point Cargo / GIAB caches at the big disk (as `actions`):

```bash
# ~/.bashrc or ~/.profile for user actions
export CARGO_TARGET_DIR=/data/cargo-target
export GIAB_CACHE_ROOT=/data/giab
```

### 5. Wire GitHub secrets (for start/stop automation)

Repo → **Settings → Secrets and variables → Actions**:

| Name | Type | Required for | Purpose |
|------|------|----------------|---------|
| `CLOUD_PROVIDER` | **Variable** | both | `hetzner` or `aws` |
| `HCLOUD_TOKEN` | Secret | Hetzner | Project API token (Read & Write) |
| `HCLOUD_SERVER_ID` | Secret | Hetzner | Numeric server id (`hcloud server list`) |
| `AWS_ACCESS_KEY_ID` / `AWS_SECRET_ACCESS_KEY` / `AWS_REGION` | Secret | AWS | IAM that can `ec2:StartInstances` + `StopInstances` only |
| `AWS_INSTANCE_ID` | Secret | AWS | `i-…` |

Optional:

| Secret | Purpose |
|--------|---------|
| `GH_RUNNER_WAIT_TOKEN` | PAT with `repo` admin (or fine-grained: Actions administration) so the start job can wait until the runner is `online` |

### 6. Smoke-test without a full autosome run

1. Manually start the VM (`./scripts/ci/cloud_vmctl.sh start` locally with env loaded).
2. Confirm runner shows **Idle** in GitHub → Settings → Actions → Runners.
3. Actions → **Genomewide validation (self-hosted)** → Run workflow with `mode=smoke`.
4. Confirm the stop job powers the VM off afterward.

---

## Cost control automation

[`scripts/ci/cloud_vmctl.sh`](../../scripts/ci/cloud_vmctl.sh) is called from the workflow:

1. **Job `start-vm`** (GitHub-hosted, free): `cloud_vmctl.sh start` → wait for runner online.
2. **Job `genomewide`** (`runs-on: [self-hosted, gatk-rs-genomewide]`): full `GIAB_MODE=autosomes` (or chosen mode).
3. **Job `stop-vm`** (`if: always()`, GitHub-hosted): `cloud_vmctl.sh stop` — **power off / stop**, **never terminate/delete**, so the disk image and `/data` cache survive.

Local dry-run:

```bash
export CLOUD_PROVIDER=hetzner
export HCLOUD_TOKEN=...
export HCLOUD_SERVER_ID=...
./scripts/ci/cloud_vmctl.sh status
./scripts/ci/cloud_vmctl.sh start
./scripts/ci/cloud_vmctl.sh stop
```

See [`scripts/ci/cloud_vmctl.env.example`](../../scripts/ci/cloud_vmctl.env.example).

---

## Cost model

Prices change; figures below are **order-of-magnitude** (EU, mid-2026). Always check your provider’s current list price before budgeting.

### Always-on storage (you pay this even when the VM is stopped)

| Item | Rough monthly cost |
|------|--------------------|
| Hetzner 500 GB Volume | ~€20–25 / mo |
| AWS 500 GB gp3 EBS | ~$40 / mo |
| Snapshots (optional) | extra |

Stopping the VM **does not** stop volume/EBS charges. Deleting the VM/volume does.

### Compute (only while powered on)

| SKU | Approx. hourly | Notes |
|-----|----------------|-------|
| Hetzner CCX43 (16 vCPU / 64 GB) | ~€0.22–€0.45 / h | Dedicated vCPU; EU regions; post-2026 list prices vary by location |
| AWS `r6i.2xlarge` on-demand | ~$0.50 / h | Region-dependent |
| AWS `r6i.2xlarge` Spot | ~$0.10–$0.25 / h | Can be interrupted — OK if you checkpoint BAM cache on EBS |

### Estimated cost **per full-autosome validation run**

Assumptions: **1 sample (HG001)**, `GIAB_MODE=autosomes`, Java + Rust HC + hap.py/RTG, wall clock **10–18 h** on CCX43-class hardware (first run longer if BAM cache cold).

| Component | Low | High |
|-----------|-----|------|
| Compute (Hetzner CCX43 @ ~€0.30/h × 10–18 h) | ~€3 | ~€8 |
| Compute (AWS Spot @ ~$0.15/h × 10–18 h) | ~$1.50 | ~$3 |
| Incremental egress / API | usually small on Hetzner | watch AWS egress if downloading BAMs from NCBI into AWS |
| **Amortized storage** (500 GB ÷ ~4 runs/month) | ~€5–10 / run | if you only run weekly |

**Ballpark total per weekly HG001 autosome run: about €8–20** (Hetzner) or **$8–25** (AWS on-demand), dominated by compute hours + shared monthly disk.

Multi-sample (`HG001,HG002,HG005`) roughly **×2–3** wall time and cost.

### Cost anti-patterns (avoid)

- Leaving the VM **running** 24/7 → hundreds of €/$ per month in compute alone.
- Using **org-wide** runners “for convenience” → security and surprise cross-repo spend.
- **Terminating** the instance each time → re-download 30× BAMs every run (more time = more compute cost).
- Triggering `genomewide-validation.yml` on every push (this workflow intentionally does **not**).

### Budget checklist

- [ ] Monthly disk budget approved (~€20–40).
- [ ] Weekly compute budget approved (~€10–60 depending on samples).
- [ ] Billing alerts set in Hetzner Console / AWS Budgets.
- [ ] Schedule kept weekly (or manual-only) — not daily.

---

## Operations cheat sheet

| Goal | Action |
|------|--------|
| Run full autosomes | Actions → **Genomewide validation (self-hosted)** → `mode=autosomes` |
| Skip paid VM this week | Disable the workflow schedule or cancel the run; keep VM stopped |
| Rotate runner token | Remove runner in GitHub UI; re-run `config.sh` on the VM |
| Emergency stop | `./scripts/ci/cloud_vmctl.sh stop` or provider console **Power off** |
| Destroy everything | Provider console: delete server **and** volume (stops storage charges) |

---

## Relation to other GIAB workflows

| Workflow | Runner | Scope |
|----------|--------|--------|
| `giab-genomewide.yml` | `ubuntu-latest` | Window / chr20–21 style (`ci-subset`) — free tier |
| `nightly-equivalence.yml` | `ubuntu-latest` | Trio E2E on sliced regions — free tier |
| **`genomewide-validation.yml`** | **`self-hosted` + `gatk-rs-genomewide`** | **`autosomes` / true WGS-scale** — **paid VM** |
| **`benchmark.yml`** | **`self-hosted` + `gatk-rs-benchmark`** | **Quiet-host PairHMM / Peak-RSS** — **separate paid VM** ([setup](PERF_BENCHMARK_HOST.md)) |
