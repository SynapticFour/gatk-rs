# GIAB multi-sample equivalence harness

External-facing runner: download GIAB truth sets, call **Java GATK4** (pinned) and **gatk-rs** on matching BAMs/intervals, score with **`gatk-rs-equiv`** (hap.py / RTG), capture wall time + peak RSS, publish a dashboard.

## What “genome-wide” means (`GIAB_MODE`)

| Mode | Meaning | Typical use |
|------|---------|-------------|
| `smoke` | Three ~50 kb windows (chr20 / chr21 / P12). P12 reads staged from **NA12878_20k** (same evidence class as P12 L* gates); chr20/21 from HG001 30×. Full-30× P12 is benchmark-host only (centromere-scale depth after positional DS). | M4 laptop / PR sanity |
| **`ci-subset` (default)** | **Full chr20 + full chr21** + **one 50 kb probe** on each other autosome (clamped to hs37d5 contig ends; chr19/22 use end-of-contig). HC runs as **matrix shards**: ~10 Mb windows of chr20/21 (`00_chr20_wNN` / `01_chr21_wNN`) plus `02_probes`, each × java/rust — stays under the GitHub-hosted **6 h** job hard cap | Nightly/weekly CI |
| `chr20-21` | Full chromosomes 20 and 21 only | Intermediate |
| `autosomes` | Full chr1–22 | Large hosts only — not for 16 GB laptops |

**`ci-subset` is not** “every base of the autosomes.” It is the practical CI definition of genome-adjacent coverage in this repo. Full autosomes require `GIAB_MODE=autosomes` and machine class beyond a MacBook Air M4.

## Scripts

```bash
# 1) Truth VCFs/BEDs (HG001/HG002/HG005) + selected stratification BEDs
./scripts/parity/giab/fetch_giab_truthsets.sh

# 2) Equivalence run (default: ci-subset, HG001)
GIAB_MODE=smoke GIAB_SAMPLES=HG001 \
  ./scripts/parity/giab/run_genomewide_equivalence.sh

# Multi-sample (heavy)
GIAB_MODE=ci-subset GIAB_SAMPLES=HG001,HG002,HG005 \
  ./scripts/parity/giab/run_genomewide_equivalence.sh

# 3) Joint-genotyping E2E (HC gVCF → CombineGVCFs → GenotypeGVCFs)
# Smoke (synthetic mini cohort; no BAM download):
TRIO_E2E_MODE=smoke ./scripts/parity/giab/run_trio_joint_genotyping_e2e.sh
# Real Ashkenazi trio (HG002/HG003/HG004) on a small interval:
TRIO_E2E_MODE=giab \
  TRIO_REFERENCE=... TRIO_INTERVAL=20:1000000-1050000 \
  TRIO_HG002_BAM=... TRIO_HG003_BAM=... TRIO_HG004_BAM=... \
  TRIO_TRUTH_VCF=... TRIO_TRUTH_BED=... \
  ./scripts/parity/giab/run_trio_joint_genotyping_e2e.sh
```

Outputs under `parity/giab/runs/<timestamp>_<mode>/`:

- `SCOPE.txt` — exact mode definition  
- `*/hc/{java,rust}.vcf` — callsets  
- `*/time/*.time.txt` — `/usr/bin/time` logs  
- `*/equiv/` — `gatk-rs-equiv` results  
- `dashboard/index.html` — publishable summary  

## CI / GitHub Pages

Workflow: `.github/workflows/giab-genomewide.yml` (weekly + `workflow_dispatch`).

- Default: `GIAB_MODE=ci-subset`, `GIAB_SAMPLES=HG001`  
- **Matrix pipeline** (GitHub-hosted 6 h hard cap): `prepare` → HC jobs per `shard × engine` (360 m each) → `finalize` (concat + hap.py/RTG)  
- Contig shards are **windowed** (`GIAB_HC_WINDOW_BP`, default **1 Mb** on hosted CI) so full chr20/21 stay under the 6 h cap **and** Rust HC fits ~16 GiB+swap (10 Mb / 2 Mb windows still OOMed dense HG001 30× shards)  

- **Reference handoff:** prepare prefers the pinned GitHub Release `giab-ref-v1` (`hs37d5.simple.fa.gz` + fai/dict via `scripts/parity/giab/fetch_hs37d5_release.sh`), caches + uploads the uncompressed FA to HC jobs (`GIAB_STAGE_REF=0`). FTP mirrors are fallback only.  
- **`GIAB_MODE=autosomes` is rejected** on this workflow — use [`genomewide-validation.yml`](../../../.github/workflows/genomewide-validation.yml) (self-hosted `gatk-rs-genomewide`)  
- Phases via `GIAB_PHASE=prepare|hc|finalize|all`; filters `GIAB_HC_SHARDS` / `GIAB_HC_ENGINES`  
- Uploads run artifacts; deploys run `dashboard/` to GitHub Pages under `/giab-ci/`
- Finalize prefers **hap.py** (Docker wrap) with RTG fallback; on green non-`smoke` runs, updates `docs/parity-site/data/history.json` and deploys the public Chart.js site

### Nightly trio joint E2E (HC → Combine → Genotype → Filter)

Workflow: `.github/workflows/nightly-equivalence.yml` (daily + `workflow_dispatch`).

```bash
./scripts/parity/giab/run_nightly_trio_equivalence.sh
```

- Regions: full chr20 + chr21 + capped hard slices (segdups / TR / alldifficult / MHC)
- BAM staging: `samtools view -L` remote slices only (no full WGS download)
- Publishes `docs/EQUIVALENCE_DASHBOARD.md` + Pages under `/equivalence/`
- Soft gate: F1 drop vs `docs/equivalence/baseline.json` opens an
  `equivalence-regression` issue (workflow does not hard-fail)  


### Public equivalence dashboard (GitHub Pages)

Static Chart.js site: [`docs/parity-site/`](../../../docs/parity-site/) →
**https://gatk-rs.github.io/gatk-rs/**  
Updated from `nightly-equivalence` + `genomewide-validation` via
`scripts/parity/giab/update_public_dashboard.py`. Scope (regions, samples,
Java GATK pin) is shown explicitly on the page.

### Full-autosome validation (paid self-hosted runner)

Workflow: [`.github/workflows/genomewide-validation.yml`](../../../.github/workflows/genomewide-validation.yml)  
Setup + **cost model**: operator-provisioned self-hosted runner for [`.github/workflows/genomewide-validation.yml`](../../../.github/workflows/genomewide-validation.yml) (internal runbook; not in the public tree)

- Runner label: `gatk-rs-genomewide` (repo-scoped self-hosted only)
- Default mode: `GIAB_MODE=autosomes` (not for GitHub-hosted 14 GB disks)
- VM is **started/stopped via API** around the job (`scripts/ci/cloud_vmctl.sh`) — compute bills only while powered on; disk still bills monthly
- Schedule: weekly Sunday 06:00 UTC + manual `workflow_dispatch` (never on push)

## Claim matrix

Signed claims live in `docs/CLAIM_MATRIX.md`.  
A green `ci-subset` (or stronger) run is required before asserting genome-adjacent GIAB equivalence; full autosome equivalence remains separate and requires the self-hosted workflow above.
