---
name: Equivalence deviation
about: Report a real Java↔Rust callset disagreement found in a pilot or local check
title: "[equivalence] "
labels: ["equivalence-deviation"]
---

## Summary

<!-- One paragraph: what disagreed (GT / missing site / FILTER / large F1 delta)? -->

## How you compared

- [ ] `scripts/pilot/compare_callsets.py` (direct Java↔Rust)
- [ ] `scripts/pilot/compare_callsets.py` + hap.py / vcfeval vs truth
- [ ] `gatk-rs-equiv`
- [ ] Other:

```bash
# Exact commands (gatk-rs + Java GATK + compare)
```

## Inputs

| Item | Value |
|------|--------|
| Reference assembly | e.g. GRCh37 / GRCh38 / hs37d5 |
| Interval(s) / BED | |
| Sample ID(s) | |
| BAM source (pipeline stage) | e.g. after BQSR / markdup |
| Java GATK version | **prefer 4.4.0.0** (see `docs/GATK_PINNED.env`) |
| gatk-rs commit / release | |

## Observed mismatch

| Class | Count / note |
|-------|----------------|
| Only-Java sites | |
| Only-Rust sites | |
| Allele / GT mismatch | |
| FILTER mismatch | |
| FORMAT soft drift only (AD/DP/GQ/PL) | |
| \|ΔF1\| vs truth (if run) | |

Paste a few concrete loci (CHROM POS REF ALT GT) — **not** whole VCFs:

```text
```

Attach or link the compare tool outputs (`REPORT.md` and `summary.json` from
`--out`; redact paths/PHI as needed).

## Why this is not just known soft drift

Confirm you read [`docs/PILOT_GUIDE.md`](../../docs/PILOT_GUIDE.md) § Expected deviations:

- [ ] GT / site presence differs (not only AD/DP/PL/QUAL)
- [ ] Or \|ΔF1\| above your agreed threshold after matching intervals + pin
- [ ] Same `-L`, same BAM, same ERC/GVCF settings on both engines

## Environment

- OS / arch:
- Threads used:
- Relevant claim scope (link [`docs/CLAIM_MATRIX.md`](../../docs/CLAIM_MATRIX.md) if known):

## Notes

<!-- Logs, Docker vs local JAR, Combine/Genotype/Filtration stage where it first appears -->
