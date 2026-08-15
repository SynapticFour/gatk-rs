# wall-losers Finalize dig (run 31857571738)

## Verdict

**Not an F1 math failure.** Finalize never reached hap.py / `gatk-rs-equiv`.

```
[giab] concat 8 Java shard VCF(s) → …/java.vcf
Checking the headers and starting positions of 8 files
Failed to open …/java.00_chr20_w09.vcf: not compressed with bgzip
exit 255
```

## Root cause

Ubuntu 24.04 **bcftools 1.19**: `bcftools concat -a` seeks via BGZF and **rejects plain `.vcf`**.
HC shards are uncompressed text. Reproduced locally:

| Command | Plain `.vcf` + contig header |
|---------|------------------------------|
| `concat -a` | fail: not compressed with bgzip |
| `concat` (no `-a`) | ok for non-overlapping shards |
| `concat --naive` | refuse: compressed only |

## Site counts (shards only; no truth F1)

| Window | Java | Rust | rust/java |
|--------|-----:|-----:|----------:|
| Σ 8 losers | 30784 | 22687 | **0.74** |

Undercall vs Java is real and belongs on the **L8 holdout / call-rate** track — separate from this concat bug. Once concat is fixed, expect the δ=0.02 equivalence gate to still fail until L8 closes; wall-losers remains a **wall** lane first.

## Fix

`giab_concat_vcfs`: ordered `bcftools concat` (no `-a`) for plain `.vcf`; keep `-a` for `.vcf.gz`/BCF; fall back to `concat_vcfs.py`.
