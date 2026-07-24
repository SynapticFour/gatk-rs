# E2E real-world parity corpus (E2E.3 / J-D01 / J-D03 / J.2.2)

Small indexed BAM + reference intervals for **L3** full HC VCF parity (`run_hc_realworld_parity.sh`).

## Bundled non-vacuous corpus (default)

Uses **p11_java_positive** so `call_region` emits a biallelic SNP (Java and Rust agree on CHROM/POS/REF/ALT).

| Asset | Path |
|-------|------|
| Reference | `parity/fixtures/p5_live_reference.fa` |
| BAM (from SAM) | `parity/fixtures/p11_java_positive.sam` → `parity/build/sam-indexed-bam/p11_java_positive.bam` |
| Interval | `chrLive:1-63` |
| Golden VCF | `parity/fixtures/hc-full-parity/e2e-real/expected/p11_java_positive_chrlive_golden.vcf` |

```bash
export PARITY_HC_REALWORLD_STRICT=1
./scripts/parity/run_hc_realworld_parity.sh
```

Refresh golden after intentional VCF emission changes:

```bash
./scripts/parity/run_hc_realworld_golden_refresh.sh
```

L3 strict mode checks CHROM/POS/REF/ALT plus QUAL/FILTER/INFO vs Docker Java (`--require-java-l3`).
Site QUAL uses GATK `GenotypeLikelihoodCalculator` log10 GLs + `AlleleFrequencyCalculator` (not legacy `2×lr` parity sums).
J-D06: `j2-vcf` row `p11_java_positive_chrlive` asserts `call-region-vcf` emits the same site.

## Legacy p5_live smoke (vacuous on call-region path)

`p5_live_case_snp` + `chrLive:1-24` remains in L2 (`j2-vcf` expects `variant_emitted false`); use p11 for L3 proof.

## NA12878 / P12 (optional, L6)

```bash
export PARITY_HC_REALWORLD_BAM=/path/to/NA12878.chr20.bam
export PARITY_HC_REALWORLD_REF=/path/to/hg38.fa
export PARITY_HC_REALWORLD_INTERVAL=chr20:10000000-10005000
export PARITY_HC_REALWORLD_GOLDEN_VCF=/path/to/java_hc_golden.vcf
```

## Gates

- L2: `j2-vcf` `call-region` row (`p5_snp_chrlive`)
- L3: `scripts/parity/run_hc_realworld_parity.sh` (golden byte + Java identity in strict CI)
- CI.1: `scripts/parity/run_hc_full_parity_ci1.sh`
