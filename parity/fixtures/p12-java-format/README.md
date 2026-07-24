# P12 Java FORMAT fixtures (L4)

Frozen **per-variant FORMAT + site QUAL + INFO AF** from the Phase E Java baseline VCF.

| Artifact | Role |
|----------|------|
| [`all_sites.tsv`](./all_sites.tsv) | Primary contract — 66 rows, one per Java-only site |
| [`sites/`](./sites/) | Per-site TSV (`<pos>_<ref>_<alt>.tsv`) for targeted tests |
| [`cluster_sites.tsv`](./cluster_sites.tsv) | Subset: `92307324`, `92307327`, `92307359` |

**Source VCF / site list:** generated under gitignored `parity/reports/` by `./scripts/parity/run_p12_realworld_na12878_20k.sh` (not checked into git).

## Regenerate

```bash
python3 scripts/parity/extract_p12_java_format_fixtures.py --per-site
```

Re-run after refreshing the Java P12 VCF (e.g. `./scripts/parity/run_p12_realworld_na12878_20k.sh`).

## Tolerances

See [`docs/CLAIM_MATRIX.md`](../../../docs/CLAIM_MATRIX.md) § P12 FORMAT.

## Tests

| Test | When |
|------|------|
| `p12_java_format_fixture_contract` | Always (`cargo test p12_java_format_fixture`) |
| `p12_format_parity` | **L4.2:** `P12_PHASE_E=1`, `P12_L4_JAVA_FORMAT` unset + BAM/FASTA, `--ignored` (algorithmic) |
| | **L4.1 harness:** add `P12_L4_JAVA_FORMAT=1` (fixture overlay) |
| `p12_cluster_format_fixture` | Always (cluster row sanity) |
