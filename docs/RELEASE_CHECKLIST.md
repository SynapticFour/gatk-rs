# Release checklist

Reusable gate before **every** git tag / GitHub Release (`vX.Y.Z`).  
Copy this file’s tables into the PR or release notes and tick only what you verified for **this** tag.  
Authority for product claims remains [`CLAIM_MATRIX.md`](CLAIM_MATRIX.md).

Do **not** treat an empty or stale equivalence dashboard as evidence.  
Quiet CI visibility and public announcement are separate decisions — this checklist covers the **technical** release bar.

---

## 0. Meta

| Item | Value |
|------|--------|
| Tag / version | |
| Commit SHA on `main` | |
| Release manager | |
| Date (UTC) | |

---

## 1. Dashboard freshness

| Check | Pass? | Evidence (path, run URL, or artifact SHA) |
|-------|-------|-------------------------------------------|
| `docs/parity-site/data/history.json` has ≥1 run **after** this release’s commit (or README still honestly states empty-by-design and no F1 is marketed) | ☐ | |
| Scope panel on the Pages site shows pin / regions / samples with any published metrics | ☐ | |
| `docs/EQUIVALENCE_DASHBOARD.md` (if linked) matches Pages / does not invent numbers | ☐ | |
| Perf claims (if any) come only from the dedicated-host fair suite, not laptop Criterion | ☐ | |

---

## 2. Claim-matrix audit

Walk [`CLAIM_MATRIX.md`](CLAIM_MATRIX.md) and README:

| Check | Pass? | Notes |
|-------|-------|--------|
| Every **Yes** claim has fresh evidence (log / CI artifact) on or after this SHA | ☐ | |
| Pending / **Not yet signed** rows are still honest (no silent promotion) | ☐ | |
| README numbers and maturity wording match the matrix (no aspirational F1) | ☐ | |
| GATK pin in `docs/GATK_PINNED.env` matches what CI and docs cite | ☐ | |
| Waivers (W-*) still named where scoped behavior is claimed | ☐ | |

---

## 3. Required gates (open points)

Tick only with reproducible evidence. Scripts live under `scripts/parity/`.

| # | Gate | Pass? | Evidence |
|---|------|-------|----------|
| 1 | `main` on origin includes release workflows + parity-site; `gh` auth works for dispatch | ☐ | |
| 2 | L2 strict green on this line (`scripts/parity/run_hc_full_parity_l2.sh` or CI equivalent) | ☐ | |
| 3 | `giab-genomewide.yml` `workflow_dispatch` with `mode=ci-subset` (HG001 at minimum) completed green; artifact SHA recorded | ☐ | |
| 4 | `nightly-equivalence.yml` (or agreed subset) completed so Pages history is non-empty **or** empty-by-design is still explicit in README | ☐ | |
| 5 | Joint / filter parity on this SHA: CombineGVCFs, GenotypeGVCFs, VariantFiltration scripts green | ☐ | |
| 6 | GitHub Pages serves `docs/parity-site` with scope-first layout (no orphaned F1-only tables) | ☐ | |
| 7 | Announcement / pins / external links deferred until this checklist is fully green (marketing is optional and separate) | ☐ | |

---

## 4. Engineering hygiene

| Check | Pass? |
|-------|-------|
| `cargo test --workspace` (or CI `pr-check` / `quality`) green on the tag SHA | ☐ |
| `python3 scripts/dev/check_doc_links.py` green | ☐ |
| Unsafe-pattern ratchets (`scripts/dev/count_unsafe_patterns.sh check`) green | ☐ |
| `CHANGELOG.md` has an entry for this version | ☐ |
| No secrets, machine-local absolute paths, or internal budget notes in files newly added for the tag | ☐ |

---

## 5. Verdict

| Question | Answer |
|----------|--------|
| Technical release (tag) allowed? | ☐ Yes / ☐ No |
| Public announcement allowed? | ☐ Yes / ☐ No (requires §1–§3 fully green) |

If **No**, leave the tag unpublished (or mark pre-release) and list blockers below.

**Blockers:**

-
-
-
