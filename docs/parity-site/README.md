# Public equivalence & performance dashboard (GitHub Pages)

Static site (HTML/CSS/JS + Chart.js CDN). No build step.

Tabs: **Equivalence** (hap.py/RTG) · **Performance** (fair HC timings vs Java
`FASTEST_AVAILABLE` on the dedicated quiet host).

| File | Role |
|------|------|
| `index.html` / `style.css` / `app.js` | UI (two tabs) |
| `data/history.json` | Equivalence run history. The live copy is on `gh-pages`, not on `main`. |
| `data/latest.json` | Most recent equivalence snapshot, published with the history file |
| `data/perf_history.json` | Fair HC comparison history (`benchmark.yml`) |
| `data/perf_latest.json` | Most recent performance snapshot |
## Update locally

```bash
python3 scripts/parity/giab/update_public_dashboard.py \
  --source nightly \
  --workflow nightly-equivalence \
  --pages-history-status absent \
  --json path/to/happy_summary.json \
  --site-dir docs/parity-site \
  --commit-sha "$(git rev-parse HEAD)"

# Preview
cd docs/parity-site && python3 -m http.server 8080
# open http://127.0.0.1:8080/
```

`--pages-history-status absent` above edits the checkout copy. It does not
read or publish GitHub Pages.

## Deploy

Three steps, kept separate on purpose:

1. **Result generation.** `nightly-equivalence.yml` scores the Ashkenazi trio.
   `giab-genomewide.yml` scores hosted GIAB modes (the usual one is
   `ci-subset`). `genomewide-validation.yml` is the self-hosted autosome path.
   Each caller passes its own workflow id (`nightly-equivalence`,
   `giab-genomewide`, or `genomewide-validation`). New history rows store that
   id. Existing rows are not relabeled.

2. **History update.** `scripts/ci/publish_parity_site.sh` fetches
   `gh-pages:data/history.json` and appends the new row. The same record id
   is not written twice, and an existing row is not replaced. The push is a
   fast-forward of that file pair onto `gh-pages`. If the push is rejected,
   the script fetches again and merges again. Workflow concurrency groups do
   not cover this branch, because those groups are per workflow file.
   The script does **not** push `docs/parity-site/data/history.json` to `main`.
   If the `gh-pages` branch is missing, or the branch has no
   `data/history.json`, the checkout file is the base. A missing checkout
   file starts an empty history. If the fetch fails or the file is not valid
   history JSON, the script stops and does not write.

3. **Site deployment.** `peaceiris/actions-gh-pages` publishes this directory
   to the `gh-pages` branch root with `keep_files: true`. Before that copy,
   the history publishers remove `data/history.json` and `data/latest.json`
   from the publish directory so the deploy cannot replace the history just
   pushed. `benchmark.yml` publishes performance JSON the same way: it removes
   those two equivalence files from its publish directory and appends
   `data/perf_*.json` via `scripts/perf/update_perf_dashboard.py`.

**Source of truth for the live equivalence history:** `data/history.json` on
the `gh-pages` branch. The copy on `main` is not that history.

Repo setting: **Settings → Pages → Deploy from branch `gh-pages` / root**.

Live URL: https://synapticfour.github.io/gatk-rs/
