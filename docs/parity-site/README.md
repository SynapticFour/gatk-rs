# Public equivalence & performance dashboard (GitHub Pages)

Static site (HTML/CSS/JS + Chart.js CDN). No build step.

Tabs: **Equivalence** (hap.py/RTG) · **Performance** (fair HC timings vs Java
`FASTEST_AVAILABLE` on the dedicated quiet host).

| File | Role |
|------|------|
| `index.html` / `style.css` / `app.js` | UI (two tabs) |
| `data/history.json` | Equivalence run history (nightly / genomewide) |
| `data/latest.json` | Most recent equivalence snapshot |
| `data/perf_history.json` | Fair HC comparison history (`benchmark.yml`) |
| `data/perf_latest.json` | Most recent performance snapshot |
## Update locally

```bash
python3 scripts/parity/giab/update_public_dashboard.py \
  --source nightly \
  --json path/to/happy_summary.json \
  --site-dir docs/parity-site \
  --commit-sha "$(git rev-parse HEAD)"

# Preview
cd docs/parity-site && python3 -m http.server 8080
# open http://127.0.0.1:8080/
```

## Deploy

`nightly-equivalence.yml` and `genomewide-validation.yml` commit equivalence
`data/*.json`. `benchmark.yml` appends `data/perf_*.json` via
`scripts/perf/update_perf_dashboard.py`. Both publish this directory to the
`gh-pages` branch root via `peaceiris/actions-gh-pages` (`keep_files: true`).

Repo setting: **Settings → Pages → Deploy from branch `gh-pages` / root**.

Live URL (org `gatk-rs`): https://gatk-rs.github.io/gatk-rs/
