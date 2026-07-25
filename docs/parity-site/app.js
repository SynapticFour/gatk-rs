/* global Chart */
(function () {
  "use strict";

  const state = {
    history: null,
    perfHistory: null,
    trendChart: null,
    barChart: null,
    speedupChart: null,
    wallChart: null,
    tab: "equivalence",
  };

  function $(id) {
    return document.getElementById(id);
  }

  function fmt(n) {
    if (n === null || n === undefined || Number.isNaN(n)) return "—";
    return Number(n).toFixed(4);
  }

  function fmtSec(n) {
    if (n === null || n === undefined || Number.isNaN(n)) return "—";
    if (n < 1) return (n * 1000).toFixed(1) + " ms";
    return Number(n).toFixed(3) + " s";
  }

  function latestRun(hist) {
    if (!hist || !hist.runs || !hist.runs.length) return null;
    return hist.runs[hist.runs.length - 1];
  }

  function escapeHtml(s) {
    return String(s)
      .replace(/&/g, "&amp;")
      .replace(/</g, "&lt;")
      .replace(/>/g, "&gt;")
      .replace(/"/g, "&quot;");
  }

  function palette(i) {
    const colors = [
      "#3dbf9c",
      "#5ec4ff",
      "#e0a35c",
      "#c792ea",
      "#82aaff",
      "#ff8b7e",
      "#c3e88d",
      "#f78c6c",
    ];
    return colors[i % colors.length];
  }

  /* ---------- tabs ---------- */
  function setTab(name) {
    state.tab = name;
    document.querySelectorAll(".tab").forEach((btn) => {
      const on = btn.dataset.tab === name;
      btn.classList.toggle("active", on);
      btn.setAttribute("aria-selected", on ? "true" : "false");
    });
    $("panel-equivalence").hidden = name !== "equivalence";
    $("panel-performance").hidden = name !== "performance";
    if (name === "performance") refreshPerf();
  }

  /* ---------- equivalence (unchanged logic) ---------- */
  function renderScope(run, meta) {
    const empty = $("scope-empty");
    const dl = $("scope-dl");
    const honesty = $("scope-honesty");
    if (!run) {
      empty.hidden = false;
      dl.hidden = true;
      honesty.hidden = true;
      return;
    }
    empty.hidden = true;
    dl.hidden = false;
    honesty.hidden = false;

    const scope = run.scope || {};
    const rows = [
      ["Workflow", run.workflow || "—"],
      ["Generated (UTC)", run.generated_utc || "—"],
      ["Commit", run.commit_sha ? run.commit_sha.slice(0, 12) : "—"],
      ["Java GATK (reference)", scope.java_gatk_version || meta.java_gatk_version || "4.4.0.0"],
      ["Java GATK Docker", scope.java_gatk_docker || meta.java_gatk_docker || "—"],
      ["Java GATK SHA", (scope.java_gatk_sha || meta.java_gatk_sha || "—").toString().slice(0, 12)],
      ["Samples", (scope.samples || []).join(", ") || "—"],
      [
        "Cohort sizes (N)",
        (scope.cohort_sizes || []).length
          ? (scope.cohort_sizes || []).join(", ")
          : "—",
      ],
      [
        "Recommended max N",
        scope.recommended_max_samples != null
          ? String(scope.recommended_max_samples)
          : "—",
      ],
      ["Regions / intervals", (scope.regions || []).join(", ") || scope.mode_description || "—"],
      ["Assembly", scope.assembly || "—"],
      ["Truth set", scope.truth || "—"],
      ["Pipeline", scope.pipeline || "—"],
      ["Engine", scope.eval_engine || "hap.py / RTG via gatk-rs-equiv"],
    ];

    dl.innerHTML = rows
      .map(
        ([k, v]) =>
          `<dt>${escapeHtml(k)}</dt><dd>${escapeHtml(String(v))}</dd>`
      )
      .join("");

    honesty.textContent =
      scope.honesty ||
      "These metrics apply only to the samples and genomic intervals listed above. " +
        "They do not certify genome-wide clinical equivalence outside that scope.";
  }

  function seriesPoints(hist, metric, vtype, engine) {
    const byRegion = new Map();
    for (const run of hist.runs || []) {
      const t = run.generated_utc || run.id || "";
      for (const row of run.metrics || []) {
        if (row.variant_type !== vtype) continue;
        if ((row.engine || "rust") !== engine) continue;
        // Cohort-scale runs: primary series key is cohort_size (N), not region.
        const region =
          row.cohort_size != null
            ? `N=${row.cohort_size}`
            : row.region || row.sample || "all";
        if (!byRegion.has(region)) byRegion.set(region, []);
        const val = row[metric];
        if (val === null || val === undefined) continue;
        byRegion.get(region).push({ t, y: Number(val), run });
      }
    }
    return byRegion;
  }

  function renderTrend(hist) {
    const metric = $("metric").value;
    const vtype = $("vtype").value;
    const engine = $("engine").value;
    const byRegion = seriesPoints(hist, metric, vtype, engine);
    const labels = (hist.runs || []).map((r) =>
      (r.generated_utc || "").replace("T", " ").replace("Z", "")
    );

    const datasets = [];
    let i = 0;
    for (const [region, pts] of byRegion) {
      const byTime = new Map(pts.map((p) => [p.t, p.y]));
      datasets.push({
        label: region,
        data: (hist.runs || []).map((r) => {
          const t = r.generated_utc || r.id || "";
          return byTime.has(t) ? byTime.get(t) : null;
        }),
        borderColor: palette(i),
        backgroundColor: palette(i),
        spanGaps: true,
        tension: 0.2,
        pointRadius: 3,
      });
      i += 1;
    }

    const ctx = $("chart-trend").getContext("2d");
    if (state.trendChart) state.trendChart.destroy();
    const yScale = { min: 0, ticks: { color: "#9aabbd" }, grid: { color: "#314050" } };
    // F1/precision/recall are unit interval; wall/RSS are unbounded.
    if (metric === "f1" || metric === "precision" || metric === "recall") {
      yScale.max = 1;
    }
    state.trendChart = new Chart(ctx, {
      type: "line",
      data: { labels, datasets },
      options: {
        responsive: true,
        maintainAspectRatio: false,
        scales: {
          y: yScale,
          x: {
            ticks: { color: "#9aabbd", maxRotation: 45 },
            grid: { color: "#1a222c" },
          },
        },
        plugins: {
          legend: { labels: { color: "#e7eef6", boxWidth: 12 } },
        },
      },
    });
  }

  function renderBars(hist) {
    const metric = $("metric").value;
    const vtype = $("vtype").value;
    const engine = $("engine").value;
    const run = latestRun(hist);
    const labels = [];
    const values = [];
    if (run) {
      const rows = (run.metrics || []).filter(
        (row) =>
          row.variant_type === vtype && (row.engine || "rust") === engine
      );
      // Sort cohort ladder by N when present.
      rows.sort((a, b) => (a.cohort_size || 0) - (b.cohort_size || 0));
      for (const row of rows) {
        labels.push(
          row.cohort_size != null
            ? `N=${row.cohort_size}`
            : row.region || row.sample || "all"
        );
        values.push(row[metric]);
      }
    }
    const ctx = $("chart-bars").getContext("2d");
    if (state.barChart) state.barChart.destroy();
    const yScale = { min: 0, ticks: { color: "#9aabbd" }, grid: { color: "#314050" } };
    if (metric === "f1" || metric === "precision" || metric === "recall") {
      yScale.max = 1;
    }
    state.barChart = new Chart(ctx, {
      type: "bar",
      data: {
        labels,
        datasets: [
          {
            label: `${engine} ${vtype} ${metric}`,
            data: values,
            backgroundColor: engine === "rust" ? "#5ec4ff99" : "#e0a35c99",
            borderColor: engine === "rust" ? "#5ec4ff" : "#e0a35c",
            borderWidth: 1,
          },
        ],
      },
      options: {
        responsive: true,
        maintainAspectRatio: false,
        scales: {
          y: yScale,
          x: {
            ticks: { color: "#9aabbd" },
            grid: { display: false },
          },
        },
        plugins: {
          legend: { display: false },
        },
      },
    });
  }

  function renderTable(hist) {
    const run = latestRun(hist);
    const wrap = $("table-wrap");
    if (!run || !(run.metrics || []).length) {
      wrap.innerHTML = '<p class="empty">No metrics in latest run.</p>';
      return;
    }
    const isCohort = (run.scope || {}).kind === "cohort_joint_scale";
    const rows = run.metrics
      .map((m) =>
        isCohort
          ? `<tr>
        <td>${escapeHtml(m.cohort_size != null ? `N=${m.cohort_size}` : m.region || "—")}</td>
        <td>${escapeHtml(m.engine || "—")}</td>
        <td>${escapeHtml(m.variant_type || "—")}</td>
        <td class="num">${fmt(m.wall_sec)}</td>
        <td class="num">${fmt(m.peak_rss_kb)}</td>
        <td class="num">${fmt(m.f1)}</td>
      </tr>`
          : `<tr>
        <td>${escapeHtml(m.region || m.sample || "—")}</td>
        <td>${escapeHtml(m.engine || "—")}</td>
        <td>${escapeHtml(m.variant_type || "—")}</td>
        <td class="num">${fmt(m.precision)}</td>
        <td class="num">${fmt(m.recall)}</td>
        <td class="num">${fmt(m.f1)}</td>
      </tr>`
      )
      .join("");
    wrap.innerHTML = isCohort
      ? `<table>
      <thead><tr>
        <th>Cohort size</th><th>Engine</th><th>Type</th>
        <th>Wall (s)</th><th>Peak RSS (KiB)</th><th>F1 (vs Java)</th>
      </tr></thead>
      <tbody>${rows}</tbody>
    </table>`
      : `<table>
      <thead><tr>
        <th>Region / sample</th><th>Engine</th><th>Type</th>
        <th>Precision</th><th>Recall</th><th>F1</th>
      </tr></thead>
      <tbody>${rows}</tbody>
    </table>`;
  }

  function refreshEq() {
    const hist = state.history;
    if (!hist) return;
    $("run-count").textContent = `${(hist.runs || []).length} run(s)`;
    if (hist.meta && hist.meta.updated_utc) {
      $("updated").textContent = ` · Updated ${hist.meta.updated_utc}`;
    }
    renderScope(latestRun(hist), hist.meta || {});
    renderTrend(hist);
    renderBars(hist);
    renderTable(hist);
  }

  /* ---------- performance ---------- */
  function refreshPerf() {
    const hist = state.perfHistory;
    const empty = $("perf-empty");
    if (!hist || !(hist.runs || []).length) {
      empty.hidden = false;
      $("perf-controls").hidden = true;
      $("perf-charts").hidden = true;
      $("perf-table-card").hidden = true;
      return;
    }
    empty.hidden = true;
    $("perf-controls").hidden = false;
    $("perf-charts").hidden = false;
    $("perf-table-card").hidden = false;
    $("perf-run-count").textContent = `${hist.runs.length} run(s)`;
    if (hist.meta && hist.meta.updated_utc) {
      $("updated").textContent = ` · Updated ${hist.meta.updated_utc}`;
    }
    renderSpeedup(hist);
    renderWallLatest(hist);
    renderPerfTable(hist);
  }

  function renderSpeedup(hist) {
    const regionFilter = $("perf-region").value;
    const labels = (hist.runs || []).map((r) => {
      const t = (r.generated_utc || "").replace("T", " ").replace("Z", "");
      const sha = (r.commit_sha || "").slice(0, 7);
      return sha ? `${t} (${sha})` : t;
    });
    const regions = ["small", "medium", "large"];
    const datasets = [];
    let i = 0;
    for (const region of regions) {
      if (regionFilter !== "all" && regionFilter !== region) continue;
      datasets.push({
        label: `${region} · vs Java FASTEST_AVAILABLE`,
        data: (hist.runs || []).map((r) => {
          const hit = (r.speedups || []).find(
            (s) =>
              s.region_size === region &&
              s.java_baseline_pair_hmm === "FASTEST_AVAILABLE"
          );
          return hit ? hit.speedup : null;
        }),
        borderColor: palette(i),
        backgroundColor: palette(i),
        spanGaps: true,
        tension: 0.2,
        pointRadius: 3,
      });
      i += 1;
    }

    const ctx = $("chart-speedup").getContext("2d");
    if (state.speedupChart) state.speedupChart.destroy();
    state.speedupChart = new Chart(ctx, {
      type: "line",
      data: { labels, datasets },
      options: {
        responsive: true,
        maintainAspectRatio: false,
        scales: {
          y: {
            title: {
              display: true,
              text: "Speedup (×)",
              color: "#9aabbd",
            },
            ticks: { color: "#9aabbd" },
            grid: { color: "#314050" },
            beginAtZero: true,
          },
          x: {
            ticks: { color: "#9aabbd", maxRotation: 45 },
            grid: { color: "#1a222c" },
          },
        },
        plugins: {
          legend: { labels: { color: "#e7eef6", boxWidth: 12 } },
          title: {
            display: true,
            text: "Java baseline = FASTEST_AVAILABLE (native AVX when loaded)",
            color: "#e0a35c",
            font: { size: 11 },
          },
        },
      },
    });
  }

  function renderWallLatest(hist) {
    const run = latestRun(hist);
    const regionFilter = $("perf-region").value;
    const labels = [];
    const values = [];
    const colors = [];
    if (run) {
      for (const c of run.cells || []) {
        if (regionFilter !== "all" && c.region_size !== regionFilter) continue;
        labels.push(`${c.region_size}/${c.config_id}`);
        values.push(c.wall_median_sec);
        colors.push(
          c.config_id && c.config_id.startsWith("rust") ? "#5ec4ff99" : "#e0a35c99"
        );
      }
    }
    const ctx = $("chart-wall").getContext("2d");
    if (state.wallChart) state.wallChart.destroy();
    state.wallChart = new Chart(ctx, {
      type: "bar",
      data: {
        labels,
        datasets: [
          {
            label: "median wall (s)",
            data: values,
            backgroundColor: colors,
            borderWidth: 1,
          },
        ],
      },
      options: {
        responsive: true,
        maintainAspectRatio: false,
        scales: {
          y: {
            title: { display: true, text: "seconds", color: "#9aabbd" },
            ticks: { color: "#9aabbd" },
            grid: { color: "#314050" },
            beginAtZero: true,
          },
          x: {
            ticks: { color: "#9aabbd", maxRotation: 60, font: { size: 9 } },
            grid: { display: false },
          },
        },
        plugins: { legend: { display: false } },
      },
    });
  }

  function renderPerfTable(hist) {
    const run = latestRun(hist);
    const wrap = $("perf-table-wrap");
    if (!run || !(run.cells || []).length) {
      wrap.innerHTML = '<p class="empty">No cells in latest run.</p>';
      return;
    }
    const host = run.host || {};
    const hostLine = host.cpu_model
      ? `<p class="host-line">Host: <code>${escapeHtml(host.cpu_model)}</code> · ${escapeHtml(
          String(host.logical_cpus || "?")
        )} CPUs · ${escapeHtml(String(host.mem_gib || "?"))} GiB · kernel <code>${escapeHtml(
          host.kernel || "?"
        )}</code></p>`
      : "";
    const baseline = `<p class="host-line">Primary speedup baseline: <strong>Java GATK <code>FASTEST_AVAILABLE</code></strong>${
      run.workflow_run_url
        ? ` · <a href="${escapeHtml(run.workflow_run_url)}">workflow run</a>`
        : ""
    }</p>`;
    const rows = run.cells
      .map(
        (c) => `<tr>
        <td>${escapeHtml(c.region_size)}</td>
        <td><code>${escapeHtml(c.config_id)}</code></td>
        <td><code>${escapeHtml(c.pair_hmm || "—")}</code></td>
        <td class="num">${fmtSec(c.wall_median_sec)}</td>
        <td class="num">${c.wall_stdev_sec != null ? fmtSec(c.wall_stdev_sec) : "—"}</td>
        <td class="num">${fmtSec(c.user_median_sec)}</td>
        <td class="num">${fmtSec(c.sys_median_sec)}</td>
        <td class="num">${
          c.peak_rss_kb_median != null
            ? (c.peak_rss_kb_median / 1024).toFixed(1) + " MiB"
            : "—"
        }</td>
        <td class="num">${
          c.energy_joules_median != null
            ? Number(c.energy_joules_median).toFixed(2) + " J"
            : "n/a"
        }</td>
      </tr>`
      )
      .join("");
    wrap.innerHTML =
      hostLine +
      baseline +
      `<table>
      <thead><tr>
        <th>Region</th><th>Config</th><th>PairHMM</th>
        <th>Wall median</th><th>Wall stdev</th><th>User</th><th>Sys</th>
        <th>Peak RSS</th><th>Energy</th>
      </tr></thead>
      <tbody>${rows}</tbody>
    </table>`;
  }

  async function main() {
    document.querySelectorAll(".tab").forEach((btn) => {
      btn.addEventListener("click", () => setTab(btn.dataset.tab));
    });

    try {
      const res = await fetch("data/history.json", { cache: "no-store" });
      if (!res.ok) throw new Error(`HTTP ${res.status}`);
      state.history = await res.json();
    } catch (err) {
      $("scope-empty").textContent =
        "Could not load data/history.json (" + err.message + ").";
    }

    try {
      const res = await fetch("data/perf_history.json", { cache: "no-store" });
      if (res.ok) state.perfHistory = await res.json();
    } catch (_) {
      /* optional until first bench */
    }

    ["metric", "vtype", "engine"].forEach((id) => {
      $(id).addEventListener("change", refreshEq);
    });
    $("perf-region").addEventListener("change", refreshPerf);

    refreshEq();
    if (state.tab === "performance") refreshPerf();
  }

  main();
})();
