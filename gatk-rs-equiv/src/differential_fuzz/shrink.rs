//! Shrink a diverging [`Scenario`] while the divergence still reproduces.

use super::runner::{evaluate_scenario, Divergence, EvalBins};
use super::scenario::Scenario;
use anyhow::Result;
use std::path::Path;

/// Try to minimize `scenario` while `pred` still reports a divergence.
pub fn shrink_scenario(
    mut scenario: Scenario,
    bins: &EvalBins,
    work_root: &Path,
    max_steps: usize,
) -> Result<(Scenario, Divergence)> {
    let mut last = evaluate_scenario(&scenario, bins, &work_root.join("shrink_0"))?
        .expect("shrink called without an initial divergence");

    let mut steps = 0usize;

    // Axis 1: fewer reads
    while steps < max_steps && scenario.n_reads > 4 {
        let mut trial = scenario.clone();
        trial.n_reads = (trial.n_reads - 2).max(4);
        steps += 1;
        if let Some(div) = evaluate_scenario(&trial, bins, &work_root.join(format!("s{steps}")))? {
            scenario = trial;
            last = div;
        } else {
            break;
        }
    }

    // Axis 2: shorter reference / reads
    while steps < max_steps && scenario.ref_len > 140 {
        let mut trial = scenario.clone();
        trial.ref_len = (trial.ref_len - 20).max(140);
        if trial.read_len >= trial.ref_len {
            trial.read_len = trial.ref_len / 2;
        }
        steps += 1;
        if let Some(div) = evaluate_scenario(&trial, bins, &work_root.join(format!("s{steps}")))? {
            scenario = trial;
            last = div;
        } else {
            break;
        }
    }

    // Axis 3: disable soft-clips / mates / plants one at a time
    for (label, apply) in [
        (
            "no_softclip",
            (|s: &mut Scenario| s.softclip_p = 0) as fn(&mut Scenario),
        ),
        ("no_mates", |s: &mut Scenario| s.mate_overlap_p = 0),
        ("no_plant_indel", |s: &mut Scenario| s.plant_indel = false),
        ("no_indels", |s: &mut Scenario| s.indel_p = 0),
        ("high_mapq", |s: &mut Scenario| {
            s.mapq_min = 60;
            s.mapq_span = 1;
        }),
    ] {
        if steps >= max_steps {
            break;
        }
        let mut trial = scenario.clone();
        apply(&mut trial);
        steps += 1;
        let dir = work_root.join(format!("s{steps}_{label}"));
        if let Some(div) = evaluate_scenario(&trial, bins, &dir)? {
            scenario = trial;
            last = div;
        }
    }

    Ok((scenario, last))
}
