"""Apply the Stage 1 and Stage 2 estimands to frozen validation overpasses."""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from analyze_built_cover import Contrast, fit_contrasts, variant_row
from analyze_landsat import load_panel, observations, write_csv
from analyze_local_comparison import DRAWS, SENSITIVITY_DRAWS, estimate_pairs, join_pair_temperatures, within_model
from audit_heldout_overpasses import SEASONS, period
from fetch_inputs import checksum
from make_validation_panel import SNAPSHOTS
from prepare_local_comparison import BLOCK, KEY, PANELS, read_design_inputs
from prepare_overpass_validation import DESIGN, DISCOVERY_YEAR, SAMPLES, verify_hashes


REGRESSION = "Stage 1 local regression"
PAIRS = "Stage 1 nearby pairs"
MODERATION = "Stage 2 contrast difference"


def local_result(frame, pairs, method, draws):
    if method == REGRESSION:
        result, _, _ = within_model(frame, draws=draws)
    else:
        if pairs.empty:
            return {"status": "no pairs"}
        result, _ = estimate_pairs(pairs, draws)
    return {k: v for k, v in result.items() if k != "coef"}


def built_result(data, target, model, draws, label):
    if data.empty:
        return {"analysis": label, "model": model, "status": "no joint support"}
    counts = data.groupby("date").band.nunique()
    eligible = data[data.date.isin(counts[counts == 2].index)]
    n = {band: len(eligible.loc[eligible.band == band, BLOCK].drop_duplicates())
         for band in ("lower", "higher")}
    if min(n.values()) < 2:
        return {"analysis": label, "model": model, "status": "fewer than two blocks in a band",
                "lower_band_blocks": n["lower"], "higher_band_blocks": n["higher"]}
    return variant_row(label, model, target, fit_contrasts(eligible, target, model=model, draws=draws))


def load_cohorts(snapshots, design, record):
    design_sha = checksum(design / "design_record.json")
    if read_design_inputs(PANELS)[2] != record["discovery_panel_sha256"]:
        raise ValueError("The discovery panel differs from the audit input")
    parts = [observations(load_panel(PANELS, "world_500"))]
    for year in record["validation_years"]:
        base = snapshots / str(year)
        p = json.loads((base / "provenance.json").read_text())
        if p["design_sha256"] != design_sha:
            raise ValueError("The validation snapshot and design do not match")
        verify_hashes(base, {"world_500.csv.gz": p["panel_sha256"]})
        dates = [d for d in record["dates"] if d.startswith(f"{year}-")]
        parts.append(observations(load_panel(base, "world_500", dates)))
    seen = pd.concat(parts, ignore_index=True)
    registry = pd.read_csv(design / "eligible_cell_dates.csv")
    rows = seen.merge(registry, on=[*KEY, "scene_id", "valid_lst_centers"], validate="one_to_one")
    if len(rows) != len(registry) or len(rows) != len(seen):
        raise ValueError("The observed samples differ from the frozen cell-date registry")
    year = rows.date.str[:4].astype(int)
    rows["year"] = year
    return {"discovery": rows[year == DISCOVERY_YEAR], "validation": rows[year != DISCOVERY_YEAR]}


def date_record(date, cohort, method, result, estimate, lo, hi):
    return {
        "cohort": cohort, "date": date, "period": period(date), "analysis": method, "status": result["status"],
        "estimate_c": result.get(estimate), "lower_c": result.get(lo), "upper_c": result.get(hi),
        **{k: result.get(k) for k in ("cell_dates", "pairs", "bootstrap_blocks", "lower_band_blocks",
                                    "higher_band_blocks", "draws_retained")},
    }


def date_and_period_tables(by_date):
    periods, balance = [], []
    for (cohort, analysis), group in by_date.groupby(["cohort", "analysis"]):
        means = []
        for label in SEASONS:
            g = group[group.period == label]
            values = g.estimate_c.dropna()
            estimate = float(values.mean()) if len(values) else None
            if estimate is not None:
                means.append(estimate)
            periods.append({
                "cohort": cohort, "analysis": analysis, "period": label, "selected_dates": len(g),
                "estimated_dates": len(values), "mean_c": estimate,
                "min_c": values.min(), "max_c": values.max(),
            })
        if means:
            balance.append({
                "cohort": cohort, "analysis": analysis, "periods_represented": len(means),
                "estimated_dates": int(group.estimate_c.notna().sum()),
                "equal_period_mean_c": float(np.mean(means)),
            })
    return pd.DataFrame(periods), pd.DataFrame(balance)


def change_table(annual):
    common = annual[annual.common]
    counts = []
    for mixed, subset in common.groupby("boundary_mixed", sort=True):
        counts.append({
            "boundary_mixed": bool(mixed), "common_cells": len(subset),
            "stable_cells": int(subset.common_stable.sum()),
            **{f"{kind}_span_over_5pp_cells": int(subset[f"{column}_span_pp"].gt(5 + 1e-6).sum())
               for kind, column in (("tree", "tree_fraction"), ("built", "built_fraction"),
                                    ("water", "esri_water_share"))},
            "water_span_over_10pp_cells": int(subset.esri_water_share_span_pp.gt(10).sum()),
        })
    return pd.DataFrame(counts)


def run(design, snapshots, output):
    record = json.loads((design / "design_record.json").read_text())
    verify_hashes(design, record["design_sha256"])
    target = Contrast(**{k: record["stage2_contrast"][k] for k in Contrast.__dataclass_fields__})
    cohorts = load_cohorts(snapshots, design, record)
    pair_ids = pd.read_csv(design / "stage1_pairs.csv", float_precision="round_trip")
    built_ids = pd.read_csv(design / "stage2_fit_sample.csv")
    local, built, dates = [], [], []
    all_pairs = []
    for cohort, observations_ in cohorts.items():
        for sample in SAMPLES:
            frame = observations_ if sample == "full" else observations_[observations_[sample]]
            ids = pair_ids[pair_ids.date.isin(frame.date.unique()) & pair_ids["sample"].eq(sample)]
            fits = built_ids[built_ids.date.isin(frame.date.unique()) & built_ids["sample"].eq(sample)]
            pairs = join_pair_temperatures(ids, frame)
            stage2 = frame.merge(fits, on=[*KEY, "year"], how="inner", validate="one_to_one")
            if len(stage2) != len(fits):
                raise ValueError("A Stage 2 support sample is missing its temperature")
            all_pairs.append(pairs.assign(cohort=cohort))
            tags = {"cohort": cohort, "sample": sample}
            for method in (REGRESSION, PAIRS):
                local.append({**tags, "analysis": method, **local_result(frame, pairs, method, DRAWS)})
            for data, model, label in (
                (stage2, "quadratic", "Primary quadratic interaction"),
                (stage2, "piecewise", "Piecewise linear interaction"),
                (stage2[stage2.shared], "quadratic", "Shared block-dates"),
            ):
                built.append({**tags, **built_result(data, target, model, DRAWS, label)})

            if sample == "full":
                for date, day in frame.groupby("date", sort=True):
                    for method in (REGRESSION, PAIRS):
                        result = local_result(day, pairs[pairs.date == date], method, SENSITIVITY_DRAWS)
                        dates.append(date_record(date, cohort, method, result, "estimate_c_per_10pp", "lower_c", "upper_c"))
                    r = built_result(stage2[stage2.date == date], target, "quadratic", SENSITIVITY_DRAWS, MODERATION)
                    dates.append(date_record(date, cohort, MODERATION, r, "higher_minus_lower_estimate_c",
                                             "higher_minus_lower_lower_c", "higher_minus_lower_upper_c"))

    by_date = pd.DataFrame(dates)
    periods, balanced = date_and_period_tables(by_date)
    annual = pd.read_csv(design / "annual_cover_comparison.csv")
    output.mkdir(parents=True, exist_ok=True)
    tables = {
        "local_estimates.csv": pd.DataFrame(local),
        "built_cover_estimates.csv": pd.DataFrame(built),
        "date_estimates.csv": by_date,
        "period_estimates.csv": periods,
        "equal_period_estimates.csv": balanced,
        "cover_stability.csv": change_table(annual),
        "paired_temperatures.csv": pd.concat(all_pairs, ignore_index=True),
    }
    for filename, table in tables.items():
        write_csv(table, output / filename)
    info = {
        "design_sha256": checksum(design / "design_record.json"), "validation_years": record["validation_years"],
        "dates": record["dates"], "common_cells": int(annual.common.sum()),
        "common_stable_cells": int(annual.common_stable.sum()),
        "result_sha256": {name: checksum(output / name) for name in tables},
    }
    (output / "analysis_record.json").write_text(json.dumps(info, indent=2) + "\n")
    print(pd.DataFrame(local)[["cohort", "sample", "analysis", "estimate_c_per_10pp", "lower_c", "upper_c"]].to_string(index=False))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--design", type=Path, default=DESIGN)
    parser.add_argument("--snapshots", type=Path, default=SNAPSHOTS)
    parser.add_argument("--output", type=Path, required=True)
    opts = parser.parse_args()
    run(opts.design, opts.snapshots, opts.output)
