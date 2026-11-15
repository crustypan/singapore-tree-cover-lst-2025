"""Estimate Stage 1 within-area associations and summarize spatial uncertainty."""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from analyze_landsat import block_coordinates, interval, load_panel, moran, observations, write_csv
from audit_coverage import ESRI_NAMES
from fetch_inputs import checksum
from make_landsat_panel import DATES, VALID_ESRI
from prepare_local_comparison import (
    BLOCK, COORDS, COVARIATES, DESIGN, KEY, PANELS, RULES, MatchRule, assign_blocks, read_design_inputs,
)


SEED = 20251116
DRAWS = 1999
SENSITIVITY_DRAWS = 999
TREE_CENTER = 0.25


def draw_totals(stats, cluster, dates, draws):
    """Resample whole geographic clusters, holding all their date records together."""
    _, spatial_id = np.unique(cluster, axis=0, return_inverse=True)
    _, date_id = np.unique(dates, return_inverse=True)
    n_blocks = int(spatial_id.max()) + 1
    n_dates = int(date_id.max()) + 1
    moments = np.zeros((n_blocks, n_dates, 1 + stats.shape[1]))
    np.add.at(moments[:, :, 0], (spatial_id, date_id), 1)
    np.add.at(moments[:, :, 1:], (spatial_id, date_id), stats)
    rng = np.random.default_rng(SEED)
    multiplicity = rng.multinomial(n_blocks, np.full(n_blocks, 1 / n_blocks), size=draws)
    flat = multiplicity @ moments.reshape(n_blocks, -1)
    samples = flat.reshape(draws, n_dates, moments.shape[2])
    samples = samples[(samples[:, :, 0] > 0).all(axis=1)]
    return np.concatenate([moments.sum(axis=0)[None, :, :], samples]), n_blocks


def within_model(data, model_size=4000, origin=0, cluster_size=None, equal_dates=True,
                 curved_tree=False, draws=DRAWS):
    if cluster_size is None:
        cluster_size = model_size
    if cluster_size % model_size:
        raise ValueError("Bootstrap blocks must contain entire fixed-effect blocks")
    data = assign_blocks(data, MatchRule(block_size=model_size, block_shift=origin))
    group_key = ["date", *BLOCK]
    group_sizes = data.groupby(group_key).cell_id.transform("size")
    records = data[group_sizes >= 2].reset_index(drop=True)
    if records.empty:
        raise ValueError("No populated date-by-block groups")
    tree = (records.tree_fraction - TREE_CENTER) * 10
    built = (records.built_fraction - 0.5) * 10
    predictors = pd.DataFrame({
        "tree": tree, "built": built, "built_squared": built ** 2, "land_share": records.land_share,
    })
    if curved_tree:
        predictors["tree_squared"] = tree ** 2
    design = pd.concat([records[group_key], predictors, records[["response_c"]]], axis=1)
    centered = design[[*predictors, "response_c"]] - design.groupby(group_key)[
        [*predictors, "response_c"]
    ].transform("mean")
    x = centered[predictors.columns].to_numpy()
    y = centered.response_c.to_numpy()
    p = x.shape[1]
    products = np.einsum("ni,nj->nij", x, x).reshape(len(x), p ** 2)
    moments, blocks = draw_totals(np.column_stack([products, x * y[:, None]]),
                                 block_coordinates(records, cluster_size, origin), records.date, draws)
    totals = moments[:, :, 1:]
    if equal_dates:
        totals = totals / moments[:, :, :1]
    totals = totals.sum(axis=1)
    gram = totals[:, :p ** 2].reshape(-1, p, p)
    rhs = totals[:, p ** 2:]
    eigen = np.linalg.eigvalsh(gram)
    ok = eigen[:, 0] > np.finfo(float).eps * p * eigen[:, -1]
    if not ok[0]:
        raise ValueError("Local regression is rank deficient")
    coefficients = np.linalg.solve(gram[ok], rhs[ok, :, None]).squeeze(axis=2)
    boot = coefficients[1:]
    if len(boot) < 0.95 * draws:
        raise ValueError(f"Only {len(boot)} of {draws} spatial draws were estimable")
    coef = coefficients[0]
    weights = np.ones(len(records))
    if equal_dates:
        weights /= records.groupby("date").cell_id.transform("size").to_numpy()
    residual = y - x @ coef
    records["fitted_c"] = records.response_c - residual
    records["residual_c"] = residual
    mse = np.average(residual ** 2, weights=weights)
    raw_within = np.average(y ** 2, weights=weights)
    stats = {
        "estimate_c_per_10pp": float(coef[0]), "lower_c": None, "upper_c": None,
        "cell_dates": len(records), "distinct_cells": int(records.cell_id.nunique()),
        "dates": int(records.date.nunique()), "bootstrap_blocks": blocks,
        "draws_requested": draws, "draws_retained": len(boot), "status": "estimated",
        "within_r2": float(1 - mse / raw_within), "rmse_c": float(np.sqrt(mse)),
        "date_by_block_groups": int(len(records[group_key].drop_duplicates())),
        "coef": dict(zip(predictors.columns, map(float, coef))),
    }
    if draws:
        stats["lower_c"], stats["upper_c"] = map(float, interval(boot[:, 0]))
    return stats, records, boot


def join_pair_temperatures(pairs, data):
    locations = data[[*KEY, *COORDS, *COVARIATES, "response_c"]]
    results = pairs.copy()
    for side in ("lower", "higher"):
        key = f"{side}_tree_cell"
        lookup = locations.rename(columns={
            c: key if c == "cell_id" else f"{side}_{c}" for c in locations.columns if c != "date"
        })
        results = results.merge(lookup, on=["date", key], how="left", validate="one_to_one")
    if results[["lower_response_c", "higher_response_c"]].isna().any().any():
        raise ValueError("A frozen pair is missing one or both eligible temperatures")
    gap = (results.higher_tree_fraction - results.lower_tree_fraction) * 100
    np.testing.assert_allclose(results.tree_gap_pp, gap, atol=1e-8, rtol=0)
    results["lst_gap_c"] = results.higher_response_c - results.lower_response_c
    results["lst_per_10pp_c"] = 10 * results.lst_gap_c / gap
    return results


def estimate_pairs(data, draws=DRAWS):
    metrics = ["lst_per_10pp_c", "lst_gap_c"]
    totals, n = draw_totals(data[metrics].to_numpy(), data[BLOCK].to_numpy(), data.date, draws)
    values = (totals[:, :, 1:] / totals[:, :, :1]).mean(axis=1)
    if len(values) - 1 < 0.95 * draws:
        raise ValueError("Too many spatial draws omit an entire study date")
    ci = interval(values[1:]) if draws else [[None, None], [None, None]]
    endpoints = pd.concat([data.lower_tree_cell, data.higher_tree_cell], ignore_index=True)
    return {
        "estimate_c_per_10pp": float(values[0, 0]), "lower_c": ci[0][0], "upper_c": ci[1][0],
        "raw_difference_c": float(values[0, 1]), "raw_lower_c": ci[0][1], "raw_upper_c": ci[1][1],
        "cell_dates": 2 * len(data), "distinct_cells": int(endpoints.nunique()),
        "dates": int(data.date.nunique()), "bootstrap_blocks": n, "pairs": len(data),
        "draws_requested": draws, "draws_retained": len(values) - 1, "status": "estimated",
    }


def matched_endpoints(primary, pairs):
    ids = pd.concat([
        pairs[["date", name]].rename(columns={name: "cell_id"})
        for name in ("lower_tree_cell", "higher_tree_cell")
    ])
    return primary.merge(ids, how="inner", on=KEY, validate="one_to_one")


def result_row(label, stats):
    return {"analysis": label, **{k: v for k, v in stats.items() if k != "coef"}}


def design_diagnostics(design, primary, pairs, output):
    composition = pd.read_csv(design / "land_composition.csv").set_index("cell_id")
    metrics = pairs[["date", *BLOCK, "distance_m", "tree_gap_pp", "built_gap_pp", "land_gap_pp"]].copy()
    for field in ("built_gap_pp", "land_gap_pp"):
        metrics["absolute_" + field] = metrics[field].abs()
    totals = metrics.groupby("date").mean(numeric_only=True).drop(columns=BLOCK).mean()
    daily = metrics.groupby("date").mean(numeric_only=True).drop(columns=BLOCK)
    daily.loc["equal_date_mean"] = totals
    write_csv(daily.reset_index(), output / "pair_balance.csv")
    rows = []
    for code in VALID_ESRI:
        fraction = composition[f"class_{code:02d}_centers"].div(composition.classified_centers)
        column = pd.DataFrame({
            "date": pairs.date,
            "lower": pairs.lower_tree_cell.map(fraction) * 100,
            "higher": pairs.higher_tree_cell.map(fraction) * 100,
        })
        mean = column.groupby("date").mean().mean()
        rows.append({
            "class": ESRI_NAMES[code], "lower_tree_cell_share_pct": mean.lower,
            "higher_tree_cell_share_pct": mean.higher, "higher_minus_lower_pp": mean.higher - mean.lower,
        })
    write_csv(pd.DataFrame(rows), output / "pair_land_composition.csv")
    design_cells = pd.read_csv(design / "cell_design.csv")
    usable = design_cells[design_cells.usable_dates > 0]
    endpoints = design_cells[design_cells.paired_dates > 0]
    a = assign_blocks(primary)
    result = []
    for label, sample in [("All usable", usable), ("Has a candidate", usable[usable.candidate_dates > 0]),
                          ("Disjoint pairs", endpoints)]:
        blocks = assign_blocks(sample)[BLOCK].drop_duplicates()
        result.append({
            "sample": label, "distinct_cells": len(sample), "area_of_cells_km2": sample.land_centers.sum() / 10000,
            "spatial_blocks": len(blocks), "mean_tree_pct": 100 * sample.tree_fraction.mean(),
            "mean_built_pct": 100 * sample.built_fraction.mean(), "mean_usable_dates": sample.usable_dates.mean(),
            "mean_land_share_pct": 100 * sample.land_share.mean(),
        })
    write_csv(pd.DataFrame(result), output / "sample_comparison.csv")
    support = a.groupby(["date", *BLOCK]).agg(
        cells=("cell_id", "size"), minimum_tree=("tree_fraction", "min"), maximum_tree=("tree_fraction", "max"),
    ).reset_index()
    write_csv(support, output / "within_block_support.csv")


def run(panels, design, output):
    output.mkdir(parents=True, exist_ok=True)
    design_record = json.loads((design / "design_record.json").read_text())
    for name, digest in design_record["design_sha256"].items():
        if checksum(design / name) != digest:
            raise ValueError(f"Frozen design changed for {name}")
    _, eligible, approved = read_design_inputs(panels)
    if approved != design_record["world_500_decompressed_sha256"]:
        raise ValueError("Temperature panel differs from the frozen design")
    world = load_panel(panels, "world_500")
    primary = observations(world)
    if not pd.MultiIndex.from_frame(primary[KEY]).sort_values().equals(
        pd.MultiIndex.from_frame(eligible[KEY]).sort_values()
    ):
        raise ValueError("Temperature observations differ from the mapped eligible rows")
    pairs = {
        name: join_pair_temperatures(pd.read_csv(design / f"pairs_{name}.csv"), primary)
        for name in RULES
    }
    regression, fitted, _ = within_model(primary)
    paired = estimate_pairs(pairs["primary"])
    matched_sample = matched_endpoints(primary, pairs["primary"])
    on_matches, _, _ = within_model(matched_sample)
    design_diagnostics(design, primary, pairs["primary"], output)
    write_csv(pd.concat([frame.assign(design=name) for name, frame in pairs.items()], ignore_index=True),
              output / "paired_temperatures.csv")
    write_csv(fitted[[*KEY, *COORDS, *BLOCK, "response_c", "fitted_c", "residual_c"]], output / "fitted_cells.csv")

    curve, curve_cells, boot = within_model(primary, curved_tree=True, draws=SENSITIVITY_DRAWS)
    shape = []
    for t in (0.10, TREE_CENTER, 0.50):
        m = 20 * (t - TREE_CENTER)
        lower, upper = interval(boot[:, 0] + m * boot[:, -1])
        close = curve_cells[curve_cells.tree_fraction.between(t - 0.05, t + 0.05)]
        shape.append({
            "tree_fraction": t, "estimate_c_per_10pp": curve["estimate_c_per_10pp"] + m * curve["coef"]["tree_squared"],
            "lower_c": lower, "upper_c": upper, "cells_within_5pp": int(close.cell_id.nunique()),
            "blocks_within_5pp": len(close[BLOCK].drop_duplicates()),
        })
    write_csv(pd.DataFrame(shape), output / "curvature.csv")

    date_rows, omitted, spatial = [], [], []
    for date in DATES:
        day = primary[primary.date == date]
        one, _, _ = within_model(day, draws=SENSITIVITY_DRAWS)
        one_pair = estimate_pairs(pairs["primary"][pairs["primary"].date == date], draws=SENSITIVITY_DRAWS)
        date_rows.append({"date": date, **result_row("Within-area regression", one)})
        date_rows.append({"date": date, **result_row("Nearby pairs", one_pair)})
        without, _, _ = within_model(primary[primary.date != date], draws=SENSITIVITY_DRAWS)
        without_pair = estimate_pairs(pairs["primary"][pairs["primary"].date != date], draws=SENSITIVITY_DRAWS)
        omitted.append({"omitted_date": date, **result_row("Within-area regression", without)})
        omitted.append({"omitted_date": date, **result_row("Nearby pairs", without_pair)})
        cells = fitted[fitted.date == date]
        east, north = cells[COORDS].to_numpy().T
        raw_i, edges = moran(cells.response_c.to_numpy(), east, north)
        resid_i, _ = moran(cells.residual_c.to_numpy(), east, north)
        spatial.append({"date": date, "used_cells": len(cells), "queen_neighbors": edges,
                        "lst_moran_i": raw_i, "residual_moran_i": resid_i})
    write_csv(pd.DataFrame(date_rows), output / "by_date.csv")
    write_csv(pd.DataFrame(omitted), output / "leave_one_date_out.csv")
    write_csv(pd.DataFrame(spatial), output / "residual_diagnostics.csv")

    variants = [result_row("Primary within-area regression", regression), result_row("Matched endpoint regression", on_matches),
                result_row("Quadratic tree term, slope at 25% trees", curve)]
    for label, records, options in [
        ("2 km comparisons and bootstrap", primary, {"model_size": 2000}),
        ("6 km comparisons and bootstrap", primary, {"model_size": 6000}),
        ("Shifted 4 km comparisons and bootstrap", primary, {"origin": 2000}),
        ("8 km nested bootstrap only", primary, {"cluster_size": 8000}),
        ("Equal weight per cell-date", primary, {"equal_dates": False}),
        ("Cells observed on all seven dates", primary[primary.groupby("cell_id").date.transform("size") == 7], {}),
        ("Esri water mask", observations(load_panel(panels, "esri_500")), {}),
        ("Cloud distance >=0.5 km and ST_QA <=4 K", observations(world, quality="4K"), {}),
        ("Cloud distance >=0.5 km and ST_QA <=2 K", observations(world, quality="2K"), {}),
    ]:
        if records.empty:
            variants.append({"analysis": label, "status": "no usable observations", "cell_dates": 0,
                             "distinct_cells": 0, "dates": 0})
            continue
        fit, _, _ = within_model(records, draws=SENSITIVITY_DRAWS, **options)
        variants.append(result_row(label, fit))
    for name, label in [
        ("primary", "Primary nearby pairs"), ("shifted", "Shifted 4 km pair blocks"),
        ("tighter", "Tighter distance and built-cover limits"),
    ]:
        fit = paired if name == "primary" else estimate_pairs(pairs[name], draws=SENSITIVITY_DRAWS)
        variants.append(result_row(label, fit))
    write_csv(pd.DataFrame(variants), output / "robustness.csv")
    summary = {
        "dates": DATES, "bootstrap_seed": SEED, "frozen_design_sha256": checksum(design / "design_record.json"),
        "within_area": regression, "nearby_pairs": paired, "regression_on_matched_cells": on_matches,
    }
    (output / "estimates.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")
    print(json.dumps(summary, indent=2, allow_nan=False))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--panels", type=Path, default=PANELS)
    parser.add_argument("--design", type=Path, default=DESIGN)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    run(args.panels, args.design, args.output)
