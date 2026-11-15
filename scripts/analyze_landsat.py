"""Fit selected-date spatial associations and check coverage, dependence, and sensitivity."""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from make_landsat_panel import DATES, GRID_VARIANTS

SEED = 20251115
MAIN_DRAWS = 1999
CHECK_DRAWS = 999
BUILT_LEVELS = (0.25, 0.50, 0.75)


def load_panel(directory, grid, dates=DATES):
    path = directory / (grid + ".csv")
    if not path.is_file():
        path = path.with_suffix(".csv.gz")
    data = pd.read_csv(path, float_precision="round_trip")
    if data.duplicated(["cell_id", "date"]).any() or tuple(sorted(data.date.unique())) != tuple(dates):
        raise ValueError(f"Invalid cell-date IDs for {grid}")
    if data.groupby("cell_id").date.size().ne(len(dates)).any():
        raise ValueError(f"Each cell should have a row for each date in {grid}")
    if not ((data.tree_fraction + data.built_fraction).dropna() <= 1 + 1e-8).all():
        raise ValueError(f"Land-cover fractions overlap for {grid}")
    if not ((data.strict_lst_centers <= data.qa4_lst_centers)
            & (data.qa4_lst_centers <= data.valid_lst_centers)).all():
        raise ValueError(f"Strict QA must be a subset of the primary QA in {grid}")
    return data


def area_rule(data, land=0.5):
    return (data.land_share >= land) & (data.classified_share >= 0.8)


def observations(data, valid=0.6, land=0.5, quality="primary"):
    share, temp = {
        "primary": ("valid_lst_share", "lst_c"),
        "2K": ("strict_lst_share", "strict_lst_c"),
        "4K": ("qa4_lst_share", "qa4_lst_c"),
    }[quality]
    kept = data.loc[area_rule(data, land) & (data[share] >= valid) & data[temp].notna()].copy()
    kept["response_c"] = kept[temp]
    return kept


def block_coordinates(data, size, shift=0):
    return np.floor((data[["easting_m", "northing_m"]].to_numpy() - shift) / size).astype(int)


def design(data, model, curved=False):
    day = pd.get_dummies(data.date, drop_first=True, dtype=float).to_numpy()
    tree = (data.tree_fraction.to_numpy() - (0.1 if curved else 0)) * 10
    built = (data.built_fraction.to_numpy() - 0.5) * 10
    terms = [np.ones(len(data)), tree]
    names = ["intercept_reference_date", "tree_at_50pct_built"]
    if model in ("additive", "interaction"):
        terms.extend([built, data.land_share.to_numpy()])
        names.extend(["built_at_reference_tree_fraction", "full_land_share"])
    if model == "interaction":
        terms.append(tree * built)
        names.append("change_in_tree_slope_per_10pp_built")
    if curved:
        # At 10% tree cover the linear term describes the local slope.
        terms.extend([tree ** 2, built ** 2])
        names.extend(["tree_squared", "built_squared"])
    names.extend("date_" + d for d in sorted(data.date.unique())[1:])
    return np.column_stack([*terms, day]), names


def interval(values):
    return np.quantile(values, [0.025, 0.975], axis=0)


def block_bootstrap(x, y, weights, blocks, draws):
    _, group = np.unique(blocks, axis=0, return_inverse=True)
    n = int(group.max()) + 1
    p = x.shape[1]
    gram, rhs = np.empty((n, p, p)), np.empty((n, p))
    for g in range(n):
        members = group == g
        xg, wg = x[members], weights[members]
        gram[g] = xg.T @ (wg[:, None] * xg)
        rhs[g] = xg.T @ (wg * y[members])
    rng = np.random.default_rng(SEED)
    # Each row contains the counts when drawing n blocks with replacement.
    repetitions = rng.multinomial(n, np.full(n, 1 / n), size=draws)
    systems = (repetitions @ gram.reshape(n, p * p)).reshape(draws, p, p)
    targets = repetitions @ rhs
    draws_ok = []
    for system, target in zip(systems, targets):
        if np.linalg.matrix_rank(system, hermitian=True) == p:
            draws_ok.append(np.linalg.solve(system, target))
    if len(draws_ok) < 0.95 * draws:
        raise ValueError(f"Only {len(draws_ok)} of {draws} spatial bootstrap samples were estimable")
    return np.asarray(draws_ok), n


def fit(data, model="interaction", draws=0, block_size=4000, shift=0, weighting="row",
        spatial_intercepts=False, curved=False):
    x, names = design(data, model, curved)
    y = data.response_c.to_numpy()
    weights = np.ones(len(data))
    if weighting == "cell":
        weights /= data.groupby("cell_id").cell_id.transform("size").to_numpy()
    root_w = np.sqrt(weights)
    if spatial_intercepts:
        # Remove block means before estimating slopes. The fixed effects themselves
        # may vary by resample, so the bootstrap uses the within-block design.
        groups = np.unique(block_coordinates(data, 4000), axis=0, return_inverse=True)[1]
        frame = pd.DataFrame(np.column_stack([y, x[:, 1:]]))
        weighted = frame.mul(weights, axis=0).groupby(groups).transform("sum")
        weight_sum = pd.Series(weights).groupby(groups).transform("sum")
        centered = (frame - weighted.div(weight_sum, axis=0)).to_numpy()
        x_fit, y_fit, fit_names = centered[:, 1:], centered[:, 0], names[1:]
    else:
        x_fit, y_fit, fit_names = x, y, names
    coef, _, rank, singular = np.linalg.lstsq(root_w[:, None] * x_fit, root_w * y_fit, rcond=None)
    if rank != x_fit.shape[1]:
        raise ValueError(f"Singular {model} model")
    residual = y_fit - x_fit @ coef
    layout = block_coordinates(data, block_size, shift)
    checked, block_count = None, len(np.unique(layout, axis=0))
    if draws:
        checked, block_count = block_bootstrap(x_fit, y_fit, weights, layout, draws)
    coefficient = dict(zip(fit_names, coef))
    idx = {name: i for i, name in enumerate(fit_names)}
    dates = sorted(data.date.unique())
    stats = {
        "rows": len(data), "cells": int(data.cell_id.nunique()), "dates": len(dates),
        "spatial_blocks": block_count, "bootstrap_draws_retained": len(checked) if draws else 0,
        "rank": int(rank), "condition_number": float(singular[0] / singular[-1]),
        "rmse_c": float(np.sqrt(np.average(residual ** 2, weights=weights))),
        "within_date_r2": None,
        "tree_fraction_at_slope": 0.10 if curved else None,
        "coefficients": coefficient,
    }
    if not spatial_intercepts:
        w_data = pd.DataFrame({"w": weights, "wy": weights * y, "date": data.date.to_numpy()})
        totals = w_data.groupby("date")[["wy", "w"]].transform("sum")
        within = y - (totals.wy / totals.w).to_numpy()
        stats["within_date_r2"] = float(1 - np.dot(weights, residual ** 2) / np.dot(weights, within ** 2))
    tree_idx = idx["tree_at_50pct_built"]
    item = "change_in_tree_slope_per_10pp_built"
    slopes = []
    if item in idx:
        b_i = idx[item]
        for b in BUILT_LEVELS:
            mult = (b - 0.5) * 10
            bounds = [None, None] if checked is None else interval(checked[:, tree_idx] + mult * checked[:, b_i])
            slopes.append({"built_fraction": b, "estimate_c": float(coef[tree_idx] + mult * coef[b_i]),
                           "lower_c": bounds[0], "upper_c": bounds[1]})
        values = [None, None] if checked is None else interval(checked[:, b_i])
        stats["interaction"] = {"estimate_c": float(coef[b_i]), "lower_c": values[0], "upper_c": values[1]}
    else:
        bounds = [None, None] if checked is None else interval(checked[:, tree_idx])
        slopes.append({"built_fraction": 0.5, "estimate_c": float(coef[tree_idx]),
                       "lower_c": bounds[0], "upper_c": bounds[1]})
    stats["slopes"] = slopes
    return stats, residual


def moran(values, x, y, grid_size=500):
    coords = list(zip(x, y))
    lookup = {point: i for i, point in enumerate(coords)}
    edges = []
    for i, (east, north) in enumerate(coords):
        for dx, dy in ((grid_size, 0), (0, grid_size), (grid_size, grid_size), (grid_size, -grid_size)):
            j = lookup.get((east + dx, north + dy))
            if j is not None:
                edges.append((i, j))
    centered = values - values.mean()
    if not edges or np.dot(centered, centered) == 0:
        return float("nan"), len(edges)
    left, right = np.array(edges).T
    return float(len(values) / len(edges) * (centered[left] @ centered[right]) / (centered @ centered)), len(edges)


def write_csv(value, path):
    value.to_csv(path, index=False, float_format="%.9f", lineterminator="\n")


def coverage_reports(data, kept, output):
    all_cells = data.drop_duplicates("cell_id").set_index("cell_id")
    all_cells["model_eligible"] = area_rule(all_cells)
    all_cells["dates_used"] = kept.groupby("cell_id").size().reindex(all_cells.index, fill_value=0)
    eligible = all_cells[all_cells.model_eligible]
    groups = [
        ("spatial eligibility failed", all_cells[~all_cells.model_eligible]),
        ("eligible, never observed", eligible[eligible.dates_used == 0]),
        ("eligible, 1–6 dates", eligible[(eligible.dates_used > 0) & (eligible.dates_used < 7)]),
        ("eligible, all 7 dates", eligible[eligible.dates_used == 7]),
    ]
    coverage = []
    for label, cells in groups:
        coverage.append({
            "group": label, "cells": len(cells), "area_km2": cells.land_centers.sum() / 10000,
            "mean_tree_fraction": cells.tree_fraction.mean(),
            "mean_built_fraction": cells.built_fraction.mean(),
            "mean_land_share": cells.land_share.mean(),
        })
    write_csv(pd.DataFrame(coverage), output / "sample_availability.csv")
    write_csv(all_cells.reset_index()[["cell_id", "easting_m", "northing_m", "model_eligible", "dates_used",
              "tree_fraction", "built_fraction", "land_share"]], output / "cell_availability.csv")
    coverage_by_cover = []
    for variable in ("tree_fraction", "built_fraction"):
        for lo in (0.0, 0.25, 0.5, 0.75):
            hi = lo + 0.25
            group = eligible[(eligible[variable] >= lo)
                             & ((eligible[variable] < hi) | ((hi == 1.0) & (eligible[variable] == hi)))]
            coverage_by_cover.append({
                "variable": variable, "from_fraction": lo, "to_fraction": hi, "cells": len(group),
                "share_observed_on_at_least_one_date": (group.dates_used > 0).mean(),
                "share_observed_on_all_seven_dates": (group.dates_used == 7).mean(),
                "mean_dates_used": group.dates_used.mean(),
            })
    write_csv(pd.DataFrame(coverage_by_cover), output / "availability_by_cover.csv")
    support = []
    used = all_cells[all_cells.dates_used > 0]
    for b in BUILT_LEVELS:
        close = used[(used.built_fraction - b).abs() <= 0.10 + 1e-12]
        support.append({
            "built_fraction_reference": b, "cells_within_10pp": len(close),
            "tree_fraction_min": close.tree_fraction.min(),
            "tree_fraction_10th_pct": close.tree_fraction.quantile(0.1),
            "tree_fraction_90th_pct": close.tree_fraction.quantile(0.9),
            "tree_fraction_max": close.tree_fraction.max(),
        })
    write_csv(pd.DataFrame(support), output / "composition_support.csv")
    return {
        "cells_intersecting_primary_land": int(len(all_cells)), "spatially_eligible_cells": int(len(eligible)),
        "eligible_with_zero_usable_dates": int((eligible.dates_used == 0).sum()),
        "usable_all_seven_dates": int((eligible.dates_used == 7).sum()),
        "total_worldcover_land_km2": float(all_cells.land_centers.sum() / 10000),
        "sampled_worldcover_land_km2": float(used.land_centers.sum() / 10000),
    }


def compact_variant(label, stats):
    slope = next(s for s in stats["slopes"] if s["built_fraction"] == 0.5)
    interaction = stats.get("interaction", {})
    return {
        "variant": label, "rows": stats["rows"], "cells": stats["cells"], "dates": stats["dates"],
        "status": "estimated",
        "tree_fraction_at_slope": stats["tree_fraction_at_slope"],
        "spatial_blocks": stats["spatial_blocks"], "draws_retained": stats["bootstrap_draws_retained"],
        "tree_at_50pct_c": slope["estimate_c"], "tree_lower_c": slope["lower_c"], "tree_upper_c": slope["upper_c"],
        "interaction_c": interaction.get("estimate_c"), "interaction_lower_c": interaction.get("lower_c"),
        "interaction_upper_c": interaction.get("upper_c"),
    }


def run(datafiles, output):
    output.mkdir(parents=True, exist_ok=True)
    panels = {name: load_panel(datafiles, name) for name in GRID_VARIANTS}
    data = panels["world_500"]
    primary = observations(data)
    sample = coverage_reports(data, primary, output)
    main, residual = fit(primary, draws=MAIN_DRAWS)
    primary["residual_c"] = residual
    base, _ = fit(primary, model="tree_only", draws=MAIN_DRAWS)
    additive, _ = fit(primary, model="additive", draws=MAIN_DRAWS)
    summary = {
        "study_dates": list(DATES), "seed": SEED, "sample": sample,
        "primary_fit": main, "tree_only_fit": base, "additive_fit": additive,
        "tree_built_pearson_cell_level": float(
            primary.drop_duplicates("cell_id")[["tree_fraction", "built_fraction"]].corr().iloc[0, 1]
        ),
    }
    (output / "estimates.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")
    write_csv(primary[["date", "cell_id", "easting_m", "northing_m", "response_c", "residual_c"]],
              output / "fitted_cells.csv")
    by_date, drop_date, record = [], [], []
    eligible = data[area_rule(data)].groupby("date").size()
    for date, cells in primary.groupby("date", sort=True):
        stats, _ = fit(cells, draws=CHECK_DRAWS)
        by_date.append(compact_variant(date, stats))
        without, _ = fit(primary[primary.date != date])
        drop_date.append(compact_variant(date, without))
        xy = cells[["easting_m", "northing_m"]].to_numpy().T
        moran_lst, links = moran(cells.response_c.to_numpy(), *xy)
        moran_res, _ = moran(cells.residual_c.to_numpy(), *xy)
        record.append({
            "date": date, "eligible_cells": int(eligible[date]), "used_cells": len(cells),
            "used_land_area_km2": cells.land_centers.sum() / 10000,
            "mean_lst_c": cells.response_c.mean(), "median_lst_c": cells.response_c.median(),
            "lst_p01_c": cells.response_c.quantile(0.01), "lst_p99_c": cells.response_c.quantile(0.99),
            "mean_tree_fraction": cells.tree_fraction.mean(), "mean_built_fraction": cells.built_fraction.mean(),
            "mean_observed_minus_fixed_tree": (cells.observed_tree_fraction - cells.tree_fraction).mean(),
            "mean_abs_observed_minus_fixed_tree": (cells.observed_tree_fraction - cells.tree_fraction).abs().mean(),
            "mean_observed_minus_fixed_built": (cells.observed_built_fraction - cells.built_fraction).mean(),
            "queen_pairs": links, "moran_lst": moran_lst, "moran_residual": moran_res,
        })
    write_csv(pd.DataFrame(by_date).rename(columns={"variant": "date"}), output / "date_estimates.csv")
    write_csv(pd.DataFrame(drop_date).rename(columns={"variant": "omitted_date"}), output / "leave_one_date_out.csv")
    write_csv(pd.DataFrame(record), output / "date_diagnostics.csv")

    strict = observations(data, quality="2K")
    relaxed = observations(data, quality="4K")
    matched = primary.merge(relaxed[["cell_id", "date"]], on=["cell_id", "date"], validate="one_to_one")
    repeat = primary.groupby("cell_id").date.transform("size")
    varying = primary.assign(tree_fraction=primary.observed_tree_fraction, built_fraction=primary.observed_built_fraction)
    trimmed = primary[primary.groupby("date").response_c.transform(
        lambda z: z.between(z.quantile(0.01), z.quantile(0.99))
    )]
    options = {
        "2 km blocks": (primary, {"block_size": 2000}),
        "6 km blocks": (primary, {"block_size": 6000}),
        "Shifted 4 km blocks": (primary, {"shift": 2000}),
        "40% valid LST area": (observations(data, valid=0.4), {}),
        "80% valid LST area": (observations(data, valid=0.8), {}),
        "Cloud distance >=0.5 km and ST_QA <=2 K": (strict, {}),
        "Cloud distance >=0.5 km and ST_QA <=4 K (exploratory)": (relaxed, {}),
        "Primary temperatures on common 4 K sample": (matched, {}),
        "Esri water mask": (observations(panels["esri_500"]), {}),
        "90% minimum land area": (observations(data, land=0.9), {}),
        "Shifted 500 m grid": (observations(panels["world_500_shift250"]), {}),
        "1 km grid": (observations(panels["world_1000"]), {}),
        "Cells with all seven dates": (primary[repeat == 7], {}),
        "Equal weight per spatial cell": (primary, {"weighting": "cell"}),
        "4 km block intercepts": (primary, {"spatial_intercepts": True}),
        "Land cover inside observed footprint": (varying, {}),
        "Both tree and built squared terms (exploratory)": (primary, {"curved": True}),
        "Omit date-specific LST tails (exploratory)": (trimmed, {}),
    }
    variants = [compact_variant("Primary", main)]
    for name, (records, settings) in options.items():
        if records.empty:
            variants.append({"variant": name, "status": "no usable observations", "rows": 0, "cells": 0, "dates": 0})
            continue
        draws = MAIN_DRAWS if name in ("2 km blocks", "6 km blocks", "Shifted 4 km blocks") else CHECK_DRAWS
        stats, _ = fit(records, draws=draws, **settings)
        variants.append(compact_variant(name, stats))
    write_csv(pd.DataFrame(variants), output / "robustness.csv")
    print(f"Primary observations  {len(primary):,} in {main['cells']:,} cells / {main['spatial_blocks']} spatial blocks")
    print("Tree slopes", json.dumps(main["slopes"]))
    print("Interaction", json.dumps(main["interaction"]))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--panels", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    run(args.panels, args.output)
