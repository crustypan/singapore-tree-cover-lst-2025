"""Estimate the tree-for-rangeland contrast on one frozen sample per comparison."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from analyze_landsat import interval, write_csv
from analyze_local_comparison import DRAWS, SEED, draw_totals
from analyze_overpass_validation import load_cohorts
from fetch_inputs import checksum
from make_validation_panel import SNAPSHOTS
from prepare_grass_tree_substitution import COVER, DESIGN, DISTANCES, GAPS, SIDES, VALIDATION
from prepare_local_comparison import BLOCK, KEY
from prepare_overpass_validation import verify_hashes


BASE = ["tree", "built", "built_squared", "other", "land"]
WATER = ["coast", "inland"]
MODELS = {
    "base": BASE, "coast": BASE + ["coast"], "inland": BASE + ["inland"], "both": BASE + WATER,
}


def join_pairs(pairs, observed, cover):
    lookup = observed[[*KEY, "response_c"]].assign(year=observed.date.str[:4].astype(int)).merge(
        cover[["year", "cell_id", *COVER, *DISTANCES]], on=["year", "cell_id"], validate="many_to_one"
    ).drop(columns="year")
    joined = pairs.copy()
    for side, key in zip(("lower", "higher"), SIDES):
        rows = lookup.rename(columns={"cell_id": key, **{
            col: f"{side}_{col}" for col in ["response_c", *COVER, *DISTANCES]
        }})
        joined = joined.merge(rows, on=["date", key], how="left", validate="one_to_one")
    cols = [f"{side}_{col}" for side in ("lower", "higher") for col in ["response_c", *COVER, *DISTANCES]]
    if len(joined) != len(pairs) or not np.isfinite(joined[cols].to_numpy()).all():
        raise ValueError("Frozen pairs differ from the available temperature or cover sample")
    for name, gap in zip(COVER, GAPS):
        np.testing.assert_allclose(100 * (joined[f"higher_{name}"] - joined[f"lower_{name}"]),
                                   pairs[gap], atol=1.1e-7, rtol=0)
    joined["lst_gap_c"] = joined.higher_response_c - joined.lower_response_c
    joined["lst_scaled_c"] = 10 * joined.lst_gap_c / joined.tree_gap_pp
    return joined


def pair_matrix(pairs, inland="inland_km", logarithm=True):
    def diff(name, transform):
        return transform(pairs[f"higher_{name}"].to_numpy()) - transform(pairs[f"lower_{name}"].to_numpy())

    def share(values):
        return 10 * values

    def squared_built(values):
        return (10 * (values - 0.50)) ** 2

    def transform_distance(km):
        return np.log1p(km) if logarithm else km

    return np.column_stack([
        diff("tree_fraction", share), diff("built_fraction", share), diff("built_fraction", squared_built),
        diff("other_fraction", share), diff("land_share", share),
        diff("coast_km", transform_distance), diff(inland, transform_distance),
    ])


def shared_fits(x, y, scaled, blocks, dates, draws):
    """Use identical pair weights and spatial draws for all nested models."""
    p = x.shape[1]
    gram_rows = np.einsum("ni,nj->nij", x, x).reshape(len(x), p ** 2)
    row_stats = np.column_stack((gram_rows, x * y[:, None], scaled))
    moments, n_blocks = draw_totals(row_stats, blocks, dates, draws)
    totals = (moments[:, :, 1:] / moments[:, :, :1]).sum(axis=1)
    gram = totals[:, :p ** 2].reshape(-1, p, p)
    rhs = totals[:, p ** 2:p ** 2 + p]
    # Full rank of the complete model guarantees rank of its nested submodels.
    eig = np.linalg.eigvalsh(gram)
    valid = eig[:, 0] > np.finfo(float).eps * p * eig[:, -1]
    if not valid[0]:
        raise ValueError("The full pair-difference design matrix is rank deficient")
    gram, rhs = gram[valid], rhs[valid]
    coef = {}
    for model, names in MODELS.items():
        cols = [(BASE + WATER).index(name) for name in names]
        coef[model] = np.linalg.solve(gram[:, cols][:, :, cols], rhs[:, cols, None])[:, :, 0]
    raw = totals[valid, -1] / moments.shape[1]
    return coef, raw, n_blocks


def summarize(pairs, draws, inland="inland_km", logarithm=True):
    x = pair_matrix(pairs, inland=inland, logarithm=logarithm)
    y = pairs.lst_gap_c.to_numpy()
    coef, raw, n_blocks = shared_fits(x, y, pairs.lst_scaled_c.to_numpy(), pairs[BLOCK].to_numpy(),
                                    pairs.date.to_numpy(), draws)
    w = 1 / pairs.groupby("date").date.transform("size").to_numpy()
    root_w = np.sqrt(w)
    pair_ids = pairs[["date", *SIDES]].to_csv(index=False, lineterminator="\n").encode()
    n_cells = int(pd.concat([pairs[s] for s in SIDES], ignore_index=True).nunique())
    fits = []

    def row(model, values, change, rmse):
        ci = interval(values[1:])
        shift_ci = interval(change[1:])
        return {
            "model": model, "estimate_c": float(values[0]), "lower_c": float(ci[0]), "upper_c": float(ci[1]),
            "change_c": float(change[0]), "change_lower_c": float(shift_ci[0]), "change_upper_c": float(shift_ci[1]),
            "rmse_c": rmse, "pairs": len(pairs), "distinct_cells": n_cells,
            "dates": int(pairs.date.nunique()), "blocks": n_blocks,
            "draws_requested": draws, "draws_retained": len(values) - 1,
        }

    if len(raw) < 2:
        raise ValueError("No estimable spatial resamples for these observations")
    for model, names in MODELS.items():
        design = x[:, [(BASE + WATER).index(name) for name in names]]
        check, _, rank, _ = np.linalg.lstsq(design * root_w[:, None], y * root_w, rcond=None)
        if rank != len(names):
            raise ValueError(f"Rank check failed in the {model} model")
        np.testing.assert_allclose(coef[model][0], check, atol=1e-9, rtol=0)
        fitted = design @ coef[model][0]
        residual = float(np.sqrt(np.average((y - fitted) ** 2, weights=w)))
        change = coef[model][:, 0] - coef["base"][:, 0]
        fits.append(row(model, coef[model][:, 0], change, residual))
    descriptive = row("paired descriptive", raw, raw - raw, None)
    for key in ("change_c", "change_lower_c", "change_upper_c"):
        descriptive[key] = None
    fits.append(descriptive)
    return fits, {
        "model_coefficients": {
            model: dict(zip(MODELS[model], map(float, values[0]))) for model, values in coef.items()
        },
        "pair_id_sha256": hashlib.sha256(pair_ids).hexdigest(),
        "weighted_design_condition_number": float(np.linalg.cond(x * root_w[:, None])),
    }


def balance(pairs):
    columns = ["distance_m", *GAPS, "lst_gap_c", "lst_scaled_c"]
    for name in ["built_gap_pp", "other_gap_pp", "land_gap_pp"]:
        key = "absolute_" + name
        columns.append(key)
    values = pairs.assign(**{f"absolute_{key}": pairs[key].abs()
                             for key in ["built_gap_pp", "other_gap_pp", "land_gap_pp"]})
    daily = values.groupby("date")[columns].mean()
    daily.insert(0, "pairs", pairs.groupby("date").size())
    daily.loc["equal_date_mean"] = daily.drop(columns="pairs").mean()
    return daily.reset_index()


def influence(pairs):
    dates = sorted(pairs.date.unique())
    rows = []
    for date in dates:
        g = pairs[pairs.date != date]
        coef, _, blocks = shared_fits(pair_matrix(g), g.lst_gap_c.to_numpy(), g.lst_scaled_c.to_numpy(),
                                     g[BLOCK].to_numpy(), g.date.to_numpy(), draws=0)
        base = float(coef["base"][0, 0])
        both = float(coef["both"][0, 0])
        rows.append({"excluded_date": date, "pairs": len(g), "blocks": blocks,
                     "base_c": base, "both_c": both, "change_c": both - base})
    return pd.DataFrame(rows)


def analyze(design, snapshots, output):
    record = json.loads((design / "design_record.json").read_text())
    verify_hashes(design, record["design_sha256"])
    verify_hashes(VALIDATION, {"design_record.json": record["validation_design_sha256"]})
    validation = json.loads((VALIDATION / "design_record.json").read_text())
    temperatures = pd.concat(load_cohorts(snapshots, VALIDATION, validation).values(), ignore_index=True)
    cover = pd.read_csv(design / "cell_covariates.csv", float_precision="round_trip")
    primary, tighter = (join_pairs(pd.read_csv(design / f"pairs_{kind}.csv"), temperatures, cover)
                        for kind in ("primary", "tighter"))
    recent = primary.date.str.startswith("2025-")
    variants = [
        ("primary", primary, {}), ("tighter", tighter, {}),
        ("2025 dates", primary[recent], {}), ("2023–24 dates", primary[~recent], {}),
        ("raw distance", primary, {"logarithm": False}),
        ("0.10 km² water", primary, {"inland": "inland_010_km"}),
    ]
    output.mkdir(parents=True, exist_ok=True)
    rows, details = [], {}
    for variant, data, kwargs in variants:
        fits, info = summarize(data, DRAWS, **kwargs)
        rows.extend({"analysis": variant, **fit} for fit in fits)
        details[variant] = info
    estimates = pd.DataFrame(rows)
    write_csv(estimates, output / "estimates.csv")
    write_csv(balance(primary), output / "pair_balance_by_date.csv")
    write_csv(influence(primary), output / "leave_one_overpass_out.csv")
    write_csv(primary[["date", *SIDES, "lower_response_c", "higher_response_c", "lst_gap_c", "lst_scaled_c"]],
              output / "paired_temperatures.csv")
    analysis = {
        "design_sha256": checksum(design / "design_record.json"), "seed": SEED, "details": details,
        "output_sha256": {name: checksum(output / name) for name in (
            "estimates.csv", "pair_balance_by_date.csv", "leave_one_overpass_out.csv", "paired_temperatures.csv"
        )},
    }
    (output / "analysis_record.json").write_text(json.dumps(analysis, indent=2) + "\n")
    print(estimates.loc[estimates.model.isin(["base", "both"])].to_string(index=False))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--design", type=Path, default=DESIGN)
    parser.add_argument("--snapshots", type=Path, default=SNAPSHOTS)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    analyze(args.design, args.snapshots, args.output)
