"""Estimate the change in a ten-point tree-cover contrast across built-cover references."""
import argparse
from dataclasses import dataclass
import json
from pathlib import Path

import numpy as np
import pandas as pd

from analyze_landsat import block_coordinates, interval, load_panel, observations, write_csv
from analyze_local_comparison import DRAWS, SEED, SENSITIVITY_DRAWS, draw_totals, join_pair_temperatures
from fetch_inputs import checksum
from prepare_built_support import DESIGN_DIR, STRATUM, selected_samples
from prepare_local_comparison import BLOCK, KEY, MatchRule, PANELS, assign_blocks, read_design_inputs


@dataclass(frozen=True)
class Contrast:
    tree_start_pct: float
    lower_built_pct: float
    higher_built_pct: float

    @property
    def tree_end_pct(self):
        return self.tree_start_pct + 10


def basis(data, contrast, model):
    t = (100 * data.tree_fraction - contrast.tree_start_pct - 5) / 10
    built_center = (contrast.lower_built_pct + contrast.higher_built_pct) / 2
    b = (100 * data.built_fraction - built_center) / 10
    x = pd.DataFrame({"tree": t, "built": b, "land_share": data.land_share})
    if model in ("additive", "quadratic"):
        x["tree_squared"] = t ** 2
        x["built_squared"] = b ** 2
    if model != "additive":
        x["tree_by_built"] = t * b
    if model == "piecewise":
        ht, hb = t.clip(lower=0), b.clip(lower=0)
        x["tree_hinge"], x["built_hinge"] = ht, hb
        x["tree_hinge_by_built"] = ht * b
        x["tree_by_built_hinge"] = t * hb
        x["tree_hinge_by_built_hinge"] = ht * hb
    return x


def contrast_matrix(contrast, model):
    rows = []
    for built in (contrast.lower_built_pct, contrast.higher_built_pct):
        predictors = basis(pd.DataFrame({
            "tree_fraction": np.array([contrast.tree_start_pct, contrast.tree_end_pct]) / 100,
            "built_fraction": [built / 100, built / 100], "land_share": [0.95, 0.95],
        }), contrast, model).to_numpy()
        rows.append(predictors[1] - predictors[0])
    return np.stack([*rows, rows[1] - rows[0]])


def fit_contrasts(data, contrast, model="quadratic", size=4000, shift=0, resample_size=None,
                  weight="date", band_intercepts=False, draws=DRAWS):
    if resample_size is None:
        resample_size = size
    if resample_size % size or (weight == "date_band" and not band_intercepts):
        raise ValueError("Resampled blocks must nest the model blocks and weights must be constant in each stratum")
    data = assign_blocks(data, MatchRule(block_size=size, block_shift=shift))
    by = [*STRATUM, "band"] if band_intercepts else STRATUM
    counts = data.groupby(by).cell_id.transform("size")
    records = data[counts >= 2].copy()
    both = records.groupby("date").band.nunique()
    records = records[records.date.isin(both[both == 2].index)].reset_index(drop=True)
    if records.empty:
        return {"status": "no usable shared dates", "cell_dates": 0, "cells": 0, "dates_used": 0}
    variables = basis(records, contrast, model)
    cols = [*variables, "response_c"]
    f = pd.concat([records[by], variables, records[["response_c"]]], axis=1)
    within = f[cols] - f.groupby(by)[cols].transform("mean")
    x, y = within[variables.columns].to_numpy(), within.response_c.to_numpy()
    p = x.shape[1]
    products = np.einsum("ni,nj->nij", x, x).reshape(len(x), p * p)
    # Requiring both arms on each retained date keeps the contrast and its uncertainty
    # tied to a common set of overpasses. The resampling unit is still the spatial block.
    by_date_band = records.date + "|" + records.band
    moments, n = draw_totals(np.column_stack([products, x * y[:, None]]),
                             block_coordinates(records, resample_size, shift), by_date_band, draws)
    date_count = int(records.date.nunique())
    moments = moments.reshape(len(moments), date_count, 2, -1)
    if weight == "date_band":
        unweighted = (moments[:, :, :, 1:] / moments[:, :, :, :1]).sum(axis=2)
    else:
        per_date = moments.sum(axis=2)
        unweighted = per_date[:, :, 1:]
        if weight == "date":
            unweighted = unweighted / per_date[:, :, :1]
    totals = unweighted.sum(axis=1)
    gram, rhs = totals[:, :p * p].reshape(-1, p, p), totals[:, p * p:]
    eig = np.linalg.eigvalsh(gram)
    full_rank = eig[:, 0] > np.finfo(float).eps * p * eig[:, -1]
    info = {
        "cell_dates": len(records), "cells": int(records.cell_id.nunique()), "dates_used": date_count,
        "date_labels": ", ".join(sorted(records.date.unique())), "bootstrap_blocks": n,
        "lower_band_blocks": len(records.loc[records.band == "lower", BLOCK].drop_duplicates()),
        "higher_band_blocks": len(records.loc[records.band == "higher", BLOCK].drop_duplicates()),
        "draws_requested": draws, "draws_with_both_bands": len(moments) - 1,
        "draws_retained": int(full_rank[1:].sum()),
        "strata": len(records[by].drop_duplicates()),
    }
    if not full_rank[0]:
        return {"status": "rank deficient", **info}
    coefs = np.linalg.solve(gram[full_rank], rhs[full_rank, :, None]).squeeze(axis=2)
    delta = coefs @ contrast_matrix(contrast, model).T
    report = []
    for i, key in enumerate(("lower_built", "higher_built", "higher_minus_lower")):
        lower, upper = (
            map(float, interval(delta[1:, i])) if draws and len(delta) > 1 and n > 1 else [None, None]
        )
        report.append((key, float(delta[0, i]), lower, upper))
    weights = np.ones(len(records))
    if weight == "date":
        weights /= records.groupby("date").cell_id.transform("size").to_numpy()
    if weight == "date_band":
        weights /= records.groupby(["date", "band"]).cell_id.transform("size").to_numpy()
    rmse = float(np.sqrt(np.average((y - x @ coefs[0]) ** 2, weights=weights)))
    return {
        "status": "estimated", **info, "within_rmse_c": rmse,
        "condition_number": float(np.sqrt(eig[0, -1] / eig[0, 0])),
        "coefficients": dict(zip(variables, map(float, coefs[0]))),
        "contrasts": {name: {"estimate_c": est, "lower_c": lo, "upper_c": hi}
                      for name, est, lo, hi in report},
    }


def variant_row(label, model, target, stats):
    row = {
        "analysis": label, "model": model, "tree_start_pct": target.tree_start_pct,
        "tree_end_pct": target.tree_end_pct, "lower_built_pct": target.lower_built_pct,
        "higher_built_pct": target.higher_built_pct,
        **{k: v for k, v in stats.items() if k not in ("contrasts", "coefficients")},
    }
    for name, v in stats.get("contrasts", {}).items():
        row.update({f"{name}_{k}": val for k, val in v.items()})
    return row


def pair_check(pairs, shared):
    if shared:
        pairs = pairs[pairs.shared_strict]
    keys = STRATUM if shared else ["date"]
    means = pairs.groupby([*keys, "band"]).lst_per_10pp_c.mean().unstack()
    if not {"lower", "higher"}.issubset(means.columns):
        return {"status": "no comparable pair comparisons", "dates": 0, "spatial_blocks": 0}
    means = means.dropna(subset=["lower", "higher"])
    if means.empty:
        return {"status": "no comparable pair comparisons", "dates": 0, "spatial_blocks": 0}
    selected = pairs.merge(means.reset_index()[keys], on=keys, how="inner", validate="many_to_one")
    by_date = means.groupby(level=0).mean() if shared else means
    lo, hi = by_date["lower"].mean(), by_date["higher"].mean()
    return {
        "status": "descriptive", "lower_built_c": float(lo), "higher_built_c": float(hi),
        "higher_minus_lower_c": float(hi - lo), "dates": len(by_date),
        "spatial_blocks": len(selected[BLOCK].drop_duplicates()), "groups": len(means),
        "lower_pairs": int(selected.band.eq("lower").sum()), "higher_pairs": int(selected.band.eq("higher").sum()),
        "lower_blocks": len(selected.loc[selected.band == "lower", BLOCK].drop_duplicates()),
        "higher_blocks": len(selected.loc[selected.band == "higher", BLOCK].drop_duplicates()),
    }


def local_sample(data, contrast, tree_radius=5., built_radius=10.):
    rows, _, _, _ = selected_samples(assign_blocks(data), contrast.tree_start_pct,
                                     contrast.lower_built_pct, contrast.higher_built_pct,
                                     tree_radius, built_radius)
    return rows


def run(panels, design, output):
    design_info = json.loads((design / "design_record.json").read_text())
    for filename, digest in design_info["design_sha256"].items():
        if checksum(design / filename) != digest:
            raise ValueError(f"Frozen design changed for {filename}")
    _, _, panel_hash = read_design_inputs(panels)
    if panel_hash != design_info["world_500_decompressed_sha256"]:
        raise ValueError("Frozen design and archived panel differ")
    target = Contrast(**{k: design_info["chosen"][k] for k in Contrast.__dataclass_fields__})
    data = load_panel(panels, "world_500")
    primary = observations(data)
    registry = pd.read_csv(design / "fit_sample.csv")
    inputs = primary.merge(registry, on=KEY, how="inner", validate="one_to_one")
    if len(inputs) != len(registry):
        raise ValueError("A frozen cell-date is missing its temperature")
    shared = inputs[inputs.shared]
    strict = join_pair_temperatures(pd.read_csv(design / "strict_pairs.csv"), primary)
    additive = fit_contrasts(inputs, target, model="additive")
    interaction = fit_contrasts(inputs, target)
    if interaction["status"] != "estimated":
        raise ValueError("Primary design has no estimable interaction")
    comparisons = [
        variant_row("Additive quadratic", "additive", target, additive),
        variant_row("Primary quadratic interaction", "quadratic", target, interaction),
    ]
    for label, model, observations_, options in [
        ("Linear product", "linear", inputs, {}),
        ("Piecewise linear interaction", "piecewise", inputs, {}),
        ("Shared block-dates", "quadratic", shared, {}),
        ("Shared block-dates with balanced bands", "quadratic", shared,
         {"weight": "date_band", "band_intercepts": True}),
        ("Narrower support windows", "quadratic", local_sample(primary, target, 2.5, 7.5), {}),
        ("Shifted 4 km local blocks", "quadratic", inputs, {"shift": 2000}),
        ("2 km local blocks", "quadratic", inputs, {"size": 2000}),
        ("6 km local blocks", "quadratic", inputs, {"size": 6000}),
        ("8 km nested bootstrap", "quadratic", inputs, {"resample_size": 8000}),
        ("Equal cell-date weights", "quadratic", inputs, {"weight": "row"}),
        ("Esri water mask", "quadratic", local_sample(observations(load_panel(panels, "esri_500")), target), {}),
        ("Cloud distance >=0.5 km and ST_QA <=4 K", "quadratic",
         local_sample(observations(data, quality="4K"), target), {}),
        ("Cloud distance >=0.5 km and ST_QA <=2 K", "quadratic",
         local_sample(observations(data, quality="2K"), target), {}),
    ]:
        estimate = fit_contrasts(observations_, target, model=model, draws=SENSITIVITY_DRAWS, **options)
        comparisons.append(variant_row(label, model, target, estimate))
    alternate = Contrast(0., target.lower_built_pct, target.higher_built_pct)
    estimate = fit_contrasts(local_sample(primary, alternate), alternate, draws=SENSITIVITY_DRAWS)
    comparisons.append(variant_row("Tree cover 0% to 10%", "quadratic", alternate, estimate))

    by_date = []
    for date in sorted(inputs.date.unique()):
        fit = fit_contrasts(inputs[inputs.date != date], target, draws=SENSITIVITY_DRAWS)
        by_date.append({"omitted_date": date, **variant_row("Primary quadratic interaction", "quadratic", target, fit)})
    by_block = []
    for e, n in inputs[BLOCK].drop_duplicates().itertuples(index=False, name=None):
        removed = inputs.block_e.eq(e) & inputs.block_n.eq(n)
        fit = fit_contrasts(inputs[~removed], target, draws=0)
        by_block.append({"omitted_block_e": e, "omitted_block_n": n, "omitted_cell_dates": int(removed.sum()),
                         **variant_row("Primary quadratic interaction", "quadratic", target, fit)})

    distribution = []
    for cohort, d in [("band-supported sample", inputs), ("shared block-dates", shared)]:
        for band, g in d.groupby("band"):
            one = g.drop_duplicates("cell_id")
            distribution.append({
                "cohort": cohort, "band": band, "cell_dates": len(g), "distinct_cells": len(one),
                "blocks": len(g[BLOCK].drop_duplicates()), "dates": int(g.date.nunique()),
                "mean_tree_pct": 100 * one.tree_fraction.mean(), "mean_built_pct": 100 * one.built_fraction.mean(),
            })
    pair_checks = {
        "same_dates": pair_check(strict, shared=False), "same_block_dates": pair_check(strict, shared=True),
    }
    output.mkdir(parents=True, exist_ok=True)
    write_csv(pd.DataFrame(comparisons), output / "robustness.csv")
    write_csv(pd.DataFrame(by_date), output / "leave_one_date_out.csv")
    write_csv(pd.DataFrame(by_block), output / "leave_one_block_out.csv")
    write_csv(pd.DataFrame(distribution), output / "sample_comparison.csv")
    write_csv(strict, output / "paired_temperatures.csv")
    result = {
        "design_sha256": checksum(design / "design_record.json"), "bootstrap_seed": SEED,
        "contrast": design_info["chosen"], "additive": additive, "interaction": interaction, "paired": pair_checks,
    }
    (output / "estimates.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--panels", default=PANELS, type=Path)
    parser.add_argument("--design", default=DESIGN_DIR, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    run(args.panels, args.design, args.output)
