"""Estimate supported five-point paired contrasts on the frozen dense-area sample."""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from analyze_grass_tree_substitution import join_pairs
from analyze_landsat import interval, write_csv
from analyze_overpass_validation import load_cohorts
from dense_swap_model import CHANGE_PP, DRAWS, SEED, adjusted_fits, assess, design, full_rank, spatial_weights
from fetch_inputs import checksum
from make_validation_panel import SNAPSHOTS
from prepare_dense_tree_swap import DESIGN, ESTIMATE_CASES, STAGE5, VALIDATION
from prepare_grass_tree_substitution import GAPS, SIDES
from prepare_local_comparison import BLOCK
from prepare_overpass_validation import verify_hashes


MODELS = ("matched", "adjusted")


def fitted_contrasts(pairs, model_support):
    x, w, blocks, ids, _ = design(pairs)
    y = pairs.lst_gap_c.to_numpy()
    normalized = y * CHANGE_PP / pairs.tree_gap_pp.to_numpy()
    weights = spatial_weights(w, ids, len(blocks))
    omit = w * (ids != np.arange(len(blocks))[:, None])
    fits, omitted = {}, {}
    if model_support["adjusted_ready"]:
        keep = full_rank(x, weights)
        if keep.sum() != model_support["adjustment_draws_retained"]:
            raise ValueError("Model rank differs from the frozen pair design")
        weights = weights[keep]
        fits["adjusted"] = adjusted_fits(x, y, np.vstack([w, weights]))
        omitted["adjusted"] = adjusted_fits(x, y, omit)
    fits["matched"] = np.r_[np.average(normalized, weights=w), (weights @ normalized) / weights.sum(axis=1)]
    omitted["matched"] = (omit @ normalized) / omit.sum(axis=1)
    return fits, omitted, blocks


def balance(pairs):
    w = 1 / pairs.groupby("date").date.transform("size").to_numpy()
    means = {f"mean_{name}": float(np.average(pairs[name], weights=w))
             for name in ("distance_m", *GAPS)}
    absolute = {f"mean_abs_{name}": float(np.average(pairs[name].abs(), weights=w)) for name in GAPS[2:]}
    built = 100 * pairs[["lower_built_fraction", "higher_built_fraction"]].to_numpy()
    return {**means, **absolute, "minimum_built_pct": float(built.min()), "maximum_built_pct": float(built.max())}


def overlap_diagnostic(temperatures):
    columns = ["case", "pair_subset", "pairs", "blocks", "weight_share", "matched_mean_c", "contribution_c"]
    identity = ["date", *SIDES]
    primary, wider = (temperatures[temperatures.case.eq(case)] for case in ("dense", "wider_balance"))
    if primary.empty or wider.empty:
        return pd.DataFrame(columns=columns)
    rows = []
    for sample, other in ((primary, wider), (wider, primary)):
        shared = pd.MultiIndex.from_frame(sample[identity]).isin(pd.MultiIndex.from_frame(other[identity]))
        w = 1 / sample.groupby("date").date.transform("size").to_numpy()
        w /= w.sum()
        for label, mask in (("shared", shared), ("only this design", ~shared)):
            share = w[mask].sum()
            contribution = w[mask] @ sample.loc[mask, "lst_c_per_5pp"].to_numpy()
            rows.append({
                "case": sample.case.iloc[0], "pair_subset": label, "pairs": int(mask.sum()),
                "blocks": len(sample.loc[mask, BLOCK].drop_duplicates()), "weight_share": share,
                "matched_mean_c": contribution / share if share > 0 else None, "contribution_c": contribution,
            })
    return pd.DataFrame(rows, columns=columns)


def analyze(dense_design, stage5, snapshots, output):
    record = json.loads((dense_design / "design_record.json").read_text())
    verify_hashes(dense_design, record["design_sha256"])
    verify_hashes(stage5, {"design_record.json": record["stage5_design_sha256"]})
    stage5_record = json.loads((stage5 / "design_record.json").read_text())
    verify_hashes(stage5, stage5_record["design_sha256"])
    verify_hashes(VALIDATION, {"design_record.json": stage5_record["validation_design_sha256"]})
    validation = json.loads((VALIDATION / "design_record.json").read_text())
    verify_hashes(VALIDATION, validation["design_sha256"])
    cover = pd.read_csv(stage5 / "cell_covariates.csv", float_precision="round_trip")
    registry = pd.read_csv(dense_design / "pairs.csv", float_precision="round_trip")
    temps = None
    if any(case["matched_ready"] for case in record["estimation"].values()):
        temps = pd.concat(load_cohorts(snapshots, VALIDATION, validation).values(), ignore_index=True)

    estimates, diagnostics, influence, fit_pairs = [], [], [], []
    for case, (comparison, band) in ESTIMATE_CASES.items():
        pairs = registry[registry.comparison.eq(comparison) & registry.built_band.eq(band)]
        expected = record["estimation"][case]
        support = assess(pairs)
        if {"comparison": comparison, "built_band": band, **support} != expected:
            raise ValueError(f"Model support no longer agrees with the outcome-free check for {case}")
        common = {"case": case, "comparison": comparison, "built_band": band, **{
            key: support.get(key) for key in ("pairs", "spatial_pairs", "dates", "blocks", "effective_blocks")
        }}
        fits = {}
        if support["matched_ready"]:
            pairs = join_pairs(pairs, temps, cover)
            fits, omitted, blocks = fitted_contrasts(pairs, support)
            fit_pairs.append(pairs[["date", *SIDES, *BLOCK, "lower_response_c", "higher_response_c", "lst_gap_c"]]
                             .assign(case=case, lst_c_per_5pp=pairs.lst_gap_c * CHANGE_PP / pairs.tree_gap_pp))
            diagnostics.append({**common, **balance(pairs)})
            for model, values in omitted.items():
                for (e, n), value in zip(blocks, values):
                    influence.append({"case": case, "model": model, "omitted_e": e, "omitted_n": n,
                                      "estimate_c_per_5pp": value})
        for model in MODELS:
            values = fits.get(model)
            lo, hi = interval(values[1:]) if values is not None else (None, None)
            estimates.append({
                **common, "model": model, "status": "estimated" if values is not None else support["reason"],
                "estimate_c_per_5pp": values[0] if values is not None else None, "lower_c": lo, "upper_c": hi,
                "draws_requested": DRAWS if values is not None else None,
                "draws_retained": len(values) - 1 if values is not None else None,
            })

    output.mkdir(parents=True, exist_ok=True)
    pair_columns = ["date", *SIDES, *BLOCK, "lower_response_c", "higher_response_c", "lst_gap_c",
                    "case", "lst_c_per_5pp"]
    temperatures = pd.concat(fit_pairs) if fit_pairs else pd.DataFrame(columns=pair_columns)
    tables = {
        "estimates.csv": pd.DataFrame(estimates), "pair_balance.csv": pd.DataFrame(diagnostics),
        "leave_one_block_out.csv": pd.DataFrame(influence), "paired_temperatures.csv": temperatures,
        "pair_overlap.csv": overlap_diagnostic(temperatures),
    }
    for name, table in tables.items():
        write_csv(table, output / name)
    details = {
        "dense_design_sha256": checksum(dense_design / "design_record.json"),
        "seed": SEED, "exchange_points": CHANGE_PP,
        "output_sha256": {name: checksum(output / name) for name in tables},
    }
    (output / "analysis_record.json").write_text(json.dumps(details, indent=2) + "\n")
    print(tables["estimates.csv"].to_string(index=False))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--design", type=Path, default=DESIGN)
    parser.add_argument("--stage5", type=Path, default=STAGE5)
    parser.add_argument("--snapshots", type=Path, default=SNAPSHOTS)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    analyze(args.design, args.stage5, args.snapshots, args.output)
