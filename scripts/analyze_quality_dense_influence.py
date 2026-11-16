"""Recompute the matched estimate and spatial interval after one deletion at a time."""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from analyze_landsat import interval, write_csv
from dense_swap_model import CHANGE_PP, DRAWS, SEED, spatial_weights
from fetch_inputs import checksum
from prepare_grass_tree_substitution import SIDES
from prepare_local_comparison import BLOCK
from prepare_overpass_validation import verify_hashes
from prepare_quality_dense_influence import DESIGN, STAGE8, date_weights


def matched_effect(rows):
    blocks, membership = np.unique(rows[BLOCK], axis=0, return_inverse=True)
    summary = {
        "pairs": len(rows), "distinct_spatial_pairs": len(rows[SIDES].drop_duplicates()),
        "dates": int(rows.date.nunique()), "blocks": len(blocks),
    }
    if summary["dates"] < 2 or summary["blocks"] < 3:
        return {**summary, "status": "fewer than two dates or three blocks"}
    w = date_weights(rows)
    temperature = CHANGE_PP * rows.lst_gap_c.to_numpy() / rows.tree_gap_pp.to_numpy()
    weights = spatial_weights(w, membership, len(blocks))
    draws = (weights @ temperature) / weights.sum(axis=1)
    lo, hi = interval(draws)
    return {
        **summary, "status": "estimated", "estimate_c_per_5pp": float(w @ temperature),
        "lower_c": float(lo), "upper_c": float(hi), "interval_includes_zero": bool(lo <= 0 <= hi),
        "draws_retained": len(draws),
    }


def influence(pairs, original):
    baseline = matched_effect(pairs)
    fields = ["estimate_c_per_5pp", "lower_c", "upper_c"]
    np.testing.assert_allclose([baseline[k] for k in fields], [original[k] for k in fields], rtol=0, atol=2e-9)
    records = [{
        "omission_type": "none", "omitted_date": None, "omitted_block_e": None, "omitted_block_n": None,
        "omitted_pair_ids": "", "omitted_pairs": 0, "original_omitted_weight": 0.,
        "original_omitted_contribution_c": 0., "change_c_per_5pp": 0., **baseline,
    }]
    for tag, group in (("date", "date"), ("block", BLOCK)):
        for key, omitted in pairs.groupby(group, sort=True):
            retained = pairs.drop(index=omitted.index)
            result = matched_effect(retained)
            change = result.get("estimate_c_per_5pp")
            e, n = key if tag == "block" else (None, None)
            records.append({
                "omission_type": tag, "omitted_date": key if tag == "date" else None,
                "omitted_block_e": e, "omitted_block_n": n,
                "omitted_pair_ids": " ".join(omitted.pair_id.astype(str)),
                "omitted_pairs": len(omitted), "original_omitted_weight": float(omitted.weight.sum()),
                "original_omitted_contribution_c": float(omitted.weighted_contribution_c.sum()),
                "change_c_per_5pp": change - baseline["estimate_c_per_5pp"] if change is not None else None,
                **result,
            })
    return pd.DataFrame(records).astype({"omitted_block_e": "Int64", "omitted_block_n": "Int64"})


def analyze(design, output):
    record = json.loads((design / "design_record.json").read_text())
    verify_hashes(design, record["design_sha256"])
    verify_hashes(STAGE8 / "design", {"design_record.json": record["stage8_design_sha256"]})
    verify_hashes(STAGE8 / "analysis", {"analysis_record.json": record["stage8_analysis_sha256"]})
    pairs = pd.read_csv(design / "weighted_pairs.csv", float_precision="round_trip")
    np.testing.assert_allclose(pairs.weight, date_weights(pairs), rtol=0, atol=5.1e-10)
    table = influence(pairs, record["stage8_baseline"])
    output.mkdir(parents=True, exist_ok=True)
    write_csv(table, output / "leave_one_out.csv")
    summary = {
        "design_sha256": checksum(design / "design_record.json"),
        "seed": SEED, "draws": DRAWS,
        "output_sha256": {"leave_one_out.csv": checksum(output / "leave_one_out.csv")},
    }
    (output / "analysis_record.json").write_text(json.dumps(summary, indent=2) + "\n")
    cols = ["omission_type", "omitted_date", "omitted_block_e", "omitted_block_n",
            "omitted_pairs", "dates", "blocks", "estimate_c_per_5pp", "lower_c", "upper_c", "change_c_per_5pp"]
    print(table[cols].to_string(index=False))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--design", type=Path, default=DESIGN)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    analyze(args.design, args.output)
