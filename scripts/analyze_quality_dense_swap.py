"""Fit supported five-point contrasts on the frozen 4 K paired sample."""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from analyze_dense_tree_swap import MODELS, fitted_contrasts
from analyze_influential_pairs import paired_changes, sample_row
from analyze_landsat import write_csv
from analyze_overpass_validation import load_cohorts
from dense_swap_model import DRAWS, SEED, assess
from fetch_inputs import checksum
from make_validation_panel import SNAPSHOTS
from prepare_grass_tree_substitution import SIDES, VALIDATION
from prepare_influential_pair_audit import IDENTITY, STAGE6, DESIGN as STAGE7
from prepare_local_comparison import BLOCK, KEY
from prepare_overpass_validation import verify_hashes
from prepare_quality_dense_swap import COUNTS, DESIGN, SAMPLES, verify_raw


TEMPERATURE = "qa4_lst_c"


def temperatures(inventory, replay):
    validation = json.loads((VALIDATION / "design_record.json").read_text())
    verify_hashes(VALIDATION, validation["design_sha256"])
    data = pd.concat(load_cohorts(SNAPSHOTS, VALIDATION, validation).values(), ignore_index=True)
    values = inventory.merge(data[[*KEY, "scene_id", *COUNTS, TEMPERATURE]],
                             on=[*KEY, "scene_id", *COUNTS], validate="one_to_one")
    if len(values) != len(inventory) or not np.isfinite(values[TEMPERATURE]).all():
        raise ValueError("The frozen quality endpoint inventory does not match the archived temperatures")
    if replay:
        verify_raw(values)
    return values[[*KEY, TEMPERATURE]]


def join_pair_outcomes(pairs, values):
    frame = pairs.copy()
    for side, field in zip(("lower", "higher"), SIDES):
        frame = frame.merge(values.rename(columns={"cell_id": field, TEMPERATURE: f"{side}_qa4_lst_c"}),
                            on=["date", field], how="left", validate="one_to_one")
    cols = ["lower_qa4_lst_c", "higher_qa4_lst_c"]
    if not np.isfinite(frame[cols].to_numpy()).all():
        raise ValueError("The QA4 pair list must retain a qualifying temperature at both endpoints")
    frame["lst_gap_c"] = frame.higher_qa4_lst_c - frame.lower_qa4_lst_c
    frame["lst_c_per_5pp"] = frame.lst_gap_c * 5 / frame.tree_gap_pp
    return frame


def analyze(design_path, output, replay):
    record = json.loads((design_path / "design_record.json").read_text())
    verify_hashes(design_path, record["design_sha256"])
    for parent, key in ((STAGE6, "stage6"), (STAGE7, "stage7"), (VALIDATION, "stage3")):
        verify_hashes(parent, {"design_record.json": record[f"{key}_design_sha256"]})
    screening = pd.read_csv(design_path / "pair_coverage.csv", float_precision="round_trip")
    endpoints = pd.read_csv(design_path / "endpoint_coverage.csv")
    output.mkdir(parents=True, exist_ok=True)
    estimates, changes, outputs = [], [], []
    if any(sample["matched_ready"] for case in record["support"].values()
           for name, sample in case.items() if name in SAMPLES):
        data = temperatures(endpoints[endpoints.passes_qa4], replay)
    else:
        data = pd.DataFrame(columns=[*KEY, TEMPERATURE])
    for case, original in screening.groupby("case", sort=True):
        kept = original[original.passes_qa4].copy()
        results = {}
        if record["support"][case]["qa4_full"]["matched_ready"]:
            kept = join_pair_outcomes(kept, data)
            outputs.append(kept)
        for label, flag in SAMPLES.items():
            selection = kept[kept[flag]]
            support = assess(selection)
            if support != record["support"][case][label]:
                raise ValueError(f"Frozen support changed for {case} {label}")
            result, coords = {}, selection[BLOCK].drop_duplicates()
            if support["matched_ready"]:
                result, _, coords = fitted_contrasts(selection, support)
            estimates.extend(sample_row(label, model, support, coords, result, case) for model in MODELS)
            results[label] = result
        first, second = (results[name] for name in SAMPLES)
        common = [model for model in MODELS if model in first and model in second]
        if common:
            changes.extend({"case": case, **row} for row in paired_changes(
                kept, kept.passes_qa4_and_outside.to_numpy(), first, second, common
            ))

    pair_columns = ["case", *IDENTITY, *BLOCK, "tree_gap_pp", "rangeland_gap_pp", "built_gap_pp",
                    "lower_qa4_lst_share", "higher_qa4_lst_share", "lower_qa4_lst_c", "higher_qa4_lst_c",
                    "lst_gap_c", "lst_c_per_5pp", "passes_qa4_and_outside"]
    change_columns = ["case", "model", "change_c_per_5pp", "lower_c", "upper_c", "paired_draws_retained",
                      "omit_with_original_date_weights_c"]
    tables = {
        "estimates.csv": pd.DataFrame(estimates),
        "omission_changes.csv": pd.DataFrame(changes, columns=change_columns),
        "qa4_paired_temperatures.csv": pd.concat(outputs, ignore_index=True)[pair_columns]
        if outputs else pd.DataFrame(columns=pair_columns),
    }
    for name, table in tables.items():
        write_csv(table, output / name)
    metadata = {
        "design_sha256": checksum(design_path / "design_record.json"),
        "draws": DRAWS, "seed": SEED,
        "output_sha256": {name: checksum(output / name) for name in tables},
    }
    (output / "analysis_record.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(tables["estimates.csv"].to_string(index=False))
    print(tables["omission_changes.csv"].to_string(index=False))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--design", type=Path, default=DESIGN)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--verify-pixels", action="store_true")
    args = parser.parse_args()
    analyze(args.design, args.output, args.verify_pixels)
