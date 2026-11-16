"""Freeze the nine Stage 8 pairs, their temperature gaps and their date weights."""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from analyze_landsat import write_csv
from dense_swap_model import CHANGE_PP
from fetch_inputs import checksum
from prepare_grass_tree_substitution import SIDES
from prepare_influential_pair_audit import IDENTITY
from prepare_local_comparison import BLOCK
from prepare_overpass_validation import verify_hashes
from prepare_quality_dense_swap import OUTPUT as STAGE8


OUTPUT = Path("results/quality_dense_influence")
DESIGN = OUTPUT / "design"
MATCHES = ("dense", "wider_balance")
BASELINE_SAMPLE = "qa4_omit_both"


def date_weights(rows):
    raw = 1 / rows.groupby("date").date.transform("size").to_numpy()
    return raw / raw.sum()


def frozen_inputs():
    design_path, analysis = STAGE8 / "design", STAGE8 / "analysis"
    record = json.loads((design_path / "design_record.json").read_text())
    run = json.loads((analysis / "analysis_record.json").read_text())
    verify_hashes(design_path, record["design_sha256"])
    verify_hashes(design_path, {"design_record.json": run["design_sha256"]})
    verify_hashes(analysis, run["output_sha256"])
    outcomes = pd.read_csv(analysis / "qa4_paired_temperatures.csv", float_precision="round_trip")
    candidates = pd.read_csv(design_path / "pair_coverage.csv")
    cover_keys = [*IDENTITY, *BLOCK]
    common = []
    for label in MATCHES:
        used = outcomes.loc[outcomes.case.eq(label) & outcomes.passes_qa4_and_outside]
        used = used.drop(columns=["case", "passes_qa4_and_outside"]).sort_values(IDENTITY).reset_index(drop=True)
        keys = candidates.loc[candidates.case.eq(label) & candidates.passes_qa4_and_outside, cover_keys]
        pd.testing.assert_frame_equal(used[cover_keys], keys.sort_values(IDENTITY).reset_index(drop=True))
        common.append(used)
    pd.testing.assert_frame_equal(*common)
    pairs = common[0]
    if (len(pairs) != 9 or pairs.duplicated(IDENTITY).any() or pairs.date.nunique() != 6
            or len(pairs[BLOCK].drop_duplicates()) != 6 or (pairs.tree_gap_pp <= 0).any()):
        raise ValueError("The frozen nine-pair cohort differs from the Stage 8 result")
    if (pairs[["lower_qa4_lst_share", "higher_qa4_lst_share"]] < .6).any().any():
        raise ValueError("A selected endpoint falls below the Stage 8 coverage floor")
    estimates = pd.read_csv(analysis / "estimates.csv")
    baseline = estimates[estimates["sample"].eq(BASELINE_SAMPLE) & estimates.model.eq("matched")]
    metrics = ["estimate_c_per_5pp", "lower_c", "upper_c"]
    np.testing.assert_array_equal(baseline[metrics].iloc[0], baseline[metrics].iloc[1])
    return pairs, baseline.iloc[0][metrics].to_dict()


def prepare(output):
    pairs, baseline = frozen_inputs()
    native = pairs.higher_qa4_lst_c - pairs.lower_qa4_lst_c
    standardized = CHANGE_PP * native / pairs.tree_gap_pp
    np.testing.assert_allclose(pairs.lst_gap_c, native, rtol=0, atol=2e-9)
    np.testing.assert_allclose(pairs.lst_c_per_5pp, standardized, rtol=0, atol=2e-9)
    pairs.insert(0, "pair_id", np.arange(1, len(pairs) + 1))
    pairs["pairs_on_date"] = pairs.groupby("date").date.transform("size")
    pairs["weight"] = date_weights(pairs)
    pairs["weighted_contribution_c"] = pairs.weight * standardized
    np.testing.assert_allclose(pairs.weight.sum(), 1, rtol=0, atol=1e-12)
    np.testing.assert_allclose(pairs.weighted_contribution_c.sum(), baseline["estimate_c_per_5pp"], rtol=0, atol=1.5e-9)
    output.mkdir(parents=True, exist_ok=True)
    write_csv(pairs, output / "weighted_pairs.csv")
    record = {
        "stage8_design_sha256": checksum(STAGE8 / "design" / "design_record.json"),
        "stage8_analysis_sha256": checksum(STAGE8 / "analysis" / "analysis_record.json"),
        "pairs": len(pairs),
        "block_ids": pairs[BLOCK].drop_duplicates().sort_values(BLOCK).values.tolist(),
        "dates": sorted(pairs.date.unique()),
        "stage8_baseline": baseline,
        "design_sha256": {"weighted_pairs.csv": checksum(output / "weighted_pairs.csv")},
    }
    (output / "design_record.json").write_text(json.dumps(record, indent=2) + "\n")
    print(pairs[["pair_id", "date", *BLOCK, *SIDES, "lst_gap_c", "lst_c_per_5pp",
                 "weight", "weighted_contribution_c"]].to_string(index=False))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    prepare(args.output)
