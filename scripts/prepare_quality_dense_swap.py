"""Freeze the dense-area pair sample using only the specified Landsat QA coverage."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from rasterio.transform import array_bounds

from analyze_influential_pairs import masks
from analyze_landsat import write_csv
from audit_heldout_overpasses import AUDIT, land_cover
from dense_swap_model import assess
from fetch_inputs import RAW, checksum
from make_influential_pair_chips import extent, grid_window, thermal_source
from make_landsat_panel import SOURCE_FILES, points_within_native
from make_validation_panel import SNAPSHOTS
from prepare_grass_tree_substitution import SIDES, VALIDATION
from prepare_influential_pair_audit import (
    DESIGN as STAGE7, DISCOVERY, IDENTITY, STAGE6, VALIDATION_MANIFEST, read_pair_inputs,
)
from prepare_local_comparison import BLOCK, KEY, PANELS
from prepare_overpass_validation import verify_hashes
from screen_landsat import land_layers


OUTPUT = Path("results/quality_dense_swap")
DESIGN = OUTPUT / "design"
COUNTS = ["land_centers", "valid_lst_centers", "qa4_lst_centers"]
SAMPLES = {"qa4_full": "passes_qa4", "qa4_omit_both": "passes_qa4_and_outside"}
FLAG_COLUMNS = ["passes_qa4", "in_omitted_blocks", "passes_qa4_and_outside"]


def read_coverage():
    columns = [*KEY, "scene_id", *COUNTS, "qa4_lst_share"]
    original = PANELS / "world_500.csv.gz"
    baseline = json.loads((PANELS / "provenance.json").read_text())
    with gzip.open(original, "rb") as data:
        digest = hashlib.sha256(data.read()).hexdigest()
    if digest != baseline["panel_sha256"]["world_500"]:
        raise ValueError("Stage 3 discovery panel contents differ from the archived inputs")
    val = json.loads((VALIDATION / "design_record.json").read_text())
    verify_hashes(VALIDATION, val["design_sha256"])
    parts = [pd.read_csv(original, usecols=columns)]
    for year in val["validation_years"]:
        folder = SNAPSHOTS / str(year)
        metadata = json.loads((folder / "provenance.json").read_text())
        verify_hashes(VALIDATION, {"design_record.json": metadata["design_sha256"]})
        verify_hashes(folder, {"world_500.csv.gz": metadata["panel_sha256"]})
        parts.append(pd.read_csv(folder / "world_500.csv.gz", usecols=columns))
    panel = pd.concat(parts, ignore_index=True)
    expected = pd.read_csv(VALIDATION / "eligible_cell_dates.csv", usecols=[*KEY, "scene_id", "valid_lst_centers"])
    return expected.merge(panel, on=[*KEY, "scene_id", "valid_lst_centers"], validate="one_to_one")


def coverage_ok(qa4_count, land_count):
    return (land_count > 0) & (5 * qa4_count >= 3 * land_count)


def endpoint_inventory(populations, coverage):
    chosen = pd.concat([
        pairs[["date", side]].rename(columns={side: "cell_id"})
        for pairs in populations.values() for side in SIDES
    ], ignore_index=True).drop_duplicates().sort_values(KEY)
    endpoints = chosen.merge(coverage, on=KEY, how="left", validate="one_to_one")
    n, m, q = (endpoints[key].to_numpy() for key in COUNTS)
    if not ((n > 0) & (q >= 0) & (q <= m) & (m <= n) & coverage_ok(m, n)).all():
        raise ValueError("A pair endpoint does not match its originally eligible coverage record")
    np.testing.assert_allclose(endpoints.qa4_lst_share, q / n, rtol=0, atol=5.1e-9)
    endpoints["passes_qa4"] = coverage_ok(q, n)
    return endpoints


def screen(pairs, endpoints, omitted):
    data = pairs.copy()
    value_cols = [*COUNTS, "qa4_lst_share", "passes_qa4"]
    for side, key in zip(("lower", "higher"), SIDES):
        values = endpoints[[*KEY, *value_cols]].rename(columns={
            "cell_id": key, **{name: f"{side}_{name}" for name in value_cols}
        })
        data = data.merge(values, on=["date", key], validate="one_to_one")
    if len(data) != len(pairs):
        raise ValueError("Pair screening must include the complete frozen pair list")
    data["passes_qa4"] = data.lower_passes_qa4 & data.higher_passes_qa4
    data["in_omitted_blocks"] = pd.MultiIndex.from_frame(data[BLOCK]).isin(pd.MultiIndex.from_frame(omitted))
    data["passes_qa4_and_outside"] = data.passes_qa4 & ~data.in_omitted_blocks
    return data


def count_support(pairs):
    return {
        "pairs": len(pairs), "spatial_pairs": len(pairs[SIDES].drop_duplicates()),
        "dates": int(pairs.date.nunique()), "blocks": len(pairs[BLOCK].drop_duplicates()),
    }


def verify_raw(endpoints):
    """Replay the earlier area weighting against the fixed native Landsat rasters."""
    sources = pd.concat([pd.read_csv(DISCOVERY), pd.read_csv(VALIDATION_MANIFEST)], ignore_index=True)
    sources = sources.assign(date=sources.datetime_utc.str[:10]).set_index("date")
    for year, panel in endpoints.groupby(endpoints.date.str[:4].astype(int)):
        reference = json.loads((AUDIT / str(year) / "audit_record.json").read_text())
        esri, _ = land_cover(year)
        verify_hashes(RAW, reference["spatial_sha256"])
        grid, land, _, _, _ = land_layers(trim_to_30m=False, esri_path=esri)
        bounds = array_bounds(*land.shape, grid)
        for date, focus in panel.groupby("date", sort=True):
            source = sources.loc[date].to_dict()
            with np.load(thermal_source(source, year, bounds), allow_pickle=False) as raw:
                if raw["urls"].tolist() != [source[key] for key in SOURCE_FILES]:
                    raise ValueError(f"Native Landsat source URLs changed for {date}")
                thermal = {key: raw[key] for key in SOURCE_FILES}
                frame = raw["source_transform"]
            for row in focus.itertuples(index=False):
                win, affine = grid_window(extent([row.cell_id]), grid)
                region = land[win]
                y, x = np.where(region)
                indices = points_within_native(affine.c + (x + .5) * affine.a, affine.f + (y + .5) * affine.e,
                                               frame, thermal["lwir11"].shape)
                (valid, qa4, _), celsius = masks(*(thermal[key][indices] for key in SOURCE_FILES))
                if (source["scene_id"], len(x), int(valid.sum()), int(qa4.sum())) != (
                        row.scene_id, row.land_centers, row.valid_lst_centers, row.qa4_lst_centers):
                    raise ValueError(f"Native QA counts differ from the archived panel for {date} {row.cell_id}")
                if hasattr(row, "qa4_lst_c"):
                    np.testing.assert_allclose(celsius[qa4].mean(), row.qa4_lst_c, rtol=0, atol=5.1e-9)
    print(f"Verified native Landsat retrievals for {len(endpoints)} distinct cell-dates.")


def prepare(output, replay):
    previous = json.loads((STAGE7 / "design_record.json").read_text())
    verify_hashes(STAGE7, previous["files_sha256"])
    verify_hashes(STAGE6, {"design_record.json": previous["stage6_design_sha256"]})
    verify_hashes(VALIDATION, {"design_record.json": previous["stage3_design_sha256"]})
    populations = read_pair_inputs(STAGE6)
    endpoints = endpoint_inventory(populations, read_coverage())
    if replay:
        verify_raw(endpoints)
    omitted = pd.DataFrame(previous["omitted_blocks"], columns=BLOCK)
    screens = pd.concat([
        screen(pairs, endpoints, omitted).assign(case=label) for label, pairs in populations.items()
    ], ignore_index=True)
    output.mkdir(parents=True, exist_ok=True)
    write_csv(endpoints, output / "endpoint_coverage.csv")
    write_csv(screens, output / "pair_coverage.csv")
    screens = pd.read_csv(output / "pair_coverage.csv", float_precision="round_trip")
    states, daily, detail = {}, [], []
    for case, rows in screens.groupby("case", sort=True):
        states[case] = {}
        for label, flag in (("original", None), *SAMPLES.items()):
            chosen = rows if flag is None else rows[rows[flag]]
            support = count_support(chosen) if flag is None else assess(chosen)
            states[case][label] = support
            detail.append({
                "case": case, "sample": label, **{key: support[key] for key in ("pairs", "spatial_pairs", "dates", "blocks")},
            })
        for date, visit in rows.groupby("date", sort=True):
            daily.append({
                "case": case, "date": date, "original_pairs": len(visit),
                **{f"{label}_pairs": int(visit[flag].sum()) for label, flag in SAMPLES.items()},
            })
    write_csv(pd.DataFrame(daily), output / "support_by_date.csv")
    write_csv(pd.DataFrame(detail), output / "support_summary.csv")
    info = {
        "stage7_design_sha256": checksum(STAGE7 / "design_record.json"),
        "stage6_design_sha256": previous["stage6_design_sha256"],
        "stage3_design_sha256": previous["stage3_design_sha256"],
        "omitted_blocks": previous["omitted_blocks"],
        "qa_rule": {"max_uncertainty_k": 4, "minimum_cloud_distance_km": .5, "minimum_land_coverage_per_endpoint": .6},
        "unique_cell_dates": len(endpoints), "unique_pair_dates": len(screens[IDENTITY].drop_duplicates()),
        "support": states,
        "design_sha256": {name: checksum(output / name) for name in (
            "endpoint_coverage.csv", "pair_coverage.csv", "support_by_date.csv", "support_summary.csv"
        )},
    }
    (output / "design_record.json").write_text(json.dumps(info, indent=2) + "\n")
    print(json.dumps(states, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--verify-pixels", action="store_true")
    args = parser.parse_args()
    prepare(args.output, args.verify_pixels)
