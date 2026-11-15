"""Freeze Stage 3 cell-date IDs, annual cover checks, and matching registries."""
import argparse
import json
from pathlib import Path

import pandas as pd

from analyze_landsat import write_csv
from audit_heldout_overpasses import AUDIT
from fetch_inputs import checksum
from prepare_built_support import DESIGN_DIR as BUILT_DESIGN, STRATUM, selected_samples
from prepare_local_comparison import BLOCK, COORDS, KEY, RULES, assign_blocks, build_matches


VALIDATION = Path("results/overpass_validation")
DESIGN = VALIDATION / "design"
VALIDATION_YEARS = (2023, 2024)
DISCOVERY_YEAR = 2025
YEARS = (*VALIDATION_YEARS, DISCOVERY_YEAR)
SAMPLES = ("full", "common", "common_stable")
COVER = ["tree_fraction", "built_fraction", "esri_water_share"]
MAX_COVER_SHIFT = 5


def verify_hashes(directory, expected):
    for name, digest in expected.items():
        if checksum(directory / name) != digest:
            raise ValueError(f"Frozen input changed at {directory / name}")


def load_audit(directory, year):
    folder = directory / str(year)
    summary = json.loads((folder / "audit_record.json").read_text())
    verify_hashes(folder, summary["result_sha256"])
    geometry = pd.read_csv(folder / "grid_land_cover.csv", float_precision="round_trip")
    cells = pd.read_csv(folder / "usable_cell_dates.csv.gz")
    scenes = pd.read_csv(folder / "scene_coverage.csv")
    use = scenes.prior_discovery if year == DISCOVERY_YEAR else scenes.scene_id.isin(summary["eligible_new_scene_ids"])
    selected = scenes.loc[use, ["date_utc", "scene_id", "overpass_sgt", "period"]].rename(columns={"date_utc": "date"})
    records = cells[cells.date.isin(selected.date)].merge(geometry, on="cell_id", validate="many_to_one")
    records = records.merge(selected, on="date", validate="many_to_one")
    if set(selected.date) != set(records.date) or records.duplicated(KEY).any():
        raise ValueError(f"Eligible cell dates are incompatible with scene selection for {year}")
    return geometry, records, selected, summary


def compare_annual_cover(grids, observations):
    keys = ["cell_id", *COORDS, "land_centers", "land_share"]
    changes = grids[DISCOVERY_YEAR][keys].copy()
    years = sorted(grids)
    for year in years:
        grid = grids[year]
        annual = grid.rename(columns={col: f"{col}_{year}" for col in grid if col not in keys})
        changes = changes.merge(annual, on=keys, how="inner", validate="one_to_one")
        if len(changes) != len(grid) or len(changes) != len(grids[DISCOVERY_YEAR]):
            raise ValueError("The WorldCover and boundary domain differs between years")
    spreads = []
    for name in COVER:
        cols = [f"{name}_{year}" for year in years]
        span = f"{name}_span_pp"
        changes[span] = 100 * (changes[cols].max(axis=1, skipna=False) - changes[cols].min(axis=1, skipna=False))
        spreads.append(span)
    availability = []
    for year in years:
        counts = observations[year].groupby("cell_id").date.nunique()
        col = f"usable_dates_{year}"
        changes[col] = changes.cell_id.map(counts).fillna(0).astype(int)
        availability.append(col)
    changes["common"] = changes[availability].gt(0).all(axis=1)
    changes["stable_cover"] = changes[spreads].le(MAX_COVER_SHIFT + 1e-6).all(axis=1)
    changes["common_stable"] = changes.common & changes.stable_cover
    changes["boundary_mixed"] = changes.land_share.lt(0.90)
    return changes


def support_summary(data, target):
    fit, arms, shared, strict = selected_samples(assign_blocks(data), target["tree_start_pct"],
                                                 target["lower_built_pct"], target["higher_built_pct"])
    counts = {"cells": int(fit.cell_id.nunique()), "cell_dates": len(fit),
              "dates_with_both_bands": int((fit.groupby("date").band.nunique() == 2).sum()),
              "shared_blocks": len({k[1:] for k in shared}), "shared_block_dates": len(shared),
              "strict_shared_blocks": len({k[1:] for k in strict}), "strict_shared_block_dates": len(strict)}
    for name, entry in arms.items():
        for field in ("endpoint_blocks", "endpoint_block_dates", "low_tight_endpoint_cells", "high_tight_endpoint_cells"):
            counts[f"{name}_{field}"] = entry[0][field]
    return fit[[*STRATUM, "cell_id", "band", "shared"]], counts


def prepare(audit, out):
    assessments = {year: load_audit(audit, year) for year in YEARS}
    validation_years = [year for year in VALIDATION_YEARS if not assessments[year][2].empty]
    if assessments[DISCOVERY_YEAR][3]["eligible_new_scene_ids"] or not validation_years:
        raise ValueError("The audited dates require a different validation scope")
    decoded = {year: assessments[year] for year in [*validation_years, DISCOVERY_YEAR]}
    geography = compare_annual_cover({year: decoded[year][0] for year in YEARS},
                                     {year: decoded[year][1] for year in YEARS})
    common_ids = geography.loc[geography.common, "cell_id"]
    stable_ids = geography.loc[geography.common_stable, "cell_id"]
    target = json.loads((BUILT_DESIGN / "design_record.json").read_text())["chosen"]
    summaries, support, pairs, members = [], [], [], []
    for year, (_, records, _, _) in decoded.items():
        records = records.assign(common=records.cell_id.isin(common_ids),
                                 common_stable=records.cell_id.isin(stable_ids))
        members.append(records[[*KEY, "scene_id", "common", "common_stable", "valid_lst_centers"]])
        for sample in SAMPLES:
            data = records if sample == "full" else records[records[sample]]
            s1_pairs, _, _, _ = build_matches(data, RULES["primary"])
            s2, counts = support_summary(data, target)
            pairs.append(s1_pairs.assign(year=year, sample=sample))
            support.append(s2.assign(year=year, sample=sample))
            summaries.append({
                "year": year, "sample": sample, "dates": int(data.date.nunique()), "cells": int(data.cell_id.nunique()),
                "cell_dates": len(data), "stage1_pairs": len(s1_pairs),
                "stage1_pair_blocks": len(s1_pairs[BLOCK].drop_duplicates()),
                **{f"stage2_{key}": value for key, value in counts.items()},
            })
    out.mkdir(parents=True, exist_ok=True)
    write_csv(geography, out / "annual_cover_comparison.csv")
    write_csv(pd.DataFrame(summaries), out / "sample_support.csv")
    write_csv(pd.concat(pairs, ignore_index=True), out / "stage1_pairs.csv")
    write_csv(pd.concat(support, ignore_index=True), out / "stage2_fit_sample.csv")
    write_csv(pd.concat(members, ignore_index=True), out / "eligible_cell_dates.csv")
    manifest = pd.concat([pd.read_csv(audit / str(y) / "eligible_new_scenes.csv") for y in validation_years],
                         ignore_index=True)
    write_csv(manifest, out / "validation_scenes.csv")
    details = {
        "validation_years": validation_years, "discovery_year": DISCOVERY_YEAR,
        "cover_stability_max_pp": MAX_COVER_SHIFT,
        "stage2_contrast": {key: target[key] for key in ("tree_start_pct", "tree_end_pct", "lower_built_pct", "higher_built_pct")},
        "scene_ids": manifest.scene_id.tolist(), "dates": manifest.datetime_utc.str[:10].tolist(),
        "discovery_panel_sha256": assessments[DISCOVERY_YEAR][3]["discovery_panel_sha256"],
        "audit_sha256": {str(year): checksum(audit / str(year) / "audit_record.json") for year in YEARS},
        "design_sha256": {p.name: checksum(p) for p in sorted(out.glob("*.csv"))},
    }
    (out / "design_record.json").write_text(json.dumps(details, indent=2) + "\n")
    print(pd.DataFrame(summaries).to_string(index=False))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--audit", type=Path, default=AUDIT)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    prepare(args.audit, args.output)
