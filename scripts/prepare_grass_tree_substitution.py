"""Freeze comparable tree-for-rangeland pair IDs without reading temperature values."""
import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import numpy as np
import pandas as pd

from analyze_landsat import write_csv
from audit_heldout_overpasses import AUDIT, BOUNDARY, land_cover
from fetch_inputs import RAW, WORLD_FILE, checksum
from grass_tree_water import ALTERNATE_AREA_KM2, MAP, MAP_SHA256, PRIMARY_AREA_KM2, distances, read_water
from prepare_local_comparison import (
    BLOCK, COORDS, DESIGN as LOCAL_DESIGN, KEY, assign_blocks, cover_composition, minimum_cost_pairs,
)
from prepare_overpass_validation import DESIGN as VALIDATION, YEARS, load_audit, verify_hashes


OUTPUT = Path("results/grass_tree_substitution")
DESIGN = OUTPUT / "design"
COVER = ["tree_fraction", "rangeland_fraction", "built_fraction", "other_fraction", "land_share"]
DISTANCES = ["coast_km", "inland_km", "inland_010_km"]
SIDES = ["lower_tree_cell", "higher_tree_cell"]
GAPS = ["tree_gap_pp", "rangeland_gap_pp", "built_gap_pp", "other_gap_pp", "land_gap_pp"]
PAIR_COLUMNS = ["date", *BLOCK, *SIDES, "distance_m", *GAPS, "cost"]
TARGET = 0.10
GAP_TOLERANCE = 0.02
LAND_TOLERANCE = 0.10
EPS = 1e-12


@dataclass(frozen=True)
class PairRule:
    distance_m: int
    built_other_gap: float
    combined_gap: float


RULES = {"primary": PairRule(2000, 0.05, 0.02), "tighter": PairRule(1500, 0.03, 0.015)}


def read_cover(geometry, esri_path):
    grid = geometry.assign(classified_centers=np.rint(geometry.classified_share * geometry.land_centers).astype(int))
    counts = cover_composition(grid, esri_path)
    if esri_path.name == "esri_impact_observatory_2025_48N_singapore.tif":
        archived = pd.read_csv(LOCAL_DESIGN / "land_composition.csv")
        pd.testing.assert_frame_equal(counts, archived)
    grid = grid[["cell_id", *COORDS, "land_centers", "land_share", "area_eligible"]].merge(
        counts, on=["cell_id", "land_centers"], validate="one_to_one"
    )
    grid = grid.loc[grid.area_eligible].drop(columns="area_eligible")
    denominator = grid.classified_centers
    if (denominator <= 0).any():
        raise ValueError("Eligible cells must contain classified study land")
    known = [grid[f"class_{code:02d}_centers"] for code in (2, 11, 7)]
    for name, count in zip(COVER[:4], [*known, denominator - sum(known)]):
        if (count < 0).any():
            raise ValueError("Invalid land-cover composition")
        grid[name] = count / denominator
    return grid


def design_inputs():
    """Read only the registered covariates and QA counts from the earlier design."""
    record = json.loads((VALIDATION / "design_record.json").read_text())
    verify_hashes(VALIDATION, record["design_sha256"])
    registry = pd.read_csv(VALIDATION / "eligible_cell_dates.csv")
    grids, dates = [], []
    for year in YEARS:
        geometry, usable, _, audit = load_audit(AUDIT, year)
        if checksum(AUDIT / str(year) / "audit_record.json") != record["audit_sha256"][str(year)]:
            raise ValueError(f"Stage 3 audit for {year} differs from the frozen design")
        esri_path, _ = land_cover(year)
        hashes = {p.name: checksum(p) for p in (esri_path, BOUNDARY, RAW / WORLD_FILE)}
        if hashes != audit["spatial_sha256"]:
            raise ValueError(f"Spatial inputs for {year} differ from the frozen audit")
        grids.append(read_cover(geometry, esri_path).assign(year=year))
        dates.append(usable[[*KEY, "scene_id", "valid_lst_centers"]])
    keys = [*KEY, "scene_id", "valid_lst_centers"]
    expected = pd.concat(dates, ignore_index=True).sort_values(KEY).reset_index(drop=True)
    pd.testing.assert_frame_equal(registry[keys].sort_values(KEY).reset_index(drop=True), expected)
    return pd.concat(grids, ignore_index=True), registry


def match_group(group, rule):
    cells = group.sort_values("cell_id", kind="stable").reset_index(drop=True)
    a, b = np.triu_indices(len(cells), k=1)
    x, y = cells[COORDS].to_numpy().T
    distances_m = np.hypot(x[a] - x[b], y[a] - y[b])
    shares = cells[COVER].to_numpy()
    lo, hi = np.where(shares[b, 0] > shares[a, 0], (a, b), (b, a))
    gaps = shares[hi] - shares[lo]
    tree, rangeland, built, other, land = gaps.T
    accept = (
        (distances_m <= rule.distance_m) & (shares[lo, 1] >= TARGET - EPS)
        & (np.abs(tree - TARGET) <= GAP_TOLERANCE + EPS)
        & (np.abs(rangeland + TARGET) <= GAP_TOLERANCE + EPS)
        & (np.abs(tree + rangeland) <= rule.combined_gap + EPS)
        & (np.abs(built) <= rule.built_other_gap + EPS)
        & (np.abs(other) <= rule.built_other_gap + EPS)
        & (np.abs(land) <= LAND_TOLERANCE + EPS)
    )
    cost = ((tree - TARGET) / GAP_TOLERANCE) ** 2 + ((rangeland + TARGET) / GAP_TOLERANCE) ** 2
    cost += (built / rule.built_other_gap) ** 2 + (other / rule.built_other_gap) ** 2
    cost += (land / LAND_TOLERANCE) ** 2 + (distances_m / rule.distance_m) ** 2
    valid = {(int(a[i]), int(b[i])): i for i in np.flatnonzero(accept)}
    edges = [(i, j, cost[k]) for (i, j), k in valid.items()]
    ids = cells.cell_id.to_numpy()
    records = []
    for pair in minimum_cost_pairs(edges):
        k = valid[pair]
        records.append((ids[lo[k]], ids[hi[k]], distances_m[k], *(100 * gaps[k]), cost[k]))
    degrees = np.bincount(np.r_[a[accept], b[accept]], minlength=len(cells))
    return records, pd.DataFrame({"cell_id": ids, "candidate_links": degrees})


def build_matches(eligible, rule):
    pairs, availability, daily = [], [], []
    for date, day in assign_blocks(eligible).groupby("date", sort=True):
        date_pairs, date_cells = [], []
        for block, group in day.groupby(BLOCK, sort=True):
            records, cells = match_group(group, rule)
            date_pairs.extend((date, *block, *record) for record in records)
            date_cells.append(cells.assign(date=date, block_e=block[0], block_n=block[1]))
        blocks = pd.concat(date_cells, ignore_index=True)
        supported = blocks.candidate_links.gt(0)
        daily.append({
            "date": date, "eligible_cells": len(day), "candidate_cells": int(supported.sum()),
            "candidate_links": int(blocks.candidate_links.sum() // 2),
            "candidate_blocks": len(blocks.loc[supported, BLOCK].drop_duplicates()),
            "pairs": len(date_pairs), "matched_cells": 2 * len(date_pairs),
        })
        pairs.extend(date_pairs)
        availability.append(blocks[[*KEY, "candidate_links"]])
    pairs = pd.DataFrame(pairs, columns=PAIR_COLUMNS)
    endpoints = pd.concat(
        [pairs[["date", side]].rename(columns={side: "cell_id"}) for side in SIDES], ignore_index=True
    ).assign(matched=True)
    support = pd.concat(availability, ignore_index=True).sort_values(KEY).merge(
        endpoints, on=KEY, how="left", validate="one_to_one"
    )
    support["matched"] = support.matched.eq(True)
    counts = pd.DataFrame(daily)
    fit_blocks = pairs.groupby("date")[BLOCK].apply(lambda g: len(g.drop_duplicates()))
    counts["matched_blocks"] = counts.date.map(fit_blocks).fillna(0).astype(int)
    return pairs, support, counts


def prepare(source, output):
    grid, registry = design_inputs()
    shore, components, inventory = read_water(source)
    cell_locations = grid[["cell_id", *COORDS]].drop_duplicates()
    if cell_locations.cell_id.duplicated().any():
        raise ValueError("A cell center changed between years")
    grid = grid.merge(distances(cell_locations, shore, components, inventory), on=["cell_id", *COORDS],
                      validate="many_to_one")
    eligible = registry[KEY].assign(year=registry.date.str[:4].astype(int)).merge(
        grid, on=["year", "cell_id"], how="left", validate="many_to_one"
    )
    if not np.isfinite(eligible[[*COVER, *DISTANCES]].to_numpy()).all():
        raise ValueError("All eligible cells must have complete composition and distance covariates")
    output.mkdir(parents=True, exist_ok=True)
    write_csv(inventory, output / "water_components.csv")
    registry_ids = eligible[["year", "cell_id"]].drop_duplicates()
    write_csv(grid.merge(registry_ids, on=["year", "cell_id"], validate="one_to_one")
              .sort_values(["year", "cell_id"]), output / "cell_covariates.csv")
    support, summaries, matching = eligible[KEY], [], {}
    for name, rule in RULES.items():
        paired, available, counts = build_matches(eligible, rule)
        write_csv(paired, output / f"pairs_{name}.csv")
        summaries.append(counts.assign(design=name))
        support = support.merge(available.rename(
            columns={"candidate_links": f"{name}_links", "matched": f"{name}_matched"}
        ), on=KEY, validate="one_to_one")
        matching[name] = {
            "rule": asdict(rule), "candidate_cells": int(available.loc[available.candidate_links > 0, "cell_id"].nunique()),
            "pairs": len(paired), "matched_cells": int(available.loc[available.matched, "cell_id"].nunique()),
            "matched_blocks": len(paired[BLOCK].drop_duplicates()),
            "matched_dates": int(paired.date.nunique()),
        }
    write_csv(pd.concat(summaries, ignore_index=True), output / "support_by_date.csv")
    write_csv(support.sort_values(KEY), output / "cell_date_support.csv")
    file_names = ["water_components.csv", "cell_covariates.csv", "support_by_date.csv", "cell_date_support.csv",
                  *[f"pairs_{name}.csv" for name in RULES]]
    record = {
        "validation_design_sha256": checksum(VALIDATION / "design_record.json"),
        "sla_source_sha256": MAP_SHA256,
        "sla_archive_sha256": checksum(source),
        "water_thresholds_km2": [PRIMARY_AREA_KM2, ALTERNATE_AREA_KM2],
        "target_tree_fraction": TARGET, "target_rangeland_fraction": -TARGET,
        "component_gap_tolerance": GAP_TOLERANCE, "land_gap_tolerance": LAND_TOLERANCE,
        "dates": sorted(eligible.date.unique()), "eligible_cell_dates": len(eligible),
        "eligible_cells": int(eligible.cell_id.nunique()), "match_block_size_m": 4000,
        "matching": matching,
        "design_sha256": {name: checksum(output / name) for name in sorted(file_names)},
    }
    (output / "design_record.json").write_text(json.dumps(record, indent=2) + "\n")
    print(json.dumps(matching, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--map", type=Path, default=MAP)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    prepare(args.map, args.output)
