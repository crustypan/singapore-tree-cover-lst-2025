"""Freeze Stage 1 eligible records and land-cover matches before temperature estimation."""
import argparse
from dataclasses import asdict, dataclass
import gzip
import hashlib
import json
from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd

from analyze_landsat import area_rule, block_coordinates, write_csv
from audit_coverage import BOUNDARY, ESRI, ESRI_CLASSES
from fetch_inputs import RAW, WORLD_FILE, checksum
from make_landsat_panel import DATES, VALID_ESRI, point_indices
from screen_landsat import land_layers


PANELS = Path("data/snapshot/landsat_panels")
BASELINE = Path("results/landsat_analysis/fitted_cells.csv")
MANIFEST = Path("data/selected_landsat_scenes.csv")
DESIGN = Path("results/local_comparison/design")
KEY = ["date", "cell_id"]
COORDS = ["easting_m", "northing_m"]
BLOCK = ["block_e", "block_n"]
COVARIATES = ["tree_fraction", "built_fraction", "land_share"]
AVAILABILITY = [
    *KEY, *COORDS, *COVARIATES, "classified_share", "valid_lst_share", "land_centers",
    "classified_centers",
]
PAIR_COLUMNS = [
    "date", *BLOCK, "lower_tree_cell", "higher_tree_cell", "distance_m",
    "tree_gap_pp", "built_gap_pp", "land_gap_pp", "cost",
]


@dataclass(frozen=True)
class MatchRule:
    block_size: int = 4000
    block_shift: int = 0
    distance: int = 2000
    built_gap: float = 0.05
    land_gap: float = 0.10
    min_tree_gap: float = 0.08
    max_tree_gap: float = 0.12


RULES = {
    "primary": MatchRule(),
    "shifted": MatchRule(block_shift=2000),
    "tighter": MatchRule(distance=1500, built_gap=0.03),
}


def read_design_inputs(directory):
    path = directory / "world_500.csv.gz"
    recorded = json.loads((directory / "provenance.json").read_text())
    with gzip.open(path, "rb") as f:
        content_hash = hashlib.sha256(f.read()).hexdigest()
    if content_hash != recorded["panel_sha256"]["world_500"]:
        raise ValueError("WorldCover panel checksum differs from the baseline snapshot")
    data = pd.read_csv(path, usecols=AVAILABILITY, float_precision="round_trip")
    observed = data[area_rule(data) & (data.valid_lst_share >= 0.6)].copy()
    keys = pd.MultiIndex.from_frame(observed[KEY]).sort_values()
    reference = pd.MultiIndex.from_frame(pd.read_csv(BASELINE, usecols=KEY)).sort_values()
    if not keys.equals(reference) or observed.duplicated(KEY).any():
        raise ValueError("QA-derived eligibility differs from baseline cell-date IDs")
    if tuple(sorted(observed.date.unique())) != DATES or observed[COVARIATES].isna().any().any():
        raise ValueError("Sample dates or land-cover covariates do not match the protocol")
    return data, observed, content_hash


def assign_blocks(data, rule=MatchRule()):
    blocks = block_coordinates(data, rule.block_size, rule.block_shift)
    return data.assign(block_e=blocks[:, 0], block_n=blocks[:, 1])


def minimum_cost_pairs(edges):
    """Maximize pair count, minimize quantized cost, and make exact ties unique."""
    edges = sorted(edges)
    graph = nx.Graph()
    # The sum of all tie-break weights fits in the low bits, below one unit of cost.
    scale = 1 << len(edges)
    for rank, (i, j, cost) in enumerate(edges):
        cost_units = round(float(cost) * 1_000_000)
        tie_break = 1 << (len(edges) - rank - 1)
        graph.add_edge(i, j, weight=-cost_units * scale + tie_break)
    pairs = nx.max_weight_matching(graph, maxcardinality=True)
    return sorted(tuple(sorted(pair)) for pair in pairs)


def match_group(group, rule):
    group = group.sort_values("cell_id", kind="stable").reset_index(drop=True)
    ids = group.cell_id.to_numpy()
    east, north = group[COORDS].to_numpy().T
    tree, built, land = group[COVARIATES].to_numpy().T
    i, j = np.triu_indices(len(group), k=1)
    distance = np.hypot(east[j] - east[i], north[j] - north[i])
    dt = np.abs(tree[j] - tree[i])
    db = np.abs(built[j] - built[i])
    dl = np.abs(land[j] - land[i])
    # The input fractions are published to eight decimals. Retain exact cutoffs
    # when the subtraction introduces a floating-point rounding error.
    eps = 1e-12
    accept = (
        (distance <= rule.distance) & (db <= rule.built_gap + eps) & (dl <= rule.land_gap + eps)
        & (dt >= rule.min_tree_gap - eps) & (dt <= rule.max_tree_gap + eps)
    )
    cost = (db / rule.built_gap) ** 2 + (dl / rule.land_gap) ** 2 + (distance / rule.distance) ** 2
    cost += ((dt - 0.10) / 0.02) ** 2
    positions = {(int(i[k]), int(j[k])): k for k in np.flatnonzero(accept)}
    edges = [(a, b, cost[k]) for (a, b), k in positions.items()]
    records = []
    for a, b in minimum_cost_pairs(edges):
        k = positions[a, b]
        lo, hi = (a, b) if tree[a] < tree[b] else (b, a)
        records.append((
            ids[lo], ids[hi], distance[k], 100 * (tree[hi] - tree[lo]),
            100 * (built[hi] - built[lo]), 100 * (land[hi] - land[lo]), cost[k],
        ))
    return records, set(ids[np.r_[i[accept], j[accept]]]), int(accept.sum())


def build_matches(observed, rule):
    pairs, counts, available = [], [], []
    assigned = assign_blocks(observed, rule)
    for date, day in assigned.groupby("date", sort=True):
        n_pairs, n_edges, n_blocks = 0, 0, 0
        day_cells = set()
        for block, group in day.groupby(BLOCK, sort=True):
            records, candidates, edges = match_group(group, rule)
            pairs.extend((date, *block, *record) for record in records)
            n_pairs += len(records)
            n_edges += edges
            n_blocks += bool(candidates)
            day_cells.update(candidates)
        counts.append({
            "date": date, "usable_cells": len(day), "candidate_cells": len(day_cells),
            "candidate_blocks": n_blocks, "candidate_links": n_edges, "pairs": n_pairs,
            "matched_cells": 2 * n_pairs,
        })
        available.extend({"date": date, "cell_id": cell} for cell in sorted(day_cells))
    frame = pd.DataFrame(pairs, columns=PAIR_COLUMNS)
    low, high = (frame[["date", column]].rename(columns={column: "cell_id"})
                 for column in ("lower_tree_cell", "higher_tree_cell"))
    matched = pd.concat([low, high], ignore_index=True)
    if matched.duplicated(KEY).any():
        raise ValueError("The design reused a cell within a date")
    return frame, pd.DataFrame(counts), pd.DataFrame(available, columns=KEY), matched


def cover_composition(panel):
    affine, world, _, _, esri = land_layers(trim_to_30m=False)
    row, col = np.nonzero(world)
    x = affine.c + (col + 0.5) * affine.a
    y = affine.f + (row + 0.5) * affine.e
    idx, centers, land = point_indices(x, y, 500, 0)
    classes = esri[row, col]
    counts = {code: np.bincount(idx[classes == code], minlength=len(land)) for code in ESRI_CLASSES}
    valid = np.sum([counts[code] for code in VALID_ESRI], axis=0)
    if not np.array_equal(np.sum(list(counts.values()), axis=0), land):
        raise ValueError("Esri codes do not account for every WorldCover land center")
    result = pd.DataFrame({
        "cell_id": [f"{e}_{n}" for e, n in centers], "land_centers": land,
        "classified_centers": valid, **{f"class_{code:02d}_centers": col for code, col in counts.items()},
    })
    check = panel.drop_duplicates("cell_id").set_index("cell_id")
    got = result.set_index("cell_id").reindex(check.index)
    for column in ("land_centers", "classified_centers"):
        np.testing.assert_array_equal(got[column], check[column])
    for code, name in [(2, "tree_fraction"), (7, "built_fraction")]:
        fraction = got[f"class_{code:02d}_centers"].div(got.classified_centers.replace(0, np.nan))
        np.testing.assert_allclose(fraction, check[name], atol=5.1e-9, rtol=0, equal_nan=True)
    return result


def prepare(panels, output):
    data, observed, content_hash = read_design_inputs(panels)
    source_hashes = {p.name: checksum(p) for p in (MANIFEST, BOUNDARY, ESRI, RAW / WORLD_FILE)}
    reference = json.loads((panels / "provenance.json").read_text())["input_sha256"]
    if source_hashes != reference:
        raise ValueError("Spatial inputs or scene IDs differ from the baseline snapshot")
    output.mkdir(parents=True, exist_ok=True)
    outputs = {}
    summary = []
    for name, rule in RULES.items():
        pairs, dates, available, matched = build_matches(observed, rule)
        write_csv(pairs, output / f"pairs_{name}.csv")
        summary.append(dates.assign(design=name))
        outputs[name] = {
            "rule": asdict(rule), "pairs": len(pairs), "matched_unique_cells": int(matched.cell_id.nunique()),
            "candidate_unique_cells": int(available.cell_id.nunique()),
            "matching_blocks": len(pairs[BLOCK].drop_duplicates()),
        }
        if name == "primary":
            cells = data.drop_duplicates("cell_id").set_index("cell_id")
            for label, rows in [("usable", observed), ("candidate", available), ("paired", matched)]:
                cells[f"{label}_dates"] = rows.groupby("cell_id").size().reindex(cells.index, fill_value=0)
            columns = [*COORDS, *COVARIATES, "classified_share", "land_centers",
                       "usable_dates", "candidate_dates", "paired_dates"]
            write_csv(cells[columns].reset_index(), output / "cell_design.csv")
    write_csv(pd.concat(summary, ignore_index=True), output / "design_by_date.csv")
    composition = cover_composition(data)
    write_csv(composition, output / "land_composition.csv")
    files = [output / name for name in (
        "cell_design.csv", "design_by_date.csv", "land_composition.csv",
        *[f"pairs_{name}.csv" for name in RULES],
    )]
    info = {
        "dates": list(DATES), "eligible_cell_dates": len(observed),
        "eligible_cells": int(observed.cell_id.nunique()), "matching": outputs,
        "world_500_decompressed_sha256": content_hash,
        "source_sha256": source_hashes,
        "design_sha256": {p.name: checksum(p) for p in sorted(files)},
    }
    (output / "design_record.json").write_text(json.dumps(info, indent=2) + "\n")
    print(json.dumps(outputs, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--panels", type=Path, default=PANELS)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    prepare(args.panels, args.output)
