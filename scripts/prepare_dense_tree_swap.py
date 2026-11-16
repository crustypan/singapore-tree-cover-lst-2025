"""Audit rangeland and local pair support before using Landsat temperatures."""
import argparse
from collections import defaultdict
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import numpy as np
import pandas as pd

from analyze_landsat import write_csv
from dense_swap_model import MIN_RESAMPLE_RETENTION, assess
from fetch_inputs import checksum
from prepare_grass_tree_substitution import COVER, DESIGN as STAGE5, EPS, PAIR_COLUMNS, VALIDATION
from prepare_local_comparison import BLOCK, COORDS, KEY, assign_blocks, minimum_cost_pairs
from prepare_overpass_validation import verify_hashes


OUTPUT = Path("results/dense_tree_swap")
DESIGN = OUTPUT / "design"
BANDS = {
    "lt50": (0., 0.50), "50to70": (0.50, 0.70), "70to80": (0.70, 0.80),
    "80to90": (0.80, 0.90), "ge90": (0.90, None), "ge80": (0.80, None),
}
STEPS = ("eligible", "donor", "opposite_gaps", "cover_balance", "within_2km", "same_block", "matched")


@dataclass(frozen=True)
class Rule:
    target: float
    margin: float
    built_other: float
    combined: float
    land: float = 0.10
    distance: int = 2000


RULES = {
    "10pp": Rule(0.10, 0.02, 0.05, 0.02),
    "5pp": Rule(0.05, 0.01, 0.025, 0.01),
    "5pp_stage5_balance": Rule(0.05, 0.01, 0.05, 0.02),
}
ESTIMATE_CASES = {
    "dense": ("5pp", "ge80"), "very_dense": ("5pp", "ge90"), "80_to_90": ("5pp", "80to90"),
    "wider_balance": ("5pp_stage5_balance", "ge80"),
}


def read_inputs(stage5):
    recorded = json.loads((stage5 / "design_record.json").read_text())
    verify_hashes(stage5, recorded["design_sha256"])
    verify_hashes(VALIDATION, {"design_record.json": recorded["validation_design_sha256"]})
    validation = json.loads((VALIDATION / "design_record.json").read_text())
    verify_hashes(VALIDATION, validation["design_sha256"])
    available = pd.read_csv(stage5 / "cell_date_support.csv", usecols=KEY).sort_values(KEY).reset_index(drop=True)
    original = pd.read_csv(VALIDATION / "eligible_cell_dates.csv", usecols=KEY).sort_values(KEY).reset_index(drop=True)
    pd.testing.assert_frame_equal(available, original)
    cells = pd.read_csv(stage5 / "cell_covariates.csv", float_precision="round_trip")
    denom = cells.classified_centers
    components = cells[["class_02_centers", "class_11_centers", "class_07_centers"]].to_numpy()
    exact = np.column_stack([components, denom - components.sum(axis=1)]) / denom.to_numpy()[:, None]
    np.testing.assert_allclose(exact, cells[COVER[:4]], atol=5.1e-10, rtol=0)
    cells[COVER[:4]] = exact
    eligible = available.assign(year=available.date.str[:4].astype(int)).merge(
        cells, on=["year", "cell_id"], how="left", validate="many_to_one"
    )
    if eligible.duplicated(KEY).any() or not np.isfinite(eligible[[*COORDS, *COVER]].to_numpy()).all():
        raise ValueError("Eligible observations and archived land-cover counts must agree")
    return assign_blocks(eligible)


def donor_funnel(day, rule):
    """Count oriented candidate links at each cutoff, then select disjoint pairs."""
    day = day.sort_values("cell_id").reset_index(drop=True)
    shares = day[COVER].to_numpy()
    coords = day[COORDS].to_numpy()
    blocks = day[BLOCK].to_numpy()
    ids = day.cell_id.to_numpy()
    donors = np.flatnonzero(shares[:, 1] >= rule.target - EPS)
    step = np.zeros(len(day), dtype=np.int8)
    step[donors] = 1

    n = len(day)
    lo, hi = np.repeat(donors, n), np.tile(np.arange(n), len(donors))
    tree = shares[hi, 0] - shares[lo, 0]
    accept = np.abs(tree - rule.target) <= rule.margin + EPS
    lo, hi = lo[accept], hi[accept]
    grass_loss = shares[lo, 1] - shares[hi, 1]
    accept = np.abs(grass_loss - rule.target) <= rule.margin + EPS
    lo, hi = lo[accept], hi[accept]
    step[lo] = 2

    gaps = shares[hi] - shares[lo]
    dt, dr, db, do, dl = gaps.T
    balanced = ((np.abs(dt + dr) <= rule.combined + EPS)
                & (np.maximum(np.abs(db), np.abs(do)) <= rule.built_other + EPS)
                & (np.abs(dl) <= rule.land + EPS))
    step[lo[balanced]] = 3
    separation = np.linalg.norm(coords[hi] - coords[lo], axis=1)
    close = balanced & (separation <= rule.distance)
    step[lo[close]] = 4
    local = close & (blocks[hi] == blocks[lo]).all(axis=1)
    step[lo[local]] = 5

    cost = ((dt - rule.target) / rule.margin) ** 2 + ((dr + rule.target) / rule.margin) ** 2
    cost += (db / rule.built_other) ** 2 + (do / rule.built_other) ** 2
    cost += (dl / rule.land) ** 2 + (separation / rule.distance) ** 2
    edges = defaultdict(list)
    positions = {}
    for k in np.flatnonzero(local):
        i, j = sorted((int(lo[k]), int(hi[k])))
        positions[i, j] = k
        edges[tuple(blocks[i])].append((i, j, cost[k]))
    pairs = []
    for block, connections in sorted(edges.items()):
        for indices in minimum_cost_pairs(connections):
            k = positions[indices]
            i, j = lo[k], hi[k]
            step[i] = len(STEPS) - 1
            pairs.append((*block, ids[i], ids[j], separation[k], *(100 * gaps[k]), cost[k]))
    return day.assign(step=step), pairs, int(local.sum())


def stratum(data, bounds):
    lower, upper = bounds
    group = data[data.built_fraction.ge(lower)]
    return group if upper is None else group[group.built_fraction.lt(upper)]


def support_table(eligible, rule):
    counts, members, links = [], [], []
    for date, day in eligible.groupby("date", sort=True):
        status, chosen, local_links = donor_funnel(day, rule)
        links.extend((date, *row) for row in chosen)
        members.append(status[["date", "cell_id", "year", *BLOCK, "step"]])
        counts.append({
            "date": date, **{name: int(status.step.ge(i).sum()) for i, name in enumerate(STEPS)},
            "local_links": local_links,
        })
    pairs = pd.DataFrame(links, columns=PAIR_COLUMNS)
    rows = pd.concat(members, ignore_index=True)
    summary = []
    for i, step in enumerate(STEPS):
        has = rows[rows.step.ge(i)]
        summary.append({
            "step": step, "cell_dates": len(has), "cell_years": len(has[["cell_id", "year"]].drop_duplicates()),
            "cells": int(has.cell_id.nunique()), "dates": int(has.date.nunique()),
            "blocks": len(has[BLOCK].drop_duplicates()),
        })
    endpoints = pd.concat([
        pairs[["date", side]].rename(columns={side: "cell_id"})
        for side in ("lower_tree_cell", "higher_tree_cell")
    ], ignore_index=True)
    if endpoints.duplicated(KEY).any() or len(pairs) != rows.step.eq(len(STEPS) - 1).sum():
        raise ValueError("Matching must be disjoint and agree with the donor registry")
    return pairs, rows, pd.DataFrame(counts), pd.DataFrame(summary)


def composition_table(cells):
    records = []
    for name, bounds in BANDS.items():
        cover = stratum(cells, bounds).rangeland_fraction
        records.append({
            "built_band": name, "cell_years": len(cover), "zero": int(cover.eq(0).sum()),
            "below_5": int(cover.lt(0.05).sum()), "at_least_5": int(cover.ge(0.05).sum()),
            "at_least_10": int(cover.ge(0.10).sum()), "rangeland_p50_pp": 100 * cover.median(),
            "rangeland_p90_pp": 100 * cover.quantile(0.9), "rangeland_p95_pp": 100 * cover.quantile(0.95),
        })
    return pd.DataFrame(records)


def prepare(stage5, output):
    eligible = read_inputs(stage5)
    output.mkdir(parents=True, exist_ok=True)
    write_csv(composition_table(eligible.drop_duplicates(["cell_id", "year"])), output / "composition_by_band.csv")
    counts, summaries, pair_tables, dense_members = [], [], [], []
    for study, rule in RULES.items():
        for band, bounds in BANDS.items():
            if study == "5pp_stage5_balance" and band not in ("ge80", "ge90", "80to90"):
                continue
            pairs, records, daily, totals = support_table(stratum(eligible, bounds), rule)
            if not pairs.empty:
                pair_tables.append(pairs.assign(comparison=study, built_band=band))
            counts.append(daily.assign(comparison=study, built_band=band))
            summaries.append(totals.assign(comparison=study, built_band=band))
            if band == "ge80":
                dense_members.append(records.assign(comparison=study, built_band=band))
    pairs = pd.concat(pair_tables, ignore_index=True)
    write_csv(pairs, output / "pairs.csv")
    write_csv(pd.concat(counts, ignore_index=True), output / "support_by_date.csv")
    write_csv(pd.concat(summaries, ignore_index=True), output / "support_by_band.csv")
    # The full registry is limited to the pooled high-built sample used on the map.
    registry = pd.concat(dense_members, ignore_index=True)
    registry["last_supported_step"] = np.asarray(STEPS)[registry.step.to_numpy()]
    registry = registry.drop(columns="step")
    registry.to_csv(output / "dense_cells.csv.gz", index=False, lineterminator="\n",
                    compression={"method": "gzip", "mtime": 0})
    pairs = pd.read_csv(output / "pairs.csv", float_precision="round_trip")
    cases = {}
    for label, (study, band) in ESTIMATE_CASES.items():
        data = pairs[pairs.comparison.eq(study) & pairs.built_band.eq(band)]
        cases[label] = {"comparison": study, "built_band": band, **assess(data)}
    record = {
        "stage5_design_sha256": checksum(stage5 / "design_record.json"),
        "comparison_rules": {key: asdict(rule) for key, rule in RULES.items()},
        "bands": BANDS, "steps": STEPS, "adjusted_draw_fraction_minimum": MIN_RESAMPLE_RETENTION,
        "estimation": cases,
        "design_sha256": {name: checksum(output / name) for name in (
            "composition_by_band.csv", "pairs.csv", "support_by_date.csv", "support_by_band.csv", "dense_cells.csv.gz"
        )},
    }
    (output / "design_record.json").write_text(json.dumps(record, indent=2) + "\n")
    print(json.dumps(cases, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage5", type=Path, default=STAGE5)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    prepare(args.stage5, args.output)
