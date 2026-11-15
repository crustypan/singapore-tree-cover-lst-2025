"""Freeze the cells and original support rules for the land-cover audit."""
import argparse
from itertools import combinations
import json
from pathlib import Path

import pandas as pd

from analyze_landsat import write_csv
from fetch_inputs import checksum
from prepare_built_support import DESIGN_DIR as BUILT_DESIGN, STRATUM, selected_samples
from prepare_local_comparison import COORDS, KEY, assign_blocks
from prepare_overpass_validation import DESIGN as VALIDATION_DESIGN, YEARS, verify_hashes


AUDIT_ROOT = Path("results/cover_support_audit")
AUDIT_DESIGN = AUDIT_ROOT / "design"
SHARES = ("tree", "built", "water")
SPAN = {"tree": "tree_fraction_span_pp", "built": "built_fraction_span_pp",
        "water": "esri_water_share_span_pp"}
EPS = 1e-6


def support_inputs():
    validation = VALIDATION_DESIGN / "design_record.json"
    built = BUILT_DESIGN / "design_record.json"
    config = json.loads(validation.read_text())
    verify_hashes(VALIDATION_DESIGN, config["design_sha256"])
    target = json.loads(built.read_text())["chosen"]
    g = pd.read_csv(VALIDATION_DESIGN / "annual_cover_comparison.csv", float_precision="round_trip")
    events = pd.read_csv(VALIDATION_DESIGN / "eligible_cell_dates.csv").query("common").copy()
    events["year"] = events.date.str[:4].astype(int)
    g = g[g.common].reset_index(drop=True)
    if set(events.cell_id) != set(g.cell_id) or events.duplicated(KEY).any():
        raise ValueError("Common geographies and the locked overpass rows disagree")
    return g, events, target, config["cover_stability_max_pp"]


def nominal_high(geometry, year, contrast):
    b = geometry[f"built_fraction_{year}"] * 100
    t = geometry[f"tree_fraction_{year}"] * 100
    return (b.sub(contrast["lower_built_pct"]).abs().le(5 + EPS)
            & t.sub(contrast["tree_end_pct"]).abs().le(2.5 + EPS)
            & geometry[f"usable_dates_{year}"].gt(0))


def esri_rows(geometry, events):
    cols = ["cell_id", *COORDS, "land_share"]
    frames = []
    for year in YEARS:
        fields = {f"{c}_{year}": c for c in ("tree_fraction", "built_fraction", "esri_water_share")}
        covariates = geometry[[*cols, *fields]].rename(columns=fields)
        frames.append(events[events.year == year].merge(covariates, on="cell_id", validate="many_to_one"))
    return pd.concat(frames, ignore_index=True)


def count_blocks(rows, target, label):
    summaries, location_rows = [], []
    for part, frame in [(str(y), rows[rows.year == y]) for y in YEARS] + [("all", rows)]:
        _, arms, shared, strict = selected_samples(assign_blocks(frame), target["tree_start_pct"],
                                                   target["lower_built_pct"], target["higher_built_pct"])
        record = {"analysis": label, "year_group": part, "cells": int(frame.cell_id.nunique()),
                  "cell_dates": len(frame), "shared_blocks": len({k[1:] for k in shared}),
                  "shared_block_dates": len(shared), "shared_strict_blocks": len({k[1:] for k in strict}),
                  "shared_strict_block_dates": len(strict)}
        for band, arm in arms.items():
            stats = arm[0]
            for name in ("endpoint_blocks", "endpoint_block_dates", "strict_blocks", "strict_block_dates",
                         "low_tight_endpoint_cells", "high_tight_endpoint_cells"):
                record[f"{band}_{name}"] = stats[name]
        summaries.append(record)
        lo, hi = arms["lower"][1], arms["higher"][1]
        for key in sorted(lo | hi):
            location_rows.append({"analysis": label, "year_group": part, **dict(zip(STRATUM, key)),
                                  "lower": key in lo, "higher": key in hi, "shared_strict": key in strict})
    return pd.DataFrame(summaries), pd.DataFrame(location_rows)


def prepare(out):
    g, events, target, threshold = support_inputs()
    cases = g.assign(**{f"target_{year}": nominal_high(g, year, target) for year in YEARS})
    cases = cases[cases[[f"target_{y}" for y in YEARS]].any(axis=1)].copy()
    cases = cases.sort_values(["northing_m", "easting_m"]).reset_index(drop=True)
    cases.insert(0, "case", [f"C{i:02d}" for i in range(1, len(cases) + 1)])
    names = ["case", "cell_id", *COORDS, "land_share", *SPAN.values()]
    for year in YEARS:
        names.extend([f"target_{year}", f"tree_fraction_{year}", f"built_fraction_{year}", f"esri_water_share_{year}"])
    names += ["common_stable", "boundary_mixed"]
    cases = cases[names]

    rows = esri_rows(g, events)
    summary, geography = [], []
    for n in range(len(SHARES) + 1):
        for combo in combinations(SHARES, n):
            mask = g[[SPAN[name] for name in combo]].le(threshold + EPS).all(axis=1)
            ids = set(g.loc[mask, "cell_id"])
            if combo == SHARES and ids != set(g.loc[g.common_stable, "cell_id"]):
                raise ValueError("The three original span checks differ from Stage 3")
            counts, sites = count_blocks(rows[rows.cell_id.isin(ids)], target, "+".join(combo) or "none")
            summary.append(counts)
            geography.append(sites)

    out.mkdir(parents=True, exist_ok=True)
    write_csv(cases, out / "critical_cells.csv")
    write_csv(pd.concat(summary, ignore_index=True), out / "esri_filter_counts.csv")
    write_csv(pd.concat(geography, ignore_index=True), out / "esri_filter_block_dates.csv")
    record = {
        "years": list(YEARS), "common_cells": len(g), "critical_cells": len(cases),
        "critical_cell_years": {str(y): int(cases[f"target_{y}"].sum()) for y in YEARS},
        "target": {k: target[k] for k in ("tree_start_pct", "tree_end_pct", "lower_built_pct", "higher_built_pct")},
        "stability_max_pp": threshold,
        "stage3_design_sha256": checksum(VALIDATION_DESIGN / "design_record.json"),
        "stage2_design_sha256": checksum(BUILT_DESIGN / "design_record.json"),
        "files_sha256": {p.name: checksum(p) for p in sorted(out.glob("*.csv"))},
    }
    (out / "design_record.json").write_text(json.dumps(record, indent=2) + "\n")
    print(json.dumps(record, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    opts = parser.parse_args()
    prepare(opts.output)
