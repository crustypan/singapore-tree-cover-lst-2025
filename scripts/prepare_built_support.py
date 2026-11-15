"""Audit land-cover support and freeze the Stage 2 samples without using LST values."""
import argparse
from itertools import combinations
import json
from pathlib import Path

import numpy as np
import pandas as pd

from analyze_landsat import write_csv
from fetch_inputs import checksum
from make_landsat_panel import DATES
from prepare_local_comparison import (
    BLOCK, COORDS, COVARIATES, KEY, PANELS, PAIR_COLUMNS,
    assign_blocks, minimum_cost_pairs, read_design_inputs,
)


DESIGN_DIR = Path("results/built_cover_moderation/design")
STRATUM = ["date", *BLOCK]
STARTS = tuple(np.arange(0, 30.1, 2.5).tolist())
BUILT_LEVELS = tuple(range(25, 91, 5))
TREE_RADIUS = 5.0
BUILT_RADIUS = 10.0
EPS = 1e-10


def group_keys(data):
    return set(data[STRATUM].itertuples(index=False, name=None))


def rows_in_groups(data, groups):
    keys = pd.DataFrame(sorted(groups), columns=STRATUM)
    return data.merge(keys, on=STRATUM, how="inner", validate="many_to_one")


def band_support(data, tree_start, built, tree_radius=TREE_RADIUS, built_radius=BUILT_RADIUS):
    tree_end = tree_start + 10
    tree = 100 * data.tree_fraction
    near_built = (100 * data.built_fraction - built).abs() < built_radius - EPS
    low = near_built & ((tree - tree_start).abs() <= tree_radius + EPS)
    high = near_built & ((tree - tree_end).abs() <= tree_radius + EPS)
    fields = [*STRATUM, "cell_id", *COORDS, *COVARIATES]
    links = data.loc[low, fields].merge(data.loc[high, fields], on=STRATUM, suffixes=("_lo", "_hi"))
    links["tree_gap_pp"] = 100 * (links.tree_fraction_hi - links.tree_fraction_lo)
    links = links[links.tree_gap_pp.between(8 - EPS, 12 + EPS)].copy()
    links["built_gap_pp"] = 100 * (links.built_fraction_hi - links.built_fraction_lo)
    links["land_gap_pp"] = 100 * (links.land_share_hi - links.land_share_lo)
    links["distance_m"] = np.hypot(links.easting_m_hi - links.easting_m_lo,
                                   links.northing_m_hi - links.northing_m_lo)
    is_strict = ((links.distance_m <= 2000) & (links.built_gap_pp.abs() <= 5 + EPS)
                 & (links.land_gap_pp.abs() <= 10 + EPS))
    strict = links[is_strict].copy()
    tighter_built = (100 * data.built_fraction - built).abs() <= 5 + EPS
    tighter_counts = [
        data.loc[tighter_built & ((tree - target).abs() <= 2.5 + EPS), "cell_id"].nunique()
        for target in (tree_start, tree_end)
    ]
    groups = group_keys(links)
    strict_groups = group_keys(strict)
    stats = {
        "tree_start_pct": tree_start, "tree_end_pct": tree_end, "built_reference_pct": built,
        "feasible_reference": bool(tree_end + built <= 97.5),
        "low_endpoint_cells": int(data.loc[low, "cell_id"].nunique()),
        "high_endpoint_cells": int(data.loc[high, "cell_id"].nunique()),
        "low_tight_endpoint_cells": int(tighter_counts[0]), "high_tight_endpoint_cells": int(tighter_counts[1]),
        "endpoint_blocks": len({k[1:] for k in groups}), "endpoint_block_dates": len(groups),
        "endpoint_dates": len({k[0] for k in groups}), "strict_blocks": len({k[1:] for k in strict_groups}),
        "strict_dates": len({k[0] for k in strict_groups}), "strict_block_dates": len(strict_groups),
        "strict_candidate_links": len(strict), "strict_low_cells": int(strict.cell_id_lo.nunique()),
        "strict_high_cells": int(strict.cell_id_hi.nunique()),
    }
    sample = rows_in_groups(data.loc[near_built & tree.between(tree_start - tree_radius - EPS,
                            tree_end + tree_radius + EPS)], groups)
    return stats, groups, strict_groups, sample, strict


def shortlist(audits):
    rows = []
    for t in STARTS:
        for lo, hi in combinations(BUILT_LEVELS, 2):
            if hi - lo < 2 * BUILT_RADIUS:
                continue
            low_stats, low, low_strict = audits[t, lo]
            high_stats, high, high_strict = audits[t, hi]
            if not high_stats["feasible_reference"]:
                continue
            both = low & high
            strict = low_strict & high_strict
            min_tight = min(s[k] for s in (low_stats, high_stats)
                            for k in ("low_tight_endpoint_cells", "high_tight_endpoint_cells"))
            rows.append({
                "tree_start_pct": t, "tree_end_pct": t + 10, "lower_built_pct": lo, "higher_built_pct": hi,
                "built_gap_pp": hi - lo, "minimum_tight_endpoint_cells": min_tight,
                "shared_blocks": len({k[1:] for k in both}), "shared_dates": len({k[0] for k in both}),
                "shared_block_dates": len(both), "shared_strict_blocks": len({k[1:] for k in strict}),
                "shared_strict_block_dates": len(strict), "shared_strict_dates": len({k[0] for k in strict}),
            })
    return pd.DataFrame(rows).sort_values(
        ["shared_blocks", "shared_block_dates", "minimum_tight_endpoint_cells", "built_gap_pp",
         "tree_start_pct", "lower_built_pct", "higher_built_pct"],
        ascending=[False, False, False, False, True, True, True],
    ).reset_index(drop=True)


def disjoint_pairs(links, label):
    rows = []
    for group, edges in links.groupby(STRATUM, sort=True):
        ids = sorted(set(edges.cell_id_lo) | set(edges.cell_id_hi))
        index = {cell: i for i, cell in enumerate(ids)}
        by_edge = {}
        graph = []
        for e in edges.itertuples(index=False):
            i, j = sorted((index[e.cell_id_lo], index[e.cell_id_hi]))
            cost = (e.built_gap_pp / 5) ** 2 + (e.land_gap_pp / 10) ** 2
            cost += (e.distance_m / 2000) ** 2 + ((e.tree_gap_pp - 10) / 2) ** 2
            by_edge[i, j] = (e, cost)
            graph.append((i, j, cost))
        for key in minimum_cost_pairs(graph):
            e, cost = by_edge[key]
            rows.append((*group, e.cell_id_lo, e.cell_id_hi, e.distance_m, e.tree_gap_pp,
                         e.built_gap_pp, e.land_gap_pp, cost, label))
    return pd.DataFrame(rows, columns=[*PAIR_COLUMNS, "band"])


def selected_samples(data, tree_start, lower, higher, tree_radius=TREE_RADIUS, built_radius=BUILT_RADIUS):
    arms = {
        name: band_support(data, tree_start, ref, tree_radius, built_radius)
        for name, ref in [("lower", lower), ("higher", higher)]
    }
    common = arms["lower"][1] & arms["higher"][1]
    common_strict = arms["lower"][2] & arms["higher"][2]
    samples = pd.concat([entry[3].assign(band=label) for label, entry in arms.items()], ignore_index=True)
    if samples.duplicated(KEY).any():
        raise ValueError("Built bands share one or more grid cell-dates")
    samples["shared"] = samples[STRATUM].apply(tuple, axis=1).isin(common)
    return samples, arms, common, common_strict


def prepare(panels, output):
    _, eligible, content_hash = read_design_inputs(panels)
    data = assign_blocks(eligible)
    support = {}
    for t in STARTS:
        for b in BUILT_LEVELS:
            support[t, b] = band_support(data, t, b)[:3]
    by_band = pd.DataFrame([entry[0] for entry in support.values()])
    candidates = shortlist(support)
    valid = candidates[(candidates.minimum_tight_endpoint_cells > 0) & (candidates.shared_blocks > 0)]
    if valid.empty:
        raise ValueError("No eligible pair of reference bands has shared geographic support")
    best = valid.iloc[0]
    t, lo, hi = (float(best[column]) for column in ("tree_start_pct", "lower_built_pct", "higher_built_pct"))
    samples, arms, common, common_strict = selected_samples(data, t, lo, hi)
    _, _, tight_common, tight_strict = selected_samples(data, t, lo, hi, 2.5, 7.5)
    pairs = pd.concat([disjoint_pairs(entry[4], name) for name, entry in arms.items()], ignore_index=True)
    pairs["shared_strict"] = pairs[STRATUM].apply(tuple, axis=1).isin(common_strict)
    unique = data.drop_duplicates("cell_id")
    edges = np.arange(0., 100.1, 5.)
    histogram, _, _ = np.histogram2d(100 * unique.built_fraction, 100 * unique.tree_fraction, bins=[edges, edges])
    density = pd.DataFrame([
        {"built_from_pct": x, "tree_from_pct": y, "cells": int(histogram[i, j])}
        for i, x in enumerate(edges[:-1]) for j, y in enumerate(edges[:-1])
    ])
    geo = set.union(*(entry[1] for entry in arms.values()))
    blocks = pd.DataFrame([
        dict(zip(STRATUM, key), lower=key in arms["lower"][1], higher=key in arms["higher"][1],
             shared=key in common, shared_strict=key in common_strict)
        for key in sorted(geo)
    ])
    output.mkdir(parents=True, exist_ok=True)
    outputs = {
        "joint_density.csv": density, "band_support.csv": by_band, "shared_support.csv": candidates,
        "selected_block_dates.csv": blocks, "fit_sample.csv": samples[[*KEY, *BLOCK, "band", "shared"]],
        "strict_pairs.csv": pairs,
    }
    for name, frame in outputs.items():
        write_csv(frame, output / name)
    record = {
        "dates": list(DATES), "all_usable_cells": int(len(unique)), "all_usable_cell_dates": len(data),
        "world_500_decompressed_sha256": content_hash, "tree_radius_pp": TREE_RADIUS,
        "built_radius_pp_exclusive": BUILT_RADIUS,
        "chosen": {k: (float(v) if k.startswith("tree") else int(v)) for k, v in best.items()},
        "band_support": {name: entry[0] for name, entry in arms.items()},
        "narrower_window_shared_blocks": len({key[1:] for key in tight_common}),
        "narrower_window_shared_strict_blocks": len({key[1:] for key in tight_strict}),
        "sample_by_band": {name: {
            "cell_dates": len(group), "distinct_cells": int(group.cell_id.nunique()),
            "shared_cell_dates": int(group.shared.sum()),
            "shared_distinct_cells": int(group.loc[group.shared, "cell_id"].nunique()),
        } for name, group in samples.groupby("band")},
        "design_sha256": {p: checksum(output / p) for p in sorted(outputs)},
    }
    (output / "design_record.json").write_text(json.dumps(record, indent=2) + "\n")
    print(json.dumps({k: v for k, v in record.items() if k != "design_sha256"}, indent=2))
    print("Leading eligible contrasts")
    print(valid.head(6).to_string(index=False))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--panels", type=Path, default=PANELS)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    prepare(args.panels, args.output)
