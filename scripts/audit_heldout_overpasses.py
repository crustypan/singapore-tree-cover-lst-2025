"""Audit 2025 held-out and earlier candidate overpasses before fitting temperatures."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timedelta, timezone
import json
from pathlib import Path
import re
from urllib.parse import urlsplit

import numpy as np
import pandas as pd
import rasterio
from rasterio.transform import Affine, array_bounds

from analyze_landsat import area_rule, write_csv
from audit_coverage import BOUNDARY, ESRI
from fetch_inputs import ESRI_URL, RAW, WORLD_FILE, atomic_target, checksum, fetch_esri_window, get_json
from make_landsat_panel import SOURCE_FILES, VALID_ESRI, counts_per_cell, fractions, point_indices, points_within_native
from prepare_built_support import DESIGN_DIR as BUILT_DESIGN, selected_samples
from prepare_local_comparison import COVARIATES, DESIGN, KEY, PANELS, RULES, assign_blocks, build_matches, read_design_inputs
from screen_landsat import FIELDS, SAS, STAC, catalog, compute_scene, land_layers, mask_values, native_window, spatial_grid


AUDIT = Path("results/overpass_validation/audit")
PRIMARY_CACHE = RAW / "validation_windows" / "primary"
ARCHIVED_SCREEN = Path("results/landsat_screening/scene_coverage.csv")
DISCOVERY_SCENES = Path("data/selected_landsat_scenes.csv")
SEASONS = ("Northeast monsoon", "First inter-monsoon", "Southwest monsoon", "Second inter-monsoon")
MANIFEST_COLUMNS = ("scene_id", "datetime_utc", "platform", "sun_elevation_deg", *SOURCE_FILES)
PRIMARY_BANDS = ("lwir11", "qa_pixel")
LOCAL_TIME = timezone(timedelta(hours=8))


def period(date):
    month = int(date[5:7])
    if month <= 3 or month == 12:
        return SEASONS[0]
    if month <= 5:
        return SEASONS[1]
    if month <= 9:
        return SEASONS[2]
    return SEASONS[3]


def land_cover(year):
    url = ESRI_URL.replace("2025", str(year))
    if year == 2023:
        url = url.replace("20231231", "20240101")
    path = ESRI if year == 2025 else RAW / f"esri_impact_observatory_{year}_48N_singapore.tif"
    fetch_esri_window(url, path)
    with rasterio.open(ESRI) as reference, rasterio.open(path) as source:
        if (source.transform, source.crs, source.shape) != (reference.transform, reference.crs, reference.shape):
            raise ValueError("Annual Esri cover grids differ from the discovery grid")
    return path, url


def read_candidates(year, frozen):
    if frozen is not None:
        candidates = pd.read_csv(frozen)
        if tuple(candidates) != MANIFEST_COLUMNS:
            raise ValueError("Candidate manifest columns have changed")
    else:
        rows = []
        for scene in catalog(year):
            props = scene["properties"]
            if props["landsat:collection_category"] != "T1":
                continue
            assets = [scene["assets"][key]["href"] for key in SOURCE_FILES]
            rows.append((scene["id"], props["datetime"], props["platform"],
                         props["view:sun_elevation"], *assets))
        candidates = pd.DataFrame(rows, columns=MANIFEST_COLUMNS)
    if (candidates.empty or not candidates.scene_id.is_unique
            or candidates.datetime_utc.str[:10].duplicated().any()
            or not candidates.datetime_utc.str.startswith(f"{year}-").all()
            or not candidates.scene_id.str.fullmatch(r"LC0[89]_L2SP_125059_\d{8}_02_T1").all()):
        raise ValueError("Expected distinct same-year Tier 1 Landsat 8/9 scenes from path 125 row 059")
    for col in SOURCE_FILES:
        for url in candidates[col]:
            asset = urlsplit(url)
            if (asset.scheme, asset.netloc, asset.query, asset.fragment) != (
                    "https", "landsateuwest.blob.core.windows.net", "", ""):
                raise ValueError("Candidate assets must use unsigned URLs from the selected Landsat mirror")
    return candidates.sort_values("datetime_utc").reset_index(drop=True)


def grid_design(transform, world, classes):
    row, col = np.nonzero(world)
    x = transform.c + (col + 0.5) * transform.a
    y = transform.f + (row + 0.5) * transform.e
    labels = classes[row, col]
    idx, centers, land = point_indices(x, y, 500, 0)
    n = len(land)
    classified = counts_per_cell(idx, np.isin(labels, VALID_ESRI), n)
    frame = pd.DataFrame({
        "cell_id": [f"{e}_{n}" for e, n in centers], "easting_m": centers[:, 0], "northing_m": centers[:, 1],
        "land_centers": land, "land_share": land / 2500, "classified_share": fractions(classified, land),
        "tree_fraction": fractions(counts_per_cell(idx, labels == 2, n), classified),
        "built_fraction": fractions(counts_per_cell(idx, labels == 7, n), classified),
        "esri_water_share": fractions(counts_per_cell(idx, labels == 1, n), land),
    })
    # Match the precision of the archived temperature panel before choosing pairs.
    frame[COVARIATES] = frame[COVARIATES].round(8)
    frame["area_eligible"] = area_rule(frame)
    return frame, (x, y, idx)


def primary_window(scene, signature, bounds, cache):
    path = cache / f"{scene['scene_id']}.npz"
    urls = [scene[key] for key in PRIMARY_BANDS]
    if not path.exists():
        dn, t = native_window(urls[0], signature, bounds, 0)
        qa, q = native_window(urls[1], signature, bounds, 1)
        if q != t or dn.shape != qa.shape:
            raise ValueError(f"Native thermal and QA grids differ for {scene['scene_id']}")
        with atomic_target(path) as tmp, tmp.open("wb") as f:
            np.savez_compressed(f, dn=dn, qa=qa, transform=tuple(t)[:6], urls=urls)
    with np.load(path, allow_pickle=False) as data:
        if data["urls"].tolist() != urls:
            raise ValueError(f"Cached URLs differ for {scene['scene_id']}")
        return data["dn"], data["qa"], Affine(*data["transform"])


def audit_scene(scene, signature, bounds, cache, scene_grid, cell_grid, points):
    dn, qa, affine = primary_window(scene, signature, bounds, cache)
    arrays = dict(zip((scene[key] for key in PRIMARY_BANDS), (dn, qa)))

    def cached_window(url, _signature, _bounds, _fill):
        return arrays[url], affine

    stamp = scene["datetime_utc"]
    date = stamp[:10]
    original = {"id": scene["scene_id"], "properties": {"datetime": stamp, "landsat:collection_category": "T1"},
                "assets": {key: {"href": scene[key]} for key in PRIMARY_BANDS}}
    row, _ = compute_scene(original, None, scene_grid, cached_window)
    metrics = {FIELDS[i]: row[i] for i in (0, 1, 3, 4, 6, 7, 8)}
    metrics["overpass_sgt"] = datetime.fromisoformat(stamp.replace("Z", "+00:00")).astimezone(LOCAL_TIME).isoformat()
    metrics["period"] = period(date)

    ok_qa, plausible, _ = mask_values(dn, qa)
    x, y, ids = points
    native = points_within_native(x, y, tuple(affine)[:6], dn.shape)
    area_passing = counts_per_cell(ids, (ok_qa & plausible)[native], len(cell_grid))
    included = cell_grid.area_eligible & (5 * area_passing >= 3 * cell_grid.land_centers)
    usable = cell_grid.loc[included].assign(date=date, valid_lst_centers=area_passing[included])
    return metrics, usable


def audit(year, frozen, output, workers):
    manifest = read_candidates(year, frozen)
    cover, cover_url = land_cover(year)
    transform, world, _, _, classes = land_layers(trim_to_30m=False, esri_path=cover)
    cell_grid, points = grid_design(transform, world, classes)
    scene_grid = spatial_grid()
    bounds = array_bounds(*world.shape, transform)
    _, old, old_sha = read_design_inputs(PANELS)
    previous_pairs = pd.read_csv(DESIGN / "pairs_primary.csv")
    pair_cells = set(previous_pairs.lower_tree_cell) | set(previous_pairs.higher_tree_cell)
    built_cells = set(pd.read_csv(BUILT_DESIGN / "fit_sample.csv").cell_id)
    contrast = json.loads((BUILT_DESIGN / "design_record.json").read_text())["chosen"]
    earlier_ids = set(pd.read_csv(DISCOVERY_SCENES).scene_id)
    old_cells = set(old.cell_id)

    cache = PRIMARY_CACHE / str(year)
    signature = (get_json(SAS)["token"] if any(not (cache / f"{s}.npz").exists() for s in manifest.scene_id)
                 else None)
    rows, availability = [], []
    eligible = []
    scenes = manifest.to_dict("records")

    def process(scene):
        return audit_scene(scene, signature, bounds, cache, scene_grid, cell_grid, points)

    with ThreadPoolExecutor(max_workers=workers) as pool:
        for result, cells in pool.map(process, scenes):
            candidate_pairs, _, _, _ = build_matches(cells, RULES["primary"])
            _, arms, shared, strict = selected_samples(assign_blocks(cells), contrast["tree_start_pct"],
                                                       contrast["lower_built_pct"], contrast["higher_built_pct"])
            is_prior = result["scene_id"] in earlier_ids
            passes = bool(result["passes_qa_and_range_60pct_44_areas"])
            result.update({
                "prior_discovery": is_prior, "usable_cells": len(cells), "prior_usable_cells": cells.cell_id.isin(old_cells).sum(),
                "prior_stage1_paired_cells": cells.cell_id.isin(pair_cells).sum(),
                "prior_stage2_model_cells": cells.cell_id.isin(built_cells).sum(),
                "stage1_disjoint_pairs": len(candidate_pairs), "stage2_lower_blocks": arms["lower"][0]["endpoint_blocks"],
                "stage2_higher_blocks": arms["higher"][0]["endpoint_blocks"], "stage2_shared_blocks": len(shared),
                "stage2_strict_shared_blocks": len(strict),
            })
            rows.append(result)
            availability.append(cells[[*KEY, "valid_lst_centers"]])
            if passes and not is_prior:
                eligible.append(result["scene_id"])
            print(result["date_utc"], result[FIELDS[6]], result[FIELDS[7]], len(cells),
                  "passes" if passes else "fails", flush=True)

    ledger = pd.DataFrame(rows)
    if year == 2025:
        prior = pd.read_csv(ARCHIVED_SCREEN).query("tier == 'T1'")
        if set(prior.scene_id) != set(manifest.scene_id):
            raise ValueError("2025 catalog has changed since the discovery freeze")
        checked = ledger.merge(prior, on="scene_id", suffixes=("", "_old"), validate="one_to_one")
        for key in (FIELDS[3], FIELDS[4], FIELDS[6], FIELDS[7], FIELDS[8]):
            np.testing.assert_allclose(checked[key], checked[f"{key}_old"], atol=1e-6, rtol=0)

    output.mkdir(parents=True, exist_ok=True)
    write_csv(manifest, output / "candidate_scenes.csv")
    write_csv(ledger, output / "scene_coverage.csv")
    write_csv(cell_grid, output / "grid_land_cover.csv")
    membership = pd.concat(availability, ignore_index=True)
    membership.to_csv(output / "usable_cell_dates.csv.gz", index=False, lineterminator="\n",
                      compression={"method": "gzip", "mtime": 0})
    write_csv(manifest[manifest.scene_id.isin(eligible)], output / "eligible_new_scenes.csv")
    counts = {}
    for season in SEASONS:
        g = ledger[ledger.period == season]
        counts[season] = {
            "candidates": len(g), "prior_discovery": int(g.prior_discovery.sum()),
            "eligible_new": int(g.scene_id.isin(eligible).sum()),
        }
    summary = {
        "year": year, "candidate_tier1_scenes": len(ledger),
        "prior_discovery_scenes": int(ledger.prior_discovery.sum()), "eligible_new_scene_ids": eligible,
        "area_eligible_cells": int(cell_grid.area_eligible.sum()), "prior_usable_cells": len(old_cells),
        "prior_paired_cells": len(pair_cells), "prior_built_model_cells": len(built_cells),
        "period_counts": counts,
        "stac_endpoint": STAC, "esri_cover_url": cover_url, "discovery_panel_sha256": old_sha,
        "spatial_sha256": {path.name: checksum(path) for path in (cover, BOUNDARY, RAW / WORLD_FILE)},
        "result_sha256": {p.name: checksum(p) for p in sorted(output.glob("*.csv*"))},
    }
    (output / "audit_record.json").write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--year", type=int, choices=(2023, 2024, 2025), required=True)
    parser.add_argument("--candidates", type=Path, help="Use a frozen candidate manifest instead of querying STAC")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()
    audit(args.year, args.candidates, args.output, args.workers)
