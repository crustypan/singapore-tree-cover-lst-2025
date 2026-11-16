"""Archive original optical and QA windows for the frozen Stage 7 cases."""
import argparse
from datetime import date, timedelta
from itertools import groupby
import json
from pathlib import Path
from urllib.parse import urlencode, urlsplit
from xml.etree import ElementTree as ET

import numpy as np
import pandas as pd
import rasterio
from rasterio.transform import Affine, array_bounds
from rasterio.warp import Resampling, transform_bounds
from rasterio.windows import from_bounds

from analyze_landsat import write_csv
from audit_heldout_overpasses import AUDIT, land_cover
from fetch_cover_audit_sentinel import SCL_SAS, clear_mask, read_window
from fetch_inputs import RAW, atomic_target, checksum, get_json, open_https
from make_landsat_panel import SOURCE_FILES, fetch_windows
from make_validation_panel import SNAPSHOTS
from prepare_influential_pair_audit import DESIGN, OPTICAL, SNAPSHOT
from prepare_local_comparison import PANELS
from prepare_overpass_validation import verify_hashes
from screen_landsat import BOUNDARY, GRID_CRS, SAS, STAC, land_layers, native_window


SENTINEL_BANDS = {"B02": 1, "B03": 2, "B04": 3, "B08": 7}
SENTINEL_HOST = "sentinel2l2a01.blob.core.windows.net"
PADDING = 4
CLOUD_WINDOW_DAYS = 30


def cell_center(cell_id):
    return tuple(map(int, cell_id.split("_")))


def extent(cell_ids, pad=0):
    xs, ys = zip(*map(cell_center, cell_ids))
    dist = 250 + pad
    return (min(xs) - dist, min(ys) - dist, max(xs) + dist, max(ys) + dist)


def grid_window(bounds, transform):
    region = from_bounds(*bounds, transform=transform)
    col, row = np.floor([region.col_off, region.row_off]).astype(int)
    end_col, end_row = np.ceil([region.col_off + region.width, region.row_off + region.height]).astype(int)
    return np.s_[row:end_row, col:end_col], transform * Affine.translation(col, row)


def save(path, **arrays):
    with atomic_target(path) as tmp, tmp.open("wb") as stream:
        np.savez_compressed(stream, **arrays)


def geography(sites, out):
    props = [feature["properties"] for feature in json.loads(BOUNDARY.read_text())["features"]]
    areas = dict(sorted((p["PLN_AREA_C"], p["PLN_AREA_N"]) for p in props))
    labels = np.array(["", *areas.values()])
    extents = {}
    for year, year_sites in sites.groupby(sites.date.str[:4].astype(int)):
        esri_path, _ = land_cover(year)
        audit = json.loads((AUDIT / str(year) / "audit_record.json").read_text())
        verify_hashes(RAW, audit["spatial_sha256"])
        transform, land, _, planning, classes = land_layers(trim_to_30m=False, esri_path=esri_path)
        extents[year] = array_bounds(*land.shape, transform)
        for site in year_sites.itertuples(index=False):
            chip = {"land": [], "classes": [], "transform": [], "planning_area": []}
            for cell_id in (site.lower_tree_cell, site.higher_tree_cell):
                crop, grid = grid_window(extent([cell_id]), transform)
                mask = land[crop]
                if mask.shape != (50, 50) or not mask.any():
                    raise ValueError(f"Invalid 10 m study area in {cell_id}")
                area = np.bincount(planning[crop][mask]).argmax()
                chip["land"].append(mask)
                chip["classes"].append(classes[crop])
                chip["transform"].append(tuple(grid)[:6])
                chip["planning_area"].append(labels[area])
            save(out / f"cover_{site.site}.npz", **chip)
    return extents


def thermal_source(row, year, grid_bounds):
    validation = year != 2025
    cache = RAW / "validation_windows" / "selected" / str(year) if validation else RAW / "landsat_selected"
    prov_file = (SNAPSHOTS / str(year) if validation else PANELS) / "provenance.json"
    provenance = json.loads(prov_file.read_text())
    source = cache / f"{row['scene_id']}.npz"
    if not source.exists():
        fetch_windows([row], grid_bounds, cache)
    key = "windows_sha256" if validation else "raw_windows_sha256"
    verify_hashes(cache, {source.name: provenance[key][row["scene_id"]]})
    return source


def landsat_windows(sources, cases, out, grid_bounds):
    signature = get_json(SAS)["token"]
    for row in sources.to_dict("records"):
        related = cases[cases.date.eq(row["date"])]
        cells = pd.concat([related.lower_tree_cell, related.higher_tree_cell]).unique()
        window = extent(cells, pad=90)
        year = int(row["date"][:4])
        source = thermal_source(row, year, grid_bounds[year])
        with np.load(source, allow_pickle=False) as originals:
            if originals["urls"].tolist() != [row[key] for key in SOURCE_FILES]:
                raise ValueError(f"Original Landsat bands differ for {row['scene_id']}")
            crop, grid = grid_window(window, Affine(*originals["source_transform"]))
            bands = {key: originals[key][crop] for key in SOURCE_FILES}
        for key in OPTICAL:
            values, transform = native_window(row[key], signature, window, 0)
            if grid != transform or values.shape != bands["lwir11"].shape:
                raise ValueError(f"Optical and thermal grids differ for {row['scene_id']}")
            bands[key] = values
        save(out / f"landsat_{row['date']}.npz", transform=tuple(grid)[:6], **bands)
        print(f"Archived Landsat window  {row['scene_id']}", flush=True)


def signed_sentinel(url, signature):
    parts = urlsplit(url)
    if (parts.scheme, parts.netloc, parts.query, parts.fragment) != ("https", SENTINEL_HOST, "", ""):
        raise ValueError("Unexpected Sentinel-2 source address")
    return url + "?" + signature


def sentinel_band(url, signature, transforms, shape, dtype):
    try:
        return np.stack([read_window("/vsicurl/" + signed_sentinel(url, signature), t, shape,
                                     Resampling.nearest, dtype) for t in transforms])
    except (rasterio.errors.RasterioError, OSError):
        raise RuntimeError(f"Failed to read Sentinel-2 source band {Path(urlsplit(url).path).name}") from None


def sentinel_candidates(day, bounds):
    delta = timedelta(days=CLOUD_WINDOW_DAYS)
    query = urlencode({
        "collections": "sentinel-2-l2a", "bbox": ",".join(f"{v:.6f}" for v in bounds),
        "datetime": f"{day - delta}T00:00:00Z/{day + delta}T23:59:59Z", "limit": 100,
    })
    items = get_json(f"{STAC}?{query}")["features"]
    rows = []
    for item in items:
        if "_T48NUG_" not in item["id"]:
            continue
        timestamp = item["properties"]["datetime"]
        offset = (date.fromisoformat(timestamp[:10]) - day).days
        if abs(offset) <= CLOUD_WINDOW_DAYS:
            rows.append({"scene_id": item["id"], "datetime_utc": timestamp, "days_from_landsat": offset,
                         **{key: item["assets"][key]["href"] for key in (*SENTINEL_BANDS, "SCL", "product-metadata")}})
    return sorted(rows, key=lambda a: (abs(a["days_from_landsat"]), a["scene_id"]))


def sentinel_calibration(url, signature):
    try:
        with open_https(signed_sentinel(url, signature)) as stream:
            root = ET.fromstring(stream.read())
    except OSError:
        raise RuntimeError("Failed to read Sentinel-2 product metadata") from None
    nodes = list(root.iter())
    quantity, = [float(e.text) for e in nodes if e.tag.rsplit("}", 1)[-1] == "BOA_QUANTIFICATION_VALUE"]
    offsets = {int(e.attrib["band_id"]): float(e.text) for e in nodes
               if e.tag.rsplit("}", 1)[-1] == "BOA_ADD_OFFSET"}
    return quantity, [offsets[i] for i in SENTINEL_BANDS.values()]


def sentinel_windows(sites, out, frozen):
    signature = get_json(SCL_SAS)["token"]
    prior = pd.read_csv(frozen).set_index("site").to_dict("index") if frozen else {}
    selected = []
    for site in sites.itertuples(index=False):
        with np.load(out / f"cover_{site.site}.npz") as geo:
            land = geo["land"]
            transforms = [Affine(*t) * Affine.translation(-PADDING, -PADDING) for t in geo["transform"]]
        shape = (land.shape[1] + 2 * PADDING, land.shape[2] + 2 * PADDING)
        bounds = transform_bounds(GRID_CRS, "EPSG:4326", *extent([site.lower_tree_cell, site.higher_tree_cell]))
        candidates = [prior[site.site]] if frozen else sentinel_candidates(date.fromisoformat(site.date), bounds)
        if not candidates:
            raise ValueError(f"No Sentinel-2 products found for site {site.site}")
        for _, same_day_offset in groupby(candidates, lambda row: abs(row["days_from_landsat"])):
            ranked = []
            for candidate in same_day_offset:
                scl = sentinel_band(candidate["SCL"], signature, transforms, shape, "uint8")
                clear = np.stack([clear_mask(p)[PADDING:-PADDING, PADDING:-PADDING] for p in scl])
                fractions = (land & clear).sum(axis=(1, 2)) / land.sum(axis=(1, 2))
                if fractions.min() >= 0.80:
                    ranked.append((fractions.min(), candidate, scl, fractions))
            if ranked:
                _, choice, scl, fractions = min(ranked, key=lambda e: (-e[0], e[1]["scene_id"]))
                break
        else:
            raise ValueError(f"No sufficiently clear Sentinel-2 image within 30 days for site {site.site}")
        quant, offsets = sentinel_calibration(choice["product-metadata"], signature)
        bands = {name: sentinel_band(choice[name], signature, transforms, shape, "uint16")
                 for name in SENTINEL_BANDS}
        save(out / f"sentinel_{site.site}.npz", scl=scl, quantification=quant, offsets=offsets,
             transform=np.array([tuple(t)[:6] for t in transforms]), **bands)
        fields = ["scene_id", "datetime_utc", "days_from_landsat", *SENTINEL_BANDS, "SCL", "product-metadata"]
        selected.append({
            "site": site.site, "landsat_reference_date": site.date, **{k: choice[k] for k in fields},
            "clear_lower": fractions[0], "clear_higher": fractions[1],
            "boa_quantification": quant, **{f"boa_offset_{key}": val for key, val in zip(SENTINEL_BANDS, offsets)},
        })
        print(f"Archived Sentinel-2 context for {site.site}  {choice['scene_id']}", flush=True)
    write_csv(pd.DataFrame(selected), out / "sentinel_sources.csv")


def create(design, output, frozen):
    record = json.loads((design / "design_record.json").read_text())
    verify_hashes(design, record["files_sha256"])
    sources = pd.read_csv(design / "landsat_sources.csv")
    cases = pd.read_csv(design / "pair_registry.csv")
    sites = pd.read_csv(design / "review_sites.csv")
    output.mkdir(parents=True, exist_ok=True)
    bounds = geography(sites, output)
    landsat_windows(sources, cases, output, bounds)
    sentinel_windows(sites, output, frozen)
    metadata = {
        "design_sha256": checksum(design / "design_record.json"), "grid_crs": GRID_CRS,
        "file_sha256": {p.name: checksum(p) for p in sorted(output.glob("*.npz"))},
    }
    metadata["file_sha256"]["sentinel_sources.csv"] = checksum(output / "sentinel_sources.csv")
    (output / "provenance.json").write_text(json.dumps(metadata, indent=2) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--design", type=Path, default=DESIGN)
    parser.add_argument("--output", type=Path, default=SNAPSHOT)
    parser.add_argument("--sentinel-selection", type=Path)
    args = parser.parse_args()
    create(args.design, args.output, args.sentinel_selection)
