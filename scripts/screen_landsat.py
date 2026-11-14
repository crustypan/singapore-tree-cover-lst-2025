"""Screen Landsat 8/9 scenes for a potential selected-date 2025 study."""
import argparse
import json
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.parse import urlencode, urlsplit

import numpy as np
import rasterio
from rasterio.features import rasterize
from rasterio.transform import Affine, array_bounds
from rasterio.warp import Resampling, reproject, transform_geom
from rasterio.windows import Window, from_bounds
from shapely.geometry import mapping, shape
from shapely.validation import make_valid

from audit_coverage import BOUNDARY, ESRI, WORLD_CLASSES, write_rows
from fetch_inputs import BBOX, RAW, READ_ENV, WORLD_FILE, atomic_target, get_json

STAC = "https://planetarycomputer.microsoft.com/api/stac/v1/search"
SAS = "https://planetarycomputer.microsoft.com/api/sas/v1/token/landsateuwest/landsat-c2"
GRID_CRS = "EPSG:32648"
SCALE = 0.00341802
OFFSET = 149.0 - 273.15


def catalog():
    query = urlencode({
        "collections": "landsat-c2-l2",
        "bbox": ",".join(map(str, BBOX)),
        "datetime": "2025-01-01T00:00:00Z/2025-12-31T23:59:59Z",
        "limit": 100,
    })
    url = f"{STAC}?{query}"
    results = []
    while url:
        page = get_json(url)
        results.extend(page["features"])
        url = next((link["href"] for link in page["links"] if link["rel"] == "next"), None)
    # The collection-items endpoint returned out-of-year records. Match the exact
    # acquisition date even when querying the dedicated STAC search endpoint.
    scenes = [
        f for f in results
        if f["properties"]["platform"] in ("landsat-8", "landsat-9")
        and f["properties"].get("landsat:correction") == "L2SP"
        and f["properties"]["datetime"].startswith("2025-")
    ]
    if not scenes or len({f["id"] for f in scenes}) != len(scenes):
        raise ValueError("Empty Landsat catalog or duplicate scene IDs")
    return sorted(scenes, key=lambda scene: scene["id"])


def land_layers(trim_to_30m):
    features = json.loads(BOUNDARY.read_text())["features"]
    names = sorted({f["properties"]["PLN_AREA_C"] for f in features})
    if len(names) != 55:
        raise ValueError("Expected 55 URA planning areas")
    idx = {name: index for index, name in enumerate(names, start=1)}
    geometries = []
    for feature in features:
        repaired = make_valid(shape(feature["geometry"]))
        geometry = transform_geom("EPSG:4326", GRID_CRS, mapping(repaired))
        geometries.append((geometry, idx[feature["properties"]["PLN_AREA_C"]]))
    with rasterio.open(ESRI) as esri:
        h, w = esri.shape
        if trim_to_30m:
            h, w = 3 * (h // 3), 3 * (w // 3)
        ds_transform = esri.transform
        esri_values = esri.read(1)[:h, :w]
    world = np.zeros((h, w), dtype=np.uint8)
    with rasterio.open(RAW / WORLD_FILE) as ds:
        reproject(rasterio.band(ds, 1), world, dst_crs=GRID_CRS, dst_transform=ds_transform,
                  resampling=Resampling.nearest, dst_nodata=0)
    ura_idx = rasterize(geometries, out_shape=(h, w), transform=ds_transform, dtype="uint8")
    in_ura = ura_idx > 0
    if not set(np.unique(world[in_ura])).issubset(WORLD_CLASSES):
        raise ValueError("WorldCover contains an unclassified area or unexpected class within URA")
    return ds_transform, in_ura & (world != 80), in_ura & (esri_values != 1), ura_idx, esri_values


def spatial_grid():
    ds_transform, world, esri, ura_idx, _ = land_layers(trim_to_30m=True)
    h, w = ura_idx.shape
    # Store the nine 10 m centers for each 30 m cell so we can credit each area to its
    # own planning area at boundaries.
    world_cells = world.reshape(h // 3, 3, w // 3, 3)
    esri_cells = esri.reshape(h // 3, 3, w // 3, 3)
    world_weight = world_cells.sum(axis=(1, 3), dtype=np.uint8)
    esri_weight = esri_cells.sum(axis=(1, 3), dtype=np.uint8)
    land_planning_areas = np.where(world, ura_idx, 0)
    land_ids = np.nonzero(land_planning_areas.ravel())[0]
    area_codes = land_planning_areas.ravel()[land_ids]
    area_totals = np.bincount(area_codes, minlength=56)[1:]
    if np.any(area_totals == 0):
        raise ValueError("An URA planning area has no WorldCover land pixels")
    subrow, subcol = np.unravel_index(land_ids, (h, w))
    cell_indices = (subrow // 3) * (w // 3) + subcol // 3
    return (
        ds_transform * Affine.scale(3), world_weight, esri_weight, area_codes,
        area_totals, cell_indices,
    )


def native_window(url, public_signature, bounds, fill, dtype="uint16"):
    location = urlsplit(url)
    if location.scheme != "https" or location.netloc != "landsateuwest.blob.core.windows.net":
        raise ValueError("Unexpected Landsat asset provider")
    try:
        with rasterio.Env(**READ_ENV), rasterio.open("/vsicurl/" + url + "?" + public_signature) as ds:
            if (ds.count != 1 or ds.dtypes[0] != dtype or ds.res != (30, 30)
                    or ds.crs is None or ds.crs.to_epsg() != 32648):
                raise ValueError(f"Unexpected Landsat raster format for {Path(location.path).name}")
            rough = from_bounds(*bounds, transform=ds.transform)
            col, row = np.floor([rough.col_off, rough.row_off]).astype(int)
            right, bottom = np.ceil([rough.col_off + rough.width, rough.row_off + rough.height]).astype(int)
            native = Window(int(col), int(row), int(right - col), int(bottom - row))
            source = ds.read(1, window=native, boundless=True, fill_value=fill)
            return source, ds.window_transform(native)
    except rasterio.errors.RasterioError:
        raise RuntimeError(f"Could not read {Path(location.path).name}") from None


def mask_values(dn, qa):
    # Bits 0–5 are fill, dilated cloud, cirrus, cloud, cloud shadow, and snow.
    valid_qa = (dn > 0) & ((qa & 63) == 0)
    temperatures = dn * SCALE + OFFSET
    plausible = (temperatures >= 15) & (temperatures <= 70)
    return valid_qa, plausible, temperatures


def masked_scene(scene, public_signature, grid_transform, extent):
    results = {}
    for key, fill in {"lwir11": 0, "qa_pixel": 1}.items():
        bounds = array_bounds(*extent, grid_transform)
        source, transform = native_window(scene["assets"][key]["href"], public_signature, bounds, fill)
        aligned = np.full(extent, fill, dtype=np.uint16)
        reproject(source, aligned, src_crs=GRID_CRS, src_transform=transform,
                  dst_crs=GRID_CRS, dst_transform=grid_transform,
                  resampling=Resampling.nearest, src_nodata=fill, dst_nodata=fill)
        results[key] = aligned
    return mask_values(results["lwir11"], results["qa_pixel"])


def native_footprints(selected, signature, output):
    transform, world_land, esri_land, ura, _ = land_layers(trim_to_30m=False)
    bounds = array_bounds(*ura.shape, transform)
    area_totals = np.bincount(ura[world_land], minlength=56)[1:]
    hits = np.zeros(ura.shape, dtype=np.uint8)
    per_scene = []
    for scene in sorted(selected, key=lambda row: row["date_utc"]):
        dn, t = native_window(scene["st_b10_source"], signature, bounds, 0)
        qc, qc_transform = native_window(scene["qa_pixel_source"], signature, bounds, 1)
        if t != qc_transform or dn.shape != qc.shape:
            raise ValueError(f"ST_B10 and QA_PIXEL have different footprints for {scene['scene_id']}")
        qa_pass, plausible, _ = mask_values(dn, qc)
        fine = np.zeros_like(hits)
        reproject((qa_pass & plausible).astype(np.uint8), fine, src_crs=GRID_CRS, src_transform=t,
                  dst_crs=GRID_CRS, dst_transform=transform, resampling=Resampling.nearest)
        hits += fine
        areas = np.bincount(ura[world_land & (fine > 0)], minlength=56)[1:]
        per_scene.append((scene["date_utc"], round(float(100 * fine[world_land].mean()), 8),
                          int((10 * areas >= 3 * area_totals).sum())))
    write_rows(output / "native_footprint_scenes.csv",
               ("date_utc", "worldcover_land_pct", "planning_areas_at_30pct"), per_scene)
    write_rows(
        output / "repeat_coverage.csv", ("min_dates", "worldcover_land_pct", "esri_land_pct"),
        [
            (i, 100 * (hits[world_land] >= i).mean(), 100 * (hits[esri_land] >= i).mean())
            for i in range(1, len(selected) + 1)
        ],
    )


def compute_scene(scene, signature, grid):
    transform, weights, esri, area_codes, area_totals, cell_indices = grid
    suitable, plausible, temperatures = masked_scene(scene, signature, transform, weights.shape)
    used = suitable & plausible

    def measures(mask):
        district = np.bincount(
            area_codes, weights=mask.ravel()[cell_indices], minlength=len(area_totals) + 1,
        )[1:]
        return (
            float(100 * weights[mask].sum() / weights.sum()),
            int((10 * district >= 3 * area_totals).sum()),
        )

    qa_pct, qa_districts = measures(suitable)
    area_pct, districts = measures(used)
    data = temperatures[suitable & (weights > 0)]
    properties = scene["properties"]
    date = properties["datetime"][:10]
    tier = properties["landsat:collection_category"]
    record = (
        scene["id"], date, tier, round(qa_pct, 6), qa_districts,
        int(0.6 <= qa_pct / 100 and qa_districts >= 44),
        round(area_pct, 6), districts, int(0.6 <= area_pct / 100 and districts >= 44),
        round(float(np.median(data)), 4) if data.size else None,
        round(float(100 * (data < 15).sum() / data.size), 4) if data.size else None,
        scene["assets"]["lwir11"]["href"], scene["assets"]["qa_pixel"]["href"],
    )
    return record, used


FIELDS = (
    "scene_id", "date_utc", "tier", "qa_worldcover_land_pct", "qa_planning_areas_at_30pct",
    "passes_qa_only_60pct_44_areas", "qa_and_range_worldcover_land_pct",
    "qa_and_range_planning_areas_at_30pct", "passes_qa_and_range_60pct_44_areas",
    "qa_only_land_lst_median_c", "qa_only_land_samples_below_15c_pct",
    "st_b10_source", "qa_pixel_source",
)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    scenes = catalog()
    signature = get_json(SAS)["token"]
    grid = spatial_grid()
    _, weights, esri, _, _, _ = grid
    union = np.zeros(weights.shape, dtype=bool)
    rows = []
    with ThreadPoolExecutor(max_workers=4) as ex:
        for record, valid in ex.map(lambda scene: compute_scene(scene, signature, grid), scenes):
            rows.append(record)
            if record[2] == "T1" and record[8]:
                union |= valid
            print(record[0], record[6], record[7], flush=True)
    with atomic_target(args.output / "scene_coverage.csv") as tmp:
        write_rows(tmp, FIELDS, sorted(rows, key=lambda row: (row[1], row[0])))
    best = [dict(zip(FIELDS, row)) for row in sorted(rows, key=lambda row: -row[6]) if row[2] == "T1"]
    summary = {
        "candidate_scenes": len(rows),
        "tier1_count": len(best),
        "tier1_passing_raw_qa_and_temperature": [row for row in best if row[FIELDS[8]]],
        "tier1_screened_union_worldcover_pct": float(100 * weights[union].sum() / weights.sum()),
        "tier1_screened_union_esri_pct": float(100 * esri[union].sum() / esri.sum()),
        "highest_coverage_tier1_scenes": best[:5],
    }
    with atomic_target(args.output / "summary.json") as tmp:
        tmp.write_text(json.dumps(summary, indent=2) + "\n")
    print("Tier 1 candidate dates", [row["date_utc"] for row in summary["tier1_passing_raw_qa_and_temperature"]])
    print("Tier 1 screened union", summary["tier1_screened_union_worldcover_pct"])
    native_footprints(summary["tier1_passing_raw_qa_and_temperature"], signature, args.output)


if __name__ == "__main__":
    main()
