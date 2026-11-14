"""Aggregate the seven selected Landsat overpasses and land-cover areas to fixed UTM grids."""
import argparse
import csv
import json
from pathlib import Path

import numpy as np
from rasterio.transform import array_bounds

from audit_coverage import BOUNDARY, ESRI, write_rows
from fetch_inputs import RAW, WORLD_FILE, atomic_target, checksum, get_json
from screen_landsat import GRID_CRS, SAS, catalog, land_layers, mask_values, native_window

SCENES = Path("results/landsat_screening/summary.json")
SOURCE_FILES = {"lwir11": (0, "uint16"), "qa_pixel": (1, "uint16"),
                "cdist": (-9999, "int16"), "qa": (-9999, "int16")}
DATES = ("2025-01-03", "2025-02-04", "2025-05-03", "2025-07-06",
         "2025-10-02", "2025-10-18", "2025-12-05")
GRID_VARIANTS = {
    "world_500": (500, 0, "WorldCover"),
    "esri_500": (500, 0, "Esri"),
    "world_500_shift250": (500, 250, "WorldCover"),
    "world_1000": (1000, 0, "WorldCover"),
}
VALID_ESRI = (2, 4, 5, 7, 8, 9, 11)


def freeze_scenes(output):
    records = json.loads(SCENES.read_text())["tier1_passing_raw_qa_and_temperature"]
    chosen = {r["date_utc"]: r for r in records}
    if sorted(chosen) != list(DATES):
        raise ValueError("The archive must contain the seven authorized dates")
    lookup = {s["id"]: s for s in catalog()}
    rows = []
    for date in DATES:
        original = chosen[date]
        scene = lookup[original["scene_id"]]
        props = scene["properties"]
        urls = {key: scene["assets"][key]["href"] for key in SOURCE_FILES}
        if (
            props["datetime"][:10] != date or props["landsat:collection_category"] != "T1"
            or urls["lwir11"] != original["st_b10_source"]
            or urls["qa_pixel"] != original["qa_pixel_source"]
        ):
            raise ValueError(f"Scene metadata changed for {date}")
        rows.append((original["scene_id"], props["datetime"], props["platform"],
                     props.get("view:sun_elevation"), *urls.values()))
    with atomic_target(output) as tmp:
        write_rows(tmp, ("scene_id", "datetime_utc", "platform", "sun_elevation_deg", *SOURCE_FILES), rows)


def read_frozen_scenes(source):
    with source.open(newline="") as f:
        rows = list(csv.DictReader(f))
    if [r["datetime_utc"][:10] for r in rows] != list(DATES):
        raise ValueError("The manifest differs from the seven authorized dates")
    if any(not r["scene_id"].endswith("_02_T1") for r in rows):
        raise ValueError("The selected scenes must be Collection 2 Tier 1")
    archived = json.loads(SCENES.read_text())["tier1_passing_raw_qa_and_temperature"]
    if {r["scene_id"] for r in rows} != {r["scene_id"] for r in archived}:
        raise ValueError("The manifest scene IDs differ from the archived screening")
    return rows


def fetch_windows(scenes, bounds, cache):
    signature = None
    for scene in scenes:
        path = cache / (scene["scene_id"] + ".npz")
        if not path.exists():
            if signature is None:
                signature = get_json(SAS)["token"]
            bands, source_transform = {}, None
            for key, (fill, dtype) in SOURCE_FILES.items():
                values, native = native_window(scene[key], signature, bounds, fill, dtype)
                if source_transform is not None and (native != source_transform or values.shape != bands["lwir11"].shape):
                    raise ValueError(f"Native raster grids differ for {scene['scene_id']}")
                source_transform = native
                bands[key] = values
            with atomic_target(path) as tmp, tmp.open("wb") as out:
                np.savez_compressed(out, source_transform=tuple(source_transform)[:6],
                                    urls=[scene[key] for key in SOURCE_FILES], **bands)
        print(f"Landsat window  {scene['scene_id']}", flush=True)


def point_indices(x, y, grid_size, shift):
    ix = np.floor((x - shift) / grid_size).astype("int32")
    iy = np.floor((y - shift) / grid_size).astype("int32")
    width = int(ix.max() - ix.min()) + 1
    codes = (iy - iy.min()) * width + ix - ix.min()
    groups, counts = np.unique(codes, return_counts=True)
    centers = np.column_stack((
        ((groups % width) + ix.min() + 0.5) * grid_size + shift,
        ((groups // width) + iy.min() + 0.5) * grid_size + shift,
    )).astype("int32")
    idx = np.full(int(codes.max()) + 1, -1, dtype="int32")
    idx[groups] = np.arange(len(groups))
    return idx[codes], centers, counts


def points_within_native(x, y, transform, shape):
    x_size, bx, x0, dy, y_size, y0 = transform
    if not (x_size == 30 and y_size == -30 and bx == dy == 0):
        raise ValueError("Expected a north-up Landsat 30 m UTM grid")
    cols = np.floor((x - x0) / x_size).astype("int32")
    rows = np.floor((y - y0) / y_size).astype("int32")
    if not ((rows >= 0) & (rows < shape[0]) & (cols >= 0) & (cols < shape[1])).all():
        raise ValueError("Landsat window excludes some of the 10 m land centers")
    return rows, cols


def counts_per_cell(indices, mask, size):
    return np.bincount(indices[mask], minlength=size)


def fractions(num, den):
    return np.divide(num, den, out=np.full(np.shape(den), np.nan), where=den > 0)


def aggregate(path, manifest, grid_size, shift, land_mask, data, cache):
    fine_transform, main_land, esri_land, esri_values = data
    row, col = np.nonzero(main_land if land_mask == "WorldCover" else esri_land)
    x = fine_transform.c + (col + 0.5) * fine_transform.a
    y = fine_transform.f + (row + 0.5) * fine_transform.e
    classes = esri_values[row, col]
    valid_class = np.isin(classes, VALID_ESRI)
    tree, built = classes == 2, classes == 7
    indices, centers, land = point_indices(x, y, grid_size, shift)
    size = len(land)
    classified = counts_per_cell(indices, valid_class, size)
    fixed_tree = fractions(counts_per_cell(indices, tree, size), classified)
    fixed_built = fractions(counts_per_cell(indices, built, size), classified)
    classified_share = fractions(classified, land)
    land_share = land / (grid_size / 10) ** 2
    ids = [f"{easting}_{northing}" for easting, northing in centers]
    records = []
    for scene in manifest:
        src = cache / (scene["scene_id"] + ".npz")
        with np.load(src, allow_pickle=False) as v:
            if v["urls"].tolist() != [scene[key] for key in SOURCE_FILES]:
                raise ValueError(f"Cached source IDs differ for {scene['scene_id']}")
            native = points_within_native(x, y, v["source_transform"], v["lwir11"].shape)
            dn, qa_pixel = v["lwir11"][native], v["qa_pixel"][native]
            cdist, stqa = v["cdist"][native], v["qa"][native]
        qa_pass, plausible, t = mask_values(dn, qa_pixel)
        passing = qa_pass & plausible
        strict = passing & (cdist >= 50) & (stqa >= 0) & (stqa <= 200)
        valid_class_obs = passing & valid_class
        usable = counts_per_cell(indices, passing, size)
        strict_usable = counts_per_cell(indices, strict, size)
        obs_classified = counts_per_cell(indices, valid_class_obs, size)
        obs_tree = fractions(counts_per_cell(indices, passing & tree, size), obs_classified)
        obs_built = fractions(counts_per_cell(indices, passing & built, size), obs_classified)
        tv = fractions(np.bincount(indices[passing], weights=t[passing], minlength=size), usable)
        ts = fractions(np.bincount(indices[strict], weights=t[strict], minlength=size), strict_usable)
        coverage = fractions(usable, land)
        strict_cover = fractions(strict_usable, land)
        if (strict_usable > usable).any() or (usable > land).any():
            raise ValueError("QA-filtered counts exceed the corresponding candidate areas")
        date = scene["datetime_utc"][:10]
        for j, (easting, northing) in enumerate(centers):
            records.append((date, scene["scene_id"], ids[j], easting, northing, land[j], classified[j],
                            land_share[j], classified_share[j], fixed_tree[j], fixed_built[j], usable[j],
                            coverage[j], tv[j], strict_usable[j], strict_cover[j], ts[j],
                            obs_classified[j], obs_tree[j], obs_built[j]))
    with atomic_target(path) as tmp, tmp.open("w", newline="") as out:
        writer = csv.writer(out, lineterminator="\n")
        writer.writerow((
            "date", "scene_id", "cell_id", "easting_m", "northing_m", "land_centers",
            "classified_centers", "land_share", "classified_share", "tree_fraction", "built_fraction",
            "valid_lst_centers", "valid_lst_share", "lst_c", "strict_lst_centers", "strict_lst_share",
            "strict_lst_c", "observed_classified_centers", "observed_tree_fraction", "observed_built_fraction",
        ))
        for record in records:
            writer.writerow([
                (f"{value:.8f}" if np.isfinite(value) else "")
                if isinstance(value, np.floating) else value for value in record
            ])
    print(f"{path.name}  {size} cells × {len(manifest)} dates", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scene-manifest", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--cache", type=Path, default=RAW / "landsat_selected")
    parser.add_argument("--freeze-scenes", action="store_true")
    args = parser.parse_args()
    if args.freeze_scenes:
        if args.scene_manifest.exists():
            raise FileExistsError(args.scene_manifest)
        freeze_scenes(args.scene_manifest)
    scene_list = read_frozen_scenes(args.scene_manifest)
    affine, world, esri_land, _, esri_values = land_layers(trim_to_30m=False)
    if (affine.a, affine.e, affine.b, affine.d) != (10.0, -10.0, 0.0, 0.0):
        raise ValueError("The area grid must be a north-up 10 m UTM grid")
    fetch_windows(scene_list, array_bounds(*world.shape, affine), args.cache)
    args.output.mkdir(parents=True, exist_ok=True)
    for key, (size, shift, mask) in GRID_VARIANTS.items():
        aggregate(args.output / (key + ".csv"), scene_list, size, shift, mask,
                  (affine, world, esri_land, esri_values), args.cache)
    provenance = {
        "crs": GRID_CRS, "input_sha256": {
            p.name: checksum(p) for p in (BOUNDARY, ESRI, RAW / WORLD_FILE, args.scene_manifest)
        }, "raw_windows_sha256": {
            s["scene_id"]: checksum(args.cache / (s["scene_id"] + ".npz")) for s in scene_list
        }, "panel_sha256": {key: checksum(args.output / (key + ".csv")) for key in GRID_VARIANTS},
    }
    with atomic_target(args.output / "provenance.json") as tmp:
        tmp.write_text(json.dumps(provenance, indent=2) + "\n")


if __name__ == "__main__":
    main()
