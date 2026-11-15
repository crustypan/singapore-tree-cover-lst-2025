"""Freeze Sentinel-2 acquisitions by local cloud coverage and cache clipped L1C data."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import date
import json
from pathlib import Path
import shutil
from tempfile import NamedTemporaryFile
from urllib.parse import quote, urlencode, urlsplit
from xml.etree import ElementTree

import numpy as np
import pandas as pd
import rasterio
from rasterio.transform import Affine, array_bounds, from_origin
from rasterio.warp import Resampling, reproject
from rasterio.windows import Window, from_bounds

from analyze_landsat import write_csv
from fetch_inputs import RAW, READ_ENV, atomic_target, checksum, get_json, open_https
from prepare_cover_support_audit import AUDIT_DESIGN, support_inputs
from prepare_overpass_validation import YEARS, verify_hashes
from screen_landsat import GRID_CRS, STAC, land_layers


CACHE = RAW / "cover_support_audit"
BANDS = ("B02", "B03", "B04", "B05", "B06", "B07", "B08", "B11", "B12")
ALL_BANDS = ("B01", "B02", "B03", "B04", "B05", "B06", "B07", "B08", "B8A", "B09", "B10", "B11", "B12")
SCL_SAS = "https://planetarycomputer.microsoft.com/api/sas/v1/token/sentinel2l2a01/sentinel2-l2"
GCS_LIST = "https://storage.googleapis.com/storage/v1/b/gcp-public-data-sentinel-2/o"
GCS_BASE = "https://storage.googleapis.com/gcp-public-data-sentinel-2/"
PAD_METERS = 300
MAX_SCENES = 5
MIN_SEPARATION_DAYS = 7


def study_grid(design, cache):
    path = cache / "geometry.npz"
    record_sha = checksum(design / "design_record.json")
    if not path.exists():
        g, _, _, _ = support_inputs()
        g = g.sort_values("cell_id").reset_index(drop=True)
        e, n = g.easting_m, g.northing_m
        offset = 250 + PAD_METERS
        west, south, east, north = e.min() - offset, n.min() - offset, e.max() + offset, n.max() + offset
        transform = from_origin(west, north, 10, 10)
        dims = ((north - south) // 10, (east - west) // 10)
        source, land, _, _, _ = land_layers(trim_to_30m=False)
        aligned = np.zeros(dims, dtype=np.uint8)
        reproject(land.astype(np.uint8), aligned, src_transform=source, src_crs=GRID_CRS,
                  dst_transform=transform, dst_crs=GRID_CRS, resampling=Resampling.nearest)
        indexes = np.full(dims, -1, dtype=np.int16)
        rows = ((north - (n + 250)) // 10).to_numpy()
        cols = (((e - 250) - west) // 10).to_numpy()
        for i, (y, x) in enumerate(zip(rows, cols)):
            indexes[y:y + 50, x:x + 50] = np.where(aligned[y:y + 50, x:x + 50] > 0, i, -1)
        np.testing.assert_array_equal(np.bincount(indexes[indexes >= 0]), g.land_centers)
        with atomic_target(path) as tmp, tmp.open("wb") as f:
            np.savez_compressed(f, transform=tuple(transform)[:6], cell_index=indexes,
                                names=g.cell_id.to_numpy(dtype="U"), rows=rows, cols=cols,
                                land_centers=g.land_centers.to_numpy(), design_sha256=record_sha)
    with np.load(path, allow_pickle=False) as data:
        if str(data["design_sha256"]) != record_sha:
            raise ValueError("Optical study grid was cached for a different design")
        return {k: data[k] for k in data.files}


def candidate_items(year):
    query = urlencode({"collections": "sentinel-2-l2a", "bbox": "103.60,1.15,104.13,1.49", "limit": 100,
                       "datetime": f"{year}-06-01T00:00:00Z/{year}-09-30T23:59:59Z"})
    url = f"{STAC}?{query}"
    candidates = {}
    while url:
        result = get_json(url)
        for item in sorted(result["features"], key=lambda x: x["id"]):
            props = item["properties"]
            identity = item["id"].split("_")
            if (identity[4] != "T48NUG" or not props["datetime"].startswith(str(year))
                    or props["eo:cloud_cover"] > 65):
                continue
            key = identity[0], identity[2]
            candidates[key] = {"scene_id": item["id"], "datetime_utc": props["datetime"],
                               "tile_cloud_pct": props["eo:cloud_cover"],
                               "scl_url": item["assets"]["SCL"]["href"]}
        url = next((link["href"] for link in result["links"] if link["rel"] == "next"), None)
    if not candidates:
        raise ValueError(f"No Sentinel-2 candidates for {year}")
    return sorted(candidates.values(), key=lambda r: r["datetime_utc"])


def read_window(source, transform, shape, method, dtype):
    with rasterio.Env(**READ_ENV), rasterio.open(source) as ds:
        if ds.crs is None or ds.crs.to_epsg() != 32648 or ds.count != 1 or ds.res[0] not in (10, 20):
            raise ValueError("Unexpected Sentinel-2 raster grid")
        rough = from_bounds(*array_bounds(*shape, transform), transform=ds.transform)
        x, y = np.floor([rough.col_off, rough.row_off]).astype(int) - 2
        end_x, end_y = np.ceil([rough.col_off + rough.width, rough.row_off + rough.height]).astype(int) + 2
        window = Window(x, y, end_x - x, end_y - y).intersection(Window(0, 0, ds.width, ds.height))
        values = ds.read(1, window=window)
        dest = np.zeros(shape, dtype=dtype)
        reproject(values, dest, src_transform=ds.window_transform(window), src_crs=ds.crs, src_nodata=0,
                  dst_transform=transform, dst_crs=GRID_CRS, dst_nodata=0, resampling=method)
    return dest


def clear_mask(scl):
    initial = np.isin(scl, [4, 5, 6])
    valid = initial.copy()
    for y, x in [(-2, 0), (2, 0), (0, -2), (0, 2), (-2, -2), (-2, 2), (2, -2), (2, 2)]:
        valid &= np.roll(initial, (y, x), axis=(0, 1))
    valid[:2] = valid[-2:] = False
    valid[:, :2] = valid[:, -2:] = False
    return valid


def get_scl(scene, cache, transform, dims, signature):
    url = scene["scl_url"]
    parts = urlsplit(url)
    if (parts.scheme, parts.netloc, parts.query) != ("https", "sentinel2l2a01.blob.core.windows.net", ""):
        raise ValueError("Unexpected unsigned Sentinel-2 SCL URL")
    if "/" in scene["scene_id"] or ".." in scene["scene_id"]:
        raise ValueError("Invalid Sentinel-2 scene ID")
    path = cache / "scl" / f"{scene['scene_id']}.npz"
    if not path.exists():
        try:
            scl = read_window("/vsicurl/" + url + "?" + signature, transform, dims, Resampling.nearest, "uint8")
        except rasterio.errors.RasterioError:
            raise RuntimeError(f"Could not read SCL for {scene['scene_id']}") from None
        with atomic_target(path) as tmp, tmp.open("wb") as f:
            np.savez_compressed(f, scl=scl, url=url)
    with np.load(path, allow_pickle=False) as data:
        if str(data["url"]) != url:
            raise ValueError(f"SCL source changed for {scene['scene_id']}")
        return data["scl"]


def l1c_sources(scene):
    platform, _, capture, orbit, tile, _ = scene["scene_id"].split("_")
    prefix = f"tiles/48/N/UG/{platform}_MSIL1C_{capture}"
    query = urlencode({"prefix": prefix, "maxResults": 500, "fields": "items(name),nextPageToken"})
    response = get_json(f"{GCS_LIST}?{query}")
    if response.get("nextPageToken"):
        raise RuntimeError(f"Too many L1C metadata hits for {scene['scene_id']}")
    products = {}
    for record in response.get("items", []):
        name = record["name"]
        if f"_{orbit}_{tile}_" not in name or ".SAFE/" not in name:
            continue
        product = name.split(".SAFE/", 1)[0]
        assets = products.setdefault(product, {})
        if name.endswith("/MTD_MSIL1C.xml"):
            assets["metadata"] = name
        elif "/IMG_DATA/" in name and name.endswith(".jp2"):
            band = Path(name).stem.rsplit("_", 1)[-1]
            assets[band] = name
    for product in sorted(products, reverse=True):
        assets = products[product]
        if not all(key in assets for key in ("metadata", *BANDS)):
            continue
        urls = {k: GCS_BASE + quote(assets[k], safe="/") for k in ("metadata", *BANDS)}
        with open_https(urls["metadata"]) as source:
            doc = ElementTree.parse(source)
        quant = float(doc.findtext(".//QUANTIFICATION_VALUE"))
        offsets = {ALL_BANDS[int(el.attrib["band_id"])]: int(el.text) for el in doc.findall(".//RADIO_ADD_OFFSET")}
        if quant != 10000 or not all(k in offsets for k in BANDS):
            raise ValueError(f"Unexpected L1C radiometry for {product}")
        return {"l1c_product": Path(product).name, "l1c_assets": urls, "quantification_value": quant,
                "offsets": {k: offsets[k] for k in BANDS}}
    return None


def choose(candidates, fractions, names, cases):
    weights = np.where(np.isin(names, cases), 5.0, 1.0)
    clear = fractions >= 0.80
    counts = np.zeros(len(names), dtype=int)
    remaining = list(range(len(candidates)))
    selected = []
    while remaining and len(selected) < MAX_SCENES:
        def rank(i):
            priority = (1.0 + (counts < 3)) * weights
            return (float(priority @ clear[i]), float(fractions[i].mean()), -candidates[i]["tile_cloud_pct"])

        best = max(remaining, key=rank)
        remaining.remove(best)
        if not clear[best].any():
            break
        scene = candidates[best]
        source = l1c_sources(scene)
        if source is None:
            continue
        counts += clear[best]
        selected.append({**scene, **source})
        day = date.fromisoformat(scene["datetime_utc"][:10])
        remaining = [i for i in remaining
                     if abs((date.fromisoformat(candidates[i]["datetime_utc"][:10]) - day).days) >= MIN_SEPARATION_DAYS]
    if len(selected) < 3:
        raise ValueError("Fewer than three locally usable matched L1C/SCL acquisitions")
    return selected


def select(design, cache, out):
    record = json.loads((design / "design_record.json").read_text())
    verify_hashes(design, record["files_sha256"])
    geo = study_grid(design, cache)
    cases = pd.read_csv(design / "critical_cells.csv").cell_id
    transform = Affine(*geo["transform"])
    dims = geo["cell_index"].shape
    positions = geo["cell_index"] >= 0
    index = geo["cell_index"][positions]
    signature = get_json(SCL_SAS)["token"]
    notes, chosen = [], []
    for year in YEARS:
        candidates = candidate_items(year)

        def coverage(scene):
            scl = get_scl(scene, cache, transform, dims, signature)
            counts = np.bincount(index[clear_mask(scl)[positions]], minlength=len(geo["names"]))
            return counts / geo["land_centers"]

        with ThreadPoolExecutor(max_workers=3) as pool:
            fractions = np.stack(list(pool.map(coverage, candidates)))
        selected = choose(candidates, fractions, geo["names"], cases)
        selected_ids = {a["scene_id"] for a in selected}
        for scene, quality in zip(candidates, fractions):
            notes.append({"year": year, **scene, "clear_share_mean": float(quality.mean()),
                          "cells_at_80pct_clear": int((quality >= 0.80).sum()),
                          "critical_at_80pct_clear": int((quality[np.isin(geo["names"], cases)] >= 0.80).sum()),
                          "selected": scene["scene_id"] in selected_ids})
        chosen.extend(selected)
        print(year, "candidate scenes", len(candidates), "selected", [s["datetime_utc"][:10] for s in selected], flush=True)
    out.mkdir(parents=True, exist_ok=True)
    write_csv(pd.DataFrame(notes), out / "sentinel_candidates.csv")
    payload = {"design_sha256": checksum(design / "design_record.json"), "scenes": chosen}
    (out / "sentinel_manifest.json").write_text(json.dumps(payload, indent=2) + "\n")


def fetch_scene(scene, cache, transform, dims):
    path = cache / "spectral" / f"{scene['scene_id']}.npz"
    if path.exists():
        return path
    signature = get_json(SCL_SAS)["token"] if not (cache / "scl" / f"{scene['scene_id']}.npz").exists() else None
    scl = get_scl(scene, cache, transform, dims, signature)
    valid = clear_mask(scl)
    dn = np.zeros((*dims, len(BANDS)), dtype=np.uint16)
    def band_values(band):
        url = scene["l1c_assets"][band]
        if not url.startswith(GCS_BASE) or urlsplit(url).query:
            raise ValueError("Expected an unsigned Sentinel-2 L1C URL")
        # Sequential HTTP reads are much faster than JP2 packet-range requests.
        with NamedTemporaryFile(suffix=".jp2") as tmp, open_https(url) as remote:
            shutil.copyfileobj(remote, tmp, 1024 * 1024)
            tmp.flush()
            size = remote.headers.get("Content-Length")
            if size is not None and tmp.tell() != int(size):
                raise ValueError(f"Incomplete L1C download for {scene['scene_id']} {band}")
            values = read_window(tmp.name, transform, dims, Resampling.bilinear, "float32")
        return values

    with ThreadPoolExecutor(max_workers=2) as pool:
        for i, values in enumerate(pool.map(band_values, BANDS)):
            valid &= values > 0
            dn[:, :, i] = np.rint(values + scene["offsets"][BANDS[i]]).clip(0, 65535).astype(np.uint16)
    with atomic_target(path) as tmp, tmp.open("wb") as f:
        np.savez(f, dn=dn, clear=valid, scene_id=scene["scene_id"], product=scene["l1c_product"])
    print("Cached original L1C reflectance and L2A cloud mask", scene["datetime_utc"][:10], flush=True)
    return path


def fetch(design, cache, manifest, year):
    geo = study_grid(design, cache)
    record = json.loads(manifest.read_text())
    if record["design_sha256"] != checksum(design / "design_record.json"):
        raise ValueError("Sentinel-2 manifest differs from the locked geographic design")
    scenes = [s for s in record["scenes"] if year is None or s["datetime_utc"].startswith(str(year))]
    for scene in scenes:
        fetch_scene(scene, cache, Affine(*geo["transform"]), geo["cell_index"].shape)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("step", choices=("select", "fetch"))
    parser.add_argument("--design", type=Path, default=AUDIT_DESIGN)
    parser.add_argument("--cache", type=Path, default=CACHE)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--year", type=int, choices=YEARS)
    opts = parser.parse_args()
    if opts.step == "select":
        select(opts.design, opts.cache, opts.output)
    else:
        fetch(opts.design, opts.cache, opts.output / "sentinel_manifest.json", opts.year)
