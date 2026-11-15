"""Fetch study-area data from publisher endpoints and retain source identifiers."""
import argparse
import csv
import datetime as dt
import hashlib
import json
import os
from concurrent.futures import ThreadPoolExecutor
from contextlib import contextmanager
from pathlib import Path
from tempfile import NamedTemporaryFile
from urllib.error import HTTPError
from urllib.parse import urlencode, urlsplit
from urllib.request import Request, urlopen

import numpy as np
import rasterio
from rasterio.warp import transform_bounds
from rasterio.windows import Window, from_bounds

RAW = Path("data/raw")
YEAR = 2025
DATES = tuple(
    (dt.date(YEAR, 1, 1) + dt.timedelta(days=day)).isoformat() for day in range(365)
)
BBOX = (103.54, 1.11, 104.16, 1.54)
URA_URL = (
    "https://api-open.data.gov.sg/v1/public/api/datasets/"
    "d_2cc750190544007400b2cfd5d7f53209/poll-download"
)
ESRI_URL = (
    "https://lulctimeseries.blob.core.windows.net/lulctimeseriesv003/lc2025/"
    "48N_20250101-20251231.tif"
)
WORLD_FILE = "ESA_WorldCover_10m_2021_v200_N00E102_Map.tif"
WORLD_URL = "https://esa-worldcover.s3.eu-central-1.amazonaws.com/v200/2021/map/" + WORLD_FILE
CMR_URL = "https://cmr.earthdata.nasa.gov/search/granules.json"
TOKEN_URL = "https://planetarycomputer.microsoft.com/api/sas/v1/token/modiseuwest/modis-061-cogs"
MODIS_COG = "https://modiseuwest.blob.core.windows.net/modis-061-cogs/MOD11A1/28/08"
READ_ENV = {
    "GDAL_DISABLE_READDIR_ON_OPEN": "EMPTY_DIR",
    "GDAL_HTTP_TIMEOUT": "30",
    "GDAL_HTTP_CONNECTTIMEOUT": "10",
    "GDAL_HTTP_MAX_RETRY": "2",
    "GDAL_PAM_ENABLED": "NO",
}


def checksum(path):
    with Path(path).open("rb") as stream:
        digest = hashlib.sha256()
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


@contextmanager
def atomic_target(destination):
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    with NamedTemporaryFile(
        dir=destination.parent, prefix=destination.stem + ".", suffix=".part", delete=False
    ) as stream:
        temp_name = Path(stream.name)
    try:
        yield temp_name
        os.replace(temp_name, destination)
    finally:
        temp_name.unlink(missing_ok=True)


def open_https(url, method="GET"):
    if urlsplit(url).scheme != "https":
        raise ValueError("The published endpoint must use HTTPS")
    return urlopen(
        Request(url, method=method, headers={"User-Agent": "Singapore-LST-2025 research"}),
        timeout=30,
    )


def get_json(url):
    with open_https(url) as reply:
        return json.load(reply)


def spatial_window(dataset):
    bounds = transform_bounds("EPSG:4326", dataset.crs, *BBOX, densify_pts=41)
    rough = from_bounds(*bounds, transform=dataset.transform)
    col, row = np.floor([rough.col_off, rough.row_off]).astype(int)
    right, bottom = np.ceil(
        [rough.col_off + rough.width, rough.row_off + rough.height]
    ).astype(int)
    if not (0 <= col < right <= dataset.width and 0 <= row < bottom <= dataset.height):
        raise ValueError("The source tile does not enclose the Singapore study window")
    return Window(int(col), int(row), int(right - col), int(bottom - row))


def fetch_esri_window(url, destination):
    if destination.is_file():
        return
    with rasterio.Env(**READ_ENV), rasterio.open("/vsicurl/" + url) as ds:
        if (
            ds.count != 1
            or ds.crs.to_epsg() != 32648
            or ds.dtypes[0] != "uint8"
            or ds.res != (10.0, 10.0)
        ):
            raise ValueError("Unexpected Esri 48N product format")
        window = spatial_window(ds)
        profile = {
            "driver": "GTiff", "width": int(window.width), "height": int(window.height),
            "count": 1, "dtype": ds.dtypes[0], "crs": ds.crs,
            "transform": ds.window_transform(window), "nodata": ds.nodata, "compress": "deflate",
        }
        with atomic_target(destination) as tmp, rasterio.open(tmp, "w", **profile) as out:
            out.write(ds.read(1, window=window), 1)
            out.write_colormap(1, ds.colormap(1))
            out.update_tags(source_url=url)


def fetch_spatial():
    boundary = RAW / "ura_planning_area_no_sea_2025.geojson"
    world = RAW / WORLD_FILE
    esri = RAW / "esri_impact_observatory_2025_48N_singapore.tif"
    if not boundary.is_file():
        signed_url = get_json(URA_URL)["data"]["url"]
        with open_https(signed_url) as reply, atomic_target(boundary) as tmp:
            contents = reply.read()
            if json.loads(contents)["type"] != "FeatureCollection":
                raise ValueError("URA did not return a GeoJSON FeatureCollection")
            tmp.write_bytes(contents)
    if not world.is_file():
        with open_https(WORLD_URL) as reply, atomic_target(world) as tmp:
            with tmp.open("wb") as out:
                for block in iter(lambda: reply.read(1024 * 1024), b""):
                    out.write(block)
            with rasterio.open(tmp) as ds:
                if ds.count != 1 or ds.crs.to_epsg() != 4326 or ds.dtypes[0] != "uint8":
                    raise ValueError("Unexpected WorldCover tile format")
    fetch_esri_window(ESRI_URL, esri)
    for path, url in [(boundary, URA_URL), (world, WORLD_URL), (esri, ESRI_URL)]:
        print(f"{path.name}  {checksum(path)}  {url}", flush=True)


def cmr_dates():
    query = urlencode({
        "collection_concept_id": "C1748058432-LPCLOUD",
        "temporal": "2025-01-01T00:00:00Z,2025-12-31T23:59:59Z",
        "producer_granule_id": "MOD11A1.A2025*.h28v08.*",
        "options[producer_granule_id][pattern]": "true",
        "page_size": 400,
    })
    entries = get_json(CMR_URL + "?" + query)["feed"]["entry"]
    days = {entry["time_start"][:10]: entry["title"] for entry in entries}
    if len(entries) != len(DATES) or set(days) != set(DATES):
        raise ValueError("CMR does not return exactly one tile per day in 2025")
    for day, granule in days.items():
        year_day = dt.date.fromisoformat(day).strftime("%Y%j")
        if not granule.startswith(f"MOD11A1.A{year_day}.h28v08.061."):
            raise ValueError(f"Unexpected product identifier for {day}")
    return sorted(days.items())


def read_modis_window(url, band):
    with rasterio.Env(**READ_ENV), rasterio.open("/vsicurl/" + url) as ds:
        expected_dtype, scale = {"LST_Day_1km": ("uint16", 0.02), "QC_Day": ("uint8", 1)}[band]
        if (
            ds.count != 1
            or ds.shape != (1200, 1200)
            or ds.dtypes[0] != expected_dtype
            or not np.isclose(ds.scales[0], scale)
            or ds.offsets[0] != 0
            or not np.allclose(ds.res, (926.625433, 926.625433), atol=0.001)
        ):
            raise ValueError(f"Unexpected MODIS format for {band}")
        window = spatial_window(ds)
        # Do not apply the raster mask to QC. The COG NoData tag is zero and zero is also best QA.
        values = ds.read(1, window=window, masked=False)
        grid = (ds.crs.to_wkt(), tuple(ds.window_transform(window))[:6], values.shape)
        return grid, values


def read_mirror_day(date, granule, public_signature):
    urls = [
        f"{MODIS_COG}/{granule.split('.')[1][1:]}/{granule}_{band}.tif"
        for band in ("LST_Day_1km", "QC_Day")
    ]
    available = []
    for url in urls:
        try:
            with open_https(url + "?" + public_signature, method="HEAD"):
                available.append(True)
        except HTTPError as error:
            if error.code != 404:
                raise RuntimeError(f"Public COG endpoint returned HTTP {error.code} for {date}") from None
            available.append(False)
        except OSError:
            raise RuntimeError(f"Could not query the public COG endpoint for {date}") from None
    if not all(available):
        return [date, granule, *urls, "missing_pair"], None
    try:
        first, dn = read_modis_window(urls[0] + "?" + public_signature, "LST_Day_1km")
        second, qc = read_modis_window(urls[1] + "?" + public_signature, "QC_Day")
    except rasterio.errors.RasterioError:
        raise RuntimeError(f"Could not read the COG window for {date}") from None
    if first != second:
        raise ValueError(f"The LST and QC grids differ for {date}")
    return [date, granule, *urls, "available"], (first, dn, qc)


def fetch_modis(workers):
    expected_days = cmr_dates()
    signature = get_json(TOKEN_URL)["token"]
    rows, observed_dates, dn, qc = [], [], [], []
    grid = None
    with ThreadPoolExecutor(max_workers=workers) as executor:
        fetches = executor.map(lambda item: read_mirror_day(*item, signature), expected_days)
        for record, observation in fetches:
            rows.append(record)
            if observation:
                actual_grid, raw_temperature, raw_qc = observation
                if grid is not None and grid != actual_grid:
                    raise ValueError(f"Inconsistent MODIS grid for {record[0]}")
                grid = actual_grid
                observed_dates.append(record[0])
                dn.append(raw_temperature)
                qc.append(raw_qc)
            if len(rows) % 31 == 0:
                print(f"Checked {len(rows)} of {len(expected_days)} dates", flush=True)
    if grid is None:
        raise RuntimeError("No published MODIS LST and QC pair could be read")
    with atomic_target(RAW / "modis_public_daily.npz") as tmp, tmp.open("wb") as out:
        np.savez_compressed(
            out, dates=observed_dates, dn=np.stack(dn), qc=np.stack(qc),
            crs=grid[0], transform=grid[1],
        )
    with atomic_target(RAW / "modis_public_daily_inventory.csv") as tmp, tmp.open("w", newline="") as out:
        writer = csv.writer(out, lineterminator="\n")
        writer.writerow(["date", "nasa_granule", "lst_url", "qc_url", "mirror_status"])
        writer.writerows(rows)
    print(f"Retrieved {len(observed_dates)} of {len(rows)} expected MODIS pairs", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("dataset", choices=("spatial", "modis"))
    parser.add_argument("--workers", type=int, choices=range(1, 9), default=6)
    args = parser.parse_args()
    if args.dataset == "spatial":
        fetch_spatial()
    else:
        fetch_modis(args.workers)


if __name__ == "__main__":
    main()
