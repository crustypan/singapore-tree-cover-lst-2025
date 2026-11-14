"""Evaluate the fixed 2025 quality and spatial coverage rules on native MODIS cells."""
import argparse
import calendar
import csv
import json
import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import rasterio
from rasterio.crs import CRS
from rasterio.features import geometry_mask
from rasterio.transform import Affine, array_bounds, xy
from rasterio.warp import Resampling, reproject, transform_geom
from shapely.geometry import mapping, shape
from shapely.ops import unary_union
from shapely.validation import make_valid

from fetch_inputs import DATES, RAW, WORLD_FILE, checksum

MIN_DAYS = 24
MIN_MONTHS = 9
MIN_QUARTER_DAYS = 2
MIN_ANNUAL_AREA = 0.8
MIN_CLASSIFIED = 0.8
WORLD_CLASSES = (10, 20, 30, 40, 50, 60, 70, 80, 90, 95, 100)
ESRI_CLASSES = (0, 1, 2, 4, 5, 7, 8, 9, 10, 11)
ESRI_NAMES = dict(zip(ESRI_CLASSES, (
    "Unclassified", "Water", "Trees", "Flooded vegetation", "Crops", "Built area",
    "Bare ground", "Snow and ice", "Clouds", "Rangeland",
)))
BOUNDARY = RAW / "ura_planning_area_no_sea_2025.geojson"
ESRI = RAW / "esri_impact_observatory_2025_48N_singapore.tif"
MIRROR = RAW / "modis_public_daily.npz"
MIRROR_INVENTORY = RAW / "modis_public_daily_inventory.csv"
ATTRIBUTION = "NASA LP DAAC · URA/data.gov.sg · ESA WorldCover · Esri / Impact Observatory"


def read_boundary(path, crs):
    geojson = json.loads(path.read_text())
    names = {f["properties"]["PLN_AREA_C"] for f in geojson["features"]}
    if geojson["type"] != "FeatureCollection" or len(names) != 55:
        raise ValueError("The URA file must contain the 55 Master Plan 2025 planning areas")
    valid = [make_valid(shape(f["geometry"])) for f in geojson["features"]]
    dissolved = unary_union(valid)
    return shape(transform_geom("EPSG:4326", crs, mapping(dissolved))), len(names)


def esri_native(boundary_path):
    with rasterio.open(ESRI) as ds:
        boundary, count = read_boundary(boundary_path, ds.crs)
        land_cover = ds.read(1)
        in_ura = geometry_mask([boundary], out_shape=ds.shape, transform=ds.transform, invert=True)
        counts = np.bincount(land_cover[in_ura], minlength=12)
        color_table = ds.colormap(1)
        pixel_area = ds.res[0] * ds.res[1]
    if set(np.nonzero(counts)[0]) - set(ESRI_CLASSES):
        raise ValueError("The Esri clip contains a class outside its published legend")
    candidate = int(counts.sum() - counts[1])
    classified = int(candidate - counts[0] - counts[10])
    return {
        "planning_areas": count,
        "crs": "EPSG:32648",
        "pixel_area_m2": pixel_area,
        "boundary_pixels": int(counts.sum()),
        "non_water_pixels": candidate,
        "classified_non_water_pixels": classified,
        "classified_non_water_pct": 100 * classified / candidate,
        "classes": [
            {"class": code, "name": name, "pixels": int(counts[code]), "rgba": color_table[code]}
            for code, name in ESRI_NAMES.items()
        ],
    }


def read_daily_mirror():
    with np.load(MIRROR, allow_pickle=False) as data:
        dates = data["dates"].tolist()
        dn, qc = data["dn"], data["qc"]
        crs, transform = CRS.from_wkt(data["crs"].item()), Affine(*data["transform"])
    with MIRROR_INVENTORY.open(newline="") as stream:
        inventory = list(csv.DictReader(stream))
    if [line["date"] for line in inventory] != list(DATES):
        raise ValueError("Mirror inventory does not contain the original 2025 calendar")
    if [line["date"] for line in inventory if line["mirror_status"] == "available"] != dates:
        raise ValueError("Mirror windows and source inventory disagree")
    if dn.shape != qc.shape or dn.ndim != 3 or dn.shape[0] != len(dates):
        raise ValueError("Invalid mirror data cube")
    return dates, dn, qc, crs, transform


def screen_modis(dn, qc):
    # Signed values let the Earth Engine export distinguish NoData from a valid zero QA flag.
    qc = qc.astype(np.int32, copy=False)
    celsius = dn.astype(np.float64) * 0.02 - 273.15
    accepted = (
        (dn > 0) & (dn <= 65535) & (qc >= 0) & (qc <= 255)
        & ((qc & 3) == 0) & (((qc >> 2) & 3) == 0)
        & (((qc >> 4) & 3) <= 1) & (((qc >> 6) & 3) <= 1)
        & (celsius >= -10) & (celsius <= 70)
    )
    return accepted, np.where(accepted, celsius, np.nan)


def temporal_eligibility(month_days):
    if month_days.shape[0] != 12:
        raise ValueError("Expected 12 monthly counts")
    quarter_days = np.stack([m.sum(axis=0) for m in np.split(month_days, 4, axis=0)])
    days = month_days.sum(axis=0)
    months = np.count_nonzero(month_days, axis=0)
    passes = (days >= MIN_DAYS) & (months >= MIN_MONTHS) & (quarter_days >= MIN_QUARTER_DAYS).all(axis=0)
    return days, months, quarter_days, passes


def check_native_window_coverage(boundary, transform, height, width):
    west, south, east, north = array_bounds(height, width, transform)
    x0, y0, x1, y1 = boundary.bounds
    if x0 < west or y0 < south or x1 > east or y1 > north:
        raise ValueError("Daily raster does not cover the whole URA boundary")


def read_earth_engine(path, mirror):
    dates_ref, dn_ref, qc_ref, crs, transform = mirror
    with rasterio.open(path) as ds:
        if ds.count == 0 or ds.count % 3 != 0 or not 1 <= ds.count // 3 <= len(DATES):
            raise ValueError("Expected up to 365 date, LST DN, QC triplets")
        if set(ds.dtypes) != {"int32"}:
            raise ValueError("Expected int32 date, LST DN, and QC bands")
        if ds.crs is None or ds.crs != crs:
            raise ValueError("Earth Engine export CRS differs from the native sinusoidal grid")
        row_col_delta = ~ds.transform * (transform.c, transform.f)
        if (
            not np.allclose(row_col_delta, np.round(row_col_delta), atol=0.001, rtol=0)
            or not np.allclose(
                (ds.transform.a, ds.transform.b, ds.transform.d, ds.transform.e),
                (transform.a, transform.b, transform.d, transform.e), atol=1e-5, rtol=0,
            )
        ):
            raise ValueError("Earth Engine export pixel alignment differs from the native grid")
        boundary, _ = read_boundary(BOUNDARY, ds.crs)
        check_native_window_coverage(boundary, ds.transform, ds.height, ds.width)
        values = ds.read()
        file_transform, shape = ds.transform, ds.shape
    codes = values[0::3, 0, 0]
    if not np.array_equal(values[0::3], np.broadcast_to(codes[:, None, None], values[0::3].shape)):
        raise ValueError("Date bands must be constant across the export")
    expected = {int(d.replace("-", "")): d for d in DATES}
    if any(code not in expected for code in codes):
        raise ValueError("Date bands must encode calendar dates in 2025")
    dates = [expected[code] for code in codes]
    if sorted(set(dates)) != dates:
        raise ValueError("Date bands must be strictly increasing")
    dn, qc = values[1::3], values[2::3]
    if not (
        np.all(((dn >= 0) & (dn <= 65535)) | (dn == -9999))
        and np.all(((qc >= 0) & (qc <= 255)) | (qc == -9999))
    ):
        raise ValueError("The two data bands must contain raw MODIS DN and QC values")
    # Restrict the comparison to the common footprint, after both grids have been aligned.
    col, row = np.round(row_col_delta).astype(int)
    row_from, col_from = max(0, row), max(0, col)
    row_to, col_to = min(shape[0], row + dn_ref.shape[1]), min(shape[1], col + dn_ref.shape[2])
    if row_from >= row_to or col_from >= col_to:
        raise ValueError("The export has no spatial overlap with the verified MODIS source")
    refs = {date: i for i, date in enumerate(dates_ref)}
    expected_dn, expected_qc, actual_dn, actual_qc = [], [], [], []
    for i, date in enumerate(dates):
        if date in refs:
            j = refs[date]
            expected_dn.append(dn_ref[j, row_from - row:row_to - row, col_from - col:col_to - col])
            expected_qc.append(qc_ref[j, row_from - row:row_to - row, col_from - col:col_to - col])
            actual_dn.append(dn[i, row_from:row_to, col_from:col_to])
            actual_qc.append(qc[i, row_from:row_to, col_from:col_to])
    if not actual_dn:
        raise ValueError("The export and mirror have no common dates to compare")
    ref_pass, ref_temp = screen_modis(np.array(expected_dn), np.array(expected_qc))
    export_pass, export_temp = screen_modis(np.array(actual_dn), np.array(actual_qc))
    if not ref_pass.any():
        raise ValueError("No accepted overlapping observations could verify the Earth Engine export")
    if not np.array_equal(ref_pass, export_pass) or not np.allclose(ref_temp, export_temp, equal_nan=True):
        raise ValueError("Earth Engine data and mirror QA results differ on overlapping native pixels")
    return dates, dn, qc, crs, file_transform, {
        "overlap_days": len(actual_dn), "accepted_cell_days_checked": int(ref_pass.sum()),
    }


def nested_land_counts(crs, native, extent, boundary, subdivisions):
    rows, cols = extent
    totals = {name: np.zeros(extent, dtype=np.int32) for name in (
        "ura", "world_land", "world_water", "world_unknown", "esri_water", "esri_valid",
        "esri_tree", "esri_built", "esri_cloud", "esri_missing",
    )}
    with rasterio.open(RAW / WORLD_FILE) as world, rasterio.open(ESRI) as esri:
        for start in range(0, rows, 4):
            stop = min(rows, start + 4)
            stripe = (stop - start) * subdivisions, cols * subdivisions
            fine = native * Affine.translation(0, start) * Affine.scale(1 / subdivisions)
            in_ura = geometry_mask([boundary], out_shape=stripe, transform=fine, invert=True)
            if not in_ura.any():
                continue
            cover = {}
            for dataset, allowed in [(world, WORLD_CLASSES), (esri, ESRI_CLASSES)]:
                values = np.zeros(stripe, dtype=np.uint8)
                reproject(
                    rasterio.band(dataset, 1), values, dst_transform=fine, dst_crs=crs,
                    resampling=Resampling.nearest, dst_nodata=0,
                )
                if set(np.unique(values[in_ura])) - set(allowed) - {0}:
                    raise ValueError(f"Unexpected classes in {Path(dataset.name).name}")
                cover[dataset.name] = values
            wc, es = cover[world.name], cover[esri.name]
            masks = {
                "ura": in_ura,
                "world_land": in_ura & (wc != 0) & (wc != 80),
                "world_water": in_ura & (wc == 80),
                "world_unknown": in_ura & (wc == 0),
                "esri_water": in_ura & (es == 1),
                "esri_valid": in_ura & np.isin(es, (2, 4, 5, 7, 8, 9, 11)),
                "esri_tree": in_ura & (es == 2),
                "esri_built": in_ura & (es == 7),
                "esri_cloud": in_ura & (es == 10),
                "esri_missing": in_ura & (es == 0),
            }
            for name, mask in masks.items():
                totals[name][start:stop] = mask.reshape(
                    stop - start, subdivisions, cols, subdivisions
                ).sum(axis=(1, 3))
    if not np.array_equal(
        totals["ura"], totals["world_land"] + totals["world_water"] + totals["world_unknown"]
    ) or not np.array_equal(
        totals["ura"],
        totals["esri_valid"] + totals["esri_water"] + totals["esri_cloud"] + totals["esri_missing"],
    ):
        raise ValueError("Land-cover classes failed to account for the entire URA footprint")
    if totals["world_unknown"].sum() > 0:
        raise ValueError("WorldCover has uncovered/unknown areas within the URA study boundary")
    return totals


def write_rows(path, fields, rows):
    with path.open("w", newline="") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(fields)
        writer.writerows(rows)


def percentage(mask, weight):
    return float(100 * weight[mask].sum() / weight.sum())


def render_quality(path, boundary, transform, esri_cloud, world_land, counts, daily, source):
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 10, "axes.spines.top": False,
        "axes.spines.right": False, "text.color": "#162a32", "axes.labelcolor": "#162a32",
        "savefig.facecolor": "#fcfcfa", "figure.facecolor": "#fcfcfa",
    })
    fig = plt.figure(figsize=(11, 7.2), layout="constrained")
    gs = fig.add_gridspec(2, 2, height_ratios=(1.4, 1))
    ax, bx = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])
    cx = fig.add_subplot(gs[1, :])
    minx, miny, maxx, maxy = boundary.bounds
    west, south, east, north = array_bounds(*world_land.shape, transform)
    extent = np.array((west - minx, east - minx, south - miny, north - miny)) / 1000
    for axis, array, title, color in [
        (ax, np.where(world_land, counts, np.nan), "QA-passing dates per native cell", "viridis"),
        (bx, esri_cloud, "Esri cloud share of candidate land", "Blues"),
    ]:
        image = axis.imshow(array, cmap=color, extent=extent, interpolation="nearest", vmin=0)
        for polygon in getattr(boundary, "geoms", (boundary,)):
            if polygon.geom_type == "Polygon":
                x, y = np.asarray(polygon.exterior.coords).T
                axis.plot((x - minx) / 1000, (y - miny) / 1000, c="#52605b", lw=0.4)
        axis.set_xlim(-1.5, (maxx - minx) / 1000 + 1.5)
        axis.set_ylim(-1.5, (maxy - miny) / 1000 + 1.5)
        axis.set_xlabel("Relative easting (km)")
        axis.set_ylabel("Relative northing (km)")
        axis.set_title(title, loc="left", fontweight="bold")
        fig.colorbar(image, ax=axis, fraction=0.035, pad=0.02)
    missing = np.array([not r[1] for r in daily], dtype=int)
    breaks = np.diff(np.r_[0, missing, 0])
    for j, (start, stop) in enumerate(zip(np.where(breaks == 1)[0], np.where(breaks == -1)[0])):
        cx.axvspan(start - 0.5, stop - 0.5, color="#f4dcd0",
                   label="Source files unavailable" if j == 0 else None)
    cx.bar(np.arange(len(daily)), [r[3] if r[1] else 0 for r in daily],
           width=1, color="#167879", linewidth=0, label="WorldCover land with passing QA")
    cx.set_ylim(0, 100)
    cx.set_xlim(-0.5, len(daily) - 0.5)
    cx.set_xticks([j for j, row in enumerate(daily) if row[0].endswith("-01")],
                  list(calendar.month_abbr)[1:])
    cx.set_ylabel("WorldCover land (%)")
    cx.set_title("Day-by-day valid land area", loc="left", fontweight="bold")
    cx.grid(axis="y", alpha=0.2)
    cx.legend(loc="upper right", frameon=False, fontsize=8)
    fig.suptitle(f"Singapore daytime surface-temperature coverage\n2025 · {source}",
                 fontsize=17, x=0.02, ha="left", fontweight="bold")
    fig.supxlabel("Data  " + ATTRIBUTION, fontsize=8, x=0.98, ha="right")
    fig.savefig(path, dpi=180)
    plt.close(fig)


def audit(args):
    mirror = read_daily_mirror()
    comparison = None
    if args.earth_engine:
        dates, dn, qc, crs, transform, comparison = read_earth_engine(args.earth_engine, mirror)
        source = "Earth Engine export"
    else:
        dates, dn, qc, crs, transform = mirror
        source = "Microsoft public mirror"
    accepted, temperature = screen_modis(dn, qc)
    months_of_dates = np.array([int(date[5:7]) for date in dates])
    month_days = np.stack([accepted[months_of_dates == m].sum(axis=0) for m in range(1, 13)])
    valid_days, valid_months, quarters, eligible = temporal_eligibility(month_days)
    boundary, _ = read_boundary(BOUNDARY, crs)
    check_native_window_coverage(boundary, transform, *dn.shape[1:])
    cover = nested_land_counts(crs, transform, dn.shape[1:], boundary, args.subdivisions)
    world_land = cover["world_land"]
    esri_land = cover["ura"] - cover["esri_water"]
    native = esri_native(BOUNDARY)
    has_land = world_land > 0
    land_area = world_land.sum() * (transform.a / args.subdivisions) ** 2 / 1e6
    pixel_area = (transform.a / args.subdivisions) ** 2
    cover_good = (esri_land > 0) & (5 * cover["esri_valid"] >= 4 * esri_land)
    pct = percentage(eligible, world_land)
    missing_days = np.array([
        calendar.monthrange(2025, month)[1] - np.count_nonzero(months_of_dates == month)
        for month in range(1, 13)
    ])
    _, _, _, potentially_eligible = temporal_eligibility(
        month_days + missing_days[:, None, None]
    )
    upper_pct = percentage(potentially_eligible, world_land)
    if pct >= 100 * MIN_ANNUAL_AREA:
        status = "pass"
    elif upper_pct < 100 * MIN_ANNUAL_AREA:
        status = "fail"
    else:
        status = "incomplete_source"
    audit_stats = {
        "year": 2025, "temperature_source": source, "annual_decision": status,
        "input_sha256": {path.name: checksum(path) for path in (
            BOUNDARY, RAW / WORLD_FILE, ESRI, MIRROR, MIRROR_INVENTORY
        )},
        "calendar_dates": len(DATES),
        "available_input_dates": len(dates),
        "dates_with_accepted_worldcover_land_samples": int(accepted[:, has_land].any(axis=1).sum()),
        "native_cell_width_m": transform.a,
        "nested_subdivisions": args.subdivisions,
        "nested_spacing_m": transform.a / args.subdivisions,
        "worldcover_land_area_km2": float(land_area),
        "worldcover_land_cells": int(has_land.sum()),
        "max_qa_valid_days_on_worldcover_land": int(valid_days[has_land].max()),
        "max_qa_valid_months_on_worldcover_land": int(valid_months[has_land].max()),
        "land_with_any_qa_pct": percentage(valid_days > 0, world_land),
        "passes_24_days_alone_pct": percentage(valid_days >= MIN_DAYS, world_land),
        "passes_9_months_alone_pct": percentage(valid_months >= MIN_MONTHS, world_land),
        "passes_2_days_each_quarter_alone_pct": percentage((quarters >= MIN_QUARTER_DAYS).all(axis=0), world_land),
        "worldcover_land_eligible_pct": pct,
        "worldcover_land_eligible_upper_bound_pct": upper_pct,
        "esri_land_eligible_pct": percentage(eligible, esri_land),
        "joint_eligibility_worldcover_land_pct": percentage(eligible & cover_good, world_land),
        "esri_candidate_land_in_classified_cells_pct": percentage(cover_good, esri_land),
        "esri_candidate_land_classified_pct": float(100 * cover["esri_valid"].sum() / esri_land.sum()),
        "esri_native_grid": native,
    }
    if args.earth_engine:
        audit_stats["input_sha256"][args.earth_engine.name] = checksum(args.earth_engine)
        audit_stats["earth_engine_overlap_audit"] = comparison
    day_index = {date: index for index, date in enumerate(dates)}
    daily = []
    for date in DATES:
        index = day_index.get(date)
        if index is None:
            daily.append((date, 0, None, None, None))
        else:
            daily.append((
                date, 1, int(accepted[index, has_land].sum()),
                round(percentage(accepted[index], world_land), 6),
                round(percentage(accepted[index], esri_land), 6),
            ))
    monthly = []
    for month in range(1, 13):
        rows = [day for day in daily if int(day[0][5:7]) == month]
        monthly.append((
            f"2025-{month:02d}", len(rows), sum(r[1] for r in rows), sum(bool(r[2]) for r in rows),
            round(percentage(month_days[month - 1] > 0, world_land), 6),
        ))
    audit_stats["best_daily_worldcover_coverage"] = sorted(
        [r for r in daily if r[2]], key=lambda r: -r[3]
    )[:5]
    land_temps = temperature[:, has_land]
    if np.isfinite(land_temps).any():
        audit_stats["qa_passed_lst_c_percentile_1_50_99"] = np.round(
            np.nanpercentile(land_temps, [1, 50, 99]), 4,
        ).tolist()
    out = args.output
    out.mkdir(parents=True, exist_ok=True)
    (out / "coverage.json").write_text(json.dumps(audit_stats, indent=2, allow_nan=False) + "\n")
    write_rows(out / "daily_coverage.csv", (
        "date", "source_pair_available", "qa_passing_worldcover_cells",
        "worldcover_land_valid_pct", "esri_land_valid_pct",
    ), daily)
    write_rows(out / "monthly_coverage.csv", (
        "month", "calendar_days", "input_days", "days_with_any_accepted_land",
        "worldcover_land_with_one_or_more_observations_pct",
    ), monthly)
    records = []
    for r, c in zip(*np.nonzero(cover["ura"])):
        x, y = xy(transform, r, c)
        tree, built, valid = [cover[f"esri_{name}"][r, c] for name in ("tree", "built", "valid")]
        records.append((
            int(r), int(c), round(x, 3), round(y, 3), round(world_land[r, c] * pixel_area, 2),
            round(esri_land[r, c] * pixel_area, 2), round(valid * pixel_area, 2),
            round(tree / valid, 6) if valid else None, round(built / valid, 6) if valid else None,
            round(cover["esri_cloud"][r, c] * pixel_area, 2), int(cover_good[r, c]),
            int(valid_days[r, c]), int(valid_months[r, c]), *[int(q[r, c]) for q in quarters],
            *[int(m[r, c]) for m in month_days], int(eligible[r, c]),
        ))
    write_rows(out / "native_cell_coverage.csv", (
        "row", "col", "x_sinusoidal_m", "y_sinusoidal_m", "worldcover_land_m2",
        "esri_candidate_land_m2", "esri_classified_land_m2", "tree_fraction_of_classified",
        "built_fraction_of_classified", "esri_cloud_m2", "passes_esri_classification_rule",
        "qa_valid_days", "qa_valid_months", "q1_days", "q2_days", "q3_days", "q4_days",
        *[f"m{m:02d}_days" for m in range(1, 13)], "passes_temperature_rule",
    ), records)
    if not args.no_figure:
        esri_cloud_pct = np.divide(
            100.0 * cover["esri_cloud"], esri_land,
            out=np.full(esri_land.shape, np.nan), where=esri_land > 0,
        )
        render_quality(out / "data_quality.png", boundary, transform,
                       esri_cloud_pct, world_land, valid_days, daily, source)
    print(json.dumps({
        "annual_decision": status, "available_dates": len(dates), "max_valid_days": audit_stats[
            "max_qa_valid_days_on_worldcover_land"
        ], "max_valid_months": audit_stats["max_qa_valid_months_on_worldcover_land"],
        "worldcover_land_eligible_pct": pct,
        "esri_native_classified_land_pct": native["classified_non_water_pct"],
    }, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--earth-engine", type=Path, help="Native raw daily export from export_earth_engine.js")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--subdivisions", type=int, choices=(50, 100, 200), default=100)
    parser.add_argument("--no-figure", action="store_true")
    args = parser.parse_args()
    mpl_cache = Path(".scratch/cache/matplotlib")
    mpl_cache.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(mpl_cache.resolve()))
    with rasterio.Env(GDAL_PAM_ENABLED="NO"):
        audit(args)


if __name__ == "__main__":
    main()
