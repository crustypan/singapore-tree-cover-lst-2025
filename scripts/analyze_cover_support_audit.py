"""Measure cell-level cover agreement and count support at the frozen references."""
import argparse
from itertools import combinations
import json
from pathlib import Path
import warnings

import numpy as np
import pandas as pd
import rasterio
from rasterio.features import shapes
from rasterio.transform import Affine
from rasterio.warp import Resampling, reproject
from shapely.geometry import shape

from analyze_landsat import write_csv
from audit_heldout_overpasses import AUDIT, land_cover
from fetch_cover_audit_sentinel import CACHE, study_grid
from fetch_inputs import RAW, checksum
from make_landsat_panel import VALID_ESRI
from prepare_built_support import STRATUM
from prepare_cover_support_audit import AUDIT_DESIGN, count_blocks, esri_rows, support_inputs
from prepare_overpass_validation import YEARS, verify_hashes
from screen_landsat import GRID_CRS


SNAPSHOT = Path("data/snapshot/cover_support_audit")
CONFIDENCE = (0.0, 0.50, 0.60, 0.70)
ESRI_FIELDS = {"tree": 2, "built": 7, "water": 1}
DW_FIELDS = {"tree": 1, "built": 6, "water": 0}
ESRI_ACCEPTED = (*VALID_ESRI, ESRI_FIELDS["water"])


def source_esri(year, geo):
    path, _ = land_cover(year)
    record = json.loads((AUDIT / str(year) / "audit_record.json").read_text())
    verify_hashes(RAW, {path.name: record["spatial_sha256"][path.name]})
    result = np.zeros(geo["cell_index"].shape, dtype=np.uint8)
    with rasterio.open(path) as ds:
        reproject(rasterio.band(ds, 1), result, dst_transform=Affine(*geo["transform"]),
                  dst_crs=GRID_CRS, dst_nodata=0, resampling=Resampling.nearest)
    return result


def pixel_cover(maps, keys, names, index, land, year, confidence, footprint, accepted, dry):
    n = len(names)
    confident = np.bincount(index[accepted], minlength=n)
    nonwater = np.bincount(index[dry], minlength=n)
    frame = pd.DataFrame({
        "cell_id": names, "year": year, "confidence": confidence, "footprint": footprint,
        "land_centers": land, "accepted_centers": confident, "common_nonwater_centers": nonwater,
        "qualified": 5 * nonwater >= 4 * land,
    })
    for source, labels, classes in zip(("esri", "dw"), maps, keys):
        for var, code in classes.items():
            mask, denom = (accepted, confident) if var == "water" else (dry, nonwater)
            values = np.bincount(index[mask & (labels == code)], minlength=n)
            frame[f"{source}_{var}_share"] = np.divide(values, denom, out=np.full(n, np.nan), where=denom > 0)
    return frame


def extraction(predictions, design, imagery, cache, out):
    frozen = json.loads((design / "design_record.json").read_text())
    verify_hashes(design, frozen["files_sha256"])
    manifest_sha = checksum(imagery / "sentinel_manifest.json")
    geo = study_grid(design, cache)
    cases = pd.read_csv(design / "critical_cells.csv")
    names, indexes = geo["names"], geo["cell_index"]
    lookup = {name: i for i, name in enumerate(names)}
    slices = []
    work = indexes >= 0
    for case in cases.itertuples():
        i = lookup[case.cell_id]
        y, x = geo["rows"][i], geo["cols"][i]
        chip = np.s_[y:y + 50, x:x + 50]
        slices.append(chip)
        work[chip] = True
    full_index = indexes[work]
    land = full_index >= 0
    index = full_index[land]
    fitted, esri, audit = {}, {}, {}
    chips = {"cases": cases.case.to_numpy(dtype="U"), "land": np.stack([indexes[chip] >= 0 for chip in slices])}
    for year in YEARS:
        pred_path = predictions / f"predictions_{year}.npz"
        info_path = predictions / f"predictions_{year}.json"
        meta = json.loads(info_path.read_text())
        if meta["prediction_sha256"] != checksum(pred_path) or meta["manifest_sha256"] != manifest_sha:
            raise ValueError(f"Inference inputs differ for {year}")
        data = np.load(pred_path, allow_pickle=False)
        np.testing.assert_array_equal(data["cell_ids"], names)
        np.testing.assert_array_equal(data["cell_index"], full_index)
        np.testing.assert_array_equal(data["cases"], chips["cases"])
        esri_raster = source_esri(year, geo)
        esri[year] = esri_raster[work][land]
        fitted[year] = (data["model_classes"][land], data["confidence"][land], data["clear_count"][land])
        for key in ("case_dn", "case_clear", "case_classes", "case_confidence", "case_count", "dates", "source_scene_ids"):
            chips[f"{key}_{year}"] = data[key]
        chips[f"esri_{year}"] = np.stack([esri_raster[chip] for chip in slices])
        audit[str(year)] = meta

    clear_all = np.logical_and.reduce([fitted[year][2] >= 2 for year in YEARS])
    esri_all = np.logical_and.reduce([np.isin(esri[year], ESRI_ACCEPTED) for year in YEARS])
    attrition = {"repeat_clear": clear_all, "esri_recognized": clear_all & esri_all}
    cell_frames = []
    for conf in CONFIDENCE:
        accepted, dry = {}, {}
        for year in YEARS:
            label, probability, count = fitted[year]
            accepted[year] = (count >= 2) & (probability >= conf) & np.isin(esri[year], ESRI_ACCEPTED)
            dry[year] = accepted[year] & (esri[year] != ESRI_FIELDS["water"]) & (label != DW_FIELDS["water"])
        all_accepted = np.logical_and.reduce(list(accepted.values()))
        all_dry = np.logical_and.reduce(list(dry.values()))
        label = "common_nonwater" if conf == 0 else f"confidence_{conf:.2f}"
        attrition[label] = all_dry
        for year in YEARS:
            pair = [esri[year], fitted[year][0]]
            for tag, valid, nonwater in [("annual", accepted[year], dry[year]), ("fixed", all_accepted, all_dry)]:
                cell_frames.append(pixel_cover(pair, [ESRI_FIELDS, DW_FIELDS], names, index, geo["land_centers"],
                                               year, conf, tag, valid, nonwater))

    out.mkdir(parents=True, exist_ok=True)
    write_csv(pd.concat(cell_frames, ignore_index=True), out / "paired_product_cover.csv")
    areas = {"cell_id": names, "land_centers": geo["land_centers"]}
    for name, mask in attrition.items():
        areas[name] = np.bincount(index[mask], minlength=len(names))
    write_csv(pd.DataFrame(areas), out / "optical_pixel_counts.csv")
    np.savez_compressed(out / "case_chips.npz", **chips)
    record = {"design_sha256": checksum(design / "design_record.json"), "manifest_sha256": manifest_sha,
              "inference": audit, "snapshot_sha256": {name: checksum(out / name) for name in (
                  "paired_product_cover.csv", "optical_pixel_counts.csv", "case_chips.npz")}}
    (out / "provenance.json").write_text(json.dumps(record, indent=2) + "\n")


def largest_patch(mask):
    if not mask.any():
        return 0
    return round(max(shape(item).area for item, _ in shapes(mask.astype(np.uint8), mask=mask, connectivity=8)))


def ndvi(data, year):
    dn = data[f"case_dn_{year}"].astype(np.float32)
    red, nir = dn[..., 2], dn[..., 3]
    keep = data[f"case_clear_{year}"] & ((red + nir) > 0)
    signal = np.divide(nir - red, nir + red, out=np.full_like(red, np.nan), where=keep)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        return {
            "n": keep.sum(axis=0), "min": np.nanmin(signal, axis=0),
            "max": np.nanmax(signal, axis=0), "median": np.nanmedian(signal, axis=0),
        }


def diagnostic_rows(chips):
    cases = chips["cases"]
    vegetation = {year: ndvi(chips, year) for year in YEARS}
    outputs = []
    for y1, y2 in combinations(YEARS, 2):
        a, b = vegetation[y1], vegetation[y2]
        for i, case in enumerate(cases):
            land = chips["land"][i]
            valid = land & (a["n"][i] >= 2) & (b["n"][i] >= 2)
            change = b["median"][i] - a["median"][i]
            gain = valid & (change >= 0.15) & (b["min"][i] > a["max"][i])
            loss = valid & (change <= -0.15) & (b["max"][i] < a["min"][i])
            item = {"case": case, "start_year": y1, "end_year": y2,
                    "comparable_land_pct": float(100 * valid.sum() / land.sum()),
                    "median_ndvi_delta": float(np.median(change[valid])) if valid.any() else None,
                    "vegetation_gain_largest_patch_100m2": largest_patch(gain),
                    "vegetation_loss_largest_patch_100m2": largest_patch(loss)}
            esri1, esri2 = chips[f"esri_{y1}"][i], chips[f"esri_{y2}"][i]
            dw1, dw2 = chips[f"case_classes_{y1}"][i], chips[f"case_classes_{y2}"][i]
            base = (valid & np.isin(esri1, ESRI_ACCEPTED) & np.isin(esri2, ESRI_ACCEPTED)
                    & (esri1 != 1) & (esri2 != 1) & (dw1 != 0) & (dw2 != 0))
            for minimum, suffix in [(0.60, ""), (0.0, "_no_score_cutoff")]:
                dry = (base & (chips[f"case_confidence_{y1}"][i] >= minimum)
                       & (chips[f"case_confidence_{y2}"][i] >= minimum))
                item["matched_nonwater_pct" + suffix] = float(100 * dry.sum() / land.sum())
                for variable in ("tree", "built"):
                    for tag, l1, l2, code in [
                        ("esri", esri1, esri2, ESRI_FIELDS[variable]),
                        ("dw", dw1, dw2, DW_FIELDS[variable]),
                    ]:
                        item[f"{tag}_{variable}_delta_pp" + suffix] = (
                            float(100 * (((l2 == code) & dry).sum() - ((l1 == code) & dry).sum()) / dry.sum())
                            if dry.any() else None
                        )
            outputs.append(item)
    return pd.DataFrame(outputs)


def transitions(chips):
    rows = []
    for y1, y2 in combinations(YEARS, 2):
        for i, case in enumerate(chips["cases"]):
            before, after = chips[f"esri_{y1}"][i], chips[f"esri_{y2}"][i]
            use = chips["land"][i] & (before != after) & np.isin(before, ESRI_ACCEPTED) & np.isin(after, ESRI_ACCEPTED)
            labels, counts = np.unique(np.stack([before[use], after[use]], axis=1), axis=0, return_counts=True)
            for (label1, label2), n in zip(labels, counts):
                rows.append({"case": case, "start_year": y1, "end_year": y2, "esri_class_from": int(label1),
                             "esri_class_to": int(label2), "centers": int(n)})
    return pd.DataFrame(rows)


def new_rows(original, covars, source):
    names = {f"{source}_{field}_share": f"{field}_fraction" for field in ("tree", "built")}
    return original.drop(columns=["tree_fraction", "built_fraction"]).merge(
        covars[["cell_id", "year", *names]].rename(columns=names),
        on=["cell_id", "year"], validate="many_to_one",
    )


def count_analysis(snapshot, design, out):
    frozen = json.loads((design / "design_record.json").read_text())
    verify_hashes(design, frozen["files_sha256"])
    record = json.loads((snapshot / "provenance.json").read_text())
    verify_hashes(snapshot, record["snapshot_sha256"])
    if record["design_sha256"] != checksum(design / "design_record.json"):
        raise ValueError("The optical snapshot does not match the frozen design")
    geometry, events, target, limit = support_inputs()
    original = esri_rows(geometry, events)
    data = pd.read_csv(snapshot / "paired_product_cover.csv", float_precision="round_trip")
    results, block_dates, agreement = [], [], []
    for (footprint, confidence), table in data.groupby(["footprint", "confidence"], sort=True):
        keep_ids = set(table.groupby("cell_id").qualified.all().loc[lambda ok: ok].index)
        paired = table[table.cell_id.isin(keep_ids)]
        reference = original[original.cell_id.isin(keep_ids)]
        comp = {name: new_rows(reference, paired, name) for name in ("esri", "dw")}
        stable_esri_ids = set(geometry.loc[geometry.common_stable, "cell_id"])
        settings = [("esri_original", reference, stable_esri_ids)]
        for name, rows in comp.items():
            vals = [f"{name}_{c}_share" for c in ("tree", "built", "water")]
            extrema = paired.groupby("cell_id")[vals]
            span = 100 * (extrema.max() - extrema.min())
            stable = set(span[span.le(limit + 1e-6).all(axis=1)].index)
            settings.append((name, rows, stable))
        joint = {}
        for label, frame, stable_ids in settings:
            for stable in (False, True):
                included = frame[frame.cell_id.isin(stable_ids)] if stable else frame
                values, sites = count_blocks(included, target, label)
                tags = {"confidence": confidence, "footprint": footprint, "stable": stable}
                results.append(values.assign(**tags))
                if not sites.empty:
                    block_dates.append(sites.assign(**tags))
                if label in ("esri", "dw"):
                    z = sites if not sites.empty else pd.DataFrame(columns=["year_group", *STRATUM, "lower",
                                                                            "higher", "shared_strict"])
                    joint[label, stable] = z
        for stable in (False, True):
            keys = ["year_group", *STRATUM]
            combined = joint["esri", stable].merge(joint["dw", stable], on=keys, suffixes=("_esri", "_dw"))
            for part in [*map(str, YEARS), "all"]:
                rows = combined[combined.year_group == part]
                supported = rows[["lower_esri", "higher_esri", "lower_dw", "higher_dw"]].all(axis=1)
                strict = rows[["shared_strict_esri", "shared_strict_dw"]].all(axis=1)
                agreement.append({
                    "year_group": part, "confidence": confidence, "footprint": footprint, "stable": stable,
                    "shared_blocks_both_products": len(rows.loc[supported, ["block_e", "block_n"]].drop_duplicates()),
                    "shared_block_dates_both_products": int(supported.sum()),
                    "strict_blocks_both_products": len(rows.loc[strict, ["block_e", "block_n"]].drop_duplicates()),
                    "strict_block_dates_both_products": int(strict.sum()),
                })
    out.mkdir(parents=True, exist_ok=True)
    write_csv(pd.concat(results, ignore_index=True), out / "matched_support_counts.csv")
    write_csv(pd.concat(block_dates, ignore_index=True), out / "matched_block_dates.csv")
    write_csv(pd.DataFrame(agreement), out / "support_agreement.csv")
    areas = pd.read_csv(snapshot / "optical_pixel_counts.csv")
    retained = []
    for col in areas.columns[2:]:
        ids = set(areas.loc[5 * areas[col] >= 4 * areas.land_centers, "cell_id"])
        totals, _ = count_blocks(original[original.cell_id.isin(ids)], target, col)
        retained.append(totals)
    write_csv(pd.concat(retained, ignore_index=True), out / "optical_filter_counts.csv")
    with np.load(snapshot / "case_chips.npz", allow_pickle=False) as chips:
        write_csv(diagnostic_rows(chips), out / "image_diagnostics.csv")
        write_csv(transitions(chips), out / "esri_label_transitions.csv")
    print(pd.concat(results).query("year_group == 'all' and confidence == 0.60")[
        ["analysis", "footprint", "stable", "cells", "lower_endpoint_blocks", "higher_endpoint_blocks", "shared_blocks"]
    ].to_string(index=False))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("step", choices=("extract", "count"))
    parser.add_argument("--design", type=Path, default=AUDIT_DESIGN)
    parser.add_argument("--imagery", type=Path, default=Path("results/cover_support_audit/imagery"))
    parser.add_argument("--predictions", type=Path, default=Path(".scratch/cover-support/inference"))
    parser.add_argument("--cache", type=Path, default=CACHE)
    parser.add_argument("--snapshot", type=Path, default=SNAPSHOT)
    parser.add_argument("--output", type=Path)
    opts = parser.parse_args()
    if opts.step == "extract":
        extraction(opts.predictions, opts.design, opts.imagery, opts.cache, opts.snapshot)
    else:
        if opts.output is None:
            parser.error("count requires --output")
        count_analysis(opts.snapshot, opts.design, opts.output)
