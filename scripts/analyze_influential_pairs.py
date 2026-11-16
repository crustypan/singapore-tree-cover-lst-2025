"""Audit the original pixel QA and compare frozen fits with two blocks omitted."""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from rasterio.transform import Affine

from analyze_dense_tree_swap import MODELS, fitted_contrasts
from analyze_landsat import interval, write_csv
from analyze_overpass_validation import load_cohorts
from dense_swap_model import DRAWS, SEED, adjusted_fits, assess, design, full_rank, spatial_weights
from fetch_cover_audit_sentinel import clear_mask
from fetch_inputs import checksum
from make_influential_pair_chips import PADDING, SENTINEL_BANDS
from make_landsat_panel import VALID_ESRI, points_within_native
from make_validation_panel import SNAPSHOTS
from prepare_dense_tree_swap import RULES as STAGE6_RULES
from prepare_grass_tree_substitution import COVER, DESIGN as STAGE5, SIDES, VALIDATION
from prepare_influential_pair_audit import DESIGN, IDENTITY, SNAPSHOT, STAGE6, read_pair_inputs
from prepare_local_comparison import BLOCK, KEY
from prepare_overpass_validation import verify_hashes
from screen_landsat import mask_values


QUALITY = ("base", "qa4", "qa2")
MIN_COVERAGE = .60
SR_SCALE, SR_OFFSET = .0000275, -.2
BLOCKS = ("full", "omit_both")


def masks(thermal, pixel, cloud_distance, uncertainty):
    clear, plausible, celsius = mask_values(thermal, pixel)
    base = clear & plausible
    far = base & (cloud_distance >= 50) & (uncertainty >= 0)
    return (base, far & (uncertainty <= 400), far & (uncertainty <= 200)), celsius


def ndvi(red, nir, clear):
    keep = clear & (red >= 0) & (nir >= 0) & ((red + nir) > 0)
    if not keep.any():
        return float("nan"), 0
    return float(np.mean((nir[keep] - red[keep]) / (nir[keep] + red[keep]))), int(keep.sum())


def shares(classes, sample):
    classified = sample & np.isin(classes, VALID_ESRI)
    if not classified.any():
        raise ValueError("Requested coverage has no classified land")
    tree, rangeland, built = [float((classes[classified] == code).mean()) for code in (2, 11, 7)]
    return np.array([tree, rangeland, built, 1 - tree - rangeland - built])


def summarize_endpoint(landsat, geometry, side):
    land = geometry["land"][side]
    classes = geometry["classes"][side][land]
    t = Affine(*geometry["transform"][side])
    y, x = np.where(land)
    points = points_within_native(t.c + (x + .5) * t.a, t.f + (y + .5) * t.e,
                                  landsat["transform"], landsat["lwir11"].shape)
    arrays = {key: landsat[key][points] for key in landsat.files if key != "transform"}
    qa, celsius = masks(arrays["lwir11"], arrays["qa_pixel"], arrays["cdist"], arrays["qa"])
    out = {"land_centers": len(classes), "planning_area": str(geometry["planning_area"][side])}
    for level, use in zip(QUALITY, qa):
        out[f"{level}_valid_centers"] = int(use.sum())
        out[f"{level}_coverage"] = float(use.mean())
        out[f"{level}_lst_c"] = float(celsius[use].mean()) if use.any() else float("nan")
    for key, use in (("full", np.ones(len(classes), dtype=bool)), ("thermal", qa[0])):
        out.update({f"{key}_{name}": float(value) for name, value in zip(COVER[:4], shares(classes, use))})
    st_uncertainty = arrays["qa"][qa[0] & (arrays["qa"] >= 0)] * .01
    st_cloud = arrays["cdist"][qa[0] & (arrays["cdist"] >= 0)] * .01
    if len(st_uncertainty) == 0 or len(st_cloud) == 0:
        raise ValueError("Original accepted pixels lack uncertainty or cloud distance")
    qa_stats = {
        "qa_k_median": np.median(st_uncertainty), "qa_k_p90": np.quantile(st_uncertainty, .9),
        "qa_k_max": st_uncertainty.max(), "cloud_km_min": st_cloud.min(), "cloud_km_median": np.median(st_cloud),
        "cloud_km_p10": np.quantile(st_cloud, .1),
    }
    out.update({key: float(value) for key, value in qa_stats.items()})
    sr = [arrays[key] for key in ("blue", "green", "red", "nir08")]
    optical = ((arrays["qa_pixel"] & 63) == 0) & ((arrays["qa_radsat"] & 2078) == 0)
    optical &= (np.stack(sr) > 0).all(axis=0)
    red, nir = [values * SR_SCALE + SR_OFFSET for values in sr[2:]]
    out["landsat_ndvi"], count = ndvi(red, nir, optical)
    out["landsat_ndvi_coverage"] = count / len(classes)
    return out


def sentinel_summary(snapshot, sites):
    records = []
    for site in sites.itertuples(index=False):
        with np.load(snapshot / f"cover_{site.site}.npz", allow_pickle=False) as geo:
            land = geo["land"]
        with np.load(snapshot / f"sentinel_{site.site}.npz", allow_pickle=False) as s2:
            offsets = dict(zip(SENTINEL_BANDS, s2["offsets"]))
            for i, key in enumerate(("lower", "higher")):
                valid = clear_mask(s2["scl"][i])[PADDING:-PADDING, PADDING:-PADDING] & land[i]
                r, ir = [s2[name][i, PADDING:-PADDING, PADDING:-PADDING] for name in ("B04", "B08")]
                valid &= (r > 0) & (ir > 0)
                red = (r + offsets["B04"]) / s2["quantification"]
                nir = (ir + offsets["B08"]) / s2["quantification"]
                value, count = ndvi(red, nir, valid)
                records.append({"site": site.site, "endpoint": key, "sentinel_ndvi": value,
                                "sentinel_ndvi_coverage": count / land[i].sum()})
    return pd.DataFrame(records)


def audit(cases, snapshot, temperatures, cover):
    panels = temperatures.set_index(KEY)
    cover = cover.set_index(["year", "cell_id"])
    eps = 5.1e-9
    rows = []
    for case in cases.to_dict("records"):
        date, site = case["date"], case["site"]
        with np.load(snapshot / f"cover_{site}.npz", allow_pickle=False) as geo:
            with np.load(snapshot / f"landsat_{date}.npz", allow_pickle=False) as landsat:
                for i, end in enumerate(("lower", "higher")):
                    cell = case[f"{end}_tree_cell"]
                    result = summarize_endpoint(landsat, geo, i)
                    ref = panels.loc[(date, cell)]
                    annual = cover.loc[(int(date[:4]), cell)]
                    for level, original in zip(QUALITY, ("valid", "qa4", "strict")):
                        assert result[f"{level}_valid_centers"] == ref[f"{original}_lst_centers"]
                        np.testing.assert_allclose(result[f"{level}_coverage"], ref[f"{original}_lst_share"],
                                                   rtol=0, atol=eps)
                        outcome = "lst_c" if level == "base" else f"{original}_lst_c"
                        np.testing.assert_allclose(result[f"{level}_lst_c"], ref[outcome], rtol=0, atol=eps,
                                                   equal_nan=True)
                    np.testing.assert_allclose([result[f"full_{key}"] for key in COVER[:4]], annual[COVER[:4]].astype(float),
                                               rtol=0, atol=eps)
                    assert result["land_centers"] == annual.land_centers
                    rows.append({
                        "site": site, "date": date, "endpoint": end, "cell_id": cell,
                        "scene_id": ref.scene_id, "coast_km": annual.coast_km, "land_share": annual.land_share, **result,
                    })
    endpoints = pd.DataFrame(rows)
    pair_metrics = []
    for case in cases.to_dict("records"):
        local = endpoints.loc[endpoints.site.eq(case["site"]) & endpoints.date.eq(case["date"])].set_index("endpoint")
        lo, hi = local.loc["lower"], local.loc["higher"]
        net = case["tree_gap_pp"] + case["rangeland_gap_pp"]
        rule = STAGE6_RULES["5pp"]
        passes = {
            "rangeland_donor": lo.full_rangeland_fraction >= rule.target - 1e-10,
            "tree_gain": abs(case["tree_gap_pp"] - 100 * rule.target) <= 100 * rule.margin + 1e-8,
            "rangeland_loss": abs(case["rangeland_gap_pp"] + 100 * rule.target) <= 100 * rule.margin + 1e-8,
            "tree_plus_rangeland": abs(net) <= 100 * rule.combined + 1e-8,
            "built": abs(case["built_gap_pp"]) <= 100 * rule.built_other + 1e-8,
            "other": abs(case["other_gap_pp"]) <= 100 * rule.built_other + 1e-8,
            "land": abs(case["land_gap_pp"]) <= 100 * rule.land + 1e-8,
            "distance": case["distance_m"] <= rule.distance,
        }
        stages = {**case, "tree_plus_rangeland_gap_pp": net,
                  "failed_primary_calipers": "|".join(key for key, ok in passes.items() if not ok)}
        for level in QUALITY:
            coverage = local[f"{level}_coverage"].min()
            compared = coverage >= MIN_COVERAGE
            difference = hi[f"{level}_lst_c"] - lo[f"{level}_lst_c"] if compared else None
            stages.update({f"{level}_eligible": compared, f"{level}_min_coverage": coverage,
                           f"{level}_gap_c": difference,
                           f"{level}_per5_c": 5 * difference / case["tree_gap_pp"] if compared else None})
        stages.update({
            "thermal_tree_gap_pp": 100 * (hi.thermal_tree_fraction - lo.thermal_tree_fraction),
            "thermal_rangeland_gap_pp": 100 * (hi.thermal_rangeland_fraction - lo.thermal_rangeland_fraction),
            "thermal_built_gap_pp": 100 * (hi.thermal_built_fraction - lo.thermal_built_fraction),
            "qa_k_median_gap": hi.qa_k_median - lo.qa_k_median,
            "landsat_ndvi_gap": hi.landsat_ndvi - lo.landsat_ndvi,
        })
        pair_metrics.append(stages)
    return endpoints, pd.DataFrame(pair_metrics)


def sample_row(label, model, status, coordinates, fitted, case):
    res = fitted.get(model)
    ci = interval(res[1:]) if res is not None else (None, None)
    return {
        "case": case, "sample": label, "model": model, "status": "estimated" if res is not None else status["reason"],
        "pairs": status["pairs"], "spatial_pairs": status["spatial_pairs"], "dates": status["dates"],
        "blocks": len(coordinates), "adjustment_draws_retained": status.get("adjustment_draws_retained"),
        "rank_after_each_block_omission": status.get("rank_after_each_block_omission"),
        "estimate_c_per_5pp": res[0] if res is not None else None, "lower_c": ci[0], "upper_c": ci[1],
        "draws_retained": len(res) - 1 if res is not None else None,
    }


def paired_changes(full, retained, before, after, models):
    x, w, blocks, ids, _ = design(full)
    remaining = full.loc[retained].reset_index(drop=True)
    z, q, _, _, _ = design(remaining)
    y = full.lst_gap_c.to_numpy()
    normalized = y * 5 / full.tree_gap_pp.to_numpy()
    mult = spatial_weights(np.ones_like(w), ids, len(blocks))
    weights_before = w * mult
    weights_after = q * mult[:, retained]
    overlap = weights_after.sum(axis=1) > 0
    result = {}
    for model in models:
        if model == "matched":
            valid = overlap
            first = np.average(normalized, axis=0, weights=w)
            second = np.average(normalized[retained], weights=q)
            aa, bb = weights_before[valid], weights_after[valid]
            boot = (bb @ normalized[retained]) / bb.sum(axis=1) - (aa @ normalized) / aa.sum(axis=1)
        else:
            valid = overlap & full_rank(x, weights_before) & full_rank(z, weights_after)
            first = adjusted_fits(x, y, w[None, :])[0]
            second = adjusted_fits(z, y[retained], q[None, :])[0]
            boot = adjusted_fits(z, y[retained], weights_after[valid]) - adjusted_fits(x, y, weights_before[valid])
        if not boot.size:
            raise ValueError("No comparable spatial resamples before and after omission")
        np.testing.assert_allclose([first, second], [before[model][0], after[model][0]], atol=1e-10, rtol=0)
        lo, hi = interval(boot)
        result[model] = {"model": model, "change_c_per_5pp": second - first, "lower_c": lo,
                         "upper_c": hi, "paired_draws_retained": len(boot)}
    result["matched"]["omit_with_original_date_weights_c"] = np.average(normalized[retained], weights=w[retained])
    return list(result.values())


def omit_blocks(populations, stage6_results, omitted):
    archive = pd.read_csv(stage6_results / "paired_temperatures.csv", float_precision="round_trip")
    recorded = pd.read_csv(stage6_results / "estimates.csv").set_index(["case", "model"])
    rows, differences, retained_all = [], [], []
    for case, pairs in populations.items():
        temperatures = archive.loc[archive.case.eq(case), [*IDENTITY, "lst_gap_c"]]
        full = pairs.merge(temperatures, on=IDENTITY, validate="one_to_one")
        if len(full) != len(pairs):
            raise ValueError("The Stage 6 temperature ledger does not match the pair list")
        mask = ~pd.MultiIndex.from_frame(full[BLOCK]).isin(pd.MultiIndex.from_frame(omitted))
        kept = full.loc[mask].reset_index(drop=True)
        kept_ids = kept[["date", *SIDES]].copy()
        retained_all.append(kept_ids.assign(case=case))
        fits = {}
        for label, selection in zip(BLOCKS, (full, kept)):
            state = assess(selection)
            if not state["matched_ready"]:
                raise ValueError("The block omission leaves too little local support for the paired contrast")
            fitted, _, coordinates = fitted_contrasts(selection, state)
            fits[label] = fitted
            for model in MODELS:
                row = sample_row(label, model, state, coordinates, fitted, case)
                if label == "full" and row["status"] == "estimated":
                    metrics = ["estimate_c_per_5pp", "lower_c", "upper_c"]
                    np.testing.assert_allclose([row[k] for k in metrics],
                                               recorded.loc[(case, model), metrics].astype(float), atol=5e-9, rtol=0)
                rows.append(row)
        common = [model for model in MODELS if all(model in fits[label] for label in BLOCKS)]
        changes = paired_changes(full, mask, fits["full"], fits["omit_both"], common)
        differences.extend({"case": case, **change} for change in changes)
    return pd.DataFrame(rows), pd.DataFrame(differences), pd.concat(retained_all, ignore_index=True)


def analyze(pair_design, snapshot, output, audit_only=False):
    record = json.loads((pair_design / "design_record.json").read_text())
    verify_hashes(pair_design, record["files_sha256"])
    verify_hashes(STAGE6, {"design_record.json": record["stage6_design_sha256"]})
    verify_hashes(VALIDATION, {"design_record.json": record["stage3_design_sha256"]})
    meta = json.loads((snapshot / "provenance.json").read_text())
    verify_hashes(pair_design, {"design_record.json": meta["design_sha256"]})
    verify_hashes(snapshot, meta["file_sha256"])
    stage6_analysis = STAGE6.parent / "analysis"
    archived = json.loads((stage6_analysis / "analysis_record.json").read_text())
    verify_hashes(STAGE6, {"design_record.json": archived["dense_design_sha256"]})
    verify_hashes(stage6_analysis, archived["output_sha256"])
    val = json.loads((VALIDATION / "design_record.json").read_text())
    verify_hashes(VALIDATION, val["design_sha256"])
    stage5 = json.loads((STAGE5 / "design_record.json").read_text())
    verify_hashes(STAGE5, stage5["design_sha256"])
    cover = pd.read_csv(STAGE5 / "cell_covariates.csv", float_precision="round_trip")
    temperatures = pd.concat(load_cohorts(SNAPSHOTS, VALIDATION, val).values(), ignore_index=True)
    cases = pd.read_csv(pair_design / "pair_registry.csv", float_precision="round_trip")
    endpoints, qaqc = audit(cases, snapshot, temperatures, cover)
    sites = pd.read_csv(pair_design / "review_sites.csv")
    s2 = sentinel_summary(snapshot, sites)
    tables = {"endpoint_quality.csv": endpoints, "pair_quality.csv": qaqc, "sentinel_context.csv": s2}
    if not audit_only:
        estimates, changes, retained = omit_blocks(read_pair_inputs(STAGE6), stage6_analysis,
                                                   pd.DataFrame(record["omitted_blocks"], columns=BLOCK))
        tables.update({"block_omission.csv": estimates, "omission_changes.csv": changes,
                       "omission_retained_pairs.csv": retained})
    output.mkdir(parents=True, exist_ok=True)
    for name, table in tables.items():
        write_csv(table, output / name)
    log = {
        "design_sha256": checksum(pair_design / "design_record.json"),
        "imagery_sha256": checksum(snapshot / "provenance.json"), "seed": SEED, "draws": DRAWS,
        "output_sha256": {name: checksum(output / name) for name in tables},
    }
    (output / "analysis_record.json").write_text(json.dumps(log, indent=2) + "\n")
    print(qaqc[["site", "date", "base_gap_c", "base_per5_c", "qa4_per5_c", "qa2_per5_c", "qa2_min_coverage"]]
          .to_string(index=False))
    if not audit_only:
        print(estimates.to_string(index=False))
        print(changes.to_string(index=False))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--design", type=Path, default=DESIGN)
    parser.add_argument("--snapshot", type=Path, default=SNAPSHOT)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--audit-only", action="store_true")
    args = parser.parse_args()
    analyze(args.design, args.snapshot, args.output, args.audit_only)
