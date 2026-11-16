"""Freeze the wider-only pair records and their original Landsat source assets."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
from urllib.parse import urlsplit

import pandas as pd

from analyze_landsat import write_csv
from fetch_inputs import checksum, get_json
from make_landsat_panel import SOURCE_FILES
from prepare_dense_tree_swap import DESIGN as STAGE6
from prepare_grass_tree_substitution import SIDES, VALIDATION
from prepare_local_comparison import BLOCK, KEY
from prepare_overpass_validation import verify_hashes


OUTPUT = Path("results/influential_pair_audit")
DESIGN = OUTPUT / "design"
SNAPSHOT = Path("data/snapshot/influential_pair_audit")
DISCOVERY = Path("data/selected_landsat_scenes.csv")
VALIDATION_MANIFEST = VALIDATION / "validation_scenes.csv"
IDENTITY = ["date", *SIDES]
OPTICAL = ("blue", "green", "red", "nir08", "qa_radsat")
CATALOG = "https://planetarycomputer.microsoft.com/api/stac/v1/collections/landsat-c2-l2/items"


def read_pair_inputs(stage6):
    record = json.loads((stage6 / "design_record.json").read_text())
    verify_hashes(stage6, record["design_sha256"])
    all_pairs = pd.read_csv(stage6 / "pairs.csv", float_precision="round_trip")
    dense = all_pairs[all_pairs.built_band.eq("ge80")]
    return {
        "dense": dense[dense.comparison.eq("5pp")],
        "wider_balance": dense[dense.comparison.eq("5pp_stage5_balance")],
    }


def focus_pairs(primary, wider):
    in_primary = pd.MultiIndex.from_frame(wider[IDENTITY]).isin(pd.MultiIndex.from_frame(primary[IDENTITY]))
    cases = wider.loc[~in_primary].copy()
    if cases.duplicated(IDENTITY).any():
        raise ValueError("Wider matching contains duplicate same-date endpoint IDs")
    counts = wider.groupby("date").size()
    cases["wider_date_weight"] = 1 / cases.date.map(counts) / len(counts)
    sites = cases.drop_duplicates(SIDES)[SIDES].copy()
    sites["site"] = [chr(65 + i) for i in range(len(sites))]
    cases = cases.merge(sites, on=SIDES, validate="many_to_one").sort_values(IDENTITY)
    if len(cases) != 9 or len(sites) != 3 or len(cases[BLOCK].drop_duplicates()) != 2:
        raise ValueError("The Stage 6 wider-only cohort changed")
    return cases


def date_manifest(cases):
    sources = pd.concat([pd.read_csv(DISCOVERY), pd.read_csv(VALIDATION_MANIFEST)], ignore_index=True)
    sources = sources.assign(date=sources.datetime_utc.str[:10])
    rows = sources[sources.date.isin(cases.date)].sort_values("date")
    if rows.date.duplicated().any() or set(rows.date) != set(cases.date):
        raise ValueError("Wider-only observations and overpass manifest disagree")
    eligible = pd.read_csv(VALIDATION / "eligible_cell_dates.csv")
    for endpoint in SIDES:
        keys = cases[["date", endpoint]].rename(columns={endpoint: "cell_id"})
        check = keys.merge(eligible[[*KEY, "scene_id"]], on=KEY, validate="one_to_one")
        if len(check) != len(keys) or not check.scene_id.eq(check.date.map(rows.set_index("date").scene_id)).all():
            raise ValueError("A focal cell-date differs from the registered original scene")
    return rows


def original_assets(row):
    item = get_json(f"{CATALOG}/{row['scene_id']}")
    if item["id"] != row["scene_id"] or item["properties"]["datetime"] != row["datetime_utc"]:
        raise ValueError("Catalog metadata differ from the frozen overpass")
    for key in [*SOURCE_FILES, *OPTICAL]:
        link = item["assets"][key]["href"]
        url = urlsplit(link)
        if url.scheme != "https" or url.netloc != "landsateuwest.blob.core.windows.net" or url.query or url.fragment:
            raise ValueError("Unexpected Landsat source address")
        if key in SOURCE_FILES and link != row[key]:
            raise ValueError(f"Frozen Landsat asset changed for {key}")
        row[key] = link
    return row


def prepare(stage6, output):
    populations = read_pair_inputs(stage6)
    cases = focus_pairs(populations["dense"], populations["wider_balance"])
    validation = json.loads((VALIDATION / "design_record.json").read_text())
    verify_hashes(VALIDATION, validation["design_sha256"])
    sources = date_manifest(cases)
    with ThreadPoolExecutor(max_workers=2) as pool:
        manifests = list(pool.map(original_assets, sources.to_dict("records")))
    # Within each repeated site, the reference date with greatest case coverage
    # uses the same calendar day for imagery when multiple sites coincide.
    shared = cases.groupby("date").size()
    review = (cases.assign(sites_on_date=cases.date.map(shared))
              .sort_values(["site", "sites_on_date", "date"], ascending=[True, False, True])
              .drop_duplicates("site")[["site", "date", *SIDES]])
    output.mkdir(parents=True, exist_ok=True)
    write_csv(cases, output / "pair_registry.csv")
    write_csv(review, output / "review_sites.csv")
    write_csv(pd.DataFrame(manifests), output / "landsat_sources.csv")
    record = {
        "stage6_design_sha256": checksum(stage6 / "design_record.json"),
        "stage3_design_sha256": checksum(VALIDATION / "design_record.json"),
        "omitted_blocks": cases[BLOCK].drop_duplicates().sort_values(BLOCK).values.tolist(),
        "wider_only_pair_observations": len(cases), "distinct_sites": int(cases.site.nunique()),
        "dates": sorted(cases.date.unique()),
        "files_sha256": {name: checksum(output / name) for name in (
            "pair_registry.csv", "review_sites.csv", "landsat_sources.csv",
        )},
    }
    (output / "design_record.json").write_text(json.dumps(record, indent=2) + "\n")
    print(json.dumps(record, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage6", type=Path, default=STAGE6)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    prepare(args.stage6, args.output)
