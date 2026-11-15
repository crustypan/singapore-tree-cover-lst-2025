"""Aggregate the frozen validation overpasses to the original 500 m grid."""
import argparse
import gzip
import json
from pathlib import Path
import shutil
from tempfile import TemporaryDirectory

import pandas as pd
from rasterio.transform import array_bounds

from analyze_landsat import load_panel, observations
from audit_heldout_overpasses import AUDIT, land_cover, read_candidates
from fetch_inputs import RAW, atomic_target, checksum
from make_landsat_panel import aggregate, fetch_windows
from prepare_local_comparison import KEY
from prepare_overpass_validation import DESIGN, verify_hashes
from screen_landsat import land_layers


SNAPSHOTS = Path("data/snapshot/overpass_validation")


def make_panel(design, audit, output):
    record = json.loads((design / "design_record.json").read_text())
    verify_hashes(design, record["design_sha256"])
    manifest = pd.read_csv(design / "validation_scenes.csv")
    if manifest.scene_id.tolist() != record["scene_ids"] or manifest.datetime_utc.str[:10].tolist() != record["dates"]:
        raise ValueError("Scene list differs from the validation design")
    registry = pd.read_csv(design / "eligible_cell_dates.csv")
    for year, rows in manifest.groupby(manifest.datetime_utc.str[:4].astype(int), sort=True):
        manifest_path = audit / str(year) / "eligible_new_scenes.csv"
        audit_path = manifest_path.with_name("audit_record.json")
        if checksum(audit_path) != record["audit_sha256"][str(year)]:
            raise ValueError("The candidate coverage record has changed")
        audit_record = json.loads(audit_path.read_text())
        verify_hashes(manifest_path.parent, audit_record["result_sha256"])
        original = read_candidates(year, manifest_path)
        if not rows.reset_index(drop=True).equals(original):
            raise ValueError("Frozen validation manifest differs from the audited candidates")

        cover, _ = land_cover(year)
        verify_hashes(RAW, audit_record["spatial_sha256"])
        transform, world, esri, _, classes = land_layers(trim_to_30m=False, esri_path=cover)
        cache = RAW / "validation_windows" / "selected" / str(year)
        scenes = original.to_dict("records")
        dates = original.datetime_utc.str[:10].tolist()
        fetch_windows(scenes, array_bounds(*world.shape, transform), cache)
        dest = output / str(year)
        dest.mkdir(parents=True, exist_ok=True)
        with TemporaryDirectory(prefix="validation-panel-") as temporary:
            folder = Path(temporary)
            csv = folder / "world_500.csv"
            aggregate(csv, scenes, 500, 0, "WorldCover", (transform, world, esri, classes), cache)
            usable = observations(load_panel(folder, "world_500", dates))
            check_columns = [*KEY, "scene_id", "valid_lst_centers"]
            actual = usable[check_columns].sort_values(KEY).reset_index(drop=True)
            expected = registry.loc[registry.date.isin(dates), check_columns].sort_values(KEY).reset_index(drop=True)
            pd.testing.assert_frame_equal(actual, expected)
            compressed = dest / "world_500.csv.gz"
            with atomic_target(compressed) as tmp, tmp.open("wb") as f, csv.open("rb") as inp:
                with gzip.GzipFile(filename="", fileobj=f, mode="wb", mtime=0) as archive:
                    shutil.copyfileobj(inp, archive)
            provenance = {
                "dates": dates, "scene_ids": original.scene_id.tolist(), "esri_land_cover_year": year,
                "design_sha256": checksum(design / "design_record.json"),
                "inputs_sha256": audit_record["spatial_sha256"],
                "windows_sha256": {s["scene_id"]: checksum(cache / f"{s['scene_id']}.npz") for s in scenes},
                "panel_sha256": checksum(compressed), "decompressed_panel_sha256": checksum(csv),
            }
            (dest / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--design", type=Path, default=DESIGN)
    parser.add_argument("--audit", type=Path, default=AUDIT)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    make_panel(args.design, args.audit, args.output)
