"""Run the published Dynamic World forward model on harmonized Sentinel-2 L1C windows."""
# The normalization constants and inference sequence are from Google's Dynamic
# World Model Runner, commit 0dd288fae5fb5291d988e5ea3c12eb91df481ec5.
# Copyright 2022 Google LLC. See third_party/dynamic_world/LICENSE and NOTICE.
import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
import tensorflow as tf


NORMALIZATION = np.array([
    [1.7417268007636313, 2.023298706048351],
    [1.7261204997060209, 2.038905204308012],
    [1.6798346251414997, 2.179592821212937],
    [1.7734969472909623, 2.2890068333026603],
    [2.289154079164943, 2.6171674549378166],
    [2.382939712192371, 2.773418590375327],
    [2.3828939530384052, 2.7578332604178284],
    [2.1952484264967844, 2.789092484314204],
    [1.554812948247501, 2.4140534947492487],
])
WEIGHTS_SHA256 = {
    "saved_model.pb": "e47c4eafa62656b688893a4d86f67f0ad9e08baaf5b5b1ce2f320208b08d4b18",
    "variables/variables.data-00000-of-00001": "c4a38299f2c0a741682cbe48810d55afdc832ff451e3d3708082751fe655fcb1",
    "variables/variables.index": "f27000c6a3d7dcd3e76e333331b108697712c734d5d5ef11b3f34571027c8436",
}
TILE = 512
HALO = 48


def sha256(path):
    with path.open("rb") as f:
        digest = hashlib.sha256()
        for block in iter(lambda: f.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def probabilities(model, image):
    resized = (np.log1p(image * 0.005) - NORMALIZATION[:, 0]) / NORMALIZATION[:, 1]
    scaled = 1.0 / (1.0 + np.exp(1.0 - 5.0 * resized))
    logits = model(tf.convert_to_tensor(scaled[None], dtype=tf.float32))
    return tf.nn.softmax(logits, axis=-1).numpy()[0]


def predict_scene(model, dn, clear, pixel, probability_sum, counts):
    height, width = pixel.shape
    for y in range(0, height, TILE):
        for x in range(0, width, TILE):
            y1, x1 = min(y + TILE, height), min(x + TILE, width)
            indexes = pixel[y:y1, x:x1]
            valid = (indexes >= 0) & clear[y:y1, x:x1]
            if not valid.any():
                continue
            top, left = max(0, y - HALO), max(0, x - HALO)
            block = dn[top:min(height, y1 + HALO), left:min(width, x1 + HALO)]
            prob = probabilities(model, block)[y - top:y1 - top, x - left:x1 - left]
            ids = indexes[valid]
            probability_sum[ids] += prob[valid]
            counts[ids] += 1


def run(year, cache, design, imagery, weights, output):
    for name, expected in WEIGHTS_SHA256.items():
        if sha256(weights / name) != expected:
            raise ValueError(f"Dynamic World weight checksum differs for {name}")
    geo = np.load(cache / "geometry.npz", allow_pickle=False)
    feed = json.loads((imagery / "sentinel_manifest.json").read_text())
    design_sha = sha256(design / "design_record.json")
    if feed["design_sha256"] != design_sha or str(geo["design_sha256"]) != design_sha:
        raise ValueError("Imagery and geometry differ from the frozen design")
    scenes = sorted([s for s in feed["scenes"] if s["datetime_utc"].startswith(str(year))],
                    key=lambda s: s["datetime_utc"])
    if len(scenes) < 3:
        raise ValueError(f"Insufficient acquisitions for {year}")
    with (design / "critical_cells.csv").open(newline="") as f:
        cases = list(csv.DictReader(f))
    names = geo["names"]
    idx = {name: i for i, name in enumerate(names)}
    slices = []
    index = geo["cell_index"]
    work = index >= 0
    for row in cases:
        cell = idx[row["cell_id"]]
        y, x = geo["rows"][cell], geo["cols"][cell]
        chip = np.s_[y:y + 50, x:x + 50]
        slices.append(chip)
        work[chip] = True
    n = int(work.sum())
    pixel = np.full(index.shape, -1, dtype=np.int32)
    pixel[work] = np.arange(n, dtype=np.int32)
    total = np.zeros((n, 9), dtype=np.float32)
    counts = np.zeros(n, dtype=np.uint8)
    case_dn = []
    case_clear = []
    source_hashes = {}
    model = tf.saved_model.load(str(weights))
    for scene in scenes:
        scene_id = scene["scene_id"]
        source = cache / "spectral" / f"{scene_id}.npz"
        source_hashes[scene_id] = sha256(source)
        with np.load(source, allow_pickle=False) as data:
            if str(data["scene_id"]) != scene_id or str(data["product"]) != scene["l1c_product"]:
                raise ValueError(f"Cached product differs for {scene_id}")
            dn, clear = data["dn"], data["clear"]
        predict_scene(model, dn, clear, pixel, total, counts)
        case_dn.append(np.stack([dn[chip][:, :, [0, 1, 2, 6]] for chip in slices]))
        case_clear.append(np.stack([clear[chip] for chip in slices]))
        print("Inferred Dynamic World", scene["datetime_utc"][:10], flush=True)

    total /= np.maximum(counts, 1)[:, None]
    labels = np.argmax(total, axis=1).astype(np.uint8)
    confidence = total.max(axis=1)
    case_pixels = np.stack([pixel[chip] for chip in slices])
    output.mkdir(parents=True, exist_ok=True)
    with (output / f"predictions_{year}.npz").open("wb") as f:
        np.savez_compressed(f, model_classes=labels, confidence=confidence, clear_count=counts,
                            cell_index=index[work], cell_ids=names,
                            cases=np.array([row["case"] for row in cases]), case_classes=labels[case_pixels],
                            case_confidence=confidence[case_pixels], case_count=counts[case_pixels],
                            case_dn=np.stack(case_dn), case_clear=np.stack(case_clear),
                            dates=np.array([s["datetime_utc"][:10] for s in scenes]),
                            source_scene_ids=np.array([s["scene_id"] for s in scenes]))
    meta = {
        "year": year, "tensorflow": tf.__version__, "halo_px": HALO, "tile_px": TILE,
        "model_sha256": WEIGHTS_SHA256, "manifest_sha256": sha256(imagery / "sentinel_manifest.json"),
        "input_sha256": source_hashes, "prediction_sha256": sha256(output / f"predictions_{year}.npz"),
    }
    (output / f"predictions_{year}.json").write_text(json.dumps(meta, indent=2) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--year", type=int, choices=(2023, 2024, 2025), required=True)
    parser.add_argument("--cache", type=Path, default=Path("data/raw/cover_support_audit"))
    parser.add_argument("--design", type=Path, default=Path("results/cover_support_audit/design"))
    parser.add_argument("--imagery", type=Path, default=Path("results/cover_support_audit/imagery"))
    parser.add_argument("--weights", type=Path, default=Path(".scratch/cover-support/model/model/forward"))
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    run(args.year, args.cache, args.design, args.imagery, args.weights, args.output)
