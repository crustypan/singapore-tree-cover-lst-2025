"""Plot the Stage 2 temperature contrasts and model and support checks."""
import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from plot_built_support import HIGH
from plot_landsat import GREEN, GRID, INK, MUTED, SOURCES, heading, style


def draw_interval(axis, y, estimate, low, high, color):
    axis.hlines(y, low, high, color=color, lw=1.65, zorder=2)
    axis.vlines([low, high], y - 0.065, y + 0.065, color=color, lw=1.15, zorder=2)
    axis.scatter(estimate, y, color=color, s=29, zorder=3)


def contrast_figure(results, output):
    summary = json.loads((results / "estimates.json").read_text())
    checks = pd.read_csv(results / "robustness.csv").set_index("analysis")
    target = summary["contrast"]
    primary = summary["interaction"]
    fig, (left, right) = plt.subplots(1, 2, figsize=(12.6, 6.3),
                                     gridspec_kw={"width_ratios": [0.94, 1.18]})
    fig.subplots_adjust(left=0.123, right=0.956, top=0.76, bottom=0.255, wspace=0.81)
    heading(fig, "Tree cover and daytime surface temperature across built cover",
            f"Singapore · 500 m cells · seven overpasses in 2025 · tree cover "
            f"{target['tree_start_pct']:g}% to {target['tree_end_pct']:g}%")
    for axis in (left, right):
        axis.axvline(0, color="#83948c", lw=1, ls="--", zorder=1)
        axis.grid(axis="x", color=GRID, linewidth=0.7)
        axis.set_axisbelow(True)
        axis.tick_params(axis="y", length=0, pad=8)
        axis.spines["left"].set_visible(False)

    lower, higher = (target[f"{band}_built_pct"] for band in ("lower", "higher"))
    left.set_title("Primary model at each built reference", loc="left", fontweight="bold", pad=10)
    left.set_yticks([1, 0], [f"{p:g}% built\n{primary[f'{band}_band_blocks']} blocks"
                            for p, band in [(lower, "lower"), (higher, "higher")]])
    left.set_ylim(-0.98, 1.98)
    left.set_xlabel("Daytime LST difference (°C)")
    for band, y, color in [("lower", 1, GREEN), ("higher", 0, HIGH)]:
        entry = primary["contrasts"][f"{band}_built"]
        draw_interval(left, y, entry["estimate_c"], entry["lower_c"], entry["upper_c"], color)
        left.annotate(f"{entry['estimate_c']:+.2f} °C".replace("-", "−"), (entry["estimate_c"], y),
                      xytext=(0, 10), textcoords="offset points", ha="center", color=color, fontsize=10)
    min_left = min(primary["contrasts"][f"{band}_built"]["lower_c"] for band in ("lower", "higher"))
    left.set_xlim(min_left - 0.22, 0.18)

    variants = [
        ("Primary quadratic interaction", "Primary interaction"),
        ("Linear product", "Linear product"),
        ("Piecewise linear interaction", "Piecewise linear"),
        ("Shared block-dates", "Shared block-dates"),
        ("Shared block-dates with balanced bands", "Shared, balanced bands"),
        ("Narrower support windows", "Narrower windows"),
    ]
    fields = [f"higher_minus_lower_{suffix}_c" for suffix in ("estimate", "lower", "upper")]
    spec = checks.loc[[name for name, _ in variants], fields].to_numpy()
    y = np.arange(len(variants) + 1)
    for i, (estimate, lo, hi) in enumerate(spec):
        draw_interval(right, i, estimate, lo, hi, GREEN if i == 0 else MUTED)
    paired = summary["paired"]["same_block_dates"]
    pair_delta = paired.get("higher_minus_lower_c")
    positions = y
    labels = [label for _, label in variants]
    if pair_delta is not None:
        right.scatter(pair_delta, y[-1], marker="D", facecolors="none",
                      edgecolors=HIGH, linewidths=1.35, s=36, zorder=3)
        labels.append("Nearby pairs, shared area")
    else:
        positions = y[:-1]
    right.set_yticks(positions, labels)
    right.set_ylim(positions[-1] + 0.65, -0.65)
    extent = np.r_[spec.ravel(), 0, pair_delta if pair_delta is not None else 0]
    right.set_xlim(extent.min() - 0.17, extent.max() + 0.16)
    right.set_title("Difference between built references", loc="left", fontweight="bold", pad=10)
    right.set_xlabel(f"Difference in tree–LST contrast (°C)\n{higher:g}% minus {lower:g}% built")
    fig.text(0.055, 0.081,
             "Bars show approximate 95% spatial-bootstrap intervals conditional on the observed dates and support.",
             color=MUTED, fontsize=9)
    if pair_delta is not None:
        blocks = paired["spatial_blocks"]
        unit = "block" if blocks == 1 else "blocks"
        fig.text(0.055, 0.049,
                 f"The diamond is a descriptive nearby-pair difference from {blocks} shared 4 km {unit}.",
                 color=MUTED, fontsize=9)
    fig.text(0.945, 0.022, SOURCES, ha="right", fontsize=8, color=MUTED)
    fig.savefig(output, dpi=180)
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    style()
    contrast_figure(args.results, args.output / "contrast_estimates.png")
