"""Draw the Stage 1 estimates and within-area residual maps."""
import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import numpy as np
import pandas as pd

from audit_coverage import BOUNDARY, read_boundary
from plot_landsat import GREEN, GRID, MUTED, SOURCES, draw_map, heading, style
from screen_landsat import GRID_CRS

PAIR = "#a6512a"


def error_bars(axis, entry, y, color, label=None):
    estimate = np.atleast_1d(entry["estimate_c_per_10pp"])
    lower = np.atleast_1d(entry["lower_c"])
    upper = np.atleast_1d(entry["upper_c"])
    axis.errorbar(
        estimate, np.atleast_1d(y), xerr=np.vstack([estimate - lower, upper - estimate]),
        color=color, fmt="o", capsize=3, lw=1.5, ms=5, label=label,
    )


def estimate_figure(results, output):
    summary = json.loads((results / "estimates.json").read_text())
    daily = pd.read_csv(results / "by_date.csv")
    fig, (left, right) = plt.subplots(
        1, 2, figsize=(11.6, 6), sharex=True, gridspec_kw={"width_ratios": [1, 1.25]},
    )
    fig.subplots_adjust(left=0.16, right=0.95, top=0.76, bottom=0.26, wspace=0.40)
    heading(fig, "Tree cover and daytime temperature within local areas",
            "Singapore · seven Landsat overpasses in 2025 · 500 m cells")
    for axis in (left, right):
        axis.axvline(0, color="#8d9a91", ls="--", lw=0.85)
        axis.grid(axis="x", lw=0.65, color=GRID, zorder=0)
        axis.set_axisbelow(True)
        axis.set_xlabel("LST difference per +10 pp tree cover (°C)")
        axis.set_xlim(min(-1.45, daily.lower_c.min() - 0.05), 0.075)
        axis.tick_params(axis="y", length=0)
    for y, method, color in (
        (0, "within_area", GREEN), (1, "nearby_pairs", PAIR), (2, "regression_on_matched_cells", "#647783"),
    ):
        error_bars(left, summary[method], y, color)
    left.set_yticks([0, 1, 2], ["Within-area\nregression", "Nearby pairs", "Regression on\npaired cells"])
    left.set_ylim(2.55, -0.55)
    left.set_title("Pooled estimates", loc="left", fontweight="bold", pad=11)
    day_labels = pd.to_datetime(daily.date.drop_duplicates()).dt.strftime("%d %b")
    positions = np.arange(len(day_labels))
    for label, color, shift in [
        ("Within-area regression", GREEN, -0.12), ("Nearby pairs", PAIR, 0.12),
    ]:
        items = daily[daily.analysis == label]
        error_bars(right, items, positions + shift, color, label)
    right.set_yticks(positions, day_labels)
    right.set_ylim(len(positions) - 0.55, -0.45)
    right.set_title("Estimates by overpass", loc="left", fontweight="bold", pad=11)
    handles, labels = right.get_legend_handles_labels()
    fig.legend(handles, labels, ncols=2, loc="lower center", bbox_to_anchor=(0.60, 0.10), frameon=False)
    fig.text(0.055, 0.057, "Bars show approximate 95% spatial-bootstrap intervals conditional on the observed dates.",
             fontsize=9, color=MUTED)
    fig.text(0.945, 0.022, SOURCES, ha="right", fontsize=8, color=MUTED)
    fig.savefig(output, dpi=180)
    plt.close(fig)


def residual_figure(results, output):
    boundary, _ = read_boundary(BOUNDARY, GRID_CRS)
    fitted = pd.read_csv(results / "fitted_cells.csv")
    diagnostics = pd.read_csv(results / "residual_diagnostics.csv").set_index("date")
    limit = float(np.quantile(np.abs(fitted.residual_c), 0.98))
    norm = TwoSlopeNorm(vmin=-limit, vcenter=0, vmax=limit)
    fig = plt.figure(figsize=(11.8, 8.7))
    panels = fig.add_gridspec(4, 3, left=0.052, right=0.96, top=0.84, bottom=0.11,
                              hspace=0.39, wspace=0.15, height_ratios=[1, 1, 1, 0.075])
    heading(fig, "Spatial structure after local adjustment",
            "Observed minus fitted daytime land surface temperature · 500 m cells")
    for i, (date, cells) in enumerate(fitted.groupby("date", sort=True)):
        axis = fig.add_subplot(panels[i // 3, i % 3])
        artist = draw_map(axis, boundary, cells, cells.residual_c, "RdBu_r", norm)
        day = pd.Timestamp(date).strftime("%d %b")
        statistic = diagnostics.loc[date, "residual_moran_i"]
        axis.set_title(f"{day}  ·  Moran’s I {statistic:.2f}", loc="left", fontsize=10.5,
                       fontweight="bold", pad=4)
    color_slot = fig.add_subplot(panels[3, :])
    bar = fig.colorbar(artist, cax=color_slot, orientation="horizontal", extend="both",
                       ticks=np.linspace(-np.floor(limit), np.floor(limit), 5))
    bar.outline.set_visible(False)
    bar.set_label("LST residual (°C)")
    fig.text(0.055, 0.039, "A common color scale saturates for the largest residuals.",
             fontsize=8.5, color=MUTED)
    fig.text(0.945, 0.018, SOURCES, ha="right", fontsize=8, color=MUTED)
    fig.savefig(output, dpi=180)
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    style()
    estimate_figure(args.results, args.output / "local_estimates.png")
    residual_figure(args.results, args.output / "local_residuals.png")
