"""Plot the frozen comparison geography and the same-pair estimates."""
import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np
import pandas as pd

from audit_coverage import BOUNDARY, read_boundary
from grass_tree_water import MAP, read_water
from plot_landsat import GREEN, GRID, INK, MUTED, draw_map, heading, style
from screen_landsat import GRID_CRS


GRAY = "#dce1dd"
GOLD = "#cc9c44"
BLUE = "#237696"
WATER_FILL = "#9fc5d6"
SOURCES = "Sources  USGS Landsat · URA · ESA WorldCover · Esri / Impact Observatory · Singapore Land Authority"


def support_figure(design, map_path, output):
    days = pd.read_csv(design / "support_by_date.csv").query("design == 'primary'")
    membership = pd.read_csv(design / "cell_date_support.csv")
    cells = (membership.assign(candidate=membership.primary_links.gt(0))
             .groupby("cell_id")[["candidate", "primary_matched"]].any().reset_index())
    covariates = pd.read_csv(design / "cell_covariates.csv")
    locations = covariates[["cell_id", "easting_m", "northing_m"]].drop_duplicates()
    cells = cells.merge(locations, on="cell_id", validate="one_to_one")
    status = cells.candidate.astype(int) + cells.primary_matched.astype(int)
    boundary, _ = read_boundary(BOUNDARY, GRID_CRS)
    _, components, inventory = read_water(map_path)
    x0, y0 = 500 * np.floor(np.array(boundary.bounds[:2]) / 500) - 500
    cmap = ListedColormap([GRAY, GOLD, GREEN])
    fig = plt.figure(figsize=(11.8, 8.4))
    grid = fig.add_gridspec(2, 1, left=0.075, right=0.965, top=0.845, bottom=0.155,
                            hspace=0.34, height_ratios=[2, 1])
    heading(fig, "Where the tree-for-rangeland comparison is supported",
            "500 m cells in Singapore · 15 daytime Landsat overpasses, 2023–25")
    ax = fig.add_subplot(grid[0])
    draw_map(ax, boundary, cells, status, cmap, BoundaryNorm(np.arange(-0.5, 3), cmap.N))
    for shape, primary in zip(components, inventory.primary):
        if primary:
            x, y = np.asarray(shape.exterior.coords).T
            ax.fill((x - x0) / 1000, (y - y0) / 1000, color=WATER_FILL, lw=0, alpha=0.85, zorder=-0.5)
    left, bottom = ax.get_xlim()[0] + 1.5, ax.get_ylim()[0] + 1.6
    ax.plot([left, left + 5], [bottom, bottom], c=INK, lw=1.3)
    ax.text(left + 2.5, bottom + 0.7, "5 km", ha="center", fontsize=8)
    handles = [
        Patch(fc=GRAY, label="Usable, no candidate link"),
        Patch(fc=GOLD, label="Candidate, never selected"),
        Patch(fc=GREEN, label="Paired on ≥1 overpass"),
        Patch(fc=WATER_FILL, label="Mapped inland water ≥0.25 km²"),
    ]
    ax.legend(handles=handles, loc="upper center", ncols=4, frameon=False, fontsize=8.5,
              bbox_to_anchor=(0.50, -0.015), handlelength=1.15)
    ax = fig.add_subplot(grid[1])
    x = np.arange(len(days))
    ax.bar(x, days.candidate_cells, color="#ead7a6", width=0.68, label="Candidate cells")
    ax.bar(x, days.matched_cells, color=GREEN, width=0.42, label="Paired cells")
    labels = pd.to_datetime(days.date).dt.strftime("%d %b\n%Y")
    ax.set_xticks(x, labels, fontsize=8)
    ax.set_ylabel("Cells")
    ax.set_ylim(0, days.candidate_cells.max() * 1.34)
    ax.set_title("Candidate and paired cells per overpass", loc="left", fontweight="bold", pad=8)
    ax.grid(axis="y", color=GRID, alpha=0.85, linewidth=0.65)
    ax.set_axisbelow(True)
    ax.legend(loc="upper right", ncols=2, frameon=False, fontsize=8.6)
    fig.text(0.075, 0.055, "Land-cover and availability rules were applied before fitting surface temperature.",
             fontsize=8.8, color=MUTED)
    fig.text(0.965, 0.026, SOURCES, ha="right", fontsize=8, color=MUTED)
    fig.savefig(output, dpi=180)
    plt.close(fig)


def interval_artist(ax, height, center, lo, hi, color):
    ax.hlines(height, lo, hi, color=color, lw=1.75)
    ax.scatter(center, height, color=color, s=28, zorder=4, edgecolors="none")


def estimate_figure(analysis, output):
    estimates = pd.read_csv(analysis / "estimates.csv")
    labels = {
        "primary": "15 overpasses",
        "tighter": "Tighter pairs",
        "2025 dates": "2025 dates",
        "2023–24 dates": "2023–24 dates",
        "raw distance": "Untransformed distances",
        "0.10 km² water": "0.10 km² inland water cutoff",
    }
    fig, axes = plt.subplots(1, 2, figsize=(11.8, 6.1), sharey=True, gridspec_kw={"width_ratios": [1.25, 1]})
    fig.subplots_adjust(left=0.235, right=0.962, top=0.774, bottom=0.247, wspace=0.19)
    heading(fig, "Tree cover, rangeland, and daytime surface temperature",
            "Fitted difference for ten percentage points of tree cover replacing ten points of rangeland")
    rows = np.arange(len(labels))
    pair_counts = estimates.query("model == 'base'").set_index("analysis").pairs
    ticks = [f"{label}\n{pair_counts[key]} pairs" for key, label in labels.items()]
    axes[0].set_yticks(rows, ticks, fontsize=9)
    axes[0].set_ylim(len(rows) - 0.5, -0.5)
    for ax in axes:
        ax.axvline(0, color="#8d9b97", lw=1)
        ax.axhspan(-0.49, 0.49, color="#edf3ef", zorder=0)
        ax.grid(axis="x", color=GRID, linewidth=0.65, alpha=0.85)
        ax.set_axisbelow(True)
        ax.tick_params(axis="y", length=0)
        ax.spines["left"].set_visible(False)
    for i, kind in enumerate(labels):
        subset = estimates[estimates.analysis == kind].set_index("model")
        for key, offset, color in (("base", -0.13, GREEN), ("both", 0.13, BLUE)):
            row = subset.loc[key]
            interval_artist(axes[0], i + offset, row.estimate_c, row.lower_c, row.upper_c, color)
        row = subset.loc["both"]
        interval_artist(axes[1], i, row.change_c, row.change_lower_c, row.change_upper_c, BLUE)
    axes[0].set_title("LST contrast", loc="left", fontweight="bold", fontsize=12)
    axes[0].set_xlabel("Difference (°C)")
    axes[1].set_title("Change after adding distances", loc="left", fontweight="bold", fontsize=12)
    axes[1].set_xlabel("With distances minus base (°C)")
    handles = [
        Line2D([], [], color=GREEN, marker="o", markersize=5, label="Base model"),
        Line2D([], [], color=BLUE, marker="o", markersize=5, label="With coast and inland water distances"),
    ]
    fig.legend(handles=handles, loc="lower left", bbox_to_anchor=(0.217, 0.086), ncols=2,
               frameon=False, fontsize=9.3)
    fig.text(0.235, 0.064, "95% spatial-block bootstrap intervals. Each row compares models on identical pairs.",
             fontsize=9, color=MUTED)
    fig.savefig(output, dpi=180)
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--design", type=Path)
    source.add_argument("--analysis", type=Path)
    parser.add_argument("--map", type=Path, default=MAP)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    style()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    if args.design:
        support_figure(args.design, args.map, args.output)
    else:
        estimate_figure(args.analysis, args.output)
