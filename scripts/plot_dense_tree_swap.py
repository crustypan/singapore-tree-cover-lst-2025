"""Draw the outcome-free support audit and supported temperature contrasts."""
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
from plot_landsat import GREEN, GRID, INK, MUTED, draw_map, heading, style
from prepare_dense_tree_swap import STAGE5, STEPS
from screen_landsat import GRID_CRS


BANDS = {
    "lt50": "<50%", "50to70": "50–<70%", "70to80": "70–<80%",
    "80to90": "80–<90%", "ge90": "90–100%", "ge80": "80–100%",
}
CASES = {
    "dense": "80–100% built · primary",
    "very_dense": "90–100% built",
    "80_to_90": "80–<90% built",
    "wider_balance": "80–100% built · wider balance limits",
}
GOLD = "#ac752a"
OCHRE = "#dfbf82"
BLUE = "#21799b"
GRAY = "#dce3df"
SOURCES = "Sources  USGS Landsat · ESA WorldCover · Esri / Impact Observatory · URA"


def support_figure(design, output):
    shares = pd.read_csv(design / "composition_by_band.csv").query("built_band != 'ge80'")
    funnel = pd.read_csv(design / "support_by_band.csv")
    fig = plt.figure(figsize=(12.7, 9.0))
    grid = fig.add_gridspec(2, 2, left=0.077, right=0.967, top=0.83, bottom=0.14,
                           hspace=0.50, wspace=0.27, height_ratios=[0.92, 1])
    heading(fig, "Rangeland and pairing support across built-cover bands",
            "Singapore · 500 m cells · annual Esri land cover and 15 daytime Landsat overpasses, 2023–25")
    ax = fig.add_subplot(grid[0, :])
    levels = np.column_stack([
        shares.zero, shares.below_5 - shares.zero, shares.at_least_5 - shares.at_least_10, shares.at_least_10,
    ]) / shares.cell_years.to_numpy()[:, None] * 100
    bx = np.arange(len(shares))
    total = np.zeros(len(shares))
    options = [(GRAY, "0%"), (OCHRE, ">0 and <5%"), ("#72ada0", "5–<10%"), (GREEN, "≥10%")]
    for layer, (color, label) in zip(levels.T, options):
        ax.bar(bx, layer, bottom=total, width=0.64, color=color, label=label)
        total += layer
    ax.set_title("Rangeland available in eligible cells", loc="left", fontweight="bold")
    ax.set_xticks(bx, [f"{BANDS[b]}\n{n:,} cell-years" for b, n in zip(shares.built_band, shares.cell_years)])
    ax.set_ylim(0, 112)
    ax.set_ylabel("Share of cell-years (%)")
    ax.set_yticks([0, 25, 50, 75, 100])
    ax.legend(loc="upper right", ncols=4, frameon=False, fontsize=9)
    ax.grid(axis="y", color=GRID, linewidth=0.65, alpha=0.75)
    ax.set_axisbelow(True)

    labels = ["Eligible", "Enough\nrangeland", "Opposite\ncover gaps", "Cover\nbalance",
              "Within\n2 km", "Same 4 km\nblock", "Disjoint\npairs"]
    for col, band in enumerate(("80to90", "ge90")):
        ax = fig.add_subplot(grid[1, col])
        for name, color, ls, label in [("10pp", GOLD, "--", "10 points"), ("5pp", GREEN, "-", "5 points")]:
            values = (funnel[funnel.comparison.eq(name) & funnel.built_band.eq(band)]
                      .set_index("step").loc[list(STEPS), "cell_dates"].to_numpy())
            ax.plot(np.arange(len(labels)), np.log10(values + 1), color=color, linestyle=ls,
                    marker="o", markersize=4, lw=1.6, label=label)
        cy = [0, 10, 100, 1000, 10000]
        ax.set_yticks(np.log10(np.array(cy) + 1), [f"{n:,}" for n in cy])
        ax.set_ylim(-0.14, 4.53)
        ax.set_xticks(np.arange(len(labels)), labels, fontsize=8.4)
        ax.set_title(f"{BANDS[band]} built cover", loc="left", fontweight="bold")
        ax.grid(axis="y", color=GRID, linewidth=0.65, alpha=0.85)
        ax.set_axisbelow(True)
        if col == 0:
            ax.set_ylabel("Lower-tree cell-dates")
        else:
            ax.legend(loc="upper right", frameon=False, fontsize=9)
    fig.text(0.077, 0.055, "Cumulative requirements. Each lower-tree cell is counted once per date at each step it reaches.",
             fontsize=9, color=MUTED)
    fig.text(0.967, 0.028, SOURCES, ha="right", fontsize=8, color=MUTED)
    fig.savefig(output, dpi=180)
    plt.close(fig)


def geography_figure(design, stage5, output):
    registry = pd.read_csv(design / "dense_cells.csv.gz")
    pair_list = pd.read_csv(design / "pairs.csv").query("built_band == 'ge80'")
    names = ["cell_id", "easting_m", "northing_m"]
    positions = pd.read_csv(stage5 / "cell_covariates.csv", usecols=names).drop_duplicates().set_index("cell_id")
    boundary, _ = read_boundary(BOUNDARY, GRID_CRS)
    x0, y0 = 500 * np.floor(np.array(boundary.bounds[:2]) / 500) - 500
    fig, axes = plt.subplots(1, 2, figsize=(12.7, 5.3))
    fig.subplots_adjust(left=0.06, right=0.972, top=0.76, bottom=0.17, wspace=0.13)
    heading(fig, "Where dense built areas have a local comparison",
            "Both endpoints have ≥80% built cover · cells and pairs shown if they qualify on any overpass")
    colors = ListedColormap([GRAY, OCHRE, GREEN])
    norm = BoundaryNorm(np.arange(-0.5, 3), colors.N)
    for ax, key, label in zip(axes, ("10pp", "5pp"), ("Ten-point benchmark", "Five-point comparison")):
        rows = registry[registry.comparison.eq(key)]
        cells = rows.assign(donor=rows.last_supported_step.ne("eligible")).groupby("cell_id").donor.any()
        pairs = pair_list[pair_list.comparison.eq(key)]
        endpoints = pd.Index(pd.concat([pairs.lower_tree_cell, pairs.higher_tree_cell]).unique())
        status = cells.astype(int).mask(cells.index.isin(endpoints), 2)
        draw_map(ax, boundary, positions.loc[cells.index], status, colors, norm)
        for lo, hi in pairs[["lower_tree_cell", "higher_tree_cell"]].drop_duplicates().itertuples(index=False):
            points = positions.loc[[lo, hi]]
            ax.plot((points.easting_m - x0) / 1000, (points.northing_m - y0) / 1000,
                    color=GREEN, lw=0.75, alpha=0.63)
        ax.set_title(f"{label}\n{len(pairs)} paired observations", loc="left", fontweight="bold",
                     fontsize=11, pad=8, linespacing=1.46)
    ax = axes[0]
    left, lower = ax.get_xlim()[0] + 1.5, ax.get_ylim()[0] + 0.7
    ax.plot([left, left + 5], [lower, lower], color=INK, lw=1.2)
    ax.text(left + 2.5, lower + 0.65, "5 km", ha="center", fontsize=8)
    fig.legend(handles=[
        Patch(fc=GRAY, label="Eligible dense cell"),
        Patch(fc=OCHRE, label="Enough rangeland on ≥1 date"),
        Patch(fc=GREEN, label="Paired on ≥1 date"),
    ], loc="lower center", bbox_to_anchor=(0.5, 0.086), frameon=False, ncols=3, fontsize=9)
    fig.text(0.967, 0.027, SOURCES, ha="right", fontsize=8, color=MUTED)
    fig.savefig(output, dpi=180)
    plt.close(fig)


def effects_figure(analysis, output):
    data = pd.read_csv(analysis / "estimates.csv")
    fig, ax = plt.subplots(figsize=(11.5, 5.6))
    fig.subplots_adjust(left=0.34, right=0.956, top=0.791, bottom=0.235)
    heading(fig, "Daytime surface temperature in supported dense-area pairs",
            "Five percentage points more tree cover with approximately five points less rangeland")
    y = np.arange(len(CASES))
    labels = []
    ax.axvline(0, color="#919e98", lw=1)
    ax.axhspan(-0.48, 0.48, color="#e9f1ec", zorder=0)
    for i, (case, label) in enumerate(CASES.items()):
        rows = data[data.case.eq(case)].set_index("model")
        n, g = rows.iloc[0][["pairs", "blocks"]]
        labels.append(f"{label}\n{n} pair observations · {g} blocks")
        if rows.estimate_c_per_5pp.isna().all():
            ax.text(0.03, i, "Sparse spatial support", transform=ax.get_yaxis_transform(),
                    va="center", color=MUTED, fontsize=9)
        for model, color, dy in [("matched", GREEN, -0.09), ("adjusted", BLUE, 0.12)]:
            row = rows.loc[model]
            if np.isfinite(row.estimate_c_per_5pp):
                height = i + dy
                ax.hlines(height, row.lower_c, row.upper_c, color=color, lw=1.8)
                ax.scatter(row.estimate_c_per_5pp, height, s=32, color=color, zorder=5, linewidths=0)
    ax.set_yticks(y, labels, fontsize=9.25)
    ax.set_ylim(len(y) - 0.5, -0.5)
    ax.set_xlabel("Higher-tree minus lower-tree surface temperature (°C per 5 points)")
    ax.grid(axis="x", color=GRID, linewidth=0.65, alpha=0.85)
    ax.set_axisbelow(True)
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    fig.legend(handles=[
        Line2D([], [], color=GREEN, marker="o", markersize=5, label="Matched contrast"),
        Line2D([], [], color=BLUE, marker="o", markersize=5, label="Adjusted for residual cover gaps"),
    ], loc="lower left", bbox_to_anchor=(0.32, 0.083), frameon=False, ncols=2, fontsize=9)
    fig.text(0.34, 0.06, "95% spatial-block resampling intervals. Each block includes its repeated dates.",
             fontsize=8.7, color=MUTED)
    fig.savefig(output, dpi=180)
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--design", type=Path)
    source.add_argument("--analysis", type=Path)
    parser.add_argument("--stage5", type=Path, default=STAGE5)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    style()
    args.output.mkdir(parents=True, exist_ok=True)
    if args.design:
        support_figure(args.design, args.output / "rangeland_and_pair_support.png")
        geography_figure(args.design, args.stage5, args.output / "dense_pair_map.png")
    else:
        effects_figure(args.analysis, args.output / "five_point_estimates.png")
