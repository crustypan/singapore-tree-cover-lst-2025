"""Plot Stage 2 cover composition and geographic support before estimation."""
import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle
import numpy as np
import pandas as pd

from audit_coverage import BOUNDARY, read_boundary
from plot_landsat import GREEN, GRID, INK, MUTED, SOURCES, draw_map, heading, style
from prepare_local_comparison import BLOCK, DESIGN as LOCAL_DESIGN
from screen_landsat import GRID_CRS


HIGH = "#a54b2d"


def joint_figure(design, output):
    info = json.loads((design / "design_record.json").read_text())
    target = info["chosen"]
    low, high = target["lower_built_pct"], target["higher_built_pct"]
    start, end = target["tree_start_pct"], target["tree_end_pct"]
    density = pd.read_csv(design / "joint_density.csv").pivot(index="tree_from_pct", columns="built_from_pct",
                                                              values="cells")
    rows = pd.read_csv(design / "band_support.csv")
    rows.loc[~rows.feasible_reference, "endpoint_blocks"] = np.nan
    counts = rows.pivot(index="tree_start_pct", columns="built_reference_pct", values="endpoint_blocks")
    fig, (ax, right) = plt.subplots(1, 2, figsize=(11.6, 6.6))
    fig.subplots_adjust(left=0.085, right=0.94, bottom=0.24, top=0.76, wspace=0.37)
    heading(fig, "Where tree and built cover can be compared",
            "Singapore · one record per usable 500 m cell · seven Landsat overpasses in 2025")
    edges = np.arange(0, 101, 5)
    values = np.ma.masked_less(density.to_numpy(), 1)
    cloud = ax.pcolormesh(edges, edges, values, shading="flat", cmap="YlGnBu",
                          norm=LogNorm(vmin=1, vmax=values.max()))
    ax.plot([0, 100], [100, 0], color="#60746e", ls="--", lw=0.9)
    ax.set(xlim=(0, 100), ylim=(0, 100), xlabel="Built cover (% of classified land)",
           ylabel="Tree cover (% of classified land)")
    ax.set_aspect("equal", adjustable="box")
    ax.set_title("Joint distribution", loc="left", fontweight="bold", pad=9)
    bar = fig.colorbar(cloud, ax=ax, fraction=0.042, pad=0.03, shrink=0.86)
    bar.set_label("Distinct cells per bin", fontsize=9)
    bar.outline.set_visible(False)
    for b, color in [(low, GREEN), (high, HIGH)]:
        ax.plot([b, b], [start, end], "o-", color=color, mec="#fcfcf8", linewidth=2.4, ms=5.2, zorder=4)
        right.scatter(b, start, s=83, facecolors="none", edgecolors=color, linewidths=2.3, zorder=3)
    right.axhline(start, color="#879892", lw=0.8, ls="--", zorder=2, alpha=0.65)
    valid = np.ma.masked_invalid(counts.to_numpy())
    blocks = right.pcolormesh(
        np.r_[counts.columns.to_numpy() - 2.5, counts.columns[-1] + 2.5],
        np.r_[counts.index.to_numpy() - 1.25, counts.index[-1] + 1.25],
        valid, shading="flat", cmap="Blues", vmin=0, vmax=valid.max(),
    )
    right.set_xticks([25, 35, 45, 55, 65, 75, 85, 90])
    right.set_yticks(np.arange(0, 31, 5))
    right.set_xlabel("Built-cover reference (%)")
    right.set_ylabel("Tree cover at start of 10-point contrast (%)")
    right.set_title("Same-date support for both endpoints", loc="left", fontweight="bold", pad=9)
    bar = fig.colorbar(blocks, ax=right, fraction=0.042, pad=0.034, shrink=0.86)
    bar.set_label("Number of 4 km blocks", fontsize=9)
    bar.outline.set_visible(False)
    fig.legend([
        Line2D([0], [0], color=GREEN, marker="o"), Line2D([0], [0], color=HIGH, marker="o"),
    ], [f"{low:g}% built", f"{high:g}% built"], ncols=2, loc="lower center",
        bbox_to_anchor=(0.50, 0.113), frameon=False)
    fig.text(0.06, 0.077,
             "Support counts use cells near each reference and with an observed tree gap of 8–12 points.",
             fontsize=8.5, color=MUTED)
    fig.text(0.935, 0.03, SOURCES, ha="right", fontsize=8, color=MUTED)
    fig.savefig(output, dpi=180)
    plt.close(fig)


def area_figure(design, output):
    target = json.loads((design / "design_record.json").read_text())["chosen"]
    cells = pd.read_csv(LOCAL_DESIGN / "cell_design.csv").query("usable_dates > 0")
    support = pd.read_csv(design / "selected_block_dates.csv")
    sample = pd.read_csv(design / "fit_sample.csv").drop_duplicates("cell_id").merge(
        cells[["cell_id", "easting_m", "northing_m"]], on="cell_id", validate="one_to_one",
    )
    boundary, _ = read_boundary(BOUNDARY, GRID_CRS)
    fig, (ax, bar) = plt.subplots(1, 2, figsize=(11.6, 6.0), gridspec_kw={"width_ratios": [1.16, 1]})
    fig.subplots_adjust(left=0.085, right=0.94, bottom=0.26, top=0.78, wspace=0.30)
    heading(fig, "Geographic support for the selected contrast",
            f"Tree cover {target['tree_start_pct']:g}% to {target['tree_end_pct']:g}% "
            f"at {target['lower_built_pct']:g}% and {target['higher_built_pct']:g}% built cover")
    draw_map(ax, boundary, cells, np.ones(len(cells)),
             matplotlib.colors.ListedColormap(["#e2e6df"]), matplotlib.colors.Normalize(0, 1))
    origin = 500 * np.floor(np.asarray(boundary.bounds[:2]) / 500) - 500
    for band, col in [("lower", GREEN), ("higher", HIGH)]:
        z = sample[sample.band == band]
        pts = (z[["easting_m", "northing_m"]].to_numpy() - origin) / 1000
        ax.scatter(*pts.T, s=4, color=col, alpha=0.70, linewidth=0, rasterized=True)
    shared = support[support.shared].drop_duplicates(BLOCK)
    for row in shared.itertuples():
        x, y = (np.array([row.block_e, row.block_n]) * 4000 - origin) / 1000
        ax.add_patch(Rectangle((x, y), 4, 4, fill=False, lw=1.05, edgecolor=INK))
    x0, y0 = ax.get_xlim()[0] + 2, ax.get_ylim()[0] + 0.7
    ax.plot([x0, x0 + 5], [y0, y0], color=INK, lw=1.3)
    ax.text(x0 + 2.5, y0 + 0.6, "5 km", ha="center", fontsize=8)
    ax.set_title("Cells in block-dates that support an endpoint contrast", loc="left",
                 fontweight="bold", fontsize=10.5, pad=6)
    per_date = support.groupby("date", sort=True)[["lower", "higher", "shared"]].sum()
    idx = np.arange(len(per_date))
    for column, color, offset in [("lower", GREEN, -0.23), ("higher", HIGH, 0), ("shared", INK, 0.23)]:
        bar.bar(idx + offset, per_date[column], width=0.22, color=color)
    bar.set_xticks(idx, pd.to_datetime(per_date.index).strftime("%d %b"), rotation=25, ha="right")
    bar.set_ylabel("4 km blocks")
    bar.set_title("Blocks with both tree endpoints, by date", loc="left", fontweight="bold", fontsize=10.5, pad=8)
    bar.grid(axis="y", linewidth=0.6, color=GRID)
    bar.set_axisbelow(True)
    fig.legend(
        [Patch(color=GREEN), Patch(color=HIGH), Patch(facecolor="none", edgecolor=INK)],
        [f"{target['lower_built_pct']:g}% built band", f"{target['higher_built_pct']:g}% built band",
         "Both bands on a common date"], ncols=3, loc="lower center", bbox_to_anchor=(0.50, 0.105),
        frameon=False, fontsize=9,
    )
    fig.text(0.935, 0.03, SOURCES, ha="right", fontsize=8, color=MUTED)
    fig.savefig(output, dpi=180)
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--design", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    style()
    joint_figure(args.design, args.output / "joint_support.png")
    area_figure(args.design, args.output / "geographic_support.png")
