"""Plot retained 4 K support and the two-block omission contrasts."""
import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle
import numpy as np
import pandas as pd

from audit_coverage import BOUNDARY, read_boundary
from plot_landsat import GREEN, GRID, INK, MUTED, draw_map, heading, style
from prepare_grass_tree_substitution import DESIGN as STAGE5, SIDES
from prepare_local_comparison import BLOCK, COORDS
from screen_landsat import GRID_CRS


CASES = (("dense", "Primary matching"), ("wider_balance", "Wider matching"))
STAGES = (("original", "Original"), ("qa4_full", "4 K screen"), ("qa4_omit_both", "4 K and omit both"))
BLUE = "#24789a"
GOLD = "#b4762c"


def support_map(ax, screen):
    passed = screen[screen.passes_qa4]
    cells = pd.concat([
        passed[[side, "in_omitted_blocks"]].rename(columns={side: "cell_id"}) for side in SIDES
    ], ignore_index=True).drop_duplicates()
    locations = pd.read_csv(STAGE5 / "cell_covariates.csv", usecols=["cell_id", *COORDS]).drop_duplicates()
    cells = cells.merge(locations, on="cell_id", validate="one_to_one")
    boundary, _ = read_boundary(BOUNDARY, GRID_CRS)
    colors = ListedColormap([GREEN, GOLD])
    draw_map(ax, boundary, cells, cells.in_omitted_blocks.astype(int), colors, BoundaryNorm([-.5, .5, 1.5], 2))
    x0, y0 = 500 * np.floor(np.array(boundary.bounds[:2]) / 500) - 500
    blocks = passed[[*BLOCK, "in_omitted_blocks"]].drop_duplicates()
    for east, north, omitted in blocks.itertuples(index=False, name=None):
        axis_coords = ((east * 4000 - x0) / 1000, (north * 4000 - y0) / 1000)
        ax.add_patch(Rectangle(axis_coords, 4, 4, fill=False, linewidth=.75, alpha=.85,
                               edgecolor=GOLD if omitted else GREEN))
    ax.set_title("Blocks with at least one 4 K pair", loc="left", fontsize=11.8, fontweight="bold", pad=6)
    return [Patch(fc=GREEN, label="Retained blocks"), Patch(fc=GOLD, label="Two omitted blocks")]


def support_figure(design, path):
    counts = pd.read_csv(design / "support_summary.csv")
    pairs = pd.read_csv(design / "pair_coverage.csv")
    daily = pd.read_csv(design / "support_by_date.csv")
    fig = plt.figure(figsize=(12.6, 8.4))
    grid = fig.add_gridspec(2, 2, left=.079, right=.963, top=.83, bottom=.145,
                            height_ratios=[1.25, .88], width_ratios=[1, 1.15], hspace=.53, wspace=.24)
    heading(fig, "How much dense-area pairing support survives the 4 K screen",
            "Both 500 m cells need ≥60% qualifying study-land coverage · minimum distance to cloud 0.5 km")
    ax = fig.add_subplot(grid[0, 0])
    x = np.arange(len(STAGES))
    for (case, label), delta, color in zip(CASES, (-.17, .17), (GREEN, BLUE)):
        table = counts[counts.case.eq(case)].set_index("sample").loc[[s[0] for s in STAGES]]
        bars = ax.bar(x + delta, table.pairs, .31, color=color, label=label)
        ax.bar_label(bars, [f"{v}" for v in table.pairs], padding=3, fontsize=9, color=INK)
    ax.set_xticks(x, [title for _, title in STAGES], fontsize=9)
    ax.set_ylabel("Pair observations")
    ax.set_ylim(0, max(55, counts.pairs.max() * 1.2))
    ax.set_title("Pair counts by design and screen", loc="left", fontsize=11.8, fontweight="bold")
    ax.grid(axis="y", color=GRID, lw=.6, alpha=.85)
    ax.set_axisbelow(True)
    ax.legend(frameon=False, fontsize=8.5, loc="upper right")

    ax = fig.add_subplot(grid[0, 1])
    handles = support_map(ax, pairs)
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(.5, -.015), ncols=2,
              frameon=False, fontsize=8.7)
    ax = fig.add_subplot(grid[1, :])
    dates = sorted(daily.date.unique())
    x = np.arange(len(dates))
    for (case, label), delta, color in zip(CASES, (-.14, .14), (GREEN, BLUE)):
        refined = daily[daily.case.eq(case)].set_index("date").reindex(dates, fill_value=0)
        ax.bar(x + delta, refined.qa4_full_pairs, color=color, alpha=.27, width=.25)
        ax.bar(x + delta, refined.qa4_omit_both_pairs, color=color, width=.25)
    ax.set_xticks(x, pd.to_datetime(dates).strftime("%d %b\n%Y"), fontsize=8.1)
    ax.set_ylabel("Pair observations")
    ax.set_title("Passing comparisons by overpass", loc="left", fontsize=11.8, fontweight="bold")
    ax.set_yticks(np.arange(0, max(4, daily.qa4_full_pairs.max()) + 1))
    ax.grid(axis="y", color=GRID, lw=.6, alpha=.85)
    ax.set_axisbelow(True)
    ax.legend(handles=[Patch(fc=GREEN, label="Primary matching"), Patch(fc=BLUE, label="Wider matching")],
              frameon=False, fontsize=8.6, ncols=2, loc="upper right")
    fig.text(.079, .064, "Pale bars show the full 4 K sample. Dark bars show pairs retained after omitting both blocks.",
             fontsize=9, color=MUTED)
    fig.text(.963, .027, "Sources  USGS Landsat · ESA WorldCover · Esri / Impact Observatory · URA",
             fontsize=8, ha="right", color=MUTED)
    fig.savefig(path, dpi=175)
    plt.close(fig)


def point(ax, y, estimate, lo, hi, color):
    ax.hlines(y, lo, hi, color=color, lw=1.7)
    ax.scatter(estimate, y, s=32, color=color, linewidths=0, zorder=3)


def temperature_figure(analysis, path):
    estimates = pd.read_csv(analysis / "estimates.csv").query("model == 'matched'")
    changes = pd.read_csv(analysis / "omission_changes.csv").query("model == 'matched'")
    fig, (before, shift) = plt.subplots(1, 2, figsize=(12, 5.4), sharey=True,
                                        gridspec_kw={"width_ratios": [1.22, 1]})
    fig.subplots_adjust(left=.195, right=.954, top=.745, bottom=.26, wspace=.19)
    heading(fig, "Five-point temperature contrasts under the 4 K screen",
            "Daytime Landsat LST · tree cover higher and rangeland approximately five points lower")
    labels = []
    for i, (key, label) in enumerate(CASES):
        rows = estimates[estimates.case.eq(key)].set_index("sample")
        labels.append(f"{label}\n{rows.loc['qa4_full', 'pairs']} to {rows.loc['qa4_omit_both', 'pairs']} pairs")
        for sample, y, color in (("qa4_full", i - .13, MUTED), ("qa4_omit_both", i + .13, GREEN)):
            row = rows.loc[sample]
            if np.isfinite(row.estimate_c_per_5pp):
                point(before, y, row.estimate_c_per_5pp, row.lower_c, row.upper_c, color)
        diff = changes[changes.case.eq(key)]
        if not diff.empty:
            row = diff.iloc[0]
            point(shift, i, row.change_c_per_5pp, row.lower_c, row.upper_c, BLUE)
    before.set_yticks([0, 1], labels, fontsize=9)
    before.set_ylim(1.51, -.51)
    for ax in (before, shift):
        ax.axvline(0, color="#8da09b", lw=.85)
        ax.grid(axis="x", color=GRID, lw=.6, alpha=.8)
        ax.set_axisbelow(True)
        ax.spines["left"].set_visible(False)
        ax.tick_params(axis="y", length=0)
    before.set_title("Matched LST contrast", loc="left", fontweight="bold", fontsize=11.8)
    before.set_xlabel("Higher-tree minus lower-tree LST (°C)")
    shift.set_title("Change after omitting the two blocks", loc="left", fontweight="bold", fontsize=11.8)
    shift.set_xlabel("Omitted minus full sample (°C)")
    handles = [
        Line2D([], [], color=MUTED, marker="o", label="Full 4 K sample"),
        Line2D([], [], color=GREEN, marker="o", label="4 K sample without both blocks"),
    ]
    fig.legend(handles=handles, loc="lower left", bbox_to_anchor=(.179, .097), frameon=False, ncols=2, fontsize=9.0)
    fig.text(.195, .061, "Intervals use 95% spatial block resampling. The change pairs the before and after draws.",
             fontsize=8.8, color=MUTED)
    fig.savefig(path, dpi=175)
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--design", type=Path, required=True)
    parser.add_argument("--analysis", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    style()
    args.output.mkdir(parents=True, exist_ok=True)
    support_figure(args.design, args.output / "quality_support.png")
    temperature_figure(args.analysis, args.output / "quality_temperature_contrasts.png")
