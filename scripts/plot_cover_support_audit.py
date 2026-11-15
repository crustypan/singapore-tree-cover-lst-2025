"""Plot the cover-filter counts and the original reflectance for each review cell."""
import argparse
from pathlib import Path
import warnings

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.patches import Patch, Rectangle
import numpy as np
import pandas as pd

from analyze_cover_support_audit import SNAPSHOT
from audit_coverage import BOUNDARY, read_boundary
from plot_landsat import GREEN, GRID, INK, MUTED, draw_map, heading, style
from prepare_cover_support_audit import AUDIT_DESIGN, support_inputs
from prepare_overpass_validation import YEARS
from screen_landsat import GRID_CRS


RUST = "#aa5836"
NEUTRAL = "#e8e8e2"
CLASS_COLORS = [NEUTRAL, "#17634c", "#a95239", "#90bbd5", "#afd18a", "#e7bc78", "#ada7bb"]
CLASS_NAMES = ["Masked", "Trees", "Built", "Water", "Other vegetation", "Bare", "Other"]
CLASSES = ListedColormap(CLASS_COLORS)
CLASS_NORM = BoundaryNorm(np.arange(-0.5, 7.5), CLASSES.N)
ATTRIBUTION = "Sources  URA/data.gov.sg · ESA WorldCover · Esri / Impact Observatory · Copernicus Sentinel-2 · Dynamic World"


def save(fig, outfile):
    fig.text(0.95, 0.020, ATTRIBUTION, ha="right", fontsize=7.8, color=MUTED)
    fig.savefig(outfile, dpi=180)
    plt.close(fig)


def plot_counts(design, analysis, output):
    a = pd.read_csv(design / "esri_filter_counts.csv").query("year_group == 'all'").set_index("analysis")
    b = pd.read_csv(analysis / "matched_support_counts.csv")
    b = b.query("year_group == 'all' and footprint == 'fixed'").set_index(["confidence", "analysis", "stable"])
    left_index = ["none", "tree", "built", "water", "tree+built+water"]
    right_index = [(0.6, "esri", False), (0.6, "dw", False), (0.0, "esri", False), (0.0, "dw", False),
                   (0.0, "dw", True)]
    axes_rows = [
        (a.loc[left_index], ["All common cells", "Stable tree share", "Stable built share", "Stable water share",
                            "All three stable"], "Esri stability conditions"),
        (b.loc[right_index], [
            "Esri\nScore ≥0.60", "Dynamic World\nScore ≥0.60", "Esri\nAll clear labels",
            "Dynamic World\nAll clear labels", "Dynamic World\nAll clear, stable shares",
        ], "Comparison on the same optical footprint"),
    ]
    fig, axes = plt.subplots(1, 2, figsize=(13.0, 7.05), sharex=True)
    fig.subplots_adjust(left=0.158, right=0.955, bottom=0.20, top=0.79, wspace=0.65)
    heading(fig, "Comparable blocks at 65% and 85% built cover",
            "Singapore · fixed 500 m cells · tree cover 2.5% to 12.5% · 15 Landsat overpasses")
    legend = [("lower_endpoint_blocks", "65% built", GREEN), ("higher_endpoint_blocks", "85% built", RUST),
              ("shared_blocks", "Both built bands on a common date", INK)]
    max_count = max(rows[["lower_endpoint_blocks", "higher_endpoint_blocks"]].max().max()
                    for rows, _, _ in axes_rows)
    for axis, (rows, labels, title) in zip(axes, axes_rows):
        y = np.arange(len(rows))
        for offset, (column, _, color) in zip((-0.235, 0, 0.235), legend):
            axis.barh(y + offset, rows[column], height=0.20, color=color, linewidth=0, zorder=2)
            for value, index in zip(rows[column], y):
                axis.text(value + 0.48, index + offset, str(value), ha="left", va="center",
                          color=color, fontsize=8)
        axis.set_yticks(y, labels)
        axis.set_ylim(len(rows) - 0.52, -0.52)
        axis.set_xlim(0, max_count + 3)
        axis.set_title(title, loc="left", fontweight="bold", pad=8, fontsize=11.5)
        axis.set_xlabel("Distinct 4 km blocks")
        axis.tick_params(axis="y", length=0)
        axis.spines["left"].set_visible(False)
        axis.grid(axis="x", color=GRID, linewidth=0.62)
        axis.set_axisbelow(True)
    fig.legend(handles=[Patch(color=color, label=label) for _, label, color in legend],
               loc="lower center", bbox_to_anchor=(0.54, 0.085), ncols=3, frameon=False)
    n = int(b.loc[(0.6, "dw", False), "cells"])
    m = int(b.loc[(0.0, "dw", False), "cells"])
    fig.text(0.06, 0.060, f"The shared footprint includes {n:,} cells at score ≥0.60 and {m:,} with all clear labels. "
             "Both use the same accepted pixels across all years.", fontsize=8.7, color=MUTED)
    save(fig, output)


def common_color(labels, source):
    conversion = ({2: 1, 7: 2, 1: 3, 4: 4, 5: 4, 11: 4, 8: 5, 9: 6} if source == "esri"
                  else {1: 1, 6: 2, 0: 3, 2: 4, 3: 4, 4: 4, 5: 4, 7: 5, 8: 6})
    out = np.zeros(labels.shape, dtype=np.uint8)
    for key, category in conversion.items():
        out[labels == key] = category
    return out


def rgb(chips, year):
    dv = chips[f"case_dn_{year}"][..., [2, 1, 0]].astype(np.float32)
    clear = chips[f"case_clear_{year}"]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        reflectance = np.nanmedian(np.where(clear[..., None], dv, np.nan), axis=0)
    view = (np.clip(np.nan_to_num(reflectance) / 3400, 0, 1)) ** 0.78
    view[clear.sum(axis=0) < 2] = [0.91, 0.91, 0.89]
    return view


def plot_chips(design, snapshot, out):
    cases = pd.read_csv(design / "critical_cells.csv")
    chips = np.load(snapshot / "case_chips.npz", allow_pickle=False)
    np.testing.assert_array_equal(chips["cases"], cases.case)
    panels = {}
    for year in YEARS:
        esri = common_color(chips[f"esri_{year}"], "esri")
        dw = common_color(chips[f"case_classes_{year}"], "dw")
        enough = chips[f"case_count_{year}"] >= 2
        dw[~enough] = 0
        drawn = CLASSES(dw)[..., :3]
        uncertain = enough & (chips[f"case_confidence_{year}"] < 0.60)
        drawn[uncertain] = 0.49 * drawn[uncertain] + 0.51
        panels[year] = (rgb(chips, year), esri, drawn)

    for offset in range(0, len(cases), 4):
        subset = cases.iloc[offset:offset + 4]
        fig = plt.figure(figsize=(13.4, 1.55 * len(subset) + 2.75))
        grid = fig.add_gridspec(len(subset), 12, left=0.059, right=0.962, bottom=0.163, top=0.802,
                               width_ratios=[1.33, 1, 1, 1, 0.17, 1, 1, 1, 0.17, 1, 1, 1],
                               hspace=0.16, wspace=0.10)
        heading(fig, f"Sentinel-2 cover audit for cells {subset.case.iloc[0]}–{subset.case.iloc[-1]}",
                "500 m windows · June–September of each year · fixed reflectance stretch")
        for r, (i, case) in enumerate(subset.iterrows()):
            label = fig.add_subplot(grid[r, 0])
            label.axis("off")
            label.text(0, 0.68, case.case, color=INK, fontsize=11, fontweight="bold")
            label.text(0, 0.47, f"{int(case.easting_m)} E\n{int(case.northing_m)} N",
                       fontsize=8, linespacing=1.3, color=MUTED)
            for year, start in zip(YEARS, (1, 5, 9)):
                for j, name in enumerate(("Sentinel-2 RGB", "Esri", "Dynamic World")):
                    ax = fig.add_subplot(grid[r, start + j])
                    color = {"cmap": CLASSES, "norm": CLASS_NORM} if j == 1 else {}
                    ax.imshow(panels[year][j][i], interpolation="nearest", **color)
                    ax.set_xticks([])
                    ax.set_yticks([])
                    ax.spines[:].set_visible(False)
                    if r == 0:
                        ax.set_title(name, color=MUTED, fontsize=8.8, pad=3)
                        if j == 1:
                            pos = ax.get_position()
                            fig.text((pos.x0 + pos.x1) / 2, pos.y1 + 0.032, str(year), ha="center",
                                     fontsize=10.5, fontweight="bold")
        legend = [Patch(facecolor=color, label=name) for color, name in zip(CLASS_COLORS, CLASS_NAMES)]
        fig.legend(handles=legend, loc="lower center", bbox_to_anchor=(0.52, 0.069), ncols=7, frameon=False,
                   fontsize=8.8)
        fig.text(0.06, 0.057, "RGB uses median clear reflectance. Pale Dynamic World colors show mean probability below 0.60.",
                 fontsize=8.3, color=MUTED)
        save(fig, out / f"{subset.case.iloc[0]}-{subset.case.iloc[-1]}.png")


def plot_sites(design, path):
    boundary, _ = read_boundary(BOUNDARY, GRID_CRS)
    g, _, _, _ = support_inputs()
    cases = pd.read_csv(design / "critical_cells.csv")
    fig, ax = plt.subplots(figsize=(10.6, 6.6))
    fig.subplots_adjust(left=0.065, right=0.96, bottom=0.10, top=0.80)
    heading(fig, "The 19 cells at the missing high-tree endpoint",
            "Nominal 65% built and 12.5% trees in at least one annual Esri classification")
    draw_map(ax, boundary, g, np.ones(len(g)),
             ListedColormap(["#e5ebe6"]), matplotlib.colors.Normalize(0, 1))
    origin = 500 * np.floor(np.asarray(boundary.bounds[:2]) / 500) - 500
    p = (cases[["easting_m", "northing_m"]].to_numpy() - origin) / 1000
    ax.scatter(*p.T, s=25, color=RUST, edgecolors="#fcfcf8", linewidths=0.5, zorder=3)
    offsets = {"C01": (-7, 7), "C02": (6, -9), "C03": (5, 8), "C04": (6, -9), "C09": (-6, 7),
               "C10": (7, -10), "C13": (-7, -6), "C14": (8, 8), "C15": (-7, 5), "C16": (4, -9),
               "C17": (8, -7), "C18": (5, 8)}
    for (x, y), name in zip(p, cases.case):
        dx, dy = offsets.get(name, (4, 6))
        ax.annotate(name, (x, y), xytext=(dx, dy), textcoords="offset points",
                    fontsize=7.9, color=INK, ha="left" if dx >= 0 else "right",
                    va="bottom" if dy >= 0 else "top")
    save(fig, path)


def plot_original_scenes(snapshot, review, out):
    rows = pd.read_csv(review).query("illustrated").reset_index(drop=True)
    if rows.empty:
        return
    data = np.load(snapshot / "case_chips.npz", allow_pickle=False)
    i = {name: p for p, name in enumerate(data["cases"])}
    fig = plt.figure(figsize=(10.8, 2.25 * len(rows) + 2.35))
    g = fig.add_gridspec(len(rows), 4, left=0.065, right=0.95, top=0.80, bottom=0.105,
                         wspace=0.24, hspace=0.34, width_ratios=[0.9, 1, 1, 1])
    heading(fig, "Original Sentinel-2 views of reviewed surface changes",
            "500 m windows · one observation per year · identical scale and outlined locations")
    for j, row in enumerate(rows.itertuples()):
        p = i[row.case]
        x0, y0, x1, y1 = row.roi_col_start, row.roi_row_start, row.roi_col_end, row.roi_row_end
        region = np.s_[int(y0):int(y1), int(x0):int(x1)]
        label = fig.add_subplot(g[j, 0])
        label.axis("off")
        label.text(0, 0.68, row.case, fontsize=13, fontweight="bold")
        label.text(0, 0.48, row.short_label.replace(" / ", "\n"), color=MUTED, fontsize=10,
                   va="top", linespacing=1.35)
        for column, year in enumerate(YEARS, start=1):
            scene_clear = data[f"case_clear_{year}"][:, p]
            local = scene_clear[:, region[0], region[1]].mean(axis=(1, 2))
            scores = local + 0.01 * scene_clear.mean(axis=(1, 2))
            scene = int(np.argmax(scores))
            a = data[f"case_dn_{year}"][scene, p][:, :, [2, 1, 0]]
            image = (np.clip(a.astype(np.float32) / 3400, 0, 1)) ** 0.78
            image[~scene_clear[scene]] = [0.91, 0.91, 0.89]
            ax = fig.add_subplot(g[j, column])
            ax.imshow(image, interpolation="nearest")
            ax.add_patch(Rectangle((x0 - 0.5, y0 - 0.5), x1 - x0, y1 - y0, lw=1.1,
                                   edgecolor="#ffe388", facecolor="none"))
            dates = data[f"dates_{year}"]
            ax.set_title(str(dates[scene]), fontsize=9.4, pad=5)
            ax.set_xticks([])
            ax.set_yticks([])
            ax.spines[:].set_visible(False)
    save(fig, out)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--design", type=Path, default=AUDIT_DESIGN)
    parser.add_argument("--snapshot", type=Path, default=SNAPSHOT)
    parser.add_argument("--analysis", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--review", type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    chips_dir = args.output / "chips"
    chips_dir.mkdir(exist_ok=True)
    style()
    plot_sites(args.design, args.output / "audit_cells.png")
    plot_counts(args.design, args.analysis, args.output / "support_counts.png")
    plot_chips(args.design, args.snapshot, chips_dir)
    if args.review is not None:
        plot_original_scenes(args.snapshot, args.review, args.output / "original_scene_examples.png")
