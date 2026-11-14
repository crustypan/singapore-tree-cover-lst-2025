"""Create reproducible maps and statistical figures for the selected overpasses."""
import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import BoundaryNorm, ListedColormap, TwoSlopeNorm
from matplotlib.ticker import MaxNLocator, PercentFormatter

from analyze_landsat import load_panel
from audit_coverage import BOUNDARY, read_boundary
from screen_landsat import GRID_CRS

INK = "#162c35"
MUTED = "#50656a"
GREEN = "#0d6a5e"
LIGHT = "#fcfcf8"
GRID = "#d9e1dc"
SOURCES = "Sources  USGS Landsat · URA/data.gov.sg · ESA WorldCover · Esri / Impact Observatory"


def style():
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 10, "axes.spines.top": False,
        "axes.spines.right": False, "text.color": INK, "axes.labelcolor": INK,
        "xtick.color": MUTED, "ytick.color": MUTED, "axes.edgecolor": "#b9c5bf",
        "savefig.facecolor": LIGHT, "figure.facecolor": LIGHT, "axes.facecolor": LIGHT,
    })


def heading(fig, title, detail, y=0.962):
    fig.text(0.055, y, title, ha="left", va="top", fontsize=22, fontweight="bold")
    fig.text(0.055, y - 0.063, detail, ha="left", va="top", fontsize=10.5, color=MUTED)


def draw_map(ax, boundary, data, val, cmap, norm):
    west, south, east, north = boundary.bounds
    x0, y0 = 500 * np.floor(np.array([west, south]) / 500) - 500
    width = int(np.ceil((east - x0) / 500)) + 1
    height = int(np.ceil((north - y0) / 500)) + 1
    raster = np.full((height, width), np.nan)
    cols = ((data.easting_m.to_numpy() - x0) // 500).astype(int)
    rows = ((data.northing_m.to_numpy() - y0) // 500).astype(int)
    inside = (cols >= 0) & (cols < raster.shape[1]) & (rows >= 0) & (rows < raster.shape[0])
    if not inside.all():
        raise ValueError("Map grid does not contain the plotted cells")
    raster[rows, cols] = val
    pm = ax.imshow(raster, origin="lower", extent=(0, width / 2, 0, height / 2),
                   cmap=cmap, norm=norm, interpolation="nearest")
    for polygon in getattr(boundary, "geoms", (boundary,)):
        if polygon.geom_type == "Polygon":
            x, y = np.asarray(polygon.exterior.coords).T
            ax.plot((x - x0) / 1000, (y - y0) / 1000, color="#4f615e", lw=0.36, alpha=0.75)
    ax.set_xlim((west - x0) / 1000 - 2, (east - x0) / 1000 + 2)
    ax.set_ylim((south - y0) / 1000 - 2, (north - y0) / 1000 + 2)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines[["left", "bottom"]].set_visible(False)
    return pm


def map_panels(folder, figure_dir):
    boundary, _ = read_boundary(BOUNDARY, GRID_CRS)
    cells = pd.read_csv(folder / "cell_availability.csv")
    fit = pd.read_csv(folder / "fitted_cells.csv")
    used = fit[fit.date == "2025-10-02"]
    eligible = cells[cells.model_eligible]
    fig, axs = plt.subplots(2, 2, figsize=(11.4, 8.2))
    fig.subplots_adjust(left=0.065, right=0.94, bottom=0.135, top=0.84, wspace=0.19, hspace=0.42)
    heading(fig, "Where Landsat measured surface temperature",
            "Singapore · 500 m cells · 7 selected daytime overpasses in 2025")
    thermal = matplotlib.colors.Normalize(*np.nanquantile(used.response_c, [0.02, 0.98]))
    counts = ListedColormap(["#d8d9d4", "#d9efea", "#bce1d9", "#8ecbbe", "#62ae9e",
                            "#3a9284", "#207a6d", "#14594f"])
    residual_scale = float(np.quantile(np.abs(used.residual_c), 0.98))
    maps = [
        (eligible, eligible.tree_fraction * 100, "Esri tree cover", "Tree share (%)", "YlGn",
         matplotlib.colors.Normalize(0, 100), [0, 25, 50, 75, 100]),
        (eligible, eligible.dates_used, "Usable overpasses per cell", "Number of dates", counts,
         BoundaryNorm(np.arange(-0.5, 8.5, 1), counts.N), list(range(8))),
        (used, used.response_c, "2 October · daytime surface temperature", "LST (°C)", "inferno",
         thermal, None),
        (used, used.residual_c, "2 October · residual after adjustment", "Observed minus fitted (°C)",
         "RdBu_r", TwoSlopeNorm(vmin=-residual_scale, vcenter=0, vmax=residual_scale), None),
    ]
    for ax, item in zip(axs.flat, maps):
        data, values, title, label, cmap, norm, ticks = item
        mapped = draw_map(ax, boundary, data, values, cmap, norm)
        ax.set_title(title, loc="left", fontsize=12, fontweight="bold", pad=2)
        bar = fig.colorbar(mapped, ax=ax, orientation="horizontal", fraction=0.061, pad=0.025,
                           ticks=ticks, extend="both" if ticks is None else "neither")
        bar.set_label(label, fontsize=9)
        bar.outline.set_visible(False)
        if ticks is None:
            bar.locator = MaxNLocator(nbins=6, integer=True)
            bar.update_ticks()
    scale_ax = axs[0, 0]
    left = scale_ax.get_xlim()[0] + 1.5
    bottom = scale_ax.get_ylim()[0] + 0.7
    scale_ax.plot([left, left + 5], [bottom, bottom], c=INK, lw=1.35)
    scale_ax.text(left + 2.5, bottom + 0.6, "5 km", ha="center", fontsize=8)
    fig.text(0.055, 0.055, "Temperature and residual color scales saturate at extreme values.",
             fontsize=8.5, color=MUTED)
    fig.text(0.945, 0.022, SOURCES, ha="right", fontsize=8, color=MUTED)
    fig.savefig(figure_dir / "landsat_maps.png", dpi=180)
    plt.close(fig)


def association(panels, folder, figure_dir):
    predictors = load_panel(panels, "world_500").drop_duplicates("cell_id").set_index("cell_id")
    points = pd.read_csv(folder / "fitted_cells.csv")
    points = points[points.date == "2025-10-02"].join(predictors[["tree_fraction"]], on="cell_id")
    result = json.loads((folder / "estimates.json").read_text())
    slopes = pd.DataFrame(result["primary_fit"]["slopes"])
    fig, axs = plt.subplots(1, 2, figsize=(11.4, 5.9), gridspec_kw={"width_ratios": (1.12, 1)})
    fig.subplots_adjust(left=0.085, right=0.94, bottom=0.225, top=0.745, wspace=0.37)
    heading(fig, "Tree cover and cooler daytime surfaces",
            "One Landsat overpass and the date-adjusted estimates across all seven")
    left, right = axs
    cells = left.hexbin(points.tree_fraction, points.response_c, gridsize=30, cmap="YlGnBu", mincnt=1,
                        bins="log", linewidths=0.1)
    col = pd.cut(points.tree_fraction, np.arange(0, 1.01, 0.1), include_lowest=True)
    means = points.groupby(col, observed=True)[["tree_fraction", "response_c"]].mean()
    left.plot(means.tree_fraction, means.response_c, marker="o", color="#a9552b",
              lw=1.7, ms=3.5, label="Mean within 10-point tree bins")
    left.xaxis.set_major_formatter(PercentFormatter(1))
    left.set_xlabel("Tree share of classified land")
    left.set_ylabel("500 m cell mean LST (°C)")
    left.set_title("2 October · observed cells", fontsize=12, loc="left", fontweight="bold", pad=10)
    left.legend(loc="upper right", frameon=False, fontsize=8)
    bar = fig.colorbar(cells, ax=left, shrink=0.8, pad=0.022)
    bar.set_label("Cell count", fontsize=9)
    bar.outline.set_visible(False)
    positions = np.arange(len(slopes))
    right.grid(axis="x", color=GRID, lw=0.6, zorder=0)
    right.axvline(0, color="#85958c", ls="--", lw=0.9)
    right.errorbar(slopes.estimate_c, positions,
                   xerr=[slopes.estimate_c - slopes.lower_c, slopes.upper_c - slopes.estimate_c],
                   fmt="o", color=GREEN, lw=2, capsize=3)
    right.set_yticks(positions, [f"{b:.0%} built" for b in slopes.built_fraction])
    right.set_ylim(-0.55, len(slopes) - 0.45)
    right.invert_yaxis()
    right.set_xlim(min(-1.12, slopes.lower_c.min() - 0.08), 0.08)
    right.set_xlabel("LST difference per +10 pp tree cover (°C)")
    right.set_title("Adjusted tree-cover slopes", fontsize=12, loc="left", fontweight="bold", pad=10)
    right.tick_params(axis="y", length=0)
    fig.text(0.055, 0.088, "Intervals use 1,999 spatial-block bootstrap draws, holding the seven dates fixed.",
             fontsize=8.5, color=MUTED)
    fig.text(0.945, 0.027, SOURCES, ha="right", fontsize=8, color=MUTED)
    fig.savefig(figure_dir / "tree_cover_association.png", dpi=180)
    plt.close(fig)


def date_estimates(folder, figure_dir):
    data = pd.read_csv(folder / "date_estimates.csv")
    omitted = pd.read_csv(folder / "leave_one_date_out.csv")
    overall = json.loads((folder / "estimates.json").read_text())["primary_fit"]
    tree = next(r for r in overall["slopes"] if r["built_fraction"] == 0.50)
    fig, axs = plt.subplots(1, 2, figsize=(11.4, 5.7), sharey=True)
    fig.subplots_adjust(left=0.155, right=0.95, bottom=0.20, top=0.75, wspace=0.27)
    heading(fig, "Association by Landsat acquisition",
            "Estimates at 50% built coverage, with each date fitted separately or omitted in turn")
    y = np.arange(len(data))
    for ax, val, lo, hi, reference, title, xlabel in [
        (axs[0], "tree_at_50pct_c", "tree_lower_c", "tree_upper_c", tree["estimate_c"],
         "Tree-cover slope", "°C per +10 pp tree cover"),
        (axs[1], "interaction_c", "interaction_lower_c", "interaction_upper_c",
         overall["interaction"]["estimate_c"], "Tree × built interaction", "Slope change per +10 pp built (°C)"),
    ]:
        ax.axvline(0, color="#a7b6af", lw=0.8)
        ax.axvline(reference, color="#a9552b", ls="--", lw=1.0, label="Pooled estimate")
        ax.grid(axis="x", color=GRID, lw=0.6)
        ax.errorbar(data[val], y, xerr=[data[val] - data[lo], data[hi] - data[val]],
                    fmt="o", color=GREEN, lw=1.5, ms=4, capsize=2, label="Date-specific estimate")
        ax.scatter(omitted[val], y + 0.19, marker="d", facecolors="none", edgecolors="#546874",
                   s=26, label="Omit this date")
        ax.set_title(title, loc="left", fontsize=12, fontweight="bold", pad=10)
        ax.set_xlabel(xlabel)
    axs[0].set_yticks(y, pd.to_datetime(data.date).dt.strftime("%d %b"))
    axs[0].set_ylim(len(data) - 0.5, -0.5)
    handles, labels = axs[1].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", bbox_to_anchor=(0.53, 0.080), ncols=3,
               frameon=False, fontsize=8.5)
    fig.text(0.055, 0.035, "Date-specific bars show 95% spatial-bootstrap intervals.",
             fontsize=8.5, color=MUTED)
    fig.savefig(figure_dir / "date_estimates.png", dpi=180)
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--panels", required=True, type=Path)
    parser.add_argument("--results", required=True, type=Path)
    parser.add_argument("--figures", required=True, type=Path)
    args = parser.parse_args()
    args.figures.mkdir(parents=True, exist_ok=True)
    style()
    map_panels(args.results, args.figures)
    association(args.panels, args.results, args.figures)
    date_estimates(args.results, args.figures)
