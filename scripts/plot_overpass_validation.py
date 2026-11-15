"""Draw validation coverage, annual cover checks, and temporal comparisons."""
import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap, Normalize
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle
import numpy as np
import pandas as pd

from analyze_overpass_validation import MODERATION, PAIRS, REGRESSION
from audit_coverage import BOUNDARY, read_boundary
from audit_heldout_overpasses import AUDIT, SEASONS
from plot_built_results import draw_interval
from plot_landsat import GREEN, GRID, INK, MUTED, SOURCES, draw_map, heading, style
from prepare_overpass_validation import DESIGN, YEARS
from screen_landsat import GRID_CRS


BLUE = "#3f638c"
GOLD = "#af682c"
GREY = "#a8b5b0"
COHORT_COLOR = {"validation": BLUE, "discovery": GREEN}


def save_figure(fig, path):
    fig.text(0.95, 0.025, SOURCES, ha="right", fontsize=8, color=MUTED)
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_candidates(audit, path):
    fig = plt.figure(figsize=(11.8, 7.35))
    spec = fig.add_gridspec(2, 3, left=0.076, right=0.96, bottom=0.11, top=0.84,
                           hspace=0.58, wspace=0.11, height_ratios=[1.3, 0.86])
    heading(fig, "Landsat coverage for additional overpasses",
            "Singapore · Tier 1 scenes · fixed land and planning-area thresholds")
    period_counts = []
    for i, year in enumerate(YEARS):
        ax = fig.add_subplot(spec[0, i])
        data = pd.read_csv(audit / str(year) / "scene_coverage.csv")
        passes = data.passes_qa_and_range_60pct_44_areas.eq(1)
        cohort = "discovery" if year == 2025 else "validation"
        ax.add_patch(Rectangle((60, 44), 40, 14, facecolor="#e4eee6", edgecolor="none", zorder=0))
        for mask, color, alpha in [(~passes, GREY, 0.72), (passes, COHORT_COLOR[cohort], 1)]:
            rows = data[mask]
            ax.scatter(rows.qa_and_range_worldcover_land_pct, rows.qa_and_range_planning_areas_at_30pct,
                       color=color, s=23, alpha=alpha, linewidths=0, zorder=3)
        ax.axvline(60, color="#7d9185", lw=0.85, ls="--")
        ax.axhline(44, color="#7d9185", lw=0.85, ls="--")
        ax.set(xlim=(-3, 100), ylim=(-2, 58), xticks=[0, 20, 40, 60, 80, 100], yticks=[0, 15, 30, 44, 55])
        ax.set_title(f"{year}  ·  {int(passes.sum())} of {len(data)} qualify", loc="left",
                     fontsize=11.5, fontweight="bold", pad=9)
        if i == 0:
            ax.set_ylabel("Planning areas with ≥30% usable land")
        else:
            ax.tick_params(axis="y", labelleft=False)
        if i == 1:
            ax.set_xlabel("Usable WorldCover land (%)")
        period_counts.append(data[passes].groupby("period").size().reindex(SEASONS, fill_value=0))

    ax = fig.add_subplot(spec[1, :])
    positions = np.arange(len(SEASONS))
    validation = pd.concat(period_counts[:-1], axis=1).sum(axis=1)
    for offset, counts, cohort in [(-0.18, validation, "validation"), (0.18, period_counts[-1], "discovery")]:
        ax.bar(positions + offset, counts, width=0.31, color=COHORT_COLOR[cohort])
    ax.set_xticks(positions, ["Dec–Mar\nNortheast monsoon", "Apr–May\nFirst inter-monsoon",
                             "Jun–Sep\nSouthwest monsoon", "Oct–Nov\nSecond inter-monsoon"])
    ax.set_ylabel("Overpasses")
    ax.set_title("Qualifying dates by calendar period", loc="left", fontweight="bold", pad=6)
    ax.set_ylim(0, max(validation.max(), period_counts[-1].max()) + 0.9)
    ax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))
    ax.grid(axis="y", color=GRID, linewidth=0.65)
    ax.set_axisbelow(True)
    ax.legend(handles=[Patch(color=BLUE, label="2023–24 validation"), Patch(color=GREEN, label="2025 discovery")],
              frameon=False, ncols=2, loc="upper right", fontsize=9)
    save_figure(fig, path)


def plot_geography(design, path):
    annual = pd.read_csv(design / "annual_cover_comparison.csv")
    boundary, _ = read_boundary(BOUNDARY, GRID_CRS)
    dates = annual.filter(regex=r"^usable_dates_\d{4}$")
    available = annual[dates.gt(0).any(axis=1)].copy()
    category = np.select([available.common_stable, available.common], [3, 2], default=1)
    fig, axes = plt.subplots(2, 2, figsize=(11.5, 8))
    fig.subplots_adjust(left=0.065, right=0.94, top=0.85, bottom=0.15, hspace=0.29, wspace=0.16)
    heading(fig, "Common geography and annual land cover",
            "Singapore · 500 m cells · share ranges across annual Esri classifications")
    colors = ListedColormap(["#d6e0e7", "#d8a85d", GREEN])
    left = axes[0, 0]
    draw_map(left, boundary, available, category, colors, BoundaryNorm([0.5, 1.5, 2.5, 3.5], 3))
    left.set_title("Geography of eligible observations", loc="left", fontweight="bold", fontsize=11)
    left.legend(handles=[
        Patch(color=GREEN, label="Common, each cover range ≤5 pp"),
        Patch(color="#d8a85d", label="Common, at least one range >5 pp"),
        Patch(color="#d6e0e7", label="Usable in some years"),
    ], frameon=False, loc="lower left", bbox_to_anchor=(0.0, -0.16), fontsize=8)
    common = annual[annual.common]
    for ax, col, title in (
        (axes[0, 1], "tree_fraction_span_pp", "Tree cover range"),
        (axes[1, 0], "built_fraction_span_pp", "Built cover range"),
        (axes[1, 1], "esri_water_share_span_pp", "Esri water share range"),
    ):
        im = draw_map(ax, boundary, common, common[col], "YlOrRd", Normalize(0, 20))
        ax.set_title(title, loc="left", fontweight="bold", fontsize=11)
    cax = fig.add_axes([0.39, 0.091, 0.22, 0.016])
    bar = fig.colorbar(im, cax=cax, orientation="horizontal", extend="max", ticks=[0, 5, 10, 15, 20])
    bar.outline.set_visible(False)
    bar.set_label("Range of shares in common cells (percentage points)", fontsize=9)
    save_figure(fig, path)


def forest_axis(ax, xlabel, title):
    ax.axvline(0, color="#879991", lw=0.9, ls="--", zorder=1)
    ax.grid(axis="x", color=GRID, lw=0.7)
    ax.set_axisbelow(True)
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.set_xlabel(xlabel)
    ax.set_title(title, loc="left", pad=8, fontweight="bold", fontsize=11.5)


def plot_pooled(design, results, path):
    record = json.loads((design / "design_record.json").read_text())
    validation_years = "–".join(map(str, record["validation_years"]))
    local = pd.read_csv(results / "local_estimates.csv")
    built = pd.read_csv(results / "built_cover_estimates.csv")
    sample_labels = {"full": "All cells", "common": "Common cells", "common_stable": "Common stable cells"}
    options = [(cohort, sample) for sample in sample_labels for cohort in ("discovery", "validation")]
    names = {"discovery": "2025", "validation": validation_years}
    fig, (left, right) = plt.subplots(1, 2, figsize=(12.15, 6.8))
    fig.subplots_adjust(left=0.153, right=0.959, top=0.79, bottom=0.222, wspace=0.60)
    heading(fig, "Local temperature contrasts on additional dates",
            "Singapore · same 500 m grid and QA rules · discovery and temporal validation")
    forest_axis(left, "Daytime LST difference per +10 tree points (°C)", "Stage 1 · local tree association")
    forest_axis(right, "Difference in tree contrast, 85% minus 65% built (°C)", "Stage 2 · variation with built cover")
    labels = []
    for y, (cohort, sample) in enumerate(options):
        color = COHORT_COLOR[cohort]
        for method, offset, shade in [(REGRESSION, -0.125, color), (PAIRS, 0.125, GOLD)]:
            row = local[local.cohort.eq(cohort) & local["sample"].eq(sample) & local.analysis.eq(method)].iloc[0]
            if row.status == "estimated":
                draw_interval(left, y + offset, row.estimate_c_per_10pp, row.lower_c, row.upper_c, shade)
        row = built[built.cohort.eq(cohort) & built["sample"].eq(sample)
                    & built.analysis.eq("Primary quadratic interaction")].iloc[0]
        if row.status == "estimated":
            draw_interval(right, y, row.higher_minus_lower_estimate_c, row.higher_minus_lower_lower_c,
                          row.higher_minus_lower_upper_c, color)
        else:
            right.annotate("Sparse support", (0.03, y), xycoords=("axes fraction", "data"),
                           va="center", color=MUTED, fontsize=9)
        labels.append(f"{names[cohort]}\n{sample_labels[sample]}")
    for axis in (left, right):
        axis.set_yticks(range(len(labels)), labels)
        axis.set_ylim(len(labels) - 0.47, -0.53)
    fig.legend(handles=[
        Line2D([0], [0], marker="o", color=GREEN, label="2025 discovery regression"),
        Line2D([0], [0], marker="o", color=BLUE, label=f"{validation_years} validation regression"),
        Line2D([0], [0], marker="o", color=GOLD, label="Nearby pairs"),
    ], ncols=3, frameon=False, loc="lower center", bbox_to_anchor=(0.52, 0.092), fontsize=9)
    fig.text(0.053, 0.065, "Bars show approximate 95% spatial-bootstrap intervals conditional on the observed dates.",
             color=MUTED, fontsize=9)
    save_figure(fig, path)


def plot_by_date(results, path):
    daily = pd.read_csv(results / "date_estimates.csv")
    dates = sorted(daily.date.unique())
    fig, (left, right) = plt.subplots(1, 2, figsize=(11.75, max(6.6, 3.0 + 0.29 * len(dates))), sharey=True)
    fig.subplots_adjust(left=0.106, right=0.951, top=0.80, bottom=0.16, wspace=0.20)
    heading(fig, "Temperature contrasts by Landsat overpass",
            "Points summarize date-to-date variation · bars quantify spatial uncertainty within each date")
    for ax, method, xlabel, title in (
        (left, REGRESSION, "LST difference per +10 tree points (°C)", "Stage 1 · local regression"),
        (right, MODERATION, "85% minus 65% built (°C)", "Stage 2 · contrast difference"),
    ):
        forest_axis(ax, xlabel, title)
        rows = daily[daily.analysis == method].set_index("date").loc[dates]
        for y, row in enumerate(rows.itertuples()):
            if pd.notna(row.estimate_c):
                if pd.notna(row.lower_c):
                    draw_interval(ax, y, row.estimate_c, row.lower_c, row.upper_c, COHORT_COLOR[row.cohort])
                else:
                    ax.scatter(row.estimate_c, y, color=COHORT_COLOR[row.cohort], s=22)
            else:
                ax.annotate("Sparse support", (0.05, y), xycoords=("axes fraction", "data"),
                            va="center", color=MUTED, fontsize=8)
        ax.set_ylim(len(dates) - 0.46, -0.54)
    left.set_yticks(range(len(dates)), pd.to_datetime(dates).strftime("%d %b %Y"))
    fig.legend(handles=[Line2D([0], [0], marker="o", color=c, label=s)
                        for c, s in [(BLUE, "Validation"), (GREEN, "Discovery")]],
               ncols=2, frameon=False, loc="lower center", bbox_to_anchor=(0.5, 0.039), fontsize=9)
    save_figure(fig, path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--audit", type=Path, default=AUDIT)
    parser.add_argument("--design", type=Path, default=DESIGN)
    parser.add_argument("--results", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    style()
    plot_candidates(args.audit, args.output / "candidate_coverage.png")
    plot_geography(args.design, args.output / "common_geography.png")
    if args.results:
        plot_pooled(args.design, args.results, args.output / "validation_estimates.png")
        plot_by_date(args.results, args.output / "by_overpass.png")
