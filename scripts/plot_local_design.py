"""Map the coverage and matched cells from the outcome-free Stage 1 design."""
import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap
import numpy as np
import pandas as pd

from audit_coverage import BOUNDARY, read_boundary
from plot_landsat import GREEN, GRID, INK, MUTED, SOURCES, draw_map, heading, style
from screen_landsat import GRID_CRS


def plot_design(design, output):
    cells = pd.read_csv(design / "cell_design.csv")
    dates = pd.read_csv(design / "design_by_date.csv").query("design == 'primary'")
    available = cells[cells.usable_dates > 0]
    boundary, _ = read_boundary(BOUNDARY, GRID_CRS)
    colors = ListedColormap(["#e4e5df", "#b6ded2", "#8ac8ba", "#60ada0", "#3d9184",
                             "#287c70", "#17675e", "#084c45"])
    norm = BoundaryNorm(np.arange(-0.5, 8.5), colors.N)
    style()
    fig = plt.figure(figsize=(11, 8.2))
    grid = fig.add_gridspec(2, 2, left=0.075, right=0.95, top=0.84, bottom=0.125,
                            hspace=0.35, wspace=0.17, height_ratios=[1.2, 1])
    heading(fig, "Where nearby grid cells can be compared",
            "Singapore · 500 m cells · seven Landsat overpasses in 2025")
    for axis, key, title in (
        (fig.add_subplot(grid[0, 0]), "usable_dates", "Usable observations"),
        (fig.add_subplot(grid[0, 1]), "paired_dates", "Same-date matched cells"),
    ):
        artist = draw_map(axis, boundary, available, available[key], colors, norm)
        axis.set_title(title, loc="left", fontweight="bold", pad=5)
        bar = fig.colorbar(artist, ax=axis, orientation="horizontal", pad=0.035, fraction=0.09,
                           ticks=np.arange(8))
        bar.set_label("Number of dates", fontsize=9)
        bar.outline.set_visible(False)
        if key == "usable_dates":
            x = axis.get_xlim()[0] + 1.5
            y = axis.get_ylim()[0] + 1
            axis.plot([x, x + 5], [y, y], color=INK, lw=1.4)
            axis.text(x + 2.5, y + 0.7, "5 km", ha="center", fontsize=8)
    ax = fig.add_subplot(grid[1, :])
    x = np.arange(len(dates))
    ax.bar(x, dates.candidate_cells, color="#d8e9e1", width=0.61, label="Cells with a candidate match")
    ax.bar(x, dates.matched_cells, color=GREEN, width=0.38, label="Cells in disjoint pairs")
    ax.set_xticks(x, pd.to_datetime(dates.date).dt.strftime("%d %b"))
    ax.set_ylabel("Number of cells")
    ax.set_title("Candidate and matched cells by date", loc="left", fontweight="bold", pad=7)
    ax.grid(axis="y", color=GRID, alpha=0.8, linewidth=0.65)
    ax.set_axisbelow(True)
    ax.set_ylim(0, dates.candidate_cells.max() * 1.20)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.11), ncols=2, frameon=False, fontsize=9)
    fig.text(0.95, 0.027, SOURCES, ha="right", fontsize=8, color=MUTED)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=180)
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--design", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    plot_design(args.design, args.output)
