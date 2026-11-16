"""Plot paired temperatures, weights, and one-deletion estimates."""
import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import PercentFormatter
import numpy as np
import pandas as pd

from plot_landsat import GREEN, GRID, INK, MUTED, heading, style


GOLD = "#b1762c"
BLUE = "#267c9b"
ZERO = "#94a39f"


def weighted_pairs_figure(design, path):
    table = pd.read_csv(design / "weighted_pairs.csv")
    fig, axes = plt.subplots(1, 3, figsize=(12.2, 6.4), sharey=True,
                             gridspec_kw={"width_ratios": [1.30, .80, 1]})
    fig.subplots_adjust(left=.198, right=.953, top=.78, bottom=.21, wspace=.29)
    heading(fig, "Temperature and weight of each retained pair",
            "Nine pairs from the Stage 8 4 K comparison after both blocks were omitted")
    i = np.arange(len(table))
    names = [f"{pair_id} · {date}" for pair_id, date in zip(table.pair_id, table.date)]
    axes[0].set_yticks(i, names, fontsize=9)
    axes[0].set_ylim(len(i) - .5, -.5)
    for ax in axes:
        ax.axvline(0, lw=.85, color=ZERO)
        ax.tick_params(axis="y", length=0)
        ax.spines["left"].set_visible(False)
        ax.grid(axis="x", color=GRID, lw=.6, alpha=.8)
        ax.set_axisbelow(True)

    ax = axes[0]
    ax.scatter(table.lst_gap_c, i - .13, s=25, marker="x", color=MUTED, linewidths=1.1)
    ax.scatter(table.lst_c_per_5pp, i + .11, s=28, color=GREEN, linewidths=0)
    ax.set_title("Temperature gap", loc="left", fontsize=11, fontweight="bold")
    ax.set_xlabel("Daytime LST difference (°C)")
    ax.set_xlim(-3.64, 3.03)

    ax = axes[1]
    weight_bars = ax.barh(i, table.weight, color=BLUE, height=.46)
    ax.bar_label(weight_bars, [f"{w:.1%}" for w in table.weight], fontsize=8.5, padding=3, color=INK)
    ax.set_title("Analysis weight", loc="left", fontsize=11, fontweight="bold")
    ax.set_xlabel("Share of total weight")
    ax.set_xlim(0, .23)
    ax.set_xticks([0, .05, .10, .15, .20])
    ax.xaxis.set_major_formatter(PercentFormatter(xmax=1, decimals=0))

    ax = axes[2]
    contribution = table.weighted_contribution_c
    signed = np.where(contribution < 0, GREEN, GOLD)
    ax.barh(i, contribution, color=signed, height=.46)
    ax.set_title("Contribution to estimate", loc="left", fontsize=11, fontweight="bold")
    ax.set_xlabel("Weight × five-point gap (°C)")

    handles = [
        Line2D([], [], linestyle="", marker="x", color=MUTED, label="Actual higher-minus-lower gap"),
        Line2D([], [], linestyle="", marker="o", color=GREEN, label="Gap scaled to five points", markersize=5),
    ]
    fig.legend(handles=handles, loc="lower left", bbox_to_anchor=(.18, .087), ncols=2, frameon=False, fontsize=8.9)
    fig.text(.198, .049, "Each overpass receives 1/6 of the total weight. The signed contributions sum to −0.85 °C.",
             fontsize=8.7, color=MUTED)
    fig.savefig(path, dpi=175)
    plt.close(fig)


def omission_figure(analysis, path):
    table = pd.read_csv(analysis / "leave_one_out.csv")
    baseline = table.loc[table.omission_type.eq("none")].iloc[0]
    fig, axes = plt.subplots(1, 2, figsize=(12.2, 6.45), sharex=True, sharey=True)
    fig.subplots_adjust(left=.173, right=.96, top=.783, bottom=.224, wspace=.45)
    heading(fig, "How the five-point estimate changes with one deletion",
            "One overpass or one 4 km block omitted at a time from the same nine-pair sample")
    for ax, kind, name in zip(axes, ("date", "block"), ("Omit one overpass", "Omit one block")):
        group = table.loc[table.omission_type.eq(kind)]
        if kind == "date":
            labels = list(group.omitted_date)
        else:
            labels = [f"({int(e)}, {int(n)})" for e, n in zip(group.omitted_block_e, group.omitted_block_n)]
        labels = [f"{label}\n{p} pairs · {b} blocks" for label, p, b in zip(labels, group.pairs, group.blocks)]
        i = np.arange(len(group))
        ax.axvline(0, lw=.85, color=ZERO)
        ax.axvline(baseline.estimate_c_per_5pp, lw=1.0, color=BLUE, linestyle="--")
        for height, row in enumerate(group.itertuples(index=False)):
            if row.status != "estimated":
                ax.text(.03, height, "Sparse support", transform=ax.get_yaxis_transform(), color=MUTED)
                continue
            color = GOLD if row.interval_includes_zero else GREEN
            ax.hlines(height, row.lower_c, row.upper_c, lw=1.7, color=color)
            ax.scatter(row.estimate_c_per_5pp, height, color=color, s=30, linewidths=0, zorder=3)
        # Separate tick labels across the two panels while preserving row positions.
        ax.set_yticks(i)
        ax.tick_params(axis="y", labelleft=False, length=0)
        for height, label in zip(i, labels):
            ax.text(-.033, height, label, transform=ax.get_yaxis_transform(), ha="right",
                    va="center", fontsize=8.6, color=MUTED, linespacing=1.28)
        ax.set_ylim(len(i) - .52, -.48)
        ax.spines["left"].set_visible(False)
        ax.set_title(name, loc="left", fontweight="bold", fontsize=11.6)
        ax.set_xlabel("Matched LST difference (°C per five points)")
        ax.grid(axis="x", color=GRID, lw=.6, alpha=.84)
        ax.set_axisbelow(True)
    handles = [
        Line2D([], [], color=BLUE, linestyle="--", label="Full sample point estimate"),
        Line2D([], [], color=GREEN, marker="o", markersize=4, label="Interval below zero"),
        Line2D([], [], color=GOLD, marker="o", markersize=4, label="Interval includes zero"),
    ]
    fig.legend(handles=handles, loc="lower left", bbox_to_anchor=(.154, .083), ncols=3, frameon=False, fontsize=8.7)
    fig.text(.173, .047, "95% spatial block resampling intervals. Weights give equal total weight to each remaining date.",
             fontsize=8.7, color=MUTED)
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
    weighted_pairs_figure(args.design, args.output / "pair_weights_and_contributions.png")
    omission_figure(args.analysis, args.output / "leave_one_out.png")
