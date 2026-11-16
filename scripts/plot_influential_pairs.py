"""Draw the original imagery, temperature-quality review, and omission contrast."""
import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap, Normalize
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np
import pandas as pd
from rasterio.transform import Affine, array_bounds

from analyze_influential_pairs import SR_OFFSET, SR_SCALE, masks
from make_influential_pair_chips import PADDING, SENTINEL_BANDS, extent
from plot_landsat import GREEN, GRID, INK, LIGHT, MUTED, heading, style
from prepare_influential_pair_audit import DESIGN, SNAPSHOT


BLUE = "#227c9b"
GOLD = "#ad7830"
COVER_COLORS = ["#f2f0e7", "#88bbc9", "#18694d", "#c6ae59", "#a99c9f", "#ceb997"]
COVER_LABELS = ["Outside study land", "Mapped water", "Trees", "Rangeland",
                "Built area", "Other / unclassified"]


def reflectance_image(rgb):
    return np.clip((rgb - .012) / .29, 0, 1) ** .63


def site_context(snapshot, sites, dest):
    s2_meta = pd.read_csv(snapshot / "sentinel_sources.csv").set_index("site")
    fig, axes = plt.subplots(len(sites), 4, figsize=(11.7, 8.6))
    fig.subplots_adjust(left=.17, right=.952, top=.81, bottom=.127, hspace=.28, wspace=.12)
    heading(fig, "Original 10 m imagery for the three wider-only pairs",
            "Sentinel-2 near the reference overpass and the annual Esri map · each square spans 500 m")
    labels = ("Lower tree · Sentinel-2", "Higher tree · Sentinel-2", "Lower tree · Esri", "Higher tree · Esri")
    for ax, name in zip(axes[0], labels):
        ax.set_title(name, loc="left", fontsize=10, fontweight="bold", pad=8)
    for row, site in zip(axes, sites.itertuples(index=False)):
        with np.load(snapshot / f"cover_{site.site}.npz", allow_pickle=False) as land:
            with np.load(snapshot / f"sentinel_{site.site}.npz", allow_pickle=False) as s2:
                offsets = dict(zip(SENTINEL_BANDS, s2["offsets"]))
                rgb = np.moveaxis(np.stack([
                    (s2[key][:, PADDING:-PADDING, PADDING:-PADDING] + offsets[key]) / s2["quantification"]
                    for key in ("B04", "B03", "B02")
                ]), 0, -1)
                for i in range(2):
                    row[i].imshow(reflectance_image(rgb[i]), interpolation="nearest")
                    original = land["classes"][i]
                    classified = np.full(original.shape, 5, dtype=np.uint8)
                    for label, code in ((1, 1), (2, 2), (3, 11), (4, 7)):
                        classified[original == code] = label
                    classified[~land["land"][i]] = 0
                    row[i + 2].imshow(classified, cmap=ListedColormap(COVER_COLORS),
                                      norm=BoundaryNorm(np.arange(-.5, 6), 6), interpolation="nearest")
                for ax in row:
                    ax.set_axis_off()
                ref = s2_meta.loc[site.site, "datetime_utc"][:10]
                row[0].text(-.13, .51, f"Pair {site.site}\n{ref}", transform=row[0].transAxes,
                            va="center", ha="right", color=INK, fontsize=10, fontweight="bold")
    fig.legend(handles=[Patch(fc=color, label=label) for color, label in zip(COVER_COLORS, COVER_LABELS)],
               loc="lower center", bbox_to_anchor=(.55, .067), frameon=False, ncols=6, fontsize=8.7)
    fig.text(.95, .028, "Sources  Copernicus Sentinel-2 · Esri / Impact Observatory · ESA WorldCover · URA",
             ha="right", fontsize=8, color=MUTED)
    fig.savefig(dest, dpi=170)
    plt.close(fig)


def show_native(axis, values, transform, location, **kwargs):
    x0, y0, x1, y1 = extent([location])
    h, w = values.shape[:2]
    left, bottom, right, top = array_bounds(h, w, Affine(*transform))
    axis.imshow(values, extent=(left - x0, right - x0, bottom - y0, top - y0),
                interpolation="nearest", **kwargs)
    axis.set_xlim(0, x1 - x0)
    axis.set_ylim(0, y1 - y0)
    axis.set_xticks([])
    axis.set_yticks([])
    for spine in axis.spines.values():
        spine.set_linewidth(.42)
        spine.set_visible(True)


def overpass_context(snapshot, cases, analysis, dest):
    values = pd.read_csv(analysis / "pair_quality.csv").set_index(["site", "date"])
    ordered = cases.sort_values(["site", "date"])
    fig, axes = plt.subplots(len(ordered), 4, figsize=(9.65, 16.8))
    fig.subplots_adjust(left=.19, right=.81, top=.884, bottom=.079, hspace=.13, wspace=.14)
    fig.text(.073, .966, "Same-overpass Landsat imagery and temperature", fontsize=20,
             fontweight="bold", va="top")
    fig.text(.073, .936, "All nine wider-only observations · 30 m delivery pixels within 500 m cells",
             fontsize=10.3, color=MUTED, va="top")
    titles = ("Lower tree · RGB", "Higher tree · RGB", "Lower tree · LST", "Higher tree · LST")
    for ax, label in zip(axes[0], titles):
        ax.set_title(label, loc="left", fontsize=9, fontweight="bold", pad=8)
    cmap = plt.get_cmap("inferno").copy()
    cmap.set_bad("#e8e4d8")
    common = Normalize(27, 49, clip=True)
    for row, pair in zip(axes, ordered.itertuples(index=False)):
        with np.load(snapshot / f"landsat_{pair.date}.npz", allow_pickle=False) as data:
            rgb = np.stack([data[key] for key in ("red", "green", "blue")], axis=-1)
            rgb = reflectance_image(rgb * SR_SCALE + SR_OFFSET)
            qa, celsius = masks(data["lwir11"], data["qa_pixel"], data["cdist"], data["qa"])
            valid_lst = np.where(qa[0], celsius, np.nan)
            cells = (pair.lower_tree_cell, pair.higher_tree_cell)
            for i, cell in enumerate(cells):
                show_native(row[i], rgb, data["transform"], cell)
                show_native(row[i + 2], valid_lst, data["transform"], cell, cmap=cmap, norm=common)
        row[0].text(-.10, .5, f"{pair.site}  {pair.date}", transform=row[0].transAxes,
                    ha="right", va="center", fontsize=8.3)
        stats = values.loc[(pair.site, pair.date)]
        row[-1].text(1.13, .51, f"{stats.base_gap_c:+.2f} °C", transform=row[-1].transAxes,
                     ha="left", va="center", fontsize=8.8, color=INK)
    cax = fig.add_axes([.246, .047, .513, .01])
    cb = fig.colorbar(plt.cm.ScalarMappable(cmap=cmap, norm=common), cax=cax,
                      orientation="horizontal", extend="both", ticks=[27, 33, 39, 45, 49])
    cb.set_label("Daytime LST (°C) · labels at right give the observed higher-minus-lower gap", fontsize=8.6)
    cb.outline.set_visible(False)
    fig.savefig(dest, dpi=165)
    plt.close(fig)


def summary_figure(analysis, dest):
    cases = pd.read_csv(analysis / "pair_quality.csv").sort_values(["site", "date"])
    effects = pd.read_csv(analysis / "block_omission.csv")
    fig = plt.figure(figsize=(12.0, 8.2))
    spec = fig.add_gridspec(2, 1, left=.111, right=.94, top=.835, bottom=.109,
                           height_ratios=[1, .96], hspace=.59)
    heading(fig, "Temperature quality and sensitivity to the two blocks",
            "All differences are higher-tree minus lower-tree surface temperature, standardized to five points")
    ax = fig.add_subplot(spec[0])
    x = np.arange(len(cases))
    for name, color, symbol, offset, label in (
        ("base", MUTED, "o", -.17, "Original QA"), ("qa4", BLUE, "s", 0, "≤4 K and ≥0.5 km from cloud"),
        ("qa2", GREEN, "^", .17, "≤2 K and ≥0.5 km from cloud"),
    ):
        ax.scatter(x + offset, cases[f"{name}_per5_c"], label=label, color=color, marker=symbol, s=32, zorder=3)
    ax.set_title("Wider-only pairs under three temperature screens", loc="left", fontweight="bold")
    ax.set_xticks(x, [f"{s}\n{d}" for s, d in zip(cases.site, cases.date)], fontsize=8.1)
    ax.set_ylabel("Temperature gap (°C)")
    ax.axhline(0, color="#93a49b", lw=.8)
    ax.grid(axis="y", color=GRID, lw=.6, alpha=.85)
    ax.legend(loc="upper center", bbox_to_anchor=(.5, -.24), frameon=False, fontsize=8.8, ncols=3)
    ax = fig.add_subplot(spec[1])
    show = effects[effects.model.eq("matched")]
    ax.axvline(0, color="#93a49b", lw=.9)
    for i, (key, label) in enumerate((("dense", "Primary matching"), ("wider_balance", "Wider matching"))):
        for s, color, height in (("full", MUTED, i-.145), ("omit_both", GREEN, i+.145)):
            row = show[show.case.eq(key) & show["sample"].eq(s)].iloc[0]
            ax.hlines(height, row.lower_c, row.upper_c, color=color, lw=1.65)
            ax.scatter(row.estimate_c_per_5pp, height, color=color, s=32, zorder=3)
    ax.set_ylim(1.48, -.48)
    ax.set_yticks([0, 1], ["Primary matching", "Wider matching"], fontsize=9)
    ax.set_xlabel("Matched LST difference (°C per five points)")
    ax.set_title("Both blocks omitted from the frozen lists", loc="left", fontweight="bold")
    ax.grid(axis="x", color=GRID, lw=.65, alpha=.85)
    ax.legend(handles=[
        Line2D([], [], color=MUTED, marker="o", label="Full sample"),
        Line2D([], [], color=GREEN, marker="o", label="Omit both blocks"),
    ], loc="upper right", frameon=False, fontsize=8.8)
    fig.text(.111, .038, "Intervals are 95% geographic block resampling intervals. Each quality screen needs ≥60% land coverage at both endpoints.",
             fontsize=8.4, color=MUTED)
    fig.savefig(dest, dpi=175)
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--snapshot", type=Path, default=SNAPSHOT)
    parser.add_argument("--design", type=Path, default=DESIGN)
    parser.add_argument("--analysis", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--audit-only", action="store_true")
    args = parser.parse_args()
    style()
    args.output.mkdir(parents=True, exist_ok=True)
    site_context(args.snapshot, pd.read_csv(args.design / "review_sites.csv"), args.output / "sentinel_site_context.png")
    overpass_context(args.snapshot, pd.read_csv(args.design / "pair_registry.csv"), args.analysis,
                     args.output / "landsat_original_pairs.png")
    if not args.audit_only:
        summary_figure(args.analysis, args.output / "temperature_quality_and_omission.png")
