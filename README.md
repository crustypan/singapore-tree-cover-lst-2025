# Tree cover and daytime land surface temperature in Singapore

This study uses seven clear-sky Landsat 8/9 overpasses in 2025 to ask whether parts of Singapore with more tree cover had cooler daytime surfaces, and whether the association varied with built-area cover. It combines 10 m land-cover classifications and Landsat surface-temperature measurements on a 500 m grid.

## Stage 1 local comparisons

Within the same Landsat date and 4 km area, a 10-percentage-point higher tree share was associated with **0.48 °C lower daytime land surface temperature** after accounting for built cover and land share. The approximate 95% spatial-bootstrap interval was **0.40–0.57 °C lower**. The regression used 14,641 cell-dates and 2,486 distinct 500 m cells. In 431 disjoint same-date pairs with comparable built and land shares, the estimated difference was 0.74 °C lower, with an interval of 0.42–1.06 °C lower. The paired cells cover a smaller and greener part of Singapore's usable sample.

[Stage 1 results and diagnostics](docs/LOCAL_COMPARISON_RESULTS.md) · [Frozen pairing protocol](docs/LOCAL_COMPARISON_PROTOCOL.md) · [Map of comparable cells](results/local_comparison/design/comparable_grid_cells.png)

## Baseline results

The main analysis uses **14,678 observations across 2,491 grid cells**. At 50% built-area cover, cells with 10 percentage points more tree cover had **0.73 °C lower daytime land surface temperature** on average within an overpass, with an approximate 95% spatial-bootstrap interval of **0.61 to 0.82 °C lower**.

The estimated tree association was stronger at higher built-area shares in the linear interaction model. Adding separate quadratic terms for tree and built-area cover made the interaction less certain. The interval included zero in the quadratic model. The negative tree association held across the dates and coverage, mask, grid and model checks.

![Observed tree cover and daytime LST on 2 October alongside the model-adjusted slopes across seven overpasses](results/landsat_analysis/figures/tree_cover_association.png)

The observations describe land surface temperature around **11:16 a.m. Singapore time** on the seven selected dates. The estimates are observational and do not establish causal cooling. They describe satellite surface temperature for the sampled overpasses. Air temperature, night conditions, and individual heat exposure require distinct measurements.

## How the sample was built

1. The feasibility study tested 2025 Terra MODIS coverage against existing annual rules. Even assuming every missing source date passed QA, at most 7.70% of land could meet the temporal rules against an 80% requirement. [Feasibility report](docs/FEASIBILITY.md)
2. The Landsat audit screened 44 Collection 2 scenes. Seven Tier 1 acquisitions passed the per-scene coverage and temperature checks. The study uses those seven acquisitions. [Scene-level checks](results/landsat_screening/scene_coverage.csv) · [Selected acquisitions and source bands](data/selected_landsat_scenes.csv)
3. The analysis aggregated the original Landsat pixel footprints to 500 m squares using land area and fit a model with date effects, tree and built fractions, their interaction, and the fraction of each square occupied by study land. Blocks of 4 km preserved nearby and repeated measurements when building intervals. [Analysis plan](docs/LANDSAT_ANALYSIS_PLAN.md)

[Read the full research report](docs/REPORT.md) · [Coverage and residual maps](results/landsat_analysis/figures/landsat_maps.png) · [Acquisition-level checks](results/landsat_analysis/figures/date_estimates.png) · [Data sources and licences](docs/DATA_SOURCES.md)

The [long-term research plan](docs/RESEARCH_ROADMAP.md) develops local comparisons, tests variation with built cover, and sets out how to validate both results on additional overpasses.

## Reproduce

The analysis ran on Linux x86-64 and Python 3.10.12. Run these commands from the repository root. Embedded snapshot tables provide a fixed model input.

```bash
python3 -m venv .venv
mkdir -p .scratch/tmp .scratch/cache/matplotlib
export TMPDIR="$PWD/.scratch/tmp"
export MPLCONFIGDIR="$PWD/.scratch/cache/matplotlib"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 GDAL_NUM_THREADS=1
.venv/bin/python -m pip install -r requirements-lock.txt
.venv/bin/python -m pip check
PYTHONPATH=scripts .venv/bin/python -B -m unittest discover -s tests -v
.venv/bin/python -B scripts/analyze_landsat.py \
  --panels data/snapshot/landsat_panels --output .scratch/reproduced-analysis
```

To regenerate the 500 m inputs from the public source rasters and draw the figures, run the following after setup. The manifest fixes all seven product IDs and band locations.

```bash
.venv/bin/python -B scripts/fetch_inputs.py spatial
.venv/bin/python -B scripts/make_landsat_panel.py \
  --scene-manifest data/selected_landsat_scenes.csv --output data/processed/landsat
.venv/bin/python -B scripts/analyze_landsat.py \
  --panels data/processed/landsat --output .scratch/reproduced-from-rasters
.venv/bin/python -B scripts/plot_landsat.py --panels data/processed/landsat \
  --results .scratch/reproduced-from-rasters --figures .scratch/reproduced-figures
```

The [data and outputs guide](docs/REPRODUCE.md) describes the archived files and the original MODIS and Landsat feasibility checks.

## Contributions

The applicant selected the project, provided the record of the earlier data work, chose the MODIS feasibility route and approved the selected-date Landsat analysis. Codex assisted with source verification, pipeline development, execution, figures, and writing.
