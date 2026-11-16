# Data and outputs

Use the environment and commands in the [README](../README.md). All code runs from the repository root. The pipeline writes new download caches and intermediate products inside the project. The result paths passed to each script can point to a new directory for comparison against the checked-in outputs.

## Stage 6 Tree-for-rangeland comparison in dense built areas

The design uses the committed Stage 5 annual class counts and the Stage 3 cell-date registry. Preparation verifies the input hashes, checks exact fractions reconstructed from integer pixel counts, and records the support and rank checks before the temperature step. The geographic figures also require the original URA boundary.

```bash
.venv/bin/python -B scripts/fetch_inputs.py spatial
.venv/bin/python -B scripts/prepare_dense_tree_swap.py --output .scratch/dense-repro/design
.venv/bin/python -B scripts/plot_dense_tree_swap.py --design .scratch/dense-repro/design \
  --output .scratch/dense-repro/figures
.venv/bin/python -B scripts/analyze_dense_tree_swap.py --design .scratch/dense-repro/design \
  --output .scratch/dense-repro/analysis
.venv/bin/python -B scripts/plot_dense_tree_swap.py --analysis .scratch/dense-repro/analysis \
  --output .scratch/dense-repro/figures
```

Use `--design results/dense_tree_swap/design` with the temperature command to start from the frozen pairs.

| File | Contents |
| --- | --- |
| [Rangeland distributions](../results/dense_tree_swap/design/composition_by_band.csv) | Cell-year denominator, percentiles and counts above the five- and ten-point thresholds by built band. |
| [Cumulative support by band](../results/dense_tree_swap/design/support_by_band.csv) and [support by date](../results/dense_tree_swap/design/support_by_date.csv) | Lower-tree cell counts at each step, including pair availability and the geography restrictions. |
| [Dense cell registry](../results/dense_tree_swap/design/dense_cells.csv.gz) and [pair registry](../results/dense_tree_swap/design/pairs.csv) | Deepest support step for each high-built cell-date, ordered endpoints for each selected pair and actual composition gaps. |
| [Frozen design](../results/dense_tree_swap/design/design_record.json) | Rules, match and model availability, block concentration and input and output checksums. |
| [Temperature contrasts](../results/dense_tree_swap/analysis/estimates.csv), [balance](../results/dense_tree_swap/analysis/pair_balance.csv) and [omit-one-block results](../results/dense_tree_swap/analysis/leave_one_block_out.csv) | Local contrasts and intervals, realized differences in cover and the estimates after omitting individual geographic blocks. |
| [Pair overlap](../results/dense_tree_swap/analysis/pair_overlap.csv) | The pair records and weighted contributions shared by or unique to the primary and wider matching designs. |
| [Paired LST values](../results/dense_tree_swap/analysis/paired_temperatures.csv) and [analysis record](../results/dense_tree_swap/analysis/analysis_record.json) | The two recorded temperatures in each fitted pair and the resampling seed and output checksums. |

[Protocol](DENSE_TREE_SWAP_PROTOCOL.md) · [Results](DENSE_TREE_SWAP_RESULTS.md)

## Stage 5 tree-for-rangeland substitution

The frozen design supplies the 15-date comparison registry, class fractions and water distances. Use the archived design to reproduce the estimates and figures. To rebuild the design from the original spatial inputs, run the preparation step as well. The preparation script verifies each annual Esri window, fetching missing windows via the Stage 3 loader. The archived SLA GeoJSON fixes the map edition for the water calculation.

```bash
.venv/bin/python -B scripts/fetch_inputs.py spatial
# Optional rebuild from the source rasters
.venv/bin/python -B scripts/prepare_grass_tree_substitution.py \
  --output .scratch/substitution-repro/design

design=results/grass_tree_substitution/design
.venv/bin/python -B scripts/analyze_grass_tree_substitution.py --design "$design" \
  --output .scratch/substitution-repro/analysis
.venv/bin/python -B scripts/plot_grass_tree_substitution.py --design "$design" \
  --output .scratch/substitution-repro/figures/comparable_grids.png
.venv/bin/python -B scripts/plot_grass_tree_substitution.py --analysis .scratch/substitution-repro/analysis \
  --output .scratch/substitution-repro/figures/temperature_contrasts.png
```

Set `design=.scratch/substitution-repro/design` to use the rebuilt design. The comparable-grids figure also reads the original URA boundary. The estimate script verifies the frozen inputs, uses all temperatures on the locked pair list, and reports the change between models from common spatial resamples.

| File | Contents |
| --- | --- |
| [Design record](../results/grass_tree_substitution/design/design_record.json) | Matched sample totals, cutoffs and SHA-256 checks. |
| [Cell-date support](../results/grass_tree_substitution/design/cell_date_support.csv) and [counts by date](../results/grass_tree_substitution/design/support_by_date.csv) | Eligible IDs, number of qualifying links, disjoint pairing flags and geographic support. |
| [Annual cell covariates](../results/grass_tree_substitution/design/cell_covariates.csv) and [water components](../results/grass_tree_substitution/design/water_components.csv) | Class counts, cover shares and distances, plus inland water component areas and source names. |
| [Primary pairs](../results/grass_tree_substitution/design/pairs_primary.csv) and [tighter pairs](../results/grass_tree_substitution/design/pairs_tighter.csv) | Date, ordered grid IDs, spatial block, center separation and cover gaps. |
| [Estimates](../results/grass_tree_substitution/analysis/estimates.csv) and [analysis record](../results/grass_tree_substitution/analysis/analysis_record.json) | Same-sample model estimates and changes, bootstrap counts, pair-list hashes and coefficients. |
| [Pair balance](../results/grass_tree_substitution/analysis/pair_balance_by_date.csv), [LST differences](../results/grass_tree_substitution/analysis/paired_temperatures.csv) and [date influence](../results/grass_tree_substitution/analysis/leave_one_overpass_out.csv) | Actual contrasts and balance per date, temperatures at the two endpoints, and the result after leaving out each overpass. |

[Protocol](GRASS_TREE_SUBSTITUTION_PROTOCOL.md) · [Results](GRASS_TREE_SUBSTITUTION_RESULTS.md)

## Stage 4 land-cover and comparison-support audit

The committed cover table and source-image chips can reproduce the support counts and figures without an Earth Engine account. The case review table records interpretation of the chips. The archived Sentinel-2 manifest fixes the optical dates and unsigned source URLs.

```bash
.venv/bin/python -B scripts/fetch_inputs.py spatial
.venv/bin/python -B scripts/prepare_cover_support_audit.py --output .scratch/cover-repro/design
.venv/bin/python -B scripts/analyze_cover_support_audit.py count --design .scratch/cover-repro/design \
  --output .scratch/cover-repro/analysis
.venv/bin/python -B scripts/plot_cover_support_audit.py --design .scratch/cover-repro/design \
  --analysis .scratch/cover-repro/analysis \
  --review results/cover_support_audit/analysis/case_review.csv --output .scratch/cover-repro/figures
```

To regenerate the optical snapshot, create a separate CPU environment and download the original forward model from the versioned repository. The inference script checks the three weight digests.

```bash
.venv/bin/python -m venv .scratch/cover-support/tf-venv
.scratch/cover-support/tf-venv/bin/python -m pip install -r requirements-dynamic-world-lock.txt
.scratch/cover-support/tf-venv/bin/python -m pip check
dw_forward=.scratch/cover-support/model/model/forward
dw_source=https://raw.githubusercontent.com/google/dynamicworld/0dd288fae5fb5291d988e5ea3c12eb91df481ec5/model/forward
mkdir -p "$dw_forward/variables"
for file in saved_model.pb variables/variables.index variables/variables.data-00000-of-00001; do
  curl --fail --silent --show-error --location "$dw_source/$file" --output "$dw_forward/$file"
done
export TF_NUM_INTEROP_THREADS=1 TF_NUM_INTRAOP_THREADS=2 TF_CPP_MIN_LOG_LEVEL=2
for year in 2023 2024 2025; do
  .venv/bin/python -B scripts/fetch_cover_audit_sentinel.py fetch --year "$year" \
    --output results/cover_support_audit/imagery
  .scratch/cover-support/tf-venv/bin/python -B scripts/run_cover_audit_dynamic_world.py --year "$year" \
    --output .scratch/cover-support/inference
done
.venv/bin/python -B scripts/analyze_cover_support_audit.py extract --snapshot .scratch/cover-repro/snapshot
```

The fetch caches clipped data from the frozen products. Run it with `select` and a scratch output directory to repeat the STAC search and local cloud ranking. The seasonal selection can change if the provider reprocesses or adds an item. [Stage 4 protocol](COVER_SUPPORT_AUDIT_PROTOCOL.md) · [Stage 4 findings](COVER_SUPPORT_AUDIT_RESULTS.md) · [Snapshot provenance](../data/snapshot/cover_support_audit/provenance.json)

## Stage 3 validation on additional overpasses

The candidate manifests pin the Tier 1 product IDs and unsigned asset URLs. The audit recomputes the original scene rule and 500 m coverage from cropped Landsat windows. It uses each year's Esri classification, with the same WorldCover domain and URA geometry. The design step freezes the usable IDs, shared geography and comparisons before the temperature fit.

```bash
.venv/bin/python -B scripts/fetch_inputs.py spatial
for year in 2023 2024 2025; do
  .venv/bin/python -B scripts/audit_heldout_overpasses.py --year "$year" \
    --candidates "results/overpass_validation/audit/$year/candidate_scenes.csv" \
    --output ".scratch/validation/audit/$year"
done
.venv/bin/python -B scripts/prepare_overpass_validation.py --audit .scratch/validation/audit \
  --output .scratch/validation/design
.venv/bin/python -B scripts/make_validation_panel.py --audit .scratch/validation/audit \
  --design .scratch/validation/design --output .scratch/validation/snapshots
.venv/bin/python -B scripts/analyze_overpass_validation.py --design .scratch/validation/design \
  --snapshots .scratch/validation/snapshots --output .scratch/validation/analysis
.venv/bin/python -B scripts/plot_overpass_validation.py --audit .scratch/validation/audit \
  --design .scratch/validation/design --results .scratch/validation/analysis --output .scratch/validation/figures
```

To start with the committed snapshots and frozen design, run `.venv/bin/python -B scripts/analyze_overpass_validation.py --output .scratch/validation/analysis`. The estimate script checks design and panel hashes and verifies the eligibility registry. The snapshot coverage fields are area weights on the same grid as the original Landsat panel.

| File | Contents |
| --- | --- |
| `results/overpass_validation/audit/<year>/scene_coverage.csv` | Scene-level screen, overpass time, calendar period, overlap with earlier samples and support for comparisons. |
| `results/overpass_validation/audit/<year>/usable_cell_dates.csv.gz` | Date, cell ID and passing area for usable cell-dates. |
| [Validation design](../results/overpass_validation/design/design_record.json) | Frozen acquisitions, covariate stability rule, source and design hashes. |
| [Annual land-cover comparison](../results/overpass_validation/design/annual_cover_comparison.csv) | Annual class shares, cover ranges and observation counts for each cell. |
| [Local estimates](../results/overpass_validation/analysis/local_estimates.csv) and [built-cover estimates](../results/overpass_validation/analysis/built_cover_estimates.csv) | Discovery and validation estimates on the full, common and common stable geographies. |
| [Date estimates](../results/overpass_validation/analysis/date_estimates.csv) and [period estimates](../results/overpass_validation/analysis/period_estimates.csv) | Date-specific spatial intervals and descriptive temporal variation. |

[Validation protocol](OVERPASS_VALIDATION_PROTOCOL.md) · [Stage 3 results](OVERPASS_VALIDATION_RESULTS.md)

## Stage 2 variation with built cover

The design uses the archived primary 500 m input panel and its accepted cell-date IDs. It reads only design columns and hashes the complete snapshot. Fit the temperature model after generating and inspecting the support maps. The spatial download supplies the boundary for the geographic map.

```bash
.venv/bin/python -B scripts/fetch_inputs.py spatial
.venv/bin/python -B scripts/prepare_built_support.py --output .scratch/built-reproduction/design
.venv/bin/python -B scripts/plot_built_support.py --design .scratch/built-reproduction/design \
  --output .scratch/built-reproduction/design
.venv/bin/python -B scripts/analyze_built_cover.py --design .scratch/built-reproduction/design \
  --output .scratch/built-reproduction/analysis
.venv/bin/python -B scripts/plot_built_results.py --results .scratch/built-reproduction/analysis \
  --output .scratch/built-reproduction/figures
```

Use `--design results/built_cover_moderation/design` in the analysis command to start from the committed design. The analysis checks the recorded hashes before reading temperatures.

| File | Contents |
| --- | --- |
| [Joint density](../results/built_cover_moderation/design/joint_density.csv) | Distinct usable cells by five-point tree and built bins. |
| [Band support](../results/built_cover_moderation/design/band_support.csv) and [shared support](../results/built_cover_moderation/design/shared_support.csv) | Counts near each endpoint and numbers of eligible spatial blocks and block-dates, with and without the nearby-pair rules. |
| [Frozen design](../results/built_cover_moderation/design/design_record.json) | Selected reference contrasts, comparable sample counts and SHA-256 checks for inputs and design tables. |
| [Fitting registry](../results/built_cover_moderation/design/fit_sample.csv) and [strict pairs](../results/built_cover_moderation/design/strict_pairs.csv) | Cell-date IDs and within-band pair IDs determined before the temperature fit. |
| [Estimates](../results/built_cover_moderation/analysis/estimates.json) and [robustness table](../results/built_cover_moderation/analysis/robustness.csv) | Finite ten-point contrasts, cross-band difference, conditional intervals, fit metrics and draw counts. |

[Protocol](BUILT_COVER_PROTOCOL.md) · [Stage 2 results](BUILT_COVER_RESULTS.md)

## Stage 1 Local comparisons

The design uses the archived 500 m input panel and the baseline usable cell-date IDs. The source rasters supply the map boundary and Esri class counts. Generate the design and its map before fitting the models.

```bash
.venv/bin/python -B scripts/fetch_inputs.py spatial
.venv/bin/python -B scripts/prepare_local_comparison.py --output .scratch/local-reproduction/design
.venv/bin/python -B scripts/plot_local_design.py --design .scratch/local-reproduction/design \
  --output .scratch/local-reproduction/design/comparable_grid_cells.png
.venv/bin/python -B scripts/analyze_local_comparison.py --design .scratch/local-reproduction/design \
  --output .scratch/local-reproduction/analysis
.venv/bin/python -B scripts/plot_local_results.py --results .scratch/local-reproduction/analysis \
  --output .scratch/local-reproduction/figures
```

The committed design contains the three lists of disjoint pairs, a table of matched dates by cell, per-date counts, counts of Esri classes per cell, and input and design SHA-256 hashes. Start at the analysis step and use `--design results/local_comparison/design` to fit from the frozen design and temperature panel. Results include the oriented temperature contrasts, sample and balance diagnostics, regression residuals, and sensitivity tables. [Protocol](LOCAL_COMPARISON_PROTOCOL.md) · [Stage 1 results](LOCAL_COMPARISON_RESULTS.md) · [Design record](../results/local_comparison/design/design_record.json)

## Landsat association study

| File | Contents |
| --- | --- |
| [Selected scene manifest](../data/selected_landsat_scenes.csv) | Seven date-time stamps and USGS product IDs with source URLs for ST_B10, QA_PIXEL, ST_CDIST and ST_QA. |
| [Input provenance](../data/snapshot/landsat_panels/provenance.json) | SHA-256 of spatial inputs, scene manifest, cropped source windows and decompressed panel CSV files. |
| `data/snapshot/landsat_panels/*.csv.gz` | Input snapshots for the primary, alternative water mask, shifted, and 1 km grids. Each has all seven dates per cell, including failed cell-dates. |
| [Model estimates](../results/landsat_analysis/estimates.json) | Sample totals, coefficients, tree slopes, interactions, fit metrics and bootstrap settings. |
| [Availability by cell](../results/landsat_analysis/cell_availability.csv) | UTM centers, spatial eligibility, land shares and accepted date counts. |
| [Fitted cells](../results/landsat_analysis/fitted_cells.csv) | The accepted cell-dates with observed temperatures and residuals from the pooled fit. |
| [Date diagnostics](../results/landsat_analysis/date_diagnostics.csv) | Usable area and sample size, mean and percentile LST, class mix and residual Moran's I. |
| [Date estimates](../results/landsat_analysis/date_estimates.csv) and [leave-one-date-out](../results/landsat_analysis/leave_one_date_out.csv) | Date-specific estimates with conditional intervals, and pooled point estimates dropping one date at a time. |
| [Robustness checks](../results/landsat_analysis/robustness.csv) | Alternative spatial and QA rules, weights, block intercepts, and exploratory curvature and trimming. |

In the panel CSVs, coordinates are UTM meters and an ID is the easting and northing of a cell center. Each candidate land center represents 100 m² before gridding. Counts are area weights, while `*_share` and `*_fraction` range from 0 to 1. Missing class fractions or temperatures have blank fields. `strict_*` fields refer to the cloud distance ≥0.5 km and ST_QA ≤2 K screen. `qa4_*` fields use the same distance with an exploratory 4 K cutoff. The original provider QA and temperature bands are part of the large cache reconstructed from the manifest.

## Initial feasibility checks

The public 2025 archives can be revisited with the source-collection commands below. To reproduce the MODIS result for the original snapshot, copy the checked-in `npz` and inventory into the raw data directory after downloading the spatial rasters.

```bash
.venv/bin/python -B scripts/fetch_inputs.py spatial
cp data/snapshot/modis_public_daily.npz data/raw/
cp results/modis_public_mirror/source_inventory.csv data/raw/modis_public_daily_inventory.csv
.venv/bin/python -B scripts/audit_coverage.py --output .scratch/reproduced-modis
.venv/bin/python -B scripts/screen_landsat.py --output .scratch/reproduced-landsat
```

Use `.venv/bin/python -B scripts/fetch_inputs.py modis` to refresh the public MODIS calendar and raster subsets. This command can retrieve different mirror coverage if the provider updates its archive. [MODIS audit](../results/modis_public_mirror/coverage.json) · [Historical methodology and thresholds](PROTOCOL.md) · [Optional Earth Engine export](MODIS_EXPORT.md)
