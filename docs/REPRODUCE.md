# Data and outputs

Use the environment and commands in the [README](../README.md). All code runs from the repository root. The pipeline writes new download caches and intermediate products inside the project. The result paths passed to each script can point to a new directory for comparison against the checked-in outputs.

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
