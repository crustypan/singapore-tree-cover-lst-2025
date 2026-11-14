# Selected-date Landsat analysis plan

## Scope and observations

The applicant chose seven 2025 Landsat 8/9 Tier 1 daytime acquisitions after reviewing the feasibility results. Use the immutable list of scene IDs with passing Tier 1 screening in the archived [coverage summary](../results/landsat_screening/summary.json). The calendar dates are 3 January, 4 February, 3 May, 6 July, 2 October, 18 October, and 5 December 2025. The analysis targets within-date spatial differences among observable parts of Singapore on those overpasses.

The main spatial unit is a 500 × 500 m square on EPSG:32648 with both edges at integer multiples of 500 m. Use the 10 m Esri grid as an area-integration grid. A center belongs to its enclosing 500 m square, the URA boundary polygon, and the original 30 m Landsat delivery pixel. Sample the native Landsat footprint without interpolation across neighboring values. The original thermal measurement has an effective resolution coarser than the 30 m delivery grid. Keep one LST mean per 500 m square and overpass.

Use the WorldCover 2021 non-water mask for the primary area. A cell is eligible if at least 50% of the full square is mapped as URA land and at least 80% of that land has a valid non-water Esri 2025 class. The tree and built fractions use just that valid classified area as their common denominator. Esri classes 1 (water), 10 (cloud), and 0 (unclassified) do not enter a classified fraction. Land-cover values are kept fixed across dates.

Convert ST_B10 to degrees Celsius as `DN × 0.00341802 + 149.0 − 273.15`. Require nonzero ST_B10, QA_PIXEL bits 0–5 clear, and LST between 15 and 70 °C inclusive. Average the passing temperature over WorldCover land centers for a cell-date and require at least 60% of its candidate land area to pass. Use the same aggregation for Esri non-water land in the alternative-mask check. The area percentages are integration weights. They do not add temperature samples.

## Model and interpretation

Fit ordinary least squares to eligible cell-date means. Give each cell-date the same weight. Use categorical date intercepts, tree fraction, built fraction centered at 50%, their interaction, and the fraction of the whole square included as study land. The land fraction accounts for an aspect of coastal mixing and boundary location. Report the interaction as the change in the tree slope for a 10-percentage-point increase in built fraction. Report tree slopes in °C per 10 percentage points of tree fraction at 25%, 50%, and 75% built coverage. Tabulate the number and range of cells around each reference value within a ±10 percentage-point window. At fixed built fraction, more trees displace a mixture of other mapped land categories.

Fit accompanying date-adjusted specifications with tree fraction alone and with the additive tree and built terms. Use the interaction model as the main specification. Check date-specific slopes and leave-one-date-out estimates, compare included and excluded land-cover distributions, and measure within-date residual autocorrelation with queen-neighbor Moran's I. Record the number of cells and repeated observations from the same location.

Resample nonoverlapping 4 km square blocks with replacement. Keep all observations from each selected block together, including repeated dates from each cell. Refit date intercepts and coefficients for 1,999 draws and use percentile intervals. Compare results with 2 km and 6 km blocks and a 2 km shift of the 4 km grid. Intervals quantify uncertainty conditional on these seven dates and the chosen spatial domain. They are approximate because neighboring blocks can remain correlated.

## Robustness checks

- Change cell-date QA area from 60% to 40% and 80%.
- Require explicit cloud distance of at least 0.5 km and ST_QA at most 2 K, and compare filtering on a common sample.
- Recompute with the Esri water mask, a 90% minimum cell land share, a 250 m shift of the 500 m grid, and an aligned 1 km grid.
- Compare a balanced seven-date sample, equal total weight per spatial cell, and a fit with 4 km spatial-block intercepts.
- Compare fixed land-cover fractions to fractions measured on the part of each cell observed on each date. Use these comparisons to diagnose potential within-cell selection.

The observations measure radiating surface temperature under usable clear-sky conditions. The analysis does not identify a causal cooling effect or quantify air temperature, nighttime conditions, or human heat exposure. Interpret regional differences with possible confounding from geography, urban morphology, shade, water, and other influences. The dates were selected for usable satellite coverage and cannot characterize the full year.
