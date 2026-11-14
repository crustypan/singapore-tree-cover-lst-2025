# Singapore tree cover and daytime land surface temperature in 2025

## Research question

Among Singapore land areas with dependable satellite observations in 2025, how does the daytime land surface temperature (LST) vary with the fraction of land classified as trees? Does the association vary with the fraction classified as built area?

The geographic boundary is the URA Master Plan 2025 Planning Area Boundary (No Sea). The study uses NASA Terra MOD11A1 Collection 6.1 for temperature feasibility and the Esri and Impact Observatory 2025 Sentinel-2 classifications for land cover. One MODIS native footprint is one spatial unit. The nominal 1 km sinusoidal grid has a spacing of about 927 m. The 10 m land-cover grid resolves the composition within these units.

## Rules carried forward from the handoff

- Convert the raw daytime LST to Celsius as `DN × 0.02 − 273.15`. Accept only nonzero DN and Celsius values from −10 to 70 inclusive.
- For `QC_Day`, require bits 0–1 equal to 0, bits 2–3 equal to 0, bits 4–5 at most 1, and bits 6–7 at most 1. Flag value zero is a valid best-quality code.
- A cell must have at least 24 valid dates, measurements in at least 9 distinct calendar months, and at least 2 valid dates in each calendar quarter.
- At least 80% of the URA land area, defined by valid ESA WorldCover 2021 non-water classes, must fall within cells meeting these time rules. Also compute coverage under the Esri 2025 water mask as a sensitivity check. Esri cloud and unclassified pixels remain identifiable in the denominator.
- A land-cover cell is usable if at least 80% of its Esri non-water candidate area has a valid class. Fractions of trees and built area use the classified non-water area of the cell as denominator. Record water, cloud and missing area explicitly.
- Leave the research question, year, scale and thresholds in place unless the applicant approves a change after seeing the measurements.

The daily calendar uses the date of the source MOD11A1 granule. Separate gaps in dated file availability from observed days on which no cell passes QA. Use the source CRS and affine transform for all grid alignment. Count and add up the monthly and quarterly observations from the accepted daily values.

## Land-area calculation

The pipeline samples categorical inputs by nearest neighbor on a finer grid nested within each MODIS cell. It rasterizes the URA boundary by fine-pixel center and sums equal-area fine-pixel counts per MODIS cell. It uses WorldCover class 80 as the primary water mask and Esri class 1 as the alternative. The land-cover audit also counts the classes directly on Esri's native 10 m UTM grid. A native scale would introduce multiple records of the same daily temperature if its temperature were copied onto smaller cells.

## Analysis after the quality decision

Summarize spatial and temporal coverage before fitting any association. Compare the included and excluded land areas by land-cover composition and geography. Use date effects to account for citywide differences between observed overpasses and a tree-by-built interaction to estimate whether the tree association changes across the observed built-area range. Express slopes in degrees Celsius per 10 percentage points of tree fraction. Report the observed joint support of the fractions, especially where a high tree fraction and a high built fraction cannot coexist.

Account for repeated observations at each cell and spatial dependence in uncertainty estimates, with spatial blocks and sensitivity to block size and placement. Check residual spatial autocorrelation, the stability of estimates across calendar periods, and the influence of coastal or mixed cells. Aggregate neighboring native cells only as a coarser-scale sensitivity check. Changes to valid-day cutoffs belong in sensitivity checks after the primary rule is evaluated.

Daytime satellite LST measures the radiating surface at Terra overpass times. Cloud screening makes observations selective. The estimates describe clear-sky overpasses and an observational association. They cannot establish causal cooling or quantify air temperature, personal heat exposure, heat-health risk or nighttime conditions.

## Starting evidence

The applicant's handoff records 89.985% annual Landsat coverage with the primary land mask after excluding the anomalous Tier 2 scene `LC09_L2SP_125059_20250127_02_T2`. The Landsat threshold was 90%. The handoff records roughly 22 valid days at the best-observed pixel in the original MODIS export and a corrected diagnosis of zero pixels meeting the temporal thresholds. The original rasters, scripts and full protocols were unavailable in the local locations searched. These figures describe earlier work and are distinct from measurements produced by this repository.

The applicant selected the subject and the subsequent MODIS feasibility check, accepted the nominal MODIS spatial scale, and supplied the history of the original Earth Engine export. Codex assisted with source checking, implementation, execution and writing. Record any later applicant work as it is performed.
