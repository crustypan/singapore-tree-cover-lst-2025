# Stage 8 protocol for the 4 K dense-area comparison

## Frozen pairs and temperature coverage

Start with the two Stage 6 pooled dense-area pair lists audited in Stage 7. The primary list contains 36 paired observations and the wider list contains 44. Both endpoints have at least 80% mapped built cover. Dates, sites, ordered endpoints and cover gaps come from the frozen lists. Every record is considered for the quality restriction. Use the two blocks identified in the Stage 7 design for the combined omission. [Stage 6 protocol](DENSE_TREE_SWAP_PROTOCOL.md) · [Stage 7 audit](INFLUENTIAL_PAIR_AUDIT_RESULTS.md)

Use the earlier Landsat 4 K screen. A native Landsat footprint must pass the original fill, cloud, shadow, snow, cirrus and plausible temperature checks. It must also have `ST_QA` at most 4 K and `ST_CDIST` at least 0.50 km. Count the 10 m study-land centers falling inside qualifying temperature footprints. Require at least 60% coverage of the study land in each of the two 500 m cells. Evaluate the exact fraction using integer counts. Use the mean temperature over the same qualifying centers at each endpoint. [USGS Level 2 product guide](https://www.usgs.gov/media/files/landsat-8-9-collection-2-level-2-science-product-guide) · [Stage 7 temperature screen](INFLUENTIAL_PAIR_AUDIT_PROTOCOL.md#image-and-temperature-review)

Reconcile the coverage counts in the Stage 3 panels with the original archived Landsat windows and study-land geometry. Record the passing pair-date IDs, distinct spatial pairs, overpasses, and 4 km blocks before fitting temperatures. From each passing list, drop every record in either focal block to form the omitted sample. Preserve the existing pair IDs and cover gaps throughout.

## Estimation and comparison

Apply the Stage 6 support criteria separately to each quality-filtered sample before and after the combined omission. A matched contrast requires at least two overpasses and three geographic blocks. For the residual cover-adjusted estimate, use the original model through the origin and require full rank with all pairs and after omitting each individual block. At least 95% of the planned spatial resamples must also retain rank.

Normalize each higher-minus-lower temperature difference to a five-point increase using the observed tree gap. Give each overpass equal total weight on each sample. Use the original block bootstrap with 1,999 draws and 95% percentile intervals. For each estimator supported both before and after dropping the two blocks, report the omitted-sample estimate minus the full-sample estimate. Use the same original block multiplicities in both estimates for the interval on the change, as in Stage 7. [Stage 6 temperature fit](DENSE_TREE_SWAP_PROTOCOL.md#five-point-temperature-contrast) · [Stage 7 omission](INFLUENTIAL_PAIR_AUDIT_PROTOCOL.md#omit-both-blocks)

The contrasts describe Landsat daytime surface temperature among the matched high-built cells and dates that retain adequate coverage under the specified screen. The attrition tables show which geography and overpasses that estimate represents.
