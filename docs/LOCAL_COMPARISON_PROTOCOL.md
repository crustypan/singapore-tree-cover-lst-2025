# Local comparison protocol for tree cover and daytime surface temperature

Stage 1 uses the seven 2025 Landsat overpasses and 500 m cells from the baseline study. It estimates differences in daytime surface temperature associated with 10 percentage points more tree cover among nearby grid cells with comparable built and land shares. The baseline outcomes were published before this follow-up. The rules below were fixed before calculating the local temperature estimates.

## Study population and pairing

Keep the baseline primary QA and spatial rules, including the WorldCover non-water mask. A cell has at least 50% mapped land and at least 80% valid Esri classification on that land. A cell-date has at least 60% of its land passing the original Landsat temperature and QA checks. Use Esri tree and built area as shares of valid classified land. Land share refers to the proportion of the entire grid square in the study land mask. Check the seven dates and eligible cell IDs against the committed baseline snapshot.

Use 4 km UTM squares with boundaries at multiples of 4,000 m in EPSG:32648. A candidate pair meets these limits.

| Quantity | Primary rule |
| --- | --- |
| Overpass and spatial block | Same acquisition date and same 4 km square |
| Distance between cell centers | At most 2,000 m |
| Absolute difference in built cover | At most 5 percentage points |
| Absolute difference in land share | At most 10 percentage points |
| Difference in tree cover | From 8 to 12 percentage points inclusive |

Form disjoint pairs within each date and block. The algorithm first finds the largest possible number of pairs. Among matchings of that size it minimizes the sum of squared differences in built share, land share, and distance, scaled by their respective limits, plus the squared deviation of the tree-cover gap from 10 points divided by 2 points squared. Round pair costs to six decimal places and break exact matching ties using sorted cell IDs. Orient each pair from lower to higher tree cover. A cell can appear again on a later date. Prepare the eligibility table and maps using coordinates, cover, and QA availability before joining Landsat temperatures. [Matching algorithm](https://networkx.org/documentation/networkx-3.4/reference/algorithms/generated/networkx.algorithms.matching.max_weight_matching.html) · [Outcome-free design principle](https://doi.org/10.1214/08-AOAS187)

Repeat matching with a 2 km shift of the 4 km block origin and, separately, with tighter center-distance and built-cover limits of 1,500 m and 3 points. Compare matched and full samples by geography, date availability, class composition, and the number of contributing blocks. The relative shares of rangeland, bare ground, flooded vegetation, and other classes help explain the land-cover tradeoff in the selected pairs.

## Temperature estimates

Fit an ordinary least squares model with a separate intercept for each date and 4 km block, linear tree cover, linear and squared built cover, and linear land share. Each date receives the same total weight among records used. Use groups with at least two eligible cells. Report the change in predicted cell-average LST in °C for 10 percentage points of tree cover at fixed built and land shares. The model uses the variation between cells on each date and within each block.

For matched pairs, compute the higher-tree cell's LST minus the lower-tree cell's LST. Report the actual difference and the result normalized to a ten-percentage-point tree difference. Give each pair equal weight within its date and average the seven date estimates equally. Refit the regression on matched endpoints to assess the role of the sample. The model and pairs describe local observed associations. A larger tree share changes the composition of remaining land cover, and local geography or land use may still contribute to the measured differences. [Spatial scale and confounding](https://doi.org/10.1214/10-STS326)

## Uncertainty and planned checks

Use 1,999 spatial block bootstrap draws for the primary regression and primary matched estimate. Keep the full date history of each sampled block together. Recompute date weights within each draw. Report percentile intervals, the number of contributing blocks, and the number of valid draws. Treat the seven dates as fixed and report date-specific estimates and leave-one-date-out results. Use 999 draws for the additional estimates. Block intervals approximate residual dependence between nearby locations. [Spatial resampling](https://doi.org/10.1214/009053606000000551)

- Check regression comparisons within 2 km, shifted 4 km, and 6 km blocks. Also resample nested 8 km units for the primary 4 km model.
- Compare equal date and equal cell-date weights, cells with valid observations on all seven dates, the Esri water mask, and the established 4 K and 2 K upper limits on reported uncertainty with at least 0.5 km cloud distance. Show any empty sample.
- Check a quadratic tree term centered at 25% tree cover and show the slopes at 10%, 25%, and 50% with nearby sample and block counts.
- Compare estimates from the two alternative pairing designs and inspect residual maps and queen-neighbor Moran's I by date.

Report the sample support, main estimates, and uncertainty for the seven observed overpasses. The interval quantifies variation over spatial blocks conditional on those dates. The geography of the sample and the remaining spatial structure matter for interpretation. The test of variation with built cover belongs to Stage 2.
