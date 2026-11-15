# Stage 5 protocol for tree cover replacing grass and shrub cover

## Question and comparison

How much cooler or warmer are nearby cells with ten percentage points more tree cover and ten points less mapped grass and shrub cover, at the same built-area share? How far does the estimate move after including distances to the coast and large inland water on those same pairs?

Use the Stage 3 registry of 15 daytime Landsat overpasses from 2023–25. Each cell uses the Esri / Impact Observatory map of its acquisition year and the fixed WorldCover and URA study domain. Compute tree (class 2), built (class 7), and rangeland (class 11) fractions from 10 m pixel counts among classified non-water study-land centers. Group the remaining recognized classes as other. Tree, built, rangeland, and other sum to one. Rangeland is the product's classification of open herbaceous and some shrub or scrub cover. The contrast compares different cells observed on the same date. [Esri service legend](https://ic.imagery1.arcgis.com/arcgis/rest/services/Sentinel2_10m_LandCover/ImageServer/rasterFunctionInfos)

## Support and sample lock

Use the Stage 1 500 m UTM grid and unshifted 4 km blocks. Orient each candidate from lower to higher tree share within one date and one 4 km block. Apply these rules to unrounded fractions.

| Requirement | Primary design | Tighter design |
| --- | ---: | ---: |
| Distance between centers | ≤2,000 m | ≤1,500 m |
| Gain in tree cover | 8–12 percentage points | 8–12 percentage points |
| Loss in rangeland cover | 8–12 percentage points | 8–12 percentage points |
| Rangeland in the lower-tree cell | ≥10% | ≥10% |
| Absolute gap in tree plus rangeland | ≤2 percentage points | ≤1.5 percentage points |
| Absolute gaps in built and in other cover | ≤5 percentage points each | ≤3 percentage points each |
| Absolute gap in land fraction of the square | ≤10 percentage points | ≤10 percentage points |

An eligible cell-date has support if it has at least one qualifying link. Report cells and blocks with a link, the observed cover range, and the number of dates on which each location is supported. Select disjoint pairs, allowing each cell to enter once per date. Maximize pair count in each date and block, then minimize the sum of squared normalized deviations in tree gain, rangeland loss, built and other cover gaps, land share, and distance. Center tree gain and rangeland loss at the ten-point target. Divide deviations by the corresponding tolerances or distance limit. Use the existing deterministic tie break. Freeze the pair IDs and distance values before opening the temperature columns. Describe the geography and finite cover ranges represented in the resulting sample.

## Fixed water distances

Use the archived Singapore Land Authority National Map Polygon GeoJSON. Project both layers from WGS 84 to UTM Zone 48N. Repair invalid geometries and use polygon parts. For the coast, dissolve Singapore's Coastal Outlines layer and take exterior rings. Exclude features labeled Malaysia. For inland water, keep Hydrographic polygons inside the land outlines (including interior holes) along with named reservoir polygons regardless of how they overlap the outlines. This retains mapped reservoirs such as Serangoon and Marina that meet the generalized coast and removes mapped harbour basins outside the land outlines. Dissolve overlaps and keep polygon components of at least 0.25 km². For a sensitivity run, lower the cutoff to 0.10 km². Measure the shortest straight-line distance from each grid center to the coast line and to a qualifying inland water polygon. A center in a water polygon has distance zero. The map edition supplies fixed geographic covariates for all 15 overpasses. [SLA National Map Polygon](https://data.gov.sg/datasets/d_29f066d67df3eae91df8a42f443863c8/view) · [PUB catchments and reservoirs](https://www.pub.gov.sg/Public/WaterLoop/OurWaterStory/Local-Catchment-Water)

## Estimation and sensitivity

Take higher-tree minus lower-tree differences of Landsat LST and all regressors. Fit a weighted least squares model through the origin, which corresponds to a model with a pair intercept before differencing. Regress on the differences in tree, built, squared centered built, other, and study-land shares. Multiply the shares by ten before taking differences, and center built share at 50% before scaling and squaring. Each pair receives the same weight within its date, with weights summing to one per date.

Rangeland is the omitted component. Holding the built and other fractions and land share fixed, a ten-point increase in tree cover is necessarily a ten-point decrease in rangeland. The tree coefficient is the fitted surface temperature difference in °C for this substitution. As a descriptive check, also report the same-date paired temperature difference standardized to the pair's actual tree gap.

Add the higher-minus-lower differences in `log(1 + distance in km)` for the coast and inland water to the base model on the identical pairs and weights. Report the estimate with each distance in turn and both together. Measure the main change as the jointly adjusted tree coefficient minus its base-model counterpart. Compute its interval directly from the paired coefficient differences. Reuse the same 4 km block bootstrap counts across models and keep every date for each sampled block. Run 1,999 draws and give percentile 95% intervals conditional on the overpasses and frozen matching. A draw omitting a date or lacking full rank in any of the compared models is excluded from all intervals in that run. Report retained draws.

Check tighter pairs, 2025 and the additional 2023–24 dates separately, untransformed distance, and the 0.10 km² water definition. Compare both models on a shared sample in each check. Record the two point estimates when omitting each date in turn. Interpret the estimates as differences in daytime satellite surface temperature in the supported sample. Unmeasured differences in buildings, shade, species, terrain, irrigation, and wind can remain within a pair. [Lahiri and Zhu on spatial resampling](https://doi.org/10.1214/009053606000000551) · [Paciorek on spatial confounding](https://doi.org/10.1214/10-STS326)
