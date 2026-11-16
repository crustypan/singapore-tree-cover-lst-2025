# Stage 6 protocol for tree cover replacing rangeland in dense built areas

## Populations and land cover

Use the 500 m cells and 15 daytime Landsat overpasses in 2023–25 registered in Stage 3. Take tree, built and rangeland shares from the annual pixel counts frozen in Stage 5. Rangeland is Esri class 11, which includes grass and some low shrubs and scrub. All four shares use the same classified, non-water study-land denominator. The study uses mapped built-area share as its density measure. [Stage 5 class definitions](GRASS_TREE_SUBSTITUTION_PROTOCOL.md) · [Esri product legend](https://ic.imagery1.arcgis.com/arcgis/rest/services/Sentinel2_10m_LandCover/ImageServer/rasterFunctionInfos)

Use five built-share bands below 50%, 50–<70%, 70–<80%, 80–<90%, and 90–100%. The dense population pools cells with at least 80% built cover, which includes the prior 85% built-cover reference. Both endpoints of any pair must belong to the analysis band in the year of their overpass. The pooled dense analysis can compare cells on opposite sides of the 90% boundary. Count rangeland availability per cell-year, keeping annual classifications separate. Count actual opportunities to compare per cell-date.

## Which requirement removes support

Repeat the Stage 5 ten-point rule within each band and the pooled dense population. Then repeat with the five-point target. For each date and population, count the lower-tree cell once at its deepest supported step, even if it has several candidate partners. The summary reports cumulative counts at each step.

1. Start with every eligible cell-date.
2. Retain possible lower-tree donors with at least the target amount of rangeland.
3. Require a same-date partner anywhere in the same built-share band with an appropriate tree gain and rangeland loss.
4. Require balance in the other components and in study-land share.
5. Require centers within 2 km.
6. Require the same unshifted 4 km UTM block.
7. Select a maximum number of disjoint pairs within each date and block and then minimize total matching cost. Each cell can appear once per date. Cost is the sum of squared deviations from the target tree gain and rangeland loss and squared differences in built, other and land shares and distance, all scaled by their respective limits. Use the existing deterministic tie break.

The counts follow this order. A cell can fail more than one requirement, and a retained partner can change as constraints are added. Keep the donor, local candidate and pair counts separate to distinguish the available land mix from the joint cover and geography requirements. The Stage 5 benchmark is re-matched with the built-band restrictions applied before matching.

| Fraction or gap | Ten-point benchmark | Five-point comparison |
| --- | ---: | ---: |
| Rangeland in the lower-tree cell | ≥10% | ≥5% |
| Tree gain and rangeland loss | 8–12 points each | 4–6 points each |
| Absolute gap in combined tree and rangeland | ≤2 points | ≤1 point |
| Absolute built and other cover gaps | ≤5 points each | ≤2.5 points each |
| Absolute study-land share gap | ≤10 points | ≤10 points |

The five-point balance sensitivity retains the 4–6 point tree and rangeland requirements. It applies the Stage 5 limits to the combined, built, and other share gaps. Geographic limits and temperature eligibility remain the same.

## Five-point temperature contrast

The primary population requires at least 80% built cover at each endpoint. Also evaluate the two dense bands and the wider-balance sensitivity. Freeze pair IDs, counts and model support before opening the cell-temperature columns.

For a sample represented on at least two dates and in at least three 4 km blocks, compute an observed paired contrast by multiplying each higher-minus-lower daytime LST difference by five divided by the actual tree gap in percentage points. Use equal total weight per overpass. Three blocks let the leave-one-block-out check retain spatial replication. The matched comparison approximates the substitution using the observed composition gaps. Its balance checks report the actual tree, rangeland and built-area differences.

Also fit a pair-difference model through the origin with tree, built, other and study-land share differences measured in five-point units. This linear adjustment accounts for residual cover gaps in the closely matched sample. An identically zero difference needs no adjustment. The coefficient on tree gives the five-point exchange with rangeland at fixed built, other, and study-land shares. Estimate this coefficient when the design has full rank on the original sample and under every deletion of a spatial block, and when at least 95% of the planned block draws retain rank. Evaluate those conditions using only the frozen differences in cover. Report descriptive support for strata too sparse for a temperature estimate.

Report 95% percentile intervals from 1,999 draws of whole 4 km blocks with replacement. Keep all date records for each sampled block and retain the original observation weights in every draw. An overpass confined to an omitted block is absent from that draw. When both models are estimable, use the same full-rank resamples for their intervals. Compute point estimates after leaving out each block to diagnose geographic influence. Spatial intervals based on a small number of blocks are approximate. Report distinct spatial pairs, dates, blocks, and the distribution of block weights with the estimates. The scope is the locally supported sample and the observed daytime overpasses. [Lahiri and Zhu on spatial resampling](https://doi.org/10.1214/009053606000000551) · [MacKinnon, Nielsen, and Webb on cluster inference](https://doi.org/10.1016/j.jeconom.2022.04.001)
