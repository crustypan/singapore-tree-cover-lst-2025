# Stage 9 influence of overpasses and blocks on the 4 K result

Use the same nine paired observations remaining under both Stage 8 matching rules after simultaneous omission of the two focal blocks. The study land at each endpoint has at least 60% coverage by temperature retrievals with reported uncertainty at most 4 K and distance from mapped cloud of at least 0.5 km. Preserve their ordered cell IDs, dates, 4 km block IDs, temperatures, and original cover gaps. [Stage 8 protocol](QUALITY_DENSE_SWAP_PROTOCOL.md) · [Stage 8 results](QUALITY_DENSE_SWAP_RESULTS.md)

## Pair temperature and weight ledger

Subtract the lower-tree endpoint temperature from the higher-tree endpoint temperature. Multiply that gap by five and divide by the observed tree-cover gap in percentage points to obtain the same five-point contrast as Stage 8. All six overpasses receive equal total weight. The weight of one observation is the reciprocal of six times the number of retained pairs on its date. Report the raw and scaled differences and each pair's contribution to the full estimate, calculated as its weight times its scaled difference. Verify that the contributions sum to the archived Stage 8 estimate.

## One deletion at a time

Start each run from all nine observations. Omit one overpass date at a time, then one geographic block at a time. A block omission removes all records in that block across all dates. On each resulting set of pairs, recalculate equal total weights by remaining date, distributing that date's weight over its remaining pairs. Refit the matched estimate and use the same 1,999 whole-block resamples and 95% percentile intervals as the Stage 8 matched result. The resampling keeps all retained date records in each resampled block together and keeps each observation's refitted weight fixed within that run.

Review the largest weights, signed contributions, changes in the point estimate, and which intervals include zero after a deletion. Report how many pairs, overpasses and blocks support each estimate. Apply the original minimum of two dates and three blocks to each remaining sample. These are sensitivity checks on the same selected observations. The spatial intervals apply to the dates and paired geography represented in each run. [Stage 6 matched estimation](DENSE_TREE_SWAP_PROTOCOL.md#five-point-temperature-contrast) · [MacKinnon, Nielsen, and Webb on cluster influence](https://doi.org/10.1016/j.jeconom.2022.04.001)
