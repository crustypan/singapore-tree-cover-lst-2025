# 2025 temperature and land-cover feasibility

## MODIS annual route

The annual coverage rule failed. The NASA CMR calendar contains 365 distinct MOD11A1.061 granules for tile h28v08. The server read 209 matching raw daytime temperature and QC pairs from the Microsoft mirror, with 156 dates lacking a pair. All of July through November was absent from the mirror snapshot. There were 24 dates in the retrieved set with at least one acceptable observation intersecting WorldCover land. The best-observed land cell had 11 such dates and the highest month count was 7.

| WorldCover land metric | Coverage |
| --- | ---: |
| At least one observed QA-passing date | 35.83% |
| All prespecified rules on the retrieved dates | 0.00% |
| Upper bound if every missing day qualifies at every pixel | 7.70% |
| Prespecified annual feasibility requirement | 80.00% |

The upper bound has a simple check. For each cell and month, add the count of missing source dates to the measured QA count. This simulates the most favorable input to the days, months and quarterly requirements, all of which are monotone in those counts. In April–June the mirror had all 91 pairs. Only 8.88% of WorldCover land met the quarterly criterion even before adding any restrictions from other quarters. Applying all three annual criteria at once yields the 7.70% upper bound. This assumes the retrieved pairs faithfully represent the corresponding NASA granules; source IDs and raw windows are archived so that this can be checked.

Among the retrieved observations passing the protocol and intersecting some WorldCover land, the 1st, 50th and 99th percentiles of LST were 25.72, 28.85 and 34.13 °C. The widest observed daily coverage was 23.70% on 23 December. The accepted observations are concentrated spatially and temporally. Coverage uncertainty here comes from the missing dates and the finite spatial grid. The upper bound is a deterministic coverage bound. There is no statistical estimate of the tree–LST association from this route.

Area calculations were repeated with fine-grid spacings of 18.53, 9.27 and 4.63 m inside the same native MODIS footprint. At all settings the upper bound rounded to 7.70%, and the area meeting the observed annual criteria was 0.00%. The QA-passing temperature observations remain at the native 926.63 m resolution. [Machine-readable audit](../results/modis_public_mirror/coverage.json) · [Monthly coverage](../results/modis_public_mirror/monthly_coverage.csv) · [Daily inventory](../results/modis_public_mirror/source_inventory.csv)

## Esri 2025 land cover

Inside the 55 URA planning areas, the Esri 10 m GeoTIFF contains 7,966,101 cell centers, of which 963,727 are water. There are 1,644,652 tree cells and 4,640,826 built-area cells. The 535,005 cells recorded as clouds in the earlier notes are class 11, rangeland. Class 10 has 1,478 cloud cells. The corrected native-grid classified share among non-water cells is 99.979%. At the tested MODIS scale, all Esri candidate land falls in cells with at least 80% classified candidate land. Classification accuracy is a separate source of uncertainty. There are 387 cells labeled snow and ice within URA, an implausible class for this location. Counting them as unclassified would still leave 99.973% of the native non-water area classified. [Class counts and palette captured in the audit](../results/modis_public_mirror/coverage.json) · [Esri publisher legend](https://ic.imagery1.arcgis.com/arcgis/rest/services/Sentinel2_10m_LandCover/ImageServer)

## Landsat candidate dates

The server screened 44 Landsat 8/9 Collection 2 L2SP scenes, including 37 Tier 1 scenes. It decoded QA_PIXEL bits 0–5 and excluded zero ST_B10. It then applied a study-specific plausibility screen of 15–70 °C, as motivated by the prior anomalous scene, and calculated area coverage from WorldCover non-water and Esri alternative masks. The screening uses nearest neighbor to align the 30 m delivery data to a grid derived from Esri’s 10 m UTM grid. Its weights count the nine 10 m land centers associated with each 30 m cell and allocate each center to its planning area. The inherited scene screen requires at least 60% of WorldCover land and at least 44 of the 55 planning areas each with at least 30% coverage.

| Tier 1 acquisition | WorldCover land with passing QA and LST range | Planning areas with at least 30% valid land |
| --- | ---: | ---: |
| 3 Jan 2025 | 84.28% | 55 |
| 4 Feb 2025 | 66.22% | 44 |
| 3 May 2025 | 62.96% | 46 |
| 6 Jul 2025 | 79.05% | 54 |
| 2 Oct 2025 | 86.00% | 55 |
| 18 Oct 2025 | 71.73% | 48 |
| 5 Dec 2025 | 84.90% | 54 |

The 27 January Tier 2 case reproduced the earlier problem. Standard QA alone leaves 77.38% of land and 52 planning areas with at least 30% coverage. Its median is −19.03 °C and 97.44% of its QA-passing land samples fall below 15 °C. Coverage falls to 1.99% after the temperature screen.

The seven Tier 1 scenes have a 90.121% spatial union on the 30 m screening grid and 90.258% under the Esri water mask. Directly assigning the same 10 m land centers to the native Landsat footprints moves the WorldCover union to 90.083%. The historical WorldCover result was 89.985% under the original rules. The old grid and code were unavailable for a like-for-like comparison. The two reprocessed estimates are on one side of the 90% threshold and the historical estimate is on the other. The annual coverage decision remains contingent on reproducing the original specification.

| Minimum number of passing dates on native Landsat footprints | WorldCover land |
| ---: | ---: |
| 1 | 90.08% |
| 3 | 89.74% |
| 5 | 79.24% |
| 7 | 34.35% |

The repeat coverage argues for retaining partially observed cells and comparing a complete-date subset as a separate sensitivity analysis if a selected-date study is authorized. Other date-level checks could include cloud distance, product-reported temperature uncertainty, and changing the land-cover mask. [Per-scene data and source IDs](../results/landsat_screening/scene_coverage.csv) · [Direct footprint results](../results/landsat_screening/native_footprint_scenes.csv) · [Repeat coverage](../results/landsat_screening/repeat_coverage.csv)

## Scope decision

The applicant selected the seven Landsat acquisitions and an approximately 500 m grid. See the [selected-date research report](REPORT.md) for the completed analysis. The annual feasibility evidence is retained as the record of the scope decision.
