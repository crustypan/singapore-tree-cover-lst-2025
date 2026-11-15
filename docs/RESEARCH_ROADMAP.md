# Research plan for tree cover and daytime surface temperature in Singapore

The study proceeds through two questions. First estimate how daytime land surface temperature varies with tree cover among nearby, comparable grid cells. Then test whether this association differs with built-area cover. An extension to additional overpasses tests how well both findings travel across dates and weather conditions once the spatial comparisons are understood. The work is planned over roughly eight to twelve weeks, with each stage completed and documented before the next analysis begins.

## Starting evidence

The [2025 baseline](REPORT.md) contains seven quality-screened Landsat 8/9 overpasses and 14,678 observations from 2,491 cells on the 500 m grid. At 50% built cover the baseline linear model associates 10 percentage points more tree cover with 0.73 °C lower surface temperature. The corresponding estimate with separate date and 4 km spatial-block intercepts is 0.64 °C lower. Residual Moran's I ranges from 0.56 to 0.75 across dates. The interval for the tree-by-built interaction includes zero once the model allows separate curvature in tree and built cover. These findings motivate two distinct analyses. [Baseline models and sensitivity estimates](../results/landsat_analysis/robustness.csv) · [Residual diagnostics](../results/landsat_analysis/date_diagnostics.csv)

Two geographic checks guide the design. Tree and built fractions have a Pearson correlation of −0.85 in the sampled cells. Within 10 percentage points of the 25%, 50%, and 75% built-cover references, respectively 4, 19, and 69 distinct cells have tree cover between 10% and 20%. Comparing all three reference slopes therefore requires an explicit check of which land-cover combinations are observed.

Within-date comparisons appear feasible at a broader set of built shares. A pilot count used the baseline fitted cell IDs and dates, cell coordinates, tree and built fractions, and land share. It counted candidate pairs within the same 4 km UTM square and date, within 2 km of each other, separated by 8–12 points of tree cover, and within 5 points of built cover and 10 points of land share. Across dates, 213 distinct cells appeared in at least one qualifying pair. Per-date counts were 112–199 cells and 22–33 occupied blocks containing at least one candidate pair. Candidate links can share cells and are a measure of possible comparisons. The cutoffs are pilot values for checking the design. [Cell table](../results/landsat_analysis/cell_availability.csv) · [Cell-date IDs](../results/landsat_analysis/fitted_cells.csv)

| Overpass | Distinct cells with a candidate match | 4 km blocks with at least one candidate pair |
| --- | ---: | ---: |
| 3 January | 189 | 30 |
| 4 February | 151 | 24 |
| 3 May | 112 | 22 |
| 6 July | 171 | 30 |
| 2 October | 199 | 33 |
| 18 October | 150 | 24 |
| 5 December | 184 | 31 |

## Stage 1 Local association on the seven Landsat dates

**Question.** How much does daytime surface temperature differ between cells in the same area and overpass when tree cover is higher and built cover is comparable?

### Define the comparisons

Use the seven acquisitions, the 500 m UTM grid, and the data-quality rules in the [baseline Landsat plan](LANDSAT_ANALYSIS_PLAN.md). Start from the versioned input panels and scene manifest. Inspect the archived land and usable-temperature fractions, the distribution of tree cover within each block and date, and the cells that can supply nearby comparisons. Map their locations and compare their land cover and date availability to the full baseline sample. Record sources and hashes, the analysis population, model settings, and the reasons for any exclusions before refitting temperatures. The seven dates have already been analyzed, so this stage is a follow-up analysis.

The comparison changes tree cover at a given built-area share. A higher tree share replaces part of the other classified land, and the mix of replaced land can differ by location. Check where that mix is dominated by grass and scrub, bare ground, or other mapped classes. These are observational contrasts in satellite surface temperature at the times of the seven overpasses.

### Estimate and check the association

Fit a regression with an intercept for each combination of date and 4 km block. Each local area gets its own baseline temperature for each overpass. Adjust for built fraction with a low-flexibility curve and include the proportion of the square covered by study land. Start with a single tree-cover slope in °C per 10 percentage points. Give the seven dates equal total weight and report date-specific estimates. Reproduce the baseline weighting by cell-date as a sensitivity analysis. Joint geography and date intercepts address temperature differences across local areas that vary between overpasses. Roads, elevation, coastline proximity, and building form can still affect comparisons within the same block. [Paciorek on spatial scale and confounding](https://doi.org/10.1214/10-STS326)

Cross-check the regression with nearby pairs. Use the pilot rule to assess overlap and freeze the final rule from covariates and availability before inspecting pairwise temperature differences. Use disjoint pairs within each date, with deterministic matching that prioritizes similar built share, similar land share, and proximity. Report the actual tree-cover difference and compare LST on the same date. If the paired and regression estimates differ, inspect the locations, composition, and coverage of their samples separately.

Resample spatial blocks with all their dates kept together. Use a fixed seed, report the number of contributing blocks and usable resamples, and check sensitivity to block size and grid origin. Compare estimates after leaving out one overpass, using cells observed on all seven dates, changing the water mask, and applying the existing stricter QA checks. Report the change in coverage and which land-cover mixes remain. Inspect residual maps and residual spatial autocorrelation by date. A low-flexibility nonlinear tree curve can show whether a single slope hides changes over the observed tree-cover range. [Lahiri and Zhu on spatial resampling](https://doi.org/10.1214/009053606000000551)

### Stage 1 deliverables

Allow about two to three weeks. Deliver a short method note, a map of covered and comparable locations, estimates for the pooled and individual dates, a pairing check, and tables showing sensitivity to coverage and model choices. Record how much of the full sample supports each analysis. Reproduce results from the committed input snapshot in the research environment and record the sample and software versions. The milestone is a reproducible estimate with stated geographic support and uncertainty, including an inconclusive result if warranted.

## Stage 2 Variation with built-area cover

**Question.** For the same increase in tree cover from a comparable starting level, is the LST contrast different at different built-area shares?

### Establish which contrasts can be estimated

Map the joint distribution of tree and built cover using distinct cells. Report the numbers of cells, overpasses, and spatial blocks in candidate ranges, including the blocks that can support both ends of a 10-percentage-point tree contrast within an overpass. Start with the 25%, 50%, and 75% references from the baseline. Choose the actual contrast levels and tree-cover interval using overlap and geographic balance before fitting an interaction. Keep the same tree-cover interval at every built-cover level being compared. If a built-cover range has sparse data, report its coverage and focus inference on the ranges with comparable cells. [Hainmueller, Mummolo, and Xu on interaction diagnostics](https://doi.org/10.1017/pan.2018.46)

### Compare the contrasts

Start with smooth additive relationships for tree and built cover and the local date-by-block controls from Stage 1. Fit a second specification that allows the tree response to vary with built cover, using limited flexibility justified by the observed support. Compare the predicted LST difference for an identical ten-point tree-cover contrast at the chosen built-cover levels. Where data permit, compare the same geographic blocks and dates or use common weights so a change in geographic mix is visible. Estimate a spatial-bootstrap interval for the **difference between the contrasts**, refitting both comparisons in each draw. Check sensitivity to plausible curvature and to omitting individual spatial blocks or dates that contribute most of the support.

### Stage 2 deliverables

Allow about three to four weeks after Stage 1. Deliver the joint coverage plot, a table of supported contrasts with uncertainty, the sensitivity to model shape, and a statement of the built-cover ranges for which the result is informative. A result that varies with the model or has wide intervals is a completed research finding when the design support and limits are documented.

## Stage 3 Validation on additional overpasses

**Question.** How consistent are the two associations across additional clear-sky dates and Singapore's monsoon and inter-monsoon periods?

Start by auditing the 30 Tier 1 scenes outside the seven prior selections for 2025. Evaluate usable area, geographic overlap with the comparison sample, time of overpass, land-cover support, and the number of dates in each weather period. Report the coverage of both a common set of locations and the full usable area. If more years or a different temperature source are needed, show the measured coverage and proposed scope before changing the year, spatial scale, or quality rules. For additional years, check that land-cover data match the measurement period and that apparent changes near coasts and development are plausible.

Once the extension has usable support and an agreed scope, lock its date selection and models, then evaluate the local tree contrast and any supported built-cover contrast on the new dates. Report date-to-date variation separately from the spatial interval conditional on each set of observed dates. Balance contributions by weather period where coverage permits and retain the clear-sky overpass definition. The expansion can improve temporal coverage and recover some previously unobserved cells. Recheck joint land-cover support on the enlarged observation set. [Meteorological Service Singapore climate periods](https://www.weather.gov.sg/climate-climate-of-singapore/) · [Ermida et al. on clear-sky sampling](https://doi.org/10.1029/2018JD029354)

Allow about three to five weeks for the coverage audit and a validation run if the data qualify. Deliver the candidate-date inventory, period and area coverage, and a comparison with the earlier results. If usable dates remain concentrated in a few periods, report that pattern and state which periods the estimates represent.

## Stage 4 Audit of the land-cover comparison

Use the Stage 3 loss of the high-tree endpoint at 65% built cover to define the review cells. Isolate the tree, built and water stability rules and inspect those cells in dated Sentinel-2 images. Reconstruct a second classification with the released Dynamic World weights. Recount the same local comparisons on a common pixel footprint and report the number of supported blocks in each band. The deliverables are the selected image inventory, case evidence, coverage diagnostics, support counts and code. [Stage 4 protocol](COVER_SUPPORT_AUDIT_PROTOCOL.md)

## Research record and handoff

Run each stage in the independent laboratory-server project. Commit the analysis decisions before fitting, then version the code, input and output manifests, figures, sensitivity results, and English findings in the dedicated GitHub repository. Keep the dated 2025 baseline available for comparison. Record actual applicant contributions when methodological choices, results review, or interpretation take place, along with AI-assisted work. A stage is ready to share once its sample, estimates, and figures can be reproduced and its conclusion matches its sampled dates and geographic support.
