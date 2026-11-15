import unittest

import pandas as pd

from analyze_overpass_validation import date_and_period_tables
from audit_heldout_overpasses import SEASONS, period
from prepare_overpass_validation import compare_annual_cover


class OverpassValidation(unittest.TestCase):
    def test_common_stable_requires_observation_and_similar_cover_in_all_years(self):
        grid = pd.DataFrame({
            "cell_id": ["unchanged", "missing", "water_shift", "built_shift"],
            "easting_m": [250, 750, 1250, 1750], "northing_m": [250] * 4,
            "land_centers": [2500] * 4, "land_share": [1.0] * 4,
            "tree_fraction": [0.10] * 4, "built_fraction": [0.60] * 4, "esri_water_share": [0.01] * 4,
        })
        grids = {
            2023: grid.assign(built_fraction=[0.55, 0.60, 0.60, 0.56]),
            2024: grid.assign(esri_water_share=[0.01, 0.01, 0.07, 0.01],
                              built_fraction=[0.60, 0.60, 0.60, 0.64]),
            2025: grid,
        }
        observations = {
            year: grid[["cell_id"]].assign(date=f"{year}-01-01") for year in grids
        }
        observations[2024] = observations[2024].query("cell_id != 'missing'")
        result = compare_annual_cover(grids, observations)
        self.assertEqual(result.loc[result.common_stable, "cell_id"].tolist(), ["unchanged"])
        self.assertEqual(result.loc[result.common, "cell_id"].tolist(), ["unchanged", "water_shift", "built_shift"])

    def test_period_balancing_weights_periods_not_acquisitions(self):
        dates = ["2024-01-01", "2024-03-13", "2024-10-23", "2024-09-13"]
        table = pd.DataFrame({"cohort": ["validation"] * 4, "analysis": ["slope"] * 4,
                              "period": [period(d) for d in dates], "estimate_c": [-1., -3., -1., None]})
        periods, balanced = date_and_period_tables(table)
        self.assertEqual(periods.period.tolist(), list(SEASONS))
        self.assertEqual(periods.selected_dates.tolist(), [2, 0, 1, 1])
        self.assertEqual(periods.estimated_dates.tolist(), [2, 0, 0, 1])
        self.assertAlmostEqual(balanced.iloc[0]["equal_period_mean_c"], -1.5)
