import unittest

import numpy as np
import pandas as pd

from analyze_local_comparison import draw_totals, estimate_pairs, within_model
from prepare_local_comparison import MatchRule, build_matches, match_group, minimum_cost_pairs


class MatchDesign(unittest.TestCase):
    def test_largest_matching_then_lowest_cost(self):
        # The lowest individual cost leaves two cells unmatched.
        edges = [(1, 2, 0), (0, 1, 2), (2, 3, 2), (0, 2, 5), (1, 3, 5)]
        self.assertEqual(minimum_cost_pairs(edges), [(0, 1), (2, 3)])
        self.assertEqual(minimum_cost_pairs(edges[::-1]), [(0, 1), (2, 3)])

    def test_exact_calipers_and_lower_to_higher_orientation(self):
        group = pd.DataFrame({
            "cell_id": ["a", "b", "c"], "tree_fraction": [0.18, 0.10, 0.181],
            "built_fraction": [0.47, 0.52, 0.571], "land_share": [0.70, 0.80, 0.801],
            "easting_m": [250, 2250, 2750], "northing_m": [250, 250, 250],
        })
        pairs, candidates, count = match_group(group, MatchRule())
        self.assertEqual(count, 1)
        self.assertEqual(candidates, {"a", "b"})
        self.assertEqual(pairs[0][:2], ("b", "a"))
        self.assertAlmostEqual(pairs[0][3], 8)

    def test_matching_never_crosses_a_block_or_a_date(self):
        obs = pd.DataFrame({
            "date": ["d1", "d1", "d2"], "cell_id": ["a", "b", "b"], "tree_fraction": [0.20, 0.10, 0.10],
            "built_fraction": [0.60, 0.60, 0.60], "land_share": [1.0, 1.0, 1.0],
            "easting_m": [3750, 4250, 4250], "northing_m": [250, 250, 250],
        })
        pairs, _, _, _ = build_matches(obs, MatchRule())
        self.assertTrue(pairs.empty)
        pairs, _, _, _ = build_matches(obs, MatchRule(block_shift=2000))
        self.assertEqual(list(pairs.date), ["d1"])


class LocalEstimation(unittest.TestCase):
    def test_within_model_handles_varying_date_by_area_baselines(self):
        rng = np.random.default_rng(802)
        panels = []
        for block in range(5):
            for date in ["2025-01-03", "2025-02-04"]:
                n = 25 if date == "2025-01-03" else 14
                t = rng.uniform(0.05, 0.30, n)
                b = rng.uniform(0.35, 0.55, n)
                land = rng.uniform(0.5, 1, n)
                bb = (b - 0.5) * 10
                y = rng.normal() * 3 - 0.65 * (t - 0.25) * 10 + bb + bb ** 2 + 3 * land
                panels.append(pd.DataFrame({
                    "date": date, "cell_id": [f"{block}_{i}" for i in range(n)],
                    "easting_m": block * 4000 + 250, "northing_m": 250,
                    "tree_fraction": t, "built_fraction": b, "land_share": land,
                    "response_c": y,
                }))
        fit, fitted, boot = within_model(pd.concat(panels), draws=19)
        self.assertAlmostEqual(fit["estimate_c_per_10pp"], -0.65)
        np.testing.assert_allclose(boot[:, 0], -0.65, atol=1e-9)
        np.testing.assert_allclose(fitted.residual_c, 0, atol=1e-9)

    def test_pairs_weight_dates_equally(self):
        pairs = pd.DataFrame({
            "date": ["d1", "d1", "d2"], "block_e": [1, 2, 2], "block_n": [0, 0, 0],
            "lower_tree_cell": ["a", "b", "b"], "higher_tree_cell": ["c", "d", "d"],
            "lst_per_10pp_c": [-2., -4., -1.], "lst_gap_c": [-2., -4., -1.],
        })
        fit = estimate_pairs(pairs, draws=0)
        self.assertAlmostEqual(fit["estimate_c_per_10pp"], -2.0)
        self.assertIsNone(fit["lower_c"])

    def test_bootstrap_preserves_all_dates_of_a_block(self):
        date = ["d1", "d2", "d1", "d2"]
        cluster = np.array([[0, 0], [0, 0], [0, 1], [0, 1]])
        value = np.array([[1.], [2.], [10.], [20.]])
        draws, n = draw_totals(value, cluster, date, 31)
        self.assertEqual(n, 2)
        np.testing.assert_array_equal(draws[:, 0, 0], draws[:, 1, 0])
        np.testing.assert_array_equal(2 * draws[:, 0, 1], draws[:, 1, 1])
