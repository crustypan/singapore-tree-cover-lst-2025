import unittest

import pandas as pd

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
