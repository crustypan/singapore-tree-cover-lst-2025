import unittest

import pandas as pd

from prepare_built_support import band_support, selected_samples


class BuiltCoverSupport(unittest.TestCase):
    def test_support_requires_both_tree_endpoints_in_same_date_block(self):
        data = pd.DataFrame({
            "date": ["a", "a", "b", "b", "c", "a"],
            "block_e": [1, 1, 1, 2, 1, 1], "block_n": [0] * 6,
            "cell_id": ["a", "b", "a", "c", "b", "d"],
            "easting_m": [4250, 4750, 4250, 8250, 4750, 5250], "northing_m": [250] * 6,
            "tree_fraction": [0.025, 0.125, 0.025, 0.125, 0.125, 0.125],
            "built_fraction": [0.65, 0.65, 0.65, 0.65, 0.65, 0.75],
            "land_share": [1.0] * 6,
        })
        stat, support, strict, sample, _ = band_support(data, 2.5, 65)
        self.assertEqual(support, {("a", 1, 0)})
        self.assertEqual(strict, support)
        self.assertEqual(set(sample.cell_id), {"a", "b"})
        self.assertEqual(stat["high_endpoint_cells"], 2)

    def test_two_bands_have_disjoint_cells_with_a_common_block_date(self):
        frame = pd.DataFrame({
            "date": ["a"] * 4, "cell_id": ["a", "b", "c", "d"],
            "block_e": [1] * 4, "block_n": [0] * 4,
            "easting_m": [4250, 4750, 4250, 4750], "northing_m": [250, 250, 750, 750],
            "tree_fraction": [0.025, 0.125, 0.025, 0.125],
            "built_fraction": [0.65, 0.65, 0.85, 0.85], "land_share": [1.0] * 4,
        })
        rows, _, common, strict = selected_samples(frame, 2.5, 65, 85)
        self.assertEqual(common, {("a", 1, 0)})
        self.assertEqual(common, strict)
        self.assertTrue(rows.shared.all())
        self.assertEqual(rows.groupby("band").size().to_dict(), {"higher": 2, "lower": 2})
