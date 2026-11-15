import unittest

import numpy as np
import pandas as pd

from analyze_built_cover import Contrast, contrast_matrix, fit_contrasts
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


class BuiltCoverContrast(unittest.TestCase):
    def test_additive_model_has_same_finite_difference_at_both_references(self):
        contrast = Contrast(2.5, 65, 85)
        c = contrast_matrix(contrast, "additive")
        np.testing.assert_allclose(c[0], c[1], atol=1e-12)
        np.testing.assert_allclose(c[2], 0, atol=1e-12)

    def test_local_model_recovers_known_joint_contrast(self):
        rng = np.random.default_rng(6731)
        panels = []
        contrast = Contrast(2.5, 65, 85)
        for block in range(4):
            for date in ["2025-01-03", "2025-05-03"]:
                for band, center in [("lower", 0.65), ("higher", 0.85)]:
                    n = 15 if band == "lower" else 28
                    tree = rng.uniform(0.01, 0.13, n)
                    built = rng.uniform(center - 0.015, center + 0.015, n)
                    land = rng.uniform(0.65, 1, n)
                    t = 10 * (tree - 0.075)
                    b = 10 * (built - 0.75)
                    y = 2 * block + (5 if date == "2025-05-03" else 0)
                    y = y - 0.50 * t + 0.25 * t ** 2 + b ** 2 - 0.12 * t * b + land
                    panels.append(pd.DataFrame({
                        "date": date, "band": band, "cell_id": [f"{block}_{band}_{i}" for i in range(n)],
                        "easting_m": 4000 * block + 250, "northing_m": 250,
                        "tree_fraction": tree, "built_fraction": built, "land_share": land, "response_c": y,
                    }))
        r = fit_contrasts(pd.concat(panels), contrast, draws=41)
        self.assertEqual(r["draws_retained"], 41)
        self.assertAlmostEqual(r["contrasts"]["lower_built"]["estimate_c"], -0.38)
        self.assertAlmostEqual(r["contrasts"]["higher_built"]["estimate_c"], -0.62)
        self.assertAlmostEqual(r["contrasts"]["higher_minus_lower"]["estimate_c"], -0.24)
        self.assertAlmostEqual(r["contrasts"]["higher_minus_lower"]["lower_c"], -0.24)
