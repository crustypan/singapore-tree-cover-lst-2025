import unittest

import numpy as np
import pandas as pd

from analyze_landsat import fit, moran
from make_landsat_panel import point_indices, points_within_native


class LandsatPipeline(unittest.TestCase):
    def test_grid_centers_and_native_pixel_alignment(self):
        x = np.array([340005, 340025, 340495, 340505])
        y = np.array([144995, 144985, 144995, 144975])
        idx, center, count = point_indices(x, y, 500, 0)
        np.testing.assert_array_equal(center, [[340250, 144750], [340750, 144750]])
        np.testing.assert_array_equal(idx, [0, 0, 0, 1])
        np.testing.assert_array_equal(count, [3, 1])
        rows, cols = points_within_native(x, y, (30, 0, 340000, 0, -30, 145000), (3, 18))
        np.testing.assert_array_equal(rows, [0, 0, 0, 0])
        np.testing.assert_array_equal(cols, [0, 0, 16, 16])

    def test_date_controls_and_tree_interaction_recover_known_slopes(self):
        rng = np.random.default_rng(123)
        n = 60
        tree = rng.uniform(0, 0.4, n)
        built = rng.uniform(0.3, 0.6, n)
        land = rng.uniform(0.6, 1, n)
        block_effect = np.repeat([0.5, 0, -1], 20)
        x = pd.DataFrame({
            "tree_fraction": tree, "built_fraction": built, "land_share": land,
            "cell_id": np.arange(n).astype(str), "easting_m": np.repeat([250, 4250, 8250], 20),
            "northing_m": np.repeat([250, 4250, 8250], 20),
        })
        t, b = 10 * tree, 10 * (built - 0.50)
        y = 31.0 - 0.7 * t + 0.2 * b - 0.04 * t * b - 1.8 * land
        data = pd.concat([x.assign(date="2025-01-03", response_c=y),
                          x.assign(date="2025-12-05", response_c=y + 6)])
        result, residual = fit(data)
        self.assertAlmostEqual(result["interaction"]["estimate_c"], -0.04)
        self.assertAlmostEqual(result["slopes"][1]["estimate_c"], -0.7)
        self.assertAlmostEqual(result["coefficients"]["date_2025-12-05"], 6)
        np.testing.assert_allclose(residual, 0, atol=1e-10)
        data["response_c"] += np.tile(block_effect, 2)
        local_result, _ = fit(data, spatial_intercepts=True)
        self.assertAlmostEqual(local_result["interaction"]["estimate_c"], -0.04)
        self.assertAlmostEqual(local_result["slopes"][1]["estimate_c"], -0.7)

    def test_moran_uses_queen_neighbors_and_omits_missing_cells(self):
        easting = np.array([250, 750, 250, 1250])
        northing = np.array([250, 250, 750, 250])
        value, edge_count = moran(np.array([1.0, 2.0, 2.0, 1.0]), easting, northing)
        self.assertEqual(edge_count, 4)
        self.assertAlmostEqual(value, -0.5)
