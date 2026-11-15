"""Exercise substitution support, fixed water geometry, and common-sample fits."""
import unittest

import numpy as np
import pandas as pd
from shapely.geometry import Polygon, box

from analyze_grass_tree_substitution import shared_fits
from grass_tree_water import distances, water_geometry
from prepare_grass_tree_substitution import COVER, RULES, match_group


class GrassTreeSubstitutionTests(unittest.TestCase):
    def test_pair_requires_rangeland_reduction_and_donor_share(self):
        cells = pd.DataFrame([
            ("z-good-donor", 0, 0, 0.10, 0.18, 0.60, 0.12, 1.00),
            ("a-target", 500, 0, 0.20, 0.08, 0.604, 0.116, 1.00),
            ("b-tree-gain-only", 0, 500, 0.20, 0.13, 0.578, 0.092, 1.00),
            ("x-no-donor", 0, 3000, 0.10, 0.09, 0.70, 0.11, 1.00),
            ("y-no-donor-target", 500, 3000, 0.20, 0.00, 0.69, 0.11, 1.00),
        ], columns=["cell_id", "easting_m", "northing_m", *COVER])
        records, degrees = match_group(cells, RULES["primary"])
        self.assertEqual(len(records), 1)
        self.assertEqual(records[0][:3], ("z-good-donor", "a-target", 500))
        np.testing.assert_allclose(records[0][3:5], [10.0, -10.0])
        self.assertEqual(degrees.loc[degrees.candidate_links > 0, "cell_id"].tolist(),
                         ["a-target", "z-good-donor"])

    def test_coast_excludes_holes_and_inland_excludes_harbour(self):
        hole = box(500, 500, 1300, 1300)
        land = Polygon(box(0, 0, 4000, 4000).exterior, [hole.exterior])
        hydro = [
            ("LARGE LAKE", hole), ("HARBOUR", box(4000, 500, 4800, 1300)),
            ("SMALL LAKE", box(2500, 2500, 2900, 2900)),
            ("COASTAL RESERVOIR", box(4000, 3000, 4800, 3800)),
        ]
        coast, components, inventory = water_geometry([land], hydro)
        self.assertEqual(set(inventory.name_largest_named_overlap),
                         {"LARGE LAKE", "SMALL LAKE", "COASTAL RESERVOIR"})
        self.assertEqual(inventory.primary.sum(), 2)
        grid = pd.DataFrame({"easting_m": [900, 2700], "northing_m": [900, 2700]})
        got = distances(grid, coast, components, inventory)
        np.testing.assert_allclose(got.coast_km, [0.9, 1.3])
        self.assertEqual(got.inland_km.iloc[0], 0)
        self.assertGreater(got.inland_km.iloc[1], 0)
        self.assertEqual(got.inland_010_km.iloc[1], 0)

    def test_full_and_base_use_identical_weighting_and_resamples(self):
        rng = np.random.default_rng(12)
        spatial_id = np.repeat(np.arange(12), 24)
        days = np.tile(np.repeat(np.arange(3), [4, 8, 12]), 12)
        x = rng.normal(size=(len(days), 7))
        x[:, 0] = 1 + 0.05 * x[:, 0]
        x[:, 5] += x[:, 0] + days
        true = np.array([-0.60, 0.20, -0.01, 0.11, 0.32, 0.08, -0.15])
        y = x @ true
        clusters = np.column_stack([spatial_id, np.zeros_like(spatial_id)])
        fits, descriptive, n = shared_fits(x, y, y / x[:, 0], clusters, days, draws=25)
        self.assertEqual(n, 12)
        self.assertEqual(len(descriptive), 26)
        root_w = 1 / np.sqrt(np.bincount(days)[days])
        base, _, _, _ = np.linalg.lstsq(x[:, :5] * root_w[:, None], y * root_w, rcond=None)
        np.testing.assert_allclose(fits["base"][0], base, atol=1e-12, rtol=0)
        np.testing.assert_allclose(fits["both"], np.tile(true, (26, 1)), atol=1e-12, rtol=0)
        np.testing.assert_allclose(fits["both"][:, 0] - fits["base"][:, 0], true[0] - fits["base"][:, 0],
                                   atol=1e-12, rtol=0)


if __name__ == "__main__":
    unittest.main()
