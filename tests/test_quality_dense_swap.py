"""Verify the inclusive area threshold and the frozen pair inclusion flags."""
import unittest

import pandas as pd

from prepare_quality_dense_swap import coverage_ok, screen


class QualityDenseSwapTests(unittest.TestCase):
    def test_each_endpoint_needs_sixty_percent_and_omission_drops_whole_blocks(self):
        endpoints = pd.DataFrame([
            ("a", 2500, 1500), ("b", 2500, 1499),
            ("c", 1500, 900), ("d", 2500, 2000),
            ("e", 2500, 2000), ("f", 2500, 2300),
        ], columns=["cell_id", "land_centers", "qa4_lst_centers"]).assign(date="2025-02-04")
        endpoints["valid_lst_centers"] = endpoints.land_centers
        endpoints["qa4_lst_share"] = endpoints.qa4_lst_centers / endpoints.land_centers
        endpoints["passes_qa4"] = coverage_ok(endpoints.qa4_lst_centers, endpoints.land_centers)
        pairs = pd.DataFrame([
            ("a", "b", 89, 36), ("c", "d", 91, 35), ("e", "f", 90, 36),
        ], columns=["lower_tree_cell", "higher_tree_cell", "block_e", "block_n"]).assign(date="2025-02-04")
        omitted = pd.DataFrame([(89, 36), (91, 35)], columns=["block_e", "block_n"])
        result = screen(pairs, endpoints, omitted)
        self.assertEqual(result.passes_qa4.tolist(), [False, True, True])
        self.assertEqual(result.passes_qa4_and_outside.tolist(), [False, False, True])
        self.assertEqual(result.lower_tree_cell.tolist(), ["a", "c", "e"])


if __name__ == "__main__":
    unittest.main()
