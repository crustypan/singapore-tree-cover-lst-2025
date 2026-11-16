"""Check the quality cutoffs and the paired two-block omission contrast."""
import unittest

import numpy as np
import pandas as pd

from analyze_dense_tree_swap import fitted_contrasts
from analyze_influential_pairs import masks, paired_changes
from dense_swap_model import COEFFICIENTS, assess


class InfluentialPairsTests(unittest.TestCase):
    def test_original_and_stricter_temperature_masks(self):
        dn = np.full(6, 48000, dtype=np.uint16)
        pixel = np.array([0, 0, 0, 0, 8, 0], dtype=np.uint16)
        distance = np.array([50, 49, 60, 100, 100, 100])
        uncertainty = np.array([200, 200, 201, 401, 200, -9999])
        (original, qa4, qa2), _ = masks(dn, pixel, distance, uncertainty)
        np.testing.assert_array_equal(original, [1, 1, 1, 1, 0, 1])
        np.testing.assert_array_equal(qa4, [1, 0, 1, 0, 0, 0])
        np.testing.assert_array_equal(qa2, [1, 0, 0, 0, 0, 0])

    def test_two_blocks_omit_only_their_rows_and_share_resampling(self):
        dates = ["2025-01-03", "2025-02-04"]
        rows = [(date, b, 0, f"lo{b}", f"hi{b}", 5, 0, 0, 0, float(b))
                for b in range(6) for date in dates]
        full = pd.DataFrame(rows, columns=["date", "block_e", "block_n", "lower_tree_cell",
                                          "higher_tree_cell", *COEFFICIENTS, "lst_gap_c"])
        retained = ~full.block_e.isin([4, 5]).to_numpy()
        reduced = full.loc[retained]
        before, _, _ = fitted_contrasts(full, assess(full))
        after, _, _ = fitted_contrasts(reduced, assess(reduced))
        shifts = paired_changes(full, retained, before, after, ["matched", "adjusted"])
        for row in shifts:
            self.assertAlmostEqual(row["change_c_per_5pp"], -1)
            self.assertLess(row["lower_c"], -1)
            self.assertGreater(row["upper_c"], -1)
            self.assertGreaterEqual(row["paired_draws_retained"], 1990)


if __name__ == "__main__":
    unittest.main()
