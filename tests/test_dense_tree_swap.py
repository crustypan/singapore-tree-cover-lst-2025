"""Check how the dense-area design counts and estimates substitution support."""
import unittest

import numpy as np
import pandas as pd

from analyze_dense_tree_swap import fitted_contrasts
from dense_swap_model import COEFFICIENTS, assess
from prepare_dense_tree_swap import RULES, donor_funnel
from prepare_grass_tree_substitution import COVER, RULES as STAGE5_RULES, build_matches
from prepare_local_comparison import assign_blocks


class DenseTreeSwapTests(unittest.TestCase):
    def test_distinguishes_donor_scarcity_calipers_and_block_boundary(self):
        day = pd.DataFrame([
            ("a", 250, 250, .05, .08, .87, .00, 1.0),
            ("b", 750, 250, .10, .03, .87, .00, 1.0),
            ("c", 3750, 3750, .12, .08, .80, .00, 1.0),
            ("d", 4250, 3750, .17, .03, .80, .00, 1.0),
            ("e", 250, 8250, .03, .08, .89, .00, 1.0),
            ("f", 750, 8250, .08, .03, .86, .03, 1.0),
            ("g", 250, 12250, .06, .04, .90, .00, 1.0),
        ], columns=["cell_id", "easting_m", "northing_m", *COVER])
        day = assign_blocks(day)
        primary, pairs, n_links = donor_funnel(day, RULES["5pp"])
        self.assertEqual(n_links, 1)
        self.assertEqual(primary.set_index("cell_id").step.to_dict(),
                         {"a": 6, "b": 0, "c": 4, "d": 0, "e": 2, "f": 0, "g": 0})
        self.assertEqual(pairs[0][2:4], ("a", "b"))
        wider, pairs, n_links = donor_funnel(day, RULES["5pp_stage5_balance"])
        self.assertEqual(n_links, 2)
        self.assertEqual(wider.loc[wider.cell_id.eq("e"), "step"].item(), 6)

        empty, _, daily = build_matches(day.assign(date="2025-10-02"), STAGE5_RULES["primary"])
        self.assertTrue(empty.empty)
        self.assertEqual(daily.matched_blocks.tolist(), [0])

    def test_supported_fit_returns_a_five_point_temperature_difference(self):
        rng = np.random.default_rng(20)
        vals = rng.uniform(-.4, .4, size=(60, 4))
        vals[:, 0] += 1
        pairs = pd.DataFrame(vals * 5, columns=COEFFICIENTS).assign(
            date=np.tile(["2025-01-01", "2025-04-01", "2025-07-01"], 20),
            block_e=np.repeat(np.arange(10), 6), block_n=0,
            lower_tree_cell=[f"l{i}" for i in range(len(vals))],
            higher_tree_cell=[f"h{i}" for i in range(len(vals))],
            lst_gap_c=vals @ np.array([-.42, .19, -.08, .035]),
        )
        support = assess(pairs)
        self.assertTrue(support["adjusted_ready"])
        fits, omission, _ = fitted_contrasts(pairs, support)
        np.testing.assert_allclose(fits["adjusted"], -.42, rtol=0, atol=1e-10)
        np.testing.assert_allclose(omission["adjusted"], -.42, rtol=0, atol=1e-10)
        self.assertEqual(len(fits["matched"]), len(fits["adjusted"]))

        pairs.loc[pairs.block_e.ne(0), "other_gap_pp"] = 0
        fragile = assess(pairs)
        self.assertTrue(fragile["matched_ready"])
        self.assertFalse(fragile["adjusted_ready"])
        self.assertFalse(fragile["rank_after_each_block_omission"])
        self.assertFalse(assess(pairs[pairs.block_e.lt(2)])["matched_ready"])


if __name__ == "__main__":
    unittest.main()
