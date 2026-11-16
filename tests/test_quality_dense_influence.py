"""Check date reweighting when a spatial omission removes part of another date."""
import unittest

import numpy as np
import pandas as pd

from analyze_quality_dense_influence import influence, matched_effect
from prepare_quality_dense_influence import date_weights


class QualityDenseInfluenceTests(unittest.TestCase):
    def test_block_deletion_reweights_pairs_and_drops_empty_dates(self):
        pairs = pd.DataFrame([
            ("2025-01-03", 1, 0, 10.0),
            ("2025-01-03", 2, 0, 1.0),
            ("2025-02-04", 1, 0, 8.0),
            ("2025-07-06", 3, 0, -3.0),
            ("2025-10-02", 4, 0, -4.0),
        ], columns=["date", "block_e", "block_n", "lst_gap_c"])
        pairs = pairs.assign(
            pair_id=np.arange(1, 6), lower_tree_cell=[f"l{i}" for i in range(5)],
            higher_tree_cell=[f"h{i}" for i in range(5)], tree_gap_pp=5.,
            weight=date_weights(pairs),
        )
        pairs["weighted_contribution_c"] = pairs.weight * pairs.lst_gap_c
        base = matched_effect(pairs)
        outcomes = influence(pairs, base)
        omission = outcomes[outcomes.omission_type.eq("block") & outcomes.omitted_block_e.eq(1)].iloc[0]
        self.assertEqual(omission.omitted_pair_ids, "1 3")
        self.assertEqual(omission.dates, 3)
        self.assertEqual(omission.blocks, 3)
        self.assertAlmostEqual(omission.original_omitted_weight, .375)
        self.assertAlmostEqual(omission.estimate_c_per_5pp, -2.)
        self.assertEqual(omission.draws_retained, 1999)


if __name__ == "__main__":
    unittest.main()
