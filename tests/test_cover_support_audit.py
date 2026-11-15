import unittest

import numpy as np

from analyze_cover_support_audit import DW_FIELDS, ESRI_FIELDS, largest_patch, pixel_cover


class CoverSupportAudit(unittest.TestCase):
    def test_cover_shares_use_same_dry_pixels_and_retain_water_diagnostics(self):
        n = 10
        index = np.zeros(n, dtype=int)
        esri = np.array([2, 2, 2, 7, 7, 7, 11, 11, 1, 7])
        dw = np.array([1, 2, 2, 6, 6, 6, 2, 6, 6, 0])
        accepted = np.ones(n, dtype=bool)
        dry = (esri != 1) & (dw != 0)
        r = pixel_cover([esri, dw], [ESRI_FIELDS, DW_FIELDS], ["cell"], index, np.array([n]),
                        2023, 0.60, "fixed", accepted, dry).iloc[0]
        self.assertTrue(r.qualified)
        self.assertEqual(r.common_nonwater_centers, 8)
        self.assertAlmostEqual(r.esri_tree_share, 3 / 8)
        self.assertAlmostEqual(r.dw_built_share, 4 / 8)
        self.assertAlmostEqual(r.esri_water_share, 1 / 10)
        self.assertAlmostEqual(r.dw_water_share, 1 / 10)

    def test_evidence_patch_connects_adjacent_change_pixels(self):
        mask = np.zeros((7, 7), dtype=bool)
        mask[:2, :2] = True
        mask[2, 2] = True
        mask[5:, 5:] = True
        self.assertEqual(largest_patch(mask), 5)
        self.assertEqual(largest_patch(np.zeros((7, 7), dtype=bool)), 0)
