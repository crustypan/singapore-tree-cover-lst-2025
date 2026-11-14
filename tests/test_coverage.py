import unittest

import numpy as np

from audit_coverage import annual_bounds, screen_modis, temporal_eligibility


class CoverageRules(unittest.TestCase):
    def test_raw_zero_qc_is_best_quality(self):
        dn = np.array([15200, 15200, 15200, 15200, 0, -9999, 25000, 13157])
        qc = np.array([0, 16 + 64, 32, -9999, 0, 0, 0, 0])
        passed, values = screen_modis(dn, qc)
        np.testing.assert_array_equal(passed, [True, True, False, False, False, False, False, False])
        np.testing.assert_allclose(values[:2], [30.85, 30.85])
        self.assertTrue(np.isnan(values[2:]).all())

    def test_every_qa_field(self):
        codes = np.array([1, 2, 3, 4, 8, 12, 16, 32, 48, 64, 128, 192, 256])
        passed, _ = screen_modis(np.full(codes.shape, 15200), codes)
        np.testing.assert_array_equal(
            passed, [False, False, False, False, False, False, True, False, False, True, False, False, False],
        )

    def test_three_independent_time_requirements(self):
        month_days = np.array([
            [3, 3, 3, 3, 3, 3, 3, 1, 1, 1, 0, 0],  # Pass at the day threshold.
            [3, 3, 3, 3, 3, 3, 3, 1, 1, 0, 1, 0],  # 24 days, only 1 in Q4.
            [3, 3, 0, 3, 3, 0, 3, 3, 0, 3, 3, 0],  # 24 days, only 8 months.
            [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1],  # 12 months, only 23 days.
        ]).T
        month_days[10, 0] = 1
        month_days[0, 0] = 2
        days, months, quarters, kept = temporal_eligibility(month_days)
        np.testing.assert_array_equal(days, [24, 24, 24, 23])
        np.testing.assert_array_equal(months, [11, 10, 8, 12])
        np.testing.assert_array_equal(quarters[:, 0], [8, 9, 5, 2])
        np.testing.assert_array_equal(kept, [True, False, False, False])

    def test_missing_dates_cannot_fix_a_measured_quarter(self):
        months = np.zeros((12, 1, 2), dtype=int)
        months[3, 0, 0] = 2
        months[4, 0, 0] = 1
        missing = np.array([1, 5, 5, 0, 0, 0, 5, 5, 5, 5, 5, 5])
        area = np.array([[1, 9]])
        self.assertEqual(annual_bounds(months, missing, area), (0.0, 10.0, "fail"))
        self.assertEqual(annual_bounds(months, missing, np.array([[9, 1]])),
                         (0.0, 90.0, "incomplete_source"))


if __name__ == "__main__":
    unittest.main()
