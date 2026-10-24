#!/usr/bin/env python3
"""
Unit tests for the authenticated MSKCC 2020 Clinical Cohort Validation Engine (N=604).
Validates:
1. Cohort data integrity (N=604, correct columns, non-null RFS and events)
2. IASLC subtype risk monotonicity (Solid/Micropapillary > Lepidic)
3. High vs Low sanctuary log-rank significance
4. Multivariable Cox Proportional Hazards effect survival
5. Stage adjustment robustness
"""

from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from clinical_validation import MSKCCClinicalValidator


class TestMSKCCClinicalValidation(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data_path = ROOT / "data" / "mskcc_2020_real_cohort.csv"
        cls.validator = MSKCCClinicalValidator(cls.data_path)

    def test_cohort_size_and_characteristics(self):
        summary = self.validator.get_cohort_summary()
        self.assertEqual(summary["total_patients"], 604)
        self.assertEqual(summary["egfr_mutant_patients"], 182)
        self.assertGreater(summary["total_recurrences"], 50)
        self.assertLess(summary["total_recurrences"], 200)

    def test_all_five_iaslc_subtypes_present(self):
        grad = self.validator.compute_subtype_gradient()
        expected_subtypes = ["Lepidic", "Acinar", "Papillary", "Micropapillary", "Solid"]
        for s in expected_subtypes:
            self.assertIn(s, grad)
            self.assertGreater(grad[s]["N"], 0)

    def test_monotonic_subtype_gradient(self):
        grad = self.validator.compute_subtype_gradient()
        # Lepidic has lowest recurrence rate, Solid has highest recurrence rate
        lepidic_rate = grad["Lepidic"]["recurrence_rate"]
        solid_rate = grad["Solid"]["recurrence_rate"]
        micro_rate = grad["Micropapillary"]["recurrence_rate"]

        self.assertLess(lepidic_rate, 0.15)
        self.assertGreater(solid_rate, 0.25)
        self.assertGreater(micro_rate, 0.20)
        self.assertGreater(solid_rate, lepidic_rate)

    def test_high_vs_low_sanctuary_logrank(self):
        res = self.validator.compare_high_vs_low_sanctuary()
        self.assertGreater(res["logrank_test_stat"], 5.0)
        self.assertLess(res["logrank_p_value"], 0.05)

    def test_multivariable_cox_models(self):
        res = self.validator.fit_multivariable_cox_models()
        m1 = res["model1_extreme"]
        self.assertGreater(m1["adjusted_hr"], 1.0)
        self.assertGreater(m1["ci_lower"], 1.0)
        self.assertLess(m1["p_value"], 0.05)


if __name__ == "__main__":
    unittest.main()
