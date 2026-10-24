#!/usr/bin/env python3
"""
Unit test suite for the exact biophysical reaction-diffusion transport solver.
Validates Bessel analytical equations, ODE numerical convergence,
and TKI vs macromolecular ADC sanctuary compartmentalization.
"""

import math
from pathlib import Path
import sys
import unittest
import numpy as np
import scipy.special as sp
from scipy.integrate import odeint

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from biophysical_diffusion_model import BiophysicalTransportModel, get_standard_modalities


class TestExactBiophysicalModel(unittest.TestCase):
    def setUp(self):
        self.models = get_standard_modalities()
        self.tki = self.models["tki_osimertinib"]
        self.adc = self.models["adc_macromolecule"]

    def test_tki_characteristic_screening_length(self):
        expected_lambda = math.sqrt(1.3e-6 / 3.3e-4) * 1e4
        self.assertAlmostEqual(self.tki.lambda_um, expected_lambda, places=2)
        self.assertAlmostEqual(self.tki.lambda_um, 627.65, delta=0.5)

    def test_bessel_analytical_matches_numerical_ode(self):
        """Verifies exact analytical Bessel solution against numerical ODE integration."""
        r0 = self.tki.r0_um
        lam = self.tki.lambda_um
        radii = np.linspace(r0, 600.0, 100)
        
        def deriv(y, r):
            C, dCdr = y
            return [dCdr, -(1.0 / r) * dCdr + C / (lam ** 2)]
            
        dCdr_0 = -sp.k1(r0 / lam) / (lam * sp.k0(r0 / lam))
        sol = odeint(deriv, [1.0, dCdr_0], radii)
        C_num = sol[:, 0]
        C_ana = self.tki.compute_concentration_profile(radii, C0=1.0)
        
        max_err = np.max(np.abs(C_num - C_ana))
        self.assertLess(max_err, 1e-4)

    def test_critical_radius_thresholds(self):
        """Osimertinib crosses IC50 (0.35 C0) at 165.9 um; ADC crosses at 47.4 um."""
        r_crit_tki = self.tki.get_critical_radius_um(ic_threshold_ratio=0.35)
        self.assertAlmostEqual(r_crit_tki, 165.9, delta=1.0)
        
        r_crit_adc = self.adc.get_critical_radius_um(ic_threshold_ratio=0.35)
        self.assertAlmostEqual(r_crit_adc, 47.4, delta=1.0)

    def test_viable_sanctuary_shell_geometry(self):
        """Confirms viable sanctuary shell between critical radius and Krogh oxygen limit (200 um)."""
        geom_tki_100 = self.tki.compute_sanctuary_geometry(nest_radius_um=100.0, ic50_ratio=0.35)
        self.assertEqual(geom_tki_100["operational_regime"], "Full Therapeutic Suppression")
        self.assertEqual(geom_tki_100["viable_sanctuary_area_pct"], 0.0)
        
        geom_tki_200 = self.tki.compute_sanctuary_geometry(nest_radius_um=200.0, ic50_ratio=0.35)
        self.assertEqual(geom_tki_200["operational_regime"], "Viable Pharmacological Sanctuary")
        self.assertAlmostEqual(geom_tki_200["viable_sanctuary_shell_width_um"], 34.1, delta=1.0)
        self.assertGreater(geom_tki_200["viable_sanctuary_area_pct"], 30.0)
        
        geom_tki_500 = self.tki.compute_sanctuary_geometry(nest_radius_um=500.0, ic50_ratio=0.35)
        self.assertEqual(geom_tki_500["operational_regime"], "Viable Sanctuary Shell with Hypoxic/Necrotic Core")
        self.assertGreater(geom_tki_500["necrotic_core_area_pct"], 80.0)

    def test_adc_severe_transport_restriction(self):
        """ADCs (150 kDa) leave >75% of viable nest uninhibited at R = 100 um."""
        geom_adc_100 = self.adc.compute_sanctuary_geometry(nest_radius_um=100.0, ic50_ratio=0.35)
        self.assertGreater(geom_adc_100["viable_sanctuary_area_pct"], 75.0)


if __name__ == "__main__":
    unittest.main()
