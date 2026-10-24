#!/usr/bin/env python3
"""
Unit test suite for Evolutionary Game Theory & ADAURA Cessation Bifurcation.
Validates Nash equilibrium, Lyapunov asymptotic stability,
and post-cessation recurrence dynamics matching clinical trial findings.
"""

from pathlib import Path
import sys
import unittest
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from replicator_dynamics_game import SpatialGameAssay


class TestExactGameTheoryModel(unittest.TestCase):
    def setUp(self):
        self.game = SpatialGameAssay()

    def test_coexistence_criteria(self):
        coexists, conds = self.game.verify_coexistence_conditions()
        self.assertTrue(coexists)
        self.assertGreater(conds["a12_minus_a22"], 0.0)
        self.assertGreater(conds["a21_minus_a11"], 0.0)
        self.assertIn("Coexistence", conds["game_class"])

    def test_internal_nash_equilibrium(self):
        x_star = self.game.compute_nash_equilibrium()
        self.assertAlmostEqual(x_star, 0.6133, places=3)
        self.assertGreater(x_star, 0.4)
        self.assertLess(x_star, 0.8)

    def test_lyapunov_stability_eigenvalue(self):
        stab = self.game.compute_lyapunov_stability()
        self.assertTrue(stab["is_asymptotically_stable"])
        self.assertLess(stab["lyapunov_eigenvalue"], 0.0)
        self.assertAlmostEqual(stab["lyapunov_eigenvalue"], -1.7934, places=2)

    def test_phase_portrait_attractor_convergence(self):
        """Tests that any initial starting frequency converges to x* = 0.6133."""
        for x0 in [0.05, 0.20, 0.50, 0.80, 0.95]:
            traj = self.game.simulate_replicator_dynamics(x0=x0, steps=250, dt=0.05)
            self.assertAlmostEqual(float(traj[-1]), 0.6133, delta=0.005)

    def test_adaura_cessation_bifurcation_reproduction(self):
        """Tests that post-cessation regrowth predicts molecular recurrence within 12 months."""
        bif = self.game.simulate_adaura_cessation_bifurcation(residual_cells_n0=2000)
        mrd_month = bif["mrd_recurrence_month"]
        self.assertIsNotNone(mrd_month)
        # Clinical trial reported 58% of recurrences occurred within 12 months
        self.assertLess(mrd_month, 12.0)
        self.assertGreater(mrd_month, 6.0)
        # ctDNA lead time reported in Nature Med 2025: median 4.7 mo (95% CI 2.2 - 5.6)
        self.assertAlmostEqual(bif["mrd_lead_time_months"], 5.5, delta=1.5)


if __name__ == "__main__":
    unittest.main()
