#!/usr/bin/env python3
"""
Comprehensive Automated Unit Test Suite for the Master Spatial Oncology Pipeline.
Validates the complete 4-pillar causal continuum:
  Pillar 1: Evolutionary Game Theory & Replicator Dynamics
  Pillar 2: Visium HD Reaction-Diffusion Biophysics (λ ≈ 624 µm)
  Pillar 3: Population Clinical Survival Validation (N = 854, HR = 19.45)
  Pillar 4: Routine Diagnostic H&E Attention-MIL AI (< 5 €, C-index ≈ 0.74)
"""

from pathlib import Path
import sys
import unittest
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from replicator_dynamics_game import SpatialGameAssay
from biophysical_diffusion_model import BiophysicalTransportModel
from clinical_survival_validation import compute_harrell_c_index, compute_logrank_statistic
from deep_attention_mil import (
    GatedAttentionLayer,
    AttentionSurvivalMIL,
    cox_negative_log_likelihood,
    project_spatial_attention_heatmap,
)
from master_pipeline import MasterSpatialOncologyPipeline


class TestMasterSpatialOncologyPipeline(unittest.TestCase):
    def setUp(self):
        self.pipeline = MasterSpatialOncologyPipeline()

    def test_pillar1_game_theory_nash_coexistence(self):
        assay = SpatialGameAssay()
        self.assertEqual(assay.A.shape, (2, 2))
        
        # Payoff matrix off-diagonal elements drive coexistence
        x_star = assay.compute_nash_equilibrium()
        self.assertGreater(x_star, 0.4)
        self.assertLess(x_star, 0.8)
        
        # Test replicator dynamics simulation
        traj = assay.simulate_replicator_dynamics(x0=0.5, steps=50)
        self.assertEqual(len(traj), 51)
        self.assertTrue(np.all(traj > 0.0) and np.all(traj < 1.0))
        self.assertAlmostEqual(traj[-1], x_star, delta=0.05)

    def test_pillar2_biophysical_diffusion_model(self):
        model = BiophysicalTransportModel(D=1.3e-6, k=3.3e-4)
        
        # Characteristic diffusion length lambda = sqrt(D/k) * 10000 µm ≈ 627.8 µm
        self.assertAlmostEqual(model.lambda_um, 627.6, delta=2.0)
        
        # Radial decay verification
        radii = np.array([10.0, 300.0, 627.6, 1200.0])
        conc = model.compute_concentration_profile(radii, C0=1.0, r0_um=10.0)
        self.assertAlmostEqual(conc[0], 1.0, places=4)
        self.assertLess(conc[-1], conc[0])
        
        # Pharmacological sanctuary threshold checks
        self.assertFalse(model.is_pharmacological_sanctuary(nest_radius_um=100.0, ic50_ratio=0.35))
        self.assertTrue(model.is_pharmacological_sanctuary(nest_radius_um=1200.0, ic50_ratio=0.35))

    def test_pillar3_clinical_survival_cohort_integrity(self):
        manifest_path = ROOT / "data" / "combined_filtered_massive_cohort.csv"
        self.assertTrue(manifest_path.exists(), f"Missing cohort manifest: {manifest_path}")
        
        df = pd.read_csv(manifest_path)
        self.assertEqual(len(df), 854)
        
        required_cols = [
            "patient_id", "cohort", "histology", "egfr_status",
            "follow_up_months", "progression_event", "os_months", "os_event"
        ]
        for col in required_cols:
            self.assertIn(col, df.columns)
            
        egfr_cases = df[df["egfr_status"] == "Mutant"]
        self.assertEqual(len(egfr_cases), 218)
        
        # Fast execution of Pillar 3
        p3_res = self.pipeline.run_pillar3_clinical_validation()
        self.assertEqual(p3_res["total_patients"], 854)
        self.assertEqual(p3_res["egfr_mutant_patients"], 218)
        self.assertLess(p3_res["logrank_p_value"], 1e-4)
        self.assertGreater(p3_res["hazard_ratio"], 1.5)
        self.assertGreater(p3_res["harrell_c_index"], 0.70)

    def test_pillar4_he_vision_mil_model(self):
        manifest_path = ROOT / "data" / "cohort_manifest_study4.csv"
        self.assertTrue(manifest_path.exists(), f"Missing Study 4 manifest: {manifest_path}")
        
        df = pd.read_csv(manifest_path)
        self.assertEqual(len(df), 570)
        
        # Gated Attention pooling
        layer = GatedAttentionLayer(in_features=64, hidden_dim=32, seed=123)
        dummy_h = np.random.randn(80, 64)
        z, attn = layer.forward(dummy_h)
        self.assertEqual(z.shape, (64,))
        self.assertEqual(attn.shape, (80,))
        self.assertAlmostEqual(float(np.sum(attn)), 1.0, places=5)
        self.assertTrue(np.all(attn >= 0.0))
        
        # Cox partial log-likelihood loss
        hazards = np.array([2.0, 0.5, -1.0])
        times = np.array([10.0, 25.0, 40.0])
        events = np.array([1, 1, 0])
        loss = cox_negative_log_likelihood(hazards, times, events)
        self.assertGreater(loss, 0.0)
        self.assertTrue(np.isfinite(loss))
        
        # Spatial attention projection
        coords = np.array([[10, 10], [50, 50], [90, 90]])
        weights = np.array([0.2, 0.5, 0.3])
        grid = project_spatial_attention_heatmap(coords, weights, grid_size=32)
        self.assertEqual(grid.shape, (32, 32))
        self.assertTrue(np.all(grid >= 0.0))
        self.assertTrue(np.all(np.isfinite(grid)))
        
        # Fast execution of Pillar 4
        p4_res = self.pipeline.run_pillar4_he_vision_mil()
        self.assertEqual(p4_res["total_wsi_slides"], 570)
        self.assertGreater(p4_res["attention_mil_c_index"], 0.65)
        self.assertLess(p4_res["cost_per_slide_eur"], 5.0)

    def test_full_master_pipeline_execution(self):
        summary = self.pipeline.execute_all()
        self.assertTrue(summary["causal_chain_verified"])
        self.assertIn("pillar1_evolutionary_game", summary)
        self.assertIn("pillar2_biophysical_diffusion", summary)
        self.assertIn("pillar3_clinical_population_validation", summary)
        self.assertIn("pillar4_he_digital_pathology_ai", summary)


if __name__ == "__main__":
    unittest.main()
