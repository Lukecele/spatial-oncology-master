#!/usr/bin/env python3
"""
Master Spatial Oncology Pipeline.
Integrates the complete multi-scale continuum:
1. Biophysical Reaction-Diffusion Transport (Bessel K0 solver)
2. Evolutionary Game Theory Replicator Dynamics (Lyapunov-stable Nash equilibrium)
3. Clinical Cohort Validation on MSKCC 2020 Real Patient Cohort (N = 604)
"""

from pathlib import Path
from typing import Dict, Any, Optional
import numpy as np

from biophysical_diffusion_model import BiophysicalTransportModel
from replicator_dynamics_game import SpatialGameAssay
from clinical_validation import MSKCCClinicalValidator


class MasterSpatialOncologyPipeline:
    """Unified pipeline orchestrating the Master Study models and validations."""

    def __init__(self, data_dir: Optional[Path] = None, results_dir: Optional[Path] = None):
        root = Path(__file__).resolve().parent.parent
        self.data_dir = Path(data_dir) if data_dir else root / "data"
        self.results_dir = Path(results_dir) if results_dir else root / "results"
        self.cohort_path = self.data_dir / "mskcc_2020_real_cohort.csv"

        self.transport_model = BiophysicalTransportModel()
        self.game_assay = SpatialGameAssay()
        self.clinical_validator = MSKCCClinicalValidator(self.cohort_path) if self.cohort_path.exists() else None

    def run_biophysical_transport(self, r_nest_um: float = 200.0) -> Dict[str, Any]:
        """Runs Pillar 1: Reaction-diffusion transport in avascular nest."""
        return self.transport_model.compute_sanctuary_geometry(nest_radius_um=r_nest_um)

    def run_evolutionary_dynamics(self) -> Dict[str, Any]:
        """Runs Pillar 2: Evolutionary game theory & ADAURA cessation dynamics."""
        x_star = self.game_assay.compute_nash_equilibrium()
        stability = self.game_assay.compute_lyapunov_stability()
        regrowth = self.game_assay.simulate_adaura_cessation_bifurcation()
        return {
            "nash_equilibrium": x_star,
            "lyapunov_eigenvalue": stability["lyapunov_eigenvalue"],
            "mrd_recurrence_month": regrowth["mrd_recurrence_month"],
            "ct_imaging_recurrence_month": regrowth["ct_imaging_recurrence_month"],
            "lead_time_months": regrowth["mrd_lead_time_months"],
        }

    def run_clinical_validation(self) -> Dict[str, Any]:
        """Runs Pillar 3: MSKCC 2020 authenticated cohort validation (N = 604)."""
        if self.clinical_validator is None:
            raise FileNotFoundError(f"Missing clinical cohort at: {self.cohort_path}")
        summary = self.clinical_validator.get_cohort_summary()
        gradient = self.clinical_validator.compute_subtype_gradient()
        logrank = self.clinical_validator.compare_high_vs_low_sanctuary()
        cox = self.clinical_validator.fit_multivariable_cox_models()
        return {
            "summary": summary,
            "gradient": gradient,
            "logrank": logrank,
            "cox_models": cox,
        }

    def execute_all(self) -> Dict[str, Any]:
        """Executes the full multi-scale pipeline and prints summary."""
        print("=" * 80)
        print("  UNIFIED SPATIAL ONCOLOGY MASTER STUDY EXECUTION")
        print("  From Subcellular Reaction-Diffusion to Clinical Recurrence Dynamics")
        print("=" * 80)

        p1 = self.run_biophysical_transport()
        print(f"[Pillar 1: Transport]    rcrit = {p1['r_critical_um']:.1f} µm | Viable Shell Δr = {p1['viable_sanctuary_shell_width_um']:.1f} µm | Sanctuary Area = {p1['viable_sanctuary_area_pct']:.1f}%")

        p2 = self.run_evolutionary_dynamics()
        print(f"[Pillar 2: Game Theory]  Nash x* = {p2['nash_equilibrium']:.4f} | Lyapunov λ = {p2['lyapunov_eigenvalue']:.4f} | ctDNA Lead Time = {p2['lead_time_months']:.1f} mo")

        if self.clinical_validator:
            p3 = self.run_clinical_validation()
            m1 = p3["cox_models"]["model1_extreme"]
            print(f"[Pillar 3: Clinical N=604] High-Sanctuary Adjusted HR = {m1['adjusted_hr']:.2f} (95% CI: {m1['ci_lower']:.2f}–{m1['ci_upper']:.2f}, p = {m1['p_value']:.2e})")

        print("=" * 80)
        return {"pillar1": p1, "pillar2": p2, "pillar3": p3 if self.clinical_validator else None}


if __name__ == "__main__":
    pipeline = MasterSpatialOncologyPipeline()
    pipeline.execute_all()
