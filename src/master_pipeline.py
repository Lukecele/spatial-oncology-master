#!/usr/bin/env python3
"""
Master Spatial Oncology Pipeline.
Integrates the complete multi-scale continuum:
1. Biophysical Reaction-Diffusion Transport (Bessel K0 solver)
2. Evolutionary Game Theory Replicator Dynamics (Lyapunov-stable Nash equilibrium)
3. Clinical Cohort Validation on MSKCC 2020 Real Patient Cohort (N = 604)
"""

import sys
from pathlib import Path
from typing import Dict, Any, Optional
import numpy as np

BASE_DIR = Path(__file__).resolve().parent
SRC_DIR = BASE_DIR / "spatial-oncology-master" / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))
if str(BASE_DIR) not in sys.path:
    sys.path.insert(0, str(BASE_DIR))

try:
    from biophysical_diffusion_model import BiophysicalTransportModel
    from replicator_dynamics_game import SpatialGameAssay
    from clinical_validation import MSKCCClinicalValidator
except ImportError:
    from CODICE_MODELLO_BIOFISICO_BESSEL import BiophysicalTransportModel
    from CODICE_TEORIA_DEI_GIOCHI_REPLICATORE import SpatialGameAssay
    from CODICE_VALIDAZIONE_CLINICA_MSKCC import MSKCCClinicalValidator


class MasterSpatialOncologyPipeline:
    """Unified pipeline orchestrating the Master Study models and validations."""

    def __init__(self, data_dir: Optional[Path] = None, results_dir: Optional[Path] = None):
        self.base_dir = BASE_DIR
        self.data_dir = Path(data_dir) if data_dir else self.base_dir
        self.results_dir = Path(results_dir) if results_dir else self.base_dir
        
        # Robust resolution of cohort path across different working directories & CI runners
        candidates = [
            self.base_dir.parent / "data" / "mskcc_2020_real_cohort.csv",
            self.base_dir / "data" / "mskcc_2020_real_cohort.csv",
            self.base_dir / "06_DATI_CLINICI_MSKCC_604_PAZIENTI.csv",
            Path.cwd() / "data" / "mskcc_2020_real_cohort.csv",
            Path.cwd() / "06_DATI_CLINICI_MSKCC_604_PAZIENTI.csv",
        ]
        if len(self.base_dir.parents) > 1:
            candidates.append(self.base_dir.parents[1] / "06_DATI_CLINICI_MSKCC_604_PAZIENTI.csv")
        self.cohort_path = next((p for p in candidates if p is not None and p.exists()), None)

        self.transport_model = BiophysicalTransportModel()
        self.game_assay = SpatialGameAssay()
        self.clinical_validator = MSKCCClinicalValidator(self.cohort_path) if self.cohort_path else None

    def run_biophysical_transport(self, r_nest_um: float = 200.0) -> Dict[str, Any]:
        """Runs Pillar 1: Reaction-diffusion transport in avascular nest."""
        return self.transport_model.compute_sanctuary_geometry(nest_radius_um=r_nest_um)

    def run_evolutionary_dynamics(self) -> Dict[str, Any]:
        """Runs Pillar 2: Evolutionary game theory & ADAURA cessation dynamics."""
        x_star = self.game_assay.compute_nash_equilibrium()
        stability = self.game_assay.compute_lyapunov_stability()
        regrowth = self.game_assay.simulate_adaura_cessation_bifurcation()
        lead_decomp = self.game_assay.compute_lead_time_decomposition()
        return {
            "nash_equilibrium": x_star,
            "lyapunov_eigenvalue": stability["lyapunov_eigenvalue"],
            "mrd_recurrence_month": regrowth["mrd_recurrence_month"],
            "ct_imaging_recurrence_month": regrowth["ct_imaging_recurrence_month"],
            "lead_time_months": regrowth["mrd_lead_time_months"],
            "analytical_lead_time_months": lead_decomp["lead_time_months"]
        }

    def run_clinical_validation(self) -> Dict[str, Any]:
        """Runs Pillar 3: MSKCC 2020 authenticated cohort validation (N = 604)."""
        if not self.clinical_validator:
            return {"error": "Clinical cohort CSV not found"}
        summary = self.clinical_validator.get_cohort_summary()
        followup = self.clinical_validator.compute_followup_and_survival_metrics()
        cox_models = self.clinical_validator.fit_multivariable_cox_models()
        return {
            "cohort_summary": summary,
            "followup_metrics": followup,
            "cox_models": cox_models,
        }

    def execute_all(self) -> Dict[str, Any]:
        """Executes the complete end-to-end analytical pipeline."""
        print("Executing Pillar 1: Reaction-Diffusion Biophysical Transport...")
        p1 = self.run_biophysical_transport()
        print("Executing Pillar 2: Evolutionary Game Theory & Post-Cessation Regrowth...")
        p2 = self.run_evolutionary_dynamics()
        print("Executing Pillar 3: Clinical Cohort Audit & Survival Modeling...")
        p3 = self.run_clinical_validation()
        return {
            "pillar1_transport": p1,
            "pillar2_evolution": p2,
            "pillar3_clinical": p3
        }


if __name__ == "__main__":
    pipeline = MasterSpatialOncologyPipeline()
    results = pipeline.execute_all()
    print("\n" + "=" * 70)
    print("  MASTER PIPELINE EXECUTION SUCCESSFUL")
    print("=" * 70)
    print(f"Pillar 1 - Sanctuary Shell Width : {results['pillar1_transport']['viable_sanctuary_shell_width_um']:.1f} um")
    print(f"Pillar 2 - Stable Nash Fraction  : {results['pillar2_evolution']['nash_equilibrium']*100:.2f}% sensitive")
    print(f"Pillar 2 - Diagnostic Lead Time  : {results['pillar2_evolution']['lead_time_months']} months")
    print(f"Pillar 3 - MSKCC Total Patients  : {results['pillar3_clinical']['cohort_summary']['total_patients']}")
    print(f"Pillar 3 - Extreme Contrast HR   : {results['pillar3_clinical']['cox_models']['model1_extreme']['adjusted_hr']:.2f} (p = {results['pillar3_clinical']['cox_models']['model1_extreme']['p_value']:.4e})")
    print(f"Pillar 3 - Full Cohort HR        : {results['pillar3_clinical']['cox_models']['model2_full_cohort']['adjusted_hr']:.2f} (p = {results['pillar3_clinical']['cox_models']['model2_full_cohort']['p_value']:.4f})")
