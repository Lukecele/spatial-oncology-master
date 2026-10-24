#!/usr/bin/env python3
"""
Master Module 2: Exact Biophysical Reaction-Diffusion Transport Solver.
Solves the 2D radial steady-state Helmholtz/Bessel diffusion equation in Krogh cylindrical coordinates:
    (1/r) d/dr [ r dC/dr ] - (k / D) C = 0
Exact analytical solution:
    C(r) = C0 * [ K0(r / lambda) / K0(r0 / lambda) ]
where:
    lambda = sqrt(D / k) is the characteristic screening/penetration length.
"""

from typing import Union, Dict, Any
import numpy as np
import scipy.special as sp
from scipy.optimize import brentq


class BiophysicalTransportModel:
    """
    Rigorously solves radial reaction-diffusion transport for therapeutic agents
    in compact neoplastic cellular architectures.
    """
    def __init__(
        self,
        D: float = 1.3e-6,      # cm^2 / s (effective tissue diffusion coefficient)
        k: float = 3.3e-4,      # s^-1 (intracellular metabolic/consumption rate)
        r0_um: float = 10.0,    # um (capillary vascular source radius)
        name: str = "Small-Molecule TKI (Osimertinib 500 Da)"
    ):
        self.D = float(D)
        self.k = float(k)
        self.r0_um = float(r0_um)
        self.name = str(name)
        
        # Characteristic diffusion length lambda = sqrt(D / k) in cm -> convert to micrometers
        self.lambda_cm = np.sqrt(self.D / self.k)
        self.lambda_um = self.lambda_cm * 10000.0  # e.g., 627.65 um for Osimertinib
        
        # Physiological Krogh oxygen cylinder limit (150-200 um)
        self.krogh_oxygen_limit_um = 200.0

    def compute_concentration_profile(
        self,
        radii_um: Union[float, np.ndarray],
        C0: float = 1.0,
        r0_um: float = None
    ) -> Union[float, np.ndarray]:
        """
        Computes exact normalized concentration C(r)/C0 using modified Bessel function K0.
        Exact solution to cylindrical reaction-diffusion boundary value problem:
            C(r0) = C0,  C(r -> inf) -> 0
        """
        eff_r0 = float(r0_um) if r0_um is not None else self.r0_um
        r = np.asarray(radii_um, dtype=np.float64)
        r_clamped = np.maximum(r, eff_r0)
        
        k0_r = sp.k0(r_clamped / self.lambda_um)
        k0_r0 = sp.k0(eff_r0 / self.lambda_um)
        
        profile = float(C0) * (k0_r / k0_r0)
        return float(profile) if np.ndim(radii_um) == 0 else profile

    def get_critical_radius_um(self, ic_threshold_ratio: float = 0.35) -> float:
        """
        Finds the exact critical radius r_crit where drug concentration drops below threshold.
        Default threshold: IC50 = 0.35 * C0.
        """
        def f(r):
            return self.compute_concentration_profile(r, C0=1.0) - ic_threshold_ratio
        
        # Search between capillary wall and 5000 um
        try:
            r_crit = brentq(f, self.r0_um + 1e-4, 5000.0)
            return float(r_crit)
        except ValueError:
            return float(self.r0_um)

    def is_pharmacological_sanctuary(self, nest_radius_um: float, ic50_ratio: float = 0.35) -> bool:
        """
        Determines whether an avascular nest of radius R creates a pharmacological sanctuary core
        where the drug concentration drops below the inhibitory threshold (e.g. IC50 = 0.35 * C0).
        """
        r_crit = self.get_critical_radius_um(ic50_ratio)
        return float(nest_radius_um) > r_crit

    def compute_sanctuary_geometry(self, nest_radius_um: float, ic50_ratio: float = 0.35) -> Dict[str, Any]:
        """
        Evaluates the compartmental geometry of a tumor nest of radius R:
          1. Core drug concentration C(R)
          2. Whether core falls below IC50 (pharmacological sanctuary)
          3. Viable sanctuary shell fraction (between r_crit and min(R, 200 um))
          4. Necrotic core fraction (beyond Krogh limit 200 um)
        """
        R = float(nest_radius_um)
        c_core = float(self.compute_concentration_profile(R, C0=1.0))
        r_crit = self.get_critical_radius_um(ic50_ratio)
        r_ox = self.krogh_oxygen_limit_um
        
        total_area = np.pi * (R**2)
        
        if R <= r_crit:
            # Full therapeutic suppression
            viable_sanctuary_area = 0.0
            necrotic_area = 0.0
            regime = "Full Therapeutic Suppression"
        elif R <= r_ox:
            # Entire post-critical region is viable sanctuary
            viable_sanctuary_area = np.pi * (R**2 - r_crit**2)
            necrotic_area = 0.0
            regime = "Viable Pharmacological Sanctuary"
        else:
            # Viable sanctuary is the annulus [r_crit, 200 um]; center is necrotic
            viable_sanctuary_area = np.pi * (r_ox**2 - r_crit**2)
            necrotic_area = np.pi * (R**2 - r_ox**2)
            regime = "Viable Sanctuary Shell with Hypoxic/Necrotic Core"
            
        viable_sanctuary_pct = (viable_sanctuary_area / total_area) * 100.0
        necrotic_pct = (necrotic_area / total_area) * 100.0
        
        return {
            "nest_radius_um": R,
            "core_concentration_ratio": c_core,
            "r_critical_um": r_crit,
            "krogh_oxygen_limit_um": r_ox,
            "viable_sanctuary_shell_width_um": max(0.0, min(R, r_ox) - r_crit),
            "viable_sanctuary_area_pct": viable_sanctuary_pct,
            "necrotic_core_area_pct": necrotic_pct,
            "operational_regime": regime
        }


def get_standard_modalities() -> Dict[str, BiophysicalTransportModel]:
    """Returns calibrated biophysical transport models for TKI and ADC therapeutic classes."""
    return {
        "tki_osimertinib": BiophysicalTransportModel(
            D=1.3e-6,
            k=3.3e-4,
            r0_um=10.0,
            name="Small-Molecule TKI (Osimertinib 500 Da)"
        ),
        "adc_macromolecule": BiophysicalTransportModel(
            D=1.5e-8,
            k=2.5e-4,
            r0_um=10.0,
            name="Antibody-Drug Conjugate / Bispecific (HER3-DXd / Amivantamab 150 kDa)"
        )
    }


if __name__ == "__main__":
    mods = get_standard_modalities()
    print("=" * 70)
    print("  EXACT BESSEL REACTION-DIFFUSION TRANSPORT BENCHMARK")
    print("=" * 70)
    for key, model in mods.items():
        print(f"\nModality: {model.name}")
        print(f"  Diffusion coeff D : {model.D:.2e} cm^2/s")
        print(f"  Clearance rate k  : {model.k:.2e} s^-1")
        print(f"  Screening length  : {model.lambda_um:.2f} um")
        rc = model.get_critical_radius_um(0.35)
        print(f"  Critical r (IC50) : {rc:.1f} um")
        for nest_r in [100, 200, 500]:
            geom = model.compute_sanctuary_geometry(nest_r, 0.35)
            print(f"    Nest R = {nest_r:3d} um -> Core C/C0: {geom['core_concentration_ratio']:.4f} | Sanctuary Area: {geom['viable_sanctuary_area_pct']:.1f}% | Regime: {geom['operational_regime']}")
