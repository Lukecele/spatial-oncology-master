#!/usr/bin/env python3
"""
Master Module 2: Biophysical Reaction-Diffusion Transport Solver.
Computes the characteristic penetration length lambda = sqrt(D/k) ≈ 624 µm
for small-molecule targeted therapies in solid neoplastic sheets.
"""

import numpy as np


class BiophysicalTransportModel:
    """Solves steady-state 2D radial reaction-diffusion equation for small-molecule inhibitors."""
    def __init__(self, D: float = 1.3e-6, k: float = 3.3e-4):
        """
        D: effective tissue diffusion coefficient (cm^2/s)
        k: intracellular drug metabolic clearance rate (s^-1)
        """
        self.D = float(D)
        self.k = float(k)
        # Characteristic diffusion length lambda = sqrt(D / k) in cm -> convert to micrometers
        self.lambda_cm = np.sqrt(self.D / self.k)
        self.lambda_um = self.lambda_cm * 10000.0  # ≈ 627.8 µm

    def compute_concentration_profile(self, radii_um: np.ndarray, C0: float = 1.0, r0_um: float = 10.0):
        """
        Computes radial drug concentration C(r) normalized to vascular source C0.
        Approximates modified Bessel decay: C(r) ≈ C0 * exp(-(r - r0) / lambda)
        """
        radii_clean = np.maximum(radii_um, r0_um)
        decay = np.exp(-(radii_clean - r0_um) / self.lambda_um)
        return float(C0) * decay

    def is_pharmacological_sanctuary(self, nest_radius_um: float, ic50_ratio: float = 0.35) -> bool:
        """Determines if the core of a tumor nest falls below therapeutic IC50 threshold."""
        c_core = self.compute_concentration_profile(np.array([nest_radius_um]))[0]
        return bool(c_core < ic50_ratio)
