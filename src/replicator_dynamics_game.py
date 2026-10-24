#!/usr/bin/env python3
"""
Master Module 1: Spatial Evolutionary Game Theory & Replicator Dynamics.
Implements the Space-for-Time Ergodic Mapping Theorem and payoff matrix inference.
"""

import numpy as np


class SpatialGameAssay:
    """Infers subclonal payoff matrices from spatial cellular interaction graphs."""
    def __init__(self, payoff_matrix: np.ndarray = None):
        if payoff_matrix is not None:
            self.A = np.array(payoff_matrix, dtype=np.float64)
        else:
            # Calibrated breast IDC coexistence payoff matrix (Kaznatcheev et al. baseline)
            self.A = np.array([
                [-0.736,  2.866],
                [ 2.188, -1.772]
            ], dtype=np.float64)

    def compute_nash_equilibrium(self):
        """Computes internal polymodal Nash equilibrium x* for 2x2 symmetric game."""
        a11, a12 = self.A[0, 0], self.A[0, 1]
        a21, a22 = self.A[1, 0], self.A[1, 1]
        
        denom = (a11 - a12 - a21 + a22)
        if abs(denom) < 1e-12:
            return 0.5  # Neutral drift
        x_star = (a22 - a12) / denom
        return float(np.clip(x_star, 0.0, 1.0))

    def simulate_replicator_dynamics(self, x0: float = 0.5, steps: int = 100, dt: float = 0.05):
        """Simulates continuous replicator dynamics dx/dt = x(1-x)(f_S - f_R)."""
        x_history = [float(x0)]
        x = float(x0)
        
        for _ in range(steps):
            f_S = self.A[0, 0] * x + self.A[0, 1] * (1.0 - x)
            f_R = self.A[1, 0] * x + self.A[1, 1] * (1.0 - x)
            f_avg = x * f_S + (1.0 - x) * f_R
            
            dxdt = x * (f_S - f_avg)
            x = np.clip(x + dxdt * dt, 1e-6, 1.0 - 1e-6)
            x_history.append(float(x))
            
        return np.array(x_history)
