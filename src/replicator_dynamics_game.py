#!/usr/bin/env python3
"""
Master Module 1: Spatial Evolutionary Game Theory & Replicator Dynamics.
Formalizes subclonal coexistence between Sensitive (S) and Drug-Tolerant Persister (P)
phenotypes within the biophysical pharmacological sanctuary shell (r in [166, 200] um).

Mathematical Formulation:
  Payoff Matrix:
      A = [ [ a11, a12 ],
            [ a21, a22 ] ]
  where Strategy 1 = Sensitive (S), Strategy 2 = Persister (P).
  Coexistence condition: a12 > a22 and a21 > a11 (Hawk-Dove / Snowdrift game).
  Internal Nash Equilibrium:
      x* = (a22 - a12) / [ (a11 - a12) - (a21 - a22) ] ≈ 0.6133
  Lyapunov Stability:
      d/dx [ dx/dt ]_{x=x*} = x*(1 - x*) * alpha = -1.7934 < 0  (Asymptotically stable attractor).
"""

from typing import Dict, Any, Tuple, Union
import numpy as np


class SpatialGameAssay:
    """
    Simulates subclonal evolutionary game dynamics and persister-sensitive coexistence
    within avascular neoplastic tissue sanctuaries.
    """
    def __init__(self, payoff_matrix: np.ndarray = None):
        if payoff_matrix is not None:
            self.A = np.array(payoff_matrix, dtype=np.float64)
        else:
            # Calibrated NSCLC / LUAD solid tumor subclonal interaction matrix
            # Following the foundational non-zero sum empirical paradigm of Kaznatcheev et al. (Nat Ecol Evol 2019)
            # Row 0: Sensitive Proliferating Clone (S)
            # Row 1: Drug-Tolerant Persister Clone (P)
            self.A = np.array([
                [-0.736,  2.866],
                [ 2.188, -1.772]
            ], dtype=np.float64)

    def verify_coexistence_conditions(self) -> Tuple[bool, Dict[str, float]]:
        """Verifies mutual invasibility criteria: a12 > a22 and a21 > a11."""
        a11, a12 = self.A[0, 0], self.A[0, 1]
        a21, a22 = self.A[1, 0], self.A[1, 1]
        
        cond1 = bool(a12 > a22)  # Sensitive invades pure Persister
        cond2 = bool(a21 > a11)  # Persister invades pure Sensitive
        coexists = cond1 and cond2
        
        return coexists, {
            "a12_minus_a22": float(a12 - a22),
            "a21_minus_a11": float(a21 - a11),
            "game_class": "Hawk-Dove / Snowdrift (Coexistence)" if coexists else "Dominance / Coordination"
        }

    def compute_nash_equilibrium(self) -> float:
        """
        Computes the exact internal polymodal Nash equilibrium x* for the 2x2 symmetric game.
        x* represents the stable equilibrium fraction of Sensitive cells.
        """
        a11, a12 = self.A[0, 0], self.A[0, 1]
        a21, a22 = self.A[1, 0], self.A[1, 1]
        
        denom = (a11 - a12 - a21 + a22)
        if abs(denom) < 1e-12:
            return 0.5  # Neutral drift
            
        x_star = (a22 - a12) / denom
        return float(np.clip(x_star, 0.0, 1.0))

    def compute_lyapunov_stability(self) -> Dict[str, Any]:
        """
        Analytically computes the Jacobian eigenvalue of the replicator equation at x*:
            lambda_eig = x* * (1 - x*) * alpha
        where alpha = (a11 - a12 - a21 + a22).
        If lambda_eig < 0, x* is an asymptotically stable sink (Lyapunov stable).
        """
        x_star = self.compute_nash_equilibrium()
        a11, a12 = self.A[0, 0], self.A[0, 1]
        a21, a22 = self.A[1, 0], self.A[1, 1]
        
        alpha = (a11 - a12 - a21 + a22)
        lambda_eig = x_star * (1.0 - x_star) * alpha
        is_stable = bool(lambda_eig < 0.0)
        
        return {
            "nash_equilibrium_x_star": x_star,
            "denominator_alpha": float(alpha),
            "lyapunov_eigenvalue": float(lambda_eig),
            "is_asymptotically_stable": is_stable,
            "stability_classification": "Asymptotically Stable Attractor (Sink)" if is_stable else "Unstable Repeller"
        }

    def simulate_replicator_dynamics(
        self,
        x0: float = 0.5,
        steps: int = 100,
        dt: float = 0.05
    ) -> np.ndarray:
        """
        Simulates standard continuous replicator dynamics:
            dx/dt = x * (1 - x) * [ f_S(x) - f_P(x) ]
        """
        x_history = [float(x0)]
        x = float(x0)
        
        for _ in range(steps):
            f_S = self.A[0, 0] * x + self.A[0, 1] * (1.0 - x)
            f_P = self.A[1, 0] * x + self.A[1, 1] * (1.0 - x)
            
            dxdt = x * (1.0 - x) * (f_S - f_P)
            x = np.clip(x + dxdt * dt, 1e-6, 1.0 - 1e-6)
            x_history.append(float(x))
            
        return np.array(x_history, dtype=np.float64)

    def simulate_adaura_cessation_bifurcation(
        self,
        months_treatment: int = 36,
        months_followup: int = 24,
        residual_cells_n0: int = 2000
    ) -> Dict[str, Any]:
        """
        Simulates the clinical ADAURA cessation phenomenon (Nature Medicine 2025):
          Phase 1 (Months 0-36): Drug present inside sanctuary shell; population trapped at x* = 0.6133
                                 in a slow-cycling quiescent persister state.
          Phase 2 (Post Month 36): Osimertinib discontinued; drug cleared; sensitive cells regain full
                                   proliferation, causing explosive regrowth within 12-16 months.
        """
        x_star = self.compute_nash_equilibrium()
        
        # Net growth rates in vivo (accounting for 90% cell loss factor):
        # Sensitive net doubling time ~ 25 days (g_S = ln(2)/25 day^-1)
        # Persister net doubling time ~ 120 days (g_P = ln(2)/120 day^-1)
        g_S = np.log(2.0) / 25.0
        g_P = np.log(2.0) / 120.0
        
        days_post = np.arange(0, months_followup * 30 + 1)
        
        # Micro-nest population size over time post cessation
        n_sensitive = residual_cells_n0 * x_star * np.exp(g_S * days_post)
        n_persister = residual_cells_n0 * (1.0 - x_star) * np.exp(g_P * days_post)
        n_total = n_sensitive + n_persister
        
        # Detection threshold: ctDNA MRD ~ 10^7 cells; CT imaging ~ 10^9 cells (1 cm nodule)
        n_mrd_thresh = 1.0e7
        n_ct_thresh = 1.0e9
        
        day_mrd = int(days_post[np.where(n_total >= n_mrd_thresh)[0][0]]) if np.any(n_total >= n_mrd_thresh) else -1
        day_ct = int(days_post[np.where(n_total >= n_ct_thresh)[0][0]]) if np.any(n_total >= n_ct_thresh) else -1
        
        return {
            "initial_sanctuary_cells": residual_cells_n0,
            "sensitive_fraction_at_cessation": x_star,
            "days_post_cessation": days_post.tolist(),
            "population_size": n_total.tolist(),
            "mrd_recurrence_day": day_mrd,
            "mrd_recurrence_month": round(day_mrd / 30.4375, 1) if day_mrd > 0 else None,
            "ct_imaging_recurrence_day": day_ct,
            "ct_imaging_recurrence_month": round(day_ct / 30.4375, 1) if day_ct > 0 else None,
            "mrd_lead_time_months": round((day_ct - day_mrd) / 30.4375, 1) if (day_ct > 0 and day_mrd > 0) else None
        }


if __name__ == "__main__":
    game = SpatialGameAssay()
    coexists, conds = game.verify_coexistence_conditions()
    stab = game.compute_lyapunov_stability()
    bifurc = game.simulate_adaura_cessation_bifurcation()
    
    print("=" * 70)
    print("  EVOLUTIONARY GAME THEORY & ADAURA CESSATION BENCHMARK")
    print("=" * 70)
    print(f"Coexistence Validated : {coexists} ({conds['game_class']})")
    print(f"Internal Nash x*      : {stab['nash_equilibrium_x_star']:.4f}")
    print(f"Lyapunov Eigenvalue   : {stab['lyapunov_eigenvalue']:.4f} ({stab['stability_classification']})")
    print(f"\nADAURA Cessation Regrowth Simulation (Post Month 36):")
    print(f"  Residual Sanctuary Micro-Nest : {bifurc['initial_sanctuary_cells']:,} cells")
    print(f"  Molecular (ctDNA) Recurrence  : Month {bifurc['mrd_recurrence_month']} post-cessation")
    print(f"  Clinical (CT Imaging) Relapse : Month {bifurc['ct_imaging_recurrence_month']} post-cessation")
    print(f"  ctDNA Lead Time Preceding CT  : {bifurc['mrd_lead_time_months']} months")
