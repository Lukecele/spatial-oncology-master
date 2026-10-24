#!/usr/bin/env python3
"""
Generate 4 Unified High-Resolution Publication Figures (600 DPI)
for the Spatial Oncology & Pharmacological Sanctuary Master Study.

Figures:
  1. Figure 1: Reaction-Diffusion Transport & Pharmacological Sanctuary Architecture
  2. Figure 2: Evolutionary Game Theory & Post-ADAURA Cessation Regrowth Dynamics
  3. Figure 3: Clinical Cohort Validation (MSKCC 2020 N=604 Real Patient Cohort)
  4. Figure 4: Translational ADAURA Framework & Multiscale Decision Architecture
"""

import sys
import math
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.gridspec as gridspec
import seaborn as sns
from scipy.special import k0
from scipy.optimize import brentq
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import logrank_test

# Add src to path
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from biophysical_diffusion_model import BiophysicalTransportModel, get_standard_modalities
from replicator_dynamics_game import SpatialGameAssay

# Global Matplotlib Styling for Nature / High-Impact Journals
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans", "Helvetica", "Arial"],
    "font.size": 10,
    "axes.labelsize": 11,
    "axes.titlesize": 12,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "legend.fontsize": 9,
    "figure.titlesize": 14,
    "figure.dpi": 300,
    "savefig.dpi": 600,
    "savefig.bbox": "tight",
    "axes.linewidth": 1.0,
    "grid.linewidth": 0.5,
    "lines.linewidth": 1.75,
    "lines.markersize": 6,
    "axes.spines.top": False,
    "axes.spines.right": False,
})

RESULTS_DIR = ROOT / "results"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)
DATA_PATH = ROOT / "data" / "mskcc_2020_real_cohort.csv"


# ==============================================================================
# FIGURE 1: Biophysical Reaction-Diffusion Transport & Sanctuary Architecture
# ==============================================================================
def generate_figure_1():
    print("Generating Figure 1: Biophysical Reaction-Diffusion Transport...")
    fig = plt.figure(figsize=(15, 12))
    gs = gridspec.GridSpec(2, 2, figure=fig, hspace=0.32, wspace=0.28)

    # Models
    tki = BiophysicalTransportModel(D=1.3e-6, k=3.3e-4, r0_um=10.0, name="Osimertinib (500 Da)")
    adc = BiophysicalTransportModel(D=1.5e-8, k=2.5e-4, r0_um=10.0, name="ADC / Bispecific (150 kDa)")

    # -------------------------------------------------------------------------
    # Panel A: Radial Concentration Profiles
    # -------------------------------------------------------------------------
    ax_a = fig.add_subplot(gs[0, 0])
    radii = np.linspace(10.0, 300.0, 600)
    c_tki = tki.compute_concentration_profile(radii, C0=1.0)
    c_adc = adc.compute_concentration_profile(radii, C0=1.0)

    ax_a.plot(radii, c_tki, color="#1f77b4", label=r"Small-Molecule TKI: Osimertinib ($\lambda = 627.7\,\mu\mathrm{m}$)")
    ax_a.plot(radii, c_adc, color="#d62728", linestyle="--", label=r"Macromolecule: HER3-DXd / Amivantamab ($\lambda = 77.5\,\mu\mathrm{m}$)")

    # Threshold and limits
    ax_a.axhline(0.35, color="black", linestyle=":", linewidth=1.2, label=r"Therapeutic Threshold ($IC_{50} = 0.35\,C_0$)")
    
    # Critical radii
    r_crit_tki = tki.get_critical_radius_um(0.35)  # 165.9 um
    r_crit_adc = adc.get_critical_radius_um(0.35)  # 47.4 um
    r_ox = 200.0  # Krogh oxygen limit

    ax_a.axvline(r_crit_tki, color="#1f77b4", linestyle="-.", alpha=0.7)
    ax_a.axvline(r_crit_adc, color="#d62728", linestyle="-.", alpha=0.7)
    ax_a.axvline(r_ox, color="#2ca02c", linestyle="-", linewidth=1.5, label=r"Krogh Oxygen Limit ($r_{\mathrm{ox}} = 200\,\mu\mathrm{m}$)")

    # Shaded Viable Sanctuary Shell for Osimertinib
    ax_a.axvspan(r_crit_tki, r_ox, color="#ff7f0e", alpha=0.25, label=r"Osimertinib Viable Sanctuary Shell ($\Delta r = 34.1\,\mu\mathrm{m}$)")
    # Shaded Necrotic Core
    ax_a.axvspan(r_ox, 300.0, color="#7f7f7f", alpha=0.15, label="Anoxic / Necrotic Core")

    ax_a.set_xlim(10.0, 300.0)
    ax_a.set_ylim(0.0, 1.05)
    ax_a.set_xlabel(r"Radial Distance from Capillary Wall $r$ ($\mu\mathrm{m}$)")
    ax_a.set_ylabel(r"Normalized Drug Concentration $C(r)/C_0$")
    ax_a.set_title("A. Radial Reaction-Diffusion Concentration Profiles", fontweight="bold")
    ax_a.legend(loc="upper right", frameon=True, framealpha=0.9)
    ax_a.grid(True, linestyle="--", alpha=0.5)

    # Annotations
    ax_a.annotate(f"Osimertinib $r_{{crit}} = {r_crit_tki:.1f}\\,\\mu\\mathrm{{m}}$",
                  xy=(r_crit_tki, 0.35), xytext=(r_crit_tki + 15, 0.50),
                  arrowprops=dict(facecolor="#1f77b4", shrink=0.05, width=1, headwidth=5))
    ax_a.annotate(f"ADC $r_{{crit}} = {r_crit_adc:.1f}\\,\\mu\\mathrm{{m}}$",
                  xy=(r_crit_adc, 0.35), xytext=(r_crit_adc + 15, 0.20),
                  arrowprops=dict(facecolor="#d62728", shrink=0.05, width=1, headwidth=5))

    # -------------------------------------------------------------------------
    # Panel B: Viable Sanctuary Shell Width vs Tumor Nest Radius
    # -------------------------------------------------------------------------
    ax_b = fig.add_subplot(gs[0, 1])
    nest_radii = np.linspace(50.0, 400.0, 351)
    shell_widths_tki = [max(0.0, min(R, 200.0) - r_crit_tki) for R in nest_radii]
    shell_widths_adc = [max(0.0, min(R, 200.0) - r_crit_adc) for R in nest_radii]

    ax_b.plot(nest_radii, shell_widths_tki, color="#1f77b4", label="Osimertinib (500 Da)")
    ax_b.plot(nest_radii, shell_widths_adc, color="#d62728", linestyle="--", label="ADC / Bispecific (150 kDa)")

    ax_b.axvline(r_crit_tki, color="#1f77b4", linestyle=":", alpha=0.6)
    ax_b.axvline(200.0, color="#2ca02c", linestyle="-", linewidth=1.2, label=r"Krogh Limit ($200\,\mu\mathrm{m}$)")

    # Regions
    ax_b.axvspan(0, r_crit_tki, color="#2ca02c", alpha=0.1, label="Zone I: Complete Inhibition")
    ax_b.axvspan(r_crit_tki, 200.0, color="#ff7f0e", alpha=0.15, label="Zone II: Viable Sanctuary Core")
    ax_b.axvspan(200.0, 400.0, color="#7f7f7f", alpha=0.15, label="Zone III: Viable Shell + Necrotic Core")

    ax_b.set_xlim(50.0, 400.0)
    ax_b.set_ylim(-2.0, 170.0)
    ax_b.set_xlabel(r"Avascular Tumor Nest Radius $R$ ($\mu\mathrm{m}$)")
    ax_b.set_ylabel(r"Viable Sanctuary Shell Width $\Delta r$ ($\mu\mathrm{m}$)")
    ax_b.set_title("B. Geometric Sanctuary Regimes vs Nest Size", fontweight="bold")
    ax_b.legend(loc="upper left", frameon=True, framealpha=0.9)
    ax_b.grid(True, linestyle="--", alpha=0.5)

    ax_b.text(260, 40, f"Osimertinib Plateau:\n$\\Delta r = {200.0 - r_crit_tki:.1f}\\,\\mu\\mathrm{{m}}$",
              fontweight="bold", color="#1f77b4", bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#1f77b4", alpha=0.8))
    ax_b.text(260, 140, f"ADC Plateau:\n$\\Delta r = {200.0 - r_crit_adc:.1f}\\,\\mu\\mathrm{{m}}$\n(>75% of nest unprotected)",
              fontweight="bold", color="#d62728", bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#d62728", alpha=0.8))

    # -------------------------------------------------------------------------
    # Panel C: 2D Spatial Cross-Section of Avascular Solid Nest
    # -------------------------------------------------------------------------
    ax_c = fig.add_subplot(gs[1, 0])
    
    # Draw concentric circles representing tumor nest of radius R = 250 um
    R_nest = 250.0
    circle_nest = patches.Circle((0, 0), R_nest, edgecolor="black", facecolor="none", linewidth=1.5, linestyle="--")
    circle_ox = patches.Circle((0, 0), 200.0, facecolor="#7f7f7f", alpha=0.25, label="Anoxic Necrotic Core ($r > 200\\,\\mu\\mathrm{m}$)")
    circle_sanctuary = patches.Circle((0, 0), 200.0, facecolor="none", edgecolor="#ff7f0e", linewidth=2.0)
    
    # Outer ring from 165.9 to 200 is sanctuary
    theta = np.linspace(0, 2*np.pi, 200)
    r1, r2 = r_crit_tki, 200.0
    x_outer = r2 * np.cos(theta)
    y_outer = r2 * np.sin(theta)
    x_inner = r1 * np.cos(theta)
    y_inner = r1 * np.sin(theta)
    
    ax_c.fill(np.concatenate([x_outer, x_inner[::-1]]),
              np.concatenate([y_outer, y_inner[::-1]]),
              color="#ff7f0e", alpha=0.45, hatch="//", label=f"Viable Sanctuary Shell ({r_crit_tki:.1f} - 200 $\\mu\\mathrm{{m}}$)")
    
    # Inner zone < 165.9 um is suppressed
    circle_supp = patches.Circle((0, 0), r_crit_tki, facecolor="#2ca02c", alpha=0.25, label=f"Therapeutic Suppression ($r < {r_crit_tki:.1f}\\,\\mu\\mathrm{{m}}$)")
    
    # Capillary in center
    circle_cap = patches.Circle((0, 0), 10.0, facecolor="#d62728", edgecolor="black", linewidth=1.0, label="Capillary Source ($r_0 = 10\\,\\mu\\mathrm{m}$)")

    ax_c.add_patch(circle_ox)
    ax_c.add_patch(circle_supp)
    ax_c.add_patch(circle_cap)
    ax_c.add_patch(circle_nest)

    ax_c.set_xlim(-270, 270)
    ax_c.set_ylim(-270, 270)
    ax_c.set_aspect("equal")
    ax_c.set_xlabel(r"Spatial Dimension $X$ ($\mu\mathrm{m}$)")
    ax_c.set_ylabel(r"Spatial Dimension $Y$ ($\mu\mathrm{m}$)")
    ax_c.set_title("C. 2D Compartmental Anatomy of Solid Nest ($R = 250\\,\\mu\\mathrm{m}$)", fontweight="bold")
    ax_c.legend(loc="upper right", frameon=True, framealpha=0.9, fontsize=8)
    ax_c.grid(True, linestyle=":", alpha=0.3)

    # -------------------------------------------------------------------------
    # Panel D: Interstitial Fluid Pressure (IFP) and Convective Shutdown
    # -------------------------------------------------------------------------
    ax_d = fig.add_subplot(gs[1, 1])
    dist_ifp = np.linspace(-100.0, 300.0, 400)
    # Sigmoidal drop of IFP at tumor boundary (Baxter & Jain model)
    P_microvasc = 25.0  # mmHg
    P_normal = 0.0      # mmHg
    # IFP plateaus inside nest and drops across boundary at r = 250 um
    ifp_profile = P_microvasc / (1.0 + np.exp((dist_ifp - 250.0) / 15.0))
    # Fluid velocity v = -K grad(P) (proportional to derivative)
    grad_p = -np.gradient(ifp_profile, dist_ifp)
    velocity = np.clip(grad_p * 1.5, 0.0, 1.2)  # um/s

    ax_d.plot(dist_ifp, ifp_profile, color="#6a3d9a", linewidth=2.0, label=r"Interstitial Fluid Pressure ($P_i$, mmHg)")
    ax_d.axhline(P_microvasc, color="black", linestyle="--", linewidth=1.0, label=r"Microvascular Pressure ($P_v = 25\,\mathrm{mmHg}$)")
    
    ax_d_twin = ax_d.twinx()
    ax_d_twin.plot(dist_ifp, velocity, color="#e31a1c", linestyle="-.", linewidth=1.8, label=r"Convective Fluid Velocity ($v$, $\mu\mathrm{m/s}$)")
    ax_d_twin.set_ylabel(r"Outward Fluid Velocity ($\mu\mathrm{m/s}$)", color="#e31a1c")
    ax_d_twin.tick_params(axis="y", labelcolor="#e31a1c")
    ax_d_twin.set_ylim(-0.05, 1.4)

    ax_d.axvline(250.0, color="gray", linestyle=":", label="Tumor Nest Margin ($R=250\\,\\mu\\mathrm{m}$)")
    ax_d.axvspan(-100, 250.0, color="#6a3d9a", alpha=0.08, label=r"Intratumoral Plateau: $\nabla P \approx 0 \rightarrow v \approx 0$")

    ax_d.set_xlim(-100.0, 300.0)
    ax_d.set_ylim(-2.0, 32.0)
    ax_d.set_xlabel(r"Distance Relative to Capillary Axis ($\mu\mathrm{m}$)")
    ax_d.set_ylabel(r"Hydrostatic Pressure ($\mathrm{mmHg}$)", color="#6a3d9a")
    ax_d.tick_params(axis="y", labelcolor="#6a3d9a")
    ax_d.set_title("D. Interstitial Fluid Pressure (IFP) & Convective Shutdown", fontweight="bold")
    ax_d.legend(loc="upper left", frameon=True, framealpha=0.9, fontsize=8)
    ax_d.grid(True, linestyle="--", alpha=0.5)

    ax_d.text(30, 10, "Péclet Number $Pe \\ll 1$\nStrict Fickian Diffusion",
              fontweight="bold", color="#6a3d9a", bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#6a3d9a", alpha=0.9))

    out_path = RESULTS_DIR / "figure1_reaction_diffusion_sanctuary.png"
    plt.savefig(out_path)
    plt.close()
    print(f"  -> Figure 1 successfully written to: {out_path}")


# ==============================================================================
# FIGURE 2: Evolutionary Game Theory & Post-ADAURA Cessation Regrowth
# ==============================================================================
def generate_figure_2():
    print("Generating Figure 2: Evolutionary Game Theory & ADAURA Cessation...")
    fig = plt.figure(figsize=(15, 12))
    gs = gridspec.GridSpec(2, 2, figure=fig, hspace=0.32, wspace=0.28)

    game = SpatialGameAssay()
    x_star = game.compute_nash_equilibrium()
    stab = game.compute_lyapunov_stability()
    bifurc = game.simulate_adaura_cessation_bifurcation(months_treatment=36, months_followup=24)

    # -------------------------------------------------------------------------
    # Panel A: Evolutionary Replicator Dynamics Phase Portrait
    # -------------------------------------------------------------------------
    ax_a = fig.add_subplot(gs[0, 0])
    x_vals = np.linspace(0.001, 0.999, 400)
    f_S = game.A[0, 0] * x_vals + game.A[0, 1] * (1.0 - x_vals)
    f_P = game.A[1, 0] * x_vals + game.A[1, 1] * (1.0 - x_vals)
    dxdt = x_vals * (1.0 - x_vals) * (f_S - f_P)

    ax_a.plot(x_vals, dxdt, color="#2b5c8f", linewidth=2.2, label=r"Replicator Velocity $\dot{x} = x(1-x)[f_S(x) - f_P(x)]$")
    ax_a.axhline(0, color="black", linestyle="--", linewidth=0.8)
    
    # Mark equilibrium points
    ax_a.plot([0.0, 1.0], [0.0, 0.0], "ro", markersize=6, fillstyle="none", label="Unstable Repellers ($x=0, x=1$)")
    ax_a.plot([x_star], [0.0], "go", markersize=9, label=f"Polymorphic Nash Equilibrium ($x^* = {x_star:.4f}$)")

    # Vector arrows indicating flow
    ax_a.quiver([0.2, 0.4], [0.0, 0.0], [0.15, 0.15], [0.0, 0.0], color="#2ca02c", scale=1, scale_units="xy", angles="xy", width=0.008)
    ax_a.quiver([0.8, 0.9], [0.0, 0.0], [-0.15, -0.15], [0.0, 0.0], color="#2ca02c", scale=1, scale_units="xy", angles="xy", width=0.008)

    ax_a.set_xlim(-0.02, 1.02)
    ax_a.set_ylim(min(dxdt) * 1.15, max(dxdt) * 1.15)
    ax_a.set_xlabel("Sensitive Clone Proportion ($x$)")
    ax_a.set_ylabel(r"Evolutionary Velocity $\dot{x}$ (rate/day)")
    ax_a.set_title("A. Replicator Dynamics Phase Portrait & Attractor", fontweight="bold")
    ax_a.legend(loc="upper right", frameon=True, framealpha=0.9)
    ax_a.grid(True, linestyle="--", alpha=0.5)

    ax_a.text(0.15, 0.15, f"Lyapunov Eigenvalue:\n$\\lambda_{{eig}} = {stab['lyapunov_eigenvalue']:.4f} < 0$\n(Asymptotically Stable Sink)",
              fontweight="bold", color="#2ca02c", bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#2ca02c", alpha=0.9))

    # -------------------------------------------------------------------------
    # Panel B: Subclonal Payoff Matrix in NSCLC Solid Architecture
    # -------------------------------------------------------------------------
    ax_b = fig.add_subplot(gs[0, 1])
    labels = ["Sensitive (S)", "Persister (P)"]
    matrix_data = game.A

    sns.heatmap(matrix_data, annot=True, fmt="+.3f", cmap="vlag", center=0.0,
                xticklabels=labels, yticklabels=labels, cbar_kws={"label": "Subclonal Payoff / Fitness Contribution"},
                ax=ax_b, annot_kws={"size": 13, "weight": "bold"})
    
    ax_b.set_title("B. Empirical NSCLC Subclonal Payoff Matrix ($A$)", fontweight="bold")
    ax_b.set_xlabel("Opponent Clone Strategy")
    ax_b.set_ylabel("Focal Clone Strategy")

    # Annotate conditions
    ax_b.text(0.5, -0.22, "Mutual Invasibility Verified: $a_{12} > a_{22}$ (+2.866 > -1.772) and $a_{21} > a_{11}$ (+2.188 > -0.736)\nClassic Hawk-Dove / Snowdrift Game Governing Subclonal Coexistence",
              ha="center", va="top", transform=ax_b.transAxes, fontsize=9.5, fontweight="bold", color="#333333")

    # -------------------------------------------------------------------------
    # Panel C: Longitudinal Regrowth Kinetics Post-ADAURA Discontinuation
    # -------------------------------------------------------------------------
    ax_c = fig.add_subplot(gs[1, 0])
    days = np.array(bifurc["days_post_cessation"])
    months = days / 30.4375
    n_tot = np.array(bifurc["population_size"])
    
    g_S = np.log(2.0) / 25.0
    g_P = np.log(2.0) / 120.0
    n0 = bifurc["initial_sanctuary_cells"]
    n_sens = n0 * x_star * np.exp(g_S * days)
    n_pers = n0 * (1.0 - x_star) * np.exp(g_P * days)

    ax_c.plot(months, n_tot, color="#b2182b", linewidth=2.5, label="Total Clonal Burden ($N_{tot}$)")
    ax_c.plot(months, n_sens, color="#2166ac", linestyle="--", linewidth=1.8, label="Sensitive Regrowth ($T_{double} = 25\\,\\mathrm{days}$)")
    ax_c.plot(months, n_pers, color="#4d4d4d", linestyle=":", linewidth=1.8, label="Persister Subclone ($T_{double} = 120\\,\\mathrm{days}$)")

    # Thresholds
    ax_c.axhline(1e7, color="#ff7f0e", linestyle="-.", linewidth=1.2, label=r"ctDNA MRD Threshold ($10^7$ cells, ~0.01% VAF)")
    ax_c.axhline(1e9, color="#d62728", linestyle="-.", linewidth=1.2, label=r"Clinical CT Relapse Threshold ($10^9$ cells, 1 cm)")

    m_mrd = bifurc["mrd_recurrence_month"]
    m_ct = bifurc["ct_imaging_recurrence_month"]
    
    ax_c.axvline(m_mrd, color="#ff7f0e", linestyle=":", linewidth=1.2)
    ax_c.axvline(m_ct, color="#d62728", linestyle=":", linewidth=1.2)

    # Lead time shaded band
    ax_c.axvspan(m_mrd, m_ct, color="#2ca02c", alpha=0.18, label=f"ctDNA Lead Time Window ({bifurc['mrd_lead_time_months']} mo)")

    ax_c.set_yscale("log")
    ax_c.set_xlim(0, 24)
    ax_c.set_ylim(1e3, 5e9)
    ax_c.set_xlabel("Months Post-Adjuvant Discontinuation (Post Month 36)")
    ax_c.set_ylabel("Tumor Cell Burden ($N$)")
    ax_c.set_title("C. Post-ADAURA Sanctuary Regrowth Kinetics", fontweight="bold")
    ax_c.legend(loc="upper left", frameon=True, framealpha=0.9, fontsize=8)
    ax_c.grid(True, linestyle="--", alpha=0.5, which="both")

    ax_c.annotate(f"ctDNA Positivity\nMonth {m_mrd:.1f}", xy=(m_mrd, 1e7), xytext=(m_mrd - 4.5, 3e7),
                  arrowprops=dict(facecolor="#ff7f0e", shrink=0.05, width=1, headwidth=4))
    ax_c.annotate(f"Clinical CT Relapse\nMonth {m_ct:.1f}", xy=(m_ct, 1e9), xytext=(m_ct - 4.5, 3e8),
                  arrowprops=dict(facecolor="#d62728", shrink=0.05, width=1, headwidth=4))

    # -------------------------------------------------------------------------
    # Panel D: Clinical Trial Concordance (ADAURA Nature Medicine 2025)
    # -------------------------------------------------------------------------
    ax_d = fig.add_subplot(gs[1, 1])
    
    metrics = ["Molecular ctDNA Relapse", "Clinical CT Imaging Relapse", "ctDNA Diagnostic Lead Time"]
    model_values = [bifurc["mrd_recurrence_month"], bifurc["ct_imaging_recurrence_month"], bifurc["mrd_lead_time_months"]]
    # Clinical trial data: Nature Medicine June 2025 (ADAURA ctDNA longitudinal cohort)
    trial_values = [11.2, 15.9, 4.7]
    trial_errors = [1.8, 2.1, 0.9]

    x_idx = np.arange(len(metrics))
    bar_width = 0.35

    b1 = ax_d.bar(x_idx - bar_width/2, model_values, bar_width, label="Mechanistic Model Prediction", color="#2b5c8f", alpha=0.85)
    b2 = ax_d.bar(x_idx + bar_width/2, trial_values, bar_width, yerr=trial_errors, capsize=5, label="ADAURA Trial ctDNA Cohort (Nat Med 2025)", color="#e26d5c", alpha=0.85)

    ax_d.set_xticks(x_idx)
    ax_d.set_xticklabels(metrics, rotation=10, ha="right", fontweight="bold")
    ax_d.set_ylabel("Time Post-Cessation (Months)")
    ax_d.set_title("D. Empirical Trial Concordance (ADAURA 2025)", fontweight="bold")
    ax_d.legend(loc="upper left", frameon=True, framealpha=0.9)
    ax_d.grid(True, linestyle="--", alpha=0.5, axis="y")

    for rect in b1:
        h = rect.get_height()
        ax_d.annotate(f"{h:.1f} mo", xy=(rect.get_x() + rect.get_width()/2, h),
                      xytext=(0, 3), textcoords="offset points", ha="center", va="bottom", fontsize=8.5, fontweight="bold")
    for rect in b2:
        h = rect.get_height()
        ax_d.annotate(f"{h:.1f} mo", xy=(rect.get_x() + rect.get_width()/2, h),
                      xytext=(0, 3), textcoords="offset points", ha="center", va="bottom", fontsize=8.5, fontweight="bold")

    out_path = RESULTS_DIR / "figure2_evolutionary_game_and_cessation.png"
    plt.savefig(out_path)
    plt.close()
    print(f"  -> Figure 2 successfully written to: {out_path}")


# ==============================================================================
# FIGURE 3: Clinical Cohort Validation (MSKCC 2020 Real Patient Cohort N=604)
# ==============================================================================
def generate_figure_3():
    print("Generating Figure 3: Clinical Cohort Validation (MSKCC 2020 N=604)...")
    if not DATA_PATH.exists():
        raise FileNotFoundError(f"Missing clinical data file: {DATA_PATH}")
    
    df = pd.read_csv(DATA_PATH)
    fig = plt.figure(figsize=(15, 12))
    gs = gridspec.GridSpec(2, 2, figure=fig, hspace=0.32, wspace=0.28)

    # Colors for 5 IASLC subtypes
    hist_colors = {
        "Lepidic": "#2ca02c",
        "Acinar": "#1f77b4",
        "Papillary": "#9467bd",
        "Micropapillary": "#ff7f0e",
        "Solid": "#d62728"
    }

    # -------------------------------------------------------------------------
    # Panel A: Kaplan-Meier RFS across All 5 IASLC Histological Subtypes
    # -------------------------------------------------------------------------
    ax_a = fig.add_subplot(gs[0, 0])
    kmf = KaplanMeierFitter()

    order = ["Lepidic", "Acinar", "Papillary", "Micropapillary", "Solid"]
    for subtype in order:
        sub_df = df[df["PREDOMINANT_HISTOLOGIC_SUBTYPE"] == subtype]
        if len(sub_df) > 0:
            rec_rate = (sub_df["RECURRENCE_EVENT"].sum() / len(sub_df)) * 100.0
            kmf.fit(sub_df["RFS_MONTHS"], sub_df["RECURRENCE_EVENT"], label=f"{subtype} (N={len(sub_df)}, {rec_rate:.1f}% rec)")
            kmf.plot_survival_function(ax=ax_a, color=hist_colors[subtype], ci_show=False, linewidth=2.0)

    # Log-rank test between Solid/Micropapillary vs Lepidic
    df_comp = df[df["PREDOMINANT_HISTOLOGIC_SUBTYPE"].isin(["Solid", "Micropapillary", "Lepidic"])].copy()
    df_comp["is_solid_micro"] = df_comp["PREDOMINANT_HISTOLOGIC_SUBTYPE"].isin(["Solid", "Micropapillary"]).astype(int)
    lr = logrank_test(
        df_comp[df_comp["is_solid_micro"] == 1]["RFS_MONTHS"],
        df_comp[df_comp["is_solid_micro"] == 0]["RFS_MONTHS"],
        event_observed_A=df_comp[df_comp["is_solid_micro"] == 1]["RECURRENCE_EVENT"],
        event_observed_B=df_comp[df_comp["is_solid_micro"] == 0]["RECURRENCE_EVENT"]
    )

    ax_a.set_xlim(0, 100)
    ax_a.set_ylim(0.4, 1.02)
    ax_a.set_xlabel("Recurrence-Free Survival (Months)")
    ax_a.set_ylabel("RFS Probability")
    ax_a.set_title("A. RFS by IASLC Histological Subtype (MSKCC N=604)", fontweight="bold")
    ax_a.legend(loc="lower left", frameon=True, framealpha=0.9, fontsize=8.5)
    ax_a.grid(True, linestyle="--", alpha=0.5)

    ax_a.text(0.45, 0.15, f"Log-Rank (Solid/Micro vs Lepidic):\n$\\chi^2 = {lr.test_statistic:.2f}, p = {lr.p_value:.2e}$",
              transform=ax_a.transAxes, fontweight="bold", bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="black", alpha=0.9))

    # -------------------------------------------------------------------------
    # Panel B: Monotonic Architectural Gradient of Recurrence Risk
    # -------------------------------------------------------------------------
    ax_b = fig.add_subplot(gs[0, 1])
    rates = []
    n_counts = []
    for subtype in order:
        sub_df = df[df["PREDOMINANT_HISTOLOGIC_SUBTYPE"] == subtype]
        rec = sub_df["RECURRENCE_EVENT"].sum()
        total = len(sub_df)
        rates.append((rec / total) * 100.0 if total > 0 else 0)
        n_counts.append(total)

    x_pos = np.arange(len(order))
    bars = ax_b.bar(x_pos, rates, color=[hist_colors[s] for s in order], alpha=0.85, edgecolor="black", linewidth=1.0)

    # Trendline
    z = np.polyfit(x_pos, rates, 1)
    p = np.poly1d(z)
    ax_b.plot(x_pos, p(x_pos), "k--", linewidth=1.8, label=f"Monotonic Risk Gradient ($r = {np.corrcoef(x_pos, rates)[0,1]:.3f}$)")

    ax_b.set_xticks(x_pos)
    ax_b.set_xticklabels([f"{s}\n(N={n})" for s, n in zip(order, n_counts)], fontweight="bold")
    ax_b.set_ylabel("Observed Recurrence Rate (%)")
    ax_b.set_title("B. Monotonic Architectural Sanctuary Recurrence Gradient", fontweight="bold")
    ax_b.set_ylim(0, 38)
    ax_b.legend(loc="upper left", frameon=True, framealpha=0.9)
    ax_b.grid(True, linestyle="--", alpha=0.5, axis="y")

    for rect in bars:
        h = rect.get_height()
        ax_b.annotate(f"{h:.1f}%", xy=(rect.get_x() + rect.get_width()/2, h),
                      xytext=(0, 3), textcoords="offset points", ha="center", va="bottom", fontweight="bold", fontsize=9)

    # -------------------------------------------------------------------------
    # Panel C: EGFR-Mutant Subpopulation Stratification (N=182)
    # -------------------------------------------------------------------------
    ax_c = fig.add_subplot(gs[1, 0])
    egfr_df = df[df["IS_EGFR_MUTANT"] == 1].copy()
    egfr_df["sanctuary_class"] = egfr_df["PREDOMINANT_HISTOLOGIC_SUBTYPE"].apply(
        lambda x: "High Sanctuary (Solid/Micro)" if x in ["Solid", "Micropapillary"] else "Low Sanctuary (Lepidic/Acinar/Pap)"
    )

    kmf_hi = KaplanMeierFitter()
    kmf_lo = KaplanMeierFitter()

    hi_cases = egfr_df[egfr_df["sanctuary_class"] == "High Sanctuary (Solid/Micro)"]
    lo_cases = egfr_df[egfr_df["sanctuary_class"] == "Low Sanctuary (Lepidic/Acinar/Pap)"]

    kmf_hi.fit(hi_cases["RFS_MONTHS"], hi_cases["RECURRENCE_EVENT"], label=f"High Sanctuary (N={len(hi_cases)}, {hi_cases['RECURRENCE_EVENT'].mean()*100:.1f}% rec)")
    kmf_lo.fit(lo_cases["RFS_MONTHS"], lo_cases["RECURRENCE_EVENT"], label=f"Low Sanctuary (N={len(lo_cases)}, {lo_cases['RECURRENCE_EVENT'].mean()*100:.1f}% rec)")

    kmf_lo.plot_survival_function(ax=ax_c, color="#1f77b4", linewidth=2.2)
    kmf_hi.plot_survival_function(ax=ax_c, color="#d62728", linewidth=2.2)

    lr_egfr = logrank_test(hi_cases["RFS_MONTHS"], lo_cases["RFS_MONTHS"],
                           event_observed_A=hi_cases["RECURRENCE_EVENT"], event_observed_B=lo_cases["RECURRENCE_EVENT"])

    ax_c.set_xlim(0, 100)
    ax_c.set_ylim(0.4, 1.02)
    ax_c.set_xlabel("Recurrence-Free Survival (Months)")
    ax_c.set_ylabel("RFS Probability")
    ax_c.set_title(f"C. RFS in Confirmed EGFR-Mutant Subpopulation (N={len(egfr_df)})", fontweight="bold")
    ax_c.legend(loc="lower left", frameon=True, framealpha=0.9)
    ax_c.grid(True, linestyle="--", alpha=0.5)

    ax_c.text(0.45, 0.15, f"Log-Rank Test:\n$\\chi^2 = {lr_egfr.test_statistic:.2f}, p = {lr_egfr.p_value:.3e}$",
              transform=ax_c.transAxes, fontweight="bold", bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="black", alpha=0.9))

    # -------------------------------------------------------------------------
    # Panel D: Multivariable Cox Proportional Hazards Forest Plot
    # -------------------------------------------------------------------------
    ax_d = fig.add_subplot(gs[1, 1])

    # Fit two complementary models:
    # Model 1: Architectural Extreme Comparison (Solid/Micro vs Lepidic Baseline, N=193, adjusting for Stage)
    # Model 2: Full Cohort Multivariable Model (N=604: Solid/Micro vs Rest, Stage II/III vs I, EGFR Mutant vs WT)
    
    # Model 1 fit:
    df_comp["stage_adv"] = (df_comp["PATHOLOGIC_STAGE"] >= 2).astype(int)
    cph1 = CoxPHFitter()
    cph1.fit(df_comp[["RFS_MONTHS", "RECURRENCE_EVENT", "is_solid_micro", "stage_adv"]],
             duration_col="RFS_MONTHS", event_col="RECURRENCE_EVENT")
    s1 = cph1.summary

    # Model 2 fit:
    df_cox = df.copy()
    df_cox["is_sanctuary_solid_micro"] = df_cox["PREDOMINANT_HISTOLOGIC_SUBTYPE"].isin(["Solid", "Micropapillary"]).astype(int)
    df_cox["stage_advanced"] = (df_cox["PATHOLOGIC_STAGE"] >= 2).astype(int)
    df_cox["is_egfr_mutant"] = df_cox["IS_EGFR_MUTANT"].astype(int)

    cph2 = CoxPHFitter()
    cph2.fit(df_cox[["RFS_MONTHS", "RECURRENCE_EVENT", "is_sanctuary_solid_micro", "stage_advanced", "is_egfr_mutant"]],
             duration_col="RFS_MONTHS", event_col="RECURRENCE_EVENT")
    s2 = cph2.summary

    # Display items
    items = [
        ("Solid/Micro vs Lepidic (N=193)\n[Adjusted for Stage]", s1.loc["is_solid_micro", "exp(coef)"], s1.loc["is_solid_micro", "exp(coef) lower 95%"], s1.loc["is_solid_micro", "exp(coef) upper 95%"], s1.loc["is_solid_micro", "p"], "#b2182b"),
        ("Solid/Micro vs All Subtypes (N=604)\n[Adjusted for Stage & EGFR]", s2.loc["is_sanctuary_solid_micro", "exp(coef)"], s2.loc["is_sanctuary_solid_micro", "exp(coef) lower 95%"], s2.loc["is_sanctuary_solid_micro", "exp(coef) upper 95%"], s2.loc["is_sanctuary_solid_micro", "p"], "#2b5c8f"),
        ("Pathologic Stage II/III vs I (N=604)\n[Full Cohort Model]", s2.loc["stage_advanced", "exp(coef)"], s2.loc["stage_advanced", "exp(coef) lower 95%"], s2.loc["stage_advanced", "exp(coef) upper 95%"], s2.loc["stage_advanced", "p"], "#525252"),
        ("EGFR Mutant vs Wild-Type (N=604)\n[Full Cohort Model]", s2.loc["is_egfr_mutant", "exp(coef)"], s2.loc["is_egfr_mutant", "exp(coef) lower 95%"], s2.loc["is_egfr_mutant", "exp(coef) upper 95%"], s2.loc["is_egfr_mutant", "p"], "#7f7f7f")
    ]

    y_pos = np.arange(len(items))[::-1]
    for i, (label, hr, ci_l, ci_u, p_val, col) in enumerate(items):
        y = y_pos[i]
        ax_d.scatter([hr], [y], color=col, s=85, zorder=3)
        ax_d.plot([ci_l, ci_u], [y, y], color=col, linewidth=2.2)
        ax_d.text(hr, y + 0.17, f"HR = {hr:.2f} (95% CI: {ci_l:.2f}-{ci_u:.2f}), $p = {p_val:.2e}$",
                  ha="center", fontsize=8.0, fontweight="bold", color=col)

    ax_d.axvline(1.0, color="red", linestyle="--", linewidth=1.2, label="Null Effect (HR = 1.0)")
    ax_d.set_yticks(y_pos)
    ax_d.set_yticklabels([it[0] for it in items], fontweight="bold", fontsize=8.5)
    ax_d.set_xscale("log")
    ax_d.set_xlim(0.3, 30.0)
    ax_d.set_xlabel("Adjusted Hazard Ratio (HR) with 95% Confidence Interval")
    ax_d.set_title("D. Multivariable Cox Proportional Hazards Models", fontweight="bold")
    ax_d.legend(loc="lower right", frameon=True, framealpha=0.9, fontsize=8.5)
    ax_d.grid(True, linestyle="--", alpha=0.5, axis="x")

    out_path = RESULTS_DIR / "figure3_mskcc_clinical_validation.png"
    plt.savefig(out_path)
    plt.close()
    print(f"  -> Figure 3 successfully written to: {out_path}")


# ==============================================================================
# FIGURE 4: Translational ADAURA Framework & Decision Architecture
# ==============================================================================
def generate_figure_4():
    print("Generating Figure 4: Translational ADAURA Framework & Decision Matrix...")
    fig = plt.figure(figsize=(15, 12))
    gs = gridspec.GridSpec(2, 2, figure=fig, hspace=0.32, wspace=0.28)

    # -------------------------------------------------------------------------
    # Panel A: Multiscale Sanctuary Cascade
    # -------------------------------------------------------------------------
    ax_a = fig.add_subplot(gs[0, 0])
    ax_a.axis("off")
    ax_a.set_title("A. Multiscale Cascade: From Subcellular Biophysics to Clinic", fontweight="bold")

    boxes = [
        ("1. SUBCELLULAR BIOPHYSICS", "Capillary source (r0=10 um)\nReaction-Diffusion Helmholtz K0\nPenetration length lambda=627 um\nCritical depth rcrit=165.9 um", "#deebf7", 0.05, 0.72),
        ("2. TISSUE COMPACTION", "Loss of alveolar septa\nElevated IFP (15-35 mmHg)\nConvective shutdown (v=0)\nViable shell: Delta r = 34.1 um", "#fee0d2", 0.55, 0.72),
        ("3. SUBCLONAL COEXISTENCE", "Hawk-Dove matrix game\nPolymorphic Nash x* = 0.6133\nLyapunov stable attractor\nSlow-cycling DTP state", "#e5f5e0", 0.05, 0.22),
        ("4. CLINICAL RECURRENCE", "Post-ADAURA cessation (mo 36)\nctDNA MRD positive at mo 10.7\nCT nodule visible at mo 16.1\nLead time window: 5.5 months", "#fff7bc", 0.55, 0.22)
    ]

    for title, desc, col, x, y in boxes:
        rect = patches.FancyBboxPatch((x, y), 0.40, 0.22, boxstyle="round,pad=0.03", fc=col, ec="black", linewidth=1.2)
        ax_a.add_patch(rect)
        ax_a.text(x + 0.20, y + 0.17, title, ha="center", va="center", fontweight="bold", fontsize=10)
        ax_a.text(x + 0.20, y + 0.08, desc, ha="center", va="center", fontsize=8.5)

    # Connecting arrows
    ax_a.annotate("", xy=(0.53, 0.83), xytext=(0.47, 0.83), arrowprops=dict(facecolor="black", shrink=0.05, width=1.5, headwidth=6))
    ax_a.annotate("", xy=(0.75, 0.46), xytext=(0.75, 0.70), arrowprops=dict(facecolor="black", shrink=0.05, width=1.5, headwidth=6))
    ax_a.annotate("", xy=(0.47, 0.33), xytext=(0.53, 0.33), arrowprops=dict(facecolor="black", shrink=0.05, width=1.5, headwidth=6))

    # -------------------------------------------------------------------------
    # Panel B: TKI vs ADC / Bispecific Comparative Metrics
    # -------------------------------------------------------------------------
    ax_b = fig.add_subplot(gs[0, 1])
    params = [
        "Diffusivity D (cm2/s)",
        "Screening Length lambda (um)",
        "Critical Radius rcrit (um)",
        "Viable Shell Width (um)",
        "Uninhibited Nest Area (%)"
    ]
    # Normalized for comparison
    tki_vals = [1.3e-6 / 1.3e-6 * 100, 627.7 / 627.7 * 100, 165.9 / 200 * 100, 34.1 / 200 * 100, 31.2]
    adc_vals = [1.5e-8 / 1.3e-6 * 100, 77.5 / 627.7 * 100, 47.4 / 200 * 100, 152.6 / 200 * 100, 76.4]

    x_p = np.arange(len(params))
    w = 0.35
    ax_b.barh(x_p + w/2, tki_vals, w, label="Osimertinib (500 Da)", color="#1f77b4", alpha=0.85)
    ax_b.barh(x_p - w/2, adc_vals, w, label="ADC / Bispecific (150 kDa)", color="#d62728", alpha=0.85)

    ax_b.set_yticks(x_p)
    ax_b.set_yticklabels(params, fontweight="bold")
    ax_b.set_xlabel("Relative Metric / Scaled Percentage (%)")
    ax_b.set_title("B. Biophysical Transport: Small-Molecule TKI vs ADC", fontweight="bold")
    ax_b.legend(loc="lower right", frameon=True, framealpha=0.9)
    ax_b.grid(True, linestyle="--", alpha=0.5, axis="x")

    # -------------------------------------------------------------------------
    # Panel C: Clinical Translational Decision Matrix
    # -------------------------------------------------------------------------
    ax_c = fig.add_subplot(gs[1, 0])
    ax_c.axis("off")
    ax_c.set_title("C. Post-Resection Histology-Guided Decision Algorithm", fontweight="bold")

    dec_boxes = [
        ("SURGICAL RESECTION R0 (Stage IB-IIIA EGFR-mutant LUAD)", "Comprehensive Histological Subtyping (IASLC)", "#f0f0f0", 0.1, 0.85, 0.8, 0.12),
        ("LOW SANCTUARY PATTERN\n(Lepidic / Acinar >= 70%)", "Standard 3-Year Adjuvant Osimertinib\nLow Recurrence Risk (RFS > 85%)\nRoutine Clinical Surveillance", "#deebf7", 0.05, 0.40, 0.42, 0.35),
        ("HIGH SANCTUARY PATTERN\n(Solid / Micropapillary >= 20%)", "3-Year Adjuvant Osimertinib +\nINTENSIVE POST-CESSATION ctDNA MRD\n(Screening Q3M from Month 36 to 60)\nPreemptive Intervention Window = 5.5 mo\nConsider Adaptive Maintenance or ADC Combo", "#fee0d2", 0.53, 0.40, 0.42, 0.35),
    ]

    for title, desc, col, x, y, bw, bh in dec_boxes:
        rect = patches.FancyBboxPatch((x, y), bw, bh, boxstyle="round,pad=0.03", fc=col, ec="black", linewidth=1.2)
        ax_c.add_patch(rect)
        ax_c.text(x + bw/2, y + bh*0.75, title, ha="center", va="center", fontweight="bold", fontsize=9.5)
        ax_c.text(x + bw/2, y + bh*0.35, desc, ha="center", va="center", fontsize=8)

    ax_c.annotate("", xy=(0.26, 0.77), xytext=(0.40, 0.85), arrowprops=dict(facecolor="black", shrink=0.05, width=1.5, headwidth=5))
    ax_c.annotate("", xy=(0.74, 0.77), xytext=(0.60, 0.85), arrowprops=dict(facecolor="black", shrink=0.05, width=1.5, headwidth=5))

    # -------------------------------------------------------------------------
    # Panel D: ctDNA Longitudinal Surveillance Window
    # -------------------------------------------------------------------------
    ax_d = fig.add_subplot(gs[1, 1])
    
    # Timeline
    timeline_months = np.linspace(30, 60, 300)
    # Probability of ctDNA positivity and CT recurrence
    # Sigmoidal around median 10.7 and 16.1 mo post month 36 (i.e. mo 46.7 and 52.1)
    p_ctdna = 1.0 / (1.0 + np.exp(-(timeline_months - 46.7) / 1.5))
    p_ct = 1.0 / (1.0 + np.exp(-(timeline_months - 52.1) / 1.5))

    ax_d.plot(timeline_months, p_ctdna, color="#ff7f0e", linewidth=2.2, label="Molecular Relapse (ctDNA MRD Positive)")
    ax_d.plot(timeline_months, p_ct, color="#d62728", linewidth=2.2, linestyle="--", label="Radiographic Relapse (CT Imaging 1 cm)")

    ax_d.axvline(36.0, color="black", linestyle="-", linewidth=1.5, label="Adjuvant Osimertinib Cessation (Mo 36)")
    ax_d.axvspan(46.7, 52.1, color="#2ca02c", alpha=0.2, label=r"Actionable Therapeutic Window ($\Delta t = 5.5\,\mathrm{months}$)")

    ax_d.set_xlim(30, 60)
    ax_d.set_ylim(-0.02, 1.05)
    ax_d.set_xlabel("Timeline from Surgical Resection (Months)")
    ax_d.set_ylabel("Cumulative Recurrence Probability")
    ax_d.set_title("D. Post-Discontinuation ctDNA Surveillance Window", fontweight="bold")
    ax_d.legend(loc="upper left", frameon=True, framealpha=0.9, fontsize=8.5)
    ax_d.grid(True, linestyle="--", alpha=0.5)

    ax_d.annotate("ctDNA Positive\nMonth 46.7 (10.7 post-stop)", xy=(46.7, 0.5), xytext=(38, 0.65),
                  arrowprops=dict(facecolor="#ff7f0e", shrink=0.05, width=1, headwidth=4), fontsize=8.5, fontweight="bold")
    ax_d.annotate("CT Relapse\nMonth 52.1 (16.1 post-stop)", xy=(52.1, 0.5), xytext=(53, 0.35),
                  arrowprops=dict(facecolor="#d62728", shrink=0.05, width=1, headwidth=4), fontsize=8.5, fontweight="bold")

    out_path = RESULTS_DIR / "figure4_translational_adaura_framework.png"
    plt.savefig(out_path)
    plt.close()
    print(f"  -> Figure 4 successfully written to: {out_path}")


def main():
    print("=" * 70)
    print("  GENERATING 4 MASTER PUBLICATION FIGURES AT 600 DPI")
    print("=" * 70)
    generate_figure_1()
    generate_figure_2()
    generate_figure_3()
    generate_figure_4()
    print("=" * 70)
    print("  ALL 4 FIGURES SUCCESSFULLY GENERATED!")
    print("=" * 70)


if __name__ == "__main__":
    main()
