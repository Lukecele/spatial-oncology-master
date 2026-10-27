#!/usr/bin/env python3
"""
Comprehensive Scientific Audit & Re-Calculation Script.
Independently verifies every biophysical, game-theoretic, and clinical survival calculation
to rigorously evaluate all points raised by Claude.
"""

import sys
from pathlib import Path
import numpy as np
import scipy.special as sp
from scipy.optimize import brentq
import pandas as pd
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import logrank_test

print("=" * 80)
print("  INDEPENDENT AUDIT & RE-CALCULATION OF THE COMPLETE RESEARCH STUDY")
print("=" * 80)

# ==============================================================================
# PILLAR 1: BESSEL TRANSPORT MODEL & SENSITIVITY PHASE SPACE
# ==============================================================================
print("\n" + "-" * 80)
print("PILLAR 1: EXACT BESSEL REACTION-DIFFUSION TRANSPORT & SENSITIVITY")
print("-" * 80)

def compute_rcrit(D, k, r0=10.0, theta=0.35):
    lam = np.sqrt(D / k) * 10000.0  # um
    def f(r):
        return sp.k0(r / lam) / sp.k0(r0 / lam) - theta
    try:
        return brentq(f, r0 + 1e-4, 5000.0)
    except Exception:
        return np.nan

D0 = 1.3e-6   # cm^2/s
k0 = 3.3e-4   # s^-1
r0 = 10.0     # um
theta0 = 0.352 # 8.8 nM / 25 nM = 0.352

rc_base = compute_rcrit(D0, k0, r0, 0.35)
rc_base_exact = compute_rcrit(D0, k0, r0, theta0)
lam_base = np.sqrt(D0 / k0) * 10000.0

print(f"Base Parameters:")
print(f"  D0 = {D0:.2e} cm^2/s, k0 = {k0:.2e} s^-1, r0 = {r0:.1f} um")
print(f"  Characteristic screening length lambda = {lam_base:.2f} um")
print(f"  r_crit (theta = 0.350) : {rc_base:.2f} um")
print(f"  r_crit (theta = 0.352) : {rc_base_exact:.2f} um")
print(f"  Sanctuary shell (r_ox = 200 um, theta = 0.350) : Delta r = {200.0 - rc_base:.2f} um")
print(f"  Sanctuary shell (r_ox = 200 um, theta = 0.352) : Delta r = {200.0 - rc_base_exact:.2f} um")

# Claude's hand calculation check:
rc_D_plus30 = compute_rcrit(D0 * 1.3, k0, r0, 0.35)
rc_k_minus30 = compute_rcrit(D0, k0 * 0.7, r0, 0.35)
rc_both = compute_rcrit(D0 * 1.3, k0 * 0.7, r0, 0.35)
rc_both_tight = compute_rcrit(D0 * 0.7, k0 * 1.3, r0, 0.35)

print(f"\nVerification of Claude's Hand Calculations (theta = 0.35):")
print(f"  D +30%           -> r_crit = {rc_D_plus30:.2f} um  (Claude said ~181 um, Diff: {rc_D_plus30 - 181:.2f} um)")
print(f"  k -30%           -> r_crit = {rc_k_minus30:.2f} um  (Claude said ~187 um, Diff: {rc_k_minus30 - 187:.2f} um)")
print(f"  D +30% & k -30%  -> r_crit = {rc_both:.2f} um  (Claude said ~204 um, Diff: {rc_both - 204:.2f} um)")
print(f"  D -30% & k +30%  -> r_crit = {rc_both_tight:.2f} um")

print(f"\nShell status under Claude's worst-case parameter combination:")
print(f"  If D +30% and k -30%: r_crit ({rc_both:.1f} um) > r_ox (200 um) -> Shell collapses (Delta r = {max(0, 200 - rc_both):.1f} um).")
print(f"  If r_ox = 150 um (Thomlinson & Gray lower limit): r_crit ({rc_base:.1f} um) > 150 um -> Baseline shell collapses.")

# Multidimensional Phase Map
print(f"\nMultidimensional Parameter Sensitivity Grid (D in [0.5, 2.5]e-6, k in [1.5, 6.0]e-4):")
D_grid = np.linspace(0.5e-6, 2.5e-6, 50)
k_grid = np.linspace(1.5e-4, 6.0e-4, 50)
total_pts = len(D_grid) * len(k_grid)

for rox in [150.0, 175.0, 200.0]:
    for th in [0.25, 0.35, 0.45]:
        viable_count = 0
        widths = []
        for D_val in D_grid:
            for k_val in k_grid:
                rc = compute_rcrit(D_val, k_val, r0, th)
                if rc < rox:
                    viable_count += 1
                    widths.append(rox - rc)
        pct = (viable_count / total_pts) * 100.0
        mean_w = np.mean(widths) if widths else 0.0
        print(f"  r_ox = {rox:5.1f} um | theta = {th:.2f} -> Shell exists in {pct:5.1f}% of grid (mean width: {mean_w:4.1f} um)")

# ADC macromolecule comparison
D_adc = 1.5e-8
k_adc = 2.5e-4
lam_adc = np.sqrt(D_adc / k_adc) * 10000.0
rc_adc = compute_rcrit(D_adc, k_adc, r0, 0.35)
print(f"\nADC Macromolecule:")
print(f"  D_adc = {D_adc:.2e} cm^2/s, k_adc = {k_adc:.2e} s^-1, lambda = {lam_adc:.2f} um")
print(f"  r_crit ADC = {rc_adc:.2f} um -> Viable sanctuary shell (r_ox=200): Delta r = {200 - rc_adc:.2f} um ({((200**2 - rc_adc**2)/(200**2))*100:.1f}% of viable cord)")


# ==============================================================================
# PILLAR 2: REPLICATOR DYNAMICS & POST-CESSATION REGROWTH
# ==============================================================================
print("\n" + "-" * 80)
print("PILLAR 2: REPLICATOR DYNAMICS & POST-CESSATION REGROWTH KINETICS")
print("-" * 80)

A = np.array([
    [-0.736,  2.866],
    [ 2.188, -1.772]
], dtype=np.float64)

a11, a12 = A[0, 0], A[0, 1]
a21, a22 = A[1, 0], A[1, 1]

alpha = a11 - a12 - a21 + a22
beta = a12 - a22
x_star = (a22 - a12) / alpha
lam_eig = x_star * (1.0 - x_star) * alpha

print(f"Payoff Matrix A:")
print(f"  [ [{a11:+.3f}, {a12:+.3f}], [{a21:+.3f}, {a22:+.3f}] ]")
print(f"Coexistence condition (Hawk-Dove / Snowdrift):")
print(f"  a12 > a22: {a12} > {a22} -> {a12 > a22}")
print(f"  a21 > a11: {a21} > {a11} -> {a21 > a11}")
print(f"Internal Nash Equilibrium x* = {x_star:.6f} ({x_star*100:.2f}% sensitive)")
print(f"Denominator alpha = {alpha:.4f}, beta = {beta:.4f}")
print(f"Jacobian eigenvalue lambda_eig = {lam_eig:.6f} (< 0 -> Asymptotically stable attractor)")

# Regrowth simulation:
g_S = np.log(2.0) / 25.0   # net doubling 25 days
g_P = np.log(2.0) / 120.0  # net doubling 120 days

# Test Claude's claim: Lead time Delta t is pure doubling arithmetic
doublings = np.log2(1e9 / 1e7) # log2(100) = 6.643856
exact_doubling_time_days = doublings * 25.0
exact_doubling_time_months = exact_doubling_time_days / 30.4375

print(f"\nClaude Lead Time Critique Analysis:")
print(f"  Ratio N_CT / N_MRD = 10^9 / 10^7 = 100")
print(f"  log2(100) = {doublings:.6f} doublings")
print(f"  doublings * 25 days = {exact_doubling_time_days:.2f} days = {exact_doubling_time_months:.2f} months")

# Simulate across varying N0
days_sim = np.arange(0, 3000)
for N0 in [200, 1000, 2000, 10000, 50000]:
    N_t = N0 * (x_star * np.exp(g_S * days_sim) + (1.0 - x_star) * np.exp(g_P * days_sim))
    d_mrd = days_sim[np.where(N_t >= 1e7)[0][0]]
    d_ct = days_sim[np.where(N_t >= 1e9)[0][0]]
    dt_days = d_ct - d_mrd
    m_mrd = d_mrd / 30.4375
    m_ct = d_ct / 30.4375
    dt_m = dt_days / 30.4375
    sens_pct_at_mrd = (N0 * x_star * np.exp(g_S * d_mrd)) / 1e7 * 100.0
    print(f"  N0 = {N0:5d} -> MRD: Month {m_mrd:4.1f} | CT: Month {m_ct:4.1f} | Lead Time: {dt_m:4.2f} months | Sens clone: {sens_pct_at_mrd:.2f}%")


# ==============================================================================
# PILLAR 3: MSKCC CLINICAL COHORT AUDIT (N = 604)
# ==============================================================================
print("\n" + "-" * 80)
print("PILLAR 3: MSKCC REAL PATIENT COHORT (N = 604)")
print("-" * 80)

csv_path = Path("/home/luca/Scrivania/PROGETTO_FINALE_ONCOLOGIA/06_DATI_CLINICI_MSKCC_604_PAZIENTI.csv")
df = pd.read_csv(csv_path)

print(f"Total Cohort Size : N = {len(df)}")
print(f"Total Recurrences : E = {int(df['RECURRENCE_EVENT'].sum())} ({df['RECURRENCE_EVENT'].mean()*100:.2f}%)")
print(f"EGFR Mutant Cases : N = {int((df['IS_EGFR_MUTANT'] == 1).sum())} ({(df['IS_EGFR_MUTANT'] == 1).mean()*100:.2f}%)")
print(f"Median Follow-Up  : {df['RFS_MONTHS'].median():.2f} months (IQR: {df['RFS_MONTHS'].quantile(0.25):.2f} - {df['RFS_MONTHS'].quantile(0.75):.2f})")
print(f"Patients at risk at 36 months: {(df['RFS_MONTHS'] >= 36).sum()}")
print(f"Patients at risk at 60 months: {(df['RFS_MONTHS'] >= 60).sum()}")

print("\nSubtype Distribution & Observed Events:")
for s in ["Lepidic", "Acinar", "Micropapillary", "Papillary", "Solid"]:
    sub = df[df["PREDOMINANT_HISTOLOGIC_SUBTYPE"] == s]
    n = len(sub)
    ev = int(sub["RECURRENCE_EVENT"].sum())
    print(f"  {s:15s}: N = {n:3d} ({n/len(df)*100:5.1f}%) | Events = {ev:2d} ({ev/n*100:5.2f}%)")

# Extreme contrast
df_ext = df[df["PREDOMINANT_HISTOLOGIC_SUBTYPE"].isin(["Solid", "Micropapillary", "Lepidic"])].copy()
df_ext["is_high"] = df_ext["PREDOMINANT_HISTOLOGIC_SUBTYPE"].isin(["Solid", "Micropapillary"]).astype(int)
df_ext["stage_adv"] = (df_ext["PATHOLOGIC_STAGE"] >= 2).astype(int)

c1 = CoxPHFitter().fit(df_ext[["RFS_MONTHS", "RECURRENCE_EVENT", "is_high", "stage_adv"]],
                       duration_col="RFS_MONTHS", event_col="RECURRENCE_EVENT")
s1 = c1.summary.loc["is_high"]
print(f"\nModel 1: Extreme Contrast (Solid/Micro vs Lepidic, N = {len(df_ext)}, Events = {int(df_ext['RECURRENCE_EVENT'].sum())}):")
print(f"  Adjusted HR = {s1['exp(coef)']:.3f} (95% CI: {s1['exp(coef) lower 95%']:.3f} - {s1['exp(coef) upper 95%']:.3f}), p = {s1['p']:.4e}")

# Full cohort
df_full = df.copy()
df_full["is_high"] = df_full["PREDOMINANT_HISTOLOGIC_SUBTYPE"].isin(["Solid", "Micropapillary"]).astype(int)
df_full["stage_adv"] = (df_full["PATHOLOGIC_STAGE"] >= 2).astype(int)

c2 = CoxPHFitter().fit(df_full[["RFS_MONTHS", "RECURRENCE_EVENT", "is_high", "stage_adv"]],
                       duration_col="RFS_MONTHS", event_col="RECURRENCE_EVENT")
s2 = c2.summary.loc["is_high"]
print(f"\nModel 2: Full Cohort (Solid/Micro vs Rest, N = {len(df_full)}, Events = {int(df_full['RECURRENCE_EVENT'].sum())}):")
print(f"  Adjusted HR = {s2['exp(coef)']:.3f} (95% CI: {s2['exp(coef) lower 95%']:.3f} - {s2['exp(coef) upper 95%']:.3f}), p = {s2['p']:.4f}")

# 5 Subtypes (Acinar Reference)
dummies = pd.get_dummies(df["PREDOMINANT_HISTOLOGIC_SUBTYPE"], prefix="sub", dtype=int)
df_sub = pd.concat([df[["RFS_MONTHS", "RECURRENCE_EVENT"]], df_full["stage_adv"],
                     dummies[["sub_Lepidic", "sub_Solid", "sub_Micropapillary", "sub_Papillary"]]], axis=1)
c3 = CoxPHFitter().fit(df_sub, duration_col="RFS_MONTHS", event_col="RECURRENCE_EVENT")
print(f"\nModel 3: 5-Subtype Multivariable Model (Acinar = Reference, Stage-Adjusted):")
for var in ["sub_Lepidic", "sub_Solid", "sub_Micropapillary", "sub_Papillary"]:
    s_row = c3.summary.loc[var]
    print(f"  {var:18s}: HR = {s_row['exp(coef)']:.3f} (95% CI: {s_row['exp(coef) lower 95%']:.3f} - {s_row['exp(coef) upper 95%']:.3f}), p = {s_row['p']:.4f}")

# EGFR mutant subset
df_egfr = df[df["IS_EGFR_MUTANT"] == 1].copy()
df_egfr["is_high"] = df_egfr["PREDOMINANT_HISTOLOGIC_SUBTYPE"].isin(["Solid", "Micropapillary"]).astype(int)
df_egfr["stage_adv"] = (df_egfr["PATHOLOGIC_STAGE"] >= 2).astype(int)

high_egfr = df_egfr[df_egfr["is_high"] == 1]
low_egfr = df_egfr[df_egfr["is_high"] == 0]
lr_egfr = logrank_test(high_egfr["RFS_MONTHS"], low_egfr["RFS_MONTHS"],
                       event_observed_A=high_egfr["RECURRENCE_EVENT"],
                       event_observed_B=low_egfr["RECURRENCE_EVENT"])
print(f"\nEGFR-Mutant Subset (N = {len(df_egfr)}, Events = {int(df_egfr['RECURRENCE_EVENT'].sum())}):")
print(f"  High sanctuary (Solid/Micro): N = {len(high_egfr)}, Events = {int(high_egfr['RECURRENCE_EVENT'].sum())} ({high_egfr['RECURRENCE_EVENT'].mean()*100:.1f}%)")
print(f"  Low sanctuary (Others):       N = {len(low_egfr)}, Events = {int(low_egfr['RECURRENCE_EVENT'].sum())} ({low_egfr['RECURRENCE_EVENT'].mean()*100:.1f}%)")
print(f"  Univariable Log-Rank p = {lr_egfr.p_value:.4f}")

print("\n" + "=" * 80)
print("  RE-CALCULATION & AUDIT EXECUTION COMPLETE")
print("=" * 80)
