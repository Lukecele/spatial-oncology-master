#!/usr/bin/env python3
"""
MSKCC 2020 Clinical Cohort Validation Engine (N = 604).
Analyzes architectural sanctuary patterns, IASLC histological subtypes,
and recurrence-free survival (RFS) in surgically resected lung adenocarcinoma.

Includes:
- Reverse Kaplan-Meier method for median follow-up estimation (Schemper & Smith 1996)
- Pairwise subtype Cox contrasts with Holm multiplicity adjustment
- Extreme biophysical boundary contrast (Solid/Micro vs Lepidic, N=193)
- Full cohort multivariable models (N=604)
- EGFR-mutant subgroup analysis demonstrating stage confounding
- Benchmark validation against the classical Freireich 6-MP clinical dataset
"""

from pathlib import Path
from typing import Dict, Any, Tuple
import numpy as np
import pandas as pd
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import logrank_test, multivariate_logrank_test


class MSKCCClinicalValidator:
    """Clinical validation engine for the MSKCC 2020 authenticated cohort."""

    def __init__(self, data_path: Path):
        self.data_path = Path(data_path)
        if not self.data_path.exists():
            raise FileNotFoundError(f"Cohort data file not found: {self.data_path}")
        self.df = pd.read_csv(self.data_path)
        self._validate_cohort_integrity()

    def _validate_cohort_integrity(self) -> None:
        """Ensures cohort matches authenticated clinical parameters."""
        expected_n = 604
        actual_n = len(self.df)
        if actual_n != expected_n:
            raise ValueError(f"Cohort size mismatch: expected {expected_n}, got {actual_n}")
        
        required_cols = [
            "patientId", "PREDOMINANT_HISTOLOGIC_SUBTYPE",
            "PATHOLOGIC_STAGE", "IS_EGFR_MUTANT",
            "RFS_MONTHS", "RECURRENCE_EVENT"
        ]
        missing = [c for c in required_cols if c not in self.df.columns]
        if missing:
            raise KeyError(f"Missing required clinical columns: {missing}")

    def get_cohort_summary(self) -> Dict[str, Any]:
        """Returns baseline characteristics of the MSKCC cohort."""
        egfr_cases = self.df[self.df["IS_EGFR_MUTANT"] == 1]
        recurrent_cases = self.df[self.df["RECURRENCE_EVENT"] == 1]
        
        subtype_counts = self.df["PREDOMINANT_HISTOLOGIC_SUBTYPE"].value_counts().to_dict()
        stage_counts = self.df["PATHOLOGIC_STAGE"].value_counts().to_dict()

        return {
            "total_patients": len(self.df),
            "egfr_mutant_patients": len(egfr_cases),
            "total_recurrences": len(recurrent_cases),
            "overall_recurrence_rate": float(len(recurrent_cases) / len(self.df)),
            "subtype_distribution": subtype_counts,
            "stage_distribution": stage_counts,
        }

    def compute_subtype_gradient(self) -> Dict[str, Dict[str, Any]]:
        """Computes observed recurrence rates and stage I distribution across the 5 IASLC subtypes."""
        order = ["Lepidic", "Acinar", "Papillary", "Micropapillary", "Solid"]
        results = {}
        for s in order:
            sub = self.df[self.df["PREDOMINANT_HISTOLOGIC_SUBTYPE"] == s]
            n = len(sub)
            rec = int(sub["RECURRENCE_EVENT"].sum()) if n > 0 else 0
            rate = float(rec / n) if n > 0 else 0.0
            st1 = int((sub["PATHOLOGIC_STAGE"] == 1).sum()) if n > 0 else 0
            st1_pct = float(st1 / n * 100.0) if n > 0 else 0.0
            results[s] = {"N": n, "recurrences": rec, "recurrence_rate": rate, "stage_1_pct": st1_pct}
        return results

    def compare_high_vs_low_sanctuary(self) -> Dict[str, Any]:
        """
        Log-rank comparison between high-sanctuary (Solid/Micropapillary)
        and baseline lepidic architecture (Extreme biophysical contrast, N=193).
        """
        sub_df = self.df[self.df["PREDOMINANT_HISTOLOGIC_SUBTYPE"].isin(["Solid", "Micropapillary", "Lepidic"])].copy()
        sub_df["is_high_sanctuary"] = sub_df["PREDOMINANT_HISTOLOGIC_SUBTYPE"].isin(["Solid", "Micropapillary"]).astype(int)

        high = sub_df[sub_df["is_high_sanctuary"] == 1]
        low = sub_df[sub_df["is_high_sanctuary"] == 0]

        lr = logrank_test(
            high["RFS_MONTHS"], low["RFS_MONTHS"],
            event_observed_A=high["RECURRENCE_EVENT"],
            event_observed_B=low["RECURRENCE_EVENT"]
        )

        return {
            "n_high_sanctuary": len(high),
            "n_low_sanctuary": len(low),
            "events_high": int(high["RECURRENCE_EVENT"].sum()),
            "events_low": int(low["RECURRENCE_EVENT"].sum()),
            "logrank_test_stat": float(lr.test_statistic),
            "logrank_p_value": float(lr.p_value)
        }

    def compute_followup_and_survival_metrics(self) -> Dict[str, Any]:
        """
        Computes follow-up duration distributions using reverse Kaplan-Meier method
        (Schemper & Smith 1996), auditing statistical power and tail reliability.
        """
        raw_median = float(self.df["RFS_MONTHS"].median())
        q25 = float(self.df["RFS_MONTHS"].quantile(0.25))
        q75 = float(self.df["RFS_MONTHS"].quantile(0.75))
        n_at_36 = int((self.df["RFS_MONTHS"] >= 36).sum())
        n_at_60 = int((self.df["RFS_MONTHS"] >= 60).sum())
        events_after_36 = int(((self.df["RFS_MONTHS"] > 36) & (self.df["RECURRENCE_EVENT"] == 1)).sum())

        # Reverse Kaplan-Meier follow-up estimator (censoring becomes event)
        kmf_rev = KaplanMeierFitter()
        kmf_rev.fit(self.df["RFS_MONTHS"], event_observed=(1 - self.df["RECURRENCE_EVENT"]))
        rev_km_median = float(kmf_rev.median_survival_time_)

        return {
            "raw_median_rfs_months": raw_median,
            "reverse_km_median_followup_months": rev_km_median,
            "iqr_rfs_followup": (q25, q75),
            "patients_at_risk_36_months": n_at_36,
            "patients_at_risk_60_months": n_at_60,
            "recurrence_events_after_36_months": events_after_36
        }

    def fit_multivariable_cox_models(self) -> Dict[str, Any]:
        """
        Fits multivariable Cox Proportional Hazards models:
        1. Pairwise contrasts vs Lepidic reference with Holm multiplicity adjustment
        2. Model 1: Extreme Boundary Comparison (N=193: Solid/Micro vs Lepidic, stage-adjusted)
        3. Model 2: Full Cohort (N=604: Solid/Micro vs Rest, stage-adjusted)
        4. Model 3: 5-Subtype Model with Acinar as reference (N=604, stage-adjusted)
        5. Subgroup: EGFR-mutant cohort analysis demonstrating stage confounding
        """
        df_full = self.df.copy()
        df_full["stage_adv"] = (df_full["PATHOLOGIC_STAGE"] >= 2).astype(int)
        df_full["is_solid_micro"] = df_full["PREDOMINANT_HISTOLOGIC_SUBTYPE"].isin(["Solid", "Micropapillary"]).astype(int)

        # Pairwise contrasts vs Lepidic reference
        subtypes = ["Solid", "Papillary", "Acinar", "Micropapillary"]
        pairwise_results = {}
        p_vals = []
        for s in subtypes:
            pair_df = df_full[df_full["PREDOMINANT_HISTOLOGIC_SUBTYPE"].isin([s, "Lepidic"])].copy()
            pair_df["is_target"] = (pair_df["PREDOMINANT_HISTOLOGIC_SUBTYPE"] == s).astype(int)
            cph = CoxPHFitter()
            cph.fit(pair_df[["RFS_MONTHS", "RECURRENCE_EVENT", "is_target", "stage_adv"]],
                    duration_col="RFS_MONTHS", event_col="RECURRENCE_EVENT")
            row = cph.summary.loc["is_target"]
            hr = float(row["exp(coef)"])
            ci_l = float(row["exp(coef) lower 95%"])
            ci_u = float(row["exp(coef) upper 95%"])
            p = float(row["p"])
            pairwise_results[s] = {"hr": hr, "ci_lower": ci_l, "ci_upper": ci_u, "p_value": p}
            p_vals.append((s, p))

        # Holm-Bonferroni correction over pairwise contrasts
        p_vals.sort(key=lambda x: x[1])
        m = len(p_vals)
        for rank, (s, p) in enumerate(p_vals):
            holm_p = min(1.0, p * (m - rank))
            pairwise_results[s]["holm_p_value"] = float(holm_p)

        # Model 1: Extreme Contrast (N=193)
        df_comp = df_full[df_full["PREDOMINANT_HISTOLOGIC_SUBTYPE"].isin(["Solid", "Micropapillary", "Lepidic"])].copy()
        cph1 = CoxPHFitter()
        cph1.fit(df_comp[["RFS_MONTHS", "RECURRENCE_EVENT", "is_solid_micro", "stage_adv"]],
                 duration_col="RFS_MONTHS", event_col="RECURRENCE_EVENT")
        s1 = cph1.summary.loc["is_solid_micro"]

        # Model 2: Full Cohort Stage-Adjusted (N=604)
        cph2a = CoxPHFitter()
        cph2a.fit(df_full[["RFS_MONTHS", "RECURRENCE_EVENT", "is_solid_micro", "stage_adv"]],
                  duration_col="RFS_MONTHS", event_col="RECURRENCE_EVENT")
        s2a = cph2a.summary.loc["is_solid_micro"]

        # Model 3: 5-Subtype with Acinar Reference (N=604)
        dummies = pd.get_dummies(df_full["PREDOMINANT_HISTOLOGIC_SUBTYPE"], prefix="sub", dtype=int)
        df_cox3 = pd.concat([df_full[["RFS_MONTHS", "RECURRENCE_EVENT", "stage_adv"]],
                             dummies[["sub_Lepidic", "sub_Solid", "sub_Micropapillary", "sub_Papillary"]]], axis=1)
        cph3 = CoxPHFitter()
        cph3.fit(df_cox3, duration_col="RFS_MONTHS", event_col="RECURRENCE_EVENT")
        s3 = cph3.summary

        # EGFR Subgroup
        df_egfr = df_full[df_full["IS_EGFR_MUTANT"] == 1].copy()
        df_egfr_ext = df_egfr[df_egfr["PREDOMINANT_HISTOLOGIC_SUBTYPE"].isin(["Solid", "Micropapillary", "Lepidic"])].copy()
        cph_egfr = CoxPHFitter()
        cph_egfr.fit(df_egfr_ext[["RFS_MONTHS", "RECURRENCE_EVENT", "is_solid_micro", "stage_adv"]],
                     duration_col="RFS_MONTHS", event_col="RECURRENCE_EVENT")
        s_egfr = cph_egfr.summary.loc["is_solid_micro"]

        model1_data = {
            "n": len(df_comp),
            "events": int(df_comp["RECURRENCE_EVENT"].sum()),
            "adjusted_hr": float(s1["exp(coef)"]),
            "ci_lower": float(s1["exp(coef) lower 95%"]),
            "ci_upper": float(s1["exp(coef) upper 95%"]),
            "p_value": float(s1["p"])
        }

        return {
            "pairwise_contrasts_vs_lepidic": pairwise_results,
            "model1_extreme": model1_data,
            "model1_extreme_contrast": model1_data,
            "model2_full_cohort": {
                "n": len(df_full),
                "events": int(df_full["RECURRENCE_EVENT"].sum()),
                "adjusted_hr": float(s2a["exp(coef)"]),
                "ci_lower": float(s2a["exp(coef) lower 95%"]),
                "ci_upper": float(s2a["exp(coef) upper 95%"]),
                "p_value": float(s2a["p"])
            },
            "model3_acinar_reference": {
                "lepidic": {"hr": float(s3.loc["sub_Lepidic", "exp(coef)"]), "ci_lower": float(s3.loc["sub_Lepidic", "exp(coef) lower 95%"]), "ci_upper": float(s3.loc["sub_Lepidic", "exp(coef) upper 95%"]), "p_value": float(s3.loc["sub_Lepidic", "p"])},
                "solid": {"hr": float(s3.loc["sub_Solid", "exp(coef)"]), "ci_lower": float(s3.loc["sub_Solid", "exp(coef) lower 95%"]), "ci_upper": float(s3.loc["sub_Solid", "exp(coef) upper 95%"]), "p_value": float(s3.loc["sub_Solid", "p"])},
                "papillary": {"hr": float(s3.loc["sub_Papillary", "exp(coef)"]), "ci_lower": float(s3.loc["sub_Papillary", "exp(coef) lower 95%"]), "ci_upper": float(s3.loc["sub_Papillary", "exp(coef) upper 95%"]), "p_value": float(s3.loc["sub_Papillary", "p"])},
                "micropapillary": {"hr": float(s3.loc["sub_Micropapillary", "exp(coef)"]), "ci_lower": float(s3.loc["sub_Micropapillary", "exp(coef) lower 95%"]), "ci_upper": float(s3.loc["sub_Micropapillary", "exp(coef) upper 95%"]), "p_value": float(s3.loc["sub_Micropapillary", "p"])},
                "stage_advanced": {"hr": float(s3.loc["stage_adv", "exp(coef)"]), "ci_lower": float(s3.loc["stage_adv", "exp(coef) lower 95%"]), "ci_upper": float(s3.loc["stage_adv", "exp(coef) upper 95%"]), "p_value": float(s3.loc["stage_adv", "p"])}
            },
            "egfr_subgroup_extreme": {
                "n": len(df_egfr_ext),
                "events": int(df_egfr_ext["RECURRENCE_EVENT"].sum()),
                "adjusted_hr": float(s_egfr["exp(coef)"]),
                "ci_lower": float(s_egfr["exp(coef) lower 95%"]),
                "ci_upper": float(s_egfr["exp(coef) upper 95%"]),
                "p_value": float(s_egfr["p"])
            }
        }


def run_freireich_benchmark_validation() -> Dict[str, float]:
    """
    Validates survival engine against the classic Freireich 6-MP leukemia benchmark
    (Gehan 1965; Freireich et al. Blood 1963).
    """
    # 6-MP vs Placebo (21 patients per arm)
    time_placebo = [1, 1, 2, 2, 3, 4, 4, 5, 5, 8, 8, 8, 8, 11, 11, 12, 12, 15, 17, 22, 23]
    event_placebo = [1] * 21
    time_6mp = [6, 6, 6, 6, 7, 9, 10, 10, 11, 13, 16, 17, 19, 20, 22, 23, 25, 32, 32, 34, 35]
    event_6mp = [1, 1, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0]

    df_bench = pd.DataFrame({
        "time": time_placebo + time_6mp,
        "event": event_placebo + event_6mp,
        "is_placebo": [1]*21 + [0]*21
    })

    lr = logrank_test(df_bench[df_bench["is_placebo"]==1]["time"],
                      df_bench[df_bench["is_placebo"]==0]["time"],
                      event_observed_A=df_bench[df_bench["is_placebo"]==1]["event"],
                      event_observed_B=df_bench[df_bench["is_placebo"]==0]["event"])

    cph = CoxPHFitter()
    cph.fit(df_bench, duration_col="time", event_col="event")
    hr_placebo = float(cph.hazard_ratios_["is_placebo"])

    return {
        "benchmark_logrank_chi2": float(lr.test_statistic),
        "benchmark_logrank_p": float(lr.p_value),
        "benchmark_efron_hr": hr_placebo
    }


if __name__ == "__main__":
    csv_file = Path(__file__).parent.parent / "data" / "mskcc_2020_real_cohort.csv"
    if not csv_file.exists():
        csv_file = Path(__file__).parent / "06_DATI_CLINICI_MSKCC_604_PAZIENTI.csv"

    print("=" * 80)
    print("  MSKCC CLINICAL VALIDATION ENGINE BENCHMARK")
    print("=" * 80)
    validator = MSKCCClinicalValidator(csv_file)
    fu = validator.compute_followup_and_survival_metrics()
    print(f"Reverse Kaplan-Meier Median Follow-Up : {fu['reverse_km_median_followup_months']:.2f} months")
    print(f"Raw Observed Median Time              : {fu['raw_median_rfs_months']:.2f} months")
    print(f"Patients at risk at 36 months         : {fu['patients_at_risk_36_months']}")
    print(f"Patients at risk at 60 months         : {fu['patients_at_risk_60_months']}")

    print("\nPairwise Cox Models vs Lepidic (Stage-Adjusted, Holm Multiplicity Correction):")
    cox = validator.fit_multivariable_cox_models()
    for s, res in cox["pairwise_contrasts_vs_lepidic"].items():
        print(f"  {s:15s}: HR = {res['hr']:.2f} (95% CI: {res['ci_lower']:.2f}-{res['ci_upper']:.2f}), p = {res['p_value']:.4f}, Holm p = {res['holm_p_value']:.4f}")

    print("\nHeadline Extreme Contrast (Solid/Micro vs Lepidic, N=193):")
    c1 = cox["model1_extreme_contrast"]
    print(f"  Adjusted HR = {c1['adjusted_hr']:.2f} (95% CI: {c1['ci_lower']:.2f}-{c1['ci_upper']:.2f}), p = {c1['p_value']:.4e}")

    print("\nFull Cohort Contrast (Solid/Micro vs All Others, N=604):")
    c2 = cox["model2_full_cohort"]
    print(f"  Adjusted HR = {c2['adjusted_hr']:.2f} (95% CI: {c2['ci_lower']:.2f}-{c2['ci_upper']:.2f}), p = {c2['p_value']:.4f}")

    print("\nEGFR-Mutant Subset Extreme Contrast (Stage Confounded):")
    ce = cox["egfr_subgroup_extreme"]
    print(f"  Adjusted HR = {ce['adjusted_hr']:.2f} (95% CI: {ce['ci_lower']:.2f}-{ce['ci_upper']:.2f}), p = {ce['p_value']:.4f}")

    print("\nClassical Statistical Benchmark Validation (Freireich 6-MP):")
    bench = run_freireich_benchmark_validation()
    print(f"  Freireich Log-Rank Chi2 = {bench['benchmark_logrank_chi2']:.2f} (Reference: 16.79)")
    print(f"  Freireich Efron HR      = {bench['benchmark_efron_hr']:.2f} (Reference: 4.82)")
