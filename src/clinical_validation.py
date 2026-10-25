#!/usr/bin/env python3
"""
MSKCC 2020 Clinical Cohort Validation Engine (N = 604).
Analyzes architectural sanctuary patterns, IASLC histological subtypes,
and recurrence-free survival (RFS) in surgically resected lung adenocarcinoma.
"""

from pathlib import Path
from typing import Dict, Any, Tuple
import numpy as np
import pandas as pd
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import logrank_test


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

    def compute_subtype_gradient(self) -> Dict[str, Dict[str, float]]:
        """Computes observed recurrence rates across the 5 IASLC histological subtypes."""
        order = ["Lepidic", "Acinar", "Papillary", "Micropapillary", "Solid"]
        results = {}
        for s in order:
            sub = self.df[self.df["PREDOMINANT_HISTOLOGIC_SUBTYPE"] == s]
            n = len(sub)
            rec = int(sub["RECURRENCE_EVENT"].sum()) if n > 0 else 0
            rate = float(rec / n) if n > 0 else 0.0
            results[s] = {"N": n, "recurrences": rec, "recurrence_rate": rate}
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
        Computes follow-up duration distributions and Kaplan-Meier metrics,
        auditing statistical power and tail reliability.
        """
        median_fu = float(self.df["RFS_MONTHS"].median())
        q25 = float(self.df["RFS_MONTHS"].quantile(0.25))
        q75 = float(self.df["RFS_MONTHS"].quantile(0.75))
        n_at_36 = int((self.df["RFS_MONTHS"] >= 36).sum())
        n_at_60 = int((self.df["RFS_MONTHS"] >= 60).sum())

        # Lepidic follow-up
        lep_df = self.df[self.df["PREDOMINANT_HISTOLOGIC_SUBTYPE"] == "Lepidic"]
        n_lep_at_60 = int((lep_df["RFS_MONTHS"] >= 60).sum())

        # Solid median RFS
        solid_df = self.df[self.df["PREDOMINANT_HISTOLOGIC_SUBTYPE"] == "Solid"]
        kmf_solid = KaplanMeierFitter()
        kmf_solid.fit(solid_df["RFS_MONTHS"], solid_df["RECURRENCE_EVENT"])
        solid_median_rfs = float(kmf_solid.median_survival_time_)
        n_solid_at_40 = int((solid_df["RFS_MONTHS"] >= 40).sum())

        return {
            "median_rfs_followup_months": median_fu,
            "iqr_rfs_followup": (q25, q75),
            "patients_at_risk_36_months": n_at_36,
            "patients_at_risk_60_months": n_at_60,
            "lepidic_patients_at_risk_60_months": n_lep_at_60,
            "solid_median_rfs_months": solid_median_rfs,
            "solid_patients_at_risk_at_median": n_solid_at_40
        }

    def fit_multivariable_cox_models(self) -> Dict[str, Any]:
        """
        Fits multivariable Cox Proportional Hazards models:
        1. Model 1: Architectural Extreme Comparison (N=193: Solid/Micro vs Lepidic, adjusted for Stage)
        2. Model 2A: Full Cohort Multivariable (N=604: Solid/Micro vs Rest, Stage II/III vs I)
        3. Model 2B: Full Cohort Multivariable (N=604: Solid/Micro vs Rest, Stage II/III vs I, EGFR status)
        4. Model 3: 5-Subtype Model with Acinar as reference (N=604, adjusted for Stage)
        """
        # Model 1: Extreme Contrast (N=193)
        df_comp = self.df[self.df["PREDOMINANT_HISTOLOGIC_SUBTYPE"].isin(["Solid", "Micropapillary", "Lepidic"])].copy()
        df_comp["is_solid_micro"] = df_comp["PREDOMINANT_HISTOLOGIC_SUBTYPE"].isin(["Solid", "Micropapillary"]).astype(int)
        df_comp["stage_adv"] = (df_comp["PATHOLOGIC_STAGE"] >= 2).astype(int)

        cph1 = CoxPHFitter()
        cph1.fit(df_comp[["RFS_MONTHS", "RECURRENCE_EVENT", "is_solid_micro", "stage_adv"]],
                 duration_col="RFS_MONTHS", event_col="RECURRENCE_EVENT")
        s1 = cph1.summary.loc["is_solid_micro"]

        # Model 2A: Full Cohort Stage-Adjusted (N=604)
        df_full = self.df.copy()
        df_full["is_solid_micro"] = df_full["PREDOMINANT_HISTOLOGIC_SUBTYPE"].isin(["Solid", "Micropapillary"]).astype(int)
        df_full["stage_adv"] = (df_full["PATHOLOGIC_STAGE"] >= 2).astype(int)

        cph2a = CoxPHFitter()
        cph2a.fit(df_full[["RFS_MONTHS", "RECURRENCE_EVENT", "is_solid_micro", "stage_adv"]],
                  duration_col="RFS_MONTHS", event_col="RECURRENCE_EVENT")
        s2a = cph2a.summary.loc["is_solid_micro"]

        # Model 2B: Full Cohort Stage + EGFR Adjusted (N=604)
        cph2b = CoxPHFitter()
        cph2b.fit(df_full[["RFS_MONTHS", "RECURRENCE_EVENT", "is_solid_micro", "stage_adv", "IS_EGFR_MUTANT"]],
                  duration_col="RFS_MONTHS", event_col="RECURRENCE_EVENT")
        s2b = cph2b.summary.loc["is_solid_micro"]

        # Model 3: 5 Subtypes with Acinar Reference (N=604)
        df_sub = pd.get_dummies(df_full["PREDOMINANT_HISTOLOGIC_SUBTYPE"], prefix="sub", dtype=int)
        df_cox3 = pd.concat([df_full[["RFS_MONTHS", "RECURRENCE_EVENT", "stage_adv"]],
                             df_sub[["sub_Lepidic", "sub_Solid", "sub_Micropapillary", "sub_Papillary"]]], axis=1)
        cph3 = CoxPHFitter()
        cph3.fit(df_cox3, duration_col="RFS_MONTHS", event_col="RECURRENCE_EVENT")
        s3 = cph3.summary

        return {
            "model1_extreme": {
                "description": "Solid/Micropapillary vs Lepidic (N=193, stage-adjusted)",
                "n": len(df_comp),
                "events": int(df_comp["RECURRENCE_EVENT"].sum()),
                "adjusted_hr": float(s1["exp(coef)"]),
                "ci_lower": float(s1["exp(coef) lower 95%"]),
                "ci_upper": float(s1["exp(coef) upper 95%"]),
                "p_value": float(s1["p"])
            },
            "model2_full_cohort": {
                "description": "Solid/Micropapillary vs Rest (N=604, stage & EGFR adjusted)",
                "n": len(df_full),
                "events": int(df_full["RECURRENCE_EVENT"].sum()),
                "adjusted_hr": float(s2b["exp(coef)"]),
                "ci_lower": float(s2b["exp(coef) lower 95%"]),
                "ci_upper": float(s2b["exp(coef) upper 95%"]),
                "p_value": float(s2b["p"])
            },
            "model2a_full_cohort_stage_only": {
                "description": "Solid/Micropapillary vs Rest (N=604, stage-adjusted only)",
                "n": len(df_full),
                "events": int(df_full["RECURRENCE_EVENT"].sum()),
                "adjusted_hr": float(s2a["exp(coef)"]),
                "ci_lower": float(s2a["exp(coef) lower 95%"]),
                "ci_upper": float(s2a["exp(coef) upper 95%"]),
                "p_value": float(s2a["p"])
            },
            "model3_subtypes_acinar_ref": {
                "description": "5-subtype multivariable model with Acinar as reference (N=604, stage-adjusted)",
                "lepidic": {
                    "hr": float(s3.loc["sub_Lepidic", "exp(coef)"]),
                    "ci_lower": float(s3.loc["sub_Lepidic", "exp(coef) lower 95%"]),
                    "ci_upper": float(s3.loc["sub_Lepidic", "exp(coef) upper 95%"]),
                    "p_value": float(s3.loc["sub_Lepidic", "p"])
                },
                "solid": {
                    "hr": float(s3.loc["sub_Solid", "exp(coef)"]),
                    "ci_lower": float(s3.loc["sub_Solid", "exp(coef) lower 95%"]),
                    "ci_upper": float(s3.loc["sub_Solid", "exp(coef) upper 95%"]),
                    "p_value": float(s3.loc["sub_Solid", "p"])
                },
                "micropapillary": {
                    "hr": float(s3.loc["sub_Micropapillary", "exp(coef)"]),
                    "ci_lower": float(s3.loc["sub_Micropapillary", "exp(coef) lower 95%"]),
                    "ci_upper": float(s3.loc["sub_Micropapillary", "exp(coef) upper 95%"]),
                    "p_value": float(s3.loc["sub_Micropapillary", "p"])
                },
                "papillary": {
                    "hr": float(s3.loc["sub_Papillary", "exp(coef)"]),
                    "ci_lower": float(s3.loc["sub_Papillary", "exp(coef) lower 95%"]),
                    "ci_upper": float(s3.loc["sub_Papillary", "exp(coef) upper 95%"]),
                    "p_value": float(s3.loc["sub_Papillary", "p"])
                }
            }
        }


if __name__ == "__main__":
    v = MSKCCClinicalValidator(Path(__file__).parent / "06_DATI_CLINICI_MSKCC_604_PAZIENTI.csv")
    print("MSKCC 2020 Validation Engine Output:")
    print("Cohort summary:", v.get_cohort_summary())
    print("\nFollow-up metrics:", v.compute_followup_and_survival_metrics())
    print("\nCox models:", v.fit_multivariable_cox_models())
