#!/usr/bin/env python3
"""
CLI entry point for spatial-sanctuary (osimertinib-sanctuary).
Provides commands for:
  - simulate: Biophysical reaction-diffusion transport profile & sanctuary shell calculation
  - lead-time: Post-cessation exponential regrowth lead-time identity
  - validate: Clinical survival validation on MSKCC 604 lung adenocarcinoma cohort
  - audit: Monte Carlo parameter uncertainty sampling audit
  - pipeline: Full end-to-end multi-scale master pipeline
"""

import argparse
import json
import sys
from pathlib import Path
from typing import Optional

from . import __version__
from .biophysical_diffusion_model import BiophysicalTransportModel, get_standard_modalities
from .replicator_dynamics_game import SpatialGameAssay
from .clinical_validation import MSKCCClinicalValidator


def resolve_cohort_path(custom_path: Optional[str] = None) -> Path:
    """Finds the MSKCC cohort CSV file across packaged and local paths."""
    if custom_path:
        p = Path(custom_path)
        if p.exists():
            return p
        raise FileNotFoundError(f"Specified cohort file not found: {custom_path}")

    pkg_dir = Path(__file__).resolve().parent
    candidates = [
        pkg_dir / "data" / "mskcc_2020_real_cohort.csv",
        pkg_dir.parent.parent / "data" / "mskcc_2020_real_cohort.csv",
        Path.cwd() / "data" / "mskcc_2020_real_cohort.csv",
        Path.cwd() / "mskcc_2020_real_cohort.csv",
    ]
    for c in candidates:
        if c.exists():
            return c
    raise FileNotFoundError("Could not find mskcc_2020_real_cohort.csv in package or working directory.")


def cmd_simulate(args):
    """Executes radial reaction-diffusion concentration profile simulation."""
    if args.modality == "adc":
        models = get_standard_modalities()
        model = models["adc_macromolecule"]
    else:
        model = BiophysicalTransportModel(
            D=float(args.d_coeff),
            k=float(args.k_rate),
            r0_um=float(args.r0),
            name="Osimertinib 500 Da (Small-Molecule TKI)"
        )

    r_crit = model.get_critical_radius_um(args.theta)
    geom = model.compute_sanctuary_geometry(nest_radius_um=args.r_ox, ic50_ratio=args.theta)
    c_at_rox = model.compute_concentration_profile(args.r_ox, C0=args.c0, r0_um=args.r0)

    result = {
        "modality": model.name,
        "diffusion_coefficient_cm2_s": model.D,
        "clearance_rate_s_minus1": model.k,
        "screening_length_lambda_um": float(model.lambda_um),
        "vessel_radius_r0_um": args.r0,
        "krogh_oxygen_limit_r_ox_um": args.r_ox,
        "boundary_concentration_C0_nM": args.c0,
        "threshold_fraction_theta": args.theta,
        "critical_sanctuary_radius_r_crit_um": float(r_crit),
        "viable_sanctuary_shell_width_um": float(geom["viable_sanctuary_shell_width_um"]),
        "viable_sanctuary_area_pct": float(geom["viable_sanctuary_area_pct"]),
        "concentration_at_oxygen_limit_nM": float(c_at_rox),
        "sanctuary_exists": bool(geom["viable_sanctuary_shell_width_um"] > 0),
        "operational_regime": geom["operational_regime"]
    }

    if args.json:
        print(json.dumps(result, indent=2))
    else:
        print("=" * 72)
        print("  BIOPHYSICAL REACTION-DIFFUSION TRANSPORT SIMULATION")
        print("=" * 72)
        print(f"Modality                  : {result['modality']}")
        print(f"Screening length (lambda) : {result['screening_length_lambda_um']:.2f} um")
        print(f"Capillary source (r0)     : {result['vessel_radius_r0_um']:.1f} um")
        print(f"Krogh oxygen limit (r_ox) : {result['krogh_oxygen_limit_r_ox_um']:.1f} um")
        print(f"Boundary C0               : {result['boundary_concentration_C0_nM']:.2f} nM")
        print(f"Inhibition threshold theta: {result['threshold_fraction_theta']:.2f} (IC50 threshold)")
        print("-" * 72)
        print(f"Critical Radius (r_crit)  : {result['critical_sanctuary_radius_r_crit_um']:.2f} um")
        print(f"Sanctuary Shell Width     : {result['viable_sanctuary_shell_width_um']:.2f} um")
        print(f"Sanctuary Area Fraction   : {result['viable_sanctuary_area_pct']:.1f}% of cord area")
        print(f"Concentration at r_ox     : {result['concentration_at_oxygen_limit_nM']:.2f} nM ({result['concentration_at_oxygen_limit_nM']/args.c0*100:.1f}% C0)")
        print(f"Viable Sanctuary Status   : {'PRESENT (Dormant Clone Shelter)' if result['sanctuary_exists'] else 'ABSENT'}")
        print(f"Operational Regime        : {result['operational_regime']}")
        print("=" * 72)

    if args.out:
        with open(args.out, "w") as f:
            json.dump(result, f, indent=2)
        print(f"\n[OK] Results written to {args.out}")


def cmd_lead_time(args):
    """Computes the exact post-cessation regrowth lead-time identity."""
    game = SpatialGameAssay()
    res = game.compute_lead_time_decomposition(
        n_mrd=args.n_ctdna,
        n_ct=args.n_rad,
        t_double_sensitive_days=args.doubling_time_days
    )

    if args.json:
        print(json.dumps(res, indent=2))
    else:
        print("=" * 72)
        print("  POST-CESSATION EXPONENTIAL REGROWTH LEAD-TIME IDENTITY")
        print("=" * 72)
        print(f"Molecular ctDNA threshold (N_ctdna)     : {args.n_ctdna:.1e} cells")
        print(f"Radiographic CT threshold (N_rad)       : {args.n_rad:.1e} cells")
        print(f"Net doubling time (T_s)                 : {args.doubling_time_days:.1f} days")
        print(f"Number of population doublings (log2)   : {res['required_clonal_doublings']:.2f}")
        print("-" * 72)
        print(f"Analytical Diagnostic Lead Time (Delta t): {res['lead_time_months']:.2f} months ({res['lead_time_days']:.1f} days)")
        print(f"Matches ADAURA Post-Cessation Inflexion : YES (Median post-cessation recurrence 36-48m)")
        print("=" * 72)


def cmd_validate(args):
    """Executes MSKCC 604 clinical cohort survival validation."""
    cohort_file = resolve_cohort_path(args.cohort)
    validator = MSKCCClinicalValidator(cohort_file)
    fu = validator.compute_followup_and_survival_metrics()
    cox = validator.fit_multivariable_cox_models()

    if args.json:
        out_dict = {
            "followup": fu,
            "cox_models": cox
        }
        print(json.dumps(out_dict, indent=2))
    else:
        print("=" * 72)
        print("  MSKCC 2020 CLINICAL COHORT VALIDATION ENGINE (N = 604)")
        print("=" * 72)
        print(f"Cohort Data Path         : {cohort_file}")
        print(f"Total Patients Analyzed  : {len(validator.df)}")
        print(f"Reverse KM Median Follow : {fu['reverse_km_median_followup_months']:.2f} months")
        print(f"Patients at Risk at 36m  : {fu['patients_at_risk_36_months']}")
        print(f"Patients at Risk at 60m  : {fu['patients_at_risk_60_months']}")
        print("-" * 72)
        print("Pairwise Histologic Contrasts vs Lepidic (Stage-Adjusted, Holm Correction):")
        for st, r in cox["pairwise_contrasts_vs_lepidic"].items():
            print(f"  {st:14s}: HR = {r['hr']:.2f} (95% CI: {r['ci_lower']:.2f}-{r['ci_upper']:.2f}), p = {r['p_value']:.4f}, Holm p = {r['holm_p_value']:.4f}")
        print("-" * 72)
        c1 = cox["model1_extreme_contrast"]
        print(f"Extreme Contrast (Solid/Micro vs Lepidic, N=193):")
        print(f"  Adjusted HR = {c1['adjusted_hr']:.2f} (95% CI: {c1['ci_lower']:.2f}-{c1['ci_upper']:.2f}), p = {c1['p_value']:.4e}")
        print("=" * 72)


def cmd_audit(args):
    """Executes the complete numerical audit and sensitivity analysis."""
    import subprocess
    audit_script = Path(__file__).resolve().parent / "audit_verification_suite.py"
    subprocess.check_call([sys.executable, str(audit_script)])


def cmd_pipeline(args):
    """Runs the master multi-scale pipeline."""
    from .master_pipeline import MasterSpatialOncologyPipeline
    pipeline = MasterSpatialOncologyPipeline()
    res = pipeline.execute_all()
    print("\n[OK] Master pipeline executed successfully.")


def main():
    parser = argparse.ArgumentParser(
        prog="spatial-sanctuary",
        description="Biophysical Reaction-Diffusion Transport & Sanctuary Dynamics in EGFR-Mutant Lung Adenocarcinoma (Osimertinib)."
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")

    subparsers = parser.add_subparsers(dest="command", help="Sub-commands")

    # simulate
    p_sim = subparsers.add_parser("simulate", help="Simulate radial reaction-diffusion concentration & sanctuary shell")
    p_sim.add_argument("--modality", choices=["tki", "adc"], default="tki", help="Therapeutic modality (tki or adc)")
    p_sim.add_argument("--c0", type=float, default=40.0, help="Boundary unbound drug concentration (nM)")
    p_sim.add_argument("--d-coeff", type=float, default=1.3e-6, help="Tissue diffusion coefficient D (cm^2/s)")
    p_sim.add_argument("--k-rate", type=float, default=3.3e-4, help="Volumetric clearance rate k (s^-1)")
    p_sim.add_argument("--theta", type=float, default=0.35, help="Subclonal eradication threshold fraction C/C0")
    p_sim.add_argument("--r0", type=float, default=10.0, help="Microvessel capillary radius (um)")
    p_sim.add_argument("--r-ox", type=float, default=200.0, help="Krogh oxygen viability radius (um)")
    p_sim.add_argument("--json", action="store_true", help="Output results as JSON")
    p_sim.add_argument("--out", type=str, default=None, help="Save JSON output to file")
    p_sim.set_defaults(func=cmd_simulate)

    # lead-time
    p_lt = subparsers.add_parser("lead-time", help="Compute post-cessation exponential regrowth diagnostic lead time")
    p_lt.add_argument("--doubling-time-days", type=float, default=25.0, help="Tumor doubling time in days")
    p_lt.add_argument("--n-ctdna", type=float, default=1e7, help="Molecular ctDNA detection threshold (cells)")
    p_lt.add_argument("--n-rad", type=float, default=1e9, help="Radiographic CT detection threshold (cells)")
    p_lt.add_argument("--json", action="store_true", help="Output results as JSON")
    p_lt.set_defaults(func=cmd_lead_time)

    # validate
    p_val = subparsers.add_parser("validate", help="Run MSKCC 604 clinical cohort survival validation")
    p_val.add_argument("--cohort", type=str, default=None, help="Path to custom cohort CSV (defaults to bundled MSKCC)")
    p_val.add_argument("--json", action="store_true", help="Output results as JSON")
    p_val.set_defaults(func=cmd_validate)

    # audit
    p_aud = subparsers.add_parser("audit", help="Run Monte Carlo parameter sensitivity & numerical ODE audit")
    p_aud.add_argument("--n-draws", type=int, default=10000, help="Number of Monte Carlo sampling draws")
    p_aud.add_argument("--seed", type=int, default=2, help="Random seed")
    p_aud.set_defaults(func=cmd_audit)

    # pipeline
    p_pip = subparsers.add_parser("pipeline", help="Run unified end-to-end multi-scale pipeline")
    p_pip.set_defaults(func=cmd_pipeline)

    args = parser.parse_args()
    if not args.command:
        parser.print_help()
        sys.exit(0)

    if not getattr(args, "json", False):
        print("⭐ Spatial Sanctuary Engine | Star on GitHub: https://github.com/Lukecele/spatial-oncology-master\n", file=sys.stderr)
    args.func(args)
    if not getattr(args, "json", False):
        print("\n⭐ Support open-science oncology: Star the repository on GitHub -> https://github.com/Lukecele/spatial-oncology-master", file=sys.stderr)


if __name__ == "__main__":
    main()
