# Conditions for Pharmacological Sanctuaries in Perivascular Tumour Cords: A Reaction–Diffusion Analysis, a Regrowth-Kinetics Identity, and a Descriptive Survival Analysis of Histologic Subtype in Resected Lung Adenocarcinoma

[![Master Study CI](https://github.com/Lukecele/spatial-oncology-master/actions/workflows/ci.yml/badge.svg)](https://github.com/Lukecele/spatial-oncology-master/actions/workflows/ci.yml)
[![Journal: Journal of Mathematical Biology](https://img.shields.io/badge/Submitted-Journal%20of%20Mathematical%20Biology%20(Springer)-007acc.svg)](https://www.springer.com/journal/285)
[![Manuscript ID: JOMB-D-26-01395](https://img.shields.io/badge/Springer%20Nature%20ID-JOMB--D--26--01395-orange.svg)](https://www.editorialmanager.com/jomb/)
[![Preprint: bioRxiv](https://img.shields.io/badge/bioRxiv-BIORXIV%2F2026%2F761610-b31b1b.svg)](https://submit.biorxiv.org/)
[![DOI: Zenodo](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.23623653-blue.svg)](https://doi.org/10.5281/zenodo.23623653)
[![License: GPL-3.0](https://img.shields.io/badge/Code_License-GPL--3.0-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.html)
[![Documentation: CC BY 4.0](https://img.shields.io/badge/Doc_License-CC_BY_4.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)
[![Python 3.10 | 3.11 | 3.12 | 3.13](https://img.shields.io/badge/python-3.10%20%7C%203.11%20%7C%203.12%20%7C%203.13-blue.svg)](https://www.python.org/downloads/)
[![Interactive Web Portal](https://img.shields.io/badge/Web%20Portal-Live%20Demo-brightgreen.svg)](https://lukecele.github.io/spatial-oncology-master/)

---

> [!NOTE]
> 🌐 **Interactive Scientific Portal:** Explore the live interactive manuscript, reaction-diffusion simulation, and survival analyses at [lukecele.github.io/spatial-oncology-master/](https://lukecele.github.io/spatial-oncology-master/).  
> **Manuscript Under Review:** Submitted to the *Journal of Mathematical Biology* (Springer Nature), Manuscript ID: **`JOMB-D-26-01395`** (*Under Consideration*).  
> **Preprint Tracking:** bioRxiv screening in progress (**`BIORXIV/2026/761610`**).  
> **Permanent Archival DOI:** Zenodo [`10.5281/zenodo.23623653`](https://doi.org/10.5281/zenodo.23623653).  
> **Author:** Luca Celebrano (School of Medicine, University of Naples Federico II & TIGEM) | ORCID: [`0009-0007-0931-3555`](https://orcid.org/0009-0007-0931-3555).

---

## Executive Summary

Targeted adjuvant therapy with the third-generation EGFR tyrosine kinase inhibitor (TKI) **osimertinib** produces profound disease-free survival benefits in resected *EGFR*-mutated non-small-cell lung cancer (NSCLC) (ADAURA trial; Wu et al. 2020, Tsuboi et al. 2023). However, clinical recurrence dynamics exhibit important temporal features: in longitudinal molecular residual disease (MRD) surveillance (Herbst et al., *Nature Medicine* 2025), **67.9% (19/28) of relapses occurred after planned therapy cessation at year 3** (months 36 to 60), with **57.9% of post-cessation events (11/19, or 39% of all relapses) concentrating within the first 12 months** (months 36–48). Crucially, these recurrences frequently lack secondary genomic gatekeeper mutations (*EGFR* C797S or *MET* amplification), retaining sensitivity to TKI rechallenge.

This repository hosts the audited, reproducible computational architecture investigating the biophysical, kinetic, and clinical foundations of microenvironmental drug distribution in lung adenocarcinoma:

1. **Biophysical Reaction-Diffusion Transport:** Exact cylindrical Helmholtz/Bessel solutions demonstrating that:
   - In **outward perivascular cords** (divergent geometry, modified Bessel $K_0$), osimertinib concentration drops below $\text{IC}_{50}$ ($\theta = 0.35$) at $r_{\text{crit}} = 165.9\,\mu\text{m}$, creating a viable pharmacological sanctuary shell of width $\Delta r = 34.1\,\mu\text{m}$ up to the Krogh oxygenation boundary ($r_{\text{ox}} \approx 200\,\mu\text{m}$).
   - In **inward avascular nests** (convergent geometry, modified Bessel $I_0$), osimertinib penetrates to $90.6\%\text{--}99.4\%$ of boundary levels throughout viable tissue, proving that avascular nests do not shelter small molecules. Sanctuaries in nests are restricted to high-molecular-weight agents such as antibody-drug conjugates (ADCs, $\lambda = 77.5\,\mu\text{m}$).
   - Monte Carlo sensitivity sampling ($N = 10,000$, seed = 2, log-normal priors) establishes that an outward sanctuary shell exists in 46.1% of parameter draws (62.8% at fixed $r_{\text{ox}} = 200\,\mu\text{m}$).

2. **Analytical Lead-Time Identity:** Mathematical proof that the diagnostic lead-time window between molecular ctDNA detection ($10^7$ cells) and radiographic CT imaging ($10^9$ cells) is an invariant algebraic identity of exponential growth:
   $$\Delta t = T_s \cdot \log_2(100) = 6.6439 \cdot T_s = 5.46\,\text{months (for } T_s = 25.0\,\text{days)}$$
   This lead time is an invariant property of cellular doubling kinetics and assay thresholds, rather than an empirical proof of spatial sanctuary physics.

3. **Descriptive Clinical Survival Analysis ($N = 604$ MSKCC Cohort):**
   - Reverse Kaplan–Meier method (Schemper & Smith 1996) determines a true median follow-up of **23.8 months** (IQR: 13.3–31.8 months; 5-group multivariate log-rank $\chi^2 = 24.52, p = 6.27 \times 10^{-5}$).
   - In stage-adjusted multivariable Cox modeling, lepidic architecture confers strong baseline protection ($HR = 0.25, 95\%\,\text{CI}: 0.09\text{--}0.70, p = 0.0084$ vs acinar).
   - In 10-fold Holm-adjusted pairwise comparisons, lepidic recurrence is significantly lower than solid ($p_{\text{adj}} = 0.0026$) and papillary ($p_{\text{adj}} = 0.0203$), while differences vs acinar ($p_{\text{adj}} = 0.0703$) and micropapillary ($p_{\text{adj}} = 0.3911$, unadjusted $p = 0.0559$) do not retain significance under strict family-wise error control.
   - The exploratory contrast comparing solid/micropapillary against lepidic ($HR = 6.75, 95\%\,\text{CI}: 2.28\text{--}19.97$) rests on only 4 lepidic events (bootstrap upper bound is unbounded), whereas across the unselected cohort the hazard ratio is $1.32$ ($p = 0.2384$).

---

## Repository Structure

```text
spatial-oncology-master/
├── data/
│   ├── mskcc_2020_real_cohort.csv      # Authenticated MSKCC clinical cohort (N=604, SHA-256 verified)
│   └── audit_clinical_cohort_report.tsv # Cohort baseline distributions and hash certificate
├── src/
│   ├── biophysical_diffusion_model.py  # Exact K0 (outward) & I0 (inward) Bessel solvers, O2 transport, Monte Carlo
│   ├── replicator_dynamics_game.py     # 2x2 matrix game, Lyapunov stability & analytical lead-time identity
│   ├── clinical_validation.py          # Reverse KM, pairwise Cox models with 10-fold Holm adjustment, Freireich benchmark
│   ├── master_pipeline.py              # End-to-end multi-scale pipeline orchestrator
│   └── audit_verification_suite.py     # Complete automated numerical audit script
├── tests/
│   ├── test_exact_biophysical_model.py # Unit tests for Bessel ODE solutions & boundary limits
│   ├── test_exact_game_theory_model.py # Unit tests for Replicator dynamics, Nash attractor & lead time
│   └── test_master_clinical_validation.py # Unit tests for MSKCC cohort statistics, log-rank & Cox models
├── figures/
│   ├── Figure1_transport_sensitivity.png # Figure 1: Bessel transport profiles & 2D sensitivity heatmap
│   ├── Figure2_regrowth_identity.png     # Figure 2: Post-cessation regrowth & lead-time identity
│   └── Figure3_cohort_cox.png            # Figure 3: MSKCC Kaplan-Meier curves & Cox models forest plot
├── manuscript/
│   ├── MANUSCRIPT_UNIFIED_MASTER_STUDY.pdf # Official 17-page vector publication PDF with MathJax
│   ├── MANUSCRIPT_UNIFIED_MASTER_STUDY.md  # Full source text in GitHub Flavored Markdown
│   ├── MANUSCRIPT_UNIFIED_MASTER_STUDY.html # High-fidelity standalone interactive web edition
│   └── 09_MANOSCRITTO_FORMATO_WORD.docx   # Publication manuscript in Word (.docx) with embedded figures
├── .github/workflows/
│   └── ci.yml                          # Continuous Integration matrix (Python 3.10, 3.11, 3.12, 3.13)
├── CITATION.cff                        # Machine-readable academic citation metadata
├── LICENSE                             # GNU General Public License v3.0 (GPL-3.0)
└── requirements.txt                    # Deterministic minimal dependencies
```

---

## Reproducibility & Continuous Integration

This repository is continuously tested across Python 3.10, 3.11, 3.12, and 3.13 via GitHub Actions.

To execute the test suite locally:

```bash
git clone https://github.com/Lukecele/spatial-oncology-master.git
cd spatial-oncology-master
pip install -r requirements.txt

# Run all 15 unit tests
python -m unittest discover -s tests -v

# Run end-to-end master pipeline
python src/master_pipeline.py

# Run complete numerical audit suite
python src/audit_verification_suite.py
```

---

## License

- **Software and Code:** Released under the [GNU General Public License v3.0](https://www.gnu.org/licenses/gpl-3.0.html) (`GPL-3.0`).
- **Manuscript, Figures, and Documentation:** Licensed under the [Creative Commons Attribution 4.0 International License](https://creativecommons.org/licenses/by/4.0/) (`CC-BY 4.0`).

---

## Citation

```bibtex
@article{celebrano2026conditions,
  title={Conditions for Pharmacological Sanctuaries in Perivascular Tumour Cords: A Reaction--Diffusion Analysis, a Regrowth-Kinetics Identity, and a Descriptive Survival Analysis of Histologic Subtype in Resected Lung Adenocarcinoma},
  author={Celebrano, Luca},
  journal={Journal of Mathematical Biology},
  year={2026},
  note={Manuscript ID: JOMB-D-26-01395. bioRxiv ID: BIORXIV/2026/761610. DOI: 10.5281/zenodo.23623653}
}
```
