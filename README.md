# Conditions for Pharmacological Sanctuaries in Perivascular Tumour Cords: A Reaction–Diffusion Analysis, a Regrowth-Kinetics Identity, and a Descriptive Survival Analysis of Histologic Subtype in Resected Lung Adenocarcinoma

[![GitHub stars](https://img.shields.io/github/stars/Lukecele/spatial-oncology-master?style=for-the-badge&logo=github&color=gold)](https://github.com/Lukecele/spatial-oncology-master)
[![Research Square DOI](https://img.shields.io/badge/Research%20Square-10.21203%2Frs.3.rs--11530312%2Fv1-007acc.svg)](https://doi.org/10.21203/rs.3.rs-11530312/v1)
[![Master Study CI](https://github.com/Lukecele/spatial-oncology-master/actions/workflows/ci.yml/badge.svg)](https://github.com/Lukecele/spatial-oncology-master/actions/workflows/ci.yml)
[![Journal: Journal of Mathematical Biology](https://img.shields.io/badge/Submitted-Journal%20of%20Mathematical%20Biology%20(Springer)-007acc.svg)](https://www.springer.com/journal/285)
[![Manuscript ID: JOMB-D-26-01395](https://img.shields.io/badge/Springer%20Nature-Editor%20Assigned%20(JOMB--D--26--01395)-brightgreen.svg)](https://www.editorialmanager.com/jomb/)
[![DOI: Zenodo](https://zenodo.org/badge/DOI/10.5281/zenodo.23623653.svg)](https://doi.org/10.5281/zenodo.23623653)
[![Bioconda Recipe](https://img.shields.io/badge/bioconda-spatial--sanctuary%20v2.2.2-brightgreen.svg?logo=conda-forge&logoColor=white)](https://anaconda.org/bioconda/spatial-sanctuary)
[![Docker Image](https://img.shields.io/badge/ghcr.io-spatial--oncology--master-blue?logo=docker&logoColor=white)](https://github.com/Lukecele/spatial-oncology-master/pkgs/container/spatial-oncology-master)
[![Substack Article](https://img.shields.io/badge/Substack-Article%20Live-FF6719?logo=substack&logoColor=white)](https://lucacelebrano.substack.com/p/why-resected-lung-cancers-relapse)
[![License: GPL-3.0](https://img.shields.io/badge/Code_License-GPL--3.0-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.html)
[![Documentation: CC BY 4.0](https://img.shields.io/badge/Doc_License-CC_BY_4.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)
[![Python >=3.9 (CI: 3.10–3.13)](https://img.shields.io/badge/python-%3E%3D3.9%20(CI%3A%203.10%E2%80%933.13)-blue.svg)](https://www.python.org/downloads/)
[![Interactive Web Portal](https://img.shields.io/badge/Web%20Portal-Live%20Demo-brightgreen.svg)](https://lukecele.github.io/spatial-oncology-master/)
[![Hugging Face Dataset](https://img.shields.io/badge/%F0%9F%A4%97%20Hugging%20Face-Dataset-ffd21e)](https://huggingface.co/datasets/Lucace/mskcc-lung-adenocarcinoma-spatial-cohort)
[![Hugging Face Spaces](https://img.shields.io/badge/%F0%9F%A4%97%20Hugging%20Face-Simulator%20%26%20Bedside%20CRI-blue)](https://huggingface.co/spaces/Lucace/spatial-oncology-simulator)
[![Expanded Study: Spatial Continuum](https://img.shields.io/badge/Expanded%20Study-spatial--continuum--oncology-darkred.svg)](https://github.com/Lukecele/spatial-continuum-oncology)
[![Preprint: medRxiv](https://img.shields.io/badge/medRxiv-Under%20Screening%20(MEDRXIV%2F2026%2F366778)-00599c.svg)](https://submit.medrxiv.org)

---

> [!NOTE]
> 📄 **Permanent Preprint:** Posted on Research Square (Springer Nature), DOI: [`10.21203/rs.3.rs-11530312/v1`](https://doi.org/10.21203/rs.3.rs-11530312/v1).  
> 🌐 **Interactive Scientific Portal:** Explore the live interactive manuscript, reaction-diffusion simulation, and survival analyses at [lukecele.github.io/spatial-oncology-master/](https://lukecele.github.io/spatial-oncology-master/).  
> **Manuscript Under Review:** Submitted to the *Journal of Mathematical Biology* (Springer Nature), Manuscript ID: **`JOMB-D-26-01395`** (*Editor Assigned*).  
> ⭐️ **Support Open Science:** If you find this computational framework or MSKCC survival analysis useful, please consider giving this repository a star on GitHub!  
> **Expanded Continuum Study:** See [Lukecele/spatial-continuum-oncology](https://github.com/Lukecele/spatial-continuum-oncology) for the grand unified biophysical reaction-diffusion theory, post-ADAURA cessation dynamics, and the SG-ADAURA clinical trial protocol.  
> **Hugging Face Dataset:** [`Lucace/mskcc-lung-adenocarcinoma-spatial-cohort`](https://huggingface.co/datasets/Lucace/mskcc-lung-adenocarcinoma-spatial-cohort) ($N = 604$ MSKCC cohort & multi-center validation).  
> **Hugging Face Space:** [`Lucace/spatial-oncology-simulator`](https://huggingface.co/spaces/Lucace/spatial-oncology-simulator) (Interactive reaction-diffusion and clinical survival simulator).  
> **Substack Publication:** [Why Resected Lung Cancers Relapse Post-Osimertinib Without Resistance Mutations](https://lucacelebrano.substack.com/p/why-resected-lung-cancers-relapse).  
> **Permanent Archival DOI:** Zenodo [`10.5281/zenodo.23623653`](https://doi.org/10.5281/zenodo.23623653).  
> **Author:** Luca Celebrano (Department of Clinical Medicine and Surgery, School of Medicine, University of Naples Federico II) | ORCID: [`0009-0007-0931-3555`](https://orcid.org/0009-0007-0931-3555).

> [!TIP]
> **Support Open Science:**  
> If you find this unified biophysical model, analytical lead-time identity, or clinical survival validation helpful in your research, please consider giving this repository a ⭐ **Star** on GitHub!

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
│   ├── mskcc_2020_real_cohort.csv      # MSKCC clinical cohort (N=604, SHA-256 verified)
│   └── audit_clinical_cohort_report.tsv # Cohort baseline distributions and hash certificate
├── src/
│   ├── spatial_sanctuary/              # Installable Python package (CLI: spatial-sanctuary)
│   ├── biophysical_diffusion_model.py  # Exact K0 (outward) & I0 (inward) Bessel solvers, O2 transport, Monte Carlo
│   ├── replicator_dynamics_game.py     # 2x2 matrix game, Lyapunov stability & analytical lead-time identity
│   ├── clinical_validation.py          # Reverse KM, pairwise Cox models with 10-fold Holm adjustment, Freireich benchmark
│   ├── master_pipeline.py              # End-to-end multi-scale pipeline orchestrator
│   └── audit_verification_suite.py     # Automated numerical audit script
├── tests/                              # Unit tests (Bessel solutions, replicator dynamics, cohort statistics, CLI)
├── scripts/                            # Demo script and publication figure generator
├── notebooks/                          # Interactive demonstration notebook
├── figures/                            # Publication figures (600 DPI)
├── manuscript/                         # Manuscript (PDF, Markdown, HTML, Word) and figures
├── docs/                               # Interactive web portal (GitHub Pages)
├── .github/workflows/                  # CI (Python 3.10-3.13) and Docker image publishing
├── CITATION.cff                        # Machine-readable citation metadata
├── Dockerfile                          # Container image (GHCR)
├── pyproject.toml                      # Package build configuration
├── LICENSE                             # GNU General Public License v3.0 (GPL-3.0)
└── requirements.txt                    # Minimal dependencies
```

---

## Reproducibility & Continuous Integration

This repository is continuously tested across Python 3.10, 3.11, 3.12, and 3.13 via GitHub Actions.

### Installation

#### Option 1: Python Package (pip)
```bash
git clone https://github.com/Lukecele/spatial-oncology-master.git
cd spatial-oncology-master
pip install -e .
```

#### Option 2: Docker Container
```bash
# Pull and run directly from GitHub Container Registry
docker run -it ghcr.io/lukecele/spatial-oncology-master:latest simulate
```

---

### Command-Line Interface (CLI)

The package installs both `spatial-sanctuary` and `osimertinib-sanctuary` executables:

```bash
# 1. Simulate radial reaction-diffusion concentration & sanctuary shell width
spatial-sanctuary simulate --c0 40.0 --lambda-um 627.65 --theta 0.35

# 2. Compute exact post-cessation regrowth lead-time identity (Delta t)
spatial-sanctuary lead-time --doubling-time-days 25.0

# 3. Validate MSKCC 604 clinical cohort survival models
spatial-sanctuary validate

# 4. Run 10,000 Monte Carlo draws parameter sensitivity audit
spatial-sanctuary audit

# 5. Run full end-to-end master pipeline
spatial-sanctuary pipeline

# Interactive zero-dependency diagnostic demo
./quickstart.sh
```

---

## License

- **Software and Code:** Released under the [GNU General Public License v3.0](https://www.gnu.org/licenses/gpl-3.0.html) (`GPL-3.0`).
- **Manuscript, Figures, and Documentation:** Licensed under the [Creative Commons Attribution 4.0 International License](https://creativecommons.org/licenses/by/4.0/) (`CC-BY 4.0`).

---

## Citation

Please cite the archived software (Zenodo). The associated manuscript is under peer review and will be added here once published.

```bibtex
@software{celebrano2026spatialsanctuary,
  title     = {Spatial Oncology Master Study: reaction-diffusion transport, evolutionary sanctuary dynamics and clinical cohort analysis in EGFR-mutant lung adenocarcinoma},
  author    = {Celebrano, Luca},
  year      = {2026},
  publisher = {Zenodo},
  doi       = {10.5281/zenodo.23623653},
  url       = {https://doi.org/10.5281/zenodo.23623653}
}
```

GitHub's "Cite this repository" button reads the same information from [`CITATION.cff`](CITATION.cff).
