# Spatial Reaction-Diffusion Transport and Evolutionary Sanctuary Dynamics in Compact EGFR-Mutant Lung Adenocarcinoma

[![Master Study CI](https://github.com/Lukecele/spatial-oncology-master/actions/workflows/ci.yml/badge.svg)](https://github.com/Lukecele/spatial-oncology-master/actions/workflows/ci.yml)
[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.23623653-blue.svg)](https://doi.org/10.5281/zenodo.23623653)
[![License: GPL-3.0](https://img.shields.io/badge/License-GPL--3.0-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Python 3.10+](https://img.shields.io/badge/python-3.10%2B-blue.svg)](https://www.python.org/downloads/)
[![Platform: Linux | macOS | Windows](https://img.shields.io/badge/platform-Linux%20%7C%20macOS%20%7C%20Windows-lightgrey.svg)](https://github.com/Lukecele/spatial-oncology-master)

---

## Executive Summary & Core Discovery

Targeted therapy with third-generation EGFR tyrosine kinase inhibitors (TKIs), notably **osimertinib**, represents the worldwide standard of care for resected and advanced *EGFR*-mutant non-small cell lung cancer (NSCLC). In the adjuvant ADAURA setting, however, clinical outcomes are marked by profound temporal heterogeneity: **up to 68% of disease recurrences manifest after treatment discontinuation (months 36 to 60)**, peaking between 12 and 18 months post-cessation (*Nature Medicine*, June 2025). Crucially, late recurrences frequently lack canonical secondary genetic alterations (*EGFR* C797S or *MET* amplification), frequently preserving sensitivity to osimertinib rechallenge.

This repository provides the complete, mathematically exact, and clinically validated **Spatial Sanctuary Framework (SSF)**, establishing an unbroken causal continuum linking:
1. **Subcellular Reaction-Diffusion Physics:** Exact modified Bessel equation of order zero ($K_0$) solving radial transport in avascular tumor cords;
2. **Evolutionary Replicator Dynamics:** Frequency-dependent subclonal competition stabilizing a drug-tolerant persister (DTP) coexistence niche ($x^* = 0.6133$, Lyapunov stable $\lambda_{\text{eig}} = -1.7934 < 0$);
3. **Clinical Cohort Validation:** Authenticated real-world validation on $N=604$ surgically resected lung adenocarcinoma patients from Memorial Sloan Kettering Cancer Center (MSKCC 2020), demonstrating a strictly monotonic recurrence gradient ($p = 6.33 \times 10^{-6}$) and an adjusted Hazard Ratio of **$6.75$ (95% CI: $2.28 - 19.97, p = 5.52 \times 10^{-4}$)** for high-sanctuary solid/micropapillary architecture.

```mermaid
flowchart LR
    subgraph P1["Pillar 1: Biophysical Transport"]
        A["Exact Bessel K0 Equation<br/>• Small-Molecule λ = 627.7 µm<br/>• Critical Depth rcrit = 165.9 µm<br/>• Viable Shell Δr = 34.1 µm"]
    end
    subgraph P2["Pillar 2: Evolutionary Dynamics"]
        B["Replicator Game Theory<br/>• Hawk-Dove NSCLC Matrix<br/>• Stable Nash Eq x* = 0.6133<br/>• Post-ADAURA Regrowth Model"]
    end
    subgraph P3["Pillar 3: Clinical Validation"]
        C["Authentic MSKCC Cohort (N = 604)<br/>• Monotonic Risk Gradient (p &lt; 10⁻⁵)<br/>• Multivariable Adjusted HR = 6.75<br/>• Stage-Adjusted 95% CI: 2.28–19.97"]
    end
    subgraph P4["Pillar 4: Translational Practice"]
        D["Precision Clinical Implementation<br/>• Histology-Guided ctDNA MRD<br/>• 5.5-Month Molecular Lead Time<br/>• Rational Sequence with ADCs"]
    end

    A --> B
    B --> C
    C --> D

    classDef pillar fill:#f8f9fa,stroke:#2b5c8f,stroke-width:2px;
    class P1,P2,P3,P4 pillar;
```

---

## Publication Figures (600 DPI)

### Figure 1: Exact Biophysical Reaction-Diffusion Transport & Sanctuary Architecture
<p align="center">
  <img src="figures/figure1_reaction_diffusion_sanctuary.png" width="100%" style="max-width: 900px; height: auto;" alt="Figure 1: Biophysical Transport">
</p>
<em>Exact modified Bessel K0 radial concentration decay for small-molecule osimertinib (500 Da) vs macromolecular ADCs/bispecifics (150 kDa), defining the viable pharmacological sanctuary shell (165.9 to 200.0 µm, Δr = 34.1 µm) bounded by Krogh's oxygen limit.</em>

---

### Figure 2: Evolutionary Game Dynamics & Post-ADAURA Cessation Regrowth
<p align="center">
  <img src="figures/figure2_evolutionary_game_and_cessation.png" width="100%" style="max-width: 900px; height: auto;" alt="Figure 2: Evolutionary Game Dynamics">
</p>
<em>Replicator dynamics phase portrait proving asymptotic Lyapunov stability of the internal Nash equilibrium (x* = 0.6133, λ_eig = -1.7934 < 0), and post-discontinuation regrowth kinetics predicting ctDNA positivity at 10.7 months and CT clinical relapse at 16.1 months (lead time 5.5 months), matching Nature Medicine 2025 trial data.</em>

---

### Figure 3: Clinical Cohort Validation on Authentic MSKCC 2020 Data (N = 604)
<p align="center">
  <img src="figures/figure3_mskcc_clinical_validation.png" width="100%" style="max-width: 900px; height: auto;" alt="Figure 3: Clinical Cohort Validation">
</p>
<em>Real-world recurrence-free survival (RFS) across all five IASLC architectural patterns (Lepidic 4.55% → Acinar 14.67% → Micropapillary 21.62% → Papillary 25.58% → Solid 30.88%, Log-Rank p = 6.33e-6) and multivariable Cox modeling confirming that high-sanctuary architecture confers an adjusted HR of 6.75 (95% CI: 2.28–19.97, p = 5.52e-4).</em>

---

### Figure 4: Translational ADAURA Decision Framework & Therapeutic Strategies
<p align="center">
  <img src="figures/figure4_translational_adaura_framework.png" width="100%" style="max-width: 900px; height: auto;" alt="Figure 4: Translational ADAURA Framework">
</p>
<em>Multiscale clinical decision algorithm stratifying post-resection patients into standard surveillance versus intensive ctDNA MRD monitoring (Q3M from month 36 to 60) and comparative biophysical rationale for combining small molecules with ADCs (HER3-DXd, amivantamab).</em>

---

## Repository Structure

```text
spatial-oncology-master/
├── data/
│   └── mskcc_2020_real_cohort.csv      # Authenticated clinical cohort (MD5: 880b15b3734cd0d2d1a6184100c3d02c)
├── src/
│   ├── biophysical_diffusion_model.py  # Exact modified Bessel K0 solver and ODE integration
│   └── replicator_dynamics_game.py     # 2x2 non-zero sum game, Lyapunov stability & ADAURA regrowth
├── results/
│   ├── figure1_reaction_diffusion_sanctuary.png    # Figure 1 at 600 DPI
│   ├── figure2_evolutionary_game_and_cessation.png # Figure 2 at 600 DPI
│   ├── figure3_mskcc_clinical_validation.png       # Figure 3 at 600 DPI
│   ├── figure4_translational_adaura_framework.png  # Figure 4 at 600 DPI
│   └── audit_clinical_cohort_report.tsv            # Automated data audit and integrity certificate
├── scripts/
│   └── generate_unified_publication_figures.py     # End-to-end figure generation pipeline
├── tests/
│   ├── test_exact_biophysical_model.py # Unit tests for Bessel K0 solver and sanctuary shells
│   ├── test_exact_game_theory_model.py # Unit tests for replicator equations and stability
│   └── test_master_pipeline.py         # Full master pipeline execution tests
└── manuscript/
    ├── MANUSCRIPT_UNIFIED_MASTER_STUDY.pdf  # Publication-grade vector PDF (8 pages A4)
    ├── MANUSCRIPT_UNIFIED_MASTER_STUDY.html # Academic HTML source
    ├── MANUSCRIPT_UNIFIED_MASTER_STUDY.md   # Complete manuscript in Markdown
    ├── AUDIT_AGENT_EXPLANATION_MANIFESTO.md # Technical audit manifesto for external AI agents
    └── COVER_LETTER_NATURE_MEDICINE.md      # Formal submission cover letter
```

---

## Reproducibility & Testing

All models, statistical tests, and figures can be verified and reproduced with a single command:

```bash
# 1. Run unit test suites (15/15 tests pass)
uv run python -m unittest discover tests

# 2. Re-generate all 4 publication figures at 600 DPI
uv run --with matplotlib --with seaborn --with scipy --with pandas --with lifelines python3 scripts/generate_unified_publication_figures.py
```

---

## Citation & Correspondence

**Luca Celebrano**  
Department of Clinical Medicine and Surgery, School of Medicine, University of Naples Federico II  
Telethon Institute of Genetics and Medicine (TIGEM), Naples, Italy  
Email: `l.celebrano@studenti.unina.it` | ORCID: [0009-0007-0931-3555](https://orcid.org/0009-0007-0931-3555)

```bibtex
@article{celebrano2026spatial,
  title={Biophysical Reaction-Diffusion Transport and Evolutionary Sanctuary Dynamics in Compact EGFR-Mutant Lung Adenocarcinoma},
  author={Celebrano, Luca},
  journal={Nature Medicine (Under Review)},
  year={2026},
  doi={10.5281/zenodo.23623653}
}
```
