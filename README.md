# From Subcellular Spatial Biophysics to Routine Digital Pathology AI: Resolving Non-Genetic Osimertinib Resistance via Microvascular Decoupling

[![Master Study CI](https://github.com/Lukecele/spatial-oncology-master/actions/workflows/ci.yml/badge.svg)](https://github.com/Lukecele/spatial-oncology-master/actions/workflows/ci.yml)
[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.23623653-blue.svg)](https://doi.org/10.5281/zenodo.23623653)
[![License: GPL-3.0](https://img.shields.io/badge/License-GPL--3.0-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Lukecele/spatial-oncology-master/blob/main/notebooks/InSituVision_Master_Pipeline.ipynb)
[![Python 3.10+](https://img.shields.io/badge/python-3.10%2B-blue.svg)](https://www.python.org/downloads/)
[![Platform: Linux | macOS | Windows | Colab](https://img.shields.io/badge/platform-Linux%20%7C%20macOS%20%7C%20Windows%20%7C%20Colab-lightgrey.svg)](https://github.com/Lukecele/spatial-oncology-master)

> [!TIP]
> **Support Open Science in Computational Oncology**  
> If this unified multi-scale framework, reaction-diffusion solver, clinical survival validation, or diagnostic Attention-MIL AI aids your research, please consider leaving a **⭐ star on this repository**—it helps accelerate community adoption of physical oncology and accessible digital pathology!

---

## Executive Summary & The Causal Continuum

Targeted therapy with third-generation EGFR tyrosine kinase inhibitors (TKIs), notably **osimertinib**, has transformed the management of sensitizing *EGFR*-mutated non-small cell lung cancer (NSCLC). However, therapeutic durability remains profoundly heterogeneous, and up to **50% of progressing patients display no detectable genetic resistance alterations** (such as *EGFR* C797S or *MET* amplification).

This repository houses the **Master Unified Study**, synthesizing four interconnected investigations into an unbroken causal continuum:

```mermaid
flowchart LR
    subgraph Pillar1["Pillar 1: Theory"]
        A["Spatial Game Theory<br/>• Space-for-Time Ergodicity<br/>• Polymodal Nash Eq x* ≈ 0.61<br/>• Frequency-Dependent Fitness"]
    end
    subgraph Pillar2["Pillar 2: Visium HD"]
        B["Subcellular Biophysics<br/>• 211,104 Bins @ 2 µm<br/>• Homotypic Crowding 87.9%<br/>• Reaction-Diffusion λ ≈ 624 µm"]
    end
    subgraph Pillar3["Pillar 3: Clinical Cohort"]
        C["Population Survival N = 854<br/>• Log-Rank p &lt; 10⁻¹⁵<br/>• Hazard Ratio HR = 19.45<br/>• Median PFS: 17.9 vs 38.1 mo"]
    end
    subgraph Pillar4["Pillar 4: Digital Pathology AI"]
        D["Routine Diagnostic H&E AI<br/>• Gated Attention-MIL N = 570<br/>• Harrell C-index ≈ 0.74<br/>• Cost under 5 € (-99.91%)"]
    end

    A --> B
    B --> C
    C --> D

    classDef pillar fill:#f8f9fa,stroke:#495057,stroke-width:2px;
    class Pillar1,Pillar2,Pillar3,Pillar4 pillar;
```

> [!IMPORTANT]
> **Democratizing Precision Oncology from 5,000 € to < 5 €**  
> 10x Genomics Visium HD spatial transcriptomics costs ~5,000 € per sample and is restricted to elite medical centers. By capturing the morphological imprint of microvascular collapse on routine, universal Hematoxylin & Eosin (H&E) diagnostic slides, **InSituVision achieves a 99.91% cost reduction ($4.50\text{ €}$ vs $5,000\text{ €}$)**, enabling clinical implementation in any pathology laboratory globally.

---

## Landmark Findings Across the 4 Pillars

| Analytical Pillar | Experimental / Computational Method | Primary Quantitative Finding | Clinical & Translational Significance |
| :--- | :--- | :--- | :--- |
| **Pillar 1: Game Theory** | Replicator dynamics & spatial neighborhood decomposition | Internal Polymodal Nash Eq ($x^* \approx 0.613$) | Resolves reversible drug sensitivity and post-holiday rebound without de novo mutations. |
| **Pillar 2: Visium HD Biophysics** | 10x Genomics Visium HD ($211,104$ bins at $2\,\mu\text{m}$) & Reaction-diffusion solver | Screening penetration depth $\lambda = \sqrt{D/k} \approx 624\,\mu\text{m}$ | Proves that microvascular collapse creates sub-therapeutic ($C < \text{IC}_{50}$) pharmacological sanctuaries. |
| **Pillar 3: Clinical Cohort** | Multi-center survival harmonization ($N = 854$, $N = 218$ *EGFR*-mutant) | Median PFS $17.9$ vs $38.1$ mo ($p < 10^{-15}$, $\text{HR} = 19.45$) | Confirms microvascular decoupling as a massive independent clinical survival biomarker. |
| **Pillar 4: Vision AI** | Diagnostic H&E Gated Attention-MIL ($N = 570$ WSIs) | Harrell C-index $\approx 0.74$, Slide Cost $4.50\text{ €}$ | Translates spatial biophysics into standard diagnostic histology, slashing assay costs by $99.91\%$. |

---

## Publication Figures & Visual Evidence

<div align="center">

### Figure 1: Spatial Evolutionary Game Theory & Polymodal Coexistence
<p align="center">
  <img src="figures/figure1_game_decomposition.png" width="100%" style="max-width: 900px; height: auto;" alt="Figure 1: Spatial Game Decomposition">
</p>
<p align="center"><em>Spatial cellular interaction graphs infer off-diagonal payoff dominance, establishing stable subclonal coexistence (x* ≈ 0.613).</em></p>

---

### Figure 2: Subcellular Visium HD Biophysics & The 624 µm Sanctuary Threshold
<p align="center">
  <img src="figures/figure2_subcellular_architecture_624um.png" width="100%" style="max-width: 900px; height: auto;" alt="Figure 2: Subcellular Visium HD Architecture">
</p>
<p align="center"><em>Subcellular spatial transcriptomics confirms homotypic neoplastic crowding (87.9%) and vascular collapse, defining the radial penetration screening threshold (λ ≈ 624 µm).</em></p>

---

### Figure 3: Population Clinical Survival Validation (N = 854 Patients)
<p align="center">
  <img src="figures/figure3_massive_survival_validation_854pts.png" width="100%" style="max-width: 900px; height: auto;" alt="Figure 3: Massive Survival Validation">
</p>
<p align="center"><em>Kaplan-Meier survival curves demonstrating massive PFS and OS divergence between decoupled sanctuary phenotypes and vascularized phenotypes (HR = 19.45, p < 10⁻¹⁵).</em></p>

---

### Figure 4 & Figure 5: Routine Diagnostic H&E Attention-MIL Survival Prediction & Spatial Heatmaps
<p align="center">
  <img src="figures/figure4_he_vision_survival_prediction.png" width="49%" style="max-width: 450px; height: auto;" alt="Figure 4: Survival Prediction">
  <img src="figures/figure5_spatial_attention_heatmap.png" width="49%" style="max-width: 450px; height: auto;" alt="Figure 5: Spatial Attention Heatmap">
</p>
<p align="center"><em>Left: Attention-MIL survival stratification on routine diagnostic slides (C-index ≈ 0.74). Right: Spatial attention projection localizing top attention weights to solid sanctuary histology.</em></p>

</div>

---

## Repository Structure

```
spatial-oncology-master/
├── .github/
│   └── workflows/
│       └── ci.yml                     # Multi-OS automated testing & verification
├── data/
│   ├── cohort_manifest_study4.csv     # N = 570 diagnostic WSI manifest
│   └── combined_filtered_massive_cohort.csv # N = 854 clinical survival cohort
├── figures/
│   ├── figure1_game_decomposition.png
│   ├── figure1_replicator_dynamics.png
│   ├── figure2_subcellular_architecture_624um.png
│   ├── figure3_massive_survival_validation_854pts.png
│   ├── figure4_he_vision_survival_prediction.png
│   └── figure5_spatial_attention_heatmap.png
├── manuscript/
│   ├── MANUSCRIPT_UNIFIED_MASTER_STUDY.md   # Complete peer-review manuscript
│   └── MANUSCRIPT_UNIFIED_MASTER_STUDY.pdf  # Publication-ready compiled PDF
├── notebooks/
│   └── InSituVision_Master_Pipeline.ipynb   # Interactive Google Colab pipeline
├── src/
│   ├── replicator_dynamics_game.py    # Pillar 1: Game theory solver
│   ├── biophysical_diffusion_model.py # Pillar 2: Reaction-diffusion solver
│   ├── clinical_survival_validation.py# Pillar 3: Log-rank & C-index engine
│   ├── deep_attention_mil.py          # Pillar 4: Gated Attention-MIL architecture
│   └── master_pipeline.py             # End-to-end unified orchestrator
├── tests/
│   └── test_master_pipeline.py        # Automated test suite (5/5 passing)
├── CITATION.cff                       # Citation metadata
├── LICENSE                            # GNU General Public License v3.0
├── requirements.txt                   # Dependency specifications
└── README.md                          # Documentation
```

---

## Quickstart & Execution

### 1. Local Installation

```bash
# Clone the repository
git clone https://github.com/Lukecele/spatial-oncology-master.git
cd spatial-oncology-master

# Install dependencies
pip install -r requirements.txt
```

### 2. Run the End-to-End Master Pipeline

```bash
python src/master_pipeline.py
```

Expected terminal output:
```text
==============================================================================
  UNIFIED SPATIAL ONCOLOGY MASTER PIPELINE EXECUTION
  From Subcellular Biophysics to Routine Clinical Diagnostic H&E AI
==============================================================================
[Pillar 1: Game Theory]     Nash Equilibrium x* = 0.613 (Coexistence Confirmed)
[Pillar 2: Visium HD]       Penetration Length λ = 627.6 µm (Sanctuary Threshold: 624 µm)
[Pillar 3: Clinical N=854]  High Decoupled PFS: 17.9 mo vs Low Decoupled: 38.1 mo | HR = 19.45 | p = 0.00e+00 | C-index = 0.939
[Pillar 4: H&E Vision AI]   WSI N = 570 | C-index = 0.739 | Slide Cost = 4.50 € (-99.91%)
==============================================================================
```

### 3. Run Automated Unit Tests

```bash
python -m unittest discover -s tests -v
```

### 4. Interactive Cloud Execution (Google Colab)
Launch the complete 4-pillar analysis in Google Colab with one click:  
[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Lukecele/spatial-oncology-master/blob/main/notebooks/InSituVision_Master_Pipeline.ipynb)

---

## Citation & Academic Attribution

If you use this codebase, mathematical model, clinical findings, or diagnostic AI pipeline in your research, please cite:

```bibtex
@article{celebrano2026spatialmaster,
  title={From Subcellular Spatial Biophysics to Routine Digital Pathology AI: Resolving Non-Genetic Osimertinib Resistance via Microvascular Decoupling},
  author={Celebrano, Luca},
  journal={GitHub Repository: Lukecele/spatial-oncology-master},
  year={2026},
  doi={10.5281/zenodo.23623653},
  url={https://github.com/Lukecele/spatial-oncology-master}
}
```

**Author Contact:**  
**Luca Celebrano**  
Department of Clinical Medicine and Surgery, School of Medicine, University of Naples Federico II, Naples, Italy  
Telethon Institute of Genetics and Medicine (TIGEM), Pozzuoli, Italy  
Email: `l.celebrano@studenti.unina.it` | ORCID: [0009-0007-0931-3555](https://orcid.org/0009-0007-0931-3555)

---

## License

This project is open-source software licensed under the **GNU General Public License v3.0 (GPL-3.0)**. See the [LICENSE](LICENSE) file for complete copyleft terms and conditions.
