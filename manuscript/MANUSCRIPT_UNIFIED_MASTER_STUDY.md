# Conditions for Pharmacological Sanctuaries in Perivascular Tumour Cords: A Reaction–Diffusion Analysis, a Regrowth-Kinetics Identity, and a Descriptive Survival Analysis of Histologic Subtype in Resected Lung Adenocarcinoma

**Luca Celebrano**  
*Department of Clinical Medicine and Surgery, School of Medicine and Surgery, University of Naples Federico II, Naples, Italy*  
*Correspondence: luca.celebrano5@gmail.com / l.celebrano@studenti.unina.it | ORCID: [0009-0006-2584-7546](https://orcid.org/0009-0006-2584-7546)*  
*Permanent Zenodo Master DOI:* [10.5281/zenodo.23623653](https://doi.org/10.5281/zenodo.23623653)  
*GitHub Code and Data Repository:* [https://github.com/Lukecele/spatial-oncology-master](https://github.com/Lukecele/spatial-oncology-master)

---

## Abstract

**Background:** After 3 years of adjuvant osimertinib in resected *EGFR*-mutant non-small cell lung cancer (NSCLC), most disease-free survival (DFS) or circulating tumour DNA (ctDNA) detection events occur after treatment discontinuation. In the exploratory molecular residual disease (MRD) analysis of the Phase III ADAURA trial (post hoc, $n = 220$), 28 of 112 patients in the osimertinib arm (25.0%) experienced a DFS or MRD event; 19 of these 28 events (67.9% $\approx 68\%$) occurred after treatment cessation, with 11 of those 19 (57.9% $\approx 58\%$) occurring within 12 months of stopping (Herbst et al., *Nature Medicine* 2025). A potential, untested biophysical hypothesis is that residual neoplastic cells reside in microenvironmental niches where drug exposure falls below therapeutic inhibitory thresholds while tissue oxygenation remains sufficient to support cellular viability.

**Methods:** We established a three-part analytical framework: (1) We solved the steady-state reaction–diffusion equation for a therapeutic agent in cylindrical geometry using modified Bessel functions for two distinct physical configurations: a perivascular cord centered on a functional microvessel (outward diffusion, $K_0$) and an avascular solid nest supplied at its perimeter (inward diffusion, $I_0$), mapping the exact parameter domains in which a sub-inhibitory but viable sanctuary shell can exist. (2) We analyzed a two-phenotype replicator game dynamic coupled to an exponential regrowth formulation, deriving an analytical identity for the model's ctDNA-to-imaging diagnostic lead time. (3) We conducted a descriptive survival analysis of the authenticated Memorial Sloan Kettering Cancer Center cohort (MSKCC 2020; $N = 604$ resected lung adenocarcinomas, 98 recurrence events) across IASLC architectural growth patterns using reverse Kaplan–Meier follow-up estimation, multivariable Cox proportional hazards modeling with Holm multiplicity adjustment, and bootstrap resamplings.

**Results:** With baseline parameters for small-molecule osimertinib (hindered tissue diffusivity $D = 1.3 \times 10^{-6}\,\text{cm}^2/\text{s}$ [Swabb et al. 1974], volumetric cellular uptake $k = 3.3 \times 10^{-4}\,\text{s}^{-1}$, characteristic screening length $\lambda = 627.6\,\mu\text{m}$, capillary radius $r_0 = 10.0\,\mu\text{m}$, inhibitory threshold $\theta = C(IC_{50})/C_0 = 0.35$), the outward perivascular solution yields an exact critical penetration depth of $r_{\text{crit}} = 165.9\,\mu\text{m}$. Assuming a physiological oxygen viability limit of $r_{\text{ox}} = 200.0\,\mu\text{m}$ (Thomlinson & Gray 1955), an outward viable pharmacological sanctuary shell emerges ($\Delta r = 34.1\,\mu\text{m}$). This shell is fragile: it collapses ($\Delta r \to 0$) under jointly unfavorable parameter variations (e.g., $D +30\%$ and $k -30\%$ yielding $r_{\text{crit}} = 201.4\,\mu\text{m}$), when $\theta < 0.310$, or when the viable cord radius contracts to $r_{\text{ox}} \le 165.9\,\mu\text{m}$. Crucially, in the inward geometry describing avascular nests, osimertinib penetrates almost completely to the center (reaching $97.5\%$ of boundary concentration at $R = 200\,\mu\text{m}$ and $90.6\%$ at $R = 400\,\mu\text{m}$); thus, avascular solid nests do not form a sanctuary for osimertinib under baseline parameters. In population dynamic modeling, the ctDNA-to-imaging diagnostic lead time satisfies the exact mathematical identity $\Delta t = T_d \cdot \log_2(100) = 6.6439 \cdot T_d$ ($5.46\,\text{months}$ for $T_d = 25\,\text{days}$), which is an invariant algebraic consequence of assumed exponential doubling between assay limits, independent of initial residual burden $N_0$. In the MSKCC cohort, reverse Kaplan–Meier median follow-up was $23.8\,\text{months}$ (IQR: 11.6–29.1). Lepidic architecture exhibited a significantly lower recurrence hazard than all other subtypes (stage-adjusted $HR = 0.25, 95\%\text{ CI: } 0.09\text{--}0.70, p = 0.0084$ relative to acinar); solid and papillary patterns showed non-significant hazard trends relative to acinar ($HR = 1.48, p = 0.133$ and $HR = 1.51, p = 0.215$), while micropapillary histology did not demonstrate increased hazard ($HR = 0.74, p = 0.468$). An exploratory boundary comparison of solid/micropapillary versus lepidic architecture ($N = 193, 33\text{ events}$) yielded an adjusted $HR = 6.75$ ($95\%\text{ CI: } 2.28\text{--}19.97$), but bootstrap resampling showed that the upper confidence limit is unbounded due to the small number of events in the lepidic reference group ($4\text{ events}$). In the full cohort ($N = 604$), the stage-adjusted hazard ratio for solid/micropapillary architecture versus all other subtypes was $1.32$ ($95\%\text{ CI: } 0.83\text{--}2.08, p = 0.239$).

**Conclusions:** A viable, sub-inhibitory pharmacological sanctuary shell around a perivascular microvessel is theoretically possible within a bounded parameter domain, but it is not an inevitable consequence of neoplastic tissue compaction, and the model demonstrates no sanctuary for osimertinib in avascular nests under physiological parameters. In clinical survival data, lepidic histology is confirmed as a strongly protective prognostic marker; however, surgical cohort data cannot validate the sanctuary hypothesis, post-cessation osimertinib kinetics, or ctDNA surveillance schedules. We formulate explicit, testable experimental predictions and refrain from unwarranted clinical recommendations.

---

## 1. Introduction

In the Phase III ADAURA trial, 3 years of daily adjuvant osimertinib significantly improved disease-free survival (DFS) and overall survival (OS) in patients with completely resected stage IB–IIIA *EGFR*-mutated non-small cell lung cancer (NSCLC) [1, 2]. However, longitudinal surveillance demonstrated that treatment discontinuation at month 36 is accompanied by an emergence of disease recurrences. In an exploratory, post hoc molecular residual disease (MRD) analysis of ADAURA ($n = 220$) [12], circulating tumour DNA (ctDNA) positivity preceded radiographic DFS relapse by a median of 4.7 months (95% CI: 2.2–5.6 months) across trial arms. Furthermore, within the osimertinib arm, 28 of 112 patients (25.0%) had a DFS or MRD event; 19 of those 28 events (67.9% $\approx 68\%$) occurred *after* treatment cessation, and 11 of those 19 post-cessation events (57.9% $\approx 58\%$) occurred within the first 12 months following cessation [12].

Post-cessation relapses in the adjuvant setting frequently retain sensitivity to EGFR tyrosine kinase inhibitor (TKI) rechallenge and lack secondary on-target genomic mutations (such as *EGFR* C797S) or off-target bypass pathway activations (such as *MET* amplification). These observations have motivated hypotheses regarding non-mutational mechanisms of persistence, including drug-tolerant persister (DTP) states [6, 7], microenvironmental sanctuary niches, and tumour dormancy. Spatial drug-delivery limitations have been proposed as one candidate mechanism: if neoplastic cellular compaction creates physical transport barriers, drug concentrations in poorly perfused regions might drop below the threshold required to eradicate sensitive clones, allowing cells to survive in an uninhibited or sub-therapeutic state while remaining viable.

The objective of this study is three-fold:
1. To formulate the steady-state reaction–diffusion transport problem in cylindrical coordinates, comparing outward diffusion (perivascular cords) against inward diffusion (avascular nests), and to establish the exact biophysical boundaries under which a viable pharmacological sanctuary can exist.
2. To examine evolutionary game dynamics between sensitive and persister phenotypes and evaluate post-cessation regrowth kinetics, proving an exact analytical identity for the model's ctDNA-to-imaging diagnostic lead time.
3. To conduct a rigorous descriptive survival analysis of architectural growth patterns in the public Memorial Sloan Kettering Cancer Center (MSKCC 2020) surgical cohort ($N = 604$), assessing what architectural features are empirically supported and identifying methodological limitations.

---

## 2. Reaction–Diffusion Transport in Neoplastic Architectures

### 2.1 Model and Analytical Solutions
Neoplastic drug delivery in cylindrical symmetry is governed by steady-state Fickian diffusion coupled to first-order cellular uptake and metabolic clearance:

$$D \nabla^2 C - k C = 0 \implies \frac{1}{r} \frac{d}{dr} \left( r \frac{dC}{dr} \right) - \frac{k}{D} C = 0$$

where $D$ is the effective tissue diffusion coefficient ($\text{cm}^2/\text{s}$), $k$ is the volumetric first-order cellular uptake and clearance rate ($\text{s}^{-1}$), and $\lambda = \sqrt{D/k}$ defines the characteristic screening penetration length. The general solution of this homogeneous Helmholtz equation is:

$$C(r) = A \cdot I_0(r/\lambda) + B \cdot K_0(r/\lambda)$$

where $I_0$ and $K_0$ are modified Bessel functions of the first and second kind of order zero, respectively. The physical boundary conditions depend fundamentally on tissue geometry:

1. **Perivascular Cord (Outward Diffusion Geometry):**  
   A central functional microvessel of radius $r_0$ maintains an unbound drug concentration $C(r_0) = C_0$. As $r \to \infty$, physical bounding demands $C(r) \to 0$. Because $I_0(r/\lambda)$ diverges exponentially ($I_0(x) \sim e^x / \sqrt{2\pi x}$), its coefficient must vanish ($A = 0$). The exact outward solution is:
   $$C_{\text{out}}(r) = C_0 \frac{K_0(r/\lambda)}{K_0(r_0/\lambda)}$$

2. **Avascular Nest (Inward Diffusion Geometry):**  
   An avascular cellular aggregate of radius $R$ is surrounded by vascularized stroma supplying an external drug concentration $C(R) = C_0$. At the nest center ($r = 0$), the drug concentration must remain finite. Because $K_0(r/\lambda)$ diverges logarithmically as $r \to 0$ ($K_0(x) \sim -\ln(x)$), its coefficient must vanish ($B = 0$). The exact inward solution is:
   $$C_{\text{in}}(r) = C_0 \frac{I_0(r/\lambda)}{I_0(R/\lambda)}$$

**The two configurations are not mathematically or physically equivalent.** Outward diffusion from a central capillary exhibits geometric divergence (dilution across expanding cross-sections $2\pi r$), causing a steep radial concentration drop. Conversely, inward diffusion toward a solid nest center exhibits geometric convergence (flux concentrating into contracting cross-sections $2\pi r$), which counteracts metabolic consumption. 

### 2.2 Physiological Parameters & Viability Limits
Defining the normalized therapeutic threshold ratio as $\theta = C(IC_{50}) / C_0$, a pharmacological sanctuary shell exists if the drug concentration drops below $\theta$ while remaining within the tissue oxygen viability limit $r_{\text{ox}}$:

$$\Delta r = \max(0, \min(R, r_{\text{ox}}) - r_{\text{crit}})$$

For the outward perivascular geometry, the critical penetration radius $r_{\text{crit}}$ is the unique root of:

$$\frac{K_0(r_{\text{crit}}/\lambda)}{K_0(r_0/\lambda)} = \theta$$

Crucially, $r_{\text{crit}}$ depends exclusively on $\lambda$, $r_0$, and $\theta$; the parameters $D$ and $k$ are not separately identifiable from spatial shell geometry alone.

Baseline parameters were calibrated as order-of-magnitude physiological estimates:
- **Small-molecule TKI (Osimertinib, 500 Da):** Effective tissue diffusivity $D = 1.3 \times 10^{-6}\,\text{cm}^2/\text{s}$ and volumetric uptake rate $k = 3.3 \times 10^{-4}\,\text{s}^{-1}$, yielding $\lambda = 627.6\,\mu\text{m}$. We emphasize that $D = 1.3 \times 10^{-6}\,\text{cm}^2/\text{s}$ reflects hindered tissue diffusion through dense extracellular matrix (Swabb, Wei & Gullino, *Cancer Res* 1974 [13]; Nugent & Jain, *Cancer Res* 1984 [14]), rather than free aqueous diffusion ($D_{\text{aq}} \approx 4\text{--}6 \times 10^{-6}\,\text{cm}^2/\text{s}$). Microvascular capillary radius is set to $r_0 = 10.0\,\mu\text{m}$.
- **Macromolecular Therapeutic (ADC / Bispecific Antibody, 150 kDa):** Hindered diffusivity $D = 1.5 \times 10^{-8}\,\text{cm}^2/\text{s}$ and endocytic clearance $k = 2.5 \times 10^{-4}\,\text{s}^{-1}$, yielding $\lambda = 77.5\,\mu\text{m}$.
- **Inhibitory Threshold Ratio ($\theta$):** Derived from steady-state clinical pharmacokinetics. For osimertinib 80 mg daily, mean steady-state total plasma concentration is $\sim 500\,\text{nM}$ ($250\text{--}300\,\text{ng/mL}$). With human plasma protein binding of 94.7% ($f_u = 5.3\%$), free vascular concentration is $C_{0,\text{free}} \approx 25.0\,\text{nM}$. In human *EGFR*-mutant cell lines, $IC_{50} \approx 8.8\,\text{nM}$ (Cross et al., *Cancer Discovery* 2014 [5]), establishing:
  $$\theta = \frac{IC_{50}}{C_{0,\text{free}}} = \frac{8.8\,\text{nM}}{25.0\,\text{nM}} = 0.352 \approx 0.35$$
- **Oxygen Viability Boundary ($r_{\text{ox}}$):** In seminal histological measurements of bronchial carcinoma cords, Thomlinson & Gray (*Br J Cancer* 1955 [4]) established that central necrosis was absent in cords under $\sim 160\,\mu\text{m}$ in radius and invariably present in cords over $200\,\mu\text{m}$, with intact viable tumour sheaths typically spanning $100\text{--}180\,\mu\text{m}$ (Worth et al., *Br J Cancer* 2022 [14]; Hirst et al., *Br J Cancer* 1991 [15]). Setting $r_{\text{ox}} = 200.0\,\mu\text{m}$ represents the upper physiological ceiling of cellular viability.

### 2.3 Transport Results & Sensitivity Mapping

In the outward perivascular configuration, baseline parameters yield:
- Osimertinib: $r_{\text{crit}} = 165.90\,\mu\text{m}$. For $r_{\text{ox}} = 200.0\,\mu\text{m}$, an outward viable pharmacological sanctuary shell emerges with width $\Delta r = 34.10\,\mu\text{m}$ (**Figure 1A**).
- Macromolecule (ADC): $r_{\text{crit}} = 47.43\,\mu\text{m}$, leaving a massive uninhibited viable shell of $\Delta r = 152.57\,\mu\text{m}$ ($>94\%$ of the cord uninhibited).

(<strong>Figure 1</strong>).

![Figure 1: Biophysical Reaction-Diffusion Transport Physics & Sensitivity Mapping](Figure1_transport_sensitivity.png)
*Figure 1. Biophysical Reaction-Diffusion Transport Physics & Sensitivity Mapping. (A) Normalized radial drug concentration profiles for small-molecule osimertinib and macromolecular ADC in outward perivascular geometry. Dotted lines indicate critical depth rcrit where drug drops below IC50 (theta = 0.35). (B) Concentration profiles in inward avascular nest geometry, demonstrating that osimertinib reaches >=90% at the center for nests up to 400 um. (C) Viable sanctuary area percentage. (D) Parameter sensitivity phase map across D, k, theta, and rox.*

#### Table 1. Sensitivity of the Osimertinib Outward Perivascular Sanctuary Shell ($\theta = 0.35, r_0 = 10\,\mu\text{m}, r_{\text{ox}} = 200\,\mu\text{m}$)

| Scenario | Penetration Length $\lambda$ ($\mu\text{m}$) | Critical Radius $r_{\text{crit}}$ ($\mu\text{m}$) | Sanctuary Shell Width $\Delta r$ ($\mu\text{m}$) | Operational State |
| :--- | :---: | :---: | :---: | :--- |
| **Baseline Calibration** | **627.6** | **165.9** | **34.1** | **Viable Sanctuary Shell** |
| $D +30\%$ | 715.6 | 180.1 | 19.9 | Shell narrows |
| $D -30\%$ | 525.1 | 148.5 | 51.5 | Shell widens |
| $k +30\%$ | 550.5 | 152.9 | 47.1 | Shell widens |
| $k -30\%$ | 750.2 | 185.5 | 14.5 | Shell narrows |
| **$D +30\%$ and $k -30\%$ Concurrently** | **855.3** | **201.4** | **0.0** | **Shell Collapses Completely** |
| $D -30\%$ and $k +30\%$ Concurrently | 460.6 | 136.9 | 63.1 | Shell expands substantially |
| Capillary radius $r_0 = 5\,\mu\text{m}$ | 627.6 | 128.3 | 71.7 | Shell expands |
| Capillary radius $r_0 = 20\,\mu\text{m}$ | 627.6 | 216.2 | 0.0 | Shell Collapses Completely |

As demonstrated in **Table 1**, the perivascular shell persists under individual $\pm 30\%$ parameter variations, but collapses completely under the concurrent combination of $+30\% D$ and $-30\% k$ ($r_{\text{crit}} = 201.4\,\mu\text{m} > r_{\text{ox}}$). Furthermore, the existence of the shell requires $\theta > 0.310$ at $r_{\text{ox}} = 200\,\mu\text{m}$; $\theta > 0.332$ at $r_{\text{ox}} = 180\,\mu\text{m}$; and collapses entirely if $r_{\text{ox}} \le 165.9\,\mu\text{m}$ (**Figure 1D**). Monte Carlo sensitivity sampling across log-normal distributions of $D, k, \theta$ and uniform $r_{\text{ox}} \in [100, 200]\,\mu\text{m}$ indicates that a viable shell exists in 46% of simulated realizations (and 63% with $r_{\text{ox}}$ fixed at $200\,\mu\text{m}$).

#### Table 2. Inward Avascular Nest Geometry: Core Drug Concentration Relative to Boundary ($C(0) / C_0$)

| Nest Radius $R$ ($\mu\text{m}$) | Osimertinib ($500\,\text{Da}, \lambda = 627.6\,\mu\text{m}$) | Macromolecule / ADC ($150\,\text{kDa}, \lambda = 77.5\,\mu\text{m}$) |
| :---: | :---: | :---: |
| $100\,\mu\text{m}$ | **0.994** (99.4%) | 0.684 (68.4%) |
| $200\,\mu\text{m}$ | **0.975** (97.5%) | **0.285** (28.5% — Sanctuary Core) |
| $300\,\mu\text{m}$ | **0.945** (94.5%) | 0.099 (9.9% — Deep Sanctuary) |
| $400\,\mu\text{m}$ | **0.906** (90.6%) | 0.032 (3.2% — Deep Sanctuary) |

**Crucial Finding for Avascular Nests:** As shown in **Table 2**, in an avascular solid nest supplied from the periphery, osimertinib concentrations at the center remain between $90.6\%$ and $99.4\%$ of the stromal boundary concentration for nests up to $400\,\mu\text{m}$ in radius. Because $C(0) \gg 0.35 C_0$, **no pharmacological sanctuary can exist for osimertinib inside avascular nests under physiological parameters**. Generating a sanctuary ($C(0) \le 0.35$) at the center of a $200\,\mu\text{m}$ avascular nest would require $\lambda \le 86.5\,\mu\text{m}$, representing an approximately 53-fold higher cellular uptake rate ($k \approx 0.017\,\text{s}^{-1}$) at fixed $D$. Therefore, our reaction-diffusion model proves that depth-limited small-molecule sanctuaries are strictly restricted to outward perivascular cord geometries, while inward avascular nests create sanctuaries exclusively for macromolecules.

---

## 3. Evolutionary Game Dynamics & Post-Cessation Regrowth Kinetics

### 3.1 Replicator Dynamics Model
To examine subclonal coexistence within the perivascular sanctuary shell, we formalized a symmetric $2 \times 2$ evolutionary game between drug-sensitive proliferating cells ($S$, fraction $x$) and drug-tolerant persister cells ($P$, fraction $1-x$). Under non-zero-sum frequency-dependent interactions (Kaznatcheev et al., *Nat Ecol Evol* 2019 [15]):

$$A = \begin{pmatrix} a_{11} & a_{12} \\ a_{21} & a_{22} \end{pmatrix} = \begin{pmatrix} -0.736 & 2.866 \\ 2.188 & -1.772 \end{pmatrix}$$

The continuous replicator equation is:

$$\frac{dx}{dt} = x(1 - x)[f_S(x) - f_P(x)] = x(1 - x)[\alpha x + \beta]$$

where $\alpha = a_{11} - a_{12} - a_{21} + a_{22} = -7.5620$ and $\beta = a_{12} - a_{22} = 4.6380$. Setting $dx/dt = 0$ yields the internal polymorphic Nash equilibrium:

$$x^* = -\frac{\beta}{\alpha} = \frac{4.6380}{7.5620} = \mathbf{0.6133} \quad (61.33\%\text{ Sensitive})$$

The Jacobian eigenvalue at $x^*$ is $\lambda_{\text{eig}} = x^*(1 - x^*)\alpha = -1.7934 < 0$, confirming that $x^* = 0.6133$ is an asymptotically stable Lyapunov attractor (**Figure 2A**). Because $a_{12} > a_{22}$ ($2.866 > -1.772$) and $a_{21} > a_{11}$ ($2.188 > -0.736$), any payoff matrix satisfying this sign pattern yields stable coexistence (Hawk–Dove / Snowdrift game). We explicitly clarify that these payoffs represent free parameters chosen to formalize this qualitative structure, rather than values fitted to empirical measurements.

(<strong>Figure 2</strong>).

![Figure 2: Evolutionary Replicator Dynamics & Post-Cessation Regrowth Identity](Figure2_regrowth_identity.png)
*Figure 2. Evolutionary Replicator Dynamics & Post-Cessation Regrowth Identity. (A) Phase portrait of the replicator equation confirming global convergence to x* = 0.6133 (eigenvalue -1.793). (B) Clonal regrowth trajectories from N0 = 2000 cells crossing ctDNA (10^7) and CT (10^9) thresholds. (C) Graphical proof of the analytical lead-time identity Delta t = Td * log2(100), showing that lead time is strictly independent of initial sanctuary load N0.*

### 3.2 Post-Cessation Regrowth & The Lead-Time Identity
Following adjuvant treatment discontinuation at month 36, drug clears from the sanctuary niche. Total cell burden expands from an initial residual micro-nest $N_0$ with initial sensitive fraction $x^* = 0.6133$:

$$N(t) = N_0 \left[ x^* \cdot 2^{t / T_s} + (1 - x^*) \cdot 2^{t / T_p} \right]$$

where $T_s = 25.0\,\text{days}$ (sensitive in vivo doubling) and $T_p = 120.0\,\text{days}$ (persister in vivo doubling). Clinical detection thresholds are defined as:
- Molecular ctDNA detection limit: $N_{\text{MRD}} = 10^7\,\text{cells}$ ($\sim 0.01\%$ variant allele fraction).
- Radiographic CT imaging limit: $N_{\text{CT}} = 10^9\,\text{cells}$ ($\sim 1\,\text{cm}^3$ nodule).

For $N_0 = 2,000\,\text{cells}$, ctDNA threshold crossing occurs at $t = 324.8\,\text{days}$ (**10.67 months**), and radiographic detection occurs at $t = 490.9\,\text{days}$ (**16.13 months**), yielding a diagnostic lead time window of $\Delta t = 166.1\,\text{days} = \mathbf{5.46\,\text{months}}$ (**Figure 2B**).

#### The Analytical Lead-Time Identity
Crucially, this lead-time window follows directly from an exact algebraic identity:

$$\Delta t = T_s \cdot \log_2\left(\frac{N_{\text{CT}}}{N_{\text{MRD}}}\right) = T_s \cdot \log_2(100) = 6.6439 \cdot T_s$$

Because sensitive cells rapidly outgrow slow-cycling persisters ($>99.9\%$ sensitive cells at $N = 10^7$), the duration required to expand 100-fold is strictly $6.6439$ doublings of the sensitive clone (**Figure 2C**). **The lead time $\Delta t = 5.46\,\text{months}$ is an invariant property of the assumed doubling time $T_s$ and assay limits; it is entirely independent of $N_0$, independent of persister kinetics, and independent of any spatial sanctuary mechanism.**

#### Table 3. Dependence of Model Regrowth Outputs on Assumed Biological Inputs

| Assumed Sensitive Doubling Time $T_s$ (days) | Resulting Lead Time $\Delta t$ (months) | Residual Inoculum $N_0$ (cells) | ctDNA Recurrence Month | Radiographic CT Relapse Month |
| :---: | :---: | :---: | :---: | :---: |
| 10 days | 2.2 months | $N_0 = 1$ | 19.7 months | 25.1 months |
| 15 days | 3.3 months | $N_0 = 100$ | 14.2 months | 19.7 months |
| 20 days | 4.4 months | $N_0 = 1,000$ | 11.5 months | 17.0 months |
| **21.5 days** | **4.7 months (ADAURA Median)** | **$N_0 = 2,000$** | **10.7 months** | **16.1 months** |
| **25.0 days** | **5.5 months** | $N_0 = 10,000$ | 8.8 months | 14.2 months |
| 30.0 days | 6.5 months | $N_0 = 100,000$ | 6.0 months | 11.5 months |

As shown in **Table 3**, the reported median lead time in ADAURA (4.7 months; 95% CI: 2.2–5.6 months [12]) corresponds mathematically to $T_s \approx 21.5\,\text{days}$ (range: 10–26 days). Therefore, agreement between our simulation and ADAURA reflects a biologically plausible choice of $T_s$, rather than an independent validation of spatial sanctuary physics. What residual sanctuary burden $N_0$ actually governs is the **calendar timing of initial relapse**: $N_0 = 200\text{--}2000$ cells reproduces the empirical onset clustering documented in ADAURA (peaking within 12–18 months post-cessation).

---

## 4. Descriptive Survival Analysis of the MSKCC Cohort ($N = 604$)

### 4.1 Cohort Characteristics & Methods
We analyzed the authenticated Memorial Sloan Kettering Cancer Center cohort (Caso et al., *J Thorac Oncol* 2020 [3]), comprising 604 consecutive surgically resected lung adenocarcinoma patients with predominant IASLC histologic subtypes, pathologic stage (I/II/III), *EGFR* mutational status, recurrence-free survival (RFS) duration, and recurrence events (98 recurrences, 16.23%). All patients were chemo/TKI-naive at resection; none received adjuvant osimertinib.

Survival analyses were performed using Kaplan–Meier estimators and multivariable Cox proportional hazards models with Efron tie handling and Wald confidence intervals. True median follow-up was computed using the **reverse Kaplan–Meier method** (Schemper & Smith, *Eff Clin Pract* 1996): **median follow-up was 23.8 months** (IQR: 11.6–29.1 months; version 1 reported 21.4 months based on raw RFS duration). A total of 77 patients remained at risk at 36 months, and only 3 beyond 60 months, with only 7 recurrence events occurring after 36 months. Consequently, this surgical cohort cannot inform recurrence kinetics occurring between months 36 and 60 post-cessation.

(<strong>Figure 3</strong>):

![Figure 3: MSKCC Clinical Cohort Validation & Survival Modeling](Figure3_cohort_cox.png)
*Figure 3. MSKCC Clinical Cohort Validation & Survival Modeling (N = 604). (A) Kaplan-Meier recurrence-free survival curves for the 5 IASLC histologic subtypes. (B) Observed crude recurrence rates across subtypes. (C) Forest plot of pairwise subtype contrasts and multivariable Cox models with Holm multiplicity adjustments. (D) Evaluation of the EGFR-mutant subset showing severe stage confounding.*

#### Table 4. Recurrence Events and Stage Distribution by Predominant Histologic Subtype

| Predominant Subtype | Total Patients ($N$) | Recurrence Events ($E$) | Crude Recurrence Rate (%) | Pathologic Stage I (%) |
| :--- | :---: | :---: | :---: | :---: |
| **Lepidic** | 88 | 4 | **4.55%** | **88.6%** |
| **Acinar** | 368 | 54 | **14.67%** | 75.5% |
| **Micropapillary** | 37 | 8 | **21.62%** | 51.4% |
| **Papillary** | 43 | 11 | **25.58%** | 76.7% |
| **Solid** | 68 | 21 | **30.88%** | 57.4% |
| **Total Cohort** | **604** | **98** | **16.23%** | **74.0%** |

Histologic subtype is significantly associated with recurrence after adjusting for pathologic stage (likelihood ratio test: $\chi^2 = 18.33$, 4 degrees of freedom, $p = 0.0011$). However, pathologic stage is unevenly distributed (**Table 4**): 88.6% of lepidic tumours presented at Stage I, compared to only 55.2% of solid and micropapillary tumours.

### 4.2 Multivariable Survival Modeling Results

#### Table 5. Stage-Adjusted Multivariable Cox Proportional Hazards Models

| Model Specification | Hazard Ratio (HR) | 95% Confidence Interval | $p$-value | Holm-Adjusted $p$-value |
| :--- | :---: | :---: | :---: | :---: |
| **Pairwise Contrasts vs Lepidic Reference:** | | | | |
| - Solid vs Lepidic | **7.82** | 2.60 – 23.55 | $0.0003$ | **0.003** |
| - Papillary vs Lepidic | **6.07** | 1.91 – 19.32 | $0.0020$ | **0.020** |
| - Acinar vs Lepidic | **3.94** | 1.41 – 10.97 | $0.0090$ | **0.070** |
| - Micropapillary vs Lepidic | **3.85** | 0.97 – 15.34 | $0.0560$ | 0.390 |
| **Composite Architectural Contrasts:** | | | | |
| - Solid/Micro vs Lepidic ($N=193$, 33 events) | **6.75** | 2.28 – 19.97 | $<0.001$ | — |
| - ...extended model (+ stage levels, sex, smoking, *EGFR*) | **7.15** | 2.40 – 21.33 | $<0.001$ | — |
| - Solid/Micro vs All Others ($N=604$, 98 events) | **1.32** | 0.83 – 2.08 | $0.239$ | — |
| - ...plus *EGFR* mutation status ($N=604$) | **1.30** | 0.82 – 2.07 | $0.258$ | — |
| **5-Subtype Model (Acinar Reference):** | | | | |
| - Lepidic | **0.25** | 0.09 – 0.70 | **0.0084** | — |
| - Solid | **1.48** | 0.89 – 2.49 | $0.1330$ | — |
| - Papillary | **1.51** | 0.79 – 2.92 | $0.2150$ | — |
| - Micropapillary | **0.74** | 0.33 – 1.66 | $0.4680$ | — |
| - Pathologic Stage II/III vs I | **4.16** | 2.73 – 6.34 | $<0.001$ | — |

As detailed in **Table 5**, the primary empirical finding from this cohort is that **lepidic histology is associated with a significantly lower recurrence hazard than any other subtype** ($HR = 0.25, p = 0.0084$ vs acinar). Conversely, the data do not demonstrate an escalating risk gradient tracking compaction: solid ($HR = 1.48, p = 0.133$) and papillary ($HR = 1.51, p = 0.215$) do not differ significantly from acinar, and micropapillary histology exhibits an adjusted $HR = 0.74$ ($p = 0.468$). Grouping solid and micropapillary patterns as a single high-sanctuary class is therefore not empirically supported by intra-cohort risk separation.

### 4.3 Robustness & Limitations of the Headline Contrast
The headline contrast comparing solid/micropapillary against lepidic architecture ($HR = 6.75, 95\%\text{ CI: } 2.28\text{--}19.97$) rests upon only 4 recurrence events in the lepidic reference group. Bootstrap resampling (500 iterations) yielded a median HR of 6.9 and a 2.5th percentile of 2.6, but the upper percentile was unbounded because numerous resamples contained zero lepidic events; the Wald upper limit of 19.97 is therefore statistically unstable. Truncating follow-up at 36 months yielded $HR = 6.11$ (95% CI: 2.04–18.31; 30 events); at 24 months, $HR = 23.8$ (95% CI: 3.14–180.3; 26 events); and at 18 months, the estimate does not exist because zero lepidic patients experienced recurrence by that time. Furthermore, this contrast was not registered in a prospective protocol and excludes 411 of 604 patients; it must be interpreted strictly as an exploratory boundary comparison.

In the confirmed *EGFR*-mutant subpopulation ($N = 182, 26\text{ events}$), recurrence occurred in 4 of 18 solid/micropapillary patients (22.2%) versus 22 of 164 other patients (13.4%; univariable log-rank $p = 0.112$). Restricting to solid/micropapillary versus lepidic ($N = 56, 6\text{ events}$) yielded an apparent univariable log-rank $p = 0.021$. However, this comparison is entirely confounded by stage: 37 of 38 lepidic tumours were Stage I (97.4%), whereas 11 of 18 solid/micropapillary tumours were Stage II/III (61.1%). After adjusting for pathologic stage, the hazard ratio collapsed to $HR = 1.10$ ($95\%\text{ CI: } 0.06\text{--}19.97, p = 0.947$). The *EGFR*-mutant subgroup is underpowered and uninformative.

---

## 5. Discussion

### 5.1 What is Supported by the Data
1. **Perivascular Cords Support Pharmacological Sanctuaries:** Within exact steady-state reaction–diffusion physics, a sub-inhibitory but oxygenated perivascular shell ($r_{\text{crit}} < r < r_{\text{ox}}$) can exist around a central microvessel. For baseline osimertinib parameters, a viable shell of width $\Delta r = 34.1\,\mu\text{m}$ exists, but collapses under moderate, jointly unfavorable parameter variations ($D +30\%$ and $k -30\%$) or in cords with oxygen limits $r_{\text{ox}} \le 165.9\,\mu\text{m}$.
2. **Avascular Nests Do Not Shelter Small-Molecule TKIs:** In inward cylindrical diffusion, osimertinib penetrates almost completely to the center ($90.6\%\text{--}99.4\%$), proving that avascular solid nests do not establish pharmacological sanctuaries for small molecules under physiological uptake rates.
3. **Lead Time is an Invariant Doubling Identity:** The ctDNA-to-imaging lead time is analytically determined by clonal doubling kinetics ($\Delta t = 6.6439 \cdot T_s$), rather than spatial sanctuary parameters.
4. **Lepidic Architecture Confers Strong Baseline Protection:** In resected lung adenocarcinoma, lepidic growth is an independent marker of favorable prognosis ($HR = 0.25, p = 0.0084$).

### 5.2 What is Not Supported & Clinical Cautions
The empirical data do not prove that tissue compaction causes non-mutational persistence, that biophysical sanctuaries explain post-ADAURA relapses, that histologic subtypes predict post-osimertinib recurrence, or that ctDNA monitoring improves patient survival. Intercepting relapses in the molecular-only phase incurs a risk of **lead-time bias**: an apparent prolongation of time-to-recurrence does not automatically confer an overall survival benefit unless prospective randomized trials demonstrate that preemptive treatment alteration alters clinical trajectory. The $\ge 20\%$ high-grade pattern cutoff defined by the IASLC Pathology Committee (Moreira et al., *J Thorac Oncol* 2020 [16]) is validated for post-surgical prognosis, not for ctDNA surveillance intervals. We therefore provide no clinical directives.

### 5.3 Testable Experimental Predictions
1. **Matrix-Assisted Laser Desorption/Ionization (MALDI) Mass-Spectrometry Imaging:** High-resolution spatial drug mapping in resected lung tumours from osimertinib-treated patients should reveal steep radial concentration gradients around perivascular cords with $r_{\text{crit}} \approx 166\,\mu\text{m}$, contrasted with flat, unhindered concentration profiles throughout avascular nests.
2. **Spatial Pharmacodynamic Gradients:** Neoplastic cells residing at distances $r > 166\,\mu\text{m}$ from functional capillaries should demonstrate persistent phospho-EGFR signaling despite therapeutic plasma concentrations.
3. **Prospective Architectural Tracking in Post-TKI Cohorts:** In prospective clinical trials with quantified architectural components, post-cessation recurrence rates should correlate specifically with the fraction of perivascular cord architecture.

---

## 6. Limitations

All reaction–diffusion parameters are literature-calibrated estimates rather than direct measurements in lung adenocarcinoma tissue. The transport model assumes steady-state conditions, idealized cylindrical geometry, non-saturable first-order uptake, and absence of convection or non-specific tissue binding. The game-theoretic formulation assumes constant payoffs and exponential clonal growth. The MSKCC cohort is observational, surgical, TKI-naive, and contains only 4 recurrence events in the lepidic reference group.

---

## Declarations

**Funding:** This research received no external grant from funding agencies in the public, commercial, or not-for-profit sectors.

**Competing Interests:** The author declares no competing financial or non-financial interests.

**Ethics Declaration:** This study analyzes publicly available, de-identified secondary clinical data from the Memorial Sloan Kettering Cancer Center cohort (cBioPortal; Caso et al., *J Thorac Oncol* 2020). Secondary analysis of anonymized public datasets requires no institutional review board approval.

**Data and Code Availability:** The authenticated clinical dataset (`06_DATI_CLINICI_MSKCC_604_PAZIENTI.csv`), analytical Bessel solvers, replicator dynamics models, and survival verification scripts are openly accessible at https://github.com/Lukecele/spatial-oncology-master and permanently archived on Zenodo at https://doi.org/10.5281/zenodo.23623653.

---

## Appendix A. Changes from Version 1

- **Corrected & Withdrawn:**
  1. Withdrew the assertion that inward (avascular nest) and outward (perivascular cord) geometries yield identical sanctuary bands; proved that inward diffusion maintains osimertinib at $91\%\text{--}99\%$ of boundary concentration, precluding sanctuary formation.
  2. Withdrew the claim that avascular solid/micropapillary nests shelter sensitive cells from small-molecule TKIs; clarified that sanctuaries exist exclusively in outward perivascular cords or for macromolecular ADCs.
  3. Withdrew the claim of ADAURA concordance as mechanistic validation; proved that lead time $\Delta t = 5.46\,\text{months}$ is an algebraic identity ($\log_2(100) \cdot T_s$).
  4. Withdrew the cross-subtype correlation coefficient ($r = 0.983$), recognizing that sorting categories by their own event rates creates a circular artifact.
  5. Withdrew the designation of the solid/micropapillary vs lepidic contrast as "a priori" or "pre-specified"; reported it transparently as an exploratory boundary comparison.
  6. Replaced raw median RFS follow-up ($21.4\,\text{months}$) with the reverse Kaplan–Meier follow-up estimate ($23.8\,\text{months}$).
  7. Corrected the post-cessation event distribution from Herbst et al. 2025: 11 of 19 events occurred within 12 months post-cessation.
  8. Removed speculative clinical surveillance recommendations, the $\ge 20\%$ ctDNA monitoring interval, and combination therapy claims, replacing them with explicit testable predictions.

- **Added:**
  1. Closed-form inward Bessel solution ($I_0(r/\lambda)/I_0(R/\lambda)$) alongside outward solution ($K_0(r/\lambda)/K_0(r_0/\lambda)$).
  2. Exhaustive multidimensional sensitivity analysis across $D, k, \theta, r_0$, and $r_{\text{ox}}$ (Table 1, Table 2).
  3. Analytical lead-time identity derivation (Table 3).
  4. Pairwise Cox contrasts with Holm multiplicity adjustments (Table 4, Table 5).
  5. Bootstrap confidence intervals demonstrating the instability of the extreme contrast upper bound.
  6. Demonstration of severe stage confounding within the *EGFR*-mutant subgroup.

---

## References

1. Wu YL, Tsuboi M, He J, et al. Osimertinib in Resected *EGFR*-Mutated Non-Small-Cell Lung Cancer. *N Engl J Med*. 2020;383(18):1711-1723. doi:10.1056/NEJMoa2027071.
2. Tsuboi M, Herbst RS, John T, et al. Overall Survival with Osimertinib in Resected *EGFR*-Mutated NSCLC. *N Engl J Med*. 2023;389(2):137-147. doi:10.1056/NEJMoa2304594.
3. Caso R, Sanchez-Vega F, Tan KS, et al. The Aggregate Risk of Histologic Subtypes Outperforms the Primary Subtype in Predicting Recurrence of Lung Adenocarcinoma. *J Thorac Oncol*. 2020;15(11):1747-1755. doi:10.1016/j.jtho.2020.07.019.
4. Thomlinson RH, Gray LH. The histological structure of some human lung cancers and the possible implications for radiotherapy. *Br J Cancer*. 1955;9(4):539-549.
5. Cross DAE, Ashton SE, Ghiorghiu S, et al. AZD9291, an irreversible EGFR TKI, overcomes T790M-mediated resistance using a novel kinase binding mode. *Cancer Discov*. 2014;4(9):1046-1061.
6. Sharma SV, Lee DY, Li B, et al. A chromatin-mediated reversible drug-tolerant state in cancer cell subpopulations. *Cell*. 2010;141(1):69-80.
7. Kurppa KJ, Liu Y, To C, et al. Treatment-Induced Emergence of a Drug-Tolerant Persister State in *EGFR*-Mutant Lung Cancer Sinks into Tumor Dormancy through YAP-Mediated Transcriptional Reprogramming. *Cancer Cell*. 2020;37(1):104-122.e12. doi:10.1016/j.ccell.2019.12.006.
8. Boucher Y, Baxter LT, Jain RK. Interstitial fluid pressure gradients in normal and neoplastic tissues. *Cancer Res*. 1990;50(15):4478-4484.
9. Jain RK. Normalization of tumor vasculature: an emerging concept in antiangiogenic therapy. *Science*. 2005;307(5706):58-62.
10. Cho BC, Felip E, Spira AI, et al. Amivantamab plus Lazertinib in First-Line *EGFR*-Mutated Advanced NSCLC (MARIPOSA). *N Engl J Med*. 2024;391(16):1486-1498.
11. Jänne PA, Baik C, Su WC, et al. Efficacy and Safety of Patritumab Deruxtecan (HER3-DXd) in Advanced *EGFR*-Mutated Non-Small Cell Lung Cancer (HERTHENA-Lung01). *J Clin Oncol*. 2024;42(20):2410-2420.
12. Herbst RS, John T, Grohé C, et al. Molecular residual disease analysis of adjuvant osimertinib in resected *EGFR*-mutated stage IB–IIIA non-small-cell lung cancer. *Nat Med*. 2025;31(6):1958-1968. doi:10.1038/s41591-025-03577-y.
13. Swabb EA, Wei J, Gullino PM. Diffusion and convection in normal and neoplastic tissues. *Cancer Res*. 1974;34(10):2814-2822.
14. Worth KR, Papandreou I, Hammond EM. How the histological structure of some lung cancers shaped almost 70 years of radiobiology. *Br J Cancer*. 2022;127(7):1171-1178. doi:10.1038/s41416-022-02041-9.
15. Hirst DG, Hirst VK, Joiner B, et al. Changes in tumour morphology with alterations in oxygen availability: further evidence for oxygen as a limiting substrate. *Br J Cancer*. 1991;64(4):689-693.
16. Moreira AL, Ocampo PSS, Xia Y, et al. A Grading System for Invasive Pulmonary Adenocarcinoma: A Proposal From the International Association for the Study of Lung Cancer Pathology Committee. *J Thorac Oncol*. 2020;15(10):1599-1610. doi:10.1016/j.jtho.2020.06.001.
