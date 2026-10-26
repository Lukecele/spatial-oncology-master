# Conditions for Pharmacological Sanctuaries in Perivascular Tumour Cords: A Reaction–Diffusion Analysis, a Regrowth-Kinetics Identity, and a Descriptive Survival Analysis of Histologic Subtype in Resected Lung Adenocarcinoma

**Luca Celebrano**$^{1,2}$  
$^{1}$Department of Clinical Medicine and Surgery, School of Medicine, University of Naples Federico II, Via Sergio Pansini 5, 80131 Naples, Italy  
$^{2}$Telethon Institute of Genetics and Medicine (TIGEM), Via Campi Flegrei 34, 80078 Pozzuoli, Naples, Italy  
*Correspondence:* luca.celebrano5@gmail.com / l.celebrano@studenti.unina.it | ORCID: [0009-0007-0931-3555](https://orcid.org/0009-0007-0931-3555)

---

### Abstract
**Background:** Adjuvant osimertinib markedly improves disease-free survival in resected *EGFR*-mutated non-small-cell lung cancer (NSCLC). However, extended follow-up in the ADAURA trial revealed that among patients who experienced disease recurrence, the majority of relapses manifested after planned treatment cessation at year 3, frequently without acquired secondary genomic resistance mutations. We examined whether biophysical drug transport limitations could establish tissue sanctuaries that shelter drug-sensitive tumor cells.

**Methods:** We analyzed steady-state reaction–diffusion transport in two canonical tumor micro-architectures: outward perivascular cords (divergent cylindrical geometry, modified Bessel $K_0$) and inward avascular nests (convergent cylindrical geometry, modified Bessel $I_0$). We derived an exact algebraic identity governing post-cessation exponential regrowth from molecular residual disease (ctDNA, $10^7$ cells) to radiologic detection (CT, $10^9$ cells). We evaluated the clinical prognostic relevance of histologic architecture in 604 resected lung adenocarcinoma patients from Memorial Sloan Kettering Cancer Center (MSKCC 2020).

**Results:** In perivascular cords, outward radial diffusion of osimertinib ($\lambda = 627.6\,\mu\text{m}$) produces a steep concentration drop beneath the $\text{IC}_{50}$ threshold ($\theta = 0.35$, corresponding to unbound $C_0 \approx 40\,\text{nM}$ and $\text{IC}_{50} \approx 14\,\text{nM}$) at $r_{\text{crit}} = 165.9\,\mu\text{m}$. Depending on physiological oxygenation limits, a viable pharmacological sanctuary shell of width $\Delta r = 34.1\,\mu\text{m}$ (for $r_{\text{ox}} = 200\,\mu\text{m}$) or $\Delta r = 14.1\,\mu\text{m}$ (for $r_{\text{ox}} = 180\,\mu\text{m}$) can persist. Under an assumed sigmoidal dose-response ($h = 1.5$), residual drug exposure in the shell ($C \approx 0.31\text{--}0.35\,C_0$) maintains partial target inhibition (~45–50%), sufficient under physiological baseline apoptosis to induce near-zero net growth (cytostasis). Conversely, in inward avascular nests, convergent geometry maintains osimertinib at $90.6\%\text{--}99.4\%$ of boundary levels throughout viable regions ($R \le 160\text{--}200\,\mu\text{m}$), precluding small-molecule sanctuaries. For macromolecular antibody–drug conjugates (ADCs, $\lambda = 77.5\,\mu\text{m}$), central concentration drops below $\theta = 0.35$ only for nest radii $R > 180\,\mu\text{m}$, an architectural domain typically constrained by hypoxia and necrosis. In Monte Carlo parameter sampling ($N = 10,000$ draws, log-normal priors, seed = 2), an outward perivascular sanctuary shell existed in 46.1% of parameter draws (62.8% when $r_{\text{ox}} = 200\,\mu\text{m}$ is fixed). The diagnostic lead time $\Delta t = T_s \cdot \log_2(100) = 5.46\,\text{months}$ (for doubling time $T_s = 25\,\text{days}$) is an analytical property of exponential expansion and assay thresholds, rather than an empirical confirmation of sanctuary existence. In the MSKCC cohort (reverse Kaplan–Meier median follow-up 23.8 months, IQR: 13.3–31.8 months), lepidic architecture was consistent with a substantially lower hazard of recurrence (4.55%, stage-adjusted $HR = 0.25$, $95\%\,\text{CI}: 0.09\text{--}0.70$, $p = 0.0084$ vs acinar). In 10-fold Holm-adjusted pairwise comparisons, lepidic recurrence was significantly lower than solid ($p_{\text{adj}} = 0.0026$) and papillary ($p_{\text{adj}} = 0.0203$), whereas differences vs acinar ($p_{\text{adj}} = 0.0703$) and micropapillary ($p_{\text{adj}} = 0.3911$, unadjusted $p = 0.0559$) did not achieve significance under family-wise error control. Grouping solid and micropapillary patterns yielded an exploratory contrast of $HR = 6.75$ ($95\%\,\text{CI}: 2.28\text{--}19.97$; Firth penalized shrinkage $HR = 2.60\text{--}3.42$) against lepidic, based on 4 lepidic events, whereas across the unselected cohort the hazard ratio was $1.32$ ($p = 0.239$).

**Conclusions:** Small-molecule pharmacological sanctuaries are geometrically restricted to outward perivascular cords and do not occur in avascular nests. Clinical cohort data demonstrate baseline protection associated with lepidic growth rather than an escalating risk across compact patterns. Surveillance strategies based on molecular lead time must be evaluated in prospective trials to rule out lead-time bias.

**Keywords:** Pharmacological sanctuary, Reaction-diffusion, Modified Bessel functions, Tumor cords, Non-small-cell lung cancer, EGFR tyrosine kinase inhibitors, ADAURA, IASLC histologic grading.

---

## 1. Introduction

Adjuvant targeted therapy with the third-generation epidermal growth factor receptor (EGFR) tyrosine kinase inhibitor (TKI) osimertinib produces profound disease-free and overall survival benefits in patients with resected stage IB–IIIA *EGFR*-mutated non-small-cell lung cancer (NSCLC) (Wu et al., *N Engl J Med* 2020 [1]; Tsuboi et al., *N Engl J Med* 2023 [2]). In the landmark ADAURA trial, three years of daily oral osimertinib suppressed disease recurrence by nearly 80% during active treatment. However, extended follow-up identified an important temporal pattern: following planned therapy completion at month 36, disease-free survival curves showed an inflection, with a substantial portion of recurrences emerging in the post-cessation observation window (months 36 to 60) (Herbst et al., *Nat Med* 2025 [12]). In an exploratory molecular residual disease (MRD) analysis of 220 enrolled patients (112 in the osimertinib arm), 28 patients in the osimertinib arm developed disease recurrence (25.0%). Among these 28 patients, 19 events (67.9%) occurred after treatment discontinuation; within this post-cessation window, 11 of 19 events (57.9%) manifested within the first 12 months post-cessation (months 36–48). Crucially, post-cessation recurrences frequently lacked canonical secondary genetic resistance alterations (such as *EGFR* T790M, C797S, or *MET* amplification), and clinical rechallenge with osimertinib in subsequent trials has demonstrated sustained disease control.

This clinical phenomenon suggests that residual cancer cells may persist through non-genetic mechanisms. While epigenetic reprogramming and drug-tolerant persister (DTP) states have been extensively characterized in vitro (Sharma et al., *Cell* 2010 [6]; Kurppa et al., *Cancer Cell* 2020 [7]), biological persisters can also undergo slow-cycling reversible quiescence. An alternative and complementary biophysical hypothesis is that microenvironmental transport barriers create spatial pharmacological sanctuaries: tissue niches where local drug exposure drops beneath cytotoxic or cytostatic thresholds while physiological nutrient and oxygen supply remains adequate to support cell survival.

In this study, we evaluate this hypothesis by:
1. Solving exact cylindrical reaction–diffusion boundary value problems across outward perivascular cords and inward avascular nests;
2. Deriving the analytical lead-time identity governing post-cessation clonal regrowth;
3. Conducting a descriptive survival audit of histologic subtype in 604 resected lung adenocarcinoma patients from Memorial Sloan Kettering Cancer Center.

---

## 2. Biophysical Reaction-Diffusion Transport Physics

### 2.1 Governing Equations & Geometric Formulations

We model steady-state drug diffusion coupled to first-order irreversible cellular uptake in cylindrical coordinates:

$$\frac{1}{r} \frac{d}{dr}\left( r \frac{dC}{dr} \right) - \frac{k}{D} C = 0$$

where $C(r)$ is local interstitial concentration, $D$ is effective tissue diffusivity, and $k$ is the volumetric uptake rate constant. The characteristic penetration length is defined by:

$$\lambda = \sqrt{\frac{D}{k}}$$

#### Geometry 1: Outward Perivascular Cord (Divergent Cylindrical Transport)
A central microvessel of radius $r_0 = 10\,\mu\text{m}$ supplies drug at constant boundary concentration $C_0$ to a surrounding cylindrical tumor cord. With boundary conditions:

$$C(r_0) = C_0, \quad \lim_{r \to \infty} C(r) = 0$$

the unique closed-form analytical solution is given by the modified Bessel function of the second kind of order zero ($K_0$):

$$C_{\text{out}}(r) = C_0 \frac{K_0(r / \lambda)}{K_0(r_0 / \lambda)}$$

The critical radius $r_{\text{crit}}$ at which local concentration falls below the threshold $\theta = C / C_0 = \text{IC}_{50} / C_0 = 0.35$ is determined implicitly by:

$$\frac{K_0(r_{\text{crit}} / \lambda)}{K_0(r_0 / \lambda)} = \theta$$

#### Geometry 2: Inward Avascular Nest (Convergent Cylindrical Transport)
A compact avascular tumor cluster of radius $R$ is supplied by surrounding microvasculature at the outer boundary $r = R$. Drug diffuses radially inward toward the center:

$$C(R) = C_0, \quad \left.\frac{dC}{dr}\right|_{r=0} = 0$$

The unique regular analytical solution is given by the modified Bessel function of the first kind of order zero ($I_0$):

$$C_{\text{in}}(r) = C_0 \frac{I_0(r / \lambda)}{I_0(R / \lambda)}$$

The minimum concentration occurs at the nest center ($r = 0$), where $I_0(0) = 1$:

$$\frac{C_{\text{in}}(0)}{C_0} = \frac{1}{I_0(R / \lambda)}$$

#### Oxygen Transport and Viability Boundaries: Cord vs Nest
Thomlinson & Gray (1955) [4] demonstrated that tumor cords surrounding central blood vessels are bounded by an outer zone of necrosis caused by oxygen consumption. With oxygen diffusivity $D_{\text{ox}} \approx 1.8 \times 10^{-5}\,\text{cm}^2/\text{s}$ and volumetric consumption $M_{\text{ox}} \approx 2.5 \times 10^{-9}\,\text{mol}/(\text{cm}^3\cdot\text{s})$ in neoplastic cords (Hirst et al. [16]; Worth et al. [15]), steady-state zero-order diffusion in a Krogh cylinder satisfies:

$$\frac{1}{r} \frac{d}{dr}\left( r \frac{dC_{\text{ox}}}{dr} \right) = \frac{M_{\text{ox}}}{D_{\text{ox}}}$$

Integrating with boundary conditions $C_{\text{ox}}(r_0) = C_{\text{vessel}}$ and zero flux at the cord outer radius $\left.\frac{dC_{\text{ox}}}{dr}\right|_{r = R_{\text{cord}}} = 0$ yields:

$$C_{\text{ox}}(r) = C_{\text{vessel}} + \frac{M_{\text{ox}}}{4 D_{\text{ox}}} (r^2 - r_0^2) - \frac{M_{\text{ox}} R_{\text{cord}}^2}{2 D_{\text{ox}}} \ln\left(\frac{r}{r_0}\right)$$

Setting $C_{\text{ox}}(R_{\text{cord}}) = 0$ at the cord radius defines the physiological limit of tissue oxygenation, established empirically and theoretically at $r_{\text{ox}} \approx 180\text{--}200\,\mu\text{m}$ (Thomlinson & Gray [4]; Hirst et al. [16]; Worth et al. [15]), corresponding to perfused capillary oxygen concentrations $C_{\text{vessel}} \approx 54\text{--}70\,\mu\text{M}$ ($pO_2 \approx 40\text{--}52\,\text{mmHg}$). Beyond $r_{\text{ox}}$, anoxia precipitates coagulative necrosis.

In inward avascular nests, oxygen diffuses inward from the surrounding stroma at $r = R_{\text{nest}}$ with boundary concentration $C_{\text{stroma}}$. Integrating with symmetry condition $dC_{\text{ox}}/dr = 0$ at $r = 0$ yields:

$$C_{\text{ox}}(r) = C_{\text{stroma}} - \frac{M_{\text{ox}}}{4 D_{\text{ox}}} (R_{\text{nest}}^2 - r^2)$$

Setting $C_{\text{ox}}(0) = 0$ gives the critical nest radius for central anoxia:

$$R_{\text{crit,ox}} = \sqrt{\frac{4 D_{\text{ox}} C_{\text{stroma}}}{M_{\text{ox}}}}$$

Because the vascular stroma surrounding avascular tumor nests in poorly perfused tumor centers exhibits low oxygen tension ($pO_2 \approx 12\text{--}20\,\text{mmHg}$, corresponding to $C_{\text{stroma}} \approx 16\text{--}27\,\mu\text{M}$) compared to direct capillary supply ($C_{\text{vessel}} \approx 54\text{--}70\,\mu\text{M}$), central anoxia in avascular nests emerges at $R_{\text{crit,ox}} \approx 150\text{--}190\,\mu\text{m}$.

### 2.2 Baseline Biophysical Parameters & Sanctuary Shell Geometry

Effective tissue diffusivity ($D = 1.3 \times 10^{-6}\,\text{cm}^2/\text{s}$) and volumetric cellular uptake constant ($k = 3.3 \times 10^{-4}\,\text{s}^{-1}$) are parameterized from classical extravascular measurements of small-molecule transport in neoplastic tissue (Swabb et al. [13]; Nugent & Jain [14]). This yields a characteristic penetration length of:

$$\lambda = \sqrt{\frac{1.3 \times 10^{-6}}{3.3 \times 10^{-4}}} = 0.06276\,\text{cm} = 627.6\,\mu\text{m}$$

In clinical pharmacokinetics, steady-state total plasma concentrations of osimertinib (80 mg once daily) reach $500\text{--}600\,\text{nM}$. Because osimertinib is ~95% plasma protein bound, the therapeutic steady-state unbound (free) interstitial concentration is $C_0 \approx 40\,\text{nM}$. Given in vitro $\text{IC}_{50} \approx 14\,\text{nM}$ for EGFR phosphorylation inhibition in mutant NSCLC cells (Cross et al. [5]), the baseline inhibitory threshold ratio is:

$$\theta = \frac{\text{IC}_{50}}{C_0} = \frac{14\,\text{nM}}{40\,\text{nM}} = 0.350$$

With $r_0 = 10\,\mu\text{m}$ and $\theta = 0.350$, modified Bessel root finding yields:

$$r_{\text{crit}} = 165.9\,\mu\text{m}$$

Because viable oxygenated tissue extends to $r_{\text{ox}}$, a perivascular pharmacological sanctuary shell exists over the interval $r \in [r_{\text{crit}}, r_{\text{ox}}]$:
- For $r_{\text{ox}} = 200\,\mu\text{m}$: $\Delta r = 200 - 165.9 = 34.1\,\mu\text{m}$
- For $r_{\text{ox}} = 180\,\mu\text{m}$: $\Delta r = 180 - 165.9 = 14.1\,\mu\text{m}$

#### Sigmoidal Dose-Response and Residual Inhibition
Within the sanctuary shell, concentration decays from $0.350\,C_0$ at $r_{\text{crit}} = 165.9\,\mu\text{m}$ to $0.310\,C_0$ at $r_{\text{ox}} = 200\,\mu\text{m}$ ($K_0(200/627.6)/K_0(10/627.6) = 0.30955$). The unbound drug concentration throughout the shell remains at $12.4\text{--}14.0\,\text{nM}$. Under an assumed standard Hill pharmacodynamic model:

$$E(C) = \frac{C^h}{\text{IC}_{50}^h + C^h}$$

with assumed Hill slope $h \approx 1.5$, local target inhibition remains at approximately $45\%\text{--}50\%$. The net cellular growth rate is governed by:

$$g_{\text{net}} = g_0 \left[ 1 - E(C) \right] - d_0$$

where $g_0$ is uninhibited proliferation and $d_0$ is baseline apoptosis and immune clearance. Net growth depends directly on the ratio $d_0 / g_0$:
- If $d_0 / g_0 \approx 0.50$, $g_{\text{net}} \approx 0$, maintaining prolonged cytostasis and near-zero net expansion during 36 months of active therapy.
- If $d_0 / g_0 \approx 0.40$, $g_{\text{net}} \approx +0.10\,g_0$, permitting slow indolent expansion.
- If $d_0 / g_0 \ge 0.55$, $g_{\text{net}} < 0$, leading to gradual clonal eradication.

#### Volumetric Scaling: From Cord Shell to Initial Residual Burden $N_0$
In a cylindrical perivascular cord of unit length $L = 1\,\text{mm} = 1000\,\mu\text{m}$, the volume of the sanctuary shell is:

$$V_{\text{shell}} = \pi \left( r_{\text{ox}}^2 - r_{\text{crit}}^2 \right) L = \pi \left( 200^2 - 165.9^2 \right) \times 1000 \approx 3.92 \times 10^7\,\mu\text{m}^3$$

Assuming a mean tumor cell diameter of $15\,\mu\text{m}$ (cell volume $\approx 1.77 \times 10^3\,\mu\text{m}^3$, corresponding to packing density $\sim 5 \times 10^8\,\text{cells/cm}^3 = 5 \times 10^{-4}\,\text{cells}/\mu\text{m}^3$), the shell contains:

$$N_{\text{cord}} \approx 3.92 \times 10^7 \times 5 \times 10^{-4} \approx 1.96 \times 10^4\,\text{cells per mm of cord}$$

For microscopic residual disease consisting of a few surviving cords following surgical resection, this physical volume corresponds directly to an initial residual cell pool of $N_0 \approx 10^3\text{--}10^5\,\text{cells}$.

---

### 2.3 Sensitivity Analysis & Parameter Bounds

#### Table 1. Outward Perivascular Cord Sensitivity (Baseline: $\theta = 0.350, r_0 = 10\,\mu\text{m}$)

| Scenario | $\lambda$ ($\mu\text{m}$) | $r_{\text{crit}}$ ($\mu\text{m}$) | Shell Width $\Delta r$ ($r_{\text{ox}} = 200\,\mu\text{m}$) | Shell Width $\Delta r$ ($r_{\text{ox}} = 180\,\mu\text{m}$) |
| :--- | :---: | :---: | :---: | :---: |
| **Baseline** | **627.6** | **165.9** | **34.1** | **14.1** |
| $D$ +30% | 715.6 | 180.1 | 19.9 | None |
| $D$ -30% | 525.1 | 148.5 | 51.5 | 31.5 |
| $k$ +30% | 550.5 | 152.9 | 47.1 | 27.1 |
| $k$ -30% | 750.2 | 185.5 | 14.5 | None |
| $D$ +30% and $k$ -30% | 855.3 | 201.4 | None | None |
| $D$ -30% and $k$ +30% | 460.6 | 136.9 | 63.1 | 43.1 |
| Vessel radius $r_0 = 5\,\mu\text{m}$ | 627.6 | 128.3 | 71.7 | 51.7 |
| Vessel radius $r_0 = 20\,\mu\text{m}$ | 627.6 | 216.2 | None | None |

As shown in **Table 1**, the sanctuary shell persists across isolated $\pm 30\%$ parameter variations at $r_{\text{ox}} = 200\,\mu\text{m}$, but at the more restrictive oxygenation boundary $r_{\text{ox}} = 180\,\mu\text{m}$, the shell collapses under $D + 30\%$ ($r_{\text{crit}} = 180.1\,\mu\text{m}$) or $k - 30\%$ ($r_{\text{crit}} = 185.5\,\mu\text{m}$) alone.

Sensitivity to the threshold ratio $\theta$:
- If $\theta = 0.333$ ($15\,\text{nM} / 45\,\text{nM}$), $r_{\text{crit}} = 179.1\,\mu\text{m}$, yielding $\Delta r = 20.9\,\mu\text{m}$ at $r_{\text{ox}} = 200\,\mu\text{m}$ and $\Delta r = 0.9\,\mu\text{m}$ at $r_{\text{ox}} = 180\,\mu\text{m}$.
- If $\theta \le 0.300$ ($15\,\text{nM} / 50\,\text{nM}$), $r_{\text{crit}} = 209.2\,\mu\text{m} > 200\,\mu\text{m}$, and the sanctuary shell collapses.

#### Table 2. Inward Avascular Nest Geometry: Central Core Concentration $C_{\text{in}}(0) / C_0 = 1 / I_0(R / \lambda)$

| Nest Outer Radius $R$ ($\mu\text{m}$) | Osimertinib ($\lambda = 627.6\,\mu\text{m}$) | Macromolecule ADC ($\lambda = 77.5\,\mu\text{m}$) | Tissue Viability Status |
| :---: | :---: | :---: | :--- |
| $100\,\mu\text{m}$ | **0.994** (99.4%) | 0.684 (68.4%) | Viable throughout ($R < R_{\text{crit,ox}}$) |
| $150\,\mu\text{m}$ | **0.986** (98.6%) | 0.459 (45.9%) | Viable throughout ($R \approx R_{\text{crit,ox}}$) |
| $200\,\mu\text{m}$ | **0.975** (97.5%) | 0.286 (28.6%) | Central hypoxia / necrosis ($R \ge R_{\text{crit,ox}}$) |
| $300\,\mu\text{m}$ | **0.945** (94.5%) | 0.099 (9.9%) | Extensive central necrosis |
| $400\,\mu\text{m}$ | **0.906** (90.6%) | 0.032 (3.2%) | Extensive central necrosis |

As shown in **Table 2**, for osimertinib ($\lambda = 627.6\,\mu\text{m}$), central core concentration remains at $\ge 90.6\%$ for any nest radius up to $400\,\mu\text{m}$, and $\ge 97.5\%$ for $R \le 200\,\mu\text{m}$. Consequently, **small-molecule osimertinib never drops beneath therapeutic thresholds in avascular nests of any viable dimension**. For macromolecular ADCs ($\lambda = 77.5\,\mu\text{m}$), central concentration falls below $\theta = 0.35$ only when $R > 180\,\mu\text{m}$. Because avascular tumor clusters typically develop central hypoxia and necrosis at $R \approx 150\text{--}190\,\mu\text{m}$, macromolecular sanctuaries within viable nest cores exist only if local stromal oxygenation sustains viability beyond $180\,\mu\text{m}$.

#### Monte Carlo Analysis of Parameter Uncertainty
Because tissue transport parameters cannot be measured non-invasively in individual patients, we evaluated the prior probability of sanctuary existence using Monte Carlo sampling ($N = 10,000$ draws, fixed pseudorandom seed = 2). Parameters were sampled from log-normal distributions centered at baseline values: $D$ and $k$ with 90% confidence spans spanning a 3-fold factor ($\sigma = \ln(3)/1.645 \approx 0.668$), $\theta$ log-normally distributed with a 2-fold factor ($\sigma = \ln(2)/1.645 \approx 0.421$), and oxygen limit $r_{\text{ox}}$ distributed uniformly between $100$ and $200\,\mu\text{m}$ (Thomlinson & Gray [4]; Worth et al. [15]).

Under this prior distribution, an outward perivascular sanctuary shell existed in **46.1% of parameter draws** (and in **62.8% of draws** if $r_{\text{ox}}$ was held fixed at $200\,\mu\text{m}$). This 46.1% figure represents the prior probability under our assumed synthetic parameter distributions, rather than an empirical patient prevalence.

---

## 3. Clonal Regrowth & The Analytical Lead-Time Identity

### 3.1 Replicator Dynamics as an Exploratory Model of Subclonal Coexistence

Within the perivascular sanctuary shell during active treatment, sub-inhibitory drug exposure allows sensitive proliferating cells ($S$) and drug-tolerant persister cells ($P$) to interact. We formalized this interaction using an illustrative continuous replicator equation (Kaznatcheev et al. [8]):

$$\frac{dx}{dt} = x(1 - x)\left[ f_S(x) - f_P(x) \right]$$

where $x$ represents the fraction of sensitive cells. To illustrate subclonal coexistence, we adopted an illustrative payoff matrix:

$$A = \begin{pmatrix} -0.736 & 2.866 \\ 2.188 & -1.772 \end{pmatrix}$$

yielding $\alpha = a_{11} - a_{12} - a_{21} + a_{22} = -7.562$ and $\beta = a_{12} - a_{22} = 4.638$. The interior fixed point is:

$$x^* = \frac{\beta}{-\alpha} = \frac{4.638}{7.562} = 0.6133$$

with Jacobian eigenvalue $\lambda_{\text{eig}} = x^*(1 - x^*)\alpha = -1.793 < 0$, confirming asymptotic Lyapunov stability. We emphasize that this game-theoretic formulation is an exploratory, uncalibrated theoretical model used to illustrate stable coexistence; the payoffs are free parameters and do not represent empirical fits from clinical trials.

### 3.2 Post-Cessation Regrowth & The Lead-Time Identity

Following adjuvant therapy cessation at month 36, drug clears rapidly from tissue. Residual tumor burden expands from an initial microscopic inoculum $N_0$ with sensitive fraction $x^* = 0.6133$:

$$N(t) = N_0 \left[ x^* \cdot 2^{t / T_s} + (1 - x^*) \cdot 2^{t / T_p} \right]$$

where $T_s = 25.0\,\text{days}$ is the doubling time of sensitive cells and $T_p = 120.0\,\text{days}$ is the doubling time of persister cells.

Because $T_s \ll T_p$, the sensitive clone rapidly expands, comprising 97.7% of the burden at $t = 150\,\text{days}$ and $>99.9\%$ by $t \approx 294\,\text{days}$. Consequently, by the time tumor burden crosses the molecular residual disease threshold ($N_1 = 10^7$ cells, corresponding approximately to ctDNA detection) and radiologic imaging detection ($N_2 = 10^9$ cells, $\approx 1\,\text{cm}^3$), growth is governed almost exclusively by the sensitive clone:

$$N(t) \approx N_0 x^* \cdot 2^{t / T_s}$$

Whenever the proliferating sensitive clone dominates the expansion between detection thresholds, the diagnostic lead time $\Delta t$ between ctDNA detection ($t_1$) and radiologic detection ($t_2$) is:

$$\Delta t = t_2 - t_1 = T_s \cdot \log_2\left( \frac{N_2}{N_1} \right) = T_s \cdot \log_2(100) \approx 6.6439 \cdot T_s$$

This result is an exact algebraic identity of exponential growth. It is strictly invariant to initial burden $N_0$, and invariant to the initial persister fraction so long as sensitive cells dominate during the threshold crossing.

#### Table 3A. The Analytical Lead-Time Identity Across Doubling Times $T_s$

| Assumed Sensitive Doubling Time $T_s$ | Diagnostic Lead Time $\Delta t$ (Days) | Diagnostic Lead Time $\Delta t$ (Months) |
| :---: | :---: | :---: |
| 10.0 days | 66.4 days | 2.18 months |
| 15.0 days | 99.7 days | 3.27 months |
| 20.0 days | 132.9 days | 4.37 months |
| **21.5 days** (matches Herbst median) | **142.8 days** | **4.70 months** |
| **25.0 days** (baseline assumed) | **166.1 days** | **5.46 months** |
| 30.0 days | 199.3 days | 6.55 months |
| 40.0 days | 265.8 days | 8.73 months |

#### Table 3B. Relapse Timing from Cessation Across Residual Inocula $N_0$ ($T_s = 25\,\text{days}$)

| Initial Inoculum $N_0$ | ctDNA Time $t_1$ ($10^7$ cells) | Imaging Time $t_2$ ($10^9$ cells) | Lead Time $\Delta t$ ($t_2 - t_1$) |
| :---: | :---: | :---: | :---: |
| 100 cells | 433.8 days (14.25 mo) | 599.9 days (19.71 mo) | 166.1 days (5.46 mo) |
| 1,000 cells | 350.7 days (11.52 mo) | 516.8 days (16.98 mo) | 166.1 days (5.46 mo) |
| **2,000 cells** | **324.8 days (10.67 mo)** | **490.9 days (16.13 mo)** | **166.1 days (5.46 mo)** |
| 10,000 cells | 267.6 days (8.79 mo) | 433.7 days (14.25 mo) | 166.1 days (5.46 mo) |
| 100,000 cells | 184.6 days (6.06 mo) | 350.7 days (11.52 mo) | 166.1 days (5.46 mo) |

As shown in **Table 3B**, varying $N_0$ shifts the calendar timing of recurrence relative to treatment cessation, but leaves the lead time $\Delta t = 5.46\,\text{months}$ invariant. The agreement between our model lead time (5.46 months) and the median lead time reported by Herbst et al. (4.7 months, 95% CI: 2.2–5.6) corresponds to an implied sensitive doubling time of $T_s \approx 21.5\,\text{days}$ (range 10–26 days). Inverting the identity to recover a plausible biological doubling time illustrates parameter consistency, rather than an independent mechanistic confirmation of sanctuary dynamics. Furthermore, ctDNA shedding varies across anatomical sites and vascularization, and earlier detection inherently introduces lead-time bias that must not be conflated with genuine therapeutic benefit.

---

## 4. Descriptive Survival Analysis of the MSKCC Cohort ($N = 604$)

### 4.1 Cohort Characteristics & Methods

We conducted a descriptive survival analysis of the authenticated Memorial Sloan Kettering Cancer Center (MSKCC 2020) lung adenocarcinoma cohort ($N = 604$ surgically resected patients; Caso et al. [3]). The cohort includes predominant IASLC histologic subtype, pathologic stage (I, II, III), *EGFR* mutation status, and recurrence-free survival (RFS). The cohort represents surgery-era patients who were chemo/TKI-naive (only 3 patients received adjuvant TKI, none osimertinib). Only the predominant histologic pattern is available, precluding assessment of secondary patterns at the $\ge 20\%$ IASLC grading threshold (Moreira et al. [17]). Adjuvant systemic chemotherapy administered in stage II–III patients represents an unmeasured clinical confounder.

Survival analysis was performed using Kaplan–Meier estimators and multivariable Cox proportional hazards models with Efron tie handling. True median follow-up was computed using the **reverse Kaplan–Meier method** (Schemper & Smith [11]): **median follow-up was 23.8 months** (IQR: 13.3–31.8 months; raw observed RFS duration median was 21.4 months). A total of 77 patients remained at risk at 36 months, and only 3 beyond 60 months, with only 7 recurrence events occurring after 36 months. Consequently, this surgical cohort cannot directly inform recurrence kinetics occurring between months 36 and 60 post-cessation.

Statistical code was validated against the benchmark Freireich leukemia dataset (Freireich et al. [9]; Gehan [10]), reproducing authenticated values ($\chi^2 = 16.79$, Efron $HR = 4.82$).

### 4.2 Cohort Recurrence Patterns

#### Table 4. Recurrence and Stage Distribution by Predominant Histologic Subtype

| Predominant Subtype | Total $N$ | Recurrence Events | Crude Recurrence % | Pathologic Stage I % |
| :--- | :---: | :---: | :---: | :---: |
| **Lepidic** | 88 | 4 | **4.55%** | **88.6%** (78/88) |
| **Acinar** | 368 | 54 | 14.67% | 75.5% (278/368) |
| **Papillary** | 43 | 11 | 25.58% | 76.7% (33/43) |
| **Micropapillary** | 37 | 8 | 21.62% | 51.4% (19/37) |
| **Solid** | 68 | 21 | 30.88% | 57.4% (39/68) |
| **Total Cohort** | 604 | 98 | 16.23% | 74.0% (447/604) |

An unadjusted log-rank test across all 5 histologic subtypes showed significant differences in recurrence-free survival (4 degrees of freedom: $\chi^2 = 24.52, p = 6.27 \times 10^{-5}$). However, stage distribution is highly skewed: 88.6% of lepidic tumors were Stage I, compared to 57.4% of solid and 51.4% of micropapillary tumors.

### 4.3 Multivariable Cox Models & Multiplicity Adjustments

#### Table 5. Multivariable Stage-Adjusted Cox Proportional Hazards Models

| Model & Comparison | Patients ($N$) | Events | Hazard Ratio ($95\%\,\text{CI}$) | Raw $p$-value | Holm-Adjusted $p$-value |
| :--- | :---: | :---: | :---: | :---: | :---: |
| **5-Subtype Model (Ref: Acinar)** | 604 | 98 | | | |
| Lepidic vs Acinar | | | 0.25 (0.09–0.70) | 0.0084 | — |
| Solid vs Acinar | | | 1.48 (0.89–2.49) | 0.1334 | — |
| Papillary vs Acinar | | | 1.51 (0.79–2.92) | 0.2149 | — |
| Micropapillary vs Acinar | | | 0.74 (0.33–1.66) | 0.4684 | — |
| Stage II/III vs Stage I | | | 4.16 (2.73–6.34) | $3.36 \times 10^{-11}$ | — |
| **Pairwise Models (10 Contrasts)** | | | | | |
| Solid vs Lepidic | 156 | 25 | 7.82 (2.60–23.55) | 0.00026 | **0.0026** |
| Papillary vs Lepidic | 131 | 15 | 6.07 (1.91–19.32) | 0.00226 | **0.0203** |
| Acinar vs Lepidic | 456 | 58 | 3.94 (1.41–10.97) | 0.00879 | 0.0703 |
| Micropapillary vs Lepidic | 125 | 12 | 3.85 (0.97–15.34) | 0.05587 | 0.3911 |
| Solid vs Acinar | 436 | 75 | 1.45 (0.87–2.44) | 0.15733 | 0.9440 |
| Papillary vs Acinar | 411 | 65 | 1.53 (0.79–2.97) | 0.20889 | 1.0000 |
| Solid vs Micropapillary | 105 | 29 | 1.72 (0.72–4.10) | 0.22416 | 1.0000 |
| Micropapillary vs Acinar | 405 | 62 | 0.67 (0.30–1.50) | 0.33304 | 1.0000 |
| Micropapillary vs Papillary | 80 | 19 | 0.68 (0.25–1.83) | 0.44273 | 1.0000 |
| Solid vs Papillary | 111 | 32 | 1.35 (0.62–2.96) | 0.45435 | 1.0000 |
| **Exploratory Composite Models** | | | | | |
| Solid/Micropapillary vs Lepidic | 193 | 33 | 6.75 (2.28–19.97) | 0.00055 | — |
| ... with Firth Penalized Likelihood ($\lambda = 0.05$) | 193 | 33 | 3.42 (1.66–7.03) | 0.0009 | — |
| ... with Firth Penalized Likelihood ($\lambda = 0.10$) | 193 | 33 | 2.60 (1.42–4.75) | 0.0019 | — |
| Solid/Micropapillary vs All Others | 604 | 98 | 1.32 (0.83–2.08) | 0.2390 | — |
| ... adjusting for EGFR, Sex, Smoking | 604 | 98 | 1.30 (0.82–2.07) | 0.2580 | — |
| **EGFR-Mutant Subgroup** | | | | | |
| Solid/Micropapillary vs Lepidic | 56 | 6 | 1.10 (0.06–19.97) | 0.9473 | — |

As detailed in **Table 5**, the 5-subtype simultaneous model (with Acinar as reference) and pairwise models are separate fits. In the 5-subtype model, lepidic histology exhibits significantly reduced hazard ($HR = 0.25, p = 0.0084$). In the complete 10-fold pairwise analysis with Holm-Bonferroni correction (Holm [18]), lepidic recurrence is significantly lower than solid ($p_{\text{adj}} = 0.0026$) and papillary ($p_{\text{adj}} = 0.0203$), but does not retain statistical significance against acinar ($p_{\text{adj}} = 0.0703$) or micropapillary ($p_{\text{adj}} = 0.3911$) under family-wise error control.

The exploratory contrast of solid/micropapillary vs lepidic ($HR = 6.75$) rests on only 4 events in the lepidic cohort. Under bootstrap resampling (500 iterations), the upper confidence limit was unbounded due to resamples containing zero lepidic events. Applying Firth-type penalized likelihood shrinkage contracts this estimate to $HR = 3.42$ ($95\%\,\text{CI}: 1.66\text{--}7.03$ at penalizer $\lambda = 0.05$) and $HR = 2.60$ ($95\%\,\text{CI}: 1.42\text{--}4.75$ at penalizer $\lambda = 0.10$). Across the full cohort, grouping solid and micropapillary yields an unadjusted $HR = 1.32$ ($p = 0.239$). In the *EGFR*-mutant subgroup ($N = 56, 6$ events), after adjusting for stage, the hazard ratio is $1.10$ ($95\%\,\text{CI}: 0.06\text{--}19.97, p = 0.947$, $z = 0.066$). This lack of significance reflects marked stage confounding: within this subgroup, 37 of 38 lepidic patients (97.4%) were Stage I, compared to only 7 of 18 solid/micropapillary patients (38.9%). The identical upper bound of 19.97 in both models is a coincidental arithmetic property of log-normal exponentiation ($0.098 + 1.96 \times 1.478 \approx 2.995$ and $1.910 + 1.96 \times 0.553 \approx 2.994$, both yielding $\exp(2.995) \approx 19.97$).

---

## 5. Discussion

### 5.1 Biophysical Transport Findings
Our reaction–diffusion analysis demonstrates a fundamental geometric asymmetry between perivascular cords and avascular nests:
1. **Perivascular Cords (Divergent Transport):** Outward diffusion from a central microvessel generates geometric dilution ($K_0$), dropping osimertinib concentration below $\text{IC}_{50}$ ($\theta = 0.350$) at $r_{\text{crit}} = 165.9\,\mu\text{m}$. Because tissue oxygenation extends to $r_{\text{ox}} \approx 180\text{--}200\,\mu\text{m}$, a viable sanctuary shell of width $\Delta r = 14.1\text{--}34.1\,\mu\text{m}$ can persist. Within this shell, partial target inhibition (~45–50%) suppresses net proliferation and maintains cytostasis.
2. **Avascular Nests (Convergent Transport):** Inward diffusion from the surrounding microvasculature produces geometric convergence ($I_0$). For osimertinib, central core concentration remains at $\ge 90.6\%$ for nest radii up to $400\,\mu\text{m}$, and $\ge 97.5\%$ in viable nests ($R \le 200\,\mu\text{m}$). Consequently, avascular nests do not form viable pharmacological sanctuaries for small-molecule TKIs. For macromolecular ADCs, central depletion below $\text{IC}_{50}$ requires $R > 180\,\mu\text{m}$, where central hypoxia and necrosis typically limit cellular viability.

### 5.2 Clinical Cohort Interpretation
The MSKCC cohort demonstrates that the primary architectural prognostic signal is **lepidic protection** ($HR = 0.25$), rather than an escalating risk across compact subtypes. Solid and papillary patterns have similar hazard ratios to acinar, and micropapillary is not worse than acinar ($HR = 0.74$). The exploratory contrast of $HR = 6.75$ reflects the exceptionally favorable baseline biology of lepidic adenocarcinoma rather than sanctuary-driven failure.

### 5.3 Testable Hypotheses
Rather than prescribing clinical management, our findings establish three testable hypotheses for future physical and translational research:
1. **Spatial Mass Spectrometry:** MALDI mass-spectrometry imaging of resected perivascular cords from treated patients should demonstrate radial osimertinib gradients dropping below $\text{IC}_{50}$ beyond $165\,\mu\text{m}$, with flat profiles across avascular nests.
2. **Pharmacodynamic Marker Retention:** Tumors with prominent perivascular cord architecture should retain phospho-EGFR signaling specifically in cord peripheries ($r > 165.9\,\mu\text{m}$).
3. **Histology-Stratified Adjuvant Trials:** In trials evaluating extended adjuvant osimertinib, post-cessation recurrence rates should correlate with the quantified proportion of perivascular cord architecture rather than compact avascular nests.

### 5.4 Limitations
All transport parameters are assumed from biophysical measurements in neoplastic tissue. The model idealizes microvasculature as isolated cylinders without inter-vessel field overlap or convection. Additionally, while plasma protein binding is parameterized at a standard baseline of ~95% ($C_0 \approx 40\,\text{nM}$), physiological inter-patient pharmacokinetic variability and albumin fluctuations can dynamically modulate the free fraction, shifting the threshold $\theta$ and either expanding or compressing the viable sanctuary shell width $\Delta r$ as quantified in our parametric sensitivity analysis. The MSKCC cohort is retrospective, surgical, lacks osimertinib exposure, provides only predominant subtype, and does not record adjuvant systemic chemotherapy.

---

## Captions for Figures

**Figure 1. Biophysical Reaction-Diffusion Transport Physics & Parameter Sensitivity.**  
(A) Outward perivascular cord geometry: radial concentration profiles for small-molecule osimertinib ($\lambda = 627.6\,\mu\text{m}$, blue line) and macromolecular ADC ($\lambda = 77.5\,\mu\text{m}$, red line) derived from modified Bessel function $K_0$. Dotted lines indicate the critical depth $r_{\text{crit}} = 165.9\,\mu\text{m}$ where osimertinib drops below the $\text{IC}_{50}$ threshold ($\theta = 0.350$). The shaded orange region illustrates the viable pharmacological sanctuary shell between $r_{\text{crit}}$ and the tissue oxygenation boundary ($r_{\text{ox}} = 200\,\mu\text{m}$, $\Delta r = 34.1\,\mu\text{m}$; $\Delta r = 14.1\,\mu\text{m}$ at $r_{\text{ox}} = 180\,\mu\text{m}$).  
(B) Inward avascular nest geometry: central core concentration $C_{\text{in}}(0) / C_0 = 1 / I_0(R / \lambda)$ as a function of nest outer radius $R$ up to $400\,\mu\text{m}$. Osimertinib remains at $\ge 90.6\%$ throughout, and $\ge 97.5\%$ within viable nests ($R \le 200\,\mu\text{m}$). ADCs drop below $\theta = 0.350$ only for $R > 180\,\mu\text{m}$. Shaded grey band indicates the critical central hypoxia/necrosis boundary ($R_{\text{crit,ox}} \approx 150\text{--}190\,\mu\text{m}$).  
(C) 2D parameter sensitivity heatmap of sanctuary shell width $\Delta r$ across tissue diffusivity $D$ and cellular uptake rate constant $k$.  
(D) Parameter sensitivity of shell width $\Delta r$ as a function of the tissue hypoxia threshold $r_{\text{ox}}$ and inhibitory threshold $\theta$.

**Figure 2. Clonal Regrowth Dynamics & The Analytical Lead-Time Identity.**  
(A) Simulated post-cessation tumor burden regrowth from residual inoculum $N_0 = 2000$ cells (sensitive doubling time $T_s = 25.0\,\text{days}$, persister doubling time $T_p = 120.0\,\text{days}$). The proliferating sensitive clone rapidly dominates, crossing the ctDNA detection threshold ($10^7$ cells) at $324.8\,\text{days}$ (10.67 months) and the CT imaging detection threshold ($10^9$ cells) at $490.9\,\text{days}$ (16.13 months).  
(B) The diagnostic lead-time identity: diagnostic lead time $\Delta t = T_s \cdot \log_2(100) = 6.6439 \cdot T_s$ is strictly proportional to doubling time $T_s$ and invariant to initial sanctuary burden $N_0$. Markers indicate the baseline model ($T_s = 25.0\,\text{days} \to \Delta t = 5.46\,\text{months}$) and the Herbst et al. (2025) median ($T_s \approx 21.5\,\text{days} \to \Delta t = 4.70\,\text{months}$).  
(C) Invariance of diagnostic lead time $\Delta t = 5.46\,\text{months}$ across initial residual burdens $N_0 \in [10^1, 10^5]$ cells, showing ctDNA detection $t_1$ (blue curve) and CT imaging detection $t_2$ (red curve).

**Figure 3. MSKCC Clinical Cohort Survival Analysis ($N = 604$).**  
(A) Kaplan-Meier recurrence-free survival curves across the 5 IASLC predominant histologic subtypes (reverse Kaplan-Meier median follow-up 23.8 months, IQR: 13.3–31.8 months; unadjusted 4-degrees-of-freedom log-rank $\chi^2 = 24.52, p = 6.27 \times 10^{-5}$).  
(B) Pathologic stage composition across histologic subtypes, showing the marked concentration of Stage I disease in lepidic adenocarcinoma (88.6%) compared to solid (57.4%) and micropapillary (51.4%) subtypes.  
(C) Forest plot of multivariable stage-adjusted Cox proportional hazards models, displaying pairwise contrasts against lepidic reference, models with Acinar reference, composite architectural models, and Firth penalized shrinkage estimates.

---

## Declarations

**Ethics Approval:** Secondary analysis of public, de-identified genomic and clinical data from cBioPortal and Caso et al. (2020) [3]. No additional ethical approval was required.  
**Consent for Publication:** Not applicable.  
**Availability of Data and Code:** The clinical dataset (`06_DATI_CLINICI_MSKCC_604_PAZIENTI.csv`) and Python verification scripts are publicly available in the GitHub repository: https://github.com/Lukecele/spatial-oncology-master and archived on Zenodo at https://doi.org/10.5281/zenodo.23623653.  
**Competing Interests:** The author declares no competing interests.  
**Funding:** No external funding was received for this study.  
**Authors' Contributions:** Luca Celebrano conceived the study, formulated the mathematical models, conducted the clinical data analysis, implemented the software pipeline, and authored the manuscript.  
**AI/LLM Declaration:** In compliance with Springer Nature and ICMJE guidelines, generative AI tools (Anthropic Claude and Google Antigravity) were utilized for code debugging, numerical cross-verification, typographic formatting, and critical editorial auditing under the direct supervision and responsibility of the author.

---

## References

1. Wu YL, Tsuboi M, He J, et al. Osimertinib in Resected *EGFR*-Mutated Non-Small-Cell Lung Cancer. *N Engl J Med*. 2020;383(18):1711-1723. doi:10.1056/NEJMoa2027071.
2. Tsuboi M, Herbst RS, John T, et al. Overall Survival with Osimertinib in Resected *EGFR*-Mutated NSCLC. *N Engl J Med*. 2023;389(2):137-147. doi:10.1056/NEJMoa2304594.
3. Caso R, Sanchez-Vega F, Tan KS, et al. The Underlying Tumor Genomics of Predominant Histologic Subtypes in Lung Adenocarcinoma. *J Thorac Oncol*. 2020;15(12):1844-1856. doi:10.1016/j.jtho.2020.08.005.
4. Thomlinson RH, Gray LH. The histological structure of some human lung cancers and the possible implications for radiotherapy. *Br J Cancer*. 1955;9(4):539-549. doi:10.1038/bjc.1955.55.
5. Cross DAE, Ashton SE, Ghiorghiu S, et al. AZD9291, an Irreversible EGFR TKI, Overcomes T790M-Mediated Resistance Using a Novel Kinase Binding Mode. *Cancer Discov*. 2014;4(9):1046-1061. doi:10.1158/2159-8290.CD-14-0337.
6. Sharma SV, Lee DY, Li B, et al. A chromatin-mediated reversible drug-tolerant state in cancer cell subpopulations. *Cell*. 2010;141(1):69-80. doi:10.1016/j.cell.2010.02.027.
7. Kurppa KJ, Liu Y, To C, et al. Treatment-Induced Tumor Dormancy through YAP-Mediated Transcriptional Reprogramming of the Apoptotic Pathway. *Cancer Cell*. 2020;37(1):104-122.e12. doi:10.1016/j.ccell.2019.12.006.
8. Kaznatcheev A, Peacock J, Basanta D, Marusyk A, Scott JG. Fibroblasts and alectinib switch the evolutionary games played by non-small cell lung cancer. *Nat Ecol Evol*. 2019;3(3):450-456. doi:10.1038/s41559-018-0768-z.
9. Freireich EJ, Gehan E, Frei E 3rd, et al. The effect of 6-mercaptopurine on the duration of steroid-induced remissions in acute leukemia: a model for evaluation of other potentially useful therapy. *Blood*. 1963;21(6):699-716.
10. Gehan EA. A generalized Wilcoxon test for comparing arbitrarily singly-censored samples. *Biometrika*. 1965;52(1/2):203-223. doi:10.1093/biomet/52.1-2.203.
11. Schemper M, Smith TL. A note on quantifying follow-up in studies of failure time. *Control Clin Trials*. 1996;17(4):343-346. doi:10.1016/0197-2456(96)00075-x.
12. Herbst RS, John T, Grohé C, et al. Molecular residual disease analysis of adjuvant osimertinib in resected EGFR-mutated stage IB-IIIA non-small-cell lung cancer. *Nat Med*. 2025;31(6):1958-1968. doi:10.1038/s41591-025-03577-y.
13. Swabb EA, Wei J, Gullino PM. Diffusion and convection in normal and neoplastic tissues. *Cancer Res*. 1974;34(10):2814-2822.
14. Nugent LJ, Jain RK. Extravascular diffusion in normal and neoplastic tissues. *Cancer Res*. 1984;44(1):238-244.
15. Worth KR, Papandreou I, Hammond EM. How the histological structure of some lung cancers shaped almost 70 years of radiobiology. *Br J Cancer*. 2023;128(3):407-412. doi:10.1038/s41416-022-02041-9.
16. Hirst DG, Hirst VK, Joiner B, et al. Changes in tumour morphology with alterations in oxygen availability: further evidence for oxygen as a limiting substrate. *Br J Cancer*. 1991;64(1):54-58. doi:10.1038/bjc.1991.238.
17. Moreira AL, Ocampo PSS, Xia Y, et al. A Grading System for Invasive Pulmonary Adenocarcinoma: A Proposal From the International Association for the Study of Lung Cancer Pathology Committee. *J Thorac Oncol*. 2020;15(10):1599-1610. doi:10.1016/j.jtho.2020.06.001.
18. Holm S. A simple sequentially rejective multiple test procedure. *Scand J Statist*. 1979;6(2):65-70.
