# Biophysical Reaction-Diffusion Transport and Evolutionary Sanctuary Dynamics in Compact EGFR-Mutant Lung Adenocarcinoma

**Luca Celebrano**  
*School of Medicine and Surgery, University of Naples Federico II, Naples, Italy*  
*Correspondence: luca.celebrano5@gmail.com | ORCID: [0009-0006-2584-7546](https://orcid.org/0009-0006-2584-7546)*  
*Permanent Zenodo Master DOI:* [10.5281/zenodo.23623653](https://doi.org/10.5281/zenodo.23623653)  
*GitHub Code and Data Repository:* [https://github.com/Lukecele/spatial-oncology-master](https://github.com/Lukecele/spatial-oncology-master)

---

## Abstract

**Background:** Third-generation epidermal growth factor receptor (EGFR) tyrosine kinase inhibitors (TKIs), notably osimertinib, constitute the global standard of care for resected and advanced *EGFR*-mutated non-small cell lung cancer (NSCLC). However, clinical outcomes are marked by substantial post-treatment heterogeneity: in prospective molecular residual disease (MRD) surveillance within the adjuvant Phase III ADAURA trial, 67.9% (~68%) of disease-free survival or circulating tumor DNA (ctDNA) detection events manifested *after* treatment discontinuation (months 36 to 60), frequently without acquiring canonical secondary resistance mutations (e.g., *EGFR* C797S or *MET* amplification). While systemic pharmacokinetics indicate extensive macroscopic organ distribution (volume of distribution $V_d \approx 918\,\text{L}$), bulk assays conceal microscopic transport barriers created by neoplastic cellular compaction.

**Methods:** We established a multi-scale biophysical and evolutionary framework integrating: (1) exact cylindrical reaction-diffusion transport physics using the modified Bessel equation of order zero ($K_0$) across both Krogh perivascular cords and avascular solid nests, accompanied by multidimensional parameter sensitivity mapping; (2) an evolutionary replicator dynamic modeling subclonal frequency-dependent competition between proliferating sensitive and drug-tolerant persister (DTP) cells, coupled to an analytical post-cessation regrowth kinetics formulation; and (3) an independent clinical cohort audit of 604 surgically resected lung adenocarcinoma patients from Memorial Sloan Kettering Cancer Center (MSKCC 2020) with pathologically classified IASLC architectural growth patterns.

**Results:** Exact Bessel solutions demonstrate that small-molecule osimertinib (500 Da, effective hindered tissue diffusivity $D = 1.3 \times 10^{-6}\,\text{cm}^2/\text{s}$ [Swabb et al. 1974], volumetric cellular uptake/clearance $k = 3.3 \times 10^{-4}\,\text{s}^{-1}$, penetration length scale $\lambda = 627.7\,\mu\text{m}$) falls below the cellular $IC_{50}$ threshold ($8.8\,\text{nM} \approx 0.35 C_{0,\text{free}}$) at an exact radial depth of $r_{\text{crit}} = 165.9\,\mu\text{m}$. Because neoplastic cords remain metabolically viable within Krogh's physiological oxygen diffusion boundary ($r_{\text{ox}} \approx 200\,\mu\text{m}$ [Thomlinson & Gray 1955]), avascular solid and micropapillary nests form a continuous **viable pharmacological sanctuary shell ($\Delta r = 34.1\,\mu\text{m}$)**. Multidimensional sensitivity analysis across $D \in [0.5, 2.5] \times 10^{-6}\,\text{cm}^2/\text{s}$ and $k \in [1.5, 6.0] \times 10^{-4}\,\text{s}^{-1}$ confirms that this sanctuary shell persists across 83.7% of the physiological parameter space at $r_{\text{ox}} = 200\,\mu\text{m}$ (mean width $41.3\,\mu\text{m}$) and 61.3% at $r_{\text{ox}} = 175\,\mu\text{m}$. Macromolecular therapeutics (antibody-drug conjugates [ADCs] / bispecifics, 150 kDa) suffer severe steric transport restriction ($r_{\text{crit}} = 47.4\,\mu\text{m}$), leaving $>94\%$ of the viable tumor cord uninhibited. Replicator game dynamics confirm that the sanctuary niche stabilizes an internal polymorphic Nash equilibrium ($x^* = 0.6133$, Lyapunov eigenvalue $\lambda_{\text{eig}} = -1.7934 < 0$). In post-cessation regrowth modeling ($N_0 = 2000$ cells, $T_{\text{double}} = 25$ days sensitive, $120$ days persister), molecular ctDNA detection occurs at 10.7 months and radiographic CT relapse at 16.1 months post-cessation. Deconstruction of the diagnostic lead time ($\Delta t = 5.46 \approx 5.4\,\text{months}$) demonstrates that it is analytically governed by exponential clonal doubling between assay limits ($\log_2(100) \times T_d$), whereas the clinical timing of initial recurrence is dictated by residual sanctuary burden ($N_0$). In the authentic MSKCC cohort ($N=604$), observed recurrence rates varied across histological patterns: Lepidic (4.55%) $\to$ Acinar (14.67%) $\to$ Micropapillary (21.62%) $\to$ Papillary (25.58%) $\to$ Solid (30.88%). In multivariable Cox modeling adjusting for pathologic stage, an extreme biophysical boundary comparison contrasting high-sanctuary solid/micropapillary architecture directly against uncompacted lepidic monolayer ($N=193$, 33 events) showed marked separation ($HR = 6.75$, 95% CI: $2.28 - 19.97, p = 5.52 \times 10^{-4}$). Across the full heterogeneous cohort ($N=604$, 98 events), stage-adjusted risk yielded $HR = 1.32$ (95% CI: $0.83 - 2.08, p = 0.238$; $HR = 1.30, p = 0.258$ including *EGFR* status), with 5-subtype multivariable analysis demonstrating that separation is predominantly driven by the pronounced protective effect of lepidic histology ($HR = 0.25, 95\%\text{ CI: } 0.09 - 0.70, p = 0.008$ relative to acinar).

**Conclusions:** Spatial architectural compaction establishes a biophysical sanctuary that shelters drug-sensitive clones from systemic drug exposure, fostering non-mutational persister niches. These findings provide a biophysical rationale for post-ADAURA recurrence kinetics, highlight the prognostic protection of uncompacted lepidic histology, and offer a quantitative framework for risk-stratified post-cessation ctDNA surveillance ($\ge 20\%$ high-grade patterns per IASLC 2020 consensus) and rational combination sequences with macromolecular modalities.

---

## 1. Introduction: The Spatial Paradigm in Targeted Therapy Relapse

Targeted inhibition of oncogenic driver kinases has transformed therapeutic paradigms across medical oncology. In advanced non-small cell lung cancer (NSCLC) harboring sensitizing *EGFR* mutations (exon 19 in-frame deletions or exon 21 L858R substitutions), third-generation EGFR tyrosine kinase inhibitors (TKIs), notably osimertinib, have established definitive superiority. In the adjuvant setting, the Phase III ADAURA trial demonstrated unprecedented disease-free and overall survival benefit for 3 years of daily osimertinib following complete surgical resection (R0) in stage IB–IIIA disease (Wu et al., *NEJM* 2020 [1]; Tsuboi et al., *NEJM* 2023 [2]). 

However, clinical outcomes following the completion of adjuvant treatment exhibit substantial heterogeneity: upon completing the 3-year treatment course, disease-free survival curves begin to decline. In prospective ctDNA monitoring within the ADAURA trial, 19 of 28 (67.9% $\approx 68\%$) disease-free survival or molecular residual disease (MRD) events occurred *after* treatment discontinuation (months 36 to 60), peaking between 12 and 18 months post-cessation (Herbst et al., *Nature Medicine* 2025 [12]).

Crucially, re-biopsy and liquid biopsy analyses reveal that late recurrences following adjuvant discontinuation frequently lack canonical on-target genomic resistance mutations (such as *EGFR* C797S) or off-target bypass activations (such as *MET* amplification). Instead, these relapses frequently preserve clonal sensitivity to EGFR inhibition, demonstrating marked objective responses to osimertinib rechallenge. This phenomenon indicates that post-discontinuation relapses do not primarily emerge from de novo genetic mutations selected under high drug exposure, but rather from persistent, non-mutational subclonal populations that survive three years of therapy in an uninhibited or sub-therapeutic spatial state.

While bulk pharmacokinetic assays demonstrate extensive systemic distribution of osimertinib (apparent volume of distribution $V_d/F \approx 918\,\text{L}$, human plasma protein binding 94.7%, steady-state free vascular concentration $C_{0,\text{free}} \approx 25.0\,\text{nM}$), bulk assays assume spatial homogeneity. Inside solid and micropapillary lung adenocarcinoma nests, the loss of alveolar architecture and elevated interstitial fluid pressure shut down microvascular perfusion, forcing drugs to penetrate exclusively via passive interstitial diffusion. Here, we present a unified biophysical, evolutionary, and clinical investigation examining the spatial sanctuary mechanism underlying non-mutational persistence and post-cessation recurrence in human lung adenocarcinoma.

---

## 2. Exact Biophysical Reaction-Diffusion Transport & Sanctuary Architecture

Neoplastic drug delivery in avascular cellular aggregates is governed by the competition between passive Fickian diffusion through the extracellular matrix and first-order metabolic uptake or intracellular sequestration. In cylindrical coordinates centered on a microvascular source (capillary radius $r_0 = 10.0\,\mu\text{m}$), the steady-state Helmholtz reaction-diffusion equation is:

$$\frac{1}{r} \frac{d}{dr} \left( r \frac{dC}{dr} \right) - \frac{k}{D} C = 0$$

where $D$ is the effective tissue diffusion coefficient ($\text{cm}^2/\text{s}$), $k$ is the volumetric first-order cellular uptake and metabolic clearance rate ($\text{s}^{-1}$), and $\lambda = \sqrt{D/k}$ defines the characteristic screening/penetration length scale.

### 2.1 Closed-Form Solution via Modified Bessel Functions
Substituting $x = r/\lambda$, the equation transforms into the canonical modified Bessel differential equation of order zero ($n=0$): $\frac{d^2C}{dx^2} + \frac{1}{x}\frac{dC}{dx} - C = 0$. The general solution is:

$$C(r) = A \cdot I_0(r/\lambda) + B \cdot K_0(r/\lambda)$$

As $r \to \infty$, physical bounding demands $C(r) \to 0$; because $I_0(x)$ diverges exponentially ($I_0(x) \sim e^x/\sqrt{2\pi x}$), the coefficient $A$ must identically equal 0. Applying the Dirichlet boundary condition at the capillary wall $C(r_0) = C_0$, we obtain the exact closed-form solution:

$$C(r) = C_0 \frac{K_0(r / \lambda)}{K_0(r_0 / \lambda)}$$

where $K_0$ is the modified Bessel function of the second kind of order zero. In earlier heuristic literature, approximate asymptotic forms ($C(r) \approx C_0 \sqrt{r_0/r} e^{-(r-r_0)/\lambda}$) were occasionally employed; however, because the relevant microenvironmental regime spans $x = r/\lambda \in [0.016, 0.5]$, asymptotic approximations incur substantial errors. Our framework strictly solves the exact Bessel equation.

From an anatomical and pathological standpoint, this cylindrical formulation captures the **Krogh perivascular tumor cord / tissue sleeve** (characteristic of perivascular cuffs and papillary fibrovascular cores), where diffusion of oxygen and small molecules proceeds **centrifugally** outward from a central functional microvessel ($r_0 = 10.0\,\mu\text{m}$). In this outward geometry, neoplastic cells proximal to the capillary ($r < 165.9\,\mu\text{m}$) experience therapeutic inhibition; cells in the intermediate annulus ($165.9 \le r \le 200.0\,\mu\text{m}$) form a viable pharmacological sanctuary; and tissue beyond the oxygen diffusion limit ($r > 200.0\,\mu\text{m}$) undergoes ischemic coagulative necrosis at the outer inter-vascular margin. In complementary avascular solid nests lacking internal microvessels, diffusion proceeds **centripetally** from peripheral stromal capillaries inward toward the center; under both topological configurations, any neoplastic cells situated at a physical distance $\Delta r \in [165.9, 200.0]\,\mu\text{m}$ from the nearest functional microvessel occupy an identical viable sanctuary band.

### 2.2 Pharmacological Calibration & Viable Sanctuary Shell
In physiological tissue, small-molecule diffusion is constrained by extracellular matrix (ECM) composition, cellular packing fraction, and interstitial tortuosity. Whereas dilute aqueous diffusivity for a 500 Da molecule is $D_{\text{aq}} \approx 4\text{--}6 \times 10^{-6}\,\text{cm}^2/\text{s}$, classical microvascular measurements by Swabb, Wei, and Gullino (*Cancer Res* 1974 [13]) and Nugent & Jain (*Cancer Res* 1984 [14]) establish that the effective tissue diffusivity in dense carcinomas is reduced by a factor of 3 to 4 ($D_{\text{tissue}} \approx 1.0\text{--}1.5 \times 10^{-6}\,\text{cm}^2/\text{s}$). We therefore calibrate $D = 1.3 \times 10^{-6}\,\text{cm}^2/\text{s}$ as the authentic tissue-hindered diffusivity for osimertinib, and $k = 3.3 \times 10^{-4}\,\text{s}^{-1}$ as the volumetric cellular consumption and metabolic uptake rate, yielding a characteristic screening penetration length $\lambda = \sqrt{D/k} = 627.7\,\mu\text{m}$.

For osimertinib at the approved standard daily dose of 80 mg orally, clinical steady-state plasma concentrations average $C_{\text{ss}} \approx 250\text{--}300\,\text{ng/mL}$ ($\approx 500\,\text{nM}$). Given a human plasma protein binding of 94.7% ($f_u = 5.3\%$), the biologically available unbound vascular concentration is $C_{0,\text{free}} \approx 25.0\,\text{nM}$. In human *EGFR*-mutant lung adenocarcinoma cell lines (PC9, H1975), the half-maximal inhibitory concentration for EGFR phosphorylation and clonal suppression is $IC_{50} \approx 8.8\,\text{nM}$ (Cross et al., *Cancer Discovery* 2014 [5]). The therapeutic threshold ratio is therefore defined directly by clinical pharmacokinetics:

$$\theta = \frac{IC_{50}}{C_{0,\text{free}}} = \frac{8.8\,\text{nM}}{25.0\,\text{nM}} = 0.352 \approx 0.35$$

At $\theta = 0.350$, solving $C(r_{\text{crit}}) = 0.350 C_0$ yields the exact critical radius: **$r_{\text{crit}} = 165.90\,\mu\text{m}$** (with $\theta = 0.352$, $r_{\text{crit}} = 164.39\,\mu\text{m}$). Neoplastic viability is concurrently constrained by oxygen diffusion. As first quantified in human lung carcinoma cords by Thomlinson & Gray (*Br J Cancer* 1955 [4]), mitochondrial consumption limits oxygen survival to $150\text{--}200\,\mu\text{m}$ from the nearest functional microvessel (with viable cords exhibiting intact cell sheaths up to $160\text{--}200\,\mu\text{m}$, beyond which coagulative necrosis is invariably observed). Setting the physiological viability ceiling at $r_{\text{ox}} = 200.0\,\mu\text{m}$, in compact nests where $R \ge 200\,\mu\text{m}$, an anatomical annulus emerges:

$$\Delta r_{\text{sanctuary}} = r_{\text{ox}} - r_{\text{crit}} = 200.0\,\mu\text{m} - 165.90\,\mu\text{m} = \mathbf{34.10\,\mu\text{m}}$$

(or $\Delta r = 35.61\,\mu\text{m}$ for $\theta = 0.352$). Within this shell, cells receive sufficient oxygen to sustain mitochondrial respiration ($pO_2 > 5\,\text{mmHg}$), but drug concentrations fall below $IC_{50}$, establishing a protected metabolic sanctuary.

For macromolecular antibody-drug conjugates (ADCs, e.g. patritumab deruxtecan 150 kDa, $5.6\,\text{mg/kg}$ Q3W) and bispecifics (amivantamab 148 kDa), hindered diffusion ($D = 1.5 \times 10^{-8}\,\text{cm}^2/\text{s}$) and rapid endocytic clearance ($k = 2.5 \times 10^{-4}\,\text{s}^{-1}$) yield $\lambda = 77.5\,\mu\text{m}$. The critical depth contracts to $r_{\text{crit}} = 47.43\,\mu\text{m}$ (or $47.16\,\mu\text{m}$ at $\theta = 0.352$). Consequently, the uninhibited viable shell spans $\Delta r = 152.57\,\mu\text{m}$, leaving over 94% of the viable tumor cord uninhibited (**Figure 1**).

### 2.3 Multidimensional Sensitivity Analysis & Phase Boundary
To test the robustness of the biophysical sanctuary shell against biological heterogeneity, we conducted an exhaustive multidimensional parameter sweep over $2,500$ parameter combinations spanning:
- Effective tissue diffusivity: $D \in [0.5, 2.5] \times 10^{-6}\,\text{cm}^2/\text{s}$
- Volumetric clearance rate: $k \in [1.5, 6.0] \times 10^{-4}\,\text{s}^{-1}$
- Inhibitory threshold ratio: $\theta \in [0.25, 0.45]$
- Oxygen diffusion boundary: $r_{\text{ox}} \in [150.0, 200.0]\,\mu\text{m}$

The viable sanctuary shell exists whenever the critical drug penetration depth is strictly less than the tissue oxygenation boundary: $\Delta r = \max(0, r_{\text{ox}} - r_{\text{crit}}) > 0$. At the standard oxygen limit $r_{\text{ox}} = 200.0\,\mu\text{m}$ and $\theta = 0.35$, the sanctuary shell persists across **83.7%** of the entire physiological grid (mean shell width $\Delta r = 41.3\,\mu\text{m}$). At an intermediate oxygen limit $r_{\text{ox}} = 175.0\,\mu\text{m}$, the shell persists across **61.3%** of parameter space (mean width $26.3\,\mu\text{m}$). If the oxygen boundary contracts to $r_{\text{ox}} = 150.0\,\mu\text{m}$ (the lower boundary in Thomlinson & Gray 1955), the shell persists across **29.0%** of the grid (mean width $17.1\,\mu\text{m}$ at $\theta = 0.35$, expanding to **99.3%** at $\theta = 0.45$). 

Analytically, if tissue diffusivity $D$ increases by 30% ($D \to 1.69 \times 10^{-6}\,\text{cm}^2/\text{s}$), $r_{\text{crit}}$ shifts outward to $180.08\,\mu\text{m}$; if cellular clearance $k$ decreases by 30% ($k \to 2.31 \times 10^{-4}\,\text{s}^{-1}$), $r_{\text{crit}}$ shifts to $185.48\,\mu\text{m}$. Only under the concurrent combination of both extremes (+30% $D$ and -30% $k$) does $r_{\text{crit}}$ reach $201.41\,\mu\text{m}$, exceeding $r_{\text{ox}} = 200\,\mu\text{m}$ and collapsing the shell ($\Delta r \to 0$). Conversely, in dense, collagen-rich solid tumors where interstitial resistance reduces diffusivity ($D \to 0.7 \times D_0$) or metabolic uptake increases ($k \to 1.3 \times k_0$), $r_{\text{crit}}$ contracts to $136.85\,\mu\text{m}$, substantially widening the viable sanctuary shell to $\Delta r = 63.15\,\mu\text{m}$.

---

## 3. Evolutionary Game Theory & Post-Cessation Regrowth Kinetics

To characterize subclonal competition within the pharmacological sanctuary shell, we formulated an evolutionary game between proliferating drug-sensitive cells ($S$, frequency $x$) and drug-tolerant persister cells ($P$, frequency $1-x$). Persister cells exhibit slow-cycling dormancy ($T_{\text{double}} \approx 120\,\text{days}$ in vivo) induced by moderate hypoxia and sub-therapeutic drug stress, mediated by epigenetic remodeling (Sharma et al., *Cell* 2010 [6]; Kurppa et al., *Cancer Cell* 2020 [7]).

### 3.1 Replicator Dynamics & Asymptotic Lyapunov Stability
Under non-zero-sum frequency-dependent interactions calibrated for NSCLC solid architectures (following the empirical game assay framework of Kaznatcheev et al., *Nat Ecol Evol* 2019 [15]), the $2 \times 2$ payoff matrix is:

$$A = \begin{pmatrix} a_{11} & a_{12} \\ a_{21} & a_{22} \end{pmatrix} = \begin{pmatrix} -0.736 & 2.866 \\ 2.188 & -1.772 \end{pmatrix}$$

where Strategy 1 = Sensitive ($S$) and Strategy 2 = Persister ($P$). Mutual invasibility criteria are strictly satisfied: $a_{12} > a_{22}$ ($2.866 > -1.772$, sensitive invades persister) and $a_{21} > a_{11}$ ($2.188 > -0.736$, persister invades sensitive), classifying the sanctuary dynamic as a Hawk-Dove / Snowdrift coexistence game.

The standard continuous replicator equation governing the sensitive fraction $x$ is:

$$\frac{dx}{dt} = x (1 - x) [f_S(x) - f_P(x)] = x (1 - x) [\alpha x + \beta]$$

where $\alpha = a_{11} - a_{12} - a_{21} + a_{22} = -7.5620$ and $\beta = a_{12} - a_{22} = 4.6380$. Setting $dx/dt = 0$ yields the internal polymorphic Nash equilibrium:

$$x^* = -\frac{\beta}{\alpha} = \frac{4.6380}{7.5620} = \mathbf{0.613330} \quad (61.33\%)$$

Linearizing the system at $x^*$, the Jacobian eigenvalue is:

$$\lambda_{\text{eig}} = \left. \frac{d}{dx} \left( \frac{dx}{dt} \right) \right|_{x=x^*} = x^* (1 - x^*) \alpha = 0.613330 \times 0.386670 \times (-7.5620) = \mathbf{-1.793376}$$

Because $\lambda_{\text{eig}} < 0$, Lyapunov's first stability theorem guarantees that $x^* = 0.6133$ is an **asymptotically stable attractor (sink)**. During 3 years of adjuvant osimertinib exposure, the subclonal population within the sanctuary shell is maintained in this stable coexistence, preventing clonal extinction while restraining macroscopic expansion.

### 3.2 Kinetic Deconstruction of Recurrence Onset vs Diagnostic Lead Time
Upon treatment cessation at month 36, drug clears from the sanctuary shell. Relieved from pharmacological suppression, the sensitive fraction ($61.3\%$) resumes exponential proliferation ($g_S = \ln(2)/25\,\text{day}^{-1}$ net in vivo doubling, accounting for cell loss), while the persister subclone cycles slowly ($g_P = \ln(2)/120\,\text{day}^{-1}$). Total cell burden expands from an initial residual micro-nest $N_0$ according to:

$$N(t) = N_0 \left[ x^* e^{g_S t} + (1 - x^*) e^{g_P t} \right]$$

Crucially, our mathematical deconstruction reveals that two distinct phenomena govern post-cessation dynamics:

1. **Analytical Invariance of Diagnostic Lead Time ($\Delta t \approx 5.4\,\text{months}$):**  
   The diagnostic window between molecular ctDNA detection ($N_{\text{MRD}} = 10^7$ cells, ~0.01% variant allele fraction) and radiographic CT relapse ($N_{\text{CT}} = 10^9$ cells, corresponding to a 1 cm solid nodule) spans a 100-fold expansion. By the time the micro-nest reaches $10^7$ cells, the sensitive clone ($T_d = 25$ days) outcompetes the persister clone ($T_d = 120$ days) to constitute $>99.9\%$ of the population. Consequently, the lead time is analytically determined by exponential doubling between assay limits:
   $$\Delta t = \log_2\left(\frac{N_{\text{CT}}}{N_{\text{MRD}}}\right) \times T_d = \log_2(100) \times 25\,\text{days} = 6.6439 \times 25 = 166.1\,\text{days} = \mathbf{5.46 \approx 5.4\,\text{months}}$$
   As verified across initial burdens from $N_0 = 200$ to $N_0 = 50,000$ cells, $\Delta t$ remains strictly invariant at $5.45\text{--}5.46$ months, representing an intrinsic property of the cellular doubling rate and assay detection thresholds.

2. **Sanctuary Burden Governs Initial Recurrence Timing:**  
   In contrast to the lead time, the **calendar timing of first molecular detection** depends directly on the residual sanctuary inoculum $N_0$:
   - For $N_0 = 20,000$ residual cells: ctDNA detection occurs at **Month 8.8** post-cessation.
   - For $N_0 = 2,000$ residual cells: ctDNA detection occurs at **Month 10.7** post-cessation (CT relapse at Month 16.1).
   - For $N_0 = 200$ residual cells: ctDNA detection occurs at **Month 13.4** post-cessation (CT relapse at Month 18.9).

We explicitly emphasize that this post-cessation regrowth simulation represents an exploratory mechanistic model calibration—parameterized by physiologically grounded cell doubling times and detection thresholds—designed to explore whether biophysical and population-dynamic mechanisms yield timescales consistent with clinical reality, rather than an independent clinical trial validation. Notably, these calibrated projections exhibit qualitative and quantitative concordance with empirical findings from the ADAURA ctDNA longitudinal surveillance study (Herbst et al., *Nature Medicine* 2025 [12]), which documented a median ctDNA lead time of 4.7 months (95% CI: 2.2–5.6 months) preceding imaging events, with post-cessation recurrences clustering predominantly between months 12 and 18 post-cessation (**Figure 2**).

---

## 4. Clinical Cohort Analysis & Transparent Survival Modeling (MSKCC 2020, $N=604$)

To evaluate histological growth patterns in human disease, we analyzed the authenticated, publicly accessible Memorial Sloan Kettering Cancer Center lung adenocarcinoma cohort (MSKCC 2020, Caso et al., *J Thorac Oncol* 2020 [3]), encompassing 604 consecutive surgically resected patients with pathologically quantified IASLC architectural growth patterns and MSK-IMPACT genomic profiling.

### 4.1 Observed Recurrence Rates Across IASLC Architectural Subtypes
Observed recurrence rates in the MSKCC cohort varied substantially across histological growth patterns (**Figure 3**):

- **Lepidic ($N=88$):** 4 recurrences (**4.55%**), Median RFS not reached.
- **Acinar ($N=368$):** 54 recurrences (**14.67%**), Median RFS not reached.
- **Micropapillary ($N=37$):** 8 recurrences (**21.62%**), Median RFS not reached.
- **Papillary ($N=43$):** 11 recurrences (**25.58%**), Median RFS not reached.
- **Solid ($N=68$):** 21 recurrences (**30.88%**), Kaplan-Meier Median RFS = **40.2 months** (95% CI: 36.5 – not reached).

Across the entire cohort ($N=604$), the median follow-up for RFS was **21.4 months** (interquartile range: 11.6 to 29.1 months). A total of 77 patients remained at risk at 36 months, and only 3 patients had documented follow-up beyond 60 months (including a single patient with lepidic histology). For the solid subtype, the Kaplan-Meier median RFS estimate of 40.2 months is supported by only 3 patients remaining at risk beyond 40 months. Consequently, long-term survival benchmarks (such as 5-year RFS) cannot be reliably extrapolated from this dataset, and we restrict our conclusions strictly to observed events and early-to-midterm recurrence dynamics.

Comparing the extreme biophysical contrast of high-sanctuary solid/micropapillary architecture ($N=105$, 29 recurrences) against the uncompacted lepidic growth pattern ($N=88$, 4 recurrences; total $N=193$), the univariable Log-Rank test demonstrated marked separation: **$\chi^2 = 20.39, p = 6.33 \times 10^{-6}$** (**Figure 3A**).

### 4.2 Multivariable Cox Proportional Hazards Modeling & Methodological Transparency
To test architectural associations while adjusting for clinical covariates, we evaluated three complementary Cox models:

1. **Model 1: Extreme Biophysical Boundary Comparison ($N=193$, 33 events):**  
   To test the biophysical hypothesis at the morphological extremes of cellular compaction (comparing avascular 3D syncytial nests against uncompacted 2D alveolar monolayers) with adjustment for pathologic stage (Stage II/III vs I), the adjusted Hazard Ratio was:
   $$\text{Adjusted } HR = \mathbf{6.75} \quad (95\%\text{ CI: } \mathbf{2.28 - 19.97}), \quad p = \mathbf{5.52 \times 10^{-4}}$$
   While this extreme comparison isolates the biophysical boundary condition, it excludes the intermediate acinar and papillary subtypes (411 patients).

2. **Model 2: Full Cohort Analysis ($N=604$, 98 events):**  
   When evaluated across the entire heterogeneous cohort ($N=604$), wherein the intermediate acinar subtype constitutes the majority (368/604, 60.9%), the prognostic hazard of solid/micropapillary architecture relative to all other subtypes combined is attenuated:
   - Adjusted for pathologic stage: $\text{Adjusted } HR = \mathbf{1.32} \quad (95\%\text{ CI: } 0.83 - 2.08, p = 0.238)$
   - Adjusted for pathologic stage and *EGFR* status: $\text{Adjusted } HR = \mathbf{1.30} \quad (95\%\text{ CI: } 0.82 - 2.07, p = 0.258)$
   In the full cohort, pathologic stage remains the predominant independent prognostic factor ($HR = 4.18$, 95% CI: $2.74 - 6.38, p = 2.85 \times 10^{-11}$).

3. **Model 3: 5-Subtype Decomposition with Acinar as Reference ($N=604$, 98 events):**  
   To dissect the source of architectural variation across the full cohort, we fitted a multivariable model treating each subtype separately with acinar histology as the reference category ($N=368$):
   - **Lepidic ($N=88$):** $HR = \mathbf{0.25} \quad (95\%\text{ CI: } 0.09 - 0.70), \quad p = \mathbf{0.0084}$
   - **Solid ($N=68$):** $HR = \mathbf{1.48} \quad (95\%\text{ CI: } 0.89 - 2.49), \quad p = \mathbf{0.133}$
   - **Papillary ($N=43$):** $HR = \mathbf{1.51} \quad (95\%\text{ CI: } 0.79 - 2.92), \quad p = \mathbf{0.215}$
   - **Micropapillary ($N=37$):** $HR = \mathbf{0.74} \quad (95\%\text{ CI: } 0.33 - 1.66), \quad p = \mathbf{0.468}$
   This decomposition reveals that the primary architectural determinant of outcome in this surgically resected cohort is the powerful protective effect of the lepidic growth pattern, whereas solid architecture exhibits a trend toward increased hazard that does not reach independent statistical significance when compared directly to acinar adenocarcinoma.

In the confirmed *EGFR*-mutant subpopulation ($N=182$, **Figure 3C**), univariable Kaplan-Meier analysis showed a trend toward separation between high-sanctuary (Solid/Micropapillary, $N=18$, 22.2% recurrence) and low-sanctuary architectures (Lepidic/Acinar/Papillary, $N=164$, 13.4% recurrence; Log-Rank $\chi^2 = 2.53, p = 0.112$; and $\chi^2 = 5.31, p = 0.021$ when comparing extreme strata Solid/Micro vs Lepidic, $N=56$). In multivariable stage-adjusted modeling within this subcohort alone, wide confidence intervals crossing unity were observed ($HR = 1.10$, 95% CI: $0.06 - 19.97, p = 0.947$). This attenuation is attributable to the fact that 78.6% of the *EGFR*-mutant patients presented with early Stage I disease, yielding only 6 recurrence events between extreme strata ($N=56$). We report these subgroup confidence intervals transparently, cautioning against overclaiming stage independence in small resected subgroups with limited event counts.

### 4.3 Study Limitations & Translational Context
Importantly, the MSKCC cohort consists of patients undergoing complete surgical resection who were chemo/TKI-naive at baseline and did not receive adjuvant osimertinib. Therefore, this dataset demonstrates the baseline prognostic value of tissue architecture in early-stage disease, rather than providing direct empirical verification of post-TKI cessation kinetics. Prospective cohorts with longitudinal post-osimertinib surveillance are required to directly evaluate the interaction between architectural sanctuary shells and post-cessation recurrence.

---

## 5. Translational ADAURA Framework & Therapeutic Decision Architecture

Integrating biophysical penetration limits, subclonal game dynamics, and clinical validation yields an actionable translational framework for post-adjuvant surveillance (**Figure 4**).

### 5.1 Histology-Guided Post-Cessation ctDNA Surveillance
Current clinical guidelines suggest uniform follow-up intervals for all resected stage IB–IIIA patients completing 3 years of adjuvant osimertinib. Our findings indicate that surgical pathology should inform post-cessation monitoring intensity:

- **Low-Sanctuary Pattern (Lepidic / Acinar Predominant, <20% High-Grade):** Associated with favorable baseline prognosis (observed recurrence rate: 4.5%–14.7%). Routine clinical and radiographic follow-up every 6 months is appropriate.
- **High-Sanctuary Pattern (Solid / Micropapillary Components $\ge 20\%$):** Grounded in the official IASLC 2020 consensus grading system (Moreira et al., *J Thorac Oncol* 2020 [16]), tumors harboring $\ge 20\%$ high-grade patterns define Grade 3 invasive adenocarcinoma. We propose that these high-sanctuary patients represent prime candidates for **quarterly (every 3 months) longitudinal ctDNA MRD surveillance from month 36 through month 60**. Because our kinetic modeling establishes an analytical lead time of $5.4\,\text{months}$, a 3-month sampling frequency satisfies the Nyquist-Shannon sampling requirement, ensuring that emerging molecular relapses are intercepted within the molecular window prior to macroscopic radiographic progression. We emphasize that this proposed schedule represents a theoretical kinetic framework that requires prospective trial validation to verify whether preemptive intervention translates into overall survival gains and to rule out lead-time bias.

### 5.2 Exploiting the Actionable ctDNA Lead Time Window
Because the biophysical sanctuary preserves drug-sensitive founder clones, molecular recurrence detected via ctDNA at ~10.7 months post-cessation represents an expanding population of *EGFR*-sensitive cells preceding macroscopic radiographic tumor formation by **5.4 months** (in our calibrated simulation). Initiating prompt **osimertinib rechallenge** during this molecular-only window (before the lesion expands to $10^9$ cells and undergoes secondary clonal diversification) represents a rational therapeutic strategy to intercept clinical relapse.

### 5.3 Rational Combination Strategies with Macromolecules (ADCs & Bispecifics)
Our transport calculations demonstrate that while small-molecule osimertinib achieves deep penetration ($r_{\text{crit}} = 165.9\,\mu\text{m}$), it leaves a viable peripheral sanctuary shell ($165.9\text{--}200.0\,\mu\text{m}$). Conversely, macromolecular ADCs (patritumab deruxtecan) and bispecific antibodies (amivantamab) deliver potent payload cytotoxicity and antibody-dependent cellular cytotoxicity (ADCC) concentrated intensely within the perivascular corona ($r \le 47.4\,\mu\text{m}$). Mechanistically, targeted cytoreduction of this perivascular cellular mantle by the macromolecule eliminates local metabolic drug consumption ($k \to 0$ across $r \le 47.4\,\mu\text{m}$) and simultaneously decompresses microvascular interstitial fluid pressure ($P_i$). This dynamic clearing expands the effective penetration length of the co-administered small-molecule TKI ($\lambda_{\text{eff}} = \sqrt{D/k_{\text{eff}}}$), shifting the TKI critical threshold ($r_{\text{crit}}$) dynamically beyond $200\,\mu\text{m}$ and collapsing the hypoxic sanctuary shell. This provides a biophysical rationale for the superior progression-free survival observed in the MARIPOSA regimen (amivantamab + lazertinib [10]) and ADC combinations [11].

---

## 6. Discussion

By linking reaction-diffusion transport physics with evolutionary population dynamics and clinical cohort survival data in 604 resected lung cancer patients, this work provides a biophysical foundation for understanding non-genetic drug tolerance. Spatial compaction acts as a physical barrier that shelters drug-sensitive cells from systemic drug levels, explaining how local treatment failure can occur despite extensive macroscopic drug biodistribution.

Importantly, our findings offer a coherent biophysical rationale for post-ADAURA recurrence patterns: the emergence of late relapses without *EGFR* C797S mutations is consistent with an evolutionary replicator dynamic released from boundary suppression upon treatment completion. Incorporating architectural subtyping and transport considerations into post-surgical risk stratification will help clinicians tailor ctDNA surveillance and deploy preemptive therapies before macroscopic disease recurrence manifests.

---

## Declarations

**Funding:** This research received no specific grant from any funding agency in the public, commercial, or not-for-profit sectors.

**Competing Interests:** The author declares no competing financial or non-financial interests.

**Ethics Declaration:** This study exclusively analyzes publicly available, de-identified clinical and genomic data from the Memorial Sloan Kettering Cancer Center cohort (cBioPortal, Caso et al., *J Thorac Oncol* 2020). In accordance with international guidelines, secondary analysis of de-identified public datasets requires no additional institutional review board approval.

**Data and Code Availability:** The authenticated clinical dataset (`06_DATI_CLINICI_MSKCC_604_PAZIENTI.csv`), exact Bessel analytical solvers, evolutionary replicator scripts, and automated verification suite are publicly available at https://github.com/Lukecele/spatial-oncology-master and archived on Zenodo at https://doi.org/10.5281/zenodo.23623653.

---

## References

1. Wu YL, Tsuboi M, He J, et al. Osimertinib in Resected *EGFR*-Mutated Non-Small-Cell Lung Cancer. *N Engl J Med*. 2020;383(18):1711-1723. doi:10.1056/NEJMoa2027071.
2. Tsuboi M, Herbst RS, John T, et al. Overall Survival with Osimertinib in Resected *EGFR*-Mutated NSCLC. *N Engl J Med*. 2023;389(2):137-147. doi:10.1056/NEJMoa2304594.
3. Caso R, Sanchez-Vega F, Tan KS, et al. The Aggregate Risk of Histologic Subtypes Outperforms the Primary Subtype in Predicting Recurrence of Lung Adenocarcinoma. *J Thorac Oncol*. 2020;15(11):1747-1755. doi:10.1016/j.jtho.2020.07.019.
4. Thomlinson RH, Gray LH. The histological structure of some human lung cancers and the possible implications for radiotherapy. *Br J Cancer*. 1955;9(4):539-549.
5. Cross DAE, Ashton SE, Ghiorghiu S, et al. AZD9291, an irreversible EGFR TKI, overcomes T790M-mediated resistance using a novel kinase binding mode. *Cancer Discov*. 2014;4(9):1046-1061.
6. Sharma SV, Lee DY, Li B, et al. A chromatin-mediated reversible drug-tolerant state in cancer cell subpopulations. *Cell*. 2010;141(1):69-80.
7. Kurppa KJ, Liu Y, To C, et al. Treatment-Induced Emergence of a Drug-Tolerant Persister State in *EGFR*-Mutant Lung Cancer. *Cancer Cell*. 2020;37(1):82-97.
8. Boucher Y, Baxter LT, Jain RK. Interstitial fluid pressure gradients in normal and neoplastic tissues. *Cancer Res*. 1990;50(15):4478-4484.
9. Jain RK. Normalization of tumor vasculature: an emerging concept in antiangiogenic therapy. *Science*. 2005;307(5706):58-62.
10. Cho BC, Felip E, Spira AI, et al. Amivantamab plus Lazertinib in First-Line *EGFR*-Mutated Advanced NSCLC (MARIPOSA). *N Engl J Med*. 2024;391(16):1486-1498.
11. Jänne PA, Baik C, Su WC, et al. Efficacy and Safety of Patritumab Deruxtecan (HER3-DXd) in Advanced *EGFR*-Mutated Non-Small Cell Lung Cancer (HERTHENA-Lung01). *J Clin Oncol*. 2024;42(20):2410-2420.
12. Herbst RS, John T, Grohé C, et al. Molecular residual disease analysis of adjuvant osimertinib in resected *EGFR*-mutated stage IB–IIIA non-small-cell lung cancer. *Nat Med*. 2025;31(6):1958-1968. doi:10.1038/s41591-025-03577-y.
13. Swabb EA, Wei J, Gullino PM. Diffusion and convection in normal and neoplastic tissues. *Cancer Res*. 1974;34(10):2814-2822.
14. Nugent LJ, Jain RK. Extravascular diffusion in normal and neoplastic tissues. *Cancer Res*. 1984;44(1):238-244.
15. Kaznatcheev A, Peacock J, Basanta D, Scott JG. Fibroblasts decelerate the evolution of resistance to EGFR tyrosine kinase inhibition in lung adenocarcinoma. *Nat Ecol Evol*. 2019;3(7):1094-1104. doi:10.1038/s41559-019-0909-7.
16. Moreira AL, Ocampo PSS, Xia Y, et al. A Grading System for Invasive Pulmonary Adenocarcinoma: A Proposal From the International Association for the Study of Lung Cancer Pathology Committee. *J Thorac Oncol*. 2020;15(10):1599-1610. doi:10.1016/j.jtho.2020.06.001.
