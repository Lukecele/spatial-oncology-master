# Biophysical Reaction-Diffusion Transport and Evolutionary Sanctuary Dynamics in Compact EGFR-Mutant Lung Adenocarcinoma

**Luca Celebrano**$^{1,2,*}$

$^{1}$Department of Clinical Medicine and Surgery, School of Medicine, University of Naples Federico II, 80131 Naples, Italy  
$^{2}$Telethon Institute of Genetics and Medicine (TIGEM), 80078 Pozzuoli, Naples, Italy  
$^*$Correspondence: `l.celebrano@studenti.unina.it` | ORCID: [0009-0007-0931-3555](https://orcid.org/0009-0007-0931-3555)  
**Article Classification:** Translational Oncology, Spatial Biophysics, Mathematical Modeling & Molecular Pathology  
**Master Repository:** [https://github.com/Lukecele/spatial-oncology-master](https://github.com/Lukecele/spatial-oncology-master)  

---

## Abstract

**Background:** Third-generation epidermal growth factor receptor (EGFR) tyrosine kinase inhibitors (TKIs), notably osimertinib, constitute the global standard of care for resected and advanced *EGFR*-mutated non-small cell lung cancer (NSCLC). However, clinical outcomes are marked by substantial heterogeneity: in exploratory molecular residual disease (MRD) analyses of the adjuvant ADAURA trial, up to 68% of disease-free survival or circulating tumor DNA (ctDNA) detection events manifested *after* treatment discontinuation (months 36 to 60), frequently without acquiring canonical secondary resistance mutations (e.g., *EGFR* C797S or *MET* amplification). While systemic pharmacokinetics indicate extensive organ distribution (volume of distribution $V_d \approx 918\,\text{L}$), bulk tissue assays conceal microscopic transport barriers created by neoplastic tissue compaction.

**Methods:** We established a multi-scale biophysical framework integrating: (1) exact cylindrical reaction-diffusion transport physics using the modified Bessel equation of order zero ($K_0$) across both Krogh perivascular cords and avascular solid nests; (2) an evolutionary replicator dynamic modeling subclonal frequency-dependent competition between proliferating sensitive and drug-tolerant persister (DTP) cells, coupled to a calibrated post-cessation regrowth simulation; and (3) an independent clinical cohort audit of 604 surgically resected lung adenocarcinoma patients from Memorial Sloan Kettering Cancer Center (MSKCC 2020) with pathologically classified IASLC architectural growth patterns.

**Results:** Exact Bessel solutions demonstrate that small-molecule osimertinib (500 Da, tissue diffusivity $D = 1.3 \times 10^{-6}\,\text{cm}^2/\text{s}$, metabolic clearance $k = 3.3 \times 10^{-4}\,\text{s}^{-1}$, penetration length scale $\lambda = 627.7\,\mu\text{m}$) falls below the cellular $IC_{50}$ threshold ($8.8\,\text{nM} \approx 0.35 C_0$) at an exact radial depth of $r_{\text{crit}} = 165.9\,\mu\text{m}$. Because neoplastic cords remain metabolically viable within Krogh's physiological oxygen diffusion boundary ($r_{\text{ox}} \approx 200\,\mu\text{m}$), avascular solid and micropapillary nests form a continuous **viable pharmacological sanctuary shell ($\Delta r = 34.1\,\mu\text{m}$)**. Macromolecular therapeutics (ADCs/bispecifics, 150 kDa) suffer severe transport restriction ($r_{\text{crit}} = 47.4\,\mu\text{m}$), leaving $>76\%$ of the viable tumor cord uninhibited. Replicator game dynamics confirm that the sanctuary niche stabilizes an internal polymorphic Nash equilibrium ($x^* = 0.6133$, Lyapunov eigenvalue $\lambda_{\text{eig}} = -1.7934 < 0$). Calibrated population dynamic simulations ($N_0 = 2000$ cells, $T_{\text{double}} = 25$ days sensitive, $120$ days persister) project molecular ctDNA detection at 10.7 months and radiographic CT relapse at 16.1 months post-cessation (lead time $\Delta t = 5.4$ months), consistent with clinical timescales documented in the prospective ADAURA ctDNA study (median lead time 4.7 months, 95% CI: 2.2–5.6). In the authentic MSKCC cohort ($N=604$), observed recurrence rates varied across histological patterns: Lepidic (4.55%) $\to$ Acinar (14.67%) $\to$ Micropapillary (21.62%) $\to$ Papillary (25.58%) $\to$ Solid (30.88%). In multivariable Cox modeling adjusting for pathologic stage, an a priori extreme biophysical contrast comparing high-sanctuary solid/micropapillary architecture directly against lepidic baseline ($N=193$) showed marked separation ($HR = 6.75$, 95% CI: $2.28 - 19.97, p = 5.52 \times 10^{-4}$). In the full cohort ($N=604$), stage-adjusted risk yielded $HR = 1.32$ (95% CI: $0.83 - 2.08, p = 0.238$; $HR = 1.30, p = 0.258$ including *EGFR* status), with 5-subtype multivariable analysis demonstrating that separation is predominantly driven by the protective effect of lepidic histology ($HR = 0.25, 95\%\text{ CI: } 0.09 - 0.70, p = 0.008$ relative to acinar).

**Conclusions:** Spatial architectural compaction establishes a biophysical sanctuary that fosters non-mutational persister niches. These findings provide a biophysical rationale for post-ADAURA recurrence kinetics, highlight the strong prognostic protection of uncompacted lepidic histology, and offer a quantitative framework for risk-stratified post-cessation ctDNA surveillance and rational sequence combinations with macromolecular modalities.

---

## 1. Introduction: The Spatial Paradigm in Targeted Therapy Relapse

Targeted inhibition of oncogenic driver kinases has transformed therapeutic paradigms across medical oncology. In advanced non-small cell lung cancer (NSCLC) harboring sensitizing *EGFR* mutations (exon 19 in-frame deletions or exon 21 L858R substitutions), third-generation EGFR tyrosine kinase inhibitors (TKIs), notably osimertinib, have established definitive superiority. In the adjuvant setting, the Phase III ADAURA trial demonstrated unprecedented disease-free and overall survival benefit for 3 years of daily osimertinib following complete surgical resection (R0) in stage IB–IIIA disease (Wu et al., *NEJM* 2020 [1]; Tsuboi et al., *NEJM* 2023 [2]). 

However, clinical outcomes following the completion of adjuvant treatment exhibit substantial heterogeneity: upon completing the 3-year treatment course, disease-free survival curves begin to drop. In prospective ctDNA monitoring within the ADAURA trial, 19 of 28 (67.9% $\approx 68\%$) disease-free survival or molecular residual disease (MRD) events occurred *after* treatment discontinuation (months 36 to 60), peaking between 12 and 18 months post-cessation (Herbst et al., *Nature Medicine* 2025 [12]).

Crucially, re-biopsy and liquid biopsy analyses reveal that late recurrences following adjuvant discontinuation frequently lack canonical on-target genomic resistance mutations (such as *EGFR* C797S) or off-target bypass activations (such as *MET* amplification). Instead, these relapses frequently preserve clonal sensitivity to EGFR inhibition, demonstrating marked objective responses to osimertinib rechallenge. This phenomenon indicates that post-discontinuation relapses do not primarily emerge from de novo genetic mutations selected under high drug exposure, but rather from persistent, non-mutational subclonal populations that survive three years of therapy in an uninhibited or sub-therapeutic spatial state.

While bulk pharmacokinetic assays demonstrate extensive systemic distribution of osimertinib (apparent volume of distribution $V_d/F \approx 918\,\text{L}$, human plasma protein binding 94.7%, steady-state free vascular concentration $C_{0,\text{free}} \approx 25.0\,\text{nM}$), bulk assays assume spatial homogeneity. Inside solid and micropapillary lung adenocarcinoma nests, the loss of alveolar architecture and elevated interstitial fluid pressure shut down microvascular perfusion, forcing drugs to penetrate exclusively via passive interstitial diffusion. Here, we present a unified biophysical, evolutionary, and clinical investigation examining the spatial sanctuary mechanism underlying non-mutational persistence and post-cessation recurrence in human lung adenocarcinoma.

---

## 2. Exact Biophysical Reaction-Diffusion Transport & Sanctuary Architecture

Neoplastic drug delivery in avascular cellular aggregates is governed by the competition between passive Fickian diffusion through the extracellular matrix and first-order metabolic uptake or intracellular sequestration. In cylindrical coordinates centered on a microvascular source (capillary radius $r_0 = 10.0\,\mu\text{m}$), the steady-state Helmholtz reaction-diffusion equation is:

$$\frac{1}{r} \frac{d}{dr} \left( r \frac{dC}{dr} \right) - \frac{k}{D} C = 0$$

where $D$ is the effective tissue diffusion coefficient ($\text{cm}^2/\text{s}$), $k$ is the volumetric cellular clearance/uptake rate ($\text{s}^{-1}$), and $\lambda = \sqrt{D/k}$ defines the characteristic penetration length scale.

### 2.1 Closed-Form Solution via Modified Bessel Functions
Substituting $x = r/\lambda$, the equation transforms into the canonical modified Bessel differential equation of order zero ($n=0$): $\frac{d^2C}{dx^2} + \frac{1}{x}\frac{dC}{dx} - C = 0$. The general solution is:

$$C(r) = A \cdot I_0(r/\lambda) + B \cdot K_0(r/\lambda)$$

As $r \to \infty$, physical bounding demands $C(r) \to 0$; because $I_0(x)$ diverges exponentially ($I_0(x) \sim e^x/\sqrt{2\pi x}$), the coefficient $A$ must identically equal 0. Applying the Dirichlet boundary condition at the capillary wall $C(r_0) = C_0$, we obtain the exact closed-form solution:

$$C(r) = C_0 \frac{K_0(r / \lambda)}{K_0(r_0 / \lambda)}$$

where $K_0$ is the modified Bessel function of the second kind of order zero. In earlier heuristic literature, approximate asymptotic forms ($C(r) \approx C_0 \sqrt{r_0/r} e^{-(r-r_0)/\lambda}$) were occasionally employed; however, because the relevant microenvironmental regime spans $x = r/\lambda \in [0.016, 0.5]$, asymptotic approximations incur substantial errors. Our framework strictly solves the exact Bessel equation.

From an anatomical and pathological standpoint, this cylindrical formulation captures the **Krogh perivascular tumor cord / tissue sleeve** (characteristic of perivascular cuffs and papillary fibrovascular cores), where diffusion of oxygen and small molecules proceeds **centrifugally** outward from a central functional microvessel ($r_0 = 10.0\,\mu\text{m}$). In this outward geometry, neoplastic cells proximal to the capillary ($r < 165.9\,\mu\text{m}$) experience therapeutic inhibition; cells in the intermediate annulus ($165.9 \le r \le 200.0\,\mu\text{m}$) form a viable pharmacological sanctuary; and tissue beyond the oxygen diffusion limit ($r > 200.0\,\mu\text{m}$) undergoes ischemic coagulative necrosis at the outer inter-vascular margin. In complementary avascular solid nests lacking internal microvessels, diffusion proceeds **centripetally** from peripheral stromal capillaries inward toward the center; under both topological configurations, any neoplastic cells situated at a physical distance $\Delta r \in [165.9, 200.0]\,\mu\text{m}$ from the nearest functional microvessel occupy an identical viable sanctuary band.

### 2.2 Pharmacological Calibration & Viable Sanctuary Shell
For osimertinib at the approved standard daily dose of 80 mg orally, clinical steady-state plasma concentrations average $C_{\text{ss}} \approx 250\text{--}300\,\text{ng/mL}$ ($\approx 500\,\text{nM}$). Given a human plasma protein binding of 94.7% ($f_u = 5.3\%$), the biologically available unbound vascular concentration is $C_{0,\text{free}} \approx 25.0\,\text{nM}$. In human *EGFR*-mutant lung adenocarcinoma cell lines (PC9, H1975), the half-maximal inhibitory concentration for EGFR phosphorylation and clonal suppression is $IC_{50} \approx 8.8\,\text{nM}$ (Cross et al., *Cancer Discovery* 2014 [5]). The therapeutic threshold ratio is therefore defined by physiology:

$$\theta = \frac{IC_{50}}{C_{0,\text{free}}} = \frac{8.8\,\text{nM}}{25.0\,\text{nM}} = 0.352 \approx 0.35$$

At $\theta = 0.350$, solving $C(r_{\text{crit}}) = 0.350 C_0$ yields the exact critical radius: **$r_{\text{crit}} = 165.90\,\mu\text{m}$** (with $\theta = 0.352$, $r_{\text{crit}} = 164.40\,\mu\text{m}$). Neoplastic viability is concurrently constrained by oxygen diffusion. As first quantified in human lung carcinoma cords by Thomlinson & Gray (*Br J Cancer* 1955 [4]), mitochondrial consumption limits oxygen survival to $150\text{--}200\,\mu\text{m}$ from the nearest functional microvessel. Beyond $r_{\text{ox}} = 200.0\,\mu\text{m}$, severe anoxia triggers coagulative necrosis. Therefore, in compact nests where $R \ge 200\,\mu\text{m}$, an anatomical annulus emerges:

$$\Delta r_{\text{sanctuary}} = r_{\text{ox}} - r_{\text{crit}} = 200.0\,\mu\text{m} - 165.90\,\mu\text{m} = \mathbf{34.10\,\mu\text{m}}$$

(or $\Delta r = 35.60\,\mu\text{m}$ for $\theta = 0.352$). Within this shell, cells receive sufficient oxygen to sustain mitochondrial respiration ($pO_2 > 5\,\text{mmHg}$), but drug concentrations fall below $IC_{50}$, establishing a protected metabolic sanctuary.

For macromolecular antibody-drug conjugates (ADCs, e.g. patritumab deruxtecan 150 kDa, $5.6\,\text{mg/kg}$ Q3W) and bispecifics (amivantamab 148 kDa), hindered diffusion ($D = 1.5 \times 10^{-8}\,\text{cm}^2/\text{s}$) and rapid endocytic clearance ($k = 2.5 \times 10^{-4}\,\text{s}^{-1}$) yield $\lambda = 77.5\,\mu\text{m}$. The critical depth contracts to $r_{\text{crit}} = 47.42\,\mu\text{m}$ (or $47.16\,\mu\text{m}$ at $\theta = 0.352$). Consequently, the uninhibited viable shell spans $\Delta r = 152.58\,\mu\text{m}$, leaving over 76% of the viable tumor cord uninhibited (**Figure 1**).

---

## 3. Evolutionary Game Theory & Calibrated ADAURA Post-Cessation Regrowth

To characterize subclonal competition within the pharmacological sanctuary shell, we formulated an evolutionary game between proliferating drug-sensitive cells ($S$, frequency $x$) and drug-tolerant persister cells ($P$, frequency $1-x$). Persister cells exhibit slow-cycling dormancy ($T_{\text{double}} \approx 120\,\text{days}$ in vivo) induced by moderate hypoxia and sub-therapeutic drug stress, mediated by epigenetic remodeling (Sharma et al., *Cell* 2010 [6]; Kurppa et al., *Cancer Cell* 2020 [7]).

### 3.1 Replicator Dynamics & Asymptotic Lyapunov Stability
Under non-zero-sum frequency-dependent interactions calibrated for NSCLC solid architectures, the $2 \times 2$ payoff matrix is:

$$A = \begin{pmatrix} a_{11} & a_{12} \\ a_{21} & a_{22} \end{pmatrix} = \begin{pmatrix} -0.736 & 2.866 \\ 2.188 & -1.772 \end{pmatrix}$$

The standard replicator equation governing the sensitive fraction $x$ is:

$$\frac{dx}{dt} = x (1 - x) [f_S(x) - f_P(x)] = x (1 - x) [\alpha x + \beta]$$

where $\alpha = a_{11} - a_{12} - a_{21} + a_{22} = -7.5620$ and $\beta = a_{12} - a_{22} = 4.6380$. Setting $dx/dt = 0$ yields the internal polymorphic Nash equilibrium:

$$x^* = -\frac{\beta}{\alpha} = \frac{4.6380}{7.5620} = \mathbf{0.613330} \quad (61.33\%)$$

Linearizing the system at $x^*$, the Jacobian eigenvalue is:

$$\lambda_{\text{eig}} = \left. \frac{d}{dx} \left( \frac{dx}{dt} \right) \right|_{x=x^*} = x^* (1 - x^*) \alpha = 0.613330 \times 0.386670 \times (-7.5620) = \mathbf{-1.793376}$$

Because $\lambda_{\text{eig}} < 0$, Lyapunov's first stability theorem guarantees that $x^* = 0.6133$ is an **asymptotically stable attractor (sink)**. During 3 years of adjuvant osimertinib exposure, the subclonal population within the sanctuary shell is maintained in this stable coexistence, preventing clonal extinction while restraining macroscopic expansion.

### 3.2 Calibrated Regrowth Kinetics Post-Adjuvant Discontinuation
Upon treatment cessation at month 36, drug clears from the sanctuary shell. Relieved from pharmacological suppression, the sensitive fraction ($61.3\%$) resumes rapid exponential proliferation ($g_S = \ln(2)/25\,\text{day}^{-1}$ net in vivo doubling), while the persister subclone cycles slowly ($g_P = \ln(2)/120\,\text{day}^{-1}$). Total cell burden expands from an initial residual micro-nest $N_0 = 2,000$ cells according to:

$$N(t) = N_0 \left[ x^* e^{g_S t} + (1 - x^*) e^{g_P t} \right]$$

The model projects:

- **Molecular ctDNA Detection Threshold ($10^7$ cells, ~0.01% VAF):** reached at day 325 post-cessation = **10.7 months** (10.68 months).
- **Radiographic CT Detection Threshold ($10^9$ cells, 1 cm lesion):** reached at day 491 post-cessation = **16.1 months** (16.13 months).
- **Diagnostic Lead Time Window:** $\Delta t = 16.13 - 10.68 = \mathbf{5.45 \approx 5.4\,\text{months}}$ (rounded to 5.5 months in preliminary summaries).

We explicitly emphasize that this post-cessation regrowth simulation represents an exploratory mechanistic model calibration—parameterized by physiologically grounded cell doubling times and detection thresholds—designed to explore whether biophysical and population-dynamic mechanisms yield timescales consistent with clinical reality, rather than an independent clinical trial validation. Notably, these calibrated projections exhibit qualitative and quantitative concordance with empirical findings from the ADAURA ctDNA longitudinal surveillance study (Herbst et al., *Nature Medicine* 2025 [12]), which documented a median ctDNA lead time of 4.7 months (95% CI: 2.2–5.6 months) preceding imaging events (**Figure 2**).

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

Across the entire cohort ($N=604$), the median follow-up for RFS was **21.4 months** (interquartile range: 11.2 to 31.9 months). A total of 77 patients remained at risk at 36 months, and only 3 patients had documented follow-up beyond 60 months (including a single patient with lepidic histology). For the solid subtype, the Kaplan-Meier median RFS estimate of 40.2 months is supported by only 3 patients remaining at risk beyond 40 months. Consequently, long-term survival benchmarks (such as 5-year RFS) cannot be reliably extrapolated from this dataset, and we restrict our conclusions strictly to observed events and early-to-midterm recurrence dynamics.

Comparing the extreme biophysical contrast of high-sanctuary solid/micropapillary architecture ($N=105$, 29 recurrences) against the uncompacted lepidic growth pattern ($N=88$, 4 recurrences; total $N=193$), the univariable Log-Rank test demonstrated marked separation: **$\chi^2 = 20.39, p = 6.33 \times 10^{-6}$** (**Figure 3A**).

### 4.2 Multivariable Cox Proportional Hazards Modeling & Methodological Transparency
To test architectural associations while adjusting for clinical covariates, we evaluated three complementary Cox models:

1. **Model 1: Extreme Biophysical Contrast ($N=193$, 33 events):**  
   Comparing high-sanctuary (Solid/Micropapillary) versus uncompacted Lepidic architecture with adjustment for pathologic stage (Stage II/III vs I), the adjusted Hazard Ratio was:
   $$\text{Adjusted } HR = \mathbf{6.75} \quad (95\%\text{ CI: } \mathbf{2.28 - 19.97}), \quad p = \mathbf{5.52 \times 10^{-4}}$$
   This comparison was pre-specified to test the biophysical hypothesis at the extremes of cellular compaction and vascular decoupling.

2. **Model 2: Full Cohort Analysis ($N=604$, 98 events):**  
   When evaluated across the entire heterogeneous cohort ($N=604$), wherein the intermediate acinar subtype constitutes the majority (368/604, 60.9%), the prognostic hazard of solid/micropapillary architecture relative to all other subtypes combined is attenuated:
   - Adjusted for pathologic stage: $\text{Adjusted } HR = \mathbf{1.32} \quad (95\%\text{ CI: } 0.83 - 2.08, p = 0.238)$
   - Adjusted for pathologic stage and *EGFR* status: $\text{Adjusted } HR = \mathbf{1.30} \quad (95\%\text{ CI: } 0.82 - 2.07, p = 0.258)$
   In the full cohort, pathologic stage remains the predominant independent prognostic factor ($HR = 4.18$, 95% CI: $2.74 - 6.38, p = 2.85 \times 10^{-11}$).

3. **Model 3: 5-Subtype Decomposition with Acinar as Reference ($N=604$, 98 events):**  
   To dissect the source of architectural variation across the full cohort, we fitted a multivariable model treating each subtype separately with acinar histology as the reference category ($N=368$):
   - **Lepidic ($N=88$):** $HR = \mathbf{0.25} \quad (95\%\text{ CI: } 0.09 - 0.70), \quad p = \mathbf{0.0084}$
   - **Solid ($N=68$):** $HR = \mathbf{1.48} \quad (95\%\text{ CI: } 0.89 - 2.49), \quad p = \mathbf{0.133}$
   - **Papillary ($N=43$):** $HR = \mathbf{1.55} \quad (95\%\text{ CI: } 0.80 - 3.00), \quad p = \mathbf{0.191}$
   - **Micropapillary ($N=37$):** $HR = \mathbf{0.72} \quad (95\%\text{ CI: } 0.32 - 1.62), \quad p = \mathbf{0.429}$
   This decomposition reveals that the primary architectural determinant of outcome in this surgically resected cohort is the powerful protective effect of the lepidic growth pattern, whereas solid architecture exhibits a trend toward increased hazard that does not reach independent statistical significance when compared directly to acinar adenocarcinoma.

In the confirmed *EGFR*-mutant subpopulation ($N=182$, **Figure 3C**), univariable Kaplan-Meier analysis showed a trend toward separation between high-sanctuary (Solid/Micropapillary, $N=18$, 22.2% recurrence) and low-sanctuary architectures (Lepidic/Acinar/Papillary, $N=164$, 13.4% recurrence; Log-Rank $\chi^2 = 2.53, p = 0.112$; and $\chi^2 = 5.31, p = 0.021$ when comparing extreme strata Solid/Micro vs Lepidic, $N=56$). In multivariable stage-adjusted modeling within this subcohort alone, wide confidence intervals crossing unity were observed ($HR = 1.10$, 95% CI: $0.06 - 19.97, p = 0.947$). This attenuation is attributable to the fact that 78.6% of the *EGFR*-mutant patients presented with early Stage I disease, yielding only 6 recurrence events between extreme strata ($N=56$). We report these subgroup confidence intervals transparently, cautioning against overclaiming stage independence in small resected subgroups with limited event counts.

### 4.3 Study Limitations & Translational Context
Importantly, the MSKCC cohort consists of patients undergoing complete surgical resection who were chemo/TKI-naive at baseline and did not receive adjuvant osimertinib. Therefore, this dataset demonstrates the baseline prognostic value of tissue architecture in early-stage disease, rather than providing direct empirical verification of post-TKI cessation kinetics. Prospective cohorts with longitudinal post-osimertinib surveillance are required to directly evaluate the interaction between architectural sanctuary shells and post-cessation recurrence.

---

## 5. Translational ADAURA Framework & Therapeutic Decision Architecture

Integrating biophysical penetration limits, subclonal game dynamics, and clinical validation yields an actionable clinical roadmap for optimizing post-adjuvant surveillance (**Figure 4**).

### 5.1 Histology-Guided Post-Cessation ctDNA Surveillance
Current clinical guidelines suggest uniform follow-up intervals for all resected stage IB–IIIA patients completing 3 years of adjuvant osimertinib. Our findings suggest that surgical pathology should inform post-cessation monitoring intensity:

- **Low-Sanctuary Pattern (Lepidic / Acinar Predominant):** Associated with favorable baseline prognosis (observed recurrence rate: 4.5%–14.7%). Routine clinical and radiographic follow-up every 6 months is appropriate.
- **High-Sanctuary Pattern (Solid / Micropapillary Components $\ge 20\%$):** Higher intrinsic recurrence risk and potential for residual microscopic sanctuaries justify **intensive longitudinal ctDNA MRD screening every 3 months starting from month 36 through month 60**.

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
