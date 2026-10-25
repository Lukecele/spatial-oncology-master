# EDITORIAL COVER LETTER

**To:** The Editor-in-Chief  
***Bulletin of Mathematical Biology*** / ***Mathematical Biosciences*** (or ***Lung Cancer***)  

**Date:** October 25, 2026  
**Subject:** Submission of Original Research Article: *"Biophysical Reaction-Diffusion Transport and Evolutionary Sanctuary Dynamics in Compact EGFR-Mutant Lung Adenocarcinoma"*  
**Corresponding Author:** Luca Celebrano (University of Naples Federico II / TIGEM)  

---

Dear Editor,

We are pleased to submit our original research manuscript titled **"Biophysical Reaction-Diffusion Transport and Evolutionary Sanctuary Dynamics in Compact EGFR-Mutant Lung Adenocarcinoma"** for consideration as an Original Research Article in your prestigious journal.

### The Unsolved Clinical and Biological Dilemma
Third-generation epidermal growth factor receptor (EGFR) tyrosine kinase inhibitors (TKIs), exemplified by osimertinib, have established definitive clinical benefit across advanced and surgically resected *EGFR*-mutated non-small cell lung cancer (NSCLC). In the landmark adjuvant ADAURA trial, 3 years of daily osimertinib dramatically improved disease-free and overall survival. However, clinical trajectories following treatment completion are marked by a perplexing paradox: up to 68% of post-surgical disease-free survival or molecular residual disease (MRD) events manifest *after* treatment discontinuation (peaking between 12 and 18 months post-cessation; Herbst et al., *Nature Medicine* 2025). 

Remarkably, genomic sequencing of late recurrences demonstrates that they frequently lack canonical secondary resistance mutations (e.g., *EGFR* C797S or *MET* amplification), retaining marked clonal sensitivity to osimertinib rechallenge. While bulk pharmacokinetic measurements indicate extensive organ biodistribution (apparent volume of distribution $V_d \approx 918\,\text{L}$), macroscopic assays overlook the microscopic transport barriers imposed by solid tissue architecture and convective vascular collapse.

### Core Scientific Contributions of This Work
In this manuscript, we present a multi-scale biophysical and evolutionary framework that bridges fundamental transport physics, population dynamics, and real-world clinical pathology:

1. **Exact Reaction-Diffusion Transport Physics (Modified Bessel $K_0$ Solutions):**  
   We model steady-state interstitial drug transport across both Krogh perivascular cords (centrifugal diffusion) and avascular solid nests (centripetal diffusion). Solving the exact Helmholtz reaction-diffusion equation with modified Bessel functions of the second kind of order zero ($K_0$), we demonstrate that osimertinib (500 Da, $\lambda = 627.7\,\mu\text{m}$) falls below the cellular $\text{IC}_{50}$ threshold ($8.8\,\text{nM} \approx 0.35 C_0$) at an exact radial depth of $r_{\text{crit}} = 165.9\,\mu\text{m}$. Because neoplastic tissue remains metabolically viable within the oxygen diffusion boundary ($r_{\text{ox}} \approx 200\,\mu\text{m}$), avascular nests establish a **viable pharmacological sanctuary shell ($\Delta r = 34.1\,\mu\text{m}$)**. In contrast, macromolecular antibody-drug conjugates (ADCs, 150 kDa) suffer severe transport restriction ($r_{\text{crit}} = 47.4\,\mu\text{m}$), leaving $>76\%$ of the viable tumor cord uninhibited.

2. **Evolutionary Game Theory & Post-Cessation Regrowth Dynamics:**  
   Formulating a non-zero-sum frequency-dependent game between proliferating drug-sensitive cells and drug-tolerant persisters (DTPs), we prove the existence of an asymptotically stable polymorphic attractor ($x^* = 0.6133$, Lyapunov eigenvalue $\lambda_{\text{eig}} = -1.7934 < 0$). Calibrated population dynamic simulations parameterized by physiological doubling times ($T_{\text{double}} = 25$ days sensitive, $120$ days persister; $N_0 = 2000$ cells) predict molecular ctDNA detection at 10.7 months and radiographic CT relapse at 16.1 months (lead time $\Delta t = 5.4$ months), quantitatively aligning with the clinical timescales observed in prospective ADAURA ctDNA monitoring.

3. **Transparent Clinical Cohort Validation (MSKCC 2020, $N = 604$):**  
   We evaluate our biophysical predictions against the authenticated, publicly accessible Memorial Sloan Kettering Cancer Center cohort of 604 surgically resected lung adenocarcinomas (Caso et al., *J Thorac Oncol* 2020). We report:
   - An a priori extreme biophysical contrast comparing high-sanctuary compact architecture (Solid/Micropapillary, $N=105$, 29 events) directly against uncompacted Lepidic tissue ($N=88$, 4 events; total $N=193$), demonstrating strong prognostic separation (Log-Rank $\chi^2 = 20.39, p = 6.33 \times 10^{-6}$; stage-adjusted multivariable $HR = 6.75$, 95% CI: $2.28 - 19.97, p = 5.52 \times 10^{-4}$).
   - The full-cohort multivariable analysis ($N=604, 98$ events), where stage-adjusted risk yields $HR = 1.32$ (95% CI: $0.83 - 2.08, p = 0.238$; $HR = 1.30, p = 0.258$ including *EGFR* status).
   - A 5-subtype multivariable model with Acinar as reference ($N=368$), revealing that the primary driver of outcome in early-stage resected disease is the pronounced protective effect of Lepidic histology ($HR = 0.25, 95\%\text{ CI } 0.09 - 0.70, p = 0.008$), while Solid architecture exhibits a non-significant trend toward increased hazard ($HR = 1.48, p = 0.133$).
   - Rigorous methodological transparency regarding cohort characteristics (predominantly Stage I, median RFS follow-up of 21.4 months, 3 patients at risk beyond 60 months; solid median RFS = 40.2 months), explicitly discussing the distinction between baseline surgical prognostic value and prospective post-TKI surveillance.

4. **Translational Decision Architecture:**  
   We synthesize these findings into an actionable framework for risk-stratified post-cessation ctDNA surveillance and provide a biophysical rationale for combining small molecules with perivascular cytoreductive macromolecules (ADCs/amivantamab) to decompress interstitial pressure, eliminate the perivascular consumption sink ($k \to 0$), and expand small-molecule penetration to collapse the sanctuary shell.

### Commitment to Open Science & Reproducibility
In adherence to the highest standards of scientific reproducibility, the complete pipeline—including the exact analytical Bessel solvers, evolutionary replicator scripts, authenticated MSKCC dataset, unit test suite, and high-resolution 600 DPI vector figures—is fully open-source and publicly accessible in our GitHub repository (`https://github.com/Lukecele/spatial-oncology-master`) and permanently archived on Zenodo (`https://doi.org/10.5281/zenodo.23623653`).

This manuscript represents original work that has not been published previously and is not under consideration elsewhere. The author has approved the submission and declares no competing interests.

Thank you very much for your time, consideration, and editorial review.

Sincerely,

**Luca Celebrano**  
Department of Clinical Medicine and Surgery, School of Medicine, University of Naples Federico II, 80131 Naples, Italy  
Telethon Institute of Genetics and Medicine (TIGEM), 80078 Pozzuoli, Naples, Italy  
Email: `l.celebrano@studenti.unina.it` | ORCID: [0009-0007-0931-3555](https://orcid.org/0009-0007-0931-3555)  
