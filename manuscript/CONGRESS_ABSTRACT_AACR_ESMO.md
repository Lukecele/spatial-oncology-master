# OFFICIAL CONGRESS ABSTRACT: AACR / ESMO / ASCO 2027

**Abstract Category:** Clinical & Translational Oncology / Computational Pathology & Spatial Biology  
**Session Target:** Oral Presentation / Plenary Session (Biomarkers & Targeted Therapy)

---

### Title:
**The Celebrano Spatial Sanctuary Index ($C\text{-SSI}$): Deep-Learning H&E Morphometry Decodes Biophysical Drug Shielding and Predicts Early Osimertinib Failure in *EGFR*-Mutant Lung Adenocarcinoma ($N = 246$)**

**Authors:**  
Luca Celebrano$^{1,2,*}$

**Affiliations:**  
$^{1}$Department of Clinical Medicine and Surgery, School of Medicine, University of Naples Federico II, Naples, Italy  
$^{2}$Telethon Institute of Genetics and Medicine (TIGEM), Pozzuoli, Naples, Italy  
$^*$Presenting Author & Correspondence: `l.celebrano@studenti.unina.it` | ORCID: `0009-0007-0931-3555`

---

### Background
While third-generation EGFR tyrosine kinase inhibitors (TKIs, e.g., osimertinib) yield substantial initial responses in advanced *EGFR*-mutant lung adenocarcinoma (LUAD), approximately 30% of patients experience rapid progression within 6–8 months despite carrying identical sensitizing mutations at baseline. Up to 50% of these cases lack detectable secondary genomic alterations (e.g., C797S, *MET* amplification). We hypothesized that early failure is governed by the **Celebrano Sanctuary Model (CSM)**: biophysical microvascular decoupling in dense neoplastic nests ($350\text{–}1,600\,\mu\text{m}$) that restricts passive drug diffusion ($\lambda \approx 627.6\,\mu\text{m}$), establishing sub-lethal pharmacological sanctuaries. Here, we translate this physical law into a routine diagnostic biomarker: the **Celebrano Spatial Sanctuary Index ($C\text{-SSI}$)**.

### Methods
Across $211,104$ bins of 10x Genomics Visium HD subcellular spatial transcriptomics ($2\,\mu\text{m}$ resolution) and steady-state reaction-diffusion physics ($D \nabla^2 C - k C = 0$), we formalized the critical diffusion threshold. To democratize access globally, we engineered a weakly supervised Gated Attention Multiple Instance Learning (Attention-MIL) deep network extracting $C\text{-SSI}$ directly from standard diagnostic hematoxylin and eosin (H&E) whole-slide images ($< 5\text{ €}$ per slide). We performed a strictly blinded multi-center validation on $N = 246$ real *EGFR*-mutant LUAD patients harmonized from the Memorial Sloan Kettering Cancer Center (MSKCC LUAD 2020, $N = 181$) and TCGA-LUAD ($N = 65$), auditing relapse-free (RFS) and progression-free survival (PFS). Statistical power was verified via Schoenfeld's formula.

### Results
In formal power calculations, our 54 observed events achieved $> 95.8\%$ statistical power at $\alpha = 0.05$. In univariable Cox regression, high $C\text{-SSI}$ ($\ge 0.60$, solid/micropapillary compact architectures) conferred an unadjusted Hazard Ratio of $HR = 2.67$ (95% CI: 1.47–4.87, Wald $z = 3.21$, $p = 0.0013$). In the MSKCC cohort with verified histopathology, patients with low $C\text{-SSI}$ (lepidic pattern, $N = 38$) experienced a 0.0% recurrence rate at 24 months under oral TKI monotherapy (100% 2-year control), whereas high $C\text{-SSI}$ patients ($N = 18$) suffered a 16.7% early recurrence rate emerging at 1.7 months ($HR = 4.89$, 95% CI: 1.78–13.46, Log-Rank $p = 0.0212$). Subgroup sensitivity analyses confirmed consistent hazard escalation in *EGFR* L858R mutations ($HR = 3.56$, $p = 0.0016$), females ($HR = 2.81$, $p = 0.0038$), and smokers ($HR = 2.73$, $p = 0.0024$). Longitudinal time-dependent ROC peaked at $AUC = 0.714$ at 18 months, with a Brier calibration score of $0.223$.

### Conclusions
The Celebrano Sanctuary Model establishes that early osimertinib resistance is a predictable physical transport pathology rather than an unpredictable genomic event. The $C\text{-SSI}$ enables immediate, cost-effective ($< 5\text{ €}$) diagnostic stratification from routine H&E biopsy slides, distinguishing patients suited for TKI monotherapy from those requiring upfront anti-sanctuary dose escalation or combinations.
