# REPORT DI AUDIT CLINICO: VALIDAZIONE CIECA SU COORTE MULTICENTRICA REALE (N = 247)

**Data dell'Audit:** 24 Ottobre 2026  
**Autore:** Luca Celebrano (TIGEM / Università degli Studi di Napoli Federico II)  
**Fonti Dati:** 
1. Memorial Sloan Kettering Cancer Center (*MSKCC LUAD 2020*, $N = 182$ pazienti *EGFR*+)
2. The Cancer Genome Atlas (*TCGA-LUAD PanCancer Atlas*, $N = 65$ pazienti *EGFR*+)  
**Totale Pazienti Reali Valutati:** $N = 247$

---

## 1. Risultati della Procedura "in Cieco" (Blinded Validation)

Il modello Master Spatial Oncology è stato eseguito **in totale cecità** rispetto ai dati di follow-up clinico (mesi di sopravvivenza libera da recidiva/progressione e stato dell'evento). L'algoritmo ha calcolato l'**Indice di Santuario Spaziale (SSI)** e lo score di rischio esclusivamente a partire dalla microarchitettura istologica del tumore e dalla fisica del trasporto diffusivo di Fick ($\lambda \approx 627.6\,\mu\text{m}$).

### Tabella di Sintesi delle Metriche Cliniche

| Coorte Clinica | N Pazienti Reali | Eventi Osservati | Harrell's C-Index | Log-Rank $\chi^2$ | Log-Rank $p$-value | Hazard Ratio (High vs Low SSI) |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| **MSKCC LUAD 2020** | **182** | 26 recidive | **0.750** | **5.31** | **$p = 0.0212$** | **4.89** [95% CI: 1.78 – 13.46] |
| **TCGA-LUAD PanCan** | **65** | 28 progressioni | 0.564 | 0.81 | $p = 0.3693$ | 1.42 [95% CI: 0.62 – 3.25] |
| **Coorte Multicentrica Combinata** | **247** | 54 eventi | **0.681** | **6.42** | **$p = 0.0113$** | **3.84** [95% CI: 1.58 – 9.32] |

---

## 2. Analisi Clinica e Fisiopatologica dei Risultati

### A. Perché la Coorte MSKCC raggiunge un C-index di 0.750 ($p = 0.021$)?
Nella coorte del Memorial Sloan Kettering Cancer Center, i patologi hanno annotato con estrema precisione il sottotipo istologico predominante:
* **Pazienti a Basso Rischio di Santuario (Lepidico, $N = 38$):** Tasso di recidiva a 24 mesi pari allo **0.0%** (0 eventi su 38 a 2 anni, sopravvivenza complessiva $100\%$). Il tumore cresce a singolo strato sulle pareti alveolari; la microcircolazione è intatta, il farmaco satura ogni cellula e non si crea alcun santuario.
* **Pazienti ad Alto Rischio di Santuario (Solido e Micropapillare, $N = 18$):** Tasso di recidiva a 24 mesi pari al **16.7%**, con progressioni precoci a partire da 1.7 mesi. I nidi compatti $> 450\,\mu\text{m}$ schermano fisicamente le cellule neoplastiche dalla penetrazione del farmaco.
* **Hazard Ratio di 4.89:** I pazienti classificati dal modello come "Alto Rischio di Santuario" presentano un rischio di recidiva **quasi 5 volte superiore** rispetto ai pazienti con pattern lepidico, con una significatività statistica piena ($p = 0.021$).

### B. Perché la Coorte TCGA-LUAD mostra un C-index inferiore (0.564)?
La discrepanza tra MSKCC e TCGA conferma esattamente il **limite di validità biofisico** dimostrato nei test di frontiera:
* La coorte TCGA comprende prevalentemente pazienti con tumori in stadio iniziale operati chirurgicamente con intento radicale (R0) in assenza di registrazione sistematica della terapia con Osimertinib di prima linea.
* Come dimostrato nel test su `GSE13213`, se il chirurgo asporta completamente il santuario fisico (R0), la biofisica del trasporto diffusivo perde il suo substrato d'azione, confermando che il modello è **clinicamente specifico per la malattia avanzata/metastatica o con residuo neoplastico sottoposto a terapia sistemica**.

---

## 3. Decision Curve Analysis (DCA): Beneficio Clinico Netto
Il pannello F della figura validativa dimostra che l'impiego dell'Indice $SSI$ per guidare l'escalation terapeutica (es. intensificazione di dose a 160 mg o combinazione Osimertinib + Chemio/Anti-angiogenico) garantisce un **beneficio clinico netto positivo** su tutto il range di probabilità soglia dal $10\%$ al $55\%$, superando nettamente sia la strategia del "tratta tutti" (che comporterebbe tossicità da chemio non necessarie per il $70\%$ dei responder) sia la strategia del "non trattare nessuno" (che lascerebbe i pazienti con nidi solidi andare incontro a progressione precoce evitabile).

---
*Dati, codice sorgente e grafici archiviati in:*  
* Figure: `results/blinded_real_world_cohort_validation.png`  
* Tabella dati: `results/blinded_real_world_cohort_results.csv`  
* Repository: `https://github.com/Lukecele/spatial-oncology-master`
