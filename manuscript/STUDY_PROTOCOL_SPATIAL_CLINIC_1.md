# PROTOCOLLO DI STUDIO CLINICO OSSERVAZIONALE RETROSPETTIVO (SPATIAL-CLINIC-1)

**Titolo del Progetto:**  
*Validazione Clinica dell'Indice di Santuario Spaziale (SSI) mediante Morfometria Computazionale su Vetrini Istologici H&E per la Stratificazione della Risposta a Inibitori Tirosin-Chinasici (Osimertinib) nell'Adenocarcinoma Polmonare Avanzato EGFR-Mutato.*

**Codice Protocollo:** `SPATIAL-CLINIC-1-LUAD`  
**Sponsor / Promotore Accademico:** Università degli Studi di Napoli Federico II / TIGEM  
**Principal Investigator (PI):** Luca Celebrano (ORCID: `0009-0007-0931-3555`)  
**Data del Protocollo:** 24 Ottobre 2026 (Versione 1.0)  
**Quadro Regolatorio:** Conforme a Good Clinical Practice (ICH-GCP E6 R2), Regolamento UE 2016/679 (GDPR), Linee Guida AIFA per gli studi osservazionali retrospettivi.

---

## 1. Razionale Scientifico e Premesse Biologiche
Nonostante l'introduzione di Osimertinib (inibitore EGFR di terza generazione) come standard di prima linea per l'adenocarcinoma polmonare non a piccole cellule (NSCLC) con mutazioni sensibilizzanti di *EGFR* (delezione esone 19 o L858R), circa il **20–30% dei pazienti manifesta progressione precoce entro 6–8 mesi** dall'avvio della terapia, a fronte di una sopravvivenza libera da progressione (PFS) mediana di 18.9 mesi (trial *FLAURA*).

Ad oggi, non esistono biomarcatori genetici predittivi affidabili al basale in grado di identificare i pazienti a rischio di progressione rapida, poiché il profilo mutazionale del DNA genomico è indistinguibile tra responder e non-responder.

Studi biofisici recenti (Celebrano et al., 2026, Zenodo DOI: `10.5281/zenodo.23623653`) hanno dimostrato che la resistenza precoce è governata dalla **microarchitettura tissutale spaziale**:
1. Nidi tumorali solidi o micropapillari compatti con diametro $> 450\,\mu\text{m}$ privi di capillarizzazione intrinseca generano una barriera diffusiva per il farmaco (Legge di Fick, lunghezza di diffusione $\lambda \approx 627\,\mu\text{m}$).
2. Nella zona centrale del nido si crea una **"finestra sub-letale di persistenza"** ($\text{IC}_{50} \le C(r) \le C_{\text{kill}}$) che protegge le cellule neoplastiche (Drug-Tolerant Persisters, DTP) favorendo l'emergere di cloni resistenti secondari.
3. Tale architettura è quantificabile mediante l'**Indice di Santuario Spaziale (SSI)** calcolato direttamente sui comuni vetrini diagnostici colorati con Ematossilina-Eosina (H&E) tramite un algoritmo di Deep Learning ad attenzione (Attention-MIL).

---

## 2. Obiettivi dello Studio

### 2.1 Endpoint Primario
* Valutare l'accuratezza predittiva (misurata mediante Harrell's Concordance Index, C-index) dell'Indice di Santuario Spaziale ($SSI$, scala $0.0\text{–}1.0$) calcolato su vetrini H&E pre-trattamento nel predire la **Sopravvivenza Libera da Progressione (PFS)** a 12 e 18 mesi in pazienti con LUAD *EGFR*-mutato trattati con Osimertinib.

### 2.2 Endpoint Secondari
1. **Separazione delle Curve di Kaplan-Meier:** Significatività statistica (Log-Rank test $p < 0.05$) tra i gruppi a Basso Rischio ($SSI < 0.35$), Intermedio ($0.35 \le SSI < 0.60$) e Alto Rischio ($SSI \ge 0.60$).
2. **Hazard Ratio (HR) Aggiustato:** Calcolo dell'Hazard Ratio multivariato mediante modello di Cox aggiustando per età, sesso, stadio clinico TNM iniziale (IIIB vs IV) e tipo di mutazione *EGFR* (Exon 19 del vs L858R).
3. **Metriche Diagnostiche di Progressione Precoce (< 8 mesi):** Calcolo di Sensibilità, Specificità, Valore Predittivo Positivo (PPV), Valore Predittivo Negativo (NPV) e Area Sotto la Curva ROC (AUC).
4. **Valutazione della Robustezza Multi-Scanner:** Verifica dell'assenza di bias legati al costruttore dello scanner digitale (Aperio vs Hamamatsu vs Philips) previa normalizzazione cromatica (algoritmo di Macenko).

---

## 3. Disegno dello Studio e Popolazione di Pazienti

### 3.1 Disegno
Studio osservazionale retrospettivo multicentrico no-profit su tessuto bioptico d'archivio e dati clinici anonimizzati.

### 3.2 Criteri di Inclusione
1. Età $\ge 18$ anni alla diagnosi.
2. Diagnosi istologica confermata di Adenocarcinoma Polmonare (LUAD).
3. Presenza documentata di mutazione sensibilizzante di *EGFR* (Exon 19 del o L858R) determinata mediante NGS o RT-PCR.
4. Trattamento sistemico con Osimertinib in prima linea a dose standard (80 mg/die).
5. Disponibilità di vetrino istologico H&E diagnostico (o scansione digitale WSI) pre-trattamento (biopsia tissutale o resezione chirurgica) con cellularità tumorale $\ge 20\%$.
6. Follow-up clinico documentato minimo di 24 mesi o fino a progressione documentata/decesso.

### 3.3 Criteri di Esclusione
1. Istotipo squamocellulare o a piccole cellule (SCLC).
2. Presenza al basale di mutazioni di resistenza note (es. *EGFR* T790M pre-trattamento, amplificazione *MET* nota prima di Osimertinib).
3. Terapia concomitante con altri agenti antineoplastici sperimentali o radioterapia concomitante sul tumore primitivo durante la prima linea.
4. Vetrini istologici gravemente danneggiati, con artefatti da piegatura o cellularità tumorale inadeguata ($< 20\%$).

---

## 4. Metodologia di Analisi e Pipeline Algoritmica

```mermaid
flowchart LR
    A["Vetrino H&E Anonimizzato (.svs / .ndpi)"] --> B["Tiling Automatico & Normalizzazione Macenko"]
    B --> C["Attention-MIL Deep Architecture"]
    C --> D["Calcolo Indice SSI (0.0 - 1.0)"]
    D --> E["Stratificazione Cieca: Basso vs Alto Rischio"]
    E --> F["Unblinding: Confronto con PFS Reale"]
```

1. **Digitalizzazione & Normalizzazione:** Le Whole Slide Images (WSI) vengono acquisite a ingrandimento 20x/40x ($0.50\text{–}0.25\,\mu\text{m/pixel}$). Vengono rimossi sfondo e frammenti acellulari, e applicata la normalizzazione cromatica (Macenko stain normalization).
2. **Estrazione delle Tessere (Tiling):** Generazione di tessere $512 \times 512$ pixel con overlap del $20\%$.
3. **Calcolo Blinded dell'Indice $SSI$:** L'algoritmo elabora le tessere senza alcuna informazione sull'outcome del paziente, stimando la percentuale di architettura solida/micropapillare e il diametro medio dei nidi compatti, restituendo l'indice $SSI$.
4. **Analisi Cieca:** Le previsioni vengono congelate e marcate con timestamp prima dell'apertura del database clinico degli eventi di progressione.

---

## 5. Aspetti Etici, Privacy e Gestione dei Dati (GDPR)

* **Pseudonimizzazione Rigorosa:** A ciascun paziente viene assegnato un codice identificativo univoco progressivo (es. `SPAT-001`, `SPAT-002`). La chiave di decifratura che associa il codice ai dati anagrafici resta custodita esclusivamente dal centro clinico promotore in ambiente protetto con crittografia AES-256.
* **Nessun Impatto Terapeutico per i Pazienti:** Trattandosi di uno studio retrospettivo su casi storici chiusi, l'analisi non interferisce in alcun modo con le decisioni cliniche né con il percorso terapeutico dei soggetti.
* **Esenzione dal Consenso Informato:** In accordo con l'Art. 110-bis del D.Lgs. 196/2003 (Codice Privacy italiano) e con le Autorizzazioni Generali del Garante per la Protezione dei Dati Personali relative alla ricerca scientifica, per i pazienti deceduti o irreperibili verrà richiesta al Comitato Etico l'autorizzazione al trattamento dei dati retrospettivi pseudonimizzati.

---

## 6. Risultati Attesi e Impatto Clinico
La validazione positiva del protocollo (C-index $\ge 0.75$, Log-Rank $p < 0.01$) permetterà di:
1. Depositare domanda di brevetto per il metodo diagnostico computazionale basato su H&E (SaMD - Software as a Medical Device).
2. Fornire agli oncologi medici uno strumento decisionale a costo quasi zero ($< 5\text{ €}$) per identificare al momento della biopsia i pazienti con "santuari fisici" attivi.
3. Proporre per tali pazienti strategie terapeutiche intensificate immediate (es. dosaggio pulsato di Osimertinib a 160 mg o combinazione precoce con anti-angiogenici/chemio come da trial *FLAURA2* / *MARIPOSA*).

---
*Documento redatto e validato in conformità con gli standard internazionali ICH-GCP.*
