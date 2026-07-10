# AKR1A1-Deficiency-in-Cancer-Metabolic-Modeling

Genome-scale metabolic modeling workflow investigating AKR1A1 deficiency as a metabolic vulnerability in renal cell carcinoma and hepatocellular carcinoma.  

AKR1A1 deficiency promotes a metabolic state associated with altered glycolytic regulation, methylglyoxal-associated glycation stress, and adaptive metabolic rewiring. To characterize these changes at the systems level, transcriptomic data from cancer cell-line models and TCGA patient tumors were integrated with genome-scale metabolic reconstructions.

The workflow enables the generation and comparison of AKR1A1-deficient and AKR1A1-expressing metabolic models, supporting the identification of altered metabolic pathways, flux distributions, and potential metabolic vulnerabilities.

This repository integrates cancer cell-line and TCGA patient tumor RNA-seq data with rFASTCORMICS/Recon3D-based genome-scale metabolic reconstruction, medium-specific constraints, FBA, FVA, flux sampling, flux-sum analysis, and IDARE/Cytoscape visualization.

---

## PUBLICATION

# Crosstalk between S-nitrosylation and glycation defines a novel metabolic vulnerability in liver and renal cancers

Chiara Pecorari, Mojca Bratina, Evelyn Gonzalez, Salvatore Rizza, Letizia Incampo, Lina Vardouli, Paola Giglio, Zsófia Márta Sztupinszki, Perrine Verdys, Mario Presti, Maria Pires Pacheco, Yuya Qiu, Trine Skov Petersen, Julie Lund Petersen, Yonglun Luo, Emmanuelle Bignon, Zoltan Szallasi, Marco Donia, Jonathan S. Stamler, Simone Cardaci, Thomas Sauter, Giuseppe Filomeni

**Status**

Accepted in *Nature Communications*

---

## SUMMARY

Metabolic reprogramming is a defining feature of cancer and contributes to therapy resistance.

This study investigates how loss of **aldo-keto reductase family 1 member A1 (AKR1A1)** rewires metabolism in:

- renal cell carcinoma (RCC)
- hepatocellular carcinoma (HCC)

AKR1A1 deficiency alters glycolytic metabolism through the crosstalk between:

- S-nitrosylation-dependent metabolic regulation
- methylglyoxal-associated glycation stress

leading to metabolic adaptations associated with NRF2 activation, chemoresistance, and tumor progression.

This repository provides the computational workflow used to reconstruct and analyze genome-scale metabolic models of AKR1A1 deficiency.

---

## METHOD OVERVIEW

The computational framework integrates transcriptomic datasets with genome-scale metabolic modeling to identify AKR1A1-dependent metabolic alterations.

The workflow includes:

## Data preprocessing

- RNA-seq data processing
- expression distribution assessment
- PCA analysis
- TPM/FPKM normalization
- gene-length correction
- AKR1A1 expression-based patient stratification

---

## Context-specific model reconstruction

Genome-scale metabolic models were reconstructed using:

- rFASTCORMICS
- Recon3D human metabolic reconstruction
- MATLAB R2019b

Models were generated from:

- TCGA patient tumor RNA-seq data
- AKR1A1-deficient cancer cell-line RNA-seq data

---

## Medium constraints

Experimental medium conditions were incorporated:

**RPMI**

- renal cancer models
- 769-P cell-line models

**DMEM**

- liver cancer models
- Huh7 cell-line models

---

## Metabolic simulations

The reconstructed models were analyzed using:

- Flux Balance Analysis (FBA)
- Flux Variability Analysis (FVA)
- FVA interval similarity analysis
- single-gene deletion simulations
- drug perturbation simulations
- flux sampling
- flux-sum analysis

Altered metabolic pathways were visualized using:

- IDARE
- Cytoscape

---

# DATASETS

## TCGA patient RNA-seq data

Patient-derived metabolic models were reconstructed using TCGA RNA-seq data.

GEO accession:

```text
GSE62944
```

Cancer cohorts analyzed:

- **KIRC** — kidney renal clear cell carcinoma
- **KIRP** — kidney renal papillary cell carcinoma
- **KICH** — kidney chromophobe carcinoma
- **LIHC** — liver hepatocellular carcinoma

For each cancer type, tumors were stratified according to AKR1A1 mRNA expression.

```text
AKR1A1 HIGH / ON = upper expression quartile

AKR1A1 LOW / OFF = lower expression quartile
```

Generated TCGA metabolic models:

```text
KIRC_ON      vs      KIRC_OFF

KIRP_ON      vs      KIRP_OFF

KICH_ON      vs      KICH_OFF

LIHC_ON      vs      LIHC_OFF
```

---

## Cancer cell-line RNA-seq data

Experimental AKR1A1-deficient models were generated from:

- 769-P renal cancer cells
- Huh7 hepatocellular carcinoma cells

Conditions:

```text
Control

AKR1A1-deficient clone sc1

AKR1A1-deficient clone sc2

AKR1A1-deficient clone sc12
```

GEO accessions:

```text
GSE310784

GSE310828
```

---

# PREREQUISITES

## Software environment

The workflow was developed and tested using:

- MATLAB R2019b
- COBRA Toolbox
- rFASTCORMICS
- IBM ILOG CPLEX Optimizer
- Recon3D metabolic reconstruction

---

# INSTALLATION

## COBRA Toolbox

Install COBRA Toolbox:

https://opencobra.github.io/cobratoolbox/latest/installation.html

Initialize in MATLAB:

```matlab
initCobraToolbox(false)
```

Check installation:

```matlab
which optimizeCbModel
which fluxVariability
```

---

# rFASTCORMICS

Install rFASTCORMICS:

https://github.com/sysbiolux/rFASTCORMICS

Add to MATLAB path:

```matlab
addpath(genpath('path_to_rFASTCORMICS'))
savepath
```

Check:

```matlab
which fastcormics_RNAseq
```

---

# Solver

CPLEX is recommended.

Activate solver:

```matlab
changeCobraSolver('ibm_cplex','all')
```

Check:

```matlab
getCobraSolver
```

---

# WORKFLOW

The repository follows the complete computational pipeline from transcriptomic preprocessing to metabolic pathway interpretation.

```text
RNA-seq datasets
        |
        v
Data preprocessing
        |
        v
Genome-scale model reconstruction
        |
        v
Medium constraints
        |
        v
Metabolic simulations
        |
        v
Flux sampling
        |
        v
Flux-sum statistics
        |
        v
IDARE/Cytoscape visualization
```

---

## 1. TCGA data preparation

Purpose:

- extract KIRC, KIRP, KICH, and LIHC datasets
- classify patients according to AKR1A1 expression
- generate HIGH/LOW groups

---

## 2. Transcriptomic preprocessing

Module:

```text
Data driver
```

Purpose:

- import TCGA and cell-line RNA-seq data
- process gene annotations
- normalize expression values
- generate rFASTCORMICS input files

---

## 3. Genome-scale model reconstruction

Module:

```text
Model driver
```

Purpose:

- reconstruct context-specific metabolic models
- generate patient-derived TCGA models
- generate AKR1A1-deficient cell-line models

---

## 4. Medium-specific constraints

Module:

```text
Medium constraint workflow
```

Purpose:

Apply experimentally defined medium constraints and validate model feasibility.

---

## 5. Metabolic analysis

Module:

```text
Universal analysis workflow
```

Includes:

- FBA
- FVA
- FVA similarity analysis
- single-gene deletion
- drug perturbation analysis

---

## 6. Flux sampling

Purpose:

Explore feasible metabolic flux distributions across biological conditions.

Comparisons:

- AKR1A1 HIGH vs LOW patient models
- Control vs AKR1A1-deficient cell-line models

---

## 7. Flux-sum statistics

Includes:

- metabolite turnover analysis
- Wilcoxon statistical testing
- signal-to-noise calculation
- pathway prioritization

---

## 8. Network visualization

Differential metabolic pathways are visualized using:

- IDARE
- Cytoscape

---

# DETAILED REPRODUCIBILITY GUIDE

Complete MATLAB execution instructions, driver descriptions, configuration options, input/output files, and validation tests are available in:

```text
docs/USER_GUIDE.md
```

---

# CITATION

If you use this workflow, please cite:

**Crosstalk between S-nitrosylation and glycation defines a novel metabolic vulnerability in liver and renal cancers**

*Nature Communications*

Associated Zenodo DOI will be provided with the archived release.
