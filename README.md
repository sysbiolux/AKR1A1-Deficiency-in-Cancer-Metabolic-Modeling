# AKR1A1-Deficiency-in-Cancer-Metabolic-Modeling

Genome-scale metabolic modeling workflow investigating AKR1A1 deficiency as a metabolic vulnerability in renal cell carcinoma (RCC) and hepatocellular carcinoma (HCC).

AKR1A1 deficiency promotes a metabolic state associated with altered glycolytic regulation, methylglyoxal-associated glycation stress, and adaptive metabolic rewiring. To characterize these changes at the systems level, transcriptomic data from cancer cell lines and TCGA patient tumors were integrated with genome-scale metabolic reconstructions.

The repository provides a reproducible computational workflow for RNA-seq preprocessing, context-specific metabolic model reconstruction, medium constraint application, metabolic simulations, flux sampling, statistical analyses, and pathway visualization.

---

# Publication

# Crosstalk between S-nitrosylation and glycation defines a novel metabolic vulnerability in liver and renal cancers

Chiara Pecorari, Mojca Bratina, Evelyn Gonzalez, Salvatore Rizza, Letizia Incampo, Lina Vardouli, Paola Giglio, Zsófia Márta Sztupinszki, Perrine Verdys, Mario Presti, Maria Pires Pacheco, Yuya Qiu, Trine Skov Petersen, Julie Lund Petersen, Yonglun Luo, Emmanuelle Bignon, Zoltan Szallasi, Marco Donia, Jonathan S. Stamler, Simone Cardaci, Thomas Sauter, Giuseppe Filomeni

**Nature Communications acepted**

---

# Project Summary

Metabolic reprogramming is a defining feature of cancer and contributes to tumor progression, therapy resistance, and adaptation to oxidative stress.

This repository investigates how AKR1A1 deficiency remodels cellular metabolism in:

- Renal cell carcinoma (RCC)
- Hepatocellular carcinoma (HCC)

using genome-scale metabolic modeling guided by transcriptomic data.

Two complementary datasets are analyzed:

- Experimental RNA-seq from AKR1A1 knockdown cancer cell lines.
- TCGA patient tumor RNA-seq stratified according to AKR1A1 expression.

The computational workflow reconstructs context-specific metabolic models and compares their metabolic behavior using multiple constraint-based modeling approaches.

---

# Computational Workflow

The repository follows the complete computational pipeline from transcriptomic data to biological interpretation.

```text
RNA-seq datasets
        │
        ▼
Data preprocessing
(driverData.m)
        │
        ▼
Context-specific model reconstruction
(driverModel_withoutO2S.m)
        │
        ▼
Medium constraints
(setMediumConstraints*.m)
        │
        ▼
Metabolic analyses
(analysis.m)
        │
        ▼
Flux sampling (HPC)
        │
        ▼
Flux Sum analysis
Flux Sampling analysis
        │
        ▼
Visualization
```

---

# Workflow Components

## 1. Data preprocessing

RNA-seq datasets are processed to generate discretized gene-expression matrices suitable for rFASTCORMICS reconstruction.

This stage includes:

- RNA-seq import
- Gene annotation processing
- Gene-length normalization
- TPM/FPKM calculation
- Expression discretization
- PCA and descriptive statistics
- AKR1A1-based patient stratification for TCGA datasets

Main script:

- `driverData.m`

---

## 2. Context-specific metabolic model reconstruction

Context-specific genome-scale metabolic models are reconstructed using:

- Recon3D
- rFASTCORMICS
- discretized RNA-seq expression

Models are generated for:

- Cell-line consensus models
- Cell-line sample-specific models
- TCGA AKR1A1 ON/OFF consensus models

Main script:

- `driverModel_withoutO2S.m`

---

## 3. Medium constraints

Experimental culture media are incorporated into the reconstructed models before metabolic simulations.

Media include:

- RPMI (769-P models)
- DMEM (Huh7 models)

Separate scripts are provided for the standard and TCGA ON/OFF workflows.

---

## 4. Metabolic analyses

The reconstructed models are analyzed using COBRA Toolbox methods including:

- Flux Balance Analysis (FBA)
- Flux Variability Analysis (FVA)
- FVA similarity analysis
- Single-gene deletion analysis
- Drug perturbation simulations
- Model similarity analyses
- Pathway activity analyses

Main script:

- `analysis.m`

---

## 5. Flux sampling

Flux sampling explores the feasible solution space of each reconstructed model using ACHR sampling.

Sampling is performed on a High Performance Computing (HPC) system.

Both workflows are supported:

- Standard cell-line models
- TCGA ON/OFF models

---

## 6. Flux statistics

Sampled flux distributions are analyzed using two complementary approaches:

### Flux Sum analysis

Evaluates metabolite turnover differences between biological conditions using:

- Mean Flux Sum
- Log2 Fold Change
- Signal-to-Noise Ratio (SNR)
- Wilcoxon rank-sum statistics

### Flux Sampling analysis

Evaluates reaction-level differences directly from sampled reaction flux distributions using the same statistical framework.

---

## 7. Visualization

Publication-quality figures are generated directly from the statistical output tables.

Visualization scripts include:

- Flux Sum figures
- Flux Sampling figures
- FVA similarity heatmaps
- Pathway-specific figures
- Combined manuscript figures

Figures are exported as:

- SVG
- PDF

---

# Repository Organization

```text
project/

├── data/
│
├── scripts/
│
├── results/
│
├── docs/
│
└── README.md
```

Detailed documentation for each component is provided within the corresponding directories.

---

# Software Requirements

The workflow was developed and tested using:

- MATLAB R2019b
- COBRA Toolbox
- rFASTCORMICS
- IBM ILOG CPLEX Optimizer
- Recon3D human metabolic reconstruction

---

# Installation

## COBRA Toolbox

Install following the official documentation:

https://opencobra.github.io/cobratoolbox/latest/installation.html

Initialize in MATLAB:

```matlab
initCobraToolbox(false)
```

Verify installation:

```matlab
which optimizeCbModel
which fluxVariability
```

---

## rFASTCORMICS

Install from:

https://github.com/sysbiolux/rFASTCORMICS

Add to the MATLAB path:

```matlab
addpath(genpath('path_to_rFASTCORMICS'))
savepath
```

Verify installation:

```matlab
which fastcormics_RNAseq
```

---

## Solver

IBM CPLEX is recommended.

Activate the solver:

```matlab
changeCobraSolver('ibm_cplex','all')
```

Verify:

```matlab
getCobraSolver
```

---

# Documentation

Additional documentation is provided within the repository:

- **data/README.md**

  Description of the datasets, GEO accessions, processed data, and directory organization.

- **scripts/README.md**

  Description of the MATLAB workflows, execution order, HPC sampling pipeline, statistical analyses, and visualization scripts.

- **docs/USER_GUIDE.md**

  Detailed reproducibility guide, execution instructions, expected inputs and outputs, and workflow validation.

---

# Citation

If you use this workflow, please cite:

**Crosstalk between S-nitrosylation and glycation defines a novel metabolic vulnerability in liver and renal cancers**

*Nature Communications*

A Zenodo DOI corresponding to the archived version of this repository will be added upon public release.

---

# License

This repository is distributed for academic research and reproducibility purposes. Please cite the associated publication when using the code or derived analyses.
