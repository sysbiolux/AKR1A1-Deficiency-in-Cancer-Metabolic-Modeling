
# Scripts

This directory contains the MATLAB scripts used throughout the computational workflow, from RNA-seq preprocessing to metabolic model reconstruction, metabolic simulations, statistical analyses, and figure generation.

The workflow supports two complementary datasets:

- Standard cell-line models
- TCGA AKR1A1 ON/OFF models

Whenever possible, the same analytical procedures are used for both datasets while preserving the original biological comparisons and output formats.

---

# Workflow

The scripts are typically executed in the following order.

## 1. Data preprocessing

Purpose:

- import RNA-seq data
- normalize expression values
- discretize expression
- prepare rFASTCORMICS inputs

Main script:

- `driverData.m`

For TCGA analyses, patients are first divided into AKR1A1 HIGH (ON) and LOW (OFF) groups before preprocessing.

---

## 2. Model reconstruction

Purpose:

- reconstruct context-specific metabolic models using rFASTCORMICS
- generate consensus and sample-specific models
- reconstruct TCGA ON/OFF models

Main script:

- `driverModel_withoutO2S.m`

---

## 3. Medium constraints

Purpose:

Apply experimentally defined medium constraints prior to metabolic simulations.

Scripts:

- `setMediumConstraints*.m`

Different media are applied to renal and liver models according to the experimental conditions.

---

## 4. Metabolic analyses

Purpose:

Analyze reconstructed models using COBRA Toolbox methods.

Main script:

- `analysis.m`

Analyses include:

- Flux Balance Analysis (FBA)
- Flux Variability Analysis (FVA)
- FVA similarity analysis
- model similarity
- pathway activity
- gene deletion
- drug perturbation analyses

---

## 5. Flux sampling

Flux sampling is performed on an HPC system using ACHR sampling.

Sampling scripts generate reaction flux distributions that are subsequently analyzed by the statistical workflows.

Separate sampling workflows are available for:

- standard cell-line models
- TCGA ON/OFF models

---

## 6. Statistical analyses

Two complementary statistical analyses are provided.

### Flux Sum analysis

Evaluates metabolite turnover differences between biological conditions using sampled flux distributions.

### Flux Sampling analysis

Evaluates reaction-level differences directly from sampled reaction distributions.

Both analyses compute descriptive statistics, effect sizes, and significance measures used for downstream interpretation.

---

## 7. Visualization

Visualization scripts generate publication-quality figures from the statistical output tables.

Outputs include:

- Flux Sum figures
- Flux Sampling figures
- FVA similarity heatmaps
- pathway-specific figures
- manuscript composite figures

Figures are exported in SVG and PDF formats.

---

# Reproducibility

The scripts preserve the original computational workflow used in the study.

Code revisions were limited to improving readability, documentation, folder organization, and relative file paths without modifying the underlying analytical methods or numerical calculations.
