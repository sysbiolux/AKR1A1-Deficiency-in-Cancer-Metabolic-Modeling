# AKR1A1 Deficiency Metabolic Modeling Pipeline

## Complete reproducibility and usage guide

This document describes how to reproduce the computational metabolic-modeling analyses associated with the study:

**Crosstalk between S-nitrosylation and glycation defines a novel metabolic vulnerability in liver and renal cancers**

The workflow integrates RNA-seq data from:

* TCGA patient tumors;
* 769-P renal cancer cells;
* Huh7 hepatocellular carcinoma cells.

The pipeline performs:

* TCGA cancer-type extraction;
* AKR1A1 HIGH/LOW patient stratification;
* RNA-seq preprocessing;
* TPM/FPKM normalization where required;
* expression discretization;
* rFASTCORMICS model reconstruction;
* medium-specific constraints;
* model feasibility testing;
* Flux Balance Analysis;
* Flux Variability Analysis;
* FVA similarity analysis;
* single-gene deletion;
* drug perturbation analysis;
* flux sampling;
* flux-sum statistics;
* IDARE/Cytoscape visualization.

The code was developed for:

```text
MATLAB R2019b
```

---

# 1. Pipeline overview

The complete analysis is organized into the following modules:

```text
RNA-seq input data
        |
        v
TCGA cancer-type extraction
        |
        v
AKR1A1 HIGH / LOW stratification
        |
        v
Expression preprocessing
        |
        v
Expression discretization
        |
        v
rFASTCORMICS model reconstruction
        |
        v
RPMI / DMEM medium constraints
        |
        v
FBA / FVA / FVA similarity
        |
        v
Gene and drug perturbation analyses
        |
        v
Flux sampling
        |
        v
Flux-sum statistics
        |
        v
IDARE / Cytoscape visualization
```

The TCGA and cell-line workflows share the same reconstruction and analysis structure, but the TCGA workflow includes the additional patient-stratification step.

---

# 2. Main MATLAB files

The core pipeline consists of the following files:

```text
pipeline_config.m
01_geo_gse62944_generate_TCGA_tables.m
02_prepare_expression_data.m
03_build_fastcormics_models.m
04_run_model_analysis.m
```

Additional project scripts include:

```text
driverData.m
driverModel.m
test_import_and_high_low.m
high_low.m
read_excel_data.m
import_TCGA_cancer_type_samples.m
TCGA_annotation_and_highlow_from_existing_files.m
run_high_low_from_existing_TCGA_files.m

setMediumConstraints_unified.m
setMediumConstraints_TCGA.m
setMediumConstraints_CellLines.m

universal_model_analysis.m
test_universal_analysis_one_model.m
FVA_similarity_Thomas.m
Generate_FVA_TCGA.m
Generate_FVA_celllines.m
Generate_FigureS1_hi.m

RUN_RAW_1500_FluxSampling.m
RUN_ONOFF_1000_FluxSampling.m
sampling_test_medium.m
sampling_test_mediumonoff.m
```

Do not rename `pipeline_config.m`, because the other scripts depend on this exact filename.

---

# 3. Recommended repository structure

```text
AKR1A1-Deficiency-in-Cancer-Metabolic-Modeling/
│
├── README.md
├── pipeline_config.m
│
├── code/
│   ├── 01_data_preparation/
│   ├── 02_model_reconstruction/
│   ├── 03_medium_constraints/
│   ├── 04_model_analysis/
│   ├── 05_flux_sampling/
│   ├── 06_flux_sum/
│   └── 07_visualization/
│
├── data/
│   ├── raw_TCGA/
│   ├── cancer_expression_TCGA/
│   ├── high_low/
│   ├── raw_cell_lines/
│   ├── model/
│   ├── medium/
│   └── drug/
│
├── results/
│   ├── 01_prepared_data/
│   ├── 02_models/
│   ├── 03_analysis/
│   ├── 04_sampling/
│   ├── 05_flux_sum/
│   └── figures/
│
└── docs/
    └── PIPELINE_USAGE.md
```

All code should use relative paths defined in `pipeline_config.m`.

Do not include local machine-specific paths in public code.

---

# 4. Software requirements

Required software:

* MATLAB R2019b;
* COBRA Toolbox;
* rFASTCORMICS;
* Recon3D reference reconstruction;
* a compatible linear-programming solver;
* IBM ILOG CPLEX is recommended.

Optional MATLAB toolboxes may include:

* Statistics and Machine Learning Toolbox;
* Bioinformatics Toolbox.

Required COBRA and rFASTCORMICS functions include:

```matlab
initCobraToolbox
changeCobraSolver
optimizeCbModel
changeObjective
changeRxnBounds
fluxVariability
singleGeneDeletion_rFASTCORMICS
DrugDeletion
FVA_similarity_Thomas
discretize_FPKM
fastcormics_RNAseq
```

Check that MATLAB can locate them:

```matlab
which initCobraToolbox
which optimizeCbModel
which changeObjective
which fluxVariability
which singleGeneDeletion_rFASTCORMICS
which DrugDeletion
which FVA_similarity_Thomas
which discretize_FPKM
which fastcormics_RNAseq
```

Each command should return a valid file path.

---

# 5. Installation

## 5.1 COBRA Toolbox

Install the COBRA Toolbox from:

```text
https://opencobra.github.io/cobratoolbox/latest/installation.html
```

Initialize it in MATLAB:

```matlab
initCobraToolbox(false)
```

Confirm installation:

```matlab
which optimizeCbModel
which fluxVariability
```

---

## 5.2 rFASTCORMICS

Install rFASTCORMICS from:

```text
https://github.com/sysbiolux/rFASTCORMICS
```

Add it to the MATLAB path:

```matlab
addpath(genpath('path_to_rFASTCORMICS'))
savepath
```

Confirm installation:

```matlab
which fastcormics_RNAseq
which discretize_FPKM
which singleGeneDeletion_rFASTCORMICS
```

---

## 5.3 IBM ILOG CPLEX

Install IBM ILOG CPLEX Optimization Studio according to the IBM instructions.

Add the MATLAB interface to the MATLAB path.

Example:

```matlab
addpath(genpath('path_to_CPLEX_matlab_interface'))
savepath
```

Select CPLEX as the COBRA solver:

```matlab
changeCobraSolver('ibm_cplex', 'all')
```

Check:

```matlab
getCobraSolver
```

A different compatible solver may be used if it has been validated for the workflow.

---

# 6. Data sources

## 6.1 TCGA patient RNA-seq data

TCGA RNA-seq data were obtained from GEO accession:

```text
GSE62944
```

Required source files:

```text
GSE62944_06_01_15_TCGA_24_CancerType_Samples.txt
GSM1536837_06_01_15_TCGA_24.tumor_Rsubread_FPKM.txt
GSM1697009_06_01_15_TCGA_24.normal_Rsubread_FPKM.txt
```

Cancer types analyzed:

```text
KICH
KIRC
KIRP
LIHC
```

Expected location:

```text
data/raw_TCGA/
```

---

## 6.2 GSE62944 memory limitation

The complete GSE62944 tumor expression matrix is very large.

In MATLAB R2019b, directly loading:

```text
GSM1536837_06_01_15_TCGA_24.tumor_Rsubread_FPKM.txt
```

with:

```matlab
readtable()
```

may produce an out-of-memory error.

The recommended approach is:

1. download the original GEO files;
2. extract each cancer type once;
3. save one expression file per cancer type;
4. use the cancer-specific files as pipeline inputs.

Expected cancer-specific files:

```text
data/cancer_expression_TCGA/KICH.txt
data/cancer_expression_TCGA/KIRC.txt
data/cancer_expression_TCGA/KIRP.txt
data/cancer_expression_TCGA/LIHC.txt
```

Example structure:

```text
Genes       TCGA_sample1     TCGA_sample2
A1BG        0.21             0.28
A1CF        0.01             0.03
AKR1A1      5.20             0.80
```

These files are equivalent to the outputs of the full GEO extraction step.

---

## 6.3 769-P cell-line RNA-seq data

GEO accession:

```text
GSE310784
```

Expected directory:

```text
data/raw_cell_lines/GSE310784/769-p/
```

Expected files:

```text
769P_project_1_normalizedCountsWithAnnotations.txt
769P_project_1_sampleAttributes.csv
```

Conditions:

```text
control
sc1
sc2
sc12
```

---

## 6.4 Huh7 cell-line RNA-seq data

GEO accession:

```text
GSE310828
```

Expected directory:

```text
data/raw_cell_lines/GSE310828/HuH7/
```

Expected files:

```text
Huh7_project_1_normalizedCountsWithAnnotations.txt
Huh7_project_1_sampleAttributes.csv
```

Conditions:

```text
control
sc1
sc2
sc12
```

---

# 7. Required reference files

Place the reference-model files in:

```text
data/model/
```

Required files:

```text
consistent_model.mat
dico2columns.mat
```

The reference model should contain:

```matlab
consistent_model
```

The project also uses the Recon3D human genome-scale metabolic reconstruction through the rFASTCORMICS workflow.

---

# 8. Required medium files

Place medium-composition files in:

```text
data/medium/
```

Required files:

```text
medium_DMEM_H.xlsx
medium_RPMI_7.xlsx
```

Medium assignments are:

| Dataset                 | Medium |
| ----------------------- | ------ |
| KICH ON/OFF             | RPMI   |
| KIRC ON/OFF             | RPMI   |
| KIRP ON/OFF             | RPMI   |
| LIHC ON/OFF             | DMEM   |
| 769-P ctrl/sc1/sc2/sc12 | RPMI   |
| Huh7 ctrl/sc1/sc2/sc12  | DMEM   |

These assignments must not be changed without biological justification.

---

# 9. Required drug file

Place the gene-drug relationship file in:

```text
data/drug/
```

Required file:

```text
GeneDrugRelationsUpdate.mat
```

Expected variable:

```matlab
GeneDrugRelations
```

Expected field:

```matlab
GeneDrugRelations.DrugName
```

---

# 10. Pipeline configuration

All project paths and execution settings are controlled by:

```text
pipeline_config.m
```

The configuration should define paths such as:

```matlab
cfg.rootDir
cfg.dataDir
cfg.rawTCGADir
cfg.cancerExpressionDir
cfg.highLowDir
cfg.preparedDir
cfg.modelDir
cfg.mediumDir
cfg.drugDir
cfg.resultsDir
cfg.modelsDir
cfg.analysisDir
cfg.samplingDir
cfg.fluxSumDir
cfg.figureDir
```

The main dataset-selection setting is:

```matlab
cfg.runMode
```

Allowed values:

```matlab
cfg.runMode = 'TCGA';
cfg.runMode = 'CELL_LINES';
cfg.runMode = 'ALL';
```

Test the configuration:

```matlab
cfg = pipeline_config;

disp(cfg.rootDir)
disp(cfg.dataDir)
disp(cfg.modelDir)
disp(cfg.modelsDir)
disp(cfg.analysisDir)
```

A valid directory should return:

```matlab
exist(cfg.modelDir, 'dir')
```

Expected result:

```text
7
```

A valid file should return:

```matlab
exist(fullfile(cfg.modelDir, 'consistent_model.mat'), 'file')
```

Expected result:

```text
2
```

---

# 11. TCGA workflow

## Step 1 — Generate cancer-specific and AKR1A1 HIGH/LOW tables

Main scripts:

```text
01_geo_gse62944_generate_TCGA_tables.m
import_TCGA_cancer_type_samples.m
high_low.m
TCGA_annotation_and_highlow_from_existing_files.m
run_high_low_from_existing_TCGA_files.m
```

The pipeline extracts:

```text
KICH
KIRC
KIRP
LIHC
```

and divides samples according to AKR1A1 expression.

Definitions:

```text
AKR1A1 HIGH / ON
=
AKR1A1 expression >= 75th percentile

AKR1A1 LOW / OFF
=
AKR1A1 expression <= 25th percentile
```

Equivalent MATLAB logic:

```matlab
thresholdHigh = prctile(AKR1A1_values, 75);
thresholdLow  = prctile(AKR1A1_values, 25);

selectON  = AKR1A1_values >= thresholdHigh;
selectOFF = AKR1A1_values <= thresholdLow;
```

Run:

```matlab
cfg = pipeline_config;
run('01_geo_gse62944_generate_TCGA_tables.m')
```

When using previously extracted cancer-specific files, run:

```matlab
run_high_low_from_existing_TCGA_files
```

Expected outputs:

```text
data/high_low/Table_ON_KICH.xlsx
data/high_low/Table_OFF_KICH.xlsx
data/high_low/KICH_anno.txt

data/high_low/Table_ON_KIRC.xlsx
data/high_low/Table_OFF_KIRC.xlsx
data/high_low/KIRC_anno.txt

data/high_low/Table_ON_KIRP.xlsx
data/high_low/Table_OFF_KIRP.xlsx
data/high_low/KIRP_anno.txt

data/high_low/Table_ON_LIHC.xlsx
data/high_low/Table_OFF_LIHC.xlsx
data/high_low/LIHC_anno.txt
```

Generated TCGA groups:

```text
KICH_ON
KICH_OFF
KIRC_ON
KIRC_OFF
KIRP_ON
KIRP_OFF
LIHC_ON
LIHC_OFF
```

---

# 12. Cell-line workflow

The cell-line workflow does not require the TCGA HIGH/LOW extraction step.

Set:

```matlab
cfg.runMode = 'CELL_LINES';
```

Then start with expression preparation.

Expected source directories:

```text
data/raw_cell_lines/GSE310784/769-p/
data/raw_cell_lines/GSE310828/HuH7/
```

---

# 13. Expression-data preparation

Main scripts:

```text
02_prepare_expression_data.m
prepare_expression_data.m
driverData.m
read_excel_data.m
```

The data-preparation workflow:

1. loads TCGA HIGH/LOW tables or cell-line RNA-seq data;
2. processes gene identifiers;
3. converts values when necessary;
4. applies TPM/FPKM normalization for cell-line data;
5. discretizes expression values;
6. saves rFASTCORMICS input files.

Run:

```matlab
cfg = pipeline_config;
run('02_prepare_expression_data.m')
```

For TCGA:

```matlab
cfg.runMode = 'TCGA';
```

For cell lines:

```matlab
cfg.runMode = 'CELL_LINES';
```

For both:

```matlab
cfg.runMode = 'ALL';
```

---

## 13.1 Cell-line normalization

The original cell-line workflow converts log-transformed values using:

```matlab
data = 2.^data - 1;
```

TPM is calculated using gene lengths:

```matlab
tempTPM = data ./ repmat(lengths, 1, size(data, 2));

TPM = 1e6 * tempTPM ./ ...
    repmat(nansum(tempTPM), size(data, 1), 1);
```

FPKM is calculated as:

```matlab
tempFPKM = 1e6 * data ./ ...
    repmat(nansum(data), size(data, 1), 1);

FPKM = tempFPKM ./ ...
    repmat(lengths, 1, size(data, 2));
```

Do not modify these transformations unless intentionally revising the preprocessing method.

---

## 13.2 Expression discretization

The workflow applies:

```matlab
discretized = discretize_FPKM( ...
    expressionMatrix, ...
    sampleNames, ...
    1);
```

The discretized matrix generally contains:

```text
 1 = expressed
 0 = intermediate or uncertain
-1 = low expression
```

The exact interpretation follows the installed rFASTCORMICS version.

---

## 13.3 Expected outputs

TCGA outputs:

```text
results/01_prepared_data/discretizedON_KICH.mat
results/01_prepared_data/discretizedOFF_KICH.mat

results/01_prepared_data/discretizedON_KIRC.mat
results/01_prepared_data/discretizedOFF_KIRC.mat

results/01_prepared_data/discretizedON_KIRP.mat
results/01_prepared_data/discretizedOFF_KIRP.mat

results/01_prepared_data/discretizedON_LIHC.mat
results/01_prepared_data/discretizedOFF_LIHC.mat
```

Cell-line outputs:

```text
results/01_prepared_data/discretized7.mat
results/01_prepared_data/discretizedH.mat
```

Expected variables include:

```matlab
discretizedon
discretizedoff
geneSymson
geneSymsoff
colnamesON
colnamesOFF
TPM
FPKM
geneSyms
```

---

## 13.4 Validate expression outputs

Example:

```matlab
cfg = pipeline_config;

S = load(fullfile( ...
    cfg.preparedDir, ...
    'discretizedON_KIRC.mat'));

fieldnames(S)
```

Check dimensions:

```matlab
size(S.discretizedon, 1) == numel(S.geneSymson)
size(S.discretizedon, 2) == numel(S.colnamesON)
```

Check discretization counts:

```matlab
sum(S.discretizedon == 1, 1)
sum(S.discretizedon == 0, 1)
sum(S.discretizedon == -1, 1)
```

For each sample, the three counts should equal the total number of genes.

---

# 14. Context-specific model reconstruction

Main scripts:

```text
03_build_fastcormics_models.m
build_fastcormics_models.m
driverModel.m
```

Required inputs:

```text
consistent_model.mat
dico2columns.mat
discretized expression files
medium files
```

The reconstruction workflow:

1. loads the reference metabolic reconstruction;
2. loads discretized expression data;
3. maps transcriptomic information to the model;
4. applies rFASTCORMICS;
5. creates context-specific models;
6. saves individual and grouped models.

Run:

```matlab
cfg = pipeline_config;
run('03_build_fastcormics_models.m')
```

The original workflow may block oxygen exchange during reconstruction:

```matlab
EX_o2s[e] = 0
```

Do not change this setting unless intentionally revising the reconstruction assumptions.

---

## 14.1 Expected TCGA models

```text
model_KICH_ON.mat
model_KICH_OFF.mat
model_KIRC_ON.mat
model_KIRC_OFF.mat
model_KIRP_ON.mat
model_KIRP_OFF.mat
model_LIHC_ON.mat
model_LIHC_OFF.mat
```

---

## 14.2 Expected cell-line models

```text
model_7_ctrl.mat
model_7_sc1.mat
model_7_sc12.mat
model_7_sc2.mat

model_H_ctrl.mat
model_H_sc1.mat
model_H_sc12.mat
model_H_sc2.mat
```

The source order is important.

Some original model collections use:

```text
ctrl
sc1
sc12
sc2
```

Do not assume:

```text
ctrl
sc1
sc2
sc12
```

unless the source file has been explicitly reordered.

---

## 14.3 Model validation

Every model should contain at least:

```matlab
model.S
model.rxns
model.mets
model.lb
model.ub
model.c
```

Check:

```matlab
size(model.S, 1) == numel(model.mets)
size(model.S, 2) == numel(model.rxns)
numel(model.lb) == numel(model.rxns)
numel(model.ub) == numel(model.rxns)
numel(model.c)  == numel(model.rxns)
```

Required objective reactions:

```matlab
any(strcmp(model.rxns, 'biomass_reaction'))
any(strcmp(model.rxns, 'DM_atp_c_'))
```

Gene-deletion analyses may additionally require:

```matlab
model.genes
model.grRules
model.rules
model.rxnGeneMat
```

---

# 15. Medium-specific constraints

Main scripts:

```text
setMediumConstraints_unified.m
setMediumConstraints_TCGA.m
setMediumConstraints_CellLines.m
```

The medium workflow:

1. loads each context-specific model;
2. assigns RPMI or DMEM;
3. applies nutrient-uptake bounds;
4. applies additional exchange constraints;
5. sets the biomass objective;
6. runs FBA;
7. runs FVA at 90% of optimum;
8. extracts exchange fluxes;
9. saves constrained models;
10. records missing reactions.

Run:

```matlab
setMediumConstraints_unified
```

---

## 15.1 Additional constraints

The workflow preserves the original constraints, including:

```text
EX_HC02161[e] closed
EX_peplys[e] closed
EX_o2s[e] closed or constrained according to the script
oxygen bounds
carbon-dioxide bounds
```

These constraints should remain unchanged unless they are deliberately revised and documented.

---

## 15.2 Required model inputs

Place individual model MAT files in:

```text
data/model/
```

A model file may contain:

```matlab
model_KICH_ON
```

or the equivalent named variable, or:

```matlab
ContextModel
```

The medium script should support the expected format used in the project.

---

## 15.3 Expected constrained outputs

Individual files:

```text
model_KICH_ON_constrained.mat
model_KICH_OFF_constrained.mat
model_KIRC_ON_constrained.mat
model_KIRC_OFF_constrained.mat
model_KIRP_ON_constrained.mat
model_KIRP_OFF_constrained.mat
model_LIHC_ON_constrained.mat
model_LIHC_OFF_constrained.mat

model_7_ctrl_constrained.mat
model_7_sc1_constrained.mat
model_7_sc12_constrained.mat
model_7_sc2_constrained.mat

model_H_ctrl_constrained.mat
model_H_sc1_constrained.mat
model_H_sc12_constrained.mat
model_H_sc2_constrained.mat
```

Grouped files:

```text
models_Kidney_RPMI.mat
models_Liver_DMEM.mat
models_769p_RPMI.mat
models_Huh7_DMEM.mat
```

Analysis file:

```text
medium_constraint_FBA_FVA_results.mat
```

Expected variables:

```matlab
FBAResults
FVAResults
ExchangeResults
missingReactions
modelList
```

---

## 15.4 Validate medium outputs

```matlab
cfg = pipeline_config;

dir(fullfile(cfg.modelsDir, '*_constrained.mat'))
dir(fullfile(cfg.analysisDir, ...
    'medium_constraint_FBA_FVA_results.mat'))
```

Inspect one model:

```matlab
S = load(fullfile( ...
    cfg.modelsDir, ...
    'model_KICH_ON_constrained.mat'));

fieldnames(S)
```

Expected fields may include:

```text
constrainedModel
mediumName
groupName
```

Inspect FBA results:

```matlab
R = load(fullfile( ...
    cfg.analysisDir, ...
    'medium_constraint_FBA_FVA_results.mat'));

R.FBAResults.model_KICH_ON.f
```

A finite numeric value indicates that FBA returned a solution.

---

# 16. Model analysis

Main scripts:

```text
04_run_model_analysis.m
run_model_analysis.m
universal_model_analysis.m
test_universal_analysis_one_model.m
Generate_FVA_TCGA.m
Generate_FVA_celllines.m
FVA_similarity_Thomas.m
```

Analyses include:

* ATP objective;
* biomass objective;
* FBA;
* FVA;
* FVA similarity;
* single-gene deletion;
* drug perturbation analysis.

Run:

```matlab
cfg = pipeline_config;
run('04_run_model_analysis.m')
```

---

# 17. Universal analysis structure

The universal workflow should receive:

```matlab
config.models
config.modelNames
config.consistentModel
config.GeneDrugRelations
config.atpObjective
config.biomassObjective
config.comparisonPairs
config.comparisonNames
config.outputFolder
config.outputPrefix
```

It should determine the model count automatically:

```matlab
numberOfModels = numel(config.models);
```

It should not hard-code eight models.

To add another dataset, create a new configuration:

```matlab
config = config_NewDataset();
```

and add it to the dataset-selection switch:

```matlab
case 'NEW_DATASET'
    config = config_NewDataset();
```

No algorithmic changes to `universal_model_analysis.m` should be required.

---

# 18. TCGA model comparisons

Recommended TCGA order:

```text
1 = KIRC_ON
2 = KIRC_OFF
3 = KIRP_ON
4 = KIRP_OFF
5 = KICH_ON
6 = KICH_OFF
7 = LIHC_ON
8 = LIHC_OFF
```

Comparison matrix:

```matlab
comparisonPairs = [ ...
    1, 2; ...
    3, 4; ...
    5, 6; ...
    7, 8];
```

Comparison names:

```text
KIRC_ON_vs_OFF
KIRP_ON_vs_OFF
KICH_ON_vs_OFF
LIHC_ON_vs_OFF
```

---

# 19. Cell-line model comparisons

Recommended cell-line analysis order:

```text
1 = 769P_ctrl
2 = 769P_sc1
3 = 769P_sc2
4 = 769P_sc12
5 = Huh7_ctrl
6 = Huh7_sc1
7 = Huh7_sc2
8 = Huh7_sc12
```

Recommended comparisons:

```matlab
comparisonPairs = [ ...
    1, 2; ...
    1, 3; ...
    1, 4; ...
    5, 6; ...
    5, 7; ...
    5, 8];
```

Comparison names:

```text
769P_ctrl_vs_sc1
769P_ctrl_vs_sc2
769P_ctrl_vs_sc12
Huh7_ctrl_vs_sc1
Huh7_ctrl_vs_sc2
Huh7_ctrl_vs_sc12
```

The configuration must explicitly map source file positions to these names.

---

# 20. Single-gene deletion analysis

The workflow evaluates:

```text
ATP maintenance
biomass production
```

Objectives:

```matlab
atpObjective = 'DM_atp_c_';
biomassObjective = 'biomass_reaction';
```

Function:

```matlab
singleGeneDeletion_rFASTCORMICS( ...
    model, ...
    'FBA', ...
    [], ...
    0, ...
    1);
```

Expected output table:

```text
geneList
grRatio_biomass
grRateKO_biomass
grRateWT_biomass
grRatio_ATP
grRateWT_ATP
grRateKO_ATP
```

Wild-type values should be expanded to vectors:

```matlab
grRateWT_biomass = ...
    ones(numel(grRatio_biomass), 1) * ...
    grRateWT_biomass;

grRateWT_ATP = ...
    ones(numel(grRatio_ATP), 1) * ...
    grRateWT_ATP;
```

---

# 21. Drug perturbation analysis

The drug list is created using:

```matlab
DrugList = unique( ...
    GeneDrugRelations.DrugName);
```

Run:

```matlab
[grRatio, grRateKO, grRateWT] = ...
    DrugDeletion( ...
        model, ...
        'FBA', ...
        DrugList);
```

A lethal drug is identified using:

```matlab
grRatio == 0
```

Expected outputs:

```matlab
DrugList
Drug_grRatio_biomass
Drug_grRateKO_biomass
Drug_grRateWT_biomass
lethalDrugs
```

TCGA and cell-line drug outputs should be stored separately.

---

# 22. Flux Variability Analysis

Run:

```matlab
[minFlux, maxFlux] = ...
    fluxVariability(model);
```

Reaction-level results should be aligned to:

```matlab
consistent_model.rxns
```

Use reaction identifiers:

```matlab
[reactionFound, referenceLocation] = ...
    ismember(model.rxns, consistent_model.rxns);
```

Do not assume identical row order between models.

---

# 23. FVA similarity analysis

Run:

```matlab
FVA_similarity_Thomas( ...
    firstMinFlux, ...
    firstMaxFlux, ...
    secondMinFlux, ...
    secondMaxFlux);
```

Subsystem names:

```matlab
uniSys = unique( ...
    consistent_model.subSystems);
```

Important:

Columns of subsystem-similarity matrices represent comparisons, not individual models.

Always save:

```matlab
comparisonPairs
comparisonNames
```

A single-model run cannot calculate between-model FVA similarity.

---

# 24. Recommended staged validation

Do not run the complete analysis immediately.

Recommended order:

```text
Stage 1: one-model FVA
Stage 2: one-model drug deletion
Stage 3: one-model single-gene deletion
Stage 4: complete one-model analysis
Stage 5: two-model FVA similarity
Stage 6: complete eight-model analysis
```

For a one-model test:

```matlab
config.runSingleGeneDeletion = false;
config.runDrugDeletion       = false;
config.runFVA                = true;
config.runFVASimilarity      = false;
```

For a two-model similarity test:

```matlab
config.runSingleGeneDeletion = false;
config.runDrugDeletion       = false;
config.runFVA                = true;
config.runFVASimilarity      = true;
```

Use:

```text
test_universal_analysis_one_model.m
```

for initial validation.

---

# 25. Flux sampling pipelines

Two separate sampling workflows are used.

Do not merge them.

## 25.1 RAW cell-line sampling

Script:

```text
RUN_RAW_1500_FluxSampling.m
```

Inputs:

```text
SamplingResults_medium_1500_model_*.mat
```

Uses:

```matlab
x.samples(:, 1:1500)
```

Generates RAW outputs ending in:

```text
_1500.xlsx
```

This workflow is used for:

```text
769-P
Huh7
```

---

## 25.2 TCGA ON/OFF sampling

Script:

```text
RUN_ONOFF_1000_FluxSampling.m
```

Inputs:

```text
samplingResults_mediumonoff_1000_model_*.mat
samplingResults_medium_KO_GLO1onoff_1000_model_*.mat
```

Uses:

```matlab
x.samples(:, 1:1000)
```

Generates ON/OFF outputs containing:

```text
5onoff
```

This workflow is used for:

```text
KIRC
KIRP
KICH
LIHC
```

Important fixes preserved in this script:

* RAW 1500 and ON/OFF 1000 remain separate;
* KO up-rank tables use `stats_KO`;
* the typo `KO5onoffonoff` was corrected to `KO5onoff`.

---

# 26. Sampling validation

Each sampling MAT file should contain:

```matlab
x.modelSampling
x.samples
```

Check:

```matlab
size(x.samples, 1) == ...
    numel(x.modelSampling.rxns)
```

TCGA and cell-line sampling files should be stored in separate folders:

```text
results/04_sampling/TCGA/
results/04_sampling/CELL_LINES/
```

Do not use identical filenames for different biological datasets in the same directory.

---

# 27. Flux-sum analysis

Flux-sum analysis quantifies metabolite turnover from sampled flux distributions.

The analysis includes:

* metabolite turnover calculation;
* comparison between conditions;
* Wilcoxon rank-sum testing;
* signal-to-noise calculation;
* pathway prioritization.

The manuscript analysis uses:

```text
500 sampled flux states per condition
```

Main outputs may include:

```text
flux-sum statistics tables
ranked metabolite tables
reaction statistics
pathway-specific output tables
```

Pathways of interest include:

```text
Pyruvate metabolism
Glycolysis
Pentose phosphate pathway
TCA cycle
Glutathione metabolism
ROS detoxification
```

---

# 28. Figure-generation scripts

Main figure scripts include:

```text
Generate_FigureS1_hi.m
Generate_FVA_TCGA.m
Generate_FVA_celllines.m
```

Recommended result-folder naming:

```text
Fig1d_fluxsum_glycolytic_metabolism
Fig1hi_FVA_similarity
Fig2ef_fluxsum_fluxsampling_pyruvate_metabolism
```

Avoid special characters such as:

```text
Log₂FC
```

in filenames.

Use:

```text
log2FC
```

inside file names when required.

---

# 29. IDARE and Cytoscape visualization

Significantly altered metabolic pathways are visualized using:

```text
IDARE
Cytoscape
```

Recommended folders:

```text
results/figures/IDARE_inputs/
results/figures/Cytoscape_sessions/
results/figures/final_figures/
```

Suggested files:

```text
IDARE_input_<comparison>.xlsx
Cytoscape_<comparison>.cys
```

---

# 30. Expected output organization

```text
results/
│
├── 01_prepared_data/
│   ├── TCGA/
│   └── CELL_LINES/
│
├── 02_models/
│   ├── reconstructed/
│   └── constrained/
│
├── 03_analysis/
│   ├── FBA/
│   ├── FVA/
│   ├── FVA_similarity/
│   ├── gene_deletion/
│   └── drug_deletion/
│
├── 04_sampling/
│   ├── TCGA/
│   └── CELL_LINES/
│
├── 05_flux_sum/
│   ├── TCGA/
│   └── CELL_LINES/
│
└── figures/
    ├── Fig1d_fluxsum_glycolytic_metabolism/
    ├── Fig1hi_FVA_similarity/
    ├── Fig2ef_fluxsum_fluxsampling_pyruvate_metabolism/
    ├── IDARE_inputs/
    └── Cytoscape_sessions/
```

---

# 31. Complete TCGA execution order

Set:

```matlab
cfg.runMode = 'TCGA';
```

Run:

```matlab
cfg = pipeline_config;

run('01_geo_gse62944_generate_TCGA_tables.m')
run('02_prepare_expression_data.m')
run('03_build_fastcormics_models.m')

setMediumConstraints_TCGA

run('04_run_model_analysis.m')

Generate_FVA_TCGA

RUN_ONOFF_1000_FluxSampling
```

Then run the TCGA flux-sum and figure-generation scripts.

---

# 32. Complete cell-line execution order

Set:

```matlab
cfg.runMode = 'CELL_LINES';
```

Run:

```matlab
cfg = pipeline_config;

run('02_prepare_expression_data.m')
run('03_build_fastcormics_models.m')

setMediumConstraints_CellLines

run('04_run_model_analysis.m')

Generate_FVA_celllines

RUN_RAW_1500_FluxSampling
```

Then run the cell-line flux-sum and figure-generation scripts.

---

# 33. Run all datasets

Set:

```matlab
cfg.runMode = 'ALL';
```

Run the modules sequentially.

The TCGA HIGH/LOW extraction must be completed before expression preparation.

The RAW 1500 and ON/OFF 1000 sampling workflows must remain separate.

---

# 34. Common errors

## Input file not found

Check:

```matlab
cfg = pipeline_config;
dir(cfg.dataDir)
```

Confirm exact filenames and extensions.

---

## Out of memory when loading GSE62944

Use the pre-extracted cancer-specific files:

```text
KICH.txt
KIRC.txt
KIRP.txt
LIHC.txt
```

Do not repeatedly load the complete tumor matrix in MATLAB R2019b.

---

## COBRA function not found

Run:

```matlab
initCobraToolbox(false)
```

Check:

```matlab
which optimizeCbModel
```

---

## Solver error

Run:

```matlab
changeCobraSolver('ibm_cplex', 'all')
```

or another installed solver.

Check:

```matlab
getCobraSolver
```

---

## Biomass objective not found

Check:

```matlab
find(strcmp(model.rxns, ...
    'biomass_reaction'))
```

Inspect similar reaction identifiers:

```matlab
model.rxns(contains( ...
    lower(model.rxns), ...
    'biomass'))
```

---

## ATP objective not found

Check:

```matlab
find(strcmp(model.rxns, ...
    'DM_atp_c_'))
```

---

## Model infeasible after medium constraints

Review:

* medium reaction identifiers;
* missing nutrient exchanges;
* oxygen bounds;
* carbon-dioxide bounds;
* blocked exchanges;
* biomass objective.

Inspect:

```matlab
R.missingReactions
R.FBAResults
```

---

## FVA similarity with one model

FVA similarity requires at least two models.

Use:

```matlab
config.runFVASimilarity = false;
```

for one-model tests.

---

## Incorrect model order

Check the source MAT-file order.

Do not assume that:

```text
element 3 = sc2
element 4 = sc12
```

because some original cell arrays use:

```text
element 3 = sc12
element 4 = sc2
```

---

## Inconsistent gene-deletion table

All columns must have equal lengths.

Expand scalar wild-type rates before creating the table.

---

# 35. Biological validation

The pipeline validates technical consistency, but it does not replace biological validation.

Before interpreting results, confirm:

* required objectives exist;
* models are feasible;
* medium constraints are biologically appropriate;
* gene associations are present;
* gene-drug relations are compatible;
* reaction identifiers map to the reference model;
* model comparisons are correctly defined;
* output differences are biologically plausible.

---

# 36. Reproducibility information

Record the following for each release:

```text
MATLAB version
COBRA Toolbox version
rFASTCORMICS version
solver name and version
operating system
Recon3D version
input dataset versions
model reconstruction parameters
medium definitions
objective reactions
number of sampling points
random seed, when applicable
Git commit
Zenodo release DOI
```

Recommended file:

```text
VERSIONS.txt
```

Example:

```text
MATLAB: R2019b
COBRA Toolbox: <version or commit>
rFASTCORMICS: <version or commit>
Solver: IBM ILOG CPLEX <version>
Reference model: Recon3D <version>
Repository release: v1.0.0
```

---

# 37. Files to include in the repository

Include:

```text
README.md
docs/PIPELINE_USAGE.md
pipeline_config.m
all MATLAB source files
small configuration files
medium definitions
model metadata
comparison metadata
VERSIONS.txt
CITATION.cff
LICENSE
```

Do not include:

```text
personal access tokens
solver license files
machine-specific absolute paths
temporary Excel files
MATLAB autosave files
unrelated workspace dumps
private reviewer links
```

Raw datasets may be referenced through GEO instead of being duplicated in GitHub.

---

# 38. Final reproducibility checklist

Before publishing the repository:

* [ ] all scripts use relative paths;
* [ ] `pipeline_config.m` works on a clean clone;
* [ ] GEO accession numbers are documented;
* [ ] all required input filenames are listed;
* [ ] MATLAB and dependency versions are recorded;
* [ ] medium assignments are documented;
* [ ] TCGA and cell-line workflows run separately;
* [ ] RAW 1500 and ON/OFF 1000 sampling remain separate;
* [ ] model order is documented;
* [ ] analysis outputs are saved in dataset-specific folders;
* [ ] figure-source folders match manuscript figure panels;
* [ ] reviewer tokens and private links are removed;
* [ ] a clean-clone test has been completed;
* [ ] the final GitHub release is archived in Zenodo.

---

# 39. Summary

The complete reproducibility workflow is:

```text
GEO RNA-seq data
        |
        v
TCGA cancer extraction
        |
        v
AKR1A1 HIGH / LOW stratification
        |
        v
Expression preprocessing
        |
        v
Expression discretization
        |
        v
rFASTCORMICS reconstruction
        |
        v
RPMI / DMEM medium constraints
        |
        v
FBA / FVA / FVA similarity
        |
        v
Gene and drug perturbation analyses
        |
        v
RAW or ON/OFF flux sampling
        |
        v
Flux-sum statistics
        |
        v
Figure generation
        |
        v
IDARE / Cytoscape visualization
```

For TCGA, begin with the GSE62944 extraction and AKR1A1 HIGH/LOW split.

For cell lines, begin with expression preparation using GSE310784 and GSE310828.

Keep TCGA and cell-line result files in separate directories throughout the workflow.
