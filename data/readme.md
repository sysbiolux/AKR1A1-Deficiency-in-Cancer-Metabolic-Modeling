## Data

The repository uses publicly available bulk RNA-seq datasets from GEO together with the Recon3D human metabolic reconstruction.

### TCGA patient tumors

TCGA RNA-seq data were obtained from the GEO SuperSeries:

- **GSE62944**
  - File used: `GSE62944_RAW.tar`
  - NCBI GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE62944 :contentReference[oaicite:0]{index=0}

Tumor cohorts included:

- KIRC — Kidney Renal Clear Cell Carcinoma
- KIRP — Kidney Renal Papillary Cell Carcinoma
- KICH — Kidney Chromophobe
- LIHC — Liver Hepatocellular Carcinoma

Patients were stratified according to **AKR1A1** expression into:

- AKR1A1 HIGH (ON) — upper quartile
- AKR1A1 LOW (OFF) — lower quartile

The script

```
split_TCGA_high_low.m
```

creates the ON/OFF patient groups used throughout the metabolic reconstruction workflow.

---

### Cancer cell-line RNA-seq

Two experimental RNA-seq datasets were used:

- **GSE310784**
  - 769-P renal carcinoma cell line (RPMI medium)

- **GSE310828**
  - Huh7 hepatocellular carcinoma cell line (DMEM medium)

NCBI GEO:
- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE310784 :contentReference[oaicite:1]{index=1}
- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE310828 :contentReference[oaicite:2]{index=2}

For each cell line, four metabolic models were reconstructed:

| Cell line | Condition | Model |
|-----------|-----------|-------|
| 769-P | Control | `model_7_ctrl` |
| 769-P | siAKR1A1_1 | `model_7_sc1` |
| 769-P | pooled siAKR1A1_1 + siAKR1A1_2 | `model_7_sc12` |
| 769-P | siAKR1A1_2 | `model_7_sc2` |
| Huh7 | Control | `model_H_ctrl` |
| Huh7 | siAKR1A1_1 | `model_H_sc1` |
| Huh7 | pooled siAKR1A1_1 + siAKR1A1_2 | `model_H_sc12` |
| Huh7 | siAKR1A1_2 | `model_H_sc2` |

---

### Experimental media

The metabolic reconstructions use the experimental culture media employed in the biological experiments:

- **RPMI** for:
  - 769-P cell-line models
  - KIRC, KIRP and KICH TCGA models

- **DMEM** for:
  - Huh7 cell-line models
  - LIHC TCGA models

Medium compositions are provided in the repository and are applied before FBA, FVA, and flux sampling.
