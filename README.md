# Metabolomics Re-analysis of Early Pregnancy Loss

## Overview

This repository contains the code and analysis workflows used to reproduce and re-run a non-targeted metabolomics and environmental exposure study investigating early pregnancy loss (EPL) and recurrent pregnancy loss (RPL).

The original study used liquid chromatography coupled with quadrupole time-of-flight high-resolution mass spectrometry (LC-QTOF HRMS) to profile maternal plasma samples collected during the first trimester of pregnancy. The goal was to identify metabolic alterations and environmental chemical exposures associated with pregnancy loss.

This work contributed to the publication:

**Ji X, Yilmaz BD, Edwards J, Nagasuri A, et al.**
*A case-control design of environmental exposures and metabolic alterations in early pregnancy loss using non-targeted analysis.*
Environmental Chemistry and Ecotoxicology (2026)

## My Role

I did not develop the original metabolomics analysis pipeline.

My contribution focused on learning, adapting, and re-running the existing workflow after sample classifications were updated during manuscript preparation.

Specifically, I:

- Learned and maintained an existing metabolomics analysis pipeline
- Updated scripts to accommodate revised sample annotations
- Re-ran analyses following sample readjudication
- Regenerated PCA visualizations and exploratory analyses
- Reproduced intermediate statistical outputs
- Performed quality-control checks on intermediate results
- Verified consistency of outputs throughout the workflow
- Documented the analysis process in a reproducible runbook
- Contributed analyses supporting the final publication


## Study Background

Early pregnancy loss (EPL) affects approximately 15% of clinically recognized pregnancies, while recurrent pregnancy loss (RPL) affects approximately 5% of reproductive-age women. Although known risk factors include chromosomal abnormalities, endocrine disorders, and anatomical abnormalities, the role of environmental exposures and maternal metabolism remains incompletely understood.

To better characterize these relationships, this study applied non-targeted analysis (NTA) to maternal plasma samples collected during the first trimester of pregnancy. By simultaneously measuring endogenous metabolites and exogenous environmental chemicals, the study aimed to identify biological pathways and environmental exposures associated with pregnancy loss.


## Study Cohort

The study included:

- 45 Early Pregnancy Loss (EPL) cases
- 33 Recurrent Pregnancy Loss (RPL) cases
- 12 First Pregnancy Loss cases
- 23 Ongoing Pregnancy controls

All participants had normal uterine anatomy and normal endocrine testing. Plasma samples were collected between 6 and 13 weeks of gestation.


## Data Generation

### Platform

- Liquid Chromatography (LC)
- Electrospray Ionization (ESI)
- Quadrupole Time-of-Flight High-Resolution Mass Spectrometry (QTOF-HRMS)

### Data Processing

Raw mass spectrometry files were processed using:

- MS-DIAL
- Peak extraction
- Peak alignment
- Chemical annotation
- Detection-frequency filtering
- Truncated-normal missing value imputation
- ComBat batch correction

Following preprocessing:

- 14,419 positive-ionization (ESI+) features were retained
- 10,823 negative-ionization (ESI−) features were retained


## Analysis Workflow

### Quality Control

- Technical replicate evaluation
- Detection-frequency filtering
- Missing-value assessment
- Batch-effect assessment
- Batch correction using ComBat

### Exploratory Analysis

- Principal Component Analysis (PCA)
- Heatmap generation
- Correlation analyses
- Feature clustering

### Statistical Analysis

- Differential abundance testing
- Fold-change analysis
- Mann-Whitney U testing
- Multiple-testing correction (Benjamini-Hochberg)

### Biological Interpretation

- Metabolite annotation
- Pathway enrichment analysis
- Metabolite-exposure correlation networks
- XGBoost feature prioritization
- SHAP interpretation


## Key Findings from the Published Study

### Strong Metabolic Separation Between Cases and Controls

More than 2,900 metabolomic features differed significantly between EPL and control samples.

Unsupervised analyses demonstrated clear separation between pregnancy loss and control metabolomic profiles.

### Altered Steroid Hormone Metabolism

Pregnancy loss samples exhibited reduced abundances of:

- Progesterone
- Corticosterone
- Additional steroid hormone intermediates

These findings suggest disruption of steroidogenesis and hormonal support pathways important for early pregnancy maintenance.

### Evidence of Altered Energy Metabolism

Pregnancy loss samples demonstrated elevated medium-chain acylcarnitines, including:

- Octanoylcarnitine (C8)
- O-decanoylcarnitine (C10)
- Cis-4-decenoylcarnitine (C10:1)

These metabolites are involved in mitochondrial fatty-acid β-oxidation and may indicate altered energy metabolism.

### Lipid Metabolism Disruption

Several glycerophospholipids and lipid-derived metabolites were significantly altered in EPL and RPL samples, suggesting disruption of membrane remodeling and fatty-acid metabolism.

### Environmental Chemical Associations

The study identified several environmental contaminants associated with pregnancy loss, including:

- Acetyl tributyl citrate
- Deethylatrazine
- Dibutyl phthalate

Correlation analyses suggested potential interactions between environmental exposures and endogenous metabolic pathways.

### Machine Learning-Based Feature Prioritization

XGBoost models successfully distinguished EPL from control samples and identified progesterone, lipid-derived metabolites, and additional metabolic features as major contributors to classification performance.


## Repository Contents

This repository contains scripts used during the re-analysis process, including:

- Data preprocessing workflows
- Sample annotation processing
- Missing-value handling
- Batch-correction analyses
- PCA visualizations
- Heatmap generation
- Correlation analyses
- Differential abundance analyses
- Quality-control workflows
- Reproducibility documentation

The repository primarily captures the stages of the workflow that were re-executed and validated following sample reclassification.


## Skills Demonstrated

- Metabolomics
- Exposomics
- Scientific Reproducibility
- High-Dimensional Data Analysis
- Quality Control
- Batch Effect Correction
- Principal Component Analysis
- Statistical Testing
- Data Visualization
- Python
- R
- Scientific Documentation


## Publication

Ji X, Yilmaz BD, Edwards J, Nagasuri A, et al.

**A case-control design of environmental exposures and metabolic alterations in early pregnancy loss using non-targeted analysis**

https://www.sciencedirect.com/science/article/pii/S2590182626000779

*Environmental Chemistry and Ecotoxicology* (2026)


## Author

**Amrita Nagasuri**

M.S. Health Data Science  
University of California, San Francisco
