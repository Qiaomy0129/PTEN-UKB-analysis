PTEN Rare Variants and Cancer Risk Analysis
Project Overview
This project analyzes the association between rare PTEN gene variants and multiple cancer risks using UK Biobank data. Through Cox proportional hazards regression, baseline characteristic analysis, sensitivity analysis, and stratified analysis, we assess the impact of PTEN rare variants on cancer incidence risk.

Project Structure
text
PTEN_Cancer_Analysis/
├── PTEN_Analysis_Full.R          # Complete analysis script
├── README.md                     # Project documentation (this file)
├── requirements.txt              # R package dependencies
├── config.R                      # Configuration file (paths and parameters)
├── results/                      # Results output directory
│   ├── cancer_model1_cox_rarenon.csv
│   ├── cancer_model2_cox_rarenon.csv
│   ├── PTEN_chi_baseline_rare_non.txt
│   ├── PTEN_chi_cancer_rare_non.txt
│   ├── cancer_model2_cox_exclude_deaths_2y.csv
│   ├── cancer_model2_cox_complete_cases.csv
│   └── session_info.txt
└── docs/                         # Documentation directory
    └── analysis_protocol.md      # Detailed analysis protocol
Quick Start
1. Environment Setup
Install R package dependencies:
r
# Method 1: Install from requirements.txt
install.packages(c("tidyverse", "ukbwranglr", "dplyr", "mice", 
                   "survival", "gttable1", "epiDisplay", "car", 
                   "AF", "rio", "forcats"))

# Method 2: Run installation script
source("install_dependencies.R")
Configure data paths:
Edit the config.R file to update the following paths:

r
# UKB data paths
UKB_MAIN_PATH <- '~/ukbdata/Main_Dataset/'

# Cancer code file paths
CANCER_CODE_PATH <- '~/share/qiaomy/Cancer/cancer.xlsx'

# PTEN variant data paths
PTEN_SAMPLE_PATH <- '~/ukbdata/WES/extractAF/PTEN_UKB_sampleDetail.txt.gz'
PTEN_AF_PATH <- '~/ukbdata/WES/extractAF/PTEN_UKB_AF.vcf'
2. Run Analysis
r
# Run complete analysis pipeline
source("PTEN_Analysis_Full.R")

# Or run step-by-step
source("config.R")
source("01_data_preparation.R")
source("02_pten_annotation.R")
source("03_cox_analysis.R")
source("04_sensitivity_analysis.R")
Analysis Pipeline
Phase 1: Data Preparation
Data Dictionary Creation: Filter UK Biobank variables

Clinical Event Extraction: Extract disease phenotypes from ICD codes

Covariate Processing: Blood pressure, BMI, demographic characteristics, etc.

Phase 2: PTEN Variant Annotation
Variant Frequency Calculation: Define rare variants based on allele frequency (AF < 0.01)

Sample Selection: Keep White ethnicity, exclude CNV-positive individuals

Case-Control Matching: PTEN rare variant carriers vs. non-carriers

Phase 3: Main Analysis
Cox Model 1: Minimal adjustment (age + sex)

Cox Model 2: Full adjustment (+ genetic principal components + lifestyle + comorbidities)

Baseline Characteristics: Generate baseline tables using gttable1

Phase 4: Sensitivity Analysis
Exclude Deaths Within 2 Years: Reduce reverse causation bias

Complete Case Analysis: Exclude missing covariates

Multiple Imputation: Handle missing data using MICE

Phase 5: Subgroup Analysis
Stratified Analysis: Stratify by sex, age, BMI, comorbidities, etc.

Interaction Testing: Assess effect modification

Analyzed Cancer Types
This project analyzes 18 cancer types:

Common Cancers: Breast, Prostate, Colorectal, Lung

Digestive System: Stomach, Oesophageal, Pancreas, Liver

Reproductive System: Cervical, Uterine, Ovarian, Testicular

Others: Bladder, Kidney, Thyroid, Laryngeal, Tongue, Brain

 Technical Details
Main R Packages
ukbwranglr: UK Biobank data processing

survival: Survival analysis

mice: Multiple imputation

gttable1: Baseline table generation

tidyverse: Data processing and visualization

Statistical Methods
Cox Proportional Hazards Models: Estimate hazard ratios (HR) and 95% confidence intervals

Multiple Testing Correction: Account for 18 cancer types

Genetic Confounding Control: Adjust for first 10 genetic principal components

Time Scale: From baseline visit to cancer diagnosis or censoring

Quality Control
Sample Exclusion: Non-White ethnicity, CNV-positive, pre-baseline diagnoses

Missing Data Handling: Multiple imputation vs. complete case analysis

Model Assumption Testing: Proportional hazards assumption, collinearity diagnostics

Expected Outputs
Main Result Files
cancer_model1_cox_rarenon.csv: Minimal adjustment Cox model results

cancer_model2_cox_rarenon.csv: Fully adjusted Cox model results

PTEN_chi_baseline_rare_non.txt: Case-control baseline characteristics comparison

PTEN_chi_cancer_rare_non.txt: Cancer incidence comparison

Sensitivity Analysis Results
cancer_model2_cox_exclude_deaths_2y.csv: Excluding deaths within 2 years

cancer_model2_cox_complete_cases.csv: Complete case analysis

 Custom Configuration
Modify Analysis Parameters
Edit in config.R:

r
# Rare variant threshold
RARE_AF_THRESHOLD <- 0.01  # or change to 0.005, etc.

# Censoring date
CENSORING_DATE <- as.Date("2023-12-31")  # Extend follow-up time

# Death exclusion window
DEATH_EXCLUSION_DAYS <- 365  # Change to 1 year
Add New Cancer Types
Add new cancer codes in the Cancer worksheet of cancer.xlsx

Add cancer name to the cancer_types vector in the script

Citations and Acknowledgments
Data Sources
UK Biobank:

PTEN variant data: From UK Biobank WES
