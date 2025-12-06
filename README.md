# ## Usage
1. Place raw UKB data in `data/raw/` (follow UKB data policy).
2. Run scripts sequentially in R:
   - 00_load_packages.R
   - 01_load_UKB_raw.R
   - 02_preprocess_covariates.R
   - 03_define_outcomes.R
   - 04_imputation_mice.R
   - 05_cox_models.R
   - 06_baseline_table.R
3. Results will be generated in `results/`.

All code includes English comments for reproducibility.
PTEN-UKB-analysis
