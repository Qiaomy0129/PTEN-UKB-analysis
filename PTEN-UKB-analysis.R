# PTEN Rare Variants and Cancer Risk Analysis - Complete Code


```r
# ============================================================================
# PTEN Rare Variants and Cancer Risk Analysis - Complete Script
# Analysis of PTEN rare variants and cancer associations in UK Biobank
# ============================================================================

# --------------------------------
# 1. LOAD PACKAGES & SETUP
# --------------------------------

# Load required libraries
library(tidyverse)
library(ukbwranglr)
library(dplyr)
library(mice)
library(survival)
library(gttable1)
library(epiDisplay)
library(car)
library(AF)
library(rio)
library(forcats)

# Set working directory
setwd("~/share/qiaomy/rare-non")

# --------------------------------
# 2. CONFIGURATION
# --------------------------------

# File paths configuration
ukb_main_path      <- '~/ukbdata/Main_Dataset/'
ukb_data_main_path <- sprintf('%s/k96511.tsv', ukb_main_path)
ukb_data_dict_path <- sprintf('%s/Data_Dictionary_Showcase.tsv', ukb_main_path)
ukb_codings_path   <- sprintf('%s/Codings_Showcase.tsv', ukb_main_path)

# Cancer code paths
cancer_code_path <- '~/share/qiaomy/Cancer/cancer.xlsx'
confounder_sheet <- 'confounder'
cancer_sheet <- 'Cancer'

# PTEN variant paths
pten_sample_path <- '~/ukbdata/WES/extractAF/PTEN_UKB_sampleDetail.txt.gz'
pten_af_path <- '~/ukbdata/WES/extractAF/PTEN_UKB_AF.vcf'
phenotype_table_path <- '~/share/phenotype_table.Rdata'

# Analysis parameters
rare_af_threshold <- 0.01
censoring_date <- as.Date("2021-04-25")
death_exclusion_days <- 730

# Create output directory
if (!dir.exists("results")) {
  dir.create("results")
}

# --------------------------------
# 3. CREATE DATA DICTIONARY
# --------------------------------

cat("Step 1: Creating data dictionary...\n")

# Load metadata
ukb_data_dict <- get_ukb_data_dict(path = ukb_data_dict_path)
ukb_codings <- get_ukb_codings(path = ukb_codings_path)

# Create data dictionary
data_dict <- make_data_dict(ukb_data_main_path, ukb_data_dict = ukb_data_dict)

# Select variables of interest
data_dict_selected <- data_dict %>%
  dplyr::filter(FieldID %in% c(
    "eid", # Participant ID
    '53', # Date of attending assessment centre
    '21022', # Age at recruitment
    '31', # Gender
    '34', # Year of birth
    '21000', # Race
    '20118', # Area
    '21001', # BMI
    '6138', # Education
    '100064', # Employment
    '738', # Income
    '20116', # Smoking status
    "20117", # Drink frequency
    '4080', # Systolic blood pressure
    '4079', # Diastolic blood pressure
    '22006', # Genetic ethnic
    '22020', # Genetic principal components (used flag)
    '22009', # Genetic principal components
    '40001','40000', # Primary death ICD-10
    '40002', # Secondary death ICD-10
    '41270','41280', # Summary HES ICD-10
    '20003', # Self-report medication
    '20002','20008', # Self-report non-cancer
    '30740', # Glucose
    '30750', # HbA1c
    '845', # Age completed education
    '20001','20006', # Self-report cancer
    '20004','20010', # Self-report operation
    '40013','40005', # Cancer register ICD-9
    '40006', # Cancer register ICD-10
    '41271','41281', # Summary HES ICD-9
    '41273','41283', # Summary HES OPCS-3
    '41272','41282' # Summary HES OPCS-4
  ))

# --------------------------------
# 4. LOAD UKB MAIN DATA
# --------------------------------

cat("Step 2: Loading UK Biobank main data...\n")

# Read selected variables from UKB main dataset
ukb_main <- read_ukb(path = ukb_data_main_path,
                     data_dict = data_dict_selected,
                     ukb_data_dict = ukb_data_dict,
                     ukb_codings = ukb_codings)
ukb_main <- as_tibble(ukb_main)

# --------------------------------
# 5. EXTRACT CLINICAL EVENTS
# --------------------------------

cat("Step 3: Extracting clinical events...\n")

# Extract clinical events (cancer, death, etc.)
clinical_events <- ukb_main %>%
  tidy_clinical_events(
    ukb_data_dict = ukb_data_dict,
    ukb_codings = ukb_codings,
    clinical_events_sources = c('summary_hes_icd10','cancer_register_icd9',
                                'cancer_register_icd10','primary_death_icd10')
  ) %>%
  dplyr::bind_rows()

# --------------------------------
# 6. LOAD DISEASE CODE LISTS
# --------------------------------

cat("Step 4: Loading disease code lists...\n")

# Load confounder codes
confounder_codes <- rio::import(cancer_code_path, sheet = confounder_sheet) %>%
  mutate(code = str_replace(code,'\\.','')) %>% 
  mutate_all(str_trim) %>% 
  as_tibble()

# Load cancer codes
disease_codes <- rio::import(cancer_code_path, sheet = cancer_sheet) %>%
  mutate(code = str_replace(code,'\\.','')) %>% 
  mutate_all(str_trim) %>% 
  as_tibble()

# --------------------------------
# 7. EXTRACT PHENOTYPES
# --------------------------------

cat("Step 5: Extracting phenotypes...\n")

# Extract confounder phenotypes
phe_confounder <- extract_phenotypes(clinical_events = clinical_events,
                                     clinical_codes = confounder_codes)

# Extract cancer phenotypes
phe_disease <- extract_phenotypes(clinical_events = clinical_events,
                                  clinical_codes = disease_codes)

# Filter only cancer registry sources
rows_to_keep <- phe_disease$source %in% c("f40006", "f40013")
phe_disease <- phe_disease[rows_to_keep, ]

# --------------------------------
# 8. EXTRACT EARLIEST DISEASE DATES
# --------------------------------

cat("Step 6: Extracting earliest disease dates...\n")

# Extract earliest date for each confounder
date_confounder <- phe_confounder %>%
  group_by(eid, disease) %>%
  summarise(min_date = min(date, na.rm = TRUE)) %>%
  pivot_wider(names_from = disease, values_from = min_date)

# Extract earliest date for each cancer
date_disease <- phe_disease %>%
  group_by(eid, disease) %>%
  summarise(min_date = min(date, na.rm = TRUE)) %>%
  pivot_wider(names_from = disease, values_from = min_date)

# --------------------------------
# 9. FILTER AND RENAME VARIABLES
# --------------------------------

cat("Step 7: Filtering and renaming variables...\n")

# Filter columns (keep only first measurement)
f <- grep(".*_0_|eid", colnames(ukb_main), value = TRUE)
f2 <- grep("f40002|f20003|f20002|f20008|f20001|f20006|f20004|f20010|f40013|f40005|f40006|f41271|f41281|f41270|f41280|f41273|f41283|f41272|f41282",
           f, invert = TRUE, value = TRUE)
f3 <- grep("f6138_0_[1-5]|f6154_0_[1-5]|f6177_0_[1-2]|f22009_0_(1[1-9]|2[0-9]|3[0-9]|40)",
           f2, invert = TRUE, value = TRUE)

ukb_finally <- ukb_main[, f3]

# Rename columns
ukb_finally <- rename(ukb_finally,
                      eid = eid,
                      sex = sex_f31_0_0,
                      birth_year = year_of_birth_f34_0_0,
                      date_attending = date_of_attending_assessment_centre_f53_0_0,
                      income = average_total_household_income_before_tax_f738_0_0,
                      qualifications = qualifications_f6138_0_0,
                      education = age_completed_full_time_education_f845_0_0,
                      DBP_1 = diastolic_blood_pressure_automated_reading_f4079_0_0,
                      DBP_2 = diastolic_blood_pressure_automated_reading_f4079_0_1,
                      SBP_1 = systolic_blood_pressure_automated_reading_f4080_0_0,
                      SBP_2 = systolic_blood_pressure_automated_reading_f4080_0_1,
                      smoke = smoking_status_f20116_0_0,
                      alcohol = alcohol_drinker_status_f20117_0_0,
                      race = ethnic_background_f21000_0_0,
                      BMI = body_mass_index_bmi_f21001_0_0,
                      age = age_at_recruitment_f21022_0_0,
                      genetic_ethnic = genetic_ethnic_grouping_f22006_0_0,
                      glu = glucose_f30740_0_0,
                      HbA1c = glycated_haemoglobin_hba1c_f30750_0_0,
                      used_in_genetic = used_in_genetic_principal_components_f22020_0_0,
                      area = home_area_population_density_urban_or_rural_f20118_0_0,
                      genetic_01 = genetic_principal_components_f22009_0_1,
                      genetic_02 = genetic_principal_components_f22009_0_2,
                      genetic_03 = genetic_principal_components_f22009_0_3,
                      genetic_04 = genetic_principal_components_f22009_0_4,
                      genetic_05 = genetic_principal_components_f22009_0_5,
                      genetic_06 = genetic_principal_components_f22009_0_6,
                      genetic_07 = genetic_principal_components_f22009_0_7,
                      genetic_08 = genetic_principal_components_f22009_0_8,
                      genetic_09 = genetic_principal_components_f22009_0_9,
                      genetic_10 = genetic_principal_components_f22009_0_10,
                      death_date = date_of_death_f40000_0_0,
                      cause_of_death = underlying_primary_cause_of_death_icd10_f40001_0_0)

# --------------------------------
# 10. DEFINE COVARIATES & DISEASE STATUS
# --------------------------------

cat("Step 8: Defining covariates and disease status...\n")

# Calculate average blood pressure
ukb_finally$SBP <- (ukb_finally$SBP_1 + ukb_finally$SBP_2) / 2
ukb_finally$DBP <- (ukb_finally$DBP_1 + ukb_finally$DBP_2) / 2

# Define hypertension
ukb_finally$hypertension <- ifelse(
  ukb_finally$eid %in% date_confounder$eid[!is.na(date_confounder$hypertension)], 'Yes',
  ifelse((!is.na(ukb_finally$SBP) & ukb_finally$SBP >= 140) | 
         (!is.na(ukb_finally$DBP) & ukb_finally$DBP >= 90), 'Yes', 'No')
)

# Define other diseases
ukb_finally$CHD <- ifelse(ukb_finally$eid %in% date_confounder$eid[!is.na(date_confounder$CHD)], 'Yes', 'No')
ukb_finally$DM <- ifelse(ukb_finally$eid %in% date_confounder$eid[!is.na(date_confounder$DM)], 'Yes',
                         ifelse((!is.na(ukb_finally$glu) & ukb_finally$glu >= 7) | 
                                (!is.na(ukb_finally$HbA1c) & ukb_finally$HbA1c >= 48), 'Yes', 'No'))
ukb_finally$stroke <- ifelse(ukb_finally$eid %in% date_confounder$eid[!is.na(date_confounder$stroke)], 'Yes', 'No')
ukb_finally$DD <- ifelse(ukb_finally$eid %in% date_confounder$eid[!is.na(date_confounder$DD)], 'Yes', 'No')

# Define cancer status
cancer_list <- c("Breast", "Prostate", "Colorectal", "Stomach", "Cervical", "Uterine", "Testicular", 
                 "Ovarian", "Bladder", "Kidney", "Pancreas", "Thyroid", "Laryngeal", "Tongue", 
                 "oesophageal", "Brain", "Liver", "Lung")

for (cancer in cancer_list) {
  ukb_finally[[cancer]] <- ifelse(ukb_finally$eid %in% date_disease$eid[!is.na(date_disease[[cancer]])], 'Yes', 'No')
}

# Join with cancer dates
ukb_finally <- left_join(ukb_finally, date_disease, by = 'eid')

# --------------------------------
# 11. RECODE COVARIATES
# --------------------------------

cat("Step 9: Recoding covariates...\n")

# Recode income
levels(ukb_finally$income)[c(1,2)] <- NA
ukb_finally$income <- ifelse(ukb_finally$income %in% c('"Less than 18,000"','"18,000 to 30,999"'), 
                             "Less than 31000", 
                             ifelse(ukb_finally$income %in% c('"31,000 to 51,999"', '"52,000 to 100,000"', '"Greater than 100,000"'),
                                    "31000 and above", NA))

# Recode race
levels(ukb_finally$race)[c(1,2)] <- NA
ukb_finally$race <- ifelse(is.na(ukb_finally$race), ukb_finally$race,
                           ifelse(ukb_finally$race %in% c("White", 'British', "Any other white background"),
                                  "White", 'other'))

# Recode education
levels(ukb_finally$qualifications)[c(1,2)] <- NA
ukb_finally$qualifications <- ifelse(ukb_finally$qualifications == 'College or University degree',
                                     'With college/university degree',
                                     'Without college/university degree')

ukb_finally <- ukb_finally %>%
  mutate(education = ifelse(!is.na(education), "Yes", "No"))

# Recode smoking, alcohol, area
levels(ukb_finally$smoke)[1] <- NA
levels(ukb_finally$alcohol)[1] <- NA
levels(ukb_finally$area)[9] <- NA

ukb_finally <- ukb_finally %>%
  mutate(area = case_when(
    area %in% c("England/Wales - Town and Fringe - sparse", 
                "England/Wales - Village - sparse", 
                "England/Wales - Hamlet and Isolated dwelling - sparse", 
                "England/Wales - Town and Fringe - less sparse", 
                "England/Wales - Village - less sparse", 
                "England/Wales - Hamlet and Isolated Dwelling - less sparse", 
                "Scotland - Accessible Small Town", 
                "Scotland - Remote Small Town", 
                "Scotland - Very Remote Small Town", 
                "Scotland - Accessible Rural", 
                "Scotland - Remote Rural", 
                "Scotland - Very Remote Rural") ~ "Rural",
    area %in% c("England/Wales - Urban - less sparse", 
                "Scotland - Large Urban Area", 
                "Scotland - Other Urban Area",
                "England/Wales - Urban - sparse") ~ "Urban",
    TRUE ~ area))

# --------------------------------
# 12. LOAD PTEN VARIANT DATA
# --------------------------------

cat("Step 10: Loading PTEN variant data...\n")

# Load PTEN sample list and allele frequency data
PTEN_samplelist <- read.delim(pten_sample_path)
colnames(PTEN_samplelist)[colnames(PTEN_samplelist) == "Sample"] <- "eid"

PTEN_UKB_AF <- read.delim(pten_af_path)

load(phenotype_table_path)

# Annotate rare vs. common variants
PTEN_UKBAF_note <- subset(PTEN_UKB_AF, select = c(3, 11, 13, 14))
PTEN_UKBAF_note$rare_non <- ifelse(PTEN_UKBAF_note$AF < rare_af_threshold, "case", "control")
colnames(PTEN_UKBAF_note)[colnames(PTEN_UKBAF_note) == "ID"] <- "Site"

cat("PTEN variant counts:\n")
print(table(PTEN_UKBAF_note$rare_non))

# Merge sample list with variant annotation
PTEN_samplelist_rarenote <- left_join(PTEN_samplelist, PTEN_UKBAF_note, by = "Site")
PTEN_samplelist_rarenote <- subset(PTEN_samplelist_rarenote, select = c("eid", "rare_non"))

# Filter duplicates: keep "case" if any, else keep "control"
filtered_PTEN_samplelist_rarenote <- PTEN_samplelist_rarenote[complete.cases(PTEN_samplelist_rarenote$rare_non), ]
filtered_PTEN_samplelist_rarenote <- filtered_PTEN_samplelist_rarenote[order(filtered_PTEN_samplelist_rarenote$eid), ]
filtered_PTEN_samplelist_rarenote <- filtered_PTEN_samplelist_rarenote %>%
  group_by(eid) %>%
  filter(if (any(rare_non == "case")) any(rare_non == "case") else all(rare_non == "control")) %>%
  distinct(eid, .keep_all = TRUE)

# --------------------------------
# 13. MERGE PHENOTYPE AND PTEN DATA
# --------------------------------

cat("Step 11: Merging phenotype and PTEN data...\n")

# Merge UKB data with PTEN variant annotation
ukb_finally1 <- left_join(ukb_finally, filtered_PTEN_samplelist_rarenote, by = "eid")

# Merge with CNV data from phenotype table
ukb_finally1 <- phenotype_table %>%
  select(eid, cnv) %>%
  full_join(ukb_finally1, by = "eid")

ukb_finally1 <- ukb_finally1 %>%
  select(-cnv, everything(), cnv)

# Recode rare_non
ukb_finally1 <- mutate(ukb_finally1, rare_non = case_when(
  rare_non == "control" ~ 0,
  rare_non == "case" ~ 1,
  is.na(rare_non) ~ 99
))

# --------------------------------
# 14. SAMPLE SELECTION
# --------------------------------

cat("Step 12: Sample selection...\n")

# Keep only carriers and non-carriers
ukb_finally2 <- ukb_finally1 %>% filter(rare_non == 1 | rare_non == 0)

# Keep only White ethnicity
ukb_finally3 <- ukb_finally2 %>% filter(race == "White")

# Remove CNV positive individuals
ukb_finally4 <- ukb_finally3 %>% filter(cnv == "N")

# --------------------------------
# 15. EXCLUDE PRE-BASELINE CANCERS
# --------------------------------

cat("Step 13: Excluding pre-baseline cancers...\n")

# Exclude cancers diagnosed before baseline visit
cancer_types <- c("Breast", "Prostate", "Colorectal", "Stomach", 
                  "Cervical", "Uterine", "Testicular", "Ovarian", "Bladder",
                  "Kidney", "Pancreas", "Thyroid", "Laryngeal", "Tongue",
                  "oesophageal", "Brain", "Liver", "Lung")

ukb_finally4$in.exclu <- 'inclu'
for (col in cancer_types) {
  col_name <- paste0(col, '.y')
  if (col_name %in% colnames(ukb_finally4)) {
    ukb_finally4$in.exclu <- ifelse(
      !is.na(ukb_finally4[[col_name]]) & 
        as.Date(ukb_finally4[[col_name]]) <= as.Date(ukb_finally4$date_attending),
      'exclu',
      ukb_finally4$in.exclu
    )
  }
}

data_ukb <- subset(ukb_finally4, in.exclu %in% 'inclu')

# --------------------------------
# 16. CALCULATE TIME-TO-EVENT
# --------------------------------

cat("Step 14: Calculating time-to-event...\n")

# Calculate time to cancer diagnosis or censoring
for (cancer in cancer_types) {
  x_col <- paste0(cancer, ".x")
  y_col <- paste0(cancer, ".y")
  time_col <- paste0(cancer, ".time")
  
  data_ukb[[time_col]] <- ifelse(data_ukb[[x_col]] %in% 'Yes', 
                                 difftime(as.Date(data_ukb[[y_col]]), 
                                         as.Date(data_ukb$date_attending), 
                                         units = "days"),
                                 difftime(censoring_date, 
                                         as.Date(data_ukb$date_attending), 
                                         units = "days"))
}

# --------------------------------
# 17. MULTIPLE IMPUTATION
# --------------------------------

cat("Step 15: Multiple imputation for missing covariates...\n")

variables_to_impute <- c("SBP", "DBP", "alcohol", "smoke",       
                         "BMI", "income", "qualifications", "area")

# Convert character variables to factors
data_ukb$alcohol <- as.factor(data_ukb$alcohol)
data_ukb$smoke <- as.factor(data_ukb$smoke)
data_ukb$income <- as.factor(data_ukb$income)
data_ukb$qualifications <- as.factor(data_ukb$qualifications)
data_ukb$area <- as.factor(data_ukb$area)

# Define imputation methods
imp_methods <- make.method(data_ukb)
imp_methods[c("SBP", "DBP", "BMI")] <- "pmm"
imp_methods[c("alcohol", "smoke", "income", "qualifications", "area")] <- "polyreg"

# Build predictor matrix
pred <- make.predictorMatrix(data_ukb)
pred[,] <- 0

pred["SBP", c("DBP", "BMI", "age", "sex")] <- 1
pred["DBP", c("SBP", "BMI", "age", "sex")] <- 1
pred["BMI", c("SBP", "DBP", "age", "sex")] <- 1

pred["alcohol", c("smoke", "income", "age", "sex")] <- 1
pred["smoke", c("alcohol", "income", "age", "sex")] <- 1

pred["income", c("qualifications", "area", "age", "sex")] <- 1
pred["qualifications", c("income", "area", "age", "sex")] <- 1
pred["area", c("income", "qualifications", "age", "sex")] <- 1

# Run MICE
imp <- mice(data_ukb,
            method = imp_methods,
            predictorMatrix = pred,
            m = 5,
            maxit = 10,
            seed = 2025)

data_imputed <- complete(imp, 1)
data_ukb1 <- data_imputed

# --------------------------------
# 18. COVARIATE FACTOR CONVERSION
# --------------------------------

cat("Step 16: Converting covariates to factors...\n")

# Factor conversion for Cox regression
data_ukb1$smoke <- factor(data_ukb1$smoke)
data_ukb1$alcohol <- factor(data_ukb1$alcohol)
data_ukb1$area <- factor(data_ukb1$area)
data_ukb1$income <- factor(data_ukb1$income)
data_ukb1$qualifications <- factor(data_ukb1$qualifications)
data_ukb1$education <- factor(data_ukb1$education)
data_ukb1$hypertension <- factor(data_ukb1$hypertension, levels = c('No', 'Yes'))
data_ukb1$CHD <- factor(data_ukb1$CHD, levels = c('No', 'Yes'))
data_ukb1$DM <- factor(data_ukb1$DM, levels = c('No', 'Yes'))
data_ukb1$stroke <- factor(data_ukb1$stroke, levels = c('No', 'Yes'))
data_ukb1$rare_non <- factor(data_ukb1$rare_non)

# Convert cancer status to factors
for (cancer in cancer_types) {
  data_ukb1[[paste0(cancer, ".x")]] <- factor(data_ukb1[[paste0(cancer, ".x")]])
}

txt_ukb <- data_ukb1

# --------------------------------
# 19. COX REGRESSION - MODEL 2 (FULL ADJUSTMENT)
# --------------------------------

cat("Step 17: Running Cox regression (Model 2 - Full adjustment)...\n")

models2 <- list()
for (cancer in cancer_types) {
  x_col <- paste0(cancer, ".x")
  event_col <- cancer
  time_col <- paste0(cancer, ".time")
  
  txt_ukb[[event_col]] <- ifelse(txt_ukb[[x_col]] == 'Yes', 1, 0)
  
  formula_str <- paste0("Surv(", time_col, ", ", event_col, ") ~ rare_non + age + sex + BMI + income + SBP + DBP + 
                        smoke + alcohol + area + hypertension + CHD + DM + stroke + qualifications + DD +
                        genetic_01 + genetic_02 + genetic_03 + genetic_04 + genetic_05 + 
                        genetic_06 + genetic_07 + genetic_08 + genetic_09 + genetic_10")
  
  formula <- as.formula(formula_str)
  models2[[cancer]] <- coxph(formula, data = txt_ukb)
}

# Export results
output_file <- "results/cancer_model2_cox_rarenon.csv"
sink(output_file)
for(cancer in cancer_types){
  cat("\nSummary for", cancer, "cancer:\n")
  summary_output <- capture.output(summary(models2[[cancer]]))
  writeLines(summary_output)
}
sink()
cat("Model 2 results saved to:", output_file, "\n")

# --------------------------------
# 20. COX REGRESSION - MODEL 1 (MINIMAL ADJUSTMENT)
# --------------------------------

cat("Step 18: Running Cox regression (Model 1 - Minimal adjustment)...\n")

models1 <- list()
for (cancer in cancer_types) {
  x_col <- paste0(cancer, ".x")
  event_col <- cancer
  time_col <- paste0(cancer, ".time")
  
  txt_ukb[[event_col]] <- ifelse(txt_ukb[[x_col]] == 'Yes', 1, 0)
  
  formula_str <- paste0("Surv(", time_col, ", ", event_col, ") ~ rare_non + age + sex")
  formula <- as.formula(formula_str)
  models1[[cancer]] <- coxph(formula, data = txt_ukb)
}

output_file <- "results/cancer_model1_cox_rarenon.csv"
sink(output_file)
for(cancer in cancer_types){
  cat("\nSummary for", cancer, "cancer:\n")
  summary_output <- capture.output(summary(models1[[cancer]]))
  writeLines(summary_output)
}
sink()
cat("Model 1 results saved to:", output_file, "\n")

# --------------------------------
# 21. BASELINE CHARACTERISTICS TABLE
# --------------------------------

cat("Step 19: Generating baseline characteristic tables...\n")

PTEN_phenotype_txt <- txt_ukb

# Baseline characteristics
cont_var <- c("genetic_01", "genetic_02", "genetic_03", "genetic_04",
              "genetic_05", "genetic_06", "genetic_07", "genetic_08",
              "genetic_09", "genetic_10", "age", "SBP", "DBP", "BMI")
cat_var <- c("sex", "area", "income", "smoke", "alcohol", "qualifications",
             "DM", "CHD", "hypertension", "stroke", "DD")
group_var <- "rare_non"

cancer_PTEN_chi <- gttable1(data = PTEN_phenotype_txt, group_var = group_var,
                            cat_var = cat_var, cont_var = cont_var,
                            out_format = 'dataframe', pDigits = 3, method = 'auto')

cancer_PTEN_chi <- cancer_PTEN_chi[, !colnames(cancer_PTEN_chi) %in% c("variable", "test_name", "var_type", 
                                                                       "var_label", "row_type", "stat_label", 
                                                                       "estimate", "parameter")]
cancer_PTEN_chi[cancer_PTEN_chi == "<br />"] <- NA

write.table(cancer_PTEN_chi, file = "results/PTEN_chi_baseline_rare_non.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("Baseline table saved to: results/PTEN_chi_baseline_rare_non.txt\n")

# --------------------------------
# 22. CANCER INCIDENCE TABLE
# --------------------------------

cat("Step 20: Generating cancer incidence table...\n")

cont_var <- c("age")
cat_var <- cancer_types

cancer_PTEN_chi2 <- gttable1(data = PTEN_phenotype_txt, group_var = group_var,
                             cat_var = cat_var, cont_var = cont_var,
                             out_format = 'dataframe', pDigits = 3, method = 'auto')

cancer_PTEN_chi2 <- cancer_PTEN_chi2[, !colnames(cancer_PTEN_chi2) %in% c("variable", "test_name", "var_type", 
                                                                          "var_label", "row_type", "stat_label", 
                                                                          "estimate", "parameter")]
cancer_PTEN_chi2[cancer_PTEN_chi2 == "<br />"] <- NA

write.table(cancer_PTEN_chi2, file = "results/PTEN_chi_cancer_rare_non.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("Cancer incidence table saved to: results/PTEN_chi_cancer_rare_non.txt\n")

# --------------------------------
# 23. SENSITIVITY ANALYSIS - EXCLUDE DEATHS WITHIN 2 YEARS
# --------------------------------

cat("Step 21: Sensitivity analysis - excluding deaths within 2 years...\n")

# Exclude deaths within 2 years after baseline
data_ukb1$death_time <- difftime(as.Date(data_ukb1$death_date), 
                                 as.Date(data_ukb1$date_attending), 
                                 units = "days")
data_ukb1$in.exclu.death <- ifelse(data_ukb1$death_time < death_exclusion_days, "exclu", "inclu")
data_ukb1$in.exclu.death[is.na(data_ukb1$in.exclu.death)] <- "inclu"

data_ukb3 <- subset(data_ukb1, in.exclu.death %in% 'inclu')

# Cox regression on subset
models2_death <- list()
for (cancer in cancer_types) {
  x_col <- paste0(cancer, ".x")
  event_col <- cancer
  time_col <- paste0(cancer, ".time")
  
  data_ukb3[[event_col]] <- ifelse(data_ukb3[[x_col]] == 'Yes', 1, 0)
  
  formula_str <- paste0("Surv(", time_col, ", ", event_col, ") ~ rare_non + age + sex + BMI + income + SBP + DBP + 
                        smoke + alcohol + area + hypertension + CHD + DM + stroke + DD + qualifications +
                        genetic_01 + genetic_02 + genetic_03 + genetic_04 + genetic_05 + 
                        genetic_06 + genetic_07 + genetic_08 + genetic_09 + genetic_10")
  
  formula <- as.formula(formula_str)
  models2_death[[cancer]] <- coxph(formula, data = data_ukb3)
}

output_file <- "results/cancer_model2_cox_exclude_deaths_2y.csv"
sink(output_file)
for(cancer in cancer_types){
  cat("\nSummary for", cancer, "cancer:\n")
  summary_output <- capture.output(summary(models2_death[[cancer]]))
  writeLines(summary_output)
}
sink()
cat("Sensitivity analysis (exclude deaths) saved to:", output_file, "\n")

# --------------------------------
# 24. SENSITIVITY ANALYSIS - COMPLETE CASES
# --------------------------------

cat("Step 22: Sensitivity analysis - complete cases...\n")

# Delete rows with missing covariates
specified_columns <- c('genetic_01','genetic_02','genetic_03','genetic_04','genetic_05',
                       'genetic_06','genetic_07','genetic_08','genetic_09','genetic_10',
                       'age','sex', 'SBP','DBP', 'BMI' ,'smoke','alcohol','area','income',
                       'DM','hypertension','CHD','stroke','DD','qualifications')

data_ukb_na <- data_ukb[complete.cases(data_ukb[, specified_columns]), ]

# Cox regression on complete cases
models2_na <- list()
for (cancer in cancer_types) {
  x_col <- paste0(cancer, ".x")
  event_col <- cancer
  time_col <- paste0(cancer, ".time")
  
  data_ukb_na[[event_col]] <- ifelse(data_ukb_na[[x_col]] == 'Yes', 1, 0)
  
  formula_str <- paste0("Surv(", time_col, ", ", event_col, ") ~ rare_non + age + sex + BMI + income + SBP + DBP + 
                        smoke + alcohol + area + hypertension + CHD + DM + stroke + DD + qualifications +
                        genetic_01 + genetic_02 + genetic_03 + genetic_04 + genetic_05 + 
                        genetic_06 + genetic_07 + genetic_08 + genetic_09 + genetic_10")
  
  formula <- as.formula(formula_str)
  models2_na[[cancer]] <- coxph(formula, data = data_ukb_na)
}

output_file <- "results/cancer_model2_cox_complete_cases.csv"
sink(output_file)
for(cancer in cancer_types){
  cat("\nSummary for", cancer, "cancer:\n")
  summary_output <- capture.output(summary(models2_na[[cancer]]))
  writeLines(summary_output)
}
sink()
cat("Sensitivity analysis (complete cases) saved to:", output_file, "\n")

# --------------------------------
# 25. STRATIFIED ANALYSIS
# --------------------------------

cat("Step 23: Running stratified analyses...\n")

# Create stratification variables
txt_ukb <- txt_ukb %>%
  mutate(
    age_group = ifelse(age < 60, "1", "2"),
    Obesity = ifelse(BMI >= 30, "Yes", "No"),
    smoke_status = ifelse(smoke == "Current", "Yes", "No"),
    alcohol_status = ifelse(alcohol == "Current", "Yes", "No")
  )

txt_ukb$sex <- as.factor(txt_ukb$sex)
txt_ukb$DD <- as.factor(txt_ukb$DD)

# Example stratified analysis for uterine cancer
cat("Example stratified analysis for Uterine cancer:\n")

# Stratification variables
strat_vars <- c("sex", "age_group", "Obesity", "smoke_status", "alcohol_status", 
                "DM", "hypertension", "CHD", "stroke", "DD")

# Run stratified analyses
for (strat_var in strat_vars) {
  cat("\nStratified by", strat_var, ":\n")
  
  results <- txt_ukb %>%
    group_by(!!sym(strat_var)) %>%
    do(
      model = coxph(Surv(Uterine.time, Uterine) ~ rare_non + age + sex + BMI + income + SBP + DBP +
                      smoke + alcohol + area + hypertension + CHD + DM + stroke + DD + qualifications +
                      genetic_01 + genetic_02 + genetic_03 + genetic_04 + genetic_05 +
                      genetic_06 + genetic_07 + genetic_08 + genetic_09 + genetic_10, 
                    data = .)
    )
  
  # Print results
  for (i in 1:nrow(results)) {
    cat("  Stratum:", results[[strat_var]][i], "\n")
    print(summary(results$model[[i]]))
  }
}

# --------------------------------
# 26. INTERACTION ANALYSIS
# --------------------------------

cat("Step 24: Running interaction analyses...\n")

# Test interactions with PTEN variant status
cat("Testing interactions for Uterine cancer:\n")

interaction_terms <- c("sex", "age_group", "Obesity", "smoke_status", "alcohol_status",
                       "DM", "hypertension", "CHD", "stroke", "DD")

for (interact_var in interaction_terms) {
  formula_str <- paste0("Surv(Uterine.time, Uterine) ~ rare_non * ", interact_var, 
                       " + age + sex + BMI + income + SBP + DBP + smoke + alcohol + area + hypertension + CHD + DM + stroke + DD + qualifications + genetic_01 + genetic_02 + genetic_03 + genetic_04 + genetic_05 + genetic_06 + genetic_07 + genetic_08 + genetic_09 + genetic_10")
  
  model <- coxph(as.formula(formula_str), data = txt_ukb)
  
  cat("\nInteraction with", interact_var, ":\n")
  print(summary(model))
}

# --------------------------------
# 27. SESSION INFO AND CLEANUP
# --------------------------------

cat("Step 25: Saving session information...\n")

# Save session info
sink("results/session_info.txt")
print(sessionInfo())
sink()

# Summary of output files
cat("\n" . strrep("=", 60) . "\n")
cat("ANALYSIS COMPLETE!\n")
cat(strrep("=", 60) . "\n\n")

cat("Output files generated:\n")
cat("1. results/cancer_model1_cox_rarenon.csv - Minimal adjustment model\n")
cat("2. results/cancer_model2_cox_rarenon.csv - Full adjustment model\n")
cat("3. results/PTEN_chi_baseline_rare_non.txt - Baseline characteristics\n")
cat("4. results/PTEN_chi_cancer_rare_non.txt - Cancer incidence\n")
cat("5. results/cancer_model2_cox_exclude_deaths_2y.csv - Sensitivity (no deaths <2y)\n")
cat("6. results/cancer_model2_cox_complete_cases.csv - Sensitivity (complete cases)\n")
cat("7. results/session_info.txt - R session information\n")

cat("\nTotal participants in final analysis:", nrow(txt_ukb), "\n")
cat("PTEN rare variant carriers:", sum(txt_ukb$rare_non == 1), "\n")
cat("Controls (non-carriers):", sum(txt_ukb$rare_non == 0), "\n")

cat("\nAnalysis completed at:", Sys.time(), "\n")
cat(strrep("=", 60) . "\n")

# End of script