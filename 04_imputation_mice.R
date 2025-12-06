# Multiple Imputation using MICE

data_ukb <- readr::read_rds('data/processed/data_ukb_preimpute.rds')

vars <- c('SBP','DBP','alcohol','smoke','BMI','income','qualifications','area')
for (v in vars[c(3:8)]) {
  data_ukb[[v]] <- as.factor(data_ukb[[v]])
}

imp_methods <- make.method(data_ukb)
imp_methods[c('SBP','DBP','BMI')] <- 'pmm'
imp_methods[c('alcohol','smoke','income','qualifications','area')] <- 'polyreg'

pred <- make.predictorMatrix(data_ukb)
pred[,] <- 0

pred['SBP', c('DBP','BMI','age','sex')] <- 1
pred['DBP', c('SBP','BMI','age','sex')] <- 1

set.seed(2025)
imp <- mice(data_ukb, method=imp_methods, predictorMatrix=pred, m=5, maxit=10)

readr::write_rds(imp, 'data/processed/imputation_mids.rds')
readr::write_rds(complete(imp,1), 'data/processed/data_ukb_imputed_example.rds')
message('04_imputation_mice done.')
