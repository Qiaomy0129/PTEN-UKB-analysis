# Cox regression models

library(survival)

imp <- readr::read_rds('data/processed/imputation_mids.rds')

formula_str <- 'Surv(Breast.time, Breast) ~ rare_non + age + sex + BMI + income + SBP + DBP + smoke + alcohol + area + hypertension + CHD + DM + stroke + qualifications'

fit <- with(imp, coxph(as.formula(formula_str), data=.))
pooled <- pool(fit)

readr::write_csv(as.data.frame(summary(pooled)), 'results/pooled_cox_Breast.csv')

message('05_cox_models done.')
