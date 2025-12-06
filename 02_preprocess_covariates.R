# Preprocess covariates from UKB main table
ukb_main <- readr::read_rds('data/processed/ukb_main.rds')

ukb_finally <- ukb_main %>%
  dplyr::select(
    eid,
    sex = sex_f31_0_0,
    birth_year = year_of_birth_f34_0_0,
    date_attending = date_of_attending_assessment_centre_f53_0_0,
    age = age_at_recruitment_f21022_0_0,
    BMI = body_mass_index_bmi_f21001_0_0,
    SBP_1 = systolic_blood_pressure_automated_reading_f4080_0_0,
    SBP_2 = systolic_blood_pressure_automated_reading_f4080_0_1,
    DBP_1 = diastolic_blood_pressure_automated_reading_f4079_0_0,
    DBP_2 = diastolic_blood_pressure_automated_reading_f4079_0_1,
    smoke = smoking_status_f20116_0_0,
    alcohol = alcohol_drinker_status_f20117_0_0,
    income = average_total_household_income_before_tax_f738_0_0,
    qualifications = qualifications_f6138_0_0,
    area = home_area_population_density_urban_or_rural_f20118_0_0
  )

ukb_finally <- ukb_finally %>%
  mutate(
    SBP = rowMeans(select(., SBP_1, SBP_2), na.rm = TRUE),
    DBP = rowMeans(select(., DBP_1, DBP_2), na.rm = TRUE)
  )

readr::write_rds(ukb_finally, 'data/processed/ukb_finally.rds')
message('02_preprocess_covariates done.')
