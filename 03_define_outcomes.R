# Define outcome variables and compute follow-up times

ukb_finally <- readr::read_rds('data/processed/ukb_finally.rds')

if (file.exists('data/processed/date_disease.rds')) {
  date_disease <- readr::read_rds('data/processed/date_disease.rds')
} else {
  message('date_disease not found, using placeholder.')
  date_disease <- tibble(eid=integer(), Breast=as.Date(character()))
}

ukb_finally <- left_join(ukb_finally, date_disease, by='eid')

cancer_types <- c('Breast','Prostate','Colorectal','Stomach','Cervical','Uterine','Testicular','Ovarian',
    'Bladder','Kidney','Pancreas','Thyroid','Laryngeal','Tongue','oesophageal','Brain','Liver','Lung')

ukb_finally$in.exclu <- 'inclu'
for (col in cancer_types) {
  y_col <- paste0(col, '.y')
  if (y_col %in% colnames(ukb_finally)) {
    ukb_finally$in.exclu <- ifelse(
      !is.na(ukb_finally[[y_col]]) &
      as.Date(ukb_finally[[y_col]]) <= as.Date(ukb_finally$date_attending),
      'exclu', ukb_finally$in.exclu)
  }
}

data_ukb <- subset(ukb_finally, in.exclu == 'inclu')

for (c in cancer_types) {
  x_col <- paste0(c, '.x')
  y_col <- paste0(c, '.y')
  time_col <- paste0(c, '.time')
  data_ukb[[time_col]] <- ifelse(
    !is.na(data_ukb[[x_col]]) & data_ukb[[x_col]] == 'Yes',
    as.numeric(difftime(as.Date(data_ukb[[y_col]]), as.Date(data_ukb$date_attending), units='days')),
    as.numeric(difftime(as.Date('2021-04-25'), as.Date(data_ukb$date_attending), units='days'))
  )
}

readr::write_rds(data_ukb, 'data/processed/data_ukb_preimpute.rds')
message('03_define_outcomes done.')
