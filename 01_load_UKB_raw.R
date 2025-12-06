# Load UKB raw data and perform initial merging
# This script does not upload raw files; paths point only to local 'data/raw/'.

raw_path <- 'data/raw/'
ukb_main_path <- file.path(raw_path, 'ukb_main.tsv')
if (!file.exists(ukb_main_path)) stop('Place ukb_main.tsv in data/raw/')

ukb_main <- readr::read_tsv(ukb_main_path)

pten_samples_path <- file.path(raw_path, 'PTEN_UKB_sampleDetail.txt.gz')
if (file.exists(pten_samples_path)) {
  PTEN_samplelist <- read.delim(pten_samples_path)
} else {
  message('PTEN sample list not found; continuing without PTEN details.')
}

dir.create('data/processed', showWarnings = FALSE)
readr::write_rds(ukb_main, 'data/processed/ukb_main.rds')
if (exists('PTEN_samplelist')) readr::write_rds(PTEN_samplelist, 'data/processed/PTEN_samplelist.rds')
message('Step 01 completed.')
