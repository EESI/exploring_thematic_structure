#!/usr/bin/env Rscript

library(readr)
library(dplyr)

DATA_DIR <- '~/AG'
SAMPS_FECAL_PATH <- file.path(DATA_DIR,'samps_fecal.txt')

ERRS <- gsub('\\_filt\\.fastq\\.gz','',list.files(file.path(DATA_DIR,'filtered')))
PROJDATA <- read_delim(file.path(DATA_DIR,'PRJEB11419.txt'),'\t')
METADATA <- readRDS(file.path(DATA_DIR,'metadata.rds'))

samps_fecal <- PROJDATA %>%
  dplyr::select(PRIMARY_ID=secondary_sample_accession,run_accession) %>%
  filter(run_accession %in% ERRS) %>%
  left_join(METADATA %>% dplyr::select(PRIMARY_ID,body_site),by='PRIMARY_ID') %>%
  filter(body_site %in% 'UBERON:feces') %>%
  dplyr::select(run_accession) 

write_lines(samps_fecal$run_accession,SAMPS_FECAL_PATH)
