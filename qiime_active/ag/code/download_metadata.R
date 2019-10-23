#!/usr/bin/env Rscript

library(readr)

DATA_DIR <- '~/AG'
DEST <- file.path(DATA_DIR,'metadata')

MAP <- read_delim(file.path(DATA_DIR,'PRJEB11419.txt'),delim='\t')

ERS <- unique(MAP$secondary_sample_accession)
ERS <- ERS[!(ERS %in% gsub('\\.xml','',list.files(DEST)))]

url <- 'http://www.ebi.ac.uk/ena/data/view/%s&display=xml'

for (ers in ERS){
  download.file(sprintf(url,ers),file.path(DEST,paste0(ers,'.xml')))
}
