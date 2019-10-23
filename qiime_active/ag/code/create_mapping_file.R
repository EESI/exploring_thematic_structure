#!/usr/bin/env Rscript

library(readr)
library(dplyr)

DATA_DIR <- '~/AG'

MAP <- tibble(`#SampleID`=list.files(file.path(DATA_DIR,'filtered_fasta'))) %>%
  mutate(BarcodeSequence='',
         LinkerPrimerSequence='',
         InputFileName=`#SampleID`,
         Description='') %>%
  mutate(`#SampleID`=gsub('^(ERR.*)_filt.fna','\\1',`#SampleID`))
                   
write_delim(MAP,file.path(DATA_DIR,'map.txt'),'\t')
