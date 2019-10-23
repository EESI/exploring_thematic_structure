#!/usr/bin/env Rscript

library(dada2)

data_dir <- '~/Qiime/Children/Processing/01-raw/PRJNA237362'
dada_dir <- '~/Dropbox/stm_microbiome/dada_active/gevers'

seqtab <- readRDS(file.path(dada_dir,'seqtab.rds'))
ref <- '~/Dada2_ref/gg_13_8_train_set_97.fa.gz'

taxa <- assignTaxonomy(seqtab,ref,outputBootstraps=TRUE,verbose=TRUE)

saveRDS(taxa,file.path(dada_dir,'tax.rds'))

