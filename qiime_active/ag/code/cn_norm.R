#!/usr/bin/env Rscript

library(stm)
library(biom)
library(readr)
library(tidyr)
library(dplyr)
library(fastICA)
library(randomForest)
library(stringr)
library(kernlab)
library(Rcpp)
library(parallel)
library(foreach)
library(ape)
library(phyloseq)
library(doParallel)
library(stm)
library(LDAvis)
library(caret)
library(glmnet)
library(ggplot2)
library(knitr)
library(gridExtra)
library(DESeq2)
library(tibble)
library(phyloseq)

source('~/Dropbox/stm_microbiome/code_active/stm_functions.R')
source('~/Dropbox/stm_microbiome/code_active/nav_froz_fxns_3.R') 
source('~/Dropbox/stm_microbiome/code_active/performance_1.R')
source('~/Dropbox/stm_microbiome/code_active/framework.R')

data_dir <- '~/Dropbox/stm_microbiome/qiime_active/ag/'

cat('Reading biom file.\n')
BIOM <- read_biom(file.path(data_dir,'picked_otus/otu_table.biom'))

cat('Extracting otu table.\n')
OTU <- as.matrix(biom_data(BIOM))
cat('Extracting tax table.\n')
TAX <- observation_metadata(BIOM)
colnames(TAX) <- c('Kingdom','Phylum','Class','Order','Family','Genus','Species')

cat('Saving otu table.\n')
saveRDS(as.data.frame(t(OTU)),file.path(data_dir,'seqtab.rds'))
saveRDS(TAX,file.path(data_dir,'seqtax.rds'))

cat('Write temp otu table biom file.\n')
BIOM_TMP <- tempfile()
write_biom(make_biom(OTU,NULL,TAX),BIOM_TMP)

cat('Normalizing by copy number.\n')
system2('normalize_by_copy_number.py',
        args=c('-i',BIOM_TMP,
               '-o',file.path(data_dir,'picked_otus/otu_table_cnnorm.biom')))

cat('Reading normalized biom file.\n')
BIOM <- read_biom(file.path(data_dir,'picked_otus/otu_table_cnnorm.biom'))
cat('Extracting normalized otu table.\n')
OTU <- round(as.matrix(biom_data(BIOM)))
cat('Extracting normalized tax table.\n')
TAX <- observation_metadata(BIOM)[,c('Kingdom','Phylum','Class','Order','Family','Genus','Species')]

cat('Saving normalized otu table.\n')
saveRDS(as.data.frame(t(OTU)),file.path(data_dir,'seqtab_cnnorm.rds'))
cat('Saving normaliezd taxa table.\n')
saveRDS(as.data.frame(TAX),file.path(data_dir,'tax.rds'))
