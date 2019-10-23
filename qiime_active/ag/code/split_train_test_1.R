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


source('~/Dropbox/stm_microbiome/code_active/stm_functions.R')
source('~/Dropbox/stm_microbiome/code_active/nav_froz_fxns_3.R') 
source('~/Dropbox/stm_microbiome/code_active/performance_1.R')
source('~/Dropbox/stm_microbiome/code_active/framework.R')
source('~/Dropbox/stm_microbiome/code_active/stm_ESTEP.R')


data_dir <- '~/Dropbox/stm_microbiome/qiime_active/ag'


seed_split <- 121


PS <- readRDS(file.path(data_dir,'data_clean.rds'))

OTU <- as(otu_table(PS),'matrix')
TAX <- as(tax_table(PS),'matrix')
META <- as(sample_data(PS),'data.frame')

taxon_ids <- data.frame(long=colnames(OTU),
                        ids=paste0('otu',1:NCOL(OTU)),
                        row.names=paste0('otu',1:NCOL(OTU)))

counts <- list(META_clean=META,
               table_clean=OTU,
               ids=taxon_ids)
colnames(counts$table_clean) <- as.character(counts$ids$ids)

if (!identical(rownames(counts$table_clean),counts$META_clean$SampleID)) stop('Please make table sample names match METAdata sample IDs!\n')

VOCAB <- as.character(counts$ids$ids)
DOCS <- lapply(1:NROW(counts$table_clean), function(i) reshape_doc(counts$table_clean[i,],VOCAB)) 
META <- counts$META_clean


cat('Total docs:',length(DOCS),'\n',
    'Total vocab:',length(VOCAB),'\n',
    'Total meta:',nrow(META),'x',ncol(META),'\n',
    'Total diet:',names(table(META$diet)),'=',c(table(META$diet)),'\n')


META$diet <- ifelse(META$diet == 'V',1,0)


set.seed(seed_split)
p <- .8
idx_train <- unlist(caret::createDataPartition(y=META$diet,times=1,p=p,list=TRUE))
meta_train <- META[idx_train,]
docs_train <- DOCS[idx_train]
meta_test <- META[-idx_train,]
docs_test <- DOCS[-idx_train]


saveRDS(list(seed=seed_split,
              META=META,
              OTU=OTU,
              VOCAB=VOCAB,
              TAX=TAX,
              p=p,
              idx_train=idx_train,
              meta_train=meta_train,
              docs_train=docs_train,
              meta_test=meta_test,
              docs_test=docs_test),
         file.path(data_dir,'train_test.rds'))