#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 1) stop('Please specify number of topics, K.', call.=FALSE)

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


K <- as.integer(args[1])

data_dir <- '~/Dropbox/stm_microbiome/qiime_active/ag'
models_dir <- file.path(data_dir,'models',K)
dir.create(models_dir,showWarnings=FALSE,recursive=TRUE)

DAT <- readRDS(file.path(data_dir,'train_test.rds'))
DOCS <- DAT$docs_train
VOCAB <- DAT$VOCAB
META <- DAT$meta_train

seed_fit <- 12

fit0 <- stm(documents=DOCS,vocab=VOCAB,K=K,
            data=META,
            max.em.its=500,init.type='Spectral',
            seed=seed_fit,
            verbose=TRUE,reportevery=25)

fit_fn <- sprintf('stm_train_k_%s.rds',K)
saveRDS(list(seed=seed_fit,
             K=K,
             OTU=DAT$OTU,
             VOCAB=DAT$VOCAB,
             TAX=DAT$TAX,
             fit=fit0,
             p=DAT$p,
             idx_train=DAT$idx_train,
             meta_train=DAT$meta_train,
             docs_train=DAT$docs_train,
             meta_test=DAT$meta_test,
             docs_test=DAT$docs_test),
        file.path(models_dir,fit_fn))