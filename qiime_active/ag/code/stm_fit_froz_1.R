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
fit_fn <- sprintf('stm_train_k_%s.rds',K)

DAT <- readRDS(file.path(models_dir,fit_fn))
DOCS <- DAT$docs_train
VOCAB <- DAT$VOCAB
META <- DAT$meta_train

fit_train <- DAT$fit
docs_train <- DAT$docs_train
docs_test <- DAT$docs_test
meta_train <- DAT$meta_train
meta_test <- DAT$meta_test

seed_fit <- DAT$seed
seed_post <- 354

froz_train <- ctm_frozen2(fit_train,docs_train,VOCAB,
                               seed=seed_fit,max.em.its=500,emtol=1e-5,avg_iters=1,
                               verbose=TRUE,
                               data=meta_train,covariate=NULL,
                               parallel=FALSE,nc=20)


froz_test <- ctm_frozen2(fit_train,docs_test,VOCAB,
                              seed=seed_fit,max.em.its=500,emtol=1e-5,avg_iters=1,
                              verbose=TRUE,
                              data=meta_test,covariate=NULL,
                              parallel=FALSE,nc=20)


set.seed(seed_post)
cat('Sampling from training posterior.\n')
theta_train <- thetaPosterior(froz_train,nsims=500,type='Global')
theta_train <- t(sapply(theta_train,colMeans))

cat('Sampling from testing posterior.\n')
theta_test <- thetaPosterior(froz_test,nsims=500,type='Global')
theta_test <- t(sapply(theta_test,colMeans))

cat('Saving results.\n')
froz_fn <- sprintf('stm_froz_k_%s.rds',K)
saveRDS(list(seeds=list(seed_fit=seed_fit,seed_post=seed_post),
             K=K,
             OTU=DAT$OTU,
             VOCAB=DAT$VOCAB,
             TAX=DAT$TAX,
             froz_train=froz_train,
             froz_test=froz_test,
             theta_train=theta_train,
             theta_test=theta_test,
             p=DAT$p,
             idx_train=DAT$idx_train,
             meta_train=DAT$meta_train,
             docs_train=DAT$docs_train,
             meta_test=DAT$meta_test,
             docs_test=DAT$docs_test),
        file.path(models_dir,froz_fn))
