#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 2) stop("Two arguments must be supplied.", call.=FALSE)

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

models_dir <- '~/Dropbox/stm_microbiome/qiime_active/ag/models'

fns <- list.files(models_dir,recursive=TRUE,full.names=TRUE)
fns <- fns[grep('/stm_k_',fns)]

seed_perm <- 89

fn <- fns[as.integer(args[1])]

dat <- readRDS(fn)

fit <- dat$fits[[as.integer(args[2])]]

K <- fit$settings$dim$K
DOCS <- dat$docs
VOCAB <- dat$vocab
META <- dat$meta

perm <- permTest(formula=~diet, fit, 'diet',
                 nruns=25,
                 documents=DOCS, vocab=VOCAB, data=META, seed=seed_perm,
                 uncertainty='Global',
                 stmverbose=TRUE)

saveRDS(list(dat=dat,
             K=fit$settings$dim$K,
             perm=perm,
             seed=seed_perm),
        file.path(models_dir,K,sprintf('perm_k_%s_%s.rds',K,names(dat$fits)[as.integer(args[2])])))

