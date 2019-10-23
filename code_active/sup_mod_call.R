#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) stop("At least one argument must be supplied (input file).n", call.=FALSE)


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


params <- expand.grid(K=c(50),
                      content=c(FALSE),
                      variable=c('DIAGNOSIS'))
params <- params %>% 
   filter(!(content == TRUE  & K == 200)) %>%
   arrange(K,variable)

param <- 1


source('~/Dropbox/stm_microbiome/code_active/stm_functions.R')
source('~/Dropbox/stm_microbiome/code_active/nav_froz_fxns_3.R') 
source('~/Dropbox/stm_microbiome/code_active/performance_1.R')
source('~/Dropbox/stm_microbiome/code_active/framework.R')




random_seed <- FALSE
K <- params[param,]$K
cn_normalize <- TRUE
content <- params[param,]$content
seq_sim <- 's97'
variable <- as.character(params[param,]$variable)
dupnum <- NULL
test <- 'dx'   
var2 <- NULL



prepare_framework_ti_rarefy(random_seed,K,cn_normalize,content,variable,seq_sim)

rm(train_docs,test_docs,train_meta,test_meta,counts,otu_taxa_xavier)


meta$PCDAI[is.na(meta$PCDAI) & meta$DIAGNOSIS == 'Not IBD'] <- 0
filt_idx <- !is.na(meta$PCDAI) 
docs <- docs[filt_idx]
meta <- meta[filt_idx,]
meta <- meta %>% mutate(PCDAI_RANK = ranker(PCDAI))



save_folder <- '~/Dropbox/stm_microbiome/data_active/Gevers_ti/supervised'
dir.create(save_folder,showWarnings=FALSE,recursive=TRUE)   




if (args[1] == 1){
   cat('Fitting content model:',variable,'\n')
   fit <- tryCatch(stm(content=as.formula(paste0('~',variable)),
                       documents=docs, vocab=vocab, K=K,
                       data=meta,
                       max.em.its=500, 
                       init.type='Spectral',
                       # ngroups=10,
                       seed=seed_fit,verbose=TRUE,reportevery=25),
                   error = function(e) e)
}

if (args[1] == 2){
   cat('Fitting prevalence model:',variable,'\n')
   fit <- tryCatch(stm(prevalence=as.formula(paste0('~',variable)),
                       documents=docs, vocab=vocab, K=K,
                       data=meta,
                       max.em.its=500, 
                       init.type='Spectral',
                       # ngroups=10,
                       seed=seed_fit,verbose=TRUE,reportevery=25),
                   error = function(e) e)
}

if (args[1] == 3){
   cat('Fitting prevalence / content model:',variable,'/', variable,'\n')
   fit <- tryCatch(stm(prevalence=as.formula(paste0('~',variable)),
                       content=as.formula(paste0('~',variable)),
                       documents=docs, vocab=vocab, K=K,
                       data=meta,
                       max.em.its=500, 
                       init.type='Spectral',
                       # ngroups=10,
                       seed=seed_fit,verbose=TRUE,reportevery=25),
                   error = function(e) e)
}

if (args[1] == 4){
   var2 <- 'PCDAI'
   cat('Fitting prevalence model:',var2,'\n')
   fit <- tryCatch(stm(prevalence=as.formula(paste0('~',var2)),
                       documents=docs, vocab=vocab, K=K,
                       data=meta,
                       max.em.its=500, 
                       init.type='Spectral',
                       # ngroups=10,
                       seed=seed_fit,verbose=TRUE,reportevery=25),
                   error = function(e) e)
}

if (args[1] == 5){
   var2 <- 'PCDAI_RANK'
   cat('Fitting prevalence model:',var2,'\n')
   fit <- tryCatch(stm(prevalence=as.formula(paste0('~',var2)),
                       documents=docs, vocab=vocab, K=K,
                       data=meta,
                       max.em.its=500, 
                       init.type='Spectral',
                       # ngroups=10,
                       seed=seed_fit,verbose=TRUE,reportevery=25),
                   error = function(e) e)
}

if (args[1] == 6){
   var2 <- 'PCDAI'
   cat('Fitting prevalence / content model:',var2,variable,'\n')
   fit <- tryCatch(stm(prevalence=as.formula(paste0('~',var2)),
                       content=as.formula(paste0('~',variable)),
                       documents=docs, vocab=vocab, K=K,
                       data=meta,
                       max.em.its=500, 
                       init.type='Spectral',
                       # ngroups=10,
                       seed=seed_fit,verbose=TRUE,reportevery=25),
                   error = function(e) e)
}

if (args[1] == 6){
   var2 <- 'PCDAI_RANK'
   cat('Fitting prevalence / content model:',var2,variable,'\n')
   fit <- tryCatch(stm(prevalence=as.formula(paste0('~',var2)),
                       content=as.formula(paste0('~',variable)),
                       documents=docs, vocab=vocab, K=K,
                       data=meta,
                       max.em.its=500, 
                       init.type='Spectral',
                       # ngroups=10,
                       seed=seed_fit,verbose=TRUE,reportevery=25),
                   error = function(e) e)
}

if (args[1] == 7){
   var2 <- 'PCDAI'
   cat('Fitting spline prevalence model:',var2,'\n')
   fit <- tryCatch(stm(prevalence=as.formula(paste0('~s(',var2,')')),
                       documents=docs, vocab=vocab, K=K,
                       data=meta,
                       max.em.its=500, 
                       init.type='Spectral',
                       # ngroups=10,
                       seed=seed_fit,verbose=TRUE,reportevery=25),
                   error = function(e) e)
}

if (args[1] == 8){
   var2 <- 'PCDAI_RANK'
   cat('Fitting spline prevalence model:',var2,'\n')
   fit <- tryCatch(stm(prevalence=as.formula(paste0('~s(',var2,')')),
                       documents=docs, vocab=vocab, K=K,
                       data=meta,
                       max.em.its=500, 
                       init.type='Spectral',
                       # ngroups=10,
                       seed=seed_fit,verbose=TRUE,reportevery=25),
                   error = function(e) e)
}

if (args[1] == 9){
   var2 <- 'PCDAI'
   cat('Fitting spline prevalence / content model:',var2,variable,'\n')
   fit <- tryCatch(stm(prevalence=as.formula(paste0('~s(',var2,')')),
                       content=as.formula(paste0('~',variable)),
                       documents=docs, vocab=vocab, K=K,
                       data=meta,
                       max.em.its=500, 
                       init.type='Spectral',
                       # ngroups=10,
                       seed=seed_fit,verbose=TRUE,reportevery=25),
                   error = function(e) e)
}

if (args[1] == 10){
   var2 <- 'PCDAI_RANK'
   cat('Fitting spline prevalence / content model:',var2,variable,'\n')
   fit <- tryCatch(stm(prevalence=as.formula(paste0('~s(',var2,')')),
                       content=as.formula(paste0('~',variable)),
                       documents=docs, vocab=vocab, K=K,
                       data=meta,
                       max.em.its=500, 
                       init.type='Spectral',
                       # ngroups=10,
                       seed=seed_fit,verbose=TRUE,reportevery=25),
                   error = function(e) e)
}

file_name <- paste0('gevers_ti_rarefied_supervised_model_type',args[1],'.rds')
saveRDS(fit,file.path(save_folder,file_name))
