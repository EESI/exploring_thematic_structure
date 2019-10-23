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


require(snow)
require(doSNOW)
ncores <- 48
cl <- makeCluster(ncores,type='SOCK')
registerDoSNOW(cl)

model <- 1
for (seed_split in c(1,12,53,74,105)){
  
  dat <- readRDS(sprintf('~/Dropbox/stm_microbiome/data_active/AG_female/stm_s97_rarefied_10000_supervised/features/seed_%s/fits_K_35_model_%s.rds',
                         seed_split,model))
    
  raw <- readRDS('~/Dropbox/stm_microbiome/data_active/AG_female/stm_s97_rarefied_10000_supervised/features/raw_kegg.rds')
  raw$kegg <- raw$kegg[rowSums(raw$kegg) != 0,]
  
  kegg_train <- t(raw$kegg[,dat$meta$meta_train$SampleID])
  meta_train <- raw$meta[dat$meta$meta_train$SampleID,]
  kegg_test <- t(raw$kegg[,dat$meta$meta_test$SampleID])
  meta_test <- raw$meta[dat$meta$meta_test$SampleID,]
  
  kegg_train <- apply(kegg_train,2,qnormalize)
  kegg_test <- apply(kegg_test,2,qnormalize)
  
  labels_train <- meta_train$SEX
  labels_test <- meta_test$SEX
  
  rf_ctrl <- trainControl(method='cv',classProbs=TRUE,
                          summaryFunction=twoClassSummary,
                          allowParallel=TRUE)
  out_rf <- train(x=kegg_train,y=labels_train,
                  method='rf',ntree=1200,tuneLength=5,metric='ROC',trControl=rf_ctrl,
                  importance=TRUE)
  
  pred_rf <- predict(out_rf,kegg_test,type='prob')[,1]
  cm_rf <- caret::confusionMatrix(ifelse(pred_rf > .5,'female','male'),labels_test)
  imp <- varImp(out_rf)
  
  save_list_rf <- list(out=out_rf,
                       pred=pred_rf,
                       cm=cm_rf,
                       imp=imp,
                       labels=list(train=labels_train,
                                   test=labels_test),
                       kegg=list(train=kegg_train,
                                 test=kegg_test),
                       meta=list(train=meta_train,
                                 test=meta_test))
  saveRDS(save_list_rf,sprintf('~/Dropbox/stm_microbiome/data_active/AG_female/stm_s97_rarefied_10000_supervised/features/seed_%s/rf_features.rds',
                               seed_split))
}

stopCluster(cl)