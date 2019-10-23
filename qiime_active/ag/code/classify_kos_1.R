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
library(doSNOW)
library(snow)
library(phyloseq)
library(stm)
library(LDAvis)
library(caret)
library(glmnet)
library(ggplot2)
library(knitr)
library(gridExtra)
library(caret)
library(pROC)


source('~/Dropbox/stm_microbiome/code_active/stm_functions.R')
source('~/Dropbox/stm_microbiome/code_active/nav_froz_fxns_3.R') 
source('~/Dropbox/stm_microbiome/code_active/performance_1.R')
source('~/Dropbox/stm_microbiome/code_active/framework.R')
source('~/Dropbox/stm_microbiome/code_active/stm_ESTEP.R')

data_dir <- '~/Dropbox/stm_microbiome/qiime_active/ag'

biom_file <- read_biom(file.path(data_dir,'picked_otus','ko_table.biom'))
ko <- as.matrix(biom_data(biom_file))
ko <- ko[rowSums(ko)>0,]
ko <- t(ko)

DAT <- readRDS(file.path(data_dir,'train_test.rds'))

META <- DAT$META

ko <- ko[rownames(META),]
idx_train <- DAT$idx_train
meta_train <- META[idx_train,]
ko_train <- ko[idx_train,]
meta_test <- META[-idx_train,]
ko_test <- ko[-idx_train,]

ko_train <- apply(ko_train,2,function(x) (x-mean(x))/(max(x)-min(x)))
ko_test <- apply(ko_test,2,function(x) (x-mean(x))/(max(x)-min(x)))
filt <- unique(c(which(is.na(ko_train),arr.ind=TRUE)[,2],which(is.na(ko_test),arr.ind=TRUE)[,2]))
ko_train <- ko_train[,-filt]
ko_test <- ko_test[,-filt]

train_labels <- factor(ifelse(meta_train$diet==1,'V','O'),levels=c('V','O'))
test_labels <- factor(ifelse(meta_test$diet==1,'V','O'),levels=c('V','O'))

ncores <- 50
cl <- makeCluster(ncores,type='SOCK')
registerDoSNOW(cl)

tr_ctrl <- trainControl(method='repeatedCV', number=10, repeats=10, 
                        classProbs=TRUE,
                        sampling='down',
                        summaryFunction=twoClassSummary,
                        allowParallel=TRUE)

param_sweep <- expand.grid(mtry=unique(floor(seq(ncol(ko_train)^.25,ncol(ko_train)^.75,length=10))))
out_rf <- train(x=ko_train,
                y=train_labels,
                method='rf',
                ntree=128,
                tuneGrid=param_sweep,
                importance=TRUE,
                proximity=TRUE,
                metric='ROC', 
                maximize=TRUE, 
                trControl=tr_ctrl)

stopCluster(cl)

pred_p <- predict(out_rf,ko_test,type='prob')
pred_Y <- factor(colnames(pred_p)[apply(pred_p,1,which.max)],levels=colnames(pred_p))
pred_out <- data.frame(obs=test_labels,
                       pred=pred_Y,
                       pred_p)
two_class_summary <- twoClassSummary(pred_out,lev=levels(pred_out$pred))
confusion_matrix <- caret::confusionMatrix(pred_Y,test_labels)

roc_curve <- roc(response=test_labels, 
                 predictor=pred_p[,1], 
                 levels=rev(unique(test_labels)))

saveRDS(list(OTU=DAT$OTU,
             ko=ko,
             VOCAB=DAT$VOCAB,
             TAX=DAT$TAX,
             ko_train=ko_train,
             ko_test=ko_test,
             train_labels=train_labels,
             test_labels=test_labels,
             tr_ctrl=tr_ctrl,
             param_sweep=param_sweep,
             out_rf=out_rf,
             pred_p=pred_p,
             pred_Y=pred_Y,
             pred_out=pred_out,
             two_class_summary=two_class_summary,
             confusion_matrix=confusion_matrix,
             roc_curve=roc_curve,
             p=DAT$p,
             idx_train=DAT$idx_train,
             meta_train=DAT$meta_train,
             meta_test=DAT$meta_test),
        file.path(data_dir,'classify_kos.rds'))
