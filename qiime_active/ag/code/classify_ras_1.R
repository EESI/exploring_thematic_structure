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

DAT <- readRDS(file.path(data_dir,'train_test.rds'))

OTU <- DAT$OTU
OTU <- OTU/rowSums(OTU)
META <- DAT$META

idx_train <- DAT$idx_train
meta_train <- META[idx_train,]
otu_train <- OTU[idx_train,]
meta_test <- META[-idx_train,]
otu_test <- OTU[-idx_train,]

otu_train <- apply(otu_train,2,function(x) (x-mean(x))/(max(x)-min(x)))
otu_test <- apply(otu_test,2,function(x) (x-mean(x))/(max(x)-min(x)))
filt <- unique(which(is.na(otu_test),arr.ind=TRUE)[,2])
otu_train <- otu_train[,-filt]
otu_test <- otu_test[,-filt]

train_labels <- factor(ifelse(meta_train$diet==1,'V','O'),levels=c('V','O'))
names(train_labels) <- meta_train$SampleID
test_labels <- factor(ifelse(meta_test$diet==1,'V','O'),levels=c('V','O'))
names(test_labels) <- meta_test$SampleID

ncores <- 60
cl <- makeCluster(ncores,type='SOCK')
registerDoSNOW(cl)

tr_ctrl <- trainControl(method='repeatedCV', number=10, repeats=10, 
                        classProbs=TRUE,
                        sampling='down',
                        summaryFunction=twoClassSummary,
                        allowParallel=TRUE,
                        verboseIter=TRUE)

param_sweep <- expand.grid(mtry=unique(floor(seq(ncol(otu_train)^.25,ncol(otu_train)^.75,length=10))))
out_rf <- train(x=otu_train,
                y=train_labels,
                method='rf',
                ntree=128,
                tuneGrid=param_sweep,
                importance=TRUE,
                proximity=TRUE,
                metric='ROC', 
                maximize=TRUE, 
                trControl=tr_ctrl,
                verbose=TRUE)

stopCluster(cl)

pred_p <- predict(out_rf,otu_test,type='prob')
pred_Y <- factor(colnames(pred_p)[apply(pred_p,1,which.max)],levels=colnames(pred_p))
pred_out <- data.frame(obs=test_labels,
                       pred=pred_Y,
                       pred_p)
two_class_summary <- twoClassSummary(pred_out,lev=levels(pred_out$pred))
confusion_matrix <- caret::confusionMatrix(pred_Y,test_labels)

roc_curve <- roc(response=test_labels, 
                 predictor=pred_p[,1], 
                 levels=rev(unique(test_labels)))

saveRDS(list(OTU=OTU,
             VOCAB=DAT$VOCAB,
             TAX=DAT$TAX,
             otu_train=otu_train,
             otu_test=otu_test,
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
        file.path(data_dir,'classify_ras.rds'))