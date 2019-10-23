#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 2) stop('Please specify number of topics, K.', call.=FALSE)

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


K <- as.integer(args[1])
fit_idx <- as.integer(args[2])

data_dir <- '~/Dropbox/stm_microbiome/qiime_active/gevers'
models_dir <- file.path(data_dir,'models',K)
froz_fn <- sprintf('stm_froz_k_%s.rds',K)
ko_fn <- file.path(models_dir,sprintf('ko_k_%s_fit%s.biom',K,fit_idx))

biom_file <- read_biom(ko_fn)
ko <- as.matrix(biom_data(biom_file))
ko <- t(ko[rowSums(ko)>0,])

DAT <- readRDS(file.path(models_dir,froz_fn))

meta_train <- DAT$meta_train
meta_test <- DAT$meta_test

theta_train <- DAT$theta_train
theta_test <- DAT$theta_test

train_labels <- factor(ifelse(meta_train$DIAGNOSIS==1,'CD','CTRL'),levels=c('CD','CTRL'))
test_labels <- factor(ifelse(meta_test$DIAGNOSIS==1,'CD','CTRL'),levels=c('CD','CTRL'))

X_train <- theta_train %*% ko
X_test <- theta_test %*% ko

X_train <- apply(X_train,2,function(x) (x-mean(x))/(max(x)-min(x)))
X_test <- apply(X_test,2,function(x) (x-mean(x))/(max(x)-min(x)))

ncores <- 5
cl <- makeCluster(ncores,type='SOCK')
registerDoSNOW(cl)

tr_ctrl <- trainControl(method='repeatedCV', number=10, repeats=10, 
                        classProbs=TRUE,
                        sampling='up',
                        summaryFunction=twoClassSummary,
                        allowParallel=TRUE)

param_sweep <- expand.grid(mtry=unique(floor(seq(ncol(X_train)^.25,ncol(X_train)^.75,length=10))))
out_rf <- train(x=X_train,
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

pred_p <- predict(out_rf,X_test,type='prob')
pred_Y <- factor(colnames(pred_p)[apply(pred_p,1,which.max)],levels=colnames(pred_p))
pred_out <- data.frame(obs=test_labels,
                       pred=pred_Y,
                       pred_p)
two_class_summary <- twoClassSummary(pred_out,lev=levels(pred_out$pred))
confusion_matrix <- caret::confusionMatrix(pred_Y,test_labels)

roc_curve <- roc(response=test_labels, 
                 predictor=pred_p[,1], 
                 levels=rev(unique(test_labels)))
#plot(ROC, col='red',lwd = 2)

classify_fn <- sprintf('classify_ko_k_%s_fit%s.rds',K,fit_idx)
saveRDS(list(K=K,
             OTU=DAT$OTU,
             VOCAB=DAT$VOCAB,
             TAX=DAT$TAX,
             theta_train=theta_train,
             theta_test=theta_test,
             ko=ko,
             X_train=X_train,
             X_test=X_test,
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
             docs_train=DAT$docs_train,
             meta_test=DAT$meta_test,
             docs_test=DAT$docs_test),
        file.path(models_dir,classify_fn))