#!/usr/bin/env Rscript

# args <- commandArgs(trailingOnly=TRUE)
# if (length(args) != 1) stop("At least one argument must be supplied (input file).n", call.=FALSE)

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

# model <- as.integer(args[1])
# seed_split <- as.integer(args[2])


require(snow)
require(doSNOW)
ncores <- 48
cl <- makeCluster(ncores,type='SOCK')
registerDoSNOW(cl)


for (seed_split in c(1,105,12,53,74)){
  for (model in 1:2){
    
    seed_general <- 1654
    set.seed(seed_general)
    
    
    dat1 <- readRDS(sprintf('~/Dropbox/stm_microbiome/data_active/AG_female/stm_s97_rarefied_10000_supervised/features/seed_%s/fits_K_35_model_%s.rds',
                           seed_split,1))
    dat2 <- readRDS(sprintf('~/Dropbox/stm_microbiome/data_active/AG_female/stm_s97_rarefied_10000_supervised/features/seed_%s/fits_K_35_model_%s.rds',
                           seed_split,2))
    
    if (!(identical(dat1$meta$meta_train$SampleID,dat2$meta$meta_train$SampleID) & 
            identical(dat1$meta$meta_test$SampleID,dat2$meta$meta_test$SampleID))){
      stop('Training and Testing set samples are different!')
    }
    
    dat <- readRDS(sprintf('~/Dropbox/stm_microbiome/data_active/AG_female/stm_s97_rarefied_10000_supervised/features/seed_%s/fits_K_35_model_%s.rds',
                           seed_split,model))
    
    fit <- dat$fits$fit
    meta <- dat$meta$meta_train
    K <- dat$fits$fit$settings$dim$K
    
    print(fit$settings$call)
    
    eff <- estimateEffect(1:K ~ FEMALE, fit, meta, uncertainty='None')
    eff_plot_female <- plot.estimateEffect(eff, 'FEMALE', model=fit, topics=1:K, method='difference',cov.value1=1,cov.value2=0)
    coef_k <- unlist(eff_plot_female$means)
    coef_k <- coef_k/sd(coef_k)
    
    biom_file <- read_biom(sprintf('~/Dropbox/stm_microbiome/data_active/AG_female/stm_s97_rarefied_10000_supervised/features/seed_%s/beta_K_35_model_%s_predicted_metagenome.biom',
                           seed_split,model))
    kegg_metadata_risk <- kegg_pw_collapse(biom_file,3)
    kegg_risk <- as.matrix(biom_data(biom_file))
    
    
    
    
    if (NCOL(kegg_risk) == K) {kegg_risk1 <- kegg_risk; kegg_risk2 <- NULL} else {kegg_risk1 <- kegg_risk[,1:K]; kegg_risk2 <- kegg_risk[,(K+1):NCOL(kegg_risk)]}
    
    
    
    
    
    kegg_risk <- kegg_risk1
    options(scipen=10000)
    abund_thresh <- 2
    abund_by_filt_x <- seq(1,1000,1)
    abund_by_filt <- sapply(abund_by_filt_x,function(i) sum(rowSums(kegg_risk) >= i))
    data.frame(x=abund_by_filt_x,y=abund_by_filt) %>%
      filter(x < 5000) %>%
      ggplot(aes(x,y)) + geom_point() + geom_vline(xintercept=abund_thresh,color='red')
    Y <- rowSums(kegg_risk > 0)+1
    X <- rowSums(kegg_risk)+1
    qplot(X,Y,geom='point',alpha=.1) + xlab('Log10 Abundance') + ylab('Prevalence') +
      scale_x_log10() + theme(legend.position='none') + geom_vline(xintercept=abund_thresh,color='red')
    
    kegg_risk <- kegg_risk[rowSums(kegg_risk) >= abund_thresh,]
    kegg_gene_names <- kegg_pw_collapse(biom_file,0)
    kegg_gene_names <- lapply(kegg_gene_names,function(x) paste0(x,collapse=' | '))
    
    DF <- data.frame(coef=coef_k) %>%
      mutate(assoc=ifelse(coef>0,'F',ifelse(coef<0,'M','ZERO')),
             coef_scaled=coef_k/sd(coef_k),
             coef_ranked=rank(coef_k))
    rownames(DF) <- colnames(kegg_risk)
    
    KEGG <- otu_table(kegg_risk,taxa_are_rows=TRUE)
    SAMP <- sample_data(DF)
    PS <- phyloseq(KEGG,SAMP)
    
    DS2 <- phyloseq_to_deseq2(PS, ~ coef_scaled) # try scaled
    DS2 <- DESeq2::DESeq(DS2, test="Wald", fitType="parametric")
    res <- DESeq2::results(DS2, cooksCutoff=FALSE, pAdjustMethod='BH')
    
    
    alpha <- 1 #.01
    res <- res[order(res$padj,decreasing=FALSE)[1:500],] ### to ensure top 500 and not more
    sigtab <- res[which(res$padj < alpha), ]
    
    kegg_sig <- kegg_risk[rownames(sigtab),]
    #rownames(kegg_sig) <- kegg_gene_names[rownames(kegg_sig)]
    
    ko_names <- rownames(kegg_sig) 
    
    #misc_filter <- !(rownames(kegg_sig) %in% c('None','hypothetical protein'))
    kegg_sig1 <- kegg_sig #kegg_sig[misc_filter,]
    ko_names1 <- ko_names #ko_names[misc_filter]
    
    ###
    ###
    ###
    
    if (!is.null(kegg_risk2)){
      kegg_risk <- kegg_risk2
      options(scipen=10000)
      abund_thresh <- 2
      abund_by_filt_x <- seq(1,1000,1)
      abund_by_filt <- sapply(abund_by_filt_x,function(i) sum(rowSums(kegg_risk) >= i))
      data.frame(x=abund_by_filt_x,y=abund_by_filt) %>%
        filter(x < 200) %>%
        ggplot(aes(x,y)) + geom_point() + geom_vline(xintercept=abund_thresh,color='red')
      Y <- rowSums(kegg_risk > 0)+1
      X <- rowSums(kegg_risk)+1
      qplot(X,Y,geom='point',alpha=.1) + xlab('Log10 Abundance') + ylab('Prevalence') +
        scale_x_log10() + theme(legend.position='none') + geom_vline(xintercept=abund_thresh,color='red')
      
      kegg_risk <- kegg_risk[rowSums(kegg_risk) >= abund_thresh,]
      kegg_gene_names <- kegg_pw_collapse(biom_file,0)
      kegg_gene_names <- lapply(kegg_gene_names,function(x) paste0(x,collapse=' | '))
      
      DF <- data.frame(coef=coef_k) %>%
        mutate(assoc=ifelse(coef>0,'F',ifelse(coef<0,'M','ZERO')),
               coef_scaled=coef_k/sd(coef_k),
               coef_ranked=rank(coef_k))
      rownames(DF) <- colnames(kegg_risk)
      
      KEGG <- otu_table(kegg_risk,taxa_are_rows=TRUE)
      SAMP <- sample_data(DF)
      PS <- phyloseq(KEGG,SAMP)
      
      DS2 <- phyloseq_to_deseq2(PS, ~ coef_scaled) ## switch to ranked
      DS2 <- DESeq2::DESeq(DS2, test="Wald", fitType="parametric")
      res <- DESeq2::results(DS2, cooksCutoff=FALSE, pAdjustMethod='BH')
      alpha <- 1 #.01
      res <- res[order(res$padj,decreasing=FALSE)[1:500],] ### to ensure top 500 and not more
      sigtab <- res[which(res$padj < alpha), ]
      
      kegg_sig <- kegg_risk[rownames(sigtab),]
      #rownames(kegg_sig) <- kegg_gene_names[rownames(kegg_sig)]
      
      ko_names <- rownames(kegg_sig) 
      
      #misc_filter <- !(rownames(kegg_sig) %in% c('None','hypothetical protein'))
      kegg_sig2 <- kegg_sig #kegg_sig[misc_filter,]
      ko_names2 <- ko_names #ko_names[misc_filter]
    }
    
    ##
    ##
    
    ko1 <- rownames(kegg_sig1)
    if (!is.null(kegg_risk2)) ko2 <- rownames(kegg_sig2) else ko2 <- NULL
    
    ko <- unique(c(ko1,ko2))
    
    
    
    
    
    
    
    
    
    
    
    # raw <- readRDS('~/Dropbox/stm_microbiome/data_active/AG_female/stm_s97_rarefied_10000_supervised/features/raw_kegg.rds')
    # raw$kegg <- raw$kegg[rowSums(raw$kegg) != 0,]
    # 
    # kegg_train <- t(raw$kegg[,dat$meta_train$SampleID])
    # meta_train <- raw$meta[dat$meta_train$SampleID,]
    # kegg_test <- t(raw$kegg[,dat$meta_test$SampleID])
    # meta_test <- raw$meta[dat$meta_test$SampleID,]
    # 
    # kegg_train <- apply(kegg_train,2,qnormalize)
    # kegg_test <- apply(kegg_test,2,qnormalize)
    # 
    # labels_train <- meta_train$SEX
    # labels_test <- meta_test$SEX
    # 
    # require(snow)
    # require(doSNOW)
    # ncores <- 48
    # cl <- makeCluster(ncores,type='SOCK')
    # registerDoSNOW(cl)
    # rf_ctrl <- trainControl(method='cv',classProbs=TRUE,
    #                         summaryFunction=twoClassSummary,
    #                         allowParallel=TRUE)
    # out_rf <- train(x=kegg_train,y=labels_train,
    #                 method='rf',ntree=1200,tuneLength=5,metric='ROC',trControl=rf_ctrl,
    #                 importance=TRUE)
    # 
    # pred_rf <- predict(out_rf,kegg_test,type='prob')[,1]
    # cm_rf <- caret::confusionMatrix(ifelse(pred_rf > .5,'female','male'),labels_test)
    # imp <- varImp(out_rf)
    # 
    # save_list_rf <- list(out=out_rf,
    #                      pred=pred_rf,
    #                      cm=cm_rf,
    #                      imp=imp,
    #                      labels=list(train=labels_train,
    #                                  test=labels_test),
    #                      kegg=list(train=kegg_train,
    #                                test=kegg_test),
    #                      meta=list(train=meta_train,
    #                                test=meta_test))
    # saveRDS(save_list_rf,sprintf('~/Dropbox/stm_microbiome/data_active/AG_female/stm_s97_rarefied_10000_supervised/features/seed_%s/rf_features.rds',
    #                              seed_split))
    
    
    
    
    
    
    
    
    
    rf_dat <- readRDS(sprintf('~/Dropbox/stm_microbiome/data_active/AG_female/stm_s97_rarefied_10000_supervised/features/seed_%s/rf_features.rds',
                              seed_split))
    
    imp <- rf_dat$imp
    kegg_train <- rf_dat$kegg$train
    kegg_test <- rf_dat$kegg$test
    labels_train <- rf_dat$labels$train
    labels_test <- rf_dat$labels$test
    
    if (!(identical(sort(rownames(kegg_train)),sort(dat1$meta$meta_train$SampleID)) &
        identical(sort(rownames(kegg_test)),sort(dat1$meta$meta_test$SampleID)))){
      stop('Raw samples are different than STM samples!')
    }
    
    
  
    
    #ko <- ko[ko %in% colnames(kegg_train)]
    rf_feat <- rownames(imp$importance)[order(imp$importance$female,decreasing=TRUE)[1:length(ko)]]
    
    
    
    
    kegg_train_feat <- kegg_train[,ko]
    kegg_test_feat <- kegg_test[,ko]
    
    
    
    rf_ctrl <- trainControl(method='cv',classProbs=TRUE,
                            summaryFunction=twoClassSummary,
                            allowParallel=TRUE)
    out_stm_feat <- train(x=kegg_train_feat,y=labels_train,
                         method='rf',ntree=1200,tuneLength=5,metric='ROC',trControl=rf_ctrl,
                         importance=FALSE)
    
    pred_stm_feat <- predict(out_stm_feat,kegg_test,type='prob')[,1]
    cm_stm_feat <- caret::confusionMatrix(ifelse(pred_stm_feat > .5,'female','male'),labels_test)
    print(cm_stm_feat)
    
    
    

    
    
    
    # kegg_table <- read_biom('~/Dropbox/stm_microbiome/data_active/AG_female/stm_s97_rarefied_10000_supervised/model/otu_table_predicted_metagenome.biom')
    # kegg <- as.matrix(biom_data(kegg_table))
    # meta <- sample_metadata(kegg_table)
    # rm(kegg_table)
    # saveRDS(list(kegg=kegg,meta=meta),'~/Dropbox/stm_microbiome/data_active/AG_female/stm_s97_rarefied_10000_supervised/features/raw_kegg.rds')
    
    
    # raw <- readRDS('~/Dropbox/stm_microbiome/data_active/AG_female/stm_s97_rarefied_10000_supervised/features/raw_kegg.rds')
    # raw$kegg <- raw$kegg[rowSums(raw$kegg) != 0,]
    # 
    # kegg_train <- t(raw$kegg[,meta_train$SampleID])
    # meta_train <- raw$meta[meta_train$SampleID,]
    # kegg_test <- t(raw$kegg[,meta_test$SampleID])
    # meta_test <- raw$meta[meta_test$SampleID,]
    # 
    # kegg_train <- apply(kegg_train,2,qnormalize)
    # kegg_test <- apply(kegg_test,2,qnormalize)
    # 
    # # kegg_train$labels <- as.factor(meta_train$SEX)
    # # kegg_test$labels <- as.factor(meta_test$SEX)
    # 
    # labels_train <- meta_train$SEX
    # labels_test <- meta_test$SEX
    # 
    # require(snow)
    # require(doSNOW)
    # ncores <- 48
    # cl <- makeCluster(ncores,type='SOCK')
    # registerDoSNOW(cl)
    # rf_ctrl <- trainControl(method='cv',classProbs=TRUE,
    #                         summaryFunction=twoClassSummary,
    #                         allowParallel=TRUE)
    # out_rf <- train(x=kegg_train,y=labels_train,
    #                 method='rf',ntree=1200,tuneLength=5,metric='ROC',trControl=rf_ctrl,
    #                 importance=TRUE)
    # 
    # pred_rf <- predict(out_rf,kegg_test,type='prob')[,1]
    # cm_rf <- caret::confusionMatrix(ifelse(pred_rf > .5,'female','male'),labels_test)
    # imp <- varImp(out_rf)
    # 
    # save_list_rf <- list(out=out_rf,
    #                      pred=pred_rf,
    #                      cm=cm_rf,
    #                      imp=imp,
    #                      labels=list(train=labels_train,
    #                                  test=labels_test),
    #                      kegg=list(train=kegg_train,
    #                                test=kegg_test),
    #                      meta=list(train=meta_train,
    #                                test=meta_test))
    # saveRDS(save_list_rf,'~/Dropbox/stm_microbiome/data_active/AG_female/stm_s97_rarefied_10000_supervised/features/rf_features.rds')
    
    
    
    
    
    
    rf_ctrl <- trainControl(method='cv',classProbs=TRUE,
                            summaryFunction=twoClassSummary,
                            allowParallel=TRUE)
    out_rf_feat <- train(x=kegg_train_feat,y=labels_train,
                    method='rf',ntree=1200,tuneLength=5,metric='ROC',trControl=rf_ctrl,
                    importance=FALSE)
    
    pred_rf_feat <- predict(out_rf_feat,kegg_test,type='prob')[,1]
    cm_rf_feat <- caret::confusionMatrix(ifelse(pred_rf_feat > .5,'female','male'),labels_test)
    print(cm_rf_feat)
    
    
    
    
    save_list <- list(rf=list(rf=out_rf_feat,stm=out_stm_feat),
                      cm=list(rf=cm_rf_feat,stm=cm_stm_feat),
                      pred=list(rf=pred_rf_feat,stm=pred_stm_feat),
                      seed=list(general=seed_general,split=seed_split))
    saveRDS(save_list,sprintf('~/Dropbox/stm_microbiome/data_active/AG_female/stm_s97_rarefied_10000_supervised/features/seed_%s/feature_perf_model_%s.rds',
                              seed_split,model))

  }
}

stopCluster(cl)