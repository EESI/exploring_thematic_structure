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
library(vegan)


params <- expand.grid(K=c(25,50,75),
                      content=c(FALSE),
                      variable=c('DIAGNOSIS'))
params <- params %>% 
   filter(!(content == TRUE  & K == 200)) %>%
   arrange(K,variable)


beta_frozens <- list()
for (param in 1:NROW(params)){
   #START LOOP FOR PARAMETERS
   #
   
   
   source('~/Dropbox/stm_microbiome/code_active/stm_functions.R')
   source('~/Dropbox/stm_microbiome/code_active/nav_froz_fxns_3.R') 
   source('~/Dropbox/stm_microbiome/code_active/performance_1.R')
   source('~/Dropbox/stm_microbiome/code_active/framework.R')
   
   
   load_fit <- TRUE
   save_fit <- FALSE
   save_output <- FALSE
   
   random_seed <- FALSE
   K <- params[param,]$K
   cn_normalize <- TRUE
   content <- params[param,]$content
   seq_sim <- 's97'
   variable <- as.character(params[param,]$variable)
   dupnum <- NULL
   test <- 'dx'   
   
   
   prepare_framework_ti_rarefy(random_seed,K,cn_normalize,content,variable,seq_sim)
   

   load_fits(file.path(save_dir,save_fit_foldername,save_fit_filename))
   
   
   train_fit <- fit_frozen[[1]]
   test_fit <- fit_frozen[[2]]
   K <- train_fit$settings$dim$K
   
   betas <- beta_prep(train_fit,test_fit,counts,Nmin,otu_taxa_xavier,vocab,
                      save_fit,save_dir,save_fit_foldername,dupnum,seed_permuted)
   
   beta_frozen_ra <- betas$beta_frozen_ra
   beta_frozen <- betas$beta_frozen
   beta_meta <- betas$beta_meta
   beta_otu <- betas$beta_otu      
   
   beta_frozens[[param]] <- beta_frozen


   rm(list = ls()[!(ls() %in% c('param','params','beta_frozens'))]) ####
} ## end save loop ####




beta_frozent25 <- beta_frozens[[1]]
shant25 <- diversity(t(beta_frozent25))
event25 <- shant25/log(specnumber(t(beta_frozent25)))

beta_frozent50 <- beta_frozens[[2]]
shant50 <- diversity(t(beta_frozent50))
event50 <- shant50/log(specnumber(t(beta_frozent50)))

beta_frozent75 <- beta_frozens[[3]]
shant75 <- diversity(t(beta_frozent75))
event75 <- shant75/log(specnumber(t(beta_frozent75)))



data.frame(Shannon=c(shant25,shant50,shant75,event25,event50,event75),
           K=rep(rep(c('T25','T50','T75'),c(25,50,75)),2),
           Metric=rep(c('Shannon Diveristy','Pielou Evenness'),each=150)) %>%
   ggplot(aes(x=K,y=Shannon,fill=K)) + geom_boxplot() + facet_wrap(~Metric,scales='free_y') +
   theme(legend.position='none')

