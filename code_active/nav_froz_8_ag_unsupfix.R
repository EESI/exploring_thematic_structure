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


### for sex_site, only the supervised stuff is right
### for some reason, the vocabs are off for the unsupervised 
### so get the output for the supervised, then rerun the fits
### for the unsupervised, but ensure the vocabs match


### NOTE I TURNED PARALLEL OFF FOR CTM_FROZEN!!!!

params <- expand.grid(K=c(50),
                      content=c(FALSE,TRUE),
                      variable=c('DIET_TYPE'))
params <- params %>% 
   filter(!(content == FALSE  & K == 25)) %>%
   arrange(K,variable)


for (param in 1:NROW(params)){
   #START LOOP FOR PARAMETERS
   #
   
   
   
   source('~/Dropbox/stm_microbiome/code_active/stm_functions.R')
   source('~/Dropbox/stm_microbiome/code_active/nav_froz_fxns_3.R') 
   source('~/Dropbox/stm_microbiome/code_active/performance_1.R')
   source('~/Dropbox/stm_microbiome/code_active/framework.R')
   
   
   load_fit <- FALSE
   save_fit <- TRUE
   save_output <- FALSE
   
   random_seed <- FALSE
   K <- params[param,]$K
   cn_normalize <- TRUE
   content <- params[param,]$content
   seq_sim <- 's97'
   variable <- as.character(params[param,]$variable)
   dupnum <- NULL
   
   
   ### Removing samples in bottom 15% percentile of total sample reads (< 975.75 reads)
   ### Capping at 10,000
   temp_dump <- tempfile()
   if (save_output) sink(temp_dump)
   prepare_framework_ag_diet_clean(random_seed,K,cn_normalize,content,variable,seq_sim)
   if (save_output) sink()
   
   
   if (load_fit){
      
      load_fits(file.path(save_dir,save_fit_foldername,save_fit_filename))
      
   }else{
      
      if (content){
         fit <- tryCatch(stm(content=as.formula(paste0('~',variable)),
                             documents=train_docs, vocab=vocab, K=K,
                             data=train_meta,
                             max.em.its=500, 
                             init.type='Spectral',
                             # ngroups=10,
                             seed=seed_fit,verbose=TRUE,reportevery=25),
                         error = function(e) e)
      }else{
         fit <- tryCatch(stm(documents=train_docs, vocab=vocab, K=K,
                             data=train_meta,
                             max.em.its=500, 
                             init.type='Spectral',
                             # ngroups=10,
                             seed=seed_fit,verbose=TRUE,reportevery=25),
                         error = function(e) e)         
      }
      
      if (inherits(fit, 'error')){
         cat('Convergence failure. Moving to next model.\n',sep='')
         rm(list = ls()[!(ls() %in% c('param','params'))]) 
         next            
      }
      
      
      if (content){
         
         fit_var_names <- c('combined',fit$settings$covariates$yvarlevels)
         fit_frozen <- list() 
         for (m in 0:(length(fit_var_names)-1)){
            
            z_covariate <- m
            cat('Exploring training set posterior for',fit_var_names[m+1],'\n')
            train_fit <- tryCatch(ctm_frozen(fit,train_docs,vocab,
                                             seed=seed_train,max.em.its=500,emtol=1e-5,avg_iters=1,
                                             verbose=TRUE,true_doc_content=TRUE,
                                             data=train_meta,covariate=z_covariate,
                                             parallel=TRUE,nc=20),
                                  error = function(e) e)
            if(exists('cl')) stopCluster(cl)
            
            cat('Exploring testing set posterior for',fit_var_names[m+1],'\n')
            test_fit <- tryCatch(ctm_frozen(fit,test_docs,vocab,
                                            seed=seed_test,max.em.its=500,emtol=1e-5,avg_iters=1,
                                            verbose=TRUE,true_doc_content=TRUE,
                                            data=test_meta,covariate=z_covariate,
                                            parallel=FALSE,nc=10),
                                 error = function(e) e)
            if(exists('cl')) stopCluster(cl)
            
            if (inherits(train_fit, 'error') | inherits(test_fit, 'error')){
               cat('Convergence failure for Z-bar. Moving to next model.\n',sep='')
               fit_frozen[[fit_var_names[m+1]]] <- list(train=NULL,test=NULL)
               next            
            }else{
               fit_frozen[[fit_var_names[m+1]]] <- list(train=train_fit,test=test_fit)
            }
            
         }
         
      }else{
         
         train_fit <- tryCatch(ctm_frozen(fit,train_docs,vocab,
                                          seed=seed_train,max.em.its=500,emtol=1e-5,avg_iters=1,
                                          verbose=TRUE,true_doc_content=TRUE,
                                          data=train_meta,covariate=NULL,
                                          parallel=TRUE,nc=20),
                               error = function(e) e)
         if(exists('cl')) stopCluster(cl)
         
         test_fit <- tryCatch(ctm_frozen(fit,test_docs,vocab,
                                         seed=seed_test,max.em.its=500,emtol=1e-5,avg_iters=1,
                                         verbose=TRUE,true_doc_content=TRUE,
                                         data=test_meta,covariate=NULL,
                                         parallel=FALSE,nc=10),
                              error = function(e) e)
         if(exists('cl')) stopCluster(cl)
         
         if (inherits(train_fit, 'error') | inherits(test_fit, 'error')){
            cat('Convergence failure for Z-bar.\n',sep='')
            fit_frozen <- list(train=NULL,test=NULL)            
         }else{
            fit_frozen <- list(train_fit,test_fit)
         }
      }
      
   }
   
   
   
   if (save_fit){
      dir.create(file.path(save_dir,save_fit_foldername),showWarnings=FALSE,recursive=TRUE)
      if (file.exists(file.path(save_dir,save_fit_foldername,save_fit_filename))){
         dupnum <- sample(1:sqrt(as.integer(Sys.time())),1)
         save_fit_filename <- paste0('dup',dupnum,'_',save_fit_filename)
      }
      saveRDS(list(fit=fit,fit_frozen=fit_frozen),
              file.path(save_dir,save_fit_foldername,save_fit_filename))
   }
   
   
   
   if (save_output) {
      output_dir <- file.path('~/Dropbox/stm_microbiome/output_active/AG_diet',save_fit_foldername)
      dir.create(output_dir,showWarnings=FALSE,recursive=TRUE)
      file.copy(temp_dump,file.path(output_dir,'info.txt'))
   }
   
   
   if (!content){
      
      if (any(unlist(lapply(fit_frozen,is.null)))) next
      
      train_fit <- fit_frozen[[1]]
      test_fit <- fit_frozen[[2]]
      K <- train_fit$settings$dim$K
      
      betas <- beta_prep(train_fit,test_fit,counts,otu_taxa_ag,vocab,
                         save_fit,save_dir,save_fit_foldername,dupnum,seed_permuted)
      
   }else{
      
      for (b in seq_along(fit_frozen)){
         
         if (any(unlist(lapply(fit_frozen[[b]],is.null)))) next
         
         train_fit <- fit_frozen[[b]][['train']]
         test_fit <- fit_frozen[[b]][['test']]
         K <- train_fit$settings$dim$K
         
         model <- str_replace_all(names(fit_frozen)[b],' ','')
         betas <- beta_prep(train_fit,test_fit,counts,otu_taxa_ag,vocab,
                            save_fit,save_dir,save_fit_foldername,dupnum,seed_permuted,model)
      }
      
   }
   
   
   
   
   
   rm(list = ls()[!(ls() %in% c('param','params'))]) ####
} ## end save loop ####
