heatmap_dump <- function(save_output,model=NULL){
   
   if (!is.null(model)) output_dir <- file.path(output_dir,model)
   dir.create(file.path(output_dir),showWarnings=FALSE)
   
   coef_k <- out$dx$lasso
   coef_sig_idx <- which(coef_k != 0)
   coef_sig <- coef_k[coef_sig_idx]
   names(coef_sig) <- coef_sig_idx
   
   if (save_output) pdf(file.path(output_dir,'zbar_heatmap.pdf'), height=25, width=25)
   plot_z_heatmap(train_fit$Z_bar,train_meta,coef_k,
                  dist1='jaccard',dist2='jaccard',
                  clust1='ward.D2',clust2='ward.D2',
                  transform='none',
                  main='Training',
                  rowclust=TRUE,variable='DIAGNOSIS')
   plot_z_heatmap(train_fit$Z_bar,train_meta,coef_k,
                  dist1='jaccard',dist2='jaccard',
                  clust1='ward.D2',clust2='ward.D2',
                  transform='none',
                  main='Training',
                  rowclust=TRUE,variable='ISOLATION_SOURCE')
   plot_z_heatmap(test_fit$Z_bar,test_meta,coef_k,
                  dist1='jaccard',dist2='jaccard',
                  clust1='ward.D2',clust2='ward.D2',
                  transform='none',
                  main='Testing',
                  rowclust=TRUE,variable='DIAGNOSIS')
   plot_z_heatmap(test_fit$Z_bar,test_meta,coef_k,
                  dist1='jaccard',dist2='jaccard',
                  clust1='ward.D2',clust2='ward.D2',
                  transform='none',
                  main='Testing',
                  rowclust=TRUE,variable='ISOLATION_SOURCE')
   
   plot_z_heatmap(train_fit$Z_bar,train_meta,coef_k,
                  dist1='jaccard',dist2='jaccard',
                  clust1='ward.D2',clust2='ward.D2',
                  transform='none',
                  main='Training',
                  rowclust=FALSE,variable='DIAGNOSIS')
   plot_z_heatmap(train_fit$Z_bar,train_meta,coef_k,
                  dist1='jaccard',dist2='jaccard',
                  clust1='ward.D2',clust2='ward.D2',
                  transform='none',
                  main='Training',
                  rowclust=FALSE,variable='ISOLATION_SOURCE')
   plot_z_heatmap(test_fit$Z_bar,test_meta,coef_k,
                  dist1='jaccard',dist2='jaccard',
                  clust1='ward.D2',clust2='ward.D2',
                  transform='none',
                  main='Testing',
                  rowclust=FALSE,variable='DIAGNOSIS')
   plot_z_heatmap(test_fit$Z_bar,test_meta,coef_k,
                  dist1='jaccard',dist2='jaccard',
                  clust1='ward.D2',clust2='ward.D2',
                  transform='none',
                  main='Testing',
                  rowclust=FALSE,variable='ISOLATION_SOURCE')
   if (save_output) dev.off()
   
   
}


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

#
# params <- expand.grid(K=c(25,50,75,125,200),
#             content=c(FALSE,TRUE),
#             variable=c('DIAGNOSIS','ISOLATION_SOURCE'))
# params <- params %>% 
#    filter(!(content == FALSE & variable == 'ISOLATION_SOURCE')) %>%
#    arrange(K,variable)

params <- expand.grid(K=c(25,50,75,125,200),
                      content=c(FALSE,TRUE),
                      variable=c('DIAGNOSIS','ISOLATION_SOURCE'))
params <- params %>% 
   filter(!(content == FALSE & variable == 'ISOLATION_SOURCE'),
          !(content == TRUE  & K == 200)) %>%
   arrange(K,variable)


for (param in 1:NROW(params)){
   #START LOOP FOR PARAMETERS
   #
   
   
   
   
   source('~/Dropbox/stm_microbiome/code_active/stm_functions.R')
   source('~/Dropbox/stm_microbiome/code_active/nav_froz_fxns_3.R') 
   source('~/Dropbox/stm_microbiome/code_active/performance_1.R')
   source('~/Dropbox/stm_microbiome/code_active/framework.R')
   
   
   load_fit <- TRUE
   save_fit <- FALSE
   save_output <- TRUE
   
   random_seed <- FALSE
   K <- params[param,]$K
   cn_normalize <- TRUE
   content <- params[param,]$content
   seq_sim <- 's97'
   variable <- as.character(params[param,]$variable)
   dupnum <- NULL
   
   
   
   temp_dump <- tempfile()
   if (save_output) sink(temp_dump)
   prepare_framework(random_seed,K,cn_normalize,content,variable,seq_sim)
   if (save_output) sink()
   
   
   
   
   if (load_fit){
      
      load_fits(file.path(save_dir,save_fit_foldername,save_fit_filename))
      
   }else{
      
      if(content){
         fit <- stm::stm(train_docs, vocab, K=K,
                         max.em.its=500, 
                         init.type='Spectral',
                         data=train_meta,content=as.formula(paste0('~',variable)),
                         seed=seed_fit,verbose=TRUE,reportevery=25)
      }else{
         fit <- stm::stm(train_docs, vocab, K=K,
                         max.em.its=500, 
                         init.type='Spectral',
                         seed=seed_fit,verbose=TRUE,reportevery=25)
      }
      
      if (content){
         
         fit_var_names <- c('combined',fit$settings$covariates$yvarlevels)
         fit_frozen <- list() 
         for (m in 0:(length(fit_var_names)-1)){
            
            z_covariate <- m
            cat('Exploring training set posterior for',fit_var_names[m+1],'\n')
            train_fit <- ctm_frozen(fit,train_docs,vocab,
                                    seed=seed_train,max.em.its=500,emtol=1e-5,avg_iters=1,
                                    verbose=TRUE,true_doc_content=TRUE,
                                    data=train_meta,covariate=z_covariate,
                                    parallel=TRUE,nc=25)
            if(exists('cl')) stopCluster(cl)
            
            cat('Exploring testing set posterior for',fit_var_names[m+1],'\n')
            test_fit <- ctm_frozen(fit,test_docs,vocab,
                                   seed=seed_test,max.em.its=500,emtol=1e-5,avg_iters=1,
                                   verbose=TRUE,true_doc_content=TRUE,
                                   data=test_meta,covariate=z_covariate,
                                   parallel=FALSE,nc=10)
            if(exists('cl')) stopCluster(cl)
            
            fit_frozen[[fit_var_names[m+1]]] <- list(train=train_fit,test=test_fit)
         }
         
      }else{
         
         train_fit <- ctm_frozen(fit,train_docs,vocab,
                                 seed=seed_train,max.em.its=500,emtol=1e-5,avg_iters=1,
                                 verbose=TRUE,true_doc_content=TRUE,
                                 data=train_meta,covariate=NULL,
                                 parallel=TRUE,nc=25)
         if(exists('cl')) stopCluster(cl)
         
         test_fit <- ctm_frozen(fit,test_docs,vocab,
                                seed=seed_test,max.em.its=500,emtol=1e-5,avg_iters=1,
                                verbose=TRUE,true_doc_content=TRUE,
                                data=test_meta,covariate=NULL,
                                parallel=FALSE,nc=10)
         if(exists('cl')) stopCluster(cl)
         
         fit_frozen <- list(train_fit,test_fit)
      }
      
      if (save_fit) {
         dir.create(file.path(save_dir,save_fit_foldername),showWarnings=FALSE)
         if (file.exists(file.path(save_dir,save_fit_foldername,save_fit_filename))){
            dupnum <- sample(1:sqrt(as.integer(Sys.time())),1)
            save_fit_filename <- paste0('dup',dupnum,'_',save_fit_filename)
         }
         saveRDS(list(fit=fit,fit_frozen=fit_frozen),
                 file.path(save_dir,save_fit_foldername,save_fit_filename))
      }
      
   }
   
   
   if (save_output) {
      output_dir <- file.path('~/Dropbox/stm_microbiome/output_active/Gevers',save_fit_foldername)
      dir.create(output_dir,showWarnings=FALSE)
      file.copy(temp_dump,file.path(output_dir,'info.txt'))
   }
   
   
   if (!content){
      train_fit <- fit_frozen[[1]]
      test_fit <- fit_frozen[[2]]
      K <- train_fit$settings$dim$K
      
      betas <- beta_prep(train_fit,test_fit,counts,otu_taxa_xavier,vocab,
                         save_fit,save_dir,save_fit_foldername,dupnum,seed_permuted)
      beta_frozen_ra <- betas$beta_frozen_ra
      beta_frozen <- betas$beta_frozen
      beta_meta <- betas$beta_meta
      beta_otu <- betas$beta_otu
      
      out <- eval_labels(save_fit,load_fit,train_fit$Z_bar,test_fit$Z_bar,train_meta,test_meta,
                         save_dir,save_fit_foldername,save_coef_filename,
                         nc=60)
      
      heatmap_dump(save_output)
      
   }else{
      
      for (b in seq_along(fit_frozen)){
         
         train_fit <- fit_frozen[[b]][['train']]
         test_fit <- fit_frozen[[b]][['test']]
         K <- train_fit$settings$dim$K
         
         model <- str_replace_all(names(fit_frozen)[b],' ','')
         betas <- beta_prep(train_fit,test_fit,counts,otu_taxa_xavier,vocab,
                            save_fit,save_dir,save_fit_foldername,dupnum,seed_permuted,model)
         beta_frozen_ra <- betas$beta_frozen_ra
         beta_frozen <- betas$beta_frozen
         beta_meta <- betas$beta_meta
         beta_otu <- betas$beta_otu
         
         cat('Evaluating predictive performance for',names(fit_frozen)[b],'\n')
         out <- eval_labels(save_fit,load_fit,train_fit$Z_bar,test_fit$Z_bar,train_meta,test_meta,
                            save_dir,save_fit_foldername,save_coef_filename,beta_type=model,
                            nc=60)
         
         heatmap_dump(save_output,model)
      }
      
   }
   
   
   
   
   
   rm(list = ls()[!(ls() %in% c('param','params','heatmap_dump'))]) ####
} ## end save loop ##
