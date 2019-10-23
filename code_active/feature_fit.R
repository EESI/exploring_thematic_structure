feature_fit <- function(num_features,out,fit,
                        train_docs,test_docs,train_meta,test_meta,
                        vocab,variable,beta=NULL){
   
   print(out$dx$imp$RF1500DS)
   
   if (is.null(beta)){
   
      topic_features <- order(out$dx$imp$RF1500DS$importance[,1],decreasing=TRUE)[1:num_features]
      fit$beta$logbeta[[3]] <- do.call('rbind',fit$beta$logbeta)[topic_features,]
      fit$settings$dim$K <- NROW(fit$beta$logbeta[[3]])
   
   }else{
      
      topic_features <- 1:NROW(beta)
      fit$beta$logbeta[[3]] <- beta
      fit$settings$dim$K <- NROW(fit$beta$logbeta[[3]])
      
   }
   
   cat('\nSelecting topics: ',topic_features,'\n\n')
   
   z_covariate <- 3
   cat('Exploring training set posterior\n')
   train_fit <- ctm_frozen(fit,train_docs,vocab,
                           seed=seed_train,max.em.its=500,emtol=1e-5,avg_iters=1,
                           verbose=TRUE,true_doc_content=TRUE,
                           data=train_meta,covariate=z_covariate,
                           parallel=TRUE,nc=25)
   if(exists('cl')) stopCluster(cl)
   
   cat('Exploring testing set posterior\n')
   test_fit <- ctm_frozen(fit,test_docs,vocab,
                          seed=seed_test,max.em.its=500,emtol=1e-5,avg_iters=1,
                          verbose=TRUE,true_doc_content=TRUE,
                          data=test_meta,covariate=z_covariate,
                          parallel=FALSE,nc=10)
   if(exists('cl')) stopCluster(cl)
   
   
   
   K <- train_fit$settings$dim$K
   
   load_fit <- FALSE
   save_fit <- FALSE
   
   betas <- beta_prep(train_fit,test_fit,counts,otu_taxa_xavier,vocab,
                      save_fit,save_dir,save_fit_foldername,dupnum,seed_permuted)
   beta_frozen_ra <- betas$beta_frozen_ra
   beta_frozen <- betas$beta_frozen
   beta_meta <- betas$beta_meta
   beta_otu <- betas$beta_otu
   
   out_new <- eval_labels(save_fit,load_fit,train_fit$Z_bar,test_fit$Z_bar,train_meta,test_meta,
                          save_dir,save_fit_foldername,save_coef_filename,beta_type=model,
                          nc=60,test='dx')
   
   colnames(train_fit$Z_bar) <- paste0('T',topic_features)
   names(out_new$dx$lasso) <- paste0('T',topic_features)
   names(out_new$dx$en1) <- paste0('T',topic_features)
   names(out_new$dx$en2) <- paste0('T',topic_features)
   rownames(out_new$dx$imp$RF1500DS$importance) <- paste0('T',topic_features)
   
   plot_z_heatmap(train_fit$Z_bar,train_meta,out_new$dx$en2,
                  dist1='jaccard',dist2='jaccard',
                  clust1='ward.D2',clust2='ward.D2',
                  transform='none',
                  main='',
                  rowclust=FALSE,variable=variable)
   
   colnames(test_fit$Z_bar) <- paste0('T',topic_features)
   
   plot_z_heatmap(test_fit$Z_bar,test_meta,out_new$dx$en2,
                           dist1='jaccard',dist2='jaccard',
                           clust1='ward.D2',clust2='ward.D2',
                           transform='none',
                           main='',
                           rowclust=FALSE,variable)
   
   return(list(out=out_new,
               z_bar_train=train_fit$Z_bar,
               train_meta=train_meta,
               z_bar_test=test_fit$Z_bar,
               test_meta=test_meta))
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


params <- expand.grid(K=c(25,50,75,125,200),
                      content=c(FALSE,TRUE),
                      variable=c('DIAGNOSIS','ISOLATION_SOURCE'))
params <- params %>% 
   filter(!(content == FALSE & variable == 'ISOLATION_SOURCE'),
          !(content == TRUE  & K == 200)) %>%
   arrange(K,variable)


param <- 2

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


prepare_framework(random_seed,K,cn_normalize,content,variable,seq_sim)


load_fits(file.path(save_dir,save_fit_foldername,save_fit_filename))

# take combined (K x M where M is the number of features, 2 for dx)
b <- 1
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

# evaluate the performance of the combined beta matrix
cat('Evaluating predictive performance for',names(fit_frozen)[b],'\n')
out <- eval_labels(save_fit,load_fit,train_fit$Z_bar,test_fit$Z_bar,train_meta,test_meta,
                   save_dir,save_fit_foldername,save_coef_filename,beta_type=model,
                   nc=60)

# take the top F topics via RF importance (tree = 1500) for predicting the target labels (dx)
# and subset the beta matrix with the indexes for these features, yielding an N x F matrix.
# Generate the topic assignments (z bar) via the subsetted beta matrix. Then, evaluluate the
# predictive performance to color code the topics via the EN output, then plot a heatmap. For
# the heatmap, the columns are ordered in terms for decreasing PCDAI score, where CD- gets a 0.
ff <- feature_fit(30,out,fit,train_docs,test_docs,train_meta,test_meta,vocab,variable)
ff$out$dx$score
ff$out$dx$imp$RF1500DS


plot_z_heatmap(ff$z_bar_train,ff$train_meta,ff$out$dx$lasso,
               dist1='jaccard',dist2='jaccard',
               clust1='ward.D2',clust2='ward.D2',
               transform='none',
               main='',
               rowclust=TRUE,variable='DIAGNOSIS')

plot_z_heatmap(ff$z_bar_test,ff$test_meta,ff$out$dx$lasso,
               dist1='jaccard',dist2='jaccard',
               clust1='ward.D2',clust2='ward.D2',
               transform='none',
               main='',
               rowclust=TRUE,variable='DIAGNOSIS')
