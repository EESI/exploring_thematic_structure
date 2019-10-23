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



beta_otu <- otu_taxa_xavier[as.character(counts$ids[vocab,'long']),]

save_folder <- '~/Dropbox/stm_microbiome/data_active/Gevers_ti/supervised'

fit_files <- list.files(save_folder)
fit_files <- fit_files[str_detect(fit_files,'model_.*.rds')]

for (i in seq_along(fit_files)){
   
   file_name <- fit_files[i]
   type_name <- str_match(file_name,'^.*(type\\d.*?).rds')[2]
   fit <- readRDS(file.path(save_folder,file_name))
   
   # gave this some thought
   # i could sample the counts from beta using the rarefyling value,
   # but I'd rather see E(beta), i.e., the average
   K <- fit$settings$dim$K
   beta <- floor(exp(do.call('rbind',fit$beta$logbeta))*Nmin)
   Ncontent <- NROW(beta)/K
   
   
   if (Ncontent==1) Ycontent <- rep('None',NROW(beta)) else Ycontent <- rep(fit$settings$covariates$yvarlevels,each=K)
   Ycontent <- str_replace(Ycontent,' ','')
   topic_names <- paste0('T',rep(1:K,Ncontent))
   
   colnames(beta) <- as.character(counts$ids[vocab,'long'])
   rownames(beta) <- paste0(topic_names,'_',Ycontent)
   beta_meta <- as.matrix(data.frame(Topic=topic_names,
                                     Content=Ycontent))
   rownames(beta_meta) <- rownames(beta)
   
   file_name_beta <- paste0(str_match(file_name,'(.*?).rds')[2],'_beta_',type_name,'.biom') 
   
   beta_biom <- make_biom(t(beta),beta_meta,beta_otu)
   write_biom(beta_biom,file.path(save_folder,file_name_beta))
   
}