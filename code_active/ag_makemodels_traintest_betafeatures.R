#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)
if (length(args) == 0) stop("At least one argument must be supplied (input file).n", call.=FALSE)

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

K <- as.integer(args[1])
rare_min <- as.integer(args[2])
model <- as.integer(args[3])
seed_split <- as.integer(args[4]) #234

seed_fit <- 6457
seed_rare <- 5346



PS <- readRDS(sprintf('~/Dropbox/stm_microbiome/AG/Processing/01-raw/otus/ag_rare_%s_ps.rds',rare_min)) ### 1000 read_min

meta <- data.frame(sample_data(PS))
otu <- as.matrix(t(otu_table(PS))@.Data)
taxa <- as.matrix(tax_table(PS)@.Data)


meta <- meta %>% filter(BODY_HABITAT %in% 'UBERON:feces', !is.na(AGE_DECADE_RANK), !is.na(SEX))
otu <- otu[meta$SampleID,]


read_min <- 2
samp_min <- .99
taxa_filter <- colSums(otu < read_min) > floor(samp_min * NROW(otu)) 

cat('Removing ', sum(taxa_filter),  ' taxa (from ',NCOL(otu),' taxa) that have fewer than ',
    read_min, ' reads in less than ',(1-samp_min)*100,'% of samples.\n',sep='')

otu <- otu[,!taxa_filter]
taxa <- taxa[colnames(otu),]
meta <- filter(meta,SampleID %in% rownames(otu)) %>% data.frame()
rownames(meta) <- meta$SampleID


taxon_ids <- data.frame(long=colnames(otu),
                        ids=paste0('otu',1:NCOL(otu)),
                        row.names=paste0('otu',1:NCOL(otu)))

counts <- list(meta_clean=meta,
               table_clean=otu,
               ids=taxon_ids)
colnames(counts$table_clean) <- as.character(counts$ids$ids)

if (!identical(rownames(counts$table_clean),counts$meta_clean$SampleID)) stop('Please make table sample names match metadata sample IDs!\n')

vocab <- as.character(counts$ids$ids)
docs <- lapply(1:NROW(counts$table_clean),function(i) reshape_doc(counts$table_clean[i,],vocab)) ### something might be wrong!!!!
meta <- counts$meta_clean



cat('Total docs:',length(docs),'\n',
    'Total vocab:',length(vocab),'\n',
    'Total Meta:',nrow(meta),'x',ncol(meta),'\n',
    'Total Sex:',names(table(meta$SEX)),'=',c(table(meta$SEX)),'\n')



meta$FEMALE <- ifelse(meta$SEX == 'female',1,0)
meta$AGE_YEARS <- as.numeric(meta$AGE_YEARS)
meta$AGE <- as.vector(scale(as.numeric(meta$AGE_YEARS)))


set.seed(seed_split)
idx_train <- caret::createDataPartition(y=with(meta,as.factor(SEX):as.factor(AGE_CAT)),times=seed_split,p=.8,list=TRUE)[[seed_split]]
otu_train <- otu[idx_train,]
meta_train <- meta[idx_train,]
docs_train <- docs[idx_train]
otu_test <- otu[-idx_train,]
meta_test <- meta[-idx_train,]
docs_test <- docs[-idx_train]



if (model == 1){
fit <- stm(documents=docs_train,vocab=vocab,K=K,
            data=meta_train,
            max.em.its=500,init.type='Spectral',
            seed=seed_fit, #645647
            verbose=TRUE,reportevery=25)
}

if (model == 2){
fit <- stm(prevalence=~FEMALE + s(AGE),
            content=~FEMALE,
            documents=docs_train,vocab=vocab,K=K,
            data=meta_train,
            max.em.its=500,init.type='Spectral',
            control=list(kappa.enet=.5),
            seed=seed_fit, #645647
            verbose=TRUE,reportevery=25)
}



dir_name <- sprintf('~/Dropbox/stm_microbiome/data_active/AG_female/stm_s97_rarefied_%s_supervised/features/seed_%s',rare_min,seed_split)
dir.create(dir_name,showWarnings=FALSE,recursive=TRUE)
fit_filename <- paste0('fits_K_',K,'_model_',model,'.rds')
saveRDS(list(fits=list(fit=fit),
             docs=list(docs_train=docs_train,docs_test=docs_test),
             vocab=vocab,
             meta=list(meta_train=meta_train,meta_test=meta_test),
             otu=list(otu_train=otu_train,otu_test=otu_test),
             taxa=taxa,
             counts=counts,
             seeds=list(seed_fit=seed_fit,seed_rare=seed_rare,seed_split=seed_split),
             filters=list(read_min=read_min,samp_min=samp_min,rare_min=rare_min)),
        file.path(dir_name,fit_filename))





# 
# 
# make_beta_biom <- function(dat,model,dir_name,override=FALSE){
#   require(biom)
#   
#   fit <- dat$fits$fit
#   K <- fit$settings$dim$K
#   Nmin <- dat$filters$rare_min
#   
#   beta <- t(floor(exp(do.call('rbind',fit$beta$logbeta))*Nmin))
#   colnames(beta) <- paste0('T',1:NCOL(beta))
#   rownames(beta) <- as.character(dat$counts$ids[dat$vocab,'long'])
#   
#   beta_meta <- data.frame(Topic=colnames(beta))
#   beta_taxa <- dat$taxa
#   
#   beta_filename <- paste0('beta_K_',K,'_model_',model,'.biom')
#   beta_biom <- make_biom(beta,beta_meta,beta_taxa)
#   
#   if (!override) if (file.exists(file.path(dir_name,beta_filename))) stop('File already exists!')
#   write_biom(beta_biom,file.path(dir_name,beta_filename))
# }
#  
# K <- 35
# cov <- 'female'
# rare_min <- 10000
# 
# for (model in 1:2){
#   for (seed_split in c(1,105,12,53,74)){
#     dir_name <- sprintf('~/Dropbox/stm_microbiome/data_active/AG_%s/stm_s97_rarefied_%s_supervised/features/seed_%s',cov,rare_min,seed_split)
#     fit_filename <- paste0('fits_K_',K,'_model_',model,'.rds')
#     
#     dat <- readRDS(file.path(dir_name,fit_filename))
#     
#     sink(file.path(dir_name,'fit_info.txt'))
#     make_beta_biom(dat,model,dir_name,override=TRUE)
#     sink()
#   }
# }