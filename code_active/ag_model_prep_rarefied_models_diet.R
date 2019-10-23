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


source('~/Dropbox/stm_microbiome/code_active/stm_functions.R')
source('~/Dropbox/stm_microbiome/code_active/nav_froz_fxns_3.R') 
source('~/Dropbox/stm_microbiome/code_active/performance_1.R')
source('~/Dropbox/stm_microbiome/code_active/framework.R')

PS <- readRDS('~/Dropbox/stm_microbiome/AG/Processing/01-raw/otus/ag_rare_ps.rds')


meta <- data.frame(sample_data(PS))
otu <- as.matrix(t(otu_table(PS))@.Data)
taxa <- as.matrix(tax_table(PS)@.Data)


meta <- meta %>% filter(BODY_HABITAT %in% 'UBERON:feces', BMI_CAT %in% c('Obese', 'Normal')) 
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




seed_fit <- 6457
seed_rare <- 5346
rare_min <- 1000

meta$OBESE <- ifelse(meta$BMI_CAT == 'Obese',1,0)
meta$DIET <- ifelse(meta$DIET_TYPE == 'V',1,0)
K <- 50
fit0 <- stm(documents=docs,vocab=vocab,K=K,
            data=meta,
            max.em.its=500,init.type='Spectral',
            seed=seed_fit, #645647
            verbose=TRUE,reportevery=25)

fit1 <- stm(prevalence=~DIET,
            documents=docs,vocab=vocab,K=K,
            data=meta,
            max.em.its=500,init.type='Spectral',
            seed=seed_fit, #645647
            verbose=TRUE,reportevery=25)

fit2 <- stm(prevalence=~OBESE,
            documents=docs,vocab=vocab,K=K,
            data=meta,
            max.em.its=500,init.type='Spectral',
            seed=seed_fit, #645647
            verbose=TRUE,reportevery=25)

fit3 <- stm(prevalence=~DIET + OBESE,
            documents=docs,vocab=vocab,K=K,
            data=meta,
            max.em.its=500,init.type='Spectral',
            seed=seed_fit, #645647
            verbose=TRUE,reportevery=25)

fit4 <- stm(prevalence=~DIET * OBESE,
            documents=docs,vocab=vocab,K=K,
            data=meta,
            max.em.its=500,init.type='Spectral',
            seed=seed_fit, #645647
            verbose=TRUE,reportevery=25)

fit5 <- stm(prevalence=~DIET,
            content=~OBESE,
            documents=docs,vocab=vocab,K=K,
            data=meta,
            max.em.its=500,init.type='Spectral',
            seed=seed_fit, #645647
            verbose=TRUE,reportevery=25)

fit6 <- stm(prevalence=~DIET + OBESE,
            content=~OBESE,
            documents=docs,vocab=vocab,K=K,
            data=meta,
            max.em.its=500,init.type='Spectral',
            seed=seed_fit, #645647
            verbose=TRUE,reportevery=25)

fit7 <- stm(prevalence=~DIET * OBESE,
            content=~OBESE, 
            documents=docs,vocab=vocab,K=K,
            data=meta,
            max.em.its=500,init.type='Spectral',
            seed=seed_fit, #645647
            verbose=TRUE,reportevery=25)


dir_name <- paste0('~/Dropbox/stm_microbiome/data_active/AG_diet/stm_s97_rarefied_supervised')
dir.create(dir_name,showWarnings=FALSE,recursive=TRUE)
fit_filename <- paste0('fits_K_',K,'.rds')
saveRDS(list(fits=list(fit1,fit2,fit3,fit4,fit5,fit6,fit7,fit0),
             docs=docs,
             vocab=vocab,
             meta=meta,
             taxa=taxa,
             counts=counts,
             seeds=list(seed_fit=seed_fit,seed_rare=seed_rare),
             filters=list(read_min=read_min,samp_min=samp_min,rare_min=rare_min)),
        file.path(dir_name,fit_filename))