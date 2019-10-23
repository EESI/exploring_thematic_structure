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



seed_fit <- 6457
seed_rare <- 5346


meta$FEMALE <- ifelse(meta$SEX == 'female',1,0)
meta$AGE_YEARS <- as.numeric(meta$AGE_YEARS)
meta$AGE <- as.vector(scale(as.numeric(meta$AGE_YEARS)))


fit0 <- stm(documents=docs,vocab=vocab,K=K,
            data=meta,
            max.em.its=500,init.type='Spectral',
            seed=seed_fit, #645647
            verbose=TRUE,reportevery=25)

fit1 <- stm(prevalence=~FEMALE,
            documents=docs,vocab=vocab,K=K,
            data=meta,
            max.em.its=500,init.type='Spectral',
            seed=seed_fit, #645647
            verbose=TRUE,reportevery=25)

fit2 <- stm(prevalence=~FEMALE,
            content=~FEMALE,
            documents=docs,vocab=vocab,K=K,
            data=meta,
            max.em.its=500,init.type='Spectral',
            control=list(kappa.enet=.5),
            seed=seed_fit, #645647
            verbose=TRUE,reportevery=25)

fit3 <- stm(prevalence=~FEMALE + s(AGE),
            content=~FEMALE,
            documents=docs,vocab=vocab,K=K,
            data=meta,
            max.em.its=500,init.type='Spectral',
            control=list(kappa.enet=.5),
            seed=seed_fit, #645647
            verbose=TRUE,reportevery=25)

# fit2 <- stm(prevalence=~FEMALE + AGE_DECADE_RANK,
#             documents=docs,vocab=vocab,K=K,
#             data=meta,
#             max.em.its=500,init.type='Spectral',
#             seed=seed_fit, #645647
#             verbose=TRUE,reportevery=25)
# 
# fit3 <- stm(prevalence=~AGE_DECADE_RANK,
#             content=~FEMALE,
#             documents=docs,vocab=vocab,K=K,
#             data=meta,
#             max.em.its=500,init.type='Spectral',
#             seed=seed_fit, #645647
#             verbose=TRUE,reportevery=25)
# 
# fit4 <- stm(prevalence=~FEMALE + AGE_DECADE_RANK,
#             content=~FEMALE,
#             documents=docs,vocab=vocab,K=K,
#             data=meta,
#             max.em.its=500,init.type='Spectral',
#             seed=seed_fit, #645647
#             verbose=TRUE,reportevery=25)

# fit5 <- stm(prevalence=~FEMALE * AGE_DECADE_RANK,
#             content=~FEMALE, 
#             documents=docs,vocab=vocab,K=K,
#             data=meta,
#             max.em.its=500,init.type='Spectral',
#             seed=seed_fit, #645647
#             verbose=TRUE,reportevery=25)


ls_var <- ls()
ls_var <- ls_var[grepl('fit[0-9]$',ls_var)]
fits <- lapply(ls_var, function(x) eval(parse(text=x)))
names(fits) <- ls_var

dir_name <- file.path(sprintf('~/Dropbox/stm_microbiome/data_active/AG_female/stm_s97_rarefied_%s_supervised',rare_min),'model')
dir.create(dir_name,showWarnings=FALSE,recursive=TRUE)
fit_filename <- paste0('fits_K_',K,'_nfits_',length(ls_var),'.rds')
saveRDS(list(fits=fits,
             docs=docs,
             vocab=vocab,
             meta=meta,
             taxa=taxa,
             counts=counts,
             seeds=list(seed_fit=seed_fit,seed_rare=seed_rare),
             filters=list(read_min=read_min,samp_min=samp_min,rare_min=rare_min)),
        file.path(dir_name,fit_filename))