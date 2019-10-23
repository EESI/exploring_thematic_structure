#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 1) stop('Please specify number of topics, K.', call.=FALSE)

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

data_dir <- '~/Dropbox/stm_microbiome/qiime_active/ag'
models_dir <- file.path(data_dir,'models',K)
dir.create(models_dir,showWarnings=FALSE,recursive=TRUE)


seed_fit <- 12



PS <- readRDS(file.path(data_dir,'data_clean.rds'))

OTU <- as(otu_table(PS),'matrix')
TAX <- as(tax_table(PS),'matrix')
META <- as(sample_data(PS),'data.frame')

taxon_ids <- data.frame(long=colnames(OTU),
                        ids=paste0('otu',1:NCOL(OTU)),
                        row.names=paste0('otu',1:NCOL(OTU)))

counts <- list(META_clean=META,
               table_clean=OTU,
               ids=taxon_ids)
colnames(counts$table_clean) <- as.character(counts$ids$ids)

if (!identical(rownames(counts$table_clean),counts$META_clean$SampleID)) stop('Please make table sample names match METAdata sample IDs!\n')

VOCAB <- as.character(counts$ids$ids)
DOCS <- lapply(1:NROW(counts$table_clean), function(i) reshape_doc(counts$table_clean[i,],VOCAB)) 
META <- counts$META_clean



cat('Total docs:',length(DOCS),'\n',
    'Total vocab:',length(VOCAB),'\n',
    'Total meta:',nrow(META),'x',ncol(META),'\n',
    'Total diet:',names(table(META$diet)),'=',c(table(META$diet)),'\n')



META$diet <- ifelse(META$diet == 'V',1,0)

fit0 <- stm(documents=DOCS,vocab=VOCAB,K=K,
            data=META,
            max.em.its=500,init.type='Spectral',
            seed=seed_fit,
            verbose=TRUE,reportevery=25)

fit1 <- stm(prevalence=~diet,
            documents=DOCS,vocab=VOCAB,K=K,
            data=META,
            max.em.its=500,init.type='Spectral',
            seed=seed_fit, 
            verbose=TRUE,reportevery=25)


ls_var <- ls()
ls_fit <- ls_var[grepl('fit[0-9]$',ls_var)]
FITS <- lapply(ls_fit, function(x) eval(parse(text=x)))
names(FITS) <- ls_fit


dat <- list(fits=FITS,
            docs=DOCS,
            vocab=VOCAB,
            meta=META,
            taxa=TAX,
            counts=counts,
            seed=seed_fit)

BETA <- vector(mode='list',length=length(ls_fit))
names(BETA) <- ls_fit
for (f in seq_along(dat$fits)){
  
  fit <- dat$fits[[f]]
  K <- fit$settings$dim$K
  
  # scaling to 10000 to counter lower probability OTUs within topic
  beta <- t(round(10000*exp(do.call('rbind',fit$beta$logbeta))))
  colnames(beta) <- paste0('T',1:NCOL(beta))
  rownames(beta) <- as.character(dat$counts$ids[dat$vocab,'long'])
  
  BETA[[f]] <- beta
  
  beta_meta <- data.frame(Topic=colnames(beta))
  beta_taxa <- dat$taxa
  
  beta_biom <- make_biom(beta,beta_meta,beta_taxa)
  
  beta_tmp <- tempfile()
  write_biom(beta_biom,beta_tmp)
  
  system2('predict_metagenomes.py',
          args=c('-i',beta_tmp,
                 '--type_of_prediction','ko',
                 '--with_confidence',
                 '-o',file.path(models_dir,sprintf('%s_k_%s_%s.biom','ko',K,names(dat$fits)[f]))))
  system2('predict_metagenomes.py',
          args=c('-i',beta_tmp,
                 '--type_of_prediction','cog',
                 '--with_confidence',
                 '-o',file.path(models_dir,sprintf('%s_k_%s_%s.biom','cog',K,names(dat$fits)[f]))))
  
}

dat$beta <- BETA

fit_fn <- sprintf('stm_k_%s.rds',K)
saveRDS(dat,file.path(models_dir,fit_fn))
