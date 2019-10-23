#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 3) stop("Three arguments must be supplied (input file).n", call.=FALSE)

library(stm,quietly=TRUE,verbose=FALSE,warn.conflicts=FALSE)
library(biom,quietly=TRUE,verbose=FALSE,warn.conflicts=FALSE)
library(readr,quietly=TRUE,verbose=FALSE,warn.conflicts=FALSE)
library(tidyr,quietly=TRUE,verbose=FALSE,warn.conflicts=FALSE)
library(dplyr,quietly=TRUE,verbose=FALSE,warn.conflicts=FALSE)
library(fastICA,quietly=TRUE,verbose=FALSE,warn.conflicts=FALSE)
library(randomForest,quietly=TRUE,verbose=FALSE,warn.conflicts=FALSE)
library(stringr,quietly=TRUE,verbose=FALSE,warn.conflicts=FALSE)
library(kernlab,quietly=TRUE,verbose=FALSE,warn.conflicts=FALSE)
library(Rcpp,quietly=TRUE,verbose=FALSE,warn.conflicts=FALSE)
library(parallel,quietly=TRUE,verbose=FALSE,warn.conflicts=FALSE)
library(foreach,quietly=TRUE,verbose=FALSE,warn.conflicts=FALSE)
library(ape,quietly=TRUE,verbose=FALSE,warn.conflicts=FALSE)
library(phyloseq,quietly=TRUE,verbose=FALSE,warn.conflicts=FALSE)
library(doParallel,quietly=TRUE,verbose=FALSE,warn.conflicts=FALSE)
library(stm,quietly=TRUE,verbose=FALSE,warn.conflicts=FALSE)
library(LDAvis,quietly=TRUE,verbose=FALSE,warn.conflicts=FALSE)
library(caret,quietly=TRUE,verbose=FALSE,warn.conflicts=FALSE)
library(glmnet,quietly=TRUE,verbose=FALSE,warn.conflicts=FALSE)
library(ggplot2,quietly=TRUE,verbose=FALSE,warn.conflicts=FALSE)
library(knitr,quietly=TRUE,verbose=FALSE,warn.conflicts=FALSE)
library(gridExtra,quietly=TRUE,verbose=FALSE,warn.conflicts=FALSE)


source('~/Dropbox/stm_microbiome/code_active/stm_functions.R')
source('~/Dropbox/stm_microbiome/code_active/nav_froz_fxns_3.R') 
source('~/Dropbox/stm_microbiome/code_active/performance_1.R')
source('~/Dropbox/stm_microbiome/code_active/framework.R')
source('~/Dropbox/stm_microbiome/code_active/stm_ESTEP.R')

K <- as.integer(args[1])
rare_min <- as.integer(args[2])
rare_method <- as.character(args[3])


PS <- readRDS(sprintf('~/Dropbox/stm_microbiome/Children/Processing/01-raw/otus/ag_%s_%s_ps.rds',rare_method,rare_min)) ### 1000 read_min

meta <- data.frame(sample_data(PS))
otu <- as.matrix(t(otu_table(PS))@.Data)
taxa <- as.matrix(tax_table(PS)@.Data)


read_min <- NA
samp_min <- NA
if (rare_method != 'deseq'){
  read_min <- 2
  samp_min <- .99
  taxa_filter <- colSums(otu < read_min) > floor(samp_min * NROW(otu)) 
  
  cat('Removing ', sum(taxa_filter),  ' taxa (from ',NCOL(otu),' taxa) that have fewer than ',
      read_min, ' reads in less than ',(1-samp_min)*100,'% of samples.\n',sep='')
  
  otu <- otu[,!taxa_filter]
}


taxa <- taxa[colnames(otu),]
meta <- filter(meta,SampleID %in% rownames(otu)) %>% as.data.frame()
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
    'Total meta:',nrow(meta),'x',ncol(meta),'\n',
    'Total diagnosis:',names(table(meta$DIAGNOSIS)),'=',c(table(meta$DIAGNOSIS)),'\n')




seed_heldout <- 56 #3653
seed_fit <- 12 #6457
seed_rare <- 5346



heldout_N <- floor(.4*length(docs)) # was .1
heldout_p <- .4 # was .5 then .8
heldout <- make.heldout(docs,vocab,N=heldout_N,proportion=heldout_p,seed=seed_heldout)



meta$DIAGNOSIS <- ifelse(meta$DIAGNOSIS == 'CD',1,0)

fit0 <- stm(documents=heldout$documents,vocab=heldout$vocab,K=K,
            data=meta,
            max.em.its=500,init.type='Random',
            seed=seed_fit, #645647
            verbose=FALSE,reportevery=25)

fit1 <- stm(prevalence=~DIAGNOSIS,
            documents=heldout$documents,vocab=heldout$vocab,K=K,
            data=meta,
            max.em.its=500,init.type='Random',
            seed=seed_fit, #645647
            verbose=FALSE,reportevery=25)

fit2 <- stm(prevalence=~DIAGNOSIS,
            content=~DIAGNOSIS,
            documents=heldout$documents,vocab=heldout$vocab,K=K,
            data=meta,
            max.em.its=500,init.type='Random',
            control=list(kappa.enet=.5),
            seed=seed_fit, #645647
            verbose=FALSE,reportevery=25)

fit3 <- stm(prevalence=~PCDAI,
            content=~DIAGNOSIS,
            documents=heldout$documents,vocab=heldout$vocab,K=K,
            data=meta,
            max.em.its=500,init.type='Random',
            control=list(kappa.enet=.5),
            seed=seed_fit, #645647
            verbose=FALSE,reportevery=25)

fit4 <- stm(prevalence=~s(PCDAI),
            content=~DIAGNOSIS,
            documents=heldout$documents,vocab=heldout$vocab,K=K,
            data=meta,
            max.em.its=500,init.type='Random',
            control=list(kappa.enet=.5),
            seed=seed_fit, #645647
            verbose=FALSE,reportevery=25)



ls_var <- ls()

ls_fit <- ls_var[grepl('fit[0-9]$',ls_var)]
fits <- lapply(ls_fit, function(x) eval(parse(text=x)))
names(fits) <- ls_fit

evals <- lapply(fits, function(x) eval.heldout(x,heldout$missing))

dir_name <- file.path(sprintf('~/Dropbox/stm_microbiome/data_active/Gevers_ti/stm_s97_%s_%s_supervised',rare_method,rare_min),'heldout')
dir.create(dir_name,showWarnings=FALSE,recursive=TRUE)
fit_filename <- paste0('heldout_K_',K,'_nfits_',length(ls_fit),'.rds')
saveRDS(list(fits=fits,
             evals=evals,
             heldout=heldout,
             docs=docs,
             vocab=vocab,
             meta=meta,
             taxa=taxa,
             counts=counts,
             seeds=list(seed_fit=seed_fit,seed_rare=seed_rare,seed_heldout=seed_heldout),
             filters=list(read_min=read_min,samp_min=samp_min,rare_min=rare_min,heldout_N=heldout_N,heldout_p=heldout_p)),
        file.path(dir_name,fit_filename))


heldout <- readRDS(file.path(dir_name,fit_filename))
print(sapply(heldout$evals,function(x) x$expected.heldout))
