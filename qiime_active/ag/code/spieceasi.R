#!/usr/bin/env Rscript

library(SpiecEasi)
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


data_dir <- '~/Dropbox/stm_microbiome/qiime_active/ag'

seed_fit <- 12


PS <- readRDS(file.path(data_dir,'data_clean.rds'))

SE <- spiec.easi(PS, method='mb', 
                     lambda.min.ratio=1e-2,nlambda=20,
                     icov.select.params=list(rep.num=50,ncores=40), 
                     verbose=TRUE)

SE_graph <- adj2igraph(SE$refit,vertex.attr=list(name=taxa_names(PS)))

saveRDS(list(PS=PS,glasso=SE,glass_graph=SE_graph),
        file.path(data_dir,'spieceasi.rds'))
