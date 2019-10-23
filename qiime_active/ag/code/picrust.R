#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)
if (length(args) == 0) stop('Please specify method.', call.=FALSE)

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
library(DESeq2)
library(tibble)
library(phyloseq)

source('~/Dropbox/stm_microbiome/code_active/stm_functions.R')
source('~/Dropbox/stm_microbiome/code_active/nav_froz_fxns_3.R') 
source('~/Dropbox/stm_microbiome/code_active/performance_1.R')
source('~/Dropbox/stm_microbiome/code_active/framework.R')

data_dir <- '~/Dropbox/stm_microbiome/qiime_active/ag/picked_otus'

if  (length(args) == 2){
  system2('predict_metagenomes.py',
          args=c('-i',file.path(data_dir,'otu_table_cnnorm.biom'),
                 '-o',file.path(data_dir,'ko_table.biom'),
                 '-t','ko'))
  
  system2('predict_metagenomes.py',
          args=c('-i',file.path(data_dir,'otu_table_cnnorm.biom'),
                 '-o',file.path(data_dir,'cog_table.biom'),
                 '-t','cog'))
}else{
  if (args[1] == 'kegg'){
    system2('predict_metagenomes.py',
            args=c('-i',file.path(data_dir,'otu_table_cnnorm.biom'),
                   '-o',file.path(data_dir,'ko_table.biom'),
                   '--with_confidence',
                   '-t','ko'))
  }
  if (args[1] == 'cog'){
    system2('predict_metagenomes.py',
            args=c('-i',file.path(data_dir,'otu_table_cnnorm.biom'),
                   '-o',file.path(data_dir,'cog_table.biom'),
                   '--with_confidence',
                   '-t','cog'))
  }
}