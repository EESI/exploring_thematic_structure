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
library(viridis)
library(DESeq2)

gm_mean <- function(x, na.rm=TRUE) exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))

source('~/Dropbox/stm_microbiome/code_active/stm_functions.R')
source('~/Dropbox/stm_microbiome/code_active/nav_froz_fxns_3.R')


data_dir <- '~/Dropbox/stm_microbiome/qiime_active/ag'

META <- read_delim(file.path(data_dir,'PRJEB11419.txt'),'\t') %>%
  select(PRIMARY_ID=secondary_sample_accession,run_accession) %>%
  left_join(readRDS(file.path(data_dir,'metadata.rds')),by='PRIMARY_ID') %>%
  dplyr::rename(SampleID=run_accession) %>%
  filter(body_site %in% c('UBERON:feces'),
         !(age_cat %in% c('baby')),
         diet_type %in% c('Omnivore','Vegan','Vegetarian')) %>%
  mutate(diet=ifelse(diet_type == 'Omnivore','O','V'))

META <- as.data.frame(META)
rownames(META) <- META$SampleID

ko_fn <- file.path(data_dir,'picked_otus/ko_table.biom')

biom_file <- read_biom(ko_fn)
ko <- as.matrix(biom_data(biom_file))
ko <- ko[rowSums(ko)>0,]

SAMPLEIDS <- intersect(colnames(ko),META$SampleID)

META <- META[SAMPLEIDS,]
ko <- ko[,SAMPLEIDS]

PS <- phyloseq(otu_table(ko,taxa_are_rows=TRUE),sample_data(META))

DS2 <- phyloseq_to_deseq2(PS, ~ diet)

DS2_OUT <- try(DESeq(DS2,test='Wald',fitType='parametric'), silent=TRUE)
if (class(DS2_OUT) == 'try-error') {
  geo_means <- apply(counts(DS2),1,gm_mean)
  DS2_OUT <- estimateSizeFactors(DS2,geoMeans=geo_means)
  DS2_OUT <- DESeq(DS2_OUT,fitType='local')
}

res <- DESeq2::results(DS2_OUT, cooksCutoff=FALSE, pAdjustMethod='BH', contrast=c('diet','V','O'))

saveRDS(list(PS=PS,DS2=DS2_OUT,res=res),file.path(data_dir,'deseq2_ko.rds'))