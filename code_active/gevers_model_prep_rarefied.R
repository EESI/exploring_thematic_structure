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
library(DESeq2)

source('~/Dropbox/stm_microbiome/code_active/stm_functions.R')
source('~/Dropbox/stm_microbiome/code_active/nav_froz_fxns_3.R') 
source('~/Dropbox/stm_microbiome/code_active/performance_1.R')
source('~/Dropbox/stm_microbiome/code_active/framework.R')

rare_method <- 'deseq' #rare #norm
seed_rare <- 5346
rare_min <- 1000 #quantile(colSums(OTU),.2)


biom_file <- read_biom('~/Dropbox/stm_microbiome/Children/Processing/01-raw/simcomp_otus_q32_sim97/sequences_filtered/otu_table_cnNormalized.biom')
data_mat_xavier <- floor(as.matrix(biom_data(biom_file)))
biom_file <- read_biom('~/Dropbox/stm_microbiome/Children/Processing/01-raw/simcomp_otus_q32_sim97/sequences_filtered/otu_table.biom')
otu_taxa_xavier <- observation_metadata(biom_file)

metadata_xavier <- read_delim('~/Dropbox/stm_microbiome/Children/Processing/01-raw/metadata.txt',delim='\t')
sample_info_xavier <- read_csv('~/Dropbox/stm_microbiome/Children/xavier_sample_info.csv')

metadata_xavier <- metadata_xavier %>%
  dplyr::select(SampleID = ends_with('SampleID'), everything()) %>%
  mutate(SAMPLE_NAME = ifelse(STRAIN == "missing", SAMPLE_NAME, STRAIN),
         SAMPLE_NAME = str_replace(SAMPLE_NAME,'-','')) %>%
  dplyr::select(SampleID,SAMPLE_NAME,ISOLATION_SOURCE)

metadata_xavier$COHORT <- ifelse(metadata_xavier$SAMPLE_NAME %in% sample_info_xavier$sample, "RISK","Other")
metadata_xavier <- metadata_xavier[metadata_xavier$SampleID %in% colnames(data_mat_xavier),] # make sure data_mat_xavier and metadata_xavier have same SRRs
metadata_xavier <- data.frame(metadata_xavier,
                              sample_info_xavier[match(metadata_xavier$SAMPLE_NAME,sample_info_xavier$sample),]) # combine metadata_xavier and sample_data

metadata_xavier <- metadata_xavier %>%
  mutate(ISOLATION_SOURCE = ifelse(ISOLATION_SOURCE == "missing", sample_location, ISOLATION_SOURCE),
         ISOLATION_SOURCE = str_to_title(ISOLATION_SOURCE),
         RACE = str_to_title(race),
         STUDY = 'Xavier',
         AGE = 'Minor',
         DIAGNOSIS = ifelse(COHORT == "Other", "CD", Diagnosis)) %>%
  dplyr::select(SampleID, SAMPLE_NAME, ISOLATION_SOURCE, DIAGNOSIS, RACE, STUDY, AGE, COHORT, PCDAI, AB_exposure)
data_mat_xavier <- data_mat_xavier[,metadata_xavier$SampleID]
metadata_risk <- metadata_xavier %>% dplyr::filter(COHORT == 'RISK',
                                                   ISOLATION_SOURCE %in% c('Terminalileum Biopsy'))
metadata_risk <- metadata_risk %>% dplyr::filter(AB_exposure == 'NoAB')

metadata_risk <- metadata_risk %>%
  mutate(PCDAI = ifelse(is.na(PCDAI) & DIAGNOSIS == 'Not IBD',0,PCDAI)) %>%
  filter(!is.na(PCDAI)) %>%
  mutate(PCDAI_RANK = ranker(PCDAI))

data_mat_risk <- data_mat_xavier[,metadata_risk$SampleID]
data_mat_risk <- data_mat_risk[rowSums(data_mat_risk) != 0,]


META <- sample_data(metadata_risk)
META$BURDEN <- metadata_risk %>%    
  mutate(burden=ifelse(PCDAI>0,1,0)) %>%
  group_by(burden) %>%
  mutate(burden_bin=ntile(PCDAI,3)) %>%
  ungroup() %>%
  mutate(burden_bin=as.factor(ifelse(DIAGNOSIS == 'CD',burden_bin,0))) %>%
  dplyr::select(burden_bin) %>%
  unlist()
rownames(META) <- META$SampleID



if (rare_method == 'norm'){
# NORMALIZATION
  
  data_mat_risk <- data_mat_risk[,colSums(data_mat_risk) > 0]
  data_mat_risk <- round(rare_min*t(t(data_mat_risk)/colSums(data_mat_risk)))
  data_mat_risk <- data_mat_risk[,colSums(data_mat_risk) > 0]
  OTU <- otu_table(data_mat_risk,taxa_are_rows=TRUE)
  
  otu_taxa_xavier <- as.matrix(otu_taxa_xavier[rownames(OTU),])
  TAXA <- tax_table(otu_taxa_xavier)
  
  PS <- phyloseq(OTU,META,TAXA)
  
}else if (rare_method == 'rare'){
# RAREFACTION
  
  OTU <- otu_table(data_mat_risk,taxa_are_rows=TRUE)
  OTU <- rarefy_even_depth(OTU,
                           rngseed=seed_rare, #6546
                           sample.size=rare_min,
                           replace=FALSE,
                           trimOTUs=TRUE,
                           verbose=TRUE)
  
  otu_taxa_xavier <- as.matrix(otu_taxa_xavier[rownames(OTU),])
  TAXA <- tax_table(otu_taxa_xavier)
  
  PS <- phyloseq(OTU,META,TAXA)
  
}else if (rare_method == 'deseq'){
# DESEQ2 VARIANCE STABLIZING NORMALIZATION

  data_mat_risk <- data_mat_risk[,colSums(data_mat_risk) > 0]
  META <- META[colnames(data_mat_risk),]
  OTU <- otu_table(data_mat_risk,taxa_are_rows=TRUE)
  
  otu_taxa_xavier <- as.matrix(otu_taxa_xavier[rownames(OTU),])
  TAXA <- tax_table(otu_taxa_xavier)
  
  PS <- phyloseq(OTU,META,TAXA)

  habdds <- phyloseq_to_deseq2(PS,~DIAGNOSIS)
  habdds_temp <- try(estimateSizeFactors(habdds), silent=TRUE)
  if (class(habdds_temp) == 'try-error') {
    geo_means <- apply(counts(habdds), 1, gm_mean)
    habdds <- estimateSizeFactors(habdds,geoMeans=geo_means)
    rm(habdds_temp)
  }else{    
    habdds <- estimateSizeFactors(habdds)
  }
  
  habdds <- estimateDispersions(habdds)
  diagvst <- getVarianceStabilizedData(habdds)
  
  diagvst[diagvst < 0] <- 0
  diagvst <- round(diagvst)
  diagvst <- diagvst[rowSums(diagvst) > 0,]
  otu_table(PS) <- otu_table(diagvst, taxa_are_rows = TRUE)
}


saveRDS(PS,sprintf('~/Dropbox/stm_microbiome/Children/Processing/01-raw/otus/ag_%s_%s_ps.rds',rare_method,rare_min))
