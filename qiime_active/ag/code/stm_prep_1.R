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
library(tibble)
library(phyloseq)

source('~/Dropbox/stm_microbiome/code_active/stm_functions.R')
source('~/Dropbox/stm_microbiome/code_active/nav_froz_fxns_3.R') 
source('~/Dropbox/stm_microbiome/code_active/performance_1.R')
source('~/Dropbox/stm_microbiome/code_active/framework.R')

data_dir <- '~/Dropbox/stm_microbiome/qiime_active/ag/'


OTU <- readRDS(file.path(data_dir,'seqtab_cnnorm.rds'))
TAX <- readRDS(file.path(data_dir,'tax.rds'))

META <- read_delim(file.path(data_dir,'PRJEB11419.txt'),'\t') %>%
  select(PRIMARY_ID=secondary_sample_accession,run_accession) %>%
  left_join(readRDS(file.path(data_dir,'metadata.rds')),by='PRIMARY_ID') %>%
  dplyr::rename(SampleID=run_accession) %>%
  filter(body_site %in% c('UBERON:feces'),
         !(age_cat %in% c('baby')),
         diet_type %in% c('Omnivore','Vegan','Vegetarian')) %>%
  mutate(diet=ifelse(diet_type == 'Omnivore','O','V'))


SAMPLEIDS <- intersect(names(which(rowSums(OTU)>=1000)),META$SampleID)

META <- META %>% filter(SampleID %in% SAMPLEIDS)
OTU <- OTU[SAMPLEIDS,]

OTUIDS <- intersect(rownames(TAX)[!is.na(TAX$Phylum)],names(which(colSums(OTU)>0)))
OTU <- OTU[,OTUIDS]


META <- as.data.frame(META)
rownames(META) <- META$SampleID
META <- sample_data(META[rownames(OTU),])
TAX <- tax_table(as.matrix(TAX[colnames(OTU),]))
OTU <- otu_table(OTU,taxa_are_rows=FALSE)
PS <- phyloseq(OTU,META,TAX)


table(sort(colSums(otu_table(PS)>0)))[1:10]
data.frame(otu=colnames(otu_table(PS)),
           abundance=colSums(otu_table(PS)),
           prevalence=colSums(otu_table(PS)>0)/nrow(otu_table(PS))) %>%
  left_join(data.frame(otu=rownames(tax_table(PS)),tax_table(PS)),by='otu') %>%
  filter(Phylum %in% names(sort(table(Phylum),decreasing=TRUE)[1:25])) %>%
  ggplot(aes(abundance,prevalence)) +
  geom_hline(yintercept=10/nrow(otu_table(PS)),color='red',linetype=2) +
  geom_point(alpha=.5) +
  facet_wrap(~Phylum) +
  scale_x_log10() +
  ggtitle('OTUs')


PS <- prune_taxa(colSums(otu_table(PS) > 0) >= 10,PS)

table(sort(colSums(otu_table(PS)>0)))[1:10]
data.frame(otu=colnames(otu_table(PS)),
           abundance=colSums(otu_table(PS)),
           prevalence=colSums(otu_table(PS)>0)/nrow(otu_table(PS))) %>%
  left_join(data.frame(otu=rownames(tax_table(PS)),tax_table(PS)),by='otu') %>%
  filter(Phylum %in% names(sort(table(Phylum),decreasing=TRUE)[1:25])) %>%
  ggplot(aes(abundance,prevalence)) +
  geom_point(alpha=.5) +
  facet_wrap(~Phylum) +
  scale_x_log10() +
  ggtitle('OTUs')

data.frame(sample=rownames(otu_table(PS)),
           abundance=rowSums(otu_table(PS)),
           prevalence=rowSums(otu_table(PS)>0)/ncol(otu_table(PS))) %>%
  ggplot(aes(abundance,prevalence)) +
  geom_point(alpha=.5) +
  scale_x_log10() +
  ggtitle('Samples')


saveRDS(PS,file.path(data_dir,'data_clean.rds'))
