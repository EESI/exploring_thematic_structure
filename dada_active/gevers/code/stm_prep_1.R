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

source('~/Dropbox/stm_microbiome/code_active/stm_functions.R')
source('~/Dropbox/stm_microbiome/code_active/nav_froz_fxns_3.R') 
source('~/Dropbox/stm_microbiome/code_active/performance_1.R')
source('~/Dropbox/stm_microbiome/code_active/framework.R')

data_dir <- '~/Dropbox/stm_microbiome/dada_active/gevers/'

OTU <- readRDS(file.path(data_dir,'seqtab.rds'))
TAX <- readRDS(file.path(data_dir,'tax.rds'))
TAX_boot <- TAX$boot
TAX <- as.data.frame(TAX$tax)

META <- read_delim('~/Dropbox/stm_microbiome/Children/Processing/01-raw/metadata.txt',delim='\t') %>%
  mutate(SAMPLE_NAME = ifelse(STRAIN == "missing", SAMPLE_NAME, STRAIN),
         SAMPLE_NAME = str_replace(SAMPLE_NAME,'-','')) %>%
  dplyr::select(SampleID=`#SampleID`,SAMPLE_NAME,ISOLATION_SOURCE) %>%
  left_join(read_csv('~/Dropbox/stm_microbiome/Children/xavier_sample_info.csv') %>% mutate(SAMPLE_NAME=sample),
            by='SAMPLE_NAME') %>%
  mutate(COHORT=ifelse(SAMPLE_NAME %in% sample,'RISK','Other'),
         ISOLATION_SOURCE=ifelse(ISOLATION_SOURCE %in% 'missing',sample_location,ISOLATION_SOURCE),
         ISOLATION_SOURCE=str_to_title(ISOLATION_SOURCE),
         ISOLATION_SOURCE=ifelse(ISOLATION_SOURCE %in% 'Terminal Ileum','Terminalileum Biopsy',ISOLATION_SOURCE),
         RACE=str_to_title(race),
         STUDY='Xavier',
         AGE='Minor',
         DIAGNOSIS=ifelse(COHORT %in% 'Other','CD',Diagnosis),
         PCDAI = ifelse(is.na(PCDAI) & DIAGNOSIS == 'Not IBD',0,PCDAI)) %>%
  dplyr::select(SampleID, SAMPLE_NAME, ISOLATION_SOURCE, DIAGNOSIS, RACE, STUDY, AGE, COHORT, PCDAI, AB_exposure)


META <- META %>%
  filter(ISOLATION_SOURCE %in% 'Terminalileum Biopsy',
         !(AB_exposure %in% 'AB'))

SAMPLEIDS <- intersect(intersect(names(which(rowSums(OTU) > 0)),META$SampleID),names(which(rowSums(OTU)>=1000)))

META <- META %>% filter(SampleID %in% SAMPLEIDS)
OTU <- OTU[SAMPLEIDS,]

OTUIDS <- intersect(names(which(!is.na(TAX$Phylum))),names(which(colSums(OTU)>0)))
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
  scale_x_log10()


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
  scale_x_log10()


saveRDS(PS,file.path(data_dir,'data_clean.rds'))
