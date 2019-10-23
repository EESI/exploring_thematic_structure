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

seed_rare <- 5346


dat <- readRDS('~/Dropbox/stm_microbiome/AG/Processing/01-raw/otus/ag_data_cnNorm.rds')
data_mat_ag <- dat$data
otu_taxa_ag <- dat$taxa
rm(dat)


metadata_ag <- read_delim('~/Dropbox/stm_microbiome/AG/Processing/01-raw/metadata.txt',delim='\t') %>%
   dplyr::select(SampleID = ends_with('SampleID'), 
                 LIVER_DISEASE,
                 BODY_PRODUCT,
                 AGE_CAT,
                 AGE_YEARS,
                 SIBO,
                 DOG,
                 ANONYMIZED_NAME,
                 IBD,
                 CAT,
                 IBS,
                 ACNE_MEDICATION_OTC,
                 DIABETES,
                 CDIFF,
                 BODY_SITE,
                 FUNGAL_OVERGROWTH,
                 FLOSSING_FREQUENCY,
                 NON_FOOD_ALLERGIES_SUN,
                 CARDIOVASCULAR_DISEASE,
                 FED_AS_INFANT,
                 ACID_REFLUX,
                 BMI,
                 TEETHBRUSHING_FREQUENCY,
                 DEPRESSION_BIPOLAR_SCHIZOPHRENIA,
                 ROOMMATES,
                 WEIGHT_CHANGE,
                 SEX,
                 ASD,
                 LUNG_DISEASE,
                 ALCOHOL_FREQUENCY,
                 SMOKING_FREQUENCY,
                 KIDNEY_DISEASE,
                 WEIGHT_KG,
                 HEIGHT_CM,
                 PROBIOTIC_FREQUENCY,
                 SLEEP_DURATION,
                 BOWEL_MOVEMENT_QUALITY,
                 MENTAL_ILLNESS_TYPE_SCHIZOPHRENIA,
                 THYROID,
                 NON_FOOD_ALLERGIES_DRUG_EG_PENICILLIN,
                 ACNE_MEDICATION,
                 ANTIBIOTIC_HISTORY,
                 BODY_HABITAT,
                 BMI_CORRECTED,
                 BMI_CAT,
                 DIET_TYPE,
                 AUTOIMMUNE) %>%
   filter(SampleID %in% colnames(data_mat_ag),
          BODY_HABITAT %in% c('UBERON:feces','UBERON:oral cavity','UBERON:skin'),
          SEX %in% c('male','female'),
          !(AGE_CAT %in% c('Unknown','baby')),
          DIET_TYPE %in% c('Vegan','Vegetarian','Omnivore')) %>%
   mutate(AGE_DECADE=cut(as.numeric(AGE_YEARS),c(seq(0,70,10),100)),
          AGE_DECADE_RANK=as.integer(as.factor(AGE_DECADE)),
          DIET_TYPE = ifelse(DIET_TYPE %in% 'Omnivore','O','V'))


data_mat_ag <- data_mat_ag[,metadata_ag$SampleID]
data_mat_ag <- data_mat_ag[rowSums(data_mat_ag) != 0,]

metadata_ag$TOTAL_READS <- colSums(data_mat_ag)


OTU <- otu_table(data_mat_ag,taxa_are_rows=TRUE)
rare_min <- 10000 # 1000
OTU <- rarefy_even_depth(OTU,
                         rngseed=seed_rare, #6546
                         sample.size=rare_min,
                         replace=FALSE,
                         trimOTUs=TRUE,
                         verbose=TRUE)

metadata_ag <- data.frame(metadata_ag)
rownames(metadata_ag) <- metadata_ag$SampleID
metadata_ag <- metadata_ag[colnames(OTU),]
META <- sample_data(metadata_ag)

otu_taxa_ag <- as.matrix(otu_taxa_ag[rownames(OTU),])
TAXA <- tax_table(otu_taxa_ag)

PS <- phyloseq(OTU,META,TAXA)


# ord <- ordinate(PS,method='NMDS',distance='bray')
# eig <- ord$values$Eigenvalues
# plot_ordination(PS,ord,color='SEX') 



saveRDS(PS,sprintf('~/Dropbox/stm_microbiome/AG/Processing/01-raw/otus/ag_rare_%s_ps.rds',rare_min))
