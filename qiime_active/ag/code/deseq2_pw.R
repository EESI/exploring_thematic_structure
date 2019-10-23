#!/usr/bin/env Rscript

library(biomformat)
library(phyloseq)
library(vegan)
library(Matrix)
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(Rcpp)
library(parallel)
library(foreach)
library(phyloseq)
library(doParallel)
library(rstan)
library(lme4)
library(DESeq2)

source('~/Dropbox/stm_microbiome/code_active/stm_functions.R')
source('~/Dropbox/stm_microbiome/code_active/nav_froz_fxns_3.R') 
source('~/Dropbox/stm_microbiome/code_active/performance_1.R')
source('~/Dropbox/stm_microbiome/code_active/framework.R')

gm_mean <- function(x, na.rm=TRUE) exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))

data_dir <- '~/Dropbox/stm_microbiome/qiime_active/ag'

PS <- readRDS('~/Dropbox/stm_microbiome/qiime_active/ag/data_clean.rds')
SAMP <- sample_data(PS)

biom_file <- read_biom('/data/sw1/Dropbox/stm_microbiome/qiime_active/ag/picked_otus/ko_table.biom')
ko <- as.matrix(biom_data(biom_file))
samp_ids <- intersect(colnames(ko),rownames(SAMP))
y <- ifelse(SAMP[samp_ids,'diet']=='V','1','0')
ko <- round(ko * 10000/max(ko)) # scale to max=10k
ko <- ko[,samp_ids]
ko <- ko[rowSums(ko)>0,]
kegg_meta_2 <- kegg_metadata(biom_file,2)[rownames(ko)]
kegg_meta_3 <- kegg_metadata(biom_file,3)[rownames(ko)]

my_kegg_order <- c('carbohydrate metabolism',
                   'amino acid metabolism',
                   'lipid metabolism',
                   'cell motility',
                   'glycan biosynthesis and metabolism',
                   'xenobiotics biodegradation and metabolism',
                   'metabolism of cofactors and vitamins',
                   'nucleotide metabolism',
                   'metabolism of terpenoids and polyketides',
                   'biosynthesis of other secondary metabolites',
                   'energy metabolism',
                   'metabolism',
                   'membrane transport',
                   'cellular processes and signaling',
                   'genetic information processing')

targets <- sapply(kegg_meta_2, function(x) which(x$pathway %in% my_kegg_order))
targets <- targets[sapply(targets,length)>0]

kegg_meta <- vector(mode='list',length=length(targets))
names(kegg_meta) <- names(targets)
for (i in seq_along(targets)){
  target <- targets[i]
  target_name <- names(target)
  target_idx <- target[[1]]
  
  temp_kegg_meta <- kegg_meta_3[[target_name]]
  temp_kegg_meta$pathway <- temp_kegg_meta$pathway[target_idx]
  
  kegg_meta[[target_name]] <- temp_kegg_meta
}


ko_colsums <- colSums(ko)
target_pws <- unique(unlist(sapply(kegg_meta,function(x) x$pathway)))
target_kos <- lapply(target_pws,function(pw) names(which(sapply(kegg_meta, function(x) any(x$pathway %in% pw)))))
names(target_kos) <- target_pws

mat <- matrix(0,length(target_pws),length(samp_ids),dimnames=list(target_pws,samp_ids))
for (pw in target_pws){
  kos <- target_kos[[pw]]
  mat[pw,samp_ids] <- colSums(ko[kos,samp_ids,drop=FALSE])
}
mat <- t(mat)
diet <- ifelse(y=='1','V','O')
diet <- as.data.frame(diet[rownames(mat),,drop=FALSE])
diet$diet <- factor(diet$diet,levels=c('V','O'))

PS <- phyloseq(otu_table(mat,taxa_are_rows=FALSE),sample_data(diet))
DS2 <- phyloseq_to_deseq2(PS, ~ diet)

DS2_OUT <- try(DESeq(DS2,test='Wald',fitType='parametric'), silent=TRUE)
if (class(DS2_OUT) == 'try-error') {
  geo_means <- apply(counts(DS2),1,gm_mean)
  DS2_OUT <- estimateSizeFactors(DS2,geoMeans=geo_means)
  DS2_OUT <- DESeq(DS2_OUT,fitType='local')
}

res <- DESeq2::results(DS2_OUT, cooksCutoff=FALSE, pAdjustMethod='BH', contrast=c('diet','V','O'))

saveRDS(list(PS=PS,DS2=DS2_OUT,res=res),file.path(data_dir,'deseq2_pw.rds'))
