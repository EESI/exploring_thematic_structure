#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 2) stop("Two arguments must be supplied.", call.=FALSE)

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
library(rstanarm)


source('~/Dropbox/stm_microbiome/code_active/stm_functions.R')
source('~/Dropbox/stm_microbiome/code_active/nav_froz_fxns_3.R') 
source('~/Dropbox/stm_microbiome/code_active/performance_1.R')
source('~/Dropbox/stm_microbiome/code_active/framework.R')

seed_stan <- 123

models_dir <- '~/Dropbox/stm_microbiome/qiime_active/ag/models'

K <- as.integer(args[1])
fitn <- as.integer(args[2])

dat_fn <- file.path(models_dir,sprintf('%s/stm_k_%s.rds',K,K))
ko_fn <- gsub('stm_k_','ko_k_',dat_fn)
ko_fn <- gsub('\\.rds$','_fit%s\\.biom',ko_fn)
ko_fn <- sprintf(ko_fn,fitn)

dat <- readRDS(dat_fn)
fit <- dat$fits[[fitn+1]]

biom_file <- read_biom(ko_fn)
ko <- as.matrix(biom_data(biom_file))

if (K > 25){
  set.seed(453)
  META <- dat$meta
  
  eff <- estimateEffect(1:K ~ diet, fit, META, uncertainty='Global')
  eff_plot <- plot.estimateEffect(eff, 'diet', model=fit, 
                                  topics=1:K, method='difference',cov.value1=1,cov.value2=0)
  topics_subset <- order(abs(unlist(eff_plot$means)),decreasing=TRUE)[1:25]
  topics_subset <- topics_subset[order(unlist(eff_plot$means)[topics_subset])]
  
  ko <- ko[,topics_subset]
  
  subset_out <- list(eff=eff_plot,topics_subset=topics_subset)
}else{
  subset_out <- NULL
}

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


j1 <- 1
j2 <- 0
ko_topics <- as.integer(gsub('T','',colnames(ko)))
gene_mat <- as.data.frame(matrix(0,sum(sapply(target_kos,length)) * length(ko_topics),4,
                                 dimnames=list(NULL,c('topic','pw','ko','count'))))
for (i in seq_along(target_kos)){
  pw <- target_kos[[i]]
  for (k in seq_along(ko_topics)){
    j2 <- j2 + length(pw)
    gene_mat[j1:j2,1] <- ko_topics[k]
    gene_mat[j1:j2,2] <- target_pws[i]
    gene_mat[j1:j2,3] <- pw
    gene_mat[j1:j2,4] <- ko[pw,k]
    j1 <- j2+1
  }
}

stan_out <- rstanarm::stan_glmer.nb(count ~ (1|pw) + (1|topic) + (1|topic:pw),
                                 link='log',
                                 data=gene_mat,
                                 chains=4,cores=4,
                                 seed=seed_stan,
                                 iter=2000)

saveRDS(list(dat=gene_mat,
             K=fit$settings$dim$K,
             stan=stan_out,
             subset=subset_out,
             seed=seed_stan),
        file.path(models_dir,sprintf('%s/stan_lev3_k_%s_fit%s.rds',K,K,fitn)))
