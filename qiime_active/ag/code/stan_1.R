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

models_dir <- '~/Dropbox/stm_microbiome/qiime_active/gevers/models'

fns <- list.files(models_dir,recursive=TRUE,full.names=TRUE)
fns <- fns[grep('/stm_k_',fns)]

dat_fn <- fns[as.integer(args[1])]

dat <- readRDS(dat_fn)
fit <- dat$fits[[as.integer(args[2])]]
K <- fit$settings$dim$K

ko_fn <- file.path(models_dir,K,sprintf('ko_k_%s_%s.biom',K,names(dat$fits)[as.integer(args[2])]))

biom_file <- read_biom(ko_fn)
ko <- as.matrix(biom_data(biom_file))
ko <- ko[rowSums(ko)>0,]
kegg_meta <- kegg_metadata(biom_file,2)[rownames(ko)]

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
                   'genetic information processing',
                   'unknown or misclassification')

ko_colsums <- colSums(ko)
target_pws <- my_kegg_order[-c(16)] 
target_kos <- lapply(target_pws,function(pw) names(which(sapply(kegg_meta, function(x) any(x$pathway %in% pw)))))
names(target_kos) <- target_pws


j1 <- 1
j2 <- 0
gene_mat <- as.data.frame(matrix(0,sum(sapply(target_kos,length)) * K,5,
                                 dimnames=list(NULL,c('topic','pw','ko','count','offset'))))
for (i in seq_along(target_kos)){
  pw <- target_kos[[i]]
  for (k in seq_len(K)){
    j2 <- j2 + length(pw)
    gene_mat[j1:j2,1] <- k
    gene_mat[j1:j2,2] <- target_pws[i]
    gene_mat[j1:j2,3] <- pw
    gene_mat[j1:j2,4] <- ko[pw,k]
    gene_mat[j1:j2,5] <- ko_colsums[k]
    j1 <- j2+1
  }
}

stan_out <- rstanarm::stan_glmer.nb(count ~ (1|pw) + (1|topic) + (1|topic:pw),
                                 link='log',
                                 data=gene_mat,
                                 chains=4,cores=4,
                                 seed=seed_stan,
                                 iter=1000)

saveRDS(list(dat=gene_mat,
             K=fit$settings$dim$K,
             stan=stan_out,
             seed=seed_stan),
        gsub('.biom','.rds',gsub('\\/ko_k_','\\/stan_k_',ko_fn)))
