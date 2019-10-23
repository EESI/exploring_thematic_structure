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

#+ fig.width=16, fig.height=16



rds_out <- readRDS('~/Qiime/Children/xavier_data_loading.rds')
otu_lookup <- rds_out$otu_lookup
data_mat_xavier <- rds_out$data_mat_xavier
metadata_xavier <- rds_out$metadata_xavier
kegg_metadata_xavier <- rds_out$kegg_metadata_xavier
kegg_mat_xavier <- rds_out$kegg_mat_xavier


counts <- collapse_table(otu_lookup,kegg_mat_xavier,metadata_xavier,
                         taxon=NULL,method='rescaled',
                         kegg=TRUE, cap=TRUE,roof=10^5,
                         filtermin=FALSE,mc=5,ps=.1,
                         filtermax=FALSE,pw=.9,pa=.9)

counts$meta <- metadata_xavier[metadata_xavier$SampleID %in% rownames(counts$table),]
counts$meta_clean <- counts$meta[!is.na(counts$meta$ISOLATION_SOURCE),]
counts$table_clean <- counts$table[counts$meta_clean$SampleID,]
rownames(counts$meta_clean) <- counts$meta_clean$SampleID

set.seed(4354)
n_oversamp <- abs(diff(table(counts$meta_clean$DIAGNOSIS)))-5
id_oversamp <- sample(counts$meta_clean$SampleID[counts$meta_clean$DIAGNOSIS == 'Not IBD'],n_oversamp,replace=TRUE)
saveRDS(id_oversamp,'~/Qiime/Children/xavier_oversamp_loading')
counts$meta_clean <- rbind(counts$meta_clean,counts$meta_clean[id_oversamp,])
counts$table_clean <- rbind(counts$table_clean,counts$table_clean[id_oversamp,])

vocab <- as.character(counts$ids$ids)
docs <- lapply(1:nrow(counts$table_clean),function(i) reshape_doc(counts$table_clean[i,],vocab))
meta <- counts$meta_clean
meta$DIAGNOSIS <- as.factor(ifelse(meta$DIAGNOSIS=='CD',1,0))
meta$ISOLATION_SOURCE <- as.factor(meta$ISOLATION_SOURCE)
meta$COHORT <- as.factor(meta$COHORT)

prep <- prepDocuments(docs,vocab,meta,lower.thresh=0) ### removes terms apparently

K <- 75
fit1 <- stm(prep$documents, prep$vocab, K=K,
            prevalence=~ISOLATION_SOURCE + DIAGNOSIS + COHORT,
            content=~DIAGNOSIS,
            gamma.prior='Pooled',
            kappa.prior='L1',
            max.em.its=75, init.type='Spectral',
            data=prep$meta)
checkResiduals(fit1,prep$documents) # 10

saveRDS(fit1,paste0('~/Qiime/Children/stm_oversamp_kegg_cohort_content_pooled_binary_l1_k',K,'_',unique_idx(),'.rds'))








counts <- collapse_table(otu_lookup,data_mat_xavier,metadata_xavier,
                         taxon=NULL,method='rescaled',
                         kegg=FALSE, cap=TRUE,roof=10^5,
                         filtermin=FALSE,mc=5,ps=.01,
                         filtermax=FALSE,pw=.95,pa=.9)

counts$meta <- metadata_xavier[metadata_xavier$SampleID %in% rownames(counts$table),]
counts$meta_clean <- counts$meta[!is.na(counts$meta$ISOLATION_SOURCE),]
counts$table_clean <- counts$table[counts$meta_clean$SampleID,]
rownames(counts$meta_clean) <- counts$meta_clean$SampleID

### using the same ids as before to oversample same 'subjects'
id_oversamp <- readRDS('~/Qiime/Children/xavier_oversamp_loading')
counts$meta_clean <- rbind(counts$meta_clean,counts$meta_clean[id_oversamp,])
counts$table_clean <- rbind(counts$table_clean,counts$table_clean[id_oversamp,])

vocab <- as.character(counts$ids$ids)
docs <- lapply(1:nrow(counts$table_clean),function(i) reshape_doc(counts$table_clean[i,],vocab))
meta <- counts$meta_clean
meta$DIAGNOSIS <- as.factor(ifelse(meta$DIAGNOSIS=='CD',1,0))
meta$ISOLATION_SOURCE <- as.factor(meta$ISOLATION_SOURCE)
meta$COHORT <- as.factor(meta$COHORT)

prep <- prepDocuments(docs,vocab,meta,lower.thresh=0) ### removes terms apparently

K <- 75
fit2 <- stm(prep$documents, prep$vocab, K=K,
            prevalence=~ISOLATION_SOURCE + DIAGNOSIS + COHORT,
            content=~DIAGNOSIS,
            gamma.prior='Pooled',
            kappa.prior='L1',
            max.em.its=75, init.type='Spectral',
            data=prep$meta)
checkResiduals(fit2,prep$documents) # 363

saveRDS(fit2,paste0('~/Qiime/Children/stm_oversamp_otus_cohort_content_pooled_binary_l1_k',K,'_',unique_idx(),'.rds'))



