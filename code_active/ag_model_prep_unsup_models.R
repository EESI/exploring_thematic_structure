#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 3) stop("At least one argument must be supplied (input file).n", call.=FALSE)

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
source('~/Dropbox/stm_microbiome/code_active/stm_ESTEP.R')


K <- as.integer(args[1])
rare_min <- as.integer(args[2])
ncores <- as.integer(args[3])
if (length(args) == 4) seed_adder <- as.integer(args[4]) else seed_adder <- 0


PS <- readRDS(sprintf('~/Dropbox/stm_microbiome/AG/Processing/01-raw/otus/ag_rare_%s_ps.rds',rare_min)) 

meta <- data.frame(sample_data(PS))
otu <- as.matrix(t(otu_table(PS))@.Data)
taxa <- as.matrix(tax_table(PS)@.Data)

meta <- meta %>% filter(BODY_HABITAT %in% 'UBERON:feces', !is.na(AGE_DECADE_RANK), !is.na(SEX))
otu <- otu[meta$SampleID,]


read_min <- 2
samp_min <- .99
taxa_filter <- colSums(otu < read_min) > floor(samp_min * NROW(otu)) 


cat('Removing ', sum(taxa_filter),  ' taxa (from ',NCOL(otu),' taxa) that have fewer than ',
    read_min, ' reads in less than ',(1-samp_min)*100,'% of samples.\n',sep='')


otu <- otu[,!taxa_filter]
taxa <- taxa[colnames(otu),]
meta <- filter(meta,SampleID %in% rownames(otu)) %>% data.frame()
rownames(meta) <- meta$SampleID


taxon_ids <- data.frame(long=colnames(otu),
                        ids=paste0('otu',1:NCOL(otu)),
                        row.names=paste0('otu',1:NCOL(otu)))

counts <- list(meta_clean=meta,
               table_clean=otu,
               ids=taxon_ids)
colnames(counts$table_clean) <- as.character(counts$ids$ids)

if (!identical(rownames(counts$table_clean),counts$meta_clean$SampleID)) stop('Please make table sample names match metadata sample IDs!\n')

vocab <- as.character(counts$ids$ids)
docs <- lapply(1:NROW(counts$table_clean),function(i) reshape_doc(counts$table_clean[i,],vocab)) ### something might be wrong!!!!
meta <- counts$meta_clean


cat('Total docs:',length(docs),'\n',
    'Total vocab:',length(vocab),'\n',
    'Total Meta:',nrow(meta),'x',ncol(meta),'\n',
    'Total Sex:',names(table(meta$SEX)),'=',c(table(meta$SEX)),'\n')


seed_fit <- 6457 + seed_adder
seed_froz <- 89 + seed_adder
seed_rare <- 5346 + seed_adder
seed_split <- 23 + seed_adder
seed_post <- 90 + seed_adder
seed_raw <- 98 + seed_adder


meta$FEMALE <- ifelse(meta$SEX == 'female',1,0)
meta$AGE_YEARS <- as.numeric(meta$AGE_YEARS)
meta$AGE <- as.vector(scale(as.numeric(meta$AGE_YEARS)))


set.seed(seed_split)
idx_train <- unlist(caret::createDataPartition(y=with(meta,as.factor(SEX):as.factor(AGE_CAT)),times=1,p=.8,list=TRUE))
otu_train <- otu[idx_train,]
meta_train <- meta[idx_train,]
docs_train <- docs[idx_train]
otu_test <- otu[-idx_train,]
meta_test <- meta[-idx_train,]
docs_test <- docs[-idx_train]



fit <- stm(documents=docs_train,vocab=vocab,K=K,
            data=meta_train,
            max.em.its=500,init.type='Spectral',
            seed=seed_fit, #645647
            verbose=TRUE,reportevery=25)

froz_train <- ctm_frozen2(fit,docs_train,vocab,
                               seed=seed_froz,max.em.its=500,emtol=1e-5,
                               verbose=TRUE,reportevery=5,
                               data=meta_train,covariate=NULL)

froz_test <- ctm_frozen2(fit,docs_test,vocab,
                              seed=seed_froz,max.em.its=500,emtol=1e-5,
                              verbose=TRUE,reportevery=5,
                              data=meta_test,covariate=NULL)



set.seed(seed_post)
nsims <- 500
theta_train <- thetaPosterior(froz_train,nsims=500,type='Local',documents=docs_train)
theta_train <- array(do.call('rbind',theta_train),dim=c(nsims,length(docs_train),K))

theta_test <- thetaPosterior(froz_test,nsims=500,type='Local',documents=docs_test)
theta_test <- array(do.call('rbind',theta_test),dim=c(nsims,length(docs_test),K))


require(snow)
require(doSNOW)

cl <- makeCluster(ncores,type='SOCK')
clusterExport(cl,c('qnormalize'))
registerDoSNOW(cl)

LABELS_TRAIN <- meta_train$FEMALE
LABELS_TEST <- meta_test$FEMALE

t1 <- Sys.time()
out <- foreach(i=seq_len(nsims),
               .verbose=FALSE,
               .errorhandling='pass',
               .packages=c('glmnet'),
               .combine='rbind') %dopar% {
                 
                 
                 TRAIN <- apply(theta_train[i,,],2,qnormalize)
                 TEST <- apply(theta_test[i,,],2,qnormalize)
                 
                 
                 en <- cv.glmnet(x=TRAIN,y=LABELS_TRAIN,
                                 family = 'binomial',alpha=.5,parallel=FALSE,
                                 standardize=TRUE,
                                 type.measure='auc')
                 
                 
                 en_pred <- ifelse(predict(en,newx=TEST,s='lambda.min')>0,1,0)
                 
                 
                 cbind(i,en_pred,LABELS_TEST)
                 
}
t2 <- Sys.time()
cat('Completed elastic net for',nsims,'simulations, taking',round(c(t2-t1),1),'seconds.\n')

cm_stm_qnorm <- caret::confusionMatrix(out[,2],out[,3])
print(cm_stm_qnorm)


t1 <- Sys.time()
out <- foreach(i=seq_len(nsims),
               .verbose=FALSE,
               .errorhandling='pass',
               .packages=c('glmnet'),
               .combine='rbind') %dopar% {
                 
                 
                 TRAIN <- apply(theta_train[i,,],2,scale)
                 TEST <- apply(theta_test[i,,],2,scale)
                 
                 
                 en <- cv.glmnet(x=TRAIN,y=LABELS_TRAIN,
                                 family = 'binomial',alpha=.5,parallel=FALSE,
                                 standardize=TRUE,
                                 type.measure='auc')
                 
                 
                 en_pred <- ifelse(predict(en,newx=TEST,s='lambda.min')>0,1,0)
                 
                 
                 cbind(i,en_pred,LABELS_TEST)
                 
               }
t2 <- Sys.time()
cat('Completed elastic net for',nsims,'simulations, taking',round(c(t2-t1),1),'seconds.\n')

cm_stm_znorm <- caret::confusionMatrix(out[,2],out[,3])
print(cm_stm_znorm)



set.seed(seed_raw)
en <- cv.glmnet(x=apply(otu_train,2,qnormalize),y=meta_train$FEMALE,
                family = 'binomial',alpha=.5,parallel=TRUE,
                standardize=TRUE,
                type.measure='auc')

en_pred <- ifelse(predict(en,newx=otu_test,s='lambda.min')>0,1,0)
cm_raw_qnorm <- caret::confusionMatrix(en_pred,meta_test$FEMALE)
print(cm_raw_qnorm)

en <- cv.glmnet(x=apply(otu_train,2,scale),y=meta_train$FEMALE,
                family = 'binomial',alpha=.5,parallel=TRUE,
                standardize=TRUE,
                type.measure='auc')

en_pred <- ifelse(predict(en,newx=otu_test,s='lambda.min')>0,1,0)
cm_raw_znorm <- caret::confusionMatrix(en_pred,meta_test$FEMALE)
print(cm_raw_znorm)



stopCluster(cl)




dir_name <- sprintf('~/Dropbox/stm_microbiome/data_active/AG_female/stm_s97_rarefied_%s_unsupervised',rare_min)
if (length(args) == 4) dir_name <- file.path(dir_name,'seeds',seed_adder)
dir.create(dir_name,showWarnings=FALSE,recursive=TRUE)
fit_filename <- paste0('fits_K_',K,'.rds')
saveRDS(list(fits=list(fit=fit,froz_train=froz_train,froz_test=froz_test),
             theta=list(theta_train=theta_train,theta_test=theta_test),
             docs=list(docs_train=docs_train,docs_test=docs_test),
             vocab=vocab,
             meta=list(meta_train=meta_train,meta_test=meta_test),
             otu=list(otu_train=otu_train,otu_test=otu_test),
             taxa=taxa,
             counts=counts,
             seeds=list(seed_fit=seed_fit,seed_rare=seed_rare,seed_froz=seed_froz,
                        seed_split=seed_split,seed_post=seed_post,seed_raw=seed_raw),
             filters=list(read_min=read_min,samp_min=samp_min,rare_min=rare_min,nsims=nsims),
             performance=list(cm_stm_qnorm=cm_stm_qnorm,cm_stm_znorm=cm_stm_znorm,
                              cm_raw_qnorm=cm_raw_qnorm,cm_raw_znorm=cm_raw_znorm)),
        file.path(dir_name,fit_filename))