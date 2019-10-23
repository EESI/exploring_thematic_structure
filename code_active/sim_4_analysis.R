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
library(snow)
library(ggrepel)
library(doSNOW)

source('~/Dropbox/stm_microbiome/code_active/stm_functions.R')
source('~/Dropbox/stm_microbiome/code_active/nav_froz_fxns_3.R') 
source('~/Dropbox/stm_microbiome/code_active/performance_1.R')
source('~/Dropbox/stm_microbiome/code_active/framework.R')

rzinb <- function(n,m,s,p){
  
  # p = probability of a zero
  # m = effect size, mean of nb
  # s = variance of nb
  
  # z = 1, presence of zero
  
  z <- ifelse(rbinom(n=n,size=1,prob=p)==1,0,1)
  w <- rnbinom(n=n,mu=m,size=s)
  
  w*z
  
}

reshape_doc <- function(doc,vocab){
  doc <- doc[doc != 0]
  idx <- match(names(doc),vocab)
  out <- rbind(idx,doc)
  rownames(out) <- c('idx','count')
  return(out)
}

preproc_ps <- function(ps){
  otu <- t(as.data.frame(otu_table(ps)))
  meta <- data.frame(sample_data(ps))
  vocab <- colnames(otu)
  docs <- lapply(1:NROW(otu),function(i) reshape_doc(otu[i,],vocab))
  
  return(list(otu=otu,meta=meta,vocab=vocab,docs=docs))
}

extract_effect <- function(out){
  
  fit <- out$fit
  K <- fit$settings$dim$K
  meta <- data.frame(sample_data(out$ps))
  
  eff <- estimateEffect(1:K ~ Effect, fit, meta, uncertainty='Global')
  eff <- plot.estimateEffect(eff, 'Effect', model=fit, topics=1:K, method='difference',cov.value1=1,cov.value2=0)
  
  p <- as.data.frame(t(rbind(sapply(eff$cis,'['),unlist(eff$means)))) %>%
    dplyr::rename(Lower=`2.5%`,Upper=`97.5%`,Mean=V3) %>%
    mutate(Topic=factor(1:K,levels=1:K),
           Sig=factor(sign(Lower) + sign(Upper),levels=c(-2,0,2))) %>%
    ggplot(aes(x=Topic,y=Mean,ymin=Lower,ymax=Upper,colour=Sig)) + 
    geom_hline(yintercept=0) + 
    geom_pointrange() + 
    theme(legend.position='none') + scale_colour_manual(limits=c(-2,0,2),values=c('blue','black','red'))
  
  sig <- sign(sapply(eff$cis,function(x) sum(sign(x))))
  
  return(list(eff=eff,sig=sig,fig=p))
  
}
entropy <- function(x){
  f <- x/sum(x)
  -sum(f * log2(f))
}
gen_ent <- function(L,model_grid,m,f,k,idx=NULL){
  if (is.null(idx)){
    sapply(L,function(x) sum(sapply(x$out[[with(model_grid,which(ps==m & model==f & K==k))]]$results$max_counts,entropy)))
  }else{
    sapply(L,function(x) sum(sapply(x$out[[with(model_grid,which(ps==m & model==f & K==k))]]$results$max_counts,entropy)[1:idx]))
  }
}


sims_dir <- '~/Dropbox/stm_microbiome/output_active/simulation/'
sims_name <- list.files(sims_dir)
sims_name <- sims_name[grepl('[0-9]',sims_name)]
sim_dir <- file.path(sims_dir,sims_name[1])


vec_n_otu <- c(500,1000,2500)
vec_n_samp <- c(50,100,500)
vec_sc_n <- c(10,15,30)
vec_sc_m <- c(5,10,25,50,100)
vec_sc_s <- c(.1,.5,1)


ps <- c('PS_EQ','PS_UNEQ','PS_RARE','PS_DESEQ2')
model <- c('unsup','prev')
model_grid <- expand.grid(ps=ps,model=model,K=c(15,25,50),stringsAsFactors=FALSE)
param_grid <- expand.grid(n_otu=vec_n_otu,n_samp=vec_n_samp,sc_n=vec_sc_n,sc_m=vec_sc_m,sc_s=vec_sc_s,stringsAsFactors = FALSE) %>% 
  arrange(n_otu,desc(sc_s),sc_n,n_samp)


N_OTU <- 500
N_SAMP <- 100
SC_S <- .1
SC_N <- 30
N_SC <- 5

param_subset_idx <- with(param_grid,which(n_otu==N_OTU & n_samp==N_SAMP & sc_s==SC_S & sc_n==SC_N))
fns_subset <- paste0('data/dat_',param_subset_idx,'.rds')
LIST_SUBSET <- vector(mode='list',length=length(fns_subset))
for (i in seq_along(fns_subset))  LIST_SUBSET[[i]] <- readRDS(file.path(sim_dir,fns_subset[i]))
param_values <- param_grid[param_subset_idx,]
param_values$idx <- rownames(param_values)





M <- LIST_SUBSET[[5]]$out[[5]]
K <- M$data$fit$settings$dim$K
OTUS <- M$data$sim$sc_otus
BETA <- exp(M$data$fit$beta$logbeta[[1]])

THRESH <- seq(0,max(BETA)/2,length=100)
SIM <- lapply(1:nrow(OTUS), function(sc) t(sapply(seq_along(THRESH), function(thr) rowSums(ifelse(BETA[,OTUS[sc,]] >= THRESH[thr],1,0))/SC_N)))
df <- data.frame(do.call('rbind',SIM),Threshold=rep(THRESH,nrow(OTUS)),SC=rep(1:nrow(OTUS),each=length(THRESH))) %>%
  gather(Topic,Proportion,-Threshold,-SC) %>%
  mutate(Topic=factor(gsub('X','',Topic),levels=1:K),
         thresh_rank=dense_rank(Threshold),
         SC=factor(SC,levels=1:nrow(OTUS))) %>%
  filter(thresh_rank==50 & SC %in% 1:3) %>%
  group_by(thresh_rank,SC) %>%
  mutate(Above=sum(Proportion>0))
  
  
df %>% 
  filter(SC %in% 1:5) %>%
  ggplot(aes(Threshold,Proportion,colour=Topic)) +
  facet_wrap(~SC) +
  geom_vline(xintercept=1/N_OTU) +
  geom_line() + 
  xlim(0,df %>% filter(Proportion>0) %>% summarise(max(Threshold)) %>% unlist()) +
  ylim(0,df %>% filter(Threshold>0) %>% summarise(max(Proportion)) %>% unlist())




df <- data.frame(model_grid,matrix(0,nrow(model_grid),length(param_subset_idx)))
colnames(df) <- c('ps','model','K',param_subset_idx)
for (i in 1:nrow(df)) df[i,4:8] <- gen_ent(LIST_SUBSET,model_grid,df[i,1],df[i,2],df[i,3])
df %>%
  gather(idx,entropy,-ps,-model,-K) %>%
  mutate(idx=gsub('X','',idx)) %>%
  left_join(param_values,by='idx') %>%
  ggplot(aes(sc_m,entropy,colour=ps,linetype=model)) + facet_grid(~K) + geom_point(size=3) + geom_line(size=1.25) +
  ylim(0,entropy(rep(1,SC_N))*N_SC*3)


model_score <- function(LIST_SUBSET,i,j){
  SS <- LIST_SUBSET[[i]]
  SSOUT <- SS$out[[j]]
  SCOTUS <- SSOUT$data$sim$sc_otus
  WTOPICS <- unlist(SSOUT$figures$eff$means)
#   scores <- sapply(1:ncol(SCOTUS), function(k) mean(colSums(exp(SSOUT$data$fit$beta$logbeta[[1]])[,SCOTUS[k,]]*WTOPICS)))
  scores <- sapply(1:nrow(SCOTUS), function(k) rowSums(exp(SSOUT$data$fit$beta$logbeta[[1]])[,SCOTUS[k,]]))
#   names(scores) <- paste0('sc',1:length(scores))
  dimnames(scores) <- list(paste0('topic',1:nrow(scores)),paste0('sc',1:ncol(scores)))
  
  norm_scores <- matrix(NA,nrow(scores),ncol(scores),dimnames=dimnames(scores))
  for (ii in 1:ncol(scores)){
    norm_scores[,ii] <- scores[,ii]/rowSums(scores[,-ii])
  }

  score <- apply(norm_scores,2,function(x) x * WTOPICS)

  colSums(score)
}

df <- data.frame(model_grid,m=rep(param_values$sc_m,each=nrow(model_grid)))
score_grid <- expand.grid(m=1:nrow(model_grid),p=1:nrow(param_values))
score_out <- t(sapply(1:nrow(score_grid), function(i) model_score(LIST_SUBSET,score_grid[i,2],score_grid[i,1])))
score_out <- data.frame(score_out,model_grid[score_grid$m,],m=param_values$sc_m[score_grid$p])
score_out %>%
  gather(sc,score,-ps,-model,-K,-m) %>%
  filter(m > 10) %>%
  mutate(sc=factor(sc,levels=paste0('sc',1:15)),
         K=factor(K,levels=c(15,25,50))) %>%
  ggplot(aes(ps,score,colour=K,shape=model)) +
  facet_grid(m~sc) + 
  geom_jitter(size=2) +
  geom_hline(yintercept=0) + 
  theme(aspect.ratio=2.5,
        axis.text.x = element_text(angle=45,hjust=1)) +
  coord_flip() 

# eq uneq rare deseq2
apply(sapply(5:8, function(k) unlist(LIST_SUBSET[[5]]$out[[k]]$figures$eff$means)),1,function(x) which.max(abs(x)))
apply(sapply(5:8, function(k) model_score(LIST_SUBSET,5,k)),1,which.max)




df <- data.frame(model_grid,matrix(0,nrow(model_grid),length(param_subset_idx)))
colnames(df) <- c('ps','model','K',param_subset_idx)
for (i in 1:nrow(df)) df[i,4:8] <- gen_ent(LIST_SUBSET,model_grid,df[i,1],df[i,2],df[i,3],10)
df %>%
  gather(idx,entropy,-ps,-model,-K) %>%
  mutate(idx=gsub('X','',idx)) %>%
  left_join(param_values,by='idx') %>%
  ggplot(aes(sc_m,entropy,colour=ps,linetype=model)) + facet_grid(~K) + geom_point(size=3) + geom_line(size=1.25) +
  ylim(0,entropy(rep(1,SC_N))*N_SC*2)


LIST_SUBSET[[5]]$out[[with(model_grid,which(model=='prev' & ps=='PS_DESEQ2' & K==15))]]$results$max_counts
LIST_SUBSET[[5]]$out[[with(model_grid,which(model=='prev' & ps=='PS_DESEQ2' & K==50))]]$results$max_counts
sum(sapply(LIST_SUBSET[[5]]$out[[with(model_grid,which(model=='unsup' & ps=='PS_DESEQ2' & K==15))]]$results$max_counts,entropy))
sum(sapply(LIST_SUBSET[[5]]$out[[with(model_grid,which(model=='prev' & ps=='PS_DESEQ2' & K==50))]]$results$max_counts,entropy))








LIST_SUBSET[[1]]$out[[model_subset_idx]]$results$max_counts
model_subset_idx <- with(model_grid,which(ps=='PS_DESEQ2' & model=='unsup' & K==50))
dot_matrix <- matrix(0,n_sc*3,model_grid[model_subset_idx,'K'])
dimnames(dot_matrix) <- list(1:nrow(dot_matrix),1:ncol(dot_matrix))
for (i in 1:nrow(dot_matrix)){
  max_counts <- LIST_SUBSET[[5]]$out[[model_subset_idx]]$results$max_counts[[i]]
  dot_matrix[i,as.integer(names(max_counts))] <- max_counts
}
max_count <- max(dot_matrix)
dot_matrix <- as.data.frame(dot_matrix)
dot_matrix$sc <- rownames(dot_matrix)
dot_matrix %>%
  gather(topic,count,-sc) %>%
  mutate(topic=as.numeric(topic),
         sc=as.numeric(sc),
         topic=factor(topic, levels=1:max(topic)),
         sc=factor(sc, levels=1:max(sc))) %>%
  ggplot(aes(topic,sc,size=count,alpha=count)) + geom_point() + scale_size(range=c(0,max_count*2))
model_subset_idx <- with(model_grid,which(ps=='PS_DESEQ2' & model=='unsup' & K==50))
dot_matrix <- matrix(0,n_sc*3,model_grid[model_subset_idx,'K'])
dimnames(dot_matrix) <- list(1:nrow(dot_matrix),1:ncol(dot_matrix))
for (i in 1:nrow(dot_matrix)){
  max_counts <- LIST_SUBSET[[3]]$out[[model_subset_idx]]$results$max_counts[[i]]
  dot_matrix[i,as.integer(names(max_counts))] <- max_counts
}
max_count <- max(dot_matrix)
dot_matrix <- as.data.frame(dot_matrix)
dot_matrix$sc <- rownames(dot_matrix)
dot_matrix %>%
  gather(topic,count,-sc) %>%
  mutate(topic=as.numeric(topic),
         sc=as.numeric(sc),
         topic=factor(topic, levels=1:max(topic)),
         sc=factor(sc, levels=1:max(sc))) %>%
  ggplot(aes(topic,sc,size=count,alpha=count)) + geom_point() + scale_size(range=c(0,max_count*2))



