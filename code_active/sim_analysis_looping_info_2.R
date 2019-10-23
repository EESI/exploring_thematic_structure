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
library(snow)
library(tibble)
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
jsd <- function(p,q){
  m <- .5*(p+q)
  sqrt(.5*(sum(p*log(p/m)) + sum(q*log(q/m))))
}
entropy <- function(x){
  x <- x[x!=0]
  f <- x/sum(x)
  -sum(f * log2(f))
}
norm10 <- function(x) (x-min(x))/(max(x)-min(x))
bcd <- function(x,y) sum(abs(x-y)/sum(x+y))




sim_dir <- '~/Dropbox/stm_microbiome/output_active/simulation/sim_1484014177'


model_grid <- read_csv(file.path(sim_dir,'data/model_grid.csv'))
param_grid <- read_csv(file.path(sim_dir,'data/param_grid.csv'))


GRID <- list(n_otu=unique(param_grid$n_otu),
             n_samp=unique(param_grid$n_samp),
             sc_n=unique(param_grid$sc_n),
             sc_samp_p=unique(param_grid$sc_samp_p),
             k=unique(model_grid$K),
             group_type <- c('group_sc','group_complete'))

figure_grid <- with(GRID,expand.grid(n_otu=n_otu,n_samp=n_samp,sc_n=sc_n,sc_samp_p=sc_samp_p,k=k,group_type=group_type,
                                     stringsAsFactors=FALSE))

N_SC <- 5
for (f in 1:nrow(figure_grid)){
  
  figure_values <- figure_grid[f,]
  
  N_OTU <- figure_values$n_otu
  N_SAMP <- figure_values$n_samp
  SC_N <- figure_values$sc_n
  SC_SAMP_P <- figure_values$sc_samp_p
  MODEL_K <- figure_values$k
  GROUP_TYPE <- figure_values$group_type
  
  figure_title <- sprintf('N_OTU=%s | N_SAMP=%s | SC_N=%s | SC_SAMP_P=%s | K=%s | Group=%s',
                          N_OTU,N_SAMP,SC_N,SC_SAMP_P,MODEL_K,GROUP_TYPE)
    
  param_subset_idx <- with(param_grid,
                           which(n_otu %in% N_OTU & n_samp %in% N_SAMP & sc_n %in% SC_N & sc_samp_p %in% SC_SAMP_P))
  fns_subset <- paste0('data/dat_',param_subset_idx,'.rds')
  LIST_SUBSET <- vector(mode='list',length=length(fns_subset))
  for (i in seq_along(fns_subset))  LIST_SUBSET[[i]] <- readRDS(file.path(sim_dir,fns_subset[i]))
  param_values <- param_grid[param_subset_idx,]
  param_values$idx <- as.integer(rownames(param_values))
  
  
  model_param_values <- with(LIST_SUBSET[[1]]$model_params,which(model %in% c('unsup') & K %in% c(MODEL_K)))
  S <- matrix(0.0,nrow(param_values),length(model_param_values),
              dimnames=list(1:nrow(param_values),LIST_SUBSET[[1]]$model_params[model_param_values,'ps']))
  TH <- MAX50 <- S
  for (p in seq_along(param_values$idx)){
    
    for (m in seq_along(model_param_values)){
      
      MM <- LIST_SUBSET[[p]]
      
      M <- MM$out[[model_param_values[m]]]
      if (class(M$data$fit) == 'try-error') {S[p,m] <- TH[p,m] <- NA; next}
      K <- M$data$fit$settings$dim$K
      VOCAB <- M$data$fit$vocab
      OTUS <- M$data$sim$sc_otus
      SAMPS <- M$data$sim$samples
      LOGBETA2 <- M$data$fit$beta$logbeta[[1]]
      #LOGBETA2[LOGBETA2==-Inf] <- min(LOGBETA2[LOGBETA2>-Inf]) ### think about this more
      BETA2 <- exp(LOGBETA2)
      colnames(BETA2) <- VOCAB
      
      
      OTU <- t(otu_table(M$data$ps))
      META <- sample_data(M$data$ps)
      SIG <- M$results$sig
      
      
      sc_df <- data.frame(sc=paste0('sc',1:(3*N_SC)),
                          sc_idx=paste0('sc',rep(1:5,3)),
                          sc_type=rep(c('1','0','b'),each=N_SC))
      
      G1 <- rownames(OTU)[rownames(OTU) %in% paste0('s',1:(nrow(OTU)/2))]
      G2 <- rownames(OTU)[!(rownames(OTU) %in% paste0('s',1:(nrow(OTU)/2)))]
      GB <- rownames(OTU)
      G <- list(G1=G1,G2=G2,GB=GB)
      
      SCORE <- sapply(1:nrow(OTUS), function(sc) {
        
        if (GROUP_TYPE == 'group_sc'){
          GROUP <- rownames(OTU)[rownames(OTU) %in% paste0('s',SAMPS[[sc]])]
        }
        if (GROUP_TYPE == 'group_complete'){
          GROUP <- G[[ifelse(sc <= 5, 1, ifelse(sc >= 10, 3, 2))]] 
        }
                
        OTU_SC <- colSums(OTU[GROUP,]) 
        OTU_SC <- OTU_SC[colnames(BETA2)]
        OTU_SC[!(names(OTU_SC) %in% paste0('otu',OTUS[sc,]))] <- 0
        OTU_SC <- OTU_SC/sum(OTU_SC)
        
        TOP_SC <- BETA2   
        TOP_SC[TOP_SC<1e-300] <- 1e-300
        
        apply(TOP_SC,1,function(p) sum(ifelse(OTU_SC==0,0,OTU_SC*log(OTU_SC/p))))
      })
      dimnames(SCORE) <- list(paste0('T',1:nrow(SCORE)),paste0('sc',1:ncol(SCORE)))
      
      THRES <- seq(min(SCORE),max(SCORE),by=.01)
      for (thres in THRES){
        col_thres <- colSums(1*(SCORE<thres))
        if (sum(col_thres) >= K){
          thres <- thres_last
          break
        }
        thres_last <- thres
      }
      
      MAX50[p,m] <- median(apply(SCORE,2,min))
      TH[p,m] <- thres
      S[p,m] <- sum(rowSums(1*(SCORE<thres)) > 1)/K
      
    }
  }
  
  colnames(MAX50) <- paste0(colnames(MAX50),'.max50')
  colnames(TH) <- paste0(colnames(TH),'.threshold')
  colnames(S) <- paste0(colnames(S),'.score')
  
  P <- param_values %>%
    left_join(data.frame(TH,idx=1:nrow(TH)),by='idx') %>%
    left_join(data.frame(S,idx=1:nrow(S)),by='idx') %>%
    left_join(data.frame(MAX50,idx=1:nrow(MAX50)),by='idx') %>%
    gather(model,score,-(n_otu:idx)) %>%
    separate(model,c('model','stat'),sep='\\.') %>%
    spread(stat,score) 
  
  
  FIG <-   P %>%
    ggplot(aes(threshold,score,shape=model,colour=model)) + 
    geom_vline(aes(xintercept=max50,colour=model),size=1,alpha=.5,linetype=2) +
    geom_point(size=5) + facet_wrap(sc_p~sc_m,ncol=4,scales='free_x') + 
    scale_colour_brewer(type='qual',palette=2) +
    theme(aspect.ratio=1) +
    ggtitle(figure_title)
  
  pdf(file.path(sim_dir,'figures',sprintf('sim_fig_%s.pdf',f)),height=10,width=10)
  print(FIG)
  dev.off()

}
  





model_param_values <- with(LIST_SUBSET[[1]]$model_params,which(model %in% c('unsup') & K %in% c(50)))
S <- matrix(0.0,nrow(param_values),length(model_param_values),
            dimnames=list(1:nrow(param_values),LIST_SUBSET[[1]]$model_params[model_param_values,'ps']))
TH <- MAX50 <- S
for (p in seq_along(param_values$idx)){
  
  for (m in seq_along(model_param_values)){
    
    MM <- LIST_SUBSET[[p]]
    
    M <- MM$out[[model_param_values[m]]]
    if (class(M$data$fit) == 'try-error') {S[p,m] <- TH[p,m] <- NA; next}
    K <- M$data$fit$settings$dim$K
    VOCAB <- M$data$fit$vocab
    OTUS <- M$data$sim$sc_otus
    SAMPS <- M$data$sim$samples
    LOGBETA2 <- M$data$fit$beta$logbeta[[1]]
    LOGBETA2[LOGBETA2==-Inf] <- min(LOGBETA2[LOGBETA2>-Inf])
    BETA2 <- exp(LOGBETA2)
    colnames(BETA2) <- VOCAB
    
    
    OTU <- t(otu_table(M$data$ps))
    META <- sample_data(M$data$ps)
    SIG <- M$results$sig
    
    
    sc_df <- data.frame(sc=paste0('sc',1:(3*N_SC)),
                        sc_idx=paste0('sc',rep(1:5,3)),
                        sc_type=rep(c('1','0','b'),each=N_SC))
    
    G1 <- rownames(OTU)[rownames(OTU) %in% paste0('s',1:(nrow(OTU)/2))]
    G2 <- rownames(OTU)[!(rownames(OTU) %in% paste0('s',1:(nrow(OTU)/2)))]
    GB <- rownames(OTU)
    G <- list(G1=G1,G2=G2,GB=GB)
    
    SCORE <- sapply(1:nrow(OTUS), function(sc) {
      #GROUP <- G[[ifelse(sc <= 5, 1, ifelse(sc >= 10, 3, 2))]]
      GROUP <- rownames(OTU)[rownames(OTU) %in% paste0('s',SAMPS[[sc]])]
      
      OTU_SC <- colSums(OTU[GROUP,]) 
      OTU_SC <- OTU_SC[colnames(BETA2)]
      OTU_SC <- OTU_SC/sum(OTU_SC)
      OTU_SC[!(names(OTU_SC) %in% paste0('otu',OTUS[sc,]))] <- 0
      
      TOP_SC <- BETA2   
      TOP_SC[TOP_SC<1e-300] <- 1e-300
      
      apply(TOP_SC,1,function(p) bcd(p,OTU_SC))
    })
    dimnames(SCORE) <- list(paste0('T',1:nrow(SCORE)),paste0('sc',1:ncol(SCORE)))
    
    THRES <- seq(min(SCORE),max(SCORE),by=.01)
    for (thres in THRES){
      col_thres <- colSums(1*(SCORE<thres))
      if (sum(col_thres) >= K){
        thres <- thres_last
        break
      }
      thres_last <- thres
    }
    
    MAX50[p,m] <- median(apply(SCORE,2,min))
    TH[p,m] <- thres
    S[p,m] <- sum(rowSums(1*(SCORE<thres)) > 1)/K
    
  }
}

colnames(MAX50) <- paste0(colnames(MAX50),'.max50')
colnames(TH) <- paste0(colnames(TH),'.threshold')
colnames(S) <- paste0(colnames(S),'.score')

P <- param_values %>%
  left_join(data.frame(TH,idx=1:nrow(TH)),by='idx') %>%
  left_join(data.frame(S,idx=1:nrow(S)),by='idx') %>%
  left_join(data.frame(MAX50,idx=1:nrow(MAX50)),by='idx') %>%
  gather(model,score,-(n_otu:idx)) %>%
  separate(model,c('model','stat'),sep='\\.') %>%
  spread(stat,score) 


P %>%
  ggplot(aes(threshold,score,shape=model,colour=model)) + 
  geom_vline(aes(xintercept=max50,colour=model),size=1,alpha=.5,linetype=2) +
  geom_point(size=5) + facet_grid(sc_p~sc_m) + 
  scale_colour_brewer(type='qual',palette=2) +
  theme(aspect.ratio=1)

