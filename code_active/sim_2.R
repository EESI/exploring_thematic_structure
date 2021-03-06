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
library(ggrepel)
library(snow)
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

entropy <- function(x,K){
  f <- x/K
  -sum(f * log2(f))
}


sim_dir <- sprintf('~/Dropbox/stm_microbiome/output_active/simulation/sim_%s',floor(as.numeric(Sys.time())))
dir.create(sim_dir,showWarnings=FALSE,recursive=TRUE)
dir.create(file.path(sim_dir,'figures'),showWarnings=FALSE,recursive=TRUE)
dir.create(file.path(sim_dir,'data'),showWarnings=FALSE,recursive=TRUE)
dir.create(file.path(sim_dir,'output'),showWarnings=FALSE,recursive=TRUE)


seed_general <- 10
seed_rare <- 4


vec_K <- c(25,50)
vec_n_otu <- c(500,1000)
vec_n_samp <- c(50,100,500)
vec_sc_n <- c(10,15,30)
vec_sc_m <- c(5,10,25,50,100)
vec_sc_s <- c(.1,.5,1)
param_grid <- expand.grid(K=vec_K,n_otu=vec_n_otu,n_samp=vec_n_samp,sc_n=vec_sc_n,sc_m=vec_sc_m,sc_s=vec_sc_s,stringsAsFactors = FALSE)


filename <- paste(names(param_grid[1,]),param_grid[1,],sep='-',collapse='_')

param_grid <- param_grid[1:60,]

nc <- 60 #ifelse(nrow(param_grid) < 20,nrow(param_grid),20)
cl <- makeCluster(nc,type='SOCK')
clusterExport(cl,c('logisticnormalcpp','gradcpp','lhoodcpp','hpbcpp'))
registerDoSNOW(cl)
SIM <- foreach(p=1:nrow(param_grid),.verbose=TRUE,.errorhandling="stop",
               .packages=c("Matrix","stm","phyloseq","ggplot2","dplyr","tidyr","DESeq2","ggrepel")) %dopar% {

  PARAMS <- param_grid[p,]
  
  K <- PARAMS$K
  rare_min <- 1000
  
  
  n_otu <- PARAMS$n_otu
  n_samp <- PARAMS$n_samp
  otu_m <- 5
  otu_s <- .1
  otu_p <- .3
  
  
  n_sc <- 5
  sc_n <- PARAMS$sc_n
  sc_p_b <- .3
  sc_p_1 <- .3
  sc_p_0 <- .3
  sc_m_b <- PARAMS$sc_m
  sc_m_1 <- PARAMS$sc_m
  sc_m_0 <- PARAMS$sc_m
  sc_s_1 <- PARAMS$sc_s
  sc_s_0 <- PARAMS$sc_s
  sc_s_b <- PARAMS$sc_s
  sc_samp_p <- .5
    
  
  set.seed(seed_general)
  
  
  n_otu_uneq <- sample(100:5000,n_samp,replace=TRUE)
  
  OTU_EQ <- replicate(n_samp,rzinb(n_otu,otu_m,otu_s,otu_p))
  
  # max(OTU_EQ)
  # colMeans(OTU_EQ==0)
  # colSums(OTU_EQ)
  # hist(log(OTU_EQ[OTU_EQ>0]),100)
  
  
  sc_samp <- c(lapply(1:n_sc,function(x) sample(1:(n_samp/2),sc_samp_p*(n_samp/2),replace=FALSE)),
               lapply(1:n_sc,function(x) sample((n_samp/2 + 1):n_samp,sc_samp_p*(n_samp/2),replace=FALSE)),
               lapply(1:n_sc,function(x) c(sample(1:(n_samp/2),sc_samp_p*(n_samp/2/2),replace=FALSE),
                                           sample((n_samp/2 + 1):n_samp,sc_samp_p*(n_samp/2/2),replace=FALSE))))
  sc_otu <- matrix(sample(1:n_otu,n_sc*sc_n*3),ncol=sc_n)
  sc_df <- data.frame(subcom=paste0('subcom',1:(3*n_sc)),
                      subcom_idx=paste0('subcom',rep(1:5,3)),
                      subcom_type=rep(c('1','0','b'),each=n_sc))
  
  for (i in 1:n_sc){
    OTU_EQ[sc_otu[i,],sc_samp[[i]]] <- rzinb(sc_n*length(sc_samp[[i]]),sc_m_1,sc_s_1,sc_p_1)
    
    OTU_EQ[sc_otu[i+n_sc,],sc_samp[[i+n_sc]]] <- rzinb(sc_n*length(sc_samp[[i+n_sc]]),sc_m_0,sc_s_0,sc_p_0)
    
    OTU_EQ[sc_otu[i+n_sc*2,],sc_samp[[i+n_sc*2]]] <- rzinb(sc_n*length(sc_samp[[i+n_sc*2]]),sc_m_b,sc_s_b,sc_p_b)
  }
  
  # max(OTU_EQ)
  # colMeans(OTU_EQ==0)
  # colSums(OTU_EQ)
  # hist(log(OTU_EQ[OTU_EQ>0]),100)
  
  
  OTU_UNEQ_counts <- sapply(1:n_samp,function(i) table(sample(1:n_otu,size=n_otu_uneq[i],replace=TRUE,prob=OTU_EQ[,i])))
  OTU_UNEQ <- matrix(0,n_otu,n_samp)
  for (i in 1:n_samp) OTU_UNEQ[as.integer(names(OTU_UNEQ_counts[[i]])),i] <- OTU_UNEQ_counts[[i]]
  
  max(OTU_UNEQ)
  colMeans(OTU_UNEQ==0)
  colSums(OTU_UNEQ)
  hist(log(OTU_UNEQ[OTU_UNEQ>0]),100)
  
  
  TAXA <- matrix(paste0('otu',1:n_otu,sep=''),ncol=1,dimnames=list(paste0('otu',1:n_otu,sep=''),'Taxon'))
  
  
  colnames(OTU_EQ) <- paste0('s',1:ncol(OTU_EQ))
  rownames(OTU_EQ) <- paste0('otu',1:nrow(OTU_EQ))
  META <- data.frame(Sample=paste0('s',1:ncol(OTU_EQ)),Effect=rep(1:0,each=n_samp/2))
  rownames(META) <- META$Sample
  PS_EQ <- phyloseq(otu_table(OTU_EQ,taxa_are_rows=TRUE),
                    sample_data(META),
                    tax_table(TAXA))
  PS_EQ <- prune_samples(sample_sums(PS_EQ)>0,PS_EQ)
  PS_EQ <- prune_taxa(taxa_sums(PS_EQ)>0,PS_EQ)
  
  
  colnames(OTU_UNEQ) <- paste0('s',1:ncol(OTU_UNEQ))
  rownames(OTU_UNEQ) <- paste0('otu',1:nrow(OTU_UNEQ))
  META <- data.frame(Sample=paste0('s',1:ncol(OTU_UNEQ)),Effect=rep(1:0,each=n_samp/2))
  rownames(META) <- META$Sample
  PS_UNEQ <- phyloseq(otu_table(OTU_UNEQ,taxa_are_rows=TRUE),
                      sample_data(META),
                      tax_table(TAXA))
  PS_UNEQ <- prune_samples(sample_sums(PS_UNEQ)>0,PS_UNEQ)
  PS_UNEQ <- prune_taxa(taxa_sums(PS_UNEQ)>0,PS_UNEQ)
  
  
  PS_RARE <- rarefy_even_depth(PS_UNEQ,
                               rngseed=seed_rare,
                               sample.size=rare_min,
                               replace=FALSE,
                               trimOTUs=TRUE,
                               verbose=FALSE)
  PS_RARE <- prune_samples(sample_sums(PS_RARE)>0,PS_RARE)
  PS_RARE <- prune_taxa(taxa_sums(PS_RARE)>0,PS_RARE)
  
  
  habdds <- phyloseq_to_deseq2(PS_UNEQ,~Effect)
  gm_mean <- function(x, na.rm=TRUE) exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  habdds_temp <- try(estimateSizeFactors(habdds), silent=TRUE)
  if (class(habdds_temp) == 'try-error') {
    geo_means <- apply(counts(habdds), 1, gm_mean)
    habdds <- estimateSizeFactors(habdds,geoMeans=geo_means)
    rm(habdds_temp)
  }else{    
    habdds <- estimateSizeFactors(habdds)
  }
  habdds <- estimateDispersions(habdds)
  diagvst <- getVarianceStabilizedData(habdds)
  diagvst[diagvst < 0] <- 0
  diagvst <- round(diagvst)
  
  PS_DESEQ2 <- PS_UNEQ
  otu_table(PS_DESEQ2) <- otu_table(diagvst, taxa_are_rows = TRUE)
  PS_DESEQ2 <- prune_samples(sample_sums(PS_DESEQ2)>0,PS_DESEQ2)
  PS_DESEQ2 <- prune_taxa(taxa_sums(PS_DESEQ2)>0,PS_DESEQ2)
  
  
  PS <- list(PS_EQ=PS_EQ,PS_UNEQ=PS_UNEQ,PS_RARE=PS_RARE,PS_DESEQ2=PS_DESEQ2)

  
  # PS_PLOT <- PS_EQ
  # ord <- ordinate(PS_PLOT,'MDS','bray')
  # plot_ordination(PS_PLOT,ord,type='sample',color='Effect')
  
  
  ps <- names(PS)
  model <- c('unsup','prev')
  model_grid <- expand.grid(ps=ps,model=model,
                            stringsAsFactors=FALSE)
  
  
  LIST_FIT <- vector(mode='list',length=nrow(model_grid))
  for (i in 1:nrow(model_grid)){
    
    param <- model_grid[i,]
    ps <- PS[[param$ps]]
    dat <- preproc_ps(ps)
    mod <- param$model
    
    if (mod == 'unsup'){
      fit <- stm(documents=dat$docs,vocab=dat$vocab,K=K,
                 data=dat$meta,
                 max.em.its=500,init.type='Spectral',
                 control=list(kappa.enet=.5),
                 verbose=FALSE,reportevery=25)
    } else if (mod == 'prev'){
      fit <- stm(prevalence=~Effect,
                 documents=dat$docs,vocab=dat$vocab,K=K,
                 data=dat$meta,
                 max.em.its=500,init.type='Spectral',
                 control=list(kappa.enet=.5),
                 verbose=FALSE,reportevery=25)
    } else if (mod == 'cont'){
      fit <- stm(content=~Effect,
                 documents=dat$docs,vocab=dat$vocab,K=K,
                 data=dat$meta,
                 max.em.its=500,init.type='Spectral',
                 control=list(kappa.enet=.5),
                 verbose=FALSE,reportevery=25)
    } else if (mod == 'prev_cont'){
      fit <- stm(prevalence=~Effect,
                 content=~Effect,
                 documents=dat$docs,vocab=dat$vocab,K=K,
                 data=dat$meta,
                 max.em.its=500,init.type='Spectral',
                 control=list(kappa.enet=.5),
                 verbose=FALSE,reportevery=25)
    }
    
    LIST_FIT[[i]] <- list(fit=fit,param=param,ps=ps)
    
  }
  
  
  
  
  LIST_OUT <- vector(mode='list',length=length(LIST_FIT))
  for (m in seq_along(LIST_FIT)){
    
    out <- LIST_FIT[[m]]
    
    fit <- out$fit
    meta <- data.frame(sample_data(out$ps))
    form <- as.formula(paste0('1:',deparse(K),' ~ Effect'))
    eff <- estimateEffect(form, fit, meta, uncertainty='Global')
    eff <- plot.estimateEffect(eff, 'Effect', model=fit, topics=1:K, method='difference',cov.value1=1,cov.value2=0)
    
    p1 <- as.data.frame(t(rbind(sapply(eff$cis,'['),unlist(eff$means)))) %>%
      dplyr::rename(Lower=`2.5%`,Upper=`97.5%`,Mean=V3) %>%
      mutate(Topic=factor(1:K,levels=1:K),
             Sig=factor(sign(Lower) + sign(Upper),levels=c(-2,0,2))) %>%
      ggplot(aes(x=Topic,y=Mean,ymin=Lower,ymax=Upper,colour=Sig)) + 
      geom_hline(yintercept=0) + 
      geom_pointrange() + 
      theme(legend.position='none') + scale_colour_manual(limits=c(-2,0,2),values=c('blue','black','red'))
    
    sig <- sign(sapply(eff$cis,function(x) sum(sign(x))))
    
    vocab <- fit$vocab
    
    logbeta <- fit$beta$logbeta
    for (i in seq_along(logbeta)) logbeta[[i]][is.infinite(logbeta[[i]])] <- -100
    for (i in seq_along(logbeta)) logbeta[[i]][logbeta[[i]] < -100] <- -100
    sig <- data.frame(topic=paste0('topic',1:K),sig=factor(sig,levels=c('-1','0','1')),stringsAsFactors=FALSE)
    
    
    ll <- t(apply(sc_otu,1,function(x) rowSums(logbeta[[1]][,vocab %in% paste0('otu',x)])))
    dimnames(ll) <- list(paste0('subcom',1:NROW(ll)),paste0('topic',1:NCOL(ll)))
    logbeta_rank <- t(apply(logbeta[[1]],1,function(x) dplyr::dense_rank(dplyr::desc(x))))
    ranks <- lapply(1:nrow(sc_otu),function(i) logbeta_rank[,vocab %in% paste0('otu',sc_otu[i,])])
    names(ranks) <- paste0('subcom',1:length(ranks))
    for (i in seq_along(ranks)) dimnames(ranks[[i]]) <- list(paste0('topic',1:NROW(ranks[[i]])),paste0('otu',sc_otu[i,]))
    
    p2 <- as.data.frame(ll) %>%
      mutate(subcom=rownames(.)) %>%
      gather(topic,logprob,-subcom) %>%
      left_join(sig,by='topic') %>%
      left_join(sc_df,by='subcom') %>%
      mutate(topic=factor(topic,levels=colnames(ll)),
             subcom=factor(subcom,levels=rownames(ll))) %>%
      ggplot() + 
      facet_grid(subcom_type~subcom_idx) +
      geom_point(aes(topic,logprob,colour=sig,size=logprob)) + 
      scale_size(range = c(0, 5)) +
      scale_colour_manual(limits=c(-1,0,1),values=c('blue','black','red'),drop=FALSE) + 
      geom_label_repel(data=. %>% group_by(subcom) %>% filter(logprob > quantile(logprob,.90)), 
                       aes(topic,logprob,label=topic,fill=sig),colour='white',size=4) +
      scale_fill_manual(limits=c(-1,0,1),values=c('blue','black','red'),drop=FALSE) +
      labs(x='',y='Log Probability') +
      ggtitle(paste(out$param,collapse=' | ')) +
      theme(axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.text.y=element_text(size=14),
            axis.title.y=element_text(size=18),
            title=element_text(size=22),
            legend.position='none') +
      ylim(-1000,0)
    
    top15 <- sapply(ranks,function(x) apply(x,1,function(y) sum(y <= sc_n)))
    max_counts <- lapply(1:nrow(sc_otu), 
                         function(i) sort(table(apply(logbeta[[1]][,vocab %in% paste0('otu',sc_otu[i,])],2,which.max)),decreasing=TRUE))
    max_counts_entropy <- sapply(max_counts,entropy,K)
    sum_max_counts_entropy <- sum(max_counts_entropy)
    
    LIST_OUT[[m]] <-   list(figures=list(p1=p1,p2=p2,eff=eff),
                            data=out,
                            results=list(top15=top15,
                                         max_counts=max_counts,
                                         max_counts_entropy=max_counts_entropy,
                                         sum_max_counts_entropy=sum_max_counts_entropy))
  }
  
  list(out=LIST_OUT,params=param_grid[p,],PARAMS=PARAMS)
  
}
stopCluster(cl)
               


saveRDS(list(figures=list(p1=p1,p2=p2),
             out=LIST_OUT,
             phyloseqs=PS,
             parameters=param_grid[i,]),
        file.path(sim_dir,'data',paste0('dat_',filename,'.rds')))

