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

N_OTU <- 500
N_SAMP <- 100
SC_P <- c(.1,.25,.50,.75)
SC_N <- 10
SC_SAMP_P <- .75
N_SC <- 5


param_subset_idx <- with(param_grid,which(n_otu %in% N_OTU & n_samp %in% N_SAMP & sc_p %in% SC_P & sc_n %in% SC_N & sc_samp_p %in% SC_SAMP_P))
fns_subset <- paste0('data/dat_',param_subset_idx,'.rds')
LIST_SUBSET <- vector(mode='list',length=length(fns_subset))
for (i in seq_along(fns_subset))  LIST_SUBSET[[i]] <- readRDS(file.path(sim_dir,fns_subset[i]))
param_values <- param_grid[param_subset_idx,]
param_values$idx <- as.integer(rownames(param_values))


model_param_values <- with(LIST_SUBSET[[1]]$model_params,which(model %in% c('unsup') & K %in% c(15)))
P_MIN <- matrix(0.0,nrow(param_values),length(model_param_values),
                dimnames=list(1:nrow(param_values),LIST_SUBSET[[1]]$model_params[model_param_values,'ps']))
E_MIN <- E_MEAN <- P_MEAN <- P_MIN
for (p in seq_along(param_values$idx)){
  
  for (m in seq_along(model_param_values)){
    
    MM <- LIST_SUBSET[[p]]
    
    M <- MM$out[[model_param_values[m]]]
    if (class(M$data$fit) == 'try-error') {P_MIN[p,m] <- P_MEAN[p,m] <- E_MEAN[p,m] <- E_MIN[p,m] <- NA; next}
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
      
      OTU_SC <- colSums(OTU[GROUP,paste0('otu',OTUS[sc,])]) 
      OTU_SC <- OTU_SC/sum(OTU_SC)
      
      TOP_SC <- BETA2[,paste0('otu',OTUS[sc,])]
      TOP_SC <- TOP_SC/rowSums(TOP_SC)
      
      TOP_SC <- TOP_SC[,names(which(!(OTU_SC==0)))]
      OTU_SC <- OTU_SC[names(which(!(OTU_SC==0)))]
      
      apply(TOP_SC,1,function(p) jsd(p,OTU_SC))
    })
    dimnames(SCORE) <- list(paste0('T',1:nrow(SCORE)),paste0('sc',1:ncol(SCORE)))
    
    NSCORE_MIN <- matrix(NA,nrow(SCORE),ncol(SCORE),dimnames=dimnames(SCORE))
    for (i in 1:NROW(NSCORE_MIN)){
      for (j in 1:NCOL(NSCORE_MIN)){
        NSCORE_MIN[i,j] <- SCORE[i,j]/min(SCORE[i,-j])
      }
    }
    
    NSCORE_MEAN <- matrix(NA,nrow(SCORE),ncol(SCORE),dimnames=dimnames(SCORE))
    for (i in 1:NROW(NSCORE_MEAN)){
      for (j in 1:NCOL(NSCORE_MEAN)){
        NSCORE_MEAN[i,j] <- SCORE[i,j]/median(SCORE[i,-j]) ### using median
      }
    }
    
    NSCORE_MIN_THRES <- 1*(NSCORE_MIN<1)
    nscore_min_prop <- 1-mean(colSums(NSCORE_MIN_THRES)==0)
    nscore_min_entr <- sum(apply(NSCORE_MIN_THRES,1,function(x) entropy(table(which(x==1)))))
    
    NSCORE_MEAN_THRES <- 1*(NSCORE_MEAN<1)
    nscore_mean_prop <- 1-mean(colSums(NSCORE_MEAN_THRES)==0)
    nscore_mean_entr <- sum(apply(NSCORE_MEAN_THRES,1,function(x) entropy(table(which(x==1)))))
    
    P_MIN[p,m] <- nscore_min_prop
    E_MIN[p,m] <- nscore_min_entr
    
    P_MEAN[p,m] <- nscore_mean_prop
    E_MEAN[p,m] <- nscore_mean_entr
  }
}

P <- param_values %>%
  left_join(data.frame(P_MIN,idx=1:nrow(P)),by='idx') %>%
  left_join(data.frame(E_MEAN,idx=1:nrow(P)),by='idx') %>%
  gather(model,score,-(n_otu:idx)) %>%
  mutate(model=gsub('\\.x','\\.pmin',model),
         model=gsub('\\.y','\\.emean',model)) %>%
  separate(model,c('model','stat'),sep='\\.') %>%
  spread(stat,score)


P %>%
  ggplot(aes(emean,pmin,shape=model,colour=model)) + geom_point(size=5) + facet_grid(sc_p~sc_m) +
  scale_x_reverse() +
  scale_colour_brewer(type='qual',palette=2)







model_param_values <- with(LIST_SUBSET[[1]]$model_params,which(model %in% c('unsup') & K %in% c(25)))
S <- matrix(0.0,nrow(param_values),length(model_param_values),
            dimnames=list(1:nrow(param_values),LIST_SUBSET[[1]]$model_params[model_param_values,'ps']))
TH <- S
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
      GROUP <- G[[ifelse(sc <= 5, 1, ifelse(sc >= 10, 3, 2))]]
      #GROUP <- rownames(OTU)[rownames(OTU) %in% paste0('s',SAMPS[[sc]])]
      
      OTU_SC <- colSums(OTU[GROUP,paste0('otu',OTUS[sc,])]) 
      OTU_SC <- OTU_SC/sum(OTU_SC)
      
      TOP_SC <- BETA2[,paste0('otu',OTUS[sc,])]
      TOP_SC <- TOP_SC/rowSums(TOP_SC)
      
      TOP_SC <- TOP_SC[,names(which(!(OTU_SC==0)))]
      OTU_SC <- OTU_SC[names(which(!(OTU_SC==0)))]
      
      apply(TOP_SC,1,function(p) jsd(p,OTU_SC))
    })
    dimnames(SCORE) <- list(paste0('T',1:nrow(SCORE)),paste0('sc',1:ncol(SCORE)))
    
    
    THRES <- seq(min(SCORE),max(SCORE),by=.01)
    for (thres in THRES){
      col_thres <- colSums(1*(SCORE<thres))
      if (all(col_thres != 0)){
        break
      }
    }
    
    TH[p,m] <- thres
    S[p,m] <- sum((rowSums(1*(SCORE<thres))))/ncol(SCORE)
    
  }
}

P <- param_values %>%
  left_join(data.frame(TH,idx=1:nrow(P)),by='idx') %>%
  left_join(data.frame(S,idx=1:nrow(P)),by='idx') %>%
  gather(model,score,-(n_otu:idx)) %>%
  mutate(model=gsub('\\.x','\\.threshold',model),
         model=gsub('\\.y','\\.score',model)) %>%
  separate(model,c('model','stat'),sep='\\.') %>%
  spread(stat,score)


P %>%
  ggplot(aes(threshold,score,shape=model,colour=model)) + geom_point(size=5) + facet_grid(sc_p~sc_m) +
  scale_colour_brewer(type='qual',palette=2)
















model_param_values <- with(LIST_SUBSET[[1]]$model_params,which(model %in% c('unsup') & K %in% c(25)))
S <- matrix(0.0,nrow(param_values),length(model_param_values),
            dimnames=list(1:nrow(param_values),LIST_SUBSET[[1]]$model_params[model_param_values,'ps']))
TH <- S
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
      OTU_SC[!(names(OTU_SC) %in% paste0('otu',OTUS[sc,]))] <- min(BETA2)
      OTU_SC <- OTU_SC/sum(OTU_SC)
      
      TOP_SC <- BETA2    
      TOP_SC <- TOP_SC/rowSums(TOP_SC)
      
      TOP_SC <- TOP_SC[,names(which(!(OTU_SC==0)))]
      OTU_SC <- OTU_SC[names(which(!(OTU_SC==0)))]
      
      apply(TOP_SC,1,function(p) jsd(p,OTU_SC))
    })
    dimnames(SCORE) <- list(paste0('T',1:nrow(SCORE)),paste0('sc',1:ncol(SCORE)))
    
    
    THRES <- seq(min(SCORE),max(SCORE),by=.01)
    for (thres in THRES){
      col_thres <- colSums(1*(SCORE<thres))
      if (all(col_thres != 0)){
        break
      }
    }
    
    TH[p,m] <- thres
    S[p,m] <- sum((rowSums(1*(SCORE<thres))))/ncol(SCORE)
    
  }
}

P <- param_values %>%
  left_join(data.frame(TH,idx=1:nrow(P)),by='idx') %>%
  left_join(data.frame(S,idx=1:nrow(P)),by='idx') %>%
  gather(model,score,-(n_otu:idx)) %>%
  mutate(model=gsub('\\.x','\\.threshold',model),
         model=gsub('\\.y','\\.score',model)) %>%
  separate(model,c('model','stat'),sep='\\.') %>%
  spread(stat,score)


P %>%
  ggplot(aes(threshold,score,shape=model,colour=model)) + geom_point(size=5) + facet_grid(sc_p~sc_m) +
  scale_colour_brewer(type='qual',palette=2)






model_param_values <- with(LIST_SUBSET[[1]]$model_params,which(model %in% c('unsup') & K %in% c(15)))
S <- matrix(0.0,nrow(param_values),length(model_param_values),
            dimnames=list(1:nrow(param_values),LIST_SUBSET[[1]]$model_params[model_param_values,'ps']))
TH <- S
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
      OTU_SC[!(names(OTU_SC) %in% paste0('otu',OTUS[sc,]))] <- min(BETA2)
      OTU_SC <- OTU_SC/sum(OTU_SC)
      
      TOP_SC <- BETA2    
      #TOP_SC <- TOP_SC/rowSums(TOP_SC)
      
      TOP_SC <- TOP_SC[,names(which(!(OTU_SC==0)))]
      OTU_SC <- OTU_SC[names(which(!(OTU_SC==0)))]
      
      apply(TOP_SC,1,function(p) jsd(p,OTU_SC))
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
    
    TH[p,m] <- thres
    S[p,m] <- sum(rowSums(1*(SCORE<thres)) != 1)/K
    
  }
}

P <- param_values %>%
  left_join(data.frame(TH,idx=1:nrow(P)),by='idx') %>%
  left_join(data.frame(S,idx=1:nrow(P)),by='idx') %>%
  gather(model,score,-(n_otu:idx)) %>%
  mutate(model=gsub('\\.x','\\.threshold',model),
         model=gsub('\\.y','\\.score',model)) %>%
  separate(model,c('model','stat'),sep='\\.') %>%
  spread(stat,score)


P %>%
  ggplot(aes(threshold,score,shape=model,colour=model)) + geom_point(size=5) + facet_grid(sc_p~sc_m) +
  scale_colour_brewer(type='qual',palette=2) +
  xlim(0,1) + ylim(0,1) +
  theme(aspect.ratio=1)






model_param_values <- with(LIST_SUBSET[[1]]$model_params,which(model %in% c('unsup') & K %in% c(15)))
S <- matrix(0.0,nrow(param_values),length(model_param_values),
            dimnames=list(1:nrow(param_values),LIST_SUBSET[[1]]$model_params[model_param_values,'ps']))
TH <- S
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
      OTU_SC[!(names(OTU_SC) %in% paste0('otu',OTUS[sc,]))] <- NA
      OTU_SC <- dense_rank(dplyr::desc(OTU_SC))
      OTU_SC[is.na(OTU_SC)] <- length(OTU_SC)
      names(OTU_SC) <- colnames(OTU)
      OTU_SC <- OTU_SC[colnames(BETA2)]
      
      TOP_SC <- t(apply(BETA2,1,function(x) dense_rank(dplyr::desc(x))))    
      dimnames(TOP_SC) <- dimnames(BETA2)
      
      apply(TOP_SC,1,function(p) cor(p,OTU_SC,method='kendall'))
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
    
    TH[p,m] <- thres
    S[p,m] <- sum(rowSums(1*(SCORE<thres)) != 1)/K
    
  }
}

P <- param_values %>%
  left_join(data.frame(TH,idx=1:nrow(P)),by='idx') %>%
  left_join(data.frame(S,idx=1:nrow(P)),by='idx') %>%
  gather(model,score,-(n_otu:idx)) %>%
  mutate(model=gsub('\\.x','\\.threshold',model),
         model=gsub('\\.y','\\.score',model)) %>%
  separate(model,c('model','stat'),sep='\\.') %>%
  spread(stat,score)


P %>%
  ggplot(aes(threshold,score,shape=model,colour=model)) + geom_point(size=5) + facet_grid(sc_p~sc_m) +
  scale_colour_brewer(type='qual',palette=2) +
  scale_x_reverse() +
  theme(aspect.ratio=1)







model_param_values <- with(LIST_SUBSET[[1]]$model_params,which(model %in% c('unsup') & K %in% c(15)))
S <- matrix(0.0,nrow(param_values),length(model_param_values),
            dimnames=list(1:nrow(param_values),LIST_SUBSET[[1]]$model_params[model_param_values,'ps']))
TH <- S
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
      
      OTU_SC <- colSums(OTU[GROUP,paste0('otu',OTUS[sc,])])
      OTU_SC <- dense_rank(dplyr::desc(OTU_SC))
      names(OTU_SC) <- paste0('otu',OTUS[sc,])

      TOP_SC <- t(apply(BETA2,1,function(x) dense_rank(dplyr::desc(x))))  
      dimnames(TOP_SC) <- dimnames(BETA2)
      TOP_SC <- TOP_SC[,paste0('otu',OTUS[sc,])]
      
      apply(TOP_SC,1,function(p) cor(p,OTU_SC,method='kendall'))
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
    
    TH[p,m] <- thres
    S[p,m] <- sum(rowSums(1*(SCORE<thres)) != 1)/K
    
  }
}

P <- param_values %>%
  left_join(data.frame(TH,idx=1:nrow(P)),by='idx') %>%
  left_join(data.frame(S,idx=1:nrow(P)),by='idx') %>%
  gather(model,score,-(n_otu:idx)) %>%
  mutate(model=gsub('\\.x','\\.threshold',model),
         model=gsub('\\.y','\\.score',model)) %>%
  separate(model,c('model','stat'),sep='\\.') %>%
  spread(stat,score)


P %>%
  ggplot(aes(threshold,score,shape=model,colour=model)) + geom_point(size=5) + facet_grid(sc_p~sc_m) +
  scale_colour_brewer(type='qual',palette=2) +
  scale_x_reverse() +
  theme(aspect.ratio=1)








model_param_values <- with(LIST_SUBSET[[1]]$model_params,which(model %in% c('unsup') & K %in% c(15)))
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
    
    q_rank <- -qnormalize(1:N_OTU)
    
    SCORE <- sapply(1:nrow(OTUS), function(sc) {
      #GROUP <- G[[ifelse(sc <= 5, 1, ifelse(sc >= 10, 3, 2))]]
      GROUP <- rownames(OTU)[rownames(OTU) %in% paste0('s',SAMPS[[sc]])]
      
      OTU_SC <- colSums(OTU[GROUP,paste0('otu',OTUS[sc,])])
      OTU_SC <- q_rank[dense_rank(dplyr::desc(OTU_SC))]
      names(OTU_SC) <- paste0('otu',OTUS[sc,])
      
      TOP_SC <- t(apply(BETA2,1,function(x) q_rank[dense_rank(dplyr::desc(x))]))  
      dimnames(TOP_SC) <- dimnames(BETA2)
      TOP_SC <- TOP_SC[,paste0('otu',OTUS[sc,])]
      
      apply(TOP_SC,1,function(p) sqrt(sum(p - OTU_SC)^2))
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
    
    MAX50[p,m] <- median(apply(SCORE,2,max))
    TH[p,m] <- thres
    S[p,m] <- sum(rowSums(1*(SCORE<thres)) != 1)/K
    
  }
}

colnames(MAX50) <- paste0(colnames(MAX50),'.max50')
colnames(TH) <- paste0(colnames(TH),'.threshold')
colnames(S) <- paste0(colnames(S),'.score')

P <- param_values %>%
  left_join(data.frame(TH,idx=1:nrow(P)),by='idx') %>%
  left_join(data.frame(S,idx=1:nrow(P)),by='idx') %>%
  left_join(data.frame(MAX50,idx=1:nrow(P)),by='idx') %>%
  gather(model,score,-(n_otu:idx)) %>%
  separate(model,c('model','stat'),sep='\\.') %>%
  spread(stat,score) 


P %>%
  ggplot(aes(threshold,score,shape=model,colour=model)) + 
  geom_point(size=5) + facet_grid(sc_p~sc_m) + 
  geom_vline(aes(xintercept=max50,colour=model),size=1,alpha=.7) +
  scale_colour_brewer(type='qual',palette=2) +
  theme(aspect.ratio=1)








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


P %>%
  ggplot(aes(threshold,score,shape=model,colour=model)) + 
  geom_vline(aes(xintercept=max50,colour=model),size=1,alpha=.5,linetype=2) +
  geom_point(size=5) + facet_grid(sc_p~sc_m) + 
  scale_colour_brewer(type='qual',palette=2) +
  theme(aspect.ratio=1)








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
      #OTU_SC[!(names(OTU_SC) %in% paste0('otu',OTUS[sc,]))] <- 0
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


P %>%
  ggplot(aes(threshold,score,shape=model,colour=model)) + 
  geom_vline(aes(xintercept=max50,colour=model),size=1,alpha=.5,linetype=2) +
  geom_point(size=5) + facet_grid(sc_p~sc_m) + 
  scale_colour_brewer(type='qual',palette=2) +
  theme(aspect.ratio=1)






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

