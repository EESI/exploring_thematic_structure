N_OTU <- 500
N_SAMP <- 100
SC_P <- c(.1,.25,.50,.75)
SC_N <- 15
SC_SAMP_P <- .75
N_SC <- 5



param_subset_idx <- with(param_grid,which(n_otu %in% N_OTU & n_samp %in% N_SAMP & sc_p %in% SC_P & sc_n %in% SC_N & sc_samp_p %in% SC_SAMP_P))
fns_subset <- paste0('data/dat_',param_subset_idx,'.rds')
LIST_SUBSET <- vector(mode='list',length=length(fns_subset))
for (i in seq_along(fns_subset))  LIST_SUBSET[[i]] <- readRDS(file.path(sim_dir,fns_subset[i]))
param_values <- param_grid[param_subset_idx,]
param_values$idx <- as.integer(rownames(param_values))


model_param_values <- with(LIST_SUBSET[[1]]$model_params,which(model %in% c('unsup') & K %in% c(25)))
P <- matrix(0.0,nrow(param_values),length(model_param_values),
            dimnames=list(1:nrow(param_values),LIST_SUBSET[[1]]$model_params[model_param_values,'ps']))
for (p in seq_along(param_values$idx)){
  
  for (m in seq_along(model_param_values)){
    
    MM <- LIST_SUBSET[[p]]
    
    M <- MM$out[[model_param_values[m]]]
    if (class(M$data$fit) == 'try-error') {P[p,m] <- NA; next}
    K <- M$data$fit$settings$dim$K
    OTUS <- M$data$sim$sc_otus
    SAMPS <- M$data$sim$samples
    LOGBETA <- M$data$fit$beta$logbeta[[1]]
    BETA <- exp(LOGBETA)
    
    
    OTU <- t(otu_table(M$data$ps))
    META <- sample_data(M$data$ps)
    
    
    sc_df <- data.frame(sc=paste0('sc',1:(3*N_SC)),
                        sc_idx=paste0('sc',rep(1:5,3)),
                        sc_type=rep(c('1','0','b'),each=N_SC))
    SIG <- M$results$sig
    VOCAB <- M$data$fit$vocab
    RANKS <- t(apply(LOGBETA,1,dense_rank))
    ZRANKS <- (RANKS-min(RANKS))/(max(RANKS)-min(RANKS)) 
    #ZZRANKS <- (RANKS-mean(RANKS))/sd(RANKS)
    dimnames(RANKS) <- list(paste0('topic',1:K),VOCAB)
    dimnames(ZRANKS) <- list(paste0('topic',1:K),VOCAB)
    #dimnames(ZZRANKS) <- list(paste0('topic',1:K),VOCAB)
    rownames(OTUS) <- paste0('sc',1:nrow(OTUS))
    ROOF <- mean(sort(unique(as.vector(ZRANKS)),decreasing = TRUE)[1:SC_N])
    
    RANKS <- t(sapply(1:nrow(OTUS), function(sc) rowMeans(RANKS[,paste0('otu',OTUS[sc,])])))
    ZRANKS <- t(sapply(1:nrow(OTUS), function(sc) rowMeans(ZRANKS[,paste0('otu',OTUS[sc,])])))
    #ZZRANKS <- t(sapply(1:nrow(OTUS), function(sc) rowMeans(ZZRANKS[,OTUS[sc,]])))
    rownames(RANKS) <- paste0('sc',1:nrow(OTUS))
    rownames(ZRANKS) <- paste0('sc',1:nrow(OTUS))
    #rownames(ZZRANKS) <- paste0('sc',1:nrow(OTUS))
    
    
    # as.data.frame(ZRANKS) %>%
    #   mutate(sc=rownames(.)) %>%
    #   gather(topic,zrank,-sc) %>%
    #   left_join(SIG,by='topic') %>%
    #   left_join(sc_df,by='sc') %>%
    #   mutate(topic=factor(topic,levels=colnames(ZRANKS)),
    #          sc=factor(sc,levels=rownames(ZRANKS)),
    #          sc_type=factor(sc_type,levels=c('1','0','b'))) %>%
    #   ggplot() + 
    #   facet_grid(sc_type~sc_idx) +
    #   geom_hline(yintercept=ROOF,size=1,linetype=2) +
    #   geom_point(aes(topic,zrank,colour=sig,size=zrank)) + 
    #   scale_size(range = c(0, 5)) +
    #   scale_colour_manual(limits=c(-1,0,1),values=c('blue','black','red'),drop=FALSE) + 
    #   geom_label_repel(data=. %>% group_by(sc) %>% filter(zrank > quantile(zrank,.90)), 
    #                    aes(topic,zrank,label=topic,fill=sig),colour='white',size=4) +
    #   scale_fill_manual(limits=c(-1,0,1),values=c('blue','black','red'),drop=FALSE) +
    #   labs(x='',y='Normalized Ranks') +
    #   ggtitle(paste(M$data$param,collapse=' | ')) +
    #   theme(axis.text.x=element_blank(),
    #         axis.ticks.x=element_blank(),
    #         axis.text.y=element_text(size=14),
    #         axis.title.y=element_text(size=18),
    #         title=element_text(size=22),
    #         legend.position='none') +
    #   ylim(0,1)
    
    
    
    #   ZZRANKS <- matrix(NA,nrow(ZRANKS),ncol(ZRANKS),dimnames=dimnames(ZRANKS))
    #   for (i in 1:NROW(ZRANKS)){
    #     for (j in 1:NCOL(ZRANKS)){
    #       ZZRANKS[i,j] <- ZRANKS[i,j]/mean(ZRANKS[-i,j])
    #     }
    #   }
    #   rowMeans(ZZRANKS)
    #   mean(rowMeans(ZZRANKS))
    #   
    
    
    
    ZRANKS <- (1/ZRANKS)^2
    ZZRANKS <- matrix(NA,nrow(ZRANKS),ncol(ZRANKS),dimnames=dimnames(ZRANKS))
    for (i in 1:NROW(ZRANKS)){
      for (j in 1:NCOL(ZRANKS)){
        ZZRANKS[i,j] <- ZRANKS[i,j]/mean(ZRANKS[-i,j])
      }
    }
    #   rowMeans(ZZRANKS)
    P[p,m] <- mean(ZZRANKS)
  }
}

P <- param_values %>%
  left_join(data.frame(P,idx=1:nrow(P)),by='idx') %>%
  gather(model,score,-(n_otu:idx))


p1 <- P %>%
  filter(model=='PS_EQ') %>%
  ggplot(aes(sc_m,score)) + facet_grid(.~sc_p) + geom_line(size=1.1) + geom_point(size=2)
p2 <- P %>%
  filter(model != 'PS_EQ') %>%
  ggplot(aes(sc_m,score,colour=model)) + facet_grid(.~sc_p) + geom_line(size=1.1) + geom_point(size=2) +
  theme(legend.position = 'bottom')
grid.arrange(p1,p2,ncol=1,top=paste0(K,' Topics'))
