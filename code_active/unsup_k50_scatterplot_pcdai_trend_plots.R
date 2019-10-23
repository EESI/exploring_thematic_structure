coef_k <- out$dx$lasso
top_ord <- order(coef_k,decreasing=TRUE)[c(1:4,(K-4+1):K)]
plot_violins(train_fit,train_meta,top_ord)

pcdai_id <- train_meta %>% 
   filter(DIAGNOSIS == 'CD',
          !is.na(PCDAI)) %>%
   arrange(PCDAI) %>%
   dplyr::select(SampleID,PCDAI)
th <- train_fit$theta
colnames(th) <- paste0('T',1:K)
rownames(th) <- train_meta$SampleID
th <- th[pcdai_id$SampleID,top_ord]

data.frame(th,PCDAI=pcdai_id$PCDAI,Rank=1:NROW(th)) %>%
   gather(Topic,Abundance,-Rank,-PCDAI) %>%
   mutate(Topic=factor(Topic,levels=colnames(th))) %>%
   ggplot(aes(Rank,Abundance,colour=PCDAI)) + facet_wrap(~Topic,nrow=2) + geom_point(alpha=.5) + scale_y_log10() +
   stat_smooth(method='lm',se=FALSE,colour='red') + xlab('PCDAI Rank') 






z <- counts$table_clean
z <- z[counts$meta_clean %>% filter(DIAGNOSIS == 'CD') %>% dplyr::select(SampleID) %>% unlist(),]
z <- z[,colSums(z) > 0]

b <- train_fit$beta$logbeta[[1]][1,]
names(b) <- paste0('otu',1:length(b))
o <- names(which(b > quantile(b,.95)))
length(o)

z <- z[,o]

df <- data.frame(z,counts$meta_clean %>% filter(DIAGNOSIS == 'CD') %>% dplyr::select(SampleID,PCDAI)) %>%   
   filter(!is.na(PCDAI)) %>%
   gather(OTU,Abundance,-SampleID,-PCDAI) %>%
   mutate(PCDAI_RANK = ranker(PCDAI)) %>%
   group_by(OTU) %>%
   mutate(Prevalence = sum(Abundance > 0)) %>%
   ungroup() %>%
   filter(Prevalence > 35) %>%
   mutate(Abundance = Abundance + 1) %>%
   rowwise() %>%
   mutate(Taxon = pretty_taxa_names(paste0(otu_taxa_xavier[as.character(counts$ids[unique(OTU),'long']),],collapse='|')),
          Facet = paste0(Taxon,' (',OTU,')'))
length(unique(df$OTU))

pdf('~/MiscOut/stm_unsup_k50_cdp.pdf', height=12, width=15)
ggplot(df, aes(PCDAI_RANK,Abundance,colour=Facet)) + 
   geom_point(data=dplyr::select(df,-Facet),aes(PCDAI_RANK,Abundance),colour='gray',alpha=.5,size=1) +
   facet_wrap(~Facet) + 
   geom_point(alpha=1,size=2,colour='black') +
   scale_y_log10() +
   stat_smooth(method='lm',se=FALSE,colour='red') +
   theme_bw() + xlab('PCDAI Rank')
dev.off()


z <- counts$table_clean
z <- z[counts$meta_clean %>% filter(DIAGNOSIS == 'CD') %>% dplyr::select(SampleID) %>% unlist(),]
z <- z[,colSums(z) > 0]

b <- fit$beta$logbeta[[1]][48,]
names(b) <- paste0('otu',1:length(b))
o <- names(which(b > quantile(b,.95)))
length(o)

z <- z[,o]

df <- data.frame(z,counts$meta_clean %>% filter(DIAGNOSIS == 'CD') %>% dplyr::select(SampleID,PCDAI)) %>%   
   filter(!is.na(PCDAI)) %>%
   gather(OTU,Abundance,-SampleID,-PCDAI) %>%
   mutate(PCDAI_RANK = ranker(PCDAI)) %>%
   group_by(OTU) %>%
   mutate(Prevalence = sum(Abundance > 0)) %>%
   ungroup() %>%
   filter(Prevalence > 35) %>%
   mutate(Abundance = Abundance + 1) %>%
   rowwise() %>%
   mutate(Taxon = pretty_taxa_names(paste0(otu_taxa_xavier[as.character(counts$ids[unique(OTU),'long']),],collapse='|')),
          Facet = paste0(Taxon,' (',OTU,')'))
length(unique(df$OTU))

pdf('~/MiscOut/stm_unsup_k50_cdn.pdf', height=12, width=15)
ggplot(df, aes(PCDAI_RANK,Abundance,colour=Facet)) + 
   geom_point(data=dplyr::select(df,-Facet),aes(PCDAI_RANK,Abundance),colour='gray',alpha=.5,size=1) +
   facet_wrap(~Facet) + 
   geom_point(alpha=1,size=2,colour='black') +
   scale_y_log10() +
   stat_smooth(method='lm',se=FALSE,colour='red') +
   theme_bw() + xlab('PCDAI Rank')
dev.off()
