kegg_barplot_dir <- gsub('heatmaps','kegg_dotplots',heatmap_dir)
dir.create(kegg_barplot_dir,showWarnings=FALSE)

pw_names <- table(unlist(kegg_metadata_risk_lev2))
pw_names <- names(pw_names[pw_names > 50])
pw_names <- pw_names[!(pw_names %in% c('None','Poorly Characterized'))]





for (pw_name in pw_names){
  
  file_name <- paste0(gsub(' ','_',pw_name),'.pdf')
  
  pw1 <- sapply(rownames(kegg_risk1), function(x) pw_name %in% kegg_metadata_risk_lev2[[x]])
  pw2 <- sapply(rownames(kegg_risk2), function(x) pw_name %in% kegg_metadata_risk_lev2[[x]])
  
  kegg1 <- kegg_risk1[pw1,]
  kegg2 <- kegg_risk2[pw2,]
  colnames(kegg1) <- 1:K
  colnames(kegg2) <- 1:K
  
  sum(rownames(kegg1) %in% names(ko_names1))
  sum(rownames(kegg2) %in% names(ko_names2))
  
  
  names(coef_k) <- 1:K
  names(coef_k_unnorm) <- 1:K
  coef_ends <- as.numeric(c(names(which(sort(coef_k,TRUE) > 0)),names(which(sort(coef_k,TRUE) < 0))))
  
  


  
#   post_theta <- thetaPosterior(fit,nsims=500,type='Global')
#   post_theta_marg <- sapply(post_theta, function(x) colMeans(x))
#   post1 <- rowMeans(post_theta_marg[,fit$settings$covariates$X[,2] > 0])
#   post0 <- rowMeans(post_theta_marg[,fit$settings$covariates$X[,2] == 0])
  
  
p11 <- rbind(data.frame(Count=colSums(kegg1),Topic=names(colSums(kegg1)),Content='CD+'),
             data.frame(Count=colSums(kegg2),Topic=names(colSums(kegg2)),Content='CD-')) %>%
  filter(Topic %in% coef_ends) %>%
  spread(Content,Count) %>%
  mutate(Topic=factor(Topic,levels=coef_ends),
         Prev=as.factor(ifelse(Topic %in% which(coef_k_unnorm > 0),'CD+','CD-')),
         FC=log2(`CD+`/`CD-`),
         Coef=coef_k_unnorm[as.character(Topic)],
         Sig=ifelse(Topic %in% sig_corr_topics,1,0)) %>%
  ggplot(aes(x=Topic,y=FC,colour=Coef,size=abs(Coef))) +
  geom_blank() +
  geom_vline(aes(xintercept = sum(Coef > 0)+.5),size=1.1,colour='black') +
  geom_hline(yintercept=0,size=1.1,colour='black') +
  geom_hline(aes(yintercept = quantile(FC,.95)),linetype=2) + 
  geom_hline(aes(yintercept = quantile(FC,.05)),linetype=2) +
  geom_point() +
  geom_label_repel(data=. %>% filter(Sig == 1),aes(x=Topic,y=FC,label=paste0('Topic ',Topic)),point.padding=unit(1,'lines')) +
#   geom_label_repel(data=. %>% filter(FC > quantile(FC,.95)),aes(x=Topic,y=FC,label=paste0('Topic ',Topic)),point.padding=unit(1,'lines')) +
#   geom_label_repel(data=. %>% filter(FC < quantile(FC,.05)),aes(x=Topic,y=FC,label=paste0('Topic ',Topic)),point.padding=unit(1,'lines')) +
  labs(x='',y='Fold Difference',title=pw_name,colour='Coefficient') +
  scale_colour_gradient2(low=scales::muted('blue'),mid='lightgray',high=scales::muted('red'),midpoint=0) +
  scale_size(range = c(4, 12)) +
  scale_x_discrete(expand=c(-1,1)) + 
  guides(size=FALSE) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_text(size=14),
        axis.title.y=element_text(size=18),
        title=element_text(size=22),
        legend.title=element_text(size=15),
        legend.text=element_text(size=10),
        aspect.ratio=2/4,
        legend.position='right')
  
  
  
  pdf(file.path(kegg_barplot_dir,file_name), height=5, width=10)
  print(p11)
  dev.off()
  
}