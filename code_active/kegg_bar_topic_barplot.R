kegg_barplot_dir <- gsub('heatmaps','kegg_barplots',heatmap_dir)
dir.create(kegg_barplot_dir,showWarnings=FALSE)

pw_names <- table(unlist(kegg_metadata_risk_lev2))
pw_names <- names(pw_names[pw_names > 120])
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
  coef_ends <- as.numeric(c(names(which(sort(coef_k,TRUE) > 0)),names(which(sort(coef_k,TRUE) < 0))))
  
  
  
  p11 <- rbind(data.frame(Count=colSums(kegg1),Topic=names(colSums(kegg1)),Content='CD+'),
            data.frame(Count=colSums(kegg2),Topic=names(colSums(kegg2)),Content='CD-')) %>%
    filter(Topic %in% coef_ends) %>%
    mutate(Topic=factor(Topic,levels=coef_ends),
           Corr=Content:as.factor(ifelse(Topic %in% which(coef_k > 0),'CD+','CD-'))) %>%
    ggplot(aes(x=Topic,y=Count,fill=Corr)) + geom_bar(stat='identity',position='dodge',colour='black') +
    geom_vline(xintercept = K/2-.5,size=1.5,alpha=.25) +
    labs(x='',title=pw_name) +
    theme(aspect.ratio=2/5,legend.position='bottom') +
    scale_fill_brewer(palette='Paired')
  
  p12 <- ggplot(data.frame(Topic=1:K,Coefficient=sort(coef_k,decreasing=TRUE),Corr=sort(coef_k,decreasing=TRUE)),
               aes(Topic,Coefficient,colour=Corr)) +
    labs(colour='Corr') + 
    geom_vline(xintercept = K/2-.5,size=1.5,alpha=.25) + 
    geom_hline(yintercept = 0,size=1.5,alpha=.25) + 
    geom_point(size=3) +
    labs(colour='',y='Standardized Coefficient') +
    theme(aspect.ratio=2/5,legend.position='bottom') +
    scale_x_continuous(breaks=1:K,labels=names(sort(coef_k,decreasing=TRUE)),expand = c(.01,.01)) +
    scale_colour_distiller(palette='RdBu')
  
  
  
  
  
  
  sd1 <- vegan::diversity(beta1,index='shannon',MARGIN=2)
  sd2 <- vegan::diversity(beta2,index='shannon',MARGIN=2)
  
  # normaliezd by diversity
  p21 <- rbind(data.frame(Count=floor(colSums(kegg1)/sd1),Topic=names(colSums(kegg1)),Content='CD+'),
              data.frame(Count=floor(colSums(kegg2)/sd2),Topic=names(colSums(kegg2)),Content='CD-')) %>%
    filter(Topic %in% coef_ends) %>%
    mutate(Topic=factor(Topic,levels=coef_ends),
           Corr=Content:as.factor(ifelse(Topic %in% which(coef_k > 0),'CD+','CD-'))) %>%
    ggplot(aes(x=Topic,y=Count,fill=Corr)) + geom_bar(stat='identity',position='dodge',colour='black') +
    geom_vline(xintercept = K/2-.5,size=1.5,alpha=.25) +
    labs(x='',title=pw_name) +
    theme(aspect.ratio=2/5,legend.position='bottom') +
    scale_fill_brewer(palette='Paired')
  
  p22 <- ggplot(data.frame(Topic=1:K,Coefficient=sort(coef_k,decreasing=TRUE),Corr=sort(coef_k,decreasing=TRUE)),
               aes(Topic,Coefficient,colour=Corr)) +
    labs(colour='Corr') + 
    geom_vline(xintercept = K/2-.5,size=1.5,alpha=.25) + 
    geom_hline(yintercept = 0,size=1.5,alpha=.25) + 
    geom_point(size=3) +
    labs(colour='',y='Standardized Coefficient') +
    theme(aspect.ratio=2/5,legend.position='bottom') +
    scale_x_continuous(breaks=1:K,labels=names(sort(coef_k,decreasing=TRUE)),expand = c(.01,.01)) +
    scale_colour_distiller(palette='RdBu')
  
  
  
  
  
  uni1 <- vegan::specnumber(beta1,MARGIN=2)
  uni2 <- vegan::specnumber(beta2,MARGIN=2)
  
  # normaliezd by unique
  p31 <- rbind(data.frame(Count=floor(colSums(kegg1)/uni1),Topic=names(colSums(kegg1)),Content='CD+'),
              data.frame(Count=floor(colSums(kegg2)/uni2),Topic=names(colSums(kegg2)),Content='CD-')) %>%
    filter(Topic %in% coef_ends) %>%
    mutate(Topic=factor(Topic,levels=coef_ends),
           Corr=Content:as.factor(ifelse(Topic %in% which(coef_k > 0),'CD+','CD-'))) %>%
    ggplot(aes(x=Topic,y=Count,fill=Corr)) + geom_bar(stat='identity',position='dodge',colour='black') +
    geom_vline(xintercept = K/2-.5,size=1.5,alpha=.25) +
    labs(x='',title=pw_name) +
    theme(aspect.ratio=2/5,legend.position='bottom') +
    scale_fill_brewer(palette='Paired')
  
  p32 <- ggplot(data.frame(Topic=1:K,Coefficient=sort(coef_k,decreasing=TRUE),Corr=sort(coef_k,decreasing=TRUE)),
               aes(Topic,Coefficient,colour=Corr)) +
    labs(colour='Corr') + 
    geom_vline(xintercept = K/2-.5,size=1.5,alpha=.25) + 
    geom_hline(yintercept = 0,size=1.5,alpha=.25) + 
    geom_point(size=3) +
    labs(colour='',y='Standardized Coefficient') +
    theme(aspect.ratio=2/5,legend.position='bottom') +
    scale_x_continuous(breaks=1:K,labels=names(sort(coef_k,decreasing=TRUE)),expand = c(.01,.01)) +
    scale_colour_distiller(palette='RdBu')
  
  
  
  
  p41 <- rbind(data.frame(Freq=colSums(kegg1)/colSums(kegg_risk1),Topic=names(colSums(kegg1)),Content='CD+'),
               data.frame(Freq=colSums(kegg2)/colSums(kegg_risk2),Topic=names(colSums(kegg2)),Content='CD-')) %>%
    filter(Topic %in% coef_ends) %>%
    mutate(Topic=factor(Topic,levels=coef_ends),
           Corr=Content:as.factor(ifelse(Topic %in% which(coef_k > 0),'CD+','CD-'))) %>%
    ggplot(aes(x=Topic,y=Freq,fill=Corr)) + geom_bar(stat='identity',position='dodge',colour='black') +
    geom_vline(xintercept = K/2-.5,size=1.5,alpha=.25) +
    labs(x='',title=pw_name) +
    theme(aspect.ratio=2/5,legend.position='bottom') +
    scale_fill_brewer(palette='Paired')
  
  p42 <- ggplot(data.frame(Topic=1:K,Coefficient=sort(coef_k,decreasing=TRUE),Corr=sort(coef_k,decreasing=TRUE)),
                aes(Topic,Coefficient,colour=Corr)) +
    labs(colour='Corr') + 
    geom_vline(xintercept = K/2-.5,size=1.5,alpha=.25) + 
    geom_hline(yintercept = 0,size=1.5,alpha=.25) + 
    geom_point(size=3) +
    labs(colour='',y='Standardized Coefficient') +
    theme(aspect.ratio=2/5,legend.position='bottom') +
    scale_x_continuous(breaks=1:K,labels=names(sort(coef_k,decreasing=TRUE)),expand = c(.01,.01)) +
    scale_colour_distiller(palette='RdBu')
  
  
  
  
  pdf(file.path(kegg_barplot_dir,file_name), height=10, width=10)
  gridExtra::grid.arrange(nrow=2,p11,p12)
  gridExtra::grid.arrange(nrow=2,p21,p22)
  gridExtra::grid.arrange(nrow=2,p31,p32)
  gridExtra::grid.arrange(nrow=2,p41,p42)
  dev.off()
  
}
