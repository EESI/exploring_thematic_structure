collapse_table <- function(otu_taxa,data_mat,metadata,
                           taxon=NULL,method='rescaled',
                           kegg=FALSE,cap=TRUE,roof=10^6,
                           filtermin=FALSE,mc=0,ps=.01,
                           filtermax=FALSE,pw=.9,pa=.9){
  
  if (is.null(dim(data_mat))){
    data_mat <- t(counts$table_clean)
    rownames(data_mat) <- counts$id[rownames(data_mat),'long']
  }
    
  if (!is.null(taxon)){
    
    otu_taxa <- otu_taxa[rownames(data_mat),]
    target_taxon <- data.frame(long=apply(otu_taxa[,1:taxon],1,function(i) paste0(i,collapse='|')),
                               lowest_name=otu_taxa[,taxon],
                               id=rownames(otu_taxa))
    
    target_length <- nrow(unique(target_taxon))
    taxon_table <- data.frame(data_mat,taxon=target_taxon$long) %>%
      group_by(taxon) %>%
      summarise_each(funs(sum)) 
    
    tax_short <- substring(str_to_lower(colnames(otu_taxa)[taxon]),1,3)
    taxon_ids <- data.frame(long=taxon_table$taxon,
                            ids=paste0(tax_short,1:nrow(taxon_table)),
                            row.names=paste0(tax_short,1:nrow(taxon_table)))
    taxon_table <- dplyr::select(taxon_table,-taxon)
    
  }else{
    
    if (kegg == TRUE) tax_short <- 'kegg' else tax_short <- 'otu'
    taxon_table <- data_mat
    taxon_ids <- data.frame(long=rownames(taxon_table),
                            ids=paste0(tax_short,1:nrow(taxon_table)),
                            row.names=paste0(tax_short,1:nrow(taxon_table)))
    
  }
  
  taxon_table <- t(taxon_table)
  colnames(taxon_table) <- rownames(taxon_ids)
  
  if (method == 'ra'){
    
    if (filtermin == TRUE) taxon_table <- filter_min(taxon_table,mc,ps)
    if (filtermax == TRUE) taxon_table <- filter_max(taxon_table,pw,pa)
    
    taxon_table <- taxon_table/rowSums(taxon_table)
    
    taxon_table <- taxon_table[rowSums(taxon_table) != 0,colSums(taxon_table) != 0]
    taxon_ids <- taxon_ids[colnames(taxon_table),]  
    
    taxon_ids$long <- as.character(taxon_ids$long)
    return(list(table=taxon_table,ids=taxon_ids))
  }
  if (method == 'rescaled'){
    
    if (filtermin == TRUE) taxon_table <- filter_min(taxon_table,mc,ps)
    if (filtermax == TRUE) taxon_table <- filter_max(taxon_table,pw,pa)
    
    taxon_table <- taxon_table/rowSums(taxon_table)
    multiplier <- 10^ceiling(abs(log(min(taxon_table[taxon_table != 0]),10)))
    
    if (cap == TRUE){
      
      cat('\nBecause multiplier ', round(multiplier,5),' > ',roof, ', capping at ',roof,'\n',sep='')
      multiplier <- roof
      
      n_otu <- ncol(taxon_table)
      taxon_table <- floor(taxon_table * multiplier)
      taxon_table <- taxon_table[,colSums(taxon_table) != 0]
      cat('\nRemoved ',n_otu-ncol(taxon_table),' out of ',n_otu,' OTUs due to capping multiplier.\n',sep='')
      
    }else{
      
      cat('\nMultiplier =', round(multiplier,5),'\n')
      taxon_table <- floor(taxon_table * multiplier)
      
    }
    
    taxon_table <- taxon_table[rowSums(taxon_table) != 0,colSums(taxon_table) != 0]
    taxon_ids <- taxon_ids[colnames(taxon_table),]    
    
    taxon_ids$long <- as.character(taxon_ids$long)
    return(list(table=taxon_table,ids=taxon_ids))
  }
  if (method == 'counts'){
    
    if (filtermin == TRUE) taxon_table <- filter_min(taxon_table,mc,ps)
    if (filtermax == TRUE) taxon_table <- filter_max(taxon_table,pw,pa)
    
    taxon_table <- taxon_table[rowSums(taxon_table) != 0,colSums(taxon_table) != 0]
    
    taxon_ids <- taxon_ids[colnames(taxon_table),]
    
    taxon_ids$long <- as.character(taxon_ids$long)
    return(list(table=taxon_table,ids=taxon_ids))
  }
  
}

filter_max <- function(x,prop_within_subjects=.9,prop_across_subjects=.9,samplerows=TRUE){
  
  if (samplerows == FALSE) x <- t(x)
  
  x_filt <- colSums(x/rowSums(x) > prop_within_subjects) > prop_across_subjects
  cat('\nRemoved ',sum(x_filt),' out of ',length(x_filt),' OTUs due to max constraints.\n',sep='')
  x <- x[,!x_filt]
  
  x <- x[rowSums(x) != 0,]
  
  return(x)
}

filter_min <- function(x,min_count=0,prop_subjects=.01,samplerows=TRUE){
  
  if (samplerows == FALSE) x <- t(x)
  
  x_sel <- colSums(x > min_count) > nrow(x)*prop_subjects
  cat('\nRemoved ',ncol(x)-sum(x_sel),' out of ',length(x_sel),' OTUs due to min constraints.\n',sep='')
  x <- x[,x_sel]
  
  x <- x[rowSums(x) != 0,]
  
  return(x)
}



run_models <- function(metadata,ra,taxinsub=.01,prop=.6,studies=c('Xavier','Morgan','AG'),disease="CD"){
  
  cat('\nStudies: ',paste0(studies,sep=' '),'\n\n',sep='')
  
  ra <- ra[,colSums(ra > 0) > nrow(ra)*taxinsub]
  
  if (length(studies) < 3){
    metadata <- subset(metadata,STUDY %in% studies)
    ra <- ra[rownames(metadata),]
  }
  
  if (disease == "IBD"){
    metadata$DIAGNOSIS[metadata$DIAGNOSIS %in% c("CD","IBD","UC")] <- disease
  }
  
  feature_ibd <- metadata %>% filter(DIAGNOSIS %in% c('Not IBD',disease)) %>% dplyr::select(SampleID,DIAGNOSIS)
  feature_ibd0 <- feature_ibd %>% filter(DIAGNOSIS == 'Not IBD') %>% dplyr::select(SampleID)
  feature_ibd1 <- feature_ibd %>% filter(DIAGNOSIS == disease) %>% dplyr::select(SampleID)
  
  N_samp <- min(nrow(feature_ibd0),nrow(feature_ibd1))
  undersamp <- c(sample(unlist(feature_ibd0),N_samp),sample(unlist(feature_ibd1),N_samp))
  
  ra_us <- ra[undersamp,]
  feature_ibd <- feature_ibd %>% filter(SampleID %in% undersamp)
  
  sets_idx <- caret::createDataPartition(feature_ibd$DIAGNOSIS,times=1,p=prop,list=FALSE)
  train_samps <- unlist(feature_ibd[sets_idx,'SampleID'])
  test_samps <- train_idx <- unlist(feature_ibd[-sets_idx,'SampleID'])
  train_table <- ra_us[train_samps,]
  train_labels <- unlist(feature_ibd %>% filter(SampleID %in% train_samps) %>% dplyr::select(DIAGNOSIS))
  test_table <- ra_us[test_samps,]
  test_labels <- unlist(feature_ibd %>% filter(SampleID %in% test_samps) %>% dplyr::select(DIAGNOSIS))
  
  print(table(train_labels))
  print(table(test_labels))
  
  cat('\n\nPrediction accuracy for test set\n')
  out_rf <- randomForest(x=train_table,y=as.factor(train_labels),
                         xtest=test_table,ytest=as.factor(test_labels),
                         ntree=100,importance=FALSE)
  cat('RF: ',mean(out_rf$test$predicted == test_labels),'\n',sep='')
  
  
  svm_train <- ksvm(train_table, train_labels,
                    type="C-svc",kernel="vanilladot",
                    C=100,scale = FALSE)
  cat('SVM (linear kernel): ',mean(predict(svm_train,test_table,type="response") == test_labels),'\n',sep='')
  
  return(list(rf=out_rf,svm=svm_train))
}

reshape_doc <- function(doc,vocab){
  doc <- doc[doc != 0]
  idx <- match(names(doc),vocab)
  return(rbind(idx,doc))
}

unique_idx <- function(){
  set.seed(as.numeric(Sys.time()))
  sample(1:10^6,1)
}

voc_gen <- function(lookup,vocab,ids,taxon,otu_id=TRUE){
  
  taxon_name <- stringr::str_to_lower(colnames(lookup)[taxon])
  taxon_names <- stringr::str_extract(lookup[ids[vocab,]$long,taxon],'(?<=__)[\\w+.-]+')
  collapse_names <- is.na(taxon_names)
  
  while (any(is.na(taxon_names))){  
    taxon <- taxon - 1
    names_update <- stringr::str_extract(lookup[ids[vocab,]$long,taxon],'(?<=__)[\\w+.-]+')
    taxon_names[is.na(taxon_names)] <- names_update[is.na(taxon_names)]
  }
  
  taxon_names[collapse_names] <- paste(taxon_names[collapse_names],taxon_name)
  if (otu_id) taxon_names <- paste0(taxon_names,' (',vocab,')')
  
  return(taxon_names)
}

print_otus_dx <- function(p,fit,ids,lookup,cov=2,stat=2,content=FALSE){
  
  cis <- p$cis
  means <- unlist(p$means)
  
  if (cov==2){
    topics_ranked <- which(sapply(cis,function(x) all(x>0)))
    topics_ranked <- topics_ranked[order(means[topics_ranked],decreasing=TRUE)]
  }else{
    topics_ranked <- which(sapply(cis,function(x) all(x<0)))
    topics_ranked <- topics_ranked[order(means[topics_ranked],decreasing=FALSE)]
  }
  
  
  if (content){
    
    K <- dim(fit$theta)[2]
    sage <- sageLabels(fit,25)
    
    for (r in topics_ranked){
      cat('\n')
      cat('\nTopic (mean):\tTaxon Family\t(Method:',names(sage$cov.betas[[cov]])[stat],')\n',sep='')
      
      tax <- lookup[as.character(ids[sage$cov.betas[[cov]][[stat]][r,],'long']),-c(1)]
      for (j in 1:nrow(tax)) cat('T',r,' (',round(means[r],4),')',':\t',paste(tax[j,],collapse=' '),'\n',sep='')
    }
    
  }else{
    
    for (r in 1:length(topics_ranked)){
      cat('\n')
      cat('\nTopic (mean):\tTaxon Family\t(Method:',method,')\n',sep='')
      tax <- lookup[as.character(ids[labelTopics(fit)[[method]][topics_ranked[r],],'long']),-c(1)]
      for (j in 1:nrow(tax)) cat('T',topics_ranked[r],' (',round(means[topics_ranked[r]],4),')',':\t',paste(tax[j,],collapse=' '),'\n',sep='')
    }
    
  }
  
}

print_otus_loc <- function(p,fit,ids,lookup,greater=TRUE,content=FALSE,method='frex'){
  
  labels <- p$uvals
  cis <- p$cis
  means <- sapply(p$means,identity)
  
  if (greater){
    topics_ranked <- apply(sapply(cis,function(x) apply(x,2,function(y) all(y>0))),1,which)
    topics_ranked <- sapply(1:length(topics_ranked), function(i) topics_ranked[[i]][order(means[i,topics_ranked[[i]]],decreasing=TRUE)])
  }else{
    topics_ranked <- apply(sapply(cis,function(x) apply(x,2,function(y) all(y<0))),1,which)
    topics_ranked <- sapply(1:length(topics_ranked), function(i) topics_ranked[[i]][order(means[i,topics_ranked[[i]]],decreasing=FALSE)])
  }
  
  if (content){
    
    K <- dim(fit$theta)[2]
    
    for (loc in 1:length(topics_ranked)){
      
      topic_int <- labelTopics(fit)$interaction[(K+1):(K*2),]
      otu_id <- topic_int[topics_ranked[[loc]],]
      
      if (is.null(dim(otu_id))){
        if (length(otu_id)==0){
          next
        }
      }else{
        if (nrow(otu_id)==0){
          next
        }
      }
      
      
      cat('\n')
      cat('\nAssociated with ',as.character(labels[loc]),'\n',sep='')
      cat('\nTopic (mean):\tTaxon Family\t(Method:',method,')\n',sep='')
      
      if (!is.null(dim(otu_id))){
        tax_id <- t(apply(otu_id,1,function(x) as.character(ids[x,'long'])))
      }else{
        tax_id <- matrix(ids[otu_id,'long'],nrow=1)
      }
      
      for (k in 1:nrow(tax_id)){
        
        top <- topics_ranked[[loc]][k]
        tax <- lookup[tax_id[k,],-c(1,7)]
        
        for (tx in 1:nrow(tax)){      
          cat('T',top,' (',round(means[loc,top],4),')',':\t',paste(tax[tx,],collapse=' '),'\n',sep='')
        }
        
      }
      
    }
    
    
  }else{
    
    for (loc in 1:length(topics_ranked)){
      
      otu_id <- labelTopics(fit)[[method]][topics_ranked[[loc]],]
      
      if (is.null(dim(otu_id))){
        if (length(otu_id)==0){
          next
        }
      }else{
        if (nrow(otu_id)==0){
          next
        }
      }
      
      
      cat('\n')
      cat('\nAssociated with ',as.character(labels[loc]),'\n',sep='')
      cat('\nTopic (mean):\tTaxon Family\t(Method:',method,')\n',sep='')
      
      if (!is.null(dim(otu_id))){
        tax_id <- t(apply(otu_id,1,function(x) ids[x,'long']))
      }else{
        tax_id <- matrix(as.character(ids[otu_id,'long']),nrow=1)
      }
      
      for (k in 1:nrow(tax_id)){
        
        top <- topics_ranked[[loc]][k]
        tax <- lookup[tax_id[k,],-c(1,7)]
        
        for (tx in 1:nrow(tax)){      
          cat('T',top,' (',round(means[loc,top],4),')',':\t',paste(tax[tx,],collapse=' '),'\n',sep='')
        }
        
      }
      
    }
    
  }
  
}

rank_topics <- function(p,cov=2){
  
  cis <- p$cis
  means <- unlist(p$means)
  
  if (cov==2){
    topics_ranked <- which(sapply(cis,function(x) all(x>0)))
    topics_ranked <- topics_ranked[order(means[topics_ranked],decreasing=TRUE)]
  }else{
    topics_ranked <- which(sapply(cis,function(x) all(x<0)))
    topics_ranked <- topics_ranked[order(means[topics_ranked],decreasing=FALSE)]
  }
  
  return(topics_ranked)
  
}


print_kegg_dx <- function(p,fit,ids,lookup,cov=2,stat=2,content=FALSE){
  
  cis <- p$cis
  means <- unlist(p$means)
  
  if (cov==2){
    topics_ranked <- which(sapply(cis,function(x) all(x>0)))
    topics_ranked <- topics_ranked[order(means[topics_ranked],decreasing=TRUE)]
  }else{
    topics_ranked <- which(sapply(cis,function(x) all(x<0)))
    topics_ranked <- topics_ranked[order(means[topics_ranked],decreasing=FALSE)]
  }
  
  
  if (content){
    
    K <- dim(fit$theta)[2]
    sage <- sageLabels(fit,25)
    
    for (r in topics_ranked){
      cat('\n')
      cat('\nTopic (mean):\tKEGG Pathway\t(Method:',names(sage$cov.betas[[cov]])[stat],')\n',sep='')
      
      kegg_id <- as.character(ids[sage$cov.betas[[cov]][[stat]][r,],'long'])
      kegg <- lookup[unique(kegg_id)]
      
      for (j in 1:length(kegg)) cat('T',r,' (',round(means[r],4),')\t',
                                    as.character(kegg_id[j]),':\t',paste(kegg[[j]],collapse=' | '),'\n',sep='')
    }
    
  }else{
    
    for (r in 1:length(topics_ranked)){
      cat('\n')
      
      cat('\nTopic (mean):\tKEGG Pathway\t(Method:',method,')\n',sep='')
      kegg_id <- as.character(ids[labelTopics(fit)[[stat]][topics_ranked[r],],'long'])
      kegg <- lookup[kegg_id]
      
      for (j in 1:length(kegg)) cat('T',topics_ranked[r],' (',round(means[topics_ranked[r]],4),')\t',
                                    as.character(kegg_id[j]),':\t',paste(kegg[[j]],collapse=' | '),'\n',sep='')
    }
    
  }
  
}



print_kegg_loc <- function(p,fit,ids,lookup,greater=TRUE,content=FALSE,method='frex'){
  
  labels <- p$uvals
  cis <- p$cis
  means <- sapply(p$means,identity)
  
  if (greater){
    topics_ranked <- apply(sapply(cis,function(x) apply(x,2,function(y) all(y>0))),1,which)
    topics_ranked <- sapply(1:length(topics_ranked), function(i) topics_ranked[[i]][order(means[i,topics_ranked[[i]]],decreasing=TRUE)])
  }else{
    topics_ranked <- apply(sapply(cis,function(x) apply(x,2,function(y) all(y<0))),1,which)
    topics_ranked <- sapply(1:length(topics_ranked), function(i) topics_ranked[[i]][order(means[i,topics_ranked[[i]]],decreasing=FALSE)])
  }
  
  if (content){
    
    K <- dim(fit$theta)[2]
    
    for (loc in 1:length(topics_ranked)){

      topic_int <- labelTopics(fit)$interaction[(K+1):(K*2),]
      otu_id <- topic_int[topics_ranked[[loc]],]
      
      if (is.null(dim(otu_id))){
        if (length(otu_id)==0){
          next
        }
      }else{
        if (nrow(otu_id)==0){
          next
        }
      }
      
      
      cat('\n')
      cat('\nAssociated with ',as.character(labels[loc]),'\n',sep='')
      cat('\nTopic (mean):\tKEGG Pathway\t(Method:',method,')\n',sep='')
      
      if (!is.null(dim(otu_id))){
        tax_id <- t(apply(otu_id,1,function(x) ids[x,'long']))
      }else{
        tax_id <- matrix(ids[otu_id,'long'],nrow=1)
      }
      
      for (k in 1:nrow(tax_id)){
        
        top <- topics_ranked[[loc]][k]
        kegg_id <- as.character(tax_id[k,])
        tax <- lookup[kegg_id]
        
        for (tx in 1:length(tax)){      
          cat('T',top,' (',round(means[loc,top],4),')\t',
              tax_id[k,tx],':\t',paste(tax[[tx]],collapse=' | '),'\n',sep='')
        }
        
      }
      
    }
    
    
  }else{
    
    for (loc in 1:length(topics_ranked)){
      
      otu_id <- labelTopics(fit)[[method]][topics_ranked[[loc]],]
      
      if (is.null(dim(otu_id))){
        if (length(otu_id)==0){
          next
        }
      }else{
        if (nrow(otu_id)==0){
          next
        }
      }
      
      cat('\n')
      cat('\nAssociated with ',as.character(labels[loc]),'\n',sep='')
      cat('\nTopic (mean):\tKEGG Pathway\t(Method:',method,')\n',sep='')
      
      if (!is.null(dim(otu_id))){
        tax_id <- t(apply(otu_id,1,function(x) ids[x,'long']))
      }else{
        tax_id <- matrix(ids[otu_id,'long'],nrow=1)
      }
      
      for (k in 1:nrow(tax_id)){
        
        top <- topics_ranked[[loc]][k]
        kegg_id <- as.character(tax_id[k,])
        tax <- lookup[kegg_id]
        
        for (tx in 1:length(tax)){      
          cat('T',top,' (',round(means[loc,top],4),')',
              tax_id[k,tx],':\t',paste(tax[[tx]],collapse=' | '),'\n',sep='')
        }
        
      }
      
    }
    
  }
  
}


plot_corr_dx <- function(pdx,corrfit,thres=.1){
  cisdx <- pdx$cis
  
  cols <- rep('gray',K)
  cols[which(sapply(cisdx,function(x) all(x>0)))] <- 'red1'
  cols[which(sapply(cisdx,function(x) all(x<0)))] <- 'blue1'
  
  sig_topics <- unique(unlist(apply(corrfit$poscor,1,function(x) which(x > thres & x != 1))))
  
  if (length(sig_topics)==0){
    print('Threshold too large')
  }else{
    plot(corrfit,vertex.color=cols[sig_topics],vertex.label.cex=.8,
         topics=sig_topics,main=labels[label])
  }
}


plot_corr_loc <- function(pdx,ploc,corrfit,label,thres=.1){
  
  cisdx <- pdx$cis
  cisloc <- ploc$cis
  labels <- ploc$uvals
  
  cols <- rep('gray',K)
  
  coldx_pos <- which(sapply(cisdx,function(x) all(x>0)))
  coldx_neg <- which(sapply(cisdx,function(x) all(x<0)))
  colloc <- apply(sapply(cisloc,function(x) apply(x,2,function(y) all(y>0))),1,which)[[label]]
  
  cols[intersect(coldx_pos,colloc)] <- 'orange1'
  cols[intersect(coldx_neg,colloc)] <- 'cyan1'
  cols[colloc[!(colloc %in% c(coldx_pos,coldx_neg))]] <- 'darkolivegreen4'
  cols[coldx_pos[!(coldx_pos %in% colloc)]] <- 'red1'
  cols[coldx_neg[!(coldx_neg %in% colloc)]] <- 'blue1'
  
  sig_topics <- unique(unlist(apply(corrfit$poscor,1,function(x) which(x > thres & x != 1))))
  
  if (length(sig_topics)==0){
    print('Threshold too large')
  }else{
    plot(corrfit,vertex.color=cols[sig_topics],vertex.label.cex=.8,
         topics=sig_topics,main=labels[label])
  }
}

plot_imp <- function(rf, n_imp=25){
  imp_acc_full <- sort(importance(rf)[,'MeanDecreaseAccuracy'],decreasing=TRUE)
  imp_gini_full <- sort(importance(rf)[,'MeanDecreaseGini'],decreasing=TRUE)
  
  imp_acc <- imp_acc_full[1:n_imp]
  imp_gini <- imp_gini_full[1:n_imp]
  
  fig <- rbind(data.frame(val=imp_acc, pw=names(imp_acc), stat='decreased accuracy'),
        data.frame(val=imp_gini,pw=names(imp_gini),stat='gini index')) %>%
    mutate(pw=as.character(pw),
           col=as.factor(ifelse(pw %in% pw[duplicated(pw)],1,0)),
           rank=row_number(),
           pw=ifelse(col==1, paste0(pw,'_',rank), pw)) %>%
    ggplot(aes(val,reorder(pw,val),colour=col)) + geom_point() + facet_wrap(~stat,scales='free') + 
    theme(legend.position='none',
          axis.text.x = element_text(size=20),
          axis.text.y = element_text(size=14),
          axis.title = element_text(size=25,face="bold"),
          strip.text.x = element_text(size = 25),
          title = element_text(size=30,face="bold")) + 
    xlab('Variable Importance') + ylab('') 
  
  print(fig)
  
  return(list(acc=imp_acc_full,gini=imp_gini_full))
}

plot_prox <- function(rf,train_labels){
  mds <- MDSplot(rf,train_labels)
  qplot(mds$points[,1],mds$points[,2],geom='point',colour=train_labels,alpha=.3) +
    labs(colour='') + xlab('Dimension 1') + ylab('Dimension 2') + ggtitle('Proximity') +
    guides(alpha=FALSE) +
    theme(axis.text.x = element_text(size=20),
          axis.text.y = element_text(size=14),
          axis.title = element_text(size=25,face="bold"),
          legend.text = element_text(size=14),
          title = element_text(size=30,face="bold"))
}

print_results <- function(out,dat_type,model,r=5){
  cat('\nBootstrapped (',nrow(out),'-fold) ',
      str_extract(deparse(substitute(out)),'[^_]*$'),
      ' ',dat_type,' prediction results (mean (error)) for test set (n=',
      out[1,1],
      ')\n',
      sep='')
  for (i in 2:ncol(out)){
    cat('\t',colnames(out)[i],':\t',
        sprintf(paste0('%.',r,'f',sep=''),round(mean(out[,i]),r)),'\t(',
        sprintf(paste0('%.',r,'f',sep=''),round(mean(out[,i])-sd(out[,i]),r)),', ',
        sprintf(paste0('%.',r,'f',sep=''),round(mean(out[,i])+sd(out[,i]),r)),')\n',
        sep='')
  }
}

clean <- function(w){
  w_split <- str_to_lower(unlist(str_split(w,' \\| ')))
  w_split <- w_split[w_split != '|']
  w_split <- str_replace(w_split,'\\[[^\\]]*\\]','')
  w_split <- str_trim(tm::removeWords(w_split,c(tm::stopwords('english'),
                                                'unclassified','none','poorly characterized',
                                                'other','others','function unknown',
                                                'general function prediction')),
                      'both')
  w_split <- str_replace_all(w_split,'\\s+','_')
  w_split <- str_replace_all(w_split,'[[:punct:]]+','_')
}

translate_kegg <- function(sage,topic,ids,metadata,cov,stat=2,print=TRUE){
  pw <- metadata[as.character(ids[sage$cov.betas[[cov]][[stat]][topic,],'long'])]
  df <- data.frame(id1=as.vector(sage$cov.betas[[cov]][[stat]][topic,]),
                   id2=ids[as.vector(sage$cov.betas[[cov]][[stat]][topic,]),'long'],
                   pw=sapply(metadata[ids[as.vector(sage$cov.betas[[cov]][[stat]][topic,]),'long']],
                             function(x) paste0(x,collapse=' | ')))
  if (print){
    cat('\nScore is set to ',names(sage$cov.betas[[cov]])[stat],'\n\n',sep='')
    apply(df, 1, function(x) cat(paste0(x,collapse='; '),'\n\n'))
  }else{
    df
  }
}

kegg_pw_collapse2 <- function(kegg,level=0){
  
  out <- list()
  
  if (level==0){
    for (i in 1:length(kegg)){
      out[[kegg[[i]]$id]] <- kegg[[i]]$description
    }
  }else{
    for (i in 1:length(kegg)){
      pw <- kegg[[i]]$pathways
      if (!('None' %in% pw)){
        out[[kegg[[i]]$id]] <- sapply(pw, function(x) x[[level]])
      }else{
        out[[kegg[[i]]$id]] <- 'None'
      }
    }
  }
  
  return(out)
}

kegg_pw_collapse <- function(biom_file,level=0,collapse=FALSE,cog=FALSE){
  
  if (cog & level==3) stop('COG only has 2 levels.')
  
  kegg <- biom_file$rows
  
  out <- list()
  
  if (level==0){
    for (i in 1:length(kegg)){
      if (cog){
        out[[kegg[[i]]$id]] <- kegg[[i]]$metadata$COG_Description
      }else{
        out[[kegg[[i]]$id]] <- kegg[[i]]$metadata$KEGG_Description
      }
    }
  }else{
    for (i in 1:length(kegg)){
      if (cog) pw <- kegg[[i]]$metadata$COG_Category else pw <- kegg[[i]]$metadata$KEGG_Pathways
      if (!('None' %in% pw)){
        if (collapse){
          out[[kegg[[i]]$id]] <- sapply(pw, function(x) paste(x[[level]],collapse=' > '))
        }else{
          out[[kegg[[i]]$id]] <- sapply(pw, function(x) x[[level]])
        }
      }else{
        out[[kegg[[i]]$id]] <- 'None'
      }
    }
  }
  
  return(out)
}

filter_subjects <- function(mat,meta,quant=.05,samp_col=TRUE){
  meta <- data.frame(meta)
  rownames(meta) <- meta$SampleID
  
  if (!samp_col) mat <- t(mat)
  samp_sums <- colSums(mat)
  q_val <- quantile(samp_sums,quant)
  mat <- mat[,samp_sums > q_val]
  mat <- mat[rowSums(mat) != 0,]
  
  meta <- meta[colnames(mat),]
  
  cat('Filtered ',length(samp_sums)-ncol(mat),' out of ',length(samp_sums),' samples with read totals less than ',q_val,'.\n',sep='')
  return(list(mat=mat,meta=meta))
}

kegg_quantile <- function(data,k,lookup,quan=.99){
  topic <- data[,k]
  words <- names(topic)[topic > quantile(topic,quan)]
  return(kegg_metadata_risk[words])  
}

# permTest <- function (formula, stmobj, treatment, nruns = 100, documents, 
#                       vocab, data, seed = NULL, stmverbose = TRUE, uncertainty = "Global") 
# {
#   if (!requireNamespace("clue", quietly = TRUE)) 
#     stop("Install the clue package to use this function")
#   settings <- stmobj$settings
#   if (!is.data.frame(data)) 
#     stop("data object must be a data.frame containing treatment variable.")
#   if (!(treatment %in% colnames(data))) 
#     stop("treatment variable must be in the data.frame")
#   if (!all(data[, treatment] %in% c(0, 1))) 
#     stop("treatment variable must be binary and coded 0 and 1.")
#   prob <- mean(data[, treatment])
#   if (is.null(seed)) {
#     seed <- floor(runif(1) * 10000000)
#     set.seed(seed)
#   }
#   else {
#     set.seed(seed)
#   }
#   settings$verbose <- stmverbose
#   settings$keepHistory <- FALSE
#   betaref <- exp(stmobj$beta$logbeta[[1]])
#   qeffects <- function(formula, mod, data, uncertainty) {
#     prep <- stm::estimateEffect(formula, stmobj = mod, metadata = data, 
#                                 uncertainty = uncertainty)
#     betas <- lapply(prep$parameters, function(x) lapply(x, 
#                                                         function(y) stm:::rmvnorm(100, y$est, y$vcov)))
#     simbetas <- lapply(betas, do.call, what = rbind)
#     parcol <- which(colnames(settings$covariates$X) == treatment)
#     par <- lapply(simbetas, function(x) quantile(x[, parcol], 
#                                                  c(0.025, 0.5, 0.975)))
#     par
#   }
#   cat("Calculating effects for reference model \n")
#   ref.effect <- qeffects(formula, stmobj, data, uncertainty)
#   tosave <- list()
#   for (i in 1:(nruns - 1)) {
#     settings$seed <- floor(runif(1) * 10000000)
#     data[, treatment] <- rbinom(n = nrow(data), size = 1, 
#                                 prob = prob)
#     termobj <- terms(formula, data = data)
#     if (attr(termobj, "response") == 1) 
#       stop("Response variables should not be included in prevalence formula.")
#     settings$covariates$X <- model.matrix(termobj, data = data)
#     cat(sprintf("Running model %i of %i \n", (i + 1), (nruns)))
#     mod <- stm:::stm.control(documents, vocab, settings = settings, model=NULL)
#     par <- qeffects(formula, mod, data, uncertainty)
#     betamod <- exp(mod$beta$logbeta[[1]])
#     align <- clue::solve_LSAP(betaref %*% t(betamod), maximum = TRUE)
#     tosave[[i]] <- par[align]
#   }
#   out <- list(ref = ref.effect, permute = tosave, variable = treatment, 
#               seed = seed)
#   class(out) <- "STMpermute"
#   return(out)
# }

permTest <- function (formula, stmobj, treatment, nruns = 100, documents, 
                      vocab, data, seed = NULL, stmverbose = TRUE, uncertainty = "Global") 
{
  if (is.null(stmobj$settings$covariates$X)){
    stmobj$settings$covariates$X <- cbind(1,data[,treatment])
    colnames(stmobj$settings$covariates$X) <- c('(Intercept)',treatment)
  }
  if (!requireNamespace("clue", quietly = TRUE)) 
    stop("Install the clue package to use this function")
  settings <- stmobj$settings
  if (!is.data.frame(data)) 
    stop("data object must be a data.frame containing treatment variable.")
  if (!(treatment %in% colnames(data))) 
    stop("treatment variable must be in the data.frame")
  if (!all(data[, treatment] %in% c(0, 1))) 
    stop("treatment variable must be binary and coded 0 and 1.")
  prob <- mean(data[, treatment])
  if (is.null(seed)) {
    seed <- floor(runif(1) * 10000000)
    set.seed(seed)
  }
  else {
    set.seed(seed)
  }
  
  settings$verbose <- stmverbose
  settings$keepHistory <- FALSE
  betaref <- exp(stmobj$beta$logbeta[[1]])
  qeffects <- function(formula, mod, data, uncertainty) {
    prep <- stm::estimateEffect(formula, stmobj = mod, metadata = data, 
                                uncertainty = uncertainty)
    betas <- lapply(prep$parameters, function(x) lapply(x, 
                                                        function(y) stm:::rmvnorm(100, y$est, y$vcov)))
    simbetas <- lapply(betas, do.call, what = rbind)
    parcol <- which(colnames(settings$covariates$X) == treatment)
    par <- lapply(simbetas, function(x) quantile(x[, parcol], 
                                                 c(0.025, 0.5, 0.975)))
    par
  }
  cat("Calculating effects for reference model \n")
  ref.effect <- qeffects(formula, stmobj, data, uncertainty)
  tosave <- list()
  for (i in seq_len(nruns - 1)) {
    settings$seed <- floor(runif(1) * 10000000)
    data[, treatment] <- rbinom(n = nrow(data), size = 1, prob = prob)
    termobj <- terms(formula, data = data)
    if (attr(termobj, "response") == 1) 
      stop("Response variables should not be included in prevalence formula.")
    settings$covariates$X <- model.matrix(termobj, data = data)
    cat(sprintf("Running model %i of %i \n", (i + 1), (nruns)))
    mod <- stm:::stm.control(documents, vocab, settings = settings, model=NULL)
    par <- qeffects(formula, mod, data, uncertainty)
    betamod <- exp(mod$beta$logbeta[[1]])
    align <- clue::solve_LSAP(betaref %*% t(betamod), maximum = TRUE)
    tosave[[i]] <- par[align]
  }
  out <- list(ref = ref.effect, permute = tosave, variable = treatment, seed = seed)
  class(out) <- "STMpermute"
  return(out)
}

kegg_metadata <- function(biom_file, level){
  
  kegg <- biom_file$rows
  
  out <- vector(mode='list',length=length(kegg))
  for (i in seq_along(kegg)){
    
    kegg_id <- kegg[[i]]$id
    
    kegg_de <- kegg[[i]]$metadata$KEGG_Description
    kegg_de <- str_to_lower(gsub('\\["(.*)"\\]','\\1',kegg_de))
    
    kegg_pw <- kegg[[i]]$metadata$KEGG_Pathways
    kegg_pw <- strsplit(kegg_pw,'\\], \\[')[[1]]
    kegg_pw <- strsplit(kegg_pw,', ')
    
    out_pw <- vector(mode='character',length=length(kegg_pw))
    for (j in seq_along(kegg_pw)){
      out_pw[[j]] <- str_to_lower(gsub('[[:punct:]]','',kegg_pw[[j]]))[level]
    }
    
    out[[i]] <- list(id=kegg_id,
                     description=kegg_de,
                     pathway=out_pw)
    
    names(out)[i] <- kegg_id
    
  }
  
  return(out)
}