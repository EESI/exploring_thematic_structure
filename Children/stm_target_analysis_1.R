rds_out <- readRDS('~/Qiime/Children/xavier_data_loading.rds')
otu_lookup <- rds_out$otu_lookup
data_mat_xavier <- rds_out$data_mat_xavier
metadata_xavier <- rds_out$metadata_xavier
kegg_metadata_xavier <- rds_out$kegg_metadata_xavier
kegg_mat_xavier <- rds_out$kegg_mat_xavier


counts <- collapse_table(otu_lookup,kegg_mat_xavier,metadata_xavier,
                         taxon=NULL,method='rescaled',
                         kegg=TRUE, cap=TRUE,
                         filtermin=FALSE,mc=5,ps=.1,
                         filtermax=FALSE,pw=.9,pa=.9)

counts$meta <- metadata_xavier[metadata_xavier$SampleID %in% rownames(counts$table),]
counts$meta_clean <- counts$meta[!is.na(counts$meta$ISOLATION_SOURCE),]
counts$table_clean <- counts$table[counts$meta_clean$SampleID,]

vocab <- as.character(counts$ids$ids)
docs <- lapply(1:nrow(counts$table_clean),function(i) reshape_doc(counts$table_clean[i,],vocab))
meta <- counts$meta_clean

prep <- prepDocuments(docs,vocab,meta,lower.thresh=0) ### removes terms apparently

f <- 'stm_kegg_cohort_content_pooled_binary_l1_k75_374728.rds'
fit <- readRDS('~/Qiime/Children/stm_kegg_cohort_content_pooled_binary_l1_k75_374728.rds')

###
###
###
uidx <- unique_idx()
dump_dir <- paste0('~/Qiime/Children/stm_out_',uidx,'/stm_kegg',sep='')
dir.create(dump_dir, showWarnings = FALSE,recursive=TRUE)
sink(file.path(dump_dir, 'stm.txt'))

cat(rep('\n',8),f,rep('\n',3))

cat('\nInteractive Topics:\n\t\thttps://dl.dropboxusercontent.com/u/907375/STM/',str_extract(f,'[^\\.]*'),'/index.html\n',sep='')

content <- !is.na(str_extract(f,'content'))
K <- dim(fit$theta)[2]

cat('\nDispersion (',K,' topics): ',checkResiduals(fit,prep$documents)$dispersion,'\n',sep='')

eff <- estimateEffect(1:K ~ as.factor(DIAGNOSIS) + as.factor(ISOLATION_SOURCE) + as.factor(COHORT), 
                      fit,
                      meta = prep$meta,
                      uncertainty = 'Global')

method <- 'frex'
eff_plot_dx <- plot(eff,'DIAGNOSIS',method='difference',cov.value1='CD',cov.value2='Not IBD',width=50,verbose.labels=0,model='frex')
cat('\n\n\nAssociated with CD\n\n')
print_kegg_dx(eff_plot_dx,fit,counts$ids,kegg_metadata_xavier,greater=TRUE,content=content,method=method)
cat('\n\n\nAssociated with Healthy\n\n')
print_kegg_dx(eff_plot_dx,fit,counts$ids,kegg_metadata_xavier,greater=FALSE,content=content,method=method)

eff_plot_loc <- plot(eff,'ISOLATION_SOURCE',method="pointestimate",width=50,model='frex',n=10)
print_kegg_loc(eff_plot_loc,fit,counts$ids,kegg_metadata_xavier,greater=TRUE,content=content,method=method)

sink()
###
###
###

for (i in 1:length(unique(prep$meta$ISOLATION_SOURCE))){
  
  corrfit <- topicCorr(fit,cutoff=.1)
  
  file_name <- paste0('corr_loc_',i,'.pdf')
  pdf(file.path(dump_dir, file_name),width=12,height=12)
  plot_corr_loc(eff_plot_dx,eff_plot_loc,corrfit,label=i)
  dev.off()
  
}


theta <- t(apply(fit$theta + 1/K,1,function(x) x/sum(x)))
phi <- t(apply(exp(fit$beta$logbeta[[1]]) + .001,1,function(x) x/sum(x)))
freq <- colSums(counts$table_clean[,prep$vocab])
len <- sapply(prep$documents,function(x) sum(x[2,]))
voc <- as.vector(sapply(kegg_metadata_xavier[counts$id[prep$vocab,]$long],function(x) x[1]))
voc[nchar(voc)==0] <- 'None'

json <- createJSON(phi=phi,
                   theta=theta,
                   doc.length=len,
                   vocab=voc,
                   term.frequency=freq,
                   reorder.topics=FALSE,
                   R=15)

serVis(json, out.dir = file.path(dump_dir,'vis'), open.browser=FALSE)


rank_size <- 50
ranked_vocab <- apply(phi[apply(theta,1,which.max),],1,
                      function(x) paste0(voc[order(x,decreasing=TRUE)[1:rank_size]],collapse='  _____  '))
prep$meta$text <- ranked_vocab

dir.create(file.path(dump_dir,'browser'), showWarnings = FALSE)
stmBrowser(fit,data=prep$meta,covariates=c('DIAGNOSIS','ISOLATION_SOURCE'),
           labeltype='topics',
           text='text',
           n=length(ranked_vocab),
           directory=file.path(dump_dir,'browser'))

dir.create(file.path(dump_dir,'viz'), showWarnings = FALSE)
stmCorrViz(fit,file_out=file.path(dump_dir,'viz/vizz.html'),
           prep$meta$text,prep$documents)










counts <- collapse_table(otu_lookup,data_mat_xavier,metadata_xavier,
                         taxon=NULL,method='rescaled',
                         kegg=FALSE, cap=TRUE,
                         filtermin=FALSE,mc=5,ps=.1,
                         filtermax=FALSE,pw=.9,pa=.9)

counts$meta <- metadata_xavier[metadata_xavier$SampleID %in% rownames(counts$table),]
counts$meta_clean <- counts$meta[!is.na(counts$meta$ISOLATION_SOURCE),]
counts$table_clean <- counts$table[counts$meta_clean$SampleID,]

vocab <- as.character(counts$ids$ids)
docs <- lapply(1:nrow(counts$table_clean),function(i) reshape_doc(counts$table_clean[i,],vocab))
meta <- counts$meta_clean

prep <- prepDocuments(docs,vocab,meta,lower.thresh=0) ### removes terms apparently

f <- 'stm_otus_cohort_content_pooled_binary_l1_k75_611210.rds'
fit <- readRDS("~/Qiime/Children/stm_otus_cohort_content_pooled_binary_l1_k75_611210.rds")

###
###
###
dump_dir <- paste0('~/Qiime/Children/stm_out_',uidx,'/stm_otu',sep='')
dir.create(dump_dir, showWarnings = FALSE,recursive=TRUE)
sink(file.path(dump_dir, 'stm.txt'))

cat(rep('\n',8),f,rep('\n',3))

cat('\nInteractive Topics:\n\t\thttps://dl.dropboxusercontent.com/u/907375/STM/',str_extract(f,'[^\\.]*'),'/index.html\n',sep='')

content <- !is.na(str_extract(f,'content'))
K <- dim(fit$theta)[2]

cat('\nDispersion (',K,' topics): ',checkResiduals(fit,prep$documents)$dispersion,'\n',sep='')

eff <- estimateEffect(1:K ~ as.factor(DIAGNOSIS) + as.factor(ISOLATION_SOURCE) + as.factor(COHORT), 
                      fit,
                      meta = prep$meta,
                      uncertainty = 'Global')

method <- 'frex'
eff_plot_dx <- plot(eff,'DIAGNOSIS',method='difference',cov.value1='CD',cov.value2='Not IBD',width=50,verbose.labels=0,model='frex')
cat('\n\n\nAssociated with CD\n\n')
print_otus_dx(eff_plot_dx,fit,counts$ids,otu_lookup,greater=TRUE,content=content,method=method)
cat('\n\n\nAssociated with Healthy\n\n')
print_otus_dx(eff_plot_dx,fit,counts$ids,otu_lookup,greater=FALSE,content=content,method=method)

eff_plot_loc <- plot(eff,'ISOLATION_SOURCE',method="pointestimate",width=50,model='frex',n=10)
print_otus_loc(eff_plot_loc,fit,counts$ids,otu_lookup,greater=TRUE,content=content,method=method)

sink()
###
###
###

for (i in 1:length(unique(prep$meta$ISOLATION_SOURCE))){
  
  corrfit <- topicCorr(fit,cutoff=.1)
  
  file_name <- paste0('corr_loc_',i,'.pdf')
  pdf(file.path(dump_dir, file_name),width=12,height=12)
  plot_corr_loc(eff_plot_dx,eff_plot_loc,corrfit,label=i)
  dev.off()
  
}


theta <- t(apply(fit$theta + 1/K,1,function(x) x/sum(x)))
phi <- t(apply(exp(fit$beta$logbeta[[1]]) + .001,1,function(x) x/sum(x)))
freq <- colSums(counts$table_clean[,prep$vocab])
len <- sapply(prep$documents,function(x) sum(x[2,]))
voc <- voc_gen(otu_lookup,prep$vocab,counts$ids,5)

json <- createJSON(phi=phi,
                   theta=theta,
                   doc.length=len,
                   vocab=voc,
                   term.frequency=freq,
                   reorder.topics=FALSE,
                   R=15)

serVis(json, out.dir = file.path(dump_dir,'vis'), open.browser=FALSE)


rank_size <- 50
ranked_vocab <- apply(phi[apply(theta,1,which.max),],1,
                      function(x) paste0(voc[order(x,decreasing=TRUE)[1:rank_size]],collapse='  _____  '))
prep$meta$text <- ranked_vocab

dir.create(file.path(dump_dir,'browser'), showWarnings = FALSE)
stmBrowser(fit,data=prep$meta,covariates=c('DIAGNOSIS','ISOLATION_SOURCE'),
           labeltype='topics',
           text='text',
           n=length(ranked_vocab),
           directory=file.path(dump_dir,'browser'))

dir.create(file.path(dump_dir,'viz'), showWarnings = FALSE)
stmCorrViz(fit,file_out=file.path(dump_dir,'viz/vizz.html'),
           prep$meta$text,prep$documents)















































rds_out <- readRDS('~/Qiime/Children/xavier_data_loading.rds')
id_oversamp <- readRDS('~/Qiime/Children/xavier_oversamp_loading')
otu_lookup <- rds_out$otu_lookup
data_mat_xavier <- rds_out$data_mat_xavier
metadata_xavier <- rds_out$metadata_xavier
kegg_metadata_xavier <- rds_out$kegg_metadata_xavier
kegg_mat_xavier <- rds_out$kegg_mat_xavier


counts <- collapse_table(otu_lookup,kegg_mat_xavier,metadata_xavier,
                         taxon=NULL,method='rescaled',
                         kegg=TRUE, cap=TRUE,roof=10^5,
                         filtermin=FALSE,mc=5,ps=.1,
                         filtermax=FALSE,pw=.9,pa=.9)

counts$meta <- metadata_xavier[metadata_xavier$SampleID %in% rownames(counts$table),]
counts$meta_clean <- counts$meta[!is.na(counts$meta$ISOLATION_SOURCE),]
counts$table_clean <- counts$table[counts$meta_clean$SampleID,]
rownames(counts$meta_clean) <- counts$meta_clean$SampleID

counts$meta_clean <- rbind(counts$meta_clean,counts$meta_clean[id_oversamp,])
counts$table_clean <- rbind(counts$table_clean,counts$table_clean[id_oversamp,])

vocab <- as.character(counts$ids$ids)
docs <- lapply(1:nrow(counts$table_clean),function(i) reshape_doc(counts$table_clean[i,],vocab))
meta <- counts$meta_clean
meta$DIAGNOSIS <- as.factor(ifelse(meta$DIAGNOSIS=='CD',1,0))
meta$ISOLATION_SOURCE <- as.factor(meta$ISOLATION_SOURCE)
meta$COHORT <- as.factor(meta$COHORT)

prep <- prepDocuments(docs,vocab,meta,lower.thresh=0) ### removes terms apparently

f <- 'stm_oversamp_kegg_cohort_content_pooled_binary_l1_k75_816225.rds'
fit <- readRDS('~/Qiime/Children/stm_oversamp_kegg_cohort_content_pooled_binary_l1_k75_816225.rds')

###
###
###
uidx <- unique_idx()
dump_dir <- paste0('~/Qiime/Children/stm_out_oversamp_',uidx,'/stm_kegg',sep='')
dir.create(dump_dir, showWarnings = FALSE,recursive=TRUE)
sink(file.path(dump_dir, 'stm.txt'))

cat(rep('\n',8),f,rep('\n',3))

cat('\nInteractive Topics:\n\t\thttps://dl.dropboxusercontent.com/u/907375/STM/',str_extract(f,'[^\\.]*'),'/index.html\n',sep='')

content <- !is.na(str_extract(f,'content'))
K <- dim(fit$theta)[2]

cat('\nDispersion (',K,' topics): ',checkResiduals(fit,prep$documents)$dispersion,'\n',sep='')

eff <- estimateEffect(1:K ~ as.factor(ISOLATION_SOURCE) + as.factor(DIAGNOSIS) + as.factor(COHORT), 
                      fit,
                      meta = prep$meta,
                      uncertainty = 'Global')

method <- 'frex'
eff_plot_dx <- plot(eff,'DIAGNOSIS',method='difference',cov.value1='1',cov.value2='0',width=50,verbose.labels=0,model='frex')
cat('\n\n\nAssociated with CD\n\n')
print_kegg_dx(eff_plot_dx,fit,counts$ids,kegg_metadata_xavier,greater=TRUE,content=content,method=method)
cat('\n\n\nAssociated with Healthy\n\n')
print_kegg_dx(eff_plot_dx,fit,counts$ids,kegg_metadata_xavier,greater=FALSE,content=content,method=method)


eff_plot_loc <- plot_effect(eff,'ISOLATION_SOURCE',method="pointestimate",width=50,model='frex',n=10)
print_kegg_loc(eff_plot_loc,fit,counts$ids,kegg_metadata_xavier,greater=TRUE,content=content,method=method)

sink()
###
###
###

for (i in 1:length(unique(prep$meta$ISOLATION_SOURCE))){
  
  corrfit <- topicCorr(fit,cutoff=.1)
  
  file_name <- paste0('corr_loc_',i,'.pdf')
  pdf(file.path(dump_dir, file_name),width=12,height=12)
  plot_corr_loc(eff_plot_dx,eff_plot_loc,corrfit,label=i)
  dev.off()
  
}


theta <- t(apply(fit$theta + 1/K,1,function(x) x/sum(x)))
phi <- t(apply(exp(fit$beta$logbeta[[1]]) + .001,1,function(x) x/sum(x)))
freq <- colSums(counts$table_clean[,prep$vocab])
len <- sapply(prep$documents,function(x) sum(x[2,]))
voc <- as.vector(sapply(kegg_metadata_xavier[counts$id[prep$vocab,]$long],function(x) x[1]))
voc[nchar(voc)==0] <- 'None'

json <- createJSON(phi=phi,
                   theta=theta,
                   doc.length=len,
                   vocab=voc,
                   term.frequency=freq,
                   reorder.topics=FALSE,
                   R=15)

serVis(json, out.dir = file.path(dump_dir,'vis'), open.browser=FALSE)


rank_size <- 50
ranked_vocab <- apply(phi[apply(theta,1,which.max),],1,
                      function(x) paste0(voc[order(x,decreasing=TRUE)[1:rank_size]],collapse='  _____  '))
prep$meta$text <- ranked_vocab

dir.create(file.path(dump_dir,'browser'), showWarnings = FALSE)
stmBrowser(fit,data=prep$meta,covariates=c('DIAGNOSIS','ISOLATION_SOURCE'),
           labeltype='topics',
           text='text',
           n=length(ranked_vocab),
           directory=file.path(dump_dir,'browser'))

dir.create(file.path(dump_dir,'viz'), showWarnings = FALSE)
stmCorrViz(fit,file_out=file.path(dump_dir,'viz/vizz.html'),
           prep$meta$text,prep$documents)










counts <- collapse_table(otu_lookup,data_mat_xavier,metadata_xavier,
                        taxon=NULL,method='rescaled',
                        kegg=FALSE, cap=TRUE,roof=10^5,
                        filtermin=FALSE,mc=5,ps=.01,
                        filtermax=FALSE,pw=.95,pa=.9)

counts$meta <- metadata_xavier[metadata_xavier$SampleID %in% rownames(counts$table),]
counts$meta_clean <- counts$meta[!is.na(counts$meta$ISOLATION_SOURCE),]
counts$table_clean <- counts$table[counts$meta_clean$SampleID,]
rownames(counts$meta_clean) <- counts$meta_clean$SampleID

counts$meta_clean <- rbind(counts$meta_clean,counts$meta_clean[id_oversamp,])
counts$table_clean <- rbind(counts$table_clean,counts$table_clean[id_oversamp,])

vocab <- as.character(counts$ids$ids)
docs <- lapply(1:nrow(counts$table_clean),function(i) reshape_doc(counts$table_clean[i,],vocab))
meta <- counts$meta_clean
meta$DIAGNOSIS <- as.factor(ifelse(meta$DIAGNOSIS=='CD',1,0))
meta$ISOLATION_SOURCE <- as.factor(meta$ISOLATION_SOURCE)
meta$COHORT <- as.factor(meta$COHORT)

prep <- prepDocuments(docs,vocab,meta,lower.thresh=0) ### removes terms apparently

f <- 'stm_oversamp_otus_cohort_content_pooled_binary_l1_k75_314644.rds'
fit <- readRDS("~/Qiime/Children/stm_oversamp_otus_cohort_content_pooled_binary_l1_k75_314644.rds")

###
###
###
dump_dir <- paste0('~/Qiime/Children/stm_out_oversamp_',uidx,'/stm_otu',sep='')
dir.create(dump_dir, showWarnings = FALSE,recursive=TRUE)
sink(file.path(dump_dir, 'stm.txt'))

cat(rep('\n',8),f,rep('\n',3))

cat('\nInteractive Topics:\n\t\thttps://dl.dropboxusercontent.com/u/907375/STM/',str_extract(f,'[^\\.]*'),'/index.html\n',sep='')

content <- !is.na(str_extract(f,'content'))
K <- dim(fit$theta)[2]

cat('\nDispersion (',K,' topics): ',checkResiduals(fit,prep$documents)$dispersion,'\n',sep='')

eff <- estimateEffect(1:K ~ as.factor(ISOLATION_SOURCE) + as.factor(DIAGNOSIS) + as.factor(COHORT), 
                      fit,
                      meta = prep$meta,
                      uncertainty = 'Global')

method <- 'frex'
eff_plot_dx <- plot(eff,'DIAGNOSIS',method='difference',cov.value1='1',cov.value2='0',width=50,verbose.labels=0,model='frex')
cat('\n\n\nAssociated with CD\n\n')
print_otus_dx(eff_plot_dx,fit,counts$ids,otu_lookup,greater=TRUE,content=content,method=method)
cat('\n\n\nAssociated with Healthy\n\n')
print_otus_dx(eff_plot_dx,fit,counts$ids,otu_lookup,greater=FALSE,content=content,method=method)

eff_plot_loc <- plot(eff,'ISOLATION_SOURCE',method="pointestimate",width=50,model='frex',n=10)
print_otus_loc(eff_plot_loc,fit,counts$ids,otu_lookup,greater=TRUE,content=content,method=method)

sink()
###
###
###

for (i in 1:length(unique(prep$meta$ISOLATION_SOURCE))){
  
  corrfit <- topicCorr(fit,cutoff=.1)
  
  file_name <- paste0('corr_loc_',i,'.pdf')
  pdf(file.path(dump_dir, file_name),width=12,height=12)
  plot_corr_loc(eff_plot_dx,eff_plot_loc,corrfit,label=i)
  dev.off()
  
}


theta <- t(apply(fit$theta + 1/K,1,function(x) x/sum(x)))
phi <- t(apply(exp(fit$beta$logbeta[[1]]) + .001,1,function(x) x/sum(x)))
freq <- colSums(counts$table_clean[,prep$vocab])
len <- sapply(prep$documents,function(x) sum(x[2,]))
voc <- voc_gen(otu_lookup,prep$vocab,counts$ids,5)

json <- createJSON(phi=phi,
                   theta=theta,
                   doc.length=len,
                   vocab=voc,
                   term.frequency=freq,
                   reorder.topics=FALSE,
                   R=15)

serVis(json, out.dir = file.path(dump_dir,'vis'), open.browser=FALSE)


rank_size <- 50
ranked_vocab <- apply(phi[apply(theta,1,which.max),],1,
                      function(x) paste0(voc[order(x,decreasing=TRUE)[1:rank_size]],collapse='  _____  '))
prep$meta$text <- ranked_vocab

dir.create(file.path(dump_dir,'browser'), showWarnings = FALSE)
stmBrowser(fit,data=prep$meta,covariates=c('DIAGNOSIS','ISOLATION_SOURCE'),
           labeltype='topics',
           text='text',
           n=length(ranked_vocab),
           directory=file.path(dump_dir,'browser'))

dir.create(file.path(dump_dir,'viz'), showWarnings = FALSE)
stmCorrViz(fit,file_out=file.path(dump_dir,'viz/vizz.html'),
           prep$meta$text,prep$documents)
























