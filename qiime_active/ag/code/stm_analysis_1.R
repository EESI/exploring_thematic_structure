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
library(viridis)


source('~/Dropbox/stm_microbiome/code_active/stm_functions.R')
source('~/Dropbox/stm_microbiome/code_active/nav_froz_fxns_3.R')


K <-  25 #50
fit_idx <- 1

models_dir <- sprintf('~/Dropbox/stm_microbiome/qiime_active/ag/models/%s',K)


fit_fn <- file.path(models_dir,sprintf('stm_k_%s.rds',K))
ko_fn <- file.path(models_dir,sprintf('ko_k_%s_fit%s.biom',K,fit_idx))
cog_fn <- file.path(models_dir,sprintf('cog_k_%s_fit%s.biom',K,fit_idx))


dat <- readRDS(fit_fn)

OTU <- dat$counts$table_clean
META <- dat$meta %>%
  mutate(age=as.integer(META$age_corrected))
VOCAB <- dat$vocab

TAX <- as.data.frame(dat$taxa) %>% 
  mutate(long=rownames(dat$taxa)) %>%
  left_join(dat$counts$ids,by='long')

fit <- dat$fits[[fit_idx+1]]

theta <- fit$theta
rownames(theta) <- META$SampleID
colnames(theta) <- paste0('T',1:K)

beta <- exp(t(fit$beta$logbeta[[1]]))
rownames(beta) <- VOCAB
colnames(beta) <- paste0('T',1:K)


eff <- estimateEffect(1:K ~ diet, fit, META, uncertainty='Global')
eff_plot <- plot.estimateEffect(eff, 'diet', model=fit, 
                                topics=1:K, method='difference',cov.value1=1,cov.value2=0)
topic_sig <- which(sapply(eff_plot$cis,function(x) sum(sign(x))) != 0)

coefs <- data.frame(coef=unlist(eff_plot$means),
                    topic=paste0('T',1:K),
                    stringsAsFactors=FALSE) %>%
  mutate(coefn=coef/sd(coef))

data.frame(log(theta),SampleID=rownames(theta)) %>%
  gather(topic,p,-SampleID) %>%
  mutate(p=ifelse(p < -20,-20,p)) %>%
  left_join(META,by='SampleID') %>%
  left_join(coefs,by='topic') %>%
  arrange(coef,desc(diet),desc(age)) %>%
  mutate(diet=as.factor(diet),
         topic=factor(topic,levels=unique(topic)),
         SampleID=factor(SampleID,levels=unique(SampleID))) %>%
  ggplot(aes(x=SampleID,y=topic,fill=p)) +
  geom_raster() +
  geom_rug(aes(color=diet),sides='b',size=.05,linetype=1,alpha=.3) +
  geom_vline(xintercept=sum(META$diet==1)+.5,colour='white',size=.4) +
  geom_hline(yintercept=sum(coefs$coef>0)+1.5,colour='white',size=.4) +
  scale_fill_viridis(name='p') +
  theme(axis.text.x = element_blank(),
        axis.ticks = element_blank()) +
  guides(color=FALSE)


top_taxa <- as.data.frame(OTU) %>%
  gather(ids,count) %>%
  left_join(TAX,by='ids') %>%
  group_by(Phylum) %>%
  summarise(count=sum(count)) %>%
  arrange(desc(count)) %>%
  select(Phylum) %>%
  unlist() %>%
  as.character()

data.frame(log(beta),ids=rownames(beta)) %>%
  gather(topic,p,-ids) %>%
  mutate(p=ifelse(p < -20,-20,p)) %>%
  left_join(TAX,by='ids') %>%
  left_join(coefs,by='topic') %>%
  mutate(Phylum=as.character(Phylum),
         Phylum=ifelse(Phylum %in% top_taxa[1:5],Phylum,'Other')) %>%
  arrange(coef,Phylum,Class,Order,Family,Genus,Species) %>%
  mutate(topic=factor(topic,levels=unique(topic)),
         ids=factor(ids,levels=unique(ids))) %>%
  ggplot(aes(x=ids,y=topic,fill=p)) +
  geom_raster() +
  facet_wrap(~Phylum,scales='free_x') +
  geom_hline(yintercept=sum(coefs$coef>0)+1.5,colour='white',size=.4) +
  scale_fill_viridis(name='p') +
  theme(axis.text.x = element_blank(),
        axis.ticks = element_blank()) 






biom_file <- read_biom(ko_fn)
ko <- as.matrix(biom_data(biom_file))
ko <- ko[rowSums(ko)>0,]
kegg_meta <- kegg_metadata(biom_file,2)[rownames(ko)]

my_kegg_order <- c('carbohydrate metabolism',
                   'amino acid metabolism',
                   'lipid metabolism',
                   'cell motility',
                   'glycan biosynthesis and metabolism',
                   'xenobiotics biodegradation and metabolism',
                   'metabolism of cofactors and vitamins',
                   'nucleotide metabolism',
                   'metabolism of terpenoids and polyketides',
                   'biosynthesis of other secondary metabolites',
                   'energy metabolism',
                   'metabolism',
                   'membrane transport',
                   'cellular processes and signaling',
                   'genetic information processing',
                   'unknown or misclassification')

ko_colsums <- colSums(ko)
target_pws <- my_kegg_order[-c(16)] 
target_kos <- lapply(target_pws,function(pw) names(which(sapply(kegg_meta, function(x) any(x$pathway %in% pw)))))
names(target_kos) <- target_pws


j1 <- 1
j2 <- 0
gene_mat <- as.data.frame(matrix(0,sum(sapply(target_kos,length)) * K,6,
                                 dimnames=list(NULL,c('topic','pw','ko','count','offset','coef'))))
for (i in seq_along(target_kos)){
  pw <- target_kos[[i]]
  for (k in seq_len(K)){
    j2 <- j2 + length(pw)
    gene_mat[j1:j2,1] <- k
    gene_mat[j1:j2,2] <- target_pws[i]
    gene_mat[j1:j2,3] <- pw
    gene_mat[j1:j2,4] <- ko[pw,k]
    gene_mat[j1:j2,5] <- ko_colsums[k]
    gene_mat[j1:j2,6] <- coefs$coef[k]
    j1 <- j2+1
  }
}

gene_mat <- gene_mat %>%
  mutate(coef=(coef-mean(coef))/sd(coef),
         log_count=log(count+1),
         ra=count/offset,
         logoffset=log(offset))

mm1 <- lmer(log_count ~  (1|pw) + (1|topic),data=gene_mat)
mm2 <- lmer(log_count ~  (1|pw) + (1|topic) + (1|topic:pw),data=gene_mat)
mm3 <- lmer(log_count ~  coef + (1|topic) + (coef|pw),data=gene_mat)

tmm2 <- lmer(log_count ~  (1|pw) + (1|topic) + (1|topic:pw),data=gene_mat,REML=FALSE)
tmm1 <- lmer(log_count ~  (1|pw) + (1|topic),data=gene_mat,REML=FALSE)

anova(tmm2,tmm1)

rmm2 <- ranef(mm2)


dotplot(ranef(mm2, condVar = TRUE), strip = FALSE)$topic
dotplot(ranef(mm2, condVar = TRUE), strip = FALSE)$pw
dotplot(ranef(mm2, condVar = TRUE), strip = FALSE)$`topic:pw`


int_coef <- rmm2$`topic:pw`
int_coef_names <- str_split(rownames(int_coef),':')
int_mat <- matrix(0,length(target_pws),K,dimnames=list(target_pws,1:K))
for (i in seq_len(nrow(int_coef))){
  k <- int_coef_names[[i]][1]
  pw <- int_coef_names[[i]][2]
  int_mat[pw,k] <- int_coef[i,1]
}

data.frame(int_mat,pw=rownames(int_mat)) %>%
  gather(topic,weight,-pw) %>%
  ggplot(aes(topic,pw,fill=weight)) +
  geom_raster() +
  scale_fill_viridis(name='weight')





library(rstanarm) 
# 
# 
# stan1 <- rstanarm::stan_glmer.nb(count ~ (1|pw) + (1|topic) + (1|topic:pw),
#                                  link='log',
#                                  data=gene_mat,
#                                  chains=4,cores=4,
#                                  seed=123,
#                                  iter=1000)


stan <- readRDS(file.path(models_dir,sprintf('stan_k_%s_fit%s.rds',K,fit_idx)))
ppd <- posterior_predict(stan$stan,draws=50)

data.frame(gene_mat,t(ppd)) %>%
  gather(sample,pred,-c(topic:logoffset)) %>%
  group_by(sample,pw) %>%
  summarise(pred=mean(pred),
            true=mean(count)) %>%
  ungroup() %>%
  gather(type,count,-pw,-sample) %>%
  ggplot(aes(pw,count,color=type,alpha=type)) +
  facet_wrap(~pw,scales='free_x') +
  geom_jitter() +
  scale_alpha_manual(values=c(.3,1)) +
  scale_color_manual(values=c('black','red'))


dotplot(ranef(stan1, condVar = TRUE), strip = FALSE)$topic
dotplot(ranef(stan1, condVar = TRUE), strip = FALSE)$pw
dotplot(ranef(stan1, condVar = TRUE), strip = FALSE)$`topic:pw`


int_coef <- ranef(stan1)$`topic:pw`
int_coef_names <- str_split(rownames(int_coef),':')
int_mat <- matrix(0,length(target_pws),K,dimnames=list(target_pws,1:K))
for (i in seq_len(nrow(int_coef))){
  k <- int_coef_names[[i]][1]
  pw <- int_coef_names[[i]][2]
  int_mat[pw,k] <- int_coef[i,1]
}

data.frame(int_mat,pw=rownames(int_mat)) %>%
  gather(topic,weight,-pw) %>%
  ggplot(aes(topic,pw,fill=weight)) +
  geom_raster() +
  scale_fill_viridis(name='weight')


















META <- dat$meta %>%
  mutate(age=ifelse(is.na(age) & diet == 'Not IBD',0,age)) %>%
  filter(!is.na(age)) %>%
  arrange(diet,desc(age))
theta <- theta[META$SampleID,]
rownames(theta) <- META$diet
colnames(theta) <- 1:NCOL(theta)



beta <- t(round(10000*exp(do.call('rbind',fit$beta$logbeta))))


taxon <- 7
taxon_names <- apply(TAX[as.character(dat$counts$ids[VOCAB,'long']),1:taxon],1,paste0,collapse='|')
taxon_table <- data.frame(beta,taxon_names) %>%
  group_by(taxon_names) %>%
  summarise_each(funs(sum)) 
taxon_names <- pretty_taxa_names(taxon_table$taxon_names)
taxon_table <- data.frame(dplyr::select(taxon_table,-taxon_names))
rownames(taxon_table) <- taxon_names
colnames(taxon_table) <- 1:NCOL(taxon_table)

filt_idx <- rowSums(taxon_table) != 0
taxon_table1 <- taxon_table[filt_idx,]




# prep data for kegg heatmap
####################################
####################################
biom_file <- read_biom(file.path(dir_name,kegg_filename))
kegg_metadata_risk_lev3 <- kegg_pw_collapse(biom_file,3)
kegg_metadata_risk_lev2 <- kegg_pw_collapse(biom_file,2)




my_kegg_order <- c('carbohydrate metabolism',
                    'amino acid metabolism',
                    'lipid metabolism',
                    'cell motility',
                    'glycan biosynthesis and metabolism',
                    'xenobiotics biodegradation and metabolism',
                    'metabolism of cofactors and vitamins',
                    'nucleotide metabolism',
                    'metabolism of terpenoids and polyketides',
                    'biosynthesis of other secondary metabolites',
                    'energy metabolism',
                    'metabolism',
                    'membrane transport',
                    'cellular processes and signaling',
                    'genetic information processing',
                    'unknown or misclassification')


kegg_risk <- as.matrix(biom_data(biom_file))

if (NCOL(kegg_risk) == K) {kegg_risk1 <- kegg_risk; kegg_risk2 <- NULL} else {kegg_risk1 <- kegg_risk[,1:K]; kegg_risk2 <- kegg_risk[,(K+1):NCOL(kegg_risk)]}


###
###
###
kegg_risk <- kegg_risk1
options(scipen=10000)
abund_thresh <- 2
abund_by_filt_x <- seq(1,1000,1)
abund_by_filt <- sapply(abund_by_filt_x,function(i) sum(rowSums(kegg_risk) >= i))
data.frame(x=abund_by_filt_x,y=abund_by_filt) %>%
  filter(x < 5000) %>%
  ggplot(aes(x,y)) + geom_point() + geom_vline(xintercept=abund_thresh,color='red')
Y <- rowSums(kegg_risk > 0)+1
X <- rowSums(kegg_risk)+1
qplot(X,Y,geom='point',alpha=.1) + xlab('Log10 Abundance') + ylab('Prevalence') +
  scale_x_log10() + theme(legend.position='none') + geom_vline(xintercept=abund_thresh,color='red')

kegg_risk <- kegg_risk[rowSums(kegg_risk) >= abund_thresh,]
kegg_gene_names <- kegg_pw_collapse(biom_file,0)
kegg_gene_names <- lapply(kegg_gene_names,function(x) paste0(x,collapse=' | '))

DF <- data.frame(coef=coef_k) %>%
  mutate(assoc=ifelse(coef>0,'CD+',ifelse(coef<0,'CD-','ZERO')),
         coef_scaled=coef_k/sd(coef_k),
         coef_ranked=rank(coef_k))
rownames(DF) <- colnames(kegg_risk)

KEGG <- otu_table(kegg_risk,taxa_are_rows=TRUE)
SAMP <- sample_data(DF)
PS <- phyloseq(KEGG,SAMP)

DS2 <- phyloseq_to_deseq2(PS, ~ coef_scaled)
DS2 <- DESeq2::DESeq(DS2, test="Wald", fitType="parametric")
res <- DESeq2::results(DS2, cooksCutoff=FALSE, pAdjustMethod='BH')
alpha <- .01
sigtab <- res[which(res$padj < alpha), ]

kegg_sig <- kegg_risk[rownames(sigtab),]
rownames(kegg_sig) <- kegg_gene_names[rownames(kegg_sig)]

ko_names <- rownames(kegg_sig) 
names(ko_names) <- rownames(sigtab)

# misc_filter <- !(rownames(kegg_sig) %in% c('None','hypothetical protein'))
# kegg_sig1 <- kegg_sig[misc_filter,]
# ko_names1 <- ko_names[misc_filter]
kegg_sig1 <- kegg_sig
ko_names1 <- ko_names

###
###
###

if (!is.null(kegg_risk2)){
  kegg_risk <- kegg_risk2
  options(scipen=10000)
  abund_thresh <- 2
  abund_by_filt_x <- seq(1,1000,1)
  abund_by_filt <- sapply(abund_by_filt_x,function(i) sum(rowSums(kegg_risk) >= i))
  data.frame(x=abund_by_filt_x,y=abund_by_filt) %>%
    filter(x < 200) %>%
    ggplot(aes(x,y)) + geom_point() + geom_vline(xintercept=abund_thresh,color='red')
  Y <- rowSums(kegg_risk > 0)+1
  X <- rowSums(kegg_risk)+1
  qplot(X,Y,geom='point',alpha=.1) + xlab('Log10 Abundance') + ylab('Prevalence') +
    scale_x_log10() + theme(legend.position='none') + geom_vline(xintercept=abund_thresh,color='red')
  
  kegg_risk <- kegg_risk[rowSums(kegg_risk) >= abund_thresh,]
  kegg_gene_names <- kegg_pw_collapse(biom_file,0)
  kegg_gene_names <- lapply(kegg_gene_names,function(x) paste0(x,collapse=' | '))
  
  DF <- data.frame(coef=coef_k) %>%
    mutate(assoc=ifelse(coef>0,'CD+',ifelse(coef<0,'CD-','ZERO')),
           coef_scaled=coef_k/sd(coef_k),
           coef_ranked=rank(coef_k))
  rownames(DF) <- colnames(kegg_risk)
  
  KEGG <- otu_table(kegg_risk,taxa_are_rows=TRUE)
  SAMP <- sample_data(DF)
  PS <- phyloseq(KEGG,SAMP)
  
  DS2 <- phyloseq_to_deseq2(PS, ~ coef_scaled)
  DS2 <- DESeq2::DESeq(DS2, test="Wald", fitType="parametric")
  res <- DESeq2::results(DS2, cooksCutoff=FALSE, pAdjustMethod='BH')
  alpha <- .01
  sigtab <- res[which(res$padj < alpha), ]
  
  kegg_sig <- kegg_risk[rownames(sigtab),]
  rownames(kegg_sig) <- kegg_gene_names[rownames(kegg_sig)]
  
  ko_names <- rownames(kegg_sig) 
  names(ko_names) <- rownames(sigtab)
  
#   misc_filter <- !(rownames(kegg_sig) %in% c('None','hypothetical protein'))
#   kegg_sig2 <- kegg_sig[misc_filter,]
#   ko_names2 <- ko_names[misc_filter]
  kegg_sig2 <- kegg_sig
  ko_names2 <- ko_names
}
####################################
####################################




# prep data for cog heatmap
####################################
####################################
biom_file <- read_biom(file.path(dir_name,cog_filename))
cog_metadata_risk <- kegg_pw_collapse(biom_file,2,cog=TRUE)
cog_risk <- as.matrix(biom_data(biom_file))

if (NCOL(cog_risk) == K) {cog_risk1 <- cog_risk; cog_risk2 <- NULL} else {cog_risk1 <- cog_risk[,1:K]; cog_risk2 <- cog_risk[,(K+1):NCOL(cog_risk)]}

###
###
###

cog_risk <- cog_risk1
options(scipen=10000)
abund_thresh <- 2
abund_by_filt_x <- seq(1,1000,1)
abund_by_filt <- sapply(abund_by_filt_x,function(i) sum(rowSums(cog_risk) >= i))
data.frame(x=abund_by_filt_x,y=abund_by_filt) %>%
  filter(x < 500) %>%
  ggplot(aes(x,y)) + geom_point() + geom_vline(xintercept=abund_thresh,color='red')
Y <- rowSums(cog_risk > 0)+1
X <- rowSums(cog_risk)+1
qplot(X,Y,geom='point',alpha=.1) + xlab('Log10 Abundance') + ylab('Prevalence') +
  scale_x_log10() + theme(legend.position='none') + geom_vline(xintercept=abund_thresh,color='red')

cog_risk <- cog_risk[rowSums(cog_risk) >= abund_thresh,]
cog_gene_names <- kegg_pw_collapse(biom_file,0,cog=TRUE)
cog_gene_names <- lapply(cog_gene_names,function(x) paste0(x,collapse=' | '))

DF <- data.frame(coef=coef_k) %>%
  mutate(assoc=ifelse(coef>0,'CD+',ifelse(coef<0,'CD-','ZERO')),
         coef_scaled=coef_k/sd(coef_k),
         coef_ranked=rank(coef_k))
rownames(DF) <- colnames(cog_risk)

COG <- otu_table(cog_risk,taxa_are_rows=TRUE)
SAMP <- sample_data(DF)
PS <- phyloseq(COG,SAMP)

DS2 <- phyloseq_to_deseq2(PS, ~ coef_ranked)
DS2 <- DESeq2::DESeq(DS2, test="Wald", fitType="parametric")
res <- DESeq2::results(DS2, cooksCutoff=FALSE, pAdjustMethod='BH')
alpha <- .01
sigtab <- res[which(res$padj < alpha), ]

cog_sig <- cog_risk[rownames(sigtab),]
rownames(cog_sig) <- cog_gene_names[rownames(cog_sig)]

cog_names <- rownames(cog_sig) 
names(cog_names) <- rownames(sigtab)

# misc_filter <- !(rownames(cog_sig) %in% c('None','hypothetical protein','Uncharacterized conserved protein','Uncharacterized conserved protein'))
# cog_sig1 <- cog_sig[misc_filter,]
# cog_names1 <- cog_names[misc_filter]
cog_sig1 <- cog_sig
cog_names1 <- cog_names

###
###
###

if (!is.null(cog_risk2)){
  cog_risk <- cog_risk2
  options(scipen=10000)
  abund_thresh <- 2
  abund_by_filt_x <- seq(1,1000,1)
  abund_by_filt <- sapply(abund_by_filt_x,function(i) sum(rowSums(cog_risk) >= i))
  data.frame(x=abund_by_filt_x,y=abund_by_filt) %>%
    filter(x < 500) %>%
    ggplot(aes(x,y)) + geom_point() + geom_vline(xintercept=abund_thresh,color='red')
  Y <- rowSums(cog_risk > 0)+1
  X <- rowSums(cog_risk)+1
  qplot(X,Y,geom='point',alpha=.1) + xlab('Log10 Abundance') + ylab('Prevalence') +
    scale_x_log10() + theme(legend.position='none') + geom_vline(xintercept=abund_thresh,color='red')
  
  cog_risk <- cog_risk[rowSums(cog_risk) >= abund_thresh,]
  cog_gene_names <- kegg_pw_collapse(biom_file,0,cog=TRUE)
  cog_gene_names <- lapply(cog_gene_names,function(x) paste0(x,collapse=' | '))
  
  DF <- data.frame(coef=coef_k) %>%
    mutate(assoc=ifelse(coef>0,'CD+',ifelse(coef<0,'CD-','ZERO')),
           coef_scaled=coef_k/sd(coef_k),
           coef_ranked=rank(coef_k))
  rownames(DF) <- colnames(cog_risk)
  
  COG <- otu_table(cog_risk,taxa_are_rows=TRUE)
  SAMP <- sample_data(DF)
  PS <- phyloseq(COG,SAMP)
  
  DS2 <- phyloseq_to_deseq2(PS, ~ coef_ranked)
  DS2 <- DESeq2::DESeq(DS2, test="Wald", fitType="parametric")
  res <- DESeq2::results(DS2, cooksCutoff=FALSE, pAdjustMethod='BH')
  alpha <- .01
  sigtab <- res[which(res$padj < alpha), ]
  
  cog_sig <- cog_risk[rownames(sigtab),]
  rownames(cog_sig) <- cog_gene_names[rownames(cog_sig)]
  
  cog_names <- rownames(cog_sig) 
  names(cog_names) <- rownames(sigtab)
  
#   misc_filter <- !(rownames(cog_sig) %in% c('None','hypothetical protein','Uncharacterized conserved protein','Uncharacterized conserved protein'))
#   cog_sig2 <- cog_sig[misc_filter,]
#   cog_names2 <- cog_names[misc_filter]   
  cog_sig2 <- cog_sig
  cog_names2 <- cog_names
}
####################################
####################################










coef_k_hm <- coef_k*5



pdf(file.path(heatmap_dir,'heatmap_theta_row.pdf'), height=20, width=20)
plot_big_heatmap3(mat,coef_k_hm,col0='black',col1='cyan',col2='orange',normtype='q',normaxis=1,rowlabs=FALSE,rowclust=FALSE,rowcols=TRUE,
                  ColSideLabs='Topic Weight',RowSideLabs='Disease',main='Samples over Topics',
                  cexCol=1.5,cexRow=.75)
dev.off()

pdf(file.path(heatmap_dir,'heatmap_theta_col.pdf'), height=20, width=20)
plot_big_heatmap3(mat,coef_k_hm,col0='black',col1='cyan',col2='orange',normtype='q',normaxis=2,rowlabs=FALSE,rowclust=FALSE,rowcols=TRUE,
                  ColSideLabs='Topic Weight',RowSideLabs='Disease',main='Samples over Topics',
                  cexCol=1.5,cexRow=.75)
dev.off()

pdf(file.path(heatmap_dir,'heatmap_theta_clust_row.pdf'), height=20, width=20)
plot_big_heatmap3(mat,coef_k_hm,col0='black',col1='cyan',col2='orange',normtype='q',normaxis=1,rowlabs=FALSE,rowclust=TRUE,rowcols=TRUE,
                  ColSideLabs='Topic Weight',RowSideLabs='Disease',main='Samples over Topics',
                  cexCol=1.5,cexRow=.75)
dev.off()

pdf(file.path(heatmap_dir,'heatmap_theta_clust_col.pdf'), height=20, width=20)
plot_big_heatmap3(mat,coef_k_hm,col0='black',col1='cyan',col2='orange',normtype='q',normaxis=2,rowlabs=FALSE,rowclust=TRUE,rowcols=TRUE,
                  ColSideLabs='Topic Weight',RowSideLabs='Disease',main='Samples over Topics',
                  cexCol=1.5,cexRow=.75)
dev.off()



mat_ends <- mat[,c(1:5,46:50)]
coef_k_hm_ends <- coef_k_hm[c(1:5,46:50)]
pdf(file.path(heatmap_dir,'heatmap_thetaends_row.pdf'), height=20, width=20)
plot_big_heatmap3(mat_ends,coef_k_hm_ends,col0='black',col1='cyan',col2='orange',normtype='q',normaxis=1,rowlabs=FALSE,rowclust=FALSE,rowcols=TRUE,
                  ColSideLabs='Topic Weight',RowSideLabs='Disease',main='Samples over Topics',
                  cexCol=1.5,cexRow=.75)
dev.off()

pdf(file.path(heatmap_dir,'heatmap_thetaends_col.pdf'), height=20, width=20)
plot_big_heatmap3(mat_ends,coef_k_hm_ends,col0='black',col1='cyan',col2='orange',normtype='q',normaxis=2,rowlabs=FALSE,rowclust=FALSE,rowcols=TRUE,
                  ColSideLabs='Topic Weight',RowSideLabs='Disease',main='Samples over Topics',
                  cexCol=1.5,cexRow=.75)
dev.off()

pdf(file.path(heatmap_dir,'heatmap_thetaends_clust_row.pdf'), height=20, width=20)
plot_big_heatmap3(mat_ends,coef_k_hm_ends,col0='black',col1='cyan',col2='orange',normtype='q',normaxis=1,rowlabs=FALSE,rowclust=TRUE,rowcols=TRUE,
                  ColSideLabs='Topic Weight',RowSideLabs='Disease',main='Samples over Topics',
                  cexCol=1.5,cexRow=.75)
dev.off()

pdf(file.path(heatmap_dir,'heatmap_thetaends_clust_col.pdf'), height=20, width=20)
plot_big_heatmap3(mat_ends,coef_k_hm_ends,col0='black',col1='cyan',col2='orange',normtype='q',normaxis=2,rowlabs=FALSE,rowclust=TRUE,rowcols=TRUE,
                  ColSideLabs='Topic Weight',RowSideLabs='Disease',main='Samples over Topics',
                  cexCol=1.5,cexRow=.75)
dev.off()





pdf(file.path(heatmap_dir,'heatmap_beta1.pdf'), height=20, width=20)
hm_beta1 <- plot_big_heatmap3(taxon_table1,coef_k_hm,col0='black',col1='cyan',col2='orange',normtype=NULL,normaxis=1,rowlabs=FALSE,rowclust=TRUE,rowcols='red',
                              dis='bray',scale='row',
                              ColSideLabs='Topic Weight',RowSideLabs='Disease',main='Topics over OTUs',
                              cexCol=1.5,cexRow=.75)
dev.off()

pdf(file.path(heatmap_dir,'heatmap_kegg1.pdf'), height=20, width=20)
plot_big_heatmap3(kegg_sig1,coef_k_hm,col0='black',col1='cyan',col2='orange',normtype=NULL,normaxis=1,rowlabs=FALSE,rowclust=TRUE,rowcols='red',
                  dis='bray',scale='row',
                  ColSideLabs='Topic Weight',RowSideLabs='Disease',main='Topics over KOs',
                  cexCol=1.5,cexRow=.75)
dev.off()

pdf(file.path(heatmap_dir,'heatmap_kegg_rownames1.pdf'), height=50, width=75)
hm_kegg1 <- plot_big_heatmap3(kegg_sig1,coef_k_hm,col0='black',col1='cyan',col2='orange',normtype=NULL,normaxis=1,rowlabs=TRUE,rowclust=TRUE,rowcols='red',
                              dis='bray',scale='row',
                              ColSideLabs='Topic Weight',RowSideLabs='Disease',main='Topics over KOs',
                              cexCol=5,cexRow=4,margins=c(20,50))
dev.off()

pdf(file.path(heatmap_dir,'heatmap_cog1.pdf'), height=20, width=20)
plot_big_heatmap3(cog_sig1,coef_k_hm,col0='black',col1='cyan',col2='orange',normtype=NULL,normaxis=1,rowlabs=FALSE,rowclust=TRUE,rowcols='red',
                  dis='bray',scale='row',
                  ColSideLabs='Topic Weight',RowSideLabs='Disease',main='Topics over COGs',
                  cexCol=1.5,cexRow=.75)
dev.off()

pdf(file.path(heatmap_dir,'heatmap_cog_rownames1.pdf'), height=60, width=75)
hm_cog1 <- plot_big_heatmap3(cog_sig1,coef_k_hm,col0='black',col1='cyan',col2='orange',normtype=NULL,normaxis=1,rowlabs=TRUE,rowclust=TRUE,rowcols='red',
                             dis='bray',scale='row',
                             ColSideLabs='Topic Weight',RowSideLabs='Disease',main='Topics over COGs',
                             cexCol=5,cexRow=4,margins=c(20,20))
dev.off()


side_names <- collapse_to_my_kegg_names(kegg_metadata_risk_lev2[names(ko_names1)], my_kegg_order)
side_names_ord <- order(side_names)
side_names <- side_names[side_names_ord]
side_idx <- as.integer(as.factor(side_names))

pdf(file.path(heatmap_dir,'heatmap_kegg_pathway1.pdf'), height=20, width=20)
plot_big_heatmap3(kegg_sig1[side_names_ord,],
                  coef_k,col0='black',col1='cyan',col2='orange',normtype=NULL,normaxis=1,rowlabs=FALSE,rowclust=FALSE,rowcols=side_names,
                  dis='bray',scale='row',
                  ColSideLabs='Topic Weight',RowSideLabs='Disease',main='Topics over KOs',
                  cexCol=1.5,cexRow=.75)
dev.off()


write.table(data.frame(PW=as.vector(side_names),
                       IDX=side_idx,
                       COLOR=rainbow(30)[side_idx]),
            quote=FALSE,
            file=file.path(heatmap_dir,'heatmap_kegg_pathway1_cckey.txt'))






if (!is.null(beta2)){
  
  pdf(file.path(heatmap_dir,'heatmap_beta2.pdf'), height=20, width=20)
  hm_beta2 <- plot_big_heatmap3(taxon_table2,coef_k_hm,col0='black',col1='cyan',col2='orange',normtype=NULL,normaxis=1,rowlabs=FALSE,rowclust=TRUE,rowcols='blue',
                                dis='bray',scale='row',
                                ColSideLabs='Topic Weight',RowSideLabs='Disease',main='Topics over OTUs',
                                cexCol=1.5,cexRow=.75)
  dev.off()
  
  pdf(file.path(heatmap_dir,'heatmap_kegg2.pdf'), height=20, width=20)
  plot_big_heatmap3(kegg_sig2,coef_k_hm,col0='black',col1='cyan',col2='orange',normtype=NULL,normaxis=1,rowlabs=FALSE,rowclust=TRUE,rowcols='blue',
                    dis='bray',scale='row',
                    ColSideLabs='Topic Weight',RowSideLabs='Disease',main='Topics over KOs',
                    cexCol=1.5,cexRow=.75)
  dev.off()
  
  pdf(file.path(heatmap_dir,'heatmap_kegg_rownames2.pdf'), height=50, width=75)
  hm_kegg2 <- plot_big_heatmap3(kegg_sig2,coef_k_hm,col0='black',col1='cyan',col2='orange',normtype=NULL,normaxis=1,rowlabs=TRUE,rowclust=TRUE,rowcols='blue',
                                dis='bray',scale='row',
                                ColSideLabs='Topic Weight',RowSideLabs='Disease',main='Topics over KOs',
                                cexCol=5,cexRow=4,margins=c(20,20))
  dev.off()
  
  pdf(file.path(heatmap_dir,'heatmap_cog2.pdf'), height=20, width=20)
  plot_big_heatmap3(cog_sig2,coef_k_hm,col0='black',col1='cyan',col2='orange',normtype=NULL,normaxis=1,rowlabs=FALSE,rowclust=TRUE,rowcols='blue',
                    dis='bray',scale='row',
                    ColSideLabs='Topic Weight',RowSideLabs='Disease',main='Topics over COGs',
                    cexCol=1.5,cexRow=.75)
  dev.off()
  
  pdf(file.path(heatmap_dir,'heatmap_cog_rownames2.pdf'), height=60, width=75)
  hm_cog2 <- plot_big_heatmap3(cog_sig2,coef_k_hm,col0='black',col1='cyan',col2='orange',normtype=NULL,normaxis=1,rowlabs=TRUE,rowclust=TRUE,rowcols='blue',
                               dis='bray',scale='row',
                               ColSideLabs='Topic Weight',RowSideLabs='Disease',main='Topics over COGs',
                               cexCol=5,cexRow=4,margins=c(20,20))
  dev.off()
  
  side_names <- collapse_to_my_kegg_names(kegg_metadata_risk_lev2[names(ko_names2)], my_kegg_order)
  side_names_ord <- order(side_names)
  side_names <- side_names[side_names_ord]
  side_idx <- as.integer(as.factor(side_names))
  
  
  pdf(file.path(heatmap_dir,'heatmap_kegg_pathway2.pdf'), height=20, width=20)
  plot_big_heatmap3(kegg_sig2[side_names_ord,],
                    coef_k,col0='black',col1='cyan',col2='orange',normtype=NULL,normaxis=1,rowlabs=FALSE,rowclust=FALSE,rowcols=side_names,
                    dis='bray',scale='row',
                    ColSideLabs='Topic Weight',RowSideLabs='Disease',main='Topics over KOs',
                    cexCol=1.5,cexRow=.75)
  dev.off()
  
  write.table(data.frame(PW=as.vector(side_names),
                         IDX=side_idx,
                         COLOR=rainbow(30)[side_idx]),
              quote=FALSE,
              file=file.path(heatmap_dir,'heatmap_kegg_pathway2_cckey.txt'))
  
}








cut_folder <- 'beta1'
hm_beta <- hm_beta1
taxon_table <- taxon_table1
cut_path <- file.path(heatmap_dir,cut_folder)
dir.create(cut_path,showWarnings=FALSE,recursive=TRUE)

tree_groups <- cutree(hm_beta,6)
table(tree_groups) # number of genes in a group
tree_groups <- tree_groups[rev(hm_beta$order)] 

sink(file.path(cut_path,'group_order.txt'))
cat(unique(tree_groups))
sink()
for (g in unique(tree_groups)){
  group_file_name <- paste0('heatmap_beta_group_',g)
  idx_groups <- rev(rev(hm_beta$order)[tree_groups == g]) # capture 4th group
  
  group_mat <- taxon_table[idx_groups,]
  h <- ifelse(floor(nrow(group_mat)/3) < 8, 8, floor(nrow(group_mat)/3))
  
  pdf(file.path(cut_path,paste0(group_file_name,'.pdf')), height=h, width=h*1.25)
  plot_big_heatmap3(group_mat,coef_k_hm,col0='black',col1='cyan',col2='orange',normtype=NULL,normaxis=1,rowlabs=TRUE,rowclust=FALSE,rowcols='red',
                    dis='bray',scale='row',
                    ColSideLabs='Topic Weight',RowSideLabs='Disease',main='Topics over OTUs',
                    cexCol=log(h,10),cexRow=log(h,7))
  dev.off()
}



cut_folder <- 'kegg1'
hm_kegg <- hm_kegg1
kegg_sig <- kegg_sig1
ko_names <- ko_names1
cut_path <- file.path(heatmap_dir,cut_folder)
dir.create(cut_path,showWarnings=FALSE,recursive=TRUE)

tree_groups <- cutree(hm_kegg,h=2)
table(tree_groups) # number of genes in a group
tree_groups <- tree_groups[rev(hm_kegg$order)] # reorder the group indexes to match the clustering indexes
# group_idx <- 4
# idx_groups <- rev(hm_kegg$order)[tree_groups == group_idx] # capture 4th group
# rownames(kegg_sig)[idx_groups] # checker

sink(file.path(cut_path,'group_order.txt'))
cat(unique(tree_groups))
sink()

sink(file.path(heatmap_dir,'heatmap_kegg_rownames1.txt'))
for (i in seq_along(hm_kegg$labels[hm_kegg$order])) cat(hm_kegg$labels[hm_kegg$order][i],'\n')
sink()

for (g in unique(tree_groups)){
  group_file_name <- paste0('heatmap_kegg_group_',g)
  idx_groups <- rev(rev(hm_kegg$order)[tree_groups == g]) # capture 4th group
  
  
  side_names <- collapse_to_my_kegg_names(kegg_metadata_risk_lev2[names(ko_names)[idx_groups]], my_kegg_order)
  side_names_ord <- order(side_names)
  
  pws <- kegg_metadata_risk_lev3[names(ko_names)[idx_groups]]
  pws2 <- kegg_metadata_risk_lev2[names(ko_names)[idx_groups]]
  
  
  sink(file.path(cut_path,paste0(group_file_name,'.txt')))
  cat('Group',g,'\n\n\n')
  cat('Gene Name [EC Identifier] (KO term): Pathways\n\n')
  for (i in seq_along(pws)){
    cat(ko_names[names(pws[i])],' (',names(pws[i]),'):  ',paste0(pws[[i]],collapse=', '),'\n',sep='')
    cat(ko_names[names(pws[i])],' (',names(pws[i]),'):  ',paste0(pws2[[i]],collapse=', '),'\n',sep='')
    cat(ko_names[names(pws[i])],' (',names(pws[i]),'):  ',paste0(side_names[[i]],collapse=', '),'\n\n',sep='')
  }
  sink()
  
  group_mat <- kegg_sig[idx_groups,]
  
  
  side_names <- side_names[side_names_ord]
  group_mat <- group_mat[side_names_ord,]
  
  
  h <- ifelse(floor(nrow(group_mat)/3) < 8, 8, floor(nrow(group_mat)/3))
  
  pdf(file.path(cut_path,paste0(group_file_name,'.pdf')), height=h, width=h*1.75)
  plot_big_heatmap3(group_mat,coef_k_hm,col0='black',col1='cyan',col2='orange',normtype=NULL,normaxis=1,rowlabs=TRUE,rowclust=FALSE,rowcols=side_names,
                    dis='bray',scale='row',
                    ColSideLabs='Topic Weight',RowSideLabs='Pathway',main='Topics over KOs',
                    cexCol=log(h,10),cexRow=log(h,7))
  dev.off()
}





cut_folder <- 'cog1'
hm_cog <- hm_cog1
cog_sig <- cog_sig1
cog_names <- cog_names1
cut_path <- file.path(heatmap_dir,cut_folder)
dir.create(cut_path,showWarnings=FALSE,recursive=TRUE)

tree_groups <- cutree(hm_cog,h=2)
table(tree_groups) # number of genes in a group
tree_groups <- tree_groups[rev(hm_cog$order)] # reorder the group indexes to match the clustering indexes
# group_idx <- 4
# idx_groups <- rev(hm_cog$order)[tree_groups == group_idx] # capture 4th group
# rownames(kegg_sig)[idx_groups] # checker

sink(file.path(cut_path,'group_order.txt'))
cat(unique(tree_groups))
sink()
for (g in unique(tree_groups)){
  group_file_name <- paste0('heatmap_cog_group_',g)
  idx_groups <- rev(rev(hm_cog$order)[tree_groups == g]) # capture 4th group
  
  pws <- cog_metadata_risk[names(cog_names)[idx_groups]]
  sink(file.path(cut_path,paste0(group_file_name,'.txt')))
  cat('Group',g,'\n\n\n')
  cat('Gene Name (COG term): Pathways\n')
  for (i in seq_along(pws)){
    cat(cog_names[names(pws[i])],' (',names(pws[i]),'):  ',paste0(pws[[i]],collapse=', '),'\n',sep='')
  }
  sink()
  
  group_mat <- cog_sig[idx_groups,]
  h <- ifelse(floor(nrow(group_mat)/3) < 8, 8, floor(nrow(group_mat)/3))
  
  pdf(file.path(cut_path,paste0(group_file_name,'.pdf')), height=h, width=h*1.85)
  plot_big_heatmap3(group_mat,coef_k_hm,col0='black',col1='cyan',col2='orange',normtype=NULL,normaxis=1,rowlabs=TRUE,rowclust=FALSE,rowcols='red',
                    dis='bray',scale='row',
                    ColSideLabs='Topic Weight',RowSideLabs='Disease',main='Topics over COGs',
                    cexCol=log(h,8),cexRow=log(h,7))
  dev.off()
}






if (!is.null(beta2)){
  
  cut_folder <- 'beta2'
  hm_beta <- hm_beta2
  taxon_table <- taxon_table2
  cut_path <- file.path(heatmap_dir,cut_folder)
  dir.create(cut_path,showWarnings=FALSE,recursive=TRUE)
  
  tree_groups <- cutree(hm_beta,6)
  table(tree_groups) # number of genes in a group
  tree_groups <- tree_groups[rev(hm_beta$order)] 
  
  sink(file.path(cut_path,'group_order.txt'))
  cat(unique(tree_groups))
  sink()
  for (g in unique(tree_groups)){
    group_file_name <- paste0('heatmap_beta_group_',g)
    idx_groups <- rev(rev(hm_beta$order)[tree_groups == g]) # capture 4th group
    
    group_mat <- taxon_table[idx_groups,]
    h <- ifelse(floor(nrow(group_mat)/3) < 8, 8, floor(nrow(group_mat)/3))
    
    pdf(file.path(cut_path,paste0(group_file_name,'.pdf')), height=h, width=15)
    plot_big_heatmap3(group_mat,coef_k_hm,col0='black',col1='cyan',col2='orange',normtype=NULL,normaxis=1,rowlabs=TRUE,rowclust=FALSE,rowcols='blue',
                      dis='bray',scale='row',
                      ColSideLabs='Topic Weight',RowSideLabs='Disease',main='Topics over OTUs',
                      cexCol=log(h,10),cexRow=log(h,7))
    dev.off()
  }
  
  
  cut_folder <- 'kegg2'
  hm_kegg <- hm_kegg2
  kegg_sig <- kegg_sig2
  ko_names <- ko_names2
  cut_path <- file.path(heatmap_dir,cut_folder)
  dir.create(cut_path,showWarnings=FALSE,recursive=TRUE)
  
  tree_groups <- cutree(hm_kegg,h=2)
  table(tree_groups) # number of genes in a group
  tree_groups <- tree_groups[rev(hm_kegg$order)] # reorder the group indexes to match the clustering indexes
  # group_idx <- 4
  # idx_groups <- rev(hm_kegg$order)[tree_groups == group_idx] # capture 4th group
  # rownames(kegg_sig)[idx_groups] # checker
  
  sink(file.path(cut_path,'group_order.txt'))
  cat(unique(tree_groups))
  sink()
  
  sink(file.path(heatmap_dir,'heatmap_kegg_rownames2.txt'))
  for (i in seq_along(hm_kegg$labels[hm_kegg$order])) cat(hm_kegg$labels[hm_kegg$order][i],'\n')
  sink()
  
  for (g in unique(tree_groups)){
    group_file_name <- paste0('heatmap_kegg_group_',g)
    idx_groups <- rev(rev(hm_kegg$order)[tree_groups == g]) # capture 4th group
    
    
    side_names <- collapse_to_my_kegg_names(kegg_metadata_risk_lev2[names(ko_names)[idx_groups]], my_kegg_order)
    side_names_ord <- order(side_names)
    
    pws <- kegg_metadata_risk_lev3[names(ko_names)[idx_groups]]
    pws2 <- kegg_metadata_risk_lev2[names(ko_names)[idx_groups]]
    
    
    sink(file.path(cut_path,paste0(group_file_name,'.txt')))
    cat('Group',g,'\n\n\n')
    cat('Gene Name [EC Identifier] (KO term): Pathways\n\n')
    for (i in seq_along(pws)){
      cat(ko_names[names(pws[i])],' (',names(pws[i]),'):  ',paste0(pws[[i]],collapse=', '),'\n',sep='')
      cat(ko_names[names(pws[i])],' (',names(pws[i]),'):  ',paste0(pws2[[i]],collapse=', '),'\n',sep='')
      cat(ko_names[names(pws[i])],' (',names(pws[i]),'):  ',paste0(side_names[[i]],collapse=', '),'\n\n',sep='')
    }
    sink()
    
    group_mat <- kegg_sig[idx_groups,]
    
    
    side_names <- side_names[side_names_ord]
    group_mat <- group_mat[side_names_ord,]
    
    
    h <- ifelse(floor(nrow(group_mat)/3) < 8, 8, floor(nrow(group_mat)/3))
    
    pdf(file.path(cut_path,paste0(group_file_name,'.pdf')), height=h, width=h*1.75)
    plot_big_heatmap3(group_mat,coef_k_hm,col0='black',col1='cyan',col2='orange',normtype=NULL,normaxis=1,rowlabs=TRUE,rowclust=FALSE,rowcols=side_names,
                      dis='bray',scale='row',
                      ColSideLabs='Topic Weight',RowSideLabs='Pathway',main='Topics over KOs',
                      cexCol=log(h,10),cexRow=log(h,7))
    dev.off()
  }
  
  
  cut_folder <- 'cog2'
  hm_cog <- hm_cog2
  cog_sig <- cog_sig2
  cog_names <- cog_names2
  cut_path <- file.path(heatmap_dir,cut_folder)
  dir.create(cut_path,showWarnings=FALSE,recursive=TRUE)
  
  tree_groups <- cutree(hm_cog,h=2)
  table(tree_groups) # number of genes in a group
  tree_groups <- tree_groups[rev(hm_cog$order)] # reorder the group indexes to match the clustering indexes
  # group_idx <- 4
  # idx_groups <- rev(hm_cog$order)[tree_groups == group_idx] # capture 4th group
  # rownames(kegg_sig)[idx_groups] # checker
  
  sink(file.path(cut_path,'group_order.txt'))
  cat(unique(tree_groups))
  sink()
  for (g in unique(tree_groups)){
    group_file_name <- paste0('heatmap_cog_group_',g)
    idx_groups <- rev(rev(hm_cog$order)[tree_groups == g]) # capture 4th group
    
    pws <- cog_metadata_risk[names(cog_names)[idx_groups]]
    sink(file.path(cut_path,paste0(group_file_name,'.txt')))
    cat('Group',g,'\n\n\n')
    cat('Gene Name (COG term): Pathways\n')
    for (i in seq_along(pws)){
      cat(cog_names[names(pws[i])],' (',names(pws[i]),'):  ',paste0(pws[[i]],collapse=', '),'\n',sep='')
    }
    sink()
    
    group_mat <- cog_sig[idx_groups,]
    h <- ifelse(floor(nrow(group_mat)/3) < 8, 8, floor(nrow(group_mat)/3))
    
    pdf(file.path(cut_path,paste0(group_file_name,'.pdf')), height=h, width=15)
    plot_big_heatmap3(group_mat,coef_k_hm,col0='black',col1='cyan',col2='orange',normtype=NULL,normaxis=1,rowlabs=TRUE,rowclust=FALSE,rowcols='blue',
                      dis='bray',scale='row',
                      ColSideLabs='Topic Weight',RowSideLabs='Disease',main='Topics over COGs',
                      cexCol=log(h,8),cexRow=log(h,7))
    dev.off()
  }
  
}















sig_corr_topics <- which(sapply(eff_plot$cis,function(x) sum(sign(x))) != 0)
pos_corr_topics <- which(sapply(eff_plot$cis,function(x) sum(sign(x))) == 2)
neg_corr_topics <- which(sapply(eff_plot$cis,function(x) sum(sign(x))) == -2)


















     

pdf(file.path(out_dir,K,mod_dir,'vio.pdf'), height=4, width=12)
data.frame(fit$theta,id=dat$meta$SampleID,dx=dat$meta$diet,age=dat$meta$age,row.names=NULL) %>%
  gather(topic,p,-dx,-age,-id) %>%
  mutate(topic=str_replace(topic,'X','T'),
         corr=as.factor(
           ifelse(topic %in% paste0('T',pos_corr_topics),1,
                      ifelse(topic %in% paste0('T',neg_corr_topics),0,NA)))) %>%
  filter(!is.na(corr)) %>%
  group_by(id,corr) %>%
  summarise(p=sum(p),
            age=unique(age),
            dx=as.factor(unique(dx))) %>%
  #ggplot(aes(p,fill=dx)) + facet_grid(corr~.) + geom_density(colour='black',alpha=.5) # sanity check
  mutate(status=ifelse(age>0,1,0)) %>%
  group_by(status) %>%
  mutate(bin=ntile(age,3)) %>%
  ungroup() %>%
  mutate(bin=as.factor(ifelse(dx==1,bin,0))) %>%
  ggplot(aes(x=dx,y=p)) + 
  geom_point(aes(colour=bin),size=1.75,alpha=.8,position=position_jitter(width=1,height=0)) + 
  geom_violin(fill=NA,size=1.75) + 
  scale_x_discrete(labels=c('CD+ Subjects','CD- Subjects')) +
  scale_y_continuous(trans='log10') +
  facet_grid(.~corr,labeller=labeller(corr=c('0'='CD- Topics','1'='CD+ Topics'))) +
  scale_colour_brewer(type='div',palette=1) +
  labs(colour='Burden',x='',y='Probability') + 
  theme_bw() +
  theme(axis.text.x=element_text(size=20),
        axis.text.y=element_text(size=14),
        axis.title.y=element_text(size=20,face='bold'),
        legend.title=element_text(size=20,face='bold'),
        legend.text=element_text(size=20),
        strip.text=element_text(size=20,face='bold'),
        aspect.ratio=2/4,
        legend.position='bottom')
dev.off()
  
  


# library(splines)
# pdf(file.path(out_dir,K,mod_dir,'vsage.pdf'), height=10, width=12)
# data.frame(fit$theta,id=dat$meta$SampleID,dx=dat$meta$diet,age=dat$meta$age,row.names=NULL) %>%
#   gather(topic,p,-dx,-age,-id) %>%
#   mutate(topic=str_replace(topic,'X','T'),
#          corr=factor(
#            ifelse(topic %in% paste0('T',pos_corr_topics),1,
#                   ifelse(topic %in% paste0('T',neg_corr_topics),0,NA)),levels=1:0)) %>%
#   filter(!is.na(corr),
#          dx == 1, 
#          !is.na(age)) %>%
#   ggplot(aes(x=age,y=p,colour=age)) + 
#   geom_point(size=1.75,alpha=.8) + 
#   stat_smooth(method='lm',formula=y~x,se=FALSE,size=2,colour='black') + 
#   scale_y_continuous(trans='log10',labels=scales::trans_format('log10', scales::math_format(10^.x))) +
#   facet_wrap(corr~topic,labeller=labeller(corr=c('0'='CD- Topics','1'='CD+ Topics'))) +
#   labs(colour='Burden',x='',y='Probability') + 
#   scale_colour_distiller(type='seq',palette=7) +
#   theme_bw() +
#   theme(axis.text.x=element_text(size=20),
#         axis.text.y=element_text(size=14),
#         axis.title.y=element_text(size=20,face='bold'),
#         legend.title=element_text(size=20,face='bold'),
#         legend.text=element_text(size=16),
#         strip.text=element_text(size=19,face='bold'),
#         aspect.ratio=1,
#         legend.position='none')
# dev.off()


library(splines)
pdf(file.path(out_dir,K,mod_dir,'vsage.pdf'), height=6, width=12)
data.frame(fit$theta,id=dat$meta$SampleID,dx=dat$meta$diet,age=dat$meta$age,row.names=NULL) %>%
  gather(topic,p,-dx,-age,-id) %>%
  mutate(topic=str_replace(topic,'X','T'),
         corr=factor(
           ifelse(topic %in% paste0('T',pos_corr_topics),1,
                  ifelse(topic %in% paste0('T',neg_corr_topics),0,NA)),levels=1:0)) %>%
  filter(!is.na(corr),
         dx == 1, 
         !is.na(age)) %>%
  ggplot(aes(x=age,y=p,colour=age)) + 
  geom_point(size=1.75,alpha=.8) + 
  stat_smooth(method='lm',formula=y~x,se=FALSE,size=2,colour='black') + 
  scale_y_continuous(trans='log10',labels=scales::trans_format('log10', scales::math_format(10^.x))) +
  facet_wrap(corr~topic,labeller=labeller(corr=c('0'='CD- Topics','1'='CD+ Topics')),nrow=2) +
  labs(colour='Burden',x='Disease Burden (age)',y='Probability') + 
  scale_colour_distiller(type='seq',palette=7) +
  theme_bw() +
  theme(axis.text.x=element_text(size=18),
        axis.text.y=element_text(size=14),
        axis.title=element_text(size=20,face='bold'),
        legend.title=element_text(size=20,face='bold'),
        legend.text=element_text(size=16),
        strip.text=element_text(size=19,face='bold'),
        aspect.ratio=1,
        legend.position='none')
dev.off()











meta$age_BIN <- meta %>% mutate(status=ifelse(age>0,1,0)) %>%
  group_by(status) %>%
  mutate(bin=ntile(age,3)) %>%
  ungroup() %>%
  mutate(bin=as.factor(ifelse(diet == 1,bin,0))) %>%
  dplyr::select(bin) %>%
  unlist()






qmat <- apply(mat,2,qnormalize)
pca <- prcomp(qmat,center=TRUE,.scale=TRUE)

df <- data.frame(x=pca$x[,1],y=pca$x[,2],diet=as.factor(meta$diet))
tpca1 <- df %>%
  ggplot(aes(x=x,y=y,colour=diet)) + 
  geom_point(alpha=1,size=1.5) +
  facet_wrap(~diet,nrow=1) +
  geom_point(data=dplyr::select(df,-diet),colour='gray',alpha=.4,size=.8) +
  theme_bw() + labs(x='PC1',y='PC2',colour='diet') + theme(legend.position='none') +
  ggtitle('PCA: diet')

df <- data.frame(x=pca$x[,1],y=pca$x[,2],age_BIN=as.factor(meta$age_BIN),diet=meta$diet) %>% filter(diet == 1)
tpca2 <- df %>%
  ggplot(aes(x=x,y=y,colour=age_BIN)) + 
  geom_point(alpha=1,size=1.5) +
  facet_wrap(~age_BIN,nrow=2) +
  geom_point(data=dplyr::select(df,-age_BIN),colour='gray',alpha=.4,size=.8) +
  theme_bw() + labs(x='PC1',y='PC2',colour='age Bin') + theme(legend.position='none') +
  ggtitle('PCA: age')


gridExtra::grid.arrange(nrow=2,tpca1,tpca2)
pdf(file.path(out_dir,K,mod_dir,'pca.pdf'), height=10, width=10); gridExtra::grid.arrange(nrow=2,tpca1,tpca2); dev.off()
# plotly::plot_ly(data.frame(x=pca$x[,1],y=pca$x[,2],z=pca$x[,3]),
#         x=x,y=y,z=z,
#         type='scatter3d',color=rownames(qmat),
#         mode='markers')



tsne <- Rtsne::Rtsne(qmat,n.comp=2)

df <- data.frame(x=tsne$Y[,1],y=tsne$Y[,2],diet=as.factor(meta$diet))
tsne1 <- df %>%
  ggplot(aes(x=x,y=y,colour=diet)) + 
  geom_point(alpha=1,size=1.5) +
  facet_wrap(~diet,nrow=1) +
  geom_point(data=dplyr::select(df,-diet),colour='gray',alpha=.4,size=.8) +
  theme_bw() + labs(x='PC1',y='PC2',colour='diet') + theme(legend.position='none') +
  ggtitle('tsne: diet')

df <- data.frame(x=tsne$Y[,1],y=tsne$Y[,2],age_BIN=as.factor(meta$age_BIN),diet=meta$diet) %>% filter(diet == 1)
tsne2 <- df %>%
  ggplot(aes(x=x,y=y,colour=age_BIN)) + 
  geom_point(alpha=1,size=1.5) +
  facet_wrap(~age_BIN,nrow=2) +
  geom_point(data=dplyr::select(df,-age_BIN),colour='gray',alpha=.4,size=.8) +
  theme_bw() + labs(x='PC1',y='PC2',colour='age Bin') + theme(legend.position='none') +
  ggtitle('tsne: age')


gridExtra::grid.arrange(nrow=2,tsne1,tsne2)
pdf(file.path(out_dir,K,mod_dir,'tsne.pdf'), height=10, width=10); gridExtra::grid.arrange(nrow=2,tsne1,tsne2); dev.off()
# plotly::plot_ly(data.frame(x=ica$S[,1],y=ica$S[,2],z=ica$S[,3]),
#         x=x,y=y,z=z,
#         type='scatter3d',color=rownames(qmat),
#         mode='markers')









pca <- prcomp(qmat[,c(pos_corr_topics,neg_corr_topics)],center=TRUE,.scale=TRUE)

df <- data.frame(x=pca$x[,1],y=pca$x[,2],diet=as.factor(meta$diet))
tpca1 <- df %>%
  ggplot(aes(x=x,y=y,colour=diet)) + 
  geom_point(alpha=1,size=1.5) +
  facet_wrap(~diet,nrow=1) +
  geom_point(data=dplyr::select(df,-diet),colour='gray',alpha=.4,size=.8) +
  theme_bw() + labs(x='PC1',y='PC2',colour='diet') + theme(legend.position='none') +
  ggtitle('PCA: diet')

df <- data.frame(x=pca$x[,1],y=pca$x[,2],age_BIN=as.factor(meta$age_BIN),diet=meta$diet) %>% filter(diet == 1)
tpca2 <- df %>%
  ggplot(aes(x=x,y=y,colour=age_BIN)) + 
  geom_point(alpha=1,size=1.5) +
  facet_wrap(~age_BIN,nrow=2) +
  geom_point(data=dplyr::select(df,-age_BIN),colour='gray',alpha=.4,size=.8) +
  theme_bw() + labs(x='PC1',y='PC2',colour='age Bin') + theme(legend.position='none') +
  ggtitle('PCA: age')


gridExtra::grid.arrange(nrow=2,tpca1,tpca2)
pdf(file.path(out_dir,K,mod_dir,'toptpca.pdf'), height=10, width=10); gridExtra::grid.arrange(nrow=2,tpca1,tpca2); dev.off()








tsne <- Rtsne::Rtsne(qmat[,c(pos_corr_topics,neg_corr_topics)],n.comp=2)

df <- data.frame(x=tsne$Y[,1],y=tsne$Y[,2],diet=as.factor(meta$diet))
tsne1 <- df %>%
  ggplot(aes(x=x,y=y,colour=diet)) + 
  geom_point(alpha=1,size=1.5) +
  facet_wrap(~diet,nrow=1) +
  geom_point(data=dplyr::select(df,-diet),colour='gray',alpha=.4,size=.8) +
  theme_bw() + labs(x='PC1',y='PC2',colour='diet') + theme(legend.position='none') +
  ggtitle('tsne: diet')

df <- data.frame(x=tsne$Y[,1],y=tsne$Y[,2],age_BIN=as.factor(meta$age_BIN),diet=meta$diet) %>% filter(diet == 1)
tsne2 <- df %>%
  ggplot(aes(x=x,y=y,colour=age_BIN)) + 
  geom_point(alpha=1,size=1.5) +
  facet_wrap(~age_BIN,nrow=2) +
  geom_point(data=dplyr::select(df,-age_BIN),colour='gray',alpha=.4,size=.8) +
  theme_bw() + labs(x='PC1',y='PC2',colour='age Bin') + theme(legend.position='none') +
  ggtitle('tsne: age')


gridExtra::grid.arrange(nrow=2,tsne1,tsne2)
pdf(file.path(out_dir,K,mod_dir,'toptsne.pdf'), height=10, width=10); gridExtra::grid.arrange(nrow=2,tsne1,tsne2); dev.off()
# plotly::plot_ly(data.frame(x=ica$S[,1],y=ica$S[,2],z=ica$S[,3]),
#         x=x,y=y,z=z,
#         type='scatter3d',color=rownames(qmat),
#         mode='markers')














perm25_fit2 <- readRDS("~/Dropbox/stm_microbiome/data_active/Gevers_ti/stm_s97_%s_1000_supervised/perm/perm_K_25_fit2.rds")


sig_K <- c(4,17,16)
sig_ref <- data.frame(do.call('rbind',perm25_fit2$perm$ref[sig_K]),
                      Topic=sig_K,
                      Sim=0) %>%
  rename(Lower=`X2.5.`,Mean=`X50.`,Upper=`X97.5.`) 

pdf(file.path(out_dir,25,'fit_2','perm25_2.pdf'), height=4, width=12)
data.frame(do.call('rbind',lapply(perm25_fit2$perm$permute[sig_K],function(x) do.call('rbind',x))),
           Topic=rep(sig_K,each=25),
           Sim=rep(1:25,3)) %>%
  rename(Lower=`X2.5.`,Mean=`X50.`,Upper=`X97.5.`) %>%
  ggplot() +
  geom_hline(yintercept=0) + 
  geom_pointrange(aes(x=Sim,y=Mean,ymin=Lower,ymax=Upper),colour='black',alpha=.75) + 
  facet_grid(.~Topic,labeller=labeller(Topic=c('4'='T4','16'='T16','17'='T17'))) +
  geom_pointrange(data=sig_ref,aes(x=0,y=Mean,ymin=Lower,ymax=Upper),colour=scales::muted('blue'),size=1.5) +
  labs(x='Simulation',y='Estimate') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_text(size=14),
        axis.title=element_text(size=20,face='bold'),
        legend.title=element_text(size=20,face='bold'),
        legend.text=element_text(size=20),
        strip.text=element_text(size=20,face='bold'),
        aspect.ratio=1)
dev.off()
















# THIS IS WHAT I DID FOR THE POSTER! (SHOULDNT BE COLSUMS; SHOULD BE COLMEANS)
# #target_pws <- my_kegg_order[-c(16)] 
# target_pws <- my_kegg_order[-c(5,8,7,9,16)]
# target_kos <- lapply(target_pws,function(pw) names(which(sapply(kegg_metadata_risk_lev2, function(x) any(str_to_lower(x) %in% pw)))))
# 
# ra_kegg_risk1 <- t(t(kegg_risk1)/colSums(kegg_risk1)) ### CONVERTED TO RA
# # ra_kegg_risk2 <- t(t(kegg_risk2)/colSums(kegg_risk2))
# 
# 
# #ra_kegg_risk1 <- round(1000*ra_kegg_risk1)
# 
# 
# target_kegg1 <- t(sapply(target_kos, function(i) colSums(ra_kegg_risk1[i,])))
# # target_kegg2 <- t(sapply(target_kos, function(i) colSums(ra_kegg_risk2[i,])))
# 
# # target_kegg1 <- t(sapply(target_kos, function(i) colSums(kegg_risk1[i,])))
# # target_kegg2 <- t(sapply(target_kos, function(i) colSums(kegg_risk2[i,])))
# 



# 
# target_pws <- my_kegg_order[-c(16)] 
# #target_pws <- my_kegg_order[-c(5,8,7,9,16)]
# target_kos <- lapply(target_pws,function(pw) names(which(sapply(kegg_metadata_risk_lev2, function(x) any(str_to_lower(x) %in% pw)))))
# target_kegg1_counts <- t(sapply(target_kos, function(i) colSums(kegg_risk1[i,])))
# target_kegg1 <- t(t(target_kegg1_counts)/colSums(target_kegg1_counts)) ### CONVERTED TO RA



kegg_risk1_colsums <- colSums(kegg_risk1)
target_pws <- my_kegg_order[-c(16)] 
#target_pws <- my_kegg_order[-c(5,8,7,9,16)]
target_kos <- lapply(target_pws,function(pw) names(which(sapply(kegg_metadata_risk_lev2, function(x) any(str_to_lower(x) %in% pw)))))
target_kegg1_counts <- t(sapply(target_kos, function(i) colSums(kegg_risk1[i,])))
target_kegg1 <- t(t(target_kegg1_counts)/kegg_risk1_colsums)


# DF <- data.frame(coef=coef_k_unnorm) %>%
#   mutate(assoc=ifelse(coef>0,'1',ifelse(coef<0,'0','UNKNOWN')),
#          coef_scaled=coef_k/sd(coef_k),
#          coef_ranked=rank(coef_k),
#          topic=row_number()) %>%
#   filter(topic %in% sig_corr_topics)
# rownames(DF) <- paste0('T',DF$topic)
# 
# rownames(target_kegg1_counts) <- target_pws
# PWS <- otu_table(target_kegg1_counts[,sig_corr_topics],taxa_are_rows=TRUE)
# SAMP <- sample_data(DF)
# PS <- phyloseq(PWS,SAMP)
# 
# DS2 <- phyloseq_to_deseq2(PS, ~ assoc)
# DS2 <- DESeq2::DESeq(DS2, test="Wald", fitType="parametric")
# res <- DESeq2::results(DS2, cooksCutoff=FALSE, pAdjustMethod='BH')
# alpha <- .05
# sigtab <- res[which(res$padj < alpha), ]
# sigtab








kegg_risk1_colsums <- colSums(kegg_risk1)
target_pws <- my_kegg_order[-c(16)] 
#target_pws <- my_kegg_order[-c(5,8,7,9,16)]
target_kos <- lapply(target_pws,function(pw) names(which(sapply(kegg_metadata_risk_lev2, function(x) any(str_to_lower(x) %in% pw)))))
names(target_kos) <- target_pws


j1 <- 1
j2 <- 0
gene_mat <- matrix(0,sum(sapply(target_kos,length)) * K,6,dimnames=list(NULL,c('topic','pw','ko','count','offset','coef')))
for (i in seq_along(target_kos)){
  pw <- target_kos[[i]]
  for (k in seq_len(K)){
    j2 <- j2 + length(pw)
    gene_mat[j1:j2,1] <- k
    gene_mat[j1:j2,2] <- target_pws[i]
    gene_mat[j1:j2,3] <- pw
    gene_mat[j1:j2,4] <- kegg_risk1[pw,k]
    gene_mat[j1:j2,5] <- kegg_risk1_colsums[k]
    gene_mat[j1:j2,6] <- coef_k_unnorm[k]
    j1 <- j2+1
  }
}
gene_mat <- as.data.frame(gene_mat)
gene_mat$coef <- scale(as.numeric(gene_mat$coef))
gene_mat$count <- as.numeric(gene_mat$count)
gene_mat$offset <- log(as.numeric(gene_mat$offset))
gene_mat$log_count <- log(gene_mat$count+1)
gene_mat$ra <- gene_mat$count/exp(gene_mat$offset)
gene_mat$log_ra <- log((gene_mat$count+1)/exp(gene_mat$offset))


mm1 <- lmer(log_count ~  (1|pw) + (1|topic),data=gene_mat)
mm2 <- lmer(log_count ~  (1|pw) + (1|topic) + (1|topic:pw),data=gene_mat)
mm3 <- lmer(log_count ~  coef + (1|topic) + (coef|pw),data=gene_mat)
rmm2 <- ranef(mm2)


rmm3 <- ranef(mm3)
dotplot(ranef(mm3, condVar = TRUE), strip = FALSE)$pw

int_coef <- rmm2$`topic:pw`
int_coef_names <- str_split(rownames(int_coef),':')
int_mat <- matrix(0,length(target_pws),K,dimnames=list(target_pws,1:K))
for (i in seq_len(nrow(int_coef))){
  k <- int_coef_names[[i]][1]
  pw <- int_coef_names[[i]][2]
  int_mat[pw,k] <- int_coef[i,1]
}



dotplot(ranef(mm1, condVar = TRUE), strip = FALSE)$`topic:pw`




target_kegg1_counts <- t(sapply(target_kos, function(i) colSums(kegg_risk1[i,])))
target_kegg1 <- t(t(target_kegg1_counts)/kegg_risk1_colsums)







# 
# 
# make_beta_biom2 <- function(dat,m,dir_name){
#   fit <- dat$fits[[m]]
#   K <- fit$settings$dim$K
#   Nmin <- dat$filters$rare_min
#   
#   #beta <- t(floor(exp(do.call('rbind',fit$beta$logbeta))*Nmin))
#   beta <- t(exp(do.call('rbind',fit$beta$logbeta)))
#   for (i in 1:ncol(beta)) beta[order(beta[,i],decreasing=TRUE)[16:nrow(beta)],i] <- 0
#   for (i in 1:ncol(beta)) beta[order(beta[,i],decreasing=TRUE)[1:15],i] <- 15:1
#   beta_filt <- rowSums(beta) != 0
#   beta <- beta[beta_filt,]
#   #beta <- t(exp(do.call('rbind',fit$beta$logbeta))*Nmin)
#   
#   colnames(beta) <- paste0('T',1:NCOL(beta))
#   rownames(beta) <- as.character(dat$counts$ids[dat$vocab,'long'])[beta_filt]
#   
#   beta_meta <- data.frame(Topic=colnames(beta))
#   beta_taxa <- dat$taxa[beta_filt,]
#   
#   beta_filename <- paste0('beta_RANK.biom')
#   beta_biom <- make_biom(beta,beta_meta,beta_taxa)
#   write_biom(beta_biom,file.path(dir_name,beta_filename))
# }
# 
# 
# beta_name <- '~/MiscOut/TestBeta'
# dir.create(beta_name,showWarnings=FALSE,recursive=TRUE)
# make_beta_biom2(dat,5,beta_name)
# 
# 
# test <- t(floor(exp(do.call('rbind',fit$beta$logbeta))*1000))
# 
# 
# 
# biom_test <- read_biom('~/MiscOut/TestBeta/beta_RANK_predicted_metagenome.biom')
# kegg_test <- as.matrix(biom_data(biom_test))
# target_kegg1 <- kegg_test[,1:K]; target_kegg2 <- kegg_test[,(K+1):NCOL(kegg_test)]
# 
# # hold <- cbind(rowSums(kegg_test1),rowSums(kegg_risk1))
# 
# 
# 
# 
# 
# target_pws <- my_kegg_order[-15]
# target_kos <- lapply(target_pws,function(pw) names(which(sapply(kegg_metadata_risk_lev2, function(x) any(str_to_lower(x) %in% pw)))))
# 
# target_kegg1_sum <- colSums(target_kegg1)
# target_kegg2_sum <- colSums(target_kegg2)
# 
# target_kegg1 <- t(sapply(target_kos, function(i) colSums(target_kegg1[i,])))
# target_kegg2 <- t(sapply(target_kos, function(i) colSums(target_kegg2[i,])))
# 
# 
# target_kegg1 <- sapply(1:25,function(i) target_kegg1[,i]/target_kegg1_sum[i])
# target_kegg2 <- sapply(1:25,function(i) target_kegg2[,i]/target_kegg2_sum[i])
# 
# # target_kegg1 <- log10(t(t(target_kegg1)/target_kegg1_sum))
# # target_kegg2 <- log10(t(t(target_kegg2)/target_kegg2_sum))
# 
# # target_kegg1 <- t(sapply(target_kos, function(i) colSums(kegg_risk1[i,])))
# # target_kegg2 <- t(sapply(target_kos, function(i) colSums(kegg_risk2[i,])))
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# target_kegg1 <- t(apply(target_kegg1,1,qnormalize))
# 
# 
# 
# colors <- scales::muted(colorRampPalette(c('purple','black','orange'))(199))
# 
# hm_coef <- coef_k_unnorm
# names(hm_coef) <- 1:K
# hm_coef <- sort(hm_coef[sig_corr_topics],decreasing=TRUE)
# 
# 
# library(RColorBrewer)
# colors <- rev(colorRampPalette(brewer.pal(11,'PRGn'))(100))
# topic_pallete <- colorRampPalette(brewer.pal(11,'RdBu'))(length(sig_corr_topics) + 4)
# 
# topic_colors <- c(topic_pallete[1:sum(hm_coef>0)],rev(rev(topic_pallete)[1:sum(hm_coef<0)]))
# 
# pdf('~/MiscOut/gev_hm_oct1.pdf', height=20, width=20)
# heatmap3::heatmap3(target_kegg1[,order(coef_k_unnorm)],
#                    
#                    labRow=target_pws,
#                    #labCol='',
#                    revC=TRUE, 
#                    
#                    Colv=NA,
#                    Rowv=NA,                  
#                    
#                    col=colors,
#                    #balanceColor=FALSE,
#                    #ColSideColors=topic_colors,
#                    
#                    scale='row')
# dev.off()
# 
# 
# 
# 
# 
# 
# 





























colors <- scales::muted(colorRampPalette(c('purple','black','orange'))(199))

hm_coef <- coef_k_unnorm
names(hm_coef) <- 1:K

topics_choose <- sig_corr_topics
#topics_choose <- c(12,5,19,34,14,43,32,47,18,48,33,23,20)

hm_coef <- sort(hm_coef[topics_choose],decreasing=TRUE)


library(RColorBrewer)
colors <- rev(colorRampPalette(brewer.pal(11,'PRGn'))(100))
topic_pallete <- colorRampPalette(brewer.pal(11,'RdBu'))(length(topics_choose) + 4)

topic_colors <- c(topic_pallete[1:sum(hm_coef>0)],rev(rev(topic_pallete)[1:sum(hm_coef<0)]))

target_pws_abr <- c('C Met', 'AA Met', 'L Met', 'Motile',
                    'Xeno','2-Mtblts',
                    'E Met','Gen Met','Mem Trans','Cell P/S','GIP')


# pdf(file.path(out_dir,K,mod_dir,'hm_pw2_kegg1.pdf'), height=4, width=8)
# heatmap3::heatmap3(target_kegg1[,as.numeric(names(hm_coef))],
#                    
#                    labRow=target_pws_abr,
#                    #labCol='',
#                    revC=TRUE, 
#                    
#                    Colv=NA,
#                    Rowv=NA,                  
#                   
#                    lasCol=2,
#                    lasRow=2,
#                    
#                    cexRow=2,
#                    cexCol=1.75,
#                    
#                    legendfun=function() showLegend(legend=c('Low','High'),col=colors[c(1,length(colors))],cex=1.5),
#                    
#                    col=colors,
# 
#                    ColSideLabs='CD-Topic Association',
#                    ColSideColors=topic_colors,
#                    
#                    scale='row')
# dev.off()


pdf(file.path(out_dir,K,mod_dir,'hm_pw2_kegg1.pdf'), height=8, width=6)
heatmap3::heatmap3(t(target_kegg1[,as.numeric(names(hm_coef))]),
                   
                   labCol=target_pws,
                   revC=TRUE, 
                   
                   Colv=NA,
                   Rowv=NA,                  
                   
                   lasCol=2,
                   lasRow=2,
                   
                   cexCol=2,
                   cexRow=1.75,
                   
                   legendfun=function() showLegend(legend=c('Low','High'),col=colors[c(1,length(colors))],cex=1.5),
                   
                   col=colors,
                   
                   RowSideLabs='CD-Topic Association',
                   RowSideColors=topic_colors,
                   
#                    highlightCell=data.frame(row=c(1:6,10:13,11:12),
#                                             col=c(rep(2,10),rep(9,2)),
#                                             color=c(rep(c('cyan','orange2'),c(6,4)),rep('pink',2)),
#                                             lwd=4,stringsAsFactors=FALSE),
                   
                   scale='col')
dev.off()
rbind(cbind(rep(1:4),2,'yellow',10),cbind(rep(8:13),2,'orange',10),
      cbind(rep(1:7),3,'yellow',10),cbind(rep(8:9),3,'orange',10),
      cbind(rep(2:3),9,'black',10))



pdf('~/MiscOut/gev_hm_oct2.pdf', height=20, width=20)
heatmap3::heatmap3(target_kegg2[,as.numeric(names(hm_coef))],
                   
                   labRow=my_kegg_order,
                   #labCol='',
                   revC=TRUE, 
                   
                   Colv=NA,
                   Rowv=NA,                  
                   
                   col=colors,
                   #balanceColor=FALSE,
                   ColSideColors=topic_colors,
                   
                   scale='column')
dev.off()









TAX <- 6
TAXA <- dat$taxa
colnames(TAXA) <- c('Kingdom','Phylum','Class','Order','Family','Genus','Species')

B1 <- beta1[,as.numeric(names(hm_coef))]
B2 <- beta2[,as.numeric(names(hm_coef))]
B <- as.data.frame(cbind(B1,B2))

TAXON <- apply(TAXA[rownames(B),1:5], 1, function(x) paste0(x,collapse='.'))
TAXONP <- pretty_taxa_names(TAXON,'\\.',rownames(B))
TAXONPU <- unique(data.frame(TAXONP=TAXONP))
TAXONPU$OTU <- rownames(TAXONPU)

BB <- as.data.frame(t(B)) %>%
  mutate(Topic=factor(rownames(.),levels=paste0('T',as.numeric(names(hm_coef))))) %>%
  filter(Topic %in% paste0('T',c(12,5,14,47,23,20))) %>% # keep only some
  gather(OTU,Abundance,-Topic) %>%
  left_join(data.frame(TAXONP,OTU=names(TAXONP)),by='OTU') %>%
  group_by(Topic,TAXONP) %>%
  summarise(Abundance=sum(Abundance)) %>%
  ungroup() %>%
  filter(Abundance > 0) %>%
  group_by(Topic) %>%
  mutate(Frequency = Abundance/sum(Abundance)) %>%
  ungroup() %>%
  left_join(TAXONPU,by='TAXONP') %>%
  left_join(data.frame(TAXA[,1:TAX],OTU=rownames(TAXA)),by='OTU')

# 
# BBF <- BB %>%
#   mutate(Phylum=str_replace(Phylum,'p__',''),
#          Genus=str_replace(Genus,'g__','')) %>%
#   group_by(Phylum) %>%
#   mutate(GroupSum=sum(Abundance)) %>%
#   ungroup() %>%
#   filter(dense_rank(desc(GroupSum)) <= 5) %>%
#   group_by(Genus) %>%
#   mutate(GroupSum=sum(Abundance)) %>%
#   ungroup() %>%
#   filter(dense_rank(desc(GroupSum)) <= 10) %>%
#   ggplot(aes(x=as.factor(1),y=Frequency,fill=Genus)) + 
#   geom_bar(colour='black',stat='identity') +
#   facet_grid(Phylum ~ Topic) +
#   labs(x='',y='Phylum') +
#   guides(fill=guide_legend(nrow=3)) +
#   theme(axis.text.y=element_text(size=9),
#         strip.text.y=element_text(size=8),
#         strip.text.x=element_text(size=14,face='bold'),
#         axis.title=element_text(size=20,face='bold'),
#         axis.text.x=element_blank(),
#         axis.ticks.x=element_blank(),
#         title=element_text(size=25),
#         legend.position='bottom',
#         legend.text=element_text(size=12,face='bold'),
#         legend.title=element_text(size=20,face='bold'),
#         aspect.ratio=1.5/1) 



BBF <- BB %>%
  mutate(Phylum=str_replace(Phylum,'p__',''),
         Class=str_replace(Class,'c__',''),
         Order=str_replace(Order,'o__',''),
         Family=str_replace(Family,'f__',''),
         Genus=str_replace(Genus,'g__','')) %>%
  mutate(Genus=ifelse(Genus != '',Genus,
                      ifelse(Family != '',paste0(Family,'*'),
                             paste0(Order,'*')))) %>%
  group_by(Phylum) %>%
  mutate(GroupSum=sum(Abundance)) %>%
  ungroup() %>%
  filter(dense_rank(desc(GroupSum)) <= 5) %>%
  group_by(Genus) %>%
  mutate(GroupSum=sum(Abundance)) %>%
  ungroup() %>%
  filter(dense_rank(desc(GroupSum)) <= 10) %>%
  ggplot(aes(x=as.factor(1),y=Frequency,fill=Genus)) + 
  geom_bar(colour='black',stat='identity') +
  facet_grid(Topic ~ Phylum) +
  labs(y='',x='Phylum') +
  guides(fill=guide_legend(nrow=3)) +
  theme(axis.text.y=element_text(size=8),
        strip.text.x=element_text(size=8,angle=90,face='bold'),
        strip.text.y=element_text(size=12,face='bold',angle=90),
        axis.title=element_text(size=18,face='bold'),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position='right',
        legend.text=element_text(size=8,face='bold'),
        legend.title=element_text(size=18,face='bold'),
        aspect.ratio=1.5/1) +
  guides(fill=guide_legend(ncol=1,keyheight=.75, #keywidth = .75,
                           label.position='bottom'))


pdf(file.path(out_dir,K,mod_dir,'taxa_beta1.pdf'), height=5, width=3.75)
BBF
dev.off()










order(coef_k_unnorm,decreasing=TRUE)
sageLabels(fit,n=5)

fit_corr <- topicCorr(fit,method='huge',.1)


which(fit_corr$cor[7,] > .1)
as.character(dat$counts$ids[fit$vocab[which(exp(fit$beta$logbeta[[1]][7,]) > .01)],'long'])
as.character(dat$counts$ids[fit$vocab[which(exp(fit$beta$logbeta[[1]][17,]) > .01)],'long'])
as.character(dat$counts$ids[fit$vocab[which(exp(fit$beta$logbeta[[1]][20,]) > .01)],'long'])

plot(fit_corr)


OTU <- dat$counts$table_clean
colnames(OTU) <- as.character(dat$counts$ids[colnames(OTU),'long'])
META <- dat$counts$meta_clean
TAXA <- dat$taxa
colnames(TAXA) <- c('Kingdom','Phylum','Class','Order','Family','Genus','Species')
PS <- phyloseq(otu_table(OTU,taxa_are_rows=FALSE),sample_data(META),tax_table(TAXA))


library(phyloseq)
## Load round 2 of American gut project
SPEZ <- spiec.easi(PS, method='mb', lambda.min.ratio=.05,
                   nlambda=10, icov.select.params=list(rep.num=50))


SPEZ_MB <- adj2igraph(SPEZ$refit, rmEmptyNodes=TRUE,  vertex.attr=list(name=taxa_names(PS)))
plot_network(SPEZ_MB, PS, type='taxa', color='Phylum',line_weight=.75)


SPEZ_cor <- as.matrix(SPEZ$beta[[SPEZ$opt.index]])
dimnames(SPEZ_cor) <- list(colnames(OTU),colnames(OTU))










target_pws <- names(which(sort(table(unlist(kegg_metadata_risk_lev3))) >= 25))
target_kos <- lapply(target_pws,function(pw) names(which(sapply(kegg_metadata_risk_lev3, function(x) any(x %in% pw)))))



target_kegg1_sum <- colSums(target_kegg1)
target_kegg2_sum <- colSums(target_kegg2)

target_kegg1 <- t(sapply(target_kos, function(i) colSums(target_kegg1[i,])))
target_kegg2 <- t(sapply(target_kos, function(i) colSums(target_kegg2[i,])))


target_kegg1 <- sapply(1:25,function(i) target_kegg1[,i]/target_kegg1_sum[i])
target_kegg2 <- sapply(1:25,function(i) target_kegg2[,i]/target_kegg2_sum[i])



# target_kegg1 <- log10(t(sapply(target_kos, function(i) colSums(kegg_risk1[i,])))+1)
# target_kegg2 <- log10(t(sapply(target_kos, function(i) colSums(kegg_risk2[i,])))+1)


row_drop1 <- which(rowSums(target_kegg1) == 0)
row_drop2 <- which(rowSums(target_kegg2) == 0)



colors <- scales::muted(colorRampPalette(c('purple','black','orange'))(199))

hm_coef <- coef_k_unnorm
names(hm_coef) <- 1:K
hm_coef <- sort(hm_coef[sig_corr_topics],decreasing=TRUE)


library(RColorBrewer)
colors <- rev(colorRampPalette(brewer.pal(11,'PRGn'))(100))
topic_pallete <- colorRampPalette(brewer.pal(11,'RdBu'))(length(sig_corr_topics) + 4)

topic_colors <- c(topic_pallete[1:sum(hm_coef>0)],rev(rev(topic_pallete)[1:sum(hm_coef<0)]))


dismat1 <- vegdist(target_kegg1[-row_drop1,as.numeric(names(hm_coef))],method='bray')
rc <- hclust(dismat1,'ward.D2')
roword <- as.dendrogram(cc)




pdf('~/MiscOut/gev_hm_oct1.pdf', height=20, width=20)
heatmap3::heatmap3(target_kegg1[-row_drop1,as.numeric(names(hm_coef))],
                   
                   labRow=target_pws[-row_drop1],
                   #labCol='',
                   revC=TRUE, 
                   
                   Colv=NA,
                   Rowv=roword,                  
                   
                   col=colors,
                   #balanceColor=FALSE,
                   ColSideColors=topic_colors,
                   
                   scale='row')
dev.off()




pdf('~/MiscOut/gev_hm_oct1.pdf', height=20, width=20)
heatmap3::heatmap3(target_kegg1[-row_drop1,order(coef_k_unnorm)],
                   
                   labRow=target_pws[-row_drop1],
                   #labCol='',
                   revC=TRUE, 
                   
                   Colv=NA,
                   Rowv=roword,                  
                   
                   col=colors,
                   #balanceColor=FALSE,
                   #ColSideColors=topic_colors,
                   
                   scale='row')
dev.off()





dismat1 <- vegdist(target_kegg2[-row_drop2,as.numeric(names(hm_coef))],method='bray')
rc <- hclust(dismat1,'ward.D2')
roword <- as.dendrogram(cc)


pdf('~/MiscOut/gev_hm_oct2.pdf', height=20, width=20)
heatmap3::heatmap3(target_kegg2[-row_drop2,as.numeric(names(hm_coef))],
                   
                   labRow=target_pws[-row_drop2],
                   #labCol='',
                   revC=TRUE, 
                   
                   Colv=NA,
                   Rowv=roword,                  
                   
                   col=colors,
                   #balanceColor=FALSE,
                   ColSideColors=topic_colors,
                   
                   scale='row')
dev.off()






target_pws <- names(which(sort(table(unlist(kegg_metadata_risk_lev3))) >= 25))
target_kos <- lapply(target_pws,function(pw) names(which(sapply(kegg_metadata_risk_lev3, function(x) any(x %in% pw)))))

test <- t(t(kegg_risk1)/colSums(kegg_risk1))
target_kegg1 <- log10(t(sapply(target_kos, function(i) colSums(test[i,])))+1)




dismat1 <- vegdist(target_kegg1[-row_drop1,order(coef_k_unnorm,decreasing=TRUE)],method='bray')
rc <- hclust(dismat1,'ward.D2')
roword <- as.dendrogram(cc)




pdf('~/MiscOut/gev_hm_oct1.pdf', height=20, width=20)
heatmap3::heatmap3(target_kegg1[-row_drop1,order(coef_k_unnorm,decreasing=TRUE)],
                   
                   labRow=target_pws[-row_drop1],
                   #labCol='',
                   revC=TRUE, 
                   
                   Colv=NA,
                   Rowv=roword,                  
                   
                   col=colors,
                   #balanceColor=FALSE,
                   #ColSideColors=topic_colors,
                   
                   scale='row')
dev.off()









eff <- estimateEffect(1:K ~ FEMALE + s(AGE), fit, meta, uncertainty='Global')

plot.estimateEffect(eff, covariate = "AGE", model = fit, 
                    method = "continuous", xlab = "AGE", moderator = "FEMALE",
                    moderator.value = 1, linecol = "blue", ylim = c(0, .15),
                    printlegend = F)
plot.estimateEffect(eff, covariate = "AGE", model = fit,
                    method = "continuous", xlab = "Days", moderator = "FEMALE",
                    moderator.value = 0, linecol = "red", add = T,
                    printlegend = F)
legend(0, .1, c("female", "male"), lwd = 2, col = c("blue", "red"))




k <- c(31)
plot.estimateEffect(eff, covariate = "AGE", model = fit, topics=k,
                    method = "continuous", xlab = "AGE", moderator = "FEMALE",
                    moderator.value = 1, linecol = "blue", ylim = c(0, .15),
                    printlegend = F)
plot.estimateEffect(eff, covariate = "AGE", model = fit, topics=k,
                    method = "continuous", xlab = "Days", moderator = "FEMALE",
                    moderator.value = 0, linecol = "red", add = T,
                    printlegend = F)
legend(0, .1, c("female", "male"), lwd = 2, col = c("blue", "red"))




labelTopics(fit,5)
labelTopics(fit,c(5,17))
plot.STM(fit,type='perspectives',topics=c(5))
plot.STM(fit,type='perspectives',topics=c(17))
plot.STM(fit,type='perspectives',topics=c(5,17))


tt <- fit$beta$logbeta[[1]][17,]
tt_idx <- which(tt > quantile(tt,.95))
tt_val <- exp(tt[tt_idx])
tt_otu <- paste0('otu',tt_idx)
tt_taxa <- pretty_taxa_names(apply(dat$taxa[dat$counts$ids[tt_otu,'long'],1:4],1,function(x) paste0(x,collapse='|')))

data.frame(Order=tt_taxa,Abundance=tt_val) %>%
  arrange(desc(Abundance)) %>%
  mutate(idx=row_number()) %>%
  ggplot(aes(idx,Abundance,fill=Order)) + geom_bar(colour='black',position='dodge',stat='identity') + xlab('')




z <- dat$counts$table_clean
z <- z[dat$counts$meta_clean %>% filter(diet == 'CD') %>% dplyr::select(SampleID) %>% unlist(),]
z <- z[,colSums(z) > 0]

b <- fit$beta$logbeta[[1]][17,]
names(b) <- paste0('otu',1:length(b))
o <- names(which(b > quantile(b,.95)))
length(o)

z <- z[,o]

df <- data.frame(z,dat$counts$meta_clean %>% filter(diet == 1) %>% dplyr::select(SampleID,age,age_RANK)) %>%
  gather(OTU,Abundance,-SampleID,-age,-age_RANK) %>%
  group_by(OTU) %>%
  mutate(Prevalence = sum(Abundance > 0)) %>%
  ungroup() %>%
  filter(Prevalence >= 35) %>%
  mutate(Abundance = Abundance + 1) %>%
  rowwise() %>%
  mutate(Taxon = pretty_taxa_names(paste0(dat$taxa[as.character(dat$counts$ids[unique(OTU),'long']),],collapse='|')),
         Facet = paste0(Taxon,' (',OTU,')'))
length(unique(df$OTU))


library(purrr)
pvals <- df %>%
  mutate(ALOG10 = log10(Abundance)) %>%
  split(.$OTU) %>%
  map(~ coef(summary(glm(Abundance ~ age, data=., family=quasipoisson(link='log'))))[2,]) %>%
  do.call('rbind',.)
pvals <- as.data.frame(pvals)
pvals$adj <- p.adjust(pvals[,4],method='BH')
pvals_sig <- rownames(pvals)[pvals$adj < .05]


pdf('~/MiscOut/testing_formanu.pdf', height=12, width=15)
df %>%
  mutate(Significant = as.factor(ifelse(OTU %in% pvals_sig,1,0))) %>%
  ggplot(aes(age,Abundance,colour=Significant)) + 
  geom_point(data=dplyr::select(df,-Facet),aes(age,Abundance),colour='gray',alpha=.5,size=1) +
  facet_wrap(~Facet) + 
  geom_point(alpha=1,size=2) +
  scale_y_log10() +
  scale_colour_manual(values=c('black','red')) + 
  stat_smooth(method='lm',se=FALSE,colour='darkgreen') +
  theme_bw() + theme(legend.position='none') + xlab('age')
dev.off()





z <- dat$counts$table_clean
z <- z[dat$counts$meta_clean %>% dplyr::select(SampleID) %>% unlist(),]
z <- z[,colSums(z) > 0]

b <- fit$beta$logbeta[[1]][17,]
names(b) <- paste0('otu',1:length(b))

o <- labelTopics(fit,c(17),n=5)$interaction

z <- z[,c(o[1,],o[2,])]

df <- data.frame(z,dat$counts$meta_clean %>% dplyr::select(SampleID,diet)) %>%
  gather(OTU,Abundance,-SampleID,-diet) %>%
  mutate(Score=ifelse(OTU %in% o[1,],'CD+','CD-')) %>%
  group_by(OTU) %>%
  mutate(Prevalence = sum(Abundance > 0)) %>%
  ungroup() %>%
  filter(Prevalence > 15) %>%
  #mutate(Abundance = Abundance + 1) %>%
  rowwise() %>%
  mutate(Taxon = pretty_taxa_names(paste0(dat$taxa[as.character(dat$counts$ids[unique(OTU),'long']),],collapse='|')),
         Facet = paste0(Taxon,' (',OTU,')'))
length(unique(df$OTU))

#pdf('~/MiscOut/testing1.pdf', height=12, width=15)
ggplot(df, aes(age_RANK,Abundance,colour=Facet)) + 
  geom_point(data=dplyr::select(df,-Facet),aes(age_RANK,Abundance),colour='gray',alpha=.5,size=1) +
  facet_wrap(~Facet) + 
  geom_point(alpha=1,size=2,colour='black') +
  scale_y_log10() +
  stat_smooth(method='lm',se=FALSE,colour='red') +
  theme_bw() + xlab('age Rank')
#dev.off()










age <- age_id$age
qplot(1:length(age_id$age),dat$counts$table_clean[age_id$SampleID,'otu22'] + 1, colour=age) + scale_y_log10() +
  xlab('age Rank') + ylab('Abundance') + ggtitle('T5: Enterobacteriaceae OTU')

dat$taxa[as.character(dat$counts$ids[c('otu356','otu270','otu21'),]$long),]
dat$taxa[as.character(dat$counts$ids[c('otu190','otu201'),]$long),]


data.frame(dat$counts$table_clean[,c('otu356','otu270','otu21','otu190','otu201')],diet=dat$counts$meta_clean$diet,age=dat$counts$meta_clean$age) %>%
  gather(OTU,Abundance,-age,-diet) %>%
  ggplot(aes(Abundance,)) + geom_point() + facet_wrap(~diet) + scale_y_log10()



data.frame(th,age=age_id$age,Rank=1:NROW(th)) %>%
  gather(Topic,Abundance,-Rank,-age) %>%
  mutate(Topic=factor(Topic,levels=colnames(th))) %>%
  filter(Topic == 'T17') %>%
  ggplot(aes(Rank,Abundance,colour=age)) + geom_point(alpha=.5) + scale_y_log10() +
  stat_smooth(method='lm',se=FALSE,colour='red') + xlab('age Rank')


th <- fit$theta[,17]
bugs <- dat$counts$table_clean[meta$SampleID,paste0('otu',which(b1 > quantile(b1,.8)))]
data.frame(T17=th,RANK=dat$meta$age_RANK,bugs=rowSums(bugs)) %>%
  ggplot(aes(RANK,bugs)) + geom_point() 





b1 <- exp(fit$beta$logbeta[[1]][17,])
paste0('otu',which(b1 > quantile(b1,.8)))


b1 <- c('otu174', 'otu11', 'otu356', 'otu138', 'otu561', 'otu502', 'otu578' )
b2 <- c('otu411', 'otu505', 'otu144', 'otu209', 'otu250', 'otu272', 'otu24')
b3 <- c('otu327', 'otu517', 'otu262', 'otu92', 'otu470', 'otu73', 'otu496')
b4 <- c('otu356','otu270','otu21','otu190','otu201')


checkBeta(fit) 
checkResiduals(fit, dat$docs)
sage <- sageLabels(fit,n=5)

# output_dir <- file.path(out_dir,K,mod_dir,'string')
# dir.create(output_dir,showWarnings=FALSE,recursive=TRUE)
# cog_dir <- file.path(dir_name,cog_filename)
# string_dump3(output_dir,cog_dir,rand=534,min_thres=.001)







filt <- rowSums(kegg_risk1) >= 20 & rowSums(kegg_risk2) >= 13
kegg1 <- kegg_risk1[filt,]
kegg2 <- kegg_risk2[filt,]
logbeta <- list(log(ra_norm(kegg1+1,FALSE)),log(ra_norm(kegg2+1,FALSE))) ### added psuedocounts
wordcounts <- cbind(rowSums(kegg1),rowSums(kegg2))

frex1 <- calcfrex_kegg(logbeta, coef_k, kegg_gene_names, n=10, w=.5)
frex2 <- calcfrex_kegg(logbeta, coef_k, kegg_gene_names, n=10, w=.5, wordcounts=apply(wordcounts,1,median))

frex3 <- lapply(seq_along(logbeta), function(i) calcfrex_kegg(logbeta[i], coef_k, kegg_gene_names, n=25, w=.5))
frex4 <- lapply(seq_along(logbeta), function(i) calcfrex_kegg(logbeta[i], coef_k, kegg_gene_names, n=25, w=.5, wordcounts=wordcounts[,i]))







rakegg1 <- ra_norm(kegg1+1,FALSE)
rakegg2 <- ra_norm(kegg2+1,FALSE)

k <- 3
n <- 15


genes_idx <- order(rakegg1[k,]-rakegg2[k,],decreasing=TRUE)
genes_idx <- c(genes_idx[1:floor(n/2)],rev(genes_idx)[1:ceiling(n/2)])

f_diff <- rakegg1[k,genes_idx] - rakegg2[k,genes_idx]
f_diff <- f_diff[sample(seq_along(f_diff))]
f1 <- sort(rakegg1[k,],decreasing=TRUE)[1:n]
f2 <- sort(rakegg2[k,],decreasing=TRUE)[1:n]


words_diff <- str_replace(unlist(kegg_gene_names[names(f_diff)]),'\\s\\[.*','')
words1 <- str_replace(unlist(kegg_gene_names[names(f1)]),'\\s\\[.*','')
words2 <- str_replace(unlist(kegg_gene_names[names(f2)]),'\\s\\[.*','')

data.frame(words=c(words_diff,words1,words2),
           f=c(f_diff,f1,f2),
           color=c(ifelse(f_diff>0,'female','male'),rep(c('female','male'),each=n)),
           group=rep(c('DIFF','female','male'),each=n)) %>%
  mutate(group=factor(group,levels=c('female','DIFF','male'))) %>%
  group_by(group) %>%
  mutate(y=sample(1:n())) %>%
  ggplot(aes(x=group,y=y,label=words,colour=color,size=f,alpha=f)) + 
  geom_text(position=position_jitter(width=.4,height=1.1),fontface='bold') + 
  scale_size(range = c(4, 10)) + scale_alpha(range=c(.3,1)) +
  theme_bw() +
  theme(legend.position='none')

