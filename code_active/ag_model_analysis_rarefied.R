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


source('~/Dropbox/stm_microbiome/code_active/stm_functions.R')
source('~/Dropbox/stm_microbiome/code_active/nav_froz_fxns_3.R')


rare_min <- 10000
K <- 35
m <- 4 ### which fit?
M <- 4


dir_name <- sprintf('~/Dropbox/stm_microbiome/data_active/AG_female/stm_s97_rarefied_%s_supervised/model',rare_min)
fit_filename <- paste0('fits_K_',K,'_nfits_',M,'.rds')
beta_filename <- paste0('beta_K_',K,'_nfits_',M,'_fit_',m,'.biom')
kegg_filename <- paste0('beta_K_',K,'_nfits_',M,'_fit_',m,'_predicted_metagenome.biom')
cog_filename <- paste0('beta_K_',K,'_nfits_',M,'_fit_',m,'_predicted_cogs.biom')

dat <- readRDS(file.path(dir_name,fit_filename))

fit <- dat$fits[[m]]
meta <- dat$meta

fit$settings$call

out_dir <- gsub('data_active','output_active',dir_name)
mod_dir <- paste0('mod_',m)
heatmap_dir <- file.path(out_dir,K,mod_dir,'heatmaps')
dir.create(heatmap_dir,showWarnings=FALSE,recursive=TRUE)

sink(file.path(out_dir,K,mod_dir,'model.txt')); fit$settings$call; cat('\nK = ',fit$settings$dim$K,'\n'); sink()


# eff <- estimateEffect(1:K ~ FEMALE + AGE_DECADE_RANK, fit, meta, uncertainty='None')
# eff <- estimateEffect(1:K ~ FEMALE + s(AGE), fit, meta, uncertainty='None')
eff <- estimateEffect(1:K ~ FEMALE, fit, meta, uncertainty='None')
eff_plot_female <- plot.estimateEffect(eff, 'FEMALE', model=fit, topics=1:K, method='difference',cov.value1=1,cov.value2=0)
eff <- estimateEffect(1:K ~ s(AGE), fit, meta, uncertainty='None')
eff_plot_age <- plot.estimateEffect(eff, 'AGE', model=fit, topics=1:K, method='continuous')

coef_k <- unlist(eff_plot_female$means)
coef_k <- coef_k/sd(coef_k)

coef_k_age <- sapply(eff_plot_age$means,function(x) max(x[20:80])-min(x[20:80]))
coef_k_age <- coef_k_age/sd(coef_k_age)

coef_k_age <- coef_k_age * coef_k
coef_k_age <- coef_k_age/sd(coef_k_age)


#eff_plot_age <- plot.estimateEffect(eff, 'AGE', model=fit, topics=1:K, method='continuous')


dat$meta$MILESTONE <- ifelse(dat$meta$AGE <= 10, 'Prepubescent',
                         ifelse(dat$meta$AGE > 10 & dat$meta$AGE <= 17, 'Pubescent',
                                ifelse(dat$meta$AGE > 17 & dat$meta$AGE <= 54, 'Adult',
                                       ifelse(dat$meta$AGE > 54 & dat$meta$AGE <= 69, 'Postmenopausal',
                                              'Elderly'))))
dat$meta$MILESTONE <- factor(dat$meta$MILESTONE,levels=c('Prepubescent','Pubescent','Adult','Postmenopausal','Elderly'))


####################################
####################################
mat_b <- fit$theta
meta <- dat$meta
rownames(mat_b) <- meta$SampleID


meta <- meta %>%
  mutate(AGE_YEARS=as.numeric(AGE_YEARS)) %>%
  arrange(desc(AGE_YEARS))
mat_b <- mat_b[meta$SampleID,]
rownames(mat_b) <- meta$AGE_DECADE_RANK
colnames(mat_b) <- 1:NCOL(mat_b)
####################################
####################################


####################################
####################################
mat_f <- fit$theta
meta <- dat$meta
rownames(mat_f) <- meta$SampleID


meta <- meta %>%
  filter(SEX=='female') %>%
  mutate(AGE_YEARS=as.numeric(AGE_YEARS)) %>%
  arrange(desc(AGE_YEARS))
mat_f <- mat_f[meta$SampleID,]
rownames(mat_f) <- meta$AGE_DECADE_RANK
colnames(mat_f) <- 1:NCOL(mat_f)
####################################
####################################



####################################
####################################
mat_m <- fit$theta
meta <- dat$meta
rownames(mat_m) <- meta$SampleID


meta <- meta %>%
  filter(SEX=='male') %>%
  mutate(AGE_YEARS=as.numeric(AGE_YEARS)) %>%
  arrange(desc(AGE_YEARS))
mat_m <- mat_m[meta$SampleID,]
rownames(mat_m) <- meta$AGE_DECADE_RANK
colnames(mat_m) <- 1:NCOL(mat_m)
####################################
####################################


####################################
####################################
mat <- fit$theta
meta <- dat$meta
rownames(mat) <- meta$SampleID


meta <- meta %>%
  mutate(AGE_YEARS=as.numeric(AGE_YEARS)) %>%
  arrange(SEX,desc(AGE_YEARS))
mat <- mat[meta$SampleID,]
rownames(mat) <- meta$SEX
colnames(mat) <- 1:NCOL(mat)
####################################
####################################



# prep data for taxon table heatmap
####################################
####################################
biom_file <- read_biom(file.path(dir_name,beta_filename))
beta <- as.matrix(biom_data(biom_file))
beta_taxa <- observation_metadata(biom_file) 
if (NCOL(beta) == K) {beta1 <- beta; beta2 <- NULL} else {beta1 <- beta[,1:K]; beta2 <- beta[,(K+1):NCOL(beta)]}


taxon <- 7 # species
beta <- beta1
taxon_names <- apply(beta_taxa[rownames(beta),1:taxon],1,paste0,collapse='|')
taxon_table <- data.frame(beta,taxon_names) %>%
  group_by(taxon_names) %>%
  summarise_each(funs(sum)) 
taxon_names <- pretty_taxa_names(taxon_table$taxon_names)
taxon_table <- data.frame(dplyr::select(taxon_table,-taxon_names))
rownames(taxon_table) <- taxon_names
colnames(taxon_table) <- 1:NCOL(taxon_table)

filt_idx <- rowSums(taxon_table) != 0
taxon_table1 <- taxon_table[filt_idx,]

if (!is.null(beta2)){
  beta <- beta2
  taxon_names <- apply(beta_taxa[rownames(beta),1:taxon],1,paste0,collapse='|')
  taxon_table <- data.frame(beta,taxon_names) %>%
    group_by(taxon_names) %>%
    summarise_each(funs(sum)) 
  taxon_names <- pretty_taxa_names(taxon_table$taxon_names)
  taxon_table <- data.frame(dplyr::select(taxon_table,-taxon_names))
  rownames(taxon_table) <- taxon_names
  colnames(taxon_table) <- 1:NCOL(taxon_table)
  
  filt_idx <- rowSums(taxon_table) != 0
  taxon_table2 <- taxon_table[filt_idx,]
}
####################################
####################################


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
  mutate(assoc=ifelse(coef>0,'F',ifelse(coef<0,'M','ZERO')),
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

misc_filter <- !(rownames(kegg_sig) %in% c('None','hypothetical protein'))
kegg_sig1 <- kegg_sig[misc_filter,]
ko_names1 <- ko_names[misc_filter]

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
    mutate(assoc=ifelse(coef>0,'F',ifelse(coef<0,'M','ZERO')),
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
  
  misc_filter <- !(rownames(kegg_sig) %in% c('None','hypothetical protein'))
  kegg_sig2 <- kegg_sig[misc_filter,]
  ko_names2 <- ko_names[misc_filter]
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
  mutate(assoc=ifelse(coef>0,'F',ifelse(coef<0,'M','ZERO')),
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

misc_filter <- !(rownames(cog_sig) %in% c('None','hypothetical protein','Uncharacterized conserved protein','Uncharacterized conserved protein'))
cog_sig1 <- cog_sig[misc_filter,]
cog_names1 <- cog_names[misc_filter]

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
    mutate(assoc=ifelse(coef>0,'F',ifelse(coef<0,'M','ZERO')),
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
  
  misc_filter <- !(rownames(cog_sig) %in% c('None','hypothetical protein','Uncharacterized conserved protein','Uncharacterized conserved protein'))
  cog_sig2 <- cog_sig[misc_filter,]
  cog_names2 <- cog_names[misc_filter]   
}
####################################
####################################














pdf(file.path(heatmap_dir,'heatmap_theta.pdf'), height=20, width=20)
plot_big_heatmap3(mat,coef_k,col1='cyan',col2='orange',normtype='q',normaxis=2,rowlabs=FALSE,rowclust=FALSE,rowcols=TRUE,
                  ColSideLabs='Topic Weight',RowSideLabs='Sex',main='Samples over Topics',
                  cexCol=1.5,cexRow=.75)
dev.off()

pdf(file.path(heatmap_dir,'heatmap_theta_age.pdf'), height=20, width=20)
plot_big_heatmap3(mat_b,coef_k_age,col1='cyan',col2='orange',normtype='q',normaxis=2,rowlabs=FALSE,rowclust=FALSE,rowcols=TRUE,
                  ColSideLabs='Topic Weight',RowSideLabs='Sex',main='Samples over Topics',
                  cexCol=1.5,cexRow=.75)
dev.off()

pdf(file.path(heatmap_dir,'heatmap_theta_age_f.pdf'), height=20, width=20)
plot_big_heatmap3(mat_f,coef_k_age,col1='cyan',col2='orange',normtype='q',normaxis=2,rowlabs=FALSE,rowclust=FALSE,rowcols=TRUE,
                  ColSideLabs='Topic Weight',RowSideLabs='Sex',main='Samples over Topics',
                  cexCol=1.5,cexRow=.75)
dev.off()

pdf(file.path(heatmap_dir,'heatmap_theta_age_m.pdf'), height=20, width=20)
plot_big_heatmap3(mat_m,coef_k_age,col1='cyan',col2='orange',normtype='q',normaxis=2,rowlabs=FALSE,rowclust=FALSE,rowcols=TRUE,
                  ColSideLabs='Topic Weight',RowSideLabs='Sex',main='Samples over Topics',
                  cexCol=1.5,cexRow=.75)
dev.off()


pdf(file.path(heatmap_dir,'heatmap_beta1.pdf'), height=20, width=20)
hm_beta1 <- plot_big_heatmap3(taxon_table1,coef_k,col1='cyan',col2='orange',normtype=NULL,normaxis=1,rowlabs=FALSE,rowclust=TRUE,rowcols='red',
                              dis='bray',scale='row',
                              ColSideLabs='Topic Weight',RowSideLabs='Sex',main='Topics over OTUs',
                              cexCol=1.5,cexRow=.75)
dev.off()

pdf(file.path(heatmap_dir,'heatmap_kegg1.pdf'), height=20, width=20)
plot_big_heatmap3(kegg_sig1,coef_k,col1='cyan',col2='orange',normtype=NULL,normaxis=1,rowlabs=FALSE,rowclust=TRUE,rowcols='red',
                  dis='bray',scale='row',
                  ColSideLabs='Topic Weight',RowSideLabs='Sex',main='Topics over KOs',
                  cexCol=1.5,cexRow=.75)
dev.off()

pdf(file.path(heatmap_dir,'heatmap_kegg_rownames1.pdf'), height=50, width=75)
hm_kegg1 <- plot_big_heatmap3(kegg_sig1,coef_k,col1='cyan',col2='orange',normtype=NULL,normaxis=1,rowlabs=TRUE,rowclust=TRUE,rowcols='red',
                              dis='bray',scale='row',
                              ColSideLabs='Topic Weight',RowSideLabs='Sex',main='Topics over KOs',
                              cexCol=5,cexRow=4,margins=c(20,50))
dev.off()

pdf(file.path(heatmap_dir,'heatmap_cog1.pdf'), height=20, width=20)
plot_big_heatmap3(cog_sig1,coef_k,col1='cyan',col2='orange',normtype=NULL,normaxis=1,rowlabs=FALSE,rowclust=TRUE,rowcols='red',
                  dis='bray',scale='row',
                  ColSideLabs='Topic Weight',RowSideLabs='Sex',main='Topics over COGs',
                  cexCol=1.5,cexRow=.75)
dev.off()

pdf(file.path(heatmap_dir,'heatmap_cog_rownames1.pdf'), height=60, width=75)
hm_cog1 <- plot_big_heatmap3(cog_sig1,coef_k,col1='cyan',col2='orange',normtype=NULL,normaxis=1,rowlabs=TRUE,rowclust=TRUE,rowcols='red',
                             dis='bray',scale='row',
                             ColSideLabs='Topic Weight',RowSideLabs='Sex',main='Topics over COGs',
                             cexCol=5,cexRow=4,margins=c(20,20))
dev.off()




side_names <- collapse_to_my_kegg_names(kegg_metadata_risk_lev2[names(ko_names1)], my_kegg_order)
side_names_ord <- order(side_names)
side_names <- side_names[side_names_ord]
side_idx <- as.integer(as.factor(side_names))

pdf(file.path(heatmap_dir,'heatmap_kegg_pathway1.pdf'), height=20, width=20)
plot_big_heatmap3(kegg_sig1[side_names_ord,],
                  coef_k,col1='cyan',col2='orange',normtype=NULL,normaxis=1,rowlabs=FALSE,rowclust=FALSE,rowcols=side_names,
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
  hm_beta2 <- plot_big_heatmap3(taxon_table2,coef_k,col1='cyan',col2='orange',normtype=NULL,normaxis=1,rowlabs=FALSE,rowclust=TRUE,rowcols='blue',
                                dis='bray',scale='row',
                                ColSideLabs='Topic Weight',RowSideLabs='Sex',main='Topics over OTUs',
                                cexCol=1.5,cexRow=.75)
  dev.off()
  
  pdf(file.path(heatmap_dir,'heatmap_kegg2.pdf'), height=20, width=20)
  plot_big_heatmap3(kegg_sig2,coef_k,col1='cyan',col2='orange',normtype=NULL,normaxis=1,rowlabs=FALSE,rowclust=TRUE,rowcols='blue',
                    dis='bray',scale='row',
                    ColSideLabs='Topic Weight',RowSideLabs='Sex',main='Topics over KOs',
                    cexCol=1.5,cexRow=.75)
  dev.off()
  
  pdf(file.path(heatmap_dir,'heatmap_kegg_rownames2.pdf'), height=50, width=75)
  hm_kegg2 <- plot_big_heatmap3(kegg_sig2,coef_k,col1='cyan',col2='orange',normtype=NULL,normaxis=1,rowlabs=TRUE,rowclust=TRUE,rowcols='blue',
                                dis='bray',scale='row',
                                ColSideLabs='Topic Weight',RowSideLabs='Sex',main='Topics over KOs',
                                cexCol=5,cexRow=4,margins=c(20,20))
  dev.off()
  
  pdf(file.path(heatmap_dir,'heatmap_cog2.pdf'), height=20, width=20)
  plot_big_heatmap3(cog_sig2,coef_k,col1='cyan',col2='orange',normtype=NULL,normaxis=1,rowlabs=FALSE,rowclust=TRUE,rowcols='blue',
                    dis='bray',scale='row',
                    ColSideLabs='Topic Weight',RowSideLabs='Sex',main='Topics over COGs',
                    cexCol=1.5,cexRow=.75)
  dev.off()
  
  pdf(file.path(heatmap_dir,'heatmap_cog_rownames2.pdf'), height=60, width=75)
  hm_cog2 <- plot_big_heatmap3(cog_sig2,coef_k,col1='cyan',col2='orange',normtype=NULL,normaxis=1,rowlabs=TRUE,rowclust=TRUE,rowcols='blue',
                               dis='bray',scale='row',
                               ColSideLabs='Topic Weight',RowSideLabs='Sex',main='Topics over COGs',
                               cexCol=5,cexRow=4,margins=c(20,20))
  dev.off()
  
  side_names <- collapse_to_my_kegg_names(kegg_metadata_risk_lev2[names(ko_names2)], my_kegg_order)
  side_names_ord <- order(side_names)
  side_names <- side_names[side_names_ord]
  side_idx <- as.integer(as.factor(side_names))
  
  
  pdf(file.path(heatmap_dir,'heatmap_kegg_pathway2.pdf'), height=20, width=20)
  plot_big_heatmap3(kegg_sig2[side_names_ord,],
                    coef_k,col1='cyan',col2='orange',normtype=NULL,normaxis=1,rowlabs=FALSE,rowclust=FALSE,rowcols=side_names,
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
  plot_big_heatmap3(group_mat,coef_k,col1='cyan',col2='orange',normtype=NULL,normaxis=1,rowlabs=TRUE,rowclust=FALSE,rowcols='red',
                    dis='bray',scale='row',
                    ColSideLabs='Topic Weight',RowSideLabs='Sex',main='Topics over OTUs',
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
  plot_big_heatmap3(group_mat,coef_k,col1='cyan',col2='orange',normtype=NULL,normaxis=1,rowlabs=TRUE,rowclust=FALSE,rowcols=side_names,
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
  plot_big_heatmap3(group_mat,coef_k,col1='cyan',col2='orange',normtype=NULL,normaxis=1,rowlabs=TRUE,rowclust=FALSE,rowcols='red',
                    dis='bray',scale='row',
                    ColSideLabs='Topic Weight',RowSideLabs='Sex',main='Topics over COGs',
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
    plot_big_heatmap3(group_mat,coef_k,col1='cyan',col2='orange',normtype=NULL,normaxis=1,rowlabs=TRUE,rowclust=FALSE,rowcols='blue',
                      dis='bray',scale='row',
                      ColSideLabs='Topic Weight',RowSideLabs='Sex',main='Topics over OTUs',
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
    plot_big_heatmap3(group_mat,coef_k,col1='cyan',col2='orange',normtype=NULL,normaxis=1,rowlabs=TRUE,rowclust=FALSE,rowcols=side_names,
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
    plot_big_heatmap3(group_mat,coef_k,col1='cyan',col2='orange',normtype=NULL,normaxis=1,rowlabs=TRUE,rowclust=FALSE,rowcols='blue',
                      dis='bray',scale='row',
                      ColSideLabs='Topic Weight',RowSideLabs='Sex',main='Topics over COGs',
                      cexCol=log(h,8),cexRow=log(h,7))
    dev.off()
  }
  
}






# ### For linear regression
# coef_sig <- which(sapply(eff_plot_female$cis,function(x) sum(sign(x))) != 0)
# coef_top <- order(abs(coef_k),decreasing=TRUE)[1:10]
# coef_top[coef_top %in% coef_sig]
# coef1 <- coef_top[coef_top %in% which(coef_k > 0)]
# coef0 <- coef_top[coef_top %in% which(coef_k < 0)]


### For nonlinear splines and the like
coef_top <- order(abs(coef_k_age),decreasing=TRUE)[1:10]
coef1 <- coef_top[coef_top %in% which(coef_k_age > 0)]
coef0 <- coef_top[coef_top %in% which(coef_k_age < 0)]




# top_ord <- order(coef_k,decreasing=TRUE)[c(1:4,(K-4+1):K)]
# vio <- plot_violins2(fit,dat$meta,top_ord)
vio1 <- plot_violins2(fit,dat$meta,coef1,1) + ggtitle('Topics Associated with Female')
vio0 <- plot_violins2(fit,dat$meta,coef0,1) + ggtitle('Topics Associated with Male')
gridExtra::grid.arrange(nrow=2,vio1,vio0)
pdf(file.path(out_dir,K,mod_dir,'vio.pdf'), height=8, width=12); gridExtra::grid.arrange(nrow=2,vio1,vio0); dev.off()




library(mgcv)
th <- fit$theta
colnames(th) <- paste0('T',1:K)
rownames(th) <- dat$meta$SampleID
th <- th[,coef1]

df <- data.frame(th,AGE_RANK=dat$meta$AGE_DECADE_RANK,AGE=as.numeric(dat$meta$AGE_YEARS),SEX=dat$meta$SEX) %>%
  gather(Topic,Abundance,-AGE_RANK,-AGE,-SEX) %>%
  mutate(Topic=factor(Topic,levels=colnames(th)))

vsage1 <- ggplot(df, aes(AGE,Abundance)) + facet_grid(SEX~Topic) + 
  geom_point(data=dplyr::select(df,-SEX),aes(AGE,Abundance),colour='gray',alpha=.5,size=.8) + 
  geom_point(alpha=.8,size=1.2,colour='royalblue4') + #scale_y_log10() +
  #stat_smooth(method='loess',se=FALSE,colour='red') + 
  xlab('AGE') +
  stat_smooth(method='gam',formula = y ~ s(x,k=-1,bs='cr'), colour='red') + 
  theme_bw() + ggtitle('Topics Associated with Female') + ylim(0,.25)

th <- fit$theta
colnames(th) <- paste0('T',1:K)
rownames(th) <- dat$meta$SampleID
th <- th[,coef0]

df <- data.frame(th,AGE_RANK=dat$meta$AGE_DECADE_RANK,AGE=as.numeric(dat$meta$AGE_YEARS),SEX=dat$meta$SEX) %>%
  gather(Topic,Abundance,-AGE_RANK,-AGE,-SEX) %>%
  mutate(Topic=factor(Topic,levels=colnames(th)))


vsage0 <- ggplot(df, aes(AGE,Abundance)) + facet_grid(SEX~Topic) + 
  geom_point(data=dplyr::select(df,-SEX),aes(AGE,Abundance),colour='gray',alpha=.5,size=.8) + 
  geom_point(alpha=.8,size=1.2,colour='royalblue4') + #scale_y_log10() +
  #stat_smooth(method='loess',se=FALSE,colour='red') + 
  xlab('AGE') +
  stat_smooth(method='gam',formula = y ~ s(x,k=-1,bs='cr'), colour='red') +
  theme_bw() + ggtitle('Topics Associated with Male')  + ylim(0,.25)

gridExtra::grid.arrange(nrow=2,vsage1,vsage0)
pdf(file.path(out_dir,K,mod_dir,'vsage.pdf'), height=8, width=12); gridExtra::grid.arrange(nrow=2,vsage1,vsage0); dev.off()






qmat <- apply(mat,2,qnormalize)
pca <- prcomp(qmat,center=TRUE,.scale=TRUE)

df <- data.frame(x=pca$x[,1],y=pca$x[,2],SEX=as.factor(meta$SEX))
tpca1 <- df %>%
  ggplot(aes(x=x,y=y,colour=SEX)) + 
  geom_point(alpha=1,size=1.5) +
  facet_wrap(~SEX,nrow=1) +
  geom_point(data=dplyr::select(df,-SEX),colour='gray',alpha=.4,size=.8) +
  theme_bw() + labs(x='PC1',y='PC2',colour='Sex') + theme(legend.position='none') +
  ggtitle('PCA: Sex')
  
df <- data.frame(x=pca$x[,1],y=pca$x[,2],AGE=as.factor(meta$AGE_DECADE_RANK),SEX=as.factor(meta$SEX)) %>% filter(SEX == 'female')
tpca2_f <- df %>%
  ggplot(aes(x=x,y=y,colour=AGE)) + 
  geom_point(alpha=1,size=1.5) +
  facet_wrap(~AGE,nrow=2) +
  geom_point(data=dplyr::select(df,-AGE),colour='gray',alpha=.4,size=.8) +
  theme_bw() + labs(x='PC1',y='PC2',colour='Age Decade') + theme(legend.position='none') +
  ggtitle('PCA: Age Decade (Female)')

df <- data.frame(x=pca$x[,1],y=pca$x[,2],AGE=as.factor(meta$AGE_DECADE_RANK),SEX=as.factor(meta$SEX)) %>% filter(SEX == 'male')
tpca2_m <- df %>%
  filter(SEX == 'male') %>%
  ggplot(aes(x=x,y=y,colour=AGE)) + 
  geom_point(alpha=1,size=1.5) +
  facet_wrap(~AGE,nrow=2) +
  geom_point(data=dplyr::select(df,-AGE),colour='gray',alpha=.4,size=.8) +
  theme_bw() + labs(x='PC1',y='PC2',colour='Age Decade') + theme(legend.position='none') +
  ggtitle('PCA: Age Decade (Male)')

gridExtra::grid.arrange(nrow=3,tpca1,tpca2_f,tpca2_m)
pdf(file.path(out_dir,K,mod_dir,'tpca.pdf'), height=10, width=10); gridExtra::grid.arrange(nrow=3,tpca1,tpca2_f,tpca2_m); dev.off()
# plotly::plot_ly(data.frame(x=pca$x[,1],y=pca$x[,2],z=pca$x[,3]),
#         x=x,y=y,z=z,
#         type='scatter3d',color=rownames(qmat),
#         mode='markers')



ica <- fastICA::fastICA(qmat,n.comp=2)

df <- data.frame(x=ica$S[,1],y=ica$S[,2],SEX=as.factor(meta$SEX))
tica1 <- df %>%
  ggplot(aes(x=x,y=y,colour=SEX)) + 
  geom_point(alpha=1,size=1.5) +
  facet_wrap(~SEX,nrow=1) +
  geom_point(data=dplyr::select(df,-SEX),colour='gray',alpha=.4,size=.8) +
  theme_bw() + labs(x='IC1',y='IC2',colour='Sex') + theme(legend.position='none') +
  ggtitle('ICA: Sex')

df <- data.frame(x=ica$S[,1],y=ica$S[,2],AGE=as.factor(meta$AGE_DECADE_RANK),SEX=as.factor(meta$SEX)) %>% filter(SEX == 'female')
tica2_f <- df %>%
  ggplot(aes(x=x,y=y,colour=AGE)) + 
  geom_point(alpha=1,size=1.5) +
  facet_wrap(~AGE,nrow=2) +
  geom_point(data=dplyr::select(df,-AGE),colour='gray',alpha=.4,size=.8) +
  theme_bw() + labs(x='IC1',y='IC2',colour='Age Decade') + theme(legend.position='none') +
  ggtitle('ICA: Age Decade (Female)')

df <- data.frame(x=ica$S[,1],y=ica$S[,2],AGE=as.factor(meta$AGE_DECADE_RANK),SEX=as.factor(meta$SEX)) %>% filter(SEX == 'male')
tica2_m <- df %>%
  ggplot(aes(x=x,y=y,colour=AGE)) + 
  geom_point(alpha=1,size=1.5) +
  facet_wrap(~AGE,nrow=2) +
  geom_point(data=dplyr::select(df,-AGE),colour='gray',alpha=.4,size=.8) +
  theme_bw() + labs(x='IC1',y='IC2',colour='Age Decade') + theme(legend.position='none') +
  ggtitle('ICA: Age Decade (Male)')


gridExtra::grid.arrange(nrow=3,tica1,tica2_f,tica2_m)
pdf(file.path(out_dir,K,mod_dir,'tica.pdf'), height=10, width=10); gridExtra::grid.arrange(nrow=3,tica1,tica2_f,tica2_m); dev.off()
# plotly::plot_ly(data.frame(x=ica$S[,1],y=ica$S[,2],z=ica$S[,3]),
#         x=x,y=y,z=z,
#         type='scatter3d',color=rownames(qmat),
#         mode='markers')









pca <- prcomp(qmat[,c(coef1,coef0)],center=TRUE,.scale=TRUE)

df <- data.frame(x=pca$x[,1],y=pca$x[,2],SEX=as.factor(meta$SEX))
tpca1 <- df %>%
  ggplot(aes(x=x,y=y,colour=SEX)) + 
  geom_point(alpha=1,size=1.5) +
  facet_wrap(~SEX,nrow=1) +
  geom_point(data=dplyr::select(df,-SEX),colour='gray',alpha=.4,size=.8) +
  theme_bw() + labs(x='PC1',y='PC2',colour='Sex') + theme(legend.position='none') +
  ggtitle('PCA: Sex')

df <- data.frame(x=pca$x[,1],y=pca$x[,2],AGE=as.factor(meta$AGE_DECADE_RANK))
tpca2 <- df %>%
  ggplot(aes(x=x,y=y,colour=AGE)) + 
  geom_point(alpha=1,size=1.5) +
  facet_wrap(~AGE,nrow=2) +
  geom_point(data=dplyr::select(df,-AGE),colour='gray',alpha=.4,size=.8) +
  theme_bw() + labs(x='PC1',y='PC2',colour='Age Decade') + theme(legend.position='none') +
  ggtitle('PCA: Age Decade')

gridExtra::grid.arrange(nrow=2,tpca1,tpca2)
pdf(file.path(out_dir,K,mod_dir,'toptpca.pdf'), height=8, width=10); gridExtra::grid.arrange(nrow=2,tpca1,tpca2); dev.off()
# plotly::plot_ly(data.frame(x=pca$x[,1],y=pca$x[,2],z=pca$x[,3]),
#         x=x,y=y,z=z,
#         type='scatter3d',color=rownames(qmat),
#         mode='markers')








ica <- fastICA::fastICA(qmat[,c(coef1,coef0)],n.comp=2)

df <- data.frame(x=ica$S[,1],y=ica$S[,2],SEX=as.factor(meta$SEX))
tica1 <- df %>%
  ggplot(aes(x=x,y=y,colour=SEX)) + 
  geom_point(alpha=1,size=1.5) +
  facet_wrap(~SEX,nrow=1) +
  geom_point(data=dplyr::select(df,-SEX),colour='gray',alpha=.4,size=.8) +
  theme_bw() + labs(x='IC1',y='IC2',colour='Sex') + theme(legend.position='none') +
  ggtitle('ICA: Sex')

df <- data.frame(x=ica$S[,1],y=ica$S[,2],AGE=as.factor(meta$AGE_DECADE_RANK))
tica2 <- df %>%
  ggplot(aes(x=x,y=y,colour=AGE)) + 
  geom_point(alpha=1,size=1.5) +
  facet_wrap(~AGE,nrow=2) +
  geom_point(data=dplyr::select(df,-AGE),colour='gray',alpha=.4,size=.8) +
  theme_bw() + labs(x='IC1',y='IC2',colour='Age Decade') + theme(legend.position='none') +
  ggtitle('ICA: Age Decade')

gridExtra::grid.arrange(nrow=2,tica1,tica2)
pdf(file.path(out_dir,K,mod_dir,'toptica.pdf'), height=8, width=10); gridExtra::grid.arrange(nrow=2,tica1,tica2); dev.off()
# plotly::plot_ly(data.frame(x=ica$S[,1],y=ica$S[,2],z=ica$S[,3]),
#         x=x,y=y,z=z,
#         type='scatter3d',color=rownames(qmat),
#         mode='markers')

































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
z <- z[dat$counts$meta_clean %>% filter(DIAGNOSIS == 'CD') %>% dplyr::select(SampleID) %>% unlist(),]
z <- z[,colSums(z) > 0]

b <- fit$beta$logbeta[[1]][17,]
names(b) <- paste0('otu',1:length(b))
o <- names(which(b > quantile(b,.95)))
length(o)

z <- z[,o]

df <- data.frame(z,dat$counts$meta_clean %>% filter(DIAGNOSIS == 'CD') %>% dplyr::select(SampleID,PCDAI,PCDAI_RANK)) %>%
  gather(OTU,Abundance,-SampleID,-PCDAI,-PCDAI_RANK) %>%
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
  map(~ coef(summary(glm(Abundance ~ PCDAI, data=., family=quasipoisson(link='log'))))[2,]) %>%
  do.call('rbind',.)
pvals <- as.data.frame(pvals)
pvals$adj <- p.adjust(pvals[,4],method='BH')
pvals_sig <- rownames(pvals)[pvals$adj < .05]


pdf('~/MiscOut/testing_formanu.pdf', height=12, width=15)
df %>%
  mutate(Significant = as.factor(ifelse(OTU %in% pvals_sig,1,0))) %>%
  ggplot(aes(PCDAI,Abundance,colour=Significant)) + 
  geom_point(data=dplyr::select(df,-Facet),aes(PCDAI,Abundance),colour='gray',alpha=.5,size=1) +
  facet_wrap(~Facet) + 
  geom_point(alpha=1,size=2) +
  scale_y_log10() +
  scale_colour_manual(values=c('black','red')) + 
  stat_smooth(method='lm',se=FALSE,colour='darkgreen') +
  theme_bw() + theme(legend.position='none') + xlab('PCDAI')
dev.off()





z <- dat$counts$table_clean
z <- z[dat$counts$meta_clean %>% dplyr::select(SampleID) %>% unlist(),]
z <- z[,colSums(z) > 0]

b <- fit$beta$logbeta[[1]][17,]
names(b) <- paste0('otu',1:length(b))

o <- labelTopics(fit,c(17),n=5)$interaction

z <- z[,c(o[1,],o[2,])]

df <- data.frame(z,dat$counts$meta_clean %>% dplyr::select(SampleID,DIAGNOSIS)) %>%
  gather(OTU,Abundance,-SampleID,-DIAGNOSIS) %>%
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
ggplot(df, aes(PCDAI_RANK,Abundance,colour=Facet)) + 
  geom_point(data=dplyr::select(df,-Facet),aes(PCDAI_RANK,Abundance),colour='gray',alpha=.5,size=1) +
  facet_wrap(~Facet) + 
  geom_point(alpha=1,size=2,colour='black') +
  scale_y_log10() +
  stat_smooth(method='lm',se=FALSE,colour='red') +
  theme_bw() + xlab('PCDAI Rank')
#dev.off()










PCDAI <- pcdai_id$PCDAI
qplot(1:length(pcdai_id$PCDAI),dat$counts$table_clean[pcdai_id$SampleID,'otu22'] + 1, colour=PCDAI) + scale_y_log10() +
  xlab('PCDAI Rank') + ylab('Abundance') + ggtitle('T5: Enterobacteriaceae OTU')

dat$taxa[as.character(dat$counts$ids[c('otu356','otu270','otu21'),]$long),]
dat$taxa[as.character(dat$counts$ids[c('otu190','otu201'),]$long),]


data.frame(dat$counts$table_clean[,c('otu356','otu270','otu21','otu190','otu201')],Diagnosis=dat$counts$meta_clean$DIAGNOSIS,PCDAI=dat$counts$meta_clean$PCDAI) %>%
  gather(OTU,Abundance,-PCDAI,-Diagnosis) %>%
  ggplot(aes(Abundance,)) + geom_point() + facet_wrap(~Diagnosis) + scale_y_log10()



data.frame(th,PCDAI=pcdai_id$PCDAI,Rank=1:NROW(th)) %>%
  gather(Topic,Abundance,-Rank,-PCDAI) %>%
  mutate(Topic=factor(Topic,levels=colnames(th))) %>%
  filter(Topic == 'T17') %>%
  ggplot(aes(Rank,Abundance,colour=PCDAI)) + geom_point(alpha=.5) + scale_y_log10() +
  stat_smooth(method='lm',se=FALSE,colour='red') + xlab('PCDAI Rank')


th <- fit$theta[,17]
bugs <- dat$counts$table_clean[meta$SampleID,paste0('otu',which(b1 > quantile(b1,.8)))]
data.frame(T17=th,RANK=dat$meta$PCDAI_RANK,bugs=rowSums(bugs)) %>%
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

