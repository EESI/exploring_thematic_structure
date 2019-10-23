library(phyloseq)
library(DESeq2)

coef_k <- out$dx$lasso

beta_table_name <- 'beta_table'
biom_file <- read_biom(file.path(save_dir,save_fit_foldername,paste0(beta_table_name,'_predicted_metagenome.biom')))
kegg_metadata_risk <- kegg_pw_collapse(biom_file,3)
kegg_risk <- as.matrix(biom_data(biom_file))
kegg_risk <- kegg_risk[rowSums(kegg_risk) >= 5,]
kegg_gene_names <- kegg_pw_collapse(biom_file,0)

DF <- data.frame(lasso=coef_k) %>%
   mutate(assoc=ifelse(lasso>0,'CD+',ifelse(lasso<0,'CD-','ZERO')),
          lasso_scaled=coef_k/sd(coef_k))
rownames(DF) <- colnames(kegg_risk)

OTU <- otu_table(kegg_risk,taxa_are_rows=TRUE)
SAMP <- sample_data(DF)
PS <- phyloseq(OTU,SAMP)

DS2 <- phyloseq_to_deseq2(PS, ~ lasso)
DS2 <- DESeq2::DESeq(DS2, test="Wald", fitType="parametric")
res <- results(DS2, cooksCutoff=FALSE, pAdjustMethod='BH')
alpha <- .1
sigtab <- res[which(res$padj < alpha), ]

kegg_sig <- log10(kegg_risk[rownames(sigtab),] + 1)
rownames(kegg_sig) <- kegg_gene_names[rownames(kegg_sig)]

pdf('~/MiscOut/tests1.pdf', height=24, width=20)
plot_big_heatmap(kegg_sig,coef_k,zbar=FALSE,col1='cyan',col2='black',col3='orange',cexCol=1.5,cexRow=.75,main='KEGG',scale='row')
dev.off()





beta_table_name <- 'beta_table'
biom_file <- read_biom(file.path(save_dir,save_fit_foldername,paste0(beta_table_name,'_predicted_metagenome.biom')))
kegg_metadata_risk <- kegg_pw_collapse(biom_file,3)
kegg_risk <- as.matrix(biom_data(biom_file))
kegg_risk <- kegg_risk[rowSums(kegg_risk) >= 5,]
kegg_gene_names <- kegg_pw_collapse(biom_file,0)

DF <- data.frame(lasso=coef_k) %>%
   mutate(assoc=ifelse(lasso>0,'CD+',ifelse(lasso<0,'CD-','ZERO')),
          lasso_scaled=coef_k/sd(coef_k))
rownames(DF) <- colnames(kegg_risk)

OTU <- otu_table(kegg_risk,taxa_are_rows=TRUE)
SAMP <- sample_data(DF)
PS <- phyloseq(OTU,SAMP)

DS2 <- phyloseq_to_deseq2(PS, ~ lasso_scaled)
DS2 <- DESeq2::DESeq(DS2, test="Wald", fitType="parametric")
res <- results(DS2, cooksCutoff=FALSE, pAdjustMethod='BH')
alpha <- .1
sigtab <- res[which(res$padj < alpha), ]

kegg_sig <- log10(kegg_risk[rownames(sigtab),] + 1)
rownames(kegg_sig) <- kegg_gene_names[rownames(kegg_sig)]

pdf('~/MiscOut/tests1.pdf', height=24, width=20)
plot_big_heatmap(kegg_sig,coef_k,zbar=FALSE,col1='cyan',col2='black',col3='orange',cexCol=1.5,cexRow=.75,main='KEGG',scale='row')
dev.off()






beta_table_name <- 'beta_table'
biom_file <- read_biom(file.path(save_dir,save_fit_foldername,paste0(beta_table_name,'_predicted_metagenome.biom')))
kegg_metadata_risk <- kegg_pw_collapse(biom_file,3)
kegg_risk <- as.matrix(biom_data(biom_file))
kegg_risk <- kegg_risk[rowSums(kegg_risk) >= 5,]
kegg_gene_names <- kegg_pw_collapse(biom_file,0)

DF <- data.frame(lasso=coef_k) %>%
   mutate(assoc=ifelse(lasso>0,'CD+',ifelse(lasso<0,'CD-','ZERO')),
          lasso_scaled=coef_k/sd(coef_k))
rownames(DF) <- colnames(kegg_risk)

OTU <- otu_table(kegg_risk,taxa_are_rows=TRUE)
SAMP <- sample_data(DF)
PS <- phyloseq(OTU,SAMP)

DS2 <- phyloseq_to_deseq2(PS, ~ lasso)
DS2 <- DESeq2::DESeq(DS2, test="Wald", fitType="parametric")
res <- results(DS2, cooksCutoff=FALSE, pAdjustMethod='BH')
alpha <- .01
sigtab <- res[which(res$padj < alpha), ]

kegg_sig <- log10(kegg_risk[rownames(sigtab),] + 1)
rownames(kegg_sig) <- kegg_gene_names[rownames(kegg_sig)]

pdf('~/MiscOut/tests1.pdf', height=24, width=20)
plot_big_heatmap(kegg_sig,coef_k,zbar=FALSE,col1='cyan',col2='black',col3='orange',cexCol=1.5,cexRow=.75,main='KEGG',scale='row')
dev.off()







beta_table_name <- 'beta_table'
biom_file <- read_biom(file.path(save_dir,save_fit_foldername,paste0(beta_table_name,'_predicted_metagenome.biom')))
kegg_metadata_risk <- kegg_pw_collapse(biom_file,3)
kegg_risk <- as.matrix(biom_data(biom_file))
kegg_risk <- kegg_risk[rowSums(kegg_risk) >= 5,]
kegg_gene_names <- kegg_pw_collapse(biom_file,0)

DF <- data.frame(lasso=coef_k) %>%
   mutate(assoc=ifelse(lasso>0,'CD+',ifelse(lasso<0,'CD-','ZERO')),
          lasso_scaled=coef_k/sd(coef_k),
          lasso_ranked=rank(coef_k))
rownames(DF) <- colnames(kegg_risk)

OTU <- otu_table(kegg_risk,taxa_are_rows=TRUE)
SAMP <- sample_data(DF)
PS <- phyloseq(OTU,SAMP)

DS2 <- phyloseq_to_deseq2(PS, ~ lasso_ranked)
DS2 <- DESeq2::DESeq(DS2, test="Wald", fitType="parametric")
res <- results(DS2, cooksCutoff=FALSE, pAdjustMethod='BH')
alpha <- .01
sigtab <- res[which(res$padj < alpha), ]

kegg_sig <- kegg_risk[rownames(sigtab),]
rownames(kegg_sig) <- kegg_gene_names[rownames(kegg_sig)]
kegg_sig <- kegg_sig[!(rownames(kegg_sig) %in% c('None','hypothetical protein')),]

pdf('~/MiscOut/tests1.pdf', height=40, width=20)
plot_big_heatmap(kegg_sig,coef_k,zbar=FALSE,col1='cyan',col2='black',col3='orange',cexCol=1.5,cexRow=.75,main='KEGG',scale='row')
dev.off()
