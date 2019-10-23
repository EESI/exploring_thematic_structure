plothm <- function(mat,coef_k,col1='cyan',col2='black',col3='orange',...){
   
   require(vegan)
   
   top_ord <- order(coef_k,decreasing=TRUE)
   mat <- mat[,top_ord]
   colnames(mat) <- top_ord
   
   K <- NCOL(mat)
   
   if (is.null(col3)) colors <- colorRampPalette(c(col1, col2))(n=100) else colors <- colorRampPalette(c(col1, col2, col3))(n=199)
   
   coefcols <- ramp_colors(coef_k)
   coefcols <- coefcols[top_ord]
   
   col_breaks <- which(sort(coef_k,decreasing=TRUE) == 0)
   col_breaks <- c(min(col_breaks)-1,max(col_breaks))
   
   labcols <- rainbow(length(unique(rownames(mat))))
   labcols <- labcols[as.integer(as.factor(rownames(mat)))]
   
   #    row_breaks <- rownames(mat)
   #    row_breaks <- na.omit(sapply(2:length(row_breaks), function(i) ifelse(row_breaks[i] != row_breaks[i-1],i,NA)))-1
   
   
   dismat1 <- vegdist(mat,method='bray')
   rc <- hclust(dismat1,'ward.D2')
   
   
   gplots::heatmap.2(as.matrix(mat),
                     density.info="histogram",
                     denscol='black',
                     labCol=top_ord,
                     labRow=rownames(mat),
                     #dendrogram = 'none',
                     dendrogram = 'row',
                     Colv=FALSE,
                     #Rowv=FALSE,
                     Rowv=as.dendrogram(rc),
                     col=colors,
                     ColSideColors=coefcols,
                     RowSideColors=labcols,
                     key=FALSE,
                     #keysize=2,
                     xlab='',ylab='',
                     trace="none",
                     tracecol='yellow',
                     #rowsep=row_breaks,
                     colsep=col_breaks,
                     margins=c(25,45),
                     ...
   )
}

plot_big_heatmap <- function(mat,coef_k,zbar,col1='black',col2='green',col3=NULL,rowlabs=TRUE,rowclust=TRUE,marg=c(25,25),...){
   
   require(vegan)
   
   top_ord <- order(coef_k,decreasing=TRUE)
   mat <- mat[,top_ord]
   colnames(mat) <- top_ord
   
   K <- NCOL(mat)
   
   if (is.null(col3)) colors <- colorRampPalette(c(col1, col2))(n=100) else colors <- colorRampPalette(c(col1, col2, col3))(n=199)
   
   coefcols <- ramp_colors(coef_k)
   coefcols <- coefcols[top_ord]
   
   col_breaks <- which(sort(coef_k,decreasing=TRUE) == 0)
   col_breaks <- c(min(col_breaks)-1,max(col_breaks))
   
   if (zbar){
      labcols <- c('red','blue')
      labcols <- labcols[as.integer(as.factor(rownames(mat)))]
      
      row_breaks <- rownames(mat)
      row_breaks <- na.omit(sapply(2:length(row_breaks), function(i) ifelse(row_breaks[i] != row_breaks[i-1],i,NA)))-1
      
      gplots::heatmap.2(as.matrix(mat),
                        density.info="histogram",
                        denscol='black',
                        labCol=top_ord,
                        labRow='',
                        dendrogram = 'none',
                        Colv=FALSE,
                        Rowv=FALSE,
                        col=colors,
                        ColSideColors=coefcols,
                        RowSideColors=labcols,
                        key=TRUE,
                        #keysize=2,
                        xlab='',ylab='',
                        trace="none",
                        tracecol='yellow',
                        rowsep=row_breaks,
                        colsep=col_breaks,
                        margins=marg,
                        ...
      )
   }else{
      
      dismat1 <- vegdist(mat,method='bray')
      rc <- hclust(dismat1,'ward.D2')
      
      if (rowclust) {
         if (rowlabs) rowlabs <- rc$labels else rowlabs <- ''
         roword <- as.dendrogram(rc)
         dendro <- 'row'
      }else{
         if (rowlabs) rowlabs <- rownames(mat) else rowlabs <- ''
         roword <- FALSE
         dendro <- 'none'
      }
      
      gplots::heatmap.2(as.matrix(mat),
                        density.info="histogram",
                        denscol='black',
                        labCol=top_ord,
                        labRow=rowlabs,
                        dendrogram = dendro,
                        Colv=FALSE,
                        Rowv=roword,
                        col=colors,
                        ColSideColors=coefcols,
                        #RowSideColors=labcols,
                        key=TRUE,
                        #keysize=2,
                        xlab='',ylab='',
                        trace="none",
                        tracecol='yellow',
                        #rowsep=row_breaks,
                        colsep=col_breaks,
                        margins=marg,
                        ...
      )
      
      return(rc)
   }
}

plot_big_heatmap2 <- function(mat,z,coef_k=NULL,top_ord=NULL,col1='black',col2='green',col3=NULL,...){
   
   require(vegan)
   
   if (z){
      if (!is.null(coef_k)){
         top_ord <- order(coef_k,decreasing=TRUE)
         mat <- mat[,top_ord]
         colnames(mat) <- top_ord
         
         coefcols <- ramp_colors(coef_k)
         coefcols <- coefcols[top_ord]
         
         col_breaks <- which(sort(coef_k,decreasing=TRUE) == 0)
         col_breaks <- c(min(col_breaks)-1,max(col_breaks))
         
         dendro <- 'none'
         cv <- FALSE
      }else{
         cdismat <- vegdist(t(mat),method='bray')
         cc <- hclust(cdismat,'ward.D2')
         top_ord <- cc$order
         
         dendro <- 'column'
         cv <- as.dendrogram(cc)
      }
   }else{
      if (!is.null(coef_k)){
         top_ord <- order(coef_k,decreasing=TRUE)
         mat <- mat[,top_ord]
         colnames(mat) <- top_ord
         
         coefcols <- ramp_colors(coef_k)
         coefcols <- coefcols[top_ord]
         
         col_breaks <- which(sort(coef_k,decreasing=TRUE) == 0)
         col_breaks <- c(min(col_breaks)-1,max(col_breaks))
      }else if(!is.null(top_ord)){
         mat <- mat[,top_ord]
         colnames(mat) <- top_ord
      }else{
         stop('Must supply either topic coeffiicents or order indexes if matrix is not zbar.')
      }
   }
   
   
   K <- NCOL(mat)
   
   if (is.null(col3)) colors <- colorRampPalette(c(col1, col2))(n=100) else colors <- colorRampPalette(c(col1, col2, col3))(n=199)
   
   
   if (z){
      labcols <- c('red','blue')
      labcols <- labcols[as.integer(as.factor(rownames(mat)))]
      
      row_breaks <- rownames(mat)
      row_breaks <- na.omit(sapply(2:length(row_breaks), function(i) ifelse(row_breaks[i] != row_breaks[i-1],i,NA)))-1
      
      if (!is.null(coef_k)){
         gplots::heatmap.2(as.matrix(mat),
                           density.info="histogram",
                           denscol='black',
                           labCol=top_ord,
                           labRow='',
                           dendrogram = dendro,
                           Colv=cv,
                           Rowv=FALSE,
                           col=colors,
                           ColSideColors=coefcols,
                           RowSideColors=labcols,
                           key=TRUE,
                           #keysize=2,
                           xlab='',ylab='',
                           trace="none",
                           tracecol='yellow',
                           rowsep=row_breaks,
                           colsep=col_breaks,
                           margins=c(25,25),
                           ...
         )
      }else{
         gplots::heatmap.2(as.matrix(mat),
                           density.info="histogram",
                           denscol='black',
                           labCol=1:K,
                           labRow='',
                           dendrogram = dendro,
                           Colv=cv,
                           Rowv=FALSE,
                           col=colors,
                           #ColSideColors=coefcols,
                           RowSideColors=labcols,
                           key=TRUE,
                           #keysize=2,
                           xlab='',ylab='',
                           trace="none",
                           tracecol='yellow',
                           rowsep=row_breaks,
                           #colsep=col_breaks,
                           margins=c(25,25),
                           ...
         )
         
         return(top_ord)
      }
   }else{
      
      rdismat <- vegdist(mat,method='bray')
      rc <- hclust(rdismat,'ward.D2')
      
      if (!is.null(coef_k)){
         gplots::heatmap.2(as.matrix(mat),
                           density.info="histogram",
                           denscol='black',
                           labCol=top_ord,
                           labRow=rc$labels,
                           #labRow=rownames(mat),
                           #dendrogram = 'none',
                           dendrogram = 'row',
                           Colv=FALSE,
                           #Rowv=FALSE,
                           Rowv=as.dendrogram(rc),
                           col=colors,
                           ColSideColors=coefcols,
                           #RowSideColors=labcols,
                           key=TRUE,
                           #keysize=2,
                           xlab='',ylab='',
                           trace="none",
                           tracecol='yellow',
                           #rowsep=row_breaks,
                           colsep=col_breaks,
                           margins=c(25,25),
                           ...
         )
      }else{
         gplots::heatmap.2(as.matrix(mat),
                           density.info="histogram",
                           denscol='black',
                           labCol=top_ord,
                           labRow=rc$labels,
                           #labRow=rownames(mat),
                           #dendrogram = 'none',
                           dendrogram = 'row',
                           Colv=FALSE,
                           #Rowv=FALSE,
                           Rowv=as.dendrogram(rc),
                           col=colors,
                           #RowSideColors=labcols,
                           key=TRUE,
                           #keysize=2,
                           xlab='',ylab='',
                           trace="none",
                           tracecol='yellow',
                           #rowsep=row_breaks,
                           margins=c(25,25),
                           ...
         )
      }
   }
}


out_dir <- '~/Dropbox/stm_microbiome/output_active/Gevers_ti'
heatmap_dir <- file.path(out_dir,save_fit_foldername,'heatmaps')
dir.create(heatmap_dir,showWarnings=FALSE,recursive=TRUE)

pdf(file.path(heatmap_dir,'heatmap_theta.pdf'), height=10, width=20)
plot_big_heatmap(mat,coef_k,zbar=TRUE,col1='cyan',col2='black',col3='orange',cexCol=1.5,main='Training',scale='column')
plot_big_heatmap(mat2,coef_k,zbar=TRUE,col1='cyan',col2='black',col3='orange',cexCol=1.5,main='Testing',scale='column')
dev.off()
pdf(file.path(heatmap_dir,'heatmap_beta.pdf'), height=20, width=20)
plot_big_heatmap(taxon_table,coef_k,zbar=FALSE,col1='cyan',col2='black',col3='orange',cexCol=1.5,cexRow=.75,main='Taxa',scale='row')
dev.off()
pdf(file.path(heatmap_dir,'heatmap_kegg.pdf'), height=20, width=20)
plot_big_heatmap(kegg_sig,coef_k,zbar=FALSE,col1='cyan',col2='black',col3='orange',rowlabs=FALSE,cexCol=1.5,cexRow=1,main='KEGG',scale='row')
dev.off()
pdf(file.path(heatmap_dir,'heatmap_kegg_rownames.pdf'), height=70, width=30)
hm_kegg <- plot_big_heatmap(kegg_sig,coef_k,zbar=FALSE,col1='cyan',col2='black',col3='orange',rowlabs=TRUE,cexCol=1.5,cexRow=1,main='KEGG',scale='row')
dev.off()
pdf(file.path(heatmap_dir,'heatmap_cog.pdf'), height=20, width=20)
plot_big_heatmap(cog_sig,coef_k,zbar=FALSE,col1='cyan',col2='black',col3='orange',rowlabs=FALSE,cexCol=1.5,cexRow=1,main='COG',scale='row')
dev.off()
pdf(file.path(heatmap_dir,'heatmap_cog_rownames.pdf'), height=70, width=30)
hm_cog <- plot_big_heatmap(cog_sig,coef_k,zbar=FALSE,col1='cyan',col2='black',col3='orange',rowlabs=TRUE,cexCol=1.5,cexRow=1,main='COG',scale='row')
dev.off()


tree_groups <- cutree(hm_kegg,h=2)
table(tree_groups) # number of genes in a group
tree_groups <- tree_groups[rev(hm_kegg$order)] # reorder the group indexes to match the clustering indexes
# group_idx <- 4
# idx_groups <- rev(hm_kegg$order)[tree_groups == group_idx] # capture 4th group
# rownames(kegg_sig)[idx_groups] # checker

sink(file.path(heatmap_dir,'group_kegg_order.txt'))
cat(unique(tree_groups))
sink()
for (g in unique(tree_groups)){
   group_file_name <- paste0('heatmap_kegg_group_',g)
   idx_groups <- rev(hm_kegg$order)[tree_groups == g] # capture 4th group
   
   pws <- kegg_metadata_risk[names(ko_names)[idx_groups]]
   sink(file.path(heatmap_dir,paste0(group_file_name,'.txt')))
   cat('Group',g,'\n\n\n')
   cat('Gene Name [EC Identifier] (KO term): Pathways')
   for (i in seq_along(pws)){
      cat(ko_names[names(pws[i])],' (',names(pws[i]),'):  ',paste0(pws[[i]],collapse=', '),'\n',sep='')
   }
   sink()
   
   group_mat <- kegg_sig[idx_groups,]
   h <- ifelse(floor(nrow(group_mat)/2) < 8, 8, floor(nrow(group_mat)/2))
   
   pdf(file.path(heatmap_dir,paste0(group_file_name,'.pdf')), height=h, width=15)
   plot_big_heatmap(group_mat,coef_k,zbar=FALSE,col1='cyan',col2='black',col3='orange',rowlabs=TRUE,rowclust=FALSE,
                    cexCol=1.5,cexRow=1,main='KEGG',scale='row')
   dev.off()
}


tree_groups <- cutree(hm_cog,h=2)
table(tree_groups) # number of genes in a group
tree_groups <- tree_groups[rev(hm_cog$order)] # reorder the group indexes to match the clustering indexes
# group_idx <- 4
# idx_groups <- rev(hm_cog$order)[tree_groups == group_idx] # capture 4th group
# rownames(kegg_sig)[idx_groups] # checker

sink(file.path(heatmap_dir,'group_cog_order.txt'))
cat(unique(tree_groups))
sink()
for (g in unique(tree_groups)){
   group_file_name <- paste0('heatmap_cog_group_',g)
   idx_groups <- rev(hm_cog$order)[tree_groups == g] # capture 4th group
   
   pws <- cog_metadata_risk[names(cog_names)[idx_groups]]
   sink(file.path(heatmap_dir,paste0(group_file_name,'.txt')))
   cat('Group',g,'\n\n\n')
   cat('Gene Name (COG term): Pathways')
   for (i in seq_along(pws)){
      cat(cog_names[names(pws[i])],' (',names(pws[i]),'):  ',paste0(pws[[i]],collapse=', '),'\n',sep='')
   }
   sink()
   
   group_mat <- cog_sig[idx_groups,]
   h <- ifelse(floor(nrow(group_mat)/2) < 8, 8, floor(nrow(group_mat)/2))
   
   pdf(file.path(heatmap_dir,paste0(group_file_name,'.pdf')), height=h, width=15)
   plot_big_heatmap(group_mat,coef_k,zbar=FALSE,col1='cyan',col2='black',col3='orange',rowlabs=TRUE,rowclust=FALSE,
                    cexCol=1.5,cexRow=1,main='KEGG',scale='row')
   dev.off()
}
















# Get coefficients from lasso
coef_k <- out$dx$lasso


# prep data for training set heatmap
####################################
####################################
mat <- train_fit$Z_bar
meta <- train_meta

filt_idx <- rowSums(mat) != 0
mat <- mat[filt_idx,]
meta <- meta[filt_idx,]
lab <- meta[,'DIAGNOSIS']
pcdai <- meta[,'PCDAI']
pcdai[is.na(meta$PCDAI) & meta$DIAGNOSIS == 'Not IBD'] <- 0
filt_idx <- !is.na(pcdai) 
mat <- mat[filt_idx,]
meta <- meta[filt_idx,]
lab <- meta[,variable]
pcdai <- meta$PCDAI
lab_ord <- order(pcdai,decreasing=TRUE)
lab <- lab[lab_ord]
mat <- mat[lab_ord,]
rownames(mat) <- lab
colnames(mat) <- 1:NCOL(mat)
####################################
####################################


# prep data for testing set heatmap
####################################
####################################
mat2 <- test_fit$Z_bar
meta2 <- test_meta

filt_idx <- rowSums(mat2) != 0
mat2 <- mat2[filt_idx,]
meta2 <- meta2[filt_idx,]
lab <- meta2[,'DIAGNOSIS']
pcdai <- meta2[,'PCDAI']
pcdai[is.na(meta2$PCDAI) & meta2$DIAGNOSIS == 'Not IBD'] <- 0
filt_idx <- !is.na(pcdai) 
mat2 <- mat2[filt_idx,]
meta2 <- meta2[filt_idx,]
lab <- meta2[,variable]
pcdai <- meta2$PCDAI
lab_ord <- order(pcdai,decreasing=TRUE)
lab <- lab[lab_ord]
mat2 <- mat2[lab_ord,]
rownames(mat2) <- lab
colnames(mat2) <- 1:NCOL(mat2)
####################################
####################################


# prep data for taxon table heatmap
####################################
####################################
taxon <- 7 # species
taxon_names <- apply(beta_otu[rownames(beta_frozen),1:taxon],1,paste0,collapse='|')
taxon_table <- data.frame(beta_frozen,taxon_names) %>%
   group_by(taxon_names) %>%
   summarise_each(funs(sum)) 
taxon_names <- pretty_taxa_names(taxon_table$taxon_names)
taxon_table <- data.frame(dplyr::select(taxon_table,-taxon_names))
beta_frozen <- t(t(beta_frozen)/colSums(beta_frozen))
rownames(taxon_table) <- taxon_names
colnames(taxon_table) <- 1:NCOL(taxon_table)

filt_idx <- rowSums(taxon_table) != 0
taxon_table <- taxon_table[filt_idx,]
####################################
####################################


# prep data for kegg heatmap
####################################
####################################
beta_table_name <- 'beta_table'
biom_file <- read_biom(file.path(save_dir,save_fit_foldername,paste0(beta_table_name,'_predicted_metagenome.biom')))
kegg_metadata_risk <- kegg_pw_collapse(biom_file,3)
kegg_risk <- as.matrix(biom_data(biom_file))

options(scipen=10000)
Y <- rowSums(kegg_risk > 0)+1
X <- rowSums(kegg_risk)+1
qplot(X,Y,geom='point',alpha=.1) + xlab('Log10 Abundance') + ylab('Prevalence') +
   scale_x_log10() + theme(legend.position='none') + geom_vline(xintercept=13,color='red')

abund_by_filt_x <- seq(1,1000,1)
abund_by_filt <- sapply(abund_by_filt_x,function(i) sum(rowSums(kegg_risk) >= i))
abund_thresh <- 13
data.frame(x=abund_by_filt_x,y=abund_by_filt) %>%
   filter(x < 500) %>%
   ggplot(aes(x,y)) + geom_point() + geom_vline(xintercept=abund_thresh,color='red')


kegg_risk <- kegg_risk[rowSums(kegg_risk) >= abund_thresh,]
kegg_gene_names <- kegg_pw_collapse(biom_file,0)
kegg_gene_names <- lapply(kegg_gene_names,function(x) paste0(x,collapse=' | '))

DF <- data.frame(lasso=coef_k) %>%
   mutate(assoc=ifelse(lasso>0,'CD+',ifelse(lasso<0,'CD-','ZERO')),
          lasso_scaled=coef_k/sd(coef_k),
          lasso_ranked=rank(coef_k))
rownames(DF) <- colnames(kegg_risk)

KEGG <- otu_table(kegg_risk,taxa_are_rows=TRUE)
SAMP <- sample_data(DF)
PS <- phyloseq(KEGG,SAMP)

DS2 <- phyloseq_to_deseq2(PS, ~ lasso_ranked)
DS2 <- DESeq2::DESeq(DS2, test="Wald", fitType="parametric")
res <- DESeq2::results(DS2, cooksCutoff=FALSE, pAdjustMethod='BH')
alpha <- .01
sigtab <- res[which(res$padj < alpha), ]

kegg_sig <- kegg_risk[rownames(sigtab),]
rownames(kegg_sig) <- kegg_gene_names[rownames(kegg_sig)]

ko_names <- rownames(kegg_sig) 
names(ko_names) <- rownames(sigtab)

misc_filter <- !(rownames(kegg_sig) %in% c('None','hypothetical protein'))
kegg_sig <- kegg_sig[misc_filter,]
ko_names <- ko_names[misc_filter]
####################################
####################################




# prep data for kegg heatmap
####################################
####################################
beta_table_name <- 'beta_table'
biom_file <- read_biom(file.path(save_dir,save_fit_foldername,paste0(beta_table_name,'_predicted_cogs.biom')))
cog_metadata_risk <- kegg_pw_collapse(biom_file,2,cog=TRUE)
cog_risk <- as.matrix(biom_data(biom_file))

options(scipen=10000)
Y <- rowSums(cog_risk > 0)+1
X <- rowSums(cog_risk)+1
qplot(X,Y,geom='point',alpha=.1) + xlab('Log10 Abundance') + ylab('Prevalence') +
   scale_x_log10() + theme(legend.position='none') + geom_vline(xintercept=13,color='red')

abund_by_filt_x <- seq(1,1000,1)
abund_by_filt <- sapply(abund_by_filt_x,function(i) sum(rowSums(cog_risk) >= i))
abund_thresh <- 9
data.frame(x=abund_by_filt_x,y=abund_by_filt) %>%
   filter(x < 500) %>%
   ggplot(aes(x,y)) + geom_point() + geom_vline(xintercept=abund_thresh,color='red')

cog_risk <- cog_risk[rowSums(cog_risk) >= abund_thresh,]
cog_gene_names <- kegg_pw_collapse(biom_file,0,cog=TRUE)
cog_gene_names <- lapply(cog_gene_names,function(x) paste0(x,collapse=' | '))

DF <- data.frame(lasso=coef_k) %>%
   mutate(assoc=ifelse(lasso>0,'CD+',ifelse(lasso<0,'CD-','ZERO')),
          lasso_scaled=coef_k/sd(coef_k),
          lasso_ranked=rank(coef_k))
rownames(DF) <- colnames(cog_risk)

COG <- otu_table(cog_risk,taxa_are_rows=TRUE)
SAMP <- sample_data(DF)
PS <- phyloseq(COG,SAMP)

DS2 <- phyloseq_to_deseq2(PS, ~ lasso_ranked)
DS2 <- DESeq2::DESeq(DS2, test="Wald", fitType="parametric")
res <- DESeq2::results(DS2, cooksCutoff=FALSE, pAdjustMethod='BH')
alpha <- .01
sigtab <- res[which(res$padj < alpha), ]

cog_sig <- cog_risk[rownames(sigtab),]
rownames(cog_sig) <- cog_gene_names[rownames(cog_sig)]

cog_names <- rownames(cog_sig) 
names(cog_names) <- rownames(sigtab)

misc_filter <- !(rownames(cog_sig) %in% c('None','hypothetical protein','Uncharacterized conserved protein','Uncharacterized conserved protein'))
cog_sig <- cog_sig[misc_filter,]
cog_names <- cog_names[misc_filter]
####################################
####################################
















# pws <- table(unlist(sapply(kegg_metadata_risk, identity)))
# kegg_table <- kegg_risk
# kegg_table_pws <- sapply(kegg_metadata_risk[rownames(kegg_table)], function(x) names(which.max(pws[x])))
# kegg_table <- data.frame(kegg_table,pw=kegg_table_pws) %>%
#    group_by(pw) %>%
#    summarise_each(funs(sum)) 
# pw_names <- as.character(kegg_table$pw)
# kegg_table <- data.frame(dplyr::select(kegg_table,-pw))
# #kegg_table <- log10(kegg_table+1)
# #kegg_table <- t(t(kegg_table)/colSums(kegg_table))
# rownames(kegg_table) <- pw_names
# 
# filt_idx <- rowSums(kegg_table) != 0
# kegg_table <- kegg_table[filt_idx,]