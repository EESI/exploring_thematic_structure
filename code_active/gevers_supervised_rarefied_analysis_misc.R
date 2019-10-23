source('~/Dropbox/stm_microbiome/code_active/stm_functions.R')
source('~/Dropbox/stm_microbiome/code_active/nav_froz_fxns_3.R') 
source('~/Dropbox/stm_microbiome/code_active/performance_1.R')
source('~/Dropbox/stm_microbiome/code_active/framework.R')

biom_file <- read_biom('~/Dropbox/stm_microbiome/Children/Processing/01-raw/simcomp_otus_q32_sim97/sequences_filtered/otu_table_cnNormalized.biom')
data_mat_xavier <- floor(as.matrix(biom_data(biom_file)))
biom_file <- read_biom('~/Dropbox/stm_microbiome/Children/Processing/01-raw/simcomp_otus_q32_sim97/sequences_filtered/otu_table.biom')
otu_taxa_xavier <- observation_metadata(biom_file)

metadata_xavier <- read_delim('~/Dropbox/stm_microbiome/Children/Processing/01-raw/metadata.txt',delim='\t')
sample_info_xavier <- read_csv('~/Dropbox/stm_microbiome/Children/xavier_sample_info.csv')

metadata_xavier <- metadata_xavier %>%
   dplyr::select(SampleID = ends_with('SampleID'), everything()) %>%
   mutate(SAMPLE_NAME = ifelse(STRAIN == "missing", SAMPLE_NAME, STRAIN),
          SAMPLE_NAME = str_replace(SAMPLE_NAME,'-','')) %>%
   dplyr::select(SampleID,SAMPLE_NAME,ISOLATION_SOURCE)

metadata_xavier$COHORT <- ifelse(metadata_xavier$SAMPLE_NAME %in% sample_info_xavier$sample, "RISK","Other")
metadata_xavier <- metadata_xavier[metadata_xavier$SampleID %in% colnames(data_mat_xavier),] # make sure data_mat_xavier and metadata_xavier have same SRRs
metadata_xavier <- data.frame(metadata_xavier,
                              sample_info_xavier[match(metadata_xavier$SAMPLE_NAME,sample_info_xavier$sample),]) # combine metadata_xavier and sample_data

metadata_xavier <- metadata_xavier %>%
   mutate(ISOLATION_SOURCE = ifelse(ISOLATION_SOURCE == "missing", sample_location, ISOLATION_SOURCE),
          ISOLATION_SOURCE = str_to_title(ISOLATION_SOURCE),
          RACE = str_to_title(race),
          STUDY = 'Xavier',
          AGE = 'Minor',
          DIAGNOSIS = ifelse(COHORT == "Other", "CD", Diagnosis)) %>%
   dplyr::select(SampleID, SAMPLE_NAME, ISOLATION_SOURCE, DIAGNOSIS, RACE, STUDY, AGE, COHORT, PCDAI, AB_exposure)
data_mat_xavier <- data_mat_xavier[,metadata_xavier$SampleID]
metadata_risk <- metadata_xavier %>% dplyr::filter(COHORT == 'RISK',
                                                   ISOLATION_SOURCE %in% c('Terminalileum Biopsy','Rectum Biopsy'))
metadata_risk <- metadata_risk %>% dplyr::filter(AB_exposure == 'NoAB')

metadata_risk <- metadata_risk %>%
   mutate(PCDAI = ifelse(is.na(PCDAI) & DIAGNOSIS == 'Not IBD',0,PCDAI)) %>%
   filter(!is.na(PCDAI)) %>%
   mutate(PCDAI_RANK = ranker(PCDAI))

data_mat_risk <- data_mat_xavier[,metadata_risk$SampleID]
data_mat_risk <- data_mat_risk[rowSums(data_mat_risk) != 0,]


OTU <- otu_table(data_mat_risk,taxa_are_rows=TRUE)
Nmin <- quantile(colSums(OTU),.2)
OTU <- rarefy_even_depth(OTU,
                         rngseed=2432, #6546
                         sample.size=Nmin,
                         replace=FALSE,
                         trimOTUs=TRUE,
                         verbose=TRUE)

data_mat_risk <- t(OTU@.Data)


read_min <- 2
samp_min <- .99
taxa_filter <- colSums(data_mat_risk < read_min) > floor(samp_min * NROW(data_mat_risk)) 

cat('Removing ', sum(taxa_filter),  ' taxa (from ',NCOL(data_mat_risk),' taxa) that have fewer than ',
    read_min, ' reads in less than ',(1-samp_min)*100,'% of samples.\n',sep='')

data_mat_risk <- data_mat_risk[,!taxa_filter]

otu_taxa_xavier <- otu_taxa_xavier[colnames(data_mat_risk),]   
rownames(metadata_risk) <- metadata_risk$SampleID
metadata_risk <- metadata_risk[rownames(data_mat_risk),]




taxon_ids <- data.frame(long=colnames(data_mat_risk),
                        ids=paste0('otu',1:NCOL(data_mat_risk)),
                        row.names=paste0('otu',1:NCOL(data_mat_risk)))

counts <- list(meta_clean=metadata_risk,
               table_clean=data_mat_risk,
               ids=taxon_ids)
colnames(counts$table_clean) <- as.character(counts$ids$ids)

if (!identical(rownames(counts$table_clean),counts$meta_clean$SampleID)) stop('Please make table sample names match metadata sample IDs!\n')

vocab <- as.character(counts$ids$ids)
docs <- lapply(1:NROW(counts$table_clean),function(i) reshape_doc(counts$table_clean[i,],vocab))
meta <- counts$meta_clean

cat('Total docs:',length(docs),'\n',
    'Total vocab:',length(vocab),'\n',
    'Total Meta:',nrow(meta),'x',ncol(meta),'\n',
    'Total DX:',names(table(meta$DIAGNOSIS)),'=',c(table(meta$DIAGNOSIS)),'\n',
    'Total IS:',names(table(meta$ISOLATION_SOURCE)),'=',c(table(meta$ISOLATION_SOURCE)),'\n')


K <- 85
fit1 <- stm(prevalence=~DIAGNOSIS,
            documents=docs,vocab=vocab,K=K,
            data=meta,
            max.em.its=500,init.type='Spectral',
            seed=7676, #645647
            verbose=TRUE,reportevery=25)

fit2 <- stm(prevalence=~PCDAI,
            documents=docs,vocab=vocab,K=K,
            data=meta,
            max.em.its=500,init.type='Spectral',
            seed=7676, #645647
            verbose=TRUE,reportevery=25)


fit3 <- stm(prevalence=~PCDAI,
            content=~ISOLATION_SOURCE,
            documents=docs,vocab=vocab,K=K,
            data=meta,
            max.em.its=500,init.type='Spectral',
            seed=7676, #645647
            verbose=TRUE,reportevery=25)


fit4 <- stm(prevalence=~PCDAI*ISOLATION_SOURCE,
            documents=docs,vocab=vocab,K=K,
            data=meta,
            max.em.its=500,init.type='Spectral',
            seed=7676, #645647
            verbose=TRUE,reportevery=25)


fit5 <- stm(prevalence=~PCDAI*ISOLATION_SOURCE,
            content=~ISOLATION_SOURCE,
            documents=docs,vocab=vocab,K=K,
            data=meta,
            max.em.its=500,init.type='Spectral',
            seed=7676, #645647
            verbose=TRUE,reportevery=25)

# 
# fit6 <- stm(prevalence=~DIAGNOSIS,
#             content=~DIAGNOSIS,
#             documents=docs,vocab=vocab,K=K,
#             data=meta,
#             max.em.its=500,init.type='Spectral',
#             seed=7676, #645647
#             verbose=TRUE,reportevery=25)
# 
# 





eff1 <- estimateEffect(1:K ~ DIAGNOSIS, fit1, meta, uncertainty='Global')
eff_plot1 <- plot.estimateEffect(eff1, 'DIAGNOSIS', model=fit1, topics=1:K, method='difference',cov.value1='CD',cov.value2='Not IBD')

coef_k <- ifelse(sapply(eff_plot1$cis,function(x) sum(sign(x))) != 0, unlist(eff_plot1$means), 0)
coef_k <- coef_k/sd(coef_k)

theta <- fit1$theta
rownames(theta) <- meta$SampleID
colnames(theta) <- paste0('T',1:K)
theta_sort <- meta %>%
   arrange(ISOLATION_SOURCE,DIAGNOSIS,desc(PCDAI)) %>%
   select(SampleID,DIAGNOSIS) 
theta <- theta[theta_sort$SampleID,]
rownames(theta) <- theta_sort$DIAGNOSIS

pdf('~/MiscOut/super_iso_dx_01.pdf', height=10, width=20)
top_ord <- plot_big_heatmap2(theta,coef_k=coef_k,z=TRUE,col1='cyan',col2='black',col3='orange',cexCol=1.5,main='',scale='column')
dev.off()





slope <- function(x){
   (x[length(x)]-x[1])/(length(x)-1)
}




eff2 <- estimateEffect(1:K ~ PCDAI, fit2, meta, uncertainty='Global')
eff_plot2 <- plot.estimateEffect(eff2, 'PCDAI', model=fit2, topics=1:K, method='continuous')

coef_k <- sapply(eff_plot2$means,slope)
coef_k <- coef_k/sd(coef_k)

theta <- fit2$theta
rownames(theta) <- meta$SampleID
colnames(theta) <- paste0('T',1:K)
theta_sort <- meta %>%
   arrange(ISOLATION_SOURCE,DIAGNOSIS,desc(PCDAI)) %>%
   select(SampleID,DIAGNOSIS) 
theta <- theta[theta_sort$SampleID,]
rownames(theta) <- theta_sort$DIAGNOSIS

pdf('~/MiscOut/super_iso_dx_02.pdf', height=10, width=20)
top_ord <- plot_big_heatmap2(theta,coef_k=coef_k,z=TRUE,col1='cyan',col2='black',col3='orange',cexCol=1.5,main='',scale='column')
dev.off()





eff3 <- estimateEffect(1:K ~ PCDAI, fit3, meta, uncertainty='Global')
eff_plot3 <- plot.estimateEffect(eff3, 'PCDAI', model=fit3, topics=1:K, method='continuous')

coef_k <- sapply(eff_plot3$means,slope)
coef_k <- coef_k/sd(coef_k)

theta <- fit3$theta
rownames(theta) <- meta$SampleID
colnames(theta) <- paste0('T',1:K)
theta_sort <- meta %>%
   arrange(ISOLATION_SOURCE,DIAGNOSIS,desc(PCDAI)) %>%
   select(SampleID,DIAGNOSIS) 
theta <- theta[theta_sort$SampleID,]
rownames(theta) <- theta_sort$DIAGNOSIS

pdf('~/MiscOut/super_iso_dx_03.pdf', height=10, width=20)
top_ord <- plot_big_heatmap2(theta,coef_k=coef_k,z=TRUE,col1='cyan',col2='black',col3='orange',cexCol=1.5,main='',scale='column')
dev.off()









qplot(fit3$beta$logbeta[[1]][1,],fit3$beta$logbeta[[2]][1,])











eff4 <- estimateEffect(1:K ~ PCDAI*ISOLATION_SOURCE, fit4, meta, uncertainty='Global')
eff_plot4 <- plot.estimateEffect(eff4, covariate='PCDAI', model=fit4, topics=1:K, method='continuous')

coef_k <- sapply(eff_plot4$means,slope)
coef_k <- coef_k/sd(coef_k)

theta <- fit4$theta
rownames(theta) <- meta$SampleID
colnames(theta) <- paste0('T',1:K)
theta_sort <- meta %>%
   arrange(ISOLATION_SOURCE,DIAGNOSIS,desc(PCDAI)) %>%
   select(SampleID,DIAGNOSIS) 
theta <- theta[theta_sort$SampleID,]
rownames(theta) <- theta_sort$DIAGNOSIS

pdf('~/MiscOut/super_iso_dx_04.pdf', height=10, width=20)
top_ord <- plot_big_heatmap2(theta,coef_k=coef_k,z=TRUE,col1='cyan',col2='black',col3='orange',cexCol=1.5,main='',scale='column')
dev.off()


topic <- 26
plot.estimateEffect(eff4, covariate='PCDAI', model=fit4, topics=topic, method='continuous',
                                 moderator='ISOLATION_SOURCE',moderator.value='Terminalileum Biopsy',linecol='yellow')
plot.estimateEffect(eff4, covariate='PCDAI', model=fit4, topics=topic, method='continuous',
                    moderator='ISOLATION_SOURCE',moderator.value='Rectum Biopsy',linecol='brown',add=TRUE)






eff5 <- estimateEffect(1:K ~ PCDAI*ISOLATION_SOURCE, fit5, meta, uncertainty='Global')
eff_plot5 <- plot.estimateEffect(eff5, covariate='PCDAI', model=fit5, topics=1:K, method='continuous')

coef_k <- sapply(eff_plot5$means,slope)
coef_k <- coef_k/sd(coef_k)

theta <- fit5$theta
rownames(theta) <- meta$SampleID
colnames(theta) <- paste0('T',1:K)
theta_sort <- meta %>%
   arrange(ISOLATION_SOURCE,DIAGNOSIS,desc(PCDAI)) %>%
   select(SampleID,DIAGNOSIS) 
theta <- theta[theta_sort$SampleID,]
rownames(theta) <- theta_sort$DIAGNOSIS

pdf('~/MiscOut/super_iso_dx_05.pdf', height=10, width=20)
top_ord <- plot_big_heatmap2(theta,coef_k=coef_k,z=TRUE,col1='cyan',col2='black',col3='orange',cexCol=1.5,main='',scale='column')
dev.off()


topic <- c(3,23)
plot.estimateEffect(eff5, covariate='PCDAI', model=fit5, topics=topic, method='continuous')
plot.estimateEffect(eff5, covariate='PCDAI', model=fit5, topics=topic, method='continuous',
                    moderator='ISOLATION_SOURCE',moderator.value='Rectum Biopsy',linecol='brown',add=TRUE)









eff6 <- estimateEffect(1:K ~ DIAGNOSIS, fit6, meta, uncertainty='Global')
eff_plot6 <- plot.estimateEffect(eff6, 'DIAGNOSIS', model=fit6, topics=1:K, method='difference',cov.value1='CD',cov.value2='Not IBD')

coef_k <- ifelse(sapply(eff_plot6$cis,function(x) sum(sign(x))) != 0, unlist(eff_plot6$means), 0)
coef_k <- coef_k/sd(coef_k)

theta <- fit6$theta
rownames(theta) <- meta$SampleID
colnames(theta) <- paste0('T',1:K)
theta_sort <- meta %>%
   arrange(DIAGNOSIS,desc(PCDAI)) %>%
   select(SampleID,DIAGNOSIS) 
theta <- theta[theta_sort$SampleID,]
rownames(theta) <- theta_sort$DIAGNOSIS

pdf('~/MiscOut/super_iso_dx_06.pdf', height=10, width=20)
top_ord <- plot_big_heatmap2(theta,coef_k=coef_k,z=TRUE,col1='cyan',col2='black',col3='orange',cexCol=1.5,main='',scale='column')
dev.off()