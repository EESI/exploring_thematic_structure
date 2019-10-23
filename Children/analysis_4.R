library(biom)
library(readr)
library(tidyr)
library(dplyr)
library(fastICA)
library(randomForest)
library(kernlab)
library(Rcpp)
library(parallel)
library(foreach)
library(ape)
library(phyloseq)
library(doParallel)

collapse_table <- function(otu_taxa,data_mat,metadata,taxon,RA=TRUE){
  otu_taxa <- otu_taxa[rownames(data_mat),]
  target_taxon <- data.frame(full=apply(otu_taxa[,1:taxon],1,function(i) paste0(i,collapse='')),genus=otu_taxa[,taxon])
  target_length <- nrow(unique(target_taxon))
  taxon_table <- data.frame(data_mat,taxon=target_taxon$full) %>%
    group_by(taxon) %>%
    summarise_each(funs(sum)) %>%
    dplyr::select(-taxon)
  taxon_table <- t(taxon_table)
  
  if (ncol(taxon_table) != target_length & any(!(metadata$SampleID == rownames(taxon_table)))){
    print("Data not properly formated.")
    break
  }
  
  if (RA == TRUE){
    taxon_table <- t(apply(taxon_table,1,function(i) i/sum(i)))
  }
  
  return(taxon_table)
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

biom_file <- read_biom('~/Qiime/Children/Processing/01-raw/otus/sequences_filtered_100nt/otu_table.biom')
otu_taxa_xavier <- observation_metadata(biom_file)
data_mat_xavier <- as.matrix(biom_data(biom_file))
metadata_xavier <- read_delim('~/Qiime/Children/Processing/01-raw/metadata.txt',delim='\t')
sample_info_xavier <- read_csv('~/Qiime/Children/xavier_sample_info.csv')

biom_file <- read_biom("~/Qiime/Morgan/Processing/01-raw/otus/sequences_filtered_100nt/otu_table.biom")
otu_taxa_morgan <- observation_metadata(biom_file)
data_mat_morgan <- as.matrix(biom_data(biom_file))
metadata_morgan <- read_delim("~/Qiime/Morgan/Processing/01-raw/metadata.txt",delim="\t",quote='\"')

biom_file <- read_biom("~/Qiime/AG/Processing/01-raw/otus/sequences_filtered_bloomed_100nt/otu_table.biom")
otu_taxa_ag <- observation_metadata(biom_file)
data_mat_ag <- as.matrix(biom_data(biom_file))
metadata_ag <- read_delim("~/Qiime/AG/Processing/01-raw/metadata.txt",delim="\t",quote='\"')

metadata_xavier <- metadata_xavier %>%
  dplyr::select(SampleID = ends_with('SampleID'), everything()) %>%
  mutate(SAMPLE_NAME = ifelse(STRAIN == "missing", SAMPLE_NAME, STRAIN),
         SAMPLE_NAME = str_replace(SAMPLE_NAME,'-','')) %>%
  dplyr::select(SampleID,SAMPLE_NAME,ISOLATION_SOURCE)

metadata_xavier$COHORT <- ifelse(metadata_xavier$SAMPLE_NAME %in% sample_info_xavier$sample, "RISK","Other")
metadata_xavier <- metadata_xavier[metadata_xavier$SampleID %in% colnames(data_mat_xavier),] # make sure data_mat_xavier and metadata_xavier have same SRRs
metadata_xavier <- data.frame(metadata_xavier,
                              sample_info_xavier[match(metadata_xavier$SAMPLE_NAME,sample_info_xavier$sample),]) # combine metadata_xavier and sample_data

### there are nearly 700 subjects from RISK (about 1000 samples), which is about 460 CH and 230 None. See this based on subject, not sample name.
### the paper says that the RISK cohort was merged with two other cohorts giving about 1750 samples. These are likely the remaining samples not 
### accounted for in the csv file (sample_info_xavier) so they should be coded as CD.

metadata_xavier <- metadata_xavier %>%
  mutate(ISOLATION_SOURCE = ifelse(ISOLATION_SOURCE == "missing", sample_location, ISOLATION_SOURCE),
         ISOLATION_SOURCE = str_to_title(ISOLATION_SOURCE),
         RACE = str_to_title(race),
         STUDY = 'Xavier',
         AGE = 'Minor',
         DIAGNOSIS = ifelse(COHORT == "Other", "CD", Diagnosis)) %>%
  dplyr::select(SampleID, SAMPLE_NAME, ISOLATION_SOURCE, DIAGNOSIS, RACE, STUDY, AGE)

data_mat_xavier <- data_mat_xavier[,metadata_xavier$SampleID]

metadata_ag <- metadata_ag %>%
  dplyr::select(SampleID = ends_with('SampleID'), everything()) %>%
  mutate(AGE = as.numeric(AGE_CORRECTED)) %>%
  filter(BODY_SITE %in% "UBERON:feces",
         AGE > 0,
         IBD %in% c("Diagnosed by a medical professional (doctor, physician assistant)", 
                    "I do not have this condition")) %>%
  dplyr::select(SampleID,RACE,AGE,IBD) %>%
  mutate(ISOLATION_SOURCE = 'Stool',
         AGE = ifelse(AGE >= 17, 'Adult',
                      ifelse(AGE < 17, 'Minor', NA)),
         DIAGNOSIS = ifelse(IBD == "Diagnosed by a medical professional (doctor, physician assistant)",'IBD','Not IBD'),
         STUDY = 'AG',
         SAMPLE_NAME = SampleID) %>%
  dplyr::select(SampleID, SAMPLE_NAME, ISOLATION_SOURCE, DIAGNOSIS, RACE, STUDY, AGE)

metadata_ag <- metadata_ag[metadata_ag$SampleID %in% colnames(data_mat_ag),] # make sure data_mat_ag and metadata_ag have same SRRs
data_mat_ag <- data_mat_ag[,metadata_ag$SampleID] # make sure data_mat_ag and metadata_ag are aligned the same
data_mat_ag <- data_mat_ag[rowSums(data_mat_ag) != 0,] # remove OTUs with no occurances

metadata_morgan <- metadata_morgan %>%
  dplyr::select(SampleID = ends_with('SampleID'), everything()) %>%
  mutate(ISOLATION_SOURCE = ifelse(ISOLATION_SOURCE == 'missing','Stool',ISOLATION_SOURCE),
         DIAGNOSIS = 'IBD',
         RACE = NA,
         STUDY = 'Morgan',
         AGE = 'Adult') %>%
  dplyr::select(SampleID, SAMPLE_NAME, ISOLATION_SOURCE, DIAGNOSIS, RACE, STUDY, AGE)

metadata_morgan <- metadata_morgan[metadata_morgan$SampleID %in% colnames(data_mat_morgan),] # make sure data_mat_morgan and metadata_morgan have same SRRs
data_mat_morgan <- data_mat_morgan[,metadata_morgan$SampleID]




otus_ag <- rownames(data_mat_ag)
otus_xavier <- rownames(data_mat_xavier)
otus_morgan <- rownames(data_mat_morgan)
samples_ag <- colnames(data_mat_ag)
samples_xavier <- colnames(data_mat_xavier)
samples_morgan <- colnames(data_mat_morgan)

otus <- unique(c(otus_ag,otus_xavier,otus_morgan))
samples <- unique(c(samples_ag,samples_xavier,samples_morgan))

merged <- matrix(0,length(otus),length(samples),dimnames=list(otus,samples))

merged[rownames(data_mat_ag),colnames(data_mat_ag)] <- data_mat_ag
merged[rownames(data_mat_xavier),colnames(data_mat_xavier)] <- data_mat_xavier
merged[rownames(data_mat_morgan),colnames(data_mat_morgan)] <- data_mat_morgan

metadata <- data.frame(rbind(dplyr::select(metadata_ag,SampleID, SAMPLE_NAME, ISOLATION_SOURCE, DIAGNOSIS, RACE, STUDY, AGE),
                  dplyr::select(metadata_xavier,SampleID, SAMPLE_NAME, ISOLATION_SOURCE, DIAGNOSIS, RACE, STUDY, AGE),
                  dplyr::select(metadata_morgan,SampleID, SAMPLE_NAME, ISOLATION_SOURCE, DIAGNOSIS, RACE, STUDY, AGE)))
rownames(metadata) <- metadata$SampleID

otu_lookup <- unique(data.frame(rbind(otu_taxa_ag,otu_taxa_xavier,otu_taxa_morgan),
                                id=c(rownames(otu_taxa_ag),rownames(otu_taxa_xavier),rownames(otu_taxa_morgan)))) %>%
  dplyr::select(-id)
otu_lookup <- otu_lookup[rownames(merged),]
colnames(otu_lookup) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")



ra <- collapse_table(otu_lookup,merged,metadata,7)
model_out <- run_models(metadata,ra,studies=c('Xavier'),disease='CD')
model_out <- run_models(metadata,ra,studies=c('Xavier','AG'),disease='IBD')
model_out <- run_models(metadata,ra,studies=c('Xavier','AG'),disease='CD')
model_out <- run_models(metadata,ra,studies=c('Xavier','Morgan'),disease='IBD')
model_out <- run_models(metadata,ra,studies=c('Xavier','AG','Morgan'),disease='IBD')



OTU <- otu_table(merged,taxa_are_rows=TRUE)
TAX <- tax_table(as.matrix(otu_lookup))
META <- sample_data(metadata)
physeq <- phyloseq(OTU,TAX)
TREE <- rtree(ntaxa(physeq), rooted = TRUE, tip.label = taxa_names(physeq))
physeq <- phyloseq(OTU,TAX,META,TREE)


# cl <- makeCluster(50)
# registerDoParallel(cl)
# 
# ord1 <- ordinate(physeq,method='NMDS')
# plot_ordination(physeq,ord1,color="STUDY",shape="DIAGNOSIS")


wh0 <- genefilter_sample(physeq, filterfun_sample(function(x) x > 5), A = 0.1 * nsamples(physeq))
physeq1 <- prune_taxa(wh0, physeq)

wh1 <- rowSums(otu_table(physeq1)) != 0
physeq1 <- prune_taxa(wh1, physeq1)

wh2 <- colSums(otu_table(physeq1)) != 0
physeq1 <- prune_samples(wh2, physeq1)

physeq1 <- transform_sample_counts(physeq1, function(x) 1e+06 * x/sum(x))

phylum_sum <- tapply(taxa_sums(physeq1), tax_table(physeq1)[, "Phylum"], sum, na.rm = TRUE)
top5phyla <- names(sort(phylum_sum, TRUE))[1:5]

physeq2 <- prune_taxa((tax_table(physeq1)[, "Phylum"] %in% top5phyla), physeq1)

ord <- ordinate(physeq2, "NMDS", "bray")

p1 <- plot_ordination(physeq1, ord, type = "taxa", color = "Phylum", title = "taxa")
print(p1)

p2 <- plot_ordination(physeq1, ord, type = "samples", color = "STUDY", shape = "DIAGNOSIS")
print(p2)

p3 <- plot_ordination(physeq1, ord, type = "samples", color = "DIAGNOSIS", shape = "STUDY") +
  facet_wrap(~STUDY) + geom_point(alpha=.1)
print(p3)


test <- metaMDS(otu_table(physeq2),k=2,distance='jaccard')

data.frame(nmds1=test$species[,1],nmds2=test$species[,2],
           dx=sample_data(physeq2)$DIAGNOSIS,
           study=sample_data(physeq2)$STUDY) %>%
  ggplot(aes(nmds1,nmds2,colour=dx)) + 
  geom_hline(yintercept=0) + geom_vline(xintercept=0) +
    geom_point(alpha=.5) + facet_wrap(~study)

data.frame(nmds1=test$species[,1],nmds2=test$species[,2],
           age=sample_data(physeq2)$AGE,
           study=sample_data(physeq2)$STUDY) %>%
  ggplot(aes(nmds1,nmds2,colour=age,alpha=age)) + 
  geom_hline(yintercept=0) + geom_vline(xintercept=0) +
  geom_point() + facet_wrap(~study)



library("doParallel")
registerDoParallel(cores=60)

ord <- ordinate(physeq2, "PCoA", distance="bray", parallel=TRUE, k=2)

p1 <- plot_ordination(physeq1, ord, type = "taxa", color = "Phylum", title = "taxa")
print(p1)

p2 <- plot_ordination(physeq1, ord, type = "samples", color = "STUDY", shape = "DIAGNOSIS")
print(p2)

p3 <- plot_ordination(physeq1, ord, type = "samples", color = "DIAGNOSIS", shape = "STUDY") +
  facet_wrap(~STUDY) + geom_point(alpha=.1)
print(p3)


ord1 <- ordinate(physeq2, "PCoA", distance="jaccard", parallel=TRUE, k=2)

p1 <- plot_ordination(physeq1, ord1, type = "taxa", color = "Phylum", title = "taxa")
print(p1)

p2 <- plot_ordination(physeq1, ord1, type = "samples", color = "STUDY", shape = "DIAGNOSIS")
print(p2)

p3 <- plot_ordination(physeq1, ord1, type = "samples", color = "DIAGNOSIS") +
  geom_hline(yintercept=0,alpha=.5) + geom_vline(xintercept=0,alpha=.5) +
  facet_wrap(~STUDY) + geom_point(alpha=.1) + 
  theme(legend.position = 'bottom',
        legend.text = element_text(size=rel(1.5)),
        strip.text = element_text(size=rel(2)))
print(p3)

p4 <- plot_ordination(physeq1, ord1, type = "samples", color = "AGE") +
  geom_hline(yintercept=0,alpha=.5) + geom_vline(xintercept=0,alpha=.5) +
  facet_wrap(~STUDY) + geom_point(alpha=.1) + 
  theme(legend.position = 'bottom',
        legend.text = element_text(size=rel(1.5)),
        strip.text = element_text(size=rel(2)))
print(p4)

p5 <- plot_ordination(physeq1, ord1, type = "samples", color = "STUDY") +
  geom_hline(yintercept=0,alpha=.5) + geom_vline(xintercept=0,alpha=.5) +
  facet_wrap(~ISOLATION_SOURCE) + geom_point(alpha=.1) + 
  theme(legend.position = 'bottom',
        legend.text = element_text(size=rel(1.5)),
        strip.text = element_text(size=rel(2)))
print(p5)
