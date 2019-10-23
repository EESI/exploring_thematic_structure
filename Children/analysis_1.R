library(biom)
library(readr)
library(tidyr)
library(dplyr)
library(fastICA)
library(randomForest)
library(kernlab)

biom_file <- read_biom('~/Qiime/Children/Processing/01-raw/otus/sequences_filtered_100nt/otu_table.biom')
otu_taxa <- observation_metadata(biom_file)
data_mat <- as.matrix(biom_data(biom_file))
metadata <- read_delim('~/Qiime/Children/Processing/01-raw/metadata.txt',delim='\t')
sample_info <- read_csv('~/Qiime/Children/xavier_sample_info.csv')

dim(data_mat)
dim(otu_taxa)
max(data_mat)

# prepare metadata
metadata <- metadata %>%
  dplyr::select(SampleID = ends_with('SampleID'), everything()) %>%
  mutate(STRAIN = str_replace(STRAIN,'-',''))
metadata <- metadata[metadata$SampleID %in% colnames(data_mat),] # make sure data_mat and metadata have same SRRs
otu_taxa <- otu_taxa[rownames(data_mat),] # make sure data_mat and otu_taxa are aligned the same
metadata <- data.frame(metadata,sample_info[match(metadata$STRAIN,sample_info$sample),]) # combine metadata and sample_data
rownames(metadata) <- metadata$SampleID
metadata <- metadata[colnames(data_mat),] # make sure data_mat and metadata are aligned the same


# group by taxonomic level and return relative abundances
collapse_table <- function(otu_taxa,data_mat,metadata,taxon,RA=TRUE){
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

ra <- collapse_table(otu_taxa,data_mat,metadata,7)


# pca <- prcomp(ra, scale.=FALSE)
# ica <- fastICA(ra,n.comp=2)
# out <- data.frame(pc1=pca$x[,1],pc2=pca$x[,2],
#            ic1=ica$S[,1],ic2=ica$S[,2],
#            location=metadata$ISOLATION_SOURCE,
#            ab=metadata$AB_exposure,
#            ibd=metadata$Diagnosis,
#            race=metadata$race)
# 
# ggplot(out, aes(x=pc1,y=pc2,colour=ibd)) + geom_point(alpha=.5) 
# ggplot(out, aes(x=ic1,y=ic2,colour=ibd)) + geom_point(alpha=.5) 

run_models <- function(metadata,ra,prop=.6){
  feature_ibd <- unlist(metadata %>% filter(Diagnosis %in% c("Not IBD","CD")) %>% dplyr::select(Diagnosis))
  feature_ibd0 <- which(feature_ibd == 'Not IBD')
  feature_ibd1 <- which(feature_ibd == 'CD')
  undersamp_idx <- c(sample(feature_ibd0,min(table(feature_ibd))),sample(feature_ibd1,min(table(feature_ibd))))
  feature_ibd <- feature_ibd[undersamp_idx]
  ra_us <- ra[undersamp_idx,]
  
  train_idx <- caret::createDataPartition(feature_ibd,times=1,p=prop,list=FALSE)
  train_table <- ra_us[train_idx,]
  train_labels <- feature_ibd[train_idx]
  test_table <- ra_us[-train_idx,]
  test_labels <- feature_ibd[-train_idx]
  
  print(table(train_labels))
  print(table(test_labels))
  
  out_rf <- randomForest(x=train_table,y=as.factor(train_labels),
                         xtest=test_table,ytest=as.factor(test_labels),
                         ntree=100,importance=FALSE)
  print((mean(out_rf$test$predicted == test_labels)))
  
  
  svm_train <- ksvm(train_table, train_labels,
                    type="C-svc",kernel="vanilladot",
                    C=100,scale = FALSE)
  print((mean(predict(svm_train,test_table,type="response") == test_labels)))
  
  return(list(rf=out_rf,svm=svm_train))
}

model_out <- run_models(metadata,ra)

biom_file <- read_biom("~/American-Gut/trunk/data/AG/AG_100nt.biom")
otu_taxa_ag <- observation_metadata(biom_file)
data_mat_ag <- as.matrix(biom_data(biom_file))
metadata_ag <- read_delim("~/American-Gut/trunk/data/AG/AG_100nt.txt",delim="\t",quote='\"')

dim(data_mat_ag)
dim(otu_taxa_ag)
max(data_mat_ag)

# # prepare metadata
# metadata_ag <- metadata_ag %>%
#   dplyr::select(SampleID = ends_with('SampleID'), everything()) %>%
#   mutate(AGE = as.numeric(AGE)) %>%
#   filter(BODY_SITE %in% "UBERON:feces",
#          AGE >= 18,
#          !is.na(AGE)) %>%
#   dplyr::select(SampleID,IBD,RACE,AGE) %>%
#   mutate(ISOLATION_SOURCE = 'Stool',
#          Diagnosis = ifelse(IBD == 'I do not have IBD','Not IBD',
#                             ifelse(IBD == "Crohn's disease","CD",
#                                    ifelse(IBD == "Ulcerative colitis","UC",NA)))) %>%
#   dplyr::select(-IBD)
# 
# metadata_ag <- metadata_ag[metadata_ag$SampleID %in% colnames(data_mat_ag),] # make sure data_mat_ag and metadata_ag have same SRRs
# otu_taxa_ag <- otu_taxa_ag[rownames(data_mat_ag),] # make sure data_mat_ag and otu_taxa_ag are aligned the same
# rownames(metadata_ag) <- metadata_ag$SampleID
# data_mat_ag <- data_mat_ag[,rownames(metadata_ag)] # make sure data_mat_ag and metadata_ag are aligned the same
# 
# 
# ra_ag <- collapse_table(otu_taxa_ag,data_mat_ag,metadata_ag,7)
# model_out <- run_models(metadata_ag,ra_ag)

# prepare metadata
metadata_ag <- metadata_ag %>%
  dplyr::select(SampleID = ends_with('SampleID'), everything()) %>%
  mutate(AGE = as.numeric(AGE)) %>%
  filter(BODY_SITE %in% "UBERON:feces",
         AGE >= 18,
         !is.na(AGE),
         IBD %in% c("Crohn's disease", "Ulcerative colitis")) %>%
  dplyr::select(SampleID,IBD,RACE,AGE) %>%
  mutate(ISOLATION_SOURCE = 'Stool',
         Diagnosis = ifelse(IBD == "Crohn's disease","CD",
                            ifelse(IBD == "Ulcerative colitis","UC",NA))) %>%
  dplyr::select(-IBD)

metadata_ag <- metadata_ag[metadata_ag$SampleID %in% colnames(data_mat_ag),] # make sure data_mat_ag and metadata_ag have same SRRs
otu_taxa_ag <- otu_taxa_ag[rownames(data_mat_ag),] # make sure data_mat_ag and otu_taxa_ag are aligned the same
rownames(metadata_ag) <- metadata_ag$SampleID
data_mat_ag <- data_mat_ag[,rownames(metadata_ag)] # make sure data_mat_ag and metadata_ag are aligned the same

added_otus <- rownames(data_mat)[!(rownames(data_mat) %in% rownames(data_mat_ag))]
merge_mat <- matrix(0,
                    nrow(data_mat_ag) + length(added_otus),
                    ncol(data_mat_ag) + ncol(data_mat))
merge_row_names <- c(rownames(data_mat_ag),
                     added_otus)
rownames(merge_mat) <- merge_row_names
merge_mat[rownames(data_mat),1:ncol(data_mat)] <- data_mat
merge_mat[rownames(data_mat_ag),(ncol(data_mat)+1):(ncol(data_mat)+ncol(data_mat_ag))] <- data_mat_ag
colnames(merge_mat) <- c(colnames(data_mat),colnames(data_mat_ag))

merged_otu_taxa <- rbind(otu_taxa_ag,otu_taxa[added_otus,])
merged_metadata <- rbind(dplyr::select(metadata_ag,-RACE,-AGE),
                         dplyr::select(metadata,SampleID,ISOLATION_SOURCE,Diagnosis))

merged_otu_taxa <- merged_otu_taxa[rownames(merge_mat),] # make sure data_mat_ag and otu_taxa_ag are aligned the same
rownames(merged_metadata) <- merged_metadata$SampleID
merge_mat <- merge_mat[,rownames(merged_metadata)] # make sure data_mat_ag and metadata_ag are aligned the same

ra_ag_merged <- collapse_table(merged_otu_taxa,merge_mat,merged_metadata,7)
model_out <- run_models(merged_metadata,ra_ag_merged)
