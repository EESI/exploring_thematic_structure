library(randomForest)
library(kernlab)
library(caret)
library(tidyr)
library(dplyr)
library(ggplot2)
library(glmnet)
library(stm)
library(parallel)
library(foreach)
library(doParallel)

source('~/Qiime/Children/stm_functions.R')

rds_out <- readRDS('~/Qiime/Children/xavier_data_loading.rds')
otu_lookup <- rds_out$otu_lookup
data_mat_xavier <- rds_out$data_mat_xavier
metadata_xavier <- rds_out$metadata_xavier
kegg_metadata_xavier <- rds_out$kegg_metadata_xavier
kegg_mat_xavier <- rds_out$kegg_mat_xavier


#id_oversamp <- readRDS('~/Qiime/Children/xavier_oversamp_loading')
# this run is to see the pred accuracy without the bootstrap method I used before
# I want to see how much over fitting
# therefore I'm going to set aside the test set first, then bootstrap within the training set


uidx <- 468727 #unique_idx()
dump_dir <- paste0('~/Qiime/Children/stm_comp_oversamp_unique_',uidx,sep='')
dump_dir_files <- paste0(dump_dir,'/kegg',sep='')
dir.create(dump_dir_files, showWarnings = FALSE,recursive=TRUE)

### KEGG oversampled model

counts <- collapse_table(otu_lookup,kegg_mat_xavier,metadata_xavier,
                         taxon=NULL,method='rescaled',
                         kegg=TRUE, cap=TRUE,roof=10^5,
                         filtermin=FALSE,mc=5,ps=.1,
                         filtermax=FALSE,pw=.9,pa=.9)

counts$meta <- metadata_xavier[metadata_xavier$SampleID %in% rownames(counts$table),]
counts$meta_clean <- counts$meta[!is.na(counts$meta$ISOLATION_SOURCE),]
counts$table_clean <- counts$table[counts$meta_clean$SampleID,]
rownames(counts$meta_clean) <- counts$meta_clean$SampleID

set.seed(343)
test_id <- c(sample(unique(counts$meta_clean$SAMPLE_NAME[counts$meta_clean$DIAGNOSIS =='Not IBD']),101),
             sample(unique(counts$meta_clean$SAMPLE_NAME[counts$meta_clean$DIAGNOSIS =='CD']),101))
test_idx <- match(test_id,counts$meta_clean$SAMPLE_NAME)
n_oversamp <- abs(diff(table(counts$meta_clean[-test_idx,]$DIAGNOSIS)))-5
id_oversamp <- sample(counts$meta_clean[-test_idx,]$SampleID[counts$meta_clean[-test_idx,]$DIAGNOSIS=='Not IBD'],
                      n_oversamp,replace=TRUE)

counts$meta_clean <- rbind(counts$meta_clean,counts$meta_clean[id_oversamp,])
counts$table_clean <- rbind(counts$table_clean,counts$table_clean[id_oversamp,])

vocab <- as.character(counts$ids$ids)
docs <- lapply(1:nrow(counts$table_clean),function(i) reshape_doc(counts$table_clean[i,],vocab))
meta <- counts$meta_clean
meta$DIAGNOSIS <- as.factor(ifelse(meta$DIAGNOSIS=='CD',1,0))
meta$ISOLATION_SOURCE <- as.factor(meta$ISOLATION_SOURCE)
meta$COHORT <- as.factor(meta$COHORT)

prep <- prepDocuments(docs,vocab,meta,lower.thresh=0) ### removes terms apparently
counts$table_clean <- counts$table_clean[prep$meta$SampleID,prep$vocab]


# set.seed(43)
# idx_train <- (1:nrow(meta))[-match(test_id,meta$SAMPLE_NAME)]
# out_rf_dx_kegg <- randomForest(x=counts$table_clean[idx_train,],y=meta$DIAGNOSIS[idx_train],
#                        xtest=counts$table_clean[-idx_train,],ytest=meta$DIAGNOSIS[-idx_train],
#                        ntree=100,importance=TRUE,proximity=TRUE)
# saveRDS(list(model=out_rf_dx_kegg,idx_train=idx_train),
#         file.path(dump_dir,paste0(deparse(substitute(out_rf_dx_kegg)),'.rds')))
# 
# set.seed(5455)
# out_rf_loc_kegg <- randomForest(x=counts$table_clean[idx_train,],y=meta$ISOLATION_SOURCE[idx_train],
#                                 xtest=counts$table_clean[-idx_train,],ytest=meta$ISOLATION_SOURCE[-idx_train],
#                                 ntree=100,importance=TRUE,proximity=TRUE)
# 
# saveRDS(list(model=out_rf_loc_kegg,idx_train=idx_train),
#         file.path(dump_dir,paste0(deparse(substitute(out_rf_loc_kegg)),'.rds')))

model <- readRDS(file.path(dump_dir,paste0(deparse(substitute(out_rf_dx_kegg)),'.rds')))        
out_rf_dx_kegg <- model$model
idx_train_dx <- model$idx_train

model <- readRDS(file.path(dump_dir,paste0(deparse(substitute(out_rf_loc_kegg)),'.rds')))        
out_rf_loc_kegg <- model$model
idx_train_loc <- model$idx_train

cf_dx <- confusionMatrix(out_rf_dx_kegg$test$predicted,meta$DIAGNOSIS[-idx_train_dx])
cf_loc <- confusionMatrix(out_rf_loc_kegg$test$predicted,meta$ISOLATION_SOURCE[-idx_train_loc])




cat('Prediction Report',rep('\n',3),sep='',file=file.path(dump_dir_files,'out.txt'))

cat('Labels: Diagnosis\nModel:\tRF\nCV:\t70/30 (preBS)\n',sep='',file=file.path(dump_dir_files,'out.txt'),append=TRUE)
cat('m (splits):\t',out_rf_dx_kegg$mtry,' (sqrt(p))\n',sep='',file=file.path(dump_dir_files,'out.txt'),append=TRUE)
cat('n (trees):\t',out_rf_dx_kegg$ntree,rep('\n',3),sep='',file=file.path(dump_dir_files,'out.txt'),append=TRUE)
cat(paste0(names(cf_dx$overall),'\t'),'\n',file=file.path(dump_dir_files,'out.txt'),append=TRUE)
cat(paste0(round(cf_dx$overall,5),'\t'),'\n',file=file.path(dump_dir_files,'out.txt'),append=TRUE)

cat(rep('\n',3),sep='',file=file.path(dump_dir_files,'out.txt'),append=TRUE)
write.table(cf_dx$table,file=file.path(dump_dir_files,'out.txt'),quote=FALSE,append=TRUE)

cat(rep('\n',3),sep='',file=file.path(dump_dir_files,'out.txt'),append=TRUE)
cat(paste0(names(cf_dx$byClass),'\t'),'\n',file=file.path(dump_dir_files,'out.txt'),append=TRUE)
cat(paste0(round(cf_dx$byClass,5),'\t'),'\n',file=file.path(dump_dir_files,'out.txt'),append=TRUE)

cat(rep('\n',10),sep='',file=file.path(dump_dir_files,'out.txt'),append=TRUE)

cat('Labels: Isolation Source\nModel:\tRF\nCV:\t70/30 (preBS)\n',sep='',file=file.path(dump_dir_files,'out.txt'),append=TRUE)
cat('m (splits):\t',out_rf_loc_kegg$mtry,' (sqrt(p))\n',sep='',file=file.path(dump_dir_files,'out.txt'),append=TRUE)
cat('n (trees):\t',out_rf_loc_kegg$ntree,rep('\n',3),sep='',file=file.path(dump_dir_files,'out.txt'),append=TRUE)
cat(paste0(names(cf_loc$overall),'\t'),'\n',file=file.path(dump_dir_files,'out.txt'),append=TRUE)
cat(paste0(round(cf_loc$overall,5),'\t'),'\n',file=file.path(dump_dir_files,'out.txt'),append=TRUE)

cat(rep('\n',3),sep='',file=file.path(dump_dir_files,'out.txt'),append=TRUE)
write.table(cf_loc$table,file=file.path(dump_dir_files,'out.txt'),quote=FALSE,append=TRUE)

cat(rep('\n',3),sep='',file=file.path(dump_dir_files,'out.txt'),append=TRUE)
cat('\t',paste0(colnames(cf_loc$byClass),'\t'),'\n',quote=FALSE,file=file.path(dump_dir_files,'out.txt'),append=TRUE)
for (r in 1:nrow(cf_loc$byClass)){
  cat(rownames(cf_loc$byClass)[r],'\t',
      paste0(cf_loc$byClass[r,],'\t'),
      '\n',file=file.path(dump_dir_files,'out.txt'),append=TRUE)
}


pdf(file.path(dump_dir_files,'imp_dx.pdf'),width=12,height=12)
imp_dx <- plot_imp(out_rf_dx_kegg)
dev.off()

pdf(file.path(dump_dir_files,'prox_dx.pdf'),width=12,height=12)
plot_prox(out_rf_dx_kegg,meta$DIAGNOSIS[idx_train_dx])
dev.off()

pdf(file.path(dump_dir_files,'imp_loc.pdf'),width=12,height=12)
imp_loc <- plot_imp(out_rf_loc_kegg)
dev.off()

pdf(file.path(dump_dir_files,'prox_loc.pdf'),width=12,height=12)
plot_prox(out_rf_loc_kegg,meta$ISOLATION_SOURCE[idx_train_loc])
dev.off()

imp <- imp_dx$acc[1:500]
write.table(
data.frame(id1=names(imp[imp > 0]),
           id2=counts$ids[names(imp[imp > 0]),'long'],
           score=imp[imp > 0],
           pw=sapply(kegg_metadata_xavier[counts$ids[names(imp[imp > 0]),'long']],
                     function(x) paste0(x,collapse=' | '))),
file=file.path(dump_dir_files,'imp_acc_dx.txt'),quote=FALSE,append=TRUE)

imp <- imp_dx$gini[1:500]
write.table(
  data.frame(id1=names(imp[imp > 0]),
             id2=counts$ids[names(imp[imp > 0]),'long'],
             score=imp[imp > 0],
             pw=sapply(kegg_metadata_xavier[counts$ids[names(imp[imp > 0]),'long']],
                       function(x) paste0(x,collapse=' | '))),
  file=file.path(dump_dir_files,'imp_gini_dx.txt'),quote=FALSE,append=TRUE)

imp <- imp_loc$acc[1:500]
write.table(
  data.frame(id1=names(imp[imp > 0]),
             id2=counts$ids[names(imp[imp > 0]),'long'],
             score=imp[imp > 0],
             pw=sapply(kegg_metadata_xavier[counts$ids[names(imp[imp > 0]),'long']],
                       function(x) paste0(x,collapse=' | '))),
  file=file.path(dump_dir_files,'imp_acc_loc.txt'),quote=FALSE,append=TRUE)

imp <- imp_loc$gini[1:500]
write.table(
  data.frame(id1=names(imp[imp > 0]),
             id2=counts$ids[names(imp[imp > 0]),'long'],
             score=imp[imp > 0],
             pw=sapply(kegg_metadata_xavier[counts$ids[names(imp[imp > 0]),'long']],
                       function(x) paste0(x,collapse=' | '))),
  file=file.path(dump_dir_files,'imp_gini_loc.txt'),quote=FALSE,append=TRUE)




if(exists('cl')) stopCluster(cl) 
ncores <- 10 # for 10 fold cv
cl <- makeCluster(ncores)
registerDoParallel(cl)
cv1 <- cv.glmnet(x=counts$table_clean[idx_train_dx,],y=as.numeric(meta$DIAGNOSIS[idx_train_dx])-1,
                 family = "binomial",alpha=1,parallel=TRUE)
if(exists('cl')) stopCluster(cl) 

pred1 <- ifelse(predict(cv1,newx=counts$table_clean[-idx_train_dx,],s='lambda.min')>0,1,0)
cf_dx_lasso <- confusionMatrix(pred1,as.numeric(meta$DIAGNOSIS[-idx_train_dx])-1)


cat('Prediction Report',rep('\n',3),sep='',file=file.path(dump_dir_files,'out_lasso.txt'))

cat('Labels: Diagnosis\nModel:\tLasso\nCV:\t70/30 (preBS)\n',sep='',file=file.path(dump_dir_files,'out_lasso.txt'),append=TRUE)
cat('lambda_min:\t',cv1$lambda.min,'\n',sep='',file=file.path(dump_dir_files,'out_lasso.txt'),append=TRUE)
cat(paste0(names(cf_dx_lasso$overall),'\t'),'\n',file=file.path(dump_dir_files,'out_lasso.txt'),append=TRUE)
cat(paste0(round(cf_dx_lasso$overall,5),'\t'),'\n',file=file.path(dump_dir_files,'out_lasso.txt'),append=TRUE)

cat(rep('\n',3),sep='',file=file.path(dump_dir_files,'out_lasso.txt'),append=TRUE)
write.table(cf_dx_lasso$table,file=file.path(dump_dir_files,'out_lasso.txt'),quote=FALSE,append=TRUE)

cat(rep('\n',3),sep='',file=file.path(dump_dir_files,'out_lasso.txt'),append=TRUE)
cat(paste0(names(cf_dx_lasso$byClass),'\t'),'\n',file=file.path(dump_dir_files,'out_lasso.txt'),append=TRUE)
cat(paste0(round(cf_dx_lasso$byClass,5),'\t'),'\n',file=file.path(dump_dir_files,'out_lasso.txt'),append=TRUE)

cat(rep('\n',10),sep='',file=file.path(dump_dir_files,'out_lasso.txt'),append=TRUE)



n_rank <- 25
idx_pos <- order(coef(cv1, s = "lambda.min"),decreasing=TRUE)[1:n_rank]
idx_neg <- order(coef(cv1, s = "lambda.min"),decreasing=FALSE)[1:n_rank]

cat('Associated with CD\n\n',sep='',file=file.path(dump_dir_files,'lasso_kegg_dx.txt'))
write.table(
  data.frame(id1=rownames(coef(cv1, s = "lambda.min"))[idx_pos],
             id2=counts$ids[rownames(coef(cv1, s = "lambda.min"))[idx_pos],'long'],
             coef=coef(cv1, s = "lambda.min")[idx_pos],
             pw=as.vector(sapply(kegg_metadata_xavier[counts$ids[rownames(coef(cv1, s = "lambda.min"))[idx_pos],'long']],
                         function(x) paste0(x,collapse=' | ')))),
  file=file.path(dump_dir_files,'lasso_kegg_dx.txt'),quote=FALSE,append=TRUE)

cat(rep('\n',5),sep='',file=file.path(dump_dir_files,'lasso_otu_dx.txt'),append=TRUE)

cat('Associated with Control\n\n',sep='',file=file.path(dump_dir_files,'lasso_kegg_dx.txt'),append=TRUE)
write.table(
  data.frame(id1=rownames(coef(cv1, s = "lambda.min"))[idx_neg],
             id2=counts$ids[rownames(coef(cv1, s = "lambda.min"))[idx_pos],'long'],
             coef=coef(cv1, s = "lambda.min")[idx_neg],
             pw=as.vector(sapply(kegg_metadata_xavier[counts$ids[rownames(coef(cv1, s = "lambda.min"))[idx_neg],'long']],
                         function(x) paste0(x,collapse=' | ')))),
  file=file.path(dump_dir_files,'lasso_kegg_dx.txt'),quote=FALSE,append=TRUE)






### OTU oversampled model



dump_dir_files <- paste0(dump_dir,'/otu',sep='')
dir.create(dump_dir_files, showWarnings = FALSE,recursive=TRUE)


counts <- collapse_table(otu_lookup,data_mat_xavier,metadata_xavier,
                         taxon=NULL,method='rescaled',
                         kegg=FALSE, cap=TRUE,roof=10^5,
                         filtermin=FALSE,mc=5,ps=.01,
                         filtermax=FALSE,pw=.95,pa=.9)

counts$meta <- metadata_xavier[metadata_xavier$SampleID %in% rownames(counts$table),]
counts$meta_clean <- counts$meta[!is.na(counts$meta$ISOLATION_SOURCE),]
counts$table_clean <- counts$table[counts$meta_clean$SampleID,]
rownames(counts$meta_clean) <- counts$meta_clean$SampleID

set.seed(343)
test_id <- c(sample(unique(counts$meta_clean$SAMPLE_NAME[counts$meta_clean$DIAGNOSIS =='Not IBD']),101),
             sample(unique(counts$meta_clean$SAMPLE_NAME[counts$meta_clean$DIAGNOSIS =='CD']),101))
test_idx <- match(test_id,counts$meta_clean$SAMPLE_NAME)
n_oversamp <- abs(diff(table(counts$meta_clean[-test_idx,]$DIAGNOSIS)))-5
id_oversamp <- sample(counts$meta_clean[-test_idx,]$SampleID[counts$meta_clean[-test_idx,]$DIAGNOSIS=='Not IBD'],
                      n_oversamp,replace=TRUE)

counts$meta_clean <- rbind(counts$meta_clean,counts$meta_clean[id_oversamp,])
counts$table_clean <- rbind(counts$table_clean,counts$table_clean[id_oversamp,])

vocab <- as.character(counts$ids$ids)
docs <- lapply(1:nrow(counts$table_clean),function(i) reshape_doc(counts$table_clean[i,],vocab))
meta <- counts$meta_clean
meta$DIAGNOSIS <- as.factor(ifelse(meta$DIAGNOSIS=='CD',1,0))
meta$ISOLATION_SOURCE <- as.factor(meta$ISOLATION_SOURCE)
meta$COHORT <- as.factor(meta$COHORT)

prep <- prepDocuments(docs,vocab,meta,lower.thresh=0) ### removes terms apparently
counts$table_clean <- counts$table_clean[prep$meta$SampleID,prep$vocab]


# set.seed(43)
# idx_train <- readRDS(file.path(dump_dir,paste0(deparse(substitute(out_rf_dx_kegg)),'.rds')))$idx_train
# out_rf_dx_otu <- randomForest(x=counts$table_clean[idx_train,],y=meta$DIAGNOSIS[idx_train],
#                                xtest=counts$table_clean[-idx_train,],ytest=meta$DIAGNOSIS[-idx_train],
#                                ntree=100,importance=TRUE,proximity=TRUE)
# 
# saveRDS(list(model=out_rf_dx_otu,idx_train=idx_train),
#         file.path(dump_dir,paste0(deparse(substitute(out_rf_dx_otu)),'.rds')))
# 
# set.seed(3455)
# idx_train <- readRDS(file.path(dump_dir,paste0(deparse(substitute(out_rf_loc_kegg)),'.rds')))$idx_train
# out_rf_loc_otu <- randomForest(x=counts$table_clean[idx_train,],y=meta$ISOLATION_SOURCE[idx_train],
#                                 xtest=counts$table_clean[-idx_train,],ytest=meta$ISOLATION_SOURCE[-idx_train],
#                                 ntree=100,importance=TRUE,proximity=TRUE)
# 
# saveRDS(list(model=out_rf_loc_otu,idx_train=idx_train),
#         file.path(dump_dir,paste0(deparse(substitute(out_rf_loc_otu)),'.rds')))


model <- readRDS(file.path(dump_dir,paste0(deparse(substitute(out_rf_dx_otu)),'.rds')))        
out_rf_dx_otu <- model$model
idx_train_dx <- model$idx_train

model <- readRDS(file.path(dump_dir,paste0(deparse(substitute(out_rf_loc_otu)),'.rds')))        
out_rf_loc_otu <- model$model
idx_train_loc <- model$idx_train




cf_dx <- confusionMatrix(out_rf_dx_otu$test$predicted,meta$DIAGNOSIS[-idx_train_dx])
cf_loc <- confusionMatrix(out_rf_loc_otu$test$predicted,meta$ISOLATION_SOURCE[-idx_train_loc])




cat('Prediction Report',rep('\n',3),sep='',file=file.path(dump_dir_files,'out.txt'))

cat('Labels: Diagnosis\nModel:\tRF\nCV:\t70/30 (preBS)\n',sep='',file=file.path(dump_dir_files,'out.txt'),append=TRUE)
cat('m (splits):\t',out_rf_dx_otu$mtry,' (sqrt(p))\n',sep='',file=file.path(dump_dir_files,'out.txt'),append=TRUE)
cat('n (trees):\t',out_rf_dx_otu$ntree,rep('\n',3),sep='',file=file.path(dump_dir_files,'out.txt'),append=TRUE)
cat(paste0(names(cf_dx$overall),'\t'),'\n',file=file.path(dump_dir_files,'out.txt'),append=TRUE)
cat(paste0(round(cf_dx$overall,5),'\t'),'\n',file=file.path(dump_dir_files,'out.txt'),append=TRUE)

cat(rep('\n',3),sep='',file=file.path(dump_dir_files,'out.txt'),append=TRUE)
write.table(cf_dx$table,file=file.path(dump_dir_files,'out.txt'),quote=FALSE,append=TRUE)

cat(rep('\n',3),sep='',file=file.path(dump_dir_files,'out.txt'),append=TRUE)
cat(paste0(names(cf_dx$byClass),'\t'),'\n',file=file.path(dump_dir_files,'out.txt'),append=TRUE)
cat(paste0(round(cf_dx$byClass,5),'\t'),'\n',file=file.path(dump_dir_files,'out.txt'),append=TRUE)

cat(rep('\n',10),sep='',file=file.path(dump_dir_files,'out.txt'),append=TRUE)

cat('Labels: Isolation Source\nModel:\tRF\nCV:\t70/30 (preBS)\n',sep='',file=file.path(dump_dir_files,'out.txt'),append=TRUE)
cat('m (splits):\t',out_rf_loc_otu$mtry,' (sqrt(p))\n',sep='',file=file.path(dump_dir_files,'out.txt'),append=TRUE)
cat('n (trees):\t',out_rf_loc_otu$ntree,rep('\n',3),sep='',file=file.path(dump_dir_files,'out.txt'),append=TRUE)
cat(paste0(names(cf_loc$overall),'\t'),'\n',file=file.path(dump_dir_files,'out.txt'),append=TRUE)
cat(paste0(round(cf_loc$overall,5),'\t'),'\n',file=file.path(dump_dir_files,'out.txt'),append=TRUE)

cat(rep('\n',3),sep='',file=file.path(dump_dir_files,'out.txt'),append=TRUE)
write.table(cf_loc$table,file=file.path(dump_dir_files,'out.txt'),quote=FALSE,append=TRUE)

cat(rep('\n',3),sep='',file=file.path(dump_dir_files,'out.txt'),append=TRUE)
cat('\t',paste0(colnames(cf_loc$byClass),'\t'),'\n',quote=FALSE,file=file.path(dump_dir_files,'out.txt'),append=TRUE)
for (r in 1:nrow(cf_loc$byClass)){
  cat(rownames(cf_loc$byClass)[r],'\t',
      paste0(cf_loc$byClass[r,],'\t'),
      '\n',file=file.path(dump_dir_files,'out.txt'),append=TRUE)
}


pdf(file.path(dump_dir_files,'imp_dx.pdf'),width=12,height=12)
imp_dx <- plot_imp(out_rf_dx_otu)
dev.off()

pdf(file.path(dump_dir_files,'prox_dx.pdf'),width=12,height=12)
plot_prox(out_rf_dx_otu,meta$DIAGNOSIS[idx_train_dx])
dev.off()

pdf(file.path(dump_dir_files,'imp_loc.pdf'),width=12,height=12)
imp_loc <- plot_imp(out_rf_loc_otu)
dev.off()

pdf(file.path(dump_dir_files,'prox_loc.pdf'),width=12,height=12)
plot_prox(out_rf_loc_otu,meta$ISOLATION_SOURCE[idx_train_loc])
dev.off()

imp <- imp_dx$acc[1:500]
write.table(
  data.frame(id1=names(imp[imp > 0]),
             id2=counts$ids[names(imp[imp > 0]),'long'],
             score=imp[imp > 0],
             taxon=apply(otu_lookup[counts$ids[names(imp[imp > 0]),'long'],],1,
                       function(x) paste0(x,collapse='; '))),
  file=file.path(dump_dir_files,'imp_acc_dx.txt'),quote=FALSE,append=TRUE)

imp <- imp_dx$gini[1:500]
write.table(
  data.frame(id1=names(imp[imp > 0]),
             id2=counts$ids[names(imp[imp > 0]),'long'],
             score=imp[imp > 0],
             taxon=apply(otu_lookup[counts$ids[names(imp[imp > 0]),'long'],],1,
                      function(x) paste0(x,collapse='; '))),
  file=file.path(dump_dir_files,'imp_gini_dx.txt'),quote=FALSE,append=TRUE)

imp <- imp_loc$acc[1:500]
write.table(
  data.frame(id1=names(imp[imp > 0]),
             id2=counts$ids[names(imp[imp > 0]),'long'],
             score=imp[imp > 0],
             taxon=apply(otu_lookup[counts$ids[names(imp[imp > 0]),'long'],],1,
                      function(x) paste0(x,collapse='; '))),
  file=file.path(dump_dir_files,'imp_acc_loc.txt'),quote=FALSE,append=TRUE)

imp <- imp_loc$gini[1:500]
write.table(
  data.frame(id1=names(imp[imp > 0]),
             id2=counts$ids[names(imp[imp > 0]),'long'],
             score=imp[imp > 0],
             taxon=apply(otu_lookup[counts$ids[names(imp[imp > 0]),'long'],],1,
                      function(x) paste0(x,collapse='; '))),
  file=file.path(dump_dir_files,'imp_gini_loc.txt'),quote=FALSE,append=TRUE)



if(exists('cl')) stopCluster(cl) 
ncores <- 10 # for 10 fold cv
cl <- makeCluster(ncores)
registerDoParallel(cl)
cv1 <- cv.glmnet(x=counts$table_clean[idx_train_dx,],y=as.numeric(meta$DIAGNOSIS[idx_train_dx])-1,
                 family = "binomial",alpha=1,parallel=TRUE)
if(exists('cl')) stopCluster(cl) 

pred1 <- ifelse(predict(cv1,newx=counts$table_clean[-idx_train_dx,],s='lambda.min')>0,1,0)
cf_dx_lasso <- confusionMatrix(pred1,as.numeric(meta$DIAGNOSIS[-idx_train_dx])-1)


cat('Prediction Report',rep('\n',3),sep='',file=file.path(dump_dir_files,'out_lasso.txt'))

cat('Labels: Diagnosis\nModel:\tLasso\nCV:\t70/30 (preBS)\n',sep='',file=file.path(dump_dir_files,'out_lasso.txt'),append=TRUE)
cat('lambda_min:\t',cv1$lambda.min,'\n',sep='',file=file.path(dump_dir_files,'out_lasso.txt'),append=TRUE)
cat(paste0(names(cf_dx_lasso$overall),'\t'),'\n',file=file.path(dump_dir_files,'out_lasso.txt'),append=TRUE)
cat(paste0(round(cf_dx_lasso$overall,5),'\t'),'\n',file=file.path(dump_dir_files,'out_lasso.txt'),append=TRUE)

cat(rep('\n',3),sep='',file=file.path(dump_dir_files,'out_lasso.txt'),append=TRUE)
write.table(cf_dx_lasso$table,file=file.path(dump_dir_files,'out_lasso.txt'),quote=FALSE,append=TRUE)

cat(rep('\n',3),sep='',file=file.path(dump_dir_files,'out_lasso.txt'),append=TRUE)
cat(paste0(names(cf_dx_lasso$byClass),'\t'),'\n',file=file.path(dump_dir_files,'out_lasso.txt'),append=TRUE)
cat(paste0(round(cf_dx_lasso$byClass,5),'\t'),'\n',file=file.path(dump_dir_files,'out_lasso.txt'),append=TRUE)

cat(rep('\n',10),sep='',file=file.path(dump_dir_files,'out_lasso.txt'),append=TRUE)


n_rank <- 25
idx_pos <- order(coef(cv1, s = "lambda.min"),decreasing=TRUE)[1:n_rank]
idx_neg <- order(coef(cv1, s = "lambda.min"),decreasing=FALSE)[1:n_rank]

cat('Associated with CD\n\n',sep='',file=file.path(dump_dir_files,'lasso_otu_dx.txt'))
write.table(
  data.frame(id1=rownames(coef(cv1, s = "lambda.min"))[idx_pos],
             coef=coef(cv1, s = "lambda.min")[idx_pos],
             taxon=apply(otu_lookup[counts$ids[rownames(coef(cv1, s = "lambda.min"))[idx_pos],'long'],],1,
                         function(x) paste0(x,collapse='; '))),
  file=file.path(dump_dir_files,'lasso_otu_dx.txt'),quote=FALSE,append=TRUE)

cat(rep('\n',5),sep='',file=file.path(dump_dir_files,'lasso_otu_dx.txt'),append=TRUE)

cat('Associated with Control\n\n',sep='',file=file.path(dump_dir_files,'lasso_otu_dx.txt'),append=TRUE)
write.table(
  data.frame(id1=rownames(coef(cv1, s = "lambda.min"))[idx_neg],
             coef=coef(cv1, s = "lambda.min")[idx_neg],
             taxon=apply(otu_lookup[counts$ids[rownames(coef(cv1, s = "lambda.min"))[idx_neg],'long'],],1,
                         function(x) paste0(x,collapse='; '))),
  file=file.path(dump_dir_files,'lasso_otu_dx.txt'),quote=FALSE,append=TRUE)

