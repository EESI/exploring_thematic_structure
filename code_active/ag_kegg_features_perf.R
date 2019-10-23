kegg_table <- read_biom('~/Dropbox/stm_microbiome/data_active/AG_female/stm_s97_rarefied_10000_supervised/model/otu_table_predicted_metagenome.biom')
kegg <- as.matrix(biom_data(kegg_table))
meta <- sample_metadata(kegg_table)
rm(kegg_table)



name_to_ko <- sapply(kegg_gene_names,identity)
ko1 <- names(name_to_ko)[match(rownames(kegg_sig1),name_to_ko)]
ko2 <- names(name_to_ko)[match(rownames(kegg_sig2),name_to_ko)]

kegg <- t(kegg)
features <- kegg[,unique(c(ko1,ko2))]



meta$FEMALE <- ifelse(meta$SEX == 'female',1,0)
meta$AGE_YEARS <- as.numeric(meta$AGE_YEARS)
meta$AGE <- as.vector(scale(as.numeric(meta$AGE_YEARS)))


set.seed(3453)
idx_train <- unlist(caret::createDataPartition(y=with(meta,as.factor(SEX):as.factor(AGE_CAT)),times=1,p=.8,list=TRUE))
features_train <- features[idx_train,]
kegg_train <- kegg[idx_train,]
meta_train <- meta[idx_train,]
features_test <- features[-idx_train,]
kegg_test <- kegg[-idx_train,]
meta_test <- meta[-idx_train,]



TRAIN <- apply(features_train,2,qnormalize)
TEST <- apply(features_test,2,qnormalize)

LABELS_TRAIN <- meta_train$FEMALE
LABELS_TEST <- meta_test$FEMALE

en <- cv.glmnet(x=TRAIN,y=LABELS_TRAIN,
                family = 'binomial',alpha=.5,parallel=FALSE,
                standardize=TRUE,
                type.measure='auc')


en_pred <- ifelse(predict(en,newx=TEST,s='lambda.min')>0,1,0)

caret::confusionMatrix(en_pred,LABELS_TEST)




TRAIN <- apply(kegg_train,2,qnormalize)
TEST <- apply(kegg_test,2,qnormalize)

LABELS_TRAIN <- meta_train$FEMALE
LABELS_TEST <- meta_test$FEMALE

en <- cv.glmnet(x=TRAIN,y=LABELS_TRAIN,
                family = 'binomial',alpha=.5,parallel=FALSE,
                standardize=TRUE,
                type.measure='auc')


en_pred <- ifelse(predict(en,newx=TEST,s='lambda.min')>0,1,0)

caret::confusionMatrix(en_pred,LABELS_TEST)