eval_diagnosis_prediction <- function(train_Zbar,test_Zbar,train_meta,test_meta,nc){
   
   require(snow)
   require(doSNOW)
   
   score <- matrix(NA,7,10,dimnames=list(NULL,c('n','acco','sens','spec','ppv','npv','prev','dr','dp','acci'))) 
   
   train_Zbar <- as.matrix(train_Zbar)
   test_Zbar <- as.matrix(test_Zbar)
   train_data <- data.frame(train_Zbar,
                            labels=as.factor(ifelse(train_meta$DIAGNOSIS == 'CD','CD','NotIBD')))
   test_data <- data.frame(test_Zbar,
                           labels=as.factor(ifelse(test_meta$DIAGNOSIS == 'CD','CD','NotIBD')))
   train_labels_binary <- ifelse(train_data$labels == 'CD',1,0)
   
   cv1 <- cv.glmnet(x=train_Zbar,y=train_labels_binary,
                    family = "binomial",alpha=1,parallel=FALSE,
                    standardize=FALSE,
                    type.measure="auc")
   cat('LASSO complete.\n')
   
   cv2 <- cv.glmnet(x=train_Zbar,y=train_labels_binary,
                    family = "binomial",alpha=.75,parallel=FALSE,
                    standardize=FALSE,
                    type.measure="auc")
   cat('Elastic net (alpha=.75) complete.\n')
   
   cv3 <- cv.glmnet(x=train_Zbar,y=train_labels_binary,
                    family = "binomial",alpha=.5,parallel=FALSE,
                    standardize=FALSE,
                    type.measure="auc")
   cat('Elastic net (alpha=.5) complete.\n')
   
   cv1_coefs <- as.vector(coef(cv1, s="lambda.min"))[-1] # drop intercept
   
   pred_dx_lasso1 <- ifelse(predict(cv1,newx=test_Zbar,s='lambda.min')>0,'CD','NotIBD')
   cf_dx_lasso1 <- caret::confusionMatrix(pred_dx_lasso1,test_data$labels)
   score[1,] <- c(length(test_data$labels),cf_dx_lasso1$overall[1],cf_dx_lasso1$byClass)
   
   pred_dx_lasso2 <- ifelse(predict(cv2,newx=test_Zbar,s='lambda.min')>0,'CD','NotIBD')
   cf_dx_lasso2 <- caret::confusionMatrix(pred_dx_lasso2,test_data$labels)
   score[2,] <- c(length(test_data$labels),cf_dx_lasso2$overall[1],cf_dx_lasso2$byClass)
   
   pred_dx_lasso3 <- ifelse(predict(cv3,newx=test_Zbar,s='lambda.min')>0,'CD','NotIBD')
   cf_dx_lasso3 <- caret::confusionMatrix(pred_dx_lasso3,test_data$labels)
   score[3,] <- c(length(test_data$labels),cf_dx_lasso3$overall[1],cf_dx_lasso3$byClass)
   
   rf_imp <- list()
   ntree <- c(1500,750,350)
   ncores <- nc
   cl <- makeCluster(ncores,type='SOCK')
   registerDoSNOW(cl)
   nmin <- min(table(train_data$labels))
   rf_ctrl <- trainControl(method='cv',classProbs=TRUE,summaryFunction=twoClassSummary,
                           allowParallel=TRUE)
   out_rf <- train(labels ~ .,
                   data=train_data,
                   method='rf',ntree=ntree[1],tuneLength=5,metric='ROC',trControl=rf_ctrl,
                   strata=train_data$labels,
                   sampsize=rep(nmin,2))
   
   pred_dx_rf <- predict(out_rf,test_data,type='prob')[,1]
   roc_dx_rf1 <- pROC::roc(response = test_data$labels, 
                           predictor = pred_dx_rf,
                           levels = rev(levels(test_data$labels)))
   cf_dx_rf1 <- caret::confusionMatrix(ifelse(pred_dx_rf > .5,'CD','NotIBD'),test_data$labels)
   rf_imp[['RF1500DS']] <- varImp(out_rf)
   score[4,] <- c(length(test_data$labels),cf_dx_rf1$overall[1],cf_dx_rf1$byClass)
   cat('Random forest (downsampled, ntree = 1500) complete.\n')
   
   out_rf <- train(labels ~ .,
                   data=train_data,
                   method='rf',ntree=ntree[2],tuneLength=5,metric='ROC',trControl=rf_ctrl,
                   strata=train_data$labels,
                   sampsize=rep(nmin,2))

   pred_dx_rf <- predict(out_rf,test_data,type='prob')[,1]
   roc_dx_rf2 <- pROC::roc(response = test_data$labels, 
                           predictor = pred_dx_rf,
                           levels = rev(levels(test_data$labels)))
   cf_dx_rf2 <- caret::confusionMatrix(ifelse(pred_dx_rf > .5,'CD','NotIBD'),test_data$labels)
   rf_imp[['RF750DS']] <- varImp(out_rf)
   score[5,] <- c(length(test_data$labels),cf_dx_rf2$overall[1],cf_dx_rf2$byClass)
   cat('Random forest (downsampled, ntree = 750) complete.\n')
   
   out_rf <- train(labels ~ .,
                   data=train_data,
                   method='rf',ntree=ntree[3],tuneLength=5,metric='ROC',trControl=rf_ctrl,
                   strata=train_data$labels,
                   sampsize=rep(nmin,2))

   pred_dx_rf <- predict(out_rf,test_data,type='prob')[,1]
   roc_dx_rf3 <- pROC::roc(response = test_data$labels, 
                           predictor = pred_dx_rf,
                           levels = rev(levels(test_data$labels)))
   cf_dx_rf3 <- caret::confusionMatrix(ifelse(pred_dx_rf > .5,'CD','NotIBD'),test_data$labels)
   rf_imp[['RF350DS']] <- varImp(out_rf)
   score[6,] <- c(length(test_data$labels),cf_dx_rf3$overall[1],cf_dx_rf3$byClass)
   cat('Random forest (downsampled, ntree = 350) complete.\n')

   out_rf <- train(labels ~ .,
                   data=train_data,
                   method='rf',ntree=ntree[3],tuneLnegth=5,metric='ROC',trControl=rf_ctrl)

   pred_dx_rf <- predict(out_rf,test_data,type='prob')[,1]
   roc_dx_rf4 <- pROC::roc(response = test_data$labels, 
                           predictor = pred_dx_rf,
                           levels = rev(levels(test_data$labels)))
   cf_dx_rf4 <- caret::confusionMatrix(ifelse(pred_dx_rf > .5,'CD','NotIBD'),test_data$labels)
   rf_imp[['RF350UB']] <- varImp(out_rf)
   score[7,] <- c(length(test_data$labels),cf_dx_rf4$overall[1],cf_dx_rf4$byClass)
   cat('Random forest (unbalanced, ntree = 350) complete.\n')
   
   annotations <- data.frame(Specificity=rep(.35,4),
                             Sensitivity=seq(.45,.1,length=4),
                             ntree=as.factor(c(ntree,ntree[3])),
                             Balancing=c(rep('Downsampled',3),'Unbalanced'),
                             AUC=paste0('ntree: ',c(ntree,ntree[3]),
                                        ' (AUC = ',round(c(roc_dx_rf1$auc[1],roc_dx_rf2$auc[1],roc_dx_rf3$auc[1],roc_dx_rf4$auc[1]),3),
                                        ')'))
   p1 <- rbind(
         data.frame(Sensitivity=roc_dx_rf1$sensitivities,
                    Specificity=roc_dx_rf1$specificities,
                    ntree=ntree[1],Balancing='Downsampled'),
         data.frame(Sensitivity=roc_dx_rf2$sensitivities,
                    Specificity=roc_dx_rf2$specificities,
                    ntree=ntree[2],Balancing='Downsampled'),
         data.frame(Sensitivity=roc_dx_rf3$sensitivities,
                    Specificity=roc_dx_rf3$specificities,
                    ntree=ntree[3],Balancing='Downsampled'),
         data.frame(Sensitivity=roc_dx_rf4$sensitivities,
                    Specificity=roc_dx_rf4$specificities,
                    ntree=ntree[3],Balancing='Unbalanced')) %>%
            mutate(ntree=as.factor(ntree)) %>%
            ggplot(aes(x=Specificity,y=Sensitivity,colour=ntree,linetype=Balancing)) +
               geom_abline(slope=1,intercept=1,colour='gray',size=2) + geom_line(size=1.5) +
               scale_x_reverse(limits=c(1,0)) + 
               ggtitle('RF Performance') +
               geom_text(data=annotations,aes(x=Specificity,y=Sensitivity,colour=ntree,linetype=Balancing,label=AUC)) +
               scale_colour_discrete(guide='none')

   row_names <- c('LASSO','EN75','EN5','RF1500DS','RF750DS','RF350DS','RF350UB')
   rownames(score) <- row_names
   
   acc <- rbind(cf_dx_lasso1$overall,cf_dx_lasso2$overall,cf_dx_lasso3$overall,
                cf_dx_rf1$overall,cf_dx_rf2$overall,cf_dx_rf3$overall,cf_dx_rf4$overall)
   rownames(acc) <- row_names

   print(p1)
   
   if(exists('cl')) stopCluster(cl)
   
   return(list(score=score,acc=acc,lasso=cv1_coefs,figure=p1,imp=rf_imp)) 
}


eval_isolation_source_prediction <- function(train_Zbar,test_Zbar,train_meta,test_meta,nc){
   
   require(snow)
   require(doSNOW)
   
   score <- list()
   
   train_Zbar <- as.matrix(train_Zbar)
   test_Zbar <- as.matrix(test_Zbar)
   train_data <- data.frame(train_Zbar,
                            labels=as.factor(str_replace_all(train_meta$ISOLATION_SOURCE,' ', '')))
   test_data <- data.frame(test_Zbar,
                           labels=as.factor(str_replace_all(test_meta$ISOLATION_SOURCE,' ', '')))
   train_labels_key <- c('Stool','RectumBiopsy', 'TerminalileumBiopsy')
   train_labels_binary <- ifelse(train_data$labels == train_labels_key[1],1,ifelse(train_data$labels == train_labels_key[2],2,3))
   
   cv1 <- cv.glmnet(x=train_Zbar,y=train_labels_binary,
                    family = 'multinomial',alpha=1,parallel=FALSE,
                    standardize=FALSE,
                    type.measure="mae",type.multinomial='ungrouped')
   cat('LASSO complete.\n')
   
   cv2 <- cv.glmnet(x=train_Zbar,y=train_labels_binary,
                    family = 'multinomial',alpha=.75,parallel=FALSE,
                    standardize=FALSE,
                    type.measure="mae",type.multinomial='ungrouped')
   cat('Elastic net (alpha=.75) complete.\n')
   
   cv3 <- cv.glmnet(x=train_Zbar,y=train_labels_binary,
                    family = 'multinomial',alpha=.5,parallel=FALSE,
                    standardize=FALSE,
                    type.measure="mae",type.multinomial='ungrouped')
   cat('Elastic net (alpha=.5) complete.\n')
   
   cv1_coefs <- coef(cv1, s="lambda.min")# drop intercept
   names(cv1_coefs) <- train_labels_key
   
   pred_dx_lasso1 <- train_labels_key[apply(predict(cv1,newx=test_Zbar,s='lambda.min'),1,which.max)]
   cf_dx_lasso1 <- caret::confusionMatrix(factor(pred_dx_lasso1,levels=train_labels_key),factor(test_data$labels,levels=train_labels_key))
   score[[1]] <- list(n=length(test_data$labels),acc=cf_dx_lasso1$overall[1],score=cf_dx_lasso1$byClass)
   
   pred_dx_lasso2 <- train_labels_key[apply(predict(cv2,newx=test_Zbar,s='lambda.min'),1,which.max)]
   cf_dx_lasso2 <- caret::confusionMatrix(factor(pred_dx_lasso2,levels=train_labels_key),factor(test_data$labels,levels=train_labels_key))
   score[[2]] <- list(n=length(test_data$labels),acc=cf_dx_lasso2$overall[1],score=cf_dx_lasso2$byClass)
   
   pred_dx_lasso3 <- train_labels_key[apply(predict(cv3,newx=test_Zbar,s='lambda.min'),1,which.max)]
   cf_dx_lasso3 <- caret::confusionMatrix(factor(pred_dx_lasso3,levels=train_labels_key),factor(test_data$labels,levels=train_labels_key))
   score[[3]] <- list(n=length(test_data$labels),acc=cf_dx_lasso3$overall[1],score=cf_dx_lasso3$byClass)
   
   
   rf_imp <- list()
   ntree <- c(1500,750,350)
   ncores <- nc
   cl <- makeCluster(ncores,type='SOCK')
   registerDoSNOW(cl)
   nmin <- table(train_data$labels)
   rf_ctrl <- trainControl(method='cv',classProbs=TRUE,summaryFunction=multiClassSummary,
                           allowParallel=TRUE)
   out_rf <- train(labels ~ .,
                   data=train_data,
                   method='rf',ntree=ntree[1],tuneLength=5,metric='logLoss',trControl=rf_ctrl,
                   strata=train_data$labels,
                   sampsize=rep(min(nmin),length(nmin)))
   
   pred_dx_rf <- predict(out_rf,test_data,type='prob')
   rf_pred1 <- as.factor(names(pred_dx_rf)[apply(pred_dx_rf,1,which.max)])
   pred_dat <- data.frame(obs=test_data$labels,
                          pred=rf_pred1,
                          pred_dx_rf)
   rf_imp[['RF1500DS']] <- varImp(out_rf)
   score[[4]] <- multiClassSummary(pred_dat,lev=levels(pred_dat$pred))
   cat('Random forest (downsampled, ntree = 1500) complete.\n')
   

   nmin <- table(train_data$labels)
   rf_ctrl <- trainControl(method='cv',classProbs=TRUE,summaryFunction=multiClassSummary,
                           allowParallel=TRUE)
   out_rf <- train(labels ~ .,
                   data=train_data,
                   method='rf',ntree=ntree[2],tuneLength=5,metric='logLoss',trControl=rf_ctrl,
                   strata=train_data$labels,
                   sampsize=rep(min(nmin),length(nmin)))

   pred_dx_rf <- predict(out_rf,test_data,type='prob')
   rf_pred2 <- as.factor(names(pred_dx_rf)[apply(pred_dx_rf,1,which.max)])
   pred_dat <- data.frame(obs=test_data$labels,
                          pred=rf_pred2,
                          pred_dx_rf)
   rf_imp[['RF750DS']] <- varImp(out_rf)
   score[[5]] <- multiClassSummary(pred_dat,lev=levels(pred_dat$pred))
   cat('Random forest (downsampled, ntree = 750) complete.\n')
   

   nmin <- table(train_data$labels)
   rf_ctrl <- trainControl(method='cv',classProbs=TRUE,summaryFunction=multiClassSummary,
                           allowParallel=TRUE)
   out_rf <- train(labels ~ .,
                   data=train_data,
                   method='rf',ntree=ntree[3],tuneLength=5,metric='logLoss',trControl=rf_ctrl,
                   strata=train_data$labels,
                   sampsize=rep(min(nmin),length(nmin)))

   pred_dx_rf <- predict(out_rf,test_data,type='prob')
   rf_pred3 <- as.factor(names(pred_dx_rf)[apply(pred_dx_rf,1,which.max)])
   pred_dat <- data.frame(obs=test_data$labels,
                          pred=rf_pred3,
                          pred_dx_rf)
   rf_imp[['RF350DS']] <- varImp(out_rf)
   score[[6]] <- multiClassSummary(pred_dat,lev=levels(pred_dat$pred))
   cat('Random forest (downsampled, ntree = 350) complete.\n')
   

   nmin <- table(train_data$labels)
   rf_ctrl <- trainControl(method='cv',classProbs=TRUE,summaryFunction=multiClassSummary,
                           allowParallel=TRUE)
   out_rf <- train(labels ~ .,
                   data=train_data,
                   method='rf',ntree=ntree[3],tuneLength=5,metric='logLoss',trControl=rf_ctrl)

   pred_dx_rf <- predict(out_rf,test_data,type='prob')
   rf_pred4 <- as.factor(names(pred_dx_rf)[apply(pred_dx_rf,1,which.max)])
   pred_dat <- data.frame(obs=test_data$labels,
                          pred=rf_pred4,
                          pred_dx_rf)
   rf_imp[['RF350UB']] <- varImp(out_rf)
   score[[7]] <- multiClassSummary(pred_dat,lev=levels(pred_dat$pred))
   cat('Random forest (unbalanced, ntree = 350) complete.\n')
   
   model_names <- c('LASSO','EN75','EN5','RF1500DS','RF750DS','RF350DS','RF350UB')
   names(score) <- model_names
   
   preds <- data.frame(pred_dx_lasso1,pred_dx_lasso2,pred_dx_lasso3,
                       rf_pred1,rf_pred2,rf_pred3,rf_pred4)
   colnames(preds) <- model_names
   preds$true <- test_data$labels
   
   if(exists('cl')) stopCluster(cl)
   
   return(list(score=score,lasso=cv1_coefs,pred=preds,imp=rf_imp)) 
}

eval_binary_prediction <- function(train_Zbar,test_Zbar,train_meta,test_meta,nc,covariate){
   
   score <- matrix(NA,7,10,dimnames=list(NULL,c('n','acco','sens','spec','ppv','npv','prev','dr','dp','acci'))) 
   
   train_Zbar <- as.matrix(train_Zbar)
   test_Zbar <- as.matrix(test_Zbar)
   train_data <- data.frame(train_Zbar,
                            labels=as.factor(ifelse(train_meta[,covariate] == 'CD','CD','NotIBD')))
   test_data <- data.frame(test_Zbar,
                           labels=as.factor(ifelse(test_meta[,covariate] == 'CD','CD','NotIBD')))
   train_labels_binary <- ifelse(train_data$labels == 'CD',1,0)
   
   cv1 <- cv.glmnet(x=train_Zbar,y=train_labels_binary,
                    family = "binomial",alpha=1,parallel=FALSE,
                    standardize=FALSE,
                    type.measure="auc")
   cat('LASSO complete.\n')
   
   cv2 <- cv.glmnet(x=train_Zbar,y=train_labels_binary,
                    family = "binomial",alpha=.75,parallel=FALSE,
                    standardize=FALSE,
                    type.measure="auc")
   cat('Elastic net (alpha=.75) complete.\n')
   
   cv3 <- cv.glmnet(x=train_Zbar,y=train_labels_binary,
                    family = "binomial",alpha=.5,parallel=FALSE,
                    standardize=FALSE,
                    type.measure="auc")
   cat('Elastic net (alpha=.5) complete.\n')
   
   cv1_coefs <- as.vector(coef(cv1, s="lambda.min"))[-1] # drop intercept
   
   pred_dx_lasso1 <- ifelse(predict(cv1,newx=test_Zbar,s='lambda.min')>0,'CD','NotIBD')
   cf_dx_lasso1 <- caret::confusionMatrix(pred_dx_lasso1,test_data$labels)
   score[1,] <- c(length(test_data$labels),cf_dx_lasso1$overall[1],cf_dx_lasso1$byClass)
   
   pred_dx_lasso2 <- ifelse(predict(cv2,newx=test_Zbar,s='lambda.min')>0,'CD','NotIBD')
   cf_dx_lasso2 <- caret::confusionMatrix(pred_dx_lasso2,test_data$labels)
   score[2,] <- c(length(test_data$labels),cf_dx_lasso2$overall[1],cf_dx_lasso2$byClass)
   
   pred_dx_lasso3 <- ifelse(predict(cv3,newx=test_Zbar,s='lambda.min')>0,'CD','NotIBD')
   cf_dx_lasso3 <- caret::confusionMatrix(pred_dx_lasso3,test_data$labels)
   score[3,] <- c(length(test_data$labels),cf_dx_lasso3$overall[1],cf_dx_lasso3$byClass)
   
   rf_imp <- list()
   ntree <- c(1500,750,350)
   ncores <- nc
   cl <- makeCluster(ncores)
   registerDoParallel(cl)
   nmin <- min(table(train_data$labels))
   rf_ctrl <- trainControl(method='cv',classProbs=TRUE,summaryFunction=twoClassSummary,
                           allowParallel=TRUE)
   out_rf <- train(labels ~ .,
                   data=train_data,
                   method='rf',ntree=ntree[1],tuneLength=5,metric='ROC',trControl=rf_ctrl,
                   strata=train_data$labels,
                   sampsize=rep(nmin,2))
   if(exists('cl')) stopCluster(cl)
   pred_dx_rf <- predict(out_rf,test_data,type='prob')[,1]
   roc_dx_rf1 <- pROC::roc(response = test_data$labels, 
                           predictor = pred_dx_rf,
                           levels = rev(levels(test_data$labels)))
   cf_dx_rf1 <- caret::confusionMatrix(ifelse(pred_dx_rf > .5,'CD','NotIBD'),test_data$labels)
   rf_imp[['RF1500DS']] <- varImp(out_rf)
   score[4,] <- c(length(test_data$labels),cf_dx_rf1$overall[1],cf_dx_rf1$byClass)
   cat('Random forest (downsampled, ntree = 1500) complete.\n')
   
   cl <- makeCluster(ncores)
   registerDoParallel(cl)
   out_rf <- train(labels ~ .,
                   data=train_data,
                   method='rf',ntree=ntree[2],tuneLength=5,metric='ROC',trControl=rf_ctrl,
                   strata=train_data$labels,
                   sampsize=rep(nmin,2))
   if(exists('cl')) stopCluster(cl)
   pred_dx_rf <- predict(out_rf,test_data,type='prob')[,1]
   roc_dx_rf2 <- pROC::roc(response = test_data$labels, 
                           predictor = pred_dx_rf,
                           levels = rev(levels(test_data$labels)))
   cf_dx_rf2 <- caret::confusionMatrix(ifelse(pred_dx_rf > .5,'CD','NotIBD'),test_data$labels)
   rf_imp[['RF750DS']] <- varImp(out_rf)
   score[5,] <- c(length(test_data$labels),cf_dx_rf2$overall[1],cf_dx_rf2$byClass)
   cat('Random forest (downsampled, ntree = 750) complete.\n')
   
   cl <- makeCluster(ncores)
   registerDoParallel(cl)
   out_rf <- train(labels ~ .,
                   data=train_data,
                   method='rf',ntree=ntree[3],tuneLength=5,metric='ROC',trControl=rf_ctrl,
                   strata=train_data$labels,
                   sampsize=rep(nmin,2))
   if(exists('cl')) stopCluster(cl)
   pred_dx_rf <- predict(out_rf,test_data,type='prob')[,1]
   roc_dx_rf3 <- pROC::roc(response = test_data$labels, 
                           predictor = pred_dx_rf,
                           levels = rev(levels(test_data$labels)))
   cf_dx_rf3 <- caret::confusionMatrix(ifelse(pred_dx_rf > .5,'CD','NotIBD'),test_data$labels)
   rf_imp[['RF350DS']] <- varImp(out_rf)
   score[6,] <- c(length(test_data$labels),cf_dx_rf3$overall[1],cf_dx_rf3$byClass)
   cat('Random forest (downsampled, ntree = 350) complete.\n')
   
   cl <- makeCluster(ncores)
   registerDoParallel(cl)
   out_rf <- train(labels ~ .,
                   data=train_data,
                   method='rf',ntree=ntree[3],tuneLnegth=5,metric='ROC',trControl=rf_ctrl)
   if(exists('cl')) stopCluster(cl)
   pred_dx_rf <- predict(out_rf,test_data,type='prob')[,1]
   roc_dx_rf4 <- pROC::roc(response = test_data$labels, 
                           predictor = pred_dx_rf,
                           levels = rev(levels(test_data$labels)))
   cf_dx_rf4 <- caret::confusionMatrix(ifelse(pred_dx_rf > .5,'CD','NotIBD'),test_data$labels)
   rf_imp[['RF350UB']] <- varImp(out_rf)
   score[7,] <- c(length(test_data$labels),cf_dx_rf4$overall[1],cf_dx_rf4$byClass)
   cat('Random forest (unbalanced, ntree = 350) complete.\n')
   
   annotations <- data.frame(Specificity=rep(.35,4),
                             Sensitivity=seq(.45,.1,length=4),
                             ntree=as.factor(c(ntree,ntree[3])),
                             Balancing=c(rep('Downsampled',3),'Unbalanced'),
                             AUC=paste0('ntree: ',c(ntree,ntree[3]),
                                        ' (AUC = ',round(c(roc_dx_rf1$auc[1],roc_dx_rf2$auc[1],roc_dx_rf3$auc[1],roc_dx_rf4$auc[1]),3),
                                        ')'))
   p1 <- rbind(
      data.frame(Sensitivity=roc_dx_rf1$sensitivities,
                 Specificity=roc_dx_rf1$specificities,
                 ntree=ntree[1],Balancing='Downsampled'),
      data.frame(Sensitivity=roc_dx_rf2$sensitivities,
                 Specificity=roc_dx_rf2$specificities,
                 ntree=ntree[2],Balancing='Downsampled'),
      data.frame(Sensitivity=roc_dx_rf3$sensitivities,
                 Specificity=roc_dx_rf3$specificities,
                 ntree=ntree[3],Balancing='Downsampled'),
      data.frame(Sensitivity=roc_dx_rf4$sensitivities,
                 Specificity=roc_dx_rf4$specificities,
                 ntree=ntree[3],Balancing='Unbalanced')) %>%
      mutate(ntree=as.factor(ntree)) %>%
      ggplot(aes(x=Specificity,y=Sensitivity,colour=ntree,linetype=Balancing)) +
      geom_abline(slope=1,intercept=1,colour='gray',size=2) + geom_line(size=1.5) +
      scale_x_reverse(limits=c(1,0)) + 
      ggtitle('RF Performance') +
      geom_text(data=annotations,aes(x=Specificity,y=Sensitivity,colour=ntree,linetype=Balancing,label=AUC)) +
      scale_colour_discrete(guide='none')
   
   row_names <- c('LASSO','EN75','EN5','RF1500DS','RF750DS','RF350DS','RF350UB')
   rownames(score) <- row_names
   
   acc <- rbind(cf_dx_lasso1$overall,cf_dx_lasso2$overall,cf_dx_lasso3$overall,
                cf_dx_rf1$overall,cf_dx_rf2$overall,cf_dx_rf3$overall,cf_dx_rf4$overall)
   rownames(acc) <- row_names
   
   print(p1)
   
   return(list(score=score,acc=acc,lasso=cv1_coefs,figure=p1,imp=rf_imp)) 
}


eval_multinomial_prediction <- function(train_Zbar,test_Zbar,train_meta,test_meta,nc,covariate){
   
   score <- list()
   
   train_y <- str_replace_all(train_y,' ','')
   train_y <- str_replace_all(train_y,'UBERON:','')
   
   test_y <- str_replace_all(test_y,' ','')
   test_y <- str_replace_all(test_y,'UBERON:','')
   
   train_Zbar <- as.matrix(train_Zbar)
   test_Zbar <- as.matrix(test_Zbar)
   train_data <- data.frame(train_Zbar,labels=as.factor(train_y))
   test_data <- data.frame(test_Zbar,labels=as.factor(test_y))
   
   train_labels_key <- unique(train_y)
   
   train_labels_binary <- ifelse(train_data$labels == train_labels_key[1],1,ifelse(train_data$labels == train_labels_key[2],2,3))
   
   cv1 <- cv.glmnet(x=train_Zbar,y=train_labels_binary,
                    family = 'multinomial',alpha=1,parallel=FALSE,
                    standardize=FALSE,
                    type.measure="mae",type.multinomial='ungrouped')
   cat('LASSO complete.\n')
   
   cv2 <- cv.glmnet(x=train_Zbar,y=train_labels_binary,
                    family = 'multinomial',alpha=.75,parallel=FALSE,
                    standardize=FALSE,
                    type.measure="mae",type.multinomial='ungrouped')
   cat('Elastic net (alpha=.75) complete.\n')
   
   cv3 <- cv.glmnet(x=train_Zbar,y=train_labels_binary,
                    family = 'multinomial',alpha=.5,parallel=FALSE,
                    standardize=FALSE,
                    type.measure="mae",type.multinomial='ungrouped')
   cat('Elastic net (alpha=.5) complete.\n')
   
   cv1_coefs <- coef(cv1, s="lambda.min")# drop intercept
   names(cv1_coefs) <- train_labels_key
   
   pred_dx_lasso1 <- train_labels_key[apply(predict(cv1,newx=test_Zbar,s='lambda.min'),1,which.max)]
   cf_dx_lasso1 <- caret::confusionMatrix(factor(pred_dx_lasso1,levels=train_labels_key),factor(test_data$labels,levels=train_labels_key))
   score[[1]] <- list(n=length(test_data$labels),acc=cf_dx_lasso1$overall[1],score=cf_dx_lasso1$byClass)
   
   pred_dx_lasso2 <- train_labels_key[apply(predict(cv2,newx=test_Zbar,s='lambda.min'),1,which.max)]
   cf_dx_lasso2 <- caret::confusionMatrix(factor(pred_dx_lasso2,levels=train_labels_key),factor(test_data$labels,levels=train_labels_key))
   score[[2]] <- list(n=length(test_data$labels),acc=cf_dx_lasso2$overall[1],score=cf_dx_lasso2$byClass)
   
   pred_dx_lasso3 <- train_labels_key[apply(predict(cv3,newx=test_Zbar,s='lambda.min'),1,which.max)]
   cf_dx_lasso3 <- caret::confusionMatrix(factor(pred_dx_lasso3,levels=train_labels_key),factor(test_data$labels,levels=train_labels_key))
   score[[3]] <- list(n=length(test_data$labels),acc=cf_dx_lasso3$overall[1],score=cf_dx_lasso3$byClass)
   
   
   rf_imp <- list()
   ntree <- c(1500,750,350)
   ncores <- nc
   cl <- makeCluster(ncores)
   registerDoParallel(cl)
   nmin <- table(train_data$labels)
   rf_ctrl <- trainControl(method='cv',classProbs=TRUE,summaryFunction=multiClassSummary,
                           allowParallel=TRUE)
   out_rf <- train(labels ~ .,
                   data=train_data,
                   method='rf',ntree=ntree[1],tuneLength=5,metric='logLoss',trControl=rf_ctrl,
                   strata=train_data$labels,
                   sampsize=rep(min(nmin),length(nmin)))
   if(exists('cl')) stopCluster(cl)
   pred_dx_rf <- predict(out_rf,test_data,type='prob')
   rf_pred1 <- as.factor(names(pred_dx_rf)[apply(pred_dx_rf,1,which.max)])
   pred_dat <- data.frame(obs=test_data$labels,
                          pred=rf_pred1,
                          pred_dx_rf)
   rf_imp[['RF1500DS']] <- varImp(out_rf)
   score[[4]] <- multiClassSummary(pred_dat,lev=levels(pred_dat$pred))
   cat('Random forest (downsampled, ntree = 1500) complete.\n')
   
   cl <- makeCluster(ncores)
   registerDoParallel(cl)
   nmin <- table(train_data$labels)
   rf_ctrl <- trainControl(method='cv',classProbs=TRUE,summaryFunction=multiClassSummary,
                           allowParallel=TRUE)
   out_rf <- train(labels ~ .,
                   data=train_data,
                   method='rf',ntree=ntree[2],tuneLength=5,metric='logLoss',trControl=rf_ctrl,
                   strata=train_data$labels,
                   sampsize=rep(min(nmin),length(nmin)))
   if(exists('cl')) stopCluster(cl)
   pred_dx_rf <- predict(out_rf,test_data,type='prob')
   rf_pred2 <- as.factor(names(pred_dx_rf)[apply(pred_dx_rf,1,which.max)])
   pred_dat <- data.frame(obs=test_data$labels,
                          pred=rf_pred2,
                          pred_dx_rf)
   rf_imp[['RF750DS']] <- varImp(out_rf)
   score[[5]] <- multiClassSummary(pred_dat,lev=levels(pred_dat$pred))
   cat('Random forest (downsampled, ntree = 750) complete.\n')
   
   cl <- makeCluster(ncores)
   registerDoParallel(cl)
   nmin <- table(train_data$labels)
   rf_ctrl <- trainControl(method='cv',classProbs=TRUE,summaryFunction=multiClassSummary,
                           allowParallel=TRUE)
   out_rf <- train(labels ~ .,
                   data=train_data,
                   method='rf',ntree=ntree[3],tuneLength=5,metric='logLoss',trControl=rf_ctrl,
                   strata=train_data$labels,
                   sampsize=rep(min(nmin),length(nmin)))
   if(exists('cl')) stopCluster(cl)
   pred_dx_rf <- predict(out_rf,test_data,type='prob')
   rf_pred3 <- as.factor(names(pred_dx_rf)[apply(pred_dx_rf,1,which.max)])
   pred_dat <- data.frame(obs=test_data$labels,
                          pred=rf_pred3,
                          pred_dx_rf)
   rf_imp[['RF350DS']] <- varImp(out_rf)
   score[[6]] <- multiClassSummary(pred_dat,lev=levels(pred_dat$pred))
   cat('Random forest (downsampled, ntree = 350) complete.\n')
   
   cl <- makeCluster(ncores)
   registerDoParallel(cl)
   nmin <- table(train_data$labels)
   rf_ctrl <- trainControl(method='cv',classProbs=TRUE,summaryFunction=multiClassSummary,
                           allowParallel=TRUE)
   out_rf <- train(labels ~ .,
                   data=train_data,
                   method='rf',ntree=ntree[3],tuneLength=5,metric='logLoss',trControl=rf_ctrl)
   if(exists('cl')) stopCluster(cl)
   pred_dx_rf <- predict(out_rf,test_data,type='prob')
   rf_pred4 <- as.factor(names(pred_dx_rf)[apply(pred_dx_rf,1,which.max)])
   pred_dat <- data.frame(obs=test_data$labels,
                          pred=rf_pred4,
                          pred_dx_rf)
   rf_imp[['RF350UB']] <- varImp(out_rf)
   score[[7]] <- multiClassSummary(pred_dat,lev=levels(pred_dat$pred))
   cat('Random forest (unbalanced, ntree = 350) complete.\n')
   
   model_names <- c('LASSO','EN75','EN5','RF1500DS','RF750DS','RF350DS','RF350UB')
   names(score) <- model_names
   
   preds <- data.frame(pred_dx_lasso1,pred_dx_lasso2,pred_dx_lasso3,
                       rf_pred1,rf_pred2,rf_pred3,rf_pred4)
   colnames(preds) <- model_names
   preds$true <- test_data$labels
   
   return(list(score=score,lasso=cv1_coefs,pred=preds,imp=rf_imp)) 
}