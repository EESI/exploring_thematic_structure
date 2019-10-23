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

#
params <- expand.grid(K=c(25,50,75,125,200),
                      content=c(FALSE,TRUE),
                      variable=c('DIAGNOSIS','ISOLATION_SOURCE'))
params <- params %>% 
   filter(!(content == FALSE & variable == 'ISOLATION_SOURCE')) %>%
   arrange(K,variable)

# 
# for (param in 1:NROW(params)){
#    #START LOOP FOR PARAMETERS
#    #
#    
#    
#    
   
   
   
   
   
   source('~/Dropbox/stm_microbiome/code_active/stm_functions.R')
   source('~/Dropbox/stm_microbiome/code_active/nav_froz_fxns_3.R') 
   source('~/Dropbox/stm_microbiome/code_active/performance_1.R')
   source('~/Dropbox/stm_microbiome/code_active/framework.R')
   
   
   load_fit <- FALSE
   save_fit <- TRUE
   save_output <- FALSE
   
   random_seed <- FALSE
   K <- params[param,]$K
   cn_normalize <- TRUE
   content <- params[param,]$content
   seq_sim <- 's97'
   variable <- as.character(params[param,]$variable)
   
   
   
   temp_dump <- tempfile()
   if (save_output) sink(temp_dump)
   prepare_framework_ag(random_seed,K,cn_normalize,content,variable,seq_sim)
   if (save_output) sink()
   
   
   if (load_fit){
      
      load_fits(file.path(save_dir,save_fit_foldername,save_fit_filename))
      
   }else{
      
      if(content){
         fit <- stm::stm(train_docs, vocab, K=K,
                         max.em.its=500, 
                         init.type='Spectral',
                         data=train_meta,content=as.formula(paste0('~',variable)),
                         ngroups=20,
                         seed=seed_fit,verbose=TRUE)
      }else{
         fit <- stm::stm(train_docs, vocab, K=K,
                         max.em.its=500, 
                         init.type='Spectral',
                         ngroups=20, # 10 worked well (6-7s per), trying 25 next
                         seed=seed_fit,verbose=TRUE)
      }
      
      
      train_fit <- ctm_frozen(fit,train_docs,vocab,
                              seed=seed_train,max.em.its=500,emtol=1e-5,avg_iters=1,
                              verbose=TRUE,true_doc_content=TRUE,
                              parallel=TRUE,nc=25)
      if(exists('cl')) stopCluster(cl)
      
      test_fit <- ctm_frozen(fit,test_docs,vocab,
                             seed=seed_train,max.em.its=500,emtol=1e-5,avg_iters=1,
                             verbose=TRUE,true_doc_content=TRUE,
                             parallel=TRUE,nc=25)
      if(exists('cl')) stopCluster(cl)
      
      if (save_fit) {
         dir.create(file.path(save_dir,save_fit_foldername),showWarnings=FALSE)
         if (file.exists(file.path(save_dir,save_fit_foldername,save_fit_filename))){
            dupnum <- sample(1:sqrt(as.integer(Sys.time())),1)
            save_fit_filename <- paste0('dup',dupnum,'_',save_fit_filename)
         }
         saveRDS(list(fit=fit,fit_frozen=list(train_fit,test_fit)),
                 file.path(save_dir,save_fit_foldername,save_fit_filename))
      }
      
   }
   
   train_fit <- fit_frozen[[1]]
   test_fit <- fit_frozen[[2]]
   train_Zbar <- train_fit$Z_bar
   test_Zbar <- test_fit$Z_bar
   
   beta_frozen <- t(exp(train_fit$beta$logbeta[[1]]))
   colnames(beta_frozen) <- paste0('T',1:K)
   rownames(beta_frozen) <- counts$ids[vocab,'long']
   beta_frozen_ra <- beta_frozen
   beta_frozen <- floor(beta_frozen*max(counts$table_clean))
   beta_meta <- data.frame(Topic=colnames(beta_frozen))
   beta_otu <- otu_taxa_xavier[rownames(beta_frozen),]
   
   if (save_fit){
      beta_biom <- make_biom(beta_frozen,beta_meta,beta_otu)
      beta_filename <- 'beta_table.biom'
      if (file.exists(file.path(save_dir,save_fit_foldername,beta_filename))){
         beta_filename <- paste0('dup',dupnum,'_',beta_filename)
      }
      write_biom(beta_biom,file.path(save_dir,save_fit_foldername,beta_filename))
      
      #### permute for null model ####
      set.seed(seed_permuted)
      beta_frozen_permuted <- t(apply(beta_frozen,1,function(x) x[sample(1:K,K,replace=FALSE)]))
      colnames(beta_frozen_permuted) <- colnames(beta_frozen)
      beta_biom <- make_biom(beta_frozen_permuted,beta_meta,beta_otu)
      beta_filename <- 'beta_table_permuted.biom'
      if (file.exists(file.path(save_dir,save_fit_foldername,beta_filename))){
         beta_filename <- paste0('dup',dupnum,'_',beta_filename)
      }
      write_biom(beta_biom,file.path(save_dir,save_fit_foldername,beta_filename))
   }
   
   beta_ppd <- t(floor(apply(train_fit$beta_sums_replicates,c(2,3),mean)))
   colnames(beta_ppd) <- paste0('T',1:K)
   rownames(beta_ppd) <- counts$ids[vocab,'long']
   beta_meta <- data.frame(Topic=colnames(beta_ppd))
   beta_otu <- otu_taxa_xavier[rownames(beta_ppd),]
   
   if (save_fit){
      beta_biom <- make_biom(beta_ppd,beta_meta,beta_otu)
      beta_filename <- 'beta_ppd_table.biom'
      if (file.exists(file.path(save_dir,save_fit_foldername,beta_filename))){
         beta_filename <- paste0('dup',dupnum,'_',beta_filename)
      }
      write_biom(beta_biom,file.path(save_dir,save_fit_foldername,beta_filename))
      
      #### permute for null model ####
      set.seed(seed_permuted)
      beta_ppd_permuted <- t(apply(beta_ppd,1,function(x) x[sample(1:K,K,replace=FALSE)]))
      colnames(beta_ppd_permuted) <- colnames(beta_ppd)
      beta_biom <- make_biom(beta_ppd_permuted,beta_meta,beta_otu)
      beta_filename <- 'beta_ppd_table_permuted.biom'
      if (file.exists(file.path(save_dir,save_fit_foldername,beta_filename))){
         beta_filename <- paste0('dup',dupnum,'_',beta_filename)
      }
      write_biom(beta_biom,file.path(save_dir,save_fit_foldername,beta_filename))
   }
   
   beta_repl_array <- train_fit$beta_sums_replicates
   beta_repl <- aperm(beta_repl_array,c(2,1,3))
   dim(beta_repl) <- c(dim(beta_repl_array)[1]*dim(beta_repl_array)[2],dim(beta_repl_array)[3])
   beta_repl <- t(beta_repl)
   colnames(beta_repl) <- paste0('R',rep(1:dim(beta_repl_array)[1],each=K),
                                 '_',
                                 'T',rep(1:dim(beta_repl_array)[2],dim(beta_repl_array)[1]))
   rownames(beta_repl) <- counts$ids[vocab,'long']
   
   beta_meta <- data_frame(ID=colnames(beta_repl)) %>%
      separate(ID,c('Topic','Replicate'))
   beta_otu <- otu_taxa_xavier[rownames(beta_repl),]
   
   if (save_fit){
      beta_biom <- make_biom(beta_repl,beta_meta,beta_otu)
      beta_filename <- 'beta_repl_table.biom'
      if (file.exists(file.path(save_dir,save_fit_foldername,beta_filename))){
         beta_filename <- paste0('dup',dupnum,'_',beta_filename)
      }
      write_biom(beta_biom,file.path(save_dir,save_fit_foldername,beta_filename))
   }
   
   
   if (save_fit){
      
      out_dx <- eval_binary_prediction(train_Zbar,test_Zbar,train_meta,test_meta,nc=60,covariate='SEX') 
      out_is <- eval_multinomial_prediction(train_Zbar,test_Zbar,train_meta,test_meta,nc=60,covariate='BODY_HABITAT') 
      
      if (file.exists(file.path(save_dir,save_fit_foldername,save_coef_filename))){
         save_coef_filename <- paste0('dup',dupnum,'_',save_coef_filename)
      }
      saveRDS(list(dx=out_dx,is=out_is),file.path(save_dir,save_fit_foldername,save_coef_filename))
      
   }else if(load_fit){
      out <- readRDS(file.path(save_dir,save_fit_foldername,save_coef_filename))
      out_dx <- out$dx
      out_is <- out$is
   }
   
   
   
   
   
   
   
   #### rm(list = ls()[!(ls() %in% c('param','params'))]) ####
   #### } ## end save loop ####
   
   output_dir <- file.path('~/Dropbox/stm_microbiome/output_active/',save_fit_foldername)
   dir.create(output_dir,showWarnings=FALSE)
   
   file.copy(temp_dump,file.path(output_dir,'info.txt'))
   
   if (save_output) pdf(file.path(output_dir,'dx1.pdf'), height=10, width=15)
   grid.table(round(out$dx$score,2))
   if (save_output) dev.off()
   if (save_output) pdf(file.path(output_dir,'dx2.pdf'), height=10, width=15)
   grid.table(round(out$dx$acc,2))
   if (save_output) dev.off()
   if (save_output) pdf(file.path(output_dir,'dx_rf_roc_imp.pdf'), height=10, width=15)
   print(out_dx$figure)
   print(
      data.frame(sapply(out_dx$imp,function(x) x[[1]][[1]]),
                 Topic=paste0('T',1:K)) %>%
         gather(Model,Importance,-Topic) %>%
         mutate(Model=factor(Model,levels=rev(c('RF1500DS','RF750DS','RF350DS','RF350UB')))) %>%
         group_by(Topic) %>%
         mutate(Labely=Importance[Model == 'RF1500DS']) %>%
         ungroup() %>%
         ggplot(aes(x=Model,y=Importance,group=Topic,colour=Topic)) + geom_line(size=2) +
         geom_text(aes(x=4.1,y=Labely,label=Topic)) + 
         ggtitle('Random Forest Importance') +
         theme(legend.position='none',
               axis.text = element_text(size = 20),
               axis.title = element_text(size = 30),
               plot.title = element_text(size = 35))
   )
   if (save_output) dev.off()
   
   if (save_output) pdf(file.path(output_dir,'is1.pdf'), height=10, width=15)
   grid.table(rbind(data.frame(Model='LASSO',round(out$is$score$LASSO$score,2)),
                    data.frame(Model='EN75',round(out$is$score$EN75$score,2)),
                    data.frame(Model='EN5',round(out$is$score$EN5$score,2))))
   if (save_output) dev.off()
   
   if (save_output) pdf(file.path(output_dir,'is2.pdf'), height=10, width=20)
   grid.table(rbind(c(Model='RF1500DS',round(out$is$score$RF1500DS,2)),
                    c(Model='RF750DS',round(out$is$score$RF750DS,2)),
                    c(Model='RF350DS',round(out$is$score$RF350DS,2)),
                    c(Model='RF350UB',round(out$is$score$RF350UB,2))))
   if (save_output) dev.off()
   if (save_output) pdf(file.path(output_dir,'is_rf_imp.pdf'), height=10, width=15)
   print(
      data.frame(sapply(out_is$imp,function(x) x[[1]][[1]]),
                 Topic=paste0('T',1:K)) %>%
         gather(Model,Importance,-Topic) %>%
         mutate(Model=factor(Model,levels=rev(c('RF1500DS','RF750DS','RF350DS','RF350UB')))) %>%
         group_by(Topic) %>%
         mutate(Labely=Importance[Model == 'RF1500DS']) %>%
         ungroup() %>%
         ggplot(aes(x=Model,y=Importance,group=Topic,colour=Topic)) + geom_line(size=2) +
         geom_text(aes(x=4.1,y=Labely,label=Topic)) + 
         ggtitle('Random Forest Importance') +
         theme(legend.position='none',
               axis.text = element_text(size = 20),
               axis.title = element_text(size = 30),
               plot.title = element_text(size = 35))
   )
   if (save_output) dev.off()
   
   
   
   coef_k <- out_dx$lasso
   coef_sig_idx <- which(coef_k != 0)
   coef_sig <- coef_k[coef_sig_idx]
   names(coef_sig) <- coef_sig_idx
   
   
   if (save_output) pdf(file.path(output_dir,'lasso_coef.pdf'), height=8, width=15)
   print(
      qplot(1:K,coef_k,label=1:K,geom='label',fill=as.factor(ifelse(coef_k > 0,'CD','NotIBD')),alpha=abs(coef_k)) + 
         scale_colour_manual(values=c('gray','darkred','darkblue')) +
         theme(legend.position='none') +
         xlab('Topics') + ylab("Coefficient")
   )
   if (save_output) dev.off()
   
   
   if (save_output) pdf(file.path(output_dir,'corr_network.pdf'), height=20, width=20)
   par(mar=c(0,0,0,0))
   plot_corr_network(fit,coef_k,thres=.1,vertex.size=5,vertex.label.cex=.5)
   if (save_output) dev.off()
   
   
   if (save_output) pdf(file.path(output_dir,'taxa_heatmap.pdf'), height=25, width=25)
   plot_taxa_heatmap(beta_frozen_ra,beta_otu,beta_meta,coef_k,taxon=5,dist='jaccard',clust='ward.D2')
   plot_taxa_heatmap(beta_frozen_ra,beta_otu,beta_meta,coef_k,taxon=6,dist='jaccard',clust='ward.D2')
   plot_taxa_heatmap(beta_frozen_ra,beta_otu,beta_meta,coef_k,taxon=7,dist='jaccard',clust='ward.D2')
   if (save_output) dev.off()
   
   
   
   
   taxa_names <- c('Kingdom','Phylum','Class','Order','Family','Genus','Species')
   facet_rows <- 8
   
   taxon <- 2
   #beta_collapse <- log10(floor(collapse_beta(beta_otu,beta_frozen_ra,beta_meta,NULL,taxon,12)*10^6) + 1)
   beta_collapse <- collapse_beta(beta_otu,beta_frozen_ra,beta_meta,NULL,taxon,12)
   
   print(
      data.frame(beta_collapse,Taxon=rownames(beta_collapse)) %>%
         #dplyr::select_(.dots = c('Taxon', paste0('T',names(coef_sig)))) %>%
         gather(Topic,Abundance,-Taxon) %>%
         mutate(Taxon = str_extract(Taxon,"[^_]*$"),
                Topic=reorder(Topic,Abundance,var),
                Correlation=as.factor(ifelse(Topic %in% paste0('T',names(coef_sig[coef_sig > 0])), 'CD+', 'CD-'))) %>%
         ggplot(aes(x=Topic,y=Abundance,fill=Taxon)) +
         geom_bar(colour='black',stat='identity') +
         facet_wrap(~Correlation,nrow=2,scales='free_x') +
         theme(legend.position='bottom',axis.text.x = element_blank()) +
         ggtitle(taxa_names[taxon])
   )
   
   if (K > 48) beta_collapse <- beta_collapse[,topic_select(coef_k,beta_collapse)]
   
   if (save_output) pdf(file.path(output_dir,'taxa_phylum.pdf'), height=18, width=24)
   print(
      data.frame(beta_collapse,Taxon=rownames(beta_collapse)) %>%
         gather(Topic,Abundance,-Taxon) %>%
         mutate(Taxon = str_extract(Taxon,"[^_]*$"),
                Taxon = reorder(Taxon,Abundance,mean),
                Topic=reorder(Topic,Abundance,mean),
                Correlation=as.factor(ifelse(Topic %in% paste0('T',names(coef_sig[coef_sig > 0])), 'CD+', 'CD-'))) %>%
         ggplot(aes(x=Taxon,y=Abundance,fill=Taxon,label=Taxon)) +
         geom_bar(colour='black',position='dodge',stat='identity') +
         geom_text(aes(y=Abundance),angle=90,hjust=-.1,vjust=.5,size=3) +
         facet_wrap(~Topic,nrow=facet_rows) +
         theme(legend.position='none',
               axis.text.x = element_blank()) +
         ggtitle(taxa_names[taxon])
   )
   print(
      data.frame(beta_collapse,Taxon=rownames(beta_collapse)) %>%
         gather(Topic,Abundance,-Taxon) %>%
         mutate(Taxon = str_extract(Taxon,"[^_]*$"),
                Taxon = reorder(Taxon,Abundance,mean),
                Topic=reorder(Topic,Abundance,mean),
                Correlation=as.factor(ifelse(Topic %in% paste0('T',names(coef_sig[coef_sig > 0])), 'CD+', 'CD-'))) %>%
         ggplot(aes(x=Topic,y=Abundance,fill=Topic,label=Topic)) +
         geom_bar(colour='black',position='dodge',stat='identity') +
         geom_text(aes(y=Abundance),angle=90,hjust=-.1,vjust=.5,size=3) +
         facet_wrap(~Taxon,nrow=facet_rows) +
         theme(legend.position='none',
               axis.text.x = element_blank()) +
         ggtitle(taxa_names[taxon])
   )
   
   taxon <- 5
   # beta_collapse <- log10(floor(collapse_beta(beta_otu,beta_frozen_ra,beta_meta,coef_sig,taxon,12)*10^6) + 1)
   beta_collapse <- collapse_beta(beta_otu,beta_frozen_ra,beta_meta,NULL,taxon,24)
   
   print(
      data.frame(beta_collapse,Taxon=rownames(beta_collapse)) %>%
         #dplyr::select_(.dots = c('Taxon', paste0('T',names(coef_sig)))) %>%
         gather(Topic,Abundance,-Taxon) %>%
         mutate(Taxon = str_extract(Taxon,"[^_]*$"),
                Topic=reorder(Topic,Abundance,var),
                Correlation=as.factor(ifelse(Topic %in% paste0('T',names(coef_sig[coef_sig > 0])), 'CD+', 'CD-'))) %>%
         ggplot(aes(x=Topic,y=Abundance,fill=Taxon)) +
         geom_bar(colour='black',stat='identity') +
         facet_wrap(~Correlation,nrow=2,scales='free_x') +
         theme(legend.position='bottom',axis.text.x = element_blank()) +
         ggtitle(taxa_names[taxon])
   )
   
   if (K > 48) beta_collapse <- beta_collapse[,topic_select(coef_k,beta_collapse)]
   
   print(
      data.frame(beta_collapse,Taxon=rownames(beta_collapse)) %>%
         gather(Topic,Abundance,-Taxon) %>%
         mutate(Taxon = str_extract(Taxon,"[^_]*$"),
                Taxon = reorder(Taxon,Abundance,mean),
                Topic=reorder(Topic,Abundance,mean),
                Correlation=as.factor(ifelse(Topic %in% paste0('T',names(coef_sig[coef_sig > 0])), 'CD+', 'CD-'))) %>%
         ggplot(aes(x=Taxon,y=Abundance,fill=Taxon,label=Taxon)) +
         geom_bar(colour='black',position='dodge',stat='identity') +
         geom_text(aes(y=Abundance),angle=90,hjust=-.1,vjust=.5,size=3) +
         facet_wrap(~Topic,nrow=facet_rows) +
         theme(legend.position='none',
               axis.text.x = element_blank()) +
         ggtitle(taxa_names[taxon])
   )
   print(
      data.frame(beta_collapse,Taxon=rownames(beta_collapse)) %>%
         gather(Topic,Abundance,-Taxon) %>%
         mutate(Taxon = str_extract(Taxon,"[^_]*$"),
                Taxon = reorder(Taxon,Abundance,mean),
                Topic=reorder(Topic,Abundance,mean),
                Correlation=as.factor(ifelse(Topic %in% paste0('T',names(coef_sig[coef_sig > 0])), 'CD+', 'CD-'))) %>%
         ggplot(aes(x=Topic,y=Abundance,fill=Topic,label=Topic)) +
         geom_bar(colour='black',position='dodge',stat='identity') +
         geom_text(aes(y=Abundance),angle=90,hjust=-.1,vjust=.5,size=3) +
         facet_wrap(~Taxon,nrow=facet_rows) +
         theme(legend.position='none',
               axis.text.x = element_blank()) +
         ggtitle(taxa_names[taxon])
   )
   if (save_output) dev.off()
   
   
   
   ufw <- read.table(file.path(save_dir,save_fit_foldername,'beta_diversity/weighted_unifrac_beta_table.txt'))
   ufuw <- read.table(file.path(save_dir,save_fit_foldername,'beta_diversity/unweighted_unifrac_beta_table.txt'))
   
   pcoa_ufw <- cmdscale(ufw,k=4)
   colnames(pcoa_ufw) <- paste0('PC',1:NCOL(pcoa_ufw))
   pcoa_ufuw <- cmdscale(ufuw,k=4)
   colnames(pcoa_ufuw) <- paste0('PC',1:NCOL(pcoa_ufuw))
   
   
   if (save_output) pdf(file.path(output_dir,'uf.pdf'), height=12, width=16)
   print(
      
      GGally::ggpairs(data.frame(pcoa_ufw,Lasso=if_else(coef_k > 0, 'CD+',if_else(coef_k < 0, 'CD-','None'))),
                      columns=1:4,
                      legends=TRUE,
                      aes(colour=Lasso))
   )
   
   ##df <- setNames(data.frame(pcoa_ufw), c("x", "y", "z"))
   ####
   ## plotly::plot_ly(df,x=x,y=y,z=z,type='scatter3d',mode='markers',
   ##         size=aa[.bincode(abs(coef_k)/max(abs(coef_k)),a)],
   ##         color=ifelse(coef_k > 0, '1', '2'))
   ####
   
   print(
      
      GGally::ggpairs(data.frame(pcoa_ufuw,Lasso=if_else(coef_k > 0, 'CD+',if_else(coef_k < 0, 'CD-','None'))),
                      columns=1:4,
                      legends=TRUE,
                      aes(colour=Lasso))
   )
   if (save_output) dev.off()
   
   ##df <- setNames(data.frame(pcoa_ufuw), c("x", "y", "z"))
   ####
   ## plotly::plot_ly(df,x=x,y=y,z=z,type='scatter3d',mode='markers',
   ##         size=aa[.bincode(abs(coef_k)/max(abs(coef_k)),a)],
   ##         color=ifelse(coef_k > 0, '1', '2'))
   ####
   
   
   
   
   
   
   
   biom_file <- read_biom(file.path(save_dir,save_fit_foldername,'beta_predicted_metagenome.biom'))
   
   kegg_metadata_risk <- kegg_pw_collapse(biom_file,3)
   kegg_risk <- as.matrix(biom_data(biom_file))
   kegg_risk_ra <- t(t(kegg_risk)/colSums(kegg_risk))
   kegg_risk_ra <- kegg_risk_ra[rowSums(kegg_risk_ra) != 0,]
   
   
   
   facet_rows <- 8
   kegg_collapse <- collapse_beta_kegg(kegg_metadata_risk,kegg_risk,NULL,method='sd',var_min=10^-5,var_len=2,top=24)
   
   print(
      data.frame(kegg_collapse,KO=rownames(kegg_collapse)) %>%
         #dplyr::select_(.dots = c('Taxon', paste0('T',names(coef_sig)))) %>%
         gather(Topic,Abundance,-KO) %>%
         mutate(KO = str_extract(KO,"[^_]*$"),
                Topic=reorder(Topic,Abundance,var),
                Correlation=as.factor(ifelse(Topic %in% paste0('T',names(coef_sig[coef_sig > 0])), 'CD+', 'CD-'))) %>%
         ggplot(aes(x=Topic,y=Abundance,fill=KO)) +
         geom_bar(colour='black',stat='identity') +
         facet_wrap(~Correlation,nrow=2,scales='free_x') +
         theme(legend.position='bottom',axis.text.x = element_blank()) +
         ggtitle('KEGG')
   )
   print(
      data.frame(kegg_collapse,KO=rownames(kegg_collapse)) %>%
         #dplyr::select_(.dots = c('Taxon', paste0('T',names(coef_sig)))) %>%
         gather(Topic,Abundance,-KO) %>%
         mutate(KO = str_extract(KO,"[^_]*$"),
                Correlation=as.factor(ifelse(Topic %in% paste0('T',names(coef_sig[coef_sig > 0])), 'CD+', 'CD-')),
                Topic=reorder(Topic,Abundance,median)) %>%
         group_by(KO,Correlation) %>%
         mutate(Topic=reorder(Topic,Abundance,mean)) %>%
         ungroup() %>%
         ggplot(aes(x=Topic,y=Abundance,fill=Correlation)) +
         geom_bar(colour='black',stat='identity',position='dodge') +
         facet_wrap(~KO,nrow=2,scales='free') +
         theme(legend.position='bottom',axis.text.x = element_blank()) +
         ggtitle('COG')
   )
   
   if (K > 48) kegg_collapse <- kegg_collapse[,topic_select(coef_k,kegg_collapse)]
   
   if (save_output) pdf(file.path(output_dir,'kegg.pdf'), height=18, width=24)
   print(
      data.frame(kegg_collapse,KO=rownames(kegg_collapse)) %>%
         gather(Topic,Abundance,-KO) %>%
         mutate(Correlation=as.factor(ifelse(Topic %in% paste0('T',names(coef_sig[coef_sig > 0])), 'CD+', 'CD-')),
                KO = reorder(KO,Abundance,mean),
                Topic=reorder(Topic,Abundance,mean)) %>%
         ggplot(aes(x=KO,y=Abundance,fill=KO,label=KO)) +
         geom_bar(colour='black',position='dodge',stat='identity') +
         geom_text(aes(y=Abundance),angle=90,hjust=-.1,vjust=.5,size=2.25) +
         facet_wrap(~Topic,nrow=facet_rows,scales='free_y') +
         theme(legend.position='none',
               axis.text.x = element_blank()) +
         ggtitle('KEGG')
   )
   print(
      data.frame(kegg_collapse,KO=rownames(kegg_collapse)) %>%
         gather(Topic,Abundance,-KO) %>%
         mutate(Correlation=as.factor(ifelse(Topic %in% paste0('T',names(coef_sig[coef_sig > 0])), 'CD+', 'CD-')),
                KO = reorder(KO,Abundance,mean),
                Topic=reorder(Topic,Abundance,mean)) %>%
         ggplot(aes(x=Topic,y=Abundance,fill=Topic,label=Topic)) +
         geom_bar(colour='black',position='dodge',stat='identity') +
         geom_text(aes(y=Abundance),angle=90,hjust=-.1,vjust=.5,size=3) +
         facet_wrap(~KO,nrow=facet_rows,scales='free_y') +
         theme(legend.position='none',
               axis.text.x = element_blank()) +
         ggtitle('KEGG')
   )
   
   kegg_collapse <- collapse_beta_kegg(kegg_metadata_risk,kegg_risk,coef_sig,method='sd',var_min=10^-5,var_len=2,top=24)
   
   print(
      data.frame(kegg_collapse,KO=rownames(kegg_collapse)) %>%
         #dplyr::select_(.dots = c('Taxon', paste0('T',names(coef_sig)))) %>%
         gather(Topic,Abundance,-KO) %>%
         mutate(KO = str_extract(KO,"[^_]*$"),
                Topic=reorder(Topic,Abundance,var),
                Correlation=as.factor(ifelse(Topic %in% paste0('T',names(coef_sig[coef_sig > 0])), 'CD+', 'CD-'))) %>%
         ggplot(aes(x=Topic,y=Abundance,fill=KO)) +
         geom_bar(colour='black',stat='identity') +
         facet_wrap(~Correlation,nrow=2,scales='free_x') +
         theme(legend.position='bottom',axis.text.x = element_blank()) +
         ggtitle('KEGG')
   )
   print(
      data.frame(kegg_collapse,KO=rownames(kegg_collapse)) %>%
         #dplyr::select_(.dots = c('Taxon', paste0('T',names(coef_sig)))) %>%
         gather(Topic,Abundance,-KO) %>%
         mutate(KO = str_extract(KO,"[^_]*$"),
                Correlation=as.factor(ifelse(Topic %in% paste0('T',names(coef_sig[coef_sig > 0])), 'CD+', 'CD-')),
                Topic=reorder(Topic,Abundance,median)) %>%
         group_by(KO,Correlation) %>%
         mutate(Topic=reorder(Topic,Abundance,mean)) %>%
         ungroup() %>%
         ggplot(aes(x=Topic,y=Abundance,fill=Correlation)) +
         geom_bar(colour='black',stat='identity',position='dodge') +
         facet_wrap(~KO,nrow=2,scales='free') +
         theme(legend.position='bottom',axis.text.x = element_blank()) +
         ggtitle('COG')
   )
   
   if (K > 48) kegg_collapse <- kegg_collapse[,topic_select(coef_k,kegg_collapse)]
   
   print(
      data.frame(kegg_collapse,KO=rownames(kegg_collapse)) %>%
         gather(Topic,Abundance,-KO) %>%
         mutate(Correlation=as.factor(ifelse(Topic %in% paste0('T',names(coef_sig[coef_sig > 0])), 'CD+', 'CD-')),
                KO = reorder(KO,Abundance,mean),
                Topic=reorder(Topic,Abundance,mean)) %>%
         ggplot(aes(x=KO,y=Abundance,fill=KO,label=KO)) +
         geom_bar(colour='black',position='dodge',stat='identity') +
         geom_text(aes(y=Abundance),angle=90,hjust=-.1,vjust=.5,size=2.25) +
         facet_wrap(~Topic,nrow=facet_rows,scales='free_y') +
         theme(legend.position='none',
               axis.text.x = element_blank()) +
         ggtitle('KEGG')
   )
   print(
      data.frame(kegg_collapse,KO=rownames(kegg_collapse)) %>%
         gather(Topic,Abundance,-KO) %>%
         mutate(Correlation=as.factor(ifelse(Topic %in% paste0('T',names(coef_sig[coef_sig > 0])), 'CD+', 'CD-')),
                KO = reorder(KO,Abundance,mean),
                Topic=reorder(Topic,Abundance,mean)) %>%
         ggplot(aes(x=Topic,y=Abundance,fill=Topic,label=Topic)) +
         geom_bar(colour='black',position='dodge',stat='identity') +
         geom_text(aes(y=Abundance),angle=90,hjust=-.1,vjust=.5,size=3) +
         facet_wrap(~KO,nrow=facet_rows,scales='free_y') +
         theme(legend.position='none',
               axis.text.x = element_blank()) +
         ggtitle('KEGG')
   )
   if (save_output) dev.off()
   
   
   
   if (save_output) pdf(file.path(output_dir,'kegg_heatmap.pdf'), height=25, width=25)
   plot_pw_heatmap(kegg_risk_ra,kegg_metadata_risk,coef_k,dist='jaccard',clust='ward.D2',var_thres=.5)
   if (save_output) dev.off()
   
   pw_list <- pw_counter(kegg_risk_ra,kegg_metadata_risk)
   pw_results <- pw_test(pw_list,kegg_risk_ra,K,fdr=TRUE,exact=FALSE,alpha=.05)
   
   if (save_output) sink(file.path(output_dir,'kegg_analysis.txt'))
   cat('Results from kegg pw analysis: wilcox test, alpha=.05, FDR corrected.\n\n\n')
   print(pw_results)
   if (save_output) sink()
   
   
   
   
   
   biom_file <- read_biom(file.path(save_dir,save_fit_foldername,'beta_predicted_cogs.biom'))
   
   cog_metadata_risk <- kegg_pw_collapse(biom_file,2,cog=TRUE)
   cog_risk <- as.matrix(biom_data(biom_file))
   cog_risk_ra <- t(t(cog_risk)/colSums(cog_risk))
   cog_risk_ra <- cog_risk_ra[rowSums(cog_risk_ra) != 0,]
   
   
   
   facet_rows <- 8
   cog_collapse <- collapse_beta_kegg(cog_metadata_risk,cog_risk,NULL,method='sd',var_min=10^-5,var_len=2,top=24)
   
   print(
      data.frame(cog_collapse,GO=rownames(cog_collapse)) %>%
         #dplyr::select_(.dots = c('Taxon', paste0('T',names(coef_sig)))) %>%
         gather(Topic,Abundance,-GO) %>%
         mutate(GO = str_extract(GO,"[^_]*$"),
                Topic=reorder(Topic,Abundance,var),
                Correlation=as.factor(ifelse(Topic %in% paste0('T',names(coef_sig[coef_sig > 0])), 'CD+', 'CD-'))) %>%
         ggplot(aes(x=Topic,y=Abundance,fill=GO)) +
         geom_bar(colour='black',stat='identity') +
         facet_wrap(~Correlation,nrow=2,scales='free_x') +
         theme(legend.position='bottom',axis.text.x = element_blank()) +
         ggtitle('COG')
   )
   print(
      data.frame(cog_collapse,GO=rownames(cog_collapse)) %>%
         #dplyr::select_(.dots = c('Taxon', paste0('T',names(coef_sig)))) %>%
         gather(Topic,Abundance,-GO) %>%
         mutate(GO = str_extract(GO,"[^_]*$"),
                Correlation=as.factor(ifelse(Topic %in% paste0('T',names(coef_sig[coef_sig > 0])), 'CD+', 'CD-')),
                Topic=reorder(Topic,Abundance,median)) %>%
         group_by(GO,Correlation) %>%
         mutate(Topic=reorder(Topic,Abundance,mean)) %>%
         ungroup() %>%
         ggplot(aes(x=Topic,y=Abundance,fill=Correlation)) +
         geom_bar(colour='black',stat='identity',position='dodge') +
         facet_wrap(~GO,nrow=2,scales='free') +
         theme(legend.position='bottom',axis.text.x = element_blank()) +
         ggtitle('COG')
   )
   
   if (K > 48) cog_collapse <- cog_collapse[,topic_select(coef_k,cog_collapse)]
   
   if (save_output) pdf(file.path(output_dir,'cog.pdf'), height=18, width=24)
   print(
      data.frame(cog_collapse,GO=rownames(cog_collapse)) %>%
         gather(Topic,Abundance,-GO) %>%
         mutate(Correlation=as.factor(ifelse(Topic %in% paste0('T',names(coef_sig[coef_sig > 0])), 'CD+', 'CD-')),
                GO = reorder(GO,Abundance,mean),
                Topic=reorder(Topic,Abundance,mean)) %>%
         ggplot(aes(x=GO,y=Abundance,fill=GO,label=GO)) +
         geom_bar(colour='black',position='dodge',stat='identity') +
         geom_text(aes(y=Abundance),angle=90,hjust=-.1,vjust=.5,size=2.25) +
         facet_wrap(~Topic,nrow=facet_rows,scales='free_y') +
         theme(legend.position='none',
               axis.text.x = element_blank()) +
         ggtitle('COG')
   )
   print(
      data.frame(cog_collapse,GO=rownames(cog_collapse)) %>%
         gather(Topic,Abundance,-GO) %>%
         mutate(Correlation=as.factor(ifelse(Topic %in% paste0('T',names(coef_sig[coef_sig > 0])), 'CD+', 'CD-')),
                GO = reorder(GO,Abundance,mean),
                Topic=reorder(Topic,Abundance,mean)) %>%
         ggplot(aes(x=Topic,y=Abundance,fill=Topic,label=Topic)) +
         geom_bar(colour='black',position='dodge',stat='identity') +
         geom_text(aes(y=Abundance),angle=90,hjust=-.1,vjust=.5,size=3) +
         facet_wrap(~GO,nrow=facet_rows,scales='free_y') +
         theme(legend.position='none',
               axis.text.x = element_blank()) +
         ggtitle('COG')
   )
   
   
   cog_collapse <- collapse_beta_kegg(cog_metadata_risk,cog_risk,coef_sig,method='sd',var_min=10^-5,var_len=2,top=24)
   
   print(
      data.frame(cog_collapse,GO=rownames(cog_collapse)) %>%
         #dplyr::select_(.dots = c('Taxon', paste0('T',names(coef_sig)))) %>%
         gather(Topic,Abundance,-GO) %>%
         mutate(GO = str_extract(GO,"[^_]*$"),
                Topic=reorder(Topic,Abundance,var),
                Correlation=as.factor(ifelse(Topic %in% paste0('T',names(coef_sig[coef_sig > 0])), 'CD+', 'CD-'))) %>%
         ggplot(aes(x=Topic,y=Abundance,fill=GO)) +
         geom_bar(colour='black',stat='identity') +
         facet_wrap(~Correlation,nrow=2,scales='free_x') +
         theme(legend.position='bottom',axis.text.x = element_blank()) +
         ggtitle('COG')
   )
   print(
      data.frame(cog_collapse,GO=rownames(cog_collapse)) %>%
         #dplyr::select_(.dots = c('Taxon', paste0('T',names(coef_sig)))) %>%
         gather(Topic,Abundance,-GO) %>%
         mutate(GO = str_extract(GO,"[^_]*$"),
                Correlation=as.factor(ifelse(Topic %in% paste0('T',names(coef_sig[coef_sig > 0])), 'CD+', 'CD-')),
                Topic=reorder(Topic,Abundance,median)) %>%
         group_by(GO,Correlation) %>%
         mutate(Topic=reorder(Topic,Abundance,mean)) %>%
         ungroup() %>%
         ggplot(aes(x=Topic,y=Abundance,fill=Correlation)) +
         geom_bar(colour='black',stat='identity',position='dodge') +
         facet_wrap(~GO,nrow=2,scales='free') +
         theme(legend.position='bottom',axis.text.x = element_blank()) +
         ggtitle('COG')
   )
   
   if (K > 48) cog_collapse <- cog_collapse[,topic_select(coef_k,cog_collapse)]
   
   print(
      data.frame(cog_collapse,GO=rownames(cog_collapse)) %>%
         gather(Topic,Abundance,-GO) %>%
         mutate(Correlation=as.factor(ifelse(Topic %in% paste0('T',names(coef_sig[coef_sig > 0])), 'CD+', 'CD-')),
                GO = reorder(GO,Abundance,mean),
                Topic=reorder(Topic,Abundance,mean)) %>%
         ggplot(aes(x=GO,y=Abundance,fill=GO,label=GO)) +
         geom_bar(colour='black',position='dodge',stat='identity') +
         geom_text(aes(y=Abundance),angle=90,hjust=-.1,vjust=.5,size=2.25) +
         facet_wrap(~Topic,nrow=facet_rows,scales='free_y') +
         theme(legend.position='none',
               axis.text.x = element_blank()) +
         ggtitle('COG')
   )
   print(
      data.frame(cog_collapse,GO=rownames(cog_collapse)) %>%
         gather(Topic,Abundance,-GO) %>%
         mutate(Correlation=as.factor(ifelse(Topic %in% paste0('T',names(coef_sig[coef_sig > 0])), 'CD+', 'CD-')),
                GO = reorder(GO,Abundance,mean),
                Topic=reorder(Topic,Abundance,mean)) %>%
         ggplot(aes(x=Topic,y=Abundance,fill=Topic,label=Topic)) +
         geom_bar(colour='black',position='dodge',stat='identity') +
         geom_text(aes(y=Abundance),angle=90,hjust=-.1,vjust=.5,size=3) +
         facet_wrap(~GO,nrow=facet_rows,scales='free_y') +
         theme(legend.position='none',
               axis.text.x = element_blank()) +
         ggtitle('COG')
   )
   if (save_output) dev.off()
   
   
   
   if (save_output) pdf(file.path(output_dir,'cog_heatmap.pdf'), height=25, width=25)
   plot_pw_heatmap(cog_risk_ra,cog_metadata_risk,coef_k,dist='jaccard',clust='ward.D2',var_thres=.1)
   if (save_output) dev.off()
   
   cog_list <- pw_counter(cog_risk_ra,cog_metadata_risk)
   cog_results <- pw_test(cog_list,cog_risk_ra,K,fdr=TRUE,exact=FALSE,alpha=.05)
   
   if (save_output)  sink(file.path(output_dir,'cog_analysis.txt'))
   cat('Results from cog pw analysis: wilcox test, alpha=.05, FDR corrected.\n\n\n')
   print(cog_results)
   if (save_output)  sink()
   
   
   
   
#    
#    
#    
#    rm(list = ls()[!(ls() %in% c('param','params'))]) ####
# } ## end save loop ####
