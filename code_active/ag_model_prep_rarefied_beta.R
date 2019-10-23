#!/usr/bin/env Rscript

make_beta_biom <- function(dat,m,dir_name,seed_permuted=534,verbose=TRUE,override=FALSE,permuted=TRUE){
  require(biom)
  
  nfits <- length(dat$fits)
  fit <- dat$fits[[m]]
  K <- fit$settings$dim$K
  Nmin <- dat$filters$rare_min
  
  beta <- t(floor(exp(do.call('rbind',fit$beta$logbeta))*Nmin))
  colnames(beta) <- paste0('T',1:NCOL(beta))
  rownames(beta) <- as.character(dat$counts$ids[dat$vocab,'long'])
  
  beta_meta <- data.frame(Topic=colnames(beta))
  beta_taxa <- dat$taxa
  
  beta_filename1 <- paste0('beta_K_',K,'_nfits_',nfits,'_fit_',m,'.biom')
  beta_biom <- make_biom(beta,beta_meta,beta_taxa)
  
  if (!override) if (file.exists(file.path(dir_name,beta_filename1))) stop('File already exists!')
  write_biom(beta_biom,file.path(dir_name,beta_filename1))
  
  if (permuted){
    set.seed(seed_permuted)
    beta_perm <- t(apply(beta,1,function(x) x[sample(1:NCOL(beta),NCOL(beta),replace=FALSE)]))
    colnames(beta_perm) <- colnames(beta)
    beta_perm <- as.matrix(beta_perm)
    
    if ((NROW(beta_perm) != NROW(beta_taxa)) | (NCOL(beta_perm) != NROW(beta_meta))) stop('Dimensions are not correct.')
    
    beta_filename2 <- paste0('beta_K_',K,'_nfits_',nfits,'_permuted_',seed_permuted,'_fit_',m,'.biom')
    beta_biom <- make_biom(beta_perm,beta_meta,beta_taxa)
    
    if (!override) if (file.exists(file.path(dir_name,beta_filename2))) stop('File already exists!')
    write_biom(beta_biom,file.path(dir_name,beta_filename2))
    
    if (verbose) {cat('\nWrote the following to',beta_filename1,'and',beta_filename2,'\n'); print(fit); print(fit$settings$call)}
  }else{
    if (verbose) {cat('\nWrote the following to',beta_filename1,'\n'); print(fit); print(fit$settings$call)}
  }
}


K <- 35
nfits <- 4
cov <- 'female'
rare_min <- 10000

dir_name <- sprintf('~/Dropbox/stm_microbiome/data_active/AG_%s/stm_s97_rarefied_%s_supervised/model',cov,rare_min)
dir.create(dir_name,showWarnings=FALSE,recursive=TRUE)
fit_filename <- paste0('fits_K_',K,'_nfits_',nfits,'.rds')

dat <- readRDS(file.path(dir_name,fit_filename))

sink(file.path(dir_name,'fit_info.txt'))
for (i in seq_along(dat$fits)) make_beta_biom(dat,i,dir_name,override=TRUE,permuted=FALSE)
sink()

# 1 fit1 prev~sex
# 2 fit2 prev~sex+adr
# 3 fit3 prev~adr | cont~sex
# 4 fit4 prev~sex+adr | cont~sex
# 5 fit4 prev~sex*adr | cont~sex
# 6 fit0 unsup







# K <- 35
# cov <- 'female'
# rare_min <- 10000
# 
# for (model in 1:2)
#   for (seed_split in c(52,550,75,1310,3110)){
#     dir_name <- sprintf('~/Dropbox/stm_microbiome/data_active/AG_%s/stm_s97_rarefied_%s_supervised/features/seed_%s',cov,rare_min,seed_split)
#     dir.create(dir_name,showWarnings=FALSE,recursive=TRUE)
#     fit_filename <- paste0('fits_K_',K,'_model_',model,'.rds')
#     
#     dat <- readRDS(file.path(dir_name,fit_filename))
#     
#     sink(file.path(dir_name,'fit_info.txt'))
#     for (i in seq_along(dat$fits)) make_beta_biom(dat,i,dir_name,override=TRUE,permuted=FALSE)
#     sink()
#   }
# }












# 
# 
# 
# require(biom)
# 
# otu <- t(dat$counts$table_clean)
# rownames(otu) <- as.character(dat$counts$ids[rownames(otu),'long'])
# meta <- dat$counts$meta_clean
# taxa <- dat$taxa
# 
# cov <- 'female'
# rare_min <- 10000
# 
# otu_filename <- 'otu_table.biom'
# otu_biom <- make_biom(otu,meta,taxa)
# 
# dir_name <- sprintf('~/Dropbox/stm_microbiome/data_active/AG_%s/stm_s97_rarefied_%s_supervised/model',cov,rare_min)
# write_biom(otu_biom,file.path(dir_name,otu_filename))