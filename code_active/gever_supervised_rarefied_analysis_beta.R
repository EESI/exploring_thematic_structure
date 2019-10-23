make_beta_biom <- function(dat,m,dir_name,seed_permuted=534){
   fit <- dat$fits[[m]]
   K <- fit$settings$dim$K
   Nmin <- dat$filters$rare_min
      
   beta <- t(floor(exp(do.call('rbind',fit$beta$logbeta))*Nmin))
   colnames(beta) <- paste0('T',1:NCOL(beta))
   rownames(beta) <- as.character(dat$counts$ids[dat$vocab,'long'])
   
   beta_meta <- data.frame(Topic=colnames(beta))
   beta_taxa <- dat$taxa
   
   beta_filename <- paste0('beta_K_',K,'_fit_',m,'.biom')
   beta_biom <- make_biom(beta,beta_meta,beta_taxa)
   write_biom(beta_biom,file.path(dir_name,beta_filename))
   
   set.seed(seed_permuted)
   beta_perm <- t(apply(beta,1,function(x) x[sample(1:NCOL(beta),NCOL(beta),replace=FALSE)]))
   colnames(beta_perm) <- colnames(beta)
   beta_perm <- as.matrix(beta_perm)
   
   if ((NROW(beta_perm) != NROW(beta_taxa)) | (NCOL(beta_perm) != NROW(beta_meta))) stop('Dimensions are not correct.')
   
   beta_filename <- paste0('beta_K_',K,'_permuted_',seed_permuted,'_fit_',m,'.biom')
   beta_biom <- make_biom(beta_perm,beta_meta,beta_taxa)
   write_biom(beta_biom,file.path(dir_name,beta_filename))
}


K <- 50
dir_name <- paste0('~/Dropbox/stm_microbiome/data_active/Gevers_ti/stm_s97_rarefied_supervised')
dir.create(dir_name,showWarnings=FALSE,recursive=TRUE)
fit_filename <- paste0('fits_K_',K,'.rds')

dat <- readRDS(file.path(dir_name,fit_filename))

# make_beta_biom(dat,1,dir_name) # prev ~ DX
# make_beta_biom(dat,2,dir_name) # prev ~ PCDAI
# make_beta_biom(dat,5,dir_name) # prev ~ PCDAI, cont ~ DX

