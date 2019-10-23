folder_name <- '~/Dropbox/stm_microbiome/data_active/AG_female/stm_s97_rarefied_10000_unsupervised'

k <- 35
file_name <- sprintf('fits_K_%s.rds',k)
fit <- readRDS(file.path(folder_name,file_name))

cats <- str_split(str_replace(names(fit$performance),'cm_',''),'_')

sink(file.path(folder_name,'performance_table.txt'))
cat('K\tModel\tNormalization\t',paste0(names(fit$performance$cm_stm_qnorm$byClass),collapse='\t'),'\n')
for (i in 3:4) cat('-',cats[[i]],sprintf('%.4f',round(fit$performance[[i]]$byClass,4)),'\n',sep='\t')
for (i in 1:2) cat(k,cats[[i]],sprintf('%.4f',round(fit$performance[[i]]$byClass,4)),'\n',sep='\t')

for (k in c(35,50,75,100,125,150)){
  
  file_name <- sprintf('fits_K_%s.rds',k)
  fit <- readRDS(file.path(folder_name,file_name))
  
  for (i in 1:2) cat(k,cats[[i]],sprintf('%.4f',round(fit$performance[[i]]$byClass,4)),'\n',sep='\t')
  
}
sink()
