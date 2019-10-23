#!/usr/bin/env Rscript

library(dada2)
library(gridExtra)
library(Biostrings)

data_dir <- '~/Qiime/Children/Processing/01-raw/PRJNA237362'
dada_dir <- '~/Dropbox/stm_microbiome/dada_active/gevers'

fqs_filt <- list.files(file.path(dada_dir,'sequences_filtered'),full.names=TRUE)
fqs_f_filt <- fqs_filt[grepl('_1.fastq.gz',fqs_filt)]
fqs_r_filt <- fqs_filt[grepl('_2.fastq.gz',fqs_filt)]

run_names_f <- sapply(strsplit(basename(fqs_f_filt), "_"), `[`, 1)
run_names_r <- sapply(strsplit(basename(fqs_r_filt), "_"), `[`, 1) 
if(!identical(run_names_f, run_names_r)) stop('Forward and reverse files do not match.')

names(fqs_f_filt) <- run_names_f
names(fqs_r_filt) <- run_names_r

set.seed(100)
run_idx_f <- sample(seq_along(fqs_f_filt))
for (i in seq(10,length(run_idx_f),by=5)){
  
  indexes <- run_idx_f[1:i]
  fqs <- readDNAStringSet(fqs_f_filt[indexes],format='fastq')
  
  cat(sprintf('%s reads over %s samples.\n',length(fqs),i))
  
  if (length(fqs) > 1e6){
    subset_idx_f <- run_idx_f[1:i]
    rm(fqs)
    break
  }
  
}

run_idx_r <- sample(seq_along(fqs_r_filt))
for (i in seq(10,length(run_idx_r),by=5)){
  
  indexes <- run_idx_r[1:i]
  fqs <- readDNAStringSet(fqs_r_filt[indexes],format='fastq')
  
  cat(sprintf('%s reads over %s samples.\n',length(fqs),i))
  
  if (length(fqs) > 1e6){
    subset_idx_r <- run_idx_r[1:i]
    rm(fqs)
    break
  }
  
}

threads <- 60L

derep_err_f <- derepFastq(fqs_f_filt[subset_idx_f],verbose=TRUE)
dd_err_f <- dada(derep_err_f,err=NULL,selfConsist=TRUE,multithread=threads)
err_f <- dd_err_f[[1]]$err_out
rm(list=c('derep_err_f','dd_err_f'))

derep_err_r <- derepFastq(fqs_r_filt[subset_idx_r],verbose=TRUE)
dd_err_r <- dada(derep_err_r,err=NULL,selfConsist=TRUE,multithread=threads)
err_r <- dd_err_r[[1]]$err_out
rm(list=c('derep_err_r','dd_err_r'))

mergers <- vector(mode='list',length(run_names_f))
names(mergers) <- run_names_f
for(run in run_names_f) {
  
  cat('Processing:', run, '\n')
  
  derep_f <- derepFastq(fqs_f_filt[[run]])
  dd_f <- dada(derep_f,err=err_f,multithread=threads)
  
  derep_r <- derepFastq(fqs_r_filt[[run]])
  dd_r <- dada(derep_r,err=err_r,multithread=threads)
  
  merger <- mergePairs(dd_f,derep_f,dd_r,derep_r)
  mergers[[run]] <- merger
  
}
rm(list=c('derep_f','derep_r'))

seqtab <- makeSequenceTable(mergers)

seqtab <- removeBimeraDenovo(seqtab,multithread=threads)

saveRDS(seqtab,file.path(dada_dir,'seqtab.rds'))