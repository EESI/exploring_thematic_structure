#!/usr/bin/env Rscript

library(dada2)
library(gridExtra)

data_dir <- '~/Qiime/Children/Processing/01-raw/PRJNA237362'
dada_dir <- '~/Dropbox/stm_microbiome/dada_active/gevers'

dir.create(file.path(dada_dir,'sequences_filtered_10'),showWarnings=FALSE,recursive=TRUE)

srrs <- list.files(data_dir,full.names=TRUE)
fqs <- vector(mode='character')
for (i in seq_along(srrs)){
  fq <- list.files(srrs[i],full.names=TRUE)
  fq <- fq[grepl('fastq.gz',fq)]
  fqs <- c(fqs,fq)
}

which(table(gsub('.*(SRR.*)\\_.*','\\1',fqs)) != 2)

fqs_f <- fqs[grepl('_1.fastq.gz',fqs)]
fqs_r <- fqs[grepl('_2.fastq.gz',fqs)]

idx <- sample(length(fqs_f),5)
for (i in idx){
  pf <- plotQualityProfile(fqs_f[i]) + ggtitle('forward')
  pr <- plotQualityProfile(fqs_r[i]) + ggtitle('reverse')
  grid.arrange(pf,pr,ncol=1)
}

fqs_f_filt <- file.path(dada_dir,'sequences_filtered_10',gsub('.*(SRR.*\\.fastq\\.gz)$','\\1',fqs_f))
fqs_r_filt <- file.path(dada_dir,'sequences_filtered_10',gsub('.*(SRR.*\\.fastq\\.gz)$','\\1',fqs_r))

for (i in seq_along(fqs_f_filt)){
  fastqPairedFilter(c(fqs_f[i],fqs_r[i]),
                    c(fqs_f_filt[i],fqs_r_filt[i]),
                    trimLeft=10,
                    truncLen=c(150,135),
                    maxN=0,maxEE=2,truncQ=2,
                    compress=TRUE)
}
