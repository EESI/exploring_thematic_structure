#!/usr/bin/env Rscript

library(dada2)
library(ggplot2)
library(gridExtra)

seq_dir <- '~/AG/fastq'
filt_dir <- '~/AG/filtered'

dir.create(filt_dir,showWarnings=FALSE,recursive=TRUE)

seqs <- list.files(seq_dir,full.names=TRUE)
seqs_filt <- gsub('\\.fastq\\.gz$','_filt\\.fastq\\.gz',seqs)
seqs_filt <- gsub('\\/fastq\\/','\\/filtered\\/',seqs_filt)

# idx <- sample(length(seqs),5)
# for (i in idx){
#   print(plotQualityProfile(seqs[i]) + geom_vline(xintercept=135,linetype=3))
# }

for (i in seq_along(seqs)){
  fastqFilter(seqs[i],seqs_filt[i],
              trimLeft=10,truncLen=135,
              maxN=0,maxEE=2,truncQ=2,
              compress=TRUE,
              verbose=TRUE)
}
