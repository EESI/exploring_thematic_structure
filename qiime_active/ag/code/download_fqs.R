#!/usr/bin/env Rscript

library(readr)

DATA_DIR <- '~/AG'

MAP <- read_delim(file.path(DATA_DIR,'PRJEB11419.txt'),delim='\t')
ERR <- gsub('\\.fastq\\.gz$','',unique(MAP$run_accession))
ERR_DL <- gsub('\\.fastq\\.gz$','',list.files(file.path(DATA_DIR,'fastq')))

MAP_TARG <-  MAP[!(ERR %in% ERR_DL),]

for (i in 1:nrow(MAP_TARG)){
  
  if (is.na(MAP_TARG$fastq_ftp[i])){
    
    if (is.na(MAP_TARG$submitted_ftp[i])){
      cat('*** No url: next.\n')
      next
    }
    
    err <- MAP_TARG$run_accession[i]
    lib <- MAP_TARG$library_name[i]
    
    cat(sprintf('Downloading %s (manually parsing).\n',err))
    
    url <- paste0('ftp://',MAP_TARG$submitted_ftp[i])
    
    FQ <- read_lines(url)
    FQ[seq(1,length(FQ),by=4)] <- paste0('@',err,'.',1:length(FQ),' ',lib,'_',1:length(FQ)-1,'/1')
    gz_temp <- gzfile(file.path(DATA_DIR,'fastq',paste0(err,'.fastq.gz')), 'w')
    writeLines(FQ, gz_temp)
    close(gz_temp)
    
  }else{
    
    err <- MAP_TARG$run_accession[i]
    
    cat(sprintf('Downloading %s.\n',err))
    
    url <- paste0('ftp://',MAP_TARG$fastq_ftp[i])
    download.file(url,file.path(DATA_DIR,'fastq',paste0(err,'.fastq.gz')),quiet=TRUE)
    
  }

}
