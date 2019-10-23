#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 1) stop('Load temp unfinished file (0/1)?', call.=FALSE)

library(XML)
library(readr)

DATA_DIR <- '~/AG'
METADATA_DIR <- file.path(DATA_DIR,'metadata')

MAP <- read_delim(file.path(DATA_DIR,'PRJEB11419.txt'),delim='\t')
ERS <- unique(MAP$secondary_sample_accession)

xmls <- list.files(METADATA_DIR,full.names=TRUE)

if (args[1] == 0){
  md_xml <- xmlParse(xmls[1],isURL=FALSE,encoding='UTF-8',
                     ignoreBlanks=FALSE,trim=TRUE,
                     addAttributeNamespaces=TRUE)
  
  FEATURES <- NULL
  FEATURES <- c(FEATURES,names(xmlToDataFrame(getNodeSet(md_xml,'//IDENTIFIERS'),homogeneous=FALSE)))
  FEATURES <- c(FEATURES,names(unlist(xmlToDataFrame(getNodeSet(md_xml,'//SAMPLE_LINK'),homogeneous=FALSE))))
  FEATURES <- c(FEATURES,names(unlist(xmlToDataFrame(getNodeSet(md_xml,'//DESCRIPTION'),homogeneous=FALSE))))
  FEATURES <- c(FEATURES,names(xmlToDataFrame(getNodeSet(md_xml,'//SAMPLE_NAME'),homogeneous=FALSE)))
  
  FEATURES <- c(FEATURES,
                xmlToDataFrame(getNodeSet(md_xml,'//SAMPLE_ATTRIBUTE'),homogeneous=FALSE)[,1])
  
  METADATA <- matrix('',length(ERS),length(FEATURES),
                     dimnames=list(ERS,FEATURES))
  
  start <- 2
}

if (args[1] == 1){
  TEMP <- readRDS(file.path(DATA_DIR,'metadata_temp.rds'))
  METADATA <- TEMP$METADATA
  
  start <- TEMP$i + 1
  
  cat(sprintf('Starting at i=%s.\n',start))
}


for (i in start:length(xmls)){
  
  cat(sprintf('%s. Extracting metadata from %s.\n',i,xmls[i]))
  
  md_xml <- xmlParse(xmls[i],isURL=FALSE,encoding='UTF-8',
                     ignoreBlanks=FALSE,trim=TRUE,
                     addAttributeNamespaces=TRUE)
  
  samp_info <- unlist(xmlToDataFrame(getNodeSet(md_xml,'//IDENTIFIERS'),homogeneous=FALSE))
  samp_id <- samp_info['PRIMARY_ID']
  METADATA[samp_id,names(samp_info)] <- samp_info
  
  samp_info <- xmlToDataFrame(getNodeSet(md_xml,'//SAMPLE_ATTRIBUTE'),homogeneous=FALSE)
  
  catch <- try(METADATA[samp_id,samp_info[,1]],silent=TRUE)
  
  if (class(catch) == 'try-error'){
    if (length(samp_info)==0){
      cat('*** No metadata: next.\n')
      next
    }
    missing_features <- samp_info[!(samp_info[,1] %in% colnames(METADATA)),1]
    
    for (j in seq_along(missing_features)){
      cat(sprintf('*** Adding %s column.\n',missing_features[j]))
      METADATA <- eval(parse(text=sprintf("cbind(METADATA,%s='')",missing_features[j]))) 
    }
  }
  
  METADATA[samp_id,samp_info[,1]] <- samp_info[,2]
  
  samp_info <- unlist(xmlToDataFrame(getNodeSet(md_xml,'//SAMPLE_LINK'),homogeneous=FALSE))
  METADATA[samp_id,names(samp_info)] <- samp_info
  
  samp_info <- unlist(xmlToDataFrame(getNodeSet(md_xml,'//DESCRIPTION'),homogeneous=FALSE))
  METADATA[samp_id,names(samp_info)] <- samp_info
  
  if (i %% 50 == 0){
    saveRDS(list(i=i,xml=xmls[i],METADATA=METADATA),file.path(DATA_DIR,'metadata_temp.rds'))
  }
}

saveRDS(as.data.frame(METADATA),file.path(DATA_DIR,'metadata.rds'))
file.remove(file.path(DATA_DIR,'metadata_temp.rds'))