prepare_framework <- function(random_seed,K,cn_normalize,content,variable,seq_sim){
   
   if (random_seed){
      ### random seeds
      seed_master <<- sample(1:sqrt(as.integer(Sys.time())),1)
      set.seed(seed_master) # replace with reference seed
      seeds <<- sample(1:999999,7)
      seed_general <<- seeds[1]
      seed_balance <<- seeds[2]
      seed_fit <<- seeds[3]
      seed_test <<- seeds[4]
      seed_train <<- seeds[5]
      seed_score <<- seeds[6]
      seed_permuted <- seeds[7]
   }else{
      ### static seeds
      ### seeds for gevers_out_1
#       seed_general <- 1231  #999
#       seed_balance <- 53453 #542
#       seed_fit <- 654       #534
#       seed_test <- 234      #14
#       seed_train <- 873     #14
#       seed_score <- 5467    #654
#       seed_permuted <- 8973
      ### seeds for gevers_out_2
#       seed_general <<- 1  
#       seed_balance <<- 2  
#       seed_fit <<- 3      
#       seed_test <<- 4     
#       seed_train <<- 5    
#       seed_score <<- 6    
#       seed_permuted <<- 7
      ### seeds for gevers_out_3
         seed_general <<- 11  
         seed_balance <<- 22  
         seed_fit <<- 33      
         seed_test <<- 44     
         seed_train <<- 55    
         seed_score <<- 66    
         seed_permuted <<- 77
   }
   if (random_seed) seed_filename <<- seed_master else seed_filename <- 'static'
   
   save_dir <<- '~/Dropbox/stm_microbiome/data_active/Gevers'
   if (content == TRUE){
      save_fit_foldername <<- paste0('stm_',seq_sim,'_unsupTheta_K_',K,'_cn_norm_',as.integer(cn_normalize),'_content_',as.integer(content),'_',variable,'_seed_',seed_filename)
      save_fit_filename <<- paste0(save_fit_foldername,'.rds')
      save_coef_filename <<- paste0(save_fit_foldername,'_coefs.rds') 
   }else{
      save_fit_foldername <<- paste0('stm_',seq_sim,'_unsupTheta_K_',K,'_cn_norm_',as.integer(cn_normalize),'_content_',as.integer(content),'_seed_',seed_filename)
      save_fit_filename <<- paste0(save_fit_foldername,'.rds')
      save_coef_filename <<- paste0(save_fit_foldername,'_coefs.rds')   
   }
   
   if (cn_normalize){
      biom_file <- read_biom('~/Dropbox/stm_microbiome/Children/Processing/01-raw/simcomp_otus_q32_sim97/sequences_filtered/otu_table_cnNormalized.biom')
      data_mat_xavier <- floor(as.matrix(biom_data(biom_file)))
      biom_file <- read_biom('~/Dropbox/stm_microbiome/Children/Processing/01-raw/simcomp_otus_q32_sim97/sequences_filtered/otu_table.biom')
      otu_taxa_xavier <<- observation_metadata(biom_file)
   }else{
      biom_file <- read_biom('~/Dropbox/stm_microbiome/Children/Processing/01-raw/simcomp_otus_q32_sim97/sequences_filtered/otu_table.biom')
      data_mat_xavier <- floor(as.matrix(biom_data(biom_file)))
      otu_taxa_xavier <<- observation_metadata(biom_file)
   }
   
   metadata_xavier <- read_delim('~/Dropbox/stm_microbiome/Children/Processing/01-raw/metadata.txt',delim='\t')
   sample_info_xavier <- read_csv('~/Dropbox/stm_microbiome/Children/xavier_sample_info.csv')
   
   metadata_xavier <- metadata_xavier %>%
      dplyr::select(SampleID = ends_with('SampleID'), everything()) %>%
      mutate(SAMPLE_NAME = ifelse(STRAIN == "missing", SAMPLE_NAME, STRAIN),
             SAMPLE_NAME = str_replace(SAMPLE_NAME,'-','')) %>%
      dplyr::select(SampleID,SAMPLE_NAME,ISOLATION_SOURCE)
   
   metadata_xavier$COHORT <- ifelse(metadata_xavier$SAMPLE_NAME %in% sample_info_xavier$sample, "RISK","Other")
   metadata_xavier <- metadata_xavier[metadata_xavier$SampleID %in% colnames(data_mat_xavier),] # make sure data_mat_xavier and metadata_xavier have same SRRs
   metadata_xavier <- data.frame(metadata_xavier,
                                 sample_info_xavier[match(metadata_xavier$SAMPLE_NAME,sample_info_xavier$sample),]) # combine metadata_xavier and sample_data
   
   ### there are nearly 700 subjects from RISK (about 1000 samples), which is about 460 CH and 230 None. See this based on subject, not sample name.
   ### the paper says that the RISK cohort was merged with two other cohorts giving about 1750 samples. These are likely the remaining samples not 
   ### accounted for in the csv file (sample_info_xavier) so they should be coded as CD.
   
   ### Keep only no AB patients with samples from rectum, stool, or TI.
   metadata_xavier <- metadata_xavier %>%
      mutate(ISOLATION_SOURCE = ifelse(ISOLATION_SOURCE == "missing", sample_location, ISOLATION_SOURCE),
             ISOLATION_SOURCE = str_to_title(ISOLATION_SOURCE),
             RACE = str_to_title(race),
             STUDY = 'Xavier',
             AGE = 'Minor',
             DIAGNOSIS = ifelse(COHORT == "Other", "CD", Diagnosis)) %>%
      dplyr::select(SampleID, SAMPLE_NAME, ISOLATION_SOURCE, DIAGNOSIS, RACE, STUDY, AGE, COHORT, PCDAI, AB_exposure)
   data_mat_xavier <- data_mat_xavier[,metadata_xavier$SampleID]
   metadata_risk <- metadata_xavier %>% dplyr::filter(COHORT == 'RISK',
                                                      ISOLATION_SOURCE %in% c('Rectum Biopsy',
                                                                              'Stool',
                                                                              'Terminalileum Biopsy'))
   metadata_risk <- metadata_risk %>% dplyr::filter(AB_exposure == 'NoAB')
   
   data_mat_risk <- data_mat_xavier[,metadata_risk$SampleID]
   data_mat_risk <- data_mat_risk[rowSums(data_mat_risk) != 0,]
   
   ### Removing samples in bottom 1% percentile of total sample reads
   data_filt_sub <- filter_subjects(data_mat_risk,metadata_risk,.01)
   data_mat_risk <- data_filt_sub$mat
   metadata_risk <- data_filt_sub$meta
   
   ### Keeping OTUs with 5 or more reads in at least 1% of samples
   counts <<- collapse_table(otu_taxa_xavier,data_mat_risk,metadata_risk,
                            taxon=NULL,method='rescaled',
                            kegg=FALSE, cap=TRUE,roof=10^5,
                            filtermin=TRUE,mc=4,ps=.01,  
                            filtermax=FALSE,pw=.9,pa=.9)
   counts$meta <<- metadata_risk[metadata_risk$SampleID %in% rownames(counts$table),]
   counts$meta_clean <<- counts$meta[!is.na(counts$meta$ISOLATION_SOURCE),]
   counts$table_clean <<- counts$table[counts$meta_clean$SampleID,]
   if (!identical(rownames(counts$table_clean),counts$meta_clean$SampleID)) stop('Please make table sample names match metadata sample IDs!\n')
   
   vocab <<- as.character(counts$ids$ids)
   docs <<- lapply(1:nrow(counts$table_clean),function(i) reshape_doc(counts$table_clean[i,],vocab))
   meta <<- counts$meta_clean
   
   ### Balanced test set, unbalanced training set
   set.seed(seed_balance)
   ### 80/20
   test_N <- floor(table(meta$DIAGNOSIS)*.2)
   
   test_idx <- c(sample(which(names(test_N)[1] == meta$DIAGNOSIS),min(test_N)),
                 sample(which(names(test_N)[2] == meta$DIAGNOSIS),min(test_N)))
   
   test_docs <<- docs[test_idx]
   test_meta <<- meta[test_idx,]
   
   train_docs <<- docs[-test_idx]
   train_meta <<- meta[-test_idx,]
   
   train_v <<- sort(unique(unlist(lapply(train_docs,function(x) x[1,]))))
   test_v <<- sort(unique(unlist(lapply(train_docs,function(x) x[1,]))))
   
   if (!isTRUE(all.equal(train_v,1:length(train_v)))) stop('There are gaps in training set vocabulary!')
   if (!all(test_v %in% train_v)) stop('Some test set words are not in training set vocabulary!')
   
   cat('\n')
   cat('Total docs:',length(docs),'\n',
       'Total vocab:',length(vocab),'\n',
       'Total Meta:',nrow(meta),'x',ncol(meta),'\n',
       'Total DX:',names(table(meta$DIAGNOSIS)),'=',c(table(meta$DIAGNOSIS)),'\n',
       'Total IS:',names(table(meta$ISOLATION_SOURCE)),'=',c(table(meta$ISOLATION_SOURCE)),'\n\n',
       'Training docs:',length(train_docs),'\n',
       'Training vocab:',length(train_v),'\n',
       'Training Meta:',nrow(train_meta),'x',ncol(train_meta),'\n',
       'Training DX:',names(table(train_meta$DIAGNOSIS)),'=',c(table(train_meta$DIAGNOSIS)),'\n',
       'Training IS:',names(table(train_meta$ISOLATION_SOURCE)),'=',c(table(train_meta$ISOLATION_SOURCE)),'\n\n',
       'Testing docs:',length(test_docs),'\n',
       'Testing vocab:',length(test_v),'\n',
       'Testing Meta:',nrow(test_meta),'x',ncol(test_meta),'\n',
       'Testing DX:',names(table(test_meta$DIAGNOSIS)),'=',c(table(test_meta$DIAGNOSIS)),'\n',
       'Testing IS:',names(table(test_meta$ISOLATION_SOURCE)),'=',c(table(test_meta$ISOLATION_SOURCE)),'\n\n'
   )
}

prepare_framework_ag <- function(random_seed,K,cn_normalize,content,variable,seq_sim){
   
   if (random_seed){
      ### random seeds
      seed_master <<- sample(1:sqrt(as.integer(Sys.time())),1)
      set.seed(seed_master) # replace with reference seed
      seeds <<- sample(1:999999,7)
      seed_general <<- seeds[1]
      seed_balance <<- seeds[2]
      seed_fit <<- seeds[3]
      seed_test <<- seeds[4]
      seed_train <<- seeds[5]
      seed_score <<- seeds[6]
      seed_permuted <- seeds[7]
   }else{
      seed_general <<- 1  
      seed_balance <<- 2  
      seed_fit <<- 3      
      seed_test <<- 4     
      seed_train <<- 5    
      seed_score <<- 6    
      seed_permuted <<- 7
   }
   if (random_seed) seed_filename <<- seed_master else seed_filename <- 'static'
   
   save_dir <<- '~/Dropbox/stm_microbiome/data_active/AG'
   if (content == TRUE){
      save_fit_foldername <<- paste0('stm_',seq_sim,'_unsupTheta_K_',K,'_cn_norm_',as.integer(cn_normalize),'_content_',as.integer(content),'_',variable,'_seed_',seed_filename)
      save_fit_filename <<- paste0(save_fit_foldername,'.rds')
      save_coef_filename <<- paste0(save_fit_foldername,'_coefs.rds') 
   }else{
      save_fit_foldername <<- paste0('stm_',seq_sim,'_unsupTheta_K_',K,'_cn_norm_',as.integer(cn_normalize),'_content_',as.integer(content),'_seed_',seed_filename)
      save_fit_filename <<- paste0(save_fit_foldername,'.rds')
      save_coef_filename <<- paste0(save_fit_foldername,'_coefs.rds')   
   }
   
   if (cn_normalize){
#       biom_file <- read_biom('~/Dropbox/stm_microbiome/AG/Processing/01-raw/otus/otu_table_cnNormalized.biom')
#       data_mat_ag <- floor(as.matrix(biom_data(biom_file)))
#       biom_file <- read_biom('~/Dropbox/stm_microbiome/AG/Processing/01-raw/otus/otu_table.biom')
#       otu_taxa_ag <<- observation_metadata(biom_file)
      dat <- readRDS('~/Dropbox/stm_microbiome/AG/Processing/01-raw/otus/ag_data_cnNorm.rds')
      data_mat_ag <- dat$data
      otu_taxa_ag <<- dat$taxa
      rm(dat)
   }else{
      biom_file <- read_biom('~/Dropbox/stm_microbiome/AG/Processing/01-raw/otus/otu_table.biom')
      data_mat_ag <- floor(as.matrix(biom_data(biom_file)))
      otu_taxa_ag <<- observation_metadata(biom_file)
   }
   
   metadata_ag <- read_delim('~/Dropbox/stm_microbiome/AG/Processing/01-raw/metadata.txt',delim='\t') %>%
   dplyr::select(SampleID = ends_with('SampleID'), 
                    LIVER_DISEASE,
                    BODY_PRODUCT,
                    AGE_CAT,
                    AGE_YEARS,
                    SIBO,
                    DOG,
                    ANONYMIZED_NAME,
                    IBD,
                    CAT,
                    IBS,
                    ACNE_MEDICATION_OTC,
                    DIABETES,
                    CDIFF,
                    BODY_SITE,
                    FUNGAL_OVERGROWTH,
                    FLOSSING_FREQUENCY,
                    NON_FOOD_ALLERGIES_SUN,
                    CARDIOVASCULAR_DISEASE,
                    FED_AS_INFANT,
                    ACID_REFLUX,
                    BMI,
                    TEETHBRUSHING_FREQUENCY,
                    DEPRESSION_BIPOLAR_SCHIZOPHRENIA,
                    ROOMMATES,
                    WEIGHT_CHANGE,
                    SEX,
                    ASD,
                    LUNG_DISEASE,
                    ALCOHOL_FREQUENCY,
                    SMOKING_FREQUENCY,
                    KIDNEY_DISEASE,
                    WEIGHT_KG,
                    HEIGHT_CM,
                    PROBIOTIC_FREQUENCY,
                    SLEEP_DURATION,
                    BOWEL_MOVEMENT_QUALITY,
                    MENTAL_ILLNESS_TYPE_SCHIZOPHRENIA,
                    THYROID,
                    NON_FOOD_ALLERGIES_DRUG_EG_PENICILLIN,
                    ACNE_MEDICATION,
                    ANTIBIOTIC_HISTORY,
                    BODY_HABITAT,
                    BMI_CORRECTED,
                    BMI_CAT,
                    DIET_TYPE,
                    AUTOIMMUNE) %>%
   filter(SampleID %in% colnames(data_mat_ag),
          BODY_HABITAT %in% c('UBERON:feces','UBERON:oral cavity','UBERON:skin'),
          SEX %in% c('male','female'))
      
   data_mat_ag <- data_mat_ag[,metadata_ag$SampleID]
   data_mat_ag <- data_mat_ag[rowSums(data_mat_ag) != 0,]
   
   ### Removing samples in bottom 1% percentile of total sample reads
   data_filt_sub <- filter_subjects(data_mat_ag,metadata_ag,.01)
   data_mat_ag <- data_filt_sub$mat
   metadata_ag <- data_filt_sub$meta
   
   ### Keeping OTUs with 5 or more reads in at least .1% of samples.
   ### For Gevers I had 1%, but after some thinking, I changed it.
   ### It should be 1% of an important target group, so in this case,
   ### 1% of skin samples (370) which amounts to:
   filter_min_floor <- min(table(metadata_ag$BODY_HABITAT))*.1/nrow(metadata_ag)

   counts <<- collapse_table(otu_taxa_ag,data_mat_ag,metadata_ag,
                             taxon=NULL,method='rescaled',
                             kegg=FALSE, cap=FALSE,roof=10^5,
                             filtermin=TRUE,mc=4,ps=filter_min_floor,  
                             filtermax=FALSE,pw=.9,pa=.9)
   counts$meta_clean <<- metadata_ag[metadata_ag$SampleID %in% rownames(counts$table),]
   counts$table_clean <<- counts$table[counts$meta_clean$SampleID,]
   if (!identical(rownames(counts$table_clean),counts$meta_clean$SampleID)) stop('Please make table sample names match metadata sample IDs!\n')
   
   vocab <<- as.character(counts$ids$ids)
   docs <<- lapply(1:nrow(counts$table_clean),function(i) reshape_doc(counts$table_clean[i,],vocab))
   meta <<- counts$meta_clean
   
   ### Balanced test set, unbalanced training set
   set.seed(seed_balance)
   ### 80/20
   test_N <- floor(table(meta$BODY_HABITAT)*.2)

   test_idx <- as.vector(sapply(1:length(test_N), function(i) sample(which(names(test_N)[i] == meta$BODY_HABITAT),min(test_N))))
   
   test_docs <<- docs[test_idx]
   test_meta <<- meta[test_idx,]
   
   train_docs <<- docs[-test_idx]
   train_meta <<- meta[-test_idx,]
   
   train_v <<- sort(unique(unlist(lapply(train_docs,function(x) x[1,]))))
   test_v <<- sort(unique(unlist(lapply(train_docs,function(x) x[1,]))))
   
   if (!isTRUE(all.equal(train_v,1:length(train_v)))) stop('There are gaps in training set vocabulary!')
   if (!all(test_v %in% train_v)) stop('Some test set words are not in training set vocabulary!')
   
   cat('\n')
   cat('Total docs:',length(docs),'\n',
       'Total vocab:',length(vocab),'\n',
       'Total Meta:',nrow(meta),'x',ncol(meta),'\n',
       'Total DX:',names(table(meta$BODY_HABITAT)),'=',c(table(meta$BODY_HABITAT)),'\n\n',
       'Training docs:',length(train_docs),'\n',
       'Training vocab:',length(train_v),'\n',
       'Training Meta:',nrow(train_meta),'x',ncol(train_meta),'\n',
       'Training DX:',names(table(train_meta$BODY_HABITAT)),'=',c(table(train_meta$BODY_HABITAT)),'\n\n',
       'Testing docs:',length(test_docs),'\n',
       'Testing vocab:',length(test_v),'\n',
       'Testing Meta:',nrow(test_meta),'x',ncol(test_meta),'\n',
       'Testing DX:',names(table(test_meta$BODY_HABITAT)),'=',c(table(test_meta$BODY_HABITAT)),'\n\n'
   )
}
