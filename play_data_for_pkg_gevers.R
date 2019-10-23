library(readr)
library(tidyr)
library(dplyr)
library(stringr)
library(phyloseq)

# source('~/Dropbox/stm_microbiome/code_active/stm_functions.R')
# source('~/Dropbox/stm_microbiome/code_active/nav_froz_fxns_3.R')
# source('~/Dropbox/stm_microbiome/code_active/performance_1.R')
# source('~/Dropbox/stm_microbiome/code_active/framework.R')

data_dir <- '~/Dropbox/stm_microbiome/qiime_active/gevers/'

BIOM <- import_biom(file.path(data_dir,'picked_otus/otu_table.biom'))
OTU <- t(as(otu_table(BIOM),'matrix'))
TAX <- as.data.frame(as(tax_table(BIOM),'matrix'))
colnames(TAX) <- c('Kingdom','Phylum','Class','Order','Family','Genus','Species')

META <- read_delim(file.path(data_dir,'metadata.txt'),delim='\t') %>%
  mutate(SAMPLE_NAME = ifelse(STRAIN == "missing", SAMPLE_NAME, STRAIN),
         SAMPLE_NAME = str_replace(SAMPLE_NAME,'-','')) %>%
  dplyr::select(SampleID=`#SampleID`,SAMPLE_NAME,ISOLATION_SOURCE) %>%
  left_join(read_csv(file.path(data_dir,'xavier_sample_info.csv')) %>% mutate(SAMPLE_NAME=sample),
            by='SAMPLE_NAME') %>%
  mutate(COHORT=ifelse(SAMPLE_NAME %in% sample,'RISK','Other'),
         ISOLATION_SOURCE=ifelse(ISOLATION_SOURCE %in% 'missing',sample_location,ISOLATION_SOURCE),
         ISOLATION_SOURCE=str_to_title(ISOLATION_SOURCE),
         ISOLATION_SOURCE=ifelse(ISOLATION_SOURCE %in% 'Terminal Ileum','Terminalileum Biopsy',ISOLATION_SOURCE),
         RACE=str_to_title(race),
         STUDY='Xavier',
         AGE='Minor',
         DIAGNOSIS=ifelse(COHORT %in% 'Other','CD',Diagnosis),
         PCDAI = ifelse(is.na(PCDAI) & DIAGNOSIS == 'Not IBD',0,PCDAI)) %>%
  dplyr::select(SampleID, SAMPLE_NAME, ISOLATION_SOURCE, DIAGNOSIS, RACE, STUDY, AGE, COHORT, PCDAI, AB_exposure)


META <- META %>%
  filter(ISOLATION_SOURCE %in% 'Terminalileum Biopsy',
         AB_exposure == 'NoAB')

SAMPLEIDS <- intersect(names(which(rowSums(OTU)>=1000)),META$SampleID)

META <- META %>% filter(SampleID %in% SAMPLEIDS)
OTU <- OTU[SAMPLEIDS,]

OTUIDS <- intersect(rownames(TAX)[!is.na(TAX$Phylum)],names(which(colSums(OTU)>0)))
OTU <- OTU[,OTUIDS]


OTU <- OTU[,names(sort(colSums(OTU),decreasing=TRUE))[1:1000]]
OTU <- OTU[rowSums(OTU) > 0,]

META <- META %>%
  group_by(SAMPLE_NAME) %>%
  filter(n() == 1) %>%
  ungroup() %>%
  select(SampleID,DIAGNOSIS,PCDAI) %>%
  group_by(DIAGNOSIS) %>%
  sample_n(100) %>%
  ungroup() %>%
  as.data.frame()

rownames(META) <- META$SampleID
OTU <- OTU[META$SampleID,]
OTU <- OTU[,colSums(OTU)>0]
TAX <- as.matrix(TAX[colnames(OTU),])


saveRDS(list(OTU=OTU,TAX=TAX,META=META),
        '~/Dropbox/stm_microbiome/playdata_for_pkg_gevers.rds')
