#!/usr/bin/env Rscript

lapply(c('ggplot2','stringr','dplyr','purrr','lubridate','readr','tidyr','scales','gridExtra','ggrepel',
         'biom','phyloseq','vegan','RColorBrewer','lme4',
         'glmnet','caret','snow','doSNOW','randomForest',
         'pROC','e1071','fastICA','Rtsne','dada2','Biostrings','DECIPHER','phangorn'),
       require,character.only=TRUE)

data_dir <- '~/Dropbox/stm_microbiome/dada_active/gevers/'

OTU <- readRDS(file.path(data_dir,'seqtab.rds'))

SEQS <- getSequences(OTU)
names(SEQS) <- SEQS 
alignment <- AlignSeqs(DNAStringSet(SEQS),anchor=NA)

phang_align <- phyDat(as.matrix(alignment),type='DNA')
dm <- dist.ml(phang_align)
treeNJ <- NJ(dm) 
tree_fit <- pml(treeNJ,data=phang_align)

fitGTR <- update(tree_fit,k=4,inv=0.2)
fitGTR <- optim.pml(fitGTR,
                    model='GTR',
                    optInv=TRUE,
                    optGamma=TRUE,
                    rearrangement='stochastic',
                    control=pml.control(trace=0))

readr::write_rds(list(tree_fit=tree_fit,fitGTR=fitGTR),
                 file.path(data_dir,'tree_data.rds'),compress='gz')

ape::write.tree(fitGTR$tree,file.path(data_dir,'tree.tree'))