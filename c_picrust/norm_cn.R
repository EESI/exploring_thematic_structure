BIOM <- '/data/sw1/Dropbox/stm_microbiome/qiime_active/gevers/picked_otus/otu_table.biom'
otu_table <- as.matrix(biom_data(read_biom(BIOM)))
otus_to_load <- rownames(otu_table)

count_table <- as.matrix(read.table('~/MiscOut/cnntable.tab.gz',sep='\t',header=FALSE,stringsAsFactors=FALSE))
rownames(count_table) <- count_table[,1]
count_table <- count_table[otus_to_load,]

norm_table <- otu_table/count_table[,2]


BIOM = '/data/sw1/Dropbox/stm_microbiome/qiime_active/gevers/picked_otus/otu_table_cnnorm.biom'
true_table <- as.matrix(biom_data(read_biom(BIOM)))



norm_table[1:5,1:5]
true_table[1:5,1:5]
round(cbind(rowSums(norm_table),rowSums(true_table)),2)
