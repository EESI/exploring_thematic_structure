#include <fstream>
#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
CharacterVector read_fasta(std::string file) {
  CharacterVector records;
  std::ifstream in(file.c_str());
  in.get(); // remove first '>'
  std::string rec;
  while(getline(in,rec,'>')){
    int newLineLoc = rec.find('\n');
    std::string header = rec.substr(0,newLineLoc);
    printf("testing: %s\n",header.c_str());
  }
  return(records);
}

/*** R
# otu_ids <- c('1042097','4129986','4426776','1135553','1052930')

# system.time(out1 <- picrust('/data/sw1/MiscOut/gg_subset.tab',otu_ids))
# system.time(out1 <- picrust('/data/sw1/MiscOut/gg_subset.tab.gz',c('1042097','4129986','4426776','1135553','1052930')))
# system.time(out2 <- picrust('/data/sw1/MiscOut/gg_subset.tab.gz',c('183970','4339831','1848900','1141646','3398712')))
# 
library(biom)
# BIOM <- '/data/sw1/Dropbox/stm_microbiome/qiime_active/gevers/picked_otus/otu_table_cnnorm.biom'
# otu_table <- as.matrix(biom_data(read_biom(BIOM)))
# otu_ids <- rownames(otu_table)

path <- '~/MiscOut/gg_subset.tab'
# path <- '/data/sw1/.local/lib/python2.7/site-packages/picrust/data/ko_13_5_precalculated.tab.gz'
system.time(out2 <- read_fasta(path))
# overlapping_otus <- intersect(out2$matches,otu_ids)
# otu_data <- otu_table[overlapping_otus,]
# rownames(out2$genome_table_out) <- out2$matches
# genome_data <- out2$genome_table_out[overlapping_otus,]
# kegg_data <- round(t(otu_data) %*% genome_data)
# 
# kegg_true <- as.matrix(biom_data(read_biom('~/Dropbox/stm_microbiome/qiime_active/gevers/picked_otus/ko_table.biom')))
*/
