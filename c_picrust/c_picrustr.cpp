#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <set>
#include <Rcpp.h>

#include "gzstream.h"

using namespace std;
using namespace Rcpp;

vector<string> split_line(const string &source, const char *delimiter = " ") {

  vector<string> results;

  size_t prev = 0;
  size_t next = 0;

  while ((next = source.find_first_of(delimiter,prev)) != string::npos) {

    if (next - prev != 0) {
      results.push_back(source.substr(prev,next - prev));
    }

    prev = next + 1;

  }

  if (prev < source.size()) {
    results.push_back(source.substr(prev));
  }

  return results;

}

List sweep_picrust(string file_path, StringVector otu_id_targets) {

  int a=0;
  int o=0;
  List out;
  string rec;
  StringVector lines;
  StringVector matches;
  StringVector genemeta;

  ifstream file(file_path.c_str()); // load in file path

  getline(file,rec,'\n'); // start with line 1, column names
  
  int id_pos = rec.find('\t');
  int id_pos_end = rec.find("metadata");
  int pimeta_pos_end = rec.find("\n");

  // extract from #OTU to picrust metadata_*
  string header_left = rec.substr(id_pos + 1,id_pos_end - 10);
  vector<string> gene_ids = split_line(header_left,"\t");
  
  // extract from picrust metadata_* to end
  string header_right = rec.substr(id_pos_end,pimeta_pos_end);
  vector<string> pimeta_ids = split_line(header_right,"\t");

  int gene_n=gene_ids.size();
  int pimeta_n=pimeta_ids.size();

  IntegerMatrix genome_table(otu_id_targets.size(),gene_n);
  NumericMatrix pimeta_table(otu_id_targets.size(),pimeta_n);

  while(getline(file,rec,'\n')) { // starting at line 2, data row 1

    if (a % 5000 == 0){
      Rcpp::checkUserInterrupt();
    }
    a += 1;
    
    // check if at gene metadata lines (final 2)
    if (rec.substr(0,8) == "metadata"){

      int line_end = rec.find('\n');
      genemeta.push_back(rec.substr(id_pos + 1,line_end));

    // otherwise, extract counts
    }else{

      // get row name
      int id_pos = rec.find('\t');
      string otu_id = rec.substr(0,id_pos);

      // look over all otu ids from otu table to check
      // if row name match otu id from otu table
      for (int j=0;j<otu_id_targets.size();j++){

        if (otu_id_targets[j] == otu_id){

          int line_end = rec.find('\n');
          string counts_start = rec.substr(id_pos + 1,line_end);
          vector<string> counts = split_line(counts_start,"\t");

          for (int k=0;k<gene_n;k++){
            genome_table(o,k) = atoi(counts[k].c_str()); // str to int
          }

          int c=0;
          for (int m=gene_n;m<counts.size();m++){
            pimeta_table(o,c) = atof(counts[m].c_str()); // str to dbl
            c += c;
          }

          o += 1;

          // for (int s=0; s<otu_id_targets.size();s++){
          //   Rcout << otu_id_targets(s) << " ";
          // }

          matches.push_back(otu_id); // record otu id of match
          otu_id_targets.erase(j); // remove otu id from targets

          // printf("\niter %d, found a match: %s\n",o,otu_id.c_str());

          break;
        }

      }

    }

  }
  file.close();

  if (o == 0){

    return(R_NilValue);

  }else{

    SubMatrix<INTSXP> genome_table_out = genome_table(Range(0,o - 1),_);
    SubMatrix<REALSXP> pimeta_table_out = pimeta_table(Range(0,o - 1),_);

    out["gene_ids"] = gene_ids;
    out["matches"] = matches;
    out["genemeta"] = genemeta;
    out["genome_table_out"] = genome_table_out;
    out["pimeta_table_out"] = pimeta_table_out;

    return(out);

  }

}

List sweep_picrust_gz(string file_path, StringVector otu_id_targets) {
  
  int a=0;
  int o=0;
  List out;
  string rec;
  StringVector lines;
  StringVector matches;
  StringVector genemeta;

  igzstream file;
  file.open(file_path.c_str());

  getline(file,rec,'\n'); // start with line 1, column names
  
  int id_pos = rec.find('\t');
  int id_pos_end = rec.find("metadata");
  int pimeta_pos_end = rec.find("\n");

  // extract from #OTU to picrust metadata_*
  string header_left = rec.substr(id_pos + 1,id_pos_end - 10);
  vector<string> gene_ids = split_line(header_left,"\t");
  
  // extract from picrust metadata_* to end
  string header_right = rec.substr(id_pos_end,pimeta_pos_end);

  vector<string> pimeta_ids = split_line(header_right,"\t");
  
  int gene_n=gene_ids.size();
  int pimeta_n=pimeta_ids.size();

  IntegerMatrix genome_table(otu_id_targets.size(),gene_n);
  NumericMatrix pimeta_table(otu_id_targets.size(),pimeta_n);
  
  while(getline(file,rec,'\n')) { // starting at line 2, data row 1
    
    if (a % 5000 == 0){
      Rcpp::checkUserInterrupt();
    }
    a += 1;
    
    // check if at gene metadata lines (final 2)
    if (rec.substr(0,8) == "metadata"){

      // printf("Testing!\n\n");

      int line_end = rec.find('\n');
      genemeta.push_back(rec.substr(id_pos + 1,line_end));

    // otherwise, extract counts
    }else{

      // get row name
      int id_pos = rec.find('\t');
      string otu_id = rec.substr(0,id_pos);

      // look over all otu ids from otu table to check
      // if row name match otu id from otu table
      for (int j=0;j<otu_id_targets.size();j++){

        if (otu_id_targets[j] == otu_id){

          int line_end = rec.find('\n');
          string counts_start = rec.substr(id_pos + 1,line_end);
          vector<string> counts = split_line(counts_start,"\t");

          for (int k=0;k<gene_n;k++){
            genome_table(o,k) = atoi(counts[k].c_str()); // str to int
          }

          int c=0;
          for (int m=gene_n;m<counts.size();m++){
            pimeta_table(o,c) = atof(counts[m].c_str()); // str to dbl
            c += c;
          }

          o += 1;

          // for (int s=0; s<otu_id_targets.size();s++){
          //   Rcout << otu_id_targets(s) << " ";
          // }

          matches.push_back(otu_id); // record otu id of match
          otu_id_targets.erase(j); // remove otu id from targets

          // printf("\niter %d, found a match: %s\n",o,otu_id.c_str());

          break;
        }

      }

    }

  }
  file.close();

  if (o == 0){

    return(R_NilValue);

  }else{

    SubMatrix<INTSXP> genome_table_out = genome_table(Range(0,o - 1),_);
    SubMatrix<REALSXP> pimeta_table_out = pimeta_table(Range(0,o - 1),_);

    out["gene_ids"] = gene_ids;
    out["matches"] = matches;
    out["genemeta"] = genemeta;
    out["genome_table_out"] = genome_table_out;
    out["pimeta_table_out"] = pimeta_table_out;

    return(out);

  }

}

// [[Rcpp::export]]

List picrust(string file_path, StringVector otu_id_targets) {

  string file_ext = file_path.substr(file_path.rfind('.'),file_path.size());

  if (file_ext == ".gz"){

    return(sweep_picrust_gz(file_path,otu_id_targets));

  }else{

    return(sweep_picrust(file_path,otu_id_targets));

  }

}


/*** R
# otu_ids <- c('1042097','4129986','4426776','1135553','1052930')

# system.time(out1 <- picrust('/data/sw1/MiscOut/gg_subset.tab',otu_ids))
# system.time(out1 <- picrust('/data/sw1/MiscOut/gg_subset.tab.gz',c('1042097','4129986','4426776','1135553','1052930')))
# system.time(out2 <- picrust('/data/sw1/MiscOut/gg_subset.tab.gz',c('183970','4339831','1848900','1141646','3398712')))



# library(biom)
# BIOM <- '/data/sw1/Dropbox/stm_microbiome/qiime_active/gevers/picked_otus/otu_table_cnnorm.biom'
# otu_table <- as.matrix(biom_data(read_biom(BIOM)))
# otu_ids <- rownames(otu_table)
# 
# path <- '~/MiscOut/gg_subset.tab'
path <- '/data/sw1/.local/lib/python2.7/site-packages/picrust/data/ko_13_5_precalculated.tab.gz'
system.time(out2 <- picrust(path,otu_ids))
# overlapping_otus <- intersect(out2$matches,otu_ids)
# otu_data <- otu_table[overlapping_otus,]
# rownames(out2$genome_table_out) <- out2$matches
# genome_data <- out2$genome_table_out[overlapping_otus,]
# kegg_data <- round(t(otu_data) %*% genome_data)
# 
# kegg_true <- as.matrix(biom_data(read_biom('~/Dropbox/stm_microbiome/qiime_active/gevers/picked_otus/ko_table.biom')))
# 
# format_gene_metadata <- function(gene_metadata){
#   
#   gene_metadata_list <- vector(mode='list')
#   for (i in seq_along(gene_metadata)){
#     
#     line <- gene_metadata[i]
#     values <- unlist(strsplit(gene_metadata[i],'\t'))
#     type <- values[1]
#     values <- values[-1]
#     
#     if (grepl('\\;',line)){
#       if (grepl('\\|',line)){
#         gene_metadata_list[[type]] <- lapply(values, function(x) lapply(strsplit(trimws(unlist(strsplit(x,'\\|'))),'\\;'),trimws))
#       }else{
#         gene_metadata_list[[type]] <- lapply(values,function(x) trimws(unlist(strsplit(x,'\\;'))))
#       }
#     }else{
#       gene_metadata_list[[list]] <- values
#     }
#   }
#   
#   return(gene_metadata_list)
#   
# }
# 
# format_gene_metadata(out2$genemeta)



*/
