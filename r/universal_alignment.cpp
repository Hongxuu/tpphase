#include <Rcpp.h>
using namespace Rcpp;

#define N_GENOME 2
CharacterVector iupac_to_char = {"-", "A", "C", "M", "G", "R", "S", "V", 
                                 "T", "W", "Y", "H", "K", "D", "B", "N"};

// [[Rcpp::export]]
List make_universal(List alignment, unsigned int for_hmm) { 
  IntegerVector reads = alignment["reads"];
  IntegerVector dim = alignment["dim"];
  unsigned int i, j, m;
  unsigned int n_alignment = dim.length()/N_GENOME;
  List universal(n_alignment);
  IntegerVector start_id(n_alignment);
  for(i = 1; i < n_alignment; ++i)
    start_id[i] = start_id[i - 1] + dim[(i - 1) * N_GENOME] * N_GENOME;
  
  unsigned int idx;
  CharacterVector seq(reads.length());
  
  if(!for_hmm)
    for(m = 0; m < reads.length(); ++m)
      seq[m] = iupac_to_char[reads[m]];
  
  for (i = 0; i < n_alignment; ++i) {
   CharacterVector uni_alignment(dim[i * N_GENOME]);
   for(j = 0; j < dim[i * N_GENOME]; ++j) {
     idx = j + start_id[i];
     if(for_hmm) {
       if(reads[idx] == 15 && reads[idx + dim[i * N_GENOME]] != 15) {
         uni_alignment[j] = "J";
       } else if(reads[idx] != 15 && reads[idx + dim[i * N_GENOME]] == 15) {
         uni_alignment[j] = "I";
       } else
         uni_alignment[j] = "M";
     } else {
       //Rcout << idx << "\t" << idx + dim[i * N_GENOME] << "\n";
       if(reads[idx] == reads[idx + dim[i * N_GENOME]]) {
         uni_alignment[j] = seq[idx];
       } else {
         uni_alignment[j] = "N";
       }
     }
   }
   universal[i] = uni_alignment;
  }
  return universal;
}

//[[Rcpp::export]]
void write_fasta(List seq, StringVector names, std::string outfile = "./out.fasta") {
  unsigned int i, j;
  FILE *fp = fopen(outfile.c_str(), "w");
  if (!fp)
    stop("Cannot open file");
  
  for(i = 0; i < seq.size(); ++i) {
    putc('>', fp);
    
    for(j = 0; j < names[i].size(); ++j)
      fprintf(fp, "%c", names[i][j]);
    fprintf(fp, "\n");
    
    CharacterVector read = seq[i];
    for(j = 0; j < read.size(); ++j)
      fprintf(fp, "%s", static_cast<const char*>(read[j]));
    fprintf(fp, "\n");
  }
  fclose(fp);
}






