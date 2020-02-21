#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List make_universal(List A_genome, List B_genome) {
  IntegerVector A_aligned = A_genome["reads"];
  IntegerVector B_aligned = B_genome["reads"];
  IntegerVector B_lengths = B_genome["dim"];
  
  unsigned int i, j;
  IntegerVector start_id(B_lengths.length());
  CharacterVector universal(B_aligned.length());
  
  for(i = 1; i < B_lengths.length(); ++i)
    start_id[i] = start_id[i - 1] + B_lengths[i - 1];
  
  for(j = 0; j < B_aligned.length(); ++j) {
    if(B_aligned[j] == 15 & A_aligned[j] != 15) {
      universal[j] = "J";
    } else if(B_aligned[j] != 15 & A_aligned[j] == 15) {
      universal[j] = "I";
    } else
      universal[j] = "M";
  }
  
  List ls = List::create(
    Named("universal_alignment") = universal,
    Named("start_id") = start_id);
  return ls;
}


