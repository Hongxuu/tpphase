#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerVector adjust_variants(List nuc_unique, IntegerVector ref_start, IntegerVector pos) {
  unsigned int i, j, count = 0;
  IntegerVector exclude_pos(pos.size());
  
  for(j = 0; j < pos.size(); ++j) {
    IntegerVector genotype(3); // count of 3 genotypes
    int base_nuc = -1;
    for(i = 0; i < nuc_unique.size(); ++i) {
      List nucs = nuc_unique[i];
      int id = pos[j] - ref_start[i];
      if(id < 0) // not covering this site
        continue;
      IntegerVector nuc = nucs[id];
      Rcout << nuc << "\t||";
      if(nuc.size() == 1) {
        if(base_nuc != nuc[0]) {
          base_nuc = nuc[0];
          genotype[0]++;
        } else
          genotype[2]++;
      } else
          genotype[1]++;
    }
    Rcout << genotype << "\n";
    double total = sum(genotype);
    double prop = genotype[1]/total;
    if(prop > 0.55 || prop < 0.45) // not a variant site
      exclude_pos[count++] = j;
  }
  exclude_pos.erase(0, count);
  
  return pos[exclude_pos];
}
