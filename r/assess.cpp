#include <Rcpp.h>
using namespace Rcpp;

#define NUM_CLASS 4
// [[Rcpp::export]] 
// snps encoding 0 1 2 3: no, homologous snp, heter snp in A, heter in B
//construct a confusion matrix for multiclass
List snp_stats (CharacterMatrix inferred_snp, IntegerVector snp_location, int hap_length, int min_ref, 
                IntegerMatrix true_hap) {
  unsigned int j, count = 0;
  IntegerVector snp_type(true_hap.ncol());
  IntegerVector snp_id(true_hap.ncol());
  for(j = 0; j < true_hap.ncol(); ++j) {
    if(true_hap(0, j) != true_hap(1, j)) {
      snp_id[count] = j;
      snp_type[count++] = 2;
    }
    else if(true_hap(2, j) != true_hap(3, j)) {
      snp_id[count] = j;
      snp_type[count++] = 3;
    }
    else if(true_hap(0, j) != true_hap(2, j) && true_hap(1, j) != true_hap(3, j)) {
      snp_id[count] = j;
      snp_type[count++] = 1;
    }
  }
  snp_id.erase(count, true_hap.ncol());
  snp_type.erase(count, true_hap.ncol());
  IntegerVector both = intersect(snp_location, snp_id);
  IntegerVector refer_match = match(both, snp_location); // 1-based
  IntegerVector true_match = match(both, snp_id); // 1-based
  IntegerVector under = setdiff(snp_id, both);
  IntegerVector over = setdiff(snp_location, both);
  
  IntegerVector inferred_snp_type(snp_location.size());
  for(j = 0; j < snp_location.size(); ++j) {
    if(inferred_snp(0, j) != inferred_snp(1, j)) {
      inferred_snp_type[j] = 2;
    }
    else if(inferred_snp(2, j) != inferred_snp(3, j)) {
      inferred_snp_type[j] = 3;
    }
    else if(inferred_snp(0, j) != inferred_snp(2, j) && inferred_snp(1, j) != inferred_snp(3, j)) {
      inferred_snp_type[j] = 1;
    }
  }
  int homo2heter = 0, heter2homo = 0, heter2heter = 0, homo2homo = 0;
  int homo2non = 0, heter2non = 0, non2heter = 0, non2homo = 0;
  CharacterVector true_snp(NUM_CLASS);
  for(j = 0; j < both.size(); ++j) {
    // Rcout << snp_type[true_match[j] - 1] << " " << inferred_snp_type[refer_match[j] - 1]<< "\n";
    if(snp_type[true_match[j] - 1] == 1 && inferred_snp_type[refer_match[j] - 1] != 1)
      homo2heter++;
    if(snp_type[true_match[j] - 1] == 1 && inferred_snp_type[refer_match[j] - 1] == 1)
      homo2homo++;
    if(snp_type[true_match[j] - 1] != 1 && inferred_snp_type[refer_match[j] - 1] != 1)
      heter2heter++;
    else if (snp_type[true_match[j] - 1] != 1 && inferred_snp_type[refer_match[j] - 1] == 1)
      heter2homo++;
  }
  if(over) {
    IntegerVector over_match = match(over, snp_location);
    for(j = 0; j < over.size(); ++j) {
      if(inferred_snp_type[over_match[j] - 1] == 1)
        non2homo++;
      else
        non2heter++;
    }
  }
  if(under) {
    IntegerVector under_match = match(under, snp_id);
    int max_ref = min_ref + hap_length;
    for(j = 0; j < under.size(); ++j) {
      if(under[j] >= min_ref && under[j] <= max_ref) {
        // Rcout << snp_type[under_match[j] - 1] << " " << inferred_snp_type[refer_match[j] - 1]<< "\n";
        if(snp_type[under_match[j] - 1] == 1)
          homo2non++;
        else
          heter2non++;
      }
    }
  }
  int non2non = hap_length - snp_location.size();
  IntegerVector homo = {homo2homo, heter2homo, non2homo};
  IntegerVector heter = {homo2heter, heter2heter, non2heter};
  IntegerVector non = {homo2non, heter2non, non2non};
  DataFrame confusion_metric = DataFrame::create(_["homo"] = homo, _["heter"] = heter, _["non"] = non);
  
  List snp_info = List::create (
    Named("confusion metric") = confusion_metric,
    Named("tsnp_id") = snp_id
    // Named("tsnp_type") = snp_type
  );
  return snp_info;
}
