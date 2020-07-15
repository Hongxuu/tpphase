#include <RcppArmadillo.h>
#include "utils.h"
using namespace Rcpp;

#define NUM_CLASS 4
CharacterVector iupac_to_char = {"-", "A", "C", "M", "G", "R", "S",
                                 "V", "T", "W", "Y", "H", "K",
                                 "D", "B", "N"};
// [[Rcpp::export]] 
List switch_err(CharacterMatrix hmm_snp, CharacterMatrix real_snp, IntegerVector both,
                IntegerVector both_refer_type, IntegerVector both_true_type) {
  unsigned int k, j, i;
  int len = hmm_snp.ncol();
  IntegerVector status(len); // 0: A->A
  for(j = 0; j < len; ++j) {
    // first determine the order
    if(both_refer_type[j] == 1 && both_true_type[j] == 1) {
      if(hmm_snp(0, j) == real_snp(2, j))
        status[j] = 1;// A->B
    } else if(both_refer_type[j] == 2 && both_true_type[j] == 2) {
      if(hmm_snp(1, j) == real_snp(1, j) && hmm_snp(0, j) == real_snp(0, j))
        status[j] = 2; // A1->A1 && A2->A2
      else 
        status[j] = 3; // A1->A2 && A2->A1
    } else if(both_refer_type[j] == 3 && both_true_type[j] == 3) {
      if(hmm_snp(2, j) == real_snp(2, j) && hmm_snp(3, j) == real_snp(3, j))
        status[j] = 4; // B1->B1 && B2->B2
      else 
        status[j] = 5; // B1->B2 && B2->B1
    } else if(both_refer_type[j] == 2 && both_true_type[j] == 3) {
      if(hmm_snp(0, j) == real_snp(2, j) && hmm_snp(1, j) == real_snp(3, j))
        status[j] = 6; // A1->B1 && A2->B2
      else
        status[j] = 7; // A1->B2 && A2->B1
    } else if(both_refer_type[j] == 3 && both_true_type[j] == 2) {
      if(real_snp(0, j) == hmm_snp(2, j) && real_snp(1, j) == hmm_snp(3, j))
        status[j] = 8; //  B1->A1 && B2->A2
      else
        status[j] = 9; // B1->A2 && B2->A1
    } else
      status[j] = 10; // homo to heter or heter to homo(wrong anyway, should we consider this?)
  }
  // Rcout << status << "\n";
  // switched A B 
  IntegerVector swAB_id(len);
  IntegerVector swAB_idx(len);
  int swi_AB = 0, non_swi_AB = 0; 
  int count = 1; // start from the second position, the first
  swAB_id[0] = both[0];
  for(j = 0, i = 1; i < len; ++i, ++j) {
    if(status[j] == 10)
      continue;
    if(status[i] == 10)
      i = i + 1;
    // switch error in A B 
    if(status[j] == 0 || status[j] == 2 || status[j] == 3
         || status[j] == 4 || status[j] == 5) {
      if(status[i] == 1 || status[i] == 6 || 
      status[i] == 7 || status[i] == 8 || status[i] == 9) {
        swAB_idx[count] = i;
        swAB_id[count++] = both[i];
        swi_AB++;
      } else
        non_swi_AB++;
    } else {
      if(status[i] == 1 || status[i] == 6 || 
         status[i] == 7 || status[i] == 8 || status[i] == 9) {
        non_swi_AB++;
      } else {
        swAB_idx[count] = i;
        swAB_id[count++] = both[i];
        swi_AB++;
      }
    }
  }
  // add the end
  swAB_idx[count] = len;
  swAB_id[count++] = both[len - 1] + 1;
  swAB_id.erase(count, len);
  swAB_idx.erase(count, len); // snp index (relative to whole hap) to seperate entire snps into blocks
  // switched heter within each correct homo region (get the switch error separately)
  // Rcout << swAB_id << "\n";
  // Rcout << swAB_idx << "\n";
  IntegerVector swi_one(count - 1);
  IntegerVector non_swi_one(count - 1); 
  IntegerVector heter_sw(len);
  int heter_c = 0;
  
  for(j = 0; j < count - 1; ++j) {
    IntegerVector inx_A(len);
    IntegerVector inx_B(len);
    int count_a = 0, count_b = 0;
    // pick A sites and B sites
    for(i = swAB_idx[j]; i < swAB_idx[j + 1]; ++i) {
      if(both_refer_type[i] == 2) //A
        inx_A[count_a++] = i;
      else if(both_refer_type[i] == 3)
        inx_B[count_b++] = i;
    }
    for(i = 0, k = 1; k < count_a; ++i, ++k) {
      // Rcout << status[inx_A[i]] <<  status[inx_A[k]] << "\n";
      if((status[inx_A[i]] == 2 && status[inx_A[k]] == 2) || 
         (status[inx_A[i]] == 3 && status[inx_A[k]] == 3) || 
         (status[inx_A[i]] == 6 && status[inx_A[k]] == 6) || 
         (status[inx_A[i]] == 7 && status[inx_A[k]] == 7)) {
        non_swi_one[j]++;
      } else {
        heter_sw[heter_c++] = both[inx_A[k]];
        swi_one[j]++;
      }
    }
    for(i = 0, k = 1; k < count_b; ++i, ++k) {
      // Rcout << status[inx_B[i]] <<  status[inx_B[k]] << "\n";
      if((status[inx_B[i]] == 4 && status[inx_B[k]] == 4) || 
         (status[inx_B[i]] == 5 && status[inx_B[k]] == 5) || 
         (status[inx_B[i]] == 8 && status[inx_B[k]] == 8) || 
         (status[inx_B[i]] == 9 && status[inx_B[k]] == 9)) {
        non_swi_one[j]++;
      } else {
        heter_sw[heter_c++] = both[inx_B[k]];
        swi_one[j]++;
      }
    }
  }
  List sw;
  if(heter_c) {
    sw = List::create (
        Named("homo_switch_id") = swAB_id,
        Named("homo_switch") = swi_AB,
        Named("heter_nswitch") = non_swi_one,
        Named("heter_switch") = swi_one,
        Named("heter_switch_id") = heter_sw[Range(0, heter_c - 1)]);
  } else {
    sw = List::create (
      Named("homo_switch_id") = swAB_id,
      Named("homo_switch") = swi_AB,
      Named("heter_nswitch") = non_swi_one,
      Named("heter_switch") = swi_one);
  }
 
  return(sw);
}

// [[Rcpp::export]] 
// snps encoding 0 1 2 3: no, homologous snp, heter snp in A, heter in B
//construct a confusion matrix for multiclass
List snp_stats (CharacterMatrix inferred_snp, IntegerVector snp_location, int hap_length, int min_ref, 
                IntegerMatrix true_hap) {
  unsigned int k, j, count = 0;
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
  std::vector<int> s = as<std::vector<int> >(both);
  std::sort(s.begin(), s.end());
  both = wrap(s);
  IntegerVector refer_match = match(both, snp_location); // 1-based
  IntegerVector true_match = match(both, snp_id); // 1-based
  IntegerVector under = setdiff(snp_id, both);
  IntegerVector over = setdiff(snp_location, both);
  
  for(j = 0; j < both.size(); ++j) {
    refer_match[j] = refer_match[j] - 1;
    true_match[j] = true_match[j] - 1;
  }
  
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
    // Rcout << snp_type[true_match[j]] << " " << inferred_snp_type[refer_match[j]]<< "\n";
    if(snp_type[true_match[j]] == 1 && inferred_snp_type[refer_match[j]] != 1)
      homo2heter++;
    if(snp_type[true_match[j]] == 1 && inferred_snp_type[refer_match[j]] == 1)
      homo2homo++;
    if(snp_type[true_match[j]] != 1 && inferred_snp_type[refer_match[j]] != 1)
      heter2heter++;
    else if (snp_type[true_match[j]] != 1 && inferred_snp_type[refer_match[j]] == 1)
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
        // Rcout << snp_type[under_match[j]] << " " << inferred_snp_type[refer_match[j]]<< "\n";
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
  // switch error get the both snps first
  int len = both.size();
  CharacterMatrix hmm_snp(NUM_CLASS, len);
  CharacterMatrix real_snp(NUM_CLASS, len);
 
  for(j = 0; j < len; ++j) {
    hmm_snp(_, j) = inferred_snp(_, refer_match[j]);
    for(k = 0; k < NUM_CLASS; ++k) {
      real_snp(k, j) = iupac_to_char[true_hap(k, snp_id[true_match[j]])]; // to character
    }
  }
  IntegerVector both_refer_type = inferred_snp_type[refer_match];
  IntegerVector both_true_type = snp_type[true_match];
  
  // List sw = switch_err(hmm_snp, real_snp, both, both_refer_type, both_true_type);
  
  List snp_info = List::create (
    Named("confusion metric") = confusion_metric,
    Named("tsnp_id") = snp_id,
    Named("hmm_snp") = hmm_snp,
    Named("real_snp") = real_snp,
    Named("both") = both,
    Named("both_refer_type") = both_refer_type,
    Named("both_true_type") = both_true_type
    // Named("tsnp_type") = snp_type
  );
  return snp_info;
}

// [[Rcpp::export]] 
// get the posterior prob for each snp location
NumericVector pp_snp(IntegerVector hs_id, List combination, 
                     List snp_loci_t, List gamma, unsigned int num_snp) {
  unsigned int t, m, j;
  NumericVector pp(num_snp);
  List pp_snp(snp_loci_t.size());
  //get the pp for each snp at each t
  int count = 0;
  for(t = 0; t < snp_loci_t.size(); ++t) {
    IntegerVector snp_t = snp_loci_t(t);
    if(snp_t[0] == -1)
      continue;
    NumericVector gam = gamma(t); 
    IntegerMatrix comb = combination(t);
    // print_intmat(comb);
    IntegerVector chosed_comb = comb(hs_id[t], _);
    // Rcout << "choosed " << chosed_comb << "\n";
    NumericVector ppt(comb.ncol());
    
    for(j = 0; j < comb.ncol(); ++j)
      for(m = 0; m < comb.nrow(); ++m)
        if(comb(m, j) == chosed_comb(j))
          ppt[j] += exp(gam[m]);
    count += snp_t.size();
    Rcout << "ppt " << ppt << "\n";
    pp_snp[t] = ppt;
  }
  
  IntegerVector snp_flat(count);
  NumericVector pp_all(count);
  count = 0;
  // average over the snp pp
  for(t = 0; t < snp_loci_t.size(); ++t) {
    IntegerVector snp_t = snp_loci_t(t);
    if(snp_t[0] == -1)
      continue;
    NumericVector ppt = pp_snp[t];
    for(j = 0; j < snp_t.size(); ++j) {
      pp_all(count) = ppt[j];
      snp_flat(count++) = snp_t[j];
    }
  }
  // Rcout << snp_flat << "\n";
  List hased_info = hash_intvec(snp_flat);
  IntegerVector value = hased_info["value"];
  List all_id = hased_info["all_id"];
  count = 0;
  for(m = 0; m < value.size(); ++m) {
    IntegerVector id = all_id[m];
    NumericVector snp_each = pp_all[id];
    pp[count++] = mean(snp_each);
  }
  return pp;
}
// 1 homo, 2 heter in first 2, 3 heter in last 2

  
  
