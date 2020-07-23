#include <RcppArmadillo.h>
#include "utils.h"
#include "data_format.h"
using namespace Rcpp;

#define NUM_CLASS 4
CharacterVector iupac_to_char = {"-", "A", "C", "M", "G", "R", "S",
                                 "V", "T", "W", "Y", "H", "K",
                                 "D", "B", "N"};
CharacterVector xy_to_c = {"A", "C", "T", "G"};

// [[Rcpp::export]]
double MEC(List dat_info, CharacterMatrix haps, IntegerVector cov_record) {
  double mec = 0;
  unsigned int i, j, k;
  int hap_min_pos = dat_info["ref_start"];
  IntegerVector index = dat_info["start_id"];
  IntegerVector length = dat_info["length"];
  IntegerVector ref_pos = dat_info["ref_pos"];
  CharacterMatrix referred_hap;
  IntegerVector missing;
  IntegerMatrix linkage;
  int len;
  int read_maxi = max(length);
  if(cov_record[0] != -1) {
    int full_len = cov_record[cov_record.size() - 1] - cov_record[0] + 1;
    IntegerVector full_record(full_len);
    full_record[0] = cov_record[0];
    for(i = 1; i < full_len; ++i)
      full_record[i] = full_record[i - 1] + 1;
    missing = setdiff(full_record, cov_record);
    IntegerVector true_match = match(missing, full_record); // 1-based
    len = missing.size();
    referred_hap = CharacterMatrix(NUM_CLASS, true_match.size());
    for(j = 0; j < len; ++j)
      true_match[j] = true_match[j] - 1;
    for(j = 0; j < len; ++j)
      referred_hap(_, j) = haps(_, true_match[j]);
  } else {
    len = haps.ncol();
    missing = IntegerVector(len);
    for(j = 0; j < len; ++j)
      missing[j] = j;
    referred_hap = haps;
  }
  linkage = linkage_info(dat_info, missing);
  unsigned int ref_j;
  // compute mec between referred_hap & reads
  for(i = 0; i < linkage.nrow(); ++i) {
    int max = read_maxi;
    unsigned int max_id = 0;
    for(k = 0; k < NUM_CLASS; ++k) {
       int hamming = 0;
        for(j = 0; j < len; ++j) {
            ref_j = missing[j] + hap_min_pos;
          if (ref_pos[index[i]] <= ref_j && ref_j <= ref_pos[index[i] + length[i] - 1]) {
            int nuc = linkage(i, j);
            if(referred_hap(k, j) == "N" && nuc == -1) {
              continue;
           } else if(referred_hap(k, j) != "N" && nuc == -1) {
              hamming++;
           } else if (xy_to_c[nuc] == referred_hap(k, j)) {
             continue;
           } else {
              hamming++;
             }
        }
      }
      if(hamming < max) {
        max = hamming;
        max_id = k;
      }
    }
    mec += max;
  }
  mec /= linkage.nrow();
  
  return(mec);
}

List switch_err(CharacterMatrix hmm_snp, CharacterMatrix real_snp, IntegerVector both,
                IntegerVector both_refer_type, IntegerVector both_true_type) {
  unsigned int k, j, i;
  int len = hmm_snp.ncol();
  IntegerVector status(len); // 0: A->A
  // Rcout << "refer " << both_refer_type << "\n";
  // Rcout << "true " << both_true_type << "\n";
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
    int tmp = i;
    for(int m = tmp; m < len; ++m) {
      if(status[m] != 10) {
        i = m;
        break;
      }
    }
      
    // Rcout <<j << " "<< status[j]  << "|| "<< i << " "<< status[i] << "\n";
    // switch error in A B 
    if(status[j] == 0 || status[j] == 2 || status[j] == 3
         || status[j] == 4 || status[j] == 5) {
      if(status[i] == 1 || status[i] == 6 || 
      status[i] == 7 || status[i] == 8 || status[i] == 9) {
        // Rcout << "switch\n";
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
        // Rcout << "switch\n";
        swAB_idx[count] = i;
        swAB_id[count++] = both[i];
        swi_AB++;
      }
    }
    i = tmp;
  }
  // add the end
  double homo_swi_err = double(swi_AB)/(non_swi_AB + swi_AB);
  swAB_idx[count] = len;
  swAB_id[count++] = both[len - 1] + 1;
  swAB_id.erase(count, len);
  swAB_idx.erase(count, len); // snp index (relative to whole hap) to seperate entire snps into blocks
  // switched heter within each correct homo region (get the switch error separately)
  // Rcout << "homo_swi " << swAB_id << "\n";
  // Rcout << "homo_swid "  << swAB_idx << "\n";
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
      if(status[i] == 10)
        continue;
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
        // Rcout << "switch\n";
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
        // Rcout << "switch\n";
        heter_sw[heter_c++] = both[inx_B[k]];
        swi_one[j]++;
      }
    }
  }
  List sw;
  if(heter_c) {
    NumericVector swi_er(count - 1);
    int c = 0;
    for(i = 0; i < count - 1; ++i) {
      if(swi_one[i] == 0 && non_swi_one[i] == 0)
        continue;
      swi_er[c++] = double(swi_one[i])/(swi_one[i] + non_swi_one[i]);
    }
    swi_er.erase(c, count - 1);
    double heter_se = median(swi_er);
    if(swi_AB) {
      sw = List::create (
        Named("homo_sw_id") = swAB_id,
        Named("homo_sw_err") = homo_swi_err,
        // Named("heter_nswitch") = non_swi_one,
        Named("heter_sw_err") = heter_se,
        Named("heter_sw_id") = heter_sw[Range(0, heter_c - 1)]);
    } else {
      sw = List::create (
        Named("heter_sw_id") = heter_sw[Range(0, heter_c - 1)],
        Named("homo_sw_err") = 0,
        Named("heter_sw_err") = heter_se);
    }
  } else if(!heter_c && swi_AB && non_swi_one) {
    sw = List::create (
      Named("homo_sw_id") = swAB_id,
      Named("homo_sw_err") = homo_swi_err,
      // Named("heter_nswitch") = non_swi_one,
      Named("heter_sw_err") = 0);
  } else {
    sw = List::create (
      Named("heter_sw_err") = 0,
      Named("homo_sw_err") = 0);
  }
 
  return(sw);
}
// [[Rcpp::export]] 
List find_snp(IntegerMatrix true_hap)  {
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
  List ls = List::create (
    Named("snp_id") = snp_id,
    Named("snp_type") = snp_type);
  return(ls);
}

List tf_table(IntegerVector inferred_snp_type, IntegerVector snp_location,
              IntegerVector snp_type, IntegerVector snp_id, int hap_length, int min_ref) {
  unsigned int j;
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
  List tf = List::create (
    Named("homo") = homo,
    Named("heter") = heter,
    Named("non") = non,
    Named("refer_match") = refer_match,
    Named("true_match") = true_match,
    Named("both") = both);
  return(tf);
}

// [[Rcpp::export]] 
List snp_stats_other (IntegerMatrix inferred_hap, int hap_length, int min_ref, IntegerMatrix true_hap) {
  unsigned int j, k;
  List true_info = find_snp(true_hap);
  IntegerVector snp_type = true_info["snp_type"];
  IntegerVector snp_id = true_info["snp_id"];
  
  List ref_info = find_snp(inferred_hap);
  IntegerVector inferred_snp_type = ref_info["snp_type"];
  IntegerVector snp_location = ref_info["snp_id"];
  
  List tf_info = tf_table(inferred_snp_type, snp_location, snp_type, snp_id, hap_length, min_ref);
  IntegerVector both = tf_info["both"];
  IntegerVector refer_match = tf_info["refer_match"];
  IntegerVector true_match = tf_info["true_match"];
  
  // switch error get the both snps first
  int len = both.size();
  CharacterMatrix hmm_snp(NUM_CLASS, len);
  CharacterMatrix real_snp(NUM_CLASS, len);
  
  for(j = 0; j < len; ++j) {
    for(k = 0; k < NUM_CLASS; ++k) {
      hmm_snp(k, j) = iupac_to_char[inferred_hap(k, snp_location[refer_match[j]])]; // to character
      real_snp(k, j) = iupac_to_char[true_hap(k, snp_id[true_match[j]])]; // to character
    }
  }
  IntegerVector both_refer_type = inferred_snp_type[refer_match];
  IntegerVector both_true_type = snp_type[true_match];
  
  List sw = switch_err(hmm_snp, real_snp, both, both_refer_type, both_true_type);
  IntegerVector homo = tf_info["homo"];
  IntegerVector heter = tf_info["heter"];
  IntegerVector non = tf_info["non"];
  DataFrame confusion_metric = DataFrame::create(_["homo"] = homo, _["heter"] = heter, _["non"] = non);
  
  List snp_info = List::create (Named("confusion metric") = confusion_metric,
                                       Named("tsnp_id") = snp_id,
                                       Named("switch") = sw,
                                       Named("both") = both);
  return snp_info;
}
// [[Rcpp::export]] 
// snps encoding 0 1 2 3: no, homologous snp, heter snp in A, heter in B
//construct a confusion matrix for multiclass
List snp_stats (CharacterMatrix inferred_snp, IntegerVector snp_location, int hap_length, int min_ref, 
                IntegerMatrix true_hap) {
  unsigned int j, k;
  List true_info = find_snp(true_hap);
  IntegerVector snp_type = true_info["snp_type"];
  IntegerVector snp_id = true_info["snp_id"];
  
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
  
  List tf_info = tf_table(inferred_snp_type, snp_location, snp_type, snp_id, hap_length, min_ref);
  IntegerVector both = tf_info["both"];
  IntegerVector refer_match = tf_info["refer_match"];
  IntegerVector true_match = tf_info["true_match"];
  
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
  
  List sw = switch_err(hmm_snp, real_snp, both, both_refer_type, both_true_type);
  IntegerVector homo = tf_info["homo"];
  IntegerVector heter = tf_info["heter"];
  IntegerVector non = tf_info["non"];
  DataFrame confusion_metric = DataFrame::create(_["homo"] = homo, _["heter"] = heter, _["non"] = non);
  
  List snp_info = List::create (
    Named("confusion metric") = confusion_metric,
    Named("tsnp_id") = snp_id,
    Named("switch") = sw,
    Named("both") = both
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
    // Rcout << "ppt " << ppt << "\n";
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

// NumericVector pp_all(IntegerVector hs_id, List combination, IntegerVector start_t, int hap_len,
//                      IntegerVector t_len, List gamma, IntegerVector ref_snp, IntegerVector t_snp) {
//   unsigned int t, m, j;
//   NumericVector pp(hap_len);
//   List pp_all(start_t.size());
//   //get the pp for each snp at each t
//   int count = 0;
//   for(t = 0; t < start_t.size(); ++t) {
//     NumericVector gam = gamma(t); 
//     IntegerMatrix comb = combination(t);
//     // print_intmat(comb);
//     IntegerVector chosed_comb = comb(hs_id[t], _);
//     // Rcout << "choosed " << chosed_comb << "\n";
//     NumericVector ppt(t_len.size());
//     
//     for(j = 0; j < comb.ncol(); ++j)
//       for(m = 0; m < comb.nrow(); ++m)
//         if(comb(m, j) == chosed_comb(j))
//           ppt[j] += exp(gam[m]);
//     count += snp_t.size();
//         // Rcout << "ppt " << ppt << "\n";
//     pp_snp[t] = ppt;
//   }
//   
//   IntegerVector snp_flat(count);
//   NumericVector pp_all(count);
//   count = 0;
//   // average over the snp pp
//   for(t = 0; t < snp_loci_t.size(); ++t) {
//     IntegerVector snp_t = snp_loci_t(t);
//     if(snp_t[0] == -1)
//       continue;
//     NumericVector ppt = pp_snp[t];
//     for(j = 0; j < snp_t.size(); ++j) {
//       pp_all(count) = ppt[j];
//       snp_flat(count++) = snp_t[j];
//     }
//   }
//   // Rcout << snp_flat << "\n";
//   List hased_info = hash_intvec(snp_flat);
//   IntegerVector value = hased_info["value"];
//   List all_id = hased_info["all_id"];
//   count = 0;
//   for(m = 0; m < value.size(); ++m) {
//     IntegerVector id = all_id[m];
//     NumericVector snp_each = pp_all[id];
//     pp[count++] = mean(snp_each);
//   }
//   return pp;
// }
  
