#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

#define NUM_CLASS 4
CharacterVector convert_char(IntegerVector nuc, int length) {
  int i = 0;
  CharacterVector nuc_char(length);
    
  for(i = 0; i < length; ++i) {
    if(nuc[i] == 0) 
      nuc_char[i] = 'A';
    if(nuc[i] == 1)
      nuc_char[i] = 'C';
    if(nuc[i] == 2)
      nuc_char[i] = 'T';
    else if(nuc[i] == 3)
      nuc_char[i] = 'G';
  }
  
  return(nuc_char);
}

arma::colvec solveC (NumericMatrix ar, NumericVector br) {
  int n = ar.nrow(), k = ar.ncol();
  
  arma::mat a(ar.begin(), n, k, false);
  arma::colvec b(br.begin(), br.size(), false);
  
  arma::colvec x = arma::solve(a, b);
  //NumericVector re;
  //re = as<NumericVector>(wrap(x)); 
  return(x);
}

// IntegerVector match_c (NumericVector ar, NumericVector br){
//   unsigned int n = ar.length(), k = br.length();
//   IntegerVector index(n);
//   unsigned int m = 0;
//   
//   for(unsigned int i = 0; i < n; ++i) 
//     for (unsigned int j = 0; j < k; ++j)
//       if(ar[i] == br[j])
//         index[m++] = i;
//   
//   return index(Range(0, m - 1));
// }

DataFrame filter_c(DataFrame x, int condition, String selection) {
  StringVector sub = x[selection];
  LogicalVector ind(sub.size());
  for (int i = 0; i < sub.size(); i++) {
    ind[i] = (sub[i] == condition);
  }
  
  // extracting each column into a vector
  IntegerVector id = x["id"];
  IntegerVector mode = x["mode"];
  IntegerVector read_pos = x["read_pos"];
  IntegerVector ref_pos = x["ref_pos"];
  IntegerVector qua = x["qua"];
  IntegerVector nuc = x["nuc"];
  IntegerVector hap_nuc = x["hap_nuc"];
  
  return DataFrame::create(Named("id") = id[ind],
                           Named("mode") = mode[ind],
                           Named("read_pos") = read_pos[ind],
                           Named("ref_pos") = ref_pos[ind],
                           Named("qua") = qua[ind],
                           Named("nuc") = nuc[ind],
                           Named("hap_nuc") = hap_nuc[ind]);
}
// Find the mutural deletions in reads and haps
//Slower than R
IntegerVector distint_tab (DataFrame x, int nrow, String selection) {
  IntegerVector id = x[selection];
  IntegerVector length(nrow);
  int unique = 0;
  for (int m = 0; m < nrow; ++m) {
    int i = id[m];
    if (!length[i-1]) {
      unique++;
    }
    ++length[i-1];
  }
  return length[Range(0, unique - 1)];
}

int top_n_map(const Rcpp::NumericVector & v)
{
  
  // Initialize a map
  std::map<double, int> Elt;
  Elt.clear();
  
  // Count each element
  for (int i = 0; i != v.size(); ++i)
    Elt[ v[i] ] += 1;
  
  // Find out how many unique elements exist... 
  int n_obs = Elt.size();
  // If the top number, n, is greater than the number of observations,
  // then drop it.  
  int n = n_obs;
  
  // Pop the last n elements as they are already sorted. 
  // Make an iterator to access map info
  std::map<double,int>::iterator itb = Elt.end();
  // Advance the end of the iterator up to 5.
  std::advance(itb, -n);
  
  // Recast for R
  NumericVector result_vals(n);
  NumericVector result_keys(n);
  
  unsigned int count = 0;
  // Start at the nth element and move to the last element in the map.
  for( std::map<double,int>::iterator it = itb; it != Elt.end(); ++it )
  {
    // Move them into split vectors
    result_keys(count) = it->first;
    result_vals(count) = it->second;
    count++;
  }
  
  int key = 0;
  for (unsigned int i = n; i --> 0;)
    if(result_vals(i) >= 4) {
      key = result_keys(i);
      break;
    }
  
  // return Rcpp::List::create(Rcpp::Named("lengths") = result_vals,
  //                           Rcpp::Named("values") = result_keys);
  return key;
}

IntegerMatrix len_hapGap(List dat_info, List hap_info) {
  List deletion = dat_info["deletion"];
  IntegerVector length = dat_info["length"];
  IntegerVector fake_length = dat_info["fake_length"];
  int n_observation = dat_info["n_observation"];
  
  IntegerVector hap_ref_pos = hap_info["deletion_pos"];
  IntegerVector hap_deletion_len = hap_info["hap_deletion_len"];
  IntegerVector strat_id = hap_info["hap_del_start_id"];
  unsigned int i, k, j, loci;
  IntegerMatrix new_len(n_observation, NUM_CLASS);
  IntegerVector cumsum_del_len(hap_deletion_len.size());
  // cumsum_del_len = cumsum(hap_deletion_len);
  std::partial_sum(hap_deletion_len.begin(), hap_deletion_len.end(), cumsum_del_len.begin());
  
  for(i = 0; i < n_observation; ++i) {
    for (k = 0; k < NUM_CLASS; ++k) {
      new_len(i, k) = length(i);
      if (strat_id(k) != -1) {
        for (j = 0; j < hap_deletion_len[k]; ++j) {
          if (k != 0) {
            loci = hap_ref_pos[j + cumsum_del_len[k - 1]];
          } else
            loci = hap_ref_pos[j];
          if(loci <= fake_length(i))
            new_len(i, k)--;
        }
      }
    }
  }
  return new_len;
}

/*
 * Viterbi Algorithm:
 * transition probablities: actgn
 * emmision probablities:  
 */
// List viterbi_hap(List par, List dat_info, unsigned int PD_LENGTH,
//                  unsigned int deletion, IntegerVector SNP, NumericMatrix transition_prob,
//                  NumericVector start_prob) {
//   NumericMatrix w_ic = par["wic"];
//   IntegerVector excluded_read = par["excluded_read"];
//   IntegerVector ref = dat_info["ref_pos"];
//   IntegerVector non_covered = dat_info["non_covered_site"];
//   NumericVector beta = par["beta"];
//   int hap_length = dat_info["ref_length_max"];
//   double del_rate = par["del_rate"]; 
//   double del_lambda = del_rate* hap_length;
//   double ins_rate = par["ins_rate"]; 
//   double ins_lambda = ins_rate* hap_length;
//   IntegerMatrix m_haplotype(NUM_CLASS, hap_length);
//   IntegerVector hap_deletion_len(NUM_CLASS);
//   unsigned int j, k, l, l1, max_id;
//   double max;
//   arma::cube e_jnk(hap_length, NONINDEL_STATE + 1, NUM_CLASS);
//   arma::cube path_knj(NUM_CLASS, NONINDEL_STATE + 1, hap_length);
//   
//   /*
//    * compute log emmision probability
//    */
//   for (k = 0; k < NUM_CLASS; ++k) {
//     for (j = 0; j < hap_length; ++j) {
//       for (l = 0; l < NONINDEL_STATE + 1; ++l) {
//         e_jnk(k, j, l) = m_haplotype_llk(j, l, k, PD_LENGTH, beta, dat_info, w_ic, excluded_read, 
//               ins_lambda, del_lambda, 1);
//       }
//     }
//   }
//   
//   /*
//    * compute log probability table to trace back the hidden state
//    */
//   double max_prob = 0;
//   for (k = 0; k < NUM_CLASS; ++k) {
//     for (j = 0; j < hap_length; ++j) {
//       for (l = 0; l < NONINDEL_STATE + 1; ++l) {
//         if (j == 0) {
//           path_knj(k, j, l) = e_jnk(k, j, l);
//         }
//         else {
//           max = -INFINITY;
//           for (l1 = 0; l1 < NONINDEL_STATE + 1; ++l1) {
//             max_prob = path_knj(k, j - 1, l1) + transition_prob(l1, l);
//             if (max_prob > max)
//               max = max_prob;
//           }
//           path_knj(k, j, l) = e_jnk(k, j, l) + max;
//         }
//       }
//     }
//   }
//   
//   /*
//    * find the path (backtrace)
//    */
//   for (k = 0; k < NUM_CLASS; ++k) {
//     for (j = hap_length; j > -1; --j) {
//       max = -INFINITY;
//       max_id = 0;
//       for (l = 0; l < NONINDEL_STATE + 1; ++l) {
//         if (path_knj(k, j, l) > max) {
//           max_id = l;
//           max = path_knj(k, j, l);
//         }
//       }
//       m_haplotype(k, j) = max_id;
//       if(m_haplotype(k, j) == 4)
//         hap_deletion_len[k]++;
//     }
//   }
//   
//   unsigned int sum = 0;
//   for(k = 0; k < NUM_CLASS; ++k)
//     sum += hap_deletion_len[k];
//   IntegerVector hap_ref_pos(sum);
//   IntegerVector strat_id(NUM_CLASS);
//   
//   sum = 0;
//   for (k = 0; k < NUM_CLASS; ++k) {
//     for (j = 0; j < hap_length; ++j) {
//       if (m_haplotype(k, j) == 4)
//         hap_ref_pos[sum++] = j;
//     }
//     if (hap_deletion_len[k] != 0) {
//       if (k != 0)
//         for (l = 0; l < k; ++l)
//           strat_id[k] += hap_deletion_len[l];
//     }
//     else
//       strat_id[k] = -1; //record the deletion starting id in vector hap_ref_pos for each hap
//   }
//   
//   List ls = List::create(
//     Named("hap") = m_haplotype,
//     Named("deletion_pos") = hap_ref_pos,
//     Named("hap_deletion_len") = hap_deletion_len,
//     Named("hap_del_start_id") = strat_id);
//   return(ls);
// }


// /*** R
// ### Slower
// rcpp_inc <- '
// using namespace Rcpp;
// using namespace arma;
// '
// src_cp <- '
// mat m1 = as<mat>(m1in);
// mat m2 = as<mat>(m2in);
// mat cp = trans(m1) * m2;
// return(wrap(cp));
// '
// library(inline)
// crossprodC <- cxxfunction(signature(m1in="numeric", m2in="numeric"), src_cp, plugin='RcppArmadillo', rcpp_inc)
// */
