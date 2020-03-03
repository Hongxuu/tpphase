#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace std;
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
List make_universal_old(List A_genome, List B_genome) {
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

// #define NUM_CLASS 4
// CharacterVector convert_char(IntegerVector nuc, int length) {
//   int i = 0;
//   CharacterVector nuc_char(length);
//   
//   for(i = 0; i < length; ++i) {
//     if(nuc[i] == 0) 
//       nuc_char[i] = 'A';
//     if(nuc[i] == 1)
//       nuc_char[i] = 'C';
//     if(nuc[i] == 2)
//       nuc_char[i] = 'T';
//     else if(nuc[i] == 3)
//       nuc_char[i] = 'G';
//   }
//   
//   return(nuc_char);
// }
// 
// arma::colvec solveC (NumericMatrix ar, NumericVector br) {
//   int n = ar.nrow(), k = ar.ncol();
//   
//   arma::mat a(ar.begin(), n, k, false);
//   arma::colvec b(br.begin(), br.size(), false);
//   
//   arma::colvec x = arma::solve(a, b);
//   //NumericVector re;
//   //re = as<NumericVector>(wrap(x)); 
//   return(x);
// }
// 
// // IntegerVector match_c (NumericVector ar, NumericVector br){
// //   unsigned int n = ar.length(), k = br.length();
// //   IntegerVector index(n);
// //   unsigned int m = 0;
// //   
// //   for(unsigned int i = 0; i < n; ++i) 
// //     for (unsigned int j = 0; j < k; ++j)
// //       if(ar[i] == br[j])
// //         index[m++] = i;
// //   
// //   return index(Range(0, m - 1));
// // }
// 
// DataFrame filter_c(DataFrame x, int condition, String selection) {
//   StringVector sub = x[selection];
//   LogicalVector ind(sub.size());
//   for (int i = 0; i < sub.size(); i++) {
//     ind[i] = (sub[i] == condition);
//   }
//   
//   // extracting each column into a vector
//   IntegerVector id = x["id"];
//   IntegerVector mode = x["mode"];
//   IntegerVector read_pos = x["read_pos"];
//   IntegerVector ref_pos = x["ref_pos"];
//   IntegerVector qua = x["qua"];
//   IntegerVector nuc = x["nuc"];
//   IntegerVector hap_nuc = x["hap_nuc"];
//   
//   return DataFrame::create(Named("id") = id[ind],
//                            Named("mode") = mode[ind],
//                                                Named("read_pos") = read_pos[ind],
//                                                                            Named("ref_pos") = ref_pos[ind],
//                                                                                                      Named("qua") = qua[ind],
//                                                                                                                        Named("nuc") = nuc[ind],
//                                                                                                                                          Named("hap_nuc") = hap_nuc[ind]);
// }
// // Find the mutural deletions in reads and haps
// //Slower than R
// IntegerVector distint_tab (DataFrame x, int nrow, String selection) {
//   IntegerVector id = x[selection];
//   IntegerVector length(nrow);
//   int unique = 0;
//   for (int m = 0; m < nrow; ++m) {
//     int i = id[m];
//     if (!length[i-1]) {
//       unique++;
//     }
//     ++length[i-1];
//   }
//   return length[Range(0, unique - 1)];
// }
// 
// 
// IntegerMatrix len_hapGap(List dat_info, List hap_info) {
//   List deletion = dat_info["deletion"];
//   IntegerVector length = dat_info["length"];
//   IntegerVector fake_length = dat_info["fake_length"];
//   int n_observation = dat_info["n_observation"];
//   
//   IntegerVector hap_ref_pos = hap_info["deletion_pos"];
//   IntegerVector hap_deletion_len = hap_info["hap_deletion_len"];
//   IntegerVector strat_id = hap_info["hap_del_start_id"];
//   unsigned int i, k, j, loci;
//   IntegerMatrix new_len(n_observation, NUM_CLASS);
//   IntegerVector cumsum_del_len(hap_deletion_len.size());
//   // cumsum_del_len = cumsum(hap_deletion_len);
//   partial_sum(hap_deletion_len.begin(), hap_deletion_len.end(), cumsum_del_len.begin());
//   
//   for(i = 0; i < n_observation; ++i) {
//     for (k = 0; k < NUM_CLASS; ++k) {
//       new_len(i, k) = length(i);
//       if (strat_id(k) != -1) {
//         for (j = 0; j < hap_deletion_len[k]; ++j) {
//           if (k != 0) {
//             loci = hap_ref_pos[j + cumsum_del_len[k - 1]];
//           } else
//             loci = hap_ref_pos[j];
//           if(loci <= fake_length(i))
//             new_len(i, k)--;
//         }
//       }
//     }
//   }
//   return new_len;
// }
// 
// 
// List dup_res(IntegerMatrix hap) {
//   unsigned int i, j, i1;
//   IntegerVector index;
//   IntegerVector key(hap.nrow());  
//   int count = 0;
//   for(i = 0; i < hap.nrow(); ++i) {
//     for(i1 = i + 1; i1 < hap.nrow(); ++i1) {
//       for(j = 0; j < hap.ncol(); ++j) {
//         if(hap[i1, j] != hap[i, j]) {
//           key(count++) = i;
//           break;
//         }
//       }
//     }
//   }
// }
// */
