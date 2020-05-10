#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <iostream>
#include <stdio.h>
#include <string>
#include <unordered_map>

using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
NumericVector comb_weight(NumericVector weight, List hash_idx) {
  int len = hash_idx.size();
  NumericVector new_weight(len);
  for (int i = 0; i < len; ++i) {
    // subset weight
    IntegerVector idx = hash_idx[i];
    NumericVector weig = weight[idx];
    for(int j = 0; j < weig.size(); ++j) {
      new_weight[i] += weig[j];
    }
  }
  return(new_weight);
}
// struct VectorHasher {
//   int operator()(const vector<int> &V) const {
//     int hash=0;
//     for(int i=0;i<V.size();i++) {
//       hash^=V[i];
//     }
//     return hash;
//   }
// };
// 
// std::unordered_map<intvec, int, VectorHasher> map;

// trans_indicator: indicate which state can transfer to which, further_limit indicate some states should not be considered
// List prepare_ini_hmm2 (unsigned int t_max, IntegerVector num_states, List trans_indicator, List further_limit) {
//   List trans_new_ind(t_max - 1);
//   // List emit(t_max);
//   unsigned int w, m, t;
//   // 
//   // IntegerVector more_limit = further_limit[0];
//   // IntegerMatrix trans_ind = trans_indicator[0];
//   // IntegerMatrix trans_id_new(num_states[0], num_states[1]);
//   // int count = 0;
//   // for (w = 0; w < trans_ind.ncol(); ++w)
//   //   if(!more_limit[w])
//   //     trans_id_new(_, count++) = trans_ind(_, w);
//   //   trans_new_ind[0] = trans_id_new;
//   
//   for(t = 0; t < t_max - 1; ++t) {
//     // Rcout << t << "\n";
//     IntegerVector more_limit_last = further_limit[t];
//     IntegerVector more_limit = further_limit[t + 1];
//     IntegerMatrix trans_ind = trans_indicator[t];
//     IntegerMatrix trans_ind_new(trans_ind.nrow(), num_states[t + 1]);
//     int count2 = 0;
//     int count = 0;
//     for (w = 0; w < trans_ind.ncol(); ++w)
//       if(!more_limit[w])
//         trans_ind_new(_, count++) = trans_ind(_, w);
//       for (m = 0; m < trans_ind.nrow(); ++m)
//         if(!more_limit_last[m])
//           trans_ind_new2(count2++, _) = trans_ind_new(m, _);
//         trans_new_ind[t] = trans_ind_new2;
//   }
//   
//   return(trans_new_ind);
// }
// 
// List trans_permit(IntegerVector num_states, List combination, List loci, int t_max) {
//   List trans_permits(t_max - 1);
//   unsigned int t, j, m, w;
//   
//   for(t = 0; t < t_max - 1; ++t) {
//     if(num_states[t + 1] != 1) {
//       IntegerMatrix trans(num_states[t], num_states[t + 1]);
//       IntegerMatrix comb_t1 = combination[t];
//       IntegerMatrix comb_t2 = combination[t + 1];
//       IntegerVector location_t1 = loci[t];
//       IntegerVector location_t2 = loci[t + 1];
//       // Rcout << location_t1 << "\n";
//       // Rcout << location_t2 << "\n";
//       // get the overlapped region, this moght be different from the overlapped states we had
//       int id = 0;
//       for(j = 0; j < location_t1.size(); ++j) 
//         if(location_t1[j] == location_t2[0]) {
//           id = j;
//           break;
//         }
//       int end = 0;
//       for(j = id; j < location_t1.size(); ++j)
//         if(location_t1[j] == location_t2[j - id])
//           end = j;
//         
//         for(m = 0; m < num_states[t]; ++m) {
//           IntegerVector hap_t1 = comb_t1(m, _);
//           for(w = 0; w < num_states[t + 1]; ++w) {
//             IntegerVector hap_t2 = comb_t2(w, _);
//             for(j = id; j < end; ++j)
//               if (hap_t1(j) != hap_t2(j - id)) {
//                 trans(m, w) = 1; // represents m cannot transfer to w
//                 break;
//               }
//           }
//         }
//         trans_permits(t) = trans;
//     } else {
//       IntegerMatrix temp(num_states[t], num_states[t + 1]);
//       trans_permits(t) = temp;
//     }
//   }
// 
//   return(trans_permits);
// }

  // List ini_hmm (unsigned int t_max,  IntegerVector num_states, List trans_indicator) {
  //   NumericVector phi(num_states[0]);
  //   List trans(t_max - 1);
  //   // List emit(t_max);
  //   unsigned int w, m, t;
  //   
  //   for(w = 0; w < num_states[0]; ++w)
  //     phi[w] = 1/double(num_states[0]);
  //   
  //   for(t = 0; t < t_max - 1; ++t) {
  //     NumericMatrix transition(num_states[t], num_states[t + 1]);
  //     // 
  //     // if(trans_indicator.isNotNull()) {
  //     IntegerMatrix trans_ind = trans_indicator[t];
  //     for (w = 0; w < num_states[t]; ++w) {
  //       int new_num = num_states[t + 1];
  //       for (m = 0; m < num_states[t + 1]; ++m)
  //         if (trans_ind(w, m)) // 1 means not the same, so cannot b transferred
  //           new_num--;
  //         // Rcout << new_num << "\n";
  //         for (m = 0; m < num_states[t + 1]; ++m)
  //           if (!trans_ind(w, m))
  //             transition(w, m) = 1/double(new_num);
  //     }
  //     // } else {
  //     // for (w = 0; w < num_states[t]; ++w)
  //     //   for (m = 0; m < num_states[t + 1]; ++m) 
  //     //     transition(w, m) = 1/double(num_states[t + 1]);
  //     // }
  //     trans(t) = transition;
  //   }
  //   // for(t = 0; t < t_max; ++t) {
  //   //   NumericVector emission(num_states[t]);
  //   //   for (w = 0; w < num_states[t]; ++w)
  //   //     emission(w) = 1/double(num_states[t]);
  //   //   emit(t) = emission;
  //   // }
  //   // 
  //   List par_hmm = List::create(
  //     Named("phi") = phi,
  //     // Named("emit") = emit
  //     Named("trans") = trans);
  //   return(par_hmm);
  // }
  // List trans_permit2(IntegerVector num_states, List combination, int t_max, IntegerVector undecided_pos, 
  //                   IntegerVector time_pos, IntegerVector p_tmax, int hap_min_pos) {
  //   List trans_permits(t_max - 1);
  //   List further_limit_col(t_max - 1);
  //   List further_limit_row(t_max - 1);
  //   IntegerVector new_num_states(t_max);
  //   unsigned int t, j, m, w;
  //   int end_t1, begin_t2, end_t2;
  //   for (t = 0; t < t_max; ++t)
  //     new_num_states[t] = num_states[t];
  //   
  //   for (t = 0; t < t_max - 1; ++t) {
  //     if(num_states[t + 1] != 1 && num_states[t] != 1) { // assume if no. state = 1, then it can transfer to all
  //       end_t1 = time_pos[t] + p_tmax[t];
  //       begin_t2 = time_pos[t + 1];
  //       end_t2 = time_pos[t + 1] + p_tmax[t + 1];
  //       if(end_t1 > end_t2)
  //         end_t1 = end_t2; // take the smaller one
  //       IntegerMatrix full_hap_t1 = combination(t);
  //       IntegerMatrix full_hap_t2 = combination(t + 1);
  //       IntegerMatrix trans(num_states[t], num_states[t + 1]);
  //       int count = 0; // number of overlapped variational sites
  //       int count1 = 0;
  //       
  //       for(j = 0; j < undecided_pos.size(); ++j)
  //         if (undecided_pos[j] + hap_min_pos < end_t1 && undecided_pos[j] >= time_pos[t]) {
  //           count1++;
  //           if(undecided_pos[j] + hap_min_pos >= begin_t2)
  //             count++;
  //         }
  //         // Rcout << t << "\t" << count << "\n";
  //       int left = count1 - count;
  //         for(m = 0; m < num_states[t]; ++m) {
  //           IntegerVector hap_t1 = full_hap_t1(m, _);
  //           // Rcout << hap_t1 << "\n";
  //           for(w = 0; w < num_states[t + 1]; ++w) {
  //             IntegerVector hap_t2 = full_hap_t2(w, _);
  //             // Rcout << hap_t2 << "\t";
  //             for(j = 0; j < count; ++j)
  //               if (hap_t2(j) != hap_t1(j + left)) {
  //                 trans(m, w) = 1; // represents m cannot transfer to w
  //                 break;
  //               }
  //           }
  //         }
  //         trans_permits(t) = trans;
  //     } else {
  //       IntegerMatrix temp(num_states[t], num_states[t + 1]);
  //       trans_permits(t) = temp;
  //     }
  //   }
  //   // label the states that should be excluded (go forward)
  //   for(t = 0; t < t_max - 1; ++t) {
  //     if(num_states[t + 1] != 1 || num_states[t] != 1) { 
  //     // if for a row, all == 1 or for a column all == 1, then these two should be excluded
  //     IntegerVector more_limits_row(num_states[t]);
  //     IntegerVector more_limits_col(num_states[t + 1]);
  //     IntegerVector more_limits_col_last;
  //     IntegerMatrix trans = trans_permits(t);
  //     int exluded_col = 0;
  //     for(w = 0; w < num_states[t + 1]; ++w) {
  //       int num_col = 0;
  //       for(m = 0; m < num_states[t]; ++m) {
  //         if(trans(m, w) == 1)
  //           num_col++;
  //       }
  //       if(num_col == num_states[t]){
  //         more_limits_col[w] = 1; // meaning no state transfer to it
  //         exluded_col++;
  //       }
  //     }
  //     for(m = 0; m < num_states[t]; ++m) {
  //       int num_row = 0;
  //       for(w = 0; w < num_states[t + 1]; ++w) {
  //         if(trans(m, w) == 1 && more_limits_col[w] != 1)
  //           num_row++;
  //       }
  //       if(num_row == num_states[t + 1] - exluded_col){
  //         more_limits_row[m] = 1; // meaning it transfer to no state, considering the excluded column 
  //       }
  //     }
  //     // found the states that never appears
  //     further_limit_col(t) = more_limits_col;
  //     further_limit_row(t) = more_limits_row;
  //   } else {
  //     IntegerVector more_limits_col(num_states[t + 1]);
  //     IntegerVector more_limits_row(num_states[t]);
  //     further_limit_col(t) = more_limits_col;
  //     further_limit_row(t) = more_limits_row;
  //   }
  //   }
  //   // go backward
  //   for(t = t_max - 1; t <= 0; --t) {
  //     if(num_states[t + 1] != 1 || num_states[t] != 1) { 
  //       // if for a row, all == 1 or for a column all == 1, then these two should be excluded
  //       IntegerVector more_limits_row = further_limit_row(t);
  //       IntegerVector more_limits_col = further_limit_col[t - 1];
  //       
  //       int index;
  //       for(w = 0; w < more_limits_col.size(); ++w) {
  //         if(more_limits_row[w] != more_limits_col[w])
  //           index = 0;
  //       }
  //      
  //       // found the states that never appears
  //       further_limit_col(t) = more_limits_col;
  //       further_limit_row(t) = more_limits_row;
  //     } 
  //   }
  //   // get new number of states
  //   List further_limit(t_max);
  //   IntegerVector more_limits_row = further_limit_row[0];
  //   for (m = 0; m < num_states[0]; ++m)
  //     if (more_limits_row[m] == 1)
  //       new_num_states[0]--;
  //   further_limit[0] = further_limit_row[0];
  //   IntegerVector more_limits_col = further_limit_col[t_max - 2];
  //   for(w = 0; w < num_states[t_max - 1]; ++w) {
  //       if(more_limits_col[w] == 1)
  //         new_num_states[t_max - 1]--;
  //   }
  //   further_limit[t_max - 1] = further_limit_col[t_max - 2];
  //   
  //   for(t = 1; t < t_max - 1; ++t) {
  //     more_limits_row = further_limit_row(t);
  //     more_limits_col = further_limit_col(t - 1);
  //     IntegerVector combine(more_limits_row.size());
  //     // find the union of these two vector
  //     for(w = 0; w < num_states[t]; ++w) {
  //       if(more_limits_col[w] == 1 || more_limits_row[w] == 1)
  //         combine[w] = 1;
  //     }
  //     further_limit[t] = combine; // indicates some states cannot appear
  //     for(w = 0; w < num_states[t]; ++w)
  //       if(combine[w] == 1)
  //         new_num_states[t]--;
  //   }
  //   
  //   List out = List::create(
  //     Named("trans_permits") = trans_permits,
  //     Named("further_limit") = further_limit, 
  //     Named("new_num_states") = new_num_states);
  //   
  //   return(out);
  // }

// vector<vector<int> > cart_product (const vector<vector<int> > &v) {
//   vector<vector<int> > s = {{}};
//   for (const auto& u : v) {
//     vector<vector<int> > r;
//     for (const auto& x : s) {
//       for (const auto y : u) {
//         r.push_back(x);
//         r.back().push_back(y);
//       }
//     }
//     s = move(r);
//   }
//   return s;
// }
// 
// IntegerMatrix call_cart_product(IntegerVector len) {
//   unsigned int row = len.size();
//   vector<vector<int> > vec(row);
//   unsigned int col, count, i, j;
//   for (i = 0; i < row; i++) {
//     count = 1;
//     col = len[i];
//     vec[i] = vector<int>(col);
//     for (j = 0; j < col; j++)
//       vec[i][j] = count++;
//   }
//   vector<vector<int> > res = cart_product(vec);
//   IntegerMatrix out(res.size(), row);
//   for(i = 0; i < res.size(); ++i)
//     for(j = 0; j < row; ++j) 
//       out(i, j) = res[i][j] - 1; //minus 1 for the index in C
//   
//   return(out);
// }
// 
// List find_combination(IntegerVector undecided_pos, IntegerVector pos_possibility, 
//                       unsigned int p_tmax, unsigned int time_pos, int hap_min_pos) {
//   //possible combination of the rest non-unique loci
//   IntegerVector location(pos_possibility.size());
//   IntegerVector location_len(pos_possibility.size());
//   unsigned int num = 0;
//   for(unsigned int i = 0; i < pos_possibility.size(); ++i)
//     if(time_pos - hap_min_pos <= undecided_pos[i] && undecided_pos[i] < time_pos + p_tmax - hap_min_pos) {
//       location(num) = undecided_pos[i];
//       location_len(num++) = pos_possibility[i];
//     }
//     IntegerMatrix combination = call_cart_product(location_len[Range(0, num - 1)]);
//     List ls = List::create(
//       Named("combination") = combination, // possible comb
//       Named("num") = num, // number of comb sites at time t
//       Named("location") = location[Range(0, num - 1)]); // undecided site at time t [here assume alignment starts from 0]
//     return(ls);
// }

// 
// List limit_comb(IntegerMatrix combination, List hidden_states, IntegerVector location,
//                 IntegerMatrix linkage_info, unsigned int num, unsigned int start_idx, unsigned int num_states) {
//   unsigned int i, j, k, idx, m;
//   unsigned int n_observation = linkage_info.nrow();
//   IntegerMatrix sub_hap(NUM_CLASS, num);
//   IntegerMatrix sub_link = linkage_info(_, Range(start_idx, start_idx + num - 1));
//   IntegerVector exclude(num_states);
//   int count, linkage_len, all_excluded;
//   linkage_len = num/2;
//   
//   for(i = 0 ; i < n_observation; ++i)
//     for (j = 0; j < num; ++j)
//      Rcout << sub_link(i, j) << "\t";
//   Rcout << "\n";
//   all_excluded = num_states;
//   while (all_excluded == num_states) {
//     all_excluded = 0;
//     Rcout << "linkage length " << linkage_len << "\n"; //actual linkage length + 1
//     for (m = 0; m < num_states; ++m) {
//       exclude(m) = 0;
//       IntegerVector comb = combination(m, _);
//       count = 0;
//       for (k = 0; k < NUM_CLASS; ++k) {
//         for (j = 0; j < num; ++j) {
//          
//           IntegerMatrix hidden = hidden_states[location[j]];
//           idx = comb[j];
//           sub_hap(k, j) = hidden(idx, k);
//           Rcout << sub_hap(k, j)  << "\t";
//         }
//         Rcout << "\n";
//         for (i = 0; i < n_observation; i++) {
//           int flag = 0;
//           for (j = 0; j < num - 1; ++j) {
//             // if (sub_link(i, j) != -1 && sub_link(i, j) != 4)
//             if (sub_hap(k, j) == sub_link(i, j) && sub_hap(k, j + 1) == sub_link(i, j + 1))
//               flag++;
//           }
//           if (flag >= linkage_len) {
//             count++;
//             break;
//           }
//         }
//       }
//       if (count != NUM_CLASS) {
//         exclude(m) = 1;
//         all_excluded++;
//       }
//     }
//     linkage_len--;
//   }
//   List out = List::create(
//     Named("num_states") = num_states - all_excluded,
//     Named("exclude") = exclude);
//   return(out);
// }

// get the overlapped variational region at each time t, notice T1: 1-30, T2: 10-20, T3: 15-40/22-45, then T3 has overlap with T1 but not T2
// List get_overlap(IntegerVector p_tmax, IntegerVector time_pos, IntegerVector num_states,
//                  IntegerVector undecided_pos, unsigned int t_max, int hap_min_pos)
// {
//   unsigned int start_t = 0;
//   // get the first t which has variation
//   for (unsigned int t = 0; t < t_max; ++t)
//     if(num_states[t] > 1) {
//       start_t = t;
//       break;
//     }
//     List overlapped(t_max);
//     List location(t_max);
//     IntegerVector overlapped_idx(t_max);
//     unsigned int begin, end, end1, min;
//     
//     for (unsigned int t = 0; t < t_max; ++t) {
//       
//       if (num_states[t] == 1) {
//         overlapped[t] = -1;
//         overlapped_idx[t] = -1;
//         location[t] = -1;
//         continue;
//       }
//       if(t == start_t) {
//         overlapped[t] = -1;
//         overlapped_idx[t] = -1;
//         begin = time_pos[start_t] - hap_min_pos;
//         end = time_pos[start_t] + p_tmax[start_t] - hap_min_pos;
//         int num = 0;
//         IntegerVector location_t(undecided_pos.size());
//         for (unsigned int m = 0; m < undecided_pos.size(); ++m)
//           if (undecided_pos[m] >= begin && undecided_pos[m] < end)
//             location_t(num++) = undecided_pos[m];
//         location[start_t] = location_t[Range(0, num - 1)];
//         continue;
//       }
//       begin = time_pos[t] - hap_min_pos;
//       end = time_pos[t] + p_tmax[t] - hap_min_pos;
//       // store location
//       int num = 0;
//       IntegerVector location_t(undecided_pos.size());
//       for (unsigned int m = 0; m < undecided_pos.size(); ++m)
//         if (undecided_pos[m] >= begin && undecided_pos[m] < end)
//           location_t(num++) = undecided_pos[m];
//       location[t] = location_t[Range(0, num - 1)];
//       int len = 0;
//       int id_t = 0;
//       int index = 0;
//       // Rcout << "location" << location_t[num - 1] << "\n";
//       for (unsigned int t1 = 0; t1 < t; ++t1) {
//         IntegerVector last_location = location[t1];
//         // begin1 = time_pos[t1] - hap_min_pos;
//         end1 = last_location[last_location.size() - 1];
//         // Rcout << end1 << "\t";
//         // end1 = time_pos[t1] + p_tmax[t1] - hap_min_pos;
//         if(begin <= end1) {
//           min = end1;
//           // minimum overlapped region
//           if(location_t[num - 1] < end1)
//             min = location_t[num - 1];
//           // find the time t which has the longest coverage
//           if(min - location_t[0] > len) {
//             // Rcout << t << " has coverage with " << t1 << "\n";
//             len = min - location_t[0];
//             id_t = min;
//             index = t1;
//           }
//         }
//       }
//       // Rcout<< "\n" << begin << "\t" << id_t << "\n";
//       int count = 0;
//       IntegerVector position(undecided_pos.size());
//       for (unsigned int m = 0; m < undecided_pos.size(); ++m) 
//         if (undecided_pos[m] >= begin && undecided_pos[m] <= id_t)
//           position(count++) = undecided_pos[m];
//       // Rcout << position << "\n";
//       if(count) { 
//         overlapped[t] = position[Range(0, count - 1)];
//         overlapped_idx[t] = index;
//       }
//     }
//     
//     List overlap = List::create(
//       Named("location") = location,
//       Named("overlapped") = overlapped,
//       Named("overlapped_id") = overlapped_idx, 
//       Named("start_t") = start_t);
//     return(overlap);
// }
// 
// List limit_comb_t0(IntegerMatrix combination, List hidden_states, IntegerVector location,
//                    IntegerMatrix linkage_info, unsigned int num, unsigned int start_idx, unsigned int num_states) {
//   unsigned int i, j, k, idx, m;
//   unsigned int n_observation = linkage_info.nrow();
//   IntegerMatrix sub_hap(NUM_CLASS, num);
//   IntegerMatrix sub_link = linkage_info(_, Range(start_idx, start_idx + num - 1));
//   IntegerVector exclude(num_states);
//   int count, linkage_len, all_excluded;
//   linkage_len = num/2; // change the linkage length to be the length appears in the read
//   int cut_off;
//   all_excluded = num_states;
//   while (all_excluded == num_states) {
//     cut_off = NUM_CLASS;
//     // Rcout << "linkage length " << linkage_len << "\n"; //actual linkage length + 1
//     while (cut_off >= 1 && all_excluded == num_states) {
//       all_excluded = 0;
//       for (m = 0; m < num_states; ++m) {
//         exclude(m) = 0;
//         IntegerVector comb = combination(m, _);
//         // Rcout << comb << "\t\t";
//         count = 0;
//         for (k = 0; k < NUM_CLASS; ++k) {
//           // Rcout << "k" << k << "\n";
//           for (j = 0; j < num; ++j) {
//             IntegerMatrix hidden = hidden_states[location[j]];
//             idx = comb[j];
//             sub_hap(k, j) = hidden(idx, k);
//             // Rcout << sub_hap(k, j) << "\t";
//           }
//           // Rcout << "read" << "\n";
//           for (i = 0; i < n_observation; i++) {
//             int flag = 0;
//             for (j = 0; j < num - 1; ++j) {
//               // Rcout << sub_link(i, j) << "\t" << sub_link(i, j + 1);
//               // if (sub_link(i, j) != -1 && sub_link(i, j) != 4)
//               if (sub_hap(k, j) == sub_link(i, j) && sub_hap(k, j + 1) == sub_link(i, j + 1))
//                 flag++;
//             }
//             // Rcout << "\n" << flag << "\n" ;
//             if (flag >= linkage_len) {
//               count++;
//               break;
//             }
//           }
//         }
//         // Rcout << "count "<< count << "\n";
//         if (count != cut_off) {
//           exclude(m) = 1;
//           all_excluded++;
//         }
//       }
//       cut_off--;
//     }
//     linkage_len--;
//   }
//   // Rcout << exclude << "\n";
//   List out = List::create(
//     Named("num_states") = num_states - all_excluded,
//     Named("exclude") = exclude);
//   return(out);
// }
// 
// 
// template <typename T>
// inline bool approx_equal_cpp(const T& lhs, const T& rhs, double tol = 0.00000001) {
//   return arma::approx_equal(lhs, rhs, "absdiff", tol);
// }
// 
// arma::mat unique_rows(const arma::mat& m) {
//   
//   arma::uvec ulmt = arma::zeros<arma::uvec>(m.n_rows);
//   
//   for (arma::uword i = 0; i < m.n_rows; i++) {
//     for (arma::uword j = i + 1; j < m.n_rows; j++) {
//       if (approx_equal_cpp(m.row(i), m.row(j))) { ulmt(j) = 1; break; }
//     }
//   }
//   
//   return m.rows(find(ulmt == 0));
//   
// }
// 
// List unique_map(const Rcpp::IntegerVector & v)
// {
//   // Initialize a map
//   std::map<double, int> Elt;
//   Elt.clear();
//   
//   // Count each element
//   for (int i = 0; i != v.size(); ++i)
//     Elt[ v[i] ] += 1;
//   
//   // Find out how many unique elements exist... 
//   int n_obs = Elt.size();
//   // If the top number, n, is greater than the number of observations,
//   // then drop it.  
//   int n = n_obs;
//   
//   // Pop the last n elements as they are already sorted. 
//   // Make an iterator to access map info
//   std::map<double,int>::iterator itb = Elt.end();
//   // Advance the end of the iterator up to 5.
//   std::advance(itb, -n);
//   
//   // Recast for R
//   NumericVector result_vals(n);
//   NumericVector result_keys(n);
//   
//   unsigned int count = 0;
//   // Start at the nth element and move to the last element in the map.
//   for (std::map<double,int>::iterator it = itb; it != Elt.end(); ++it) {
//     // Move them into split vectors
//     result_keys(count) = it->first;
//     result_vals(count) = it->second;
//     count++;
//   }
//   return List::create(Named("lengths") = result_vals,
//                       Named("values") = result_keys);
// }
// 
// IntegerMatrix unique_overlap(IntegerVector overlapped, IntegerVector exclude_last, IntegerMatrix overlap_comb, IntegerVector overlap_loci, 
//                             unsigned int overlap_new_states, unsigned int overlap_num_states) {
//   unsigned int i, m;
//   // unsigned int n_observation = linkage_info.nrow();
//   //decide the first t 
//   unsigned int overlap_len = overlapped.size();
//   // limit the space based on the last limition first
//   int start_overlap = 0;
//   for(i = 0; i < overlap_loci.size(); ++i)
//     if (overlap_loci[i] == overlapped[0]){
//       start_overlap = i;
//       break;
//     }
//     
//     //find the unique states
//     IntegerMatrix comb_last(overlap_new_states, overlap_len);
//     unsigned int count = 0;
//     for (m = 0; m < overlap_num_states; ++m)
//       if(!exclude_last[m]) {
//         IntegerVector tmp = overlap_comb(m, _);
//         comb_last(count++, _) = tmp[Range(start_overlap, start_overlap + overlap_len - 1)];
//       }
//     arma::mat out = unique_rows(as<arma::mat>(comb_last));
//     return(wrap(out));
// }
// 
// IntegerMatrix new_combination(List hmm_info, IntegerVector location, IntegerVector overlapped, IntegerVector exclude_last, IntegerMatrix overlap_comb, 
//                      IntegerVector overlap_loci, IntegerMatrix linkage_info, unsigned int overlap_new_states, unsigned int overlap_num_states) {
//   
//   IntegerMatrix first_comb = unique_overlap(overlapped, exclude_last, overlap_comb, overlap_loci, 
//                                             overlap_new_states, overlap_num_states);
//   //find combinatio of the rest position(include 1 overlap to make sure the linkage)
//   List hidden_states = hmm_info["hidden_states"];
//   IntegerVector pos_possibility = hmm_info["pos_possibility"];
//   IntegerVector undecided_pos = hmm_info["undecided_pos"];
//   unsigned int i, j, m, w, k;
//   int overlap_len = overlapped.size();
//   
//   IntegerVector left_loci = location[Range(overlap_len - 1, location.size() - 1)];
//   int count = 0;
//   int len = left_loci.size();
//   IntegerVector left_possible(len);
//   for(i = 0; i < undecided_pos.size(); ++i) {
//     if (undecided_pos[i] >= left_loci[0] && undecided_pos[i] <= left_loci[len - 1]) {
//       left_possible[count++] = pos_possibility[i];
//     }
//   }
//   
//   IntegerMatrix combination = call_cart_product(left_possible[Range(0, count - 1)]);
//   // get the appeared possiblilities at the overlapped position
//   IntegerVector last_col = first_comb(_, first_comb.ncol() - 1);
//   List first_uni = unique_map(last_col); // start might not from 0 (e.g. ailgnment starts from 2)
//   IntegerVector n1 = first_uni["lengths"];
//   IntegerVector allowed = first_uni["values"];
//   // IntegerVector allowed = unique(last_col);
//   IntegerVector first_col = combination(_, 0);
//   IntegerVector exist = unique(first_col);
//   List exclude_info(2);
//   int flag = 0;
//   int num = 0;
//   IntegerMatrix new_combination(combination.nrow(), combination.ncol());
//   unsigned int start_idx = 0;
//   for(i = 0; i < undecided_pos.size(); ++i)
//     if (undecided_pos[i] == left_loci[0]) {
//       start_idx = i;
//       break;
//     }
//   Rcout <<  "exists: " << exist << "\n";
//   Rcout << "allowed: " << n1 << "\n";
//   Rcout << allowed << "\n";
//   for(m = 0; m < exist.size(); ++m)
//     for(w = 0; w < allowed.size(); ++w)
//       if(exist[m] == allowed[w]) {
//         num++;
//         }
//   if(num != exist.size())
//     flag = 1;
//   num = 0;
//   if(flag) {
//     Rcout << "different unique at same col" << "\n";
//     for(m = 0; m < combination.nrow(); ++m)
//       for(w = 0; w < allowed.size(); ++w)
//         if(allowed(w) == first_col(m)) {
//           new_combination(num++, _) = combination(m, _);
//         }
//     exclude_info = limit_comb_t0(new_combination(Range(0, num - 1), _), hidden_states, left_loci, linkage_info, combination.ncol(), start_idx, num);
//   } else {
//     exclude_info = limit_comb_t0(combination, hidden_states, left_loci, linkage_info, combination.ncol(), start_idx, combination.nrow());
//   }
// 
//   // // Now give the limited combination
//   int num_states = exclude_info["num_states"];
//   IntegerVector exclude = exclude_info["exclude"];
//   Rcout << exclude << "\n";
//   Rcout << "left states:" << num_states << "\n";
//   Rcout << "first states " << first_comb.nrow() << "\n";
//   IntegerMatrix next_comb(num_states, combination.ncol());
//   count = 0;
//   // Now combine first and second part
//   if(flag) {
//     for(m = 0; m < num; ++m) {
//       if(!exclude[m])
//         next_comb(count++, _) = new_combination(m, _);
//     }
//   } else {
//     for(m = 0; m < combination.nrow(); ++m) {
//       if(!exclude[m])
//         next_comb(count++, _) = combination(m, _);
//     }
//   }
//   for(m = 0; m < count; ++m) {
//     for(k = 0; k < next_comb.ncol(); ++k) {
//       Rcout << next_comb(m, k) << "\t";
//     }
//     Rcout << "\n";}
//   IntegerVector new_col = next_comb(_, 0);
//   List second_uni = unique_map(new_col); // start might not from 0 (e.g. ailgnment starts from 2)
//   IntegerVector n2 = second_uni["lengths"];
//   IntegerVector possible = second_uni["values"];
//   int all = 0;
// 
//   for(m = 0; m < n2.size(); ++m)
//     all += n2[m] * n1[m];
//   IntegerMatrix final_comb(all, location.size());
//   Rcout << all << "\n";
//   all = 0;
//   
//   for(k = 0; k < last_col.size(); ++k)
//     for(w = 0; w < count; ++w) {
//       if(new_col(w) == last_col(k)) {
//         for(j = 0; j < overlap_len; ++j) {
//           Rcout << first_comb(k, j) << "\t";
//           final_comb(all, j) = first_comb(k, j);
//         }
//         for(i = 1; i < len; ++i) {
//           Rcout << next_comb(w, i) << "\t";
//           final_comb(all, i + overlap_len - 1) = next_comb(w, i);
//         }
//         all++;
//         Rcout << "\n";
//       }
//     }
//  
//   return(final_comb);
// }
// function to do left join by the last column of first and first of second

// IntegerMatrix left_join(IntegerMatrix first, IntegerMatrix second, int nrow, int ncol) {
//   IntegerMatrix final_comb(nrow, ncol);
//   int all = 0;
//   for(int k = 0; k < first.nrow(); ++k) {
//     for(int w = 0; w < second.nrow(); ++w) {
//       if(second(w, 0) == first(k, first.ncol() - 1)) {
//           for(int j = 0; j < first.ncol(); ++j) {
//           Rcout << first(k, j) << "\t";
//            final_comb(all, j) = first(k, j);
//           }
//           for(int i = 1; i < second.ncol(); ++i) {
//           Rcout << second(w, i) << "\t";
//           final_comb(all, i + first.ncol() - 1) = second(w, i);
//           }
//           all++;
//           Rcout << "\n";
//         }
//     }
//   }
//   return(final_comb);
// }
// 
// 
// IntegerMatrix fill_all_hap(List hidden_states, unsigned int hap_length, IntegerVector n_row) {
//   unsigned int j, k;
//   IntegerMatrix haplotype(NUM_CLASS, hap_length);
//   IntegerVector nuc_j(NUM_CLASS);
//   for(j = 0; j < hap_length; ++j)
//     if(n_row(j) == 1) { // fill the loci with only 1 possibility, can be optimized by using the part has been filled and add/trim
//       nuc_j = hidden_states[j];
//       for(k = 0; k < NUM_CLASS; ++k)
//         haplotype(k, j) = nuc_j[k];  
//     }
//     return haplotype;
// }
// 
// IntegerMatrix make_hap(List hidden_states, IntegerMatrix haplotype, IntegerVector location, unsigned int p_tmax,
//                        IntegerVector combination, unsigned int time_pos, unsigned int num, int hap_min_pos) {
//   unsigned int j, k, idx;
//   
//   for(j = 0; j < num; ++j) {
//     IntegerMatrix hidden = hidden_states[location[j]];
//     //Rcout << hidden << "\n";
//     idx = combination[j];
//     //Rcout << idx << "\n";
//     for(k = 0; k < NUM_CLASS; ++k)
//       haplotype(k, location[j]) = hidden(idx, k);
//   }
//   return(haplotype(_, Range(time_pos - hap_min_pos, time_pos + p_tmax - hap_min_pos - 1)));
// }

// 
// List full_hap_new (List hmm_info, IntegerMatrix linkage_info, unsigned int hap_length, int hap_min_pos) {
//   List hidden_states = hmm_info["hidden_states"];
//   IntegerVector num_states = hmm_info["num_states"];
//   IntegerVector time_pos = hmm_info["time_pos"];
//   IntegerVector p_tmax = hmm_info["p_tmax"];
//   IntegerVector n_row = hmm_info["n_row"];
//   IntegerVector pos_possibility = hmm_info["pos_possibility"];
//   IntegerVector undecided_pos = hmm_info["undecided_pos"];
//   unsigned int t_max = hmm_info["t_max"];
//   unsigned int t, m;
//   List full_hap(t_max);
//   List comb(t_max);
//   IntegerMatrix hap = fill_all_hap(hidden_states, hap_length, n_row);
//   IntegerVector new_num_states(t_max);
//   List overlap_info = get_overlap(p_tmax, time_pos, num_states, undecided_pos, t_max, hap_min_pos);
//   List overlapped = overlap_info["overlapped"];
//   IntegerVector overlapped_id = overlap_info["overlapped_id"];
//   int start_t = overlap_info["start_t"];
//   List loci = overlap_info["location"];
//   
//   //start t info
//   List comb_info_t0 = find_combination(undecided_pos, pos_possibility, p_tmax[start_t], time_pos[start_t], hap_min_pos);
//   IntegerVector location = comb_info_t0["location"];
//   IntegerMatrix combination = comb_info_t0["combination"];
//   unsigned int num = comb_info_t0["num"];
//   List t0 = limit_comb_t0(combination, hidden_states, location, linkage_info, num, 0, num_states[start_t]);
//   IntegerVector exclude = t0["exclude"];
//   IntegerVector exclude_last = exclude;
//   IntegerMatrix comb_in = combination; 
//   // get the states  
//   for(t = 0; t < t_max; ++t) {
//     List full_hap_t;
//       if(num_states[t] != 1) {
//         int count = 0;
//         if(t == start_t) {
//           // decide start_t first
//           new_num_states[t] = t0["num_states"];
//           IntegerMatrix new_comb(new_num_states[t], combination.ncol());
//           full_hap_t = List(new_num_states[t]);
//           for(m = 0; m < num_states[t]; ++m)
//             if(!exclude[m]) {
//               IntegerMatrix haplotype = make_hap(hidden_states, hap, location, p_tmax[t], combination(m, _), time_pos[t], num, hap_min_pos);
//               new_comb(count, _) = combination(m, _);
//               full_hap_t(count++) = haplotype;
//             }
//           comb[t] = new_comb;
//           Rcout << "start done" << "\n";
//         }
//         else {
//           int identical = 0;
//           int last_t = overlapped_id[t];
//           IntegerVector overlapped_t = overlapped[t];
//           IntegerVector loci_lastt = loci[last_t];
//           IntegerVector loci_currt = loci[t];
//           int old_state;
//           Rcout << t << "\t" << last_t << "\n";
//           Rcout << "overlapped: " << overlapped_t << "\n";
//           Rcout << "loci_lastt: " << loci_lastt << "\n";
//           Rcout << "loci_currt: " << loci_currt << "\n";
//           
//           if(loci_lastt[0] <= loci_currt[0] && loci_lastt[loci_lastt.size() - 1] >= loci_currt[loci_currt.size() - 1]) {
//             if(loci_lastt.size() > loci_currt.size()) { // if current is in its overlap
//               // get the unique overlapped combination from the last t
//               if(last_t == start_t) {
//                 exclude_last = IntegerVector(exclude.size());
//                 exclude_last = exclude;
//                 old_state = num_states[last_t];
//                 comb_in = combination;
//               }
//               else {
//                 exclude_last = IntegerVector(new_num_states[last_t]);
//                 for(m = 0; m < new_num_states[last_t]; ++m)
//                   exclude_last[m] = 0;
//                 old_state = new_num_states[last_t];
//                 IntegerMatrix tmp = comb[last_t];
//                 comb_in = IntegerMatrix(tmp.nrow(), tmp.ncol());
//                 comb_in = tmp;
//               }
//               Rcout << "last t:\t" << new_num_states[last_t] << "\told\t" << old_state << "\n";
//               IntegerMatrix new_comb = unique_overlap(overlapped_t, exclude_last, comb_in, loci_lastt, new_num_states[last_t], old_state);
//               new_num_states[t] = new_comb.nrow();
//               comb[t] = new_comb;
//               full_hap_t = List(new_num_states[t]);
//               for(m = 0; m < new_num_states[t]; ++m) {
//                 IntegerMatrix haplotype = make_hap(hidden_states, hap, loci_currt, p_tmax[t], new_comb(m, _), time_pos[t], loci_currt.size(), hap_min_pos);
//                 full_hap_t(count++) = haplotype;
//               }
//             } else if (loci_lastt.size() == loci_currt.size()) {
//               for(m = 0; m < loci_lastt.size(); ++m)
//                 if(loci_lastt[m] != loci_currt[m]) {
//                   identical = 1;
//                   break;
//                 }
//                 if(!identical) {
//                   Rcout << "same sites\n";
//                   comb[t] = comb[last_t];
//                   new_num_states[t] = new_num_states[last_t];
//                   full_hap_t = List(new_num_states[t]);
//                   full_hap_t = full_hap[last_t];
//                 }
//             }
//           } 
//           else {
//             if(last_t == start_t) {
//               exclude_last = IntegerVector(exclude.size());
//               exclude_last = exclude;
//               old_state = num_states[last_t];
//               comb_in = combination;
//             }
//             else {
//               exclude_last = IntegerVector(new_num_states[last_t]);
//               for(m = 0; m < new_num_states[last_t]; ++m)
//                 exclude_last[m] = 0;
//               old_state = new_num_states[last_t];
//               IntegerMatrix tmp = comb[last_t];
//               comb_in = IntegerMatrix(tmp.nrow(), tmp.ncol());
//               comb_in = tmp;
//             }
//             IntegerMatrix new_comb = new_combination(hmm_info, loci_currt, overlapped_t, exclude_last, comb_in,
//                               loci_lastt, linkage_info, new_num_states[last_t], old_state);
//             new_num_states[t] = new_comb.nrow();
//             comb[t] = new_comb;
//             full_hap_t = List(new_num_states[t]);
//             for(m = 0; m < new_num_states[t]; ++m) {
//               IntegerMatrix haplotype = make_hap(hidden_states, hap, loci_currt, p_tmax[t], new_comb(m, _), time_pos[t], loci_currt.size(), hap_min_pos);
//               full_hap_t(count++) = haplotype;
//             }
//           }
//         }
//       } 
//       else {
//         full_hap_t = List(1);
//         full_hap_t[0] = hap(_, Range(time_pos[t] - hap_min_pos, time_pos[t] + p_tmax[t] - hap_min_pos - 1));
//         new_num_states[t] = 1;
//         comb[t] = -1;
//       }
//       Rcout << "new no. states: " << new_num_states[t] << "\n";
//       full_hap[t] = full_hap_t;
//   }
//   
//   List out = List::create(
//     Named("full_hap") = full_hap,
//     Named("new_num_states") = new_num_states, 
//     Named("combination") = comb);
//   
//   return(out);
// }
// List limit_comb_rest(IntegerVector overlapped, IntegerVector exclude_last, IntegerMatrix overlap_comb, IntegerVector overlap_loci, 
//                      IntegerMatrix combination, List hidden_states, IntegerVector location, IntegerMatrix linkage_info, 
//                      unsigned int num, unsigned int start_idx, unsigned int num_states, unsigned int overlap_num_states) {
//   unsigned int i, j, k, idx, m, w;
//   unsigned int n_observation = linkage_info.nrow();
//   //decide the first t 
//   unsigned int overlap_len = overlapped.size();
//   // limit the space based on the last limition first
//   int start_overlap = 0;
//   
//   for(i = 0; i < overlap_loci.size(); ++i) {
//     if (overlap_loci[i] == overlapped[0]){
//       start_overlap = i;
//       break;
//     }
//   }
//   Rcout << start_overlap << "\n";
//   IntegerVector exclude(num_states);
//   for (m = 0; m < overlap_num_states; ++m)
//     if(exclude_last[m]) {
//       IntegerVector comb_last(overlap_len);
//       IntegerVector tmp = overlap_comb(m, _);
//       comb_last = tmp[Range(start_overlap, start_overlap + overlap_len - 1)];
//       Rcout << comb_last << "\n";
//       for(w = 0; w < num_states; ++w) {
//         exclude[w] = 1;
//         IntegerVector comb(overlap_len);
//         IntegerVector tmp2 = combination(w, _);
//         comb = tmp2[Range(0, overlap_len - 1)];
//         Rcout << comb << "\t";
//         // if the overlapped combination is the same, then exclude it as well, so this makes sure there is always one former stage transfer to the next stage
//         for(i = 0; i < overlap_len; ++i)
//           if(comb[i] != comb_last[i]){
//             exclude[w] = 0;
//             break;
//           }
//       }
//     }
//     Rcout << "first stage done" << "\n";
//     // further examine linkage info, if do not meet the linkage, remove as well
//     IntegerMatrix sub_hap(NUM_CLASS, num);
//     IntegerMatrix sub_link = linkage_info(_, Range(start_idx, start_idx + num - 1));
//     int count, linkage_len, all_excluded;
//     linkage_len = num/2; // change the linkage length to be the length appears in the read
//     
//     all_excluded = num_states;
//     while (all_excluded == num_states) {
//       all_excluded = 0;
//       // Rcout << "linkage length " << linkage_len << "\n"; //actual linkage length + 1
//       for (m = 0; m < num_states; ++m) {
//         if(!exclude[m]) {
//           IntegerVector comb = combination(m, _);
//           // Rcout << comb << "\t\t";
//           count = 0;
//           for (k = 0; k < NUM_CLASS; ++k) {
//             for (j = 0; j < num; ++j) {
//               IntegerMatrix hidden = hidden_states[location[j]];
//               idx = comb[j];
//               sub_hap(k, j) = hidden(idx, k);
//             }
//             for (i = 0; i < n_observation; i++) {
//               int flag = 0;
//               for (j = 0; j < num - 1; ++j) {
//                 // if (sub_link(i, j) != -1 && sub_link(i, j) != 4)
//                 if (sub_hap(k, j) == sub_link(i, j) && sub_hap(k, j + 1) == sub_link(i, j + 1))
//                   flag++;
//               }
//               if (flag >= linkage_len) {
//                 count++;
//                 break;
//               }
//             }
//           }
//           if (count != NUM_CLASS) {
//             exclude(m) = 1;
//             all_excluded++;
//           }
//         } else {
//           all_excluded++;
//           }
//       }
//       linkage_len--;
//     }
//     List out = List::create(
//       Named("num_states") = num_states - all_excluded,
//       Named("exclude") = exclude);
//     return(out);
// }


// List trans_permit2(IntegerVector num_states, List hap_info, int t_max, IntegerVector undecided_pos, 
//                   IntegerVector time_pos, IntegerVector p_tmax, int hap_min_pos) {
//   List trans_permits(t_max - 1);
//   unsigned int t, j, m, w, k;
//   int end_t1, begin_t2, end_t2;
//   IntegerVector position(undecided_pos.size());
//   for (t = 0; t < t_max - 1; ++t) {
//     if(num_states[t + 1] != 1) {
//       end_t1 = time_pos[t] + p_tmax[t];
//       begin_t2 = time_pos[t + 1];
//       end_t2 = time_pos[t + 1] + p_tmax[t + 1];
//     
//       if(end_t1 > end_t2)
//         end_t1 = end_t2; // take the smaller one
//       List full_hap_t1 = hap_info(t);
//       List full_hap_t2 = hap_info(t + 1);
//       IntegerMatrix trans(num_states[t], num_states[t + 1]);
//       int count = 0; // number of overlapped variational sites
//       int count1 = 0;
//       
//       for(j = 0; j < undecided_pos.size(); ++j)
//         if (undecided_pos[j] + hap_min_pos < end_t1 && undecided_pos[j] >= time_pos[t]) {
//           count1++;
//           if (undecided_pos[j] + hap_min_pos >= begin_t2)
//             position[count++] = undecided_pos(j) + hap_min_pos;
//         }
//         for(j = 0; j < count; ++j){
//           Rcout << position[j] - begin_t2 << "\t" << position[j] - time_pos[t] << "\n";
//         }
//         for(m = 0; m < num_states[t]; ++m) {
//           IntegerMatrix hap_t1 = full_hap_t1(m);
//           for(w = 0; w < num_states[t + 1]; ++w) {
//             IntegerMatrix hap_t2 = full_hap_t2(w);
//             for(j = 0; j < count; ++j){
//               // Rcout << position[j] - begin_t2 << "\t" << position[j] - time_pos[t] << "\n";
//               for(k = 0; k < NUM_CLASS; ++k)
//                 if (hap_t2(position[j] - begin_t2, k) != hap_t1(position[j] - time_pos[t], k)) {
//                   
//                   trans(m, w) = 1; // represents m cannot transfer to w
//                   break;
//                 }}
//           }
//         }
//         trans_permits(t) = trans;
//     } else {
//       IntegerMatrix temp(num_states[t], num_states[t + 1]);
//       trans_permits(t) = temp;
//     }
//   }
//   return(trans_permits);
// }



