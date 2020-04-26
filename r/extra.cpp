#include <RcppArmadillo.h>
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <ctype.h>
#include <limits.h>
#include <unistd.h>
#include <math.h>
#include <iostream>
#include <vector>
#include <algorithm>

using namespace Rcpp;
using namespace std;
// [[Rcpp::depends(RcppArmadillo)]]
#define NUM_CLASS 4
#define MLOGIT_CLASS 4


// [[Rcpp::export]]
// trans_indicator: indicate which state can transfer to which, further_limit indicate some states should not be considered
List prepare_ini_hmm (unsigned int t_max, IntegerVector num_states, List trans_indicator, List further_limit) {
  List trans_new_ind(t_max - 1);
  // List emit(t_max);
  unsigned int w, m, t;
  // 
  // IntegerVector more_limit = further_limit[0];
  // IntegerMatrix trans_ind = trans_indicator[0];
  // IntegerMatrix trans_id_new(num_states[0], num_states[1]);
  // int count = 0;
  // for (w = 0; w < trans_ind.ncol(); ++w)
  //   if(!more_limit[w])
  //     trans_id_new(_, count++) = trans_ind(_, w);
  //   trans_new_ind[0] = trans_id_new;
  
  for(t = 0; t < t_max - 1; ++t) {
    // Rcout << t << "\n";
    IntegerVector more_limit_last = further_limit[t];
    IntegerVector more_limit = further_limit[t + 1];
    IntegerMatrix trans_ind = trans_indicator[t];
    IntegerMatrix trans_ind_new(trans_ind.nrow(), num_states[t + 1]);
    int count2 = 0;
    int count = 0;
    for (w = 0; w < trans_ind.ncol(); ++w)
      if(!more_limit[w])
        trans_ind_new(_, count++) = trans_ind(_, w);
      IntegerMatrix trans_ind_new2(num_states[t], num_states[t + 1]);
      for (m = 0; m < trans_ind.nrow(); ++m)
        if(!more_limit_last[m])
          trans_ind_new2(count2++, _) = trans_ind_new(m, _);
        trans_new_ind[t] = trans_ind_new2;
  }
  
  return(trans_new_ind);
}
  
  // [[Rcpp::export]]
  List trans_permit2(IntegerVector num_states, List combination, int t_max, IntegerVector undecided_pos, 
                    IntegerVector time_pos, IntegerVector p_tmax, int hap_min_pos) {
    List trans_permits(t_max - 1);
    List further_limit_col(t_max - 1);
    List further_limit_row(t_max - 1);
    IntegerVector new_num_states(t_max);
    unsigned int t, j, m, w;
    int end_t1, begin_t2, end_t2;
    for (t = 0; t < t_max; ++t)
      new_num_states[t] = num_states[t];
    
    for (t = 0; t < t_max - 1; ++t) {
      if(num_states[t + 1] != 1 && num_states[t] != 1) { // assume if no. state = 1, then it can transfer to all
        end_t1 = time_pos[t] + p_tmax[t];
        begin_t2 = time_pos[t + 1];
        end_t2 = time_pos[t + 1] + p_tmax[t + 1];
        if(end_t1 > end_t2)
          end_t1 = end_t2; // take the smaller one
        IntegerMatrix full_hap_t1 = combination(t);
        IntegerMatrix full_hap_t2 = combination(t + 1);
        IntegerMatrix trans(num_states[t], num_states[t + 1]);
        int count = 0; // number of overlapped variational sites
        int count1 = 0;
        
        for(j = 0; j < undecided_pos.size(); ++j)
          if (undecided_pos[j] + hap_min_pos < end_t1 && undecided_pos[j] >= time_pos[t]) {
            count1++;
            if(undecided_pos[j] + hap_min_pos >= begin_t2)
              count++;
          }
          // Rcout << t << "\t" << count << "\n";
        int left = count1 - count;
          for(m = 0; m < num_states[t]; ++m) {
            IntegerVector hap_t1 = full_hap_t1(m, _);
            // Rcout << hap_t1 << "\n";
            for(w = 0; w < num_states[t + 1]; ++w) {
              IntegerVector hap_t2 = full_hap_t2(w, _);
              // Rcout << hap_t2 << "\t";
              for(j = 0; j < count; ++j)
                if (hap_t2(j) != hap_t1(j + left)) {
                  trans(m, w) = 1; // represents m cannot transfer to w
                  break;
                }
            }
          }
          trans_permits(t) = trans;
      } else {
        IntegerMatrix temp(num_states[t], num_states[t + 1]);
        trans_permits(t) = temp;
      }
    }
    // label the states that should be excluded (go forward)
    for(t = 0; t < t_max - 1; ++t) {
      if(num_states[t + 1] != 1 || num_states[t] != 1) { 
      // if for a row, all == 1 or for a column all == 1, then these two should be excluded
      IntegerVector more_limits_row(num_states[t]);
      IntegerVector more_limits_col(num_states[t + 1]);
      IntegerVector more_limits_col_last;
     
      int exluded_col = 0;
      for(w = 0; w < num_states[t + 1]; ++w) {
        int num_col = 0;
        for(m = 0; m < num_states[t]; ++m) {
          if(trans(m, w) == 1)
            num_col++;
        }
        if(num_col == num_states[t]){
          more_limits_col[w] = 1; // meaning no state transfer to it
          exluded_col++;
        }
      }
      for(m = 0; m < num_states[t]; ++m) {
        int num_row = 0;
        for(w = 0; w < num_states[t + 1]; ++w) {
          if(trans(m, w) == 1 && more_limits_col[w] != 1)
            num_row++;
        }
        if(num_row == num_states[t + 1] - exluded_col){
          more_limits_row[m] = 1; // meaning it transfer to no state, considering the excluded column 
        }
      }
      // found the states that never appears
      further_limit_col(t) = more_limits_col;
      further_limit_row(t) = more_limits_row;
    } else {
      IntegerVector more_limits_col(num_states[t + 1]);
      IntegerVector more_limits_row(num_states[t]);
      further_limit_col(t) = more_limits_col;
      further_limit_row(t) = more_limits_row;
    }
    }
    // go backward
    for(t = t_max - 1; t <= 0; --t) {
      if(num_states[t + 1] != 1 || num_states[t] != 1) { 
        // if for a row, all == 1 or for a column all == 1, then these two should be excluded
        IntegerVector more_limits_row = further_limit_row(t);
        IntegerVector more_limits_col = further_limit_col[t - 1];
        
        int index;
        for(w = 0; w < more_limits_col.size(); ++w) {
          if(more_limits_row[w] != more_limits_col[w])
            index = 0;
        }
       
        // found the states that never appears
        further_limit_col(t) = more_limits_col;
        further_limit_row(t) = more_limits_row;
      } 
    }
    // get new number of states
    List further_limit(t_max);
    IntegerVector more_limits_row = further_limit_row[0];
    for (m = 0; m < num_states[0]; ++m)
      if (more_limits_row[m] == 1)
        new_num_states[0]--;
    further_limit[0] = further_limit_row[0];
    IntegerVector more_limits_col = further_limit_col[t_max - 2];
    for(w = 0; w < num_states[t_max - 1]; ++w) {
        if(more_limits_col[w] == 1)
          new_num_states[t_max - 1]--;
    }
    further_limit[t_max - 1] = further_limit_col[t_max - 2];
    
    for(t = 1; t < t_max - 1; ++t) {
      more_limits_row = further_limit_row(t);
      more_limits_col = further_limit_col(t - 1);
      IntegerVector combine(more_limits_row.size());
      // find the union of these two vector
      for(w = 0; w < num_states[t]; ++w) {
        if(more_limits_col[w] == 1 || more_limits_row[w] == 1)
          combine[w] = 1;
      }
      further_limit[t] = combine; // indicates some states cannot appear
      for(w = 0; w < num_states[t]; ++w)
        if(combine[w] == 1)
          new_num_states[t]--;
    }
    
    List out = List::create(
      Named("trans_permits") = trans_permits,
      Named("further_limit") = further_limit, 
      Named("new_num_states") = new_num_states);
    
    return(out);
  }

vector<vector<int> > cart_product (const vector<vector<int> > &v) {
  vector<vector<int> > s = {{}};
  for (const auto& u : v) {
    vector<vector<int> > r;
    for (const auto& x : s) {
      for (const auto y : u) {
        r.push_back(x);
        r.back().push_back(y);
      }
    }
    s = move(r);
  }
  return s;
}

IntegerMatrix call_cart_product(IntegerVector len) {
  unsigned int row = len.size();
  vector<vector<int> > vec(row);
  unsigned int col, count, i, j;
  for (i = 0; i < row; i++) {
    count = 1;
    col = len[i];
    vec[i] = vector<int>(col);
    for (j = 0; j < col; j++)
      vec[i][j] = count++;
  }
  vector<vector<int> > res = cart_product(vec);
  IntegerMatrix out(res.size(), row);
  for(i = 0; i < res.size(); ++i)
    for(j = 0; j < row; ++j) 
      out(i, j) = res[i][j] - 1; //minus 1 for the index in C
  
  return(out);
}

List find_combination(IntegerVector undecided_pos, IntegerVector pos_possibility, 
                      unsigned int p_tmax, unsigned int time_pos, int hap_min_pos) {
  //possible combination of the rest non-unique loci
  IntegerVector location(pos_possibility.size());
  IntegerVector location_len(pos_possibility.size());
  unsigned int num = 0;
  for(unsigned int i = 0; i < pos_possibility.size(); ++i)
    if(time_pos - hap_min_pos <= undecided_pos[i] && undecided_pos[i] < time_pos + p_tmax - hap_min_pos) {
      location(num) = undecided_pos[i];
      location_len(num++) = pos_possibility[i];
    }
    IntegerMatrix combination = call_cart_product(location_len[Range(0, num - 1)]);
    List ls = List::create(
      Named("combination") = combination, // possible comb
      Named("num") = num, // number of comb sites at time t
      Named("location") = location[Range(0, num - 1)]); // undecided site at time t [here assume alignment starts from 0]
    return(ls);
}

List limit_comb(IntegerMatrix combination, List hidden_states, IntegerVector location,
                IntegerMatrix linkage_info, unsigned int num, unsigned int start_idx, unsigned int num_states) {
  unsigned int i, j, k, idx, m;
  unsigned int n_observation = linkage_info.nrow();
  IntegerMatrix sub_hap(NUM_CLASS, num);
  IntegerMatrix sub_link = linkage_info(_, Range(start_idx, start_idx + num - 1));
  IntegerVector exclude(num_states);
  int count, linkage_len, all_excluded;
  linkage_len = num/2;
  
  for(i = 0 ; i < n_observation; ++i)
    for (j = 0; j < num; ++j)
     Rcout << sub_link(i, j) << "\t";
  Rcout << "\n";
  all_excluded = num_states;
  while (all_excluded == num_states) {
    all_excluded = 0;
    Rcout << "linkage length " << linkage_len << "\n"; //actual linkage length + 1
    for (m = 0; m < num_states; ++m) {
      exclude(m) = 0;
      IntegerVector comb = combination(m, _);
      count = 0;
      for (k = 0; k < NUM_CLASS; ++k) {
        for (j = 0; j < num; ++j) {
         
          IntegerMatrix hidden = hidden_states[location[j]];
          idx = comb[j];
          sub_hap(k, j) = hidden(idx, k);
          Rcout << sub_hap(k, j)  << "\t";
        }
        Rcout << "\n";
        for (i = 0; i < n_observation; i++) {
          int flag = 0;
          for (j = 0; j < num - 1; ++j) {
            // if (sub_link(i, j) != -1 && sub_link(i, j) != 4)
            if (sub_hap(k, j) == sub_link(i, j) && sub_hap(k, j + 1) == sub_link(i, j + 1))
              flag++;
          }
          if (flag >= linkage_len) {
            count++;
            break;
          }
        }
      }
      if (count != NUM_CLASS) {
        exclude(m) = 1;
        all_excluded++;
      }
    }
    linkage_len--;
  }
  List out = List::create(
    Named("num_states") = num_states - all_excluded,
    Named("exclude") = exclude);
  return(out);
}
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



