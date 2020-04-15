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

IntegerMatrix linkage_info(List dat_info, IntegerVector undecided_pos) {
  IntegerMatrix ref_index = dat_info["ref_idx"];
  IntegerVector obs = dat_info["nuc"];
  IntegerVector index = dat_info["start_id"];
  IntegerVector ref_pos = dat_info["ref_pos"];
  IntegerVector length = dat_info["length"];
  int hap_min_pos = dat_info["ref_start"];
  int n_observation = dat_info["n_observation"];
  IntegerMatrix link(n_observation, undecided_pos.size());
  
  unsigned int i, j;
  int idx;
  
  for (j = 0; j < undecided_pos.size(); ++j) {
    unsigned int ref_j = undecided_pos[j] + hap_min_pos;
    for (i = 0; i < n_observation; i++) {
      if (ref_pos[index[i]] <= ref_j && ref_j < ref_pos[index[i] + length[i] - 1]) {
        idx = ref_index(i, ref_j - hap_min_pos); // read pos, start from 0
        if (idx != -1)
          link(i, j) = obs[index[i] + idx];
        else 
          link(i, j) = 4; // meaning deletion
      } else 
        link(i, j) = -1; // meaning not covered
    }
  }
  return(link);
}

// [[Rcpp::export]]
IntegerVector limit_comb(IntegerMatrix combination, List hidden_states, IntegerVector location,
                         IntegerMatrix linkage_info, unsigned int num, unsigned int start_idx, unsigned int num_states) {
  unsigned int i, j, k, idx, m;
  unsigned int n_observation = linkage_info.nrow();
  IntegerMatrix sub_hap(NUM_CLASS, num);
  IntegerMatrix sub_link = linkage_info(_, Range(start_idx, start_idx + num - 1));
  IntegerVector exclude(num_states);
  int count, linkage_len, all_excluded;
  linkage_len = num/2;
  
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
          Rcout << sub_hap(k, j) << "\t";
        }
        Rcout << "reads\n";
        for (i = 0; i < n_observation; i++) {
          int flag = 0;
          for (j = 0; j < num - 1; ++j) {
            Rcout << sub_link(i, j) << "\t" << sub_link(i, j + 1);
            // if (sub_link(i, j) != -1 && sub_link(i, j) != 4)
            if (sub_hap(k, j) == sub_link(i, j) && sub_hap(k, j + 1) == sub_link(i, j + 1))
              flag++;
          }
          Rcout << "\n";
          if (flag == linkage_len) {
            count++;
            break;
          }
        }
      }
      if (count != NUM_CLASS) {
        exclude(m) = 1;
        all_excluded++;
        Rcout << m << "th possible combination is excluded\n";
      }
    }
    linkage_len--;
  }
  return(exclude);
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

// [[Rcpp::export]]
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
