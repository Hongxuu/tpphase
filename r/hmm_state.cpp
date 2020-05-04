#include <RcppArmadillo.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "utils.h"
#include "hmm_state.h"
using namespace Rcpp;
using namespace std;
#define NUM_CLASS 4
#define LINKAGE_LEN 2

// [[Rcpp::depends(RcppArmadillo)]]
class Permutation {
public:
  vector<vector<int> > permuteUnique(vector<int>& nums) {
    vector<vector<int> > res;
    sort(nums.begin(), nums.end());
    res.push_back(nums);
    while (next_permutation(nums.begin(), nums.end())) {
      res.push_back(nums);
    }
    return res;
  }
};

IntegerMatrix Twopossible(IntegerVector a) {
  IntegerVector temp = {a[0], a[1], a[0], a[1], a[1], a[0], a[1], a[0]};
  temp.attr("dim") = Dimension(2, 4);
  return(as<IntegerMatrix>(temp));
}

IntegerMatrix Fourpossible(IntegerVector small, int big) {
  IntegerVector temp = {small[0], small[1], big, big, small[1], small[0], big, big,
                        big, big, small[0], small[1], big, big, small[1], small[0]};
  temp.attr("dim") = Dimension(4, 4);
  return(as<IntegerMatrix>(temp));
}

IntegerMatrix call_permute(vector<int> a) {
  Permutation per;
  vector<vector<int> > res =  per.permuteUnique(a);
  IntegerMatrix permutation(res.size(), res[0].size());
  for (int i = 0; i < res.size(); i++)
    for (int j = 0; j < res[i].size(); j++) 
      permutation(i, j) = res[i][j];
  
  return(permutation);
}

IntegerMatrix call_permute_N(vector<int> a, unsigned int genome_A) {
  Permutation per;
  vector<vector<int> > res =  per.permuteUnique(a);
  IntegerMatrix permutation(res.size(), res[0].size());
  for (int i = 0; i < res.size(); i++)
    for (int j = 0; j < res[i].size(); j++) 
      permutation(i, j) = res[i][j];
  IntegerVector del(res.size() * 2, -1);
  del.attr("dim") = Dimension(res.size(), 2);
  arma::mat m1 = as<arma::mat>(permutation);
  arma::mat m2 = as<arma::mat>(del);
  arma::mat out;
  if(genome_A)
    out = join_rows(m2, m1);
  else
    out = join_rows(m1, m2);
  
  return(wrap(out));
}

List aux_noN_S2(IntegerVector sum_site, IntegerVector hap_site) {
  
  int n_row;
  double sum;
  IntegerMatrix temp;
  IntegerMatrix temp2;
  
  sum = sum_site[1] + sum_site[0];
  if (sum_site[0]/sum <= 0.38) { //1st one appears 3 times
    temp = call_permute({hap_site[0], hap_site[1], hap_site[1], hap_site[1]});
    temp2 = Twopossible(hap_site);
    arma::mat m1 = as<arma::mat>(temp);
    arma::mat m2 = as<arma::mat>(temp2);
    arma::mat out = join_cols(m1, m2);
    temp = wrap(out);
    n_row = 6;
  } else if (sum_site[0]/sum >= 0.62) {//2nd one appears 3 times
    temp = call_permute({hap_site[0], hap_site[0], hap_site[0], hap_site[1]});
    temp2 = Twopossible(hap_site);
    arma::mat m1 = as<arma::mat>(temp);
    arma::mat m2 = as<arma::mat>(temp2);
    arma::mat out = join_cols(m1, m2);
    temp = wrap(out);
    n_row = 6;
  } else {
    temp = Twopossible(hap_site);
    n_row = 2;
  }
  
  List ls = List::create(
    Named("n_row") = n_row,
    Named("temp") = temp);
  return ls; 
}

List aux_noN_S3(IntegerVector sum_site, IntegerVector hap_site) {
  int n_row;
  double sum;
  IntegerMatrix temp;
  List from_2(2);
  
  sum = sum_site[2] + sum_site[1] + sum_site[0];
  double p0 = sum_site[0]/sum;
  double p1 = sum_site[1]/sum;
  double p2 = sum_site[2]/sum;
  IntegerVector small(2);
  IntegerVector small_count(2);
  // different ceriterion to reduce the possibilities, basically rule out things like ATAG
  if (((sum_site[2] + sum_site[1])/sum_site[0] <= 1.2 && (sum_site[2] + sum_site[1])/sum_site[0] >= 0.8) ||
      (sum_site[0]/sum >= 0.45)) {
    small = {hap_site[1], hap_site[2]};
    temp = Fourpossible(small, hap_site[0]);
    n_row = 4;
  } else if (((sum_site[1] + sum_site[0])/sum_site[2] <= 1.2 && (sum_site[0] + sum_site[1])/sum_site[2] >= 0.8) ||
    (sum_site[2]/sum >= 0.45)) {
    small = {hap_site[1], hap_site[0]};
    temp = Fourpossible(small, hap_site[2]);
    n_row = 4;
  } else if (((sum_site[2] + sum_site[0])/sum_site[1] <= 1.2 && (sum_site[0] + sum_site[2])/sum_site[1] >= 0.8) ||
    (sum_site[1]/sum >= 0.45)) {
    small = {hap_site[2], hap_site[0]};
    temp = Fourpossible(small, hap_site[1]);
    n_row = 4;
  } else if (p0 <= 0.25) {
    small = {hap_site[1], hap_site[2]};
    small_count = {sum_site[1], sum_site[2]};
    List from_2 = aux_noN_S2(small_count, small);
    IntegerMatrix ttemp = from_2["temp"];
    temp = ttemp;
    n_row = from_2["n_row"];
  } else if (p1 <= 0.25) {
    small = {hap_site[0], hap_site[2]};
    small_count = {sum_site[0], sum_site[2]};
    List from_2 = aux_noN_S2(small_count, small);
    IntegerMatrix ttemp = from_2["temp"];
    temp = ttemp;
    n_row = from_2["n_row"];
  } else if (p2 <= 0.25) {
    small = {hap_site[0], hap_site[1]};
    small_count = {sum_site[0], sum_site[1]};
    List from_2 = aux_noN_S2(small_count, small);
    IntegerMatrix ttemp = from_2["temp"];
    temp = ttemp;
    n_row = from_2["n_row"];
  } else { //hopefully this won't happen, when three roughly equal happens
    small = {hap_site[1], hap_site[2]};
    temp = Fourpossible(small, hap_site[0]);
    arma::mat m1 = as<arma::mat>(temp);
    small = {hap_site[1], hap_site[0]};
    temp = Fourpossible(small, hap_site[2]);
    arma::mat m2 = as<arma::mat>(temp);
    small = {hap_site[2], hap_site[0]};
    temp = Fourpossible(small, hap_site[1]);
    arma::mat m3 = as<arma::mat>(temp);
    arma::mat out = join_cols(m1, m2, m3);
    temp = wrap(out);
    n_row = 12;
  }
  // also include the situation of 2 possibles
  int min_id = which_min(sum_site);
  IntegerVector sum_site2(2);
  IntegerVector hap_site2(2);
  int num = 0;
  for(unsigned int i = 0; i < 3; ++i) {
    if(i != min_id) {
      sum_site2[num] = sum_site[i];
      hap_site2[num] = hap_site[i];
      num++;
    }
  }
  List list_2pos = aux_noN_S2(sum_site2, hap_site2);
  IntegerMatrix temp2 = list_2pos["temp"];
  int more_row = list_2pos["n_row"];
  n_row += more_row;
  arma::mat m1 = as<arma::mat>(temp);
  arma::mat m2 = as<arma::mat>(temp2);
  arma::mat out = join_cols(m1, m2);
  temp = wrap(out);
  
  List ls = List::create(
    Named("n_row") = n_row,
    Named("temp") = temp);
  return ls; 
}
/*
 * determin the number of hidden states site by site
 */
List sbs_state(unsigned int num, unsigned int ref_j, IntegerVector hap_site, IntegerVector sum_site, 
               CharacterVector uni_alignment) {
  unsigned int l, m;
  double sum;
  unsigned int n_row;
  List haplotype(1);
  //record possible hidden states[only suitable for ployploids]
  if (num == 1) {
    IntegerVector temp = {hap_site[0], hap_site[0], hap_site[0], hap_site[0]};
    haplotype[0] = temp;
    n_row = 1;
  }
  else if (num == 2) {
    if(hap_site[0] == -1) {
      IntegerVector temp(NUM_CLASS);
      if(uni_alignment[ref_j] == "I")
        temp = {-1, -1, hap_site[1], hap_site[1]};
      else if (uni_alignment[ref_j] == "J")
        temp = {hap_site[1], hap_site[1], -1, -1};
      else
        temp = {hap_site[1], hap_site[1], hap_site[1], hap_site[1]};
      haplotype[0] = temp;
      n_row = 1;
    } else { // if N not appears
      List out = aux_noN_S2(sum_site, hap_site);
      haplotype[0] = out["temp"];
      n_row = out["n_row"];
      // undecided_pos(more_than1) = j; // relative to the haplotype currently try to infer (alignment start from 0)
      // pos_possibility(more_than1++) = n_row;
    }
  }
  else if (num == 3) {
    sum = sum_site[2] + sum_site[1] + sum_site[0];
    if (hap_site[0] == -1) {
      // if N appears, 4 possibilities
      if(sum_site[0]/sum >= 0.45) {
        IntegerVector temp(2 * NUM_CLASS);
        if (uni_alignment[ref_j] == "I")
          temp = {-1, -1,  -1, -1, hap_site[1], hap_site[2], //J in universal alignment
                  hap_site[2], hap_site[1]};
        else if(uni_alignment[ref_j] == "J")
          temp = {hap_site[1], hap_site[2], hap_site[2],hap_site[1],  //I in universal alignment
                  -1,  -1,  -1, -1};
        temp.attr("dim") = Dimension(2, 4); // by column
        haplotype[0] = temp;
        n_row = 2;
      } else {
        IntegerVector sum_site_cleaned = {sum_site[1], sum_site[2]};
        IntegerVector hap_site_cleaned = {hap_site[1], hap_site[2]};
        List out = aux_noN_S2(sum_site_cleaned, hap_site_cleaned);
        haplotype[0] = out["temp"];
        n_row = out["n_row"];
      }
    } else {
      List out = aux_noN_S3(sum_site, hap_site);
      haplotype[0] = out["temp"];
      n_row = out["n_row"];
    }
    // undecided_pos(more_than1) = j;
    // pos_possibility(more_than1++) = n_row[j];
  }
  else if(num == 4) {
    sum = sum_site[3] + sum_site[2] + sum_site[1] + sum_site[0];
    if(hap_site[0] == -1) {
      if(sum_site[0]/sum >= 0.45) {
        arma::mat temp;
        IntegerMatrix inner_tmp(2, NUM_CLASS);
        for(l = 0; l < num; ++l)
          for(m = l + 1; m < num; ++m) {
            if(uni_alignment[ref_j] == "I") //deletion in B genome
              inner_tmp = call_permute_N({hap_site[l], hap_site[m]}, 0);
            else if(uni_alignment[ref_j] == "J")
              inner_tmp = call_permute_N({hap_site[l], hap_site[m]}, 1);
            temp = join_cols(temp, as<arma::mat>(inner_tmp));
          }
          haplotype[0] = wrap(temp);
        n_row = 2;
      } else {
        IntegerVector sum_site_cleaned = {sum_site[1], sum_site[2], sum_site[3]};
        IntegerVector hap_site_cleaned = {hap_site[1], hap_site[2], hap_site[3]};
        List out = aux_noN_S3(sum_site_cleaned, hap_site_cleaned);
        haplotype[0] = out["temp"];
        n_row = out["n_row"];
      }
    } else {
      IntegerMatrix temp = call_permute({hap_site[0], hap_site[1], hap_site[2], hap_site[3]});
      haplotype[0] = temp;
      n_row = 24;
    }
    // undecided_pos(more_than1) = j;
    // pos_possibility(more_than1++) = n_row[j];
  }
  List ls = List::create(
    Named("n_row") = n_row,
    Named("haplotype") = haplotype);
  
  return(ls);
}

// limit combination based on linkage information
List limit_comb_t0(IntegerMatrix combination, List hidden_states, IntegerVector location,
                   IntegerMatrix linkage_info, unsigned int num, unsigned int start_idx, unsigned int num_states) {
  unsigned int i, j, k, idx, m;
  unsigned int n_observation = linkage_info.nrow();
  IntegerMatrix sub_hap(NUM_CLASS, num);
  IntegerMatrix sub_link = linkage_info(_, Range(start_idx, start_idx + num - 1));
  IntegerVector exclude(num_states);
  int count, linkage_len, all_excluded;
  linkage_len = num/2; // change the linkage length to be the length appears in the read
  int cut_off;
  all_excluded = num_states;
  while (all_excluded == num_states) {
    cut_off = NUM_CLASS;
    // Rcout << "linkage length " << linkage_len << "\n"; //actual linkage length + 1
    while (cut_off >= 1 && all_excluded == num_states) {
      all_excluded = 0;
      for (m = 0; m < num_states; ++m) {
        exclude(m) = 0;
        IntegerVector comb = combination(m, _);
        // Rcout << comb << "\t\t";
        count = 0;
        for (k = 0; k < NUM_CLASS; ++k) {
          // Rcout << "k" << k << "\n";
          for (j = 0; j < num; ++j) {
            IntegerMatrix hidden = hidden_states[location[j]];
            idx = comb[j];
            sub_hap(k, j) = hidden(idx, k);
            // Rcout << sub_hap(k, j) << "\t";
          }
          // Rcout << "read" << "\n";
          for (i = 0; i < n_observation; i++) {
            int flag = 0;
            for (j = 0; j < num - 1; ++j) {
              // Rcout << sub_link(i, j) << "\t" << sub_link(i, j + 1);
              // if (sub_link(i, j) != -1 && sub_link(i, j) != 4)
              if (sub_hap(k, j) == sub_link(i, j) && sub_hap(k, j + 1) == sub_link(i, j + 1))
                flag++;
            }
            // Rcout << "\n" << flag << "\n" ;
            if (flag >= linkage_len) {
              count++;
              break;
            }
          }
        }
        // Rcout << "count "<< count << "\n";
        if (count != cut_off) {
          exclude(m) = 1;
          all_excluded++;
        }
      }
      cut_off--;
    }
    linkage_len--;
  }
  // Rcout << exclude << "\n";
  List out = List::create(
    Named("num_states") = num_states - all_excluded,
    Named("exclude") = exclude);
  return(out);
}
 // trans_indicator: indicate which state can transfer to which, further_limit indicate some states should not be considered
 List prepare_ini_hmm (unsigned int t_max, IntegerVector num_states, List trans_indicator, List further_limit) {
   List trans_new_ind(t_max - 1);
   unsigned int w, m, t;
   
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

template <typename T>
inline bool approx_equal_cpp(const T& lhs, const T& rhs, double tol = 0.00000001) {
  return arma::approx_equal(lhs, rhs, "absdiff", tol);
}

arma::mat unique_rows(const arma::mat& m) {
  arma::uvec ulmt = arma::zeros<arma::uvec>(m.n_rows);
  for (arma::uword i = 0; i < m.n_rows; i++)
    for (arma::uword j = i + 1; j < m.n_rows; j++)
      if (approx_equal_cpp(m.row(i), m.row(j))) { ulmt(j) = 1; break; }
  
  return m.rows(find(ulmt == 0));
}

// get the unique rows for overlapped region
IntegerMatrix unique_overlap(IntegerVector overlapped, IntegerVector exclude_last, IntegerMatrix overlap_comb, IntegerVector overlap_loci, 
                             unsigned int overlap_new_states, unsigned int overlap_num_states) {
  unsigned int i, m;
  // unsigned int n_observation = linkage_info.nrow();
  //decide the first t 
  unsigned int overlap_len = overlapped.size();
  // limit the space based on the last limition first
  int start_overlap = 0;
  for(i = 0; i < overlap_loci.size(); ++i)
    if (overlap_loci[i] == overlapped[0]) {
      start_overlap = i;
      break;
    }
    
    //find the unique states
    IntegerMatrix comb_last(overlap_new_states, overlap_len);
    unsigned int count = 0;
    for (m = 0; m < overlap_num_states; ++m)
      if(!exclude_last[m]) {
        IntegerVector tmp = overlap_comb(m, _);
        comb_last(count++, _) = tmp[Range(start_overlap, start_overlap + overlap_len - 1)];
      }
    arma::mat out = unique_rows(as<arma::mat>(comb_last));
    return(wrap(out));
}

List get_overlap(IntegerVector p_tmax, IntegerVector time_pos, IntegerVector num_states,
                 IntegerVector undecided_pos, unsigned int t_max, int hap_min_pos)
{
  unsigned int start_t = 0;
  // get the first t which has variation
  for (unsigned int t = 0; t < t_max; ++t)
    if(num_states[t] > 1) {
      start_t = t;
      break;
    }
    List overlapped(t_max);
    List location(t_max);
    IntegerVector overlapped_idx(t_max);
    unsigned int begin, end, end1, min;
    
    for (unsigned int t = 0; t < t_max; ++t) {
      
      if (num_states[t] == 1) {
        overlapped[t] = -1;
        overlapped_idx[t] = -1;
        location[t] = -1;
        continue;
      }
      if(t == start_t) {
        overlapped[t] = -1;
        overlapped_idx[t] = -1;
        begin = time_pos[start_t] - hap_min_pos;
        end = time_pos[start_t] + p_tmax[start_t] - hap_min_pos;
        int num = 0;
        IntegerVector location_t(undecided_pos.size());
        for (unsigned int m = 0; m < undecided_pos.size(); ++m)
          if (undecided_pos[m] >= begin && undecided_pos[m] < end)
            location_t(num++) = undecided_pos[m];
          location[start_t] = location_t[Range(0, num - 1)];
          continue;
      }
      begin = time_pos[t] - hap_min_pos;
      end = time_pos[t] + p_tmax[t] - hap_min_pos;
      // store location
      int num = 0;
      IntegerVector location_t(undecided_pos.size());
      for (unsigned int m = 0; m < undecided_pos.size(); ++m)
        if (undecided_pos[m] >= begin && undecided_pos[m] < end)
          location_t(num++) = undecided_pos[m];
        location[t] = location_t[Range(0, num - 1)];
        int len = 0;
        int id_t = 0;
        int index = 0;
        // Rcout << "location" << location_t[num - 1] << "\n";
        for (unsigned int t1 = 0; t1 < t; ++t1) {
          IntegerVector last_location = location[t1];
          // begin1 = time_pos[t1] - hap_min_pos;
          end1 = last_location[last_location.size() - 1];
          // Rcout << end1 << "\t";
          // end1 = time_pos[t1] + p_tmax[t1] - hap_min_pos;
          if(begin <= end1) {
            min = end1;
            // minimum overlapped region
            if(location_t[num - 1] < end1)
              min = location_t[num - 1];
            // find the time t which has the longest coverage
            if(min - location_t[0] > len) {
              // Rcout << t << " has coverage with " << t1 << "\n";
              len = min - location_t[0];
              id_t = min;
              index = t1;
            }
          }
        }
        // Rcout<< "\n" << begin << "\t" << id_t << "\n";
        int count = 0;
        IntegerVector position(undecided_pos.size());
        for (unsigned int m = 0; m < undecided_pos.size(); ++m) 
          if (undecided_pos[m] >= begin && undecided_pos[m] <= id_t)
            position(count++) = undecided_pos[m];
          // Rcout << position << "\n";
          if(count) { 
            overlapped[t] = position[Range(0, count - 1)];
            overlapped_idx[t] = index;
          }
    }
    
    List overlap = List::create(
      Named("location") = location,
      Named("overlapped") = overlapped,
      Named("overlapped_id") = overlapped_idx, 
      Named("start_t") = start_t);
    return(overlap);
}

IntegerMatrix new_combination(List hmm_info, IntegerVector location, IntegerVector overlapped, IntegerVector exclude_last, IntegerMatrix overlap_comb, 
                              IntegerVector overlap_loci, IntegerMatrix linkage_info, unsigned int overlap_new_states, unsigned int overlap_num_states) {
  
  IntegerMatrix first_comb = unique_overlap(overlapped, exclude_last, overlap_comb, overlap_loci, 
                                            overlap_new_states, overlap_num_states);
  //find combinatio of the rest position(include 1 overlap to make sure the linkage)
  List hidden_states = hmm_info["hidden_states"];
  IntegerVector pos_possibility = hmm_info["pos_possibility"];
  IntegerVector undecided_pos = hmm_info["undecided_pos"];
  unsigned int i, j, m, w, k;
  int overlap_len = overlapped.size();
  
  IntegerVector left_loci = location[Range(overlap_len - 1, location.size() - 1)];
  int count = 0;
  int len = left_loci.size();
  IntegerVector left_possible(len);
  for(i = 0; i < undecided_pos.size(); ++i) {
    if (undecided_pos[i] >= left_loci[0] && undecided_pos[i] <= left_loci[len - 1]) {
      left_possible[count++] = pos_possibility[i];
    }
  }
  
  IntegerMatrix combination = call_cart_product(left_possible[Range(0, count - 1)]);
  // get the appeared possiblilities at the overlapped position
  IntegerVector last_col = first_comb(_, first_comb.ncol() - 1);
  List first_uni = unique_map(last_col); // start might not from 0 (e.g. ailgnment starts from 2)
  IntegerVector n1 = first_uni["lengths"];
  IntegerVector allowed = first_uni["values"];
  // IntegerVector allowed = unique(last_col);
  IntegerVector first_col = combination(_, 0);
  IntegerVector exist = unique(first_col);
  List exclude_info(2);
  int flag = 0;
  int num = 0;
  IntegerMatrix new_combination(combination.nrow(), combination.ncol());
  unsigned int start_idx = 0;
  for(i = 0; i < undecided_pos.size(); ++i)
    if (undecided_pos[i] == left_loci[0]) {
      start_idx = i;
      break;
    }
    
    // Rcout << allowed << "\n";
  for(m = 0; m < exist.size(); ++m)
    for(w = 0; w < allowed.size(); ++w)
      if(exist[m] != allowed[w]) {
        flag = 1;
        break;
      }
      
  if(flag) {
    for(m = 0; m < combination.nrow(); ++m)
      for(w = 0; w < allowed.size(); ++w)
        if(allowed(w) == combination(m, 0))
          new_combination(num++, _) = combination(m, _);
    exclude_info = limit_comb_t0(new_combination(Range(0, num - 1), _), hidden_states, left_loci, linkage_info, combination.ncol(), start_idx, num);
  } else {
    exclude_info = limit_comb_t0(combination, hidden_states, left_loci, linkage_info, combination.ncol(), start_idx, combination.nrow());
  }
        
        // // Now give the limited combination
  int num_states = exclude_info["num_states"];
  IntegerVector exclude = exclude_info["exclude"];
  // Rcout << num_states << "\n";
  
  IntegerMatrix next_comb(num_states, combination.ncol());
  count = 0;
        // Now combine first and second part
  if(flag) {
    for(m = 0; m < num; ++m) {
      if(!exclude[m])
          next_comb(count++, _) = new_combination(m, _);
      }
  } else {
      for(m = 0; m < combination.nrow(); ++m) {
        if(!exclude[m])
          next_comb(count++, _) = combination(m, _);
      }
  }
        
  IntegerVector new_col = next_comb(_, 0);
  List second_uni = unique_map(new_col); // start might not from 0 (e.g. ailgnment starts from 2)
  IntegerVector n2 = second_uni["lengths"];
  IntegerVector possible = second_uni["values"];
  int all = 0;
  
  for(m = 0; m < n2.size(); ++m)
    all += n2[m] * n1[m];
  IntegerMatrix final_comb(all, location.size());
  
  all = 0;
        
  for(k = 0; k < last_col.size(); ++k)
    for(w = 0; w < count; ++w) {
      if(new_col(w) == last_col(k)) {
        for(j = 0; j < overlap_len; ++j)
          final_comb(all, j) = first_comb(k, j);
        for(i = len - 1; i < len; ++i)
          final_comb(all++, i + overlap_len - len + 1) = next_comb(w, i);
      }
    }
        
  return(final_comb);
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
/*
 * determine the hidden state by each time t, only use this on the first time t, then, the rest non_overlapped still use site by site
 */
// List tbt_state(unsigned int num, unsigned int ref_j, IntegerVector hap_site, IntegerVector sum_site, 
//                List hmm_info, List dat_info) {
//   
//   List n_in_t = hmm_info["n_in_t"];
//   IntegerVector n_t = hmm_info["n_t"];
//   IntegerVector length = dat_info["length"];
//   IntegerVector obs = dat_info["nuc"];
//   IntegerVector obs_index = dat_info["id"];
//   IntegerVector time_pos = hmm_info["time_pos"];
//   IntegerVector p_tmax = hmm_info["p_tmax"];
//   IntegerVector index = dat_info["start_id"];
//   int hap_max_pos = dat_info["ref_length_max"];
//   int hap_min_pos = dat_info["ref_start"];
//   int hap_length = hap_max_pos - hap_min_pos;
//   unsigned int t_max = hmm_info["t_max"];
//   unsigned int qua_in, obs_in;
//   // hash the reads at the first time position
//   std::string outfile = "./out.fasta";
//   FILE *fp = fopen(outfile.c_str(), "w");
//   if (!fp)
//     stop("Cannot open file");
//   
//   for (unsigned int i = 0; i < n_t[0]; ++i) {
//     int ind = n_in_t[i];
//     for (unsigned int j = 0; j < length[ind]; ++j) {
//       // qua_in = qua[index[ind] + j];
//       obs_in = obs[index[ind] + j];
//     }
//   }
//   
// }





