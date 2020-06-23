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
#define MIN_COVERGE 4
// [[Rcpp::depends(RcppArmadillo)]]

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
    arma::Mat<int> m1 = as<arma::Mat<int>>(temp);
    arma::Mat<int> m2 = as<arma::Mat<int>>(temp2);
    arma::Mat<int> m3 = join_cols(m1, m2);
    if(sum_site[0] == 1 && sum >= 4.0) { // 0 is too few
      arma::Mat<int> m4 = {hap_site[1], hap_site[1], hap_site[1], hap_site[1]};
      arma::Mat<int> out = join_cols(m4, m3);
      temp = wrap(out);
      n_row = 7;
    } else {
      temp = wrap(m3);
      n_row = 6;
    }
  } else if (sum_site[0]/sum >= 0.62) {//2nd one appears 3 times
    temp = call_permute({hap_site[0], hap_site[0], hap_site[0], hap_site[1]});
    temp2 = Twopossible(hap_site);
    arma::Mat<int> m1 = as<arma::Mat<int>>(temp);
    arma::Mat<int> m2 = as<arma::Mat<int>>(temp2);
    arma::Mat<int> m3 = join_cols(m1, m2);
    if(sum_site[1] == 1 && sum >= 4.0) {
      arma::Mat<int> m4 = {hap_site[0], hap_site[0], hap_site[0], hap_site[0]};
      arma::Mat<int> out = join_cols(m4, m3);
      temp = wrap(out);
      n_row = 7;
    } else {
      temp = wrap(m3);
      n_row = 6;
    }
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
  // different rules to reduce the possibilities, basically rule out things like ATAG
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
    arma::Mat<int> m1 = as<arma::Mat<int>>(temp);
    small = {hap_site[1], hap_site[0]};
    temp = Fourpossible(small, hap_site[2]);
    arma::Mat<int> m2 = as<arma::Mat<int>>(temp);
    small = {hap_site[2], hap_site[0]};
    temp = Fourpossible(small, hap_site[1]);
    arma::Mat<int> m3 = as<arma::Mat<int>>(temp);
    arma::Mat<int> out = join_cols(m1, m2, m3);
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
  arma::Mat<int> m1 = as<arma::Mat<int>>(temp);
  arma::Mat<int> m2 = as<arma::Mat<int>>(temp2);
  arma::Mat<int> out = join_cols(m1, m2);
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
        arma::Mat<int> temp;
        IntegerMatrix inner_tmp(2, NUM_CLASS);
        for(l = 0; l < num; ++l)
          for(m = l + 1; m < num; ++m) {
            if(uni_alignment[ref_j] == "I") //deletion in B genome
              inner_tmp = call_permute_N({hap_site[l], hap_site[m]}, 0);
            else if(uni_alignment[ref_j] == "J")
              inner_tmp = call_permute_N({hap_site[l], hap_site[m]}, 1);
            temp = join_cols(temp, as<arma::Mat<int>>(inner_tmp));
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

// remake linkage (based on observed data)
IntegerMatrix remake_linkage(IntegerMatrix sub_link, unsigned int num) {
  unsigned int i, j, k, i1;
  arma::mat sub_uni = unique_rows(as<arma::mat>(sub_link));
  IntegerMatrix link_uni = wrap(sub_uni);
  List new_link(link_uni.nrow());
  unsigned int total_row = 0;
  for (i = 0; i < link_uni.nrow(); ++i) {
    List nuc_info = unique_map(link_uni(i, _));
    IntegerVector nuc = nuc_info["values"];
    IntegerVector nuc_count = nuc_info["lengths"];
    if (nuc.size() == 1 && nuc[0] == -1) // skip the read does not cover any site
      continue;
    // if(nuc_count[0] == num - 1 || nuc_count[0] == 1) 
    //   continue;
    if(nuc[0] != -1) { // read covers all site
      total_row++;
      new_link[i] = link_uni(i, _);
      continue;
    }
    IntegerVector idx(num);
    for(j = 0; j < num; ++j)
      if(link_uni(i, j) == -1)
        idx(j) = 1; // indicate -1 is here
      
      int move_out = 0;
      // if this read is contained in others
      for (i1 = 0; i1 < link_uni.nrow(); ++i1) {
        if (i1 == i)
          continue;
        int count = 0;
        // List nuc_info = unique_map(link_uni(i1, _));
        // IntegerVector nuc1 = nuc_info["values"];
        // IntegerVector nuc_count1 = nuc_info["lengths"];
        // if(nuc_count1[0] >= nuc_count[0])
        //   continue;
        for(j = 0; j < num; ++j)
          if(!idx(j))
            if(link_uni(i, j) == link_uni(i1, j))
              count++;
            if(count == num - nuc_count[0]) {
              // Rcout << i << "move" << "\n";
              move_out = 1;
              break;
            }
      }
      if(move_out)
        continue;
      List missing(num);
      IntegerVector flag(num);
      int in_row_num = 0;
      for (j = 0; j < num; ++j) {
        if (link_uni(i, j) == -1) {
          List nuc_col = unique_map(link_uni(_, j));
          IntegerVector nuc_unique = nuc_col["values"];
          // Rcout << nuc_unique << "\n";
          if(nuc_unique[0] == -1)
            missing[j] = nuc_unique[Range(1, nuc_unique.size() - 1)];
          else
            missing[j] = nuc_unique;
          in_row_num++;
        } else
          flag[j] = 1;
      }
      
      int add_row = 0;
      if(in_row_num != 1) {
        IntegerMatrix missing_rows = comb_element(missing, flag, in_row_num);
        add_row = missing_rows.nrow();
        IntegerMatrix new_link_i(add_row, num);
        int count = 0;
        for (j = 0; j < num; ++j) { // make fake reads with missing linkage info
          if(flag[j] != 1)
            new_link_i(_, j) = missing_rows(_, count++);
          else
            for (k = 0; k < add_row; ++k)
              new_link_i(k, j) = link_uni(i, j); // repeat the non-missing ones
        }
        new_link[i] = new_link_i;
      } else {
        IntegerMatrix new_link_i;
        IntegerVector tmp;
        for (j = 0; j < num; ++j)
          if(flag[j] != 1) {
            tmp = missing[j];
            add_row = tmp.size();
            new_link_i = IntegerMatrix(add_row, num);
          }
          for (j = 0; j < num; ++j) { // make fake reads with missing linkage info
            if(flag[j] != 1) {
              new_link_i(_, j) = tmp;
            } else {
              for (k = 0; k < add_row; ++k)
                new_link_i(k, j) = link_uni(i, j);
            }
          }
          new_link[i] = new_link_i;
      }
      total_row += add_row;
  }
  
  IntegerMatrix new_link_out(total_row, num);
  total_row = 0;
  for(i = 0; i < link_uni.nrow(); ++i) {
    if(new_link[i] == R_NilValue)
      continue;
    IntegerVector tmp = new_link[i];
    tmp.attr("dim") = Dimension(tmp.size()/num, num);
    IntegerMatrix new_link_i = as<IntegerMatrix>(tmp);
    for (k = 0; k < tmp.size()/num; ++k)
      new_link_out(total_row++, _) = new_link_i(k, _);
  }
  // finally, remove duplcated rows
  arma::mat new_linkage = unique_rows(as<arma::mat>(new_link_out));
  IntegerMatrix out = wrap(new_linkage);
  return(out);
}

IntegerVector best_branch(IntegerMatrix link, List transition, NumericVector initial, 
                          List possi_nuc, int i) {
  unsigned int j, k, l, m , w;
  List comb_in(link.ncol());
  List llk_in(link.ncol());
  int id = 0;
  // possible or determined nuc at each position
  for(j = 0; j < link.ncol(); ++j) {
    IntegerVector nuc = possi_nuc[j];
    if(link(i, j) != -1) {
      for(l = 0 ; l < nuc.size(); ++l)
        if(link(i, j) == nuc[l]) {
          id = l;
          break;
        }
        comb_in(j) = link(i, j);
    } else {
      comb_in(j) = nuc;
    }
  }
  IntegerVector flag(link.ncol());
  // IntegerMatrix poss_reads = comb_element(comb_in, flag, link.ncol());
  // get state likelihood
  for(j = 0; j < link.ncol(); ++j) {
    // Rcout << j << "\t";
    IntegerVector nuc = possi_nuc[j];
    if(link(i, j) != -1) {
      for(l = 0 ; l < nuc.size(); ++l)
        if(link(i, j) == nuc[l]) {
          id = l;
          // Rcout << "nuc " << nuc[l] << "\t";
          break;
        }
        if(j == 0) {
          llk_in(j) = initial[id];
          // Rcout << "ini " << initial[id] << "\t";
        } else {
          NumericMatrix trans = transition[j - 1];
          IntegerVector nuc2 = possi_nuc[j - 1];
          int id1 = 0;
          if(link(i, j - 1) != -1) {
            for(l = 0 ; l < nuc2.size(); ++l)
              if(link(i, j - 1) == nuc2[l]) {
                // Rcout << "nuc(j-1) " << nuc2[l] << "\t";
                id1 = l;
                break;
              }
              llk_in(j) = trans(id1, id);
              // Rcout << "trans " << trans(id1, id) << "\n";
          } else{
            // Rcout << "trans: all\n";
            llk_in(j) = trans(_, id);
          }
        }
    } else {
      if(j == 0)
        llk_in[j] = initial;
      else {
        NumericMatrix trans = transition[j - 1];
        IntegerVector nuc2 = possi_nuc[j - 1];
        if(link(i, j - 1) != -1) {
          for(l = 0 ; l < nuc2.size(); ++l)
            if(link(i, j) == nuc2[l]) {
              // Rcout << "nuc(j-1) " << nuc2[l] << "\n";
              id = l;
              break;
            }
            llk_in(j) = trans(id, _);
        } else {
          llk_in(j) = trans;
        }
      }
    }
  }
  IntegerVector hidden_state(llk_in.size());
  List path(llk_in.size());
  List backptr(llk_in.size() - 1);
  int b_next = 0;
  for(k = 0; k < llk_in.size(); ++k) {
    // Rcout << k << "\n";
    IntegerVector nuc = comb_in(k);
    NumericVector path_t(nuc.size());
    IntegerVector backptr_t(nuc.size());
    if(k == 0) {
      NumericVector trans = llk_in(k);
      for(l = 0; l < nuc.size(); ++l) {
        path_t(l) = trans(l);
      }
      // Rcout << path_t << "\n";
    } else {
      NumericVector tran = llk_in(k);
      int len = tran.size();
      int nrow = len/nuc.size();
      tran.attr("dim") = Dimension(nrow, nuc.size());
      NumericMatrix trans = as<NumericMatrix>(tran);
      NumericVector path_last = path[k - 1];
      for(m = 0; m < trans.ncol(); ++m) {
        double max = -INFINITY;
        int max_id = 0;
        for(w = 0; w < trans.nrow(); ++w) {
          double max_prob = path_last(w) + trans(w, m);
          if (max_prob > max) {
            max = max_prob;
            max_id = w;
          }
        }
        path_t(m) = max;
        backptr_t[m] = max_id;
      }
      backptr(k - 1) = backptr_t;
      // Rcout << path_t << "\n";
      // Rcout << backptr_t << "\n";
    }
    path(k) = path_t;
  }
  // Rcout << "decode\n";
  k = llk_in.size() - 1;
  double max = -INFINITY;
  IntegerVector nuc = comb_in(k);
  NumericVector path_t = path(k);
  for(m = 0; m < nuc.size(); ++m) {
    if (path_t(m) > max) {
      b_next = m;
      max = path_t(m);
    }
  }
  hidden_state(k) = nuc[b_next];
  
  while (k--) {
    // Rcout << k << "\n";
    IntegerVector nuc = comb_in(k);
    IntegerVector backptr_t = backptr(k);
    // Rcout << backptr_t << "\n";
    b_next = backptr_t[b_next];
    hidden_state(k) = nuc[b_next];
  }
  
  // get the combination of reads and corresponding transition matrix
  // List ls = List::create( // possible comb
  //   Named("poss_reads") = poss_reads,
  //   Named("comb_in") = comb_in,
  //   Named("llk_in") = llk_in, 
  //   Named("read") = hidden_state);
  
  return(hidden_state);
}

// use markov chain to make the linkage
IntegerMatrix mc_linkage(IntegerMatrix sub_link, int num) {
  unsigned int i, j, k, l;
  //remove non-covered reads
  NumericVector initial;
  IntegerMatrix link_pre(sub_link.nrow(), sub_link.ncol());
  int count = 0;
  for(i = 0; i < sub_link.nrow(); ++i) {
    IntegerVector read = sub_link(i, _);
    int rowsum = sum(read);
    if(rowsum == -num)
      continue;
    link_pre(count++, _) = sub_link(i, _);
  }
  IntegerMatrix link = link_pre(Range(0, count - 1), _);
  List transition(link.ncol() - 1);
  List possi_nuc(link.ncol());
  for(j = 0 ; j < link.ncol() - 1; ++j) {
    List nuc_info = unique_map(link(_, j));
    IntegerVector nuc = nuc_info["values"];
    IntegerVector nuc_count = nuc_info["lengths"];
    // Rcout << j << "\t" << nuc << "\t" << nuc_count << "\n";
    int state1 = nuc.size();
    int start = 0;
    possi_nuc[j] = nuc;
    if(nuc[0] == -1) {
      state1 = nuc.size() - 1;
      start = 1;
      possi_nuc[j] = nuc[Range(1, nuc.size() - 1)];
    }
    if(j == 0) {
      double total = count;
      initial = IntegerVector(state1);
      if(nuc[0] == -1) {
        for(k = 1; k < nuc.size(); ++k)
          initial[k - 1] = log(nuc_count[k]/(total - nuc_count[0])); // in the ascending order
      }
      else {
        for(k = 0; k < nuc.size(); ++k)
          initial[k] = log(nuc_count[k]/total);
      }
    }
    // unique rows and the count
    List nuc1_info = unique_map(link(_, j + 1));
    IntegerVector nuc1 = nuc1_info["values"];
    int state2 = nuc1.size();
    possi_nuc[j + 1] = nuc1;
    if(nuc1[0] == -1) {
      possi_nuc[j + 1] = nuc1[Range(1, nuc1.size() - 1)];
      state2 = nuc1.size() - 1;
    }
    
    IntegerMatrix link_unique(link.nrow(), 2);
    int unique_ct = 0;
    for(k = 0; k < link.nrow(); ++k)
      if(link(k, j + 1) != -1 && link(k, j) != -1) {
        link_unique(unique_ct, 0) = link(k, j );
        link_unique(unique_ct++, 1) = link(k, j + 1);
      }
      // for(k = 0; k < link.nrow(); ++k)
      //   for(l = 0; l < 2; ++l)
      //     Rcout << link_unique(k, j) << "\t";
      IntegerMatrix in_link = link_unique(Range(0, unique_ct - 1), _);
      List unique_row = hash_mat(in_link);
      List all_id = unique_row["all_id"];
      IntegerVector idx = unique_row["idx"];
      NumericMatrix trans(state1, state2);
      for(k = 0; k < state1; ++k)
        for(i = 0; i < state2; ++i)
          trans(k, i) = R_NegInf;
      // Rcout << state1 << "\t" << state2 << "\n";
      List nuc_info_uni = unique_map(in_link(_, 0));
      IntegerVector nuc_count_uni = nuc_info_uni["lengths"];
      IntegerVector nuc_uni = nuc_info_uni["values"];
      
      for(k = start; k < nuc.size(); ++k) {
        // Rcout<< "\n" << nuc[k] << "\n";
        int nuc_id = 0;
        for(l = 0; l < nuc_uni.size(); ++l)
          if(nuc_uni[l] == nuc[k]) {
            nuc_id = l;
            break;
          }
          for(i = 0; i < idx.size(); ++i) {
            double de = nuc_count_uni[nuc_id];
            IntegerVector sub_read = in_link(idx[i], _);
            // Rcout << "unique sub " << sub_read << "\n";
            // IntegerVector sub_read = ordered_read(i, _);
            // if(sub_read[1] == -1) {
            //   de--;
            //   continue;
            // }
            if(nuc[k] != sub_read[0])
              continue;
            IntegerVector sub_id = all_id[i];
            // Rcout << "uniques " << sub_id << "\n";
            double nu = sub_id.size();
            int col_id = 0;
            int flag = 0;
            if(state2 < nuc1.size())
              flag = 1;
            for(l = 0; l < nuc1.size(); ++l)
              if(nuc1[l] != -1 && nuc1[l] == sub_read[1])
                col_id = l - flag;
              // Rcout <<"col_id " << col_id << "\n" ;
              // if(nuc[k] == sub_read[0]) {
              // Rcout <<"de " << de << " nu " << nu<< "\n" ;
              trans(k - start, col_id) = log(nu/de); // log likelihood
              // }
          }
      }
      transition[j] = trans;
  }
  // 
  // // now use MC to impute the missing nuc
  arma::mat uniqu_link = unique_rows(as<arma::mat>(link));
  IntegerMatrix sub_uni_link = wrap(uniqu_link);
  IntegerMatrix mc_reads(sub_uni_link.nrow(), sub_uni_link.ncol());
  for(i = 0; i < sub_uni_link.nrow(); ++i) {
    int flag = 0;
    for(j = 0; j < sub_uni_link.ncol(); ++j) {
      if(sub_uni_link(i, j) == -1) {
        flag = 1;
        break;
      }
    }
    if(!flag) {
      mc_reads(i, _) = sub_uni_link(i, _);
      continue;
    }
    IntegerVector tmp = sub_uni_link(i, _);
    // Rcout << tmp << "\n";
    mc_reads(i, _) = best_branch(sub_uni_link, transition, initial, possi_nuc, i);
  }
  
  arma::mat uniqu = unique_rows(as<arma::mat>(mc_reads));
  IntegerMatrix uni = wrap(uniqu);
  // List ls = List::create(
  //   Named("transition") = transition, // possible comb
  //   Named("initial") = initial,
  //   Named("possi_nuc") = possi_nuc,
  //   Named("link") = sub_uni_link,
  //   Named("reads") = uni);
  // 
  return(uni);
}
// limit combination based on linkage information. For the states not start a new sequence(overlap, 
// keep the last overlapped combination and then make the rest)
List limit_comb_t0(IntegerMatrix combination, List hidden_states, IntegerVector location,
                   IntegerMatrix linkage_info, unsigned int num, unsigned int start_idx, 
                   unsigned int num_states, unsigned int use_MC) {
  unsigned int i, j, k, idx, m;
  IntegerMatrix sub_hap(NUM_CLASS, num);
  IntegerMatrix old_sub_link = linkage_info(_, Range(start_idx, start_idx + num - 1));
  IntegerVector exclude(num_states);
  int count, linkage_len, all_excluded;
  linkage_len = num - 1; 
  // int cut_off;
  // all_excluded = num_states;
  //remake the linkage
  IntegerMatrix sub_link;
  if(!use_MC)
    sub_link = remake_linkage(old_sub_link, num);
  else
    sub_link = mc_linkage(old_sub_link, num);
  unsigned int n_observation = sub_link.nrow();
  // while (all_excluded == num_states) {
  //   cut_off = NUM_CLASS;
  //   // Rcout << "linkage length " << linkage_len << "\n"; //actual linkage length + 1
  //   while (cut_off >= 2 && all_excluded == num_states) {
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
          // Rcout << "\n read " << "\n";
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
        if (count != NUM_CLASS) {
          exclude(m) = 1;
          all_excluded++;
        }
      }
  //     cut_off--; //TODO:how to make sure include more relaiable possibles
  //   }
  //   linkage_len--;
  // }
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

//fill haplotype at the rest sites (with variation) at time t
IntegerMatrix make_hap(List hidden_states, IntegerMatrix haplotype, IntegerVector location, unsigned int p_tmax,
                       IntegerVector combination, unsigned int time_pos, unsigned int num, int hap_min_pos) {
  unsigned int j, k, idx;
  
  for(j = 0; j < num; ++j) {
    IntegerMatrix hidden = hidden_states[location[j]];
    //Rcout << hidden << "\n";
    idx = combination[j];
    //Rcout << idx << "\n";
    for(k = 0; k < NUM_CLASS; ++k)
      haplotype(k, location[j]) = hidden(idx, k);
  }
  return(haplotype(_, Range(time_pos - hap_min_pos, time_pos + p_tmax - hap_min_pos - 1)));
}

/*
 * Return the combination of sites with variation within one time t[NOTE:time_pos might not start from 0, but undecided_pos starts from 0]
 */
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

IntegerMatrix new_combination(List hmm_info, IntegerVector location, IntegerVector overlapped, IntegerVector exclude_last, IntegerMatrix overlap_comb, 
                              IntegerVector overlap_loci, IntegerMatrix linkage_info, unsigned int overlap_new_states, unsigned int overlap_num_states,
                              unsigned int use_MC) {
  
  IntegerMatrix first_comb = unique_overlap(overlapped, exclude_last, overlap_comb, overlap_loci, 
                                            overlap_new_states, overlap_num_states);
  //find combination of the rest position(include 1 overlap to make sure the linkage)
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
  // get the coombination of the non-overlapped location (include last overlapped)
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
  
  // Rcout << "exist: " << exist << "\n";
  // Rcout << "allowed: " << allowed << "\n";
  // Rcout << "num_comb: " << combination.nrow() << "\n";
  if(!setequal(exist, allowed))
    flag = 1;
      
  if(flag) {  // if the exist contains more possibility than allowed, only keeps the allowed          
    for(m = 0; m < combination.nrow(); ++m)
      for(w = 0; w < allowed.size(); ++w)
        if(allowed(w) == combination(m, 0))
          new_combination(num++, _) = combination(m, _);
    exclude_info = limit_comb_t0(new_combination(Range(0, num - 1), _), hidden_states, left_loci, linkage_info, combination.ncol(), start_idx, num, use_MC);
  } else {
    exclude_info = limit_comb_t0(combination, hidden_states, left_loci, linkage_info, combination.ncol(), start_idx, combination.nrow(), use_MC);
  }
        
  // Now give the limited combination
  int num_states = exclude_info["num_states"];
  IntegerVector exclude = exclude_info["exclude"];
  // Rcout << "exclude "<< exclude << "\n";
  IntegerMatrix next_comb(num_states, combination.ncol());
  count = 0;
  // Now combine first and second part[make sure the connection states appear]
  if(flag) {
    for(m = 0; m < num; ++m)
      if(!exclude[m])
          next_comb(count++, _) = new_combination(m, _);
  } else {
      for(m = 0; m < combination.nrow(); ++m)
        if(!exclude[m])
          next_comb(count++, _) = combination(m, _);
  }
  // if next_comb does not contain one of the states in allowed, add it back (use the one w/ smallest index)
  // this will introduce more states not shown in the reads linkage, but to keep the trans works, have to...
  IntegerVector new_exist = next_comb(_, 0);
  if(!setequal(new_exist, allowed)) {
    IntegerVector diff = setdiff(allowed, new_exist);
    Rcout << "add more possibilities " << diff << "\n";;
    IntegerMatrix extra(diff.size(), combination.ncol());
    arma::Mat<int> m1 = as<arma::Mat<int>>(next_comb);
    next_comb = IntegerMatrix(num_states + diff.size(), combination.ncol());
    count = 0;
    if(flag) {
      for(w = 0; w < diff.size(); ++w)
        for(m = 0; m < new_combination.nrow(); ++m) {
          IntegerVector tmp = new_combination(m, _);
          if(new_combination(m, 0) == diff[w]) {
            extra(count++, _) = new_combination(m, _);
            break;
          }
        }
    } else {
      for(w = 0; w < diff.size(); ++w)
        for(m = 0; m < combination.nrow(); ++m) {
          IntegerVector tmp = combination(m, _);
          if(combination(m, 0) == diff[w]) {
            extra(count++, _) = combination(m, _);
            break;
          }
        }
    }
    arma::Mat<int> m2 = as<arma::Mat<int>>(extra);
    m1.insert_rows(1, m2);
    next_comb = wrap(m1);
  }
 
  IntegerVector new_col = next_comb(_, 0);
  List second_uni = unique_map(new_col); // start might not from 0 (e.g. ailgnment starts from 2)
  IntegerVector n2 = second_uni["lengths"];
  // IntegerVector possible = second_uni["values"];
  int all = 0;
  for(m = 0; m < n2.size(); ++m)
    all += n2[m] * n1[m];
  IntegerMatrix final_comb(all, location.size());
  
  all = 0;
  for(k = 0; k < last_col.size(); ++k)
    for(w = 0; w < new_col.size(); ++w) {
      if(new_col(w) == last_col(k)) {
        for(j = 0; j < overlap_len; ++j)
          final_comb(all, j) = first_comb(k, j);
        for(i = 1; i < len; ++i)
          final_comb(all, i + overlap_len - 1) = next_comb(w, i);
        all++;
      }
    }
        
  return(final_comb);
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





