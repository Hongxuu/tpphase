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

#include "data_format.h"
#include "baumwelch.h"
#include "hmm_state.h"
using namespace Rcpp;
using namespace std;

#define MLOGIT_CLASS 4
#define NUM_CLASS 4

// [[Rcpp::depends(RcppArmadillo)]]
List ini_hmm (unsigned int t_max,  IntegerVector num_states, Rcpp::Nullable<Rcpp::List> trans_indicator = R_NilValue) {
  NumericVector phi(num_states[0]);
  List trans(t_max - 1);
  // List emit(t_max);
  unsigned int w, m, t;
  
  for(w = 0; w < num_states[0]; ++w)
    phi[w] = 1/double(num_states[0]);
  
  for(t = 0; t < t_max - 1; ++t) {
    NumericMatrix transition(num_states[t], num_states[t + 1]);
    // 
    // if(trans_indicator.isNotNull()) {
    //   IntegerMatrix trans_ind = trans_indicator[t];
    //   for (w = 0; w < num_states[t]; ++w) {
    //     int new_num = num_states[t + 1];
    //     for (m = 0; m < num_states[t + 1]; ++m)
    //       if (trans_ind(w, m)) // 1 means not the same, so cannot b transferred
    //         new_num--;
    //       for (m = 0; m < num_states[t + 1]; ++m) 
    //         if (!trans_ind(w, m)) 
    //           transition(w, m) = 1/double(new_num);
    //   }
    // } else {
      for (w = 0; w < num_states[t]; ++w)
        for (m = 0; m < num_states[t + 1]; ++m) 
          transition(w, m) = 1/double(num_states[t + 1]);
    // }
    trans(t) = transition;
  }
  // for(t = 0; t < t_max; ++t) {
  //   NumericVector emission(num_states[t]);
  //   for (w = 0; w < num_states[t]; ++w)
  //     emission(w) = 1/double(num_states[t]);
  //   emit(t) = emission;
  // }
  // 
  List par_hmm = List::create(
    Named("phi") = phi,
    // Named("emit") = emit
    Named("trans") = trans);
  return(par_hmm);
}

List forward_backward(List par_hmm, unsigned int t_max, IntegerVector num_states)
{
  NumericVector phi = par_hmm["phi"];
  List trans = par_hmm["trans"];
  List emit = par_hmm["emit"];
  
  List alpha(t_max);
  List beta_wt(t_max);
  unsigned int w, t, m;
  
  for (t = 0; t < t_max; ++t) {
    NumericVector alp(num_states[t]);
    NumericVector emission = emit(t);
    NumericVector alp_last;
    if (t == 0)
      for (w = 0; w < num_states[t]; ++w)
        alp(w) = phi[w] * emission(w);
    else {
      alp_last = alpha(t - 1);
      NumericMatrix transition = trans(t - 1);
      for (w = 0; w < num_states[t]; ++w) {
        for (m = 0; m < num_states[t - 1]; ++m)
          alp(w) += alp_last(m) * transition(m, w);
        alp(w) = alp(w) * emission(w); // will this cause underflow?
      }
    }
    alpha(t) = alp;
  }
  
  for (t = t_max; t --> 0;) {
    NumericVector beta(num_states[t]);
    NumericVector beta_after;
    
    if(t == t_max - 1) {
      for (w = 0; w < num_states[t]; ++w)
        beta(w) = 1.0;
    } else {
      NumericVector emission = emit(t + 1);
      NumericMatrix transition = trans(t);
      beta_after = beta_wt(t + 1);
      for (w = 0; w < num_states[t]; ++w) {
        for (m = 0; m < num_states[t + 1]; ++m)
          beta(w) += beta_after(m) * transition(w, m) * emission(m);
      }
    }
    beta_wt(t) = beta;
  }
  
  List gamma(t_max);
  double sum = 0;
  for (t = 0; t < t_max; ++t) {
    NumericVector beta = beta_wt(t);
    NumericVector alp = alpha(t); // maybe change this to a long vector instead a list?
    sum = 0;
    NumericVector gam(num_states[t]);
    for (w = 0; w < num_states[t]; ++w) {
      gam(w) = alp(w) * beta(w);
      sum += gam(w);
    }
    for (w = 0; w < num_states[t]; ++w)
      gam(w) /= sum;
    gamma(t) = gam;
  }
  
  List xi(t_max - 1);
  //arma::Cube<double> xi(total_state, total_state, t_max);
  for (t = 0; t < t_max - 1; ++t) {
    NumericVector emission = emit(t + 1);
    NumericMatrix transition = trans(t);
    NumericVector beta = beta_wt(t + 1);
    NumericVector alp = alpha(t);
    NumericMatrix x(num_states[t], num_states[t + 1]);
    sum = 0;
    for (w = 0; w < num_states[t]; ++w)
      for (m = 0; m < num_states[t + 1]; ++m) {
        // x(w, m) = log(alp(w)) + log(transition(w, m)) + log(beta(m)) + log(emission(m)); // avoid underflow?
        // x(w, m) = exp(x(w, m));
        x(w, m) = alp(w) * transition(w, m) * beta(m) * emission(m); // not logging it since some transition are 0
        sum += x(w, m);
      }
      for (w = 0; w < num_states[t]; ++w)
        for (m = 0; m < num_states[t + 1]; ++m)
          x(w, m) /= sum;
    xi(t) = x;
  }
  List par_hmm_out = List::create(
    Named("phi") = phi,
    Named("trans") = trans,
    Named("emit") = emit,
    Named("alpha") = alpha,
    Named("beta_wt") = beta_wt,
    Named("gamma") = gamma,
    Named("xi") = xi);
  return(par_hmm_out);
}

NumericVector update_phi(List par_hmm, unsigned int num_states) {
  List gamma = par_hmm["gamma"];
  NumericVector phi(num_states);
  NumericVector gam = gamma[0];
  
  for (unsigned int w = 0; w < num_states; ++w)
    phi[w] = gam(w);
  return phi;
}

List update_trans(List par_hmm, unsigned int t_max, IntegerVector num_states) {
  List gamma = par_hmm["gamma"];
  List xi = par_hmm["xi"];
  List trans(t_max - 1);
  unsigned int t, m, w;
  for (t = 0; t < t_max - 1; ++t) {
    NumericMatrix transition(num_states[t], num_states[t + 1]);
    NumericVector gam = gamma[t];
    NumericMatrix x = xi[t];
    for (w = 0; w < num_states[t]; ++w)
      for (m = 0; m < num_states[t + 1]; ++m)
        transition(w, m) = x(w, m)/gam(w);
    trans[t] = transition;
  }
  return(trans);
}

double site_likelihood (unsigned int i, unsigned int K, NumericMatrix beta, unsigned int PD_LENGTH, 
                        List dat_info, IntegerMatrix haplotype, unsigned int time_pos) {
  unsigned int j, l;
  double qua_in, read_pos_in, ref_pos_in;
  double tail, sum;
  double xb = 0;
  
  IntegerVector qua = dat_info["qua"];
  IntegerVector obs = dat_info["nuc"];
  IntegerVector ref_pos = dat_info["ref_pos"];
  IntegerVector read_pos = dat_info["read_pos"];
  IntegerVector length = dat_info["length"];
  IntegerVector index = dat_info["start_id"];
  IntegerMatrix ref_index = dat_info["ref_idx"];
  
  NumericVector hap_nuc(MLOGIT_CLASS);
  NumericVector hnuc_qua(MLOGIT_CLASS);
  NumericVector pred_beta(MLOGIT_CLASS - 1);
  double read_llk = 0.;
  arma::vec predictor(PD_LENGTH);
  //NumericVector site_llk(length[i]);
  
  /* non-indel positions for both reads and haplotypes */
  for(j = 0; j < length[i]; ++j) {
    
    //Rprintf("j %d, position %d\n", j, index[i] + j);
    read_pos_in = read_pos[index[i] + j];
    qua_in = qua[index[i] + j];
    ref_pos_in = ref_pos[index[i] + j]; // this is different since the sub-hap is counted 0 from the time t
    
    for (l = 0; l < MLOGIT_CLASS; ++l) {
      hap_nuc[l] = 0;
      hnuc_qua[l] = 0;
    }
    //Rcout << "haplotype " << haplotype(K, ref_pos_in) << "\t";
    
    if(haplotype(K, ref_pos_in - time_pos) == 1) {
      hap_nuc[0] = 1;
      hnuc_qua[0] = qua_in;
    } else if(haplotype(K, ref_pos_in - time_pos) == 3) {
      hap_nuc[1] = 1;
      hnuc_qua[1] = qua_in;
    } else if(haplotype(K, ref_pos_in - time_pos) == 2) {
      hap_nuc[3] = 1;
      hnuc_qua[3] = qua_in;
    }
    
    predictor = {1, read_pos_in, ref_pos_in, qua_in, hap_nuc[0], hap_nuc[1],
                 hap_nuc[3], hnuc_qua[0], hnuc_qua[1], hnuc_qua[3]};
    
    arma::mat beta_ar = as<arma::mat>(beta);
    arma::vec pb = beta_ar.t() * predictor;
    pred_beta = as<NumericVector>(wrap(pb));
    
    //Rcout << "predictor : " << predictor.t() << "\n";
    //Rcout << "beta_ar : " << beta_ar << "\n";
    //Rcout << "pred_beta : " << pred_beta << "\n";
    
    sum = 0.0;
    for (l = 0; l < MLOGIT_CLASS - 1; ++l)
      sum += exp(pred_beta[l]);
    tail = log(1/(1 + sum));
    
    if(obs[index[i] + j] == 1) {
      xb = pred_beta[0];
    } else if(obs[index[i] + j] == 3) {
      xb = pred_beta[1];
    } else if(obs[index[i] + j] == 2) {
      xb = pred_beta[2];
    } else if(obs[index[i] + j] == 0) {
      xb = 0;
    }
    read_llk += xb + tail; /* Notice here use log likelihood, not likelihood */
  }
  //Rcout << "read_llk: " << read_llk << "\t";
  return read_llk;
} /* site_likelihood */


/*
 * Fill the positions without variation
 */
// [[Rcpp::export]]
IntegerMatrix fill_all_hap(List hidden_states, unsigned int hap_length, IntegerVector n_row) {
  unsigned int j, k;
  IntegerMatrix haplotype(NUM_CLASS, hap_length);
  IntegerVector nuc_j(NUM_CLASS);
  for(j = 0; j < hap_length; ++j)
    if(n_row(j) == 1) { // fill the loci with only 1 possibility, can be optimized by using the part has been filled and add/trim
      nuc_j = hidden_states[j];
      for(k = 0; k < NUM_CLASS; ++k)
        haplotype(k, j) = nuc_j[k];  
    }
  return haplotype;
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

// [[Rcpp::export]]
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
      } 
      else 
        link(i, j) = -1; // meaning not covered
    }
  }
  return(link);
}

// [[Rcpp::export]]
List full_hap_new (List hmm_info, IntegerMatrix linkage_info, unsigned int hap_length, int hap_min_pos) {
  List hidden_states = hmm_info["hidden_states"];
  IntegerVector num_states = hmm_info["num_states"];
  IntegerVector time_pos = hmm_info["time_pos"];
  IntegerVector p_tmax = hmm_info["p_tmax"];
  IntegerVector n_row = hmm_info["n_row"];
  IntegerVector pos_possibility = hmm_info["pos_possibility"];
  IntegerVector undecided_pos = hmm_info["undecided_pos"];
  unsigned int t_max = hmm_info["t_max"];
  unsigned int t, m;
  List full_hap(t_max);
  List comb(t_max);
  IntegerMatrix hap = fill_all_hap(hidden_states, hap_length, n_row);
  IntegerVector new_num_states(t_max);
  List overlap_info = get_overlap(p_tmax, time_pos, num_states, undecided_pos, t_max, hap_min_pos);
  List overlapped = overlap_info["overlapped"];
  IntegerVector overlapped_id = overlap_info["overlapped_id"];
  int start_t = overlap_info["start_t"];
  List loci = overlap_info["location"];
  
  //start t info
  List comb_info_t0 = find_combination(undecided_pos, pos_possibility, p_tmax[start_t], time_pos[start_t], hap_min_pos);
  IntegerVector location = comb_info_t0["location"];
  IntegerMatrix combination = comb_info_t0["combination"];
  unsigned int num = comb_info_t0["num"];
  List t0 = limit_comb_t0(combination, hidden_states, location, linkage_info, num, 0, num_states[start_t]);
  IntegerVector exclude = t0["exclude"];
  IntegerVector exclude_last = exclude;
  IntegerMatrix comb_in = combination; 
  // get the states  
  for(t = 0; t < t_max; ++t) {
    List full_hap_t;
    if(num_states[t] != 1) {
      int count = 0;
      if(t == start_t) {
        // decide start_t first
        new_num_states[t] = t0["num_states"];
        IntegerMatrix new_comb(new_num_states[t], combination.ncol());
        full_hap_t = List(new_num_states[t]);
        for(m = 0; m < num_states[t]; ++m)
          if(!exclude[m]) {
            IntegerMatrix haplotype = make_hap(hidden_states, hap, location, p_tmax[t], combination(m, _), time_pos[t], num, hap_min_pos);
            new_comb(count, _) = combination(m, _);
            full_hap_t(count++) = haplotype;
          }
          comb[t] = new_comb;
          Rcout << "start done" << "\n";
      }
      else {
        int identical = 0;
        int last_t = overlapped_id[t];
        IntegerVector overlapped_t = overlapped[t];
        IntegerVector loci_lastt = loci[last_t];
        IntegerVector loci_currt = loci[t];
        int old_state;
        Rcout << t << "\t" << last_t << "\n";
        Rcout << "overlapped: " << overlapped_t << "\n";
        Rcout << "loci_lastt: " << loci_lastt << "\n";
        Rcout << "loci_currt: " << loci_currt << "\n";
        
        if(loci_lastt[0] <= loci_currt[0] && loci_lastt[loci_lastt.size() - 1] >= loci_currt[loci_currt.size() - 1]) {
          if(loci_lastt.size() > loci_currt.size()) { // if current is in its overlap
            // get the unique overlapped combination from the last t
            if(last_t == start_t) {
              exclude_last = IntegerVector(exclude.size());
              exclude_last = exclude;
              old_state = num_states[last_t];
              comb_in = combination;
            }
            else {
              exclude_last = IntegerVector(new_num_states[last_t]);
              for(m = 0; m < new_num_states[last_t]; ++m)
                exclude_last[m] = 0;
              old_state = new_num_states[last_t];
              IntegerMatrix tmp = comb[last_t];
              comb_in = IntegerMatrix(tmp.nrow(), tmp.ncol());
              comb_in = tmp;
            }
            Rcout << "last t:\t" << new_num_states[last_t] << "\told\t" << old_state << "\n";
            IntegerMatrix new_comb = unique_overlap(overlapped_t, exclude_last, comb_in, loci_lastt, new_num_states[last_t], old_state);
            new_num_states[t] = new_comb.nrow();
            comb[t] = new_comb;
            full_hap_t = List(new_num_states[t]);
            for(m = 0; m < new_num_states[t]; ++m) {
              IntegerMatrix haplotype = make_hap(hidden_states, hap, loci_currt, p_tmax[t], new_comb(m, _), time_pos[t], loci_currt.size(), hap_min_pos);
              full_hap_t(count++) = haplotype;
            }
          } else if (loci_lastt.size() == loci_currt.size()) {
            for(m = 0; m < loci_lastt.size(); ++m)
              if(loci_lastt[m] != loci_currt[m]) {
                identical = 1;
                break;
              }
              if(!identical) {
                Rcout << "same sites\n";
                comb[t] = comb[last_t];
                new_num_states[t] = new_num_states[last_t];
                full_hap_t = List(new_num_states[t]);
                full_hap_t = full_hap[last_t];
              }
          }
        } 
        else {
          if(last_t == start_t) {
            exclude_last = IntegerVector(exclude.size());
            exclude_last = exclude;
            old_state = num_states[last_t];
            comb_in = combination;
          }
          else {
            exclude_last = IntegerVector(new_num_states[last_t]);
            for(m = 0; m < new_num_states[last_t]; ++m)
              exclude_last[m] = 0;
            old_state = new_num_states[last_t];
            IntegerMatrix tmp = comb[last_t];
            comb_in = IntegerMatrix(tmp.nrow(), tmp.ncol());
            comb_in = tmp;
          }
          IntegerMatrix new_comb = new_combination(hmm_info, loci_currt, overlapped_t, exclude_last, comb_in,
                                                   loci_lastt, linkage_info, new_num_states[last_t], old_state);
          new_num_states[t] = new_comb.nrow();
          comb[t] = new_comb;
          full_hap_t = List(new_num_states[t]);
          for(m = 0; m < new_num_states[t]; ++m) {
            IntegerMatrix haplotype = make_hap(hidden_states, hap, loci_currt, p_tmax[t], new_comb(m, _), time_pos[t], loci_currt.size(), hap_min_pos);
            full_hap_t(count++) = haplotype;
          }
        }
      }
    } 
    else {
      full_hap_t = List(1);
      full_hap_t[0] = hap(_, Range(time_pos[t] - hap_min_pos, time_pos[t] + p_tmax[t] - hap_min_pos - 1));
      new_num_states[t] = 1;
      comb[t] = -1;
    }
    Rcout << "new no. states: " << new_num_states[t] << "\n";
    full_hap[t] = full_hap_t;
  }
  
  List out = List::create(
    Named("full_hap") = full_hap,
    Named("new_num_states") = new_num_states, 
    Named("combination") = comb);
  
  return(out);
}

// function to return all the possible haps at each time t, further reduce the number of hidden states as well
// [[Rcpp::export]]
List full_hap(List hmm_info, IntegerMatrix linkage_info, unsigned int hap_length, int hap_min_pos) {
  List hidden_states = hmm_info["hidden_states"];
  IntegerVector num_states = hmm_info["num_states"];
  IntegerVector time_pos = hmm_info["time_pos"];
  IntegerVector p_tmax = hmm_info["p_tmax"];
  IntegerVector n_row = hmm_info["n_row"];
  IntegerVector pos_possibility = hmm_info["pos_possibility"];
  IntegerVector undecided_pos = hmm_info["undecided_pos"];
  unsigned int t_max = hmm_info["t_max"];
  unsigned int t, m, j, start_idx = 0;
  List full_hap(t_max);
  List comb(t_max);
  IntegerMatrix hap = fill_all_hap(hidden_states, hap_length, n_row);
  IntegerVector new_num_states(t_max);
  for(t = 0; t < t_max; ++t) {
    // find the start id for linkage_info
    unsigned int id = time_pos[t] - hap_min_pos;
    for(j = 0; j < undecided_pos.size(); ++j) {
      start_idx = j;
      if(undecided_pos[j] >= id)
        break;
    }
    // give h_t, each t has many possible combinations
    List full_hap_t(num_states(t));
    if(num_states[t] != 1) {
      // if the current covered the same site as the last one, then directly use the last combined information
      
      List comb_info = find_combination(undecided_pos, pos_possibility, p_tmax[t], time_pos[t], hap_min_pos);
      IntegerVector location = comb_info["location"];
      // Rcout << t << "\t" << start_idx << "\t" << location << "\n";
      IntegerMatrix combination = comb_info["combination"];
      unsigned int num = comb_info["num"];
      List exclude_info = limit_comb_t0(combination, hidden_states, location, linkage_info, num, start_idx, num_states[t]);
      IntegerVector exclude = exclude_info["exclude"];
      new_num_states[t] = exclude_info["num_states"];
      IntegerMatrix new_comb(new_num_states[t], combination.ncol());
      int count = 0;
      for(m = 0; m < num_states[t]; ++m)
        if(!exclude[m]) {
          IntegerMatrix haplotype = make_hap(hidden_states, hap, location, p_tmax[t], combination(m, _), time_pos[t], num, hap_min_pos);
          new_comb(count, _) = combination(m, _);
          full_hap_t(count++) = haplotype;
        }
      comb[t] = new_comb;
    } else {
      full_hap_t[0] = hap(_, Range(time_pos[t] - hap_min_pos, time_pos[t] + p_tmax[t] - hap_min_pos - 1));
      new_num_states[t] = num_states[t];
    }
    full_hap[t] = full_hap_t[Range(0, new_num_states[t] - 1)];
  }
  
  List out = List::create(
    Named("full_hap") = full_hap,
    Named("new_num_states") = new_num_states, 
    Named("combination") = comb);
  
  return(out);
}

List compute_emit(List hmm_info, List dat_info, List hap_info, NumericMatrix beta, NumericVector eta, int PD_LENGTH) 
{
  IntegerVector num_states = hmm_info["num_states"];
  IntegerVector n_t = hmm_info["n_t"];
  List n_in_t = hmm_info["n_in_t"];
  IntegerVector time_pos = hmm_info["time_pos"];
  IntegerVector p_tmax = hmm_info["p_tmax"];
  unsigned int t_max = hmm_info["t_max"];
  unsigned int i, k, m, t;
  double sum_emit_prob;
  List w_ic(t_max);
  List emit(t_max);
  NumericVector weight_llk(NUM_CLASS);
  // compute emission based on the initial value of eta and beta
  for(t = 0; t < t_max; ++t) {
    // give h_t, each t has many possible combinations
    List full_hap_t = hap_info(t);
    List w_icm(num_states[t]);
    NumericVector emission(num_states[t]);
    sum_emit_prob = 0;
    IntegerVector idx = n_in_t[t];
    for(m = 0; m < num_states[t]; ++m) {
      //Rcout << m;
      NumericVector read_likelihood(n_t[t]);
      NumericMatrix read_class_llk(n_t[t], NUM_CLASS);
      IntegerMatrix haplotype = full_hap_t(m);
      NumericMatrix w_tic(n_t[t], NUM_CLASS);
      for (i = 0; i < n_t[t]; ++i) {
        unsigned int id = idx[i];
        for (k = 0; k < NUM_CLASS; ++k) {
          read_class_llk(i, k) = site_likelihood(id, k, beta, PD_LENGTH, dat_info, haplotype, time_pos[t]);
          weight_llk[k] = eta[k] * exp(read_class_llk(i, k));
          read_likelihood[i] += weight_llk[k];
        }
        for (k = 0; k < NUM_CLASS; ++k)
          w_tic(i, k) = weight_llk[k]/read_likelihood[i];
        emission[m] += log(read_likelihood[i]);
      }
      emission[m] = exp(emission[m]);
      sum_emit_prob += emission[m];
      w_icm(m) = w_tic;
      // Rcout << w_tic;
    }
    w_ic(t) = w_icm;
    for(m = 0; m < num_states[t]; ++m) 
      emission[m] = emission[m]/sum_emit_prob;
    if(num_states[t] == 1)
      emission[0] = 1;
    emit(t) = emission;
  }
  List par_hmm_out = List::create(
    Named("emit") = emit,
    Named("w_ic") = w_ic);
  return(par_hmm_out);
}

NumericVector update_eta(List w_ic, List gamma, IntegerVector num_states, IntegerVector n_t, unsigned int t_max, unsigned int n_observation) {
  // update eta where P(Ht|R, theta) = gamma
  double inner;
  unsigned int i, k, m, t;
  NumericVector eta_new(NUM_CLASS);
  
  for(k = 0; k < NUM_CLASS; ++k) {
    for(t = 0; t < t_max; ++t) {
      List w_icm = w_ic(t);
      NumericVector gam = gamma(t); 
      inner = 0;
      for(m = 0; m < num_states[t]; ++m) {
        NumericMatrix w_tic = w_icm(m);
        double sum_wic = 0;
        for (i = 0; i < n_t[t]; ++i) {
          sum_wic += w_tic(i, k);
        }
        inner += gam(m) * sum_wic;
      }
      eta_new[k] += inner;
    }
  }
  for(k = 0; k < NUM_CLASS; ++k) 
    eta_new[k] = eta_new[k]/n_observation;
  
  return(eta_new);
}

List sub_sample(List hmm_info, List dat_info) {
  unsigned int i, j, t;
  IntegerVector qua = dat_info["qua"];
  IntegerVector obs = dat_info["nuc"];
  IntegerVector obs_index = dat_info["id"];
  IntegerVector ref_pos = dat_info["ref_pos"];
  IntegerVector read_pos = dat_info["read_pos"];
  IntegerVector n_t = hmm_info["n_t"];
  IntegerVector index = dat_info["start_id"];
  IntegerVector length = dat_info["length"];
  IntegerVector hap_len = hmm_info["p_tmax"];
  unsigned int t_max = hmm_info["t_max"];
  List n_in_t = hmm_info["n_in_t"];
  
  List dat_out(t_max);
  unsigned int id, len, total;
  for(t = 0; t < t_max; ++t) {
    len = 0;
    List dat_info_t(7);
    IntegerVector idx = n_in_t[t];
    for(i = 0; i < n_t[t]; ++i) {
      id = idx[i];
      len += length[id];
    }
    IntegerVector qua_out(len);
    IntegerVector obs_out(len);
    IntegerVector obs_index_out(len);
    IntegerVector ref_pos_out(len);
    IntegerVector read_pos_out(len);
    total = 0;
    for(i = 0; i < n_t[t]; ++i) {
      id = idx[i];
      for(j = 0; j < length[id]; ++j) {
        qua_out(total) =  qua[index[id] + j];
        obs_out(total) = obs[index[id] + j];
        obs_index_out(total) = obs_index[index[id] + j];
        ref_pos_out(total) = ref_pos[index[id] + j];
        read_pos_out(total++) = read_pos[index[id] + j];
      }
    }
    dat_info_t["total"] = total;
    dat_info_t["ref_length_max"] = hap_len[t];
    dat_info_t["qua"] = qua_out;
    dat_info_t["nuc"] = obs_out;
    dat_info_t["id"] = obs_index_out;
    dat_info_t["ref_pos"] = ref_pos_out;
    dat_info_t["read_pos"] = read_pos_out;
    dat_out(t) = dat_info_t;
  }
  return(dat_out);
}

/*
 * format the data for mnlogit, note here ref_pos should be shifted according the strating t
 */
// [[Rcpp::export]]
DataFrame format_data2(List hmm_info, List d_info, List hap_info) {
  IntegerVector num_states = hmm_info["num_states"];
  IntegerVector time_pos = hmm_info["time_pos"];
  unsigned int t_max = hmm_info["t_max"];
  unsigned int t, m, i;
  
  List subsample = sub_sample(hmm_info, d_info);
  
  unsigned int total = 0;
  for(t = 0; t < t_max; ++t) {
    List data_info = subsample(t);
    unsigned int tt = data_info["total"];
    total += num_states[t] * tt * MLOGIT_CLASS * NUM_CLASS;
  }
  
  IntegerVector r_ref_pos(total);
  IntegerVector r_read_pos(total);
  IntegerVector r_qua(total);
  IntegerVector r_obs(total);
  IntegerVector r_hap_nuc(total);
  IntegerVector r_mode(total);
  IntegerVector r_id(total);
  
  unsigned int count = 0;
  for(t = 0; t < t_max; ++t) {
    List full_hap_t = hap_info(t);
    List data_info = subsample(t);
    unsigned int num = data_info["total"];
    unsigned int len = num * MLOGIT_CLASS * NUM_CLASS;
    for(m = 0; m < num_states[t]; ++m) {
      IntegerMatrix haplotype = full_hap_t(m);
      DataFrame df = format_data(data_info, haplotype, time_pos[t]);
      IntegerVector ref_pos = df["ref_pos"];
      IntegerVector read_pos = df["read_pos"];
      IntegerVector qua = df["qua"];
      IntegerVector obs = df["nuc"];
      IntegerVector hap_nuc = df["hap_nuc"];
      IntegerVector mode = df["mode"];
      IntegerVector id = df["id"];
      for(i = 0; i < len; ++i) {
        r_ref_pos(count) = ref_pos(i);
        r_read_pos(count) = read_pos(i);
        r_qua(count) = qua(i);
        r_obs(count) = obs(i);
        r_hap_nuc(count) = hap_nuc(i);
        r_mode(count) = mode(i);
        r_id(count++) = id(i);
      }
    }
  }
  
  DataFrame df_new = DataFrame::create(
    Named("id") = r_id,
    Named("mode") = r_mode,
    Named("read_pos") = r_read_pos,
    Named("ref_pos") = r_ref_pos,
    Named("qua") = r_qua,
    Named("nuc") = r_obs,
    Named("hap_nuc") = r_hap_nuc);
  return(df_new);
}

NumericVector make_weight(List wic, List gamma, List hmm_info, List dat_info) {
  IntegerVector num_states = hmm_info["num_states"];
  IntegerVector n_t = hmm_info["n_t"];
  IntegerVector length = dat_info["length"];
  List n_in_t = hmm_info["n_in_t"];
  unsigned int t_max = hmm_info["t_max"];
  
  unsigned int t, m, i, k, j;
  unsigned int id;
  unsigned int num, count = 0;
  IntegerVector n_set;
  for(t = 0; t < t_max; ++t) {
    num = 0;
    n_set = n_in_t[t];
    for (i = 0; i < n_t[t]; ++i) { 
      id = n_set[i];
      num += length[id];
    }
    count += num * NUM_CLASS * num_states[t];
  }
  NumericVector weight(count);
  count = 0;
  for(t = 0; t < t_max; ++t) {
    List w_icm = wic(t);
    NumericVector gam = gamma(t);
    n_set = n_in_t[t];
    for(m = 0; m < num_states[t]; ++m) {
      NumericMatrix w_tic = w_icm(m);
      for (i = 0; i < n_t[t]; ++i) { // todo: maybe better to make a read set for each t
        id = n_set[i];
        for(j = 0; j < length[id]; ++j)
          for (k = 0; k < NUM_CLASS; ++k)
            weight(count++) = w_tic(i, k) * gam(m);
      }
    }
  }
  return(weight);
}

// [[Rcpp::export]]
List baum_welch_init(List hmm_info, List data_info, List hap_info, int PD_LENGTH, List par, 
                     Nullable<List> trans_indicator_new = R_NilValue)
{
  IntegerVector num_states = hmm_info["num_states"];
  IntegerVector n_t = hmm_info["n_t"];
  NumericMatrix beta = par["beta"];
  NumericVector eta = par["eta"];
  unsigned int t_max = hmm_info["t_max"];
  
  List par_hmm;
  /* initialize hmm par */
  par_hmm = ini_hmm(t_max, num_states, trans_indicator_new);
  NumericVector phi = par_hmm["phi"];
  List trans = par_hmm["trans"];
  List for_emit = compute_emit(hmm_info, data_info, hap_info, beta, eta, PD_LENGTH);
  par_hmm["emit"] = for_emit["emit"];
  List w_ic = for_emit["w_ic"];
  List par_hmm_bf = forward_backward(par_hmm, t_max, num_states);
  List gamma = par_hmm_bf["gamma"];
  /* prepare weight for beta */
  NumericVector weight = make_weight(w_ic, gamma, hmm_info, data_info);
  
  //store the parmaeters for calling mnlogit
  List par_aux = List::create(
    Named("beta") = beta,
    Named("w_ic") = w_ic,
    Named("weight") = weight);
  
  List ls = List::create(
    Named("par_hmm_bf") = par_hmm_bf,
    Named("par_aux") = par_aux,
    Named("par_hmm") = par_hmm);
  
  return(ls);
}
// [[Rcpp::export]]
List baum_welch_iter(List hmm_info, List par_hmm, List data_info, List hap_info, NumericMatrix beta, int PD_LENGTH)
{
  List par_aux = par_hmm["par_aux"];
  List par_hmm_bf = par_hmm["par_hmm_bf"];
  IntegerVector num_states = hmm_info["num_states"];
  IntegerVector n_t = hmm_info["n_t"];
  List w_ic = par_aux["w_ic"];
  List gamma = par_hmm_bf["gamma"];
  unsigned int t_max = hmm_info["t_max"];
  unsigned int n_observation = data_info["n_observation"];
  List par_hmm_new;
  /* update eta */
  NumericVector eta_new = update_eta(w_ic, gamma, num_states, n_t, t_max, n_observation);
  
  /* update emit */
  List for_emit_new = compute_emit(hmm_info, data_info, hap_info, beta, eta_new, PD_LENGTH);
  par_hmm_new["emit"] = for_emit_new["emit"];
  List wic_new = for_emit_new["w_ic"];
  
  /* update hmm par (except for emit) */
  NumericVector phi_new = update_phi(par_hmm_bf, num_states[0]);
  List trans_new = update_trans(par_hmm_bf, t_max, num_states);
  par_hmm_new["phi"] = phi_new;
  par_hmm_new["trans"] = trans_new;
  
  List par_hmm_bf_new = forward_backward(par_hmm_new, t_max, num_states);
  List gamma_new = par_hmm_bf_new["gamma"];
  
  /* prepare weight for beta */
  NumericVector weight = make_weight(wic_new, gamma_new, hmm_info, data_info);
  
  List par_aux_out = List::create(
    Named("beta") = beta,
    Named("eta") = eta_new,
    Named("w_ic") = wic_new,
    Named("weight") = weight);
  
  List ls = List::create(
    Named("par_aux") = par_aux_out,
    Named("par_hmm") = par_hmm_new,
    Named("par_hmm_bf") = par_hmm_bf_new);
  
  return(ls);
}

// find which states at t can transfer to the next states at t+1
// undecided_pos starts from 0; but time_pos starts from int hap_min_pos = dat_info["ref_start"];

// [[Rcpp::export]]
List trans_permit(IntegerVector num_states, List combination, int t_max, IntegerVector undecided_pos, 
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
    if(num_states[t + 1] != 1) {
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
          // Rcout << "\n";
        }
        // if for a row, all == 1 or for a column all == 1, then these two should be excluded
        IntegerVector more_limits_row(num_states[t]);
        IntegerVector more_limits_col(num_states[t + 1]);
        trans_permits(t) = trans;
        for(w = 0; w < num_states[t + 1]; ++w) {
          int num_col = 0;
          for(m = 0; m < num_states[t]; ++m) {
            if(trans(m, w) == 1)
              num_col++;
          }
          if(num_col == num_states[t]){
            more_limits_col[w] = 1; 
          }// meaning no states transfer to it, there must be some overlaps??? what if not?
        }
        for(m = 0; m < num_states[t]; ++m) {
          int num_row = 0;
          for(w = 0; w < num_states[t + 1]; ++w) {
            if(trans(m, w) == 1)
              num_row++;
          }
          if(num_row == num_states[t + 1]){
            more_limits_row[m] = 1; 
          }
        }
        // found the states that never appears
        further_limit_col(t) = more_limits_col;
        further_limit_row(t) = more_limits_row;
    } else {
      IntegerVector more_limits_col(num_states[t + 1]);
      IntegerVector more_limits_row(num_states[t]);
      IntegerMatrix temp(num_states[t], num_states[t + 1]);
      trans_permits(t) = temp;
      further_limit_col(t) = more_limits_col;
      further_limit_row(t) = more_limits_row;
    }
  }
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

// [[Rcpp::export]]
List final_exclude(List full_hap, List further_limit, int t_max, IntegerVector num_states) {
  unsigned int t, m;
  List full_hap_new(t_max);
  
  for (t = 0; t < t_max; ++t) {
    List full_hap_t = full_hap[t];
    List full_hap_t_new(num_states[t]);
    IntegerVector more_limits = further_limit[t];
    int count = 0;
    for (m = 0; m < more_limits.size(); ++m) {
      if (!more_limits[m]) {
        IntegerMatrix haplotype = full_hap_t(m);
        full_hap_t_new[count++] = haplotype;
      }
    }
    full_hap_new[t] = full_hap_t_new;
  }
  
  return(full_hap_new);
}







