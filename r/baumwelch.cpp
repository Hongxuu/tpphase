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
using namespace Rcpp;
using namespace std;

#define MLOGIT_CLASS 4
#define NUM_CLASS 4

// [[Rcpp::depends(RcppArmadillo)]]


List ini_hmm (unsigned int t_max,  IntegerVector num_states) {
  NumericVector phi(num_states[0]);
  List trans(t_max - 1);
  // List emit(t_max);
  unsigned int w, m, t;
  
  for(w = 0; w < num_states[0]; ++w)
    phi[w] = 1/double(num_states[0]);
  
  for(t = 0; t < t_max - 1; ++t) {
    NumericMatrix transition(num_states[t], num_states[t + 1]);
    for (w = 0; w < num_states[t]; ++w) {
      for (m = 0; m < num_states[t + 1]; ++m)
        transition(w, m) = 1/double(num_states[t + 1]);
    }
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
    if(t == 0)
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
        x(w, m) = log(alp(w)) + log(transition(w, m)) + log(beta(m)) + log(emission(m)); // avoid underflow?
        x(w, m) = exp(x(w, m));
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
    ref_pos_in = ref_pos[index[i] + j] - time_pos; // this is different since the sub-hap are counted 0 from the time t
    
    for (l = 0; l < MLOGIT_CLASS; ++l) {
      hap_nuc[l] = 0;
      hnuc_qua[l] = 0;
    }
    //Rcout << "haplotype " << haplotype(K, ref_pos_in) << "\t";
    
    if(haplotype(K, ref_pos_in) == 1) {
      hap_nuc[0] = 1;
      hnuc_qua[0] = qua_in;
    } else if(haplotype(K, ref_pos_in) == 3) {
      hap_nuc[1] = 1;
      hnuc_qua[1] = qua_in;
    } else if(haplotype(K, ref_pos_in) == 2) {
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
 * Fill the positions without variation
 */
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
 * Return the combination of sites with variation within one time t
 */
List find_combination(IntegerVector undecided_pos, IntegerVector pos_possibility, 
                      unsigned int p_tmax, unsigned int time_pos) {
  //possible combination of the rest non-unique loci
  IntegerVector location(pos_possibility.size());
  IntegerVector location_len(pos_possibility.size());
  unsigned int num = 0;
  for(unsigned int i = 0; i < pos_possibility.size(); ++i)
    if(time_pos <= undecided_pos[i] && undecided_pos[i] < time_pos + p_tmax) {
      location(num) = undecided_pos[i];
      location_len(num++) = pos_possibility[i];
    }
    IntegerMatrix combination = call_cart_product(location_len[Range(0, num - 1)]);
    List ls = List::create(
      Named("combination") = combination,
      Named("num") = num,
      Named("location") = location[Range(0, num - 1)]);
    return(ls);
}

//fill haplotype at the rest sites (with variation)
IntegerMatrix make_hap(List hidden_states, IntegerMatrix haplotype, IntegerVector location, unsigned int p_tmax,
                       IntegerVector combination, unsigned int time_pos, unsigned int num) {
  unsigned int j, k, idx;
  
  for(j = 0; j < num; ++j) {
    IntegerMatrix hidden = hidden_states[location[j]];
    //Rcout << hidden << "\n";
    idx = combination[j];
    //Rcout << idx << "\n";
    for(k = 0; k < NUM_CLASS; ++k)
      haplotype(k, location[j]) = hidden(idx, k);
  }
  return(haplotype(_, Range(time_pos, time_pos + p_tmax - 1)));
}

// [[Rcpp::export]]
List full_hap(List hmm_info, unsigned int hap_length) {
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
  IntegerMatrix hap = fill_all_hap(hidden_states, hap_length, n_row);
  for(t = 0; t < t_max; ++t) {
    List full_hap_t(num_states(t));
    // give h_t, each t has many possible combinations
    if(num_states[t] != 1) {
      List comb_info = find_combination(undecided_pos, pos_possibility, p_tmax[t], time_pos[t]);
      IntegerVector location = comb_info["location"];
      //Rcout << location << "\n";
      IntegerMatrix combination = comb_info["combination"];
      unsigned int num = comb_info["num"];
      for(m = 0; m < num_states[t]; ++m) {
        IntegerMatrix haplotype = make_hap(hidden_states, hap, location, p_tmax[t], combination(m, _), time_pos[t], num);
        full_hap_t(m) = haplotype;
      }
    } else {
      full_hap_t[0] = hap(_, Range(time_pos[t], time_pos[t] + p_tmax[t] - 1));
    }
    full_hap(t) = full_hap_t;
  }
  return(full_hap);
}

List compute_emit(List hmm_info, List dat_info, List hap_info, NumericMatrix beta, NumericVector eta, int PD_LENGTH) 
{
  List hidden_states = hmm_info["hidden_states"];
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
  for(t = 0; t < t_max; ++t){
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
      for (i = 0; i < n_t[t]; ++i) { // todo: maybe better to make a read set for each t
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
    List dat_info_t;
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
    unsigned int num = data_info["total"];
    total += num_states[t] * num * MLOGIT_CLASS * NUM_CLASS;
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
      //Rcout << m << "\n";
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
List baum_welch_init(List hmm_info, List data_info, List hap_info, List par, int PD_LENGTH)
{
  IntegerVector num_states = hmm_info["num_states"];
  IntegerVector n_t = hmm_info["n_t"];
  NumericMatrix beta = par["beta"];
  NumericVector eta = par["eta"];
  unsigned int t_max = hmm_info["t_max"];
  
  List par_hmm;
  
  /* initialize hmm par */
  par_hmm = ini_hmm(t_max, num_states);
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