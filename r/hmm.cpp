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

#define MLOGIT_CLASS 4
#define NUM_CLASS 4

// [[Rcpp::depends(RcppArmadillo)]]

enum {	NO_ERROR /*!< no error */ };
int foward_backward(List par_hmm, unsigned int t_max,  IntegerVector num_states);
int compute_emit(IntegerVector n_t, NumericVector read_likelihood, 
                 unsigned int t_max,  IntegerVector num_states);
int compute_trans(List par_hmm, unsigned int t_max,  IntegerVector num_states);

List ini_hmm (unsigned int t_max,  IntegerVector num_states) {
  NumericVector phi(num_states[0], 1/num_states[0]);
  List trans(t_max - 1);
  List emit(t_max);
  unsigned int w, m, t;
  for(t = 0; t < t_max - 1; ++t) {
    NumericMatrix transition(num_states[t], num_states[t + 1]);
    for (w = 0; w < num_states[t]; ++w) {
      for (m = 0; m < num_states[t + 1]; ++m)
        transition(w, m) = 1/num_states[t + 1];
    }
    trans(t) = transition;
  }
  for(t = 0; t < t_max; ++t) {
    NumericVector emission(num_states[t]);
    for (w = 0; w < num_states[t]; ++w)
      emission(w) = 1/num_states[t];
    emit(t) = emission;
  }
  
  List par_hmm = List::create(
    Named("phi") = phi,
    Named("trans") = trans,
    Named("emit") = emit);
  return(par_hmm);
}

List foward_backward(List par_hmm, unsigned int t_max, IntegerVector num_states)
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
    NumericVector transition = trans(t);
    NumericVector alp_last;
    if(t == 0)
      for (w = 0; w < num_states[t]; ++w)
        alp(w) = phi[w] * emission(w);
    else {
      alp_last = alpha(t - 1);
      for (w = 0; w < num_states[t]; ++w) {
        for (m = 0; m < num_states[t - 1]; ++m)
          alp(w) += alp_last(m) * transition(m, w);
        alp(w) = alp(w) * emission(w, t);
      }
    }
  }
  alpha(t) = alp;
}
for (t = t_max - 1; t >= 0; --t) {
  NumericVector beta(num_states[t]);
  NumericVector emission = emit(t);
  NumericVector transition = trans(t);
  NumericVector beta_after;
  if(t == t_max - 1)
    for (w = 0; w < num_states[t]; ++w)
      beta(w) = 1;
  else {
    beta_after = beta_wt(t + 1);
    for (w = 0; w < num_states[t]; ++w) {
      for (m = 0; m < num_states[t + 1]; ++m) 
        beta(w) += beta_after(w, t + 1) * transition(w, m) * emission(m, t + 1);
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
  List xi(t_max);
  //arma::Cube<double> xi(total_state, total_state, t_max);
  for (t = 0; t < t_max - 1; ++t) {
    NumericVector emission = emit(t);
    NumericVector transition = trans(t);
    NumericVector beta = beta_wt(t);
    NumericVector alp = alpha(t + 1);
    NumericMatrix x(num_states[t], num_states[t + 1]);
    sum = 0;
    for (w = 0; w < num_states[t]; ++w)
      for (m = 0; m < num_states[t + 1]; ++m) {
        x(w, m) = alp(t) * transition(w, m) * emission(m, t + 1) * emission(m, t + 1);
        sum += x(w, m);
      }
      for (w = 0; w < num_states[t]; ++w)
        for (m = 0; m < num_states[t + 1]; ++m)
          x(w, m) /= sum;
    xi(t) = x;
  }
  List par_hmm = List::create(
    Named("phi") = phi,
    Named("trans") = trans,
    Named("emit") = emit,
    Named("alpha") = alpha,
    Named("beta_wt") = beta_wt,
    Named("gamma") = gamma,
    Named("xi") = xi);
  return(par_hmm);
}
NumericVector update_phi(List par_hmm, unsigned int num_states) {
  List gamma = par_hmm["gamma"];
  NumericVector phi(num_states);
  NumericVector gam = gamma[0];
  
  for (w = 0; w < num_states; ++w)
    phi[w] = gam(w);
  return phi;
}

List update_trans(List par_hmm, unsigned int t_max,  IntegerVector num_states) {
  List gamma = par_hmm["gamma"];
  List xi = par_hmm["xi"];
  List trans(t_max);
  
  for (t = 0; t < t_max; ++t) {
    NumericMatrix transition(num_states[t], num_states[t + 1]);
    NumericVector gam = gamma[t];
    NumericVector x = xi[t];
    for (w = 0; w < num_states[t]; ++w)
      for (m = 0; m < num_states[t + 1]; ++m) {
        transition(w, m) = x(w, m)/gam(w);
      }
      trans[t] = transition;
  }
  return(trans);
}

NumericMatrix update_emit(unsigned int n_observation, NumericVector read_likelihood, 
                          unsigned int t_max,  IntegerVector num_states, IntegerVector which_t) {
  unsigned int t, w, i;
  List emit(t_max);
  for (t = 0; t < t_max; ++t) {
    NumericVector emission(num_states[t]);
    for (w = 0; w < num_states[t]; ++w) {
      for (i = 0; i < n_observation; ++i) {
        if(which_t[i] == t)
          emission(w) *= read_likelihood[i]; //no indel penalty now
      }
    }
    emit(t) = emission;
  }
  return(emit);
}

double site_likelihood (unsigned int i, unsigned int K, NumericMatrix beta, unsigned int PD_LENGTH, double ins_rate, double del_rate, 
                        List dat_info, IntegerMatrix haplotype, unsigned int del_count, unsigned int ins_count)
{
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
    if(haplotype(K, ref_pos_in) == 4)
      continue;
    qua_in = qua[index[i] + j];
    ref_pos_in = ref_pos[index[i] + j];
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
  /* relative deletion penalty */
  int hap_length = dat_info["ref_length_max"];
  if(del_count != 0)
    read_llk += R::dpois(del_count, hap_length * del_rate, true);
  //Rcout << "deletion llk" << R::dpois(del_count, hap_length * del_rate, true) << "\t";
  /* relative insertion penalty */
  if(ins_count != 0)
    read_llk += R::dpois(ins_count, hap_length * ins_rate, true);
  //Rcout << "ins llk" << R::dpois(ins_count, hap_length * ins_rate, true) << "\n";
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

List full_hap(List hmm_info, unsigned int hap_length) {
  List hidden_states = hmm_info["hidden_states"];
  IntegerVector num_states = hmm_info["num_states"];
  IntegerVector time_pos = hmm_info["time_pos"];
  IntegerVector p_tmax = hmm_info["p_tmax"];
  IntegerVector n_row = hmm_info["n_row"];
  IntegerVector pos_possibility = hmm_info["pos_possibility"];
  IntegerVector undecided_pos = hmm_info["undecided_pos"];
  unsigned int t_max = 1;
  unsigned int t, m;
  List full_hap(t_max);
  IntegerMatrix haplotype = fill_all_hap(hidden_states, hap_length, n_row);
  for(t = 0; t < t_max; ++t){
    List full_hap_t(num_states(t));
    IntegerMatrix hap_t = haplotype(_, Range(time_pos[t], time_pos[t] + p_tmax[t] - 1));
    // give h_t, each t has many possible combinations
    List comb_info = find_combination(undecided_pos, pos_possibility, p_tmax[t], time_pos[t]);
    IntegerVector location = comb_info["location"];
    IntegerMatrix combination = comb_info["combination"];
    unsigned int num = comb_info["num"];
    for(m = 0; m < num_states[t]; ++m) {
      IntegerMatrix haplotype = make_hap(hidden_states, hap_t, location, p_tmax[t], combination(m, _), time_pos[t], num);
      full_hap_t(m) = haplotype;
    }
    full_hap(t) = full_hap_t;
  }
  return full_hap;
}

List baum_welch(List hmm_info, List data_info，List par)
{
  List hidden_states = hmm_info["hidden_states"];
  IntegerVector num_states = hmm_info["num_states"];
  IntegerVector n_t = hmm_info["n_t"];
  IntegerVector which_t = hmm_info["which_t"];
  IntegerVector time_pos = hmm_info["time_pos"];
  IntegerVector p_tmax = hmm_info["p_tmax"];
  unsigned int t_max = time_pos.length();
  unsigned int n_observation = dat_info["n_observation"];
  unsigned int i, k, l, m, t;
  
  List par_hmm(7);
  /* initialize hmm par */
  par_hmm = ini_hmm(t_max, num_states);
  
  NumericVector eta = par["eta"];
  NumericVector phi = par_hmm["phi"];
  List trans = par_hmm["trans"];
  List emit = par_hmm["emit"];
  
  List w_tic(t_max);
  NumericVector mixture_prop_updated(NUM_CLASS);
  NumericMatrix read_class_llk(n_observation, NUM_CLASS);
  NumericMatrix read_likelihood(n_observation);
  // compute likelihood for emission prob (beta should be estimated by mnlogit outside cpp)
  for(t = 0; t < t_max; ++t){
    // give h_t, each t has many possible combinations
    for(m = 0; m < num_states[t]; ++m) {
      NumericMatrix haplotype = ;
      NumericMatrix w_ic(n_t[t], NUM_CLASS);
      for (i = 0; i < n_observation; ++i) {
        if(which_t[i] == t)
          for (k = 0; k < NUM_CLASS; ++k) {
            read_class_llk(i, k) = site_likelihood(i, k, beta, PD_LENGTH, ins_rate, del_rate, dat_info, 
                           haplotype, del_count(i, k), ins_count(i, k));
            weight_llk[k] = eta[k] * exp(read_class_llk(i, k));
            read_likelihood[i] += weight_llk[k];
          }
          for (k = 0; k < NUM_CLASS; ++k)
            w_ic(i, k) = weight_llk[k]/read_likelihood[i];
      }
    }
    
    w_tic(t) = w_ic(i, k);
  }
  
  par_hmm = foward_backward(par_hmm, t_max, num_states);
  trans = update_trans(par_hmm, t_max, num_states);
  emit = update_emit(n_t, read_likelihood, t_max, num_states);
  
  // update eta where P(Ht|R, theta) = gamma
}


//create list t of list haplotypes segments



