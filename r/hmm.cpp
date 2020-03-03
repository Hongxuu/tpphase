#include <RcppArmadillo.h>
#include <stdio.h>
#include <stdlib.h>

#include <stddef.h>
#include <ctype.h>
#include <limits.h>
#include <unistd.h>
#include <math.h>
using namespace Rcpp;

#define MLOGIT_CLASS 4
#define NUM_CLASS 4

// [[Rcpp::depends(RcppArmadillo)]]

enum {	NO_ERROR /*!< no error */ };
int foward_backward(List par_hmm, unsigned int t_max, unsigned int num_states);
int compute_emit(IntegerVector n_t, NumericVector read_likelihood, 
                 unsigned int t_max, unsigned int num_states);
int compute_trans(List par_hmm, unsigned int t_max, unsigned int num_states);

List ini_hmm (unsigned int t_max, unsigned int num_states) {
  NumericVector phi(num_states, 1/num_states);
  NumericMatrix trans(num_states, num_states);
  NumericMatrix emit(num_states, t_max);
  unsigned int w, m, t;
  for (w = 0; w < num_states; ++w) {
    for (m = 0; m < num_states; ++m) {
      trans(w, m) = 1/num_states;
    }
    for (t = 0; t < t_max; ++t) {
      emit(w, t) = 1/t_max;
    }
  }
  List par_hmm = List::create(
    Named("phi") = phi,
    Named("trans") = trans,
    Named("emit") = emit);
  return(par_hmm);
}

List foward_backward(List par_hmm, unsigned int t_max, unsigned int num_states)
{
  NumericVector phi = par_hmm["phi"];
  NumericMatrix trans = par_hmm["trans"];
  NumericMatrix emit = par_hmm["emit"];
  NumericMatrix alpha(num_states, t_max);
  NumericMatrix beta_wt(num_states, t_max);
  unsigned int w, t, m;
  for (w = 0; w < num_states; ++w) {
    alpha(w, 0) = phi[w] * emit(w, 0);
    for (t = 1; t < t_max; ++t) {
      for (m = 0; m < num_states; ++m) {
        alpha(w, t) += alpha(m, t - 1) * trans(m, w);
      }
      alpha(w, t) = emit(w, t);
    }
  }
  for (w = 0; w < num_states; ++w) {
    beta_wt(w, t_max - 1) = 1;
    for (t = t_max - 2; t >= 0; --t) {
      for (m = 0; m < num_states; ++m) {
        beta_wt(w, t) += beta_wt(w, t + 1) * trans(w, m) * emit(m, t + 1);
      }
    }
  }
  NumericMatrix gamma(num_states, t_max);
  arma::Cube<double> xi(num_states, num_states, t_max);
  double sum = 0;
  for (t = 0; t < t_max; ++t) {
    sum = 0;
    for (w = 0; w < num_states; ++w) {
      gamma(w, t) = alpha(w, t) * beta_wt(w, t);
      sum += gamma(w, t);
    }
    for (w = 0; w < num_states; ++w)
      gamma(w, t) /= sum;
    
  }
  for (t = 0; t < t_max; ++t) {
    sum = 0;
    for (w = 0; w < num_states; ++w)
      for (m = 0; m < num_states; ++m) {
        xi(w, m, t) = alpha(w, t) * trans(w, m) * emit(m, t + 1) * emit(m, t + 1);
        sum += xi(w, m, t);
      }
      for (w = 0; w < num_states; ++w)
        for (m = 0; m < num_states; ++m)
          xi(w, m, t) /= sum;
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

NumericMatrix update_trans(List par_hmm, unsigned int t_max, unsigned int num_states) {
  NumericMatrix gamma = par_hmm["gamma"];
  arma::Cube<double> xi = par_hmm["xi"];
  NumericMatrix trans(num_states, num_states);
  double sum_de, sum_nu;
  for (w = 0; w < num_states; ++w)
    phi[w] = gamma(w, 0);
  
  for (w = 0; w < num_states; ++w)
    for (m = 0; m < num_states; ++m) {
      sum_de = sum_nu = 0;
      for (t = 0; t < t_max; ++t) {
        sum_nu += xi(w, m, t);
        sum_de += gamma(w, t);
      }
      trans(w, m) = sum_nu/sum_de;
    }
  
  return(trans);
}

NumericMatrix update_emit(IntegerVector n_t, NumericVector read_likelihood, 
                unsigned int t_max, unsigned int num_states) {
  unsigned int t, w, i;
  NumericMatrix emit(num_states, t_max);
  for (t = 0; t < t_max; ++t) {
    for (w = 0; w < num_states; ++w) {
      for (i = 0; i < n_t[t]; ++i) {
        emit(w, t) *= read_likelihood[i]; //no indel penalty now
      }
    }
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

List baum_welch(List hmm_info, List data_info， List par)
{
  IntegerVector n_t = hmm_info["n_t"];
  IntegerVector time_pos = hmm_info["time_pos"];
  IntegerVector p_tmax = hmm_info["p_tmax"];
  int n_observation = dat_info["n_observation"];
  unsigned int t_max = time_pos.length();
  unsigned int num_states = hmm_info["num_states"];
  unsigned int i, k, l, m, t;
  
  List par_hmm(7);
  NumericMatrix emit(num_states, t_max);
  NumericMatrix trans(num_states, num_states);
  
  /* initialize hmm par */
  par_hmm = ini_hmm(t_max, num_states);
  
  NumericVector eta = par["eta"];
  NumericVector phi = par_hmm["phi"];
  NumericMatrix trans = par_hmm["trans"];
  NumericMatrix emit = par_hmm["emit"];
  
  List w_tic(t_max);
  NumericVector mixture_prop_updated(NUM_CLASS);
  NumericMatrix read_class_llk(n_observation, NUM_CLASS);
  NumericMatrix read_likelihood(n_observation);
  // estimate eta for emission prob (beta should be estimated by mnlogit outside cpp)
  
  for(t = 0; t < t_max; ++t){
    // give h_t
    NumericMatrix haplotype_t(NUM_CLASS, p_tmax[t]);
    NumericMatrix w_ic(n_t[t], NUM_CLASS);
    for (i = 0; i < n_t[t]; ++i) {
      for (k = 0; k < NUM_CLASS; ++k) {
        read_class_llk(i, k) = site_likelihood(i, k, beta, PD_LENGTH, ins_rate, del_rate, dat_info, 
                       haplotype_t, del_count(i, k), ins_count(i, k));
        weight_llk[k] = eta[k] * exp(read_class_llk(i, k));
        read_likelihood[i] += weight_llk[k];
      }
      for (k = 0; k < NUM_CLASS; ++k)
        w_ic(i, k) = weight_llk[k]/read_likelihood[i];
    }
    w_tic(t) = w_ic(i, k);
  }
  //compute P(H|R, theta) or with t?
  unsigned int total; // number of possible combination of H
  arma::Cube<double> H_marginal(total, );
  for (m = 0; m < total; ++m) {
    //for each m, we have a H, then we know the order of num_states
    H_marginal(m) = phi(w) * 
  
  }
  
  
  par_hmm = foward_backward(par_hmm, t_max, num_states);
  trans = update_trans(par_hmm, t_max, num_states);
  emit = update_emit(n_t, read_likelihood, t_max, num_states);
  
  
}







