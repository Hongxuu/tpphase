#include <RcppArmadillo.h>
#include <stdio.h>
#include <stdlib.h>

#include <stddef.h>
#include <ctype.h>
#include <limits.h>
#include <unistd.h>
#include <math.h>

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

#define MLOGIT_CLASS 4
#define NUM_CLASS 4

double site_likelihood (unsigned int i, unsigned int K, NumericMatrix beta, unsigned int PD_LENGTH, double rate, 
                        List dat_info, IntegerMatrix haplotype, unsigned int del_count, unsigned int ins_count);

double site_likelihood (unsigned int i, unsigned int K, NumericMatrix beta, unsigned int PD_LENGTH, double rate, 
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
//if(!N_in) {
    predictor = {1, read_pos_in, ref_pos_in, qua_in, hap_nuc[0], hap_nuc[1], 
                 hap_nuc[3], hnuc_qua[0], hnuc_qua[1], hnuc_qua[3]};
    // } else {
    //   predictor = {1, read_pos_in, ref_pos_in, qua_in, hap_nuc[0], hap_nuc[1], 
    //                          hap_nuc[2], hap_nuc[3], hnuc_qua[0], hnuc_qua[1], hnuc_qua[2], hnuc_qua[3]};
    // }
    
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
    //Rcout << "xb + tail: " << xb + tail << "\t";
    read_llk += xb + tail; /* Notice here use log likelihood, not likelihood */
  }
  
  /* relative deletion penalty */
  int hap_length = dat_info["ref_length_max"];
  if(del_count != 0)
    read_llk += R::dpois(del_count, hap_length * rate_del, true);
    //Rcout << "deletion llk" << R::dpois(del_count, hap_length * rate, true) << "\t";
  /* relative insertion penalty */
  if(ins_count != 0)
    read_llk += R::dpois(ins_count, hap_length * rate_ins, true);
  
  return read_llk;
} /* site_likelihood */


// [[Rcpp::export]]
List em_eta (List par, List dat_info, List hap_info, IntegerMatrix haplotype,
             unsigned int PD_LENGTH) {
  
  int n_observation = dat_info["n_observation"];
  int hap_length = dat_info["ref_length_max"];
  double full_llk;
  unsigned int i, k, l, m;
  
  NumericMatrix w_ic(n_observation, NUM_CLASS);
  NumericVector eta = par["eta"];
  NumericVector mixture_prop_updated(NUM_CLASS);
  IntegerVector excluded_read = par["excluded_read"];
  NumericMatrix read_class_llk(n_observation, NUM_CLASS);
  NumericMatrix sum_weight(n_observation);
  NumericMatrix beta = par["beta"];
  double del_rate = par["del_rate"];
  double ins_rate = par["ins_rate"];
  // NumericVector m_hap_llk(NUM_CLASS * hap_length * MLOGIT_CLASS);
  // m_hap_llk.attr("dim") = Dimension(NUM_CLASS, hap_length, MLOGIT_CLASS);
  
  /* deletion */
  List deletion = dat_info["deletion"];
  IntegerVector del_ref_pos = deletion["del_ref_pos"];
  IntegerVector del_length_all = deletion["del_length_all"];
  //IntegerVector del_length = deletion["del_length"];
  IntegerVector del_flag = deletion["del_flag"];
  IntegerVector del_strat_id = deletion["del_strat_id"];
  IntegerVector hap_deletion_len = hap_info["hap_deletion_len"];
  IntegerVector hap_deletion_pos = hap_info["deletion_pos"];
  IntegerVector hap_del_start_id = hap_info["hap_del_start_id"];
  IntegerMatrix del_count(n_observation, NUM_CLASS);
  IntegerMatrix ins_count(n_observation, NUM_CLASS);
  IntegerMatrix indel_both(n_observation, NUM_CLASS);
  
  /* e step */
  for (i = 0; i < n_observation; ++i) {
    // Exclude the read has -inf likelihood
    if (excluded_read[i] == 1)
      continue;
    
    IntegerVector read_del_pos(del_length_all[i]);
    if (del_flag[i])
      for (l = 0; l < del_length_all[i]; ++l) 
        read_del_pos[l] = del_ref_pos[del_strat_id[i] + l];
    
    //Rcout << "del_length: " << del_length_all[i] << "read_del_pos: " << read_del_pos << "\n";
    /* reset eta %*% llk at each class for each read */
    NumericVector weight_llk(NUM_CLASS);
    /* compute the likelihood for each read under each haplotype */
    //sum_weight = 0.;
    for (k = 0; k < NUM_CLASS; ++k) {
      // find how many deletion positions are the NOT same between read and haplotype 
      // (both read and haplotype need to contain deletion)
      if (del_flag[i]) {
        if (hap_del_start_id[k] != -1) {
          IntegerVector hap_del_pos(hap_deletion_len[k]);
          for (l = 0; l < hap_deletion_len[k]; ++l) 
            hap_del_pos[l] = hap_deletion_pos[hap_del_start_id[k] + l];
          //Rcout<< "del_length: " << hap_deletion_len[k] << "hap_del_pos: " << hap_del_pos << "\n";
          for (m = 0; m < del_length_all[i]; ++m)
            for (l = 0; l < hap_deletion_len[k]; ++l) 
              if (read_del_pos[m] == hap_del_pos[l])
                indel_both(i, k)++;
          del_count(i, k) = del_length_all[i] - indel_both(i, k);
          ins_count(i, k) = hap_deletion_len[k] - indel_both(i, k);
        } 
      }
      read_class_llk(i, k) = site_likelihood(i, k, beta, PD_LENGTH, rate, dat_info, haplotype, del_count(i, k), ins_count(i, k));
      weight_llk[k] = eta[k] * exp(read_class_llk(i, k));
      //Rprintf("\n read %d class %d llk %f; log wei_lk %.10lf; del count %d\n", i, k, read_class_llk(i, k), log(weight_llk[k]), del_count(i, k));
      sum_weight[i] += weight_llk[k];
    }
    //Rprintf("sum_weight %.30lf\n", sum_weight[i]);
    for (k = 0; k < NUM_CLASS; ++k)
      w_ic(i, k) = weight_llk[k]/sum_weight[i];
    
    /* some reads may not be explained by either of the haplotype */
    if (exp(read_class_llk(i, 0)) == 0 || isnan(w_ic(i, 0)))
      excluded_read[i] = 1;
  }
  //Rcout << "w_ic : " << w_ic << "\n";
  //Rcout << "read class llk : " << read_class_llk << "\n";
  
  full_llk = 0.;
  int count = 0;
  //double mixture_llk = 0.;
  /* Record the full likelihood for terminating EM
   Some of the likelihood is 0, exclude those reads for now */
  for (i = 0; i < n_observation; ++i) {
    //mixture_llk = 0.;
    if (excluded_read[i] == 1)
    {
      count++;
      continue;
    }
    full_llk += log(sum_weight[i]);
  }
  
  /* update eta, still exclude nan */
  for (i = 0; i < n_observation; ++i) {
      if (excluded_read[i] == 1)
        continue;
      for (k = 0; k < NUM_CLASS; ++k)
      mixture_prop_updated[k] += w_ic(i, k);
  }
  /* assume the excluded reads do not come from any of the haplotypes */
  for (k = 0; k < NUM_CLASS; ++k)
    mixture_prop_updated[k] = mixture_prop_updated[k]/(n_observation - count);
  
  /* update rate */
  double wic_di = 0, wic_pc = 0, wic_ii = 0;
  for (i = 0; i < n_observation; ++i) {
    if (excluded_read[i] == 1)
      continue;
    for (k = 0; k < NUM_CLASS; ++k) {
      if(del_count(i, k) != 0)
        wic_di += w_ic(i, k) * del_count(i, k);
      if(ins_count(i, k) != 0)
        wic_ii += w_ic(i, k) * ins_count(i, k);
      wic_pc += w_ic(i, k) * hap_length;
    }
  }
  del_rate = wic_di/wic_pc;
  ins_rate = wic_ii/wic_pc;
  
  k = 0;
  IntegerVector excluded_id(count);
  for (i = 0; i < n_observation; ++i)
    if (excluded_read[i] == 1)
      excluded_id[k++] = i + 1; // Plus one for the use of R (index starts from 1)
  
  Rcout << "full loglik : " << full_llk << "\n";
  Rcout << "Mixture proportions : " << mixture_prop_updated << "\n";

  List ls_par = List::create(
      Named("mixture_prop") = mixture_prop_updated,
      Named("w_ic") = w_ic, 
      Named("beta") = beta,
      Named("del_rate") = del_rate,
      Named("ins_rate") = ins_rate,
      Named("excluded_read") = excluded_read);
  
  List ls = List::create(
    Named("full_llk") = full_llk,
    Named("param") = ls_par,
    Named("excluded_id") = excluded_id);
  
  return ls;
}
