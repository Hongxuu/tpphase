#include <RcppArmadillo.h>
#include <stdio.h>
#include <stdlib.h>

#include <stddef.h>
#include <ctype.h>
#include <limits.h>
#include <unistd.h>
#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

#define MLOGIT_CLASS 4
#define NUM_CLASS 4
#define NONINDEL_STATE 4
#define HIDDEN_STATE 2

// A C G T N: 0 1 2 3 4

double m_haplotype_llk (unsigned int j, unsigned int l, unsigned int K, unsigned int PD_LENGTH, NumericVector beta,
                        List dat_info, NumericMatrix w_ic, IntegerVector excluded_read, double ins_lambda, double del_lambda,
                        unsigned int viterbi);
  
double m_haplotype_llk (unsigned int j, unsigned int l, unsigned int K, unsigned int PD_LENGTH, NumericVector beta, 
                        List dat_info, NumericMatrix w_ic, IntegerVector excluded_read, double ins_lambda, double del_lambda,
                        unsigned int viterbi)
{
  IntegerMatrix ref_index = dat_info["ref_idx"];
  IntegerVector qua = dat_info["qua"];
  IntegerVector obs = dat_info["nuc"];
  IntegerVector ref_pos = dat_info["ref_pos"];
  IntegerVector read_pos = dat_info["read_pos"];
  IntegerVector length = dat_info["length"];
  IntegerVector index = dat_info["start_id"];
  int n_observation = dat_info["n_observation"];
  
  NumericVector hap_nuc(MLOGIT_CLASS);
  NumericVector hnuc_qua(MLOGIT_CLASS);
  NumericVector pred_beta(MLOGIT_CLASS - 1);
  arma::vec predictor(PD_LENGTH);
  double qua_in, read_pos_in, ref_pos_in;
  double tail, sum;
  double xb = 0;
  double likelihood = 0;
  unsigned int i, k, c;
  
  for (i = 0; i < n_observation; ++i) {
    //Rprintf("idx position %d\n", index[i] + ref_index(i, j));
    if (excluded_read[i] == 1)
      continue;
    
    if (ref_index(i, j) != -1) {
      qua_in = qua[index[i] + ref_index(i, j)];
      read_pos_in = read_pos[index[i] + ref_index(i, j)];
      ref_pos_in = ref_pos[index[i] + ref_index(i, j)];
      
      for (k = 0; k < MLOGIT_CLASS; ++k) {
        hap_nuc[k] = 0;
        hnuc_qua[k] = 0;
      }
      
      if (l == 1) {
        hap_nuc[0] = 1;
        hnuc_qua[0] = qua_in;
      } else if (l == 3) {
        hap_nuc[1] = 1;
        hnuc_qua[1] = qua_in;
      } else if (l == 2) {
        hap_nuc[3] = 1;
        hnuc_qua[3] = qua_in;
      }
      
      predictor = {1, read_pos_in, ref_pos_in, qua_in, hap_nuc[0], hap_nuc[1], 
                   hap_nuc[3], hnuc_qua[0], hnuc_qua[1], hnuc_qua[3]};
      
      arma::mat beta_ar = as<arma::mat>(beta);
      arma::vec pb = beta_ar.t() * predictor;
      pred_beta = as<NumericVector>(wrap(pb)); 
      
      //Rcout << "predictor : " << predictor.t() << "\n";
      // Rcout << "pred_beta : " << pb << "\t";
      sum = 0.0;
      for (c = 0; c < MLOGIT_CLASS - 1; ++c)
        sum += exp(pred_beta[c]);
      tail = log(1/(1 + sum));
      
      if(obs[index[i] + ref_index(i, j)] == 1) {
        xb = pred_beta[0];
      } else if(obs[index[i] + ref_index(i, j)] == 3) {
        xb = pred_beta[1];
      } else if(obs[index[i] + ref_index(i, j)] == 2) {
        xb = pred_beta[2];
      } else if(obs[index[i] + ref_index(i, j)] == 0) {
        xb = 0;
      }
      if (viterbi)
        likelihood += log(w_ic(i, K)) + xb + tail;
      else
        likelihood += w_ic(i, K) * exp(xb + tail);
      // if(i == 0)
      //   Rcout << "normal add " << w_ic(i, K) * exp(xb + tail) << "\t";
    }
    
      if(ref_index(i, j) == -1) {
        if((j >= ref_pos[index[i]]) && (j <= ref_pos[index[i] + length[i] - 1])) {
          if (viterbi)
            likelihood += log(w_ic(i, K)) + R::dpois(1, del_lambda, true);
          else
            likelihood += w_ic(i, K) * R::dpois(1, del_lambda, false);
        }
      }
      // if(i == 0)
      //   Rcout << "del add " << w_ic(i, K) * R::dpois(1, del_lambda, false) << "\t";
    // if(l == 4 && ref_index(i, j) != -1) {
    //   if (viterbi)
    //     likelihood += log(w_ic(i, K)) + R::dpois(1, ins_lambda, true);
    //   else
    //   likelihood += w_ic(i, K) * R::dpois(1, ins_lambda, false);
    // }
    
    // if(K == 0 && l == 0) {
    //   Rcout << "predictor : " << predictor.t() << "\n";
    //   Rcout << "nuc" << obs[index[i] + ref_index(i, j)] << "\t";
    // }
    // Rcout << likelihood << "\t";
  }
  return likelihood;
}/* m_haplotype_llk */

/*
 * Consider deletion as a state: M step (M1) to update the deletions conditional on the haplotypes 
 * (allowing, I guess, to drop some parts of the haplotype?),
 * go back to E step, then take another M step (M2) to update the haplotypes conditional on the deletions.
 */

// [[Rcpp::export]]
IntegerMatrix m_hap (List par, List dat_info, IntegerMatrix haplotype, unsigned int PD_LENGTH,
            IntegerMatrix hidden_state, IntegerVector SNP) {
  NumericMatrix w_ic = par["wic"];
  IntegerVector excluded_read = par["excluded_read"];
  IntegerVector ref = dat_info["ref_pos"];
  NumericVector beta = par["beta"];
  int hap_length = dat_info["ref_length_max"];
  double del_rate = par["del_rate"]; 
  double del_lambda = del_rate* hap_length;
  double ins_rate = par["ins_rate"]; 
  double ins_lambda = ins_rate* hap_length;
  IntegerMatrix m_haplotype(NUM_CLASS, hap_length);
  IntegerVector hap_deletion_len(NUM_CLASS);
  unsigned int j, k, l, max_id;
  double max;
  double llk;
  
  // Rcpp::Nullable<IntegerVector> SNP = R_NilValue
  // IntegerVector snp(hap_length);
  // if (SNP.isNull())
  //   snp = 0;
  // else
  //   snp = clone(SNP);
  
  /* pick out the updated sites for each haplotype */
  for (k = 0; k < NUM_CLASS; ++k) {
     //Rprintf("k: %d\n", k);
    for (j = 0; j < hap_length; ++j) {
       //Rprintf("j: %d\n", j);
      if(SNP[j] == 0) {
        m_haplotype(k, j) = haplotype(k, j);
        continue;
      }
      // Update other sites condition on indels
        if (hidden_state(k, j) == 1) {
          m_haplotype(k, j) = 4;
          continue;
        }
        max = 0;
        max_id = 0;
        for (l = 0; l < NONINDEL_STATE; ++l) {
          llk = m_haplotype_llk(j, l, k, PD_LENGTH, beta, dat_info, w_ic, excluded_read, 
                                ins_lambda, del_lambda, 0);
          // Rprintf("\n neu likelihood %d %f\n", l, log(llk));
          if (llk > max) {
            max = llk;
            max_id = l;
          }
        }
        if(max_id != haplotype(k, j))
          Rprintf(" Cluster %d site %d update from %d to %d\n", k, j, haplotype(k, j), max_id);
        m_haplotype(k, j) = max_id;
    }
  }
    return(m_haplotype);
}

// [[Rcpp::export]]
/*
 * Viterbi Algorithm:
 * hidden: MN
 */
List viterbi_MN(List par, List dat_info) {
  NumericMatrix w_ic = par["wic"];
  IntegerMatrix ref_index = dat_info["ref_idx"];
  IntegerVector excluded_read = par["excluded_read"];
  List deletion = dat_info["deletion"];
  IntegerVector del_strat_id = deletion["del_strat_id"];
  IntegerVector del_length_all = deletion["del_length_all"];
  IntegerVector del_ref_pos = deletion["del_ref_pos"];
  IntegerVector del_flag = deletion["del_flag"];
  IntegerVector ref_pos = dat_info["ref_pos"];
  IntegerVector length = dat_info["length"];
  IntegerVector index = dat_info["start_id"];
  int del_total = deletion["del_total"];
  int nondel_total = dat_info["total"];
  double del_rate = par["del_rate"]; 
  double ins_rate = par["ins_rate"]; 
  unsigned int hap_length = dat_info["ref_length_max"];
  unsigned int n_observation = dat_info["n_observation"];
  
  arma::Cube<double> path_jnk(hap_length, HIDDEN_STATE, NUM_CLASS);
  arma::Cube<double> transition_prob(HIDDEN_STATE, HIDDEN_STATE, NUM_CLASS);
  IntegerVector m_to_n(n_observation);
  IntegerVector n_to_n(n_observation);
  IntegerMatrix hidden_state(NUM_CLASS, hap_length);
  IntegerVector hap_deletion_len(NUM_CLASS);
  
  double del_lambda = del_rate * hap_length;
  double ins_lambda = ins_rate * hap_length;
  double max;
  int total = del_total + nondel_total;
  unsigned int i, j, j1, k, l, l1, max_id;
  /*
   * log transition prob, 0: state M [NOTICE: some reads may be aligned to the middle of a genome, when update haplotypes, 
   * remember to count that out from deletion]
   * [TODO: no need to compute mton ntom everytime]
   */
  for (k = 0; k < NUM_CLASS; ++k) {
    //Rprintf("%d\n", k);
    for (i = 0; i < n_observation; ++i) {
      //Rprintf("%d\t", i);
      if (excluded_read[i] || del_flag[i] == 0)
        continue;
      m_to_n(i) = 1;
      n_to_n(i) = 0;
      for (j = 0, j1 = 1; j1 < del_length_all[i]; ++j, ++j1) {
        if ((del_ref_pos(j1 + del_strat_id(i)) - del_ref_pos(j + del_strat_id(i))) > 1)
          m_to_n(i)++;
        if ((del_ref_pos(j1 + del_strat_id(i)) - del_ref_pos(j + del_strat_id(i))) == 1)
          n_to_n(i)++;
      }
      //Rcout << m_to_n(i) << " " << n_to_n(i) << "\n";
      transition_prob(0, 1, k) += m_to_n(i) * w_ic(i, k);
      transition_prob(1, 1, k) += n_to_n(i) * w_ic(i, k);
    }
    //Rcout << transition_prob(0, 1, k) << " " << transition_prob(1, 1, k) << "\n";
    transition_prob(0, 1, k) = transition_prob(0, 1, k)/(total);
    transition_prob(0, 0, k) = log(1 - transition_prob(0, 1, k));
    transition_prob(0, 1, k) = log(transition_prob(0, 1, k));
    transition_prob(1, 1, k) = transition_prob(1, 1, k)/(total);
    transition_prob(1, 0, k) = log(1 - transition_prob(1, 1, k));
    transition_prob(1, 1, k) = log(transition_prob(1, 1, k));
  }
   // Rprintf("transition prob:\n");
   // Rcout << transition_prob << "\n";
  //log emission prob for the entire datasets(treat them independently)
  NumericMatrix emmis_prob(hap_length, HIDDEN_STATE);
  NumericMatrix sin_emmis(HIDDEN_STATE, HIDDEN_STATE);
  sin_emmis(0, 1) = R::dpois(1, del_lambda, true); //m to deletion, 0 -> 1
  sin_emmis(0, 0) = log(1 - R::dpois(1, del_lambda, false));
  sin_emmis(1, 1) = log(1 - R::dpois(1, ins_lambda, false)); //n to deletion
  sin_emmis(1, 0) = R::dpois(1, ins_lambda, true);

  for (j = 0; j < hap_length; ++j)
    for(i = 0; i < n_observation; ++i) {
      // if indel position
        if (ref_index(i, j) == -1) {
          if((j >= ref_pos[index[i]]) && (j <= ref_pos[index[i] + length[i] - 1])) {
          emmis_prob(j, 0) += sin_emmis(0, 1);
          emmis_prob(j, 1) += sin_emmis(1, 1);
          }
        }
        else {
          emmis_prob(j, 0) += sin_emmis(0, 0);
          emmis_prob(j, 1) += sin_emmis(1, 0);
        }
    }
     // Rprintf("emmission prob:\n");
     // Rcout << emmis_prob << "\n";
  /*
   * compute log probability table to trace back the hidden state: l = 0, state M
   */
  double max_prob = 0;
  for (k = 0; k < NUM_CLASS; ++k)
    for (j = 0; j < hap_length; ++j)
      for (l = 0; l < HIDDEN_STATE; ++l) {
        if (j == 0) {
          path_jnk(j, l, k) = emmis_prob(j, l);
        } else {
          max = -INFINITY;
          for (l1 = 0; l1 < HIDDEN_STATE; ++l1) {
            max_prob = path_jnk(j - 1, l1, k) + transition_prob(l1, l, k);
            if (max_prob > max)
              max = max_prob;
          }
          path_jnk(j, l, k) = emmis_prob(j, l) + max;
        }
      }
   // Rprintf("path_jnk prob:\n");
   // Rcout << path_jnk << "\n";  
  /*
   * find the path (backtrace)
   */
  for (k = 0; k < NUM_CLASS; ++k)
    for (j = hap_length; j --> 0;) {
      max = -INFINITY;
      max_id = 0;
      for (l = 0; l < HIDDEN_STATE; ++l)
        if (path_jnk(j, l, k) > max) {
          max_id = l;
          max = path_jnk(j, l, k);
        }
      hidden_state(k, j) = max_id;
      if(hidden_state(k, j) == 1) {
        // Rprintf("%d %d deletion state\n", k, j);
        hap_deletion_len[k]++;
      }
    }
  unsigned int sum = 0;
  for(k = 0; k < NUM_CLASS; ++k)
    sum += hap_deletion_len[k];
  IntegerVector hap_ref_pos(sum);
  IntegerVector strat_id(NUM_CLASS);
  
  sum = 0;
  for (k = 0; k < NUM_CLASS; ++k) {
    for (j = 0; j < hap_length; ++j)
      if (hidden_state(k, j) == 1)
        hap_ref_pos[sum++] = j;
    if (hap_deletion_len[k] != 0) {
      if (k != 0)
        for (l = 0; l < k; ++l)
          strat_id[k] += hap_deletion_len[l];
    } else
      strat_id[k] = -1; //record the deletion starting id in vector hap_ref_pos for each hap
  }
  
  List ls = List::create(
    Named("hidden") = hidden_state,
    Named("deletion_pos") = hap_ref_pos,
    Named("hap_deletion_len") = hap_deletion_len,
    Named("hap_del_start_id") = strat_id);
  return(ls);
}


