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
#define HAP_STATE 5

// A C G T N: 0 1 2 3 4

double m_haplotype_llk (unsigned int j, unsigned int l, unsigned int K, unsigned int PD_LENGTH,
                        NumericVector beta, List dat_info, NumericMatrix w_ic, IntegerVector excluded_read, double lambda);

double m_haplotype_llk (unsigned int j, unsigned int l, unsigned int K, unsigned int PD_LENGTH,
                        NumericVector beta, List dat_info, NumericMatrix w_ic, IntegerVector excluded_read, double lambda)
{
  unsigned int i, k, c;
  double qua_in, read_pos_in, ref_pos_in;
  //unsigned int index = 0;
  double tail, sum;
  double xb = 0;
  double likelihood;
  
  int n_observation = dat_info["n_observation"];
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
  arma::vec predictor(PD_LENGTH);
  
  likelihood = 0;
  for (i = 0; i < n_observation; ++i) {
    //Rprintf("idx position %d\n", index[i] + ref_index(i, j));
    if (excluded_read[i] == 1)
      continue;
    
    if (l != 4 && ref_index(i, j) != -1) {
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
      
      //if(!N_in) {
        predictor = {1, read_pos_in, ref_pos_in, qua_in, hap_nuc[0], hap_nuc[1], 
                     hap_nuc[3], hnuc_qua[0], hnuc_qua[1], hnuc_qua[3]};
      // } else {
      //   predictor = {1, read_pos_in, ref_pos_in, qua_in, hap_nuc[0], hap_nuc[1], 
      //                hap_nuc[2], hap_nuc[3], hnuc_qua[0], hnuc_qua[1], hnuc_qua[2], hnuc_qua[3]};
      // }
      
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
      likelihood += w_ic(i, K) * exp(xb + tail);
      //Rcout << "normal add " << w_ic(i, K) * exp(xb + tail) << "\t";
    }
    if(l == 4 && ref_index(i, j) == -1) {
      likelihood += w_ic(i, K) * R::dpois(1, lambda, false);
      //Rcout << "del add " << w_ic(i, K) * R::dpois(1, lambda, false) << "\t";
    }
    
    // if(K == 0 && l == 0) {
    //   Rcout << "predictor : " << predictor.t() << "\n";
    //   Rcout << "nuc" << obs[index[i] + ref_index(i, j)] << "\t";
    // }
    //Rcout << likelihood << "\t";
  }
  return likelihood;
}/* m_haplotype_llk */

// [[Rcpp::export]]
// Consider deletion as a state
List m_hap (List par, List dat_info, IntegerMatrix haplotype, unsigned int PD_LENGTH, unsigned int N_in,
                     IntegerVector SNP) {
  NumericMatrix w_ic = par["wic"];
  IntegerVector excluded_read = par["excluded_read"];
  IntegerVector ref = dat_info["ref_pos"];
  IntegerVector non_covered = dat_info["non_covered_site"];
  NumericVector beta = par["beta"];
  int hap_length = dat_info["ref_length_max"];
  double rate = par["rate"]; 
  double lambda = rate* hap_length;
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
      if(non_covered[j] == 1 || SNP[j] == 0) {
        m_haplotype(k, j) = haplotype(k, j);
        //Rprintf("%d %d!\n", non_covered[j], SNP[j]);
        continue;
      }
      max = 0;
      max_id = 0;
      for (l = 0; l < HAP_STATE; ++l) {
        llk = m_haplotype_llk(j, l, k, PD_LENGTH, beta, dat_info, w_ic, excluded_read, lambda);
        //Rprintf("\n neu likelihood %d %f\n", l, log(llk));
        if (llk > max) {
          max = llk;
          max_id = l;
        }
      }
      if(max_id != haplotype(k, j)) {
        Rprintf(" Cluster %d site %d update from %d to %d\n", k, j, haplotype(k, j), max_id);
      }
      if (max_id == 4) {
        hap_deletion_len[k]++;
      }
      m_haplotype(k, j) = max_id;
    }
  }
  unsigned int sum = 0;
  for(k = 0; k < NUM_CLASS; ++k)
    sum += hap_deletion_len[k];
  IntegerVector hap_ref_pos(sum);
  IntegerVector strat_id(NUM_CLASS);
  
  sum = 0;
  for (k = 0; k < NUM_CLASS; ++k) {
    for (j = 0; j < hap_length; ++j) {
      if (m_haplotype(k, j) == 4)
        hap_ref_pos[sum++] = j;
    }
    if (hap_deletion_len[k] != 0) {
      if(k != 0)
        for (l = 0; l < k; ++l)
          strat_id[k] += hap_deletion_len[l];
    }
    else
      strat_id[k] = -1; //record the deletion starting id in vector hap_ref_pos for each hap
  }
  List ls = List::create(
    Named("hap") = m_haplotype,
    Named("deletion_pos") = hap_ref_pos,
    Named("hap_deletion_len") = hap_deletion_len,
    Named("hap_del_start_id") = strat_id);
  return(ls);
}
