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
#define PD_LENGTH 10

double m_haplotype_llk (unsigned int j, unsigned int l, unsigned int K,
                        List par, List dat_info, NumericMatrix w_ic);

double m_haplotype_llk (unsigned int j, unsigned int l, unsigned int K, 
                        NumericVector beta, List dat_info, NumericMatrix w_ic,
                        IntegerVector excluded_read)
{
  unsigned int i, k, c;
  double qua_in, read_pos_in, ref_pos_in;
  //unsigned int index = 0;
  double tail, sum;
  double xb = 0;
  double likelihood = 0.;
  
  int n_observation = dat_info["n_observation"];
  IntegerVector qua = dat_info["qua"];
  IntegerVector obs = dat_info["nuc"];
  IntegerVector ref_pos = dat_info["ref_pos"];
  IntegerVector read_pos = dat_info["read_pos"];
  IntegerVector length = dat_info["length"];
  IntegerVector index = dat_info["start_id"];
  IntegerMatrix ref_index = dat_info["ref_idx"];
  
  NumericVector hap_nuc(MLOGIT_CLASS - 1);
  NumericVector hnuc_qua(MLOGIT_CLASS - 1);
  NumericVector pred_beta(MLOGIT_CLASS - 1);
  
  for (i = 0; i < n_observation; ++i) {
    // if (i != 0)
    //   index += length[i - 1];
    //Rprintf("idx position %d\n", index[i] + ref_index(i, j));
    if (ref_index(i, j) == -1 || excluded_read[i] == 1)
      continue;
    
    qua_in = qua[index[i] + ref_index(i, j)];
    read_pos_in = read_pos[index[i] + ref_index(i, j)];
    ref_pos_in = ref_pos[index[i] + ref_index(i, j)];
    
    for (k = 0; k < MLOGIT_CLASS - 1; ++k) {
      hap_nuc[k] = 0;
      hnuc_qua[k] = 0;
      pred_beta[k] = 0;
    }
    
    if (l == 1) {
      hap_nuc[0] = 1;
      hnuc_qua[0] = qua_in;
    } else if (l == 3) {
      hap_nuc[1] = 1;
      hnuc_qua[1] = qua_in;
    } else if (l == 2) {
      hap_nuc[2] = 1;
      hnuc_qua[2] = qua_in;
    }
    
    arma::vec predictor = {1, read_pos_in, ref_pos_in, qua_in, hap_nuc[0], hap_nuc[1], 
                           hap_nuc[2], hnuc_qua[0], hnuc_qua[1], hnuc_qua[2]};
    
    // cblas_dgemv(CblasRowMajor, CblasNoTrans, MLOGIT_CLASS - 1, PD_LENGTH, 1, beta, PD_LENGTH, 
    //             predictor, 1, 0, pred_beta, 1);
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
    // if(K == 0 && l == 0) {
    //   Rcout << "predictor : " << predictor.t() << "\n";
    //   Rcout << "nuc" << obs[index[i] + ref_index(i, j)] << "\t";
    // }
    //Rcout << "xb : " << xb << "\n";
    /* Still exclude nan */
    //if (!isnan(w_ic(i, K))) {
    likelihood += w_ic(i, K) * exp(xb + tail);
    //}
  }
  return likelihood;
}/* m_haplotype_llk */
    
    
// [[Rcpp::export]]
IntegerMatrix m_hap (List par, List dat_info, IntegerMatrix haplotype, Rcpp::Nullable<IntegerVector> SNP = R_NilValue) {
  NumericMatrix w_ic = par["wic"];
  IntegerVector excluded_read = par["excluded_read"];
  IntegerVector ref = dat_info["ref_pos"];
  IntegerVector non_covered = dat_info["non_covered_site"];
  NumericVector beta = par["beta"];
  int hap_length = dat_info["ref_length_max"];
  
  IntegerMatrix m_haplotype(NUM_CLASS, hap_length);
  
  unsigned int j, k, l, max_id;
  double max;
  double llk;
  
  IntegerVector snp(hap_length);
  if (SNP.isNull())
    snp = 0;
  else
    snp = clone(SNP);
  
  /* pick out the updated sites for each haplotype */
  for (k = 0; k < NUM_CLASS; ++k) {
    
    // Rprintf("k: %d\n", k);
    for (j = 0; j < hap_length; ++j) {
      if(non_covered[j] == 1 || snp[j] == 0) {
        m_haplotype(k, j) = haplotype(k, j);
        //Rprintf("%d is not covered by any of the reads!\n", j);
        continue;
      }
      
      // Rprintf("j: %d\n", j);
      // For future reference, only update variation sites
      // if (excluded_site(k, j) == 1)
      //   continue;
      max = 0;
      max_id = 0;
      //llk = 0;
      for (l = 0; l < MLOGIT_CLASS; ++l) {
        llk = m_haplotype_llk(j, l, k, beta, dat_info, w_ic, excluded_read);
        // Rprintf("neu likelihood %d %f\t\t", l, llk);
        if (llk > max) {
          max = llk;
          max_id = l;
        }
      }
      if(max_id != haplotype(k, j)) {
        Rprintf(" Cluster %d site %d update from %d to %d\n", k, j, haplotype(k, j), max_id);
      } 
      m_haplotype(k, j) = max_id;
    }
  }
  return(m_haplotype);
}
