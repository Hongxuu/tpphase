#ifndef UTILS_H
#define UTILS_H

#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

List unique_map(const Rcpp::IntegerVector & v);
List hash_mat(IntegerMatrix x);
IntegerMatrix ss(IntegerMatrix X_, IntegerVector ind_);
arma::mat unique_rows(const arma::mat& m);
IntegerVector find_max(List ls, int n_obs);
#endif