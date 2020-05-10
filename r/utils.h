#ifndef UTILS_H
#define UTILS_H

#include <Rcpp.h>
using namespace Rcpp;

List unique_map(const Rcpp::IntegerVector & v);
List hash_mat(IntegerMatrix x);
IntegerMatrix ss(IntegerMatrix X_, IntegerVector ind_);
#endif