#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
CharacterVector convert_char(IntegerVector nuc, int length) {
  int i = 0;
  CharacterVector nuc_char(length);
    
  for(i = 0; i < length; ++i) {
    if(nuc[i] == 0) 
      nuc_char[i] = 'A';
    if(nuc[i] == 1)
      nuc_char[i] = 'C';
    if(nuc[i] == 2)
      nuc_char[i] = 'T';
    else if(nuc[i] == 3)
      nuc_char[i] = 'G';
  }
  
  return(nuc_char);
}

// [[Rcpp::export]]
arma::colvec solveC (NumericMatrix ar, NumericVector br) {
  int n = ar.nrow(), k = ar.ncol();
  
  arma::mat a(ar.begin(), n, k, false);
  arma::colvec b(br.begin(), br.size(), false);
  
  arma::colvec x = arma::solve(a, b);
  //NumericVector re;
  //re = as<NumericVector>(wrap(x)); 
  return(x);
}

// [[Rcpp::export]]
IntegerVector match_c (NumericVector ar, NumericVector br){
  unsigned int n = ar.length(), k = br.length();
  IntegerVector index(n);
  unsigned int m = 0;
  
  for(unsigned int i = 0; i < n; ++i) 
    for (unsigned int j = 0; j < k; ++j)
      if(ar[i] == br[j])
        index[m++] = i;
  
  return index(Range(0, m - 1));
}

// [[Rcpp::export]]

// Find the mutural deletions in reads and haps


/*** R
### Slower
rcpp_inc <- '
using namespace Rcpp;
using namespace arma;
'
src_cp <- '
mat m1 = as<mat>(m1in);
mat m2 = as<mat>(m2in);
mat cp = trans(m1) * m2;
return(wrap(cp));
'
library(inline)
crossprodC <- cxxfunction(signature(m1in="numeric", m2in="numeric"), src_cp, plugin='RcppArmadillo', rcpp_inc)
*/
