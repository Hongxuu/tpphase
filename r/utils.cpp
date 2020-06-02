#include <iostream>
#include <stdio.h>
#include <string>
#include <unordered_map>
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace std;
// hash a vector 
List unique_map(const Rcpp::IntegerVector & v)
{
  // Initialize a map
  std::map<double, int> Elt;
  Elt.clear();
  
  // Count each element
  for (int i = 0; i != v.size(); ++i)
    Elt[ v[i] ] += 1;
  
  // Find out how many unique elements exist... 
  int n_obs = Elt.size();
  // If the top number, n, is greater than the number of observations,
  // then drop it.  
  int n = n_obs;
  
  // Pop the last n elements as they are already sorted. 
  // Make an iterator to access map info
  std::map<double,int>::iterator itb = Elt.end();
  // Advance the end of the iterator up to 5.
  std::advance(itb, -n);
  
  // Recast for R
  NumericVector result_vals(n);
  NumericVector result_keys(n);
  
  unsigned int count = 0;
  // Start at the nth element and move to the last element in the map.
  for (std::map<double,int>::iterator it = itb; it != Elt.end(); ++it) {
    // Move them into split vectors
    result_keys(count) = it->first;
    result_vals(count) = it->second;
    count++;
  }
  return List::create(Named("lengths") = result_vals,
                      Named("values") = result_keys);
}

// hash a matrix by row
List hash_mat(IntegerMatrix x) {
  int n = x.nrow() ;
  int nc = x.ncol() ;
  std::vector<string> hashes(n) ;
  // arma::Mat<int> X = as<arma::Mat<int>>(x);
  for (int i = 0; i < n; i++) {
    string s = "";  
    for(int j = 0; j < nc; j++)  
      s += to_string(x(i,j));  
    hashes[i] = s;
  }
  
  std::unordered_map<string, vector<int>> map;
  for (int i = 0; i < n; i++)
    map[hashes[i]].push_back(i);
  
  int nres = map.size();
  IntegerVector idx(nres);
  List all_id(nres);
  
  int i = 0; 
  for (auto itr = map.begin(); itr != map.end(); ++itr) { 
    idx[i] = itr->second[0];
    all_id[i++] = wrap(itr->second);
  } 
  return List::create( _["all_id"] = all_id, _["idx"] = idx );
}

IntegerVector find_max(List ls, int n_obs) {
  unsigned int i, j, k = 0;
  IntegerVector max_id(n_obs);
  for(i = 0; i < ls.size(); ++i) {
    NumericMatrix mat = ls[i];
    for(j = 0; j < mat.nrow(); ++j) {
      NumericVector vec = mat(j, _);
      int id = arma::index_max(as<arma::vec>(vec));
      max_id[k++] = id;
    }
  }
  return(max_id);
}
// maybe a slower version
// List hash_mat(IntegerMatrix x) {
//   int n = x.nrow() ;
//   int nc = x.ncol() ;
//   std::vector<string> hashes(n) ;
//   // arma::Mat<int> X = as<arma::Mat<int>>(x);
//   for (int i = 0; i < n; i++) {
//     string s = "";  
//     for(int j = 0; j < nc; j++)  
//       s += to_string(x(i,j));  
//     hashes[i] = s;
//   }
//   
//   using Pair = std::pair<int, vector<int>>;
//   std::unordered_map<string, Pair> map_counts;
//   for (int i = 0; i < n; i++) {
//     Pair& p = map_counts[hashes[i]];
//     if(p.first == 0) {
//       p.first = i;
//     }
//     p.second.push_back(i);
//   }
//   
//   int nres = map_counts.size();
//   IntegerVector idx(nres);
//   List all_id(nres);
//   auto it = map_counts.begin();
//   for(int i = 0; i < nres; i++, ++it) {
//     idx[i] = it->second.first;
//     all_id[i] = wrap(it->second.second);
//   }
//   
//   return List::create( _["all_id"] = all_id, _["idx"] = idx );
// }

// subset matrix by row index
IntegerMatrix ss(IntegerMatrix X_, IntegerVector ind_) {
  
  int n = X_.nrow(), k = X_.ncol();
  arma::Mat<int> X(X_.begin(), n, k, false);
  arma::uvec ind = as<arma::uvec>(ind_);
  arma::Mat<int> submat = X.rows(ind);
  
  return wrap(submat);
}

template <typename T>
inline bool approx_equal_cpp(const T& lhs, const T& rhs, double tol = 0.00000001) {
  return arma::approx_equal(lhs, rhs, "absdiff", tol);
}

arma::mat unique_rows(const arma::mat& m) {
  arma::uvec ulmt = arma::zeros<arma::uvec>(m.n_rows);
  for (arma::uword i = 0; i < m.n_rows; i++)
    for (arma::uword j = i + 1; j < m.n_rows; j++)
      if (approx_equal_cpp(m.row(i), m.row(j))) { ulmt(j) = 1; break; }
      
      return m.rows(find(ulmt == 0));
}
