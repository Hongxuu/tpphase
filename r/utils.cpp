#include <iostream>
#include <stdio.h>
#include <string>
#include <unordered_map>
#include <RcppArmadillo.h>

#include "utils.h"
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

// hash a matrix by row, output the index of unique and the same rows
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

List hash_intvec(IntegerVector x) {
  int n = x.size();
  std::map<int, vector<int>> map;
  for (int i = 0; i < n; i++)
    map[x[i]].push_back(i);
  int nres = map.size();
  IntegerVector value(nres);
  List all_id(nres);
  
  int i = 0; 
  for (auto itr = map.begin(); itr != map.end(); ++itr) { 
    value[i] = itr->first;
    all_id[i++] = wrap(itr->second);
  }
  
  return List::create(_["all_id"] = all_id, _["value"] = value);
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

// order hased mat by the second column
// IntegerMatrix subset(idx.size(), 2);
// IntegerMatrix ordered_read(idx.size(), 2);
// IntegerVector order(idx.size());
// int summ = 0;
// 
// for(k = 0; k < nuc1.size(); ++k)
//   for(i = 0; i < idx.size(); ++i) {
//     int pre_read = link(idx[i], j + 1);
//     if(pre_read == nuc1[k])
//       order[i] = summ++;
//   }
// 
// Rcout << order<< "\n";
// for(i = 0; i < idx.size(); ++i) {
//   IntegerVector tmp = link(idx[i], _);
//   subset(i, _) = tmp[Range(j, j + 1)];
// }
// for(i = 0; i < idx.size(); ++i)
//   ordered_read(order[i], _) = subset(i, _);

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

IntegerMatrix comb_element(List len, IntegerVector flag, unsigned int row) {
  vector<vector<int> > vec;
  int col, i, j;
  for (i = 0; i < len.size(); i++) {
    if(flag[i])
      continue;
    IntegerVector row_vec = len[i];
    col = row_vec.size();
    vector<int> v1(col);
    for (j = 0; j < col; j++) {
      v1[j] = row_vec[j];
    }
    vec.push_back(v1);
  }
  // for (int i = 0; i < vec.size(); i++) { 
  //   for (int j = 0; j < vec[i].size(); j++) 
  //     Rcout << vec[i][j] << " "; 
  //   Rcout << "\n"; 
  // }
  vector<vector<int> > res = cart_product(vec);
  IntegerMatrix out(res.size(), row);
  for(i = 0; i < res.size(); ++i)
    for(j = 0; j < row; ++j)
      out(i, j) = res[i][j];
  return(out);
}

IntegerMatrix call_permute(vector<int> a) {
  Permutation per;
  vector<vector<int> > res =  per.permuteUnique(a);
  IntegerMatrix permutation(res.size(), res[0].size());
  for (int i = 0; i < res.size(); i++)
    for (int j = 0; j < res[i].size(); j++) 
      permutation(i, j) = res[i][j];
  
  return(permutation);
}

// sort matrix by every column
IntegerMatrix sort_mat(IntegerMatrix mat, int nrow, int ncol) {
  int i, j;
  IntegerMatrix sorted(nrow, ncol);
  vector<int> vec(nrow);
  for (i = 0; i < nrow; i++) {
    for (j = 0; j < ncol; j++) {
      vec[i] += pow(10, ncol - 1 - j) * mat(i, j);
    }
  }
  vector<pair<int, int> > vp; 
  for (int i = 0; i < nrow; ++i) { 
    vp.push_back(make_pair(vec[i], i)); 
  } 
  
  // Sorting pair vector 
  sort(vp.begin(), vp.end()); 
  
  for(i = 0; i < nrow; ++i) 
    sorted(i, _) = mat(vp[i].second, _);
  return sorted;
}

// matrix to vector
IntegerVector matrix2vec(IntegerMatrix m, const bool byrow) {
  if (byrow){
    m = transpose(m);
  }
  IntegerVector x = IntegerVector(m);
  x.attr("dim") = R_NilValue;
  return(x);
}

void print_intmat(IntegerMatrix m) {
  for(int i = 0; i < m.nrow(); ++i) {
    for(int j = 0; j < m.ncol(); ++j) 
      Rcout << m(i, j) << "\t";
    Rcout << "\n";
  }
}