#ifndef UTILS_H
#define UTILS_H

#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace std;

class Permutation {
public:
  vector<vector<int> > permuteUnique(vector<int>& nums) {
    vector<vector<int> > res;
    sort(nums.begin(), nums.end());
    res.push_back(nums);
    while (next_permutation(nums.begin(), nums.end())) {
      res.push_back(nums);
    }
    return res;
  }
};
// struct node  
// { 
//   int data; 
//   double prob;
//   std::vector<node *>child;
// }; 
// 
// /* newNode() allocates a new node with the given data and NULL left and  
//  right pointers. */
// struct node* newNode(int data, double prob) 
// { 
//   // Allocate memory for new node  
//   struct node *temp = new node;
//   // Assign data to this node 
//   temp->data = data; 
//   temp->prob = prob;
//   return(temp); 
// } 
IntegerMatrix call_cart_product(IntegerVector len);
IntegerMatrix comb_element(List len, IntegerVector flag, unsigned int row);
List unique_map(const Rcpp::IntegerVector & v);
List hash_mat(IntegerMatrix x);
IntegerMatrix ss(IntegerMatrix X_, IntegerVector ind_);
arma::mat unique_rows(const arma::mat& m);
IntegerVector find_max(List ls, int n_obs);
IntegerMatrix call_permute(vector<int> a);
#endif