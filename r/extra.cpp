#include <RcppArmadillo.h>
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <ctype.h>
#include <limits.h>
#include <unistd.h>
#include <math.h>
#include <iostream>
#include <vector>
#include <algorithm>

#include "data_format.h"
using namespace Rcpp;
using namespace std;
// [[Rcpp::depends(RcppArmadillo)]]
#define NUM_CLASS 4
#define MLOGIT_CLASS 4

// [[Rcpp::export]]
List make_universal_old(List A_genome, List B_genome) {
  IntegerVector A_aligned = A_genome["reads"];
  IntegerVector B_aligned = B_genome["reads"];
  IntegerVector B_lengths = B_genome["dim"];

  unsigned int i, j;
  IntegerVector start_id(B_lengths.length());
  CharacterVector universal(B_aligned.length());

  for(i = 1; i < B_lengths.length(); ++i)
    start_id[i] = start_id[i - 1] + B_lengths[i - 1];

  for(j = 0; j < B_aligned.length(); ++j) {
    if(B_aligned[j] == 15 & A_aligned[j] != 15) {
      universal[j] = "J";
    } else if(B_aligned[j] != 15 & A_aligned[j] == 15) {
      universal[j] = "I";
    } else
      universal[j] = "M";
  }

  List ls = List::create(
    Named("universal_alignment") = universal,
    Named("start_id") = start_id);
  return ls;
}

