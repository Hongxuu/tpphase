#ifndef BAUMWELCH_H
#define BAUMWELCH_H

#include <Rcpp.h>
using namespace Rcpp;
List find_combination(IntegerVector undecided_pos, IntegerVector pos_possibility, 
                      unsigned int p_tmax, unsigned int time_pos, int hap_min_pos);
IntegerMatrix fill_all_hap(List hidden_states, unsigned int hap_length, IntegerVector n_row);
List get_overlap(IntegerVector p_tmax, IntegerVector time_pos, IntegerVector num_states,
                 IntegerVector undecided_pos, unsigned int t_max, int hap_min_pos);
#endif
