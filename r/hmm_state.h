#ifndef HMM_STATE_H
#define HMM_STATE_H

#include <Rcpp.h>
using namespace Rcpp;
List prepare_ini_hmm (unsigned int t_max, IntegerVector num_states, List trans_indicator, List further_limit);
List sbs_state(unsigned int num, unsigned int ref_j, IntegerVector hap_site, IntegerVector sum_site, 
               CharacterVector uni_alignment);
List limit_comb(IntegerMatrix combination, List hidden_states, IntegerVector location,
                         IntegerMatrix linkage_info, unsigned int num, unsigned int start_idx, unsigned int num_states);
#endif
