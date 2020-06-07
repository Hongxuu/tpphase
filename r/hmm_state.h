#ifndef HMM_STATE_H
#define HMM_STATE_H

#include <Rcpp.h>

using namespace Rcpp;

IntegerMatrix mc_linkage(IntegerMatrix sub_link, int num);
IntegerVector best_branch(IntegerMatrix link, List transition, NumericVector initial, 
                          List possi_nuc, int i);
IntegerMatrix remake_linkage(IntegerMatrix sub_link, unsigned int num);
List prepare_ini_hmm (unsigned int t_max, IntegerVector num_states, List trans_indicator, List further_limit);
List sbs_state(unsigned int num, unsigned int ref_j, IntegerVector hap_site, IntegerVector sum_site, 
               CharacterVector uni_alignment);
List limit_comb_t0(IntegerMatrix combination, List hidden_states, IntegerVector location,
                   IntegerMatrix linkage_info, unsigned int num, unsigned int start_idx, 
                   unsigned int num_states, unsigned int use_MC);
List find_combination(IntegerVector undecided_pos, IntegerVector pos_possibility, 
                      unsigned int p_tmax, unsigned int time_pos, int hap_min_pos);
IntegerMatrix make_hap(List hidden_states, IntegerMatrix haplotype, IntegerVector location, unsigned int p_tmax,
                       IntegerVector combination, unsigned int time_pos, unsigned int num, int hap_min_pos);
IntegerMatrix unique_overlap(IntegerVector overlapped, IntegerVector exclude_last, IntegerMatrix overlap_comb, IntegerVector overlap_loci, 
                             unsigned int overlap_new_states, unsigned int overlap_num_states);
IntegerMatrix new_combination(List hmm_info, IntegerVector location, IntegerVector overlapped, IntegerVector exclude_last, IntegerMatrix overlap_comb, 
                              IntegerVector overlap_loci, IntegerMatrix linkage_info, unsigned int overlap_new_states, unsigned int overlap_num_states,
                              unsigned int use_MC);
#endif
