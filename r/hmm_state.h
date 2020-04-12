#ifndef DATA_FORMAT_H
#define DATA_FORMAT_H

#include <Rcpp.h>
using namespace Rcpp;
List sbs_state(unsigned int num, unsigned int ref_j, IntegerVector hap_site, IntegerVector sum_site, 
               CharacterVector uni_alignment);
#endif
