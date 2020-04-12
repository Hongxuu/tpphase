#include <RcppArmadillo.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

using namespace Rcpp;
using namespace std;
#define NUM_CLASS 4

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

IntegerMatrix Twopossible(IntegerVector a) {
  IntegerVector temp = {a[0], a[1], a[0], a[1], a[1], a[0], a[1], a[0]};
  temp.attr("dim") = Dimension(2, 4);
  return(as<IntegerMatrix>(temp));
}

IntegerMatrix Fourpossible(IntegerVector small, int big) {
  IntegerVector temp = {small[0], small[1], big, big, small[1], small[0], big, big,
                        big, big, small[0], small[1], big, big, small[1], small[0]};
  temp.attr("dim") = Dimension(4, 4);
  return(as<IntegerMatrix>(temp));
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

IntegerMatrix call_permute_N(vector<int> a, unsigned int genome_A) {
  Permutation per;
  vector<vector<int> > res =  per.permuteUnique(a);
  IntegerMatrix permutation(res.size(), res[0].size());
  for (int i = 0; i < res.size(); i++)
    for (int j = 0; j < res[i].size(); j++) 
      permutation(i, j) = res[i][j];
  IntegerVector del(res.size() * 2, -1);
  del.attr("dim") = Dimension(res.size(), 2);
  arma::mat m1 = as<arma::mat>(permutation);
  arma::mat m2 = as<arma::mat>(del);
  arma::mat out;
  if(genome_A)
    out = join_rows(m2, m1);
  else
    out = join_rows(m1, m2);
  
  return(wrap(out));
}

List aux_noN_S2(IntegerVector sum_site, IntegerVector hap_site) {
  
  int n_row;
  double sum;
  IntegerMatrix temp;
  
  sum = sum_site[1] + sum_site[0];
  if (sum_site[0]/sum <= 0.38) { //1st one appears 3 times
    temp = call_permute({hap_site[0], hap_site[1], hap_site[1], hap_site[1]});
    n_row = 4;
  } else if (sum_site[0]/sum >= 0.62) {//2nd one appears 3 times
    temp = call_permute({hap_site[0], hap_site[0], hap_site[0], hap_site[1]});
    n_row = 4;
  } else { // the possibility of occuring ACAC is small? so elimiate this possibility???
    temp = Twopossible(hap_site);
    n_row = 2;
  }
  
  List ls = List::create(
    Named("n_row") = n_row,
    Named("temp") = temp);
  return ls; 
}

List aux_noN_S3(IntegerVector sum_site, IntegerVector hap_site) {
  int n_row;
  double sum;
  IntegerMatrix temp;
  List from_2(2);
  
  sum = sum_site[2] + sum_site[1] + sum_site[0];
  double p0 = sum_site[0]/sum;
  double p1 = sum_site[1]/sum;
  double p2 = sum_site[2]/sum;
  IntegerVector small(2);
  IntegerVector small_count(2);
  // different ceriterion to reduce the possibilities, basically rule out things like ATAG
  if (((sum_site[2] + sum_site[1])/sum_site[0] <= 1.2 && (sum_site[2] + sum_site[1])/sum_site[0] >= 0.8) ||
      (sum_site[0]/sum >= 0.45)) {
    small = {hap_site[1], hap_site[2]};
    temp = Fourpossible(small, hap_site[0]);
    n_row = 4;
  } else if (((sum_site[1] + sum_site[0])/sum_site[2] <= 1.2 && (sum_site[0] + sum_site[1])/sum_site[2] >= 0.8) ||
    (sum_site[2]/sum >= 0.45)) {
    small = {hap_site[1], hap_site[0]};
    temp = Fourpossible(small, hap_site[2]);
    n_row = 4;
  } else if (((sum_site[2] + sum_site[0])/sum_site[1] <= 1.2 && (sum_site[0] + sum_site[2])/sum_site[1] >= 0.8) ||
    (sum_site[1]/sum >= 0.45)) {
    small = {hap_site[2], hap_site[0]};
    temp = Fourpossible(small, hap_site[1]);
    n_row = 4;
  } else if (p0 <= 0.25) {
    small = {hap_site[1], hap_site[2]};
    small_count = {sum_site[1], sum_site[2]};
    List from_2 = aux_noN_S2(small_count, small);
    IntegerMatrix ttemp = from_2["temp"];
    temp = ttemp;
    n_row = from_2["n_row"];
  } else if (p1 <= 0.25) {
    small = {hap_site[0], hap_site[2]};
    small_count = {sum_site[0], sum_site[2]};
    List from_2 = aux_noN_S2(small_count, small);
    IntegerMatrix ttemp = from_2["temp"];
    temp = ttemp;
    n_row = from_2["n_row"];
  } else if (p2 <= 0.25) {
    small = {hap_site[0], hap_site[1]};
    small_count = {sum_site[0], sum_site[1]};
    List from_2 = aux_noN_S2(small_count, small);
    IntegerMatrix ttemp = from_2["temp"];
    temp = ttemp;
    n_row = from_2["n_row"];
  } else { //hopefully this won't happen, when three roughly equal happens
    small = {hap_site[1], hap_site[2]};
    temp = Fourpossible(small, hap_site[0]);
    arma::mat m1 = as<arma::mat>(temp);
    small = {hap_site[1], hap_site[0]};
    temp = Fourpossible(small, hap_site[2]);
    arma::mat m2 = as<arma::mat>(temp);
    small = {hap_site[2], hap_site[0]};
    temp = Fourpossible(small, hap_site[1]);
    arma::mat m3 = as<arma::mat>(temp);
    arma::mat out = join_cols(m1, m2, m3);
    temp = wrap(out);
    n_row = 12;
  }
  List ls = List::create(
    Named("n_row") = n_row,
    Named("temp") = temp);
  return ls; 
}
/*
 * determin the number of hidden states site by site
 */
List sbs_state(unsigned int num, unsigned int ref_j, IntegerVector hap_site, IntegerVector sum_site, 
               CharacterVector uni_alignment) {
  unsigned int l, m;
  double sum;
  unsigned int n_row;
  List haplotype(1);
  //record possible hidden states[only suitable for ployploids]
  if (num == 1) {
    IntegerVector temp = {hap_site[0], hap_site[0], hap_site[0], hap_site[0]};
    haplotype[0] = temp;
    n_row = 1;
  }
  else if (num == 2) {
    if(hap_site[0] == -1) {
      IntegerVector temp(NUM_CLASS);
      if(uni_alignment[ref_j] == "I")
        temp = {-1, -1, hap_site[1], hap_site[1]};
      else if (uni_alignment[ref_j] == "J")
        temp = {hap_site[1], hap_site[1], -1, -1};
      else
        temp = {hap_site[1], hap_site[1], hap_site[1], hap_site[1]};
      haplotype[0] = temp;
      n_row = 1;
    } else { // if N not appears
      List out = aux_noN_S2(sum_site, hap_site);
      haplotype[0] = out["temp"];
      n_row = out["n_row"];
      // undecided_pos(more_than1) = j; // relative to the haplotype currently try to infer (alignment start from 0)
      // pos_possibility(more_than1++) = n_row;
    }
  }
  else if (num == 3) {
    sum = sum_site[2] + sum_site[1] + sum_site[0];
    if (hap_site[0] == -1) {
      // if N appears, 4 possibilities
      if(sum_site[0]/sum >= 0.45) {
        IntegerVector temp(2 * NUM_CLASS);
        if (uni_alignment[ref_j] == "I")
          temp = {-1, -1,  -1, -1, hap_site[1], hap_site[2], //J in universal alignment
                  hap_site[2], hap_site[1]};
        else if(uni_alignment[ref_j] == "J")
          temp = {hap_site[1], hap_site[2], hap_site[2],hap_site[1],  //I in universal alignment
                  -1,  -1,  -1, -1};
        temp.attr("dim") = Dimension(2, 4); // by column
        haplotype[0] = temp;
        n_row = 2;
      } else {
        IntegerVector sum_site_cleaned = {sum_site[1], sum_site[2]};
        IntegerVector hap_site_cleaned = {hap_site[1], hap_site[2]};
        List out = aux_noN_S2(sum_site_cleaned, hap_site_cleaned);
        haplotype[0] = out["temp"];
        n_row = out["n_row"];
      }
    } else {
      List out = aux_noN_S3(sum_site, hap_site);
      haplotype[0] = out["temp"];
      n_row = out["n_row"];
    }
    // undecided_pos(more_than1) = j;
    // pos_possibility(more_than1++) = n_row[j];
  }
  else if(num == 4) {
    sum = sum_site[3] + sum_site[2] + sum_site[1] + sum_site[0];
    if(hap_site[0] == -1) {
      if(sum_site[0]/sum >= 0.45) {
        arma::mat temp;
        IntegerMatrix inner_tmp(2, NUM_CLASS);
        for(l = 0; l < num; ++l)
          for(m = l + 1; m < num; ++m) {
            if(uni_alignment[ref_j] == "I") //deletion in B genome
              inner_tmp = call_permute_N({hap_site[l], hap_site[m]}, 0);
            else if(uni_alignment[ref_j] == "J")
              inner_tmp = call_permute_N({hap_site[l], hap_site[m]}, 1);
            temp = join_cols(temp, as<arma::mat>(inner_tmp));
          }
          haplotype[0] = wrap(temp);
        n_row = 2;
      } else {
        IntegerVector sum_site_cleaned = {sum_site[1], sum_site[2], sum_site[3]};
        IntegerVector hap_site_cleaned = {hap_site[1], hap_site[2], hap_site[3]};
        List out = aux_noN_S3(sum_site_cleaned, hap_site_cleaned);
        haplotype[0] = out["temp"];
        n_row = out["n_row"];
      }
    } else {
      IntegerMatrix temp = call_permute({hap_site[0], hap_site[1], hap_site[2], hap_site[3]});
      haplotype[0] = temp;
      n_row = 24;
    }
    // undecided_pos(more_than1) = j;
    // pos_possibility(more_than1++) = n_row[j];
  }
  List ls = List::create(
    Named("n_row") = n_row,
    Named("haplotype") = haplotype);
  
  return(ls);
}

/*
 * determine the hidden state by each time t, only use this on the first time t, then, the rest non_overlapped still use site by site
 */
List tbt_state(unsigned int num, unsigned int ref_j, IntegerVector hap_site, IntegerVector sum_site, 
               List hmm_info, List dat_info) {
  
  List n_in_t = hmm_info["n_in_t"];
  IntegerVector n_t = hmm_info["n_t"];
  IntegerVector length = dat_info["length"];
  IntegerVector obs = dat_info["nuc"];
  IntegerVector obs_index = dat_info["id"];
  IntegerVector time_pos = hmm_info["time_pos"];
  IntegerVector p_tmax = hmm_info["p_tmax"];
  IntegerVector index = dat_info["start_id"];
  int hap_max_pos = dat_info["ref_length_max"];
  int hap_min_pos = dat_info["ref_start"];
  int hap_length = hap_max_pos - hap_min_pos;
  unsigned int t_max = hmm_info["t_max"];
  unsigned int qua_in, obs_in;
  // hash the reads at the first time position
  std::string outfile = "./out.fasta";
  FILE *fp = fopen(outfile.c_str(), "w");
  if (!fp)
    stop("Cannot open file");
  
  for (unsigned int i = 0; i < n_t[0]; ++i) {
    int ind = n_in_t[i];
    for (unsigned int j = 0; j < length[ind]; ++j) {
      // qua_in = qua[index[ind] + j];
      obs_in = obs[index[ind] + j];
    }
  }
  
}




