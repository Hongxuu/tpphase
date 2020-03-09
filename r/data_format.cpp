#include <RcppArmadillo.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

using namespace Rcpp;
using namespace std;

#define MLOGIT_CLASS 4
#define NUM_CLASS 4

// [[Rcpp::depends(RcppArmadillo)]]

//https://stackoverflow.com/questions/51110244/in-rcpp-how-to-get-a-user-defined-structure-from-c-into-r
// namespace Rcpp {
//   template <>
//   SEXP wrap(const par& x);
// }
// 
// 
// 
// namespace Rcpp {
//   template <>
//   SEXP wrap(const par& x) {
//     Rcpp::NumericVector beta;
//     return Rcpp::wrap(Rcpp::Named("beta") = Rcpp::wrap(x.beta));
//   };
// }

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

typedef unsigned char xy_t;
/**
 * Convert char to xy.
 *
 * @param c	ASCII char
 * @return	xy_t
 */
inline xy_t char_to_xy(char c)
{
  if (c == 'A' || c == 'a')
    return 'A' >> 1 & 3L;
  else if (c == 'C' || c == 'c')
    return 'C' >> 1 & 3L;
  else if (c == 'G' || c == 'g')
    return 'G' >> 1 & 3L;
  else if (c == 'T' || c == 't' || c == 'U' || c == 'u')
    return 'T' >> 1 & 3L;
  else
    return 1 << 7;	/* non-nuc */
} /* char_to_xy */

unsigned char const xy_to_char[MLOGIT_CLASS] = {'A', 'C', 'T', 'G'};

// [[Rcpp::export]]
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

// [[Rcpp::export]] 
int uni_sum(List unique_map, unsigned int cut_off)
{
  NumericVector result_vals = unique_map["lengths"];
  NumericVector result_keys = unique_map["values"];
  int n = result_keys.size();
  int cumsum = 0;
  int key = 0;
  for (unsigned int i = 0; i < n; ++i) {
    cumsum += result_vals(i);
    if(cumsum >= cut_off) {
      key = result_keys(i);
      break;
    }
  }
  return key;
}

// [[Rcpp::export]]
int top_n_map(List unique_map) 
{
  NumericVector result_vals = unique_map["lengths"];
  NumericVector result_keys = unique_map["values"];
  int n = result_keys.size();
  int key = 0;
  for (unsigned int i = n; i --> 0;)
    if(result_vals(i) >= NUM_CLASS) {
      key = result_keys(i);
      break;
    }
    return key;
}

// [[Rcpp::export]]
List read_data(std::string path, unsigned int old_v) {
  int j, k, l;
  unsigned int i, m, count = 0, count_del = 0, count_ins = 0;
  char c;
  char str[100];
  int n_observation = 0;
  
  FILE* fp = fopen(path.c_str(), "r");
  if(!fp) {
    Rcpp::stop("File opening failed");
  }
  int temp_n = -1;
  while (fgets(str, sizeof(str), fp)) {
    std::sscanf(str, "%d %d %d %d %c", &i, &j, &k, &l, &c);
    /* exclude -1 in ref_pos and qua */
    if (k == -1) {
      count_ins++;
      continue;
    }
    else if (l == -1)
      count_del++;
    else {
      if (temp_n != i) {
        temp_n = i;
        n_observation++;
      }
      count++;
    }
  }
  
  unsigned int del_num = 0;
  IntegerVector del_obs_index(count_del);
  IntegerVector del_ref_pos(count_del);
  IntegerVector del_read_pos(count_del);
  IntegerVector del_id(n_observation);
  IntegerVector del_flag(n_observation);
  
  unsigned int ins_num = 0;
  IntegerVector ins_obs_index(count_ins);
  IntegerVector ins_id(n_observation);
  IntegerVector ins_read_pos(count_ins);
  IntegerVector ins_flag(n_observation);
  
  IntegerVector qua(count);
  IntegerVector obs(count);
  StringVector obs_str(count);
  IntegerVector obs_index(count);
  IntegerVector ref_pos(count);
  IntegerVector read_pos(count);
  int total;
  
  rewind(fp);
  
  count = 0;
  count_del = 0;
  count_ins = 0;
  while (fgets(str, sizeof(str), fp)) {
    std::sscanf(str, "%d %d %d %d %c", &i, &j, &k, &l, &c);
    /* exclude -1 in ref_pos and qua */
    if (k == -1) {
      ins_obs_index[count_ins] = i;
      ins_read_pos[count_ins] = j;
      ins_flag[i-1] = 1;
      count_ins++;
    }
    else if (l == -1) {
      del_obs_index[count_del] = i;
      del_read_pos[count_del] = j;
      del_ref_pos[count_del] = k;
      del_flag[i-1] = 1;
      count_del++;
    }
    else {
      obs_index[count] = i;
      read_pos[count] = j;
      ref_pos[count] = k;
      qua[count] = l;
      obs[count] = c;
      count++;
    }
  }
  fclose(fp);
  fp = NULL;
  
  total = count;
  IntegerVector length(n_observation);
  
  //insertion
  for (m = 0; m < count_ins; ++m) 
    if (m < count_ins - 1 && ins_obs_index[m] != ins_obs_index[m + 1])
      ins_id[ins_num++] = ins_obs_index[m];
  if(ins_id[ins_num- 1] != ins_obs_index[count_ins - 1])
    ins_id[ins_num++] = ins_obs_index[count_ins - 1];
    
  IntegerVector ins_count(ins_num);
  for (m = 0; m < ins_num; ++m)
    ins_count[m] = 1;
  i = 0;
  for (m = 0; m < count_ins; ++m) {
    if (m < count_ins - 1 && ins_obs_index[m] == ins_obs_index[m + 1])
      ins_count[i]++;
    else
      i++;
  }
  
  IntegerVector ins_length_all(n_observation);
  for (m = 0; m < ins_num; ++m)
        ins_length_all[ins_id[m] - 1] = ins_count[m];
        
  // deletion
  for (m = 0; m < count_del; ++m)
    if (m < count_del - 1 && del_obs_index[m] != del_obs_index[m + 1])
      del_id[del_num++] = del_obs_index[m];
  if (del_id[del_num - 1] != del_obs_index[count_del - 1])
      del_id[del_num++] = del_obs_index[count_del - 1];
  
  IntegerVector del_count(del_num);
  for (m = 0; m < del_num; ++m)
    del_count[m] = 1;
  i = 0;
  for (m = 0; m < count_del; ++m) {
    if (m < count_del - 1 && del_obs_index[m] == del_obs_index[m + 1])
      del_count[i]++;
    else
      i++;
  }
  
  IntegerVector del_length_all(n_observation);
  for (m = 0; m < del_num; ++m)
        del_length_all[del_id[m] - 1] = del_count[m];
  
  IntegerVector del_strat_id(n_observation);
  for (i = 1; i < n_observation; ++i)
    if(del_flag[i] == 1)
      for (j = 0; j < i; ++j)
        del_strat_id[i] += del_length_all[j];
    else
      del_strat_id[i] = -1;
  
  // non indel
  for (m = 0; m < count; ++m) {
    i = obs_index[m];
    // if (!length[i-1]) {
    //   n_observation++;
    // }
    ++length[i-1];
    obs[m] = char_to_xy(obs[m]);
  }
  
  /* index of read in non-indel loci */
  IntegerVector index(n_observation);
  for (m = 1; m < n_observation; ++m)
    index[m] = index[m - 1] + length[m - 1];
  
  // true length (include insertion) and fake length (include deletion) (some reads' alignment does not start from 0, other one also count that in)
  IntegerVector true_length(n_observation);
  IntegerVector fake_length(n_observation);
  for (i = 0; i < n_observation; ++i) {
    true_length[i] = length[i] + ins_length_all[i];
    fake_length[i] = length[i] + del_length_all[i];
    //fake_length[i] = length[i] + del_length_all[i] + ref_pos[index[i]];
  }
  
  int max_len = 0;
  int over_hapmax = 0; // indicate if the length of reads is more than the hap_max
  if(old_v == 1) {
    /* find the longest reference position && appears more than no. of classes (4) */
    List uni_map = unique_map(ref_pos);
    max_len = top_n_map(uni_map);
    for(i = 0; i < total; ++i)
      if (max_len < ref_pos[i]) {
        over_hapmax = 1;
        break;
      }
  } else {
    for(i = 0; i < total; ++i)
      if (max_len < ref_pos[i])
        max_len = ref_pos[i];
  }
  max_len = max_len + 1;
  
  IntegerMatrix ref_index(n_observation, max_len);
  for(i = 0; i < n_observation; ++i)
     for(m = 0; m < max_len; ++m)
      ref_index(i, m) = -1;
    
  /* Find the index of ref position under different read of every j */
  for(i = 0; i < n_observation; ++i)
    for(m = 0; m < max_len; ++m)
      for(j = 0; j < length[i]; ++j)
        if(ref_pos[index[i] + j] == m) {
          ref_index(i, m) = j;
          break;
        }
  
  // IntegerVector non_covered_site(max_len);
  // unsigned int num;
  // for (m = 0; m < max_len; ++m) {
  //   num = 0;
  //     for(i = 0; i < n_observation; ++i)
  //       if (ref_index(i, m) == -1)
  //         num++;
  //     if (num == n_observation)
  //       non_covered_site[m] = 1;
  // }
        
  List del = List::create(
    Named("del_id") = del_id[Range(0, del_num - 1)],
    Named("del_id_all") = del_obs_index,
    Named("del_flag") = del_flag,
    Named("del_read_pos") = del_read_pos,
    Named("del_ref_pos") = del_ref_pos,
    Named("del_length") = del_count, //no. of deletion in each read
    Named("del_num") = del_num,  //no. of reads have deletion
    Named("del_total") = count_del,
    Named("del_length_all") = del_length_all,
    Named("del_strat_id") = del_strat_id);
  
  List ins = List::create(
    Named("ins_id") = ins_id[Range(0, ins_num - 1)],
    Named("ins_read_pos") = ins_read_pos,
    Named("ins_length") = ins_count,
    Named("ins_num") = ins_num,
    Named("ins_flag") = ins_flag,
    Named("ins_length_all") = ins_length_all);
  
  List ls = List::create(
    Named("id") = obs_index,
    Named("read_pos") = read_pos,
    Named("ref_pos") = ref_pos,
    Named("nuc") = obs,
    Named("qua") = qua,
    Named("n_observation") = n_observation,
    Named("length") = length, 
    Named("true_length") = true_length, 
    Named("fake_length") = fake_length, 
    Named("total") = total,
    Named("start_id") = index,
    Named("ref_length_max") = max_len, 
    Named("ref_idx") = ref_index, 
    Named("over_hapmax") = over_hapmax,
    // Named("non_covered_site") = non_covered_site,
    Named("deletion") = del,
    Named("insertion") = ins);
  
  return ls;
}

//prepare data for mnlogit, mnlogit only takes data without indels in read or in haplotypes
// [[Rcpp::export]]
DataFrame format_data(List dat_info, IntegerMatrix haplotype) {
  unsigned int i, k, l;
  int input_arr[] = {0, 1, 2, 3};
  size_t input_arr_sz = sizeof input_arr / sizeof *input_arr;
  
  int total = dat_info["total"];
  int hap_length = dat_info["ref_length_max"];
  IntegerVector qua = dat_info["qua"];
  IntegerVector obs = dat_info["nuc"];
  IntegerVector obs_index = dat_info["id"];
  IntegerVector ref_pos = dat_info["ref_pos"];
  IntegerVector read_pos = dat_info["read_pos"];

  IntegerVector r_ref_pos(total * MLOGIT_CLASS * NUM_CLASS);
  IntegerVector r_read_pos(total * MLOGIT_CLASS * NUM_CLASS);
  IntegerVector r_qua(total * MLOGIT_CLASS * NUM_CLASS);
  IntegerVector r_obs(total * MLOGIT_CLASS * NUM_CLASS);
  IntegerVector r_hap_nuc(total * MLOGIT_CLASS * NUM_CLASS);
  IntegerVector mode(total * MLOGIT_CLASS * NUM_CLASS);
  IntegerVector id(total * MLOGIT_CLASS * NUM_CLASS);

  /* pick out the haplotypes accordingg to ref_pos */
  for (i = 0; i < total; ++i)
    for (k = 0; k < NUM_CLASS; ++k)
      for (l = 0; l < MLOGIT_CLASS; ++l) {
        // the reference position might be longer than the sampled haplotypes
        if(ref_pos[i] > hap_length - 1)
          r_hap_nuc[i * MLOGIT_CLASS * NUM_CLASS + MLOGIT_CLASS * k + l] = -1;
        else
          r_hap_nuc[i * MLOGIT_CLASS * NUM_CLASS + MLOGIT_CLASS * k + l] = haplotype(k, ref_pos[i]);
      }
      
  /* repeat the data for mnlogit */
  for (k = 0; k < total; ++k)
    for (i = 0; i < MLOGIT_CLASS * NUM_CLASS; ++i) {
      r_ref_pos[k * MLOGIT_CLASS * NUM_CLASS + i] = ref_pos[k];
      r_read_pos[k * MLOGIT_CLASS * NUM_CLASS + i] = read_pos[k];
      r_qua[k * MLOGIT_CLASS * NUM_CLASS + i] = qua[k];
      id[k * MLOGIT_CLASS * NUM_CLASS + i] = obs_index[k];
    }

    for (i = 0; i < total * MLOGIT_CLASS * NUM_CLASS; ++i)
      r_obs[i] = input_arr[i % input_arr_sz];

  for (i = 0; i < total; ++i)
    for (k = 0; k < MLOGIT_CLASS * NUM_CLASS; ++k)
      if (r_obs[k + MLOGIT_CLASS * NUM_CLASS * i] == obs[i])
        mode[k + MLOGIT_CLASS * NUM_CLASS * i] = 1;

  DataFrame df_new = DataFrame::create(
    Named("id") = id,
    Named("mode") = mode,
    Named("read_pos") = r_read_pos,
    Named("ref_pos") = r_ref_pos,
    Named("qua") = r_qua,
    Named("nuc") = r_obs,
    Named("hap_nuc") = r_hap_nuc);

  return(df_new);
}

// [[Rcpp::export]]
List sample_hap (List dat_info, IntegerVector start, IntegerVector idx, IntegerVector hap_deletion_len) {
  List deletion = dat_info["deletion"];
  IntegerVector del_flag = deletion["del_flag"];
  IntegerVector del_id_all = deletion["del_id_all"];
  IntegerVector del_ref_pos = deletion["del_ref_pos"];
  IntegerVector ref_pos = dat_info["ref_pos"];
  IntegerVector obs = dat_info["nuc"];
  IntegerVector length = dat_info["length"];
  int del_total = deletion["del_total"];
  int hap_length = dat_info["ref_length_max"];
  unsigned int i, j, m;
  unsigned int sum = 0;
  unsigned int count = 0;
  IntegerMatrix hap_nuc(NUM_CLASS, hap_length);
  
  for(i = 0; i < NUM_CLASS; ++i)
    sum += hap_deletion_len[i];
  IntegerVector hap_ref_pos(sum);
  IntegerVector strat_id(NUM_CLASS);
  //Rcout << start << "\n";
  for (i = 0; i < NUM_CLASS; ++i) {
    //Rprintf("\n%d\n", i);
    for (j = 0; j < length[idx[i] - 1]; ++j)  {
      //Rprintf("%d\t", obs[start[i] + j]);
      hap_nuc(i, ref_pos[start[i] + j]) = obs[start[i] + j]; //idx start from 1!!!
    }
    if (del_flag[idx[i] - 1])
      for (m = 0; m < del_total; ++m)
        if (del_id_all[m] == idx[i]) {
          hap_nuc(i, del_ref_pos[m]) = 4;
          hap_ref_pos[count++] = del_ref_pos[m];
        }
  }
  
  for (i = 0; i < NUM_CLASS; ++i) {
    if (hap_deletion_len[i] != 0) {
      if(i != 0)
        for (j = 0; j < i; ++j)
          strat_id[i] += hap_deletion_len[j];
    }
    else
      strat_id[i] = -1; //record the deletion starting id in vector hap_ref_pos for each hap
  }
  
  List ls = List::create(
    Named("hap") = hap_nuc,
    Named("deletion_pos") = hap_ref_pos,
    Named("hap_deletion_len") = hap_deletion_len,
    Named("hap_del_start_id") = strat_id);
  return ls;
}

// [[Rcpp::export]] 

/* NOT WORKING??? */
CharacterVector to_char(IntegerVector nuc, int length) {
  int i = 0;
  CharacterVector nuc_char(length);
  for(i = 0; i < length; ++i) {
    nuc_char[i] = xy_to_char[nuc[i]];
  }
  return(nuc_char);
}

// [[Rcpp::export]] 
IntegerVector find_snp (CharacterMatrix hap) {
  unsigned int k, j, flag, count = 0;
  IntegerVector snp_idx(hap.ncol());
  for(j = 0; j < hap.ncol(); ++j) {
    flag = 0;
    for(k = 1; k < hap.nrow(); ++k) {
      if(hap(k, j) != hap(0, j))
        flag = 1;
    }
    if(flag)
      snp_idx(count++) = j;
  }
  return snp_idx[Range(0, count - 1)];
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

List aux_noN(IntegerVector sum_site, IntegerVector hap_site, unsigned int which_stage) {
  List out(2);
  int n_row, sum;
  IntegerMatrix temp;
  if (which_stage == 2) {
    sum = sum_site[1] + sum_site[0];
    if(sum_site[0]/sum <= 0.35) { //1st one appears 3 times
      temp = call_permute({hap_site[0], hap_site[1], hap_site[1], hap_site[1]});
      n_row = 4;
    } if(sum_site[0]/sum >= 0.65) {//2nd one appears 3 times
      temp = call_permute({hap_site[0], hap_site[0], hap_site[0], hap_site[1]});
      n_row = 4;
    } else {
      temp = call_permute({hap_site[0], hap_site[0], hap_site[1], hap_site[1]});
      n_row = 6;
    }
  }
  if(which_stage == 3) {
    sum = sum_site[2] + sum_site[1] + sum_site[0];
    if(sum_site[0]/sum >= 0.45) {
      temp = call_permute({hap_site[0], hap_site[0], hap_site[1], hap_site[2]});
    } if(sum_site[1]/sum >= 0.45) {
      temp = call_permute({hap_site[0], hap_site[1], hap_site[1], hap_site[2]});
    } else {
      temp = call_permute({hap_site[0], hap_site[1], hap_site[2], hap_site[2]});
    }
    n_row = 12;
  }
  List ls = List::create(
    Named("n_row") = n_row,
    Named("temp") = temp);
  return ls; 
}


// [[Rcpp::export]] 
List hmm_info(List dat_info, double cut_off, CharacterVector uni_alignment) {
  unsigned int i, j, t, l, m;
  int n_observation = dat_info["n_observation"];
  unsigned int total = dat_info["total"];
  unsigned int hap_length = dat_info["ref_length_max"];
  IntegerMatrix ref_index = dat_info["ref_idx"];
  IntegerVector ref_pos = dat_info["ref_pos"];
  IntegerVector index = dat_info["start_id"];
  IntegerVector fake_length = dat_info["fake_length"];
  IntegerVector length = dat_info["length"];
  IntegerVector obs = dat_info["nuc"];
  /* Find the number of reads with alignment start from each t (hash) */
  IntegerVector read_start(n_observation);
  for(i = 0; i < n_observation; ++i)
    read_start(i) = ref_pos[index[i]];
  List start_t = unique_map(read_start);
  IntegerVector n_t = start_t["lengths"];
  IntegerVector time_pos = start_t["values"];
  
  // get the coverage of each site (including -)
  List deletion = dat_info["deletion"];
  IntegerVector del_ref_pos = deletion["del_ref_pos"];
  unsigned int del_total = deletion["del_total"];
  unsigned int enuma = total + del_total;
  IntegerVector pos(enuma);
  for(i = 0; i < total; ++i)
    pos(i) = ref_pos(i);
  for(i = 0; i < del_total; ++i)
    pos(i + total) = del_ref_pos(i);
  List coverage_stat = unique_map(pos);
  IntegerVector coverage = coverage_stat["lengths"];
  NumericVector prop(hap_length);
  for(i = 0; i < hap_length; ++i)
    prop[i] = coverage[i] * cut_off;
  
  /* Find the max p in each t and record the observations in each set t*/
  IntegerVector which_t(n_observation);
  unsigned int t_max = time_pos.length();
  IntegerVector p_tmax(t_max);
  unsigned int max;
  for(t = 0; t < t_max; ++t) {
    max = 0;
    for(i = 0; i < n_observation; ++i)
      if(ref_pos[index[i]] == time_pos[t]) {
        which_t[i] = t;
        if(max < fake_length[i])
          max = fake_length[i];
      }
    p_tmax[t] = max;
  }
  
  unsigned int count, num, more_than1 = 0;
  IntegerVector nuc_j(n_observation);
  List nuc_unique(hap_length);
  List nuc_count(hap_length);
  List nuc(2);
  NumericVector keys(MLOGIT_CLASS + 1);
  NumericVector vals(MLOGIT_CLASS + 1);
  List haplotype(hap_length);
  IntegerVector n_row(hap_length);
  IntegerVector pos_possibility(hap_length);
  IntegerVector undecided_pos(hap_length);
  for (j = 0; j < hap_length; ++j) {
    count = 0;
    for (i = 0; i < n_observation; ++i)
      if((j >= ref_pos[index[i]]) && (j <= ref_pos[index[i] + length[i] - 1])) {
        if(ref_index(i, j) == -1) {
          nuc_j(count++) = -1;
        } else {
          nuc_j(count++) = obs[index[i] + ref_index(i, j)];
        }
      }
    /* get the nuc table at each site */  
    nuc = unique_map(nuc_j[Range(0, count - 1)]);
    IntegerVector key = nuc["values"];
    IntegerVector val = nuc["lengths"];
    num = 0;
    /* remove some unlikely occurred nuc */
    for(t = 0; t < val.length(); ++t) {
      if(val[t] >= prop(j)) {
        keys(num) = key(t);
        vals(num++) = val(t);
      }
    }
    nuc_unique[j] = keys[Range(0, num - 1)];
    nuc_count[j] = vals[Range(0, num - 1)];
    IntegerVector hap_site = nuc_unique[j];
    IntegerVector sum_site = nuc_count[j];
    unsigned int sum;
    //record possible hidden states[only suitable for ployploids]
    if (num == 1) { //TODO: there is a potential problem that at one site, most of the nuc are -1, but the haplotype cannot have 4 gaps
      IntegerVector temp = {hap_site[0], hap_site[0], hap_site[0], hap_site[0]};
      haplotype(j) = temp;
      n_row[j] = 1;
    } 
    else if (num == 2) {
      if(hap_site[0] == -1) { // if N appears, 2 possibilities, but based on universal, 1 possible
        IntegerVector temp(NUM_CLASS);
        if(uni_alignment[j] == "I") 
          temp = {-1, -1, hap_site[1], hap_site[1]};
        if(uni_alignment[j] == "J")
          temp = {hap_site[1], hap_site[1], -1, -1};
        else
          temp = {hap_site[1], hap_site[1], hap_site[1], hap_site[1]};
        haplotype(j) = temp;
        n_row[j] = 1;
      } else { // if N not appears
        List out = aux_noN(sum_site, hap_site, 2);
        haplotype(j) = out["temp"];
        n_row[j] = out["n_row"];
        undecided_pos(more_than1) = j;
        pos_possibility(more_than1++) = n_row[j];
      }
    } 
    else if (num == 3) {
      sum = sum_site[2] + sum_site[1] + sum_site[0];
      if (hap_site[0] == -1) {
        // if N appears, 4 possibilities
        if(sum_site[0]/sum >= 0.45) {
          IntegerVector temp(2 * NUM_CLASS);
          if(uni_alignment[j] == "I")
            temp = {-1, -1, hap_site[1], hap_site[2], //J in universal alignment
                    -1, -1, hap_site[2], hap_site[1]};
          if(uni_alignment[j] == "J")
            temp = {hap_site[1], hap_site[2], -1, -1,  //I in universal alignment                                       
                    hap_site[2], hap_site[1], -1, -1};
          temp.attr("dim") = Dimension(2, 4);
          haplotype(j) = temp;
          n_row[j] = 2;
        } else {
          IntegerVector sum_site_cleaned = {sum_site[1], sum_site[2]};
          IntegerVector hap_site_cleaned = {hap_site[1], hap_site[2]};
          List out = aux_noN(sum_site_cleaned, hap_site_cleaned, 2);
          haplotype(j) = out["temp"];
          n_row[j] = out["n_row"];
        }
      } else {
        List out = aux_noN(sum_site, hap_site, 3);
        haplotype(j) = out["temp"];
        n_row[j] = out["n_row"];
      }
      undecided_pos(more_than1) = j;
      pos_possibility(more_than1++) = n_row[j];
    }
    else if(num == 4) {
      sum = sum_site[3] + sum_site[2] + sum_site[1] + sum_site[0];
      if(hap_site[0] == -1) {
        if(sum_site[0]/sum >= 0.45) {
          arma::mat temp;
          IntegerMatrix inner_tmp(2, NUM_CLASS);
          for(l = 0; l < num; ++l)
            for(m = l + 1; m < num; ++m) {
              if(uni_alignment[j] == "I") //deletion in B genome
                inner_tmp = call_permute_N({hap_site[l], hap_site[m]}, 0);
              if(uni_alignment[j] == "J")
                inner_tmp = call_permute_N({hap_site[l], hap_site[m]}, 1);
              temp = join_cols(temp, as<arma::mat>(inner_tmp));
            }
            haplotype(j) = wrap(temp);
          n_row[j] = 2;
        } else {
          IntegerVector sum_site_cleaned = {sum_site[1], sum_site[2], sum_site[3]};
          IntegerVector hap_site_cleaned = {hap_site[1], hap_site[2], hap_site[3]};
          List out = aux_noN(sum_site_cleaned, hap_site_cleaned, 3);
          haplotype(j) = out["temp"];
          n_row[j] = out["n_row"];
        }
      } else {
        IntegerMatrix temp = call_permute({hap_site[0], hap_site[1], hap_site[2], hap_site[3]});
        haplotype(j) = temp;
        n_row[j] = 24;
      }
      undecided_pos(more_than1) = j;
      pos_possibility(more_than1++) = n_row[j];
    }
    // else if(num == 5) {//which is not likely
    //   arma::mat temp;
    //   IntegerMatrix inner_tmp(2, NUM_CLASS);
    //   for(l = 0; l < num; ++l)
    //     for(m = l + 1; m < num; ++m) {
    //       if(uni_alignment[j] == "I") //deletion in B genome
    //         inner_tmp = call_permute_N({hap_site[l], hap_site[m]}, 0);
    //       if(uni_alignment[j] == "J")
    //         inner_tmp = call_permute_N({hap_site[l], hap_site[m]}, 1);
    //       temp = join_cols(temp, as<arma::mat>(inner_tmp));
    //     }
    //   haplotype(j) = wrap(temp);
    //   n_row[j] = 2;
    // }
  }
  // find the number of hidden states at each t
  IntegerVector num_states(t_max, 1);
  for(t = 0; t < t_max; ++t)
    for(j = time_pos[t]; j < time_pos[t] + p_tmax[t]; ++j)
      num_states[t] *= n_row[j];
  
 // compute the number of hidden states
  List ls = List::create(
    Named("undecided_pos") = undecided_pos[Range(0, more_than1 - 1)],
    Named("pos_possibility") = pos_possibility[Range(0, more_than1 - 1)], //number of possible snp site
    Named("n_row") = n_row, // number of possible genome at each site
    Named("num_states") = num_states,// number of states at each t
    Named("hidden_states") = haplotype,
    Named("nuc_unique") = nuc_unique, // unique nuc
    Named("nuc_count") = nuc_count,
    Named("p_tmax") = p_tmax, // max length at each t
    Named("n_t") = n_t,// no. emission reads at each time t
    Named("time_pos") = time_pos, // emission start position
    Named("which_t") = which_t, // indicate reads in which t 
    Named("t_max") = t_max); //biggest t
  
  return ls;
}

// 
// IntegerMatrix len_hapGap(List dat_info, List hap_info) {
//   IntegerVector length = dat_info["length"];
//   IntegerVector ref_index = dat_info["ref_idx"];
//   int n_observation = dat_info["n_observation"];
//   
//   IntegerVector hap_ref_pos = hap_info["deletion_pos"];
//   IntegerVector hap_deletion_len = hap_info["hap_deletion_len"];
//   IntegerVector strat_id = hap_info["hap_del_start_id"];
//   unsigned int i, k, j, loci;
//   IntegerMatrix new_len(n_observation, NUM_CLASS);
//   IntegerVector cumsum_del_len(hap_deletion_len.size());
//   // cumsum_del_len = cumsum(hap_deletion_len);
//   std::partial_sum(hap_deletion_len.begin(), hap_deletion_len.end(), cumsum_del_len.begin());
//   
//   for(i = 0; i < n_observation; ++i)
//     for (k = 0; k < NUM_CLASS; ++k) {
//       new_len(i, k) = length(i);
//       if (strat_id(k) != -1)
//         for (j = 0; j < hap_deletion_len[k]; ++j) {
//           if (k != 0) {
//             loci = hap_ref_pos[j + cumsum_del_len[k - 1]];
//           } else
//             loci = hap_ref_pos[j];
//           if(ref_index(i, loci) != -1)
//             new_len(i, k)--;
//         }
//     }
//   return new_len;
// }












