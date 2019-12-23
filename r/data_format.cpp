#include <Rcpp.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

using namespace Rcpp;

#define MLOGIT_CLASS 4
#define NUM_CLASS 4
#define PD_LENGTH 10

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

// 
// par init_par(NumericVector beta) {
//   par *par;
//   par->beta = NULL;
//   par->eta = NULL;
//   par->w_ic = NULL;
//   MAKE_1DOUBLEARRAY(par->eta, 4);
//   MAKE_1DOUBLEARRAY(par->beta, (MLOGIT_CLASS - 1) * PD_LENGTH);
//   
//   unsigned int i;
//   double sum = 0.;
//   
//   for (i = 0; i < 4; ++i) {
//     par->eta[i] = (double) rand() / RAND_MAX;
//     sum += par->eta[i];
//   }
//   for(i = 0; i < beta.length(); ++i) {
//     par->beta[i] = beta[i];
//   }
//   return par;
// }
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

// FIND UNIQUE IF THE SEQUENCE IS ORDERED
// int uniq(int *a, unsigned int len);
// int uniq(int *a, unsigned int len)
// {
//   int i, j;
//   for (i = j = 0; i < len; i++)
//     if (a[i] != a[j]) a[++j] = a[i];
//     return j + 1;
// }

// FIND UNIQUE UNDER MY CASE
// List uniq(IntegerVector a, unsigned int len) {
//   IntegerVector unique(len);
//   unsigned int uniq_len;
//   for (unsigned int m = 0; m < len; ++m) 
//     if (m < len - 1 && a[m] != a[m + 1])
//       unique[uniq_len++] = a[m];
//     if(unique[uniq_len] != a[len])
//       unique[uniq_len++] = a[len];
//     
//   return ls; 
// }

// [[Rcpp::export]]
List read_data(std::string path) {
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
  if(ins_id[ins_num] != ins_obs_index[count_ins])
    ins_id[ins_num++] = ins_obs_index[count_ins];
    
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
  for (i = 0; i < n_observation; ++i)
    if (ins_flag[i] == 1)
      for (m = 0; m < ins_num; ++m)
        if (ins_id[m] == i + 1)
          ins_length_all[i] = ins_count[m];
        
  // deletion
  for (m = 0; m < count_del; ++m)
    if (m < count_del - 1 && del_obs_index[m] != del_obs_index[m + 1])
      del_id[del_num++] = del_obs_index[m];
  if(del_id[del_num] != del_obs_index[count_del])
      del_id[del_num++] = del_obs_index[count_del];
  
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
  for (i = 0; i < n_observation; ++i)
    if (del_flag[i] == 1)
      for (m = 0; m < del_num; ++m)
        if (del_id[m] == i + 1)
          del_length_all[i] = del_count[m];
  
  // non indel
  for (m = 0; m < count; ++m) {
    i = obs_index[m];
    // if (!length[i-1]) {
    //   n_observation++;
    // }
    ++length[i-1];
    obs[m] = char_to_xy(obs[m]);
  }
  
  // true length
  IntegerVector true_length(n_observation);
  IntegerVector fake_length(n_observation);
  for (i = 0; i < n_observation; ++i) {
    true_length[i] = length[i] + ins_length_all[i];
    fake_length[i] = length[i] + del_length_all[i];
  }
  
  /* index of read */
  IntegerVector index(n_observation);
  for (m = 1; m < n_observation; ++m)
      index[m] = index[m - 1] + length[m - 1];
  
  int max_len = 0;
  for(i = 0; i < total; ++i)
    if(ref_pos[i] > max_len)
      max_len = ref_pos[i]; // Actual length need to plus 1
    
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
  
  IntegerVector non_covered_site(max_len);
  unsigned int num;
  for (m = 0; m < max_len; ++m) {
    num = 0;
      for(i = 0; i < n_observation; ++i)
        if (ref_index(i, m) == -1)
          num++;
      if (num == n_observation)
        non_covered_site[m] = 1;
  }
        
  List del = List::create(
    Named("del_id") = del_id[Range(0, del_num - 1)],
    Named("del_id_all") = del_obs_index,
    Named("del_flag") = del_flag,
    Named("del_read_pos") = del_read_pos,
    Named("del_ref_pos") = del_ref_pos,
    Named("del_length") = del_count, //no. of deletion in each read that contains deletion
    Named("del_num") = del_num,  //no. of reads have deletion
    Named("del_total") = count_del,
    Named("del_length_all") = del_length_all);
  
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
    Named("non_covered_site") = non_covered_site,
    Named("deletion") = del,
    Named("insertion") = ins);
  
  return ls;
}

// [[Rcpp::export]]
DataFrame format_data(List dat_info, IntegerMatrix haplotype) {
  unsigned int i, k, l;
  int input_arr[] = {0, 1, 2, 3};
  size_t input_arr_sz = sizeof input_arr / sizeof *input_arr;
  
  int total = dat_info["total"];
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
      for (l = 0; l < MLOGIT_CLASS; ++l)
        r_hap_nuc[i * MLOGIT_CLASS * NUM_CLASS + MLOGIT_CLASS * k + l] = haplotype(k, ref_pos[i]);

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
IntegerMatrix sample_hap (List dat_info, IntegerVector start, IntegerVector idx) {
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
  
  IntegerMatrix hap_nuc(NUM_CLASS, hap_length);
  
  for (i = 0; i < NUM_CLASS; ++i) {
    for (j = 0; j < length[idx[i]]; ++j)
      hap_nuc(i, ref_pos[start[i] + j]) = obs[start[i] + j];
    if (del_flag[idx[i]] == 1)
      for (m = 0; m < del_total; ++m)
        if (del_id_all[m] == idx[i])
          hap_nuc(i, del_ref_pos[m]) = -1;
  }
  return hap_nuc;
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
















