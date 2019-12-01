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
// [[Rcpp::export]]
List read_data(std::string path) {
  int j, k, l;
  unsigned int i, m, count = 0;
  char c;
  char str[100];
  FILE* fp = fopen(path.c_str(), "r");
  if(!fp) {
    Rcpp::stop("File opening failed");
  }
  while (fgets(str, sizeof(str), fp)) {
    std::sscanf(str, "%d %d %d %d %c", &i, &j, &k, &l, &c);
    /* exclude -1 in ref_pos and qua */
    if (k == -1 || l == -1) {
      continue;
    }
    count++;
  }
  
  IntegerVector qua(count);
  IntegerVector obs(count);
  IntegerVector obs_index(count);
  IntegerVector ref_pos(count);
  IntegerVector read_pos(count);
  int total;
  int n_observation = 0;
  
  rewind(fp);
  
  count = 0;
  while (fgets(str, sizeof(str), fp)) {
    std::sscanf(str, "%d %d %d %d %c", &i, &j, &k, &l, &c);
    /* exclude -1 in ref_pos and qua */
    if (k == -1 || l == -1) {
      continue;
    }
    obs_index[count] = i;
    read_pos[count] = j;
    ref_pos[count] = k;
    qua[count] = l;
    obs[count] = c;
    count++;
  }
  total = count;
  
  IntegerVector length(total);
  
  for (m = 0; m < count; ++m) {
    i = obs_index[m];
    if (!length[i-1]) {
      n_observation++;
    }
    ++length[i-1];
    obs[m] = char_to_xy(obs[m]);
  }
  fclose(fp);
  fp = NULL;
  
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
        
  // DataFrame df = DataFrame::create(
  //   Named("id") = obs_index,
  //   Named("read_pos") = read_pos,
  //   Named("ref_pos") = ref_pos,
  //   Named("nuc") = obs,
  //   Named("qua") = qua);
  
  List ls = List::create(
    Named("id") = obs_index,
    Named("read_pos") = read_pos,
    Named("ref_pos") = ref_pos,
    Named("nuc") = obs,
    Named("qua") = qua,
    Named("n_observation") = n_observation,
    Named("length") = length[Range(0, n_observation - 1)], 
    Named("total") = total,
    Named("start_id") = index,
    Named("ref_length_max") = max_len, 
    Named("ref_idx") = ref_index, 
    Named("non_covered_site") = non_covered_site);
  
  return ls;
}


// NumericVector read_fastq(std::string path){
//   unsigned int i;
//   FILE* fp = path.c_str();
//   fastq_data *fdata;		/* fasta data */
//   fastq_options fop = {.read_encoding = XY_ENCODING};	/* fastq file options */
//   read_fastq(fp, fdata, &fop);
//   
//   NumericVector haplotype(NUM_CLASS * fdata->n_max_length);
//   for(i = 0; i < NUM_CLASS * fdata->n_max_length; ++i)
//     haplotype[i] = fdata->reads[i];
//   
//   haplotype.attr("dim") = Dimension(NUM_CLASS, fdata->n_max_length);
//   
//   return haplotype;
// }

// [[Rcpp::export]]
DataFrame fromat_data(List dat_info, IntegerMatrix haplotype) {
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

/* NOT WORKING??? */
CharacterVector to_char(IntegerVector nuc, int length) {
  int i = 0;
  CharacterVector nuc_char(length);
  
  for(i = 0; i < length; ++i) {
    nuc_char[i] = xy_to_char[nuc[i]];
  }
  return(nuc_char);
}
















