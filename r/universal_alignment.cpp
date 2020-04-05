#include <Rcpp.h>
using namespace Rcpp;

#define N_GENOME 2
CharacterVector iupac_to_char = {"-", "A", "C", "M", "G", "R", "S", "V", 
                                 "T", "W", "Y", "H", "K", "D", "B", "N"};

// [[Rcpp::export]]
List make_universal(List alignment, unsigned int for_hmm, int ref_idx = -1) { 
  IntegerVector reads = alignment["reads"];
  IntegerVector dim = alignment["dim"];
  
  if(dim.length() == 2) // if only give one alignment
    dim[0] = dim[1];
  unsigned int i, j, m;
  unsigned int n_alignment = dim.length()/N_GENOME;
  List universal(n_alignment);
  IntegerVector start_id(n_alignment);
  for(i = 1; i < n_alignment; ++i)
    start_id[i] = start_id[i - 1] + dim[(i - 1) * N_GENOME] * N_GENOME;
  
  unsigned int idx;
  CharacterVector seq(reads.length());
  
  if(!for_hmm)
    for(m = 0; m < reads.length(); ++m)
      seq[m] = iupac_to_char[reads[m]];
  
  for (i = 0; i < n_alignment; ++i) {
    if(i != ref_idx && for_hmm)
      continue;
    
    CharacterVector uni_alignment(dim[i * N_GENOME]);
    for(j = 0; j < dim[i * N_GENOME]; ++j) {
      idx = j + start_id[i];
      if(for_hmm) {
        if(reads[idx] == 0 && reads[idx + dim[i * N_GENOME]] != 0) {
          uni_alignment[j] = "J";
        } else if(reads[idx] != 0 && reads[idx + dim[i * N_GENOME]] == 0) {
          uni_alignment[j] = "I";
        } else
          uni_alignment[j] = "M";
      } else {
        // Rcout << idx << "\t" << idx + dim[i * N_GENOME] << "\n";
        if(reads[idx] == reads[idx + dim[i * N_GENOME]]) {
          uni_alignment[j] = seq[idx];
        } else {
          uni_alignment[j] = "N";
        }
      }
    }
    universal[i] = uni_alignment;
  }
  return universal;
}

//[[Rcpp::export]]
void write_fasta(List seq, StringVector names, std::string outfile = "./out.fasta") {
  unsigned int i, j;
  FILE *fp = fopen(outfile.c_str(), "w");
  if (!fp)
    stop("Cannot open file");
  
  for(i = 0; i < seq.size(); ++i) {
    putc('>', fp);
    for(j = 0; j < names[i].size(); ++j)
      fprintf(fp, "%c", names[i][j]);
    fprintf(fp, "\n");
    
    CharacterVector read = seq[i];
    for(j = 0; j < read.size(); ++j)
      fprintf(fp, "%s", static_cast<const char*>(read[j]));
    fprintf(fp, "\n");
  }
  fclose(fp);
}

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

// whole genome data

// chr2 <- read_fasta("../../data/peanut_consensus/chr21to1.fasta")
//   uni <- make_universal(alignment = chr2, for_hmm = 0)
//   write_fasta(seq = uni, names = paste0("chr2.", seq(1:length(uni))) , outfile = "./uni.fasta")
// comman_path = "~/data/peanut/"
// primary_dirs = list.files(comman_path)
//   for(i in 1:length(primary_dirs)) {
//     path <- paste0(comman_path, "ch", i, "/", "ch", i)
//     seq <- read_fasta(datafile = paste0(path, "_1to1.fasta"))
//     uni <- make_universal(alignment = seq, for_hmm = 0)
//     write_fasta(seq = uni, names = paste0("ch", i, ".", seq(1:length(uni))), outfile = paste0(path, "_uni.fasta"))
//   }

// targeted data
// A <- read_fasta("../../data/peanut_consensus/A_aligned_target.fasta")
// B <- read_fasta("../../data/peanut_consensus/B_aligned_target.fasta")
// universial <- make_universal_old(A, B)
// U28 <- universial$universal_alignment[(universial$start_id[8]+1):universial$start_id[9]]
