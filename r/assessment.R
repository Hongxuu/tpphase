### read the truth
sourceCpp("./r/assess.cpp")

confusion_matric <- function(res = NULL, truth_file) {
  truth <- read_fasta(truth_file)
  true_geno <- truth$reads
  snps <- res$snps
  snp_location <- res$snp_location
  hap_length <- ncol(res$haplotypes)
  start_pos <- res$start_pos
  snp_call <- snp_stats(res$snps, snp_location - 1, hap_length, start_pos, true_geno)

  return(snp_call)
}
