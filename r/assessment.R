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


ref6.mid.snp <- confusion_matric(res = ref6.mid, truth_file = "../../data/hmm/WGS/simu/ref6/indiv0.fsa")