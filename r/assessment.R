### read the truth
sourceCpp("./r/assess.cpp")

error_rates <- function(res = NULL, truth_file) {
  truth <- read_fasta(truth_file)
  true_geno <- truth$reads
  snps <- res$snps
  snp_location <- res$snp_location
  hap_length <- ncol(res$haplotypes$hap_final)
  start_pos <- res$start_pos
  snp_call <- snp_stats(res$snps, snp_location - 1, hap_length, start_pos, true_geno)
  heter_fp <- snp_call$`confusion metric`[2, 2]/sum(snp_call$`confusion metric`[2, ])
  homo_fp <- snp_call$`confusion metric`[1, 1]/sum(snp_call$`confusion metric`[1, ])
  return(snp_call)
}

h0.005_3.0 <- read_rds("~/Documents/hmm_res0")
h0.005_4.0 <- read_rds("~/Downloads/hmm_res0")
h0.005_3.0.snp <- error_rates(res = h0.005_3.0, truth_file = "../../data/hmm/WGS/simu/homr0.005/indiv0.fsa")
h0.005_4.0.snp <- error_rates(res = h0.005_4.0, truth_file = "../../data/hmm/WGS/homr0.005/indiv0.fsa")

pp <- pp_snp(ref1.snp$haplotypes$chosed_state, ref1.snp$combination, ref1.snp$loci, 
             ref1.snp$par$par_hmm_bf$gamma, length(ref1.snp$undecided_pos))
datafile = "../../data/hmm/WGS/simu/homr0.005/cov3/out0.txt"
alignment = "../../data/hmm/WGS/simu/homr0.005/ref.fsa"
ref1.snp = altragenotype(datafile = datafile, alignment = alignment)
ref1 <- error_rates(res = ref1.snp, truth_file = "../../data/hmm/WGS/simu/homr0.005/indiv0.fsa")


