### read the truth
sourceCpp("./r/assess.cpp")

error_rates <- function(res, truth_file) {
  truth <- read_fasta(truth_file)
  true_geno <- truth$reads
  
  if(class(res) != "list") {
    a <- list()
    a[[1]] <- res
    res <- a
  }
  value <- list()
  for(i in 1:length(res)) {
    individual <- res[[i]]
    snps <- individual$snps
    snp_location <- individual$snp_location
    hap_length <- ncol(individual$haplotypes$hap_final)
    start_pos <- individual$start_pos
    snp_call <- snp_stats(individual$snps, snp_location - 1, hap_length, start_pos, true_geno)
    heter_fp <- (sum(snp_call$`confusion metric`[2, ]) - snp_call$`confusion metric`[2, 2])/sum(snp_call$`confusion metric`[2, ])
    homo_fp <- (sum(snp_call$`confusion metric`[1, ]) - snp_call$`confusion metric`[1, 1])/sum(snp_call$`confusion metric`[1, ])
    heter_swe <- snp_call$switch$heter_sw_err
    homo_swe <- snp_call$switch$homo_sw_err
    stats <- list("heter_fp" = heter_fp, "homo_fp" = homo_fp,
                  "heter_swe" = heter_swe, "homo_swe" = homo_swe)
    stats <- append(stats, snp_call[-3])
    if(heter_swe) {
      heter_msw <- mean(diff(snp_call$switch$heter_sw_id))
      stats <- append(stats, list("heter_msw" = heter_msw))
    }
    if(heter_swe) {
      homo_swe <- mean(diff(snp_call$switch$homo_sw_id))
      stats <- append(stats, list("heter_msw" = homo_swe))
    }
    value[[i]] <- stats
  }
  if(class(res) != "list") {
    a <- value[[1]]
    value <- a
  }
  return(value)
}
parent_path <- "~/Documents/tmp/"
res_all <- list()
for(j in c(0:2)) {
  res_ind <- list()
  count = 1
  for(i in c(3, 4)) {
    coverage_path <- paste0(parent_path, "cov", i)
    res_file <- paste0(coverage_path, "/hmm_res", j)
    res_ind[[count]] <- read_rds(res_file)
    count = count + 1
  }
  res_all[[j + 1]] <- res_ind
}

ind0 <- error_rates(res = res_all[[1]], truth_file = "../../data/hmm/WGS/simu/homr0.005/indiv0.fsa")
for(i in 1:length(ind0)) {
  res <- ind0[[i]]
}



