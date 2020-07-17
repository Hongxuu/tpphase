### read the truth
library(reshape2)
sourceCpp("./r/assess.cpp")

error_rates <- function(res, truth_file) {
  truth <- read_fasta(truth_file)
  true_geno <- truth$reads
  
  if(length(res) == 8) {
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
    # stats <- append(stats, snp_call[-3])
    if(heter_swe) {
      heter_msw <- mean(diff(c(0, snp_call$switch$heter_sw_id)))
      stats <- append(stats, list("heter_msw" = heter_msw))
    }
    if(homo_swe) {
      homo_msw <- mean(diff(c(0, snp_call$switch$homo_sw_id)))
      stats <- append(stats, list("homo_msw" = homo_msw))
    }
    value[[i]] <- stats
  }
  if(class(res) != "list") {
    a <- value[[1]]
    value <- a
  }
  return(value)
}
parent_path <- "../../../../peanut_simu/homr0.02/"
res_all <- list()
covergae <- c(3, 4, 8, 12, 16)
individual <- c(0:49)
individual0.02 <- individual[-c(5, 11, 12, 13, 19, 20, 21, 32, 33, 24, 26, 40, 47, 25, 
                                27, 34, 35, 36, 41, 48)]
for(j in individual0.02) {
  res_ind <- list()
  count = 1
  for(i in covergae) {
    coverage_path <- paste0(parent_path, "cov", i)
    res_file <- paste0(coverage_path, "/hmm_res", j)
    res_ind[[count]] <- read_rds(res_file)
    count = count + 1
  }
  res_all[[j + 1]] <- res_ind
}
individual0.005 <- individual[-c(10, 11, 12, 13, 19, 36, 42, 45, 48, 49, 50)]
individual0.01 <- individual[-c(15, 25, 41, 7, 5, 6, 8, 27, 32, 43, 47, 48, 50)]
get_res <- function(individual, parent_path, res_all) {
  summary <- data.frame()
  for(j in individual) {
    truth_file <- paste0(parent_path, "indiv", j, ".fsa")
    tmp <- error_rates(res = res_all[[j + 1]], truth_file) %>% 
      bind_rows() %>% 
      add_column(coverage = covergae, individual = j)
    summary <- rbind(summary, tmp)
  }
  summary$coverage <- as.factor(summary$coverage)
  summary[, -ncol(summary)] %>% melt(id.vars = "coverage") %>% 
    ggplot(aes(coverage, value)) +
    geom_boxplot() + 
    facet_wrap(~variable, scales = "free")
}
get_res(individual0.005, parent_path, res_all)

get_res(individual0.02, parent_path, res_all)



