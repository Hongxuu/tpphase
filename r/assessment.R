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

####### get the simulation results
parent_path <- "../../../../peanut_simu/homr0.005/"
covergae <- c(3, 4, 8, 12, 16)
individual <- c(0:49)

get_res <- function(parent_path, covergae, individual) {
  res_all <- list()
  n_ind = 1
  for(j in individual) {
    res_ind <- list()
    count = 1
    for(i in covergae) {
      coverage_path <- paste0(parent_path, "cov", i)
      res_file <- paste0(coverage_path, "/hmm_res", j)
      res_ind[[count]] <- read_rds(res_file)
      count = count + 1
    }
    res_all[[n_ind]] <- res_ind
    n_ind = n_ind + 1
  }
  return(res_all)
}
res_all.0.005 <- get_res(parent_path, covergae, individual)

###### err rates

get_err <- function(individual, parent_path, res_all, covergae) {
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

get_err(individual, parent_path, res_all, covergae)
get_err(individual0.005, parent_path, res_all.0.005, covergae)

######## pp for snps
get_pp <- function(individual, res_all, covergae) {
  pp_roc <- data.frame()
  ind = 1
  for(j in individual) {
    ind_res <- res_all[[ind]]
    cov = 1
    for(i in covergae) {
      pp <- pp_snp(ind_res[[cov]]$haplotypes$chosed_state, ind_res[[cov]]$combination, 
                 ind_res[[cov]]$loci, ind_res[[cov]]$par$par_hmm_bf$gamma, ncol(ind_res[[cov]]$snps))
      cov = cov + 1
      pp_roc <- rbind(pp_roc, data.frame("pp_all" = pp, coverage = i, individual = j)) 
    }
    ind = ind + 1
  }
  return(pp_roc)
}
pp_1 <- get_pp(individual, res_all.0.005, covergae)

library(precrec)	# auc from here
library(ROCR)		# plots from here

#### heterozygotes ####
d <- read.table("roc_het.txt", header=T)
boxplot(PP ~ Truth, data = d)	# great separation
em <- evalmod(scores = d$PP, labels = d$Truth, mode = "rocprc")
d.hc <- d[d$CovA>10 & d$CovB>10,]
em <- evalmod(scores = d.hc$PP, labels = d.hc$Truth, mode = "rocprc")

## make plots
pred <- prediction(d$PP, d$Truth)
perf <- performance(pred, "tpr", "fpr")
plot(perf, colorize=T)
perf <- performance(pred, "prec", "rec")
plot(perf, colorize=T)
