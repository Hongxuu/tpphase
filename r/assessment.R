

### fasta from other methods

get_res_other <- function(parent_path, covergae, individual, name) {
  res_all <- list()
  n_ind = 1
  for(j in individual) {
    res_ind <- list()
    count = 1
    for(i in covergae) {
      coverage_path <- paste0(parent_path, "cov", i, "/", name, "_res")
      res_file <- paste0(coverage_path, "/sim", j, "_", name, ".fa")
      res_ind[[count]] <- res_file
      count = count + 1
    }
    res_all[[n_ind]] <- res_ind
    n_ind = n_ind + 1
  }
  return(res_all)
}


error_rates <- function(res, truth_file, is_hmm, fdr = TRUE) {
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
    if(is_hmm) {
      snp_location <- individual$snp_location
      hap_length <- ncol(individual$haplotypes$hap_final)
      start_pos <- individual$start_pos
      snp_call <- snp_stats(individual$snps, snp_location - 1, hap_length, start_pos, true_geno)
    } else {
      cat(individual, "\n");
      inferred <- read_fasta(individual)
      inferred_geno <- inferred$reads
      snp_call <- snp_stats_other(inferred_hap = inferred_geno, hap_length = ncol(inferred_geno), min_ref = 0, true_hap = true_geno)
    }
    heter_swe <- snp_call$switch$heter_sw_err
    homo_swe <- snp_call$switch$homo_sw_err
    a <- sum(snp_call$`confusion metric`[, 2]) - snp_call$`confusion metric`[2, 2]
    b <- sum(snp_call$`confusion metric`[2, ]) - snp_call$`confusion metric`[2, 2]
    c <- sum(snp_call$`confusion metric`[, 1]) - snp_call$`confusion metric`[1, 1]
    d <- sum(snp_call$`confusion metric`[1, ]) - snp_call$`confusion metric`[1, 1]
    if(fdr) {
      heter_fp <- (a)/(sum(snp_call$`confusion metric`[, 2]))
      homo_fp <- (c)/(sum(snp_call$`confusion metric`[, 1]))
      stats <- list("heter_fdr" = heter_fp, "homo_fdr" = homo_fp,
                    "heter_swe" = heter_swe, "homo_swe" = homo_swe)
    } else {
      heter_fp <- (a)/(a + b)
      homo_fp <- (c)/(c + d)
      stats <- list("heter_fp" = heter_fp, "homo_fp" = homo_fp,
                    "heter_swe" = heter_swe, "homo_swe" = homo_swe)
    }
    
    # stats <- append(stats, snp_call[-3])
    if(heter_swe != 0) {
      heter_msw <- mean(diff(c(0, snp_call$switch$heter_sw_id)))
      stats <- append(stats, list("heter_msw" = heter_msw))
    } else
      stats <- append(stats, list("heter_msw" = NA))
    if(homo_swe != 0) {
      homo_msw <- mean(diff(c(0, snp_call$switch$homo_sw_id)))
      stats <- append(stats, list("homo_msw" = homo_msw))
    } else
      stats <- append(stats, list("homo_msw" = NA))
    value[[i]] <- stats
  }
  if(length(res) == 8) {
    a <- value[[1]]
    value <- a
  }
  return(value)
}

####### get the simulation results


get_res <- function(parent_path, covergae, individual, name) {
  res_all <- list()
  n_ind = 1
  for(j in individual) {
    res_ind <- list()
    count = 1
    for(i in covergae) {
      coverage_path <- paste0(parent_path, "cov", i, "/", name)
      res_file <- paste0(coverage_path, "/", name, j)
      res_ind[[count]] <- read_rds(res_file)
      count = count + 1
    }
    res_all[[n_ind]] <- res_ind
    n_ind = n_ind + 1
  }
  return(res_all)
}

###### err rates

get_err <- function(individual, parent_path, res_all, covergae, is_hmm, verbose = 0) {
  summary <- data.frame()
  for(j in individual) {
    truth_file <- paste0(parent_path, "indiv", j, ".fsa")
    tmp <- error_rates(res = res_all[[j + 1]], truth_file, is_hmm) %>% 
      bind_rows() %>% 
      add_column(coverage = covergae, individual = j)
    print(tmp, "\n")
    summary <- rbind(summary, tmp)
  }
  summary$coverage <- as.factor(summary$coverage)
  df <- summary[, -ncol(summary)] %>% melt(id.vars = "coverage") 
  if(verbose) {
    df %>% 
      ggplot(aes(coverage, value)) +
      geom_boxplot() + 
      facet_wrap(~variable, scales = "free")
  }
  
  return(df)
}

################ pp for snps
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
