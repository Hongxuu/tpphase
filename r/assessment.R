

### fasta from other methods

get_res_other <- function(parent_path, covergae, individual, name, is_pair = 0) {
  res_all <- list()
  n_ind = 1
  for(j in individual) {
    res_ind <- list()
    count = 1
    for(i in covergae) {
      if(is_pair)
        coverage_path <- paste0(parent_path, "pair/cov", i, "/", name, "_res")
      else
        coverage_path <- paste0(parent_path, "pair/cov", i, "/", name, "_res")
      res_file <- paste0(coverage_path, "/sim", j, "_", name, ".fa")
      res_ind[[count]] <- res_file
      count = count + 1
    }
    res_all[[n_ind]] <- res_ind
    n_ind = n_ind + 1
  }
  return(res_all)
}

error_rates <- function(res, truth_file, datfile, is_hmm, fdr = TRUE, old = 1) {
  truth <- read_fasta(truth_file)
  true_geno <- truth$reads
  
  if(length(res) == 9) {
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
      if(old) {
        snp_call <- snp_stats(individual$snps, snp_location - 1, hap_length, start_pos, true_geno)
      } else {
        datafile <- datfile[[i]]
        dat_info <- read_data(datafile, old_v = 0)
        snp_call <- sw_hmm(inferred_snp = individual$snps, snp_location = individual$snp_location - 1, 
             hap_length = ncol(individual$haplotypes$hap_final), individual$start_pos, true_geno, dat_info)
        }
    } else {
      cat(individual, "\n");
      inferred <- read_fasta(individual)
      inferred_geno <- inferred$reads
      if(old) {
        snp_call <- snp_stats_other(inferred_hap = inferred_geno, hap_length = ncol(inferred_geno),
                                    min_ref = 0, true_hap = true_geno)
      } else {
        datafile <- datfile[[i]]
        dat_info <- read_data(datafile, old_v = 0)
        snp_call <- sw_other(dat_info = dat_info, inferred_hap = inferred_geno, hap_length = ncol(inferred_geno),
               min_ref = 0, true_hap = true_geno)
        }
    }
    heter_swe <- snp_call$switch$heter_sw_err
    homo_swe <- snp_call$switch$homo_sw_err
    if(!old) {
      len <- homo_swe %>% na.omit() %>% length()
      if(len != 0)
        homo_swe <- homo_swe %>% na.omit() %>% mean()
      else
        homo_swe <- NA
      len <- heter_swe %>% na.omit() %>% length()
      if(len != 0)
        heter_swe <- heter_swe %>% na.omit() %>% mean()
      else
        heter_swe <- NA
    }
    
    a <- sum(snp_call$`confusion metric`[, 2]) - snp_call$`confusion metric`[2, 2]
    b <- sum(snp_call$`confusion metric`[2, ]) - snp_call$`confusion metric`[2, 2]
    c <- sum(snp_call$`confusion metric`[, 1]) - snp_call$`confusion metric`[1, 1]
    d <- sum(snp_call$`confusion metric`[1, ]) - snp_call$`confusion metric`[1, 1]
    if(fdr) {
      heter_fp <- (a)/(sum(snp_call$`confusion metric`[, 2]))
      homo_fp <- (c)/(sum(snp_call$`confusion metric`[, 1]))
      stats <- list("heter_fdr" = heter_fp, "homeo_fdr" = homo_fp,
                    "heter_swe" = heter_swe, "homeo_swe" = homo_swe)
    } else {
      heter_fp <- (a)/(a + b)
      homo_fp <- (c)/(c + d)
      stats <- list("heter_fp" = heter_fp, "homeo_fp" = homo_fp,
                    "heter_swe" = heter_swe, "homeo_swe" = homo_swe)
    }
    
    # stats <- append(stats, snp_call[-3])
    if(old) {
      if(!is.na(heter_swe)) {
        heter_msw <- mean(diff(c(0, snp_call$switch$heter_sw_id)))
        stats <- append(stats, list("heter_msw" = heter_msw))
      } else {
        stats <- append(stats, list("heter_msw" = NA))
      }
      if(!is.na(homo_swe)) {
        homo_msw <- mean(diff(c(0, snp_call$switch$homo_sw_id)))
        stats <- append(stats, list("homeo_msw" = homo_msw))
      } else {
        stats <- append(stats, list("homeo_msw" = NA))
      }
    }
    
    value[[i]] <- stats
  }
  if(length(res) == 8) {
    a <- value[[1]]
    value <- a
  }
  return(value)
}

####### get the simulation results
# res_file = "../../../../peanut_simu/homr0.005/cov16/hmm_res/hmm_res45"
# truth_file = "../../../../peanut_simu/homr0.005/indiv45.fsa"
# read_rds(res_file) -> individual
# gatk <- gatk.0.005[[16]][[5]]


get_res <- function(parent_path, covergae, individual, name, is_pair = 0) {
  res_all <- list()
  n_ind = 1
  for(j in individual) {
    res_ind <- list()
    count = 1
    for(i in covergae) {
      if(is_pair)
        coverage_path <- paste0(parent_path, "pair/cov", i, "/", name)
      else
        coverage_path <- paste0(parent_path, "single/cov", i, "/", name)
      res_file <- paste0(coverage_path, "/", name, j)
      res_ind[[count]] <- read_rds(res_file)
      count = count + 1
    }
    res_all[[n_ind]] <- res_ind
    n_ind = n_ind + 1
  }
  return(res_all)
}

################ pp for snps
# get_pp <- function(individual, res_all, covergae) {
#   pp_roc <- data.frame()
#   ind = 1
#   for(j in individual) {
#     ind_res <- res_all[[ind]]
#     cov = 1
#     for(i in covergae) {
#       pp <- pp_snp(ind_res[[cov]]$haplotypes$chosed_state, ind_res[[cov]]$combination, 
#                  ind_res[[cov]]$loci, ind_res[[cov]]$par$par_hmm_bf$gamma, ncol(ind_res[[cov]]$snps))
#       cov = cov + 1
#       pp_roc <- rbind(pp_roc, data.frame("pp_all" = pp, coverage = i, individual = j)) 
#     }
#     ind = ind + 1
#   }
#   return(pp_roc)
# }

############
iu_to_char_r <- function(x) {
  as.character(c("1" = "A", "8" = "T", "2" = "C", "4" = "G", "16" = "N")[as.character(x)])
}

get_mec <- function(datfile, res, is_hmm) {
  if(length(res) == 9) {
    a <- list()
    a[[1]] <- res
    res <- a
  }
  value <- c()
  for(i in 1:length(res)) {
    individual <- res[[i]]
    datafile <- datfile[[i]]
    dat_info <- read_data(datafile, old_v = 0)
    if(!is_hmm) {
      inferred <- read_fasta(individual)
      haps <- matrix(iu_to_char_r(inferred$reads), nrow = 4)
      cov_record = -1;
    } else {
      individual <- readRDS(individual)
      haps = individual$haplotypes$hap_final
      haps <- matrix(to_char_r(haps), nrow = 4)
      cov_record = individual$cov_record - dat_info$ref_start
    }
    mec <- MEC(dat_info, haps, cov_record)
    value[i] <- mec
  }
  return(value)
}

## all error metric
get_err <- function(individual, parent_path, res_all, covergae, is_hmm, 
                    datfile_name = NULL, verbose = 0, old = 0, compute_mec = 1, is_pair = 0) {
  if(!is.null(datfile_name)) {
    datfile <- list()
    resfile <- list()
    n_ind = 1
    for(j in individual) {
      d_ind <- list()
      res_ind <- list()
      count = 1
      for(i in covergae) {
        if(is_pair) {
          coverage_path <- paste0(parent_path, "pair/cov", i, "/", datfile_name, "_res")
        } else {
          coverage_path <- paste0(parent_path, "single/cov", i, "/", datfile_name, "_res")
        }
        data_f <- paste0(coverage_path, "/out", j, datfile_name, ".txt")
        if(datfile_name == "hmm") {
          res_f <- paste0(coverage_path, "/", datfile_name, "_res", j)
          res_ind[[count]] <- res_f
        }
        d_ind[[count]] <- data_f
        count = count + 1
      }
      datfile[[n_ind]] <- d_ind
      if(datfile_name == "hmm")
        resfile[[n_ind]] <- res_ind
      n_ind = n_ind + 1
    }
    if(datfile_name == "gatk")
      resfile = res_all
  }
  
  summary <- data.frame()
  for(j in individual) {
    truth_file <- paste0(parent_path, "indiv", j, ".fsa")
   
    tmp <- error_rates(res = res_all[[j + 1]], truth_file, datfile = datfile[[j + 1]], is_hmm, old = old) %>% 
      bind_rows() %>% 
      add_column(coverage = covergae, individual = j)
    if(compute_mec) {
      mec <- get_mec(datfile = datfile[[j + 1]], res = resfile[[j + 1]], is_hmm)
      tmp <- tmp %>% add_column("mec" = mec, .before = 1)
    }
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


load_res <- function(parent_path, covergae, individual, datfile_name, true_sw, is_pair, compute_mec = 1) {
  old = true_sw
  if(datfile_name == "hmm") {
    name = "hmm_res"
    is_hmm = 1
    res_all <- get_res(parent_path = parent_path, covergae = covergae, individual = individual, name = name)
    
  } else if(datfile_name == "gatk"){
    name = "gatk"
    is_hmm = 0
    res_all <- get_res_other(parent_path = parent_path, covergae = covergae, individual = individual, name = name)
  }
  res <- get_err(individual = individual, parent_path = parent_path, res_all, covergae = covergae, 
                 datfile_name = datfile_name, is_hmm = is_hmm, old = old, is_pair = is_pair, compute_mec = compute_mec)
  return(res)
}



