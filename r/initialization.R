
#' @description Random initialize parameters. Uniform(0, 1) for now to initialize etas, but for betas
#' use mnlogit to initialize since it is important to get roughly good betas to compute llk. According to 
#' Hongxu interaction term associated with position are not that significant compared with other terms, so 
#' exclude them for now, so we have 10 beta under each category
#' 
#' @param n_class number of class (A1 A2 D1 D2 here).
#' @param num_cat number of category (for multinomial logistic regression ATCG here)
#' @param num_beta number of betas
#' @param seed
#' @return parameters for EM

ini_par <- function(dat, n_observation, formula, n_class, weight_id, num_cat, ncores) {
  par <- list()
  #set.seed(seed)
  #par$eta <- runif(n_class, 0, 1)
  #par$eta <- par$eta/sum(par$eta)
  par$eta <- rep(1/num_cat, n_class)
  par$wic <- matrix(1/num_cat, nrow = n_observation, ncol = n_class)
  par$ins_rate <- 1e-5
  par$del_rate <- 1e-5
  par$excluded_read <- rep(0, n_observation)
  Mpar <- Mstep(dat, read_length = NULL, par, weight_id = weight_id, formula, num_cat, ncores, weights = FALSE)
  par$beta <- Mpar$beta
  return(par)
}


ini_hap <- function(d, init, ampliclust_command, fastq_file, ac_outfile, n_class, 
                    hap_length, fasta_file, deletion_cut) {
  hap_infom <- list()
  if(init == "ampliclust") {
    hapinit <- call_ampliclust(ampliclust_command, fastq_file, ac_outfile)
    hap_infom$hap <- hapinit[1:n_class, 1:hap_length]
    hap_infom$hap_deletion_len <- 0
  }
  
  if(init == "in_file") {
    hapinit <- read_fasta(fasta_file)[1:n_class, 1:hap_length]
    hap_infom$hap_deletion_len <- 0
    hap_infom$hap <- hapinit
  }
  
  if(init == "random") {
    uni_map <- unique_map(v = d$deletion$del_length_all)
    if(uni_map$lengths[1] >= deletion_cut) {
      deletion_num <- uni_map$values[1]
    } else 
      deletion_num <- uni_sum(uni_map, cut_off = deletion_cut)
    samp <- which(d$fake_length == hap_length & 
                    d$length >= (hap_length - median(uni_map$values)) & 
                    d$deletion$del_length_all <= deletion_num)
    if(length(samp) < n_class * 3)
      stop("Not enough sample with the same length as the haplotypes to infer, 
           adjust the deletion_cut to be larger!\n")
    samp_id <- sample(samp, n_class)
    cat("Selected reads: ", samp_id, "\n");
    start <- d$start_id[samp_id] #index is right in R!
    hap_deletion_len <- d$deletion$del_length_all[samp_id]
    hap_infom <- sample_hap(d, start, samp_id, hap_deletion_len)
  }
  return(hap_infom)
}
