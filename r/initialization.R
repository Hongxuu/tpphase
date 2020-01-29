
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

ini_par <- function(dat, n_observation, formula, n_class, weight_id, num_cat, seed, ncores) {
  par <- list()
  #set.seed(seed)
  #par$eta <- runif(n_class, 0, 1)
  #par$eta <- par$eta/sum(par$eta)
  par$eta <- rep(0.25, n_class)
  par$wic <- matrix(0.25, nrow = n_observation, ncol = n_class)
  par$ins_rate <- 1e-5
  par$del_rate <- 1e-5
  par$excluded_read <- rep(0, n_observation)
  Mpar <- Mstep(dat, read_length = NULL, par, weight_id = weight_id, formula, num_cat, ncores, weights = FALSE)
  par$beta <- Mpar$beta
  return(par)
}


ini_hap <- function(init, ampliclust_command, fastq_file, ac_outfile, n_class, 
                    hap_length, fasta_file, seed, deletion_cut) {
  if(init == "ampliclust") {
    hap_info <- list()
    hapinit <- call_ampliclust(ampliclust_command, fastq_file, ac_outfile)
    hap_info$hap <- hapinit[1:n_class, 1:hap_length]
    hap_info$hap_deletion_len <- 0
  }
  
  if(init == "in_file") {
    hap_info <- list()
    hapinit <- read_fasta(fasta_file)[1:n_class, 1:hap_length]
    hap_info$hap_deletion_len <- 0
    hap_info$hap <- hapinit
  }
  
  if(init == "random") {
    set.seed(seed)
    uni_map <- unique_map(v = d$deletion$del_length_all)
    deletion_num <- uni_sum(uni_map, cut_off = deletion_cut)
    samp <- which(d$fake_length == hap_length & d$deletion$del_length_all <= deletion_num)
    if(length(samp) < n_class)
      stop("Not enough sample with the same length as the haplotypes to infer, 
           adjust the deletion_cut to be larger!")
    samp_id <- sample(samp, n_class)
    start <- d$start_id[samp_id] #index is right in R!
    hap_deletion_len <- d$deletion$del_length_all[samp_id]
    hap_info <- sample_hap(d, start, samp_id, hap_deletion_len)
  }
  return(hap_info)
}
