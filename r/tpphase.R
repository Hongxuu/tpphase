tpphase <- function(dat_info, hap_info, par, hap, old_hap, tol, id, weight_id, 
                    data, formula, read_length, ncores, snp, max) {
  full_llk <- rep(0, max)
  CE_llk_iter <- rep(0, max)
  haps <- list()
  resu <- list()
  for (m in 1:max) {
    cat("Iteartion", m, "\n")
    cat("Infer hidden states\n")
    hap_info <- viterbi_MN(par = par, dat_info = dat_info)
    hidden_state <- hap_info$hidden
    
    res <- em_eta(par = par, dat_info = dat_info, hap_info = hap_info, haplotype = hap, PD_LENGTH = nrow(par$beta))
    resu[[m]] <- res
    if(length(res$excluded_id) != 0)
      cat(res$excluded_id, "don't (doesn't) belong to any of the haplotypes\n")
    full_llk[m] <- res$full_llk
    
    #if(abs(par$eta - res$param$mixture_prop) < tol && abs(par$beta - res$param$beta) < tol)
    if(m > 1)
      if(abs(full_llk[m] - full_llk[m-1]) < tol | 
         (abs(par$eta - res$param$mixture_prop) < tol & 
         abs(par$beta - res$param$beta) < tol))
        break;
    
    if(any(old_hap != hap)) {
      weight_id <- NULL
      data <- format_data(dat_info = dat_info, haplotype = hap)
      if(dat_info$over_hapmax)
        data <- data %>% filter(hap_nuc != -1)
      if(sum(hap_info$hap_deletion_len) != 0) {
        data_rm <- data %>% filter(mode == 1) 
        weight_id <- which(data_rm$hap_nuc == 4)
        data <- data %>% filter(hap_nuc != 4) # mnlogit only takes data without indels in read or in haplotypes
      }
      data$nuc <- to_char_r(data$nuc)
      data$hap_nuc <- to_char_r(data$hap_nuc)
      id <- data["id"]
    } 
    tmp <- m_beta(res = res, id = id, weight_id = weight_id, data = data, formula = formula, 
                  reads_lengths = read_length, ncores)
    par <- tmp$par
    CE_llk_iter[m] <- tmp$CEllk
    
    old_hap <- hap
    cat("Update non-indel positions\n")
    hap <- m_hap(par = par, dat_info = dat_info, PD_LENGTH = nrow(par$beta),
                 haplotype = hap, hidden_state = hidden_state, SNP = snp)
    hap_info$hap <- hap
    haps[[m]] <- hap_info
  }
  r <- list(resu = resu[[m]], haps = haps[[m - 1]])
  return(r)
}