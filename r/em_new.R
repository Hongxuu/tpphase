library(Rcpp)
library(RcppArmadillo)
library(tidyverse)
library(mnlogit)
library(foreach)
library(doParallel)
library(Formula)

source("./read_data.R")
source("~/Documents/Karin/GWAS/code/r/modified_mnlogit.R")
source("~/Documents/Karin/GWAS/code/r/EM.R")
source("~/Documents/Karin/GWAS/code/r/formula.R")
source("~/Documents/Karin/GWAS/code/r/newton.R")
source("~/Documents/Karin/GWAS/code/r/likelihood.R")
sourceCpp("~/Documents/Karin/GWAS/code/r/data_format.cpp")
sourceCpp("~/Documents/Karin/GWAS/code/r/mnlogit.cpp")
sourceCpp("~/Documents/Karin/GWAS/code/r/e_step.cpp")
sourceCpp("~/Documents/Karin/GWAS/code/r/m_hap.cpp")


#' @description 
#' @usage final <- tpphase(samfile = "../../data/tetraploid/308-TAN-B.sam", ref_name = "Aipa71:5269634_P3", 
#' ampliclust_command = "../amplici/run_ampliclust", output = "308TAN_B_P3.txt")
#' @return assignments and haplotypes, etc

tpphase <- function(samfile = NULL, ref_name = NULL, init = "ampliclust", FastaFile, ampliclust_command,
                    fastq_file = "res.fastq", datafile = "res.txt", output = NULL, n_class = 4, 
                    num_cat = 4, seed = 0, max = 50, tol = 1e-06, ncores = 2, ini_iter = 5) {
  
  registerDoParallel(cores = ncores)
  
  if(is.null(samfile) == FALSE)
    sam <- read_sam(samfile, ref_name, fastq_file, datafile)
  
  d <- read_data(datafile)
  read_length <- d$length
  hap_length <- d$ref_length_max
  
  if(init == "ampliclust") {
    hapinit <- call_ampliclust(ampliclust_command, fastq_file)
    hapinit <- hapinit[1:n_class, 1:hap_length]
  }
    
  if(init == "in_file")
    hapinit <- read_fasta(FastaFile)[1:n_class, 1:hap_length]
  
  if(init == "random") { ## length shoud be the same as the longest read, if only one longest not applicable
    set.seed(seed)
    hapinit <- readFastq(fastq_file)
    a <- sread(hapinit)[sample(which(hapinit@sread@ranges@width == hap_length), n_class)] %>% as.data.frame()
    reads <- plyr::ldply(apply(a, MARGIN = 1, FUN = function(x) strsplit(x, "")) %>% flatten, rbind)
    reads <- t(t(reads) %>% na.omit)[, -1]
    ncol <- ncol(reads)
    reads_num <- to_xy_r(reads)
    hapinit <- matrix(reads_num, ncol) %>% t()
    #hap <- read_fastq(fastq_file) // TODO: read with different length 
    #hap <- hap$reads[sample(1:nrow(hap$read), n_class), ]
  }
  
  data <- fromat_data(dat_info = d, haplotype = hapinit)
  data$nuc <- to_char_r(data$nuc)
  data$hap_nuc <- to_char_r(data$hap_nuc)
  id <- data["id"]
  data <- data[, !names(data) %in% c("id")]
  
  par <- list()
  par <- ini(dat = data, n_observation = d$n_observation, seed = seed)
  old_hap <- hapinit
  hap <- hapinit
  
  full_llk <- rep(0, max)
  CE_llk_iter <- rep(0, max)
  haps <- list()
  resu <- list()
  
  data <- cbind(id, data)
  ### Iteration
  for (m in 1:max) {
    cat("iteartion", m, "\n")
    #if(flag == 0) {
      res <- em_eta(par = par, dat_info = d, haplotype = hap)
      resu[[m]] <- res
    #} else {
      #flag = 0
      #resu[[m]] <- res
    #}
    
    if(length(res$excluded_id) != 0)
      cat(res$excluded_id, "don't belong to any of the haplotypes\n")
    full_llk[m] <- res$full_llk
    #if(abs(par$eta - res$param$mixture_prop) < tol && abs(par$beta - res$param$beta) < tol)
    if(m > 1)
      if(abs(full_llk[m] - full_llk[m-1]) < tol && 
         abs(par$eta - res$param$mixture_prop) < tol && 
         abs(par$beta - res$param$beta) < tol)
        break;
    
    if(any(old_hap != hap)) {
      data <- fromat_data(dat_info = d, haplotype = hap)
      data$nuc <- to_char_r(data$nuc)
      data$hap_nuc <- to_char_r(data$hap_nuc)
    } 
    
    tmp <- m_beta(res = res, d = d, id = id, data = data, reads_lengths = read_length, ncores)
    par <- tmp$par
    old_hap <- hap
    hap <- m_hap(par, d, haplotype = old_hap)
    haps[[m]] <- hap
    CE_llk_iter[m] <- tmp$CEllk
  }
  
  final_res <- list()
  final_res$full_llk <- resu[[m]]$full_llk
  final_res$mixture_prop <- resu[[m]]$param$mixture_prop
  final_res$assignments <- apply(resu[[m]]$param$w_ic, MARGIN = 1, FUN = which.max)
  final_res$logistic_coeff <- resu[[m]]$param$beta
  final_res$haplotypes <- matrix(to_char_r(haps[[m-1]]), nrow = n_class)
  
  if(is.null(output) == FALSE)
    fnlist(final_res, output)
  
  return(final_res)
}

# old_id <- 0
# flag <- 0
# init_llk <- rep(0, ini_iter)
# CE_llk <- rep(0, ini_iter)
# change <- rep(0, ini_iter)
# 
# for (n in 1:max) {
#   res <- em_eta(par = par, dat_info = d, haplotype = hap)
#   init_llk[n] <- res$full_llk
#   if(all(res$excluded_id == old_id) && (n >= ini_iter)) {
#     flag <- 1
#     break;
#   } else if (any(res$excluded_id != old_id)){
#     change[n] <- 1
#   }
#   old_id <- res$excluded_id
#   use_d = 1
#   if(all(old_hap == hap)) ## Maybe it is slower than format data (write a C comparision func)
#     use_d = 0
#   if(use_d) {
#     data <- fromat_data(dat_info = d, haplotype = hap)
#     data$nuc <- as.character(c("0" = "A", "2" = "T", "1" = "C", "3" = "G")[as.character(data$nuc)])
#     data$hap_nuc <- as.character(c("0" = "A", "2" = "T", "1" = "C", "3" = "G")[as.character(data$hap_nuc)])
#   } else {
#     data <- cbind(id, data)
#   }
#   
#   tmp <- m_beta(res = res, d = d, id = id, data = data, reads_lengths = read_length, ncores)
#   par <- tmp$par
#   print(sum(log(par$eta)))
#   old_hap <- hap
#   hap <- m_hap(par, d, haplotype = old_hap)
#   CE_llk[n] <- tmp$CEllk
# }