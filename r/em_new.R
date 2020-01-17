library(Rcpp)
library(RcppArmadillo)
library(tidyverse)
library(foreach)
library(mnlogit)
library(doParallel)
library(Formula)

source("./read_data.R")
source("../r/data_process.R")
source("../r/modified_mnlogit.R")
source("../r/m_beta.R")
source("../r/initialization.R")
source("../r/formula.R")
source("../r/newton.R")
source("../r/likelihood.R")
sourceCpp("../r/data_format.cpp")
sourceCpp("../r/mnlogit.cpp")
sourceCpp("../r/e_step.cpp")
sourceCpp("../r/m_hap.cpp")

#' @description The data used filter soft and hard clip
#' @param samfile Input sam file
#' @param ref_name Reference name in sam file
#' @param init Initialization method (Options: "ampliclust", "in_file", "random", Default: "ampliclust")
#' @param FastaFile Input a fasta file if initialization is "in_file"
#' @param ampliclust_command "ampliclust" command (indicate path as well)
#' @param snp Vector or file to indicate variation sites
#' @param fastq_file fastq file output from sam file
#' @param datafile data file output from sam file
#' @param output Indicate if write results in a file
#' @param max Max iteration
#' @param n_class Number of classes in the mixture model
#' @param num_cat Number of categories in the logistic regression
#' @param ncores Number of cores to register
#' @param seed Seed
#' @param tol Convergence tolarance
#' 
#' @usage final <- tpphase(samfile = "./308-TAN-B.sam", ref_name = "Aipa71:5269634_P3", 
#' ampliclust_command = "../amplici/run_ampliclust", output = "308TAN_B_P3.txt")
#' @return assignments and haplotypes, etc

tpphase <- function(samfile = NULL, ref_name = NULL, init = "ampliclust", FastaFile, ampliclust_command, deletion_num = 2,
                    fastq_file = "./res.fastq", datafile = "./res.txt", ac_outfile = "./init", snp = NULL, output = NULL, 
                    formula = formula(mode~1|read_pos + ref_pos + qua + hap_nuc + qua:hap_nuc),
                    n_class = 4, num_cat = 4, seed = 0, max = 50, tol = 1e-06, ncores = 2, n_predictor = 10) {
  
  registerDoParallel(cores = ncores)
  
  d <- read_data(datafile)
  read_length <- d$length
  hap_length <- d$ref_length_max
  
  if(is.null(snp) == TRUE)
    snp <- rep(1, hap_length)
  if(typeof(snp) == "character")
    snp <- read.delim(snp, header = FALSE, sep = " ") %>% as.integer
  else
    snp <- as.integer(snp)
  
  if(is.null(samfile) == FALSE)
    sam <- read_sam(samfile, ref_name, fastq_file, datafile)
  
  if(init == "ampliclust") {
    hapinit <- call_ampliclust(ampliclust_command, fastq_file, ac_outfile)
    hapinit <- hapinit[1:n_class, 1:hap_length]
  }
    
  if(init == "in_file")
    hapinit <- read_fasta(FastaFile)[1:n_class, 1:hap_length]
  
  if(init == "random") {
    set.seed(seed)
    samp <- which(d$fake_length == hap_length & d$deletion$del_length_all <= deletion_num)
    if(length(samp) < n_class)
      stop("Not enough sample with the same length as the haplotypes to infer, 
           adjust the deletion_num to be more than 2!")
    samp_id <- sample(samp, n_class)
    start <- d$start_id[samp_id] #index is right in R!
    hap_deletion_len <- d$deletion$del_length_all[samp_id]
    hap_info <- sample_hap(d, start, samp_id, hap_deletion_len)
    hapinit <- hap_info$hap
    # if (any(d$deletion$del_flag[samp_id] == 1))
    #   N_in = 1
  }
  
  data <- format_data(dat_info = d, haplotype = hapinit)
  data$nuc <- to_char_r(data$nuc)
  data$hap_nuc <- to_char_r(data$hap_nuc)
  id <- data["id"]
  data <- data[, !names(data) %in% c("id")]
  
  par <- list()
  par <- ini(dat = data, n_observation = d$n_observation, formula = formula, seed = seed, n_class = n_class, num_cat = num_cat)
  hap <- hapinit
  
  full_llk <- rep(0, max)
  CE_llk_iter <- rep(0, max)
  haps <- list()
  resu <- list()
  
  data <- cbind(id, data)
  ### Iteration
  for (m in 1:10) {
    cat("iteartion", m, "\n")
    if(nrow(par$beta) > n_predictor) {
      N_in = 1
    } else 
      N_in = 0
    
    #sink(paste0("~/Documents/debug",m,".txt"))
    res <- em_eta(par = par, dat_info = d, haplotype = hap, hap_info = hap_info, PD_LENGTH = nrow(par$beta), N_in = N_in)
    #sink()
    par$wic <- res$param$w_ic
    par$rate <- res$param$rate
    par$excluded_read <- res$param$excluded_read
    resu[[m]] <- res
    if(length(res$excluded_id) != 0)
      cat(res$excluded_id, "don't (doesn't) belong to any of the haplotypes\n")
    full_llk[m] <- res$full_llk
    
    old_hap <- hap
    #sink(paste0("~/Documents/debug_hap",m,".txt"))
    hap_info <- m_hap(par = par, dat_info = d, PD_LENGTH = nrow(par$beta), N_in = N_in, haplotype = old_hap, SNP = snp)
    hap <- hap_info$hap
    haps[[m]] <- hap_info
    #sink()
    if(any(old_hap != hap)) {
      data <- format_data(dat_info = d, haplotype = hap)
      data$nuc <- to_char_r(data$nuc)
      data$hap_nuc <- to_char_r(data$hap_nuc)
    } 
    
    if (nrow(res$param$beta) == nrow(par$beta)) {
      change = 1
    } else 
      change = 0
      
    tmp <- m_beta(res = res, d = d, id = id, data = data, formula = formula, change = change, reads_lengths = read_length, ncores)
    par <- tmp$par
    CE_llk_iter[m] <- tmp$CEllk
    
    #if(abs(par$eta - res$param$mixture_prop) < tol && abs(par$beta - res$param$beta) < tol)
    if(m > 1)
      if(abs(full_llk[m] - full_llk[m-1]) < tol && 
         abs(par$eta - res$param$mixture_prop) < tol && 
         abs(par$beta - res$param$beta) < tol)
        break;
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

final <- tpphase(samfile = "../../../data/peanut_consensus/308-TAN-consensus.sam", ref_name = "Adur313:17300_Adur313:17334_P30", init = "random",
                 ampliclust_command = "../../amplici/run_ampliclust", fastq_file = "../../../data/tpphase_res_consensus/308TAN/resp30.fastq",
                 datafile = "../../../data/tpphase_res_consensus/308TAN/resp30.txt",
                 ac_outfile = "../../../data/tpphase_res_consensus/initp30", 
                 output = "../../../data/tpphase_res_consensus/308TAN_p30.txt")


# hapinit <- readFastq(fastq_file)
# samp <- sample(which(hapinit@sread@ranges@width == hap_length), n_class) # id starts from 1 in data
# a <- sread(hapinit)[samp] %>% as.data.frame()
# reads <- plyr::ldply(apply(a, MARGIN = 1, FUN = function(x) strsplit(x, "")) %>% flatten, rbind)
# reads <- t(t(reads) %>% na.omit)[, -1]
# ncol <- ncol(reads)
# reads_num <- to_xy_r(reads)
# hapinit <- matrix(reads_num, ncol) %>% t()
