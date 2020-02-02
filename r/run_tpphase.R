library(Rcpp)
library(RcppArmadillo)
library(tidyverse)
library(foreach)
library(mnlogit)
library(doParallel)
library(Formula)

source("./src/read_data.R")
source("./r/data_process.R")
source("./r/tpphase.R")
source("./r/modified_mnlogit.R")
source("./r/m_beta.R")
source("./r/initialization.R")
source("./r/formula.R")
source("./r/newton.R")
source("./r/likelihood.R")
sourceCpp("./r/data_format.cpp")
sourceCpp("./r/mnlogit.cpp")
sourceCpp("./r/e_step.cpp")
sourceCpp("./r/m_hap.cpp")

#' @description The data used filter soft and hard clip
#' @param samfile Input sam file
#' @param ref_name Reference name in sam file
#' @param init Initialization method (Options: "ampliclust", "in_file", "random", Default: "ampliclust")
#' @param fasta_file Input a fasta file if initialization is "in_file"
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

run_tpphase <- function(samfile = NULL, ref_name = NULL, init = "random", fasta_file, ampliclust_command, deletion_cut = 15,
                    fastq_file = "./res.fastq", datafile = "./res.txt", ac_outfile = "./init", snp = NULL, output = NULL, 
                    formula = mode~1|read_pos + ref_pos + qua + hap_nuc + qua:hap_nuc, n_initialization = 1,
                    n_class = 4, num_cat = 4, seed = 0, max = 50, tol = 1e-06, ncores = 2) {
  
  registerDoParallel(cores = ncores)
  
  ## read the data
  if(is.null(samfile) == FALSE)
    sam <- read_sam(samfile, ref_name, fastq_file, datafile)
  
  d <- read_data(datafile)
  hap_length <- d$ref_length_max
  read_length <- d$length
  ## read in non-snps sites
  if(is.null(snp) == TRUE) {
    snp <- rep(1, hap_length)
  } else if(typeof(snp) == "character") {
    snp <- read.delim(snp, header = FALSE, sep = " ") %>% as.integer
  } else {
    snp <- as.integer(snp)
  }
  
  set.seed(seed)
  best_llk <- -Inf
  for (i in 1:n_initialization) {
    ## initialize haplotype
    hap_info <- ini_hap(d, init, ampliclust_command, fastq_file, ac_outfile, n_class, 
                        hap_length, fasta_file, deletion_cut)
    hapinit <- hap_info$hap
    ## prepare data
    data <- format_data(dat_info = d, haplotype = hapinit)
    if(d$over_hapmax)
      data <- data %>% filter(hap_nuc != -1)
    weight_id <- NULL
    if(sum(hap_info$hap_deletion_len) != 0) {
      data_rm <- data %>% filter(mode == 1) 
      weight_id <- which(data_rm$hap_nuc == 4)
      data <- data %>% filter(hap_nuc != 4) # mnlogit only takes data without indels in read or in haplotypes
      #read_length <- len_hapGap(dat_info = d, hap_info = hap_info)
      #read_length <- collect(data %>% filter(mode == 1) %>% count(id) %>% select(n))[[1]]/4
    }
    data$nuc <- to_char_r(data$nuc)
    data$hap_nuc <- to_char_r(data$hap_nuc)
    id <- data["id"]
    data <- data[, !names(data) %in% c("id")]
    
    ## initialize parameters
    par <- list()
    par <- ini_par(dat = data, n_observation = d$n_observation, formula = formula, 
                   n_class = n_class, num_cat = num_cat, ncores = ncores, weight_id = weight_id)
    
    hap <- hapinit
    old_hap <- hap
    data <- cbind(id, data)
    ## Iteration
    results <- tpphase(dat_info = d, hap_info = hap_info, par = par, hap = hap, old_hap = old_hap, tol = tol, 
                       id = id, weight_id = weight_id, data = data, formula = formula, read_length = read_length, 
                       ncores = ncores, snp = snp, max = max)
    cat("Log likelihood in the", n_initialization, "th", "initialization:", results$resu$full_llk, "\n")
    if(best_llk < results$resu$full_llk) {
      best_llk <- results$resu$full_llk
      final_res <- dereplicate_res(resu = results$resu, haps = results$haps, n_class = n_class)
    }
  }
  if(is.null(output) == FALSE)
    fnlist(final_res, output)
  
  return(final_res)
}

final <- run_tpphase(samfile = NULL, ref_name = NULL, 
                     init = "random", deletion_cut = 15, fastq_file = "../../../data/tpphase_res_consensus/308TAN/resp5.fastq",
                     datafile = "../../data/tpphase_res_consensus/308TAN/resp5.txt", snp = NULL,
                     output = "../../data/tpphase_res_consensus/308TAN/308TAN_p5.txt", 
                     formula = mode~1|read_pos + ref_pos + qua + hap_nuc + qua:hap_nuc, n_initialization = 1,
                     n_class = 4, num_cat = 4, seed = 6, max = 50, tol = 1e-06, ncores = 2)

"../"
# hapinit <- readFastq(fastq_file)
# samp <- sample(which(hapinit@sread@ranges@width == hap_length), n_class) # id starts from 1 in data
# a <- sread(hapinit)[samp] %>% as.data.frame()
# reads <- plyr::ldply(apply(a, MARGIN = 1, FUN = function(x) strsplit(x, "")) %>% flatten, rbind)
# reads <- t(t(reads) %>% na.omit)[, -1]
# ncol <- ncol(reads)
# reads_num <- to_xy_r(reads)
# hapinit <- matrix(reads_num, ncol) %>% t()
# ampliclust_command = "../../amplici/run_ampliclust"
# ac_outfile = "../../../data/tpphase_res_consensus/initp30"
