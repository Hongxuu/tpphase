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
source("./r/initialization.R")
source("./r/formula.R")
source("./r/newton.R")
source("./r/m_beta.R")
source("./r/likelihood.R")
sourceCpp("./r/data_format.cpp")
sourceCpp("./r/mnlogit.cpp")
sourceCpp("./r/baumwelch.cpp")
#sourceCpp("./r/extra.cpp")
sourceCpp("./r/initialization.cpp")

## prepare data
datafile = "./test_28.txt"
A <- read_fasta("../../data/peanut_consensus/A_aligned_target.fasta")
B <- read_fasta("../../data/peanut_consensus/B_aligned_target.fasta")
universial <- make_universal_old(A, B)
U28 <- universial$universal_alignment[(universial$start_id[8]+1):universial$start_id[9]]
old_version = 0
dat_info <- read_data(datafile, old_v = old_version)
HMM <- hmm_info(dat_info = dat_info, cut_off = 0.45, uni_alignment = U28) ## cut-off too high TODO: dataset

altragnotype <- function(samfile = NULL, ref_name = NULL, 
         fastq_file = "./res.fastq", datafile = "./res.txt", output = NULL, 
         formula = mode~1|read_pos + ref_pos + qua + hap_nuc + qua:hap_nuc, n_initialization = 1,
         n_class = 4, num_cat = 4, seed = 0, max = 20, tol = 1e-06, ncores = 2, old_version = 0)  {
  registerDoParallel(cores = ncores)  
}
## initialize hap
hap_info <- sample_hap2(HMM, dat_info$ref_length_max)
hapinit <- hap_info$haplotype
data <- format_data(dat_info = dat_info, haplotype = hapinit)
weight_id <- NULL
if(hap_info$gap_in) {
  data_rm <- data %>% filter(mode == 1) 
  weight_id <- which(data_rm$hap_nuc == -1)
  data <- data %>% filter(hap_nuc != -1) # mnlogit only takes data without indels in read or in haplotypes
}
data$nuc <- to_char_r(data$nuc)
data$hap_nuc <- to_char_r(data$hap_nuc)
id <- data["id"]
data <- data[, !names(data) %in% c("id")]
par <- list()
par <- ini_par(dat = data, n_observation = dat_info$n_observation, formula = formula, old_version = old_version,
               n_class = n_class, num_cat = num_cat, ncores = ncores, weight_id = weight_id)

## baum-welch
bw <- baum_welch(hmm_info = HMM, data_info = dat_info, hap_info = hap_full, par = par, PD_LENGTH = nrow(par$beta))

###### estimate beta
weight_id <- NULL
data <- format_data2(hmm_info = HMM, d_info = dat_info, hap_info = hap_full)
if(hap_info$gap_in) {
  data_rm <- data %>% filter(mode == 1) 
  weight_id <- which(data_rm$hap_nuc == -1)
  data <- data %>% filter(hap_nuc != -1) # mnlogit only takes data without indels in read or in haplotypes
}
data$nuc <- to_char_r(data$nuc)
data$hap_nuc <- to_char_r(data$hap_nuc)
id <- data["id"]

tmp <- m_beta(res = bw, id = id, weight_id = weight_id, data = data, formula = formula, 
              reads_lengths = read_length, ncores = ncores, old_version = 0, weight = bw$param$weight)










