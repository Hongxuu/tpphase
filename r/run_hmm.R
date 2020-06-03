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
sourceCpp("./r/universal_alignment.cpp")
sourceCpp("./r/initialization.cpp")
sourceCpp("./r/viterbi.cpp")
cut_off = 0.1
formula = mode~1|read_pos + ref_pos + qua + hap_nuc + qua:hap_nuc
n_class = 4
num_cat = 4
seed = 0
tol = 1e-06
ncores = 2
old_version = 0
altragenotype <- function(ref_name = NULL, alignment = NULL, ref_delim = ".", datafile = NULL, output = NULL, cut_off = 0.1,
                          formula = mode~1|read_pos + ref_pos + qua + hap_nuc + qua:hap_nuc, max_iter = 50, res_file = NULL,
                          n_class = 4, num_cat = 4, seed = 0, tol = 1e-06, ncores = 2, old_version = 0, eta = rep(0.25, 4))  {
  registerDoParallel(cores = ncores)  
  ## make universial reference
  align <- read_fasta(alignment)
  if(nrow(align$reads) != 1) { #only read in one reference pair
    universial <- make_universal(alignment = align, for_hmm = 1, ref_idx = 0)
    universial <- universial %>% unlist()
  } else { ## read in many references, take the one we want
    ref_in <- strsplit(ref_name, ref_delim, fixed = TRUE) %>% unlist()
    ref_index <- ref_in[2] %>% as.integer() - 1 ## index in C
    universial <- make_universal(alignment = align, for_hmm = 1, ref_idx = ref_index)
    universial <- Filter(Negate(is.null), universial) %>% unlist()
  }
  rm(align)
  
  ## prepare data
  dat_info <- read_data(datafile, old_v = old_version)
  HMM <- hmm_info(dat_info = dat_info, cut_off = cut_off, uni_alignment = universial)
  
  ########################## baum-welch (iterate until converge)
  
  ## initialization
  ##### use linkage info to limit some unlikelily happened transition
  hap_length <- dat_info$ref_length_max - dat_info$ref_start
  linkage_in <- linkage_info(dat_info = dat_info, undecided_pos = HMM$undecided_pos)
  overlap_info <- get_overlap(HMM$p_tmax, HMM$time_pos, HMM$num_states, HMM$undecided_pos, HMM$t_max, dat_info$ref_start)
  hap_full_info <- full_hap_new(HMM, linkage_in, overlap_info, hap_length, dat_info$ref_start)
  hap_full <- hap_full_info$full_hap
  HMM$num_states <- hap_full_info$new_num_states
  
  ## initialize hap
  set.seed(seed)
  hap_info <- sample_hap2(hmm_info = HMM, hap_length = hap_length, hap_min_pos = dat_info$ref_start)
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
  # 
  ###indicate which transfer could happen
  trans_indicator <- trans_permit(num_states = HMM$num_states, combination = hap_full_info$combination, 
                                  loci = overlap_info$location, t_max = HMM$t_max)
  ### start initializing
  ###### estimate beta
  weight_id <- NULL
  par$eta <- eta
  cat("initialization: \n");
  data_new <- format_data2(hmm_info = HMM, d_info = dat_info, hap_info = hap_full)
  data <- data_new$df_new
  bw <- baum_welch_init(hmm_info = HMM, data_info = dat_info, hap_info = hap_full, par = par, 
                        PD_LENGTH = nrow(par$beta), trans_indicator = trans_indicator, hash_idx = data_new$idx)
  rm(hap_full_info)
  rm(trans_indicator)
  if(hap_info$gap_in) {
    data_rm <- data %>% filter(mode == 1) 
    weight_id <- which(data_rm$hap_nuc == -1)
    data <- data %>% filter(hap_nuc != -1) # mnlogit only takes data without indels in read or in haplotypes
  }
  data$nuc <- to_char_r(data$nuc)
  data$hap_nuc <- to_char_r(data$hap_nuc)
  id <- data["id"]
  for (m in (1:max_iter)) {
    cat("iter: ", m, "\n");
    par_hmm_old <- bw$par_hmm
    phi_old <- bw$par_hmm$phi
    tmp <- m_beta(res = bw$par_aux, id = id, weight_id = weight_id, data = data, formula = formula, 
                  reads_lengths = read_length, ncores = ncores, old_version = 0, weight = bw$par_aux$weight)
    # init[[m]] <- bw
    ## estimation other parameters
    bw1 <- baum_welch_iter(hmm_info = HMM, par_hmm = bw, data_info = dat_info, hap_info = hap_full, 
                          beta = tmp$par$beta, PD_LENGTH = nrow(par$beta), hash_idx = data_new$idx)
    
    if ((all(abs(exp(bw$par_hmm$phi) - exp(phi_old)) < tol) == TRUE) &&
        all(compare_par(new = bw$par_hmm, old = par_hmm_old, name = "emit", tol) == TRUE) &&
        all(compare_par(new = bw$par_hmm, old = par_hmm_old, name = "trans", tol) == TRUE))
      break;
  }
  
  ### viterbi decoding
  cat("viterbi decoding\n");
  res <- list()
  hap <- viterbi(hmm_info = HMM, dat_info = dat_info, hap_info = hap_full, par_hmm = bw$par_hmm)
  haplotypes <- matrix(to_char_r(hap), nrow = n_class)
  # idx <- duplicated(haplotypes)
  # derepliacte_h <- haplotypes[!idx, ]
  snp_location <- HMM$undecided_pos + 1
  snps <- haplotypes[, snp_location]
  snp_location <- snp_location + dat_info$ref_start
  
  res$no_reads_t <- HMM$n_t
  res$time_pos <- HMM$time_pos
  res$haplotypes <- haplotypes
  res$snps <- snps
  res$snp_location <- snp_location
  if(!is.null(res_file))
    fnlist(res, fil = res_file)
  return(res)
}


a<- altragenotype(datafile = "../../data/tpphase/WGS/simu/L_SNP/high_cov/out1.txt",
              alignment = "../../data/tpphase/WGS/simu/L_SNP/ref.fsa", max_iter = 1)

