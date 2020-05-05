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

##### targeted data

## read the data
samfile = "../../data/tpphase_res_consensus/WGS/test.sam"
ref_name = "target.2"
fastq_file = "./res.fastq"
datafile = "./res.txt"
alignment = "../../data/tpphase_res_consensus/WGS/new.fasta"

datafile = "../../data/tpphase_res_consensus/WGS/out.txt"
#######
formula = mode~1|read_pos + ref_pos + qua + hap_nuc + qua:hap_nuc
n_class = 4
num_cat = 4
seed = 0
max = 20
tol = 1e-06
ncores = 2
ref_delim = "."
old_version = 0
altragenotype <- function(samfile = NULL, ref_name = NULL, alignment = NULL, ref_delim = ".",
                          fastq_file = "./res.fastq", datafile = "./res.txt", output = NULL, 
                          formula = mode~1|read_pos + ref_pos + qua + hap_nuc + qua:hap_nuc, max_iter = 1,
                          n_class = 4, num_cat = 4, seed = 0, max = 20, tol = 1e-06, ncores = 2, old_version = 0)  {
  registerDoParallel(cores = ncores)  
  ## read the data
  if(is.null(samfile) == FALSE)
    sam <- read_sam(samfile, ref_name, fastq_file, datafile)
  
  ## make universial reference
  align <- read_fasta(alignment)
  ref_in <- strsplit(ref_name, ref_delim, fixed = TRUE) %>% unlist()
  if(length(align$dim) == 2) { #only read in one reference
    universial <- make_universal(alignment = align, for_hmm = 1, ref_idx = 0)
  } else { ## read in many reference, take the one we want
    ref_index <- ref_in[2] %>% as.integer() - 1 ## index in C
    universial <- make_universal(alignment = align, for_hmm = 1, ref_idx = ref_index)
  }
  universial <- Filter(Negate(is.null), universial) %>% unlist()
  rm(align)
  
  ## prepare data
  dat_info <- read_data(datafile, old_v = old_version)
  HMM <- hmm_info(dat_info = dat_info, cut_off = 0.15, uni_alignment = universial)
  
  ########################## baum-welch (iterate until converge)
  
  ## initialization
  ##### use linkage info to limit some unlikelily happened transition
  hap_length <- dat_info$ref_length_max - dat_info$ref_start
  linkage_info <- linkage_info(dat_info = dat_info, undecided_pos = HMM$undecided_pos)
  overlap_info <- get_overlap(HMM$p_tmax, HMM$time_pos, HMM$num_states, HMM$undecided_pos, HMM$t_max, dat_info$ref_start)
  hap_full_info <- full_hap_new(HMM, linkage_info, overlap_info, hap_length, dat_info$ref_start)
  # hap_full_info <- full_hap(hmm_info = HMM, linkage_info = linkage_info, hap_length = hap_length, 
  #                           hap_min_pos = dat_info$ref_start)
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
  ###some transition could not happen, can be set to null
  trans_indicator <- trans_permit(num_states = HMM$num_states, combination = hap_full_info$combination, 
                                  loci = overlap_info$location, t_max = HMM$t_max)
  # trans_indicator_new = NULL;
  # trans_indicator <- trans_permit(HMM$num_states, hap_full_info$combination, HMM$t_max, HMM$undecided_pos, 
  #              HMM$time_pos, HMM$p_tmax, dat_info$ref_start)
  # HMM$num_states <- trans_indicator$new_num_states
  # hap_full <- final_exclude(full_hap = hap_full, further_limit = trans_indicator$further_limit, 
  #                                t_max = HMM$t_max, num_states = HMM$num_states)
  # trans_indicator_new <- prepare_ini_hmm(HMM$t_max, HMM$num_states, 
  #                                        trans_indicator$trans_permits, trans_indicator$further_limit);
  
  ### start initializing
  bw <- baum_welch_init(hmm_info = HMM, data_info = dat_info, hap_info = hap_full, 
                        par = par, PD_LENGTH = nrow(par$beta), trans_indicator_new = trans_indicator)
  rm(hap_full_info)
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
  for (m in (1:max_iter)) {
    par_hmm_old <- bw$par_hmm
    phi_old <- bw$par_hmm$phi
    
    tmp <- m_beta(res = bw$par_aux, id = id, weight_id = weight_id, data = data, formula = formula, 
                  reads_lengths = read_length, ncores = ncores, old_version = 0, weight = bw$par_aux$weight)
    
    ## estimation other parameters
    bw <- baum_welch_iter(hmm_info = HMM, par_hmm = bw, data_info = dat_info, hap_info = hap_full, 
                          beta = tmp$par$beta, PD_LENGTH = nrow(par$beta))
    
    if (abs(bw$par_hmm$phi - phi_old) < tol & 
       compare_par(new = bw$par_hmm, old = par_hmm_old, name = "emit") &
       compare_par(new = bw$par_hmm, old = par_hmm_old, name = "trans"))
      break;
  }
  
  ### viterbi decoding
  res <- list()
  hap <- viterbi(hmm_info = HMM, dat_info = dat_info, hap_info = hap_full, par_hmm = bw$par_hmm)
  haplotypes <- matrix(to_char_r(hap), nrow = n_class)
  idx <- duplicated(haplotypes)
  derepliacte_h <- haplotypes[!idx, ]
  snp_location <- HMM$undecided_pos + 1
  snps <- derepliacte_h[, snp_location]
  
  res$no_reads_t <- HMM$n_t
  res$time_pos <- HMM$time_pos
  res$haplotypes <- derepliacte_h
  res$snps <- snps
  res$snp_location <- snp_location
  fnlist(res, fil = "./test.res")
}

# if(trans_indicator_new.isNotNull()){
#   trans_permits = trans_indicator["trans_permits"];
#   further_limit = trans_indicator["further_limit"];
#   trans_indicator_new = prepare_ini_hmm(t_max, num_states, trans_permits, further_limit);
# }
sourceCpp("./r/extra.cpp")
comb_info1 = find_combination(HMM$undecided_pos, HMM$pos_possibility, HMM$p_tmax[1], HMM$time_pos[1], dat_info$ref_start)
comb_info2 = find_combination(HMM$undecided_pos, HMM$pos_possibility, HMM$p_tmax[2], HMM$time_pos[2], dat_info$ref_start)
comb_info3 = find_combination(HMM$undecided_pos, HMM$pos_possibility, HMM$p_tmax[3], HMM$time_pos[3], dat_info$ref_start)
comb_info4 = find_combination(HMM$undecided_pos, HMM$pos_possibility, HMM$p_tmax[4], HMM$time_pos[4], dat_info$ref_start)
comb_info5 = find_combination(HMM$undecided_pos, HMM$pos_possibility, HMM$p_tmax[5], HMM$time_pos[5], dat_info$ref_start)
comb_info6 = find_combination(HMM$undecided_pos, HMM$pos_possibility, HMM$p_tmax[6], HMM$time_pos[6], dat_info$ref_start)
comb_info7 = find_combination(HMM$undecided_pos, HMM$pos_possibility, HMM$p_tmax[7], HMM$time_pos[7], dat_info$ref_start)
comb_info8 = find_combination(HMM$undecided_pos, HMM$pos_possibility, HMM$p_tmax[8], HMM$time_pos[8], dat_info$ref_start)


t0 <- limit_comb_t0(comb_info1$combination, HMM$hidden_states, comb_info1$location, linkage_info = linkage_info, 
                    comb_info1$num, 0, HMM$num_states[1])
t1 <- t0
t2 <- unique_overlap(overlapped = overlap_info$overlapped[[3]], exclude_last = t0$exclude, overlap_comb = comb_info1$combination, 
                     overlap_loci = overlap_info$location[[1]], overlap_new_states = t0$num_states, overlap_num_states = HMM$num_states[1])

t3 <- new_combination(HMM, location = overlap_info$location[[4]], overlapped = overlap_info$overlapped[[4]], exclude_last = t0$exclude, 
                     overlap_comb = comb_info1$combination, 
                overlap_loci = overlap_info$location[[1]], linkage_info = linkage_info,
                overlap_new_states = t0$num_states, overlap_num_states = HMM$num_states[1])

t4 <- unique_overlap(overlapped = overlap_info$overlapped[[5]], exclude_last = t0$exclude, overlap_comb = comb_info1$combination, 
                     overlap_loci = overlap_info$location[[1]], overlap_new_states = t0$num_states, overlap_num_states = HMM$num_states[1])
exclude <- rep(0, nrow(t3))
t5 <- new_combination(HMM, location = overlap_info$location[[6]], overlapped = overlap_info$overlapped[[6]], exclude_last = exclude, 
                      overlap_comb = t3, overlap_loci = overlap_info$location[[4]], linkage_info = linkage_info,
                      overlap_new_states = nrow(t3), overlap_num_states = nrow(t3))
t6 <- unique_overlap(overlapped = overlap_info$overlapped[[7]], exclude_last = t0$exclude, overlap_comb = comb_info1$combination, 
                     overlap_loci = overlap_info$location[[1]], overlap_new_states = t0$num_states, overlap_num_states = HMM$num_states[1])
t11 <- hap_full_info$combination[[11]]
exclude <- rep(0, nrow(t11))
t7 <- hap_full_info$combination[[8]]
t12 <- new_combination(HMM, location = overlap_info$location[[12]], overlapped = overlap_info$overlapped[[12]], exclude_last = exclude, 
                      overlap_comb = t11, overlap_loci = overlap_info$location[[11]], linkage_info = linkage_info,
                      overlap_new_states = nrow(t11), overlap_num_states = nrow(t11))
exclude <- rep(0, nrow(t7))
first_comb = unique_overlap(overlap_info$overlapped[[11]], exclude, t7, overlap_info$location[[11]], 
                            nrow(t7),  nrow(t7))
t11 <- new_combination(HMM, location = overlap_info$location[[11]], overlapped = overlap_info$overlapped[[11]], exclude_last = exclude, 
                       overlap_comb = t7, overlap_loci = overlap_info$location[[8]], linkage_info = linkage_info,
                       overlap_new_states = nrow(t7), overlap_num_states = nrow(t7))
tmp <- new_combination(HMM, location = overlap_info$location[[11]], overlapped = overlap_info$overlapped[[11]], exclude_last = rep(0, nrow(hap_full_info$combination[[8]])), 
                       overlap_comb = hap_full_info$combination[[8]], overlap_loci = overlap_info$location[[8]], linkage_info = linkage_info,
                       overlap_new_states = nrow(hap_full_info$combination[[8]]), 
                       overlap_num_states =  nrow(hap_full_info$combination[[8]]))

unique_overlap(overlapped = overlap_info$overlapped[[11]], exclude_last = rep(0, nrow(hap_full_info$combination[[8]])), 
               overlap_comb = hap_full_info$combination[[8]], 
               overlap_loci = overlap_info$location[[8]], overlap_new_states = nrow(hap_full_info$combination[[8]]), 
               overlap_num_states =  nrow(hap_full_info$combination[[8]]))



