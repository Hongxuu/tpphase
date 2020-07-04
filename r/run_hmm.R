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

use_MC = 1

seed = 0

datafile = "../../data/hmm/WGS/simu/ref4/low_cov/out.txt"
alignment = "../../data/hmm/WGS/simu/ref4/ref.fsa"
call_aln(ref_nameA = "Genome_A:0-2000", ref_nameB = "Genome_B:0-2000",
         ref_fsa = alignment,
         ref_sam = "../../data/hmm/WGS/simu/ref4/ref.sam",
         alnA = "../../data/hmm/WGS/simu/ref4/high_cov/aln0A.sam",
         alnB = "../../data/hmm/WGS/simu/ref4/high_cov/aln0B.sam",
         out_file = datafile)

opts <- new.env()
assign("cut_off", 0.1, envir = opts)
assign("three_hap", 0.62, envir = opts)
assign("third_nuc", 0, envir = opts)
assign("two_hap", 0.45, envir = opts)
assign("left_range", 1.2, envir = opts)

setOpt <- function(...) {
  opts <- getOpt()
  args <- list(...)
  if(length(args)==1 && is.list(args[[1]])) {
    args <- args[[1]]
  }
  for(opnm in names(args)) {
    if(opnm %in% names(opts)) { # class() OK here, since all opts are simple objects with single classes
      if(class(getOpt(opnm)) == class(args[[opnm]])) {
        assign(opnm, args[[opnm]], envir=opts)
      } else {
        warning(paste0(opnm, " not set, value provided has different class (", class(args[[opnm]]), 
                       ") then current option value (", class(getOpt(opnm)), ")"))
      }
    } else {
      warning(opnm, " is not a valid option.")
    }
  }
}

getOpt <- function(option = NULL) {
  if(is.null(option)) option <- ls(opts)
  
  if(!all(option %in% ls(opts))) {
    warning("Tried to get an invalid option: ", option[!(option %in% ls(opts))])
    option <- option[option %in% ls(opts)]
  }
  
  ropts <- lapply(option, function(x) get(x, envir=opts))
  names(ropts) <- option
  if(length(ropts) == 1) ropts <- ropts[[1]]  # If just one option requested, return it alone
  return(ropts)
}

altragenotype <- function(datafile = NULL, alignment = NULL, ref_name = NULL, res_file = NULL, ref_delim = ".",
                          formula = mode~1|read_pos + ref_pos + qua + hap_nuc + qua:hap_nuc, max_iter = 50, 
                          seed = 0, tol = 1e-05, use_MC = 1, ncores = 2, ...)  {
  registerDoParallel(cores = ncores)  
  call <- sys.call(1)
  # Read in default opts and then replace with any that were passed in to the function
  opts <- getOpt()
  args <- list(...)
  for(opnm in names(args)) {
    if(opnm %in% names(opts)) {
      opts[[opnm]] <- args[[opnm]]
    } else {
      warning(opnm, " is not a valid option.")
    }
  }
  ## make universal reference
  cat("preparing data: \n");
  align <- read_fasta(alignment)
  if(nrow(align$reads) != 1) { #only read in one reference pair
    universal <- make_universal(alignment = align, for_hmm = 1, ref_idx = 0)
    universal <- universal %>% unlist()
  } else { ## read in many references, take the one we want
    ref_in <- strsplit(ref_name, ref_delim, fixed = TRUE) %>% unlist()
    ref_index <- ref_in[2] %>% as.integer() - 1 ## index in C
    universal <- make_universal(alignment = align, for_hmm = 1, ref_idx = ref_index)
    universal <- Filter(Negate(is.null), universal) %>% unlist()
  }
  rm(align)
  
  ## prepare data
  genotype_target = 0
  dat_info <- read_data(datafile, old_v = genotype_target)
  HMM <- hmm_info(dat_info = dat_info, uni_alignment = universal, opt = opts)
  
  ########################## baum-welch (iterate until converge)
  
  ## initialization
  ##### use linkage info to limit some unlikely happened transition
  hap_length <- dat_info$ref_length_max - dat_info$ref_start
  linkage_in <- linkage_info(dat_info = dat_info, undecided_pos = HMM$undecided_pos)
  overlap_info <- get_overlap(HMM$p_tmax, HMM$time_pos, HMM$num_states, HMM$undecided_pos, HMM$t_max, dat_info$ref_start)
  hap_full_info <- full_hap_new(HMM, linkage_in, overlap_info, hap_length, dat_info$ref_start, use_MC = use_MC)
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
  n_class = 4 
  num_cat = 4
  par <- ini_par(dat = data, n_observation = dat_info$n_observation, formula = formula, old_version = genotype_target,
                 n_class = n_class, num_cat = num_cat, ncores = ncores, weight_id = weight_id)
  
  ###indicate which transfer could happen
  trans_indicator <- trans_permit(num_states = HMM$num_states, combination = hap_full_info$combination, 
                                  loci = overlap_info$location, t_max = HMM$t_max)
  ### start initializing
  weight_id <- NULL
  cat("initialization: \n");
  data_new <- format_data2(hmm_info = HMM, d_info = dat_info, hap_info = hap_full)
  data <- data_new$df_new
  # par$eta <- eta
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
    cat("iter: ", m, "\n")
    full_llk <- bw$par_hmm_bf$full_llk
    par_hmm_old <- bw$par_hmm
    phi_old <- bw$par_hmm$phi
    start.time <- Sys.time()
    tmp <- m_beta(res = bw$par_aux, id = id, weight_id = weight_id, data = data, formula = formula, 
                  reads_lengths = read_length, ncores = ncores, old_version = genotype_target, weight = bw$par_aux$weight)
    end.time <- Sys.time()
    cat("mnlogit time: ", end.time - start.time, "\n")
    ## estimation other parameters
    start.time <- Sys.time()
    bw <- baum_welch_iter(hmm_info = HMM, par_hmm = bw, data_info = dat_info, hap_info = hap_full, 
                          beta = tmp$par$beta, PD_LENGTH = nrow(par$beta), hash_idx = data_new$idx)
    end.time <- Sys.time()
    cat("bw time: ", end.time - start.time, "\n")
    cat(bw$par_aux$eta, "\n")
    if (((all(abs(exp(bw$par_hmm$phi) - exp(phi_old)) < tol) == TRUE) &&
        all(compare_par(new = bw$par_hmm, old = par_hmm_old, name = "emit", tol) == TRUE) &&
        all(compare_par(new = bw$par_hmm, old = par_hmm_old, name = "trans", tol) == TRUE)) ||
        abs(bw$par_hmm_bf$full_llk - full_llk) < tol)
      break;
  }
  
  ### viterbi decoding
  cat("viterbi decoding\n");
  res <- list()
  hap <- viterbi(hmm_info = HMM, dat_info = dat_info, hap_info = hap_full, par_hmm = bw$par_hmm)
  haplotypes <- matrix(to_char_r(hap$hap_final), nrow = n_class)
  # ass = find_ass(selected_hap = hap$chosed_state, n_in_t = HMM$n_in_t, wic = bw$par_aux$w_ic, 
  #                t_max = HMM$t_max, n_obs = dat_info$n_observation)
  # idx <- duplicated(haplotypes)
  # derepliacte_h <- haplotypes[!idx, ]
  snp_location <- HMM$undecided_pos + 1
  ## check if the snps do contain snp
  id <- which(apply(haplotypes[, snp_location], 2, function(x) length(unique(x))) == 1)
  if(length(id) != 0) {
    snp_location <- snp_location[-id] }
  snps <- haplotypes[, snp_location]
  snp_location <- snp_location + dat_info$ref_start
  
  ## choose the weight of selected haplotypes
  res$full_llk <- bw$par_hmm_bf$full_llk
  res$haplotypes <- haplotypes
  res$snps <- snps
  res$snp_location <- snp_location
  # res$par <- bw
  if(!is.null(res_file))
    fnlist(res, fil = res_file)
  return(res)
}

altragenotype(datafile = "../../data/hmm/WGS/simu/ref2/LOW_SNP/low_cov/out.txt", 
              alignment = "../../data/hmm/WGS/simu/ref2/LOW_SNP/ref.fsa",
              res_file = "./ref2_cov50.txt") -> short_low






