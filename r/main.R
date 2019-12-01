library(tidyverse)
#library(tidyverse, lib="~/local/R_libs/")
library(ShortRead)
library(mnlogit)
library(foreach)
library(doParallel)

source("~/R/data_clean.R")
source("~/R/EM.R")

registerDoParallel(cores=2)

main <- function(dat_path, hap_path , n_class = 4, num_cat = 4, seed = 0, ncores = 2, max = 10000, filter = TRUE) {
  ### Process data
  all_dat <- data_clean(dat_path, hap_path)
  dat_short <- all_dat$dat_long 
  dat <- all_dat$lldat
  n_observation <- max(dat_short$idx)
  read_length <- as.numeric(table(dat_short$idx))/(n_class)
  
  ### Initialization
  par <- list()
  par <- ini(dat = dat, seed = seed)
  
  print("prepare done!")
  
  ### E step and estimate new haplotypes
  res <- list()
  
  full_llk <- rep(0, max)
  for (m in 1:max) {
    res <- Estep(n_observation, read_length, par, dat_short, dat, n_class, num_cat)
    
    print("Estep")
    full_llk[m] <- res$full_llk
    if(m != 1 && full_llk[m] < full_llk[m-1]) 
      break;
    ### M step: mnlogit
    par <- res$par
    par$beta <- Mstep(dat = res$dat, read_length = read_length, par = par, ncores)
    print("Mstep")
  }
  return(res)
}

res <- main(dat_path = "~/data/short.txt", hap_path = "~/data/test_hap.fa")

sink("~/res.txt")
print(res)
sink()







