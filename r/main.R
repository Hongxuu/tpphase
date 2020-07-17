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
source("./r/run_hmm.R")
sourceCpp("./r/data_format.cpp")
sourceCpp("./r/mnlogit.cpp")
sourceCpp("./r/baumwelch.cpp")
sourceCpp("./r/universal_alignment.cpp")
sourceCpp("./r/initialization.cpp")
sourceCpp("./r/viterbi.cpp")

parent_folder = "../../data/hmm/WGS/"
ref_nameA = "Genome_A:0-5000"
ref_nameB = "Genome_B:0-5000"
registerDoParallel(cores = 7)

m = 1
datfile <- list()
ref_alignment <- list()
res_file <- list()
for(i in c(0.02, 0.01, 0.005)) {
  hr = paste0(parent_folder , "homr", i)
  alignment = paste0(hr, "/ref.fsa")
  ref_sam = paste0(hr, "/ref.sam")
  for(j in c(3, 4, 8, 12, 16)) {
    for(l in c(0:49)) {
      alnA = paste0(hr, "/cov", j, "/aln", l, "A.sam")
      alnB = paste0(hr, "/cov", j, "/aln", l, "B.sam")
      datafile = paste0(hr, "/cov", j, "/out", l, ".txt")
      call_aln(ref_nameA = ref_nameA,
               ref_nameB = ref_nameB,
               ref_fsa = alignment,
               ref_sam = ref_sam,
               alnA = alnA,
               alnB = alnB,
               out_file = datafile)
      datfile[[m]] = datafile
      ref_alignment[[m]] = alignment
      res_file[[m]] = paste0(hr, "/cov", j, "/hmm_res", l)
      m = m + 1
    }
  }   
}

foreach(m=1:length(datfile)) %dopar% {
  altragenotype(datafile = datfile[[m]], 
                alignment = ref_alignment[[m]], 
                res_file = res_file[[m]])
}
datafile = "../../data/hmm/WGS/simu/homr0.01/cov3/out14.txt"
alignment = "../../data/hmm/WGS/simu/homr0.01/ref.fsa"
altragenotype(datafile = datfile, 
              alignment = alignment) -> a
