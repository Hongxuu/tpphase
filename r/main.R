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

parent_folder = "../data/peanut/peanut_simu/"
ref_nameA = "Genome_A:0-5000"
ref_nameB = "Genome_B:0-5000"
registerDoParallel(cores = 16)

m = 1
datfile <- list()
ref_alignment <- list()
res_file <- list()
uni_file <- list()
for(i in c(0.007)) {
  hr = paste0(parent_folder , "ale", i)
  alignment = paste0(hr, "/ref.fsa")
  ref_sam = paste0(hr, "/ref.sam")
  for(j in c(16)) {
    for(l in c(33)) {
      alnA = paste0(hr, "/pair/cov", j, "/aln", l, "A.sam")
      alnB = paste0(hr, "/pair/cov", j, "/aln", l, "B.sam")
      datafile = paste0(hr, "/pair/cov", j, "/hmm_res", "/out", l, ".txt")
      uni_geno = paste0(hr, "/pair/cov", j, "/hmm_res", "/uni", l, ".fa")
      call_aln(ref_nameA = ref_nameA,
               ref_nameB = ref_nameB,
               ref_fsa = alignment,
               ref_sam = ref_sam,
               alnA = alnA,
               alnB = alnB,
               out_file = datafile,
               uni_geno_file = uni_geno)
      datfile[[m]] = datafile
      uni_file[[m]] = uni_geno
      ref_alignment[[m]] = alignment
      res_file[[m]] = paste0(hr, "/pair/cov", j, "/hmm_res", l)
      m = m + 1
    }
  }   
}

foreach(m=1:length(datfile)) %dopar% {
  altragenotype(datafile = datfile[[m]], 
                alignment = uni_file[[m]], 
                res_file = res_file[[m]])
}
read_sam(samfile = "../../../../peanut_simu/heter0.005/cov4/gatk_res/sim0test.sam")


datafile = "../../../../peanut_simu/homr0.005/pair/cov16/hmm_res/out0hmm.txt"
alignment = "../../../../peanut_simu/homr0.005/pair/cov16/hmm_res/uni0.fa"

alignment = "../alleotetraploid/data/roshan/real/target_homeo.fasta"
ref_sam = "../alleotetraploid/data/roshan/real/combine.sam" 
alnA = "../alleotetraploid/data/roshan/real/Tifguard_subset_A_WGS.sam"
alnB = "../alleotetraploid/data/roshan/real/Tifguard_subset_B_WGS.sam"
datafile = "../../out1.txt"
uni_geno = "../../uni1.fa"
res <- altraphase(datafile = "../../../../out11.txt", 
                  alignment = "../../../../uni11.fa", max_iter = 2)
alnA = "../alleotetraploid/data/hmm/WGS/real/high_cov_genotypes/SRR8361736-A_subset.sam"
alnB = "../alleotetraploid/data/hmm/WGS/real/high_cov_genotypes/SRR8361736-B_subset.sam"
#######################real data

call_aln(ref_nameA = "aradu.V14167.gnm2.chr04:119340965-119341964",
         ref_nameB = "araip.K30076.gnm2.B04:141199721-141200721",
         ref_fsa = alignment,
         ref_sam = ref_sam,
         alnA = alnA,
         alnB = alnB,
         out_file = datafile,
         uni_geno_file = uni_geno)



