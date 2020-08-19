### read the truth
library(reshape2)
sourceCpp("./r/assess.cpp")
source("./r/assessment.R")

covergae <- c(4, 8, 12, 16)
individual <- c(0:49)

### get the out data from sam for gatk
# samfile <- list()
# outfile <- list()
# for(i in c(0.005)) {
#   hr = paste0("../../../../peanut_simu/" , "heter", i)
#   for(j in covergae) {
#     for(l in individual) {
#       samfile = paste0(hr, "/cov", j, "/gatk_res/sim", l, "test.sam")
#       outfile = paste0(hr, "/cov", j,  "/gatk_res/out", l, "gatk.txt")
#       read_sam(samfile, datafile = outfile)
#     }
#   }
# }

##################

hmm_res.0.01.pair <- load_res(parent_path = "../../../../peanut_simu/homr0.01/", covergae, individual, datfile_name = "hmm", 
                              true_sw = 1, is_pair = 1, compute_mec = 1)
hmm_res.0.005.pair <- load_res(parent_path = "../../../../peanut_simu/homr0.005/", covergae, individual, datfile_name = "hmm", 
                               true_sw = 1, is_pair = 1, compute_mec = 1)
hmm_res.0.008.pair <- load_res(parent_path = "../../../../peanut_simu/homr0.008/", covergae, individual, datfile_name = "hmm", 
                               true_sw = 1, is_pair = 1, compute_mec = 1)

hmm_res.0.01.sig <- load_res(parent_path = "../../../../peanut_simu/homr0.01/", covergae, individual, datfile_name = "hmm", 
                              true_sw = 1, is_pair = 0, compute_mec = 1)
hmm_res.0.005.sig <- load_res(parent_path = "../../../../peanut_simu/homr0.005/", covergae, individual, datfile_name = "hmm", 
                               true_sw = 1, is_pair = 0, compute_mec = 1)
hmm_res.0.008.sig <- load_res(parent_path = "../../../../peanut_simu/homr0.008/", covergae, individual, datfile_name = "hmm", 
                               true_sw = 1, is_pair = 0, compute_mec = 1)

gatk_res.0.01.pair <- load_res(parent_path = "../../../../peanut_simu/homr0.01/", covergae, individual, datfile_name = "gatk", 
                              true_sw = 1, is_pair = 1, compute_mec = 1)
gatk_res.0.005.pair <- load_res(parent_path = "../../../../peanut_simu/homr0.005/", covergae, individual, datfile_name = "gatk", 
                               true_sw = 1, is_pair = 1, compute_mec = 1)
gatk_res.0.008.pair <- load_res(parent_path = "../../../../peanut_simu/homr0.008/", covergae, individual, datfile_name = "gatk", 
                               true_sw = 1, is_pair = 1, compute_mec = 1)


######compare alleic rate 
hmm_res.0.003.p <- load_res(parent_path = "../../../../peanut_simu/ale0.003/", covergae, individual, datfile_name = "hmm", 
                              true_sw = 1, is_pair = 1, compute_mec = 1)
hmm_res.0.005.p <- load_res(parent_path = "../../../../peanut_simu/ale0.005/", covergae, individual, datfile_name = "hmm", 
                              true_sw = 1, is_pair = 1, compute_mec = 1)
hmm_res.0.007.p <- load_res(parent_path = "../../../../peanut_simu/ale0.007/", covergae, individual, datfile_name = "hmm", 
                              true_sw = 1, is_pair = 1, compute_mec = 1)

hmm_res.0.003.sig <- load_res(parent_path = "../../../../peanut_simu/ale0.003/", covergae, individual, datfile_name = "hmm", 
                             true_sw = 1, is_pair = 0, compute_mec = 1)
hmm_res.0.005.sig <- load_res(parent_path = "../../../../peanut_simu/ale0.005/", covergae, individual, datfile_name = "hmm", 
                              true_sw = 1, is_pair = 0, compute_mec = 1)
hmm_res.0.007.sig <- load_res(parent_path = "../../../../peanut_simu/ale0.007/", covergae, individual, datfile_name = "hmm", 
                              true_sw = 1, is_pair = 0, compute_mec = 1)

hmm0.005.sig <- load_res(parent_path = "../../../../peanut_simu/hom_ale0.005/", covergae, individual, datfile_name = "hmm", 
                              true_sw = 1, is_pair = 0, compute_mec = 1)
hmm0.003.sig <- load_res(parent_path = "../../../../peanut_simu/homr0.005/", covergae, individual, datfile_name = "hmm", 
                              true_sw = 1, is_pair = 0, compute_mec = 1)

all_res.p <- rbind(hmm_res.0.003.p, hmm_res.0.005.p) %>% 
  add_column(heter_rate = rep(c(0.003, 0.005), each = nrow(hmm_res.0.003.p))) %>% 
  filter(variable %in% c("homeo_fdr", "homeo_swe", "heter_fdr", "heter_swe", "mec")) %>% 
  add_column(method = "hmm")

all_res2 <- rbind(hmm0.003.sig, hmm0.005.sig) %>% 
  add_column(heter_rate = rep(c(0.003, 0.005), each = nrow(hmm0.003.sig))) %>% 
  filter(variable %in% c("homeo_fdr", "homeo_swe", "heter_fdr", "heter_swe", "mec")) %>% 
  add_column(method = "hmm")

all_res <- rbind(hmm_res.0.03.sig, hmm_res.0.005.sig) %>% 
  add_column(heter_rate = rep(c(0.003, 0.005), each = nrow(hmm_res.0.03.sig))) %>% 
  filter(variable %in% c("homeo_fdr", "homeo_swe", "heter_fdr", "heter_swe", "mec")) %>% 
  add_column(method = "hmm")

res_0.005 <- rbind(hmm_res.0.005.sig, hmm_res.0.005.pair) %>% 
  add_column(method = rep(c("single", "pair"), c(nrow(hmm_res.0.005.sig), nrow(hmm_res.0.005.pair)))) %>% 
  add_column("homeo_rate" = 0.005)
res_0.008 <- rbind(hmm_res.0.008.sig, hmm_res.0.008.pair) %>% 
  add_column(method = rep(c("single", "pair"), c(nrow(hmm_res.0.008.sig), nrow(hmm_res.0.008.pair)))) %>% 
  add_column("homeo_rate" = 0.008)
res_0.01 <- rbind(hmm_res.0.01.sig, hmm_res.0.01.pair) %>% 
  add_column(method = rep(c("single", "pair"), c(nrow(hmm_res.0.01.sig), nrow(hmm_res.0.01.pair)))) %>% 
  add_column("homeo_rate" = 0.01)

hmm_homo0.005 <- rbind(hmm_res.0.005, hmm_res.heter0.005) %>% 
  add_column(heter_rate = rep(c(0.003, 0.005), each = nrow(hmm_res.0.005))) %>% 
  filter(variable %in% c("heter_fdr", "heter_swe")) %>% 
  add_column(method = "hmm")
gatk_homo0.005 <- rbind(gatk_res.0.005, gatk_res.heter0.005) %>% 
  add_column(heter_rate = rep(c(0.003, 0.005), each = nrow(gatk_res.0.005))) %>% 
  filter(variable %in% c("heter_fdr", "heter_swe")) %>% 
  add_column(method = "gatk+hapcut2")

res_0.005 <- rbind(hmm_res.0.005, gatk_res.0.005) %>% 
  add_column(method = rep(c("hmm", "gatk+hapcut2"), c(nrow(hmm_res.0.005), nrow(gatk_res.0.005)))) %>% 
  add_column("homeo_rate" = 0.005)
res_0.008 <- rbind(hmm_res.0.008, gatk_res.0.008) %>% 
  add_column(method = rep(c("hmm", "gatk+hapcut2"), c(nrow(hmm_res.0.005), nrow(gatk_res.0.005)))) %>% 
  add_column("homeo_rate" = 0.008)
res_0.01 <- rbind(hmm_res.0.01, gatk_res.0.01) %>% 
  add_column(method = rep(c("hmm", "gatk+hapcut2"), c(nrow(hmm_res.0.005), nrow(gatk_res.0.005)))) %>% 
  add_column("homeo_rate" = 0.01)

all_res <- rbind(res_0.005, res_0.008, res_0.01) %>% 
  `colnames<-`(c("coverage", "variable", "error rate", "method", "homeo_rate"))

res_heter <- all_res %>% 
  filter(variable %in% c("heter_fdr", "heter_swe")) 
  
res_mec <- all_res %>% 
  filter(variable == "mec" & `error rate` != 0)

res_homeo <- all_res %>% 
  filter(variable %in% c("homeo_fdr", "homeo_swe", "heter_fdr", "heter_swe") & method == "hmm")

res_mec2 <- all_res2 %>% 
  filter(variable == "mec" & value != 0)

res_homeo2 <- all_res2 %>% 
  filter(variable %in% c("homeo_fdr", "homeo_swe", "heter_fdr", "heter_swe") & method == "hmm")

res_mec.p <- all_res.p %>% 
  filter(variable == "mec" & value != 0)

res_homeo.p <- all_res.p %>% 
  filter(variable %in% c("homeo_fdr", "homeo_swe", "heter_fdr", "heter_swe") & method == "hmm")

rbind(res_heter, res_mec) %>% 
  ggplot(aes(coverage, `error rate`)) +
  geom_boxplot(aes(fill = method)) + 
  facet_grid(variable~homeo_rate, scales = "free") + 
  ylab("value")

dat <- rbind(res_homeo, res_mec)
dat$heter_rate <- as.factor(dat$heter_rate)
dat2 <- rbind(res_homeo2, res_mec2)
dat2$heter_rate <- as.factor(dat2$heter_rate)
dat.p <- rbind(res_homeo.p, res_mec.p)
dat.p$heter_rate <- as.factor(dat.p$heter_rate)

dat %>% 
  ggplot(aes(coverage, value)) +
  geom_boxplot(aes(fill = heter_rate)) + 
  facet_wrap(variable~., scales = "free", ncol = 4) + 
  ylab("value")

dat2 %>% 
  ggplot(aes(coverage, value)) +
  geom_boxplot(aes(fill = heter_rate)) + 
  facet_wrap(variable~., scales = "free", ncol = 4) + 
  ylab("value")

dat.p %>% 
  ggplot(aes(coverage, value)) +
  geom_boxplot(aes(fill = heter_rate)) + 
  facet_wrap(variable~., scales = "free", ncol = 4) + 
  ylab("value")

rbind(hmm_homo0.005, gatk_homo0.005) %>% 
  ggplot(aes(coverage, value)) +
  geom_boxplot(aes(fill = method)) + 
  facet_grid(variable~heter_rate, scales = "free") + 
  ylab("value")
###########MEC



# pp_1 <- get_pp(individual, res_all.0.005, covergae)