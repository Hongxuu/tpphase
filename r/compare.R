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
res_all.heter0.005 <- get_res(parent_path = "../../../../peanut_simu/heter0.005/", covergae, individual, name = "hmm_res")
res_all.0.005 <- get_res(parent_path = "../../../../peanut_simu/homr0.005/", covergae, individual, name = "hmm_res")
res_all.0.008 <- get_res(parent_path = "../../../../peanut_simu/homr0.008/", covergae, individual, name = "hmm_res")
res_all.0.01 <- get_res(parent_path = "../../../../peanut_simu/homr0.01/", covergae, individual, name = "hmm_res")

hmm_res.heter0.005 <- get_err(individual, parent_path = "../../../../peanut_simu/heter0.005/", 
                              res_all.heter0.005, covergae, datfile_name = "hmm", is_hmm = 1, compute_mec = 0)
hmm_res.0.01 <- get_err(individual, parent_path = "../../../../peanut_simu/homr0.01/", 
                        res_all.0.01, covergae, datfile_name = "hmm", is_hmm = 1, compute_mec = 0)
hmm_res.0.005 <- get_err(individual, parent_path = "../../../../peanut_simu/homr0.005/", 
                         res_all.0.005, covergae , datfile_name = "hmm", is_hmm = 1, compute_mec = 0)
hmm_res.0.008 <- get_err(individual, parent_path = "../../../../peanut_simu/homr0.008/", 
                         res_all.0.008, covergae, datfile_name = "hmm", is_hmm = 1, compute_mec = 0)



############## original calculation
hmm_res2.0.01 <- get_err(individual, parent_path = "../../../../peanut_simu/homr0.01/", 
                        res_all.0.01, covergae, is_hmm = 1, old = 1)
hmm_res2.0.005 <- get_err(individual, parent_path = "../../../../peanut_simu/homr0.005/", 
                         res_all.0.005, covergae, is_hmm = 1, old = 1)
hmm_res2.0.008 <- get_err(individual, parent_path = "../../../../peanut_simu/homr0.008/", 
                         res_all.0.008, covergae, is_hmm = 1, old = 1)

gatk.heter0.005 <- get_res_other(parent_path = "../../../../peanut_simu/heter0.005/", covergae, individual, name = "gatk")
gatk_res.heter0.005 <- get_err(individual, parent_path = "../../../../peanut_simu/heter0.005/",
                               gatk.heter0.005, covergae, datfile_name = "gatk", is_hmm = 0, compute_mec = 0)

gatk_res2.heter0.005 <- get_err(individual = individual, parent_path = "../../../../peanut_simu/heter0.005/",
                                res_all = gatk.heter0.005, covergae = covergae, old = 1, is_hmm = 0)

gatk.0.005 <- get_res_other(parent_path = "../../../../peanut_simu/homr0.005/", covergae, individual, name = "gatk")
gatk_res.0.005 <- get_err(individual, parent_path = "../../../../peanut_simu/homr0.005/",
                          gatk.0.005, covergae, datfile_name = "gatk", is_hmm = 0, compute_mec = 0)
gatk_res2.0.005 <- get_err(individual = individual, parent_path = "../../../../peanut_simu/homr0.005/",
                                res_all = gatk.0.005, covergae = covergae, old = 1, is_hmm = 0)

gatk.0.008 <- get_res_other(parent_path = "../../../../peanut_simu/homr0.008/", covergae, individual, name = "gatk")
gatk_res.0.008 <- get_err(individual, parent_path = "../../../../peanut_simu/homr0.008/", 
                          gatk.0.008, covergae, datfile_name = "gatk", is_hmm = 0, compute_mec = 0)
gatk_res2.0.008 <- get_err(individual = individual, parent_path = "../../../../peanut_simu/homr0.008/",
                           res_all = gatk.0.008, covergae = covergae, old = 1, is_hmm = 0)

gatk.0.01 <- get_res_other(parent_path = "../../../../peanut_simu/homr0.01/", covergae, individual, name = "gatk")
gatk_res.0.01 <- get_err(individual, parent_path = "../../../../peanut_simu/homr0.01/",
                         gatk.0.01, covergae, datfile_name = "gatk", is_hmm = 0, compute_mec = 0)
gatk_res2.0.01 <- get_err(individual = individual, parent_path = "../../../../peanut_simu/homr0.01/",
                           res_all = gatk.0.01, covergae = covergae, old = 1, is_hmm = 0)


hmm_homo0.005 <- rbind(hmm_res.0.005, hmm_res.heter0.005) %>% 
  add_column(heter_rate = rep(c(0.003, 0.005), each = nrow(hmm_res.0.005))) %>% 
  filter(variable %in% c("heter_fdr", "heter_swe")) %>% 
  add_column(method = "hmm")
gatk_homo0.005 <- rbind(gatk_res.0.005, gatk_res.heter0.005) %>% 
  add_column(heter_rate = rep(c(0.003, 0.005), each = nrow(gatk_res.0.005))) %>% 
  filter(variable %in% c("heter_fdr", "heter_swe")) %>% 
  add_column(method = "gatk+hapcut2")

res_0.005 <- rbind(hmm_res.0.005, gatk_res.0.005) %>% 
  add_column(method = rep(c("hmm", "gatk+hapcut2"), each = nrow(gatk_res.0.005))) %>% 
  add_column("homeo_rate" = 0.005)
res_0.008 <- rbind(hmm_res.0.008, gatk_res.0.008) %>% 
  add_column(method = rep(c("hmm", "gatk+hapcut2"), each = nrow(gatk_res.0.008))) %>% 
  add_column("homeo_rate" = 0.008)
res_0.01 <- rbind(hmm_res.0.01, gatk_res.0.01) %>% 
  add_column(method = rep(c("hmm", "gatk+hapcut2"), each = nrow(gatk_res.0.01))) %>% 
  add_column("homeo_rate" = 0.01)

all_res <- rbind(res_0.005, res_0.008, res_0.01) %>% 
  `colnames<-`(c("coverage", "variable", "error rate", "method", "homeo_rate"))

res2_0.005 <- rbind(hmm_res2.0.005, gatk_res2.0.005) %>% 
  add_column(method = rep(c("hmm", "gatk+hapcut2"), each = nrow(gatk_res2.0.005))) %>% 
  add_column("homeo_rate" = 0.005)
res2_0.008 <- rbind(hmm_res2.0.008, gatk_res2.0.008) %>% 
  add_column(method = rep(c("hmm", "gatk+hapcut2"), each = nrow(gatk_res2.0.008))) %>% 
  add_column("homeo_rate" = 0.008)
res2_0.01 <- rbind(hmm_res2.0.01, gatk_res2.0.01) %>% 
  add_column(method = rep(c("hmm", "gatk+hapcut2"), each = nrow(gatk_res2.0.01))) %>% 
  add_column("homeo_rate" = 0.01)

all_res2 <- rbind(res2_0.005, res2_0.008, res2_0.01) %>% 
  `colnames<-`(c("coverage", "variable", "error rate", "method", "homeo_rate"))

res_heter <- all_res %>% 
  filter(variable %in% c("heter_fdr", "heter_swe")) 

res_heter2 <- all_res2 %>% 
  filter(variable %in% c("heter_fdr", "heter_swe"))
  
res_mec <- all_res %>% 
  filter(variable == "mec" & `error rate` != 0)

res_homeo <- all_res %>% 
  filter(variable %in% c("homeo_fdr", "homeo_swe", "heter_fdr", "heter_swe") & method == "hmm")


rbind(res_heter, res_mec) %>% 
  ggplot(aes(coverage, `error rate`)) +
  geom_boxplot(aes(fill = method)) + 
  facet_grid(variable~homeo_rate, scales = "free") + 
  ylab("value")

rbind(res_homeo, res_mec) %>% 
  ggplot(aes(coverage, `error rate`)) +
  geom_boxplot() + 
  facet_wrap(homeo_rate~variable, scales = "free", ncol = 5) + 
  ylab("value")

rbind(hmm_homo0.005, gatk_homo0.005) %>% 
  ggplot(aes(coverage, value)) +
  geom_boxplot(aes(fill = method)) + 
  facet_grid(variable~heter_rate, scales = "free") + 
  ylab("value")
###########MEC



# pp_1 <- get_pp(individual, res_all.0.005, covergae)