### read the truth
library(reshape2)
sourceCpp("./r/assess.cpp")
source("./r/assessment.R")

covergae <- c(4, 8, 12, 16)
individual <- c(0:49)

### get the out data from sam for gatk
# samfile <- list()
# outfile <- list()
# for(i in c(0.008, 0.005, 0.01)) {
#   hr = paste0("../../../../peanut_simu/" , "homr", i)
#   for(j in covergae) {
#     for(l in individual) {
#       samfile = paste0(hr, "/cov", j, "/gatk_res/sim", l, "test.sam")
#       outfile = paste0(hr, "/cov", j,  "/gatk_res/out", l, "gatk.txt")
#       read_sam(samfile, datafile = outfile)
#     }
#   }   
# }

##################

res_all.0.005 <- get_res(parent_path = "../../../../peanut_simu/homr0.005/", covergae, individual, name = "hmm_res")
res_all.0.008 <- get_res(parent_path = "../../../../peanut_simu/homr0.008/", covergae, individual, name = "hmm_res")
res_all.0.01 <- get_res(parent_path = "../../../../peanut_simu/homr0.01/", covergae, individual, name = "hmm_res")

hmm_res.0.01 <- get_err(individual, parent_path = "../../../../peanut_simu/homr0.01/", res_all.0.01, covergae, is_hmm = 1)
hmm_res.0.005 <- get_err(individual, parent_path = "../../../../peanut_simu/homr0.005/", res_all.0.005, covergae, is_hmm = 1)
hmm_res.0.008 <- get_err(individual, parent_path = "../../../../peanut_simu/homr0.008/", res_all.0.008, covergae, is_hmm = 1)

gatk.0.005 <- get_res_other(parent_path = "../../../../peanut_simu/homr0.005/", covergae, individual, name = "gatk")
gatk_res.0.005 <- get_err(individual, parent_path= "../../../../peanut_simu/homr0.005/", gatk.0.005, covergae, is_hmm = 0)
gatk.0.008 <- get_res_other(parent_path = "../../../../peanut_simu/homr0.008/", covergae, individual, name = "gatk")
gatk_res.0.008 <- get_err(individual, parent_path= "../../../../peanut_simu/homr0.008/", gatk.0.008, covergae, is_hmm = 0)
gatk.0.01 <- get_res_other(parent_path = "../../../../peanut_simu/homr0.01/", covergae, individual, name = "gatk")
gatk_res.0.01 <- get_err(individual, parent_path= "../../../../peanut_simu/homr0.01/", gatk.0.01, covergae, is_hmm = 0)

res_0.005 <- rbind(hmm_res.0.005, gatk_res.0.005) %>% 
  add_column(method = rep(c("hmm", "gatk+hapcut2"), each = nrow(gatk_res.0.005))) %>% 
  add_column("homeo_rate" = 0.005)
res_0.008 <- rbind(hmm_res.0.008, gatk_res.0.008) %>% 
  add_column(method = rep(c("hmm", "gatk+hapcut2"), each = nrow(gatk_res.0.008))) %>% 
  add_column("homeo_rate" = 0.008)
res_0.01 <- rbind(hmm_res.0.01, gatk_res.0.01) %>% 
  add_column(method = rep(c("hmm", "gatk+hapcut2"), each = nrow(gatk_res.0.01))) %>% 
  add_column("homeo_rate" = 0.01)

res_heter <- rbind(res_0.005, res_0.008, res_0.01) %>% 
  `colnames<-`(c("coverage", "variable", "error rate", "method", "homeo_rate")) %>% 
  filter(variable %in% c("heter_fdr", "heter_swe"))

res_homeo <- rbind(res_0.005, res_0.008, res_0.01) %>% 
  `colnames<-`(c("coverage", "variable", "error rate", "method", "homeo_rate")) %>% 
  filter(variable %in% c("homeo_fdr", "homeo_swe", "heter_fdr", "heter_swe") & method == "hmm")

res_heter %>% 
  ggplot(aes(coverage, `error rate`)) +
  geom_boxplot(aes(fill = method)) + 
  facet_wrap(homeo_rate~variable, scales = "free")

res_homeo %>% 
  ggplot(aes(coverage, `error rate`)) +
  geom_boxplot() + 
  facet_grid(homeo_rate~variable, scales = "free")

###########MEC



# pp_1 <- get_pp(individual, res_all.0.005, covergae)