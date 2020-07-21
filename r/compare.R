### read the truth
library(reshape2)
sourceCpp("./r/assess.cpp")
source("./r/assessment.R")

covergae <- c(3, 4, 8, 12, 16)
individual <- c(0:49)

res_all.0.0052 <- get_res(parent_path = "../../../../peanut_simu/homr0.005/", covergae, individual, name = "hmm2_res")

res_all.0.005 <- get_res(parent_path = "../../../../peanut_simu/homr0.005/", covergae, individual, name = "hmm_res")
res_all.0.01 <- get_res(parent_path = "../../../../peanut_simu/homr0.01/", covergae, individual, name = "hmm_res")

hmm_res.0.01 <- get_err(individual, parent_path = "../../../../peanut_simu/homr0.01/", res_all.0.01, covergae, is_hmm = 1)
hmm_res.0.005 <- get_err(individual, parent_path = "../../../../peanut_simu/homr0.005/", res_all.0.005, covergae, is_hmm = 1)
hmm_res.0.0052 <- get_err(individual, parent_path = "../../../../peanut_simu/homr0.005/", res_all.0.0052, covergae, is_hmm = 1)
# 
# roshan <- get_res_other(parent_path, covergae, individual, name = "roshan")
# roshan_res <- get_err(individual, parent_path, roshan, covergae, is_hmm = 0)
gatk.0.005 <- get_res_other(parent_path = "../../../../peanut_simu/homr0.005/", covergae, individual, name = "gatk")
gatk_res.0.005 <- get_err(individual, parent_path= "../../../../peanut_simu/homr0.005/", gatk.0.005, covergae, is_hmm = 0)
gatk.0.01 <- get_res_other(parent_path = "../../../../peanut_simu/homr0.01/", covergae, individual, name = "gatk")
gatk_res.0.01 <- get_err(individual, parent_path= "../../../../peanut_simu/homr0.01/", gatk.0.01, covergae, is_hmm = 0)

res_0.005 <- rbind(hmm_res.0.0052, gatk_res.0.005) %>% 
  add_column(method = rep(c("hmm", "gatk+hapcut2"), each = nrow(gatk_res.0.005))) %>% 
  add_column("homeo_rate" = 0.005)
res_0.01 <- rbind(hmm_res.0.01, gatk_res.0.01) %>% 
  add_column(method = rep(c("hmm", "gatk+hapcut2"), each = nrow(gatk_res.0.01))) %>% 
  add_column("homeo_rate" = 0.01)

res_heter <- rbind(res_0.005, res_0.01) %>% 
  filter(variable == c("heter_fdr", "heter_swe") & covergae != 3)

res_homeo <- rbind(res_0.005, res_0.01) %>% 
  filter(variable == c("homo_fdr", "homo_swe") & covergae != 3 & method == "hmm")

res_heter %>% 
  ggplot(aes(coverage, value)) +
  geom_boxplot(aes(fill = method)) + 
  facet_wrap(homeo_rate~variable, scales = "free")

res_homeo %>% 
  ggplot(aes(coverage, value)) +
  geom_boxplot() + 
  facet_grid(homeo_rate~variable, scales = "free")


# pp_1 <- get_pp(individual, res_all.0.005, covergae)