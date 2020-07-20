### read the truth
library(reshape2)
sourceCpp("./r/assess.cpp")

parent_path <- "../../../../peanut_simu/homr0.005/"
covergae <- c(3, 4, 8, 12, 16)
individual <- c(0:49)

res_all.0.005 <- get_res(parent_path = "../../../../peanut_simu/homr0.005/", covergae, individual, name = "hmm_res")

res_all.0.005 <- get_res(parent_path = "../../../../peanut_simu/homr0.005/", covergae, individual)
res_all.0.01 <- get_res(parent_path = "../../../../peanut_simu/homr0.01/", covergae, individual)

hmm_res <- get_err(individual, parent_path = "../../../../peanut_simu/homr0.01/", res_all.0.01, covergae, is_hmm = 1)
hmm_res <- get_err(individual, parent_path, res_all.0.005, covergae, is_hmm = 1)


roshan <- get_res_other(parent_path, covergae, individual, name = "roshan")
roshan_res <- get_err(individual, parent_path, roshan, covergae, is_hmm = 0)

gatk <- get_res_other(parent_path, covergae, individual, name = "gatk")
gatk_res <- get_err(individual, parent_path, gatk, covergae, is_hmm = 0)

res_0.005 <- rbind(hmm_res, roshan_res, gatk_res) %>% add_column(alg = rep(c("hmm", "roshan", "gatk"), each = nrow(gatk_res)))
res_0.01 <- rbind(hmm_res, roshan_res, gatk_res) %>% add_column(alg = rep(c("hmm", "roshan", "gatk"), each = nrow(gatk_res)))
res_0.005 %>% 
  ggplot(aes(coverage, value)) +
  geom_boxplot(aes(fill = alg)) + 
  facet_wrap(~variable, scales = "free")

 

# pp_1 <- get_pp(individual, res_all.0.005, covergae)